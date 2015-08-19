// Copyright 2015 Cloudera Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef IMPALA_FLE_ENCODING_H
#define IMPALA_FLE_ENCODING_H

#include <immintrin.h>
#include <math.h>
#include <bitset>

#include <boost/bind.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/function.hpp>

#include "common/compiler-util.h"
#include "util/bit-stream-utils.inline.h"
#include "util/bit-util.h"

using namespace std;
using namespace boost;

namespace impala {

class FleDecoder {
 public:
  FleDecoder(uint8_t* buffer, int buffer_len, int bit_width)
    : bit_width_(bit_width) {
    DCHECK_GE(bit_width_, 0);
    DCHECK_LE(bit_width_, 64);
    buffer_ = reinterpret_cast<uint64_t*>(buffer);
    buffer_end_ = buffer_;
    buffer_end_guard_ = reinterpret_cast<uint64_t*>(buffer + buffer_len);
    count_ = 64;
    memset(&current_value_, 0, 64 * 4);
    pcurrent_value8_ = reinterpret_cast<uint8_t*>(current_value_);
    pcurrent_value16_ = reinterpret_cast<uint16_t*>(current_value_);

    switch(bit_width) {
      case 1:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_1), this);
        break;
      case 2:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_2), this);
        break;
      case 3:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_3), this);
        break;
      case 4:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_4), this);
        break;
      case 5:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_5), this);
        break;
      case 6:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_6), this);
        break;
      case 7:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_7), this);
        break;
      case 8:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_8), this);
        break;
      case 9:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_9), this);
        break;
      case 10:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_10), this);
        break;
      case 11:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_11), this);
        break;
      case 12:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_12), this);
        break;
      case 13:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_13), this);
        break;
      case 14:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_14), this);
        break;
      case 15:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_15), this);
        break;
      case 16:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_16), this);
        break;
      case 17:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_17), this);
        break;
      case 18:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_18), this);
        break;
      case 19:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_19), this);
        break;
      case 20:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_20), this);
        break;
      case 21:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_21), this);
        break;
      case 22:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_22), this);
        break;
      case 23:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_23), this);
        break;
      case 24:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_24), this);
        break;
      case 25:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_25), this);
        break;
      case 26:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_26), this);
        break;
      case 27:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_27), this);
        break;
      case 28:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_28), this);
        break;
      case 29:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_29), this);
        break;
      case 30:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_30), this);
        break;
      case 31:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_31), this);
        break;
      case 32:
        unpack_function_ = bind(mem_fn(&FleDecoder::Unpack_32), this);
        break;
    }

    clr_mask = _mm256_set1_epi64x(0x0102040810204080);
    for (int i = 0; i < 8; ++i) {
      bitv[i] = _mm256_set1_epi8(0x01 << i);
    }
    shf0_mask = _mm256_setr_epi64x(0x0707070707070707, 0x0606060606060606,
        0x0505050505050505, 0x0404040404040404);
    shf1_mask = _mm256_setr_epi64x(0x0303030303030303, 0x0202020202020202,
        0x0101010101010101, 0x0000000000000000);
    shf2_mask = _mm256_setr_epi64x(0x0f0f0f0f0f0f0f0f, 0x0e0e0e0e0e0e0e0e,
        0x0d0d0d0d0d0d0d0d, 0x0c0c0c0c0c0c0c0c);
    shf3_mask = _mm256_setr_epi64x(0x0b0b0b0b0b0b0b0b, 0x0a0a0a0a0a0a0a0a,
        0x0909090909090909, 0x0808080808080808);
    shf_ret_mask = _mm256_setr_epi64x(0x8003800280018000, 0x8007800680058004,
        0x800b800a80098008, 0x800f800e800d800c);
    shf_ret1_mask = _mm256_setr_epi64x(0x0380028001800080, 0x0780068005800480,
        0x0b800a8009800880, 0x0f800e800d800c80);

  }

  FleDecoder() {}

  template<typename T>
  bool Get(T* val);

  template<typename T>
  bool Get(T* val, int skip_rows);

  bool Skip(int skip_rows);

  void Unpack_1();
  void Unpack_2();
  void Unpack_3();
  void Unpack_4();
  void Unpack_5();
  void Unpack_6();
  void Unpack_7();
  void Unpack_8();
  void Unpack_9();
  void Unpack_10();
  void Unpack_11();
  void Unpack_12();
  void Unpack_13();
  void Unpack_14();
  void Unpack_15();
  void Unpack_16();
  void Unpack_17();
  void Unpack_18();
  void Unpack_19();
  void Unpack_20();
  void Unpack_21();
  void Unpack_22();
  void Unpack_23();
  void Unpack_24();
  void Unpack_25();
  void Unpack_26();
  void Unpack_27();
  void Unpack_28();
  void Unpack_29();
  void Unpack_30();
  void Unpack_31();
  void Unpack_32();

  void Unpack_16_32();

  void Eq(int64_t num_rows, dynamic_bitset<>& skip_bitset, uint64_t value);
  void Lt(int64_t num_rows, dynamic_bitset<>& skip_bitset, uint64_t value);
  void Le(int64_t num_rows, dynamic_bitset<>& skip_bitset, uint64_t value);
  void Gt(int64_t num_rows, dynamic_bitset<>& skip_bitset, uint64_t value);
  void Ge(int64_t num_rows, dynamic_bitset<>& skip_bitset, uint64_t value);
  void In(int64_t num_rows, dynamic_bitset<>& skip_bitset, vector<uint64_t>& values);
  void Clear() {}

  int bit_width() { return bit_width_; }
 private:
  typedef boost::function<void()> UnpackFunction;
  UnpackFunction unpack_function_;
  __m256i clr_mask;
  __m256i bitv[8];
  __m256i shf0_mask;
  __m256i shf1_mask;
  __m256i shf2_mask;
  __m256i shf3_mask;
  __m256i shf_ret_mask;
  __m256i shf_ret1_mask;
  uint64_t* buffer_;
  uint64_t* buffer_end_;
  uint64_t* buffer_end_guard_;
  int bit_width_;
  uint32_t current_value_[64];
  uint8_t* pcurrent_value8_;
  uint16_t* pcurrent_value16_;
  int count_;
};

class FleEncoder {
 public:
  FleEncoder(uint8_t* buffer, int buffer_len, int bit_width)
    : bit_width_(bit_width) {
    DCHECK_GE(bit_width_, 1);
    DCHECK_LE(bit_width_, 64);
    DCHECK_GE(buffer_len % 8, 0);
    buffer_ = reinterpret_cast<uint64_t*>(buffer);
    buffer_end_ = buffer_;
    buffer_end_guard_ = buffer_ + buffer_len / 8 - bit_width;
    count_ = 0;
    pcurrent_value8_ = reinterpret_cast<uint8_t*>(current_value_);
    pcurrent_value16_ = reinterpret_cast<uint16_t*>(current_value_);

    fa[1] = &FleEncoder::Pack_1;
    fa[2] = &FleEncoder::Pack_2;
    fa[3] = &FleEncoder::Pack_3;
    fa[4] = &FleEncoder::Pack_4;
    fa[5] = &FleEncoder::Pack_5;
    fa[6] = &FleEncoder::Pack_6;
    fa[7] = &FleEncoder::Pack_7;
    fa[8] = &FleEncoder::Pack_8;
    fa[9] = &FleEncoder::Pack_9;
    fa[10] = &FleEncoder::Pack_10;
    fa[11] = &FleEncoder::Pack_11;
    fa[12] = &FleEncoder::Pack_12;
    fa[13] = &FleEncoder::Pack_13;
    fa[14] = &FleEncoder::Pack_14;
    fa[15] = &FleEncoder::Pack_15;
    fa[16] = &FleEncoder::Pack_16;

    fa[17] = &FleEncoder::Pack_16_32;
    fa[18] = &FleEncoder::Pack_16_32;
    fa[19] = &FleEncoder::Pack_16_32;
    fa[20] = &FleEncoder::Pack_16_32;
    fa[21] = &FleEncoder::Pack_16_32;
    fa[22] = &FleEncoder::Pack_16_32;
    fa[23] = &FleEncoder::Pack_16_32;
    fa[24] = &FleEncoder::Pack_16_32;
    fa[25] = &FleEncoder::Pack_16_32;
    fa[26] = &FleEncoder::Pack_16_32;
    fa[27] = &FleEncoder::Pack_16_32;
    fa[28] = &FleEncoder::Pack_16_32;
    fa[29] = &FleEncoder::Pack_16_32;
    fa[30] = &FleEncoder::Pack_16_32;
    fa[31] = &FleEncoder::Pack_16_32;
    fa[32] = &FleEncoder::Pack_16_32;

    for (int i = 0; i < 8; ++i) {
      bitv[i] = _mm256_set1_epi8(0x01 << i);
    }

    clr_mask = _mm256_set1_epi64x(0x00ff00ff00ff00ff);
    clr32_mask = _mm256_set1_epi64x(0x0000ffff0000ffff);
  }

  bool Put(uint64_t value);

  int Flush();

  void Clear() {
    count_ = 0;
    buffer_end_ = buffer_;
  }

  uint8_t* buffer() { return reinterpret_cast<uint8_t*>(buffer_); }
  int len() { return (buffer_end_ - buffer_) * 8;};

  void Pack_1();
  void Pack_2();
  void Pack_3();
  void Pack_4();
  void Pack_5();
  void Pack_6();
  void Pack_7();
  void Pack_8();
  void Pack_9();
  void Pack_10();
  void Pack_11();
  void Pack_12();
  void Pack_13();
  void Pack_14();
  void Pack_15();
  void Pack_16();
  void Pack_16_32();

 private:
  typedef void(impala::FleEncoder::*PPACK)();
  PPACK fa[33];
  __m256i bitv[8];
  __m256i clr_mask;
  __m256i clr32_mask;
  uint64_t* buffer_;
  uint64_t* buffer_end_;
  uint64_t* buffer_end_guard_;
  uint32_t current_value_[64];
  uint8_t* pcurrent_value8_;
  uint16_t* pcurrent_value16_;
  int bit_width_;
  int count_;
};

template<typename T>
inline bool FleDecoder::Get(T* val, int skip_rows) {
  if (UNLIKELY(count_ == 64)) {
    count_ = skip_rows & (64 - 1);
    int skip_64 = (skip_rows & ~(64 - 1)) >> 6;
    buffer_end_ += bit_width_ * skip_64;
    if (buffer_end_ >= buffer_end_guard_) return false;
    unpack_function_();
    buffer_end_ += bit_width_;
  } else {
    skip_rows += count_;
    count_ = skip_rows & (64 - 1);
    if (skip_rows >= 64) {
      int skip_64 = ((skip_rows & ~(64 - 1)) >> 6) - 1;
      buffer_end_ += bit_width_ * skip_64;
      if (buffer_end_ >= buffer_end_guard_) return false;
      unpack_function_();
      buffer_end_ += bit_width_;
    }
  }

  if (bit_width_ <= 8) {
    *val = pcurrent_value8_[count_];
  } else if (bit_width_ <= 16) {
    *val = pcurrent_value16_[count_];
  } else {
    *val = current_value_[count_];
  }

  if (*val == 0) {
    *val = 0;
  }
  ++count_;

  return true;
}

inline bool FleDecoder::Skip(int skip_rows) {
  if (UNLIKELY(count_ == 64)) {
    count_ = skip_rows & (64 - 1);
    int skip_64 = (skip_rows & ~(64 - 1)) >> 6;
    buffer_end_ += bit_width_ * skip_64;
    if (buffer_end_ >= buffer_end_guard_) return false;
    unpack_function_();
    buffer_end_ += bit_width_;
  } else {
    skip_rows += count_;
    count_ = skip_rows & (64 - 1);
    if (skip_rows >= 64) {
      int skip_64 = ((skip_rows & ~(64 - 1)) >> 6) - 1;
      buffer_end_ += bit_width_ * skip_64;
      if (buffer_end_ >= buffer_end_guard_) return false;
      unpack_function_();
      buffer_end_ += bit_width_;
    }
  }

  return true;
}

template<typename T>
inline bool FleDecoder::Get(T* val) {
  if (UNLIKELY(count_ == 64)) {
    if (buffer_end_ >= buffer_end_guard_) return false;
    count_ = 0;
//    for (int i = 0; i < bit_width_; ++i) {
//      uint64_t value = buffer_end_[i];
//      uint64_t set1 = 0x01ULL << i;
//      uint64_t index = _tzcnt_u64(value);
//      while (index != 64) {
//        current_value_[63 - index] |= set1;
//        value = _blsr_u64(value);
//        index = _tzcnt_u64(value);
//      }

//      uint64_t set1 = 0x01ULL << i;
//      for (int l = 0; l < 64; ++l) {
//        if (buffer_end_[i] & (0x01ULL << l)) current_value_[63 - l] |= set1;
//      }
//    }

//    memset(current_value_, 0, 64 * 8);
//    for (int i = 0; i < bit_width_; ++i) {
//      for (int l = 0; l < 64; ++l) {
//        current_value_[63 - l] |= ((buffer_end_[i] >> l) & 0x01ULL) << i;
//      }
//    }

//    memset(current_value_, 0, 64 * 4);
    unpack_function_();
//    (this->*fa[bit_width_])();
/*
    if (bit_width_ == 1) {
      __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
      __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
      buf = _mm256_inserti128_si256(buf, buf_hi, 1);
      __m256i bit0_mask, bit, ret;
      __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
      __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
          0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);

      for (int k = 0; k < 8; ++k) {
        bit0_mask = _mm256_set1_epi32(7 - k);
        bit = _mm256_shuffle_epi8(buf, bit0_mask);
        bit = _mm256_srlv_epi32(bit, rshift);
        ret = _mm256_and_si256(bit, clr_mask);
        _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
      }
    } else if (bit_width_ & 0x01) {
      int index = bit_width_ / 2;
      __m256i buf[16];
      for (int i = 0; i < index + 1; ++i) {
        buf[i] = _mm256_loadu_si256((__m256i const*)(buffer_end_ + i * 2));
        __m128i buf_hi = _mm256_extracti128_si256(buf[i], 0);
        buf[i] = _mm256_inserti128_si256(buf[i], buf_hi, 1);
      }

      __m256i bit0_mask, bit, bit1_mask, ret;
      __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
      __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
          0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);

      for (int k = 0; k < 8; ++k) {
        bit0_mask = _mm256_set1_epi32(7 - k);
        bit1_mask = _mm256_set1_epi32(15 - k);
        bit = _mm256_shuffle_epi8(buf[0], bit0_mask);
        bit = _mm256_srlv_epi32(bit, rshift);
        ret = _mm256_and_si256(bit, clr_mask);
        bit = _mm256_shuffle_epi8(buf[0], bit1_mask);
        bit = _mm256_srlv_epi32(bit, rshift);
        bit = _mm256_and_si256(bit, clr_mask);
        bit = _mm256_slli_epi32(bit, 1);
        ret = _mm256_or_si256(ret, bit);
        for (int l = 1; l < index; ++l) {
          bit = _mm256_shuffle_epi8(buf[l], bit0_mask);
          bit = _mm256_srlv_epi32(bit, rshift);
          bit = _mm256_and_si256(bit, clr_mask);
          bit = _mm256_slli_epi32(bit, 2 * l);
          ret = _mm256_or_si256(ret, bit);
          bit = _mm256_shuffle_epi8(buf[l], bit1_mask);
          bit = _mm256_srlv_epi32(bit, rshift);
          bit = _mm256_and_si256(bit, clr_mask);
          bit = _mm256_slli_epi32(bit, 2 * l + 1);
          ret = _mm256_or_si256(ret, bit);
        }
        bit = _mm256_shuffle_epi8(buf[index], bit0_mask);
        bit = _mm256_srlv_epi32(bit, rshift);
        bit = _mm256_and_si256(bit, clr_mask);
        bit = _mm256_slli_epi32(bit, bit_width_ - 1);
        ret = _mm256_or_si256(ret, bit);

        _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
      }
    } else {
      __m256i buf[16];
      for (int i = 0; i < bit_width_ / 2; ++i) {
        buf[i] = _mm256_loadu_si256((__m256i const*)(buffer_end_ + i * 2));
        __m128i buf_hi = _mm256_extracti128_si256(buf[i], 0);
        buf[i] = _mm256_inserti128_si256(buf[i], buf_hi, 1);
      }

      __m256i bit0_mask, bit, bit1_mask, ret;
      __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
      __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
          0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);

      for (int k = 0; k < 8; ++k) {
        bit0_mask = _mm256_set1_epi32(7 - k);
        bit1_mask = _mm256_set1_epi32(15 - k);
        bit = _mm256_shuffle_epi8(buf[0], bit0_mask);
        bit = _mm256_srlv_epi32(bit, rshift);
        ret = _mm256_and_si256(bit, clr_mask);
        bit = _mm256_shuffle_epi8(buf[0], bit1_mask);
        bit = _mm256_srlv_epi32(bit, rshift);
        bit = _mm256_and_si256(bit, clr_mask);
        bit = _mm256_slli_epi32(bit, 1);
        ret = _mm256_or_si256(ret, bit);
        for (int l = 1; l < bit_width_ / 2; ++l) {
          bit = _mm256_shuffle_epi8(buf[l], bit0_mask);
          bit = _mm256_srlv_epi32(bit, rshift);
          bit = _mm256_and_si256(bit, clr_mask);
          bit = _mm256_slli_epi32(bit, 2 * l);
          ret = _mm256_or_si256(ret, bit);
          bit = _mm256_shuffle_epi8(buf[l], bit1_mask);
          bit = _mm256_srlv_epi32(bit, rshift);
          bit = _mm256_and_si256(bit, clr_mask);
          bit = _mm256_slli_epi32(bit, 2 * l + 1);
          ret = _mm256_or_si256(ret, bit);
        }

        _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
      }
    }
*/
    buffer_end_ += bit_width_;
  }

//  *val = current_value_[count_];
//  current_value_[count_] = 0;
//  ++count_;

//  int shift = 63 - count_;
//  for (int i = 0; i < bit_width_; ++i) {
//    *val |= ((buffer_end_[i] >> shift) & 0x01ULL) << i;
//  }
//  ++count_;

  if (bit_width_ <= 8) {
    *val = pcurrent_value8_[count_];
  } else if (bit_width_ <= 16) {
    *val = pcurrent_value16_[count_];
  } else {
    *val = current_value_[count_];
  }

  if (*val == 0) {
    *val = 0;
  }

//  *val = current_value_[count_];
  ++count_;

  return true;
}

inline void FleDecoder::Unpack_1() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i bit, ret_mask, ret;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_), ret);
  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_ + 32), ret);
}

inline void FleDecoder::Unpack_2() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i bit, ret_mask, tmp_ret, ret;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_), ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_ + 32), ret);
}

inline void FleDecoder::Unpack_3() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m256i buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i bit, ret_mask, tmp_ret, ret;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_), ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_ + 32), ret);
}


inline void FleDecoder::Unpack_4() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m256i buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i bit, ret_mask, tmp_ret, ret;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_), ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_ + 32), ret);
}

inline void FleDecoder::Unpack_5() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m256i buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);
  __m256i bit, ret_mask, tmp_ret, ret;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_), ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_ + 32), ret);
}


inline void FleDecoder::Unpack_6() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m256i buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);
  __m256i bit, ret_mask, tmp_ret, ret;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_), ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_ + 32), ret);
}

inline void FleDecoder::Unpack_7() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m256i buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  __m256i buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);
  __m256i bit, ret_mask, tmp_ret, ret;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_), ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_ + 32), ret);
}

inline void FleDecoder::Unpack_8() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m256i buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  __m256i buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);
  __m256i bit, ret_mask, tmp_ret, ret;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);

  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_), ret);


  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);

  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value8_ + 32), ret);
}

inline void FleDecoder::Unpack_9() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m256i buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  __m256i buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);
  __m256i buf4 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf4 = _mm256_permute4x64_epi64(buf4, 0x44);
  __m256i bit, ret_mask, tmp_ret, tmp_ret1, ret, ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);

  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf4, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret1 = _mm256_and_si256(bitv[0], ret_mask);

  tmp_ret = _mm256_unpacklo_epi8(ret, ret1);
  tmp_ret1 = _mm256_unpackhi_epi8(ret, ret1);
  ret = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret1 = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_), ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 16), ret1);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);

  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf4, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret1 = _mm256_and_si256(bitv[0], ret_mask);

  tmp_ret = _mm256_unpacklo_epi8(ret, ret1);
  tmp_ret1 = _mm256_unpackhi_epi8(ret, ret1);
  ret = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret1 = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 32), ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 48), ret1);
}

inline void FleDecoder::Unpack_10() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m256i buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  __m256i buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);
  __m256i buf4 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf4 = _mm256_permute4x64_epi64(buf4, 0x44);
  __m256i bit, ret_mask, tmp_ret, tmp_ret1, ret, ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf4, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret1 = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf4, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret, ret1);
  tmp_ret1 = _mm256_unpackhi_epi8(ret, ret1);
  ret = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret1 = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_), ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 16), ret1);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf4, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret1 = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf4, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret, ret1);
  tmp_ret1 = _mm256_unpackhi_epi8(ret, ret1);
  ret = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret1 = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 32), ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 48), ret1);
}

inline void FleDecoder::Unpack_11() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m256i buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  __m256i buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);
  __m256i buf4 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  __m256i buf5 = _mm256_permute4x64_epi64(buf4, 0xee);
  buf4 = _mm256_permute4x64_epi64(buf4, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1, ret, ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf4, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret1 = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf4, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret, ret1);
  tmp_ret1 = _mm256_unpackhi_epi8(ret, ret1);
  ret = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret1 = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_), ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 16), ret1);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf4, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret1 = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf4, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret, ret1);
  tmp_ret1 = _mm256_unpackhi_epi8(ret, ret1);
  ret = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret1 = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 32), ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 48), ret1);
}

inline void FleDecoder::Unpack_12() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m256i buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  __m256i buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);
  __m256i buf4 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  __m256i buf5 = _mm256_permute4x64_epi64(buf4, 0xee);
  buf4 = _mm256_permute4x64_epi64(buf4, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1, ret, ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf4, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret1 = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf4, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret, ret1);
  tmp_ret1 = _mm256_unpackhi_epi8(ret, ret1);
  ret = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret1 = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_), ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 16), ret1);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf4, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret1 = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf4, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret, ret1);
  tmp_ret1 = _mm256_unpackhi_epi8(ret, ret1);
  ret = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret1 = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 32), ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 48), ret1);
}

inline void FleDecoder::Unpack_13() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m256i buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  __m256i buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);
  __m256i buf4 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  __m256i buf5 = _mm256_permute4x64_epi64(buf4, 0xee);
  buf4 = _mm256_permute4x64_epi64(buf4, 0x44);
  __m256i buf6 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf6 = _mm256_permute4x64_epi64(buf6, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1, ret, ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf4, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret1 = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf4, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf6, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret, ret1);
  tmp_ret1 = _mm256_unpackhi_epi8(ret, ret1);
  ret = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret1 = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_), ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 16), ret1);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf4, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret1 = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf4, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf6, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret, ret1);
  tmp_ret1 = _mm256_unpackhi_epi8(ret, ret1);
  ret = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret1 = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 32), ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 48), ret1);
}



inline void FleDecoder::Unpack_14() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m256i buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  __m256i buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);
  __m256i buf4 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  __m256i buf5 = _mm256_permute4x64_epi64(buf4, 0xee);
  buf4 = _mm256_permute4x64_epi64(buf4, 0x44);
  __m256i buf6 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf6 = _mm256_permute4x64_epi64(buf6, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1, ret, ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf4, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret1 = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf4, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf6, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf6, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret, ret1);
  tmp_ret1 = _mm256_unpackhi_epi8(ret, ret1);
  ret = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret1 = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_), ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 16), ret1);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf4, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret1 = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf4, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf6, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf6, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret, ret1);
  tmp_ret1 = _mm256_unpackhi_epi8(ret, ret1);
  ret = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret1 = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 32), ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 48), ret1);
}


inline void FleDecoder::Unpack_15() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m256i buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  __m256i buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);
  __m256i buf4 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  __m256i buf5 = _mm256_permute4x64_epi64(buf4, 0xee);
  buf4 = _mm256_permute4x64_epi64(buf4, 0x44);
  __m256i buf6 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  __m256i buf7 = _mm256_permute4x64_epi64(buf6, 0xee);
  buf6 = _mm256_permute4x64_epi64(buf6, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1, ret, ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf4, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret1 = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf4, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf6, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf6, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf7, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret, ret1);
  tmp_ret1 = _mm256_unpackhi_epi8(ret, ret1);
  ret = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret1 = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_), ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 16), ret1);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret = _mm256_or_si256(ret, tmp_ret);

  bit = _mm256_shuffle_epi8(buf4, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret1 = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf4, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf5, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf6, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf6, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);
  bit = _mm256_shuffle_epi8(buf7, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret1 = _mm256_or_si256(ret1, tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret, ret1);
  tmp_ret1 = _mm256_unpackhi_epi8(ret, ret1);
  ret = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret1 = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 32), ret);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 48), ret1);
}


inline void FleDecoder::Unpack_16() {
  __m256i ret[4], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_), ret[0]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 16), ret[2]);

  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 32), ret[1]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(pcurrent_value16_ + 48), ret[3]);
}

inline void FleDecoder::Unpack_17() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf = _mm256_permute4x64_epi64(buf, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);

  ret[6] = _mm256_setzero_si256();
  ret[7] = _mm256_setzero_si256();

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_18() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf = _mm256_permute4x64_epi64(buf, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  ret[6] = _mm256_setzero_si256();
  ret[7] = _mm256_setzero_si256();

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_19() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  ret[6] = _mm256_setzero_si256();
  ret[7] = _mm256_setzero_si256();

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_20() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  ret[6] = _mm256_setzero_si256();
  ret[7] = _mm256_setzero_si256();

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_21() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 20));
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  ret[6] = _mm256_setzero_si256();
  ret[7] = _mm256_setzero_si256();

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_22() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 20));
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  ret[6] = _mm256_setzero_si256();
  ret[7] = _mm256_setzero_si256();

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_23() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 20));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  ret[6] = _mm256_setzero_si256();
  ret[7] = _mm256_setzero_si256();

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_24() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 20));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  ret[6] = _mm256_setzero_si256();
  ret[7] = _mm256_setzero_si256();

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_25() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 20));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 24));
  buf = _mm256_permute4x64_epi64(buf, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[6] = _mm256_and_si256(bitv[0], ret_mask);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[7] = _mm256_and_si256(bitv[0], ret_mask);

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_26() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 20));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 24));
  buf = _mm256_permute4x64_epi64(buf, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[6] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[7] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_27() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 20));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 24));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[6] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[7] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_28() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 20));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 24));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[6] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[7] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_29() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 20));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 24));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 28));
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[6] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[7] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_30() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 20));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 24));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 28));
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[6] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[7] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_31() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 20));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 24));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 28));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[6] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[7] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_32() {
  __m256i ret[8], buf, buf1, buf2, buf3;
  buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  __m256i bit, ret_mask, tmp_ret, tmp_ret1;

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[0] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[0] = _mm256_or_si256(ret[0], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[1] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[1] = _mm256_or_si256(ret[1], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 12));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[2] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[2] = _mm256_or_si256(ret[2], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[3] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[3] = _mm256_or_si256(ret[3], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 16));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 20));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[4] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[4] = _mm256_or_si256(ret[4], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[5] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[5] = _mm256_or_si256(ret[5], tmp_ret);

  buf = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 24));
  buf1 = _mm256_permute4x64_epi64(buf, 0xee);
  buf = _mm256_permute4x64_epi64(buf, 0x44);
  buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 28));
  buf3 = _mm256_permute4x64_epi64(buf2, 0xee);
  buf2 = _mm256_permute4x64_epi64(buf2, 0x44);

  bit = _mm256_shuffle_epi8(buf, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[6] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf0_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf2_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[6] = _mm256_or_si256(ret[6], tmp_ret);

  bit = _mm256_shuffle_epi8(buf, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  ret[7] = _mm256_and_si256(bitv[0], ret_mask);
  bit = _mm256_shuffle_epi8(buf, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[1], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[2], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf1, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[3], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[4], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf2, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[5], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf1_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[6], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);
  bit = _mm256_shuffle_epi8(buf3, shf3_mask);
  bit = _mm256_and_si256(bit, clr_mask);
  ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
  tmp_ret = _mm256_and_si256(bitv[7], ret_mask);
  ret[7] = _mm256_or_si256(ret[7], tmp_ret);

  tmp_ret = _mm256_unpacklo_epi8(ret[0], ret[2]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[0], ret[2]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[1], ret[3]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[1], ret[3]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[4], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[4], ret[6]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  tmp_ret = _mm256_unpacklo_epi8(ret[5], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi8(ret[5], ret[7]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);

  tmp_ret = _mm256_unpacklo_epi16(ret[0], ret[4]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[0], ret[4]);
  ret[0] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_), ret[0]);
  ret[4] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 8), ret[4]);
  tmp_ret = _mm256_unpacklo_epi16(ret[2], ret[6]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[2], ret[6]);
  ret[2] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 16), ret[2]);
  ret[6] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 24), ret[6]);
  tmp_ret = _mm256_unpacklo_epi16(ret[1], ret[5]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[1], ret[5]);
  ret[1] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 32), ret[1]);
  ret[5] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 40), ret[5]);
  tmp_ret = _mm256_unpacklo_epi16(ret[3], ret[7]);
  tmp_ret1 = _mm256_unpackhi_epi16(ret[3], ret[7]);
  ret[3] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x20);
  _mm256_storeu_si256((__m256i*)(current_value_ + 48), ret[3]);
  ret[7] = _mm256_permute2x128_si256(tmp_ret, tmp_ret1, 0x31);
  _mm256_storeu_si256((__m256i*)(current_value_ + 56), ret[7]);
}

inline void FleDecoder::Unpack_16_32() {
  memset(current_value_, 0, 64 * 4);

  if (bit_width_ & 0x01) {
    int index = bit_width_ / 2;
    __m256i buf[16];
    for (int i = 0; i < index + 1; ++i) {
      buf[i] = _mm256_loadu_si256((__m256i const*)(buffer_end_ + i * 2));
      __m128i buf_hi = _mm256_extracti128_si256(buf[i], 0);
      buf[i] = _mm256_inserti128_si256(buf[i], buf_hi, 1);
    }

    __m256i bit0_mask, bit, bit1_mask, ret;
    __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
        0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);

    for (int k = 0; k < 8; ++k) {
      bit0_mask = _mm256_set1_epi32(7 - k);
      bit1_mask = _mm256_set1_epi32(15 - k);
      bit = _mm256_shuffle_epi8(buf[0], bit0_mask);
      bit = _mm256_srlv_epi32(bit, rshift);
      ret = _mm256_and_si256(bit, clr_mask);
      bit = _mm256_shuffle_epi8(buf[0], bit1_mask);
      bit = _mm256_srlv_epi32(bit, rshift);
      bit = _mm256_and_si256(bit, clr_mask);
      bit = _mm256_slli_epi32(bit, 1);
      ret = _mm256_or_si256(ret, bit);
      for (int l = 1; l < index; ++l) {
        bit = _mm256_shuffle_epi8(buf[l], bit0_mask);
        bit = _mm256_srlv_epi32(bit, rshift);
        bit = _mm256_and_si256(bit, clr_mask);
        bit = _mm256_slli_epi32(bit, 2 * l);
        ret = _mm256_or_si256(ret, bit);
        bit = _mm256_shuffle_epi8(buf[l], bit1_mask);
        bit = _mm256_srlv_epi32(bit, rshift);
        bit = _mm256_and_si256(bit, clr_mask);
        bit = _mm256_slli_epi32(bit, 2 * l + 1);
        ret = _mm256_or_si256(ret, bit);
      }
      bit = _mm256_shuffle_epi8(buf[index], bit0_mask);
      bit = _mm256_srlv_epi32(bit, rshift);
      bit = _mm256_and_si256(bit, clr_mask);
      bit = _mm256_slli_epi32(bit, bit_width_ - 1);
      ret = _mm256_or_si256(ret, bit);

      _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
    }
  } else {
    __m256i buf[16];
    for (int i = 0; i < bit_width_ / 2; ++i) {
      buf[i] = _mm256_loadu_si256((__m256i const*)(buffer_end_ + i * 2));
      __m128i buf_hi = _mm256_extracti128_si256(buf[i], 0);
      buf[i] = _mm256_inserti128_si256(buf[i], buf_hi, 1);
    }

    __m256i bit0_mask, bit, bit1_mask, ret;
    __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
        0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);

    for (int k = 0; k < 8; ++k) {
      bit0_mask = _mm256_set1_epi32(7 - k);
      bit1_mask = _mm256_set1_epi32(15 - k);
      bit = _mm256_shuffle_epi8(buf[0], bit0_mask);
      bit = _mm256_srlv_epi32(bit, rshift);
      ret = _mm256_and_si256(bit, clr_mask);
      bit = _mm256_shuffle_epi8(buf[0], bit1_mask);
      bit = _mm256_srlv_epi32(bit, rshift);
      bit = _mm256_and_si256(bit, clr_mask);
      bit = _mm256_slli_epi32(bit, 1);
      ret = _mm256_or_si256(ret, bit);
      for (int l = 1; l < bit_width_ / 2; ++l) {
        bit = _mm256_shuffle_epi8(buf[l], bit0_mask);
        bit = _mm256_srlv_epi32(bit, rshift);
        bit = _mm256_and_si256(bit, clr_mask);
        bit = _mm256_slli_epi32(bit, 2 * l);
        ret = _mm256_or_si256(ret, bit);
        bit = _mm256_shuffle_epi8(buf[l], bit1_mask);
        bit = _mm256_srlv_epi32(bit, rshift);
        bit = _mm256_and_si256(bit, clr_mask);
        bit = _mm256_slli_epi32(bit, 2 * l + 1);
        ret = _mm256_or_si256(ret, bit);
      }

      _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
    }
  }
}


/*
inline void FleDecoder::Unpack_1() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
  buf = _mm256_inserti128_si256(buf, buf_hi, 1);
  __m256i bit0_mask, bit, ret;

  for (int k = 0; k < 8; ++k) {
    bit0_mask = _mm256_set1_epi32(7 - k);
    bit = _mm256_shuffle_epi8(buf, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    ret = _mm256_and_si256(bit, clr_mask);
    _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
  }
}

inline void FleDecoder::Unpack_2() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
  buf = _mm256_inserti128_si256(buf, buf_hi, 1);
  __m256i bit0_mask, bit1_mask, bit, ret;

  for (int k = 0; k < 8; ++k) {
    bit0_mask = _mm256_set1_epi32(7 - k);
    bit1_mask = _mm256_set1_epi32(15 - k);
    bit = _mm256_shuffle_epi8(buf, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    ret = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_shuffle_epi8(buf, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 1);
    ret = _mm256_or_si256(ret, bit);
    _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
  }
}

inline void FleDecoder::Unpack_3() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
  buf = _mm256_inserti128_si256(buf, buf_hi, 1);
  __m256i buf1 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 2));
  buf_hi = _mm256_extracti128_si256(buf1, 0);
  buf1 = _mm256_inserti128_si256(buf1, buf_hi, 1);
  __m256i bit0_mask, bit1_mask, bit, ret;

  for (int k = 0; k < 8; ++k) {
    bit0_mask = _mm256_set1_epi32(7 - k);
    bit1_mask = _mm256_set1_epi32(15 - k);
    bit = _mm256_shuffle_epi8(buf, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    ret = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_shuffle_epi8(buf, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 1);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 2);
    ret = _mm256_or_si256(ret, bit);
    _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
  }
}

inline void FleDecoder::Unpack_4() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
  buf = _mm256_inserti128_si256(buf, buf_hi, 1);
  __m256i buf1 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 2));
  buf_hi = _mm256_extracti128_si256(buf1, 0);
  buf1 = _mm256_inserti128_si256(buf1, buf_hi, 1);
  __m256i bit0_mask, bit1_mask, bit, ret;

  for (int k = 0; k < 8; ++k) {
    bit0_mask = _mm256_set1_epi32(7 - k);
    bit1_mask = _mm256_set1_epi32(15 - k);
    bit = _mm256_shuffle_epi8(buf, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    ret = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_shuffle_epi8(buf, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 1);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 2);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 3);
    ret = _mm256_or_si256(ret, bit);
    _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
  }
}

inline void FleDecoder::Unpack_5() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
  buf = _mm256_inserti128_si256(buf, buf_hi, 1);
  __m256i buf1 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 2));
  buf_hi = _mm256_extracti128_si256(buf1, 0);
  buf1 = _mm256_inserti128_si256(buf1, buf_hi, 1);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf_hi = _mm256_extracti128_si256(buf2, 0);
  buf2 = _mm256_inserti128_si256(buf2, buf_hi, 1);
  __m256i bit0_mask, bit1_mask, bit, ret;

  for (int k = 0; k < 8; ++k) {
    bit0_mask = _mm256_set1_epi32(7 - k);
    bit1_mask = _mm256_set1_epi32(15 - k);
    bit = _mm256_shuffle_epi8(buf, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    ret = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_shuffle_epi8(buf, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 1);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 2);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 3);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf2, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 4);
    ret = _mm256_or_si256(ret, bit);
    _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
  }
}

inline void FleDecoder::Unpack_6() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
  buf = _mm256_inserti128_si256(buf, buf_hi, 1);
  __m256i buf1 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 2));
  buf_hi = _mm256_extracti128_si256(buf1, 0);
  buf1 = _mm256_inserti128_si256(buf1, buf_hi, 1);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf_hi = _mm256_extracti128_si256(buf2, 0);
  buf2 = _mm256_inserti128_si256(buf2, buf_hi, 1);
  __m256i bit0_mask, bit1_mask, bit, ret;

  for (int k = 0; k < 8; ++k) {
    bit0_mask = _mm256_set1_epi32(7 - k);
    bit1_mask = _mm256_set1_epi32(15 - k);
    bit = _mm256_shuffle_epi8(buf, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    ret = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_shuffle_epi8(buf, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 1);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 2);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 3);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf2, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 4);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf2, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 5);
    ret = _mm256_or_si256(ret, bit);
    _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
  }
}

inline void FleDecoder::Unpack_7() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
  buf = _mm256_inserti128_si256(buf, buf_hi, 1);
  __m256i buf1 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 2));
  buf_hi = _mm256_extracti128_si256(buf1, 0);
  buf1 = _mm256_inserti128_si256(buf1, buf_hi, 1);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf_hi = _mm256_extracti128_si256(buf2, 0);
  buf2 = _mm256_inserti128_si256(buf2, buf_hi, 1);
  __m256i buf3 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 6));
  buf_hi = _mm256_extracti128_si256(buf3, 0);
  buf3 = _mm256_inserti128_si256(buf3, buf_hi, 1);
  __m256i bit0_mask, bit1_mask, bit, ret;

  for (int k = 0; k < 8; ++k) {
    bit0_mask = _mm256_set1_epi32(7 - k);
    bit1_mask = _mm256_set1_epi32(15 - k);
    bit = _mm256_shuffle_epi8(buf, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    ret = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_shuffle_epi8(buf, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 1);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 2);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 3);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf2, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 4);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf2, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 5);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf3, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 6);
    ret = _mm256_or_si256(ret, bit);
    _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
  }

}

inline void FleDecoder::Unpack_8() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
  buf = _mm256_inserti128_si256(buf, buf_hi, 1);
  __m256i buf1 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 2));
  buf_hi = _mm256_extracti128_si256(buf1, 0);
  buf1 = _mm256_inserti128_si256(buf1, buf_hi, 1);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf_hi = _mm256_extracti128_si256(buf2, 0);
  buf2 = _mm256_inserti128_si256(buf2, buf_hi, 1);
  __m256i buf3 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 6));
  buf_hi = _mm256_extracti128_si256(buf3, 0);
  buf3 = _mm256_inserti128_si256(buf3, buf_hi, 1);
  __m256i bit0_mask, bit1_mask, bit, ret;

  for (int k = 0; k < 8; ++k) {
    bit0_mask = _mm256_set1_epi32(7 - k);
    bit1_mask = _mm256_set1_epi32(15 - k);
    bit = _mm256_shuffle_epi8(buf, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    ret = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_shuffle_epi8(buf, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 1);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 2);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 3);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf2, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 4);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf2, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 5);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf3, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 6);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf3, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 7);
    ret = _mm256_or_si256(ret, bit);

    _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
  }
}

inline void FleDecoder::Unpack_9() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
  buf = _mm256_inserti128_si256(buf, buf_hi, 1);
  __m256i buf1 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 2));
  buf_hi = _mm256_extracti128_si256(buf1, 0);
  buf1 = _mm256_inserti128_si256(buf1, buf_hi, 1);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf_hi = _mm256_extracti128_si256(buf2, 0);
  buf2 = _mm256_inserti128_si256(buf2, buf_hi, 1);
  __m256i buf3 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 6));
  buf_hi = _mm256_extracti128_si256(buf3, 0);
  buf3 = _mm256_inserti128_si256(buf3, buf_hi, 1);
  __m256i buf4 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf_hi = _mm256_extracti128_si256(buf4, 0);
  buf4 = _mm256_inserti128_si256(buf4, buf_hi, 1);
  __m256i bit0_mask, bit1_mask, bit, ret;

  for (int k = 0; k < 8; ++k) {
    bit0_mask = _mm256_set1_epi32(7 - k);
    bit1_mask = _mm256_set1_epi32(15 - k);
    bit = _mm256_shuffle_epi8(buf, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    ret = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_shuffle_epi8(buf, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 1);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 2);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 3);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf2, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 4);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf2, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 5);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf3, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 6);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf3, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 7);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf4, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 8);
    ret = _mm256_or_si256(ret, bit);
    _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
  }

}

inline void FleDecoder::Unpack_10() {
  __m256i buf = _mm256_loadu_si256((__m256i const*)buffer_end_);
  __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
  buf = _mm256_inserti128_si256(buf, buf_hi, 1);
  __m256i buf1 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 2));
  buf_hi = _mm256_extracti128_si256(buf1, 0);
  buf1 = _mm256_inserti128_si256(buf1, buf_hi, 1);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 4));
  buf_hi = _mm256_extracti128_si256(buf2, 0);
  buf2 = _mm256_inserti128_si256(buf2, buf_hi, 1);
  __m256i buf3 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 6));
  buf_hi = _mm256_extracti128_si256(buf3, 0);
  buf3 = _mm256_inserti128_si256(buf3, buf_hi, 1);
  __m256i buf4 = _mm256_loadu_si256((__m256i const*)(buffer_end_ + 8));
  buf_hi = _mm256_extracti128_si256(buf4, 0);
  buf4 = _mm256_inserti128_si256(buf4, buf_hi, 1);
  __m256i bit0_mask, bit1_mask, bit, ret;

  for (int k = 0; k < 8; ++k) {
    bit0_mask = _mm256_set1_epi32(7 - k);
    bit1_mask = _mm256_set1_epi32(15 - k);
    bit = _mm256_shuffle_epi8(buf, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    ret = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_shuffle_epi8(buf, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 1);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 2);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf1, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 3);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf2, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 4);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf2, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 5);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf3, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 6);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf3, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 7);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf4, bit0_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 8);
    ret = _mm256_or_si256(ret, bit);
    bit = _mm256_shuffle_epi8(buf4, bit1_mask);
    bit = _mm256_srlv_epi32(bit, rshift);
    bit = _mm256_and_si256(bit, clr_mask);
    bit = _mm256_slli_epi32(bit, 9);
    ret = _mm256_or_si256(ret, bit);
    _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
  }

}


inline void FleDecoder::Unpack_11_32() {
  if (bit_width_ & 0x01) {
    int index = bit_width_ / 2;
    __m256i buf[16];
    for (int i = 0; i < index + 1; ++i) {
      buf[i] = _mm256_loadu_si256((__m256i const*)(buffer_end_ + i * 2));
      __m128i buf_hi = _mm256_extracti128_si256(buf[i], 0);
      buf[i] = _mm256_inserti128_si256(buf[i], buf_hi, 1);
    }

    __m256i bit0_mask, bit, bit1_mask, ret;
    __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
        0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);

    for (int k = 0; k < 8; ++k) {
      bit0_mask = _mm256_set1_epi32(7 - k);
      bit1_mask = _mm256_set1_epi32(15 - k);
      bit = _mm256_shuffle_epi8(buf[0], bit0_mask);
      bit = _mm256_srlv_epi32(bit, rshift);
      ret = _mm256_and_si256(bit, clr_mask);
      bit = _mm256_shuffle_epi8(buf[0], bit1_mask);
      bit = _mm256_srlv_epi32(bit, rshift);
      bit = _mm256_and_si256(bit, clr_mask);
      bit = _mm256_slli_epi32(bit, 1);
      ret = _mm256_or_si256(ret, bit);
      for (int l = 1; l < index; ++l) {
        bit = _mm256_shuffle_epi8(buf[l], bit0_mask);
        bit = _mm256_srlv_epi32(bit, rshift);
        bit = _mm256_and_si256(bit, clr_mask);
        bit = _mm256_slli_epi32(bit, 2 * l);
        ret = _mm256_or_si256(ret, bit);
        bit = _mm256_shuffle_epi8(buf[l], bit1_mask);
        bit = _mm256_srlv_epi32(bit, rshift);
        bit = _mm256_and_si256(bit, clr_mask);
        bit = _mm256_slli_epi32(bit, 2 * l + 1);
        ret = _mm256_or_si256(ret, bit);
      }
      bit = _mm256_shuffle_epi8(buf[index], bit0_mask);
      bit = _mm256_srlv_epi32(bit, rshift);
      bit = _mm256_and_si256(bit, clr_mask);
      bit = _mm256_slli_epi32(bit, bit_width_ - 1);
      ret = _mm256_or_si256(ret, bit);

      _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
    }
  } else {
    __m256i buf[16];
    for (int i = 0; i < bit_width_ / 2; ++i) {
      buf[i] = _mm256_loadu_si256((__m256i const*)(buffer_end_ + i * 2));
      __m128i buf_hi = _mm256_extracti128_si256(buf[i], 0);
      buf[i] = _mm256_inserti128_si256(buf[i], buf_hi, 1);
    }

    __m256i bit0_mask, bit, bit1_mask, ret;
    __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
        0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);

    for (int k = 0; k < 8; ++k) {
      bit0_mask = _mm256_set1_epi32(7 - k);
      bit1_mask = _mm256_set1_epi32(15 - k);
      bit = _mm256_shuffle_epi8(buf[0], bit0_mask);
      bit = _mm256_srlv_epi32(bit, rshift);
      ret = _mm256_and_si256(bit, clr_mask);
      bit = _mm256_shuffle_epi8(buf[0], bit1_mask);
      bit = _mm256_srlv_epi32(bit, rshift);
      bit = _mm256_and_si256(bit, clr_mask);
      bit = _mm256_slli_epi32(bit, 1);
      ret = _mm256_or_si256(ret, bit);
      for (int l = 1; l < bit_width_ / 2; ++l) {
        bit = _mm256_shuffle_epi8(buf[l], bit0_mask);
        bit = _mm256_srlv_epi32(bit, rshift);
        bit = _mm256_and_si256(bit, clr_mask);
        bit = _mm256_slli_epi32(bit, 2 * l);
        ret = _mm256_or_si256(ret, bit);
        bit = _mm256_shuffle_epi8(buf[l], bit1_mask);
        bit = _mm256_srlv_epi32(bit, rshift);
        bit = _mm256_and_si256(bit, clr_mask);
        bit = _mm256_slli_epi32(bit, 2 * l + 1);
        ret = _mm256_or_si256(ret, bit);
      }

      _mm256_storeu_si256((__m256i*)(current_value_ + 8 * k), ret);
    }
  }
}

*/

inline void FleDecoder::Eq(int64_t num_rows, dynamic_bitset<>& skip_bitset,
    uint64_t value) {
  vector<uint64_t> C;
  for (int i = 0; i < bit_width_; ++i) {
    if (value & 0x01 << i) {
      C.push_back(~0x0);
    } else {
      C.push_back(0x0);
    }
  }

  for (int i = count_; i != 64 && num_rows > 0; ++i, --num_rows) {
    if (bit_width_ <= 8) {
      skip_bitset.push_back(pcurrent_value8_[i] == value);
    } else if (bit_width_ <= 16) {
      skip_bitset.push_back(pcurrent_value16_[i] == value);
    } else {
      skip_bitset.push_back(current_value_[i] == value);
    }
  }

  uint64_t j = 0;
  while (num_rows > 0) {
    uint64_t Meq = ~0x0;
    for (int i = 0; i < bit_width_; ++i, ++j) {
      Meq = Meq & ~(buffer_end_[j] ^ C[i]);
    }
    Meq = ((Meq >> 1) & 0x5555555555555555) | ((Meq & 0x5555555555555555) << 1);
    // swap consecutiMeqe pairs
    Meq = ((Meq >> 2) & 0x3333333333333333) | ((Meq & 0x3333333333333333) << 2);
    // swap nibbles ...
    Meq = ((Meq >> 4) & 0x0F0F0F0F0F0F0F0F) | ((Meq & 0x0F0F0F0F0F0F0F0F) << 4);
    // swap bytes
    Meq = ((Meq >> 8) & 0x00FF00FF00FF00FF) | ((Meq & 0x00FF00FF00FF00FF) << 8);
    // swap 2-byte long pairs
    Meq = ((Meq >> 16) & 0x0000FFFF0000FFFF) | ((Meq & 0x0000FFFF0000FFFF) << 16);
    Meq = ( Meq >> 32                      ) | ( Meq                       << 32);

    if (num_rows < 64) {
      std::bitset<64> tmp_bitset(Meq);
      for (int i = 0; i < num_rows; ++i) {
        skip_bitset.push_back(tmp_bitset[i]);
      }
      break;
    }
    skip_bitset.append(Meq);
    num_rows -= 64;
  }
}

inline void FleDecoder::Lt(int64_t num_rows, dynamic_bitset<>& skip_bitset,
    uint64_t value) {
  vector<uint64_t> C;
  for (int i = 0; i < bit_width_; ++i) {
    if (value & 0x01 << i) {
      C.push_back(~0x0);
    } else {
      C.push_back(0x0);
    }
  }

  for (int i = count_; i != 64 && num_rows > 0; ++i, --num_rows) {
    if (bit_width_ <= 8) {
      skip_bitset.push_back(pcurrent_value8_[i] < value);
    } else if (bit_width_ <= 16) {
      skip_bitset.push_back(pcurrent_value16_[i] < value);
    } else {
      skip_bitset.push_back(current_value_[i] < value);
    }
  }

  uint64_t start_idx = 0;
  while (num_rows > 0) {
    uint64_t Mlt = 0x0;
    uint64_t Meq = ~0x0;
    int i = bit_width_ - 1;
    int j = start_idx + bit_width_ - 1;
    for (; i >= 0; --i, --j) {
      Mlt = Mlt | (Meq & C[i] & ~buffer_end_[j]);
      Meq = Meq & ~(buffer_end_[j] ^ C[i]);
    }
    start_idx += bit_width_;

    Mlt = ((Mlt >> 1) & 0x5555555555555555) | ((Mlt & 0x5555555555555555) << 1);
    // swap consecutiMlte pairs
    Mlt = ((Mlt >> 2) & 0x3333333333333333) | ((Mlt & 0x3333333333333333) << 2);
    // swap nibbles ...
    Mlt = ((Mlt >> 4) & 0x0F0F0F0F0F0F0F0F) | ((Mlt & 0x0F0F0F0F0F0F0F0F) << 4);
    // swap bytes
    Mlt = ((Mlt >> 8) & 0x00FF00FF00FF00FF) | ((Mlt & 0x00FF00FF00FF00FF) << 8);
    // swap 2-byte long pairs
    Mlt = ((Mlt >> 16) & 0x0000FFFF0000FFFF) | ((Mlt & 0x0000FFFF0000FFFF) << 16);
    Mlt = ( Mlt >> 32                      ) | ( Mlt                       << 32);

    if (num_rows < 64) {
      std::bitset<64> tmp_bitset(Mlt);
      for (int i = 0; i < num_rows; ++i) {
        skip_bitset.push_back(tmp_bitset[i]);
      }
      break;
    }
    skip_bitset.append(Mlt);
    num_rows -= 64;
  }
}

inline void FleDecoder::Le(int64_t num_rows, dynamic_bitset<>& skip_bitset,
    uint64_t value) {
  vector<uint64_t> C;
  for (int i = 0; i < bit_width_; ++i) {
    if (value & 0x01 << i) {
      C.push_back(~0x0);
    } else {
      C.push_back(0x0);
    }
  }

  for (int i = count_; i != 64 && num_rows > 0; ++i, --num_rows) {
    if (bit_width_ <= 8) {
      skip_bitset.push_back(pcurrent_value8_[i] <= value);
    } else if (bit_width_ <= 16) {
      skip_bitset.push_back(pcurrent_value16_[i] <= value);
    } else {
      skip_bitset.push_back(current_value_[i] <= value);
    }
  }

  uint64_t start_idx = 0;
  while (num_rows > 0) {
    uint64_t Mlt = 0x0;
    uint64_t Meq = ~0x0;
    int i = bit_width_ - 1;
    int j = start_idx + bit_width_ - 1;
    for (; i >= 0; --i, --j) {
      Mlt = Mlt | (Meq & C[i] & ~buffer_end_[j]);
      Meq = Meq & ~(buffer_end_[j] ^ C[i]);
    }
    start_idx += bit_width_;
    uint64_t Mle = Mlt|Meq;
    Mle = ((Mle >> 1) & 0x5555555555555555) | ((Mle & 0x5555555555555555) << 1);
    // swap consecutiMlee pairs
    Mle = ((Mle >> 2) & 0x3333333333333333) | ((Mle & 0x3333333333333333) << 2);
    // swap nibbles ...
    Mle = ((Mle >> 4) & 0x0F0F0F0F0F0F0F0F) | ((Mle & 0x0F0F0F0F0F0F0F0F) << 4);
    // swap bytes
    Mle = ((Mle >> 8) & 0x00FF00FF00FF00FF) | ((Mle & 0x00FF00FF00FF00FF) << 8);
    // swap 2-byte long pairs
    Mle = ((Mle >> 16) & 0x0000FFFF0000FFFF) | ((Mle & 0x0000FFFF0000FFFF) << 16);
    Mle = ( Mle >> 32                      ) | ( Mle                       << 32);

    if (num_rows < 64) {
      std::bitset<64> tmp_bitset(Mle);
      for (int i = 0; i < num_rows; ++i) {
        skip_bitset.push_back(tmp_bitset[i]);
      }
      break;
    }
    skip_bitset.append(Mle);
    num_rows -= 64;
  }
}

inline void FleDecoder::Gt(int64_t num_rows, dynamic_bitset<>& skip_bitset,
    uint64_t value) {
  vector<uint64_t> C;
  for (int i = 0; i < bit_width_; ++i) {
    if (value & 0x01 << i) {
      C.push_back(~0x0);
    } else {
      C.push_back(0x0);
    }
  }

  for (int i = count_; i != 64 && num_rows > 0; ++i, --num_rows) {
    if (bit_width_ <= 8) {
      skip_bitset.push_back(pcurrent_value8_[i] > value);
    } else if (bit_width_ <= 16) {
      skip_bitset.push_back(pcurrent_value16_[i] > value);
    } else {
      skip_bitset.push_back(current_value_[i] > value);
    }
  }

  uint64_t start_idx = 0;
  while (num_rows > 0) {
    uint64_t Mgt = 0x0;
    uint64_t Meq = ~0x0;
    int i = bit_width_ - 1;
    int j = start_idx + bit_width_ - 1;
    for (; i >= 0; --i, --j) {
      Mgt = Mgt | (Meq & ~C[i] & buffer_end_[j]);
      Meq = Meq & ~(buffer_end_[j] ^ C[i]);
    }
    start_idx += bit_width_;

    Mgt = ((Mgt >> 1) & 0x5555555555555555) | ((Mgt & 0x5555555555555555) << 1);
    // swap consecutiMgte pairs
    Mgt = ((Mgt >> 2) & 0x3333333333333333) | ((Mgt & 0x3333333333333333) << 2);
    // swap nibbles ...
    Mgt = ((Mgt >> 4) & 0x0F0F0F0F0F0F0F0F) | ((Mgt & 0x0F0F0F0F0F0F0F0F) << 4);
    // swap bytes
    Mgt = ((Mgt >> 8) & 0x00FF00FF00FF00FF) | ((Mgt & 0x00FF00FF00FF00FF) << 8);
    // swap 2-byte long pairs
    Mgt = ((Mgt >> 16) & 0x0000FFFF0000FFFF) | ((Mgt & 0x0000FFFF0000FFFF) << 16);
    Mgt = ( Mgt >> 32                      ) | ( Mgt                       << 32);

    if (num_rows < 64) {
      std::bitset<64> tmp_bitset(Mgt);
      for (int i = 0; i < num_rows; ++i) {
        skip_bitset.push_back(tmp_bitset[i]);
      }
      break;
    }
    skip_bitset.append(Mgt);
    num_rows -= 64;
  }
}

inline void FleDecoder::Ge(int64_t num_rows, dynamic_bitset<>& skip_bitset,
    uint64_t value) {
  vector<uint64_t> C;
  for (int i = 0; i < bit_width_; ++i) {
    if (value & 0x01 << i) {
      C.push_back(~0x0);
    } else {
      C.push_back(0x0);
    }
  }

  for (int i = count_; i != 64 && num_rows > 0; ++i, --num_rows) {
    if (bit_width_ <= 8) {
      skip_bitset.push_back(pcurrent_value8_[i] >= value);
    } else if (bit_width_ <= 16) {
      skip_bitset.push_back(pcurrent_value16_[i] >= value);
    } else {
      skip_bitset.push_back(current_value_[i] >= value);
    }
  }

  uint64_t start_idx = 0;
  while (num_rows > 0) {
    uint64_t Mgt = 0x0;
    uint64_t Meq = ~0x0;
    int i = bit_width_ - 1;
    int j = start_idx + bit_width_ - 1;
    for (; i >= 0; --i, --j) {
      Mgt = Mgt | (Meq & ~C[i] & buffer_end_[j]);
      Meq = Meq & ~(buffer_end_[j] ^ C[i]);
    }
    start_idx += bit_width_;
    uint64_t Mge = Mgt|Meq;
    Mge = ((Mge >> 1) & 0x5555555555555555) | ((Mge & 0x5555555555555555) << 1);
    // swap consecutiMgee pairs
    Mge = ((Mge >> 2) & 0x3333333333333333) | ((Mge & 0x3333333333333333) << 2);
    // swap nibbles ...
    Mge = ((Mge >> 4) & 0x0F0F0F0F0F0F0F0F) | ((Mge & 0x0F0F0F0F0F0F0F0F) << 4);
    // swap bytes
    Mge = ((Mge >> 8) & 0x00FF00FF00FF00FF) | ((Mge & 0x00FF00FF00FF00FF) << 8);
    // swap 2-byte long pairs
    Mge = ((Mge >> 16) & 0x0000FFFF0000FFFF) | ((Mge & 0x0000FFFF0000FFFF) << 16);
    Mge = ( Mge >> 32                      ) | ( Mge                       << 32);

    if (num_rows < 64) {
      std::bitset<64> tmp_bitset(Mge);
      for (int i = 0; i < num_rows; ++i) {
        skip_bitset.push_back(tmp_bitset[i]);
      }
      break;
    }
    skip_bitset.append(Mge);
    num_rows -= 64;
  }
}

inline void FleDecoder::In(int64_t num_rows, dynamic_bitset<>& skip_bitset,
    vector<uint64_t>& values) {
  int values_size = values.size();
  vector<vector<uint64_t> > VC;
  for (int j = 0; j < values_size; ++j) {
    vector<uint64_t> C;
    for (int i = 0; i < bit_width_; ++i) {
      if (values[j] & 0x01 << i) {
        C.push_back(~0x0);
      } else {
        C.push_back(0x0);
      }
    }
    VC.push_back(C);
  }

  bool flag = false;
  if (bit_width_ <= 8) {
    for (int i = count_; i != 64 && num_rows > 0; ++i, --num_rows) {
      flag = false;
      for (int j = 0; j < values_size; ++j) {
        flag = flag || pcurrent_value8_[i] == values[j];
      }
      skip_bitset.push_back(flag);
    }
  } else if (bit_width_ <= 16) {
    for (int i = count_; i != 64 && num_rows > 0; ++i, --num_rows) {
      flag = false;
      for (int j = 0; j < values_size; ++j) {
        flag = flag || pcurrent_value16_[i] == values[j];
      }
      skip_bitset.push_back(flag);
    }
  } else {
    for (int i = count_; i != 64 && num_rows > 0; ++i, --num_rows) {
      flag = false;
      for (int j = 0; j < values_size; ++j) {
        flag = flag || current_value_[i] == values[j];
      }
      skip_bitset.push_back(flag);
    }
  }

  uint64_t j = 0;
  while (num_rows > 0) {
    uint64_t start_j = j;
    uint64_t Many_eq = 0x0;
    for (int v = 0; v < values_size; ++v) {
      j = start_j;
      uint64_t Meq = ~0x0;
      for (int i = 0; i < bit_width_; ++i, ++j) {
        Meq = Meq & ~(buffer_end_[j] ^ VC[v][i]);
      }
      Many_eq |= Meq;
    }

    Many_eq = ((Many_eq >> 1) & 0x5555555555555555) | ((Many_eq & 0x5555555555555555) << 1);
    // swap consecutiMany_eqe pairs
    Many_eq = ((Many_eq >> 2) & 0x3333333333333333) | ((Many_eq & 0x3333333333333333) << 2);
    // swap nibbles ...
    Many_eq = ((Many_eq >> 4) & 0x0F0F0F0F0F0F0F0F) | ((Many_eq & 0x0F0F0F0F0F0F0F0F) << 4);
    // swap bytes
    Many_eq = ((Many_eq >> 8) & 0x00FF00FF00FF00FF) | ((Many_eq & 0x00FF00FF00FF00FF) << 8);
    // swap 2-byte long pairs
    Many_eq = ((Many_eq >> 16) & 0x0000FFFF0000FFFF) | ((Many_eq & 0x0000FFFF0000FFFF) << 16);
    Many_eq = ( Many_eq >> 32                      ) | ( Many_eq                       << 32);

    if (num_rows < 64) {
      std::bitset<64> tmp_bitset(Many_eq);
      for (int i = 0; i < num_rows; ++i) {
        skip_bitset.push_back(tmp_bitset[i]);
      }
      break;
    }
    skip_bitset.append(Many_eq);
    num_rows -= 64;
  }
}

inline bool FleEncoder::Put(uint64_t value) {
  DCHECK(bit_width_ == 64 || value < (1LL << bit_width_));
  if (UNLIKELY(count_ == 64)) {
    if (buffer_end_ >= buffer_end_guard_) return false;

//    for (int i = 0; i < bit_width_; ++i) {
//      buffer_end_[i] = 0;
//    }

    (this->*fa[bit_width_])();

    buffer_end_ += bit_width_;
    count_ = 0;
  }

//  uint64_t set1 = 0x01ULL << (63 - count_);
//  uint64_t index = _tzcnt_u64(value);
//  while (index != 64) {
//    buffer_end_[index] |= set1;
//    value = _blsr_u64(value);
//    index = _tzcnt_u64(value);
//  }

//  for (int l = 0; l < bit_width_; ++l) {
//    buffer_end_[l] |= ((value >> l) & 0x01ULL) << (63 - count_);
//  }

//  uint64_t set1 = 0x01ULL << (63 - count_);
//  for (int l = 0; l < bit_width_; ++l) {
//    if (value & (0x01ULL << l)) buffer_end_[l] |= set1;
//  }


//  for (int i = 0; i < 64 -  _lzcnt_u64(value); ++i) {
//  for (int i = 0; i < bit_width_; ++i) {
//    *current_value_[i] |= ((value >> i) & 0x01) << (63 - count_);
//  }


  if (bit_width_ <= 8) {
    pcurrent_value8_[63 - count_] = value;
  } else if (bit_width_ <= 16) {
    pcurrent_value16_[63 - count_] = value;
  } else {
    current_value_[63 - count_] = value;
//    current_value_[count_] = value;
  }

  ++count_;
  return true;
}

inline void FleEncoder::Pack_1() {
  __m256i dat, ret;
  dat = _mm256_loadu_si256((__m256i const*)pcurrent_value8_);
  ret = _mm256_cmpeq_epi8(dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  dat = _mm256_loadu_si256((__m256i const*)(pcurrent_value8_ + 32));
  ret = _mm256_cmpeq_epi8(dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);
}

inline void FleEncoder::Pack_2() {
  __m256i dat, dat1, tmp_dat, ret;
  dat = _mm256_loadu_si256((__m256i const*)pcurrent_value8_);
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value8_ + 32));

  tmp_dat = _mm256_and_si256(dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);
}

inline void FleEncoder::Pack_3() {
  __m256i dat, dat1, tmp_dat, ret;
  dat = _mm256_loadu_si256((__m256i const*)pcurrent_value8_);
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value8_ + 32));

  tmp_dat = _mm256_and_si256(dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);
}

inline void FleEncoder::Pack_4() {
  __m256i dat, dat1, tmp_dat, ret;
  dat = _mm256_loadu_si256((__m256i const*)pcurrent_value8_);
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value8_ + 32));

  tmp_dat = _mm256_and_si256(dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[6] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[7] = _mm256_movemask_epi8(ret);
}

inline void FleEncoder::Pack_5() {
  __m256i dat, dat1, tmp_dat, ret;
  dat = _mm256_loadu_si256((__m256i const*)pcurrent_value8_);
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value8_ + 32));

  tmp_dat = _mm256_and_si256(dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[6] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[7] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[8] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[9] = _mm256_movemask_epi8(ret);
}

inline void FleEncoder::Pack_6() {
  __m256i dat, dat1, tmp_dat, ret;
  dat = _mm256_loadu_si256((__m256i const*)pcurrent_value8_);
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value8_ + 32));

  tmp_dat = _mm256_and_si256(dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[6] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[7] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[8] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[9] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[10] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[11] = _mm256_movemask_epi8(ret);
}

inline void FleEncoder::Pack_7() {
  __m256i dat, dat1, tmp_dat, ret;
  dat = _mm256_loadu_si256((__m256i const*)pcurrent_value8_);
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value8_ + 32));

  tmp_dat = _mm256_and_si256(dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[6] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[7] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[8] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[9] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[10] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[11] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[12] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[13] = _mm256_movemask_epi8(ret);
}

inline void FleEncoder::Pack_8() {
  __m256i dat, dat1, tmp_dat, ret;
  dat = _mm256_loadu_si256((__m256i const*)pcurrent_value8_);
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value8_ + 32));

  tmp_dat = _mm256_and_si256(dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[6] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[7] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[8] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[9] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[10] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[11] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[12] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[13] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[14] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat1, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[15] = _mm256_movemask_epi8(ret);
}

inline void FleEncoder::Pack_9() {
  __m256i dat, dat1, dat2, dat3, pack_dat, pack_dat1, tmp_dat, tmp_dat1, ret;
  dat = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_));
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 16));
  tmp_dat = _mm256_and_si256(dat, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  dat2 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 32));
  dat3 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 48));
  tmp_dat = _mm256_and_si256(dat2, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat3, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[6] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[7] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[8] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[9] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[10] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[11] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[12] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[13] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[14] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[15] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_srli_si256(dat, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat1, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  tmp_dat = _mm256_srli_si256(dat2, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat3, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[16] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[17] = _mm256_movemask_epi8(ret);
}


inline void FleEncoder::Pack_10() {
  __m256i dat, dat1, dat2, dat3, pack_dat, pack_dat1, tmp_dat, tmp_dat1, ret;
  dat = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_));
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 16));
  tmp_dat = _mm256_and_si256(dat, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  dat2 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 32));
  dat3 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 48));
  tmp_dat = _mm256_and_si256(dat2, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat3, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[6] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[7] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[8] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[9] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[10] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[11] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[12] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[13] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[14] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[15] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_srli_si256(dat, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat1, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  tmp_dat = _mm256_srli_si256(dat2, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat3, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[16] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[17] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[18] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[19] = _mm256_movemask_epi8(ret);
}

inline void FleEncoder::Pack_11() {
  __m256i dat, dat1, dat2, dat3, tmp_dat, tmp_dat1, pack_dat, pack_dat1, ret;
  dat = _mm256_loadu_si256((__m256i const*)pcurrent_value16_);
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 16));
  tmp_dat = _mm256_and_si256(dat, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  dat2 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 32));
  dat3 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 48));
  tmp_dat = _mm256_and_si256(dat2, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat3, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[6] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[7] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[8] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[9] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[10] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[11] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[12] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[13] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[14] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[15] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_srli_si256(dat, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat1, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  tmp_dat = _mm256_srli_si256(dat2, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat3, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[16] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[17] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[18] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[19] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[20] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[21] = _mm256_movemask_epi8(ret);
}


inline void FleEncoder::Pack_12() {
  __m256i dat, dat1, dat2, dat3, tmp_dat, tmp_dat1, pack_dat, pack_dat1, ret;
  dat = _mm256_loadu_si256((__m256i const*)pcurrent_value16_);
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 16));
  tmp_dat = _mm256_and_si256(dat, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  dat2 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 32));
  dat3 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 48));
  tmp_dat = _mm256_and_si256(dat2, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat3, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[6] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[7] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[8] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[9] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[10] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[11] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[12] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[13] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[14] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[15] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_srli_si256(dat, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat1, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  tmp_dat = _mm256_srli_si256(dat2, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat3, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[16] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[17] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[18] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[19] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[20] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[21] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[22] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[23] = _mm256_movemask_epi8(ret);
}

inline void FleEncoder::Pack_13() {
  __m256i dat, dat1, dat2, dat3, tmp_dat, tmp_dat1, pack_dat, pack_dat1, ret;
  dat = _mm256_loadu_si256((__m256i const*)pcurrent_value16_);
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 16));
  tmp_dat = _mm256_and_si256(dat, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  dat2 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 32));
  dat3 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 48));
  tmp_dat = _mm256_and_si256(dat2, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat3, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[6] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[7] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[8] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[9] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[10] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[11] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[12] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[13] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[14] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[15] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_srli_si256(dat, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat1, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  tmp_dat = _mm256_srli_si256(dat2, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat3, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[16] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[17] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[18] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[19] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[20] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[21] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[22] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[23] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[24] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[25] = _mm256_movemask_epi8(ret);
}

inline void FleEncoder::Pack_14() {
  __m256i dat, dat1, dat2, dat3, tmp_dat, tmp_dat1, pack_dat, pack_dat1, ret;
  dat = _mm256_loadu_si256((__m256i const*)pcurrent_value16_);
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 16));
  tmp_dat = _mm256_and_si256(dat, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  dat2 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 32));
  dat3 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 48));
  tmp_dat = _mm256_and_si256(dat2, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat3, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[6] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[7] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[8] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[9] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[10] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[11] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[12] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[13] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[14] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[15] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_srli_si256(dat, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat1, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  tmp_dat = _mm256_srli_si256(dat2, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat3, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[16] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[17] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[18] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[19] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[20] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[21] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[22] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[23] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[24] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[25] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[26] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[27] = _mm256_movemask_epi8(ret);
}

inline void FleEncoder::Pack_15() {
  __m256i dat, dat1, dat2, dat3, tmp_dat, tmp_dat1, pack_dat, pack_dat1, ret;
  dat = _mm256_loadu_si256((__m256i const*)pcurrent_value16_);
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 16));
  tmp_dat = _mm256_and_si256(dat, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  dat2 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 32));
  dat3 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 48));
  tmp_dat = _mm256_and_si256(dat2, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat3, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[6] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[7] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[8] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[9] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[10] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[11] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[12] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[13] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[14] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[15] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_srli_si256(dat, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat1, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  tmp_dat = _mm256_srli_si256(dat2, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat3, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[16] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[17] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[18] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[19] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[20] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[21] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[22] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[23] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[24] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[25] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[26] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[27] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[28] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[29] = _mm256_movemask_epi8(ret);
}

inline void FleEncoder::Pack_16() {
  __m256i dat, dat1, dat2, dat3, tmp_dat, tmp_dat1, pack_dat, pack_dat1, ret;
  dat = _mm256_loadu_si256((__m256i const*)pcurrent_value16_);
  dat1 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 16));
  tmp_dat = _mm256_and_si256(dat, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  dat2 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 32));
  dat3 = _mm256_loadu_si256((__m256i const*)(pcurrent_value16_ + 48));
  tmp_dat = _mm256_and_si256(dat2, clr_mask);
  tmp_dat1 = _mm256_and_si256(dat3, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[6] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[7] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[8] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[9] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[10] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[11] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[12] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[13] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[14] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[15] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_srli_si256(dat, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat1, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat = _mm256_permute4x64_epi64(pack_dat, 0xd8);
  tmp_dat = _mm256_srli_si256(dat2, 1);
  tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
  tmp_dat1 = _mm256_srli_si256(dat3, 1);
  tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
  pack_dat1 = _mm256_packus_epi16(tmp_dat, tmp_dat1);
  pack_dat1 = _mm256_permute4x64_epi64(pack_dat1, 0xd8);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[16] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[17] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[18] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[19] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[20] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[21] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[22] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[23] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[24] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[25] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[26] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[27] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[28] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[29] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(pack_dat, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[30] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(pack_dat1, bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[31] = _mm256_movemask_epi8(ret);
}

inline void FleEncoder::Pack_16_32() {
//  memset(buffer_end_, 0, bit_width_ * 8);
//  for (int i = 0; i < count_; ++i) {
//    for (int l = 0; l < bit_width_; ++l) {
//      buffer_end_[l] |= ((current_value_[i] >> l) & 0x01ULL) << (63 - i);
//    }
//  }


  __m256i dat[8], tmp_dat, tmp_dat1, pack_dat[8], ret;
  for (int l = 0; l < 8; ++l) {
    dat[l] =  _mm256_loadu_si256((__m256i const*)(current_value_ + 8 * l));
  }
  for (int l = 0; l < 4; ++l) {
    tmp_dat = _mm256_and_si256(dat[2 * l], clr32_mask);
    tmp_dat1 = _mm256_and_si256(dat[2 * l + 1], clr32_mask);
    pack_dat[l] = _mm256_packus_epi32(tmp_dat, tmp_dat1);
    pack_dat[l] = _mm256_permute4x64_epi64(pack_dat[l], 0xd8);
    tmp_dat = _mm256_srli_si256(dat[2 * l], 2);
    tmp_dat = _mm256_and_si256(tmp_dat, clr32_mask);
    tmp_dat1 = _mm256_srli_si256(dat[2 * l + 1], 2);
    tmp_dat1 = _mm256_and_si256(tmp_dat1, clr32_mask);
    pack_dat[4 + l] = _mm256_packus_epi32(tmp_dat, tmp_dat1);
    pack_dat[4 + l] = _mm256_permute4x64_epi64(pack_dat[4 + l], 0xd8);
  }
  for (int l = 0; l < 4; ++l) {
    tmp_dat = _mm256_and_si256(pack_dat[2 * l], clr_mask);
    tmp_dat1 = _mm256_and_si256(pack_dat[2 * l + 1], clr_mask);
    dat[l] = _mm256_packus_epi16(tmp_dat, tmp_dat1);
    dat[l] = _mm256_permute4x64_epi64(dat[l], 0xd8);
    tmp_dat = _mm256_srli_si256(pack_dat[2 * l], 1);
    tmp_dat = _mm256_and_si256(tmp_dat, clr_mask);
    tmp_dat1 = _mm256_srli_si256(pack_dat[2 * l + 1], 1);
    tmp_dat1 = _mm256_and_si256(tmp_dat1, clr_mask);
    dat[4 + l] = _mm256_packus_epi16(tmp_dat, tmp_dat1);
    dat[4 + l] = _mm256_permute4x64_epi64(dat[4 + l], 0xd8);
  }

// bit 0-7
  tmp_dat = _mm256_and_si256(dat[0], bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[0] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[1], bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[1] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(dat[0], bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[2] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[1], bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[3] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(dat[0], bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[4] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[1], bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[5] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(dat[0], bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[6] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[1], bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[7] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(dat[0], bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[8] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[1], bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[9] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(dat[0], bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[10] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[1], bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[11] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(dat[0], bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[12] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[1], bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[13] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(dat[0], bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[14] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[1], bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[15] = _mm256_movemask_epi8(ret);

// bit 8-15
  tmp_dat = _mm256_and_si256(dat[4], bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[16] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[5], bitv[0]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
  reinterpret_cast<uint32_t*>(buffer_end_)[17] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(dat[4], bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[18] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[5], bitv[1]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
  reinterpret_cast<uint32_t*>(buffer_end_)[19] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(dat[4], bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[20] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[5], bitv[2]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
  reinterpret_cast<uint32_t*>(buffer_end_)[21] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(dat[4], bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[22] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[5], bitv[3]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
  reinterpret_cast<uint32_t*>(buffer_end_)[23] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(dat[4], bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[24] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[5], bitv[4]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
  reinterpret_cast<uint32_t*>(buffer_end_)[25] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(dat[4], bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[26] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[5], bitv[5]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
  reinterpret_cast<uint32_t*>(buffer_end_)[27] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(dat[4], bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[28] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[5], bitv[6]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
  reinterpret_cast<uint32_t*>(buffer_end_)[29] = _mm256_movemask_epi8(ret);

  tmp_dat = _mm256_and_si256(dat[4], bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[30] = _mm256_movemask_epi8(ret);
  tmp_dat = _mm256_and_si256(dat[5], bitv[7]);
  ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
  reinterpret_cast<uint32_t*>(buffer_end_)[31] = _mm256_movemask_epi8(ret);

// bit 16-23
  if (bit_width_ >= 24) {
    tmp_dat = _mm256_and_si256(dat[2], bitv[0]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
    reinterpret_cast<uint32_t*>(buffer_end_)[32] = _mm256_movemask_epi8(ret);
    tmp_dat = _mm256_and_si256(dat[3], bitv[0]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[0]);
    reinterpret_cast<uint32_t*>(buffer_end_)[33] = _mm256_movemask_epi8(ret);

    tmp_dat = _mm256_and_si256(dat[2], bitv[1]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
    reinterpret_cast<uint32_t*>(buffer_end_)[34] = _mm256_movemask_epi8(ret);
    tmp_dat = _mm256_and_si256(dat[3], bitv[1]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[1]);
    reinterpret_cast<uint32_t*>(buffer_end_)[35] = _mm256_movemask_epi8(ret);

    tmp_dat = _mm256_and_si256(dat[2], bitv[2]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
    reinterpret_cast<uint32_t*>(buffer_end_)[36] = _mm256_movemask_epi8(ret);
    tmp_dat = _mm256_and_si256(dat[3], bitv[2]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[2]);
    reinterpret_cast<uint32_t*>(buffer_end_)[37] = _mm256_movemask_epi8(ret);

    tmp_dat = _mm256_and_si256(dat[2], bitv[3]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
    reinterpret_cast<uint32_t*>(buffer_end_)[38] = _mm256_movemask_epi8(ret);
    tmp_dat = _mm256_and_si256(dat[3], bitv[3]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[3]);
    reinterpret_cast<uint32_t*>(buffer_end_)[39] = _mm256_movemask_epi8(ret);

    tmp_dat = _mm256_and_si256(dat[2], bitv[4]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
    reinterpret_cast<uint32_t*>(buffer_end_)[40] = _mm256_movemask_epi8(ret);
    tmp_dat = _mm256_and_si256(dat[3], bitv[4]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[4]);
    reinterpret_cast<uint32_t*>(buffer_end_)[41] = _mm256_movemask_epi8(ret);

    tmp_dat = _mm256_and_si256(dat[2], bitv[5]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
    reinterpret_cast<uint32_t*>(buffer_end_)[42] = _mm256_movemask_epi8(ret);
    tmp_dat = _mm256_and_si256(dat[3], bitv[5]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[5]);
    reinterpret_cast<uint32_t*>(buffer_end_)[43] = _mm256_movemask_epi8(ret);

    tmp_dat = _mm256_and_si256(dat[2], bitv[6]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
    reinterpret_cast<uint32_t*>(buffer_end_)[44] = _mm256_movemask_epi8(ret);
    tmp_dat = _mm256_and_si256(dat[3], bitv[6]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[6]);
    reinterpret_cast<uint32_t*>(buffer_end_)[45] = _mm256_movemask_epi8(ret);

    tmp_dat = _mm256_and_si256(dat[2], bitv[7]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
    reinterpret_cast<uint32_t*>(buffer_end_)[46] = _mm256_movemask_epi8(ret);
    tmp_dat = _mm256_and_si256(dat[3], bitv[7]);
    ret = _mm256_cmpeq_epi8(tmp_dat, bitv[7]);
    reinterpret_cast<uint32_t*>(buffer_end_)[47] = _mm256_movemask_epi8(ret);

    int idx = 48;
    for (int l = 0; l < bit_width_ - 24; ++l) {
      tmp_dat = _mm256_and_si256(dat[6], bitv[l]);
      ret = _mm256_cmpeq_epi8(tmp_dat, bitv[l]);
      reinterpret_cast<uint32_t*>(buffer_end_)[idx] = _mm256_movemask_epi8(ret);
      tmp_dat = _mm256_and_si256(dat[7], bitv[l]);
      ret = _mm256_cmpeq_epi8(tmp_dat, bitv[l]);
      reinterpret_cast<uint32_t*>(buffer_end_)[idx + 1] = _mm256_movemask_epi8(ret);
      idx += 2;
    }
  } else {
    int idx = 32;
    for (int l = 0; l < bit_width_ - 16; ++l) {
      tmp_dat = _mm256_and_si256(dat[2], bitv[l]);
      ret = _mm256_cmpeq_epi8(tmp_dat, bitv[l]);
      reinterpret_cast<uint32_t*>(buffer_end_)[idx] = _mm256_movemask_epi8(ret);
      tmp_dat = _mm256_and_si256(dat[3], bitv[l]);
      ret = _mm256_cmpeq_epi8(tmp_dat, bitv[l]);
      reinterpret_cast<uint32_t*>(buffer_end_)[idx + 1] = _mm256_movemask_epi8(ret);
      idx += 2;
    }
  }

}


inline int FleEncoder::Flush() {
  (this->*fa[bit_width_])();

  buffer_end_ += bit_width_;
  count_ = 0;
  return (buffer_end_ - buffer_) * 8;
}

}

#endif

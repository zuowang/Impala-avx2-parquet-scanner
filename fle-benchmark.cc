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

#include <iostream>
#include <sstream>
#include <immintrin.h>

#include "runtime/mem-pool.h"
#include "runtime/mem-tracker.h"
#include "experiments/bit-stream-utils.8byte.inline.h"
#include "util/benchmark.h"
#include "util/bit-stream-utils.inline.h"
#include "util/cpu-info.h"
#include "util/fle-encoding.h"
#include "util/rle-encoding.h"

#include "common/names.h"

// Benchmark to measure how quickly we can do bit encoding and decoding.

using namespace impala;

const int BUFFER_LEN = 64 * 4096;

struct TestData {
  uint64_t* buffer;
  uint64_t* buffer_end;
  uint32_t* buffer32_end;
  int buffer_len;
  int bit_width;
  uint64_t* array;
  uint64_t* array_end;
  uint32_t* array32_end;
  uint16_t* array16_end;
  uint8_t* array8_end;
  int num_values;
  int max_value;
  MemPool* pool;
  bool result;
};
/*
class BaseEncoder {
 public:
virtual inline void encode(uint64_t* buffer, int count) = 0;
virtual ~BaseEncoder() {}
static BaseEncoder** encoders;
};
BaseEncoder** BaseEncoder::encoders = NULL;

template<uint64_t val>
class Encoder : public BaseEncoder {
 public:
inline void encode(uint64_t* buffer, int count) {
  uint64_t v = 0x01ULL << (63 - count);
  if (val & (0x01 << 0)) buffer[0] |= v;
  if (val & (0x01 << 1)) buffer[1] |= v;
  if (val & (0x01 << 2)) buffer[2] |= v;
  if (val & (0x01 << 3)) buffer[3] |= v;
  if (val & (0x01 << 4)) buffer[4] |= v;
  if (val & (0x01 << 5)) buffer[5] |= v;
  if (val & (0x01 << 6)) buffer[6] |= v;
  if (val & (0x01 << 7)) buffer[7] |= v;
  if (val & (0x01 << 8)) buffer[8] |= v;
  if (val & (0x01 << 9)) buffer[9] |= v;
  if (val & (0x01 << 10)) buffer[10] |= v;
  if (val & (0x01 << 11)) buffer[11] |= v;
  if (val & (0x01 << 12)) buffer[12] |= v;
  if (val & (0x01 << 13)) buffer[13] |= v;
  if (val & (0x01 << 14)) buffer[14] |= v;
  if (val & (0x01 << 15)) buffer[15] |= v;
  if (val & (0x01 << 16)) buffer[16] |= v;
  if (val & (0x01 << 17)) buffer[17] |= v;
  if (val & (0x01 << 18)) buffer[18] |= v;
  if (val & (0x01 << 19)) buffer[19] |= v;
  if (val & (0x01 << 20)) buffer[20] |= v;
  if (val & (0x01 << 21)) buffer[21] |= v;
  if (val & (0x01 << 22)) buffer[22] |= v;
  if (val & (0x01 << 23)) buffer[23] |= v;
  if (val & (0x01 << 24)) buffer[24] |= v;
  if (val & (0x01 << 25)) buffer[25] |= v;
  if (val & (0x01 << 26)) buffer[26] |= v;
  if (val & (0x01 << 27)) buffer[27] |= v;
  if (val & (0x01 << 28)) buffer[28] |= v;
  if (val & (0x01 << 29)) buffer[29] |= v;
  if (val & (0x01 << 30)) buffer[30] |= v;
  if (val & (0x01 << 31)) buffer[31] |= v;
}

};

template<uint64_t val>
class GenEncoder : public GenEncoder<val - 1> {
 public:
GenEncoder(int n = val) : GenEncoder<val - 1>(n) {
  BaseEncoder::encoders[val] = new Encoder<val>();
}
~GenEncoder() {
  delete BaseEncoder::encoders[val];
}
};

template<>
class GenEncoder<0> {
 public:
GenEncoder<0>(int n = 0) {
  BaseEncoder::encoders = new BaseEncoder*[n + 1];
  BaseEncoder::encoders[0] = new Encoder<0>();
}
~GenEncoder<0>() {
  delete BaseEncoder::encoders[0];
  delete [] BaseEncoder::encoders;
  BaseEncoder::encoders = NULL;
}
};

void TestFleTableEncode(int batch_size, void* d) {
  TestData* data = reinterpret_cast<TestData*>(d);
  for (int i = 0; i < batch_size; ++i) {
    // Unroll this to focus more on Put performance.
    if (data->bit_width <= 8) {
      for (int j = 0; j < data->num_values; j += 64) {
        for (int k = 0; k < 64; ++k) {
          BaseEncoder::encoders[data->array_end[k]]->encode(data->buffer_end, k);
        }
        data->buffer_end += data->bit_width;
        data->array_end += 64;
      }
    } else if (data->bit_width <= 16) {
      for (int j = 0; j < data->num_values; j += 64) {
        for (int k = 0; k < 64; ++k) {
          BaseEncoder::encoders[data->array_end[k] & 0xff]->encode(data->buffer_end, k);
          BaseEncoder::encoders[(data->array_end[k] >> 8) & 0xff]->encode(data->buffer_end + 8, k);
        }
        data->buffer_end += data->bit_width;
        data->array_end += 64;
      }
    }

    data->buffer_end = data->buffer;
    data->array_end = data->array;
  }
}
*/

const uint64_t SHIFTMASK[64] = {
0x1,
0x2,
0x4,
0x8,
0x10,
0x20,
0x40,
0x80,
0x100,
0x200,
0x400,
0x800,
0x1000,
0x2000,
0x4000,
0x8000,
0x10000,
0x20000,
0x40000,
0x80000,
0x100000,
0x200000,
0x400000,
0x800000,
0x1000000,
0x2000000,
0x4000000,
0x8000000,
0x10000000,
0x20000000,
0x40000000,
0x80000000,
0x100000000,
0x200000000,
0x400000000,
0x800000000,
0x1000000000,
0x2000000000,
0x4000000000,
0x8000000000,
0x10000000000,
0x20000000000,
0x40000000000,
0x80000000000,
0x100000000000,
0x200000000000,
0x400000000000,
0x800000000000,
0x1000000000000,
0x2000000000000,
0x4000000000000,
0x8000000000000,
0x10000000000000,
0x20000000000000,
0x40000000000000,
0x80000000000000,
0x100000000000000,
0x200000000000000,
0x400000000000000,
0x800000000000000,
0x1000000000000000,
0x2000000000000000,
0x4000000000000000,
0x8000000000000000
};


void TestFleWithAVXV2Encode(int batch_size, void* d) {
  TestData* data = reinterpret_cast<TestData*>(d);
  for (int i = 0; i < batch_size; ++i) {
    // Unroll this to focus more on Put performance.
    if (data->bit_width == 1) {
    __m256i bit = _mm256_set1_epi64x(0x0101010101010101);
    for (int j = 0; j < data->num_values; j += 64) {
      __m256i dat = _mm256_loadu_si256((__m256i const*)data->array8_end);
      __m256i ret = _mm256_cmpeq_epi8(dat, bit);
      data->buffer32_end[0] = _mm256_movemask_epi8(ret);
      dat = _mm256_loadu_si256((__m256i const*)(data->array8_end + 32));
      ret = _mm256_cmpeq_epi8(dat, bit);
      data->buffer32_end[1] = _mm256_movemask_epi8(ret);
      data->buffer32_end += 2;
      data->array8_end += 64;
    }
    data->buffer32_end = reinterpret_cast<uint32_t*>(data->buffer);
    data->array8_end = reinterpret_cast<uint8_t*>(data->array);
    } else if (data->bit_width <= 8) {
    __m256i bitv[8];
    for (int l = 0; l < 8; ++l) {
      bitv[l] = _mm256_set1_epi64x(0x0101010101010101 << l);
    }
    for (int j = 0; j < data->num_values; j += 64) {
      __m256i dat, dat1, tmp_dat, ret;
      dat = _mm256_loadu_si256((__m256i const*)data->array8_end);
      dat1 = _mm256_loadu_si256((__m256i const*)(data->array8_end + 32));
      for (int l = 0; l < data->bit_width; ++l) {
        tmp_dat = _mm256_and_si256(dat, bitv[l]);
        ret = _mm256_cmpeq_epi8(tmp_dat, bitv[l]);
        data->buffer32_end[2 * l] = _mm256_movemask_epi8(ret);
        tmp_dat = _mm256_and_si256(dat1, bitv[l]);
        ret = _mm256_cmpeq_epi8(tmp_dat, bitv[l]);
        data->buffer32_end[2 * l + 1] = _mm256_movemask_epi8(ret);
      }
      data->buffer32_end += 2 * data->bit_width;
      data->array8_end += 64;
    }
    data->buffer32_end = reinterpret_cast<uint32_t*>(data->buffer);
    data->array8_end = reinterpret_cast<uint8_t*>(data->array);
    } else if (data->bit_width <= 16) {
    __m256i bitv[8];
    for (int l = 0; l < 8; ++l) {
      bitv[l] = _mm256_set1_epi64x(0x0101010101010101 << l);
    }
    for (int j = 0; j < data->num_values; j += 64) {
      __m256i dat, dat1, dat2, dat3, tmp_dat, ret, cmp_dat, cmp_dat1;
      dat = _mm256_loadu_si256((__m256i const*)data->array16_end);
      dat1 = _mm256_loadu_si256((__m256i const*)(data->array16_end + 16));
      cmp_dat = _mm256_packus_epi16(dat, dat1);
      dat2 = _mm256_loadu_si256((__m256i const*)(data->array16_end + 32));
      dat3 = _mm256_loadu_si256((__m256i const*)(data->array16_end + 48));
      cmp_dat1 = _mm256_packus_epi16(dat2, dat3);
      for (int l = 0; l < 8; ++l) {
        tmp_dat = _mm256_and_si256(cmp_dat, bitv[l]);
        ret = _mm256_cmpeq_epi8(tmp_dat, bitv[l]);
        data->buffer32_end[2 * l] = _mm256_movemask_epi8(ret);
        tmp_dat = _mm256_and_si256(cmp_dat1, bitv[l]);
        ret = _mm256_cmpeq_epi8(tmp_dat, bitv[l]);
        data->buffer32_end[2 * l + 1] = _mm256_movemask_epi8(ret);
      }
      data->buffer32_end += 16;
      dat = _mm256_srli_si256(dat, 1);
      dat1 = _mm256_srli_si256(dat1, 1);
      cmp_dat = _mm256_packus_epi16(dat, dat1);

      dat2 = _mm256_srli_si256(dat2, 1);
      dat3 = _mm256_srli_si256(dat3, 1);
      cmp_dat1 = _mm256_packus_epi16(dat2, dat3);
      for (int l = 0; l < data->bit_width - 8; ++l) {
        tmp_dat = _mm256_and_si256(cmp_dat, bitv[l]);
        ret = _mm256_cmpeq_epi8(tmp_dat, bitv[l]);
        data->buffer32_end[2 * l] = _mm256_movemask_epi8(ret);
        tmp_dat = _mm256_and_si256(cmp_dat1, bitv[l]);
        ret = _mm256_cmpeq_epi8(tmp_dat, bitv[l]);
        data->buffer32_end[2 * l + 1] = _mm256_movemask_epi8(ret);
      }
      data->buffer32_end += 2 * (data->bit_width - 8);
      data->array16_end += 64;
    }
    data->buffer32_end = reinterpret_cast<uint32_t*>(data->buffer);
    data->array16_end = reinterpret_cast<uint16_t*>(data->array);
    } else {
    __m256i bitv[8];
    for (int l = 0; l < 8; ++l) {
      bitv[l] = _mm256_set1_epi64x(0x0101010101010101 << l);
    }
    for (int j = 0; j < data->num_values; j += 64) {
      __m256i tmp_dat, ret;
      __m256i dat[8];
      for (int l = 0; l < 8; ++l) {
        dat[i] =  _mm256_loadu_si256((__m256i const*)(data->array32_end + 8 * i));
      }
      __m256i cmp_dat[8];
      for (int l = 0; l < 4; ++l) {
        cmp_dat[l] = _mm256_packus_epi32(dat[2 * l], dat[2 * l + 1]);
        dat[2 * l] = _mm256_srli_si256(dat[2 * l], 2);
        dat[2 * l + 1] = _mm256_srli_si256(dat[2 * l + 1], 2);
        cmp_dat[4 + l] = _mm256_packus_epi32(dat[2 * l], dat[2 * l + 1]);
      }
      for (int l = 0; l < 4; ++l) {
        dat[l] = _mm256_packus_epi16(cmp_dat[2 * l], cmp_dat[2 * l + 1]);
        cmp_dat[2 * l] = _mm256_srli_si256(cmp_dat[2 * l], 1);
        cmp_dat[2 * l + 1] = _mm256_srli_si256(cmp_dat[2 * l + 1], 1);
        dat[4 + l] = _mm256_packus_epi16(cmp_dat[2 * l], cmp_dat[2 * l + 1]);
      }

      for (int l = 0; l < data->bit_width; ++l) {
        int idx = l % 8;
        int idx2 = (l / 8) * 2;
        if (idx2>=7) {printf("%d\n", data->bit_width);exit(-1);}
        tmp_dat = _mm256_and_si256(dat[idx2], bitv[idx]);
        ret = _mm256_cmpeq_epi8(tmp_dat, bitv[idx]);
        data->buffer32_end[2 * l] = _mm256_movemask_epi8(ret);
//        tmp_dat = _mm256_and_si256(dat[idx2 + 1], bitv[idx]);
        ret = _mm256_cmpeq_epi8(tmp_dat, bitv[idx]);
        data->buffer32_end[2 * l + 1] = _mm256_movemask_epi8(ret);
      }
      data->buffer32_end += 2 * data->bit_width;
      data->array32_end += 64;
    }
    data->buffer32_end = reinterpret_cast<uint32_t*>(data->buffer);
    data->array32_end = reinterpret_cast<uint32_t*>(data->array);
    }
  }
}

//void __attribute__((optimize("O1"))) TestFleWithAVXEncode(int batch_size, void* d) {
//  TestData* data = reinterpret_cast<TestData*>(d);
//  for (int i = 0; i < batch_size; ++i) {
//    // Unroll this to focus more on Put performance.
//    for (int j = 0; j < data->num_values; j += 64) {
//      for (int k = 0; k < 64; ++k) {
////        uint64_t set1 = 0x01ULL << (63 - k);
//        uint64_t set1 = SHIFTMASK[63 - k];
//        uint64_t value = data->array_end[k];
//        int index = _lzcnt_u64(value);
//        while (index != 64) {
//          data->buffer_end[index] |= set1;
//          value = _blsr_u64(value);
//          index = _tzcnt_u64(value);
//        }
//      }
//      data->buffer_end += data->bit_width;
//      data->array_end += 64;
//    }
//    data->buffer_end = data->buffer;
//    data->array_end = data->array;
//  }
//}

void TestFleEncode(int batch_size, void* d) {
  TestData* data = reinterpret_cast<TestData*>(d);
  for (int i = 0; i < batch_size; ++i) {
    // Unroll this to focus more on Put performance.
    for (int j = 0; j < data->num_values; j += 64) {
      for (int k = 0; k < 64; ++k) {
//        uint64_t set1 = 0x01ULL << (63 - k);
        uint64_t set1 = SHIFTMASK[63 - k];
        for (int l = 0; l < data->bit_width; ++l) {
          if (data->array_end[k] && (0x01ULL << l)) data->buffer_end[l] |= set1;
        }
      }
      data->buffer_end += data->bit_width;
      data->array_end += 64;
    }
    data->buffer_end = data->buffer;
    data->array_end = data->array;
  }
}

void TestNativeEncode(int batch_size, void* d) {
  TestData* data = reinterpret_cast<TestData*>(d);
  for (int i = 0; i < batch_size; ++i) {
    // Unroll this to focus more on Put performance.
    for (int j = 0; j < data->num_values; j += 64) {
      for (int k = 0; k < 64; ++k) {
        for (int l = 0; l < data->bit_width; ++l) {
          data->buffer_end[l] |= ((data->array_end[k] >> l) & 0x01ULL) << (63 - k);
        }
      }
      data->buffer_end += data->bit_width;
      data->array_end += 64;
    }
    data->buffer_end = data->buffer;
    data->array_end = data->array;
  }
}

//void __attribute__((optimize("O1"))) TestFleWithAVXDecode(int batch_size, void* d) {
//  TestData* data = reinterpret_cast<TestData*>(d);
//  for (int i = 0; i < batch_size; ++i) {
//    // Unroll this to focus more on Put performance.
//    for (int j = 0; j < data->num_values; j += 64) {
//      for (int k = 0; k < data->bit_width; ++k) {
//        uint64_t val = data->buffer_end[k];
////        uint64_t set1 = 0x01ULL << k;
//        uint64_t set1 = SHIFTMASK[k];
//        int index = _tzcnt_u64(val);
//        while (index != 64) {
//          data->array_end[index] |= set1;
//          val = _blsr_u64(val);
//          index = _tzcnt_u64(val);
//        }
//      }
//      data->buffer_end += data->bit_width;
//      data->array_end += 64;
//    }
//    data->buffer_end = data->buffer;
//    data->array_end = data->array;
//  }
//}

void TestFleDecode(int batch_size, void* d) {
  TestData* data = reinterpret_cast<TestData*>(d);
  for (int i = 0; i < batch_size; ++i) {
    // Unroll this to focus more on Put performance.
    for (int j = 0; j < data->num_values; j += 64) {
      for (int k = 0; k < data->bit_width; ++k) {
//        uint64_t set1 = 0x01ULL << k;
        uint64_t set1 = SHIFTMASK[k];
        for (int l = 0; l < 64; ++l) {
//          if (data->buffer_end[k] & (0x01 << l)) data->array_end[l] |= set1;
            if (data->buffer_end[k] & SHIFTMASK[l]) data->array_end[l] |= set1;
        }
      }
      data->buffer_end += data->bit_width;
      data->array_end += 64;
    }
    data->buffer_end = data->buffer;
    data->array_end = data->array;
  }
}

void TestNativeDecode(int batch_size, void* d) {
  TestData* data = reinterpret_cast<TestData*>(d);
  for (int i = 0; i < batch_size; ++i) {
    // Unroll this to focus more on Put performance.
    for (int j = 0; j < data->num_values; j += 64) {
      for (int k = 0; k < data->bit_width; ++k) {
        for (int l = 0; l < 64; ++l) {
          data->array_end[l] |= ((data->buffer_end[k] >> l) & 0x01ULL) << k;
        }
      }
      data->buffer_end += data->bit_width;
      data->array_end += 64;
    }
    data->buffer_end = data->buffer;
    data->array_end = data->array;
  }
}
/*
void TestFleWithAVXV2Decode(int batch_size, void* d) {
  TestData* data = reinterpret_cast<TestData*>(d);
  for (int i = 0; i < batch_size; ++i) {
    // Unroll this to focus more on Put performance.
    if (data->bit_width == 1) {
    for (int j = 0; j < data->num_values; j += 64) {
      for (int l = 0; l < 64; ++l) {
        data->array_end[l] |= ((data->buffer_end[0] >> l) & 0x01ULL);
      }

      data->buffer_end += data->bit_width;
      data->array_end += 64;
    }
    } else if (data->bit_width == 4) {
    for (int j = 0; j < data->num_values; j += 64) {
      for (int k = 0; k < data->bit_width; ++k) {
        for (int l = 0; l < 64; ++l) {
          data->array_end[l] |= ((data->buffer_end[k] >> l) & 0x01ULL) << k;
        }
      }
      data->buffer_end += data->bit_width;
      data->array_end += 64;
    }
    } else {
    for (int j = 0; j < data->num_values; j += 64) {
      for (int k = 0; k < data->bit_width; ++k) {
        for (int l = 0; l < 64; ++l) {
          data->array_end[l] |= ((data->buffer_end[k] >> l) & 0x01ULL) << k;
        }
      }
      data->buffer_end += data->bit_width;
      data->array_end += 64;
    }
    }

    data->buffer_end = data->buffer;
    data->array_end = data->array;
  }
}
*/
void TestFleWithAVXV3Decode(int batch_size, void* d) {
  TestData* data = reinterpret_cast<TestData*>(d);
  for (int i = 0; i < batch_size; ++i) {
    // Unroll this to focus more on Put performance.
    if (data->bit_width == 1) {
    __m256i bit0_mask, bit0, ret;
    __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
        0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);
    for (int j = 0; j < data->num_values; j += 128) {
      __m256i buf = _mm256_loadu_si256((__m256i const*)data->buffer_end);
      __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
      buf = _mm256_inserti128_si256(buf, buf_hi, 1);

      for (int k = 0; k < 8; ++k) {
        bit0_mask = _mm256_set1_epi32(7 - k);
        bit0 = _mm256_shuffle_epi8(buf, bit0_mask);
        bit0 = _mm256_srlv_epi32(bit0, rshift);
        ret = _mm256_and_si256(bit0, clr_mask);
        _mm256_storeu_si256((__m256i*)(data->array32_end + 8 * k), ret);
      }
      data->array32_end += 64;
      for (int k = 0; k < 8; ++k) {
        bit0_mask = _mm256_set1_epi32(15 - k);
        bit0 = _mm256_shuffle_epi8(buf, bit0_mask);
        bit0 = _mm256_srlv_epi32(bit0, rshift);
        ret = _mm256_and_si256(bit0, clr_mask);
        _mm256_storeu_si256((__m256i*)(data->array32_end + 8 * k), ret);
      }
      data->array32_end += 64;
      data->buffer_end += 2;
    }
    data->buffer_end = data->buffer;
    data->array32_end = reinterpret_cast<uint32_t*>(data->array);
    } else if (data->bit_width == 2) {
    __m256i bit0_mask, bit0, bit1_mask, bit1, ret;
    __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
        0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);
    for (int j = 0; j < data->num_values; j += 64) {
      __m256i buf = _mm256_loadu_si256((__m256i const*)data->buffer_end);
      __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
      buf = _mm256_inserti128_si256(buf, buf_hi, 1);
      for (int k = 0; k < 8; ++k) {
        bit0_mask = _mm256_set1_epi32(7 - k);
        bit0 = _mm256_shuffle_epi8(buf, bit0_mask);
        bit0 = _mm256_srlv_epi32(bit0, rshift);
        ret = _mm256_and_si256(bit0, clr_mask);
        bit1_mask = _mm256_set1_epi32(15 - k);
        bit1 = _mm256_shuffle_epi8(buf, bit1_mask);
        bit1 = _mm256_srlv_epi32(bit1, rshift);
        bit1 = _mm256_and_si256(bit1, clr_mask);
        bit1 = _mm256_slli_epi32(bit1, 1);
        ret = _mm256_or_si256(ret, bit1);
        _mm256_storeu_si256((__m256i*)(data->array32_end + 8 * k), ret);
      }
      data->array32_end += 64;
    }
    data->buffer_end = data->buffer;
    data->array32_end = reinterpret_cast<uint32_t*>(data->array);
    } else if (data->bit_width == 3) {
    __m256i bit0_mask, bit0, bit1_mask, bit1, bit2, ret;
    __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
        0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);
    for (int j = 0; j < data->num_values; j += 128) {
      __m256i buf = _mm256_loadu_si256((__m256i const*)data->buffer);
      __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
      buf = _mm256_inserti128_si256(buf, buf_hi, 1);
      __m256i buf23 = _mm256_loadu_si256((__m256i const*)(data->buffer + 2));
      buf_hi = _mm256_extracti128_si256(buf23, 0);
      buf23 = _mm256_inserti128_si256(buf23, buf_hi, 1);
      __m256i buf45 = _mm256_loadu_si256((__m256i const*)(data->buffer + 4));
      buf_hi = _mm256_extracti128_si256(buf45, 0);
      buf45 = _mm256_inserti128_si256(buf45, buf_hi, 1);
      for (int k = 0; k < 8; ++k) {
        bit0_mask = _mm256_set1_epi32(7 - k);
        bit0 = _mm256_shuffle_epi8(buf, bit0_mask);
        bit0 = _mm256_srlv_epi32(bit0, rshift);
        ret = _mm256_and_si256(bit0, clr_mask);
        bit1_mask = _mm256_set1_epi32(15 - k);
        bit1 = _mm256_shuffle_epi8(buf, bit1_mask);
        bit1 = _mm256_srlv_epi32(bit1, rshift);
        bit1 = _mm256_and_si256(bit1, clr_mask);
        bit1 = _mm256_slli_epi32(bit1, 1);
        ret = _mm256_or_si256(ret, bit1);
//        bit2_mask = _mm256_set1_epi32(7 - k);
        bit2 = _mm256_shuffle_epi8(buf23, bit0_mask);
        bit2 = _mm256_srlv_epi32(bit2, rshift);
        bit2 = _mm256_and_si256(bit2, clr_mask);
        bit2 = _mm256_slli_epi32(bit2, 2);
        ret = _mm256_or_si256(ret, bit2);
        _mm256_storeu_si256((__m256i*)(data->array32_end + 8 * k), ret);
      }
      data->array32_end += 64;
      for (int k = 0; k < 8; ++k) {
        bit0_mask = _mm256_set1_epi32(15 - k);
        bit0 = _mm256_shuffle_epi8(buf23, bit0_mask);
        bit0 = _mm256_srlv_epi32(bit0, rshift);
        ret = _mm256_and_si256(bit0, clr_mask);
        bit1_mask = _mm256_set1_epi32(7 - k);
        bit1 = _mm256_shuffle_epi8(buf45, bit1_mask);
        bit1 = _mm256_srlv_epi32(bit1, rshift);
        bit1 = _mm256_and_si256(bit1, clr_mask);
        bit1 = _mm256_slli_epi32(bit1, 1);
        ret = _mm256_or_si256(ret, bit1);
//        bit2_mask = _mm256_set1_epi32(15 - k);
        bit2 = _mm256_shuffle_epi8(buf45, bit0_mask);
        bit2 = _mm256_srlv_epi32(bit2, rshift);
        bit2 = _mm256_and_si256(bit2, clr_mask);
        bit2 = _mm256_slli_epi32(bit2, 2);
        ret = _mm256_or_si256(ret, bit2);
        _mm256_storeu_si256((__m256i*)(data->array32_end + 8 * k), ret);
      }
      data->array32_end += 64;
    }
    data->buffer_end = data->buffer;
    data->array32_end = reinterpret_cast<uint32_t*>(data->array);
    } else if (data->bit_width == 4) {
    __m256i bit0_mask, bit0, bit1_mask, bit1, bit2, bit3, ret;
    __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
        0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);
    for (int j = 0; j < data->num_values; j += 64) {
      __m256i buf = _mm256_loadu_si256((__m256i const*)data->buffer);
      __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
      buf = _mm256_inserti128_si256(buf, buf_hi, 1);
      __m256i buf23 = _mm256_loadu_si256((__m256i const*)(data->buffer + 2));
      buf_hi = _mm256_extracti128_si256(buf23, 0);
      buf23 = _mm256_inserti128_si256(buf23, buf_hi, 1);
      for (int k = 0; k < 8; ++k) {
        bit0_mask = _mm256_set1_epi32(7 - k);
        bit0 = _mm256_shuffle_epi8(buf, bit0_mask);
        bit0 = _mm256_srlv_epi32(bit0, rshift);
        ret = _mm256_and_si256(bit0, clr_mask);
        bit1_mask = _mm256_set1_epi32(15 - k);
        bit1 = _mm256_shuffle_epi8(buf, bit1_mask);
        bit1 = _mm256_srlv_epi32(bit1, rshift);
        bit1 = _mm256_and_si256(bit1, clr_mask);
        bit1 = _mm256_slli_epi32(bit1, 1);
        ret = _mm256_or_si256(ret, bit1);
        bit2 = _mm256_shuffle_epi8(buf23, bit0_mask);
        bit2 = _mm256_srlv_epi32(bit2, rshift);
        bit2 = _mm256_and_si256(bit2, clr_mask);
        bit2 = _mm256_slli_epi32(bit2, 2);
        ret = _mm256_or_si256(ret, bit2);
        bit3 = _mm256_shuffle_epi8(buf23, bit1_mask);
        bit3 = _mm256_srlv_epi32(bit3, rshift);
        bit3 = _mm256_and_si256(bit3, clr_mask);
        bit3 = _mm256_slli_epi32(bit3, 3);
        ret = _mm256_or_si256(ret, bit3);
        _mm256_storeu_si256((__m256i*)(data->array32_end + 8 * k), ret);
      }
      data->array32_end += 64;
    }
    data->buffer_end = data->buffer;
    data->array32_end = reinterpret_cast<uint32_t*>(data->array);
    } else if (data->bit_width == 5) {
    __m256i bit0_mask, bit0, bit1_mask, bit1, bit2, bit3, bit4, ret;
    __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
        0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);
    for (int j = 0; j < data->num_values; j += 128) {
      __m256i buf = _mm256_loadu_si256((__m256i const*)data->buffer);
      __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
      buf = _mm256_inserti128_si256(buf, buf_hi, 1);
      __m256i buf23 = _mm256_loadu_si256((__m256i const*)(data->buffer + 2));
      buf_hi = _mm256_extracti128_si256(buf23, 0);
      buf23 = _mm256_inserti128_si256(buf23, buf_hi, 1);
      __m256i buf45 = _mm256_loadu_si256((__m256i const*)(data->buffer + 4));
      buf_hi = _mm256_extracti128_si256(buf45, 0);
      buf45 = _mm256_inserti128_si256(buf45, buf_hi, 1);
      for (int k = 0; k < 8; ++k) {
        bit0_mask = _mm256_set1_epi32(7 - k);
        bit0 = _mm256_shuffle_epi8(buf, bit0_mask);
        bit0 = _mm256_srlv_epi32(bit0, rshift);
        ret = _mm256_and_si256(bit0, clr_mask);
        bit1_mask = _mm256_set1_epi32(15 - k);
        bit1 = _mm256_shuffle_epi8(buf, bit1_mask);
        bit1 = _mm256_srlv_epi32(bit1, rshift);
        bit1 = _mm256_and_si256(bit1, clr_mask);
        bit1 = _mm256_slli_epi32(bit1, 1);
        ret = _mm256_or_si256(ret, bit1);
        bit2 = _mm256_shuffle_epi8(buf23, bit0_mask);
        bit2 = _mm256_srlv_epi32(bit2, rshift);
        bit2 = _mm256_and_si256(bit2, clr_mask);
        bit2 = _mm256_slli_epi32(bit2, 2);
        ret = _mm256_or_si256(ret, bit2);
        bit3 = _mm256_shuffle_epi8(buf23, bit1_mask);
        bit3 = _mm256_srlv_epi32(bit3, rshift);
        bit3 = _mm256_and_si256(bit3, clr_mask);
        bit3 = _mm256_slli_epi32(bit3, 3);
        ret = _mm256_or_si256(ret, bit3);
        bit4 = _mm256_shuffle_epi8(buf45, bit0_mask);
        bit4 = _mm256_srlv_epi32(bit4, rshift);
        bit4 = _mm256_and_si256(bit4, clr_mask);
        bit4 = _mm256_slli_epi32(bit4, 4);
        ret = _mm256_or_si256(ret, bit4);
        _mm256_storeu_si256((__m256i*)(data->array32_end + 8 * k), ret);
      }
      data->array32_end += 64;

      __m256i buf67 = _mm256_loadu_si256((__m256i const*)(data->buffer + 6));
      buf_hi = _mm256_extracti128_si256(buf67, 0);
      buf67 = _mm256_inserti128_si256(buf67, buf_hi, 1);
      __m256i buf89 = _mm256_loadu_si256((__m256i const*)(data->buffer + 8));
      buf_hi = _mm256_extracti128_si256(buf89, 0);
      buf89 = _mm256_inserti128_si256(buf89, buf_hi, 1);
      for (int k = 0; k < 8; ++k) {
        bit0_mask = _mm256_set1_epi32(7 - k);
        bit1_mask = _mm256_set1_epi32(15 - k);
        bit0 = _mm256_shuffle_epi8(buf45, bit1_mask);
        bit0 = _mm256_srlv_epi32(bit0, rshift);
        ret = _mm256_and_si256(bit0, clr_mask);
        bit1 = _mm256_shuffle_epi8(buf67, bit0_mask);
        bit1 = _mm256_srlv_epi32(bit1, rshift);
        bit1 = _mm256_and_si256(bit1, clr_mask);
        bit1 = _mm256_slli_epi32(bit1, 1);
        ret = _mm256_or_si256(ret, bit1);
        bit2 = _mm256_shuffle_epi8(buf67, bit1_mask);
        bit2 = _mm256_srlv_epi32(bit2, rshift);
        bit2 = _mm256_and_si256(bit2, clr_mask);
        bit2 = _mm256_slli_epi32(bit2, 2);
        ret = _mm256_or_si256(ret, bit2);
        bit3 = _mm256_shuffle_epi8(buf89, bit0_mask);
        bit3 = _mm256_srlv_epi32(bit3, rshift);
        bit3 = _mm256_and_si256(bit3, clr_mask);
        bit3 = _mm256_slli_epi32(bit3, 3);
        ret = _mm256_or_si256(ret, bit3);
        bit4 = _mm256_shuffle_epi8(buf89, bit1_mask);
        bit4 = _mm256_srlv_epi32(bit4, rshift);
        bit4 = _mm256_and_si256(bit4, clr_mask);
        bit4 = _mm256_slli_epi32(bit4, 4);
        ret = _mm256_or_si256(ret, bit4);
        _mm256_storeu_si256((__m256i*)(data->array32_end + 8 * k), ret);
      }
      data->array32_end += 64;
    }
    data->buffer_end = data->buffer;
    data->array32_end = reinterpret_cast<uint32_t*>(data->array);
    } else if (data->bit_width == 6) {
    __m256i bit0_mask, bit0, bit1_mask, bit1, bit2, bit3, bit4, bit5, ret;
    __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
        0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);
    for (int j = 0; j < data->num_values; j += 64) {
      __m256i buf = _mm256_loadu_si256((__m256i const*)data->buffer);
      __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
      buf = _mm256_inserti128_si256(buf, buf_hi, 1);
      __m256i buf23 = _mm256_loadu_si256((__m256i const*)(data->buffer + 2));
      buf_hi = _mm256_extracti128_si256(buf23, 0);
      buf23 = _mm256_inserti128_si256(buf23, buf_hi, 1);
      __m256i buf45 = _mm256_loadu_si256((__m256i const*)(data->buffer + 4));
      buf_hi = _mm256_extracti128_si256(buf45, 0);
      buf45 = _mm256_inserti128_si256(buf45, buf_hi, 1);

      for (int k = 0; k < 8; ++k) {
        bit0_mask = _mm256_set1_epi32(7 - k);
        bit0 = _mm256_shuffle_epi8(buf, bit0_mask);
        bit0 = _mm256_srlv_epi32(bit0, rshift);
        ret = _mm256_and_si256(bit0, clr_mask);
        bit1_mask = _mm256_set1_epi32(15 - k);
        bit1 = _mm256_shuffle_epi8(buf, bit1_mask);
        bit1 = _mm256_srlv_epi32(bit1, rshift);
        bit1 = _mm256_and_si256(bit1, clr_mask);
        bit1 = _mm256_slli_epi32(bit1, 1);
        ret = _mm256_or_si256(ret, bit1);
        bit2 = _mm256_shuffle_epi8(buf23, bit0_mask);
        bit2 = _mm256_srlv_epi32(bit2, rshift);
        bit2 = _mm256_and_si256(bit2, clr_mask);
        bit2 = _mm256_slli_epi32(bit2, 2);
        ret = _mm256_or_si256(ret, bit2);
        bit3 = _mm256_shuffle_epi8(buf23, bit1_mask);
        bit3 = _mm256_srlv_epi32(bit3, rshift);
        bit3 = _mm256_and_si256(bit3, clr_mask);
        bit3 = _mm256_slli_epi32(bit3, 3);
        ret = _mm256_or_si256(ret, bit3);
        bit4 = _mm256_shuffle_epi8(buf45, bit0_mask);
        bit4 = _mm256_srlv_epi32(bit4, rshift);
        bit4 = _mm256_and_si256(bit4, clr_mask);
        bit4 = _mm256_slli_epi32(bit4, 4);
        ret = _mm256_or_si256(ret, bit4);
        bit5 = _mm256_shuffle_epi8(buf45, bit1_mask);
        bit5 = _mm256_srlv_epi32(bit5, rshift);
        bit5 = _mm256_and_si256(bit5, clr_mask);
        bit5 = _mm256_slli_epi32(bit5, 5);
        ret = _mm256_or_si256(ret, bit5);
        _mm256_storeu_si256((__m256i*)(data->array32_end + 8 * k), ret);
      }
      data->array32_end += 64;
    }
    data->buffer_end = data->buffer;
    data->array32_end = reinterpret_cast<uint32_t*>(data->array);
    } else if (data->bit_width == 10) {
    __m256i bit0_mask, bit0, bit1_mask, bit1, bit2, bit3, bit4, bit5, bit6, bit7, bit8, bit9, ret;
    __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
        0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);
    for (int j = 0; j < data->num_values; j += 64) {
      __m256i buf = _mm256_loadu_si256((__m256i const*)data->buffer);
      __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
      buf = _mm256_inserti128_si256(buf, buf_hi, 1);
      __m256i buf23 = _mm256_loadu_si256((__m256i const*)(data->buffer + 2));
      buf_hi = _mm256_extracti128_si256(buf23, 0);
      buf23 = _mm256_inserti128_si256(buf23, buf_hi, 1);
      __m256i buf45 = _mm256_loadu_si256((__m256i const*)(data->buffer + 4));
      buf_hi = _mm256_extracti128_si256(buf45, 0);
      buf45 = _mm256_inserti128_si256(buf45, buf_hi, 1);
      __m256i buf67 = _mm256_loadu_si256((__m256i const*)(data->buffer + 4));
      buf_hi = _mm256_extracti128_si256(buf67, 0);
      buf67 = _mm256_inserti128_si256(buf67, buf_hi, 1);
      __m256i buf89 = _mm256_loadu_si256((__m256i const*)(data->buffer + 4));
      buf_hi = _mm256_extracti128_si256(buf89, 0);
      buf89 = _mm256_inserti128_si256(buf89, buf_hi, 1);
      for (int k = 0; k < 8; ++k) {
        bit0_mask = _mm256_set1_epi32(7 - k);
        bit0 = _mm256_shuffle_epi8(buf, bit0_mask);
        bit0 = _mm256_srlv_epi32(bit0, rshift);
        ret = _mm256_and_si256(bit0, clr_mask);
        bit1_mask = _mm256_set1_epi32(15 - k);
        bit1 = _mm256_shuffle_epi8(buf, bit1_mask);
        bit1 = _mm256_srlv_epi32(bit1, rshift);
        bit1 = _mm256_and_si256(bit1, clr_mask);
        bit1 = _mm256_slli_epi32(bit1, 1);
        ret = _mm256_or_si256(ret, bit1);
        bit2 = _mm256_shuffle_epi8(buf23, bit0_mask);
        bit2 = _mm256_srlv_epi32(bit2, rshift);
        bit2 = _mm256_and_si256(bit2, clr_mask);
        bit2 = _mm256_slli_epi32(bit2, 2);
        ret = _mm256_or_si256(ret, bit2);
        bit3 = _mm256_shuffle_epi8(buf23, bit1_mask);
        bit3 = _mm256_srlv_epi32(bit3, rshift);
        bit3 = _mm256_and_si256(bit3, clr_mask);
        bit3 = _mm256_slli_epi32(bit3, 3);
        ret = _mm256_or_si256(ret, bit3);
        bit4 = _mm256_shuffle_epi8(buf45, bit0_mask);
        bit4 = _mm256_srlv_epi32(bit4, rshift);
        bit4 = _mm256_and_si256(bit4, clr_mask);
        bit4 = _mm256_slli_epi32(bit4, 4);
        ret = _mm256_or_si256(ret, bit4);
        bit5 = _mm256_shuffle_epi8(buf45, bit1_mask);
        bit5 = _mm256_srlv_epi32(bit5, rshift);
        bit5 = _mm256_and_si256(bit5, clr_mask);
        bit5 = _mm256_slli_epi32(bit5, 5);
        ret = _mm256_or_si256(ret, bit5);
        bit6 = _mm256_shuffle_epi8(buf67, bit0_mask);
        bit6 = _mm256_srlv_epi32(bit6, rshift);
        bit6 = _mm256_and_si256(bit6, clr_mask);
        bit6 = _mm256_slli_epi32(bit6, 6);
        ret = _mm256_or_si256(ret, bit6);
        bit7 = _mm256_shuffle_epi8(buf67, bit1_mask);
        bit7 = _mm256_srlv_epi32(bit7, rshift);
        bit7 = _mm256_and_si256(bit7, clr_mask);
        bit7 = _mm256_slli_epi32(bit7, 7);
        ret = _mm256_or_si256(ret, bit7);
        bit8 = _mm256_shuffle_epi8(buf89, bit0_mask);
        bit8 = _mm256_srlv_epi32(bit8, rshift);
        bit8 = _mm256_and_si256(bit8, clr_mask);
        bit8 = _mm256_slli_epi32(bit8, 8);
        ret = _mm256_or_si256(ret, bit8);
        bit9 = _mm256_shuffle_epi8(buf89, bit1_mask);
        bit9 = _mm256_srlv_epi32(bit9, rshift);
        bit9 = _mm256_and_si256(bit9, clr_mask);
        bit9 = _mm256_slli_epi32(bit9, 9);
        ret = _mm256_or_si256(ret, bit9);
        _mm256_storeu_si256((__m256i*)(data->array32_end + 8 * k), ret);
      }
      data->array32_end += 64;
    }
    data->buffer_end = data->buffer;
    data->array32_end = reinterpret_cast<uint32_t*>(data->array);
    } else if (data->bit_width == 32) {
    __m256i bit0_mask, bit, bit1_mask, ret;
    __m256i rshift = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256i clr_mask = _mm256_set_epi32(0x00000001, 0x00000001, 0x00000001,
        0x00000001, 0x00000001, 0x00000001, 0x00000001, 0x00000001);
    for (int j = 0; j < data->num_values; j += 64) {
      __m256i buf[16];
      for (int i = 0; i < 16; ++i) {
        buf[i] = _mm256_loadu_si256((__m256i const*)(data->buffer + i * 2));
        __m128i buf_hi = _mm256_extracti128_si256(buf[i], 0);
        buf[i] = _mm256_inserti128_si256(buf[i], buf_hi, 1);
      }
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
        for (int l = 1; l < 16; ++l) {
          bit = _mm256_shuffle_epi8(buf[l], bit0_mask);
          bit = _mm256_srlv_epi32(bit, rshift);
          bit = _mm256_and_si256(bit, clr_mask);
          bit = _mm256_slli_epi32(bit, 2 * l);
          ret = _mm256_or_si256(ret, bit);
          bit = _mm256_shuffle_epi8(buf[l], bit1_mask);
          bit = _mm256_srlv_epi32(bit, rshift);
          bit = _mm256_and_si256(bit, clr_mask);
          bit = _mm256_slli_epi32(bit, 2 * l + 1);
        }

        _mm256_storeu_si256((__m256i*)(data->array32_end + 8 * k), ret);
      }
      data->array32_end += 64;
    }
    data->buffer_end = data->buffer;
    data->array32_end = reinterpret_cast<uint32_t*>(data->array);
    } else {
    // Unroll this to focus more on Put performance.
    for (int j = 0; j < data->num_values; j += 64) {
      for (int k = 0; k < data->bit_width; ++k) {
        for (int l = 0; l < 64; ++l) {
          data->array_end[l] |= ((data->buffer_end[k] >> l) & 0x01ULL) << k;
        }
      }
      data->buffer_end += data->bit_width;
      data->array_end += 64;
    }
    data->buffer_end = data->buffer;
    data->array_end = data->array;
    }
  }
}

void TestFleWithAVXV4Decode(int batch_size, void* d) {
  TestData* data = reinterpret_cast<TestData*>(d);
  for (int i = 0; i < batch_size; ++i) {
    // Unroll this to focus more on Put performance.
    if (data->bit_width == 1) {
      __m256i clr_mask = _mm256_set1_epi64x(0x0102040810204080);
      __m256i bit0 = _mm256_set1_epi8(0x01);
      __m256i shf0_mask = _mm256_setr_epi64x(0x0707070707070707, 0x0606060606060606,
          0x0505050505050505, 0x0404040404040404);
      __m256i shf1_mask = _mm256_setr_epi64x(0x0303030303030303, 0x0202020202020202,
          0x0101010101010101, 0x0000000000000000);
      __m256i shf2_mask = _mm256_setr_epi64x(0x0f0f0f0f0f0f0f0f, 0x0e0e0e0e0e0e0e0e,
          0x0d0d0d0d0d0d0d0d, 0x0c0c0c0c0c0c0c0c);
      __m256i shf3_mask = _mm256_setr_epi64x(0x0b0b0b0b0b0b0b0b, 0x0a0a0a0a0a0a0a0a,
          0x0909090909090909, 0x0808080808080808);
      for (int j = 0; j < data->num_values; j += 128) {
        __m256i buf = _mm256_loadu_si256((__m256i const*)data->buffer_end);
        __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
        buf = _mm256_inserti128_si256(buf, buf_hi, 1);
        __m256i bit, ret_mask, ret;

        bit = _mm256_shuffle_epi8(buf, shf0_mask);
        bit = _mm256_and_si256(bit, clr_mask);
        ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
        ret = _mm256_and_si256(bit0, ret_mask);
        _mm256_storeu_si256((__m256i*)(data->array8_end), ret);
        data->array8_end += 32;
        bit = _mm256_shuffle_epi8(buf, shf1_mask);
        bit = _mm256_and_si256(bit, clr_mask);
        ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
        ret = _mm256_and_si256(bit0, ret_mask);
        _mm256_storeu_si256((__m256i*)(data->array8_end), ret);
        data->array8_end += 32;
        bit = _mm256_shuffle_epi8(buf, shf2_mask);
        bit = _mm256_and_si256(bit, clr_mask);
        ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
        ret = _mm256_and_si256(bit0, ret_mask);
        _mm256_storeu_si256((__m256i*)(data->array8_end), ret);
        data->array8_end += 32;
        bit = _mm256_shuffle_epi8(buf, shf3_mask);
        bit = _mm256_and_si256(bit, clr_mask);
        ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
        ret = _mm256_and_si256(bit0, ret_mask);
        _mm256_storeu_si256((__m256i*)(data->array8_end), ret);
        data->array8_end += 32;
        data->buffer_end += 2;
      }
      data->buffer_end = data->buffer;
      data->array8_end = reinterpret_cast<uint8_t*>(data->array);
    } else if (data->bit_width == 2) {
      __m256i clr_mask = _mm256_set1_epi64x(0x0102040810204080);
      __m256i bit0 = _mm256_set1_epi8(0x01);
      __m256i bit1 = _mm256_set1_epi8(0x02);
      __m256i shf0_mask = _mm256_setr_epi64x(0x0707070707070707, 0x0606060606060606,
          0x0505050505050505, 0x0404040404040404);
      __m256i shf1_mask = _mm256_setr_epi64x(0x0303030303030303, 0x0202020202020202,
          0x0101010101010101, 0x0000000000000000);
      __m256i shf2_mask = _mm256_setr_epi64x(0x0f0f0f0f0f0f0f0f, 0x0e0e0e0e0e0e0e0e,
          0x0d0d0d0d0d0d0d0d, 0x0c0c0c0c0c0c0c0c);
      __m256i shf3_mask = _mm256_setr_epi64x(0x0b0b0b0b0b0b0b0b, 0x0a0a0a0a0a0a0a0a,
          0x0909090909090909, 0x0808080808080808);
      for (int j = 0; j < data->num_values; j += 64) {
        __m256i buf = _mm256_loadu_si256((__m256i const*)data->buffer_end);
        __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
        buf = _mm256_inserti128_si256(buf, buf_hi, 1);
        __m256i bit, ret_mask, tmp_ret, ret;

        bit = _mm256_shuffle_epi8(buf, shf0_mask);
        bit = _mm256_and_si256(bit, clr_mask);
        ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
        ret = _mm256_and_si256(bit0, ret_mask);
        bit = _mm256_shuffle_epi8(buf, shf2_mask);
        bit = _mm256_and_si256(bit, clr_mask);
        ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
        tmp_ret = _mm256_and_si256(bit1, ret_mask);
        ret = _mm256_or_si256(ret, tmp_ret);
        _mm256_storeu_si256((__m256i*)(data->array8_end), ret);
        data->array8_end += 32;

        bit = _mm256_shuffle_epi8(buf, shf1_mask);
        bit = _mm256_and_si256(bit, clr_mask);
        ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
        ret = _mm256_and_si256(bit0, ret_mask);
        bit = _mm256_shuffle_epi8(buf, shf3_mask);
        bit = _mm256_and_si256(bit, clr_mask);
        ret_mask = _mm256_cmpeq_epi8(bit, clr_mask);
        tmp_ret = _mm256_and_si256(bit1, ret_mask);
        ret = _mm256_or_si256(ret, tmp_ret);
        _mm256_storeu_si256((__m256i*)(data->array8_end), ret);
        data->array8_end += 32;
        data->buffer_end += 2;
      }
      data->buffer_end = data->buffer;
      data->array8_end = reinterpret_cast<uint8_t*>(data->array);
    } else if (data->bit_width == 9) {
      __m256i clr_mask = _mm256_set1_epi64x(0x0102040810204080);
      __m256i shf0_mask = _mm256_setr_epi64x(0x0707070707070707, 0x0606060606060606,
          0x0505050505050505, 0x0404040404040404);
      __m256i shf1_mask = _mm256_setr_epi64x(0x0303030303030303, 0x0202020202020202,
          0x0101010101010101, 0x0000000000000000);
      __m256i shf2_mask = _mm256_setr_epi64x(0x0f0f0f0f0f0f0f0f, 0x0e0e0e0e0e0e0e0e,
          0x0d0d0d0d0d0d0d0d, 0x0c0c0c0c0c0c0c0c);
      __m256i shf3_mask = _mm256_setr_epi64x(0x0b0b0b0b0b0b0b0b, 0x0a0a0a0a0a0a0a0a,
          0x0909090909090909, 0x0808080808080808);

      __m256i shf_ret_mask = _mm256_setr_epi64x(0x8003800280018000, 0x8007800680058004,
          0x800b800a80098008, 0x800f800e800d800c);
      __m256i shf_ret1_mask = _mm256_setr_epi64x(0x0380028001800080, 0x0780068005800480,
          0x0b800a8009800880, 0x0f800e800d800c80);
     __m256i bitv[8];
    for (int l = 0; l < 8; ++l) {
      bitv[l] = _mm256_set1_epi8(0x01 << l);
    }


      for (int j = 0; j < data->num_values; j += 64) {
        __m256i buf = _mm256_loadu_si256((__m256i const*)data->buffer_end);
        __m128i buf_hi = _mm256_extracti128_si256(buf, 0);
        buf = _mm256_inserti128_si256(buf, buf_hi, 1);
  __m256i buf1 = _mm256_loadu_si256((__m256i const*)(data->buffer_end + 2));
  buf_hi = _mm256_extracti128_si256(buf1, 0);
  buf1 = _mm256_inserti128_si256(buf1, buf_hi, 1);
  __m256i buf2 = _mm256_loadu_si256((__m256i const*)(data->buffer_end + 4));
  buf_hi = _mm256_extracti128_si256(buf2, 0);
  buf2 = _mm256_inserti128_si256(buf2, buf_hi, 1);
  __m256i buf3 = _mm256_loadu_si256((__m256i const*)(data->buffer_end + 6));
  buf_hi = _mm256_extracti128_si256(buf3, 0);
  buf3 = _mm256_inserti128_si256(buf3, buf_hi, 1);
  __m256i buf4 = _mm256_loadu_si256((__m256i const*)(data->buffer_end + 8));
  buf_hi = _mm256_extracti128_si256(buf4, 0);
  buf4 = _mm256_inserti128_si256(buf4, buf_hi, 1);

  __m256i bit, ret_mask, tmp_ret, ret, ret1, half_ret, half_ret1, shf_ret;
  __m128i ret_hi, ret1_hi;


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
  ret_hi = _mm256_extracti128_si256(ret, 0);
  half_ret = _mm256_inserti128_si256(ret, ret_hi, 1);
  shf_ret = _mm256_shuffle_epi8(half_ret, shf_ret_mask);

  ret1_hi = _mm256_extracti128_si256(ret1, 0);
  half_ret1 = _mm256_inserti128_si256(ret1, ret1_hi, 1);
  tmp_ret = _mm256_shuffle_epi8(half_ret1, shf_ret1_mask);
  shf_ret =  _mm256_or_si256(shf_ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(data->array16_end), shf_ret);

  data->array16_end += 16;

  ret_hi = _mm256_extracti128_si256(ret, 1);
  half_ret = _mm256_inserti128_si256(ret, ret_hi, 0);
  shf_ret = _mm256_shuffle_epi8(half_ret, shf_ret_mask);

  ret1_hi = _mm256_extracti128_si256(ret1, 1);
  half_ret1 = _mm256_inserti128_si256(ret1, ret1_hi, 0);
  tmp_ret = _mm256_shuffle_epi8(half_ret1, shf_ret1_mask);
  shf_ret =  _mm256_or_si256(shf_ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(data->array16_end), shf_ret);

  data->array16_end += 16;


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
  ret_hi = _mm256_extracti128_si256(ret, 0);
  half_ret = _mm256_inserti128_si256(ret, ret_hi, 1);
  shf_ret = _mm256_shuffle_epi8(half_ret, shf_ret_mask);

  ret1_hi = _mm256_extracti128_si256(ret1, 0);
  half_ret1 = _mm256_inserti128_si256(ret1, ret1_hi, 1);
  tmp_ret = _mm256_shuffle_epi8(half_ret1, shf_ret1_mask);
  shf_ret =  _mm256_or_si256(shf_ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(data->array16_end), shf_ret);

  data->array16_end += 16;

  ret_hi = _mm256_extracti128_si256(ret, 1);
  half_ret = _mm256_inserti128_si256(ret, ret_hi, 0);
  shf_ret = _mm256_shuffle_epi8(half_ret, shf_ret_mask);

  ret1_hi = _mm256_extracti128_si256(ret1, 1);
  half_ret1 = _mm256_inserti128_si256(ret1, ret1_hi, 0);
  tmp_ret = _mm256_shuffle_epi8(half_ret1, shf_ret1_mask);
  shf_ret =  _mm256_or_si256(shf_ret, tmp_ret);
  _mm256_storeu_si256((__m256i*)(data->array16_end), shf_ret);

        data->array16_end += 16;
        data->buffer_end += 9;
      }
      data->buffer_end = data->buffer;
      data->array16_end = reinterpret_cast<uint16_t*>(data->array);
    } else {
    // Unroll this to focus more on Put performance.
    for (int j = 0; j < data->num_values; j += 64) {
      for (int k = 0; k < data->bit_width; ++k) {
        for (int l = 0; l < 64; ++l) {
          data->array_end[l] |= ((data->buffer_end[k] >> l) & 0x01ULL) << k;
        }
      }
      data->buffer_end += data->bit_width;
      data->array_end += 64;
    }
    data->buffer_end = data->buffer;
    data->array_end = data->array;
    }
  }
}


int main(int argc, char** argv) {
  CpuInfo::Init();

  MemTracker tracker;
  MemPool pool(&tracker);

  int num_values = 4096;
  int max_bits = 32;

//  GenEncoder<((0x01 << 8) - 1)> genEncoder;

  srand((int)time(0));
  Benchmark encode_suite("encode");
  TestData data[max_bits];
  for (int i = 0; i < max_bits; ++i) {
    data[i].buffer = new uint64_t[BUFFER_LEN];
    data[i].buffer_end = data[i].buffer;
    data[i].buffer_len = BUFFER_LEN;
    data[i].buffer32_end = reinterpret_cast<uint32_t*>(data[i].buffer);
    data[i].array = new uint64_t[num_values];
    data[i].array_end = data[i].array;
    data[i].array32_end = reinterpret_cast<uint32_t*>(data[i].array);
    data[i].array16_end = reinterpret_cast<uint16_t*>(data[i].array);
    data[i].array8_end = reinterpret_cast<uint8_t*>(data[i].array);
    data[i].num_values = num_values;
    data[i].bit_width = i + 1;
    data[i].max_value = 1 << i;
    data[i].pool = &pool;

    uint64_t max_value = 0x01ULL << data[i].bit_width;
    for (int j = 0; j < num_values; ++j) {
      data[i].array[j] = rand() % (max_value);
    }

    stringstream suffix;
    suffix << " " << (i+1) << "-Bit";

    stringstream name;
//    name << "\"Native" << suffix.str() << "\"";
//    int baseline =
//        encode_suite.AddBenchmark(name.str(), TestNativeEncode, &data[i], -1);

//    name.str("");
//    name << "\"Fle" << suffix.str() << "\"";
//    encode_suite.AddBenchmark(name.str(), TestFleEncode, &data[i], baseline);

    name.str("");
    name << "\"Fle With AVX" << suffix.str() << "\"";
//    int baseline = encode_suite.AddBenchmark(name.str(), TestFleWithAVXEncode, &data[i], baseline);
    int baseline = encode_suite.AddBenchmark(name.str(), TestFleWithAVXV2Encode, &data[i], -1);

    name.str("");
    name << "\"Fle With AVX V2" << suffix.str() << "\"";
    encode_suite.AddBenchmark(name.str(), TestFleWithAVXV2Encode, &data[i], baseline);
  }
  cout << encode_suite.Measure() << endl;

/*
  Benchmark decode_suite("decode");
  for (int i = 0; i < max_bits; ++i) {
    stringstream suffix;
    suffix << " " << (i+1) << "-Bit";

    stringstream name;
    name << "\"Native" << suffix.str() << "\"";
    int baseline =
        decode_suite.AddBenchmark(name.str(), TestNativeDecode, &data[i], -1);

    name.str("");
    name << "\"Fle" << suffix.str() << "\"";
    decode_suite.AddBenchmark(name.str(), TestFleDecode, &data[i], baseline);

    name.str("");
    name << "\"Fle With AVX" << suffix.str() << "\"";
    decode_suite.AddBenchmark(name.str(), TestFleWithAVXDecode, &data[i], baseline);

    name.str("");
    name << "\"Fle With AVX V3" << suffix.str() << "\"";
    decode_suite.AddBenchmark(name.str(), TestFleWithAVXV3Decode, &data[i], baseline);

    name.str("");
    name << "\"Fle With AVX V4" << suffix.str() << "\"";
    decode_suite.AddBenchmark(name.str(), TestFleWithAVXV4Decode, &data[i], baseline);

  }
  cout << decode_suite.Measure() << endl;
*/

  return 0;
}

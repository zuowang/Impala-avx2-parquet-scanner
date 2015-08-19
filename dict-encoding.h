// Copyright 2012 Cloudera Inc.
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

#ifndef IMPALA_UTIL_DICT_ENCODING_H
#define IMPALA_UTIL_DICT_ENCODING_H

#include <map>

#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/unordered_map.hpp>

#include "exec/parquet-common.h"
#include "runtime/mem-pool.h"
#include "runtime/string-value.h"
#include "util/bit-util.h"
#include "util/fle-encoding.h"
#include "util/rle-encoding.h"
#include "util/runtime-profile.h"

using namespace impala_udf;

namespace impala {

// See the dictionary encoding section of https://github.com/Parquet/parquet-format.
// This class supports dictionary encoding of all Impala types.
// The encoding supports streaming encoding. Values are encoded as they are added while
// the dictionary is being constructed. At any time, the buffered values can be
// written out with the current dictionary size. More values can then be added to
// the encoder, including new dictionary entries.
// TODO: if the dictionary was made to be ordered, the dictionary would compress better.
// Add this to the spec as future improvement.

// Base class for encoders. This is convenient so users can have a type that
// abstracts over the actual dictionary type.
// Note: it does not provide a virtual Put(). Users are expected to know the subclass
// type when using Put().
// TODO: once we can easily remove virtual calls with codegen, this interface can
// rely less on templating and be easier to follow. The type should be passed in
// as an argument rather than template argument.
class DictEncoderBase {
 public:
  virtual ~DictEncoderBase() {
    DCHECK(buffered_indices_.empty());
  }

  // Writes out the encoded dictionary to buffer. buffer must be preallocated to
  // dict_encoded_size() bytes.
  virtual void WriteDict(uint8_t* buffer) = 0;

  // The number of entries in the dictionary.
  virtual int num_entries() const = 0;

  // Clears all the indices (but leaves the dictionary).
  void ClearIndices() { buffered_indices_.clear(); }

  // Returns a conservative estimate of the number of bytes needed to encode the buffered
  // indices. Used to size the buffer passed to WriteData().
//  int EstimatedDataEncodedSize() {
//    return 1 + RleEncoder::MaxBufferSize(bit_width(), buffered_indices_.size());
//  }

  // The minimum bit width required to encode the currently buffered indices.
  int bit_width() const {
    if (UNLIKELY(num_entries() == 0)) return 0;
    if (UNLIKELY(num_entries() == 1)) return 1;
    return BitUtil::Log2(num_entries());
  }

  // Writes out any buffered indices to buffer preceded by the bit width of this data.
  // Returns the number of bytes written.
  // If the supplied buffer is not big enough, returns -1.
  // buffer must be preallocated with buffer_len bytes. Use EstimatedDataEncodedSize()
  // to size buffer.
  int WriteData(uint8_t* buffer, int buffer_len);

  int WriteData(uint8_t* buffer, int buffer_len, vector<int>& node_indices);

  int dict_encoded_size() { return dict_encoded_size_; }

 protected:
  DictEncoderBase(MemPool* pool)
    : dict_encoded_size_(0), pool_(pool) {
    count0 = 0;
  }

  // Indices that have not yet be written out by WriteData().
  std::vector<int> buffered_indices_;

  std::vector<int> to_sorted_indice_;

  // The number of bytes needed to encode the dictionary.
  int dict_encoded_size_;

  // Pool to store StringValue data. Not owned.
  MemPool* pool_;

  int count0;
};

template<typename T>
class DictEncoder : public DictEncoderBase {
 public:
  DictEncoder(MemPool* pool, int encoded_value_size) :
      DictEncoderBase(pool), buckets_(HASH_TABLE_SIZE, Node::INVALID_INDEX),
      encoded_value_size_(encoded_value_size) { }

  // Encode value. Returns the number of bytes added to the dictionary page length
  // (will be 0 if this value is already in the dictionary) or -1 if the dictionary is
  // full (in which case the caller should give up on dictionary encoding). Note that
  // this does not actually write any data, just buffers the value's index to be
  // written later.
  int Put(const T& value);

  int Put(const T& value, vector<int>& node_indices);

  virtual void WriteDict(uint8_t* buffer);

  virtual int num_entries() const { return nodes_.size(); }

 private:
  // Size of the table. Must be a power of 2.
  enum { HASH_TABLE_SIZE = 1 << 16 };

  // Dictates an upper bound on the capacity of the hash table.
  typedef uint16_t NodeIndex;

  // Hash table mapping value to dictionary index (i.e. the number used to encode this
  // value in the data). Each table entry is a index into the nodes_ vector (giving the
  // first node of a chain for this bucket) or Node::INVALID_INDEX for an empty bucket.
  std::vector<NodeIndex> buckets_;

  // Node in the chained hash table.
  struct Node {
    Node(const T& v, const NodeIndex& n) : value(v), next(n) { }

    // The dictionary value.
    T value;

    // Index into nodes_ for the next Node in the chain. INVALID_INDEX indicates end.
    NodeIndex next;

    // The maximum number of values in the dictionary.  Chosen to be around 60% of
    // HASH_TABLE_SIZE to limit the expected length of the chains.
    enum { INVALID_INDEX = 40000 };
  };

  // The nodes of the hash table. Ordered by dictionary index (and so also represents
  // the reverse mapping from encoded index to value).
  std::vector<Node> nodes_;

  // Size of each encoded dictionary value. -1 for variable-length types.
  int encoded_value_size_;

  // Hash function for mapping a value to a bucket.
  inline uint32_t Hash(const T& value) const;

  // Adds value to the hash table and updates dict_encoded_size_. Returns the
  // number of bytes added to dict_encoded_size_.
  // bucket gives a pointer to the location (i.e. chain) to add the value
  // so that the hash for value doesn't need to be recomputed.
  int AddToTable(const T& value, NodeIndex* bucket);

  static bool NodeValLess(const Node& i, const Node& j);
};

// Decoder class for dictionary encoded data. This class does not allocate any
// buffers. The input buffers (dictionary buffer and RLE buffer) must be maintained
// by the caller and valid as long as this object is.
class DictDecoderBase {
 public:
  // The rle encoded indices into the dictionary.
  void SetData(uint8_t* buffer, int buffer_len) {
    DCHECK_GT(buffer_len, 0);
    uint8_t bit_width = *buffer;
    DCHECK_GE(bit_width, 0);
    ++buffer;
    --buffer_len;
    data_decoder_.reset(new FleDecoder(buffer, buffer_len, bit_width));
  }

  virtual ~DictDecoderBase() {}

  virtual int num_entries() const = 0;

 protected:
  boost::scoped_ptr<FleDecoder> data_decoder_;
};

template<typename T>
class DictDecoder : public DictDecoderBase {
 public:
  // The input buffer containing the dictionary.  'dict_len' is the byte length
  // of dict_buffer.
  // For string data, the decoder returns StringValues with data directly from
  // dict_buffer (i.e. no copies).
  // fixed_len_size is the size that must be passed to decode fixed-length
  // dictionary values (values stored using FIXED_LEN_BYTE_ARRAY).
  DictDecoder(uint8_t* dict_buffer, int dict_len, int fixed_len_size);

  virtual int num_entries() const { return dict_.size(); }

  // Returns the next value.  Returns false if the data is invalid.
  // For StringValues, this does not make a copy of the data.  Instead,
  // the string data is from the dictionary buffer passed into the c'tor.
  bool GetValue(T* value);

  bool GetValue(T* value, int skip_rows);

  bool SkipValue(int skip_rows);

  inline void Eq(int64_t num_rows, dynamic_bitset<>& skip_bitset, T& val);
  inline void Gt(int64_t num_rows, dynamic_bitset<>& skip_bitset, T& val);
  inline void Lt(int64_t num_rows, dynamic_bitset<>& skip_bitset, T& val);
  inline void Ge(int64_t num_rows, dynamic_bitset<>& skip_bitset, T& val);
  inline void Le(int64_t num_rows, dynamic_bitset<>& skip_bitset, T& val);
  inline void In(int64_t num_rows, dynamic_bitset<>& skip_bitset, vector<T>& vals);
 private:
  std::vector<T> dict_;
};

template<typename T>
inline int DictEncoder<T>::Put(const T& value, vector<int>& node_indices) {
  NodeIndex* bucket = &buckets_[Hash(value) & (HASH_TABLE_SIZE - 1)];
  NodeIndex i = *bucket;
  // Look for the value in the dictionary.
  while (i != Node::INVALID_INDEX) {
    const Node* n = &nodes_[i];
    if (LIKELY(n->value == value)) {
      // Value already in dictionary.
      node_indices.push_back(i);
      return 0;
    }
    i = n->next;
  }
  // Value not found. Add it to the dictionary if there's space.
  i = nodes_.size();
  if (UNLIKELY(i >= Node::INVALID_INDEX)) return -1;
  node_indices.push_back(i);
  return AddToTable(value, bucket);
}

template<typename T>
inline int DictEncoder<T>::Put(const T& value) {
  NodeIndex* bucket = &buckets_[Hash(value) & (HASH_TABLE_SIZE - 1)];
  NodeIndex i = *bucket;
  // Look for the value in the dictionary.
  while (i != Node::INVALID_INDEX) {
    const Node* n = &nodes_[i];
    if (LIKELY(n->value == value)) {
      // Value already in dictionary.
      buffered_indices_.push_back(i);
      return 0;
    }
    i = n->next;
  }
  // Value not found. Add it to the dictionary if there's space.
  i = nodes_.size();
  if (UNLIKELY(i >= Node::INVALID_INDEX)) return -1;
  buffered_indices_.push_back(i);
  return AddToTable(value, bucket);
}

template<typename T>
inline uint32_t DictEncoder<T>::Hash(const T& value) const {
  return HashUtil::Hash(&value, sizeof(value), 0);
}

template<>
inline uint32_t DictEncoder<StringValue>::Hash(const StringValue& value) const {
  return HashUtil::Hash(value.ptr, value.len, 0);
}

template<typename T>
inline int DictEncoder<T>::AddToTable(const T& value, NodeIndex* bucket) {
  DCHECK_GT(encoded_value_size_, 0);
  // Prepend the new node to this bucket's chain.
  nodes_.push_back(Node(value, *bucket));
  *bucket = nodes_.size() - 1;
  dict_encoded_size_ += encoded_value_size_;
  return encoded_value_size_;
}

template<>
inline int DictEncoder<StringValue>::AddToTable(const StringValue& value,
    NodeIndex* bucket) {
  char* ptr_copy = reinterpret_cast<char*>(pool_->Allocate(value.len));
  memcpy(ptr_copy, value.ptr, value.len);
  StringValue sv(ptr_copy, value.len);
  // Prepend the new node to this bucket's chain.
  nodes_.push_back(Node(sv, *bucket));
  *bucket = nodes_.size() - 1;
  int bytes_added = ParquetPlainEncoder::ByteSize(sv);
  dict_encoded_size_ += bytes_added;
  return bytes_added;
}

template<typename T>
inline bool DictDecoder<T>::GetValue(T* value) {
  DCHECK(data_decoder_.get() != NULL);
  int index;
  bool result = data_decoder_->Get(&index);
  if (!result) return false;
  if (index >= dict_.size()) return false;
  *value = dict_[index];
  return true;
}

template<typename T>
inline bool DictDecoder<T>::GetValue(T* value, int skip_rows) {
  DCHECK(data_decoder_.get() != NULL);
  int index;
  bool result = data_decoder_->Get(&index, skip_rows);
  if (!result) return false;
  if (index >= dict_.size()) return false;
  *value = dict_[index];
  return true;
}

template<typename T>
inline bool DictDecoder<T>::SkipValue(int skip_rows) {
  DCHECK(data_decoder_.get() != NULL);
  return data_decoder_->Skip(skip_rows);
}

template<>
inline bool DictDecoder<Decimal16Value>::GetValue(Decimal16Value* value) {
  DCHECK(data_decoder_.get() != NULL);
  int index;
  bool result = data_decoder_->Get(&index);
  if (!result) return false;
  if (index >= dict_.size()) return false;
  // Workaround for IMPALA-959. Use memcpy instead of '=' so addresses
  // do not need to be 16 byte aligned.
  uint8_t* addr = reinterpret_cast<uint8_t*>(&dict_[0]);
  addr = addr + index * sizeof(*value);
  memcpy(value, addr, sizeof(*value));
  return true;
}

template<>
inline bool DictDecoder<Decimal16Value>::GetValue(Decimal16Value* value,
    int skip_rows) {
  DCHECK(data_decoder_.get() != NULL);
  int index;
  bool result = data_decoder_->Get(&index, skip_rows);
  if (!result) return false;
  if (index >= dict_.size()) return false;
  // Workaround for IMPALA-959. Use memcpy instead of '=' so addresses
  // do not need to be 16 byte aligned.
  uint8_t* addr = reinterpret_cast<uint8_t*>(&dict_[0]);
  addr = addr + index * sizeof(*value);
  memcpy(value, addr, sizeof(*value));
  return true;
}

template<typename T>
inline bool DictEncoder<T>::NodeValLess(const Node& i, const Node& j) {
  return i.value < j.value;
}

template <>
inline bool DictEncoder<StringVal>::NodeValLess(const Node& i, const Node& j) {
  int n = min(i.value.len, j.value.len);
  int result = memcmp(&i.value.ptr[0], &j.value.ptr[0], n);
  if (result == 0) return i.value.len < j.value.len;
  return result < 0;
}

template <>
inline bool DictEncoder<DecimalVal>::NodeValLess(const Node& i, const Node& j) {
  return i.value.val16 < j.value.val16;
}

template <>
inline bool DictEncoder<TimestampVal>::NodeValLess(const Node& i, const Node& j) {
  if (i.value.date == j.value.date) return i.value.time_of_day < j.value.time_of_day;
  else return i.value.date < j.value.date;
}

template<typename T>
inline void DictEncoder<T>::WriteDict(uint8_t* buffer) {
  for (int i = 0; i < nodes_.size(); ++i) {
    nodes_[i].next = i;
  }
  sort(nodes_.begin(), nodes_.end(), DictEncoder<T>::NodeValLess);
  to_sorted_indice_.resize(nodes_.size());
  for (int i = 0; i < nodes_.size(); ++i) {
    to_sorted_indice_[nodes_[i].next] = i;
  }
  BOOST_FOREACH(const Node& node, nodes_) {
    buffer += ParquetPlainEncoder::Encode(buffer, encoded_value_size_, node.value);
  }
}

inline int DictEncoderBase::WriteData(uint8_t* buffer, int buffer_len) {
  for (int i = 0; i < buffered_indices_.size(); ++i) {
    buffered_indices_[i] = to_sorted_indice_[buffered_indices_[i]];
  }
  // Write bit width in first byte
  *buffer = bit_width();
  ++buffer;
  --buffer_len;

  FleEncoder encoder(buffer, buffer_len, bit_width());
  BOOST_FOREACH(int index, buffered_indices_) {
    if (!encoder.Put(index)) return -1;
  }
  encoder.Flush();
  return 1 + encoder.len();
}

inline int DictEncoderBase::WriteData(uint8_t* buffer, int buffer_len, vector<int>& node_indices) {
  // Write bit width in first byte
  //*buffer = bit_width();
  //++buffer;
  --buffer_len;

  int max_indice = -1;
  FleEncoder encoder(buffer + 1, buffer_len, bit_width());
  BOOST_FOREACH(int index, node_indices) {
    if (!encoder.Put(to_sorted_indice_[index])) return -1;
    if (index > max_indice) max_indice = index;
  }
  // Write bit width in first byte
  if (UNLIKELY(max_indice == -1)) {
    *buffer = 0;
  } else if (UNLIKELY(max_indice == 0)) {
    *buffer = 1;
  } else {
    *buffer = BitUtil::Log2(max_indice + 1);
  }
  encoder.Flush();
  return 1 + encoder.len();
}

template<typename T>
inline DictDecoder<T>::DictDecoder(uint8_t* dict_buffer, int dict_len,
    int fixed_len_size) {
  uint8_t* end = dict_buffer + dict_len;
  while (dict_buffer < end) {
    T value;
    dict_buffer +=
        ParquetPlainEncoder::Decode(dict_buffer, fixed_len_size, &value);
    dict_.push_back(value);
  }
}

template<typename T>
inline void DictDecoder<T>::Eq(int64_t num_rows, dynamic_bitset<>& skip_bitset,
    T& val) {
  typename vector<T>::iterator it = std::lower_bound(dict_.begin(), dict_.end(), val);
  if (it == dict_.end() || val < *it) {
    skip_bitset.resize(num_rows, false);
  } else {
    data_decoder_->Eq(num_rows, skip_bitset, std::distance(dict_.begin(), it));
  }
}

template<typename T>
inline void DictDecoder<T>::Gt(int64_t num_rows, dynamic_bitset<>& skip_bitset,
    T& val) {
  if (dict_.back() <= val) {
    skip_bitset.resize(num_rows, false);
  } else if (dict_[0] > val) {
    skip_bitset.resize(num_rows, true);
  } else {
    typename vector<T>::iterator it = std::upper_bound(dict_.begin(), dict_.end(), val);
    data_decoder_->Ge(num_rows, skip_bitset, std::distance(dict_.begin(), it));
  }
}

template<typename T>
inline void DictDecoder<T>::Lt(int64_t num_rows, dynamic_bitset<>& skip_bitset,
    T& val) {
  if (dict_[0] >= val) {
    skip_bitset.resize(num_rows, false);
  } else if (dict_.back() < val) {
    skip_bitset.resize(num_rows, true);
  } else {
    typename vector<T>::iterator it = std::lower_bound(dict_.begin(), dict_.end(), val);
    data_decoder_->Lt(num_rows, skip_bitset, std::distance(dict_.begin(), it));
  }
}

template<typename T>
inline void DictDecoder<T>::Ge(int64_t num_rows, dynamic_bitset<>& skip_bitset,
    T& val) {
  if (dict_.back() < val) {
    skip_bitset.resize(num_rows, false);
  } else if (dict_[0] >= val) {
    skip_bitset.resize(num_rows, true);
  } else {
    typename vector<T>::iterator it = std::lower_bound(dict_.begin(), dict_.end(), val);
    data_decoder_->Ge(num_rows, skip_bitset, std::distance(dict_.begin(), it));
  }
}

template<typename T>
inline void DictDecoder<T>::Le(int64_t num_rows, dynamic_bitset<>& skip_bitset, T& val) {
  if (dict_[0] > val) {
    skip_bitset.resize(num_rows, false);
  } else if (dict_.back() <= val) {
    skip_bitset.resize(num_rows, true);
  } else {
    typename vector<T>::iterator it = std::upper_bound(dict_.begin(), dict_.end(), val);
    data_decoder_->Lt(num_rows, skip_bitset, std::distance(dict_.begin(), it));
  }
}

template<typename T>
inline void DictDecoder<T>::In(int64_t num_rows, dynamic_bitset<>& skip_bitset,
    vector<T>& vals) {
  vector<uint64_t> tmp_vals;
  int vals_size = vals.size();
  for (int i = 0; i < vals_size; ++i) {
    typename vector<T>::iterator it = std::lower_bound(dict_.begin(), dict_.end(), vals[i]);
    if (it == dict_.end() || vals[i] < *it) {
      continue;
    } else {
      tmp_vals.push_back(std::distance(dict_.begin(), it));
    }
  }
  if (tmp_vals.size() == 0) {
    skip_bitset.resize(num_rows, false);
  } else {
    data_decoder_->In(num_rows, skip_bitset, tmp_vals);
  }
}

}
#endif

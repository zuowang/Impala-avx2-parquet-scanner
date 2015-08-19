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


#ifndef IMPALA_EXPRS_SIMPLE_PREDICATES_H_
#define IMPALA_EXPRS_SIMPLE_PREDICATES_H_

#include "exec/hdfs-parquet-scanner.h"

#include <string>
#include <boost/dynamic_bitset.hpp>

using namespace impala_udf;

namespace impala {

class HdfsParquetScanner;

class SimplePredicate {
 public:
  virtual void GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
      boost::dynamic_bitset<>& skip_bitset) = 0;
  virtual ~SimplePredicate() { }
};

class AndOperate : public SimplePredicate {
 public:
  virtual void GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
      boost::dynamic_bitset<>& skip_bitset);

  AndOperate(SimplePredicate* child0, SimplePredicate* child1)
  : child0_(child0),
    child1_(child1) {
  }

 private:
  SimplePredicate* child0_;
  SimplePredicate* child1_;
};

class OrOperate : public SimplePredicate {
 public:
  virtual void GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
      boost::dynamic_bitset<>& skip_bitset);

  OrOperate(SimplePredicate* child0, SimplePredicate* child1)
  : child0_(child0),
    child1_(child1) {
  }

 private:
  SimplePredicate* child0_;
  SimplePredicate* child1_;
};

template <typename T>
class EqOperate : public SimplePredicate {
 public:
  virtual void GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
      boost::dynamic_bitset<>& skip_bitset);

  EqOperate(int idx, T val) : idx_(idx), val_(val) { }

 private:
  int idx_;
  T val_;
};

template <typename T>
class LtOperate : public SimplePredicate {
 public:
  virtual void GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
      boost::dynamic_bitset<>& skip_bitset);

  LtOperate(int idx, T val) : idx_(idx), val_(val) { }

 private:
  int idx_;
  T val_;
};

template <typename T>
class LeOperate : public SimplePredicate {
 public:
  virtual void GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
      boost::dynamic_bitset<>& skip_bitset);

  LeOperate(int idx, T val) : idx_(idx), val_(val) { }

 private:
  int idx_;
  T val_;
};

template <typename T>
class GtOperate : public SimplePredicate {
 public:
  virtual void GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
      boost::dynamic_bitset<>& skip_bitset);

  GtOperate(int idx, T val) : idx_(idx), val_(val) { }

 private:
  int idx_;
  T val_;
};

template <typename T>
class GeOperate : public SimplePredicate {
 public:
  virtual void GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
      boost::dynamic_bitset<>& skip_bitset);

  GeOperate(int idx, T val) : idx_(idx), val_(val) { }

 private:
  int idx_;
  T val_;
};

template <typename T>
class InOperate : public SimplePredicate {
 public:
  virtual void GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
      boost::dynamic_bitset<>& skip_bitset);

  InOperate(int idx, vector<T> vals) : idx_(idx), vals_(vals.begin(), vals.end()) { }

 private:
  int idx_;
  vector<T> vals_;
};

inline void AndOperate::GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
    boost::dynamic_bitset<>& skip_bitset) {
  DCHECK(child0_);
  DCHECK(child1_);
  child0_->GetBitset(scanner, num_rows, skip_bitset);
  boost::dynamic_bitset<> tmp_bitset;
  child1_->GetBitset(scanner, num_rows, tmp_bitset);
  skip_bitset &= tmp_bitset;
}

inline void OrOperate::GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
    boost::dynamic_bitset<>& skip_bitset) {
  DCHECK(child0_);
  DCHECK(child1_);
  child0_->GetBitset(scanner, num_rows, skip_bitset);
  boost::dynamic_bitset<> tmp_bitset;
  child1_->GetBitset(scanner, num_rows, tmp_bitset);
  skip_bitset |= tmp_bitset;
}

template <typename T>
inline void EqOperate<T>::GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
    boost::dynamic_bitset<>& skip_bitset) {
  DCHECK(scanner);
  scanner->Eq(idx_, num_rows, skip_bitset, val_);
}

template <typename T>
inline void LtOperate<T>::GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
    boost::dynamic_bitset<>& skip_bitset) {
  DCHECK(scanner);
  scanner->Lt(idx_, num_rows, skip_bitset, val_);
}

template <typename T>
inline void LeOperate<T>::GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
    boost::dynamic_bitset<>& skip_bitset) {
  DCHECK(scanner);
  scanner->Le(idx_, num_rows, skip_bitset, val_);
}

template <typename T>
inline void GtOperate<T>::GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
    boost::dynamic_bitset<>& skip_bitset) {
  DCHECK(scanner);
  scanner->Gt(idx_, num_rows, skip_bitset, val_);
}

template <typename T>
inline void GeOperate<T>::GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
    boost::dynamic_bitset<>& skip_bitset) {
  DCHECK(scanner);
  scanner->Ge(idx_, num_rows, skip_bitset, val_);
}

template <typename T>
inline void InOperate<T>::GetBitset(HdfsParquetScanner* scanner, int64_t num_rows,
    boost::dynamic_bitset<>& skip_bitset) {
  DCHECK(scanner);
  scanner->In(idx_, num_rows, skip_bitset, vals_);
}

}

#endif

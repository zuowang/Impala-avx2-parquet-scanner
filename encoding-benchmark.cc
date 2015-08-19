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

// Benchmark to measure performance of fle and rle.

using namespace impala;

const int BUFFER_LEN = 64 * 4096;

struct TestData {
  uint64_t* buffer;
  uint64_t* buffer_end;
  int buffer_len;
  int bit_width;
  uint64_t* array;
  int num_values;
  int max_value;
  int rle_encode_len;
  int fle_encode_len;
  MemPool* pool;
  bool result;
};

void TestFleDecoder(int batch_size, void* d) {
  TestData* data = reinterpret_cast<TestData*>(d);
  for (int i = 0; i < batch_size; ++i) {
    // Unroll this to focus more on Put performance.
    FleDecoder decoder(reinterpret_cast<uint8_t*>(data->buffer), data->buffer_len,
        data->bit_width, data->num_values);
    for (int j = 0; j < data->num_values; j += 8) {
      uint64_t val;
      decoder.Get(&val);
      decoder.Get(&val);
      decoder.Get(&val);
      decoder.Get(&val);
      decoder.Get(&val);
      decoder.Get(&val);
      decoder.Get(&val);
      decoder.Get(&val);
    }
  }
}

void TestRleDecoder(int batch_size, void* d) {
  TestData* data = reinterpret_cast<TestData*>(d);
  for (int i = 0; i < batch_size; ++i) {
    // Unroll this to focus more on Put performance.
    RleDecoder decoder(reinterpret_cast<uint8_t*>(data->buffer), data->buffer_len,
        data->bit_width);
    for (int j = 0; j < data->num_values; j += 8) {
      uint64_t val;
      decoder.Get(&val);
      decoder.Get(&val);
      decoder.Get(&val);
      decoder.Get(&val);
      decoder.Get(&val);
      decoder.Get(&val);
      decoder.Get(&val);
      decoder.Get(&val);
    }
  }
}

void TestFleEncoder(int batch_size, void* d) {
  TestData* data = reinterpret_cast<TestData*>(d);
  for (int i = 0; i < batch_size; ++i) {
    // Unroll this to focus more on Put performance.
    FleEncoder encoder(reinterpret_cast<uint8_t*>(data->buffer), data->buffer_len,
        data->bit_width);
    for (int j = 0; j < data->num_values; j += 8) {
      encoder.Put(data->array[j + 0]);
      encoder.Put(data->array[j + 1]);
      encoder.Put(data->array[j + 2]);
      encoder.Put(data->array[j + 3]);
      encoder.Put(data->array[j + 4]);
      encoder.Put(data->array[j + 5]);
      encoder.Put(data->array[j + 6]);
      encoder.Put(data->array[j + 7]);
    }
    data->rle_encode_len = encoder.len();
  }
}

void TestRleEncoder(int batch_size, void* d) {
  TestData* data = reinterpret_cast<TestData*>(d);
  for (int i = 0; i < batch_size; ++i) {
    // Unroll this to focus more on Put performance.
    RleEncoder encoder(reinterpret_cast<uint8_t*>(data->buffer), data->buffer_len,
        data->bit_width);
    for (int j = 0; j < data->num_values; j += 8) {
      encoder.Put(data->array[j + 0]);
      encoder.Put(data->array[j + 1]);
      encoder.Put(data->array[j + 2]);
      encoder.Put(data->array[j + 3]);
      encoder.Put(data->array[j + 4]);
      encoder.Put(data->array[j + 5]);
      encoder.Put(data->array[j + 6]);
      encoder.Put(data->array[j + 7]);
    }
    data->fle_encode_len = encoder.len();
  }
}

int main(int argc, char** argv) {
  CpuInfo::Init();

  MemTracker tracker;
  MemPool pool(&tracker);

  int num_values = 4096;
  int max_bits = 32;

  srand((int)time(0));
  Benchmark encoder_suite("Encoder: Fle vs Rle");
  TestData data[max_bits];
  for (int i = 0; i < max_bits; ++i) {
    data[i].buffer = new uint64_t[BUFFER_LEN];
    data[i].buffer_end = data[i].buffer;
    data[i].buffer_len = BUFFER_LEN;
    data[i].array = new uint64_t[num_values];
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
    name << "\"Rle Encoder" << suffix.str() << "\"";
    int baseline = encoder_suite.AddBenchmark(name.str(), TestRleEncoder, &data[i], -1);

    name.str("");
    name << "\"Fle Encoder" << suffix.str() << "\"";
    encoder_suite.AddBenchmark(name.str(), TestFleEncoder, &data[i], baseline);
  }
  cout << encoder_suite.Measure() << endl;

  Benchmark decoder_suite("Decoder: Fle vs Rle");
  for (int i = 0; i < max_bits; ++i) {
    stringstream suffix;
    suffix << " " << (i+1) << "-Bit" << " compression ratio rle/fle:" << fixed <<
        data[i].rle_encode_len * 1.0/data[i].fle_encode_len;

    stringstream name;
    name.str("");
    name << "\"Rle Decoder" << suffix.str() << "\"";
    int baseline = decoder_suite.AddBenchmark(name.str(), TestRleDecoder, &data[i], -1);

    name.str("");
    name << "\"Fle Decoder" << suffix.str() << "\"";
    decoder_suite.AddBenchmark(name.str(), TestFleDecoder, &data[i], baseline);
  }
  cout << decoder_suite.Measure() << endl;
  return 0;
}

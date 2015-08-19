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

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include <boost/utility.hpp>
#include <gtest/gtest.h>
#include <math.h>

#include "common/init.h"
#include "util/fle-encoding.h"
#include "util/bit-stream-utils.inline.h"


using namespace std;

namespace impala {

const int MAX_WIDTH = 32;

TEST(BitArray, TestBool) {
  const int len = 8;
  uint8_t buffer[len];

  BitWriter writer(buffer, len);

  // Write alternating 0's and 1's
  for (int i = 0; i < 8; ++i) {
    bool result = writer.PutValue(i % 2, 1);
    EXPECT_TRUE(result);
  }
  writer.Flush();
  EXPECT_EQ((int)buffer[0], BOOST_BINARY(1 0 1 0 1 0 1 0));

  // Write 00110011
  for (int i = 0; i < 8; ++i) {
    bool result = false;
    switch (i) {
      case 0:
      case 1:
      case 4:
      case 5:
        result = writer.PutValue(false, 1);
        break;
      default:
        result = writer.PutValue(true, 1);
        break;
    }
    EXPECT_TRUE(result);
  }
  writer.Flush();

  // Validate the exact bit value
  EXPECT_EQ((int)buffer[0], BOOST_BINARY(1 0 1 0 1 0 1 0));
  EXPECT_EQ((int)buffer[1], BOOST_BINARY(1 1 0 0 1 1 0 0));

  // Use the reader and validate
  BitReader reader(buffer, len);
  for (int i = 0; i < 8; ++i) {
    bool val;
    bool result = reader.GetValue(1, &val);
    EXPECT_TRUE(result);
    EXPECT_EQ(val, i % 2);
  }

  for (int i = 0; i < 8; ++i) {
    bool val;
    bool result = reader.GetValue(1, &val);
    EXPECT_TRUE(result);
    switch (i) {
      case 0:
      case 1:
      case 4:
      case 5:
        EXPECT_EQ(val, false);
        break;
      default:
        EXPECT_EQ(val, true);
        break;
    }
  }
}

// Writes 'num_vals' values with width 'bit_width' and reads them back.
void TestBitArrayValues(int bit_width, int num_vals) {
  const int len = BitUtil::Ceil(bit_width * num_vals, 8);
  const uint64_t mod = bit_width == 64? 1 : 1LL << bit_width;

  uint8_t buffer[len];
  BitWriter writer(buffer, len);
  for (int i = 0; i < num_vals; ++i) {
    bool result = writer.PutValue(i % mod, bit_width);
    EXPECT_TRUE(result);
  }
  writer.Flush();
  EXPECT_EQ(writer.bytes_written(), len);

  BitReader reader(buffer, len);
  for (int i = 0; i < num_vals; ++i) {
    int64_t val;
    bool result = reader.GetValue(bit_width, &val);
    EXPECT_TRUE(result);
    EXPECT_EQ(val, i % mod);
  }
  EXPECT_EQ(reader.bytes_left(), 0);
}

TEST(BitArray, TestValues) {
  for (int width = 0; width <= MAX_WIDTH; ++width) {
    TestBitArrayValues(width, 1);
    TestBitArrayValues(width, 2);
    // Don't write too many values
    TestBitArrayValues(width, (width < 12) ? (1 << width) : 4096);
    TestBitArrayValues(width, 1024);
  }
}

// Test some mixed values
TEST(BitArray, TestMixed) {
  const int len = 1024;
  uint8_t buffer[len];
  bool parity = true;

  BitWriter writer(buffer, len);
  for (int i = 0; i < len; ++i) {
    bool result;
    if (i % 2 == 0) {
      result = writer.PutValue(parity, 1);
      parity = !parity;
    } else {
      result = writer.PutValue(i, 10);
    }
    EXPECT_TRUE(result);
  }
  writer.Flush();

  parity = true;
  BitReader reader(buffer, len);
  for (int i = 0; i < len; ++i) {
    bool result;
    if (i % 2 == 0) {
      bool val;
      result = reader.GetValue(1, &val);
      EXPECT_EQ(val, parity);
      parity = !parity;
    } else {
      int val;
      result = reader.GetValue(10, &val);
      EXPECT_EQ(val, i);
    }
    EXPECT_TRUE(result);
  }
}

// Validates encoding of values by encoding and decoding them.  If
// expected_encoding != NULL, also validates that the encoded buffer is
// exactly 'expected_encoding'.
// if expected_len is not -1, it will validate the encoded size is correct.
void ValidateFle(const vector<int>& values, int bit_width,
                 uint8_t* expected_encoding, int expected_len) {
  const int len = 64 * 1024;
  uint8_t buffer[len];
  EXPECT_LE(expected_len, len);

  FleEncoder encoder(buffer, len, bit_width);
  for (int i = 0; i < values.size(); ++i) {
    bool result = encoder.Put(values[i]);
    EXPECT_TRUE(result);
  }
  int encoded_len = encoder.Flush();

  if (expected_len != -1) {
    EXPECT_EQ(encoded_len, expected_len);
  }
  if (expected_encoding != NULL) {
    EXPECT_TRUE(memcmp(buffer, expected_encoding, expected_len) == 0);
  }

  // Verify read
  FleDecoder decoder(buffer, len, bit_width);
  for (int i = 0; i < values.size(); ++i) {
    uint64_t val;
    bool result = decoder.Get(&val);
    EXPECT_TRUE(result);
    EXPECT_EQ(values[i], val);
  }
}

TEST(Fle, SpecificSequences) {
//  const int len = 1024;
//  uint8_t expected_buffer[len];
  vector<int> values;

  // Test 50 0' followed by 50 1's
  values.resize(100);
  for (int i = 0; i < 50; ++i) {
    values[i] = 0;
  }
  for (int i = 50; i < 100; ++i) {
    values[i] = 1;
  }

  // expected_buffer valid for bit width <= 1 byte
//  *(reinterpret_cast<int64_t*>(&expected_buffer[0])) = 0x3fff;
//  *(reinterpret_cast<int64_t*>(&expected_buffer[8])) = 0xfffffffff0000000;
  ValidateFle(values, 1, NULL, 16);
//  *(reinterpret_cast<int64_t*>(&expected_buffer[0])) = 0x3fff;
//  *(reinterpret_cast<int64_t*>(&expected_buffer[8])) = 0x0;
//  *(reinterpret_cast<int64_t*>(&expected_buffer[16])) = 0xfffffffff0000000;
//  *(reinterpret_cast<int64_t*>(&expected_buffer[24])) = 0x0;
  ValidateFle(values, 2, NULL, 32);

  for (int width = 2; width <= MAX_WIDTH; ++width) {
    ValidateFle(values, width, NULL, -1);
  }

  // Test 100 0's and 1's alternating
  for (int i = 0; i < 100; ++i) {
    values[i] = i % 2;
  }

//  *(reinterpret_cast<int64_t*>(&expected_buffer[0])) = 0x5555555555555555;
//  *(reinterpret_cast<int64_t*>(&expected_buffer[8])) = 0x5555555550000000;

  ValidateFle(values, 1, NULL, 16);
  for (int width = 2; width <= MAX_WIDTH; ++width) {
    ValidateFle(values, width, NULL, -1);
  }
}

// ValidateFle on 'num_vals' values with width 'bit_width'. If 'value' != -1, that value
// is used, otherwise alternating values are used.
void TestFleValues(int bit_width, int num_vals, int value = -1) {
  const uint64_t mod = (bit_width == 64) ? 1 : 1LL << bit_width;
  vector<int> values;
  for (int v = 0; v < num_vals; ++v) {
    values.push_back((value != -1) ? value : (v % mod));
  }
  ValidateFle(values, bit_width, NULL, -1);
}

TEST(Fle, TestValues) {
  for (int width = 1; width <= MAX_WIDTH; ++width) {
    TestFleValues(width, 1);
    TestFleValues(width, 1024);
    TestFleValues(width, 1024, 0);
    TestFleValues(width, 1024, 1);
  }
}

TEST(Fle, TestRandomValues) {
  srand((int)time(0));
  for (int width = 1; width <= MAX_WIDTH; ++width) {
    const uint64_t mod = (width == 64) ? 1 : 1LL << width;
    vector<int> values;
    for (int v = 0; v < 1024; ++v) {
      values.push_back(rand() % mod);
    }

    ValidateFle(values, width, NULL, -1);
  }
}
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  impala::InitCommonRuntime(argc, argv, true);
  return RUN_ALL_TESTS();
}

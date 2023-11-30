#include <gtest/gtest.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <random>

#include "hmm/internal/hmmintrin.hpp"

TEST(hmmTest, abs)
{
  std::mt19937                          gen(98);
  std::uniform_real_distribution<float> dis(-1000.f, +1000.f);
  for (std::size_t i = 0; i < 1000; ++i)
  {
    float value = dis(gen);
    EXPECT_FLOAT_EQ(std::fabs(value), hmm::fops::abs(value));
  }
}

TEST(hmmTest, f32hex)
{
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0xFF800000), -std::numeric_limits<float>::infinity());
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0xC0000000), -2.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0xBF800000), -1.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0xBF000000), -0.5f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x80000000), -0.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x00000000), +0.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x3C8EFA35), M_PI / 180.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x3DCCCCCD), 0.1f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x3E9A209B), log10f(2.f));
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x3EBC5AB2), 1.f / M_E);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x3EDE5BD9), 1.f / M_LN10);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x3F000000), 0.5f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x3F317218), M_LN2);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x3F3504F3), M_SQRT1_2);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x3F690570), 1.f / logf(3.f));
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x3F800000), 1.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x3F8C9F54), log(3.f));
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x3FB504F3), M_SQRT2);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x3FB8AA3B), 1.f / M_LN2);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x3FDDB3D7), sqrtf(3.f));
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x40000000), 2.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x40135D8E), M_LN10);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x402DF854), M_E);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x40400000), 3.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x40490FDB), M_PI);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x404A62C2), sqrtf(10.f));
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x40549A78), log2f(10.f));
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x40800000), 4.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x40A00000), 5.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x40C00000), 6.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x40C90FDB), 2.f * M_PI);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x40E00000), 7.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x41000000), 8.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x41100000), 9.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x41200000), 10.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x41300000), 11.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x41400000), 12.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x41500000), 13.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x41600000), 14.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x41700000), 15.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x41800000), 16.f);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x42652EE1), 180.f * M_1_PI);
  EXPECT_FLOAT_EQ(hmm::fops::hex_to_f32(0x7F800000), std::numeric_limits<float>::infinity());
}

TEST(hmmTest, hexf32)
{
  EXPECT_EQ(0xFF800000, hmm::fops::f32_to_hex(-std::numeric_limits<float>::infinity()));
  EXPECT_EQ(0xC0000000, hmm::fops::f32_to_hex(-2.f));
  EXPECT_EQ(0xBF800000, hmm::fops::f32_to_hex(-1.f));
  EXPECT_EQ(0xBF000000, hmm::fops::f32_to_hex(-0.5f));
  EXPECT_EQ(0x80000000, hmm::fops::f32_to_hex(-0.f));
  EXPECT_EQ(0x00000000, hmm::fops::f32_to_hex(+0.f));
  EXPECT_EQ(0x3C8EFA35, hmm::fops::f32_to_hex(M_PI / 180.f));
  EXPECT_EQ(0x3DCCCCCD, hmm::fops::f32_to_hex(0.1f));
  EXPECT_EQ(0x3E9A209B, hmm::fops::f32_to_hex(log10f(2.f)));
  EXPECT_EQ(0x3EBC5AB2, hmm::fops::f32_to_hex(1.f / M_E));
  EXPECT_EQ(0x3EDE5BD9, hmm::fops::f32_to_hex(1.f / M_LN10));
  EXPECT_EQ(0x3F000000, hmm::fops::f32_to_hex(0.5f));
  EXPECT_EQ(0x3F317218, hmm::fops::f32_to_hex(M_LN2));
  EXPECT_EQ(0x3F3504F3, hmm::fops::f32_to_hex(M_SQRT1_2));
  EXPECT_EQ(0x3F690570, hmm::fops::f32_to_hex(1.f / logf(3.f)));
  EXPECT_EQ(0x3F800000, hmm::fops::f32_to_hex(1.f));
  EXPECT_EQ(0x3F8C9F54, hmm::fops::f32_to_hex(log(3.f)));
  EXPECT_EQ(0x3FB504F3, hmm::fops::f32_to_hex(M_SQRT2));
  EXPECT_EQ(0x3FB8AA3B, hmm::fops::f32_to_hex(1.f / M_LN2));
  EXPECT_EQ(0x3FDDB3D7, hmm::fops::f32_to_hex(sqrtf(3.f)));
  EXPECT_EQ(0x40000000, hmm::fops::f32_to_hex(2.f));
  EXPECT_EQ(0x40135D8E, hmm::fops::f32_to_hex(M_LN10));
  EXPECT_EQ(0x402DF854, hmm::fops::f32_to_hex(M_E));
  EXPECT_EQ(0x40400000, hmm::fops::f32_to_hex(3.f));
  EXPECT_EQ(0x40490FDB, hmm::fops::f32_to_hex(M_PI));
  EXPECT_EQ(0x404A62C2, hmm::fops::f32_to_hex(sqrtf(10.f)));
  EXPECT_EQ(0x40549A78, hmm::fops::f32_to_hex(log2f(10.f)));
  EXPECT_EQ(0x40800000, hmm::fops::f32_to_hex(4.f));
  EXPECT_EQ(0x40A00000, hmm::fops::f32_to_hex(5.f));
  EXPECT_EQ(0x40C00000, hmm::fops::f32_to_hex(6.f));
  EXPECT_EQ(0x40C90FDB, hmm::fops::f32_to_hex(2.f * M_PI));
  EXPECT_EQ(0x40E00000, hmm::fops::f32_to_hex(7.f));
  EXPECT_EQ(0x41000000, hmm::fops::f32_to_hex(8.f));
  EXPECT_EQ(0x41100000, hmm::fops::f32_to_hex(9.f));
  EXPECT_EQ(0x41200000, hmm::fops::f32_to_hex(10.f));
  EXPECT_EQ(0x41300000, hmm::fops::f32_to_hex(11.f));
  EXPECT_EQ(0x41400000, hmm::fops::f32_to_hex(12.f));
  EXPECT_EQ(0x41500000, hmm::fops::f32_to_hex(13.f));
  EXPECT_EQ(0x41600000, hmm::fops::f32_to_hex(14.f));
  EXPECT_EQ(0x41700000, hmm::fops::f32_to_hex(15.f));
  EXPECT_EQ(0x41800000, hmm::fops::f32_to_hex(16.f));
  EXPECT_EQ(0x42652EE1, hmm::fops::f32_to_hex(180.f * M_1_PI));
  EXPECT_EQ(0x7F800000, hmm::fops::f32_to_hex(std::numeric_limits<float>::infinity()));
}

TEST(hmmTest, compare)
{
  float f1 = 0.1f * 9.f;
  float f2 = 0.9f;

  // 0.1f * 9.f == 0.9f yields false
  EXPECT_NE(f1, f2);
  EXPECT_FLOAT_EQ(f1, f2);

  EXPECT_FALSE(hmm::fops::cmplt(f1, f2));
  EXPECT_TRUE(hmm::fops::cmple(f1, f2));
  EXPECT_FALSE(hmm::fops::cmpgt(f1, f2));
  EXPECT_TRUE(hmm::fops::cmpge(f1, f2));
  EXPECT_TRUE(hmm::fops::cmpeq(f1, f2));
}
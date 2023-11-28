#include <gtest/gtest.h>

#include "hmm/audio/Windows.hpp"

TEST(hmmTest, amplitude_correction_factors)
{
  const std::size_t len = 16777216;
  EXPECT_FLOAT_EQ(1.f, hmm::Rectangular::correct(len));
  EXPECT_FLOAT_EQ(2.f, hmm::Hann::correct(len));
}

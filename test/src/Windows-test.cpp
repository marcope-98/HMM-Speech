#include <gtest/gtest.h>

#include "hmm/audio/Windows.hpp"

TEST(hmmTest, amplitude_correction_factors)
{
  const std::size_t len = 1 << 16;
  const float       eps = 1e-4;
  EXPECT_NEAR(1.f, hmm::Rectangular::correct(len), eps);
  EXPECT_NEAR(2.f, hmm::Hann::correct(len), eps);
  EXPECT_NEAR(1.8519f, hmm::Hamming::correct(len), eps);
  EXPECT_NEAR(2.f, hmm::Barthann::correct(len), eps);
  EXPECT_NEAR(2.f, hmm::Bartlett::correct(len), eps);
  EXPECT_NEAR(2.3810f, hmm::Blackman::correct(len), eps);
  EXPECT_NEAR(2.7875f, hmm::BlackmanHarris::correct(len), eps);
  EXPECT_NEAR(2.7505f, hmm::Nuttall::correct(len), eps);
  EXPECT_NEAR(4.6387f, hmm::Flattop::correct(len), eps);
}

#include <gtest/gtest.h>

#include "fft-test.hpp"

#include "hmm/audio/fft.hpp"
#include "hmm/audio/windows.hpp"
#include "hmm/internal/hmmintrin.hpp"

void compare_ctfft(const float *const src, const std::size_t &Nfft)
{
  std::complex<float> *reference = hmm::fft::dft(src, Nfft);
  std::complex<float> *to_test   = hmm::fft::cooley_tukey(src, Nfft);

  for (std::size_t i = 0; i < Nfft; ++i)
  {
    ASSERT_FLOAT_EQ(reference[i].real(), to_test[i].real()) << i;
    ASSERT_FLOAT_EQ(reference[i].real(), to_test[i].real()) << i;
  }
  delete[] reference;
  delete[] to_test;
}

TEST(hmmTest, round_up_to_multiple_of)
{
  EXPECT_EQ(0, hmm::fft::round_up_to_multiple_of(0, 256));
  EXPECT_EQ(256, hmm::fft::round_up_to_multiple_of(1, 256));
  EXPECT_EQ(512, hmm::fft::round_up_to_multiple_of(257, 256));
  EXPECT_EQ(768, hmm::fft::round_up_to_multiple_of(513, 256));
  EXPECT_EQ(1280, hmm::fft::round_up_to_multiple_of(1025, 256));
  EXPECT_EQ(2304, hmm::fft::round_up_to_multiple_of(2049, 256));
  EXPECT_EQ(4352, hmm::fft::round_up_to_multiple_of(4097, 256));
}

TEST(hmmTest, clp2)
{
  for (std::size_t i = 0; i < 31; ++i)
    EXPECT_EQ((1u << (i + 1)), hmm::fft::clp2((1u << i) + 1));
}

TEST(hmmTest, amplitude_cf)
{
  const std::size_t len = 16777216;
  EXPECT_FLOAT_EQ(1.f, hmm::fft::amplitude_cf(&hmm::windows::rectangular, len));
  EXPECT_FLOAT_EQ(2.f, hmm::fft::amplitude_cf(&hmm::windows::hann, len));
}

TEST(hmmTest, cooley_tukey)
{
  const std::size_t Nfft = 4096;
  float *           src  = new float[Nfft]();

  // zeros test
  fft_test::fill(src, Nfft, &fft_test::zeros);
  compare_ctfft(src, Nfft);

  // ones test
  fft_test::fill(src, Nfft, &fft_test::ones);
  compare_ctfft(src, Nfft);

  // alternate test
  fft_test::fill(src, Nfft, &fft_test::alternate);
  compare_ctfft(src, Nfft);

  // cosine1 test
  fft_test::fill(src, Nfft, &fft_test::cosine1);
  compare_ctfft(src, Nfft);

  // cosine2 test
  fft_test::fill(src, Nfft, &fft_test::cosine2);
  compare_ctfft(src, Nfft);
}

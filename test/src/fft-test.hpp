#ifndef TEST_FFT_TEST_HPP_
#define TEST_FFT_TEST_HPP_

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <cstdint>

using namespace std::complex_literals;

struct fft_test
{
  static float ones(const std::size_t &n, const std::size_t &N)
  {
    (void)sizeof(n);
    (void)sizeof(N);
    return 1.f;
  }

  static float zeros(const std::size_t &n, const std::size_t &N)
  {
    (void)sizeof(n);
    (void)sizeof(N);
    return 0.f;
  }

  static float alternate(const std::size_t &n, const std::size_t &N)
  {
    (void)sizeof(N);
    return (n % 2 == 0) ? 1.f : -1.f;
  }

  static float cosine1(const std::size_t &n, const std::size_t &N)
  {
    return cosf(8.f * 2.f * M_PI * float(n) / float(N));
  }

  static float cosine2(const std::size_t &n, const std::size_t &N)
  {
    return cosf((43.f / 7.f) * 2.f * M_PI * float(n) / float(N));
  }

  static void fill(float *const src, const std::size_t &Nfft,
                   float (*func)(const std::size_t &, const std::size_t &))
  {
    for (std::size_t i = 0; i < Nfft; ++i)
      src[i] = func(i, Nfft);
  }
};

#endif
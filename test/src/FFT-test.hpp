#ifndef TEST_FFT_TEST_HPP_
#define TEST_FFT_TEST_HPP_

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <cstdint>

using namespace std::complex_literals;
struct FFT_test
{
  static float sign(const float &value) { return value > 0.f ? 1.f : -1.f; }

  static std::complex<float> impulse_response(const std::size_t &n, const std::size_t &N)
  {
    (void)sizeof(N);
    return float(n == 1);
  }

  static std::complex<float> complex_sinusoid(const std::size_t &n, const std::size_t &N)
  {
    const float arg = 2.f * M_PI * float(n) / float(N);
    return cosf(arg) + sinf(arg) * 1.if;
  }

  static std::complex<float> rectangle(const std::size_t &n, const std::size_t &N)
  {
    const float arg = 2.f * M_PI * float(n) / float(N);
    return sign(cosf(arg)) + sign(sinf(arg)) * 1.if;
  }

  static void fill(std::complex<float> *const src, const std::size_t &Nfft,
                   std::complex<float> (*func)(const std::size_t &, const std::size_t &))
  {
    for (std::size_t i = 0; i < Nfft; ++i)
      src[i] = func(i, Nfft);
  }
};

#endif // TEST_FFT_TEST_HPP_
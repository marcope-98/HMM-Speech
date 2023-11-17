#ifndef HMM_AUDIO_WINDOWS_HPP_
#define HMM_AUDIO_WINDOWS_HPP_

#define _USE_MATH_DEFINES
#include <cmath>

#include "hmm/internal/hmmintrin.hpp"

namespace hmm
{
  struct windows
  {
    static float rectangular(const std::size_t &n, const std::size_t &N)
    {
      (void)sizeof(n);
      (void)sizeof(N);
      return 1.f;
    }

    static float hann(const std::size_t &n, const std::size_t &N)
    {
      return 0.5f * (1.f - cosf(2.f * M_PI * float(n) / float(N)));
    }

    static float hamming(const std::size_t &n, const std::size_t &N)
    {
      return 0.54f - 0.46f * cosf(2.f * M_PI * float(n) / float(N));
    }

    static float barthann(const std::size_t &n, const std::size_t &N)
    {
      const float arg = float(n) / float(N) - 0.5f;
      return 0.62f - 0.48f * fops::abs(arg) + 0.38f * cosf(2.f * M_PI * arg);
    }

    static float bartlett(const std::size_t &n, const std::size_t &N)
    {
      const float arg = 2.f * float(n) / float(N);
      if (fops::cmple(arg, 1.f))
        return arg;
      else
        return 2.f - arg;
    }

    static float blackman(const std::size_t &n, const std::size_t &N)
    {
      const float arg = 2.f * M_PI * float(n) / float(N);
      return 0.42f -
             0.5f * cosf(arg) +
             0.08f * cosf(2.f * arg);
    }

    static float blackmanharris(const std::size_t &n, const std::size_t &N)
    {
      const float arg = 2.f * M_PI * float(n) / float(N);
      return 0.35875f -
             0.48829f * cosf(arg) +
             0.14128f * cosf(2.f * arg) -
             0.01168f * cosf(3.f * arg);
    }

    static float nuttall(const std::size_t &n, const std::size_t &N)
    {
      const float arg = 2.f * M_PI * float(n) / float(N);
      return 0.3635819f -
             0.4891775 * cosf(arg) +
             0.1365995f * cosf(2.f * arg) -
             0.0106411f * cosf(3.f * arg);
    }

    static float flattop(const std::size_t &n, const std::size_t &N)
    {
      const float arg = 2.f * M_PI * float(n) / float(N);
      return 0.21557895f -
             0.41663158f * cosf(arg) +
             0.277263158f * cosf(2.f * arg) -
             0.083578947f * cosf(3.f * arg) +
             0.006947368f * cosf(4.f * arg);
    }
  };
} // namespace hmm

#endif // HMM_AUDIO_WINDOWS_HPP_
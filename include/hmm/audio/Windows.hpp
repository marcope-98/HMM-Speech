#ifndef HMM_AUDIO_WINDOWS_HPP_
#define HMM_AUDIO_WINDOWS_HPP_

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdint>
#include <iostream>

#include "hmm/internal/Ptr.hpp"
#include "hmm/internal/hmmintrin.hpp"

namespace hmm
{
  struct Rectangular
  {
    auto         execute(Ptr p) { return p; }
    static float correct(const std::size_t &N)
    {
      (void)sizeof(N);
      return 1.f;
    }
  };

  struct Hann
  {
    auto execute(Ptr p)
    {
      float *     dst = (float *)p.data;
      const float arg = 2.f * M_PI / float(p.size - 1);
      for (std::size_t i = 0; i < p.size; ++i)
        dst[i] *= 0.5f * (1.f - cosf(arg * float(i)));
      return p;
    }

    static float correct(const std::size_t &N)
    {
      float       Aw  = 0.f;
      const float arg = 2.f * M_PI / float(N - 1);
      for (std::size_t i = 0; i < N; ++i)
        Aw += 0.5f * (1.f - cosf(arg * float(i)));
      return float(N) / Aw;
    }
  };

  struct Hamming
  {
    auto execute(Ptr p)
    {
      float *     dst = (float *)p.data;
      const float arg = 2.f * M_PI / float(p.size - 1);
      for (std::size_t i = 0; i < p.size; ++i)
        dst[i] *= 0.54f - 0.46f * cosf(arg * float(i));
      return p;
    }

    static float correct(const std::size_t &N)
    {
      float       Aw  = 0.f;
      const float arg = 2.f * M_PI / float(N - 1);
      for (std::size_t i = 0; i < N; ++i)
        Aw += 0.54f - 0.46f * cosf(arg * float(i));
      return float(N) / Aw;
    }
  };

  struct Barthann
  {
    auto execute(Ptr p)
    {
      float *     dst = (float *)p.data;
      const float pi2 = 2.f * M_PI;
      for (std::size_t i = 0; i < p.size; ++i)
      {
        float arg = float(i) / float(p.size - 1) - 0.5f;
        dst[i] *= 0.62f - 0.48f * fops::abs(arg) + 0.38f * cosf(pi2 * arg);
      }
      return p;
    }

    static float correct(const std::size_t &N)
    {
      float       Aw  = 0.f;
      const float pi2 = 2.f * M_PI;
      for (std::size_t i = 0; i < N; ++i)
      {
        float arg = float(i) / float(N - 1) - 0.5f;
        Aw += 0.62f - 0.48f * fops::abs(arg) + 0.38f * cosf(pi2 * arg);
      }
      return float(N) / Aw;
    }
  };

  struct Bartlett
  {
    auto execute(Ptr p)
    {
      float *     dst = (float *)p.data;
      const float arg = 2.f / float(p.size - 1);
      for (std::size_t i = 0; i < (p.size - 1) / 2; ++i)
        dst[i] *= arg * float(i);
      for (std::size_t i = (p.size - 1) / 2; i < p.size; ++i)
        dst[i] *= 2.f - arg * float(i);
      return p;
    }

    static float correct(const std::size_t &N)
    {
      float       Aw  = 0.f;
      const float arg = 2.f / float(N - 1);
      for (std::size_t i = 0; i < (N - 1) / 2; ++i)
        Aw += arg * float(i);
      for (std::size_t i = (N - 1) / 2; i < N; ++i)
        Aw += 2.f - arg * float(i);
      return float(N) / Aw;
    }
  };

  struct Blackman
  {
    auto execute(Ptr p)
    {
      float *     dst = (float *)p.data;
      const float arg = 2.f * M_PI / float(p.size - 1);
      for (std::size_t i = 0; i < p.size; ++i)
      {
        float temp = arg * float(i);
        dst[i] *= 0.42f -
                  0.5f * cosf(temp) +
                  0.08f * cosf(2.f * temp);
      }
      return p;
    }

    static float correct(const std::size_t &N)
    {
      float       Aw  = 0.f;
      const float arg = 2.f * M_PI / float(N - 1);
      for (std::size_t i = 0; i < N; ++i)
      {
        float temp = arg * float(i);
        Aw += 0.42f -
              0.5f * cosf(temp) +
              0.08f * cosf(2.f * temp);
      }
      return float(N) / Aw;
    }
  };

  struct BlackmanHarris
  {
    auto execute(Ptr p)
    {
      float *     dst = (float *)p.data;
      const float arg = 2.f * M_PI / float(p.size - 1);
      for (std::size_t i = 0; i < p.size; ++i)
      {
        float temp = arg * float(i);
        dst[i] *= 0.35875f -
                  0.48829f * cosf(temp) +
                  0.14128f * cosf(2.f * temp) -
                  0.01168f * cosf(3.f * temp);
      }
      return p;
    }

    static float correct(const std::size_t &N)
    {
      float       Aw  = 0.f;
      const float arg = 2.f * M_PI / float(N - 1);
      for (std::size_t i = 0; i < N; ++i)
      {
        float temp = arg * float(i);
        Aw += 0.35875f -
              0.48829f * cosf(temp) +
              0.14128f * cosf(2.f * temp) -
              0.01168f * cosf(3.f * temp);
      }
      return float(N) / Aw;
    }
  };

  struct Nuttall
  {
    auto execute(Ptr p)
    {
      float *     dst = (float *)p.data;
      const float arg = 2.f * M_PI / float(p.size - 1);
      for (std::size_t i = 0; i < p.size; ++i)
      {
        float temp = arg * float(i);
        dst[i] *= 0.3635819f -
                  0.4891775 * cosf(temp) +
                  0.1365995f * cosf(2.f * temp) -
                  0.0106411f * cosf(3.f * temp);
      }
      return p;
    }

    static float correct(const std::size_t &N)
    {
      float       Aw  = 0.f;
      const float arg = 2.f * M_PI / float(N - 1);
      for (std::size_t i = 0; i < N; ++i)
      {
        float temp = arg * float(i);
        Aw += 0.3635819f -
              0.4891775 * cosf(temp) +
              0.1365995f * cosf(2.f * temp) -
              0.0106411f * cosf(3.f * temp);
      }
      return float(N) / Aw;
    }
  };

  struct Flattop
  {
    auto execute(Ptr p)
    {
      float *     dst = (float *)p.data;
      const float arg = 2.f * M_PI / float(p.size - 1);
      for (std::size_t i = 0; i < p.size; ++i)
      {
        float temp = arg * float(i);
        dst[i] *= 0.21557895f -
                  0.41663158f * cosf(temp) +
                  0.277263158f * cosf(2.f * temp) -
                  0.083578947f * cosf(3.f * temp) +
                  0.006947368f * cosf(4.f * temp);
      }
      return p;
    }

    static float correct(const std::size_t &N)
    {
      float       Aw  = 0.f;
      const float arg = 2.f * M_PI / float(N - 1);
      for (std::size_t i = 0; i < N; ++i)
      {
        float temp = arg * float(i);
        Aw += 0.21557895f -
              0.41663158f * cosf(temp) +
              0.277263158f * cosf(2.f * temp) -
              0.083578947f * cosf(3.f * temp) +
              0.006947368f * cosf(4.f * temp);
      }
      return float(N) / Aw;
    }
  };

} // namespace hmm

#endif // HMM_AUDIO_WINDOWS_HPP_
#ifndef HMM_INTERNAL_FFT_HPP_
#define HMM_INTERNAL_FFT_HPP_

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdint>
#include <utility>

#include <complex>

#include "hmm/internal/Ptr.hpp"

using namespace std::complex_literals;

namespace hmm
{
  struct ToComplex
  {
    auto execute(Ptr p)
    {
      const std::size_t size = p.size;
      // check if there is enough space to avoid allocation
      if (size * sizeof(std::complex<float>) <= p.capacity)
      {
        float *              src = (float *)p.data;
        std::complex<float> *dst = (std::complex<float> *)p.data;

        for (std::size_t i = size - 1; i != 0; --i)
          dst[i] = src[i];
        dst[0] = src[0];
      }
      else
      {
        // allocate
        void *temp = malloc(sizeof(std::complex<float>) * size);

        float *              src = (float *)p.data;
        std::complex<float> *dst = (std::complex<float> *)temp;

        // copy
        for (std::size_t i = 0; i < size; ++i)
          dst[i] = src[i];

        // swap
        std::swap(temp, p.data);
        free(temp);

        p.capacity = sizeof(std::complex<float>) * size;
      }

      return p;
    }
  };

  struct DFT
  {
    auto execute(Ptr p)
    {
      // We cannot avoid memory allocation for these method
      const std::size_t size = p.size;
      void *            temp = calloc(size, sizeof(std::complex<float>));

      // dft
      std::complex<float> *src      = (std::complex<float> *)p.data;
      std::complex<float> *dst      = (std::complex<float> *)temp;
      const float          constant = -2.f * M_PI / float(size);

      for (std::size_t k = 0; k < size; ++k)
      {
        float factor = constant * float(k);
        for (std::size_t n = 0; n < size; ++n)
          dst[k] += src[n] * std::exp(factor * float(n) * 1.if);
      }

      // swap
      std::swap(temp, p.data);
      free(temp);
      return p;
    }
  };

  struct CooleyTukey
  {
    auto execute(Ptr p)
    {
      std::size_t         k     = p.size;
      const float         theta = M_PI / float(p.size);
      std::complex<float> phi   = std::complex<float>(cosf(theta), -sinf(theta));
      std::complex<float> T;
      std::size_t         n;

      std::complex<float> *src = (std::complex<float> *)p.data;
      while (k > 1)
      {
        n = k;
        k >>= 1;
        phi *= phi;
        T = 1.f;
        for (std::size_t l = 0; l < k; ++l)
        {
          for (std::size_t a = l; a < p.size; a += n)
          {
            std::size_t         b = a + k;
            std::complex<float> t = src[a] - src[b];
            src[a] += src[b];
            src[b] = t * T;
          }
          T *= phi;
        }
      }

      // Bit Reverse
      const std::uint32_t m = (std::uint32_t)log2f(p.size);
      for (std::uint32_t a = 0u; a < (std::uint32_t)(p.size); ++a)
      {
        std::uint32_t b = a;
        b               = (((b & 0xaaaaaaaau) >> 1) | ((b & 0x55555555u) << 1));
        b               = (((b & 0xccccccccu) >> 2) | ((b & 0x33333333u) << 2));
        b               = (((b & 0xf0f0f0f0u) >> 4) | ((b & 0x0f0f0f0fu) << 4));
        b               = (((b & 0xff00ff00u) >> 8) | ((b & 0x00ff00ffu) << 8));
        b               = ((b >> 16u) | (b << 16u)) >> (32u - m);

        if (b > a)
        {
          std::complex<float> t = src[a];
          src[a]                = src[b];
          src[b]                = t;
        }
      }
      return p;
    }
  };

  struct Magnitude
  {
    auto execute(Ptr p)
    {
      const std::size_t    size = p.size;
      std::complex<float> *src  = (std::complex<float> *)p.data;
      float *              dst  = (float *)p.data;
      for (std::size_t i = 0; i < size; ++i)
        dst[i] = std::abs(src[i]);
      return p;
    }
  };

  struct Phase
  {
    auto execute(Ptr p)
    {
      const std::size_t    size = p.size;
      std::complex<float> *src  = (std::complex<float> *)p.data;
      float *              dst  = (float *)p.data;
      for (std::size_t i = 0; i < size; ++i)
        dst[i] = std::arg(src[i]);
      return p;
    }
  };
} // namespace hmm

#endif // HMM_INTERNAL_FFT_HPP_
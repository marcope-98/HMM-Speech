#ifndef HMM_FFT_HPP_
#define HMM_FFT_HPP_

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <cstdint>

using namespace std::complex_literals;

namespace hmm
{
  struct fft
  {
    static float *zero_padding(const float *const src, const std::size_t &prevsize, const size_t &size)
    {
      float *dst = new float[size]();
      // TODO: optimize
      for (std::size_t i = 0; i < prevsize; ++i)
        dst[i] = src[i];
      return dst;
    }

    static std::size_t round_up_to_multiple_of(const std::size_t &value, const std::size_t &mul) { return (value + (mul - 1)) & ~(mul - 1); }

    // TODO: what about energy correction factor ?
    static void apply_window(float *const src, const std::size_t &window_length,
                             float (*window)(const std::size_t &, const std::size_t &))
    {
      for (std::size_t i = 0; i < window_length; ++i)
        src[i] *= window(i, window_length - 1);
    }

    static std::size_t clp2(std::size_t x)
    {
      // https://www.bitweenie.com/listings/fft-zero-padding/
      // The FFT resolution should at least support the same resolution as your waveform frequency resolution.    Additionally, some highly-efficient implementations of the FFT require that the number of FFT points be a power of two.
      x -= 1;
      x |= x >> 1;
      x |= x >> 2;
      x |= x >> 4;
      x |= x >> 8;
      x |= x >> 16;
      return x + 1;
    }

    static float mean(const float *const src, const std::size_t &N)
    {
      float sum = 0.f;
      // TODO: to optimize
      for (std::size_t i = 0; i < N; ++i)
        sum += src[i];
      return sum / float(N);
    }

    static std::complex<float> *cooley_tukey(const float *const src, const std::size_t &Nfft,
                                             std::complex<float> *out = nullptr)
    {
      std::complex<float> *res{nullptr};
      if (out == nullptr)
        res = new std::complex<float>[Nfft]();
      else
        res = out;

      std::size_t         k     = Nfft;
      const float         theta = M_PI / float(Nfft);
      std::complex<float> phi   = std::complex<float>(cosf(theta), -sinf(theta));
      std::complex<float> T;
      std::size_t         n;
      // TODO: optimize
      for (std::size_t i = 0; i < Nfft; ++i)
        res[i] = src[i];

      while (k > 1)
      {
        n = k;
        k >>= 1;
        phi *= phi;
        T = 1.f;
        for (std::size_t l = 0; l < k; ++l)
        {
          for (std::size_t a = l; a < Nfft; a += n)
          {
            std::size_t         b = a + k;
            std::complex<float> t = res[a] - res[b];
            res[a] += res[b];
            res[b] = t * T;
          }
          T *= phi;
        }
      }

      // Bit reverse
      const std::uint32_t m = (std::size_t)log2f(Nfft);
      for (std::uint32_t a = 0u; a < (std::uint32_t)(Nfft); ++a)
      {
        std::uint32_t b = a;
        b               = (((b & 0xaaaaaaaau) >> 1) | ((b & 0x55555555u) << 1));
        b               = (((b & 0xccccccccu) >> 2) | ((b & 0x33333333u) << 2));
        b               = (((b & 0xf0f0f0f0u) >> 4) | ((b & 0x0f0f0f0fu) << 4));
        b               = (((b & 0xff00ff00u) >> 8) | ((b & 0x00ff00ffu) << 8));
        b               = ((b >> 16u) | (b << 16u)) >> (32u - m);

        if (b > a)
        {
          std::complex<float> t = res[a];
          res[a]                = res[b];
          res[b]                = t;
        }
      }

      return res;
    }

    static std::complex<float> *dft(const float *const src, const std::size_t &Nfft, std::complex<float> *out = nullptr)
    {
      std::complex<float> *res;
      // NOTE: this function assumes that src has at least Nfft elements
      if (out == nullptr)
        res = new std::complex<float>[Nfft]();
      else
        res = out;
      const float constant = -2.f * M_PI / float(Nfft);
      for (std::size_t k = 0; k < Nfft; ++k)
      {
        float factor = constant * float(k);
        for (std::size_t n = 0; n < Nfft; ++n)
          res[k] += src[n] * std::exp(factor * n * 1.if);
      }
      return res;
    }

    static float amplitude_cf(float (*window)(const std::size_t &, const std::size_t &),
                              const std::size_t &size)
    {
      float Aw = 0.f;
      for (std::size_t i = 0; i < size; ++i)
        Aw += window(i, size - 1);
      return float(size) / Aw;
    }

    static float energy_cf(float (*window)(const std::size_t &, const std::size_t &),
                           const std::size_t &size)
    {
      float Ew = 0.f;
      for (std::size_t i = 0; i < size; ++i)
      {
        float w = window(i, size - 1);
        Ew += w * w;
      }
      return sqrtf(float(size) / Ew);
    }
  };
} // namespace hmm

#endif // HMM_FFT_HPP_
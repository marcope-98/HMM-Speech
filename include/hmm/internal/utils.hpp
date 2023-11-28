#ifndef HMM_INTERNAL_UTILS_HPP_
#define HMM_INTERNAL_UTILS_HPP_

#include <cstdint>

namespace hmm
{
  struct utils
  {
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
  };
} // namespace hmm
#endif // HMM_INTERNAL_UTILS_HPP_
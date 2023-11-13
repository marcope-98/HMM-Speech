#ifndef HMM_SIMD_HPP_
#define HMM_SIMD_HPP_

#include <cstdint>

#include <immintrin.h>

namespace hmm
{
  struct simd
  {
    static std::size_t align(const std::size_t &value, const std::size_t &bytes) { return (value + (bytes - 1)) & ~(bytes - 1); }
    static void        copy(const float *const src, const std::size_t &size, float *const dst)
    {
      for (std::size_t i = 0; i < size; i += 4)
        _mm_store_ps(dst + i, _mm_load_ps(src + i));
    }

    static void copy(const std::uint8_t *const src, const std::size_t &size, std::uint8_t *const dst)
    {
      for (std::size_t i = 0; i < size; i += 16)
        _mm_store_si128((__m128i *)(dst + i), _mm_load_si128((__m128i *)(src + i)));
    }
  };
} // namespace hmm
#endif // HMM_SIMD_HPP_
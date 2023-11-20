#ifndef HMM_INTERNAL_HMMINTRIN_HPP_
#define HMM_INTERNAL_HMMINTRIN_HPP_

#include <cstdint>
#include <limits>

namespace hmm
{
  struct fops
  {
    inline static std::uint32_t f32_to_hex(const float &f32) { return *(std::uint32_t *)&f32; }
    inline static float         hex_to_f32(const std::uint32_t &hex) { return *(float *)&hex; }

    inline static float abs(const float &value)
    {
      std::uint32_t u32 = f32_to_hex(value) & 0x7FFFFFFFU;
      return hex_to_f32(u32);
    }

    inline static float max(const float &a, const float &b) { return a > b ? a : b; }
    inline static float min(const float &a, const float &b) { return a < b ? a : b; }

    inline static bool cmpgt(const float &a, const float &b)
    {
      const float epsilon = std::numeric_limits<float>::epsilon();
      const float fabs_a  = abs(a);
      const float fabs_b  = abs(b);
      return (a - b) > (max(fabs_a, fabs_b) * epsilon);
    }

    inline static bool cmplt(const float &a, const float &b)
    {
      const float epsilon = std::numeric_limits<float>::epsilon();
      const float fabs_a  = abs(a);
      const float fabs_b  = abs(b);
      return (b - a) > (max(fabs_a, fabs_b) * epsilon);
    }
    inline static bool cmpge(const float &a, const float &b) { return !cmplt(a, b); }
    inline static bool cmple(const float &a, const float &b) { return !cmpgt(a, b); }
    inline static bool cmpeq(const float &a, const float &b) { return !(cmplt(a, b) || cmpgt(a, b)); }
  };
} // namespace hmm

#endif // HMM_INTERNAL_HMMINTRIN_HPP_
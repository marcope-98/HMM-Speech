#ifndef HMM_INTERNAL_HMMINTRIN_HPP_
#define HMM_INTERNAL_HMMINTRIN_HPP_

#include <cstdint>
#include <limits>

namespace hmm
{
  struct fops
  {
  private:
    union cvt
    {
      float         f32;
      std::uint32_t u32{0u};
      cvt(){};
      explicit cvt(const float &value) : f32{value} {}
    };

  public:
    static float f32_to_hex(const float &f32)
    {
      cvt temp(f32);
      return temp.u32;
    }

    static float hex_to_f32(const std::uint32_t &hex)
    {
      cvt temp;
      temp.u32 = hex;
      return temp.f32;
    }

    static float abs(const float &value)
    {
      cvt temp(value);
      temp.u32 &= (0x7FFFFFFFU);
      return temp.f32;
    }

    static bool cmpgt(const float &a, const float &b)
    {
      const float epsilon = std::numeric_limits<float>::epsilon();
      const float fabs_a  = abs(a);
      const float fabs_b  = abs(b);
      return (a - b) > ((fabs_a < fabs_b ? fabs_b : fabs_a) * epsilon);
    }

    static bool cmplt(const float &a, const float &b)
    {
      const float epsilon = std::numeric_limits<float>::epsilon();
      const float fabs_a  = abs(a);
      const float fabs_b  = abs(b);
      return (b - a) > ((fabs_a < fabs_b ? fabs_b : fabs_a) * epsilon);
    }
    static bool cmpge(const float &a, const float &b) { return !cmplt(a, b); }
    static bool cmple(const float &a, const float &b) { return !cmpgt(a, b); }
    static bool cmpeq(const float &a, const float &b) { return !(cmplt(a, b) || cmpgt(a, b)); }
  };
} // namespace hmm

#endif // HMM_INTERNAL_HMMINTRIN_HPP_
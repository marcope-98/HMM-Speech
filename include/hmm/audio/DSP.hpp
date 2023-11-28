#ifndef HMM_AUDIO_DSP_HPP_
#define HMM_AUDIO_DSP_HPP_

#include <cstdint>
#include <utility>

#include "hmm/internal/Ptr.hpp"
#include "hmm/internal/utils.hpp"

namespace hmm
{
  struct MeanRemover
  {
    auto execute(Ptr p)
    {
      /* 
        NOTE: in reality the quantity to subtract is (mean(src.*window)/mean(window))
              only under these conditions the fft[0] is exactly zero
              still this is a fair approximation
      */
      float *     data = (float *)p.data;
      std::size_t size = p.size;
      // compute mean
      float mean{0.f};
      for (std::size_t i = 0; i < size; ++i)
        mean += data[i];
      mean /= float(size);

      // subtract mean
      for (std::size_t i = 0; i < size; ++i)
        data[i] -= mean;

      return p;
    }
  };

  struct ZeroPadding
  {
    auto execute(Ptr p)
    {
      const std::size_t size = utils::clp2(p.size);
      float *           src  = (float *)p.data;
      if (size * sizeof(float) <= p.capacity)
      {
        // make sure the remaining size is zero padded
        for (std::size_t i = p.size; i < size; ++i)
          src[i] = 0.f;
        // adjust size to be power of two
        p.size = size;
        // return pointer
        return p;
      }

      // reallocate
      void *temp = malloc(size * sizeof(float));

      // copy
      float *dst = (float *)temp;
      for (std::size_t i = 0; i < p.size; ++i)
        dst[i] = src[i];
      for (std::size_t i = p.size; i < size; ++i)
        dst[i] = 0.f;

      // swap
      std::swap(temp, p.data);
      p.size     = size;
      p.capacity = size * sizeof(float);
      free(temp);

      return p;
    }
  };

  struct Preemphesis
  {
    auto execute(Ptr p)
    {
      const std::size_t size = p.size;
      float *           dst  = (float *)p.data;
      for (std::size_t i = size - 1; i > 0; --i)
        dst[i] -= 0.95f * dst[i - 1];
      return p;
    }
  };

} // namespace hmm

#endif // HMM_AUDIO_DSP_HPP_
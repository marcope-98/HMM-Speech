#ifndef HMM_AUDIO_AUDIO_HPP_
#define HMM_AUDIO_AUDIO_HPP_
#include <cassert>
#include <complex>
#include <cstdint>
#include <utility>

using namespace std::complex_literals;
#include "hmm/simd.hpp"

namespace hmm
{
  class Audio
  {
  public:
    Audio() = default;
    explicit Audio(const char *filename);
    Audio(const Audio &other) { *this = other; }
    Audio(Audio &&other) noexcept { *this = std::move(other); }
    ~Audio();

    Audio &operator=(const Audio &other);
    Audio &operator=(Audio &&other) noexcept;
    float  operator()(const std::size_t &sample, const std::size_t &channel) const { return this->p_data[channel + sample * this->d_channels]; }

  public:
    std::size_t size() const { return this->d_samples * this->d_channels; }
    std::size_t capacity() const { return simd::align(this->size(), 4); }

    float at(const std::size_t &sample, const std::size_t &channel) const
    {
      assert(channel < this->d_channels);
      assert(sample < this->d_samples);
      return this->p_data[channel + sample * this->d_channels];
    }

    void preemphesis();
    void plot() const;
    void fft_plot(float (*window)(float, float)) const;

  private:
    std::complex<float> *dft(const float *const src, const std::size_t &N) const;
    std::size_t          clp2(std::size_t x) const;

    void deallocate_all();
    void reallocate(const std::size_t &size);

  private:
    std::size_t d_samples{0};
    std::size_t d_channels{0};
    std::size_t d_sampling_frequency{0};
    float *     p_data{nullptr};
  };
} // namespace hmm
#endif
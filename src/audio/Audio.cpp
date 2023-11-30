#include "hmm/audio/Audio.hpp"

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <iostream>

// hmm/io
#ifndef FPM_WAV_IMPLEMENTATION
#define FPM_WAV_IMPLEMENTATION
#include "hmm/io/fpm_wav.h"
#endif
#include "hmm/io/graphics.hpp"

// hmm/audio
#include "hmm/audio/DSP.hpp"
#include "hmm/audio/FFT.hpp"
#include "hmm/audio/Windows.hpp"

// hmm/internal
#include "hmm/internal/Pipeline.hpp"
#include "hmm/internal/Ptr.hpp"

hmm::Audio::Audio(const char *filename) : Audio()
{
  this->p_data = fpm_wav_load(filename,
                              &this->d_samples,
                              &this->d_channels,
                              &this->d_sampling_frequency);
  if (this->p_data == NULL)
  {
    std::cerr << "[ERROR] could not read filename " << filename << "\n";
    exit(1);
  }
}

hmm::Audio &hmm::Audio::operator=(const Audio &other)
{
  if (this->size() != other.size()) this->reallocate(other.size());
  this->d_samples            = other.d_samples;
  this->d_channels           = other.d_channels;
  this->d_sampling_frequency = other.d_sampling_frequency;

  // TODO: optimize
  for (std::size_t i = 0; i < this->size(); ++i)
    this->p_data[i] = other.p_data[i];
  return *this;
}

hmm::Audio &hmm::Audio::operator=(Audio &&other) noexcept
{
  std::swap(this->p_data, other.p_data);
  std::swap(this->d_samples, other.d_samples);
  std::swap(this->d_channels, other.d_channels);
  std::swap(this->d_sampling_frequency, other.d_sampling_frequency);
  return *this;
}

hmm::Audio::~Audio() { this->deallocate_all(); }

void hmm::Audio::deallocate_all()
{
  delete[] this->p_data;
  this->d_samples            = 0;
  this->d_channels           = 0;
  this->d_sampling_frequency = 0;
}

void hmm::Audio::reallocate(const std::size_t &size)
{
  this->deallocate_all();
  this->p_data = new float[size]();
}

void hmm::Audio::plot() const
{
  const float delta = 1.f / float(this->d_sampling_frequency);
  float *     time  = graphics::generate_axis(0.f, delta, this->d_samples);

  graphics::open_window();
  graphics::plot(time, this->p_data, this->d_samples, "Time [s]", "Amplitude [-]");
  graphics::wait();

  delete[] time;
}

void hmm::Audio::fft() const
{
  Ptr ptr;
  ptr.allocate<float>(this->d_samples);
  ptr.copy<float>(this->p_data, this->d_samples);

  const float Aw = Hann::correct(ptr.size) / float(ptr.size);

  // MeanRemover -> Hann -> ZeroPadding -> ToComplex -> CooleyTukey -> Magnitude
  Pipeline<MeanRemover, Hann, ZeroPadding, ToComplex, CooleyTukey, Magnitude> pipeline;
  // Execute Pipeline
  ptr = pipeline.execute(ptr);

  // Prepare output
  const std::size_t half_spectrum = ptr.size / 2 + 1;
  float *           magnitude     = ptr.cast<float>();

  const float delta = float(this->d_sampling_frequency) / float(ptr.size);
  float *     freq  = graphics::generate_axis(0.f, delta, half_spectrum);

  magnitude[0] *= Aw;
  for (std::size_t i = 1; i < half_spectrum; ++i)
    magnitude[i] *= 2.f * Aw;

  // plot
  graphics::open_window();
  graphics::plot(freq, magnitude, half_spectrum, "Frequency [Hz]", "Magnitude [-]");
  graphics::wait();

  // clean up
  delete[] freq;
  ptr.deallocate();
}

void hmm::Audio::spectrogram(const std::size_t &length) const
{
  const std::size_t overlap = length / 2; // NOTE: we only support 50% overlap
  /*
    let 
    n = this->d_samples
    d = overlap
    width = ceil(n / d)            - 1
          = floor((n + d - 1) / d) - 1
    since unsigned integer division rounds towards zero it is equivalent to floor
  */
  const std::size_t width  = ((this->d_samples + overlap - 1) / overlap) - 1;
  const std::size_t height = utils::clp2(length);

  // we basically allocate enough space to zero pad and convert to complex
  Ptr ptr;
  ptr.allocate<std::complex<float>>(width * height);

  // copy overlapping data into ptr
  std::size_t stride = height * sizeof(std::complex<float>) / sizeof(float); // to get next block
  float *     src    = this->p_data;
  float *     dst    = ptr.cast<float>();
  for (std::size_t col = 0; col < width - 1; ++col)
  {
    for (std::size_t row = 0; row < height; ++row)
      dst[row] = src[row];
    dst += stride;
    src += overlap;
  }
  // handle special case where this->d_samples % overlap != 0
  // this extends the last block to reach the size of length specified as input
  float *end1 = this->p_data + this->d_samples;
  float *end2 = dst + length;
  // clang-format off
  while (src != end1) *dst++ = *src++;
  while (dst != end2) *dst++ = 0.f;
  // clang-format on

  Pipeline<ZeroPadding, Hamming, ToComplex, CooleyTukey, Magnitude> pipeline;
  Ptr                                                               temp;
  temp.capacity           = height * sizeof(std::complex<float>);
  stride                  = height;
  std::complex<float> *it = ptr.cast<std::complex<float>>();
  for (std::size_t col = 0; col < width; ++col)
  {
    temp.data = it;
    temp.size = length; // this will change every time so we need to reassign it at each iteration
    temp      = pipeline.execute(temp);
    it += stride;
  }
  const float Aw = Hamming::correct(length) / float(length);

  // Graphics:
  const std::size_t half_spectrum = height / 2 + 1;
  float *           z             = new float[half_spectrum * width]();
  stride                          = height * sizeof(std::complex<float>) / sizeof(float);
  src                             = ptr.cast<float>();
  for (std::size_t col = 0; col < width; ++col)
  {
    z[col] = 20.f * log10f(Aw * src[0]);
    for (std::size_t i = 0; i < half_spectrum; ++i)
    {
      std::size_t index_dst = i * width + col;
      z[index_dst]          = 20.f * log10f(2.f * Aw * src[i]);
    }
    src += stride;
  }

  // src = ptr.cast<float>();
  // for (std::size_t col = 0; col < width; ++col)
  // {
  //   for (std::size_t i = 0; i < half_spectrum; ++i)
  //     std::cout << src[i] << "\n";
  //   src += stride;
  // }

  ptr.deallocate();

  float delta;
  delta    = float(this->d_samples) / float(width * this->d_sampling_frequency);
  float *x = graphics::generate_axis(0.5f * delta, delta, width);
  delta    = float(this->d_sampling_frequency) / float(height);
  float *y = graphics::generate_axis(0.f, delta, half_spectrum);

  graphics::open_window();
  graphics::matrix(x, y, z,
                   half_spectrum, width,
                   "Time [Block #]", "Frequency [Hz]");
  graphics::wait();
  delete[] x;
  delete[] y;
  delete[] z;
}

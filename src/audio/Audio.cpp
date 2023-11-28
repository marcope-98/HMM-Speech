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
#if 0
  this->d_channels           = 1;
  this->d_sampling_frequency = 1000;
  this->d_samples            = 1500;
  this->p_data               = new float[this->d_samples];
  const float dt             = 1.f / float(this->d_sampling_frequency);
  float       t              = 0.f;
  for (std::size_t i = 0; i < this->d_samples; ++i)
  {
    this->p_data[i] = 0.8f + 0.7f * sinf(2.f * M_PI * 50.f * t) + sinf(2.f * M_PI * 120.f * t);
    t += dt;
  }
#else
  this->p_data = fpm_wav_load(filename,
                              &this->d_samples,
                              &this->d_channels,
                              &this->d_sampling_frequency);
  if (this->p_data == NULL)
  {
    std::cerr << "[ERROR] could not read filename " << filename << "\n";
    exit(1);
  }
#endif
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
  ptr.data     = malloc(this->d_samples * sizeof(float));
  ptr.size     = this->d_samples;
  ptr.capacity = this->d_samples * sizeof(float);

  float *src = (float *)ptr.data;
  for (std::size_t i = 0; i < this->d_samples; ++i)
    src[i] = this->p_data[i];

  const float Aw = Hann::correct(ptr.size) / float(ptr.size);

  // MeanRemover -> Hann -> ZeroPadding -> ToComplex -> CooleyTukey -> Magnitude
  Pipeline<MeanRemover, Hann, ZeroPadding, ToComplex, CooleyTukey, Magnitude> pipeline;
  // Execute Pipeline
  ptr = pipeline.execute(ptr);

  // Prepare output
  const std::size_t half_spectrum = ptr.size / 2 + 1;
  float *           magnitude     = (float *)ptr.data;

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
  free(ptr.data);
}

void hmm::Audio::spectrogram() const
{
  const std::size_t window_length = 256;                                                            // NOTE: Length of each time block
  const std::size_t overlap       = window_length / 2;                                              // NOTE: 50% overlap
  const std::size_t Nfft          = utils::round_up_to_multiple_of(this->d_samples, window_length); // NOTE: Need to do this manually
  const std::size_t blocks        = (Nfft / overlap) - 1;                                           // NOTE: this assumes 50% overlap
  const float       Aw            = Hann::correct(window_length) / float(window_length);            // NOTE: Hann window amplitude correction factor

  // Initialize pipeline: Hann -> ToComplex -> FFT -> Magnitude
  Pipeline<Hann, ToComplex, CooleyTukey, Magnitude> pipeline;

  // main ptr of pipeline
  Ptr ptr;
  ptr.size     = blocks * window_length;
  ptr.capacity = ptr.size * sizeof(std::complex<float>); // make sure there's enough space for std::complex<float> data type (avoid reallocation)
  ptr.data     = malloc(ptr.capacity);

  // copy overlapping data into ptr
  std::size_t stride = window_length * sizeof(std::complex<float>) / sizeof(float);
  float *     src    = this->p_data;
  float *     dst    = (float *)ptr.data;
  for (std::size_t block = 0; block < blocks - 1; ++block)
  {
    for (std::size_t sample = 0; sample < window_length; ++sample)
      dst[sample] = src[sample];
    dst += stride;
    src += overlap;
  }
  // handle special case where Nfft != this->d_samples
  for (std::size_t sample = 0; sample < Nfft - this->d_samples; ++sample)
    dst[sample] = src[sample];
  for (std::size_t sample = Nfft - this->d_samples; sample < window_length; ++sample)
    dst[sample] = 0.f;

  // Pipeline: create temporary Ptr for each block
  Ptr temp;
  temp.size                             = window_length;
  temp.capacity                         = window_length * sizeof(std::complex<float>);
  stride                                = window_length;
  std::complex<float> *ptr_data_complex = (std::complex<float> *)ptr.data;
  for (std::size_t block = 0; block < blocks; ++block)
  {
    temp.data = ptr_data_complex;
    temp      = pipeline.execute(temp);
    ptr_data_complex += stride;
  }

  // Graphics:
  const std::size_t half_spectrum = window_length / 2 + 1;
  float *           z             = new float[half_spectrum * blocks]();
  stride                          = window_length * sizeof(std::complex<float>) / sizeof(float);
  src                             = (float *)ptr.data;
  for (std::size_t block = 0; block < blocks; ++block)
  {
    z[block] = 20.f * log10f(Aw * src[0]);
    for (std::size_t i = 0; i < half_spectrum; ++i)
    {
      std::size_t index_dst = i * blocks + block;
      z[index_dst]          = 20.f * log10f(2.f * Aw * src[i]);
    }
    src += stride;
  }

  free(ptr.data);

  float delta;
  delta    = float(this->d_samples) / float(blocks * this->d_sampling_frequency);
  float *x = graphics::generate_axis(0.5f * delta, delta, blocks);
  delta    = float(this->d_sampling_frequency) / float(window_length);
  float *y = graphics::generate_axis(0.f, delta, half_spectrum);

  graphics::open_window();
  graphics::matrix(x, y, z,
                   half_spectrum, blocks,
                   "Time [Block #]", "Frequency [Hz]");
  graphics::wait();
  delete[] x;
  delete[] y;
  delete[] z;
}

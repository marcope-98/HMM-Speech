#include "hmm/audio/Audio.hpp"

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <iostream>

#include "hmm/audio/fft.hpp"
#include "hmm/io/graphics.hpp"

#ifndef FPM_WAV_IMPLEMENTATION
#define FPM_WAV_IMPLEMENTATION
#include "hmm/io/fpm_wav.h"
#endif

hmm::Audio::Audio(const char *filename) : Audio()
{
#if 1
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

void hmm::Audio::preemphesis()
{
  fft::preemphesis(this->p_data, this->d_samples);
}

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
  float *time = new float[this->d_samples]();

  const float dt = 1.f / float(this->d_sampling_frequency);
  float       t  = 0.f;
  for (std::size_t i = 0; i < this->d_samples; ++i)
  {
    t += dt;
    time[i] = t;
  }

  graphics::open_window();
  graphics::plot(time, this->p_data, this->d_samples, "Time [s]", "Amplitude [-]");
  graphics::wait();

  delete[] time;
}

// TODO: consider making this a templated function to accept a lambda instead
void hmm::Audio::fft_plot(float (*window)(const std::size_t &, const std::size_t &)) const
{
  // next power of 2
  const std::size_t Nfft = fft::clp2(this->d_samples);
  // zero pad to next power of 2
  float *src = fft::zero_padding(this->p_data, this->d_samples, Nfft);
  // remove mean of signal:
  /* 
  NOTE: in reality the quantity to subtract is (mean(src.*window)/mean(window))
        only under these conditions the fft[0] is exactly zero
        still this is a fair approximation
  */
  float mean = fft::mean(src, this->d_samples);
  // TODO: to optimize
  for (std::size_t i = 0; i < this->d_samples; ++i)
    src[i] -= mean;

  // apply window to signal
  fft::apply_window(src, this->d_samples, window);
  // compute amplitude correction factor
  float Aw = fft::amplitude_cf(window, this->d_samples);
  Aw /= float(this->d_samples);

  // Cooley-Tukey fft algorithm
  std::complex<float> *temp = fft::cooley_tukey(src, Nfft);
  delete[] src;

  // plot dft magnitude
  const std::size_t half_spectrum = Nfft / 2 + 1;
  float *           freq          = new float[half_spectrum]();
  float *           magnitude     = new float[half_spectrum]();
  const float       delta         = float(this->d_sampling_frequency) / float(Nfft);
  float             frequency     = 0.f;

  freq[0]      = frequency;
  magnitude[0] = Aw * std::abs(temp[0]);
  for (std::size_t i = 1; i < half_spectrum; ++i)
  {
    frequency += delta;
    freq[i]      = frequency;
    magnitude[i] = 2.f * Aw * std::abs(temp[i]);
  }
  delete[] temp;

  graphics::open_window();
  graphics::plot(freq, magnitude, half_spectrum, "Frequency [Hz]", "Magnitude [-]");
  graphics::wait();

  delete[] freq;
  delete[] magnitude;
}

void hmm::Audio::spectrogram(float (*window)(const std::size_t &, const std::size_t &)) const
{
  /*
  To be precise assuming 50% overlap I need: 
   - a number of samples is at least window size
   - a number of samples multiple of overlap
  */
  const std::size_t window_length = 256;
  const std::size_t overlap       = window_length / 2; // 50% overlap
  const std::size_t Nfft          = fft::round_up_to_multiple_of(this->d_samples, window_length);
  float *           src           = fft::zero_padding(this->p_data, this->d_samples, Nfft);

  const std::size_t    blocks = (Nfft / overlap) - 1; // NOTE: This assumes 50% overlap
  std::complex<float> *spectr = new std::complex<float>[window_length * blocks]();
  float *              tmp    = new float[window_length]();

  std::size_t offset = 0;
  const float Aw     = fft::amplitude_cf(window, window_length) / float(window_length);

  for (std::size_t block = 0; block < blocks; ++block)
  {
    // copy src to tmp
    for (std::size_t i = 0; i < window_length; ++i)
      tmp[i] = src[i + offset];
    // apply window to tmp
    fft::apply_window(tmp, window_length, window);
    // compute fft on tmp
    fft::cooley_tukey(tmp, window_length, spectr + block * window_length);
    offset += overlap;
  }
  delete[] tmp;
  delete[] src;

  // compute spectrogram + transpose
  const std::size_t half_spectrum = window_length / 2 + 1;
  float *           z             = new float[half_spectrum * blocks]();
  for (std::size_t block = 0; block < blocks; ++block)
  {
    z[block] = 20.f * log10f(Aw * std::abs(spectr[block * window_length]));
    for (std::size_t i = 1; i < half_spectrum; ++i)
    {
      std::size_t index_src = i + block * window_length;
      std::size_t index_dst = i * blocks + block;
      z[index_dst]          = 20.f * log10f(2.f * Aw * std::abs(spectr[index_src]));
    }
  }
  delete[] spectr;

  float *     x  = new float[blocks]();
  const float dt = float(this->d_samples) / float(blocks * this->d_sampling_frequency);
  float       t  = dt * 0.5f;
  for (std::size_t i = 0; i < blocks; ++i)
  {
    x[i] = t;
    t += dt;
  }

  float *     y         = new float[half_spectrum]();
  const float delta     = float(this->d_sampling_frequency) / float(window_length);
  float       frequency = 0.f;
  for (std::size_t i = 0; i < half_spectrum; ++i)
  {
    y[i] = frequency;
    frequency += delta;
  }

  graphics::open_window();
  graphics::matrix(x, y, z,
                   half_spectrum, blocks,
                   "Time [Block #]", "Frequency [Hz]");
  graphics::wait();
  delete[] x;
  delete[] y;
  delete[] z;
}

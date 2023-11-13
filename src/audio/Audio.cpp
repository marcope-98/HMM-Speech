#include "hmm/audio/Audio.hpp"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include "hmm/graphics.hpp"

#ifndef FPM_WAV_IMPLEMENTATION
#define FPM_WAV_IMPLEMENTATION
#include "hmm/io/fpm_wav.h"
#endif

hmm::Audio::Audio(const char *filename) : Audio()
{
#if 0
  this->d_channels           = 1;
  this->d_sampling_frequency = 1000;
  this->d_samples            = 1500;
  this->p_data               = new float[this->d_samples]();
  float t                    = 0.f;
  float dt                   = 1.f / float(this->d_sampling_frequency);
  for (std::size_t i = 0; i < this->d_samples; ++i)
  {
    this->p_data[i] = 0.8f +
                      0.7f * sinf(2.f * M_PI * 50.f * t) +
                      1.f * sinf(2.f * M_PI * 120.f * t);
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
  if (this->capacity() != other.capacity()) this->reallocate(other.capacity());
  this->d_samples            = other.d_samples;
  this->d_channels           = other.d_channels;
  this->d_sampling_frequency = other.d_sampling_frequency;
  simd::copy(other.p_data, this->size(), this->p_data);
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
  for (std::size_t i = this->d_samples - 1; i > 0; --i)
    this->p_data[i] -= 0.95f * this->p_data[i - 1];
}

void hmm::Audio::deallocate_all()
{
  delete[] this->p_data;
  this->d_samples            = 0;
  this->d_channels           = 0;
  this->d_sampling_frequency = 0;
}

void hmm::Audio::reallocate(const std::size_t &capacity)
{
  this->deallocate_all();
  this->p_data = new float[capacity]();
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

// TODO: move this in a separate file
std::size_t hmm::Audio::clp2(std::size_t x) const
{
  x -= 1;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;
  return x + 1;
}

// TODO: this should be in a separate header
void hmm::Audio::fft_plot(float (*window)(float, float)) const
{
  // NOTE: zero padding
  std::size_t N_padded = this->clp2(this->d_samples); // next power of 2
  float *     src      = new float[N_padded]();
  for (std::size_t i = 0; i < this->d_samples; ++i)
    src[i] = this->p_data[i];

  // NOTE: Windowing + computation of amplitude correcting factor
  // TODO: make it a separate function call
  float Aw = 0.f; // Amplitude correcting factor
  for (std::size_t i = 0; i < this->d_samples; ++i)
  {
    // TODO: value is too generic, find a better naming
    float value = window(float(i), float(this->d_samples - 1));
    Aw += value;
    src[i] *= value;
  }
  Aw = 1.f / Aw;

  // NOTE: Discrete Fourier Transform
  // TODO: separate this in different translation unit
  // TODO: consider including window() as input parameter
  std::complex<float> *temp = this->dft(src, this->d_samples);

  delete[] src;

  // plot dft magnitude
  const std::size_t half_spectrum = N_padded / 2 + 1;
  float *           freq          = new float[half_spectrum]();
  float *           magnitude     = new float[half_spectrum]();
  const float       delta         = float(this->d_sampling_frequency) / float(N_padded);

  // TODO: instead of doing evil bool->float conversion consider initializing the first element and loop trough the rest
  for (std::size_t i = 0; i < half_spectrum; ++i)
  {
    freq[i] = float(i) * delta;
    // NOTE: (2.f - float(i == 0)) => double every magnitude except for the mean
    magnitude[i] = Aw * (2.f - float(i == 0)) * std::abs(temp[i]);
  }
  delete[] temp;

  graphics::open_window();
  graphics::plot(freq, magnitude, half_spectrum, "Frequency [Hz]", "Magnitude [-]");
  graphics::wait();
  delete[] freq;
  delete[] magnitude;
}

// TODO: this should be in a separate header
std::complex<float> *hmm::Audio::dft(const float *const src, const std::size_t &N) const
{
  const std::size_t    N_padded = this->clp2(N);
  std::complex<float> *answer   = new std::complex<float>[N_padded]();

  // exp(- i 2 pi k n / N)
  const float constant = -2.f * M_PI / float(N);
  for (std::size_t k = 0; k < N_padded; ++k)
  {
    float factor = constant * float(k);
    // answer[k]    = sum(0.f, 0.f);
    for (std::size_t n = 0; n < N; ++n)
      answer[k] += src[n] * std::exp(factor * n * 1.if);
  }
  return answer;
}

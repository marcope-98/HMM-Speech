#include "hmm/audio/Audio.hpp"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include <RInside.h>
#include <Rcpp.h>

#ifndef FPM_WAV_IMPLEMENTATION
#define FPM_WAV_IMPLEMENTATION
#include "hmm/io/fpm_wav.h"
#endif

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
  static RInside R;
  Rcpp::Function x11("x11");
  Rcpp::Function plot("plot");
  x11();
  Rcpp::NumericVector x;
  Rcpp::NumericVector y;
  const float         increment = 1.f / float(this->d_sampling_frequency);
  float               time      = 0.f;
  for (std::size_t i = 0; i < this->d_samples; ++i)
  {
    time += increment;
    x.push_back(time);
    y.push_back(this->p_data[i]);
  }

  plot(Rcpp::Named("x")    = x,
       Rcpp::Named("y")    = y,
       Rcpp::Named("type") = "l",
       Rcpp::Named("xlab") = "Time [s]",
       Rcpp::Named("ylab") = "Amplitude [-]");
  R.parseEval("while(names(dev.cur()) != 'null device') Sys.sleep(1)");
}

std::size_t hmm::Audio::clp2(std::size_t x)
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

void hmm::Audio::fft_plot(float (*window)(float, float)) const
{
  // zero padding
  std::size_t N_padded = clp2(this->d_samples); // next power of 2
  float *     src      = new float[N_padded]();
  for (std::size_t i = 0; i < this->d_samples; ++i)
    src[i] = this->p_data[i];

  float Aw = 0.f; // Amplitude correcting factor
  for (std::size_t i = 0; i < this->d_samples; ++i)
  {
    float value = window(float(i), float(this->d_samples - 1));
    Aw += value;
    src[i] *= value;
  }
  // Aw = 1.f / mean(window) = N / sum_i=0^N window(i, N - 1);
  Aw = float(this->d_samples) / Aw;

  // compute zero-padded fft with window
  std::complex<float> *temp = this->dft(src, this->d_samples);

  // plot dft magnitude
  static RInside      R;
  Rcpp::Function      x11("x11");
  Rcpp::Function      plot("plot");
  Rcpp::NumericVector x;
  Rcpp::NumericVector y;
  float               delta = float(this->d_sampling_frequency) / float(N_padded);
  for (std::size_t i = 0; i < N_padded / 2 + 1; ++i)
  {
    x.push_back(float(i) * delta);
    y.push_back(Aw * (2.f - float(i == 0)) * std::abs(temp[i]) / float(this->d_samples));
  }

  delete[] temp;
  delete[] src;

  x11();
  plot(Rcpp::Named("x")    = x,
       Rcpp::Named("y")    = y,
       Rcpp::Named("type") = "l",
       Rcpp::Named("xlab") = "Frequency [Hz]",
       Rcpp::Named("ylab") = "log(Magnitude) [-]");
  R.parseEval("while(names(dev.cur()) != 'null device') Sys.sleep(1)");
}

std::complex<float> *hmm::Audio::dft(const float *const src, const std::size_t &N) const
{
  const std::size_t    N_padded = clp2(N);
  std::complex<float> *answer   = new std::complex<float>[N_padded]();

  // exp(- i 2 pi k n / N)
  const float constant = -2.f * M_PI / float(N);
  for (std::size_t k = 0; k < N_padded; ++k)
  {
    float factor = constant * float(k);
    answer[k]    = sum(0.f, 0.f);
    for (std::size_t n = 0; n < N; ++n)
      answer[k] += src[n] * std::exp(factor * n * 1.if);
  }
  return answer;
}

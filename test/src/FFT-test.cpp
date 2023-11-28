#include <gtest/gtest.h>

#include "FFT-test.hpp"

#include "hmm/audio/FFT.hpp"
#include "hmm/internal/Pipeline.hpp"
#include "hmm/internal/Ptr.hpp"

TEST(hmmTest, ToComplex)
{
  hmm::Ptr                      ptr;
  hmm::Pipeline<hmm::ToComplex> pipeline;
  const std::size_t             size = 1000;
  ptr.size                           = size;
  ptr.data                           = malloc(size * sizeof(float));

  float *src = (float *)ptr.data;
  for (std::size_t i = 0; i < size; ++i)
    src[i] = 20.f;

  ptr                       = pipeline.execute(ptr);
  std::complex<float> *temp = (std::complex<float> *)ptr.data;
  for (std::size_t i = 0; i < size; ++i)
  {
    ASSERT_FLOAT_EQ(temp[i].real(), 20.f);
    ASSERT_FLOAT_EQ(temp[i].imag(), 0.f);
  }

  // ToComplex no allocation
  src = (float *)ptr.data;
  for (std::size_t i = 0; i < ptr.size; ++i)
    src[i] = 0.f;
  for (std::size_t i = 0; i < size; ++i)
    src[i] = 20.f;

  ptr  = pipeline.execute(ptr);
  temp = (std::complex<float> *)ptr.data;
  ASSERT_EQ(ptr.size, 1000);
  ASSERT_EQ(ptr.capacity, 1000 * sizeof(std::complex<float>));
  ASSERT_TRUE((float *)(ptr.data) == src);
  for (std::size_t i = 0; i < size; ++i)
  {
    ASSERT_FLOAT_EQ(temp[i].real(), 20.f);
    ASSERT_FLOAT_EQ(temp[i].imag(), 0.f);
  }

  free(ptr.data);
}

TEST(hmmTest, CooleyTukey)
{
  hmm::Ptr                        ptr1, ptr2;
  hmm::Pipeline<hmm::DFT>         pipeline_dft;
  hmm::Pipeline<hmm::CooleyTukey> pipeline_ct;

  const std::size_t size = 4096;
  const float       N    = float(size);
  ptr1.size              = size;
  ptr1.data              = malloc(size * sizeof(std::complex<float>));
  ptr2.size              = size;
  ptr2.data              = malloc(size * sizeof(std::complex<float>));

  std::complex<float> *ref = (std::complex<float> *)ptr1.data;
  std::complex<float> *src = (std::complex<float> *)ptr2.data;
  // Test 1: impulse response
  FFT_test::fill(ref, size, &FFT_test::impulse_response);
  FFT_test::fill(src, size, &FFT_test::impulse_response);
  ptr1 = pipeline_dft.execute(ptr1);
  ptr2 = pipeline_ct.execute(ptr2);
  ref  = (std::complex<float> *)ptr1.data;
  src  = (std::complex<float> *)ptr2.data;

  for (std::size_t i = 0; i < size; ++i)
    ASSERT_NEAR(std::abs(ref[i]) / N, std::abs(src[i]) / N, 1e-3);

  // Test 2: complex sinusoid
  FFT_test::fill(ref, size, &FFT_test::complex_sinusoid);
  FFT_test::fill(src, size, &FFT_test::complex_sinusoid);
  ptr1 = pipeline_dft.execute(ptr1);
  ptr2 = pipeline_ct.execute(ptr2);
  ref  = (std::complex<float> *)ptr1.data;
  src  = (std::complex<float> *)ptr2.data;

  for (std::size_t i = 0; i < size; ++i)
    ASSERT_NEAR(std::abs(ref[i]) / N, std::abs(src[i]) / N, 1e-3);

  // Test 3: rectangle
  FFT_test::fill(ref, size, &FFT_test::rectangle);
  FFT_test::fill(src, size, &FFT_test::rectangle);
  ptr1 = pipeline_dft.execute(ptr1);
  ptr2 = pipeline_ct.execute(ptr2);
  ref  = (std::complex<float> *)ptr1.data;
  src  = (std::complex<float> *)ptr2.data;

  for (std::size_t i = 0; i < size; ++i)
    ASSERT_NEAR(std::abs(ref[i]) / N, std::abs(src[i]) / N, 1e-3);

  free(ptr1.data);
  free(ptr2.data);
}
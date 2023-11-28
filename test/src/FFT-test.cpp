#include <gtest/gtest.h>

#include "FFT-test.hpp"

#include "hmm/audio/FFT.hpp"
#include "hmm/internal/Pipeline.hpp"
#include "hmm/internal/Ptr.hpp"

TEST(hmmTest, ToComplex)
{
  const std::size_t             size = 1000;
  hmm::Pipeline<hmm::ToComplex> pipeline;
  hmm::Ptr                      ptr;
  ptr.allocate<float>(size);

  float *src = ptr.cast<float>();
  for (std::size_t i = 0; i < size; ++i)
    src[i] = 20.f;

  ptr = pipeline.execute(ptr);
  ASSERT_EQ(ptr.size, 1000);
  ASSERT_EQ(ptr.capacity, 1000 * sizeof(std::complex<float>));
  std::complex<float> *temp = ptr.cast<std::complex<float>>();
  for (std::size_t i = 0; i < size; ++i)
  {
    ASSERT_FLOAT_EQ(temp[i].real(), 20.f);
    ASSERT_FLOAT_EQ(temp[i].imag(), 0.f);
  }

  // ToComplex no allocation
  src = ptr.cast<float>();
  for (std::size_t i = 0; i < ptr.size; ++i)
    src[i] = 0.f;
  for (std::size_t i = 0; i < size; ++i)
    src[i] = 20.f;

  ptr  = pipeline.execute(ptr);
  temp = ptr.cast<std::complex<float>>();
  ASSERT_EQ(ptr.size, 1000);
  ASSERT_EQ(ptr.capacity, 1000 * sizeof(std::complex<float>));
  ASSERT_TRUE(ptr.cast<float>() == src);
  for (std::size_t i = 0; i < size; ++i)
  {
    ASSERT_FLOAT_EQ(temp[i].real(), 20.f);
    ASSERT_FLOAT_EQ(temp[i].imag(), 0.f);
  }
  ptr.deallocate();
}

TEST(hmmTest, CooleyTukey)
{
  hmm::Ptr                        ptr1, ptr2;
  hmm::Pipeline<hmm::DFT>         pipeline_dft;
  hmm::Pipeline<hmm::CooleyTukey> pipeline_ct;

  const std::size_t size = 4096;
  const float       N    = float(size);
  ptr1.allocate<std::complex<float>>(size);
  ptr2.allocate<std::complex<float>>(size);
  std::complex<float> *ref = ptr1.cast<std::complex<float>>();
  std::complex<float> *src = ptr2.cast<std::complex<float>>();

  // Test 1: impulse response
  FFT_test::fill(ref, size, &FFT_test::impulse_response);
  FFT_test::fill(src, size, &FFT_test::impulse_response);
  ptr1 = pipeline_dft.execute(ptr1);
  ptr2 = pipeline_ct.execute(ptr2);
  ref  = ptr1.cast<std::complex<float>>();
  for (std::size_t i = 0; i < size; ++i)
    ASSERT_NEAR(std::abs(ref[i]) / N, std::abs(src[i]) / N, 1e-3);

  // Test 2: complex sinusoid
  FFT_test::fill(ref, size, &FFT_test::complex_sinusoid);
  FFT_test::fill(src, size, &FFT_test::complex_sinusoid);
  ptr1 = pipeline_dft.execute(ptr1);
  ptr2 = pipeline_ct.execute(ptr2);
  ref  = ptr1.cast<std::complex<float>>();

  for (std::size_t i = 0; i < size; ++i)
    ASSERT_NEAR(std::abs(ref[i]) / N, std::abs(src[i]) / N, 1e-3);

  // Test 3: rectangle
  FFT_test::fill(ref, size, &FFT_test::rectangle);
  FFT_test::fill(src, size, &FFT_test::rectangle);
  ptr1 = pipeline_dft.execute(ptr1);
  ptr2 = pipeline_ct.execute(ptr2);
  ref  = ptr1.cast<std::complex<float>>();

  for (std::size_t i = 0; i < size; ++i)
    ASSERT_NEAR(std::abs(ref[i]) / N, std::abs(src[i]) / N, 1e-3);

  ptr1.deallocate();
  ptr2.deallocate();
}
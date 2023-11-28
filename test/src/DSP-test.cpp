#include <gtest/gtest.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdint>

#include "hmm/audio/DSP.hpp"
#include "hmm/internal/Pipeline.hpp"
#include "hmm/internal/Ptr.hpp"

TEST(hmmTest, MeanRemover)
{
  // initialize
  const std::size_t               size = 1000;
  hmm::Pipeline<hmm::MeanRemover> pipeline;
  hmm::Ptr                        input;
  input.allocate<float>(size);
  float *src = input.cast<float>();

  // create sinewave for tests later
  float *     sinewave = new float[size];
  const float dt       = 1.f / float(size);
  float       t        = 0.f;
  for (std::size_t i = 0; i < size; ++i)
  {
    sinewave[i] = sinf(2.f * M_PI * t);
    t += dt;
  }

  // Test 1: constant
  for (std::size_t i = 0; i < size; ++i)
    src[i] = 20.f;
  input = pipeline.execute(input);
  for (std::size_t i = 0; i < size; ++i)
    ASSERT_FLOAT_EQ(src[i], 0.f);

  // Test 2: sinewave
  for (std::size_t i = 0; i < size; ++i)
    src[i] = sinewave[i];
  input = pipeline.execute(input);
  for (std::size_t i = 0; i < size; ++i)
    ASSERT_NEAR(src[i], sinewave[i], 1e-4);

  // Test3 : sinewave + constant
  for (std::size_t i = 0; i < size; ++i)
    src[i] = 20.f + sinewave[i];
  input = pipeline.execute(input);
  for (std::size_t i = 0; i < size; ++i)
    ASSERT_NEAR(src[i], sinewave[i], 1e-4);

  // clean up
  input.deallocate();
  delete[] sinewave;
}

TEST(hmmTest, ZeroPadding)
{
  // initialize
  const std::size_t               size = 1000;
  hmm::Pipeline<hmm::ZeroPadding> pipeline;
  hmm::Ptr                        input;
  input.allocate<float>(size);
  float *src = input.cast<float>();

  const float constant = 20.f;
  for (std::size_t i = 0; i < size; ++i)
    src[i] = constant;
  input = pipeline.execute(input);

  src = input.cast<float>();

  ASSERT_EQ(input.size, 1024);
  ASSERT_EQ(input.capacity, 1024 * sizeof(float));
  // first n elements unchanged
  for (std::size_t i = 0; i < size; ++i)
    ASSERT_FLOAT_EQ(src[i], constant);
  // remaining elements are 0
  for (std::size_t i = size; i < input.size; ++i)
    ASSERT_FLOAT_EQ(src[i], 0.f);

  // zero padding without allocation (p.data is the same pointer)
  for (std::size_t i = 0; i < input.size; ++i)
    src[i] = constant;

  input.size = 1000;
  input      = pipeline.execute(input);
  ASSERT_EQ(input.size, 1024);
  ASSERT_EQ(input.capacity, 1024 * sizeof(float));
  ASSERT_TRUE(input.cast<float>() == src);
  for (std::size_t i = 0; i < size; ++i)
    ASSERT_FLOAT_EQ(src[i], constant);
  for (std::size_t i = size; i < input.size; ++i)
    ASSERT_FLOAT_EQ(src[i], 0.f);

  // clean up
  input.deallocate();
}

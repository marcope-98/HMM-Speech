#include <gtest/gtest.h>

#include "hmm/internal/utils.hpp"

TEST(hmmTest, round_up_to_multiple_of)
{
  EXPECT_EQ(0, hmm::utils::round_up_to_multiple_of(0, 256));
  EXPECT_EQ(256, hmm::utils::round_up_to_multiple_of(1, 256));
  EXPECT_EQ(512, hmm::utils::round_up_to_multiple_of(257, 256));
  EXPECT_EQ(768, hmm::utils::round_up_to_multiple_of(513, 256));
  EXPECT_EQ(1280, hmm::utils::round_up_to_multiple_of(1025, 256));
  EXPECT_EQ(2304, hmm::utils::round_up_to_multiple_of(2049, 256));
  EXPECT_EQ(4352, hmm::utils::round_up_to_multiple_of(4097, 256));
}

TEST(hmmTest, clp2)
{
  for (std::size_t i = 0; i < 31; ++i)
    EXPECT_EQ((1u << (i + 1)), hmm::utils::clp2((1u << i) + 1));
}

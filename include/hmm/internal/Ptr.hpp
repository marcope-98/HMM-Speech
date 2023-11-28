#ifndef HMM_INTERNAL_PTR_HPP_
#define HMM_INTERNAL_PTR_HPP_

#include <cstdint>

namespace hmm
{
  struct Ptr
  {
    void *      data{nullptr};
    std::size_t size{0};     // this is the number of elements
    std::size_t capacity{0}; // this is in bytes
  };
} // namespace hmm

#endif // HMM_INTERNAL_PTR_HPP_
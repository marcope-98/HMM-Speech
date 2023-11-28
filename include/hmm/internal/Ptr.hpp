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

    Ptr()                     = default;
    Ptr(const Ptr &other)     = default;
    Ptr(Ptr &&other) noexcept = default;
    ~Ptr()                    = default;
    Ptr &operator=(const Ptr &other) = default;
    Ptr &operator=(Ptr &&other) noexcept = default;

    template<typename T>
    void allocate(const std::size_t &size)
    {
      this->size     = size;
      this->capacity = size * sizeof(T);
      this->data     = malloc(this->capacity);
    }

    void deallocate()
    {
      free(this->data);
      this->size     = 0;
      this->capacity = 0;
    }

    template<typename T>
    void reallocate(const std::size_t &size)
    {
      free(this->data);
      this->size     = size;
      this->capacity = size * sizeof(T);
      this->data     = malloc(this->capacity);
    }

    // clang-format off
    template<typename T>
    T *cast() { return (T *)this->data; }
    // clang-format on

    template<typename T>
    void copy(const T *const src, const std::size_t &size)
    {
      if (size * sizeof(T) > this->capacity) this->reallocate<T>(size);
      if (size > this->size) this->size = size;
      T *dst = this->cast<T>();
      for (std::size_t i = 0; i < size; ++i)
        dst[i] = src[i];
    }
  };
} // namespace hmm

#endif // HMM_INTERNAL_PTR_HPP_
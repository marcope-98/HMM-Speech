#ifndef HMM_INTERNAL_PIPELINE_HPP_
#define HMM_INTERNAL_PIPELINE_HPP_

#include <utility>

#include "hmm/internal/Ptr.hpp"

namespace hmm
{
  template<typename... T>
  struct Pipeline : T...
  {
    Pipeline() = default;
    explicit Pipeline(T... t) : T(std::move(t))... {}

    auto execute(Ptr p)
    {
      ((p = static_cast<T &>(*this).execute(std::move(p))), ...);
      return p;
    }
  };

} // namespace hmm

#endif // HMM_INTERNAL_PIPELINE_HPP_
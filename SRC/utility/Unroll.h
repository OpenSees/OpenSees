//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
#pragma once
#include <utility>
#include <cstddef>

namespace OpenSees {
namespace
{
  template<std::size_t N>
  struct num { static constexpr std::size_t value = N; };

  // Helper that unpacks the index_sequence into calls to func
  template <class F, std::size_t... Is>
  [[gnu::always_inline]] inline constexpr void 
  repeat_impl(F func, std::index_sequence<Is...>) noexcept
  {
    (func(num<Is>{}), ...);
  }

    // Helper that unpacks the index_sequence into calls to func, offset by Start
    template <std::size_t Start, class F, std::size_t... Is>
    [[gnu::always_inline]] inline constexpr void 
    unroll_impl(F func, std::index_sequence<Is...>) noexcept
    {
        (func(num<Start + Is>{}), ...);
    }
}

#if 0 //  __cplusplus >= 202002L
template <std::size_t N, class F>
[[gnu::always_inline]] inline constexpr void
Repeat(F func) noexcept
{
  // Lambda with templated parameter pack (C++20 feature)
  ([]<std::size_t... Is>(F func, std::index_sequence<Is...>){
      (func(num<Is>{}), ...);
  })(func, std::make_index_sequence<N>{});
}
#else
template <std::size_t N, class F>
[[gnu::always_inline]] inline constexpr void
Repeat(F func) noexcept
{
  repeat_impl(func, std::make_index_sequence<N>{});
}
#endif 

template <std::size_t Start, std::size_t Stop, class F>
inline constexpr void
Unroll(F func) noexcept
{
  static_assert(Stop >= Start, "Stop must be greater than or equal to Start");
  constexpr std::size_t N = Stop - Start;
  unroll_impl<Start>(func, std::make_index_sequence<N>{});
}

} // namespace OpenSees
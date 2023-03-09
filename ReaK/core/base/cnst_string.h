/**
 * \file cnst_string.h
 *
 * This library defines a class to build compile-time linked-lists of string literals.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
 *
 *    THIS SOFTWARE IS DISTRIBUTED UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE v3 (GPLv3).
 *
 *    This file is part of ReaK.
 *
 *    ReaK is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    ReaK is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with ReaK (as LICENSE in the root folder).
 *    If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef REAK_CORE_BASE_CNST_STRING_H_
#define REAK_CORE_BASE_CNST_STRING_H_

#include "ReaK/core/base/defs.h"

#include <array>
#include <cstdint>
#include <string>
#include <string_view>

/** Main namespace for ReaK */
namespace ReaK {

constexpr std::uint64_t fnv_prime = 1099511628211U;
constexpr std::uint64_t fnv_offset = 14695981039346656037U;

constexpr std::uint64_t fnv_1a_hash_impl(const std::string_view& text,
                                         unsigned int i) {
  return (i == 0 ? 0
                 : (i == 1 ? ((fnv_offset ^ text[0]) * fnv_prime)
                           : ((fnv_1a_hash_impl(text, i - 1) ^ text[i - 1]) *
                              fnv_prime)));
}

constexpr std::size_t fnv_1a_hash(const std::string_view& text) {
  return fnv_1a_hash_impl(text, text.size());
}

namespace detail {
namespace {

template <std::string_view const&... Strs>
struct concat_constexpr_strings {
  // Join all strings into a single std::array of chars
  static constexpr auto impl() noexcept {
    constexpr std::size_t len = (Strs.size() + ... + 0);
    std::array<char, len> arr{};
    auto append = [i = 0, &arr](auto const& s) mutable {
      for (auto c : s) {
        arr[i++] = c;
      }
    };
    (append(Strs), ...);
    return arr;
  }
  // Give the joined string static storage
  static constexpr auto arr = impl();
  static constexpr std::string_view value{arr.data(), arr.size()};
};

}  // namespace
}  // namespace detail

template <std::string_view const&... Strs>
static constexpr auto ct_concat_v =
    detail::concat_constexpr_strings<Strs...>::value;

namespace detail {
namespace {

template <unsigned int N>
struct ct_itoa_impl {
  static constexpr auto impl() noexcept {
    constexpr unsigned int len = []() {
      unsigned int len = 0;
      for (auto n = N; n != 0; len++, n /= 10) {}
      return len;
    }();
    std::array<char, len> arr{};
    auto ptr = arr.data() + arr.size();
    for (auto n = N; n != 0; n /= 10) {
      *--ptr = "0123456789"[n % 10];
    }
    return arr;
  }

  // Give the joined string static storage
  static constexpr auto arr = impl();
  static constexpr std::string_view value{arr.data(), arr.size()};
};

}  // namespace
}  // namespace detail

template <unsigned int N>
static constexpr auto ct_itoa_v = detail::ct_itoa_impl<N>().value;

}  // namespace ReaK

#endif  // REAK_CORE_BASE_CNST_STRING_H_

/**
 * \file endian_conversions.h
 *
 * This library provides some utility functions for doing endianness conversions.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date June 2014
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

#ifndef REAK_CORE_BASE_ENDIAN_CONVERSIONS_H_
#define REAK_CORE_BASE_ENDIAN_CONVERSIONS_H_

#include <cstdint>
#include <bit>
#include <cstring>
#include <type_traits>

namespace ReaK {

// NOLINTBEGIN
#if defined(_MSC_VER)
#define RK_BSWAP_U16(X) _byteswap_ushort(X)
#define RK_BSWAP_U32(X) _byteswap_ulong(X)
#define RK_BSWAP_U64(X) _byteswap_uint64(X)
#define RK_BSWAP_I16(X) std::bit_cast<const std::int16_t>(_byteswap_ushort(std::bit_cast<const std::uint16_t>(X)))
#define RK_BSWAP_I32(X) std::bit_cast<const std::int32_t>(_byteswap_ulong(std::bit_cast<const std::uint32_t>(X)))
#define RK_BSWAP_I64(X) std::bit_cast<const std::int64_t>(_byteswap_uint64(std::bit_cast<const std::uint64_t>(X)))
#else
#define RK_BSWAP_U16(X) __builtin_bswap16(X)
#define RK_BSWAP_U32(X) __builtin_bswap32(X)
#define RK_BSWAP_U64(X) __builtin_bswap64(X)
#define RK_BSWAP_I16(X) __builtin_bswap16(X)
#define RK_BSWAP_I32(X) __builtin_bswap32(X)
#define RK_BSWAP_I64(X) __builtin_bswap64(X)
#endif
// NOLINTEND

template <std::endian En, typename T>
void to_endian(T& v) {
  if constexpr (std::endian::native == En || sizeof(T) == 1) {
    return;
  } else {
    if constexpr (std::is_integral_v<T>) {
      if constexpr (sizeof(T) == 2) {
        if constexpr (std::is_signed_v<T>) {
          v = RK_BSWAP_I16(v);
        } else {
          v = RK_BSWAP_U16(v);
        }
      } else if constexpr (sizeof(T) == 4) {
        if constexpr (std::is_signed_v<T>) {
          v = RK_BSWAP_I32(v);
        } else {
          v = RK_BSWAP_U32(v);
        }
      } else if constexpr (sizeof(T) == 8) {
        if constexpr (std::is_signed_v<T>) {
          v = RK_BSWAP_I64(v);
        } else {
          v = RK_BSWAP_U64(v);
        }
      } else {
        static_assert(sizeof(T) == 2 || sizeof(T) == 4 || sizeof(T) == 8);
      }
    } else {
      using U = std::conditional_t<sizeof(T) == 2, std::uint16_t, std::conditional_t<sizeof(T) == 4, std::uint32_t, std::uint64_t>>;
      U u{};
      static_assert(sizeof(u) == sizeof(T));
      std::memcpy(&u, &v, sizeof(u));
      if constexpr (sizeof(T) == 2) {
        u = RK_BSWAP_U16(u);
      } else if constexpr (sizeof(T) == 4) {
        u = RK_BSWAP_U32(u);
      } else if constexpr (sizeof(T) == 8) {
        u = RK_BSWAP_U64(u);
      } else {
        static_assert(sizeof(T) == 2 || sizeof(T) == 4 || sizeof(T) == 8);
      }
      std::memcpy(&v, &u, sizeof(u));
    }
  }
}

template <std::endian En, typename T>
void from_endian(T& v) {
  // Conversion from or to is the same code.
  to_endian<En, T>(v);
}

// Overloaded versions for different primitive types:

template <typename T>
void ntoh_any(T& v) {
  static_assert(std::is_integral_v<T> || std::is_floating_point_v<T>);
  to_endian<std::endian::big>(v);
}
template <typename T>
void hton_any(T& v) {
  static_assert(std::is_integral_v<T> || std::is_floating_point_v<T>);
  from_endian<std::endian::big>(v);
}

}  // namespace ReaK

#endif  // REAK_CORE_BASE_ENDIAN_CONVERSIONS_H_

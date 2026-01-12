/**
 * \file typed_primitives.h
 *
 * This library associates type information to primitive "built-in" types of C++.
 * This allows built-in types to be integrated to the ReaK::rtti system.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date april 2011
 */

/*
 *    Copyright 2011 Sven Mikael Persson
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

#ifndef REAK_CORE_RTTI_TYPED_PRIMITIVES_H_
#define REAK_CORE_RTTI_TYPED_PRIMITIVES_H_

#include <cstdint>
#include <memory>

#include "ReaK/core/rtti/so_type.h"

namespace ReaK::rtti {

template <>
struct get_type_id<std::int32_t> {
  static constexpr std::uint32_t id = 0x00000001;
  static constexpr auto type_name = std::string_view{"int"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = std::int32_t;
  using load_type = std::int32_t&;
};

template <>
struct get_type_id<std::int64_t> {
  static constexpr std::uint32_t id = 0x00000001;
  static constexpr auto type_name = std::string_view{"int"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = std::int64_t;
  using load_type = std::int64_t&;
};

template <>
struct get_type_id<std::uint32_t> {
  static constexpr std::uint32_t id = 0x00000002;
  static constexpr auto type_name = std::string_view{"unsigned int"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = std::uint32_t;
  using load_type = std::uint32_t&;
};

template <>
struct get_type_id<std::uint64_t> {
  static constexpr std::uint32_t id = 0x00000002;
  static constexpr auto type_name = std::string_view{"unsigned int"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = std::uint64_t;
  using load_type = std::uint64_t&;
};

template <>
struct get_type_id<char> {
  static constexpr std::uint32_t id = 0x00000031;
  static constexpr auto type_name = std::string_view{"char"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = char;
  using load_type = char&;
};

template <>
struct get_type_id<unsigned char> {
  static constexpr std::uint32_t id = 0x00000032;
  static constexpr auto type_name = std::string_view{"unsigned char"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = unsigned char;
  using load_type = unsigned char&;
};

template <>
struct get_type_id<float> {
  static constexpr std::uint32_t id = 0x00000003;
  static constexpr auto type_name = std::string_view{"float"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = float;
  using load_type = float&;
};

template <>
struct get_type_id<double> {
  static constexpr std::uint32_t id = 0x00000004;
  static constexpr auto type_name = std::string_view{"double"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = double;
  using load_type = double&;
};

template <>
struct get_type_id<bool> {
  static constexpr std::uint32_t id = 0x00000005;
  static constexpr auto type_name = std::string_view{"bool"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = bool;
  using load_type = bool&;
};

template <>
struct get_type_id<std::string> {
  static constexpr std::uint32_t id = 0x00000006;
  static constexpr auto type_name = std::string_view{"string"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::string&;
  using load_type = std::string&;
};

template <typename Tail>
struct get_type_info<std::string, Tail> {
  using type = so_type_details::type_id<std::string, typename Tail::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::string>::type_name,
                  get_type_name_tail<Tail>::value>;
};

template <typename T>
struct get_type_id<std::shared_ptr<T>> {
  static constexpr std::uint32_t id = get_type_id<T>::id;
  static constexpr auto type_name = std::string_view{"shared_ptr"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::shared_ptr<T>&;
  using load_type = std::shared_ptr<T>&;
};

template <typename T>
struct get_type_id<std::weak_ptr<T>> {
  static constexpr std::uint32_t id = get_type_id<T>::id;
  static constexpr auto type_name = std::string_view{"weak_ptr"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::weak_ptr<T>&;
  using load_type = std::weak_ptr<T>&;
};

template <typename T>
struct get_type_id<std::unique_ptr<T>> {
  static constexpr std::uint32_t id = get_type_id<T>::id;
  static constexpr auto type_name = std::string_view{"unique_ptr"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::unique_ptr<T>&;
  using load_type = std::unique_ptr<T>&;
};

}  // namespace ReaK::rtti

#endif  // REAK_CORE_RTTI_TYPED_PRIMITIVES_H_

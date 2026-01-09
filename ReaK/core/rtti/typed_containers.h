/**
 * \file typed_containers.h
 *
 * This library associates type information to STL containers types of C++ standard libraries.
 * This allows STL containers to be integrated to the ReaK::rtti system.
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

#ifndef REAK_CORE_RTTI_TYPED_CONTAINERS_H_
#define REAK_CORE_RTTI_TYPED_CONTAINERS_H_

#include "ReaK/core/rtti/so_type.h"

#include <vector>

namespace ReaK::rtti {

template <typename T>
struct get_type_id<std::vector<T>> {
  static constexpr unsigned int id = 0x00000008;
  static constexpr auto type_name = std::string_view{"std::vector"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::vector<T>&;
  using load_type = std::vector<T>&;
};

template <typename T, typename Tail>
struct get_type_info<std::vector<T>, Tail> {
  using type = so_type_details::type_id<std::vector<T>,
                                        typename get_type_info<T, Tail>::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::vector<T>>::type_name, lsl_left_bracket,
                  get_type_id<T>::type_name, lsl_right_bracket,
                  get_type_name_tail<Tail>::value>;
};
}  // namespace ReaK::rtti

#include <list>

namespace ReaK::rtti {

template <typename T>
struct get_type_id<std::list<T>> {
  static constexpr unsigned int id = 0x00000009;
  static constexpr auto type_name = std::string_view{"std::list"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::list<T>&;
  using load_type = std::list<T>&;
};

template <typename T, typename Tail>
struct get_type_info<std::list<T>, Tail> {
  using type = so_type_details::type_id<std::list<T>,
                                        typename get_type_info<T, Tail>::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::list<T>>::type_name, lsl_left_bracket,
                  get_type_id<T>::type_name, lsl_right_bracket,
                  get_type_name_tail<Tail>::value>;
};
}  // namespace ReaK::rtti

#include <map>

namespace ReaK::rtti {

template <typename Key, typename T>
struct get_type_id<std::map<Key, T>> {
  static constexpr unsigned int id = 0x0000000A;
  static constexpr auto type_name = std::string_view{"std::map"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::map<Key, T>&;
  using load_type = std::map<Key, T>&;
};

template <typename Key, typename T, typename Tail>
struct get_type_info<std::map<Key, T>, Tail> {
  using type = so_type_details::type_id<
      std::map<Key, T>,
      typename get_type_info<Key, get_type_info<T, Tail>>::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::map<Key, T>>::type_name, lsl_left_bracket,
                  get_type_id<Key>::type_name, lsl_comma,
                  get_type_id<T>::type_name, lsl_right_bracket,
                  get_type_name_tail<Tail>::value>;
};

template <typename Key, typename T>
struct get_type_id<std::multimap<Key, T>> {
  static constexpr unsigned int id = 0x0000000D;
  static constexpr auto type_name = std::string_view{"std::multimap"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::multimap<Key, T>&;
  using load_type = std::multimap<Key, T>&;
};

template <typename Key, typename T, typename Tail>
struct get_type_info<std::multimap<Key, T>, Tail> {
  using type = so_type_details::type_id<
      std::multimap<Key, T>,
      typename get_type_info<Key, get_type_info<T, Tail>>::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::multimap<Key, T>>::type_name,
                  lsl_left_bracket, get_type_id<Key>::type_name, lsl_comma,
                  get_type_id<T>::type_name, lsl_right_bracket,
                  get_type_name_tail<Tail>::value>;
};
}  // namespace ReaK::rtti

#include <set>

namespace ReaK::rtti {

template <typename T>
struct get_type_id<std::set<T>> {
  static constexpr unsigned int id = 0x0000000B;
  static constexpr auto type_name = std::string_view{"std::set"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::set<T>&;
  using load_type = std::set<T>&;
};

template <typename T, typename Tail>
struct get_type_info<std::set<T>, Tail> {
  using type = so_type_details::type_id<std::set<T>,
                                        typename get_type_info<T, Tail>::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::set<T>>::type_name, lsl_left_bracket,
                  get_type_id<T>::type_name, lsl_right_bracket,
                  get_type_name_tail<Tail>::value>;
};

template <typename T>
struct get_type_id<std::multiset<T>> {
  static constexpr unsigned int id = 0x0000000E;
  static constexpr auto type_name = std::string_view{"std::multiset"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::multiset<T>&;
  using load_type = std::multiset<T>&;
};

template <typename T, typename Tail>
struct get_type_info<std::multiset<T>, Tail> {
  using type = so_type_details::type_id<std::multiset<T>,
                                        typename get_type_info<T, Tail>::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::multiset<T>>::type_name, lsl_left_bracket,
                  get_type_id<T>::type_name, lsl_right_bracket,
                  get_type_name_tail<Tail>::value>;
};
}  // namespace ReaK::rtti

#include <utility>

namespace ReaK::rtti {

template <typename T1, typename T2>
struct get_type_id<std::pair<T1, T2>> {
  static constexpr unsigned int id = 0x0000000C;
  static constexpr auto type_name = std::string_view{"std::pair"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::pair<T1, T2>&;
  using load_type = std::pair<T1, T2>&;
};

template <typename T1, typename T2, typename Tail>
struct get_type_info<std::pair<T1, T2>, Tail> {
  using type = so_type_details::type_id<
      std::pair<T1, T2>,
      typename get_type_info<T1, get_type_info<T2, Tail>>::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::pair<T1, T2>>::type_name, lsl_left_bracket,
                  get_type_id<T1>::type_name, lsl_comma,
                  get_type_id<T2>::type_name, lsl_right_bracket,
                  get_type_name_tail<Tail>::value>;
};
}  // namespace ReaK::rtti

#include <forward_list>

namespace ReaK::rtti {

template <typename T>
struct get_type_id<std::forward_list<T>> {
  static constexpr unsigned int id = 0x00000040;
  static constexpr auto type_name = std::string_view{"std::forward_list"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::forward_list<T>&;
  using load_type = std::forward_list<T>&;
};

template <typename T, typename Tail>
struct get_type_info<std::forward_list<T>, Tail> {
  using type = so_type_details::type_id<std::forward_list<T>,
                                        typename get_type_info<T, Tail>::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::forward_list<T>>::type_name,
                  lsl_left_bracket, get_type_id<T>::type_name,
                  lsl_right_bracket, get_type_name_tail<Tail>::value>;
};
}  // namespace ReaK::rtti

#include <array>

namespace ReaK::rtti {

template <typename T, std::size_t N>
struct get_type_id<std::array<T, N>> {
  static constexpr unsigned int id = 0x00000041;
  static constexpr auto type_name = std::string_view{"std::array"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::array<T, N>&;
  using load_type = std::array<T, N>&;
};

template <typename T, std::size_t N, typename Tail>
struct get_type_info<std::array<T, N>, Tail> {
  using type = so_type_details::type_id<
      std::array<T, N>,
      typename get_type_info<
          T,
          get_type_info<std::integral_constant<unsigned int, N>, Tail>>::type>;
  static constexpr auto type_name = ct_concat_v<
      get_type_id<std::array<T, N>>::type_name, lsl_left_bracket,
      get_type_id<T>::type_name, lsl_comma,
      get_type_id<std::integral_constant<unsigned int, N>>::type_name,
      lsl_right_bracket, get_type_name_tail<Tail>::value>;
};
}  // namespace ReaK::rtti

#include <tuple>

namespace ReaK::rtti {

template <typename... T>
struct get_type_id<std::tuple<T...>> {
  static constexpr unsigned int id = 0x00000042;
  static constexpr auto type_name = std::string_view{"std::tuple"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::tuple<T...>&;
  using load_type = std::tuple<T...>&;
};

template <typename Tail, typename... T>
struct get_type_info<std::tuple<T...>, Tail> {
  using type =
      so_type_details::type_id<std::tuple<T...>,
                               typename so_type_details::get_type_info_seq<
                                   T...>::template with_tail<Tail>::type::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::tuple<T...>>::type_name, lsl_left_bracket,
                  so_type_details::get_type_info_seq<T...>::type_name,
                  lsl_right_bracket, get_type_name_tail<Tail>::value>;
};
}  // namespace ReaK::rtti

#include <unordered_map>

namespace ReaK::rtti {

template <typename Key, typename T>
struct get_type_id<std::unordered_map<Key, T>> {
  static constexpr unsigned int id = 0x00000044;
  static constexpr auto type_name = std::string_view{"std::unordered_map"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::unordered_map<Key, T>&;
  using load_type = std::unordered_map<Key, T>&;
};

template <typename Key, typename T, typename Tail>
struct get_type_info<std::unordered_map<Key, T>, Tail> {
  using type = so_type_details::type_id<
      std::unordered_map<Key, T>,
      typename get_type_info<Key, get_type_info<T, Tail>>::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::unordered_map<Key, T>>::type_name,
                  lsl_left_bracket, get_type_id<Key>::type_name, lsl_comma,
                  get_type_id<T>::type_name, lsl_right_bracket,
                  get_type_name_tail<Tail>::value>;
};

template <typename Key, typename T>
struct get_type_id<std::unordered_multimap<Key, T>> {
  static constexpr unsigned int id = 0x00000046;
  static constexpr auto type_name = std::string_view{"std::unordered_multimap"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::unordered_multimap<Key, T>&;
  using load_type = std::unordered_multimap<Key, T>&;
};

template <typename Key, typename T, typename Tail>
struct get_type_info<std::unordered_multimap<Key, T>, Tail> {
  using type = so_type_details::type_id<
      std::unordered_multimap<Key, T>,
      typename get_type_info<Key, get_type_info<T, Tail>>::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::unordered_multimap<Key, T>>::type_name,
                  lsl_left_bracket, get_type_id<Key>::type_name, lsl_comma,
                  get_type_id<T>::type_name, lsl_right_bracket,
                  get_type_name_tail<Tail>::value>;
};
}  // namespace ReaK::rtti

#include <unordered_set>

namespace ReaK::rtti {

template <typename T>
struct get_type_id<std::unordered_set<T>> {
  static constexpr unsigned int id = 0x00000045;
  static constexpr auto type_name = std::string_view{"std::unordered_set"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::unordered_set<T>&;
  using load_type = std::unordered_set<T>&;
};

template <typename T, typename Tail>
struct get_type_info<std::unordered_set<T>, Tail> {
  using type = so_type_details::type_id<std::unordered_set<T>,
                                        typename get_type_info<T, Tail>::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::unordered_set<T>>::type_name,
                  lsl_left_bracket, get_type_id<T>::type_name,
                  lsl_right_bracket, get_type_name_tail<Tail>::value>;
};

template <typename T>
struct get_type_id<std::unordered_multiset<T>> {
  static constexpr unsigned int id = 0x00000047;
  static constexpr auto type_name = std::string_view{"std::unordered_multiset"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const std::unordered_multiset<T>&;
  using load_type = std::unordered_multiset<T>&;
};

template <typename T, typename Tail>
struct get_type_info<std::unordered_multiset<T>, Tail> {
  using type = so_type_details::type_id<std::unordered_multiset<T>,
                                        typename get_type_info<T, Tail>::type>;
  static constexpr auto type_name =
      ct_concat_v<get_type_id<std::unordered_multiset<T>>::type_name,
                  lsl_left_bracket, get_type_id<T>::type_name,
                  lsl_right_bracket, get_type_name_tail<Tail>::value>;
};

}  // namespace ReaK::rtti

#endif  // REAK_CORE_RTTI_TYPED_CONTAINERS_H_

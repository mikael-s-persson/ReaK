/*
 * \file rpc_function_helpers.hpp
 *
 * This library defines some underlying details for the implementation of the 
 * remote procedure call (rpc) function class template from ReaK. None of the elements within this header 
 * are meant to be used in user-side code, only library code.
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

#ifndef REAK_DETAIL_RPC_FUNCTION_HELPERS_HPP
#define REAK_DETAIL_RPC_FUNCTION_HELPERS_HPP

#include <functional>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include "ReaK/core/serialization/archiver.hpp"

namespace ReaK::rpc::detail {

template <int Dummy>
std::string concat_arg_types() {
  return "";
}

template <typename T>
std::string arg_type_to_string() {
  if constexpr (std::is_volatile_v<T>) {
    return "volatile" + arg_type_to_string<std::remove_volatile_t<T>>();
  } else if constexpr (!std::is_reference_v<T> && std::is_const_v<T>) {
    return "const" + arg_type_to_string<std::remove_const_t<T>>();
  } else if constexpr (std::is_lvalue_reference_v<T>) {
    return arg_type_to_string<std::remove_reference_t<T>>() + "&";
  } else if constexpr (std::is_rvalue_reference_v<T>) {
    return arg_type_to_string<std::remove_reference_t<T>>() + "&&";
  } else {
    return std::string(rtti::get_type_info<T>::type_name);
  }
}

template <int Dummy, typename Arg1, typename... Args>
std::string concat_arg_types() {
  return arg_type_to_string<Arg1>() + ", " + concat_arg_types<Dummy, Args...>();
}

template <int ByteCount>
struct fnv_1a_numbers {
  static const std::uint32_t offset = 2166136261;
  static const std::uint32_t prime = 16777619;
};

template <>
struct fnv_1a_numbers<8> {
  static const std::uint64_t offset = 14695981039346656037UL;
  static const std::uint64_t prime = 1099511628211UL;
};

std::uint32_t get_fnv_1a_hash(const std::string& s) {
  std::uint32_t result = fnv_1a_numbers<4>::offset;
  union {
    char c;
    unsigned char uc;
  } c_to_uc;
  for (char i : s) {
    c_to_uc.c = i;
    result ^= c_to_uc.uc;
    result *= fnv_1a_numbers<4>::prime;
  }
  return result;
}

template <typename... Args>
struct get_func_input_tuple {
  using type = std::tuple<typename std::remove_cv<
      typename std::remove_reference<Args>::type>::type...>;
};

template <int Idx>
struct tuple_saver_impl;

template <>
struct tuple_saver_impl<0> {
  template <typename Tuple>
  static void apply(serialization::oarchive& /*unused*/,
                    const Tuple& /*unused*/){};
};

template <int Idx>
struct tuple_saver_impl {
  template <typename Tuple>
  static void apply(serialization::oarchive& lhs, const Tuple& rhs) {
    tuple_saver_impl<Idx - 1>::apply(lhs, rhs);
    lhs << std::get<Idx - 1>(rhs);
  }
};

template <int Idx>
struct tuple_loader_impl;

template <>
struct tuple_loader_impl<0> {
  template <typename Tuple>
  static void apply(serialization::iarchive& /*unused*/, Tuple& /*unused*/) {}
};

template <int Idx>
struct tuple_loader_impl {
  template <typename Tuple>
  static void apply(serialization::iarchive& lhs, Tuple& rhs) {
    tuple_loader_impl<Idx - 1>::apply(lhs, rhs);
    lhs >> std::get<Idx - 1>(rhs);
  }
};

template <typename T>
struct generic_return_type {
  T value;
  template <typename Func, typename... Args>
  void apply(const Func& f, Args&... args) {
    value = f(args...);
  }
  void save_to(serialization::oarchive& out) { out << value; }
  void load_from(serialization::iarchive& in) { in >> value; }
  T take_value() { return std::move(value); };
};

template <>
struct generic_return_type<void> {
  template <typename Func, typename... Args>
  void apply(const Func& f, Args&... args) {
    f(args...);
  }
  void save_to(serialization::oarchive& /*unused*/) {}
  void load_from(serialization::iarchive& /*unused*/) {}
  void take_value() {}
};

template <int Idx>
struct input_tuple_unroller;

template <>
struct input_tuple_unroller<0> {
  template <typename Func, typename Tuple, typename... Args>
  static generic_return_type<typename Func::result_type> apply(
      const Func& f, Tuple& /*unused*/, Args&... args) {
    generic_return_type<typename Func::result_type> result;
    result.apply(f, args...);
    return result;
  }
};

template <int Idx>
struct input_tuple_unroller {
  template <typename Func, typename Tuple, typename... Args>
  static generic_return_type<typename Func::result_type> apply(const Func& f,
                                                               Tuple& t,
                                                               Args&... args) {
    return input_tuple_unroller<Idx - 1>::apply(f, t, std::get<Idx - 1>(t),
                                                args...);
  }
};

template <int Idx, typename Tuple>
void save_output_from_input(serialization::oarchive& out, Tuple& inputs) {
  /*Finished all the args.*/
}

template <int Idx, typename Tuple, typename Arg1, typename... Args>
void save_output_from_input(serialization::oarchive& out, Tuple& inputs) {
  if constexpr (std::is_lvalue_reference_v<Arg1> &&
                !std::is_const_v<std::remove_reference_t<Arg1>>) {
    out << std::get<Idx>(inputs);
  }
  save_output_from_input<Idx + 1, Tuple, Args...>(out, inputs);
}

template <typename... LArgs>
struct input_saver {
  static void apply(serialization::oarchive& out){};
};

template <typename LArg1, typename... LArgs>
struct input_saver<LArg1, LArgs...> {
  template <typename RArg1, typename... RArgs>
  static void apply(serialization::oarchive& out, RArg1&& a1, RArgs&&... args) {
    out << LArg1(std::forward<RArg1>(a1));
    input_saver<LArgs...>::apply(out, std::forward<RArgs>(args)...);
  }
};

template <typename... LArgs>
struct output_loader {
  static void apply(serialization::iarchive& in){};
};

template <typename LArg1, typename... LArgs>
struct output_loader<const LArg1&, LArgs...> {
  template <typename RArg1, typename... RArgs>
  static void apply(serialization::iarchive& in, RArg1&& a1, RArgs&&... args) {
    output_loader<LArgs...>::apply(in, std::forward<RArgs>(args)...);
  }
};

template <typename LArg1, typename... LArgs>
struct output_loader<LArg1&, LArgs...> {
  template <typename RArg1, typename... RArgs>
  static void apply(serialization::iarchive& in, RArg1&& a1, RArgs&&... args) {
    LArg1 tmp_a1;
    in >> tmp_a1;
    a1 = std::move(tmp_a1);
    output_loader<LArgs...>::apply(in, std::forward<RArgs>(args)...);
  }
};

template <typename LArg1, typename... LArgs>
struct output_loader<LArg1, LArgs...> {
  template <typename RArg1, typename... RArgs>
  static void apply(serialization::iarchive& in, RArg1&& a1, RArgs&&... args) {
    output_loader<LArgs...>::apply(in, std::forward<RArgs>(args)...);
  }
};

}  // namespace ReaK::rpc::detail

#endif

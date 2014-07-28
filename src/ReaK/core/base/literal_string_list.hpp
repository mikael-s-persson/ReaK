/**
 * \file literal_string_list.hpp
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

#ifndef REAK_LITERAL_STRING_LIST_HPP
#define REAK_LITERAL_STRING_LIST_HPP

#include "defs.hpp"

#include <cstddef>
#include <cstring>

#include <array>
#include <tuple>
#include <string>

/** Main namespace for ReaK */
namespace ReaK {


BOOST_CONSTEXPR_OR_CONST std::size_t fnv_prime  = (sizeof(std::size_t) == 8 ? 1099511628211u : 16777619u);
BOOST_CONSTEXPR_OR_CONST std::size_t fnv_offset = (sizeof(std::size_t) == 8 ? 14695981039346656037u : 2166136261u);

BOOST_CONSTEXPR std::size_t fnv_1a_hash_impl(const char* text_ptr, unsigned int i) {
  return (i == 0 ? 
    ((fnv_offset ^ text_ptr[0]) * fnv_prime) : 
    ((fnv_1a_hash_impl(text_ptr, i-1) ^ text_ptr[i]) * fnv_prime));
};

template <unsigned int N>
BOOST_CONSTEXPR std::size_t fnv_1a_hash(const char (&text_ptr)[N]) {
  return fnv_1a_hash_impl(text_ptr, N-2);
};


template <unsigned int... Ns>
struct literal_string_array {
  
  BOOST_STATIC_CONSTEXPR unsigned int count = sizeof...(Ns);
  
  std::tuple< std::array<char, Ns>... > data;
  
  BOOST_CONSTEXPR literal_string_array(const std::array<char, Ns>&... aStr) : data(aStr...) { };
  
  
  template <unsigned int Id>
  BOOST_CONSTEXPR typename std::enable_if< (Id == 0), 
  std::size_t >::type size_impl() const {
    return std::get<0>(data).size() - 1;
  };
  template <unsigned int Id>
  BOOST_CONSTEXPR typename std::enable_if< (Id > 0), 
  std::size_t >::type size_impl() const {
    return std::get<Id>(data).size() - 1 + size_impl<Id-1>();
  };
  BOOST_CONSTEXPR std::size_t size() const {
    return size_impl<count-1>();
  };
  
  
  template <unsigned int Id>
  BOOST_CONSTEXPR typename std::enable_if< (Id == count), 
  char >::type ind_operator_impl(unsigned int i) const {
    return '\0';
  };
  template <unsigned int Id>
  BOOST_CONSTEXPR typename std::enable_if< (Id < count), 
  char >::type ind_operator_impl(unsigned int i) const {
    return ((i + 1 < std::get<Id>(data).size()) ? std::get<Id>(data)[i] : 
      ind_operator_impl<Id+1>(i + 1 - std::get<Id>(data).size()));
  };
  BOOST_CONSTEXPR char operator[](unsigned int i) const {
    return ind_operator_impl<0>(i);
  };
  
  
  template <unsigned int Id>
  BOOST_CONSTEXPR typename std::enable_if< (Id == 0), 
  std::size_t >::type hash_impl(unsigned int i) const {
    return (i <= 1 ? 
      ((fnv_offset ^ std::get<Id>(data)[0]) * fnv_prime) : 
      ((hash_impl<Id>(i-1) ^ std::get<Id>(data)[i-1]) * fnv_prime));
  };
  template <unsigned int Id>
  BOOST_CONSTEXPR typename std::enable_if< (Id > 0), 
  std::size_t >::type hash_impl(unsigned int i) const {
    return (i <= 1 ? 
      ((hash_impl<Id-1>(std::get<Id-1>(data).size()-1) ^ std::get<Id>(data)[0]) * fnv_prime) : 
      ((hash_impl<Id>(i-1) ^ std::get<Id>(data)[i-1]) * fnv_prime));
  };
  BOOST_CONSTEXPR std::size_t hash() const {
    return hash_impl<count-1>(std::get<count-1>(data).size()-1);
  };
  
  
  template <unsigned int Id>
  typename std::enable_if< (Id == count), 
  void >::type to_string_impl(std::string& , unsigned int ) const { };
  template <unsigned int Id>
  typename std::enable_if< (Id < count), 
  void >::type to_string_impl(std::string& result, unsigned int i) const {
    for(unsigned int j = 0; j < std::get<Id>(data).size()-1; ++j)
      result[i++] = std::get<Id>(data)[j];
    to_string_impl<Id+1>(result, i);
  };
  std::string to_string() const {
    std::string result(size_impl<count-1>(),' ');
    to_string_impl<0>(result, 0);
    return result;
  };
  
};

#define RK_LSA(LITERAL_STR) ::ReaK::literal_string_array<sizeof(LITERAL_STR)>{std::array<char,sizeof(LITERAL_STR)>{LITERAL_STR}}


namespace detail {

namespace {

template <unsigned int ...>
struct index_seq { };

template <unsigned int N, unsigned int ...S>
struct gen_index_seq : gen_index_seq<N-1, N-1, S...> { };

template <unsigned int ...S>
struct gen_index_seq<0, S...> {
  typedef index_seq<S...> type;
};

};

};


template <unsigned int... Ns, unsigned int... NIds, 
          unsigned int... Ms, unsigned int... MIds>
BOOST_CONSTEXPR literal_string_array<Ns..., Ms...> add_literal_arrays_impl(
  const literal_string_array<Ns...>& lhs, detail::index_seq<NIds...>,
  const literal_string_array<Ms...>& rhs, detail::index_seq<MIds...>) {
  return literal_string_array<Ns..., Ms...>(
    std::get<NIds>(lhs.data)..., std::get<MIds>(rhs.data)...);
};

template <unsigned int... Ns, unsigned int... Ms>
BOOST_CONSTEXPR literal_string_array<Ns..., Ms...> operator+(const literal_string_array<Ns...>& lhs, 
                                                       const literal_string_array<Ms...>& rhs) {
  return add_literal_arrays_impl(lhs, typename detail::gen_index_seq<sizeof...(Ns)>::type(), 
                                 rhs, typename detail::gen_index_seq<sizeof...(Ms)>::type());
};


namespace detail {

namespace {

BOOST_CONSTEXPR const char str_digits[] = "0123456789";

template <unsigned int...>
struct digit_seq { };

template <unsigned int N, unsigned int... S>
struct gen_digit_seq : gen_digit_seq<N / 10, N % 10, S...> { };

template <unsigned int... S>
struct gen_digit_seq<0, S...> {
  typedef digit_seq<S...,10> type;
  BOOST_STATIC_CONSTEXPR unsigned int count = sizeof...(S) + 1;
};

};

};

template <unsigned int N>
struct ct_itoa {
  template <unsigned int... S>
  static BOOST_CONSTEXPR std::array<char, detail::gen_digit_seq<N>::count > text_impl(detail::digit_seq<S...>) {
    return std::array<char, detail::gen_digit_seq<N>::count >{detail::str_digits[S]...};
  };
  BOOST_STATIC_CONSTEXPR literal_string_array< detail::gen_digit_seq<N>::count > text{text_impl(typename detail::gen_digit_seq<N>::type())};
};


};

#endif





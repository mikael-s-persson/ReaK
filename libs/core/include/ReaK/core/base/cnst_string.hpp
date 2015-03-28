/**
 * \file cnst_string.hpp
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

#ifndef REAK_CNST_STRING_HPP
#define REAK_CNST_STRING_HPP

#include "defs.hpp"

#if( !defined( BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX ) && !defined( BOOST_NO_CXX11_HDR_ARRAY ) \
     && !defined( BOOST_NO_CXX11_AUTO_DECLARATIONS ) )

#include <string>
#include <array>

#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES
#include "index_sequence.hpp"
#endif

/** Main namespace for ReaK */
namespace ReaK {


BOOST_CONSTEXPR_OR_CONST std::size_t fnv_prime = ( sizeof( std::size_t ) == 8 ? 1099511628211u : 16777619u );
BOOST_CONSTEXPR_OR_CONST std::size_t fnv_offset = ( sizeof( std::size_t ) == 8 ? 14695981039346656037u : 2166136261u );

BOOST_CONSTEXPR std::size_t fnv_1a_hash_impl( const char* text_ptr, unsigned int i ) {
  return ( i == 0 ? ( ( fnv_offset ^ text_ptr[0] ) * fnv_prime )
                  : ( ( fnv_1a_hash_impl( text_ptr, i - 1 ) ^ text_ptr[i] ) * fnv_prime ) );
};

template < unsigned int N >
BOOST_CONSTEXPR std::size_t fnv_1a_hash( const char ( &text_ptr )[N] ) {
  return fnv_1a_hash_impl( text_ptr, N - 2 );
};


template < std::size_t N >
struct cnst_string {

  std::array< char, N > data;

  BOOST_CONSTEXPR cnst_string( const std::array< char, N >& aStr ) : data( aStr ){};

  BOOST_CONSTEXPR std::size_t size() const { return N - 1; };

  BOOST_CONSTEXPR char operator[]( std::size_t i ) const { return ( i < N - 1 ? data[i] : '\0' ); };

  BOOST_CONSTEXPR std::size_t fnv_1a_hash_impl( unsigned int i ) const {
    return ( i == 0 ? ( ( fnv_offset ^ data[0] ) * fnv_prime )
                    : ( ( fnv_1a_hash_impl( i - 1 ) ^ data[i] ) * fnv_prime ) );
  };
  BOOST_CONSTEXPR std::size_t hash() const { return fnv_1a_hash_impl( N - 2 ); };

  std::string to_string() const { return std::string( data.data() ); };
};

#define RK_LSA( LITERAL_STR )                                 \
  ::ReaK::cnst_string< sizeof( LITERAL_STR ) > {              \
    std::array< char, sizeof( LITERAL_STR ) > { LITERAL_STR } \
  }

#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES

namespace detail {
namespace {

template < std::size_t N, std::size_t... NIds, std::size_t M, std::size_t... MIds >
BOOST_CONSTEXPR cnst_string< N + M - 1 >
  concat_char_arrays_impl( const std::array< char, N >& lhs, index_sequence< NIds... >,
                           const std::array< char, M >& rhs, index_sequence< MIds... > ) {
  return cnst_string< N + M - 1 >{std::array< char, N + M - 1 >{{lhs[NIds]..., rhs[MIds]..., '\0'}}};
};
};
};

template < std::size_t N, std::size_t M >
BOOST_CONSTEXPR cnst_string< N + M - 1 > operator+( const cnst_string< N >& lhs, const cnst_string< M >& rhs ) {
  typedef typename make_index_sequence< N - 1 >::type LIndices;
  typedef typename make_index_sequence< M - 1 >::type RIndices;
  return detail::concat_char_arrays_impl( lhs.data, LIndices(), rhs.data, RIndices() );
};


namespace detail {
namespace {

BOOST_CONSTEXPR const char str_digits[] = "0123456789";

template < std::size_t... >
struct digit_seq {};

template < std::size_t N, std::size_t... S >
struct gen_digit_seq : gen_digit_seq< N / 10, N % 10, S... > {};

template < std::size_t... S >
struct gen_digit_seq< 0, S... > {
  typedef digit_seq< S..., 10 > type;
  BOOST_STATIC_CONSTEXPR std::size_t count = sizeof...(S)+1;
};
};
};

template < std::size_t N >
struct ct_itoa {
  template < std::size_t... S >
  static BOOST_CONSTEXPR cnst_string< detail::gen_digit_seq< N >::count > text_impl( detail::digit_seq< S... > ) {
    return cnst_string< detail::gen_digit_seq< N >::count >{
      std::array< char, detail::gen_digit_seq< N >::count >{detail::str_digits[S]...}};
  };
  BOOST_STATIC_CONSTEXPR cnst_string< detail::gen_digit_seq< N >::count > text{
    text_impl( typename detail::gen_digit_seq< N >::type() )};
};

#endif
};

#else

#pragma message( \
  "Warning: The 'cnst_string.hpp' header (from ReaK library) was included, but this compiler does not support all the necessary C++11 features: unified initialization, std::array, auto declarations, and constexpr. cnst_string will be disabled, which could cause errors!" )

#endif

#endif

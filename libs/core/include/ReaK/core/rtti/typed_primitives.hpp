/**
 * \file typed_primitives.hpp
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

#ifndef REAK_TYPED_PRIMITIVES_HPP
#define REAK_TYPED_PRIMITIVES_HPP

#include "so_type.hpp"
#include <cstdint>

namespace ReaK {

namespace rtti {


template <>
struct get_type_id< std::int32_t > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000001 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "int" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "int"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef std::int32_t save_type;
  typedef std::int32_t& load_type;
};

template <>
struct get_type_id< std::int64_t > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000001 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "int" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "int"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef std::int64_t save_type;
  typedef std::int64_t& load_type;
};

template <>
struct get_type_id< std::uint32_t > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000002 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "unsigned int" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "unsigned int"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef std::uint32_t save_type;
  typedef std::uint32_t& load_type;
};

template <>
struct get_type_id< std::uint64_t > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000002 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "unsigned int" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "unsigned int"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef std::uint64_t save_type;
  typedef std::uint64_t& load_type;
};

template <>
struct get_type_id< char > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000031 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "char" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "char"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef char save_type;
  typedef char& load_type;
};

template <>
struct get_type_id< unsigned char > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000032 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "unsigned char" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "unsigned char"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef unsigned char save_type;
  typedef unsigned char& load_type;
};

template <>
struct get_type_id< float > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000003 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "float" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "float"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef float save_type;
  typedef float& load_type;
};

template <>
struct get_type_id< double > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000004 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "double" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "double"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef double save_type;
  typedef double& load_type;
};

template <>
struct get_type_id< bool > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000005 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "bool" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "bool"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef bool save_type;
  typedef bool& load_type;
};

template <>
struct get_type_id< std::string > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000006 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "string" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "string"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef const std::string& save_type;
  typedef std::string& load_type;
};

template < typename Tail >
struct get_type_info< std::string, Tail > {
  typedef type_id< std::string, typename Tail::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::string >::type_name + get_type_name_tail< Tail >::value;
#else
  static std::string type_name() {
    std::string result = get_type_id< std::string >::type_name();
    result += get_type_name_tail< Tail >::value();
    return result; // NRVO
  };
#endif
};


template < typename T >
struct get_type_id< shared_ptr< T > > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = get_type_id< T >::ID );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "shared_ptr" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "shared_ptr"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef const shared_ptr< T >& save_type;
  typedef shared_ptr< T >& load_type;
};

template < typename T >
struct get_type_id< weak_ptr< T > > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = get_type_id< T >::ID );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "weak_ptr" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "weak_ptr"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef const weak_ptr< T >& save_type;
  typedef weak_ptr< T >& load_type;
};


#ifndef BOOST_NO_CXX11_SMART_PTR

template < typename T >
struct get_type_id< unique_ptr< T > > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = get_type_id< T >::ID );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "unique_ptr" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "unique_ptr"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef const unique_ptr< T >& save_type;
  typedef unique_ptr< T >& load_type;
};

#endif
};
};


#endif

/**
 * \file disc_headers.hpp
 *
 * This library declares helpers to generate and parse the headers of the ReaK.DISC protocol.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2015
 */

/*
 *    Copyright 2015 Sven Mikael Persson
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


#ifndef REAK_DISC_HEADERS_HPP
#define REAK_DISC_HEADERS_HPP

#include <ReaK/core/serialization/archiver.hpp>
#include <ReaK/core/base/serializable.hpp>

#include <ReaK/core/rpc/version.hpp>

#include <string>
#include <tuple>
#include <iostream>

namespace ReaK {

namespace rpc {

namespace detail {


enum disc_header_type { disc_newport = 0, disc_giveport, disc_advertise, disc_unadvertise, disc_echo };

extern const char* disc_header_type_to_str[];

std::string generate_disc_header( disc_header_type h, msg_format fmt, std::size_t msg_size );

std::tuple< disc_header_type, msg_format, std::size_t > parse_disc_header( std::istream& hdr_in );


struct disc_basic_header : serializable {
  std::string name;
  std::string addr;

  disc_basic_header( const std::string& aName = "", const std::string& aAddr = "127.0.0.0" )
      : name( aName ), addr( aAddr ){};


  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    A& RK_SERIAL_SAVE_WITH_NAME( name ) & RK_SERIAL_SAVE_WITH_NAME( addr );
  };
  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    A& RK_SERIAL_LOAD_WITH_NAME( name ) & RK_SERIAL_LOAD_WITH_NAME( addr );
  };

  RK_RTTI_REGISTER_CLASS_1BASE( disc_basic_header, 1, serializable )
};

inline bool operator<( const disc_basic_header& lhs, const disc_basic_header& rhs ) {
  if( lhs.name < rhs.name )
    return true;
  if( lhs.name == rhs.name ) {
    if( lhs.addr < rhs.addr )
      return true;
  };
  return false;
};


struct disc_with_port_header : disc_basic_header {
  unsigned int port;

  disc_with_port_header( const std::string& aName = "", const std::string& aAddr = "127.0.0.0", unsigned int aPort = 0 )
      : disc_basic_header( aName, aAddr ), port( aPort ){};


  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    disc_basic_header::save( A, 1 );
    A& RK_SERIAL_SAVE_WITH_NAME( port );
  };
  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    disc_basic_header::load( A, 1 );
    A& RK_SERIAL_LOAD_WITH_NAME( port );
  };

  RK_RTTI_REGISTER_CLASS_1BASE( disc_with_port_header, 1, disc_basic_header )
};
};
};


namespace rtti {

#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS

#define RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS( VARIABLE, CLASSID )                                               \
  template <>                                                                                                          \
  struct get_type_id< ::ReaK::rpc::detail::VARIABLE > {                                                                \
    BOOST_STATIC_CONSTANT( unsigned int, ID = CLASSID );                                                               \
    BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( BOOST_PP_STRINGIZE( VARIABLE ) );                                  \
    static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };                                               \
                                                                                                                       \
    typedef const serializable& save_type;                                                                             \
    typedef serializable& load_type;                                                                                   \
  };                                                                                                                   \
                                                                                                                       \
  template < typename Tail >                                                                                           \
  struct get_type_info< ::ReaK::rpc::detail::VARIABLE, Tail > {                                                        \
    typedef type_id< ::ReaK::rpc::detail::VARIABLE, typename Tail::type > type;                                        \
    BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< std::string >::type_name + get_type_name_tail< Tail >::value; \
  };

#else

#define RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS( VARIABLE, CLASSID )                      \
  template <>                                                                                 \
  struct get_type_id< ::ReaK::rpc::detail::VARIABLE > {                                       \
    BOOST_STATIC_CONSTANT( unsigned int, ID = CLASSID );                                      \
    static const char* type_name() BOOST_NOEXCEPT { return BOOST_PP_STRINGIZE( VARIABLE ); }; \
    static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };                      \
                                                                                              \
    typedef const serializable& save_type;                                                    \
    typedef serializable& load_type;                                                          \
  };                                                                                          \
                                                                                              \
  template < typename Tail >                                                                  \
  struct get_type_info< ::ReaK::rpc::detail::VARIABLE, Tail > {                               \
    typedef type_id< ::ReaK::rpc::detail::VARIABLE, typename Tail::type > type;               \
    static std::string type_name() {                                                          \
      std::string result = get_type_id< ::ReaK::rpc::detail::VARIABLE >::type_name();         \
      result += get_type_name_tail< Tail >::value();                                          \
      return result;                                                                          \
    };                                                                                        \
  };

#endif

RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS( disc_basic_header, 0xC1600001 )
RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS( disc_with_port_header, 0xC1600002 )

#undef RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS
};
};


#endif

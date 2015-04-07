/**
 * \file rpc_headers.hpp
 *
 * This library declares helpers to generate and parse the headers of the ReaK.RPC protocol.
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


#ifndef REAK_RPC_HEADERS_HPP
#define REAK_RPC_HEADERS_HPP

#include <ReaK/core/serialization/archiver.hpp>
#include <ReaK/core/base/serializable.hpp>

#include <ReaK/core/rpc/version.hpp>

#include <string>
#include <tuple>
#include <iostream>

namespace ReaK {

namespace rpc {

namespace detail {


enum rpc_header_type { rpc_call = 0, rpc_return, rpc_exception, rpc_unrecognized, rpc_startportserver };

extern const char* rpc_header_type_to_str[];

std::string generate_rpc_header( rpc_header_type h, msg_format fmt, std::size_t msg_size );

std::tuple< rpc_header_type, msg_format, std::size_t > parse_rpc_header( std::istream& hdr_in );


struct rpc_basic_header : serializable {
  std::string name;
  std::size_t params_hash;
  std::size_t call_seq;

  rpc_basic_header( const std::string& aName = "", std::size_t aParamsHash = 0, std::size_t aCallSeq = 0 )
      : name( aName ), params_hash( aParamsHash ), call_seq( aCallSeq ){};


  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    A& RK_SERIAL_SAVE_WITH_NAME( name ) & RK_SERIAL_SAVE_WITH_NAME( params_hash )
      & RK_SERIAL_SAVE_WITH_NAME( call_seq );
  };
  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    A& RK_SERIAL_LOAD_WITH_NAME( name ) & RK_SERIAL_LOAD_WITH_NAME( params_hash )
      & RK_SERIAL_LOAD_WITH_NAME( call_seq );
  };

  RK_RTTI_REGISTER_CLASS_1BASE( rpc_basic_header, 1, serializable )
};

inline bool operator<( const rpc_basic_header& lhs, const rpc_basic_header& rhs ) {
  if( lhs.name < rhs.name )
    return true;
  if( lhs.name == rhs.name ) {
    if( lhs.params_hash < rhs.params_hash )
      return true;
    if( ( lhs.params_hash == rhs.params_hash ) && ( lhs.call_seq < rhs.call_seq ) )
      return true;
  };
  return false;
};


struct rpc_call_header : rpc_basic_header {
  unsigned int reply_port;

  rpc_call_header( const std::string& aName = "", std::size_t aParamsHash = 0, std::size_t aCallSeq = 0,
                   unsigned int aReplyPort = 0 )
      : rpc_basic_header( aName, aParamsHash, aCallSeq ), reply_port( aReplyPort ){};


  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    rpc_basic_header::save( A, 1 );
    A& RK_SERIAL_SAVE_WITH_NAME( reply_port );
  };
  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    rpc_basic_header::load( A, 1 );
    A& RK_SERIAL_LOAD_WITH_NAME( reply_port );
  };

  RK_RTTI_REGISTER_CLASS_1BASE( rpc_call_header, 1, rpc_basic_header )
};

// CALL {
//   string name
//   uint32 params_hash
//   uint32 call_seq
//   uint32 reply_port
//   [input_params]
// }


struct rpc_return_header : rpc_basic_header {

  rpc_return_header( const std::string& aName = "", std::size_t aParamsHash = 0, std::size_t aCallSeq = 0 )
      : rpc_basic_header( aName, aParamsHash, aCallSeq ){};


  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const { rpc_basic_header::save( A, 1 ); };
  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) { rpc_basic_header::load( A, 1 ); };

  RK_RTTI_REGISTER_CLASS_1BASE( rpc_return_header, 1, rpc_basic_header )
};

// RETURN {
//   string name
//   uint32 params_hash
//   uint32 call_seq
//   [return]
//   [output_params]
// }


struct rpc_exception_header : rpc_basic_header {
  std::string except_type;
  std::string except_msg;

  rpc_exception_header( const std::string& aName = "", std::size_t aParamsHash = 0, std::size_t aCallSeq = 0,
                        const std::string& aExceptType = "", const std::string& aExceptMsg = "" )
      : rpc_basic_header( aName, aParamsHash, aCallSeq ), except_type( aExceptType ), except_msg( aExceptMsg ){};


  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    rpc_basic_header::save( A, 1 );
    A& RK_SERIAL_SAVE_WITH_NAME( except_type ) & RK_SERIAL_SAVE_WITH_NAME( except_msg );
  };
  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    rpc_basic_header::load( A, 1 );
    A& RK_SERIAL_LOAD_WITH_NAME( except_type ) & RK_SERIAL_LOAD_WITH_NAME( except_msg );
  };

  RK_RTTI_REGISTER_CLASS_1BASE( rpc_exception_header, 1, rpc_basic_header )
};

// EXCEPTION {
//   string name
//   uint32 params_hash
//   uint32 call_seq
//   string except_type
//   string except_msg
// }

struct rpc_unrecognized_header : rpc_basic_header {

  rpc_unrecognized_header( const std::string& aName = "", std::size_t aParamsHash = 0, std::size_t aCallSeq = 0 )
      : rpc_basic_header( aName, aParamsHash, aCallSeq ){};


  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const { rpc_basic_header::save( A, 1 ); };
  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) { rpc_basic_header::load( A, 1 ); };

  RK_RTTI_REGISTER_CLASS_1BASE( rpc_unrecognized_header, 1, rpc_basic_header )
};

// UNRECOGNIZED {
//   string name
//   uint32 params_hash
//   uint32 call_seq
// }


struct rpc_startportserver_header : serializable {

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {};
  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ){};

  RK_RTTI_REGISTER_CLASS_1BASE( rpc_startportserver_header, 1, serializable )
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

RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS( rpc_basic_header, 0xC1600001 )
RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS( rpc_call_header, 0xC1600002 )
RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS( rpc_return_header, 0xC1600003 )
RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS( rpc_exception_header, 0xC1600004 )
RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS( rpc_unrecognized_header, 0xC1600005 )

#undef RK_RTTI_MAKE_TYPE_INFO_FOR_RPC_DETAIL_CLASS
};
};


#endif

/**
 * \file endian_conversions.hpp
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

#ifndef REAK_ENDIAN_CONVERSIONS_HPP
#define REAK_ENDIAN_CONVERSIONS_HPP

#include "defs.hpp"

#include <boost/cstdint.hpp>

#ifdef WIN32
#include <winsock2.h>
#else
#include <netinet/in.h>
#endif

namespace ReaK {

  
union short_to_ushort {
  boost::int16_t  i16;
  boost::uint16_t ui16;
};


union float_to_ulong {
  float    f;
  boost::uint32_t ui32;
};

union long_to_ulong {
  boost::int32_t  i32;
  boost::uint32_t ui32;
};

union double_to_ulong {
  double   d;
  boost::uint32_t ui32[2];
};

union llong_to_ulong {
#ifndef BOOST_NO_INT64_T
  boost::int64_t  i64;
  boost::uint64_t ui64;
#endif
  boost::uint32_t ui32[2];
};


template <typename UnionT>
void ntoh_1ui16(UnionT& value) {
#if RK_BYTE_ORDER != RK_ORDER_BIG_ENDIAN
  value.ui16 = ntohs(value.ui16);
#endif
};

template <typename UnionT>
void hton_1ui16(UnionT& value) {
#if RK_BYTE_ORDER != RK_ORDER_BIG_ENDIAN
  value.ui16 = htons(value.ui16);
#endif
};

template <typename UnionT>
void ntoh_1ui32(UnionT& value) {
#if RK_BYTE_ORDER != RK_ORDER_BIG_ENDIAN
  value.ui32 = ntohl(value.ui32);
#endif
};

template <typename UnionT>
void hton_1ui32(UnionT& value) {
#if RK_BYTE_ORDER != RK_ORDER_BIG_ENDIAN
  value.ui32 = htonl(value.ui32);
#endif
};

template <typename UnionT>
void ntoh_2ui32(UnionT& value) {
#if RK_BYTE_ORDER == RK_ORDER_LITTLE_ENDIAN
  uint32_t tmp = ntohl(value.ui32[0]);
  value.ui32[0] = ntohl(value.ui32[1]);
  value.ui32[1] = tmp;
#endif
  // NOTE: for 64-bit values, there is no point in supporting PDP-endianness, as 64-bit values are not supported by PDP platforms.
};

template <typename UnionT>
void hton_2ui32(UnionT& value) {
#if RK_BYTE_ORDER == RK_ORDER_LITTLE_ENDIAN
  uint32_t tmp = htonl(value.ui32[0]);
  value.ui32[0] = htonl(value.ui32[1]);
  value.ui32[1] = tmp;
#endif
  // NOTE: for 64-bit values, there is no point in supporting PDP-endianness, as 64-bit values are not supported by PDP platforms.
};


// Overloaded versions for different primitive types:

inline void ntoh_any(boost::uint8_t&) { };
inline void ntoh_any(boost::int8_t&) { };

inline void hton_any(boost::uint8_t&) { };
inline void hton_any(boost::int8_t&) { };

inline void hton_any(boost::uint16_t& s) {
  short_to_ushort tmp; tmp.ui16 = s;
  hton_1ui16(tmp);
  s = tmp.ui16;
};
inline void hton_any(boost::int16_t& s) { 
  short_to_ushort tmp; tmp.i16 = s;
  hton_1ui16(tmp);
  s = tmp.i16;
};

inline void ntoh_any(boost::uint16_t& s) {
  short_to_ushort tmp; tmp.ui16 = s;
  ntoh_1ui16(tmp);
  s = tmp.ui16;
};
inline void ntoh_any(boost::int16_t& s) { 
  short_to_ushort tmp; tmp.i16 = s;
  ntoh_1ui16(tmp);
  s = tmp.i16;
};

inline void hton_any(boost::uint32_t& i) {
  long_to_ulong tmp; tmp.ui32 = i;
  hton_1ui32(tmp);
  i = tmp.ui32;
};
inline void hton_any(boost::int32_t& i) { 
  long_to_ulong tmp; tmp.i32 = i;
  hton_1ui32(tmp);
  i = tmp.i32;
};

inline void ntoh_any(boost::uint32_t& i) {
  long_to_ulong tmp; tmp.ui32 = i;
  ntoh_1ui32(tmp);
  i = tmp.ui32;
};
inline void ntoh_any(boost::int32_t& i) { 
  long_to_ulong tmp; tmp.i32 = i;
  ntoh_1ui32(tmp);
  i = tmp.i32;
};

#ifndef BOOST_NO_INT64_T
inline void hton_any(boost::uint64_t& i) {
  llong_to_ulong tmp; tmp.ui64 = i;
  hton_2ui32(tmp);
  i = tmp.ui64;
};
inline void hton_any(boost::int64_t& i) { 
  llong_to_ulong tmp; tmp.i64 = i;
  hton_2ui32(tmp);
  i = tmp.i64;
};

inline void ntoh_any(boost::uint64_t& i) {
  llong_to_ulong tmp; tmp.ui64 = i;
  ntoh_2ui32(tmp);
  i = tmp.ui64;
};
inline void ntoh_any(boost::int64_t& i) { 
  llong_to_ulong tmp; tmp.i64 = i;
  ntoh_2ui32(tmp);
  i = tmp.i64;
};
#endif



inline void hton_any(float& f) { 
  float_to_ulong tmp; tmp.f = f;
  hton_1ui32(tmp);
  f = tmp.f;
};

inline void ntoh_any(float& f) { 
  float_to_ulong tmp; tmp.f = f;
  ntoh_1ui32(tmp);
  f = tmp.f;
};

inline void hton_any(double& d) { 
  double_to_ulong tmp; tmp.d = d;
  hton_2ui32(tmp);
  d = tmp.d;
};

inline void ntoh_any(double& d) { 
  double_to_ulong tmp; tmp.d = d;
  ntoh_2ui32(tmp);
  d = tmp.d;
};


};


#endif


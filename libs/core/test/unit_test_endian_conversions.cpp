
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

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/endian_conversions.hpp>

#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE endian_conversions
#include <boost/test/unit_test.hpp>

union wrap_uint8_t {
  boost::uint8_t u;
  boost::int8_t i;
  unsigned char c[1];
};

union wrap_uint16_t {
  boost::uint16_t u;
  boost::int16_t i;
  unsigned char c[2];
};

union wrap_uint32_t {
  boost::uint32_t u;
  boost::int32_t i;
  unsigned char c[4];
};

#ifndef BOOST_NO_INT64_T
union wrap_uint64_t {
  boost::uint64_t u;
  boost::int64_t i;
  unsigned char c[8];
};
#endif

BOOST_AUTO_TEST_CASE(endian_conversions_test) {
  using namespace ReaK;

  // Tests for little-endian systems:
  if (1000 != htonl(1000)) {

    wrap_uint8_t u8;
    u8.u = 0x0A;
    hton_any(u8.u);
    if (u8.c[0] != 0x0A)
      BOOST_ERROR("uint8 conversion to network order failed!");
    u8.i = 0x0A;
    hton_any(u8.i);
    if (u8.c[0] != 0x0A)
      BOOST_ERROR("int8 conversion to network order failed!");

    wrap_uint16_t u16;
    u16.u = 0x0A0B;
    hton_any(u16.u);
    if ((u16.c[0] != 0x0A) || (u16.c[1] != 0x0B))
      BOOST_ERROR("uint16 conversion to network order failed!");
    u16.i = 0x0A0B;
    hton_any(u16.i);
    if ((u16.c[0] != 0x0A) || (u16.c[1] != 0x0B))
      BOOST_ERROR("int16 conversion to network order failed!");

    wrap_uint32_t u32;
    u32.u = 0x0A0B0C0D;
    hton_any(u32.u);
    if ((u32.c[0] != 0x0A) || (u32.c[1] != 0x0B) || (u32.c[2] != 0x0C) ||
        (u32.c[3] != 0x0D))
      BOOST_ERROR("uint32 conversion to network order failed!");
    u32.i = 0x0A0B0C0D;
    hton_any(u32.i);
    if ((u32.c[0] != 0x0A) || (u32.c[1] != 0x0B) || (u32.c[2] != 0x0C) ||
        (u32.c[3] != 0x0D))
      BOOST_ERROR("int32 conversion to network order failed!");

#ifndef BOOST_NO_INT64_T
    wrap_uint64_t u64;
    u64.u = 0x0A0B0C0D0E0F;
    hton_any(u64.u);
    if ((u64.c[0] != 0x00) || (u64.c[1] != 0x00) || (u64.c[2] != 0x0A) ||
        (u64.c[3] != 0x0B) || (u64.c[4] != 0x0C) || (u64.c[5] != 0x0D) ||
        (u64.c[6] != 0x0E) || (u64.c[7] != 0x0F))
      BOOST_ERROR("uint64 conversion to network order failed!");
    u64.i = 0x0A0B0C0D0E0F;
    hton_any(u64.i);
    if ((u64.c[0] != 0x00) || (u64.c[1] != 0x00) || (u64.c[2] != 0x0A) ||
        (u64.c[3] != 0x0B) || (u64.c[4] != 0x0C) || (u64.c[5] != 0x0D) ||
        (u64.c[6] != 0x0E) || (u64.c[7] != 0x0F))
      BOOST_ERROR("int64 conversion to network order failed!");
#endif
  };

  boost::uint8_t u8 = 0x0A;
  ntoh_any(u8);
  hton_any(u8);
  BOOST_CHECK_EQUAL(u8, 0x0A);

  boost::int8_t i8 = 0x0A;
  ntoh_any(i8);
  hton_any(i8);
  BOOST_CHECK_EQUAL(i8, 0x0A);

  boost::uint16_t u16 = 0x0A0B;
  hton_any(u16);
  ntoh_any(u16);
  BOOST_CHECK_EQUAL(u16, 0x0A0B);

  boost::int16_t i16 = 0x0A0B;
  hton_any(i16);
  ntoh_any(i16);
  BOOST_CHECK_EQUAL(i16, 0x0A0B);

  boost::uint32_t u32 = 0x0A0B0C0D;
  hton_any(u32);
  ntoh_any(u32);
  BOOST_CHECK_EQUAL(u32, 0x0A0B0C0D);

  boost::int32_t i32 = 0x0A0B0C0D;
  hton_any(i32);
  ntoh_any(i32);
  BOOST_CHECK_EQUAL(i32, 0x0A0B0C0D);

#ifndef BOOST_NO_INT64_T
  boost::uint64_t u64 = 0x0A0B0C0D0E0F;
  hton_any(u64);
  ntoh_any(u64);
  BOOST_CHECK_EQUAL(u64, 0x0A0B0C0D0E0F);

  boost::int64_t i64 = 0x0A0B0C0D0E0F;
  hton_any(i64);
  ntoh_any(i64);
  BOOST_CHECK_EQUAL(i64, 0x0A0B0C0D0E0F);
#endif

  float f(3.14159);
  hton_any(f);
  ntoh_any(f);
  BOOST_CHECK_EQUAL(f, float(3.14159));

  double d(3.14159);
  hton_any(d);
  ntoh_any(d);
  BOOST_CHECK_EQUAL(d, double(3.14159));
};

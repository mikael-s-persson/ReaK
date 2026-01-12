
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

#include "ReaK/core/base/endian_conversions.h"
#include <netinet/in.h>

#include <cstdint>

#include "gtest/gtest.h"

namespace ReaK {
namespace {

union wrap_uint8_t {
  std::uint8_t u;
  std::int8_t i;
  unsigned char c[1];  // NOLINT
};

union wrap_uint16_t {
  std::uint16_t u;
  std::int16_t i;
  unsigned char c[2];  // NOLINT
};

union wrap_uint32_t {
  std::uint32_t u;
  std::int32_t i;
  unsigned char c[4];  // NOLINT
};

union wrap_uint64_t {
  std::uint64_t u;
  std::int64_t i;
  unsigned char c[8];  // NOLINT
};

TEST(EndianConversionTests, AllCases) {
  // Tests for other-endian systems:
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunreachable-code"
  constexpr auto other_en = (std::endian::native == std::endian::little
                                 ? std::endian::big
                                 : std::endian::little);
#pragma GCC diagnostic pop

  // NOLINTBEGIN(cppcoreguidelines-pro-type-union-access)
  wrap_uint8_t v8{};
  v8.u = 0x0A;
  to_endian<other_en>(v8.u);
  EXPECT_EQ(v8.c[0], 0x0A) << "uint8 conversion to other endianness failed!";

  v8.i = 0x0A;
  to_endian<other_en>(v8.i);
  EXPECT_EQ(v8.c[0], 0x0A) << "int8 conversion to other endianness failed!";

  wrap_uint16_t v16{};
  v16.u = 0x0A0B;
  to_endian<other_en>(v16.u);
  EXPECT_EQ(v16.c[0], 0x0A) << "uint16 conversion to other endianness failed!";
  EXPECT_EQ(v16.c[1], 0x0B) << "uint16 conversion to other endianness failed!";

  v16.i = 0x0A0B;
  to_endian<other_en>(v16.i);
  EXPECT_EQ(v16.c[0], 0x0A) << "int16 conversion to other endianness failed!";
  EXPECT_EQ(v16.c[1], 0x0B) << "int16 conversion to other endianness failed!";

  wrap_uint32_t v32{};
  v32.u = 0x0A0B0C0D;
  to_endian<other_en>(v32.u);
  EXPECT_EQ(v32.c[0], 0x0A) << "uint32 conversion to other endianness failed!";
  EXPECT_EQ(v32.c[1], 0x0B) << "uint32 conversion to other endianness failed!";
  EXPECT_EQ(v32.c[2], 0x0C) << "uint32 conversion to other endianness failed!";
  EXPECT_EQ(v32.c[3], 0x0D) << "uint32 conversion to other endianness failed!";

  v32.i = 0x0A0B0C0D;
  to_endian<other_en>(v32.i);
  EXPECT_EQ(v32.c[0], 0x0A) << "int32 conversion to other endianness failed!";
  EXPECT_EQ(v32.c[1], 0x0B) << "int32 conversion to other endianness failed!";
  EXPECT_EQ(v32.c[2], 0x0C) << "int32 conversion to other endianness failed!";
  EXPECT_EQ(v32.c[3], 0x0D) << "int32 conversion to other endianness failed!";

  wrap_uint64_t v64{};
  v64.u = 0x0A0B0C0D0E0FULL;
  to_endian<other_en>(v64.u);
  EXPECT_EQ(v64.c[0], 0x00) << "uint64 conversion to other endianness failed!";
  EXPECT_EQ(v64.c[1], 0x00) << "uint64 conversion to other endianness failed!";
  EXPECT_EQ(v64.c[2], 0x0A) << "uint64 conversion to other endianness failed!";
  EXPECT_EQ(v64.c[3], 0x0B) << "uint64 conversion to other endianness failed!";
  EXPECT_EQ(v64.c[4], 0x0C) << "uint64 conversion to other endianness failed!";
  EXPECT_EQ(v64.c[5], 0x0D) << "uint64 conversion to other endianness failed!";
  EXPECT_EQ(v64.c[6], 0x0E) << "uint64 conversion to other endianness failed!";
  EXPECT_EQ(v64.c[7], 0x0F) << "uint64 conversion to other endianness failed!";

  v64.i = 0x0A0B0C0D0E0FLL;
  to_endian<other_en>(v64.i);
  EXPECT_EQ(v64.c[0], 0x00) << "int64 conversion to other endianness failed!";
  EXPECT_EQ(v64.c[1], 0x00) << "int64 conversion to other endianness failed!";
  EXPECT_EQ(v64.c[2], 0x0A) << "int64 conversion to other endianness failed!";
  EXPECT_EQ(v64.c[3], 0x0B) << "int64 conversion to other endianness failed!";
  EXPECT_EQ(v64.c[4], 0x0C) << "int64 conversion to other endianness failed!";
  EXPECT_EQ(v64.c[5], 0x0D) << "int64 conversion to other endianness failed!";
  EXPECT_EQ(v64.c[6], 0x0E) << "int64 conversion to other endianness failed!";
  EXPECT_EQ(v64.c[7], 0x0F) << "int64 conversion to other endianness failed!";
  // NOLINTEND(cppcoreguidelines-pro-type-union-access)

  std::uint8_t u8 = 0x0A;
  to_endian<other_en>(u8);
  from_endian<other_en>(u8);
  EXPECT_EQ(u8, 0x0A);

  std::int8_t i8 = 0x0A;
  to_endian<other_en>(i8);
  from_endian<other_en>(i8);
  EXPECT_EQ(i8, 0x0A);

  std::uint16_t u16 = 0x0A0B;
  to_endian<other_en>(u16);
  from_endian<other_en>(u16);
  EXPECT_EQ(u16, 0x0A0B);

  std::int16_t i16 = 0x0A0B;
  to_endian<other_en>(i16);
  from_endian<other_en>(i16);
  EXPECT_EQ(i16, 0x0A0B);

  std::uint32_t u32 = 0x0A0B0C0D;
  to_endian<other_en>(u32);
  from_endian<other_en>(u32);
  EXPECT_EQ(u32, 0x0A0B0C0D);

  std::int32_t i32 = 0x0A0B0C0D;
  to_endian<other_en>(i32);
  from_endian<other_en>(i32);
  EXPECT_EQ(i32, 0x0A0B0C0D);

  std::uint64_t u64 = 0x0A0B0C0D0E0FULL;
  to_endian<other_en>(u64);
  from_endian<other_en>(u64);
  EXPECT_EQ(u64, 0x0A0B0C0D0E0F);

  std::int64_t i64 = 0x0A0B0C0D0E0FLL;
  to_endian<other_en>(i64);
  from_endian<other_en>(i64);
  EXPECT_EQ(i64, 0x0A0B0C0D0E0F);

  float f(3.1F);
  to_endian<other_en>(f);
  from_endian<other_en>(f);
  EXPECT_EQ(f, 3.1F);

  double d(3.1);
  to_endian<other_en>(d);
  from_endian<other_en>(d);
  EXPECT_EQ(d, 3.1);
}

}  // namespace
}  // namespace ReaK

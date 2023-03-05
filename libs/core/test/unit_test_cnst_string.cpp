
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
#include "ReaK/core/base/cnst_string.hpp"

#include "ReaK/core/base/defs.hpp"

#include "gtest/gtest.h"

namespace ReaK {
namespace {

TEST(CnstStringTest, AllCases) {
  constexpr auto s1 = std::string_view{"Hello World!"};

  EXPECT_EQ(s1[0], 'H') << "Initialization from a string literal has failed!";
  EXPECT_EQ(s1[6], 'W') << "Initialization from a string literal has failed!";
  EXPECT_EQ(s1[12], '\0') << "Initialization from a string literal has failed!";

  EXPECT_EQ(s1.size(), 12);

  static constexpr auto hello = std::string_view{"Hello "};
  static constexpr auto world = std::string_view{"World!"};
  static constexpr auto nice_to_see_you = std::string_view{" Nice to see you!"};
  constexpr auto s3 = ct_concat_v<hello, world, nice_to_see_you>;
  const std::string_view compare_s3 = "Hello World! Nice to see you!";

  EXPECT_EQ(s3[0], compare_s3[0]) << "Initialization from a concatenation of "
                                     "wrapped string literals has failed!";
  EXPECT_EQ(s3[6], compare_s3[6]) << "Initialization from a concatenation of "
                                     "wrapped string literals has failed!";
  EXPECT_EQ(s3[13], compare_s3[13]) << "Initialization from a concatenation of "
                                       "wrapped string literals has failed!";

  EXPECT_EQ(s3.size(), compare_s3.size());

  constexpr std::size_t s3_hash = fnv_1a_hash(s3);
  EXPECT_EQ(s3_hash, fnv_1a_hash(compare_s3));

  constexpr auto h_54582 = ct_itoa_v<54582>;

  EXPECT_EQ(h_54582[0], '5')
      << "Initialization from a integer conversion has failed!";
  EXPECT_EQ(h_54582[1], '4')
      << "Initialization from a integer conversion has failed!";
  EXPECT_EQ(h_54582[2], '5')
      << "Initialization from a integer conversion has failed!";
  EXPECT_EQ(h_54582[3], '8')
      << "Initialization from a integer conversion has failed!";
  EXPECT_EQ(h_54582[4], '2')
      << "Initialization from a integer conversion has failed!";

  std::string s3_str{s3};
  EXPECT_EQ(s3_str, compare_s3);
}

}  // namespace
}  // namespace ReaK

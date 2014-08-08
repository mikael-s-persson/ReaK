
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

#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE cnst_string
#include <boost/test/unit_test.hpp>

#if (!defined(BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX) && !defined(BOOST_NO_CXX11_HDR_ARRAY) && !defined(BOOST_NO_CXX11_AUTO_DECLARATIONS))

#include <ReaK/core/base/cnst_string.hpp>


BOOST_AUTO_TEST_CASE( cnst_string_test )
{
  using namespace ReaK;
  
  BOOST_CONSTEXPR auto s1 = RK_LSA("Hello World!");
  
  if((s1[0] != 'H') || (s1[6] != 'W') || (s1[12] != '\0'))
    BOOST_ERROR("Initialization from a string literal has failed!");
  
  BOOST_CHECK_EQUAL( s1.size(), 12 );
  
  BOOST_CONSTEXPR auto s3 = RK_LSA("Hello ") + RK_LSA("World!") + RK_LSA(" Nice to see you!");
  const char compare_s3[] = "Hello World! Nice to see you!";
  
  if((s3[0] != compare_s3[0]) || (s3[6] != compare_s3[6]) || (s3[13] != compare_s3[13]) || (s3[sizeof(compare_s3)-1] != '\0'))
    BOOST_ERROR("Initialization from a concatenation of wrapped string literals has failed!");
  
  BOOST_CHECK_EQUAL( s3.size(), (sizeof(compare_s3)-1) );
  
  BOOST_CONSTEXPR unsigned int s3_size = s3.size();
  BOOST_CHECK_EQUAL( s3_size, (sizeof(compare_s3)-1) );
  
  BOOST_CONSTEXPR char s3_c13 = s3[13];
  BOOST_CHECK_EQUAL( s3_c13, compare_s3[13] );
  
  BOOST_CONSTEXPR std::size_t s3_hash = s3.hash();
  BOOST_CHECK_EQUAL( s3_hash, fnv_1a_hash(compare_s3));
  
#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES
  
  BOOST_CONSTEXPR auto h_54582 = ct_itoa<54582>::text;
  
  if((h_54582[0] != '5') || (h_54582[1] != '4') || (h_54582[2] != '5') || 
     (h_54582[3] != '8') || (h_54582[4] != '2') || (h_54582[5] != '\0'))
    BOOST_ERROR("Initialization from a integer conversion has failed!");
  
#endif
  
  std::string s3_str = s3.to_string();
  BOOST_CHECK_EQUAL( s3_str, compare_s3 );
  
};

#else

BOOST_AUTO_TEST_CASE( cnst_string_test )
{
  
};

#endif













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

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/lin_alg/mat_alg.hpp>

#include <iostream>
#include <fstream>
#include <cstdio>

#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE mat_alg
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


BOOST_AUTO_TEST_CASE( mat_operator_tests )
{
  using namespace ReaK;
  using std::fabs;
  
  mat<double,mat_structure::rectangular> m1234(1.0,2.0,3.0,4.0);
  BOOST_CHECK( ( ( fabs(m1234(0,0) - 1.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m1234(0,1) - 2.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m1234(1,0) - 3.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m1234(1,1) - 4.0) < std::numeric_limits<double>::epsilon() ) ) );
  
  m1234 = transpose(m1234);
  BOOST_CHECK( ( ( fabs(m1234(0,0) - 1.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m1234(0,1) - 3.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m1234(1,0) - 2.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m1234(1,1) - 4.0) < std::numeric_limits<double>::epsilon() ) ) );
  
  m1234 = transpose_move(m1234);
  BOOST_CHECK( ( ( fabs(m1234(0,0) - 1.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m1234(0,1) - 2.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m1234(1,0) - 3.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m1234(1,1) - 4.0) < std::numeric_limits<double>::epsilon() ) ) );
  
  mat<double,mat_structure::rectangular> m_ident(2,2,true);
  BOOST_CHECK( ( ( fabs(m_ident(0,0) - 1.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m_ident(0,1)) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m_ident(1,0)) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m_ident(1,1) - 1.0) < std::numeric_limits<double>::epsilon() ) ) );
  
  BOOST_CHECK( ( is_identity_mat(m_ident,std::numeric_limits<double>::epsilon()) ) );
  
  mat_identity<double>::type m_ident2(2);
  BOOST_CHECK( ( is_identity_mat(m_ident2,std::numeric_limits<double>::epsilon()) ) );
  
  mat_null<double>::type m_zeroes(2,2);
  BOOST_CHECK( ( is_null_mat(m_zeroes,std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_null_mat(m_zeroes * m_ident2,std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_null_mat(m_ident2 * m_zeroes,std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_null_mat(m_zeroes * m1234,std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_null_mat(m1234 * m_zeroes,std::numeric_limits<double>::epsilon()) ) );
  
  mat<double,mat_structure::rectangular> m1234_inv(-2.0, 1.0, 1.5, -0.5);
  
  //****************  Rectangular vs. Rectangular  **************************
  BOOST_CHECK( ( is_identity_mat(m1234_inv * m1234,std::numeric_limits<double>::epsilon()) ) );
  
  BOOST_CHECK( ( ( is_identity_mat(0.5 * (m1234_inv * m1234 + m1234_inv * m1234),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat((m1234_inv * m1234 + m1234_inv * m1234) * 0.5,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat((m1234_inv * m1234 - m1234_inv * m1234) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_null_mat(m1234_inv * m1234 - m1234_inv * m1234,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_null_mat(m1234_inv * m1234 + (-(m1234_inv * m1234)),std::numeric_limits<double>::epsilon()) ) ) );
  
  //****************  Square vs. Square  **************************
  mat<double,mat_structure::square> m1234_sqr(1.0,2.0,3.0,4.0);
  mat<double,mat_structure::square> m1234_inv_sqr(-2.0,1.0,1.5,-0.5);
  BOOST_CHECK( ( is_identity_mat(m1234_inv_sqr * m1234_sqr,std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( ( is_identity_mat(0.5 * (m1234_inv_sqr * m1234_sqr + m1234_inv_sqr * m1234_sqr),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat((m1234_inv_sqr * m1234_sqr + m1234_inv_sqr * m1234_sqr) * 0.5,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat((m1234_inv_sqr * m1234_sqr - m1234_inv_sqr * m1234_sqr) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_null_mat(m1234_inv_sqr * m1234_sqr - m1234_inv_sqr * m1234_sqr,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_null_mat(m1234_inv_sqr * m1234_sqr + (-(m1234_inv_sqr * m1234_sqr)),std::numeric_limits<double>::epsilon()) ) ) );
  
  //****************  Square vs. Rectangular  **************************
  BOOST_CHECK( ( ( is_identity_mat(m1234_inv * ((m1234_inv_sqr * m1234_sqr) * m1234),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m1234_inv * (m1234 * (m1234_inv_sqr * m1234_sqr)),std::numeric_limits<double>::epsilon()) ) ) );
  BOOST_CHECK( ( ( is_identity_mat(0.5 * (m1234 * m1234_inv + (m1234_inv_sqr * m1234_sqr)),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(0.5 * ((m1234_inv_sqr * m1234_sqr) + m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat((m1234_inv_sqr * m1234_sqr) + (m1234 * m1234_inv - (m1234_inv_sqr * m1234_sqr)),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat((m1234_inv_sqr * m1234_sqr) + ((m1234_inv_sqr * m1234_sqr) - m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) ) );
  
  //****************  Symmetric vs. Symmetric  **************************
  mat<double,mat_structure::symmetric> m123_sym(1.0,2.0,3.0);
  mat<double,mat_structure::symmetric> m123_inv_sym(-3.0,2.0,-1.0);
  BOOST_CHECK( ( is_identity_mat(m123_inv_sym * m123_sym,std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( ( is_identity_mat(0.5 * (m123_inv_sym * m123_sym + m123_inv_sym * m123_sym),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat((m123_inv_sym * m123_sym + m123_inv_sym * m123_sym) * 0.5,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat((m123_inv_sym * m123_sym - m123_inv_sym * m123_sym) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_null_mat(m123_inv_sym * m123_sym - m123_inv_sym * m123_sym,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_null_mat(m123_inv_sym * m123_sym + (-(m123_inv_sym * m123_sym)),std::numeric_limits<double>::epsilon()) ) ) ); 
  
  //****************  Symmetric vs. Rectangular  **************************
  BOOST_CHECK( ( ( is_identity_mat(m1234_inv * ((m123_sym * m123_inv_sym) * m1234),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m1234_inv * (m1234 * (m123_sym * m123_inv_sym)),std::numeric_limits<double>::epsilon()) ) ) );
  BOOST_CHECK( ( ( is_identity_mat(0.5 * (m1234 * m1234_inv + (m123_sym * m123_inv_sym)),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(0.5 * ((m123_sym * m123_inv_sym) + m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat((m123_sym * m123_inv_sym) + (m1234 * m1234_inv - (m123_sym * m123_inv_sym)),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat((m123_sym * m123_inv_sym) + ((m123_sym * m123_inv_sym) - m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) ) );
  
  //****************  Symmetric vs. Square  **************************
  BOOST_CHECK( ( ( is_identity_mat(m1234_inv_sqr * ((m123_sym * m123_inv_sym) * m1234_sqr),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m1234_inv_sqr * (m1234_sqr * (m123_sym * m123_inv_sym)),std::numeric_limits<double>::epsilon()) ) ) );
  BOOST_CHECK( ( ( is_identity_mat(0.5 * (m1234_sqr * m1234_inv_sqr + (m123_sym * m123_inv_sym)),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(0.5 * ((m123_sym * m123_inv_sym) + m1234_sqr * m1234_inv_sqr),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat((m123_sym * m123_inv_sym) + (m1234_sqr * m1234_inv_sqr - (m123_sym * m123_inv_sym)),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat((m123_sym * m123_inv_sym) + ((m123_sym * m123_inv_sym) - m1234_sqr * m1234_inv_sqr),std::numeric_limits<double>::epsilon()) ) ) );
  
  //****************  Rectangular vs. Diagonal  **************************
  mat<double,mat_structure::diagonal> m_ident_diag(2,1.0);
  BOOST_CHECK( ( ( is_identity_mat(m1234_inv * (m_ident_diag * m1234),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m1234_inv * (m1234 * m_ident_diag),std::numeric_limits<double>::epsilon()) ) ) );
  BOOST_CHECK( ( ( is_identity_mat(0.5 * (m1234 * m1234_inv + m_ident_diag),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(0.5 * (m_ident_diag + m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m_ident_diag + (m1234 * m1234_inv - m_ident_diag),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m_ident_diag + (m_ident_diag - m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) ) );
  
  //****************  Rectangular vs. Scalar  **************************
  mat<double,mat_structure::scalar> m_ident_scalar(2,1.0);
  BOOST_CHECK( ( ( is_identity_mat(m1234_inv * (m_ident_scalar * m1234),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m1234_inv * (m1234 * m_ident_scalar),std::numeric_limits<double>::epsilon()) ) ) );
  BOOST_CHECK( ( ( is_identity_mat(0.5 * (m1234 * m1234_inv + m_ident_scalar),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(0.5 * (m_ident_scalar + m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m_ident_scalar + (m1234 * m1234_inv - m_ident_scalar),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m_ident_scalar + (m_ident_scalar - m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) ) );
  
  //****************  Rectangular vs. Identity  **************************
  BOOST_CHECK( ( ( is_identity_mat(m1234_inv * (m_ident2 * m1234),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m1234_inv * (m1234 * m_ident2),std::numeric_limits<double>::epsilon()) ) ) );
  BOOST_CHECK( ( ( is_identity_mat(0.5 * (m1234 * m1234_inv + m_ident2),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(0.5 * (m_ident2 + m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m_ident2 + (m1234 * m1234_inv - m_ident2),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m_ident2 + (m_ident2 - m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) ) );
  
  //****************  Rectangular vs. Null  **************************
  BOOST_CHECK( ( ( is_null_mat(m1234 * m_zeroes,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_null_mat(m_zeroes * m1234,std::numeric_limits<double>::epsilon()) ) ) );
  BOOST_CHECK( ( ( is_identity_mat(m_zeroes - m1234_inv * (-m1234),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m1234_inv * m1234 - m_zeroes,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m1234_inv * m1234 + m_zeroes,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m_zeroes + m1234_inv * m1234,std::numeric_limits<double>::epsilon()) ) ) );
  
  //****************  Identity vs. Null  **************************
  BOOST_CHECK( ( ( is_identity_mat((m_zeroes * m_ident2) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat((m_ident2 * m_zeroes) + m_ident2,std::numeric_limits<double>::epsilon()) ) ) );
  BOOST_CHECK( ( ( is_identity_mat((m_zeroes + m_zeroes) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m_ident2 - (m_zeroes + m_zeroes),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m_ident2 + (m_zeroes - m_zeroes),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(-((m_zeroes - m_zeroes) - m_ident2),std::numeric_limits<double>::epsilon()) ) ) );
  
  //****************  Identity vs. Identity  **************************
  BOOST_CHECK( ( ( is_identity_mat((m_ident2 + m_ident2) * 0.5,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(0.5 * (m_ident2 + m_ident2),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m_ident2 * m_ident2,std::numeric_limits<double>::epsilon()) ) ) );
  
  BOOST_CHECK( ( ( is_identity_mat((m_ident2 - m_ident2) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m_ident2 + (m_ident2 - m_ident2),std::numeric_limits<double>::epsilon()) ) ) );
  
  //****************  Identity vs. Null  **************************
  BOOST_CHECK( ( ( is_identity_mat((m_zeroes * m_ident2) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat((m_ident2 * m_zeroes) + m_ident2,std::numeric_limits<double>::epsilon()) ) ) );
  BOOST_CHECK( ( ( is_identity_mat((m_zeroes + m_zeroes) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m_ident2 - (m_zeroes + m_zeroes),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(m_ident2 + (m_zeroes - m_zeroes),std::numeric_limits<double>::epsilon()) ) &&
                 ( is_identity_mat(-((m_zeroes - m_zeroes) - m_ident2),std::numeric_limits<double>::epsilon()) ) ) );
  
  mat<double,mat_structure::rectangular> m4321_orig(4.0,2.0,3.0,1.0);
  double f4321[] = {4.0,3.0,2.0,1.0};
  std::vector<double> v4321(f4321,f4321 + 4);
  mat<double> m4321(v4321,2,2);
  BOOST_CHECK( ( is_null_mat(m4321 - m4321_orig, std::numeric_limits<double>::epsilon()) ) );
  
  mat<double> m4321_cpy(m4321);
  BOOST_CHECK( ( is_null_mat(m4321_cpy - m4321_orig, std::numeric_limits<double>::epsilon()) ) );
  
  m4321_cpy = m4321;
  BOOST_CHECK( ( is_null_mat(m4321_cpy - m4321_orig, std::numeric_limits<double>::epsilon()) ) );
  
  m4321_cpy += mat<double>(2,2,true);
  BOOST_CHECK( ( is_identity_mat(m4321_cpy - m4321_orig, std::numeric_limits<double>::epsilon()) ) );
  
  m4321_cpy += mat_identity<double>::type(2);
  BOOST_CHECK( ( is_identity_mat(m4321_cpy - m4321_orig - mat<double>(2,2,true), std::numeric_limits<double>::epsilon()) ) );
  
  m4321_cpy -= mat<double>(2,2,true);
  BOOST_CHECK( ( is_identity_mat(m4321_cpy - m4321_orig, std::numeric_limits<double>::epsilon()) ) );
  
  m4321_cpy -= mat_identity<double>::type(2);
  BOOST_CHECK( ( is_null_mat(m4321_cpy - m4321_orig, std::numeric_limits<double>::epsilon()) ) );
  
  m4321_cpy *= 2.0;
  BOOST_CHECK( ( is_null_mat(m4321_cpy - m4321_orig - m4321_orig, std::numeric_limits<double>::epsilon()) ) );
  
  m4321_cpy = m4321_cpy * double(0.5);
  BOOST_CHECK( ( is_null_mat(m4321_cpy - m4321_orig, std::numeric_limits<double>::epsilon()) ) );
  
  m4321_cpy = m4321_cpy + mat<double>(2,2,true);
  BOOST_CHECK( ( is_identity_mat(m4321_cpy - m4321_orig, std::numeric_limits<double>::epsilon()) ) );
  
  m4321_cpy = m4321_cpy + mat_identity<double>::type(2);
  BOOST_CHECK( ( is_identity_mat(m4321_cpy - m4321_orig - mat<double>(2,2,true), std::numeric_limits<double>::epsilon()) ) );
  
  m4321_cpy = m4321_cpy - mat<double>(2,2,true);
  BOOST_CHECK( ( is_identity_mat(m4321_cpy - m4321_orig, std::numeric_limits<double>::epsilon()) ) );
  
  m4321_cpy = m4321_cpy - mat_identity<double>::type(2);
  BOOST_CHECK( ( is_null_mat(m4321_cpy - m4321_orig, std::numeric_limits<double>::epsilon()) ) );
  BOOST_CHECK( ( is_null_mat((-m4321_cpy) + m4321_orig, std::numeric_limits<double>::epsilon()) ) );
  
  m4321_cpy = m4321_cpy * mat_identity<double>::type(2);
  BOOST_CHECK( ( is_null_mat(m4321_cpy - m4321_orig, std::numeric_limits<double>::epsilon()) ) );
  
  BOOST_CHECK( (m4321_cpy == m4321) );
  BOOST_CHECK( (m4321_cpy != mat<double>(2,2,true)) );
  BOOST_CHECK( (m4321_cpy != mat_identity<double>::type(2)) );
  
  BOOST_CHECK( ( is_null_mat(m4321_cpy - m4321_orig, std::numeric_limits<double>::epsilon()) ) );
  vect_n<double> v85 = m4321_cpy * vect_n<double>(2,1.0);
  BOOST_CHECK( ( ( fabs(v85[0] - 6.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(v85[1] - 4.0) < std::numeric_limits<double>::epsilon() ) ) );
  vect_n<double> v104 =  vect_n<double>(2,1.0) * m4321_cpy;
  BOOST_CHECK( ( ( fabs(v104[0] - 7.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(v104[1] - 3.0) < std::numeric_limits<double>::epsilon() ) ) );
  
  mat<double,mat_structure::rectangular> m_ident3(3,3,true);
  set_block(m_ident3,m4321_cpy,1,1);
  BOOST_CHECK( ( ( fabs(m_ident3(0,0) - 1.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m_ident3(0,1) - 0.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m_ident3(0,2) - 0.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m_ident3(1,0) - 0.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m_ident3(1,1) - 4.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m_ident3(1,2) - 2.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m_ident3(2,0) - 0.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m_ident3(2,1) - 3.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m_ident3(2,2) - 1.0) < std::numeric_limits<double>::epsilon() ) ) );
  
  m4321_cpy = get_block(m_ident3,1,1,2,2);
  BOOST_CHECK( ( ( fabs(m4321_cpy(0,0) - 4.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_cpy(0,1) - 2.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_cpy(1,0) - 3.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_cpy(1,1) - 1.0) < std::numeric_limits<double>::epsilon() ) ) );
  mat<double,mat_structure::rectangular> m4321_trans( (transpose(m4321_cpy)) );
  BOOST_CHECK( ( ( fabs(m4321_trans(0,0) - 4.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_trans(0,1) - 3.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_trans(1,0) - 2.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_trans(1,1) - 1.0) < std::numeric_limits<double>::epsilon() ) ) );
  mat<double,mat_structure::symmetric> m4321_sym = mat<double,mat_structure::symmetric>(m4321_cpy);
  BOOST_CHECK( ( ( fabs(m4321_sym(0,0) - 4.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_sym(0,1) - 2.5) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_sym(1,0) - 2.5) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_sym(1,1) - 1.0) < std::numeric_limits<double>::epsilon() ) ) );
  mat<double,mat_structure::skew_symmetric> m4321_skew = mat<double,mat_structure::skew_symmetric>(m4321_cpy);
  const mat<double,mat_structure::skew_symmetric>& m4321_skew_ref = m4321_skew;
  BOOST_CHECK( ( ( fabs(m4321_skew_ref(0,0) - 0.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_skew_ref(0,1) + 0.5) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_skew_ref(1,0) - 0.5) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_skew_ref(1,1) - 0.0) < std::numeric_limits<double>::epsilon() ) ) );
  mat<double,mat_structure::rectangular> m4321_twice( (double(2.0) * m4321_cpy) );
  BOOST_CHECK( ( ( fabs(m4321_twice(0,0) - 8.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_twice(0,1) - 4.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_twice(1,0) - 6.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m4321_twice(1,1) - 2.0) < std::numeric_limits<double>::epsilon() ) ) );
  
  mat<double,mat_structure::diagonal> m123(vect_n<double>(1.0,2.0,3.0));
  const mat<double,mat_structure::diagonal>& m123_ref = m123;
  BOOST_CHECK( ( ( fabs(m123_ref(0,0) - 1.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_ref(1,1) - 2.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_ref(2,2) - 3.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_ref(0,1)) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_ref(0,2)) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_ref(1,0)) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_ref(1,2)) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_ref(2,0)) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_ref(2,1)) < std::numeric_limits<double>::epsilon() ) ) );
  
  mat<double,mat_structure::skew_symmetric> m123_skew = mat<double,mat_structure::skew_symmetric>(vect_n<double>(1.0,2.0,3.0));
  const mat<double,mat_structure::skew_symmetric>& m123_skew_ref = m123_skew;
  BOOST_CHECK( ( ( fabs(m123_skew_ref(0,0)) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_skew_ref(1,1)) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_skew_ref(2,2)) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_skew_ref(0,1) + 3.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_skew_ref(0,2) - 2.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_skew_ref(1,0) - 3.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_skew_ref(1,2) + 1.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_skew_ref(2,0) + 2.0) < std::numeric_limits<double>::epsilon() ) &&
                 ( fabs(m123_skew_ref(2,1) - 1.0) < std::numeric_limits<double>::epsilon() ) ) );
  
};



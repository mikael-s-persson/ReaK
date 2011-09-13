
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

#include "base/defs.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>

#include "mat_alg.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979324
#define M_PI_2 1.57079632679489662
#endif

int main() {

using namespace ReaK;

  unsigned int passed = 0;

  try {

    if(true){
      RK_NOTICE(2,"/*********************************************/");
      RK_NOTICE(2,"/**** PRIMITVE VARIABLE-SIZE MATRIX TESTS ****/");
      RK_NOTICE(2,"/*********************************************/");

      mat<double,mat_structure::rectangular> m1234(1.0,2.0,3.0,4.0);
      if( ( std::fabs(m1234(0,0) - 1.0) < std::numeric_limits<double>::epsilon() ) &&
	  ( std::fabs(m1234(0,1) - 2.0) < std::numeric_limits<double>::epsilon() ) &&
	  ( std::fabs(m1234(1,0) - 3.0) < std::numeric_limits<double>::epsilon() ) &&
	  ( std::fabs(m1234(1,1) - 4.0) < std::numeric_limits<double>::epsilon() ) )
	++passed;
      else
	RK_ERROR("mat<double,mat_structure::rectangular>::mat(a11,a12,a21,a22) test did not pass!");
      { std::stringstream ss;
        ss << m1234;
	if( ss.str() == "((1; 2); (3; 4))" )
	  ++passed;
	else
	  RK_ERROR("operator<<(ostream&,mat<double,mat_structure::rectangular>) test did not pass!");
      };
      m1234 = transpose(m1234);
      if( ( std::fabs(m1234(0,0) - 1.0) < std::numeric_limits<double>::epsilon() ) &&
	  ( std::fabs(m1234(0,1) - 3.0) < std::numeric_limits<double>::epsilon() ) &&
	  ( std::fabs(m1234(1,0) - 2.0) < std::numeric_limits<double>::epsilon() ) &&
	  ( std::fabs(m1234(1,1) - 4.0) < std::numeric_limits<double>::epsilon() ) )
	++passed;
      else
	RK_ERROR("transpose(mat<double,mat_structure::rectangular>) test did not pass!");
      m1234 = transpose_move(m1234);
      if( ( std::fabs(m1234(0,0) - 1.0) < std::numeric_limits<double>::epsilon() ) &&
	  ( std::fabs(m1234(0,1) - 2.0) < std::numeric_limits<double>::epsilon() ) &&
	  ( std::fabs(m1234(1,0) - 3.0) < std::numeric_limits<double>::epsilon() ) &&
	  ( std::fabs(m1234(1,1) - 4.0) < std::numeric_limits<double>::epsilon() ) )
	++passed;
      else
	RK_ERROR("transpose_move(mat<double,mat_structure::rectangular>) test did not pass!");
      mat<double,mat_structure::rectangular> m_ident(2,2,true);
      if( ( std::fabs(m_ident(0,0) - 1.0) < std::numeric_limits<double>::epsilon() ) &&
	  ( std::fabs(m_ident(0,1)) < std::numeric_limits<double>::epsilon() ) &&
	  ( std::fabs(m_ident(1,0)) < std::numeric_limits<double>::epsilon() ) &&
	  ( std::fabs(m_ident(1,1) - 1.0) < std::numeric_limits<double>::epsilon() ) )
	++passed;
      else
	RK_ERROR("mat<double,mat_structure::rectangular>::mat(M,N,bool) test did not pass!");
      if( is_identity_mat(m_ident,std::numeric_limits<double>::epsilon()) )
	++passed;
      else
	RK_ERROR("is_identity_mat(M,tol) test did not pass!");
      mat_identity<double>::type m_ident2(2);
      if( is_identity_mat(m_ident2,std::numeric_limits<double>::epsilon()) )
	++passed;
      else
	RK_ERROR("mat<double,identity>(N) test did not pass!");
      
      mat_null<double>::type m_zeroes(2,2);
      if( is_null_mat(m_zeroes,std::numeric_limits<double>::epsilon()) )
	++passed;
      else
	RK_ERROR("is_null_mat(M,tol) test did not pass!");
      if( is_null_mat(m_zeroes * m_ident2,std::numeric_limits<double>::epsilon()) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,null>,mat<T,identity>) test did not pass!");
      if( is_null_mat(m_ident2 * m_zeroes,std::numeric_limits<double>::epsilon()) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,identity>,mat<T,null>) test did not pass!");
      if( is_null_mat(m_zeroes * m1234,std::numeric_limits<double>::epsilon()) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,null>,mat<T,rect>) test did not pass!");
      if( is_null_mat(m1234 * m_zeroes,std::numeric_limits<double>::epsilon()) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,rect>,mat<T,null>) test did not pass!");
      
      mat<double,mat_structure::rectangular> m1234_inv(-2.0, 1.0, 1.5, -0.5);
      
      //****************  Rectangular vs. Rectangular  **************************
      if( is_identity_mat(m1234_inv * m1234,std::numeric_limits<double>::epsilon()) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,rect>,mat<T,rect>) test did not pass!");
      if( ( is_identity_mat(0.5 * (m1234_inv * m1234 + m1234_inv * m1234),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat((m1234_inv * m1234 + m1234_inv * m1234) * 0.5,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat((m1234_inv * m1234 - m1234_inv * m1234) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_null_mat(m1234_inv * m1234 - m1234_inv * m1234,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_null_mat(m1234_inv * m1234 + (-(m1234_inv * m1234)),std::numeric_limits<double>::epsilon()) ) ) 
	++passed;
      else
	RK_ERROR("operator+-(mat<T,rect>,mat<T,rect>) and operator*(Scalar,mat<T,rect>) test did not pass!");
      
      //****************  Square vs. Square  **************************
      mat<double,mat_structure::square> m1234_sqr(1.0,2.0,3.0,4.0);
      mat<double,mat_structure::square> m1234_inv_sqr(-2.0,1.0,1.5,-0.5);
      if( is_identity_mat(m1234_inv_sqr * m1234_sqr,std::numeric_limits<double>::epsilon()) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,square>,mat<T,square>) test did not pass!");
      if( ( is_identity_mat(0.5 * (m1234_inv_sqr * m1234_sqr + m1234_inv_sqr * m1234_sqr),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat((m1234_inv_sqr * m1234_sqr + m1234_inv_sqr * m1234_sqr) * 0.5,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat((m1234_inv_sqr * m1234_sqr - m1234_inv_sqr * m1234_sqr) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_null_mat(m1234_inv_sqr * m1234_sqr - m1234_inv_sqr * m1234_sqr,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_null_mat(m1234_inv_sqr * m1234_sqr + (-(m1234_inv_sqr * m1234_sqr)),std::numeric_limits<double>::epsilon()) ) ) 
	++passed;
      else
	RK_ERROR("operator+-(mat<T,square>,mat<T,square>) and operator*(Scalar,mat<T,square>) test did not pass!");
      
      //****************  Square vs. Rectangular  **************************
      if( ( is_identity_mat(m1234_inv * ((m1234_inv_sqr * m1234_sqr) * m1234),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m1234_inv * (m1234 * (m1234_inv_sqr * m1234_sqr)),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,square>,mat<T,rect>) test did not pass!");
      if( ( is_identity_mat(0.5 * (m1234 * m1234_inv + (m1234_inv_sqr * m1234_sqr)),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(0.5 * ((m1234_inv_sqr * m1234_sqr) + m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat((m1234_inv_sqr * m1234_sqr) + (m1234 * m1234_inv - (m1234_inv_sqr * m1234_sqr)),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat((m1234_inv_sqr * m1234_sqr) + ((m1234_inv_sqr * m1234_sqr) - m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator-+(mat<T,rect>,mat<T,square>) test did not pass!");
      
      //****************  Symmetric vs. Symmetric  **************************
      mat<double,mat_structure::symmetric> m123_sym(1.0,2.0,3.0);
      mat<double,mat_structure::symmetric> m123_inv_sym(-3.0,2.0,-1.0);
      if( is_identity_mat(m123_inv_sym * m123_sym,std::numeric_limits<double>::epsilon()) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,square>,mat<T,square>) test did not pass!");
      if( ( is_identity_mat(0.5 * (m123_inv_sym * m123_sym + m123_inv_sym * m123_sym),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat((m123_inv_sym * m123_sym + m123_inv_sym * m123_sym) * 0.5,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat((m123_inv_sym * m123_sym - m123_inv_sym * m123_sym) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_null_mat(m123_inv_sym * m123_sym - m123_inv_sym * m123_sym,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_null_mat(m123_inv_sym * m123_sym + (-(m123_inv_sym * m123_sym)),std::numeric_limits<double>::epsilon()) ) ) 
	++passed;
      else
	RK_ERROR("operator+-(mat<T,square>,mat<T,square>) and operator*(Scalar,mat<T,square>) test did not pass!");
      
      //****************  Symmetric vs. Rectangular  **************************
      if( ( is_identity_mat(m1234_inv * ((m123_sym * m123_inv_sym) * m1234),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m1234_inv * (m1234 * (m123_sym * m123_inv_sym)),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,symmetric>,mat<T,rect>) test did not pass!");
      if( ( is_identity_mat(0.5 * (m1234 * m1234_inv + (m123_sym * m123_inv_sym)),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(0.5 * ((m123_sym * m123_inv_sym) + m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat((m123_sym * m123_inv_sym) + (m1234 * m1234_inv - (m123_sym * m123_inv_sym)),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat((m123_sym * m123_inv_sym) + ((m123_sym * m123_inv_sym) - m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator-+(mat<T,rect>,mat<T,symmetric>) test did not pass!");
      
      //****************  Symmetric vs. Square  **************************
      if( ( is_identity_mat(m1234_inv_sqr * ((m123_sym * m123_inv_sym) * m1234_sqr),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m1234_inv_sqr * (m1234_sqr * (m123_sym * m123_inv_sym)),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,symmetric>,mat<T,square>) test did not pass!");
      if( ( is_identity_mat(0.5 * (m1234_sqr * m1234_inv_sqr + (m123_sym * m123_inv_sym)),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(0.5 * ((m123_sym * m123_inv_sym) + m1234_sqr * m1234_inv_sqr),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat((m123_sym * m123_inv_sym) + (m1234_sqr * m1234_inv_sqr - (m123_sym * m123_inv_sym)),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat((m123_sym * m123_inv_sym) + ((m123_sym * m123_inv_sym) - m1234_sqr * m1234_inv_sqr),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator-+(mat<T,square>,mat<T,symmetric>) test did not pass!");
      
      //****************  Rectangular vs. Diagonal  **************************
      mat<double,mat_structure::diagonal> m_ident_diag(2,1.0);
      if( ( is_identity_mat(m1234_inv * (m_ident_diag * m1234),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m1234_inv * (m1234 * m_ident_diag),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,diagonal>,mat<T,rect>) test did not pass!");
      if( ( is_identity_mat(0.5 * (m1234 * m1234_inv + m_ident_diag),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(0.5 * (m_ident_diag + m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m_ident_diag + (m1234 * m1234_inv - m_ident_diag),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m_ident_diag + (m_ident_diag - m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator-+(mat<T,rect>,mat<T,diagonal>) test did not pass!");
      
      //****************  Rectangular vs. Scalar  **************************
      mat<double,mat_structure::scalar> m_ident_scalar(2,1.0);
      if( ( is_identity_mat(m1234_inv * (m_ident_scalar * m1234),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m1234_inv * (m1234 * m_ident_scalar),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,scalar>,mat<T,rect>) test did not pass!");
      if( ( is_identity_mat(0.5 * (m1234 * m1234_inv + m_ident_scalar),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(0.5 * (m_ident_scalar + m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m_ident_scalar + (m1234 * m1234_inv - m_ident_scalar),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m_ident_scalar + (m_ident_scalar - m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator-+(mat<T,rect>,mat<T,scalar>) test did not pass!");
      
      //****************  Rectangular vs. Identity  **************************
      if( ( is_identity_mat(m1234_inv * (m_ident2 * m1234),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m1234_inv * (m1234 * m_ident2),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,identity>,mat<T,rect>) test did not pass!");
      if( ( is_identity_mat(0.5 * (m1234 * m1234_inv + m_ident2),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(0.5 * (m_ident2 + m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m_ident2 + (m1234 * m1234_inv - m_ident2),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m_ident2 + (m_ident2 - m1234 * m1234_inv),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator-+(mat<T,rect>,mat<T,identity>) test did not pass!");
      
      //****************  Rectangular vs. Null  **************************
      if( ( is_null_mat(m1234 * m_zeroes,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_null_mat(m_zeroes * m1234,std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,null>,mat<T,rect>) test did not pass!");
      if( ( is_identity_mat(m_zeroes - m1234_inv * (-m1234),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m1234_inv * m1234 - m_zeroes,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m1234_inv * m1234 + m_zeroes,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m_zeroes + m1234_inv * m1234,std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator-+(mat<T,null>,mat<T,rect>) test did not pass!");
      
      //****************  Identity vs. Null  **************************
      if( ( is_identity_mat((m_zeroes * m_ident2) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat((m_ident2 * m_zeroes) + m_ident2,std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,null>,mat<T,identity>) test did not pass!");
      if( ( is_identity_mat((m_zeroes + m_zeroes) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m_ident2 - (m_zeroes + m_zeroes),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m_ident2 + (m_zeroes - m_zeroes),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(-((m_zeroes - m_zeroes) - m_ident2),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator-+(mat<T,null>,mat<T,identity>) test did not pass!");
      
      //****************  Identity vs. Identity  **************************
      if( ( is_identity_mat((m_ident2 + m_ident2) * 0.5,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(0.5 * (m_ident2 + m_ident2),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m_ident2 * m_ident2,std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator+(mat<T,identity>,mat<T,identity>) and operator*(mat<T,identity>,mat<T,identity>) test did not pass!");
      if( ( is_identity_mat((m_ident2 - m_ident2) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m_ident2 + (m_ident2 - m_ident2),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator-(mat<T,identity>,mat<T,identity>) test did not pass!");
      
      //****************  Identity vs. Null  **************************
      if( ( is_identity_mat((m_zeroes * m_ident2) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat((m_ident2 * m_zeroes) + m_ident2,std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator*(mat<T,null>,mat<T,identity>) test did not pass!");
      if( ( is_identity_mat((m_zeroes + m_zeroes) + m_ident2,std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m_ident2 - (m_zeroes + m_zeroes),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(m_ident2 + (m_zeroes - m_zeroes),std::numeric_limits<double>::epsilon()) ) &&
	  ( is_identity_mat(-((m_zeroes - m_zeroes) - m_ident2),std::numeric_limits<double>::epsilon()) ) )
	++passed;
      else
	RK_ERROR("operator-+(mat<T,null>,mat<T,identity>) test did not pass!");
      
      
#if 1
      double f4321[] = {4.0,3.0,2.0,1.0};
      std::vector<double> v4321(f4321,f4321 + 4);
      mat<double> m4321(v4321,2,2);
      RK_NOTICE(2,"m4321 = " << m4321);

      mat<double> m4321_cpy(m4321);
      RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

      m4321_cpy = m4321;
      RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

      m4321_cpy += mat<double>(2,2,true);
      RK_NOTICE(2,"m4321 copy + identity = " << m4321_cpy);

      m4321_cpy += mat_identity<double>::type(2);
      RK_NOTICE(2,"m4321 copy + identity = " << m4321_cpy);

      m4321_cpy -= mat<double>(2,2,true);
      RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

      m4321_cpy -= mat_identity<double>::type(2);
      RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

      m4321_cpy *= 2.0;
      RK_NOTICE(2,"m4321 copy twice = " << m4321_cpy);

      m4321_cpy = m4321_cpy * double(0.5);
      RK_NOTICE(2,"m4321 copy = " << m4321_cpy);

      m4321_cpy = m4321_cpy + mat<double>(2,2,true);
      RK_NOTICE(2,"m4321 copy + identity = " << m4321_cpy);

      m4321_cpy = m4321_cpy + mat_identity<double>::type(2);
      RK_NOTICE(2,"m4321 copy + identity = " << m4321_cpy);

      m4321_cpy = m4321_cpy - mat<double>(2,2,true);
      RK_NOTICE(2,"m4321 copy = " << m4321_cpy);
      m4321_cpy = m4321_cpy - mat_identity<double>::type(2);
      RK_NOTICE(2,"m4321 copy = " << m4321_cpy);
      RK_NOTICE(2,"m4321 copy negated = " << -m4321_cpy);

      m4321_cpy = m4321_cpy * mat_identity<double>::type(2);
      RK_NOTICE(2,"m4321 copy * identity = " << m4321_cpy);

      if(m4321_cpy == m4321)
        RK_NOTICE(2,"m4321 copy is equal to m4321");
      if(m4321_cpy != mat<double>(2,2,true))
        RK_NOTICE(2,"m4321 copy is not equal to identity");
      if(m4321_cpy != mat_identity<double>::type(2))
        RK_NOTICE(2,"m4321 copy is not equal to identity");

      RK_NOTICE(2,"m4321 copy = " << m4321_cpy);
      vect_n<double> v85 = m4321_cpy * vect_n<double>(2,1.0);
      RK_NOTICE(2,"m4321 copy * vect(1.0,1.0) = " << v85);
      vect_n<double> v104 =  vect_n<double>(2,1.0) * m4321_cpy;
      RK_NOTICE(2,"vect(1.0,1.0) * m4321 copy = " << v104);

      mat<double,mat_structure::rectangular> m_ident3(3,3,true);
      set_block(m_ident3,m4321_cpy,1,1);
      RK_NOTICE(2,"identity 3x3 with sub matrix m4321 is " << m_ident3);

      m4321_cpy = get_block(m_ident3,1,1,2,2);
      RK_NOTICE(2,"m4321 copy = " << m4321_cpy);
      RK_NOTICE(2,"m4321 copy transpose = " << transpose(m4321_cpy));
      RK_NOTICE(2,"m4321 copy sym part = " << (mat<double,mat_structure::symmetric>(m4321_cpy)));
      RK_NOTICE(2,"m4321 copy skew part = " << (mat<double,mat_structure::skew_symmetric>(m4321_cpy)));
      RK_NOTICE(2,"m4321 copy twice = " << double(2.0) * m4321_cpy);

      mat<double,mat_structure::diagonal> m123(vect_n<double>(1.0,2.0,3.0));
      RK_NOTICE(2,"m123 diagonal = " << m123);
      //mat_cm<double> m_blk_diag(block_diag_mat(mat_cm<double>(1.0,2.0,3.0,4.0),mat_cm<double>(1,1,true)));
      //RK_NOTICE(2,"m_blk_diag block diagonal = " << m_blk_diag);
      //mat_cm<double> m_blk(block_mat(mat_cm<double>(1.0,2.0,3.0,4.0),m4321.getSubMat(0,0,2,1),m4321.getSubMat(1,0,1,2),mat_cm<double>(1,1,true)));
      //RK_NOTICE(2,"m_blk block matrix = " << m_blk);
      RK_NOTICE(2,"vect(1,2,3) skew matrix = " << (mat<double,mat_structure::skew_symmetric>(vect_n<double>(1.0,2.0,3.0))));
#endif
      RK_NOTICE(2,"/!!!!!! CONGRATULATIONS! SECTION PASSED !!!!!!/");
    };
  
  } catch(std::exception& e) {
    RK_ERROR("An exception has occurred during the math_gen test: '" << e.what() << "'");
  } catch(...) {
    RK_ERROR("An unexpected and unidentified exception has occurred during the math_gen test.");
  };
  
  RK_NOTICE(2,"There were " << passed << " successful tests passed on the lin_alg::mat_alg library, out of 45 possible successes.");
  
  return 0;
};



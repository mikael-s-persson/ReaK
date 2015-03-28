
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
#include <ReaK/math/lin_alg/mat_norms.hpp>
#include <ReaK/math/kinetostatics/rotations.hpp>

#include <iostream>
#include <fstream>
#include <cstdio>
#include <sstream>


#ifndef M_PI
#define M_PI 3.14159265358979324
#define M_PI_2 1.57079632679489662
#endif


#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE rotations
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


BOOST_AUTO_TEST_CASE( rotations_2D_tests ) {
  using namespace ReaK;
  const double rel_tol = 10.0 * std::numeric_limits< double >::epsilon();

  rot_mat_2D< double > r_ident;
  std::stringstream ss_r_ident;
  ss_r_ident << r_ident;
  BOOST_CHECK( ss_r_ident.str() == "(angle = 0)" );
  BOOST_CHECK_SMALL( r_ident.getAngle(), rel_tol );
  BOOST_CHECK_CLOSE( r_ident( 0, 0 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( r_ident( 0, 1 ), rel_tol );
  BOOST_CHECK_SMALL( r_ident( 1, 0 ), rel_tol );
  BOOST_CHECK_CLOSE( r_ident( 1, 1 ), 1.0, rel_tol );
  rot_mat_2D< double > r_45deg( 0.25f * double( M_PI ) );
  BOOST_CHECK_CLOSE( r_45deg.getAngle(), 0.25f * double( M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( r_45deg( 0, 0 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( r_45deg( 0, 1 ), -std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( r_45deg( 1, 0 ), std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( r_45deg( 1, 1 ), std::sqrt( 0.5 ), rel_tol );
  rot_mat_2D< double > r_45deg_cpy( r_45deg );
  BOOST_CHECK_CLOSE( r_45deg_cpy.getAngle(), 0.25f * double( M_PI ), rel_tol );
  mat< double, mat_structure::square > m_45deg = r_45deg.getMat();
  BOOST_CHECK_CLOSE( m_45deg( 0, 0 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( m_45deg( 0, 1 ), -std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( m_45deg( 1, 0 ), std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( m_45deg( 1, 1 ), std::sqrt( 0.5 ), rel_tol );
  r_45deg_cpy.setAngle( 0.0 );
  BOOST_CHECK_SMALL( r_45deg_cpy.getAngle(), rel_tol );
  r_45deg_cpy = r_45deg;
  rot_mat_2D< double > r_90deg = r_45deg * r_45deg_cpy;
  BOOST_CHECK_CLOSE( r_90deg.getAngle(), 0.5f * double( M_PI ), 10.0 * rel_tol );
  r_90deg *= r_45deg;
  BOOST_CHECK_CLOSE( r_90deg.getAngle(), 0.75f * double( M_PI ), rel_tol );
  vect< double, 2 > v1 = r_45deg * vect< double, 2 >( 1.0, 1.0 );
  BOOST_CHECK_SMALL( v1[0], rel_tol );
  BOOST_CHECK_CLOSE( v1[1], std::sqrt( 2.0 ), 10.0 * rel_tol );
  v1 = vect< double, 2 >( 1.0, 1.0 ) * r_45deg;
  BOOST_CHECK_SMALL( v1[1], rel_tol );
  BOOST_CHECK_CLOSE( v1[0], std::sqrt( 2.0 ), 10.0 * rel_tol );
  r_90deg *= invert( r_45deg );
  BOOST_CHECK_CLOSE( r_90deg.getAngle(), 0.5f * double( M_PI ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( trace( r_45deg_cpy ), std::sqrt( 2.0 ), rel_tol );
  BOOST_CHECK_CLOSE( determinant( r_45deg_cpy ), 1.0, rel_tol );
  const mat< double, mat_structure::symmetric > msym_45deg = r_45deg_cpy.getSymPart();
  BOOST_CHECK_CLOSE( msym_45deg( 0, 0 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_SMALL( msym_45deg( 0, 1 ), rel_tol );
  BOOST_CHECK_SMALL( msym_45deg( 1, 0 ), rel_tol );
  BOOST_CHECK_CLOSE( msym_45deg( 1, 1 ), std::sqrt( 0.5 ), rel_tol );
  const mat< double, mat_structure::skew_symmetric > mskw_45deg = r_45deg_cpy.getSkewSymPart();
  BOOST_CHECK_SMALL( mskw_45deg( 0, 0 ), rel_tol );
  BOOST_CHECK_CLOSE( mskw_45deg( 0, 1 ), -std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( mskw_45deg( 1, 0 ), std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_SMALL( mskw_45deg( 1, 1 ), rel_tol );
  BOOST_CHECK_SMALL( ( r_45deg * invert( r_45deg ) ).getAngle(), rel_tol );

  BOOST_CHECK( r_ident == r_ident );
  BOOST_CHECK( r_45deg != r_ident );

  mat< double, mat_structure::rectangular > m23( 2, 3 );
  m23( 0, 0 ) = 1.0;
  m23( 1, 0 ) = 0.0;
  m23( 0, 1 ) = 0.0;
  m23( 1, 1 ) = 1.0;
  m23( 0, 2 ) = 1.0;
  m23( 1, 2 ) = 1.0;
  mat< double, mat_structure::rectangular > m23_r45( r_45deg * m23 );
  BOOST_CHECK_CLOSE( m23_r45( 0, 0 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( m23_r45( 0, 1 ), -std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( m23_r45( 1, 0 ), std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( m23_r45( 1, 1 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_SMALL( m23_r45( 0, 2 ), rel_tol );
  BOOST_CHECK_CLOSE( m23_r45( 1, 2 ), std::sqrt( 2.0 ), 10.0 * rel_tol );

  mat< double, mat_structure::rectangular > m32( 3, 2 );
  m32( 0, 0 ) = 1.0;
  m32( 1, 0 ) = 0.0;
  m32( 2, 0 ) = 1.0;
  m32( 0, 1 ) = 0.0;
  m32( 1, 1 ) = 1.0;
  m32( 2, 1 ) = 1.0;
  mat< double, mat_structure::rectangular > m32_r45( m32 * r_45deg );
  BOOST_CHECK_CLOSE( m32_r45( 0, 0 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( m32_r45( 0, 1 ), -std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( m32_r45( 1, 0 ), std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( m32_r45( 1, 1 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_SMALL( m32_r45( 2, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( m32_r45( 2, 0 ), std::sqrt( 2.0 ), 10.0 * rel_tol );

  trans_mat_2D< double > t_ident;
  std::stringstream ss_t_ident;
  ss_t_ident << t_ident;
  BOOST_CHECK( ( ss_t_ident.str() == "(angle = 0; translation = (0; 0))" ) );
  BOOST_CHECK_CLOSE( t_ident( 0, 0 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( t_ident( 0, 1 ), rel_tol );
  BOOST_CHECK_SMALL( t_ident( 0, 2 ), rel_tol );
  BOOST_CHECK_SMALL( t_ident( 1, 0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 1, 1 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( t_ident( 1, 2 ), rel_tol );
  BOOST_CHECK_SMALL( t_ident( 2, 0 ), rel_tol );
  BOOST_CHECK_SMALL( t_ident( 2, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 2, 2 ), 1.0, rel_tol );

  trans_mat_2D< double > t_45deg_11( 0.25f * double( M_PI ), vect< double, 2 >( 1.0, 1.0 ) );
  BOOST_CHECK_CLOSE( t_45deg_11( 0, 0 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45deg_11( 0, 1 ), -std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( t_45deg_11( 0, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45deg_11( 1, 0 ), std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( t_45deg_11( 1, 1 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45deg_11( 1, 2 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( t_45deg_11( 2, 0 ), rel_tol );
  BOOST_CHECK_SMALL( t_45deg_11( 2, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45deg_11( 2, 2 ), 1.0, rel_tol );

  mat< double, mat_structure::square > m_45deg_11 = t_45deg_11.getMat();
  BOOST_CHECK_CLOSE( m_45deg_11( 0, 0 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( m_45deg_11( 0, 1 ), -std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( m_45deg_11( 0, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( m_45deg_11( 1, 0 ), std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( m_45deg_11( 1, 1 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( m_45deg_11( 1, 2 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( m_45deg_11( 2, 0 ), rel_tol );
  BOOST_CHECK_SMALL( m_45deg_11( 2, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( m_45deg_11( 2, 2 ), 1.0, rel_tol );

  trans_mat_2D< double > t_45deg_11_cpy( t_45deg_11 );
  BOOST_CHECK_CLOSE( t_45deg_11_cpy.getAngle(), 0.25 * double( M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45deg_11_cpy.getTranslation()[0], 1.0, rel_tol );
  t_45deg_11_cpy.setAngle( 0.0 );
  BOOST_CHECK_SMALL( t_45deg_11_cpy.getAngle(), rel_tol );
  t_45deg_11_cpy.setRotMat( r_45deg );
  BOOST_CHECK_CLOSE( t_45deg_11_cpy.getAngle(), 0.25 * double( M_PI ), rel_tol );
  t_45deg_11_cpy.setTranslation( vect< double, 2 >( -1.0, -1.0 ) );
  BOOST_CHECK_CLOSE( t_45deg_11_cpy.getTranslation()[0], -1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45deg_11_cpy.getTranslation()[1], -1.0, rel_tol );
  t_45deg_11_cpy.setTranslation( vect< double, 2 >( 1.0, 1.0 ) );
  BOOST_CHECK_CLOSE( trace( t_45deg_11_cpy ), 1.0 + std::sqrt( 2.0 ), rel_tol );
  BOOST_CHECK_CLOSE( determinant( t_45deg_11_cpy ), 1.0, rel_tol );
  const mat< double, mat_structure::symmetric > msym_45deg_11_cpy( t_45deg_11_cpy.getSymPart() );
  BOOST_CHECK_CLOSE( msym_45deg_11_cpy( 0, 0 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_SMALL( msym_45deg_11_cpy( 0, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( msym_45deg_11_cpy( 0, 2 ), 0.5, rel_tol );
  BOOST_CHECK_SMALL( msym_45deg_11_cpy( 1, 0 ), rel_tol );
  BOOST_CHECK_CLOSE( msym_45deg_11_cpy( 1, 1 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( msym_45deg_11_cpy( 1, 2 ), 0.5, rel_tol );
  BOOST_CHECK_CLOSE( msym_45deg_11_cpy( 2, 0 ), 0.5, rel_tol );
  BOOST_CHECK_CLOSE( msym_45deg_11_cpy( 2, 1 ), 0.5, rel_tol );
  BOOST_CHECK_CLOSE( msym_45deg_11_cpy( 2, 2 ), 1.0, rel_tol );
  const mat< double, mat_structure::skew_symmetric > mskw_45deg_11_cpy( t_45deg_11_cpy.getSkewSymPart() );
  BOOST_CHECK_SMALL( mskw_45deg_11_cpy( 0, 0 ), rel_tol );
  BOOST_CHECK_CLOSE( mskw_45deg_11_cpy( 0, 1 ), -std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( mskw_45deg_11_cpy( 0, 2 ), 0.5, rel_tol );
  BOOST_CHECK_CLOSE( mskw_45deg_11_cpy( 1, 0 ), std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_SMALL( mskw_45deg_11_cpy( 1, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( mskw_45deg_11_cpy( 1, 2 ), 0.5, rel_tol );
  BOOST_CHECK_CLOSE( mskw_45deg_11_cpy( 2, 0 ), -0.5, rel_tol );
  BOOST_CHECK_CLOSE( mskw_45deg_11_cpy( 2, 1 ), -0.5, rel_tol );
  BOOST_CHECK_SMALL( mskw_45deg_11_cpy( 2, 2 ), rel_tol );
  BOOST_CHECK_CLOSE( invert( t_45deg_11_cpy ).getAngle(), -0.25 * double( M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( invert( t_45deg_11_cpy ).getTranslation()[0], -std::sqrt( 2.0 ), 10.0 * rel_tol );
  mat< double, mat_structure::square > mt_45deg_11_cpy( transpose( t_45deg_11_cpy ) );
  BOOST_CHECK_CLOSE( mt_45deg_11_cpy( 0, 0 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( mt_45deg_11_cpy( 0, 1 ), std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_SMALL( mt_45deg_11_cpy( 0, 2 ), rel_tol );
  BOOST_CHECK_CLOSE( mt_45deg_11_cpy( 1, 0 ), -std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( mt_45deg_11_cpy( 1, 1 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_SMALL( mt_45deg_11_cpy( 1, 2 ), rel_tol );
  BOOST_CHECK_CLOSE( mt_45deg_11_cpy( 2, 0 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( mt_45deg_11_cpy( 2, 1 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( mt_45deg_11_cpy( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( ( invert( t_45deg_11_cpy ) * t_45deg_11_cpy ).getAngle(), rel_tol );
  BOOST_CHECK_SMALL( ( invert( t_45deg_11_cpy ) * t_45deg_11_cpy ).getTranslation()[0], rel_tol );
  t_ident *= t_45deg_11_cpy;
  BOOST_CHECK_CLOSE( t_ident.getAngle(), 0.25 * double( M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_ident.getTranslation()[0], 1.0, rel_tol );
  t_ident *= invert( t_45deg_11_cpy );
  BOOST_CHECK_SMALL( t_ident.getAngle(), rel_tol );
  BOOST_CHECK_SMALL( t_ident.getTranslation()[0], rel_tol );

  vect< double, 2 > v2 = t_45deg_11 * vect< double, 2 >( 1.0, 1.0 );
  BOOST_CHECK_CLOSE( v2[0], 1.0, rel_tol );
  BOOST_CHECK_CLOSE( v2[1], 1.0 + std::sqrt( 2.0 ), rel_tol );

  vect< double, 3 > v3 = t_45deg_11 * vect< double, 3 >( 1.0, 1.0, 1.0 );
  BOOST_CHECK_CLOSE( v3[0], 1.0, rel_tol );
  BOOST_CHECK_CLOSE( v3[1], 1.0 + std::sqrt( 2.0 ), rel_tol );
  BOOST_CHECK_CLOSE( v3[2], 1.0, rel_tol );
  v2 = t_45deg_11.rotate( vect< double, 2 >( 1.0, 1.0 ) );
  BOOST_CHECK_SMALL( v2[0], rel_tol );
  BOOST_CHECK_CLOSE( v2[1], std::sqrt( 2.0 ), 10.0 * rel_tol );

  mat< double, mat_structure::rectangular > m23_t( 2, 3 );
  m23_t( 0, 0 ) = 1.0;
  m23_t( 1, 0 ) = 0.0;
  m23_t( 0, 1 ) = 0.0;
  m23_t( 1, 1 ) = 1.0;
  m23_t( 0, 2 ) = 1.0;
  m23_t( 1, 2 ) = 1.0;
  mat< double, mat_structure::rectangular > m23_r45_11( m23_t * t_45deg_11 );
  BOOST_CHECK_CLOSE( m23_r45_11( 0, 0 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( m23_r45_11( 0, 1 ), -std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( m23_r45_11( 0, 2 ), 2.0, rel_tol );
  BOOST_CHECK_CLOSE( m23_r45_11( 1, 0 ), std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( m23_r45_11( 1, 1 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( m23_r45_11( 1, 2 ), 2.0, rel_tol );

  mat< double, mat_structure::rectangular > m32_t( 3, 2 );
  m32_t( 0, 0 ) = 1.0;
  m32_t( 1, 0 ) = 0.0;
  m32_t( 2, 0 ) = 1.0;
  m32_t( 0, 1 ) = 0.0;
  m32_t( 1, 1 ) = 1.0;
  m32_t( 2, 1 ) = 1.0;
  mat< double, mat_structure::rectangular > m32_r45_11( t_45deg_11 * m32_t );
  BOOST_CHECK_CLOSE( m32_r45_11( 0, 0 ), std::sqrt( 0.5 ) + 1.0, rel_tol );
  BOOST_CHECK_CLOSE( m32_r45_11( 0, 1 ), -std::sqrt( 0.5 ) + 1.0, 100.0 * rel_tol );
  BOOST_CHECK_CLOSE( m32_r45_11( 1, 0 ), std::sqrt( 0.5 ) + 1.0, rel_tol );
  BOOST_CHECK_CLOSE( m32_r45_11( 1, 1 ), std::sqrt( 0.5 ) + 1.0, rel_tol );
  BOOST_CHECK_CLOSE( m32_r45_11( 2, 0 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( m32_r45_11( 2, 1 ), 1.0, rel_tol );
};


BOOST_AUTO_TEST_CASE( rotations_3D_tests ) {
  using namespace ReaK;
  const double rel_tol = 10.0 * std::numeric_limits< double >::epsilon();

  rot_mat_3D< double > r_ident;
  BOOST_CHECK( ( is_diagonal( r_ident, rel_tol ) ) );
  BOOST_CHECK_CLOSE( elem_norm_2( r_ident ), std::sqrt( 3.0 ), rel_tol );
  double r_45z_a[] = {double( std::cos( 0.25 * M_PI ) ),  double( std::sin( 0.25 * M_PI ) ), 0.0,
                      double( -std::sin( 0.25 * M_PI ) ), double( std::cos( 0.25 * M_PI ) ), 0.0,
                      0.0,                                0.0,                               1.0};
  rot_mat_3D< double > r_45z( r_45z_a );
  BOOST_CHECK_CLOSE( r_45z( 0, 0 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( r_45z( 0, 1 ), -std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_SMALL( r_45z( 0, 2 ), rel_tol );
  BOOST_CHECK_CLOSE( r_45z( 1, 0 ), std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( r_45z( 1, 1 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_SMALL( r_45z( 1, 2 ), rel_tol );
  BOOST_CHECK_SMALL( r_45z( 2, 0 ), rel_tol );
  BOOST_CHECK_SMALL( r_45z( 2, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( r_45z( 2, 2 ), 1.0, rel_tol );
  rot_mat_3D< double > r_45z_cpy( r_45z );
  BOOST_CHECK_CLOSE( r_45z_cpy( 0, 0 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( r_45z_cpy( 0, 1 ), -std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_SMALL( r_45z_cpy( 0, 2 ), rel_tol );
  BOOST_CHECK_CLOSE( r_45z_cpy( 1, 0 ), std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( r_45z_cpy( 1, 1 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_SMALL( r_45z_cpy( 1, 2 ), rel_tol );
  BOOST_CHECK_SMALL( r_45z_cpy( 2, 0 ), rel_tol );
  BOOST_CHECK_SMALL( r_45z_cpy( 2, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( r_45z_cpy( 2, 2 ), 1.0, rel_tol );
  r_45z_cpy = r_45z;
  BOOST_CHECK_CLOSE( r_45z_cpy( 0, 0 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( r_45z_cpy( 0, 1 ), -std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_SMALL( r_45z_cpy( 0, 2 ), rel_tol );
  BOOST_CHECK_CLOSE( r_45z_cpy( 1, 0 ), std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( r_45z_cpy( 1, 1 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_SMALL( r_45z_cpy( 1, 2 ), rel_tol );
  BOOST_CHECK_SMALL( r_45z_cpy( 2, 0 ), rel_tol );
  BOOST_CHECK_SMALL( r_45z_cpy( 2, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( r_45z_cpy( 2, 2 ), 1.0, rel_tol );
  r_ident *= r_45z;
  BOOST_CHECK_CLOSE( r_ident( 0, 0 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_CLOSE( r_ident( 0, 1 ), -std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_SMALL( r_ident( 0, 2 ), rel_tol );
  BOOST_CHECK_CLOSE( r_ident( 1, 0 ), std::sqrt( 0.5 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( r_ident( 1, 1 ), std::sqrt( 0.5 ), rel_tol );
  BOOST_CHECK_SMALL( r_ident( 1, 2 ), rel_tol );
  BOOST_CHECK_SMALL( r_ident( 2, 0 ), rel_tol );
  BOOST_CHECK_SMALL( r_ident( 2, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( r_ident( 2, 2 ), 1.0, rel_tol );
  r_ident *= invert( r_45z );
  BOOST_CHECK( is_diagonal( r_ident, rel_tol ) );
  BOOST_CHECK_CLOSE( elem_norm_2( r_ident ), std::sqrt( 3.0 ), rel_tol );
  BOOST_CHECK_CLOSE( trace( r_45z_cpy ), std::sqrt( 2.0 ) + 1.0, rel_tol );
  BOOST_CHECK_CLOSE( determinant( r_45z_cpy ), 1.0, rel_tol );
  mat< double, mat_structure::symmetric > msym_45z( r_45z_cpy.getSymPart() );
  BOOST_CHECK( is_diagonal( msym_45z, rel_tol ) );
  BOOST_CHECK_CLOSE( elem_norm_2( msym_45z ), std::sqrt( 2.0 ), rel_tol );
  mat< double, mat_structure::skew_symmetric > mskw_45z( r_45z_cpy.getSkewSymPart() );
  BOOST_CHECK_SMALL( elem_norm_2( mskw_45z + transpose( mskw_45z ) ), rel_tol );
  BOOST_CHECK_SMALL( elem_norm_2( r_45z_cpy.getMat() - r_45z_cpy ), rel_tol );
  BOOST_CHECK_SMALL( elem_norm_2( invert( r_45z_cpy ).getMat() - transpose( r_45z_cpy ) ), rel_tol );
  BOOST_CHECK( is_diagonal( ( invert( r_45z_cpy ) * r_45z_cpy ), rel_tol ) );
  BOOST_CHECK_CLOSE( elem_norm_2( ( invert( r_45z_cpy ) * r_45z_cpy ) ), std::sqrt( 3.0 ), rel_tol );

  BOOST_CHECK( ( r_ident == r_ident ) );
  BOOST_CHECK( ( r_ident != r_45z ) );
  vect< double, 3 > v1( 1.0, 1.0, 2.0 );
  BOOST_CHECK_SMALL( norm_2( ( r_45z * v1 ) - vect< double, 3 >( 0.0, std::sqrt( 2.0 ), 2.0 ) ), 10.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( ( v1 * r_45z ) - vect< double, 3 >( std::sqrt( 2.0 ), 0.0, 2.0 ) ), 10.0 * rel_tol );

  quaternion< double > q_45z( r_45z );
  BOOST_CHECK_CLOSE( q_45z[0], std::cos( 0.125 * M_PI ), rel_tol );
  BOOST_CHECK_SMALL( q_45z[1], rel_tol );
  BOOST_CHECK_SMALL( q_45z[2], rel_tol );
  BOOST_CHECK_CLOSE( q_45z[3], std::sin( 0.125 * M_PI ), 10.0 * rel_tol );
  BOOST_CHECK( is_diagonal( ( r_45z * invert( q_45z ) ), rel_tol ) );
  BOOST_CHECK_CLOSE( elem_norm_2( r_45z * invert( q_45z ) ), std::sqrt( 3.0 ), rel_tol );
  quaternion< double > q_ident;
  BOOST_CHECK_CLOSE( q_ident[0], 1.0, rel_tol );
  BOOST_CHECK_SMALL( q_ident[1], rel_tol );
  BOOST_CHECK_SMALL( q_ident[2], rel_tol );
  BOOST_CHECK_SMALL( q_ident[3], rel_tol );
  quaternion< double > q_45z_cpy( q_45z );
  BOOST_CHECK_CLOSE( q_45z_cpy[0], std::cos( 0.125 * M_PI ), rel_tol );
  BOOST_CHECK_SMALL( q_45z_cpy[1], rel_tol );
  BOOST_CHECK_SMALL( q_45z_cpy[2], rel_tol );
  BOOST_CHECK_CLOSE( q_45z_cpy[3], std::sin( 0.125 * M_PI ), 10.0 * rel_tol );
  q_45z_cpy = q_45z;
  BOOST_CHECK_CLOSE( q_45z_cpy[0], std::cos( 0.125 * M_PI ), rel_tol );
  BOOST_CHECK_SMALL( q_45z_cpy[1], rel_tol );
  BOOST_CHECK_SMALL( q_45z_cpy[2], rel_tol );
  BOOST_CHECK_CLOSE( q_45z_cpy[3], std::sin( 0.125 * M_PI ), 10.0 * rel_tol );
  q_ident = q_45z * q_45z;
  BOOST_CHECK_CLOSE( q_ident[0], std::cos( 0.25 * M_PI ), 10.0 * rel_tol );
  BOOST_CHECK_SMALL( q_ident[1], rel_tol );
  BOOST_CHECK_SMALL( q_ident[2], rel_tol );
  BOOST_CHECK_CLOSE( q_ident[3], std::sin( 0.25 * M_PI ), rel_tol );
  q_ident *= invert( q_45z );
  BOOST_CHECK_CLOSE( q_ident[0], std::cos( 0.125 * M_PI ), 10.0 * rel_tol );
  BOOST_CHECK_SMALL( q_ident[1], rel_tol );
  BOOST_CHECK_SMALL( q_ident[2], rel_tol );
  BOOST_CHECK_CLOSE( q_ident[3], std::sin( 0.125 * M_PI ), 10.0 * rel_tol );
  q_ident *= invert( q_45z );
  BOOST_CHECK_CLOSE( q_ident[0], 1.0, 100.0 * rel_tol );
  BOOST_CHECK_SMALL( q_ident[1], rel_tol );
  BOOST_CHECK_SMALL( q_ident[2], rel_tol );
  BOOST_CHECK_SMALL( q_ident[3], rel_tol );
  mat< double, mat_structure::square > rm_90z( q_45z * r_45z );
  BOOST_CHECK_SMALL( rm_90z( 0, 0 ), rel_tol );
  BOOST_CHECK_CLOSE( rm_90z( 0, 1 ), -1.0, rel_tol );
  BOOST_CHECK_SMALL( rm_90z( 0, 2 ), rel_tol );
  BOOST_CHECK_CLOSE( rm_90z( 1, 0 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( rm_90z( 1, 1 ), rel_tol );
  BOOST_CHECK_SMALL( rm_90z( 1, 2 ), rel_tol );
  BOOST_CHECK_SMALL( rm_90z( 2, 0 ), rel_tol );
  BOOST_CHECK_SMALL( rm_90z( 2, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( rm_90z( 2, 2 ), 1.0, rel_tol );

  BOOST_CHECK_CLOSE( trace( q_45z_cpy ), 1.0 + std::sqrt( 2.0 ), rel_tol );
  BOOST_CHECK_CLOSE( determinant( q_45z_cpy ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( elem_norm_2( q_45z_cpy.getSymPart() - r_45z.getSymPart() ), rel_tol );
  BOOST_CHECK_SMALL( elem_norm_2( q_45z_cpy.getSkewSymPart() - r_45z.getSkewSymPart() ), rel_tol );
  BOOST_CHECK_CLOSE( ( q_45z_cpy * invert( q_45z_cpy ) )[0], 1.0, rel_tol );
  BOOST_CHECK_CLOSE( ( q_45z_cpy * transpose( q_45z_cpy ) )[0], 1.0, rel_tol );

  BOOST_CHECK( ( q_ident == q_ident ) );
  BOOST_CHECK( ( q_ident != q_45z ) );
  BOOST_CHECK_SMALL( norm_2( ( q_45z * v1 ) - vect< double, 3 >( 0.0, std::sqrt( 2.0 ), 2.0 ) ), 2.0 * rel_tol );

  euler_angles_TB< double > e_45z( 0.25 * M_PI, 0.0, 0.0 );
  quaternion< double > q_e_90z( q_45z * e_45z );
  BOOST_CHECK_CLOSE( q_e_90z[0], std::cos( 0.25 * M_PI ), 10.0 * rel_tol );
  BOOST_CHECK_SMALL( q_e_90z[1], rel_tol );
  BOOST_CHECK_SMALL( q_e_90z[2], rel_tol );
  BOOST_CHECK_CLOSE( q_e_90z[3], std::sin( 0.25 * M_PI ), rel_tol );
  axis_angle< double > a_45z( 0.25 * M_PI, vect< double, 3 >( 0.0, 0.0, 1.0 ) );
  quaternion< double > q_a_90z( q_45z * a_45z );
  BOOST_CHECK_CLOSE( q_a_90z[0], std::cos( 0.25 * M_PI ), 10.0 * rel_tol );
  BOOST_CHECK_SMALL( q_a_90z[1], rel_tol );
  BOOST_CHECK_SMALL( q_a_90z[2], rel_tol );
  BOOST_CHECK_CLOSE( q_a_90z[3], std::sin( 0.25 * M_PI ), rel_tol );

  euler_angles_TB< double > e_ident;
  BOOST_CHECK_SMALL( e_ident.yaw(), rel_tol );
  BOOST_CHECK_SMALL( e_ident.pitch(), rel_tol );
  BOOST_CHECK_SMALL( e_ident.roll(), rel_tol );
  euler_angles_TB< double > e_45z_cpy( e_45z );
  BOOST_CHECK_CLOSE( e_45z_cpy.yaw(), 0.25 * M_PI, rel_tol );
  BOOST_CHECK_SMALL( e_45z_cpy.pitch(), rel_tol );
  BOOST_CHECK_SMALL( e_45z_cpy.roll(), rel_tol );
  euler_angles_TB< double > e_45z_r( r_45z );
  BOOST_CHECK_CLOSE( e_45z_cpy.yaw(), 0.25 * M_PI, rel_tol );
  BOOST_CHECK_SMALL( e_45z_cpy.pitch(), rel_tol );
  BOOST_CHECK_SMALL( e_45z_cpy.roll(), rel_tol );
  euler_angles_TB< double > e_45z_q( q_45z );
  BOOST_CHECK_CLOSE( e_45z_cpy.yaw(), 0.25 * M_PI, rel_tol );
  BOOST_CHECK_SMALL( e_45z_cpy.pitch(), rel_tol );
  BOOST_CHECK_SMALL( e_45z_cpy.roll(), rel_tol );
  BOOST_CHECK_SMALL( elem_norm_2( e_45z_q.getRotMat().getMat() - r_45z.getMat() ), rel_tol );
  BOOST_CHECK_SMALL( elem_norm_2( e_45z_q.getMat() - r_45z.getMat() ), rel_tol );
  e_45z_cpy = e_45z;
  BOOST_CHECK_CLOSE( e_45z_cpy.yaw(), 0.25 * M_PI, rel_tol );
  BOOST_CHECK_SMALL( e_45z_cpy.pitch(), rel_tol );
  BOOST_CHECK_SMALL( e_45z_cpy.roll(), rel_tol );
  e_45z_cpy = r_45z;
  BOOST_CHECK_CLOSE( e_45z_cpy.yaw(), 0.25 * M_PI, rel_tol );
  BOOST_CHECK_SMALL( e_45z_cpy.pitch(), rel_tol );
  BOOST_CHECK_SMALL( e_45z_cpy.roll(), rel_tol );
  e_45z_cpy = q_45z;
  BOOST_CHECK_CLOSE( e_45z_cpy.yaw(), 0.25 * M_PI, rel_tol );
  BOOST_CHECK_SMALL( e_45z_cpy.pitch(), rel_tol );
  BOOST_CHECK_SMALL( e_45z_cpy.roll(), rel_tol );
  e_45z_cpy = a_45z;
  BOOST_CHECK_CLOSE( e_45z_cpy.yaw(), 0.25 * M_PI, 10.0 * rel_tol );
  BOOST_CHECK_SMALL( e_45z_cpy.pitch(), rel_tol );
  BOOST_CHECK_SMALL( e_45z_cpy.roll(), rel_tol );

  BOOST_CHECK( ( e_ident == e_ident ) );
  BOOST_CHECK( ( e_ident != e_45z ) );

  BOOST_CHECK_CLOSE( trace( e_45z_cpy ), 1.0 + std::sqrt( 2.0 ), rel_tol );
  BOOST_CHECK_CLOSE( determinant( e_45z_cpy ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( elem_norm_2( e_45z_cpy.getSymPart() - r_45z.getSymPart() ), rel_tol );
  BOOST_CHECK_SMALL( elem_norm_2( e_45z_cpy.getSkewSymPart() - r_45z.getSkewSymPart() ), rel_tol );
  BOOST_CHECK_CLOSE( ( q_45z_cpy * invert( e_45z_cpy ) )[0], 1.0, rel_tol );
  BOOST_CHECK_CLOSE( ( q_45z_cpy * transpose( e_45z_cpy ) )[0], 1.0, rel_tol );

  axis_angle< double > a_ident;
  BOOST_CHECK_SMALL( a_ident.angle(), rel_tol );
  axis_angle< double > a_45z_cpy( a_45z );
  BOOST_CHECK_CLOSE( a_45z_cpy.angle(), 0.25 * M_PI, rel_tol );
  BOOST_CHECK_SMALL( norm_2( a_45z_cpy.axis() - vect< double, 3 >( 0.0, 0.0, 1.0 ) ), rel_tol );
  axis_angle< double > a_45z_r( r_45z );
  BOOST_CHECK_CLOSE( a_45z_r.angle(), 0.25 * M_PI, 10.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( a_45z_r.axis() - vect< double, 3 >( 0.0, 0.0, 1.0 ) ), rel_tol );
  axis_angle< double > a_45z_q( q_45z );
  BOOST_CHECK_CLOSE( a_45z_q.angle(), 0.25 * M_PI, 10.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( a_45z_q.axis() - vect< double, 3 >( 0.0, 0.0, 1.0 ) ), rel_tol );
  axis_angle< double > a_45z_e( e_45z );
  BOOST_CHECK_CLOSE( a_45z_e.angle(), 0.25 * M_PI, 10.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( a_45z_e.axis() - vect< double, 3 >( 0.0, 0.0, 1.0 ) ), rel_tol );
  BOOST_CHECK_SMALL( elem_norm_2( a_45z_q.getRotMat().getMat() - r_45z.getMat() ), rel_tol );
  BOOST_CHECK_SMALL( elem_norm_2( a_45z_q.getMat() - r_45z.getMat() ), rel_tol );
  a_45z_cpy = a_45z;
  BOOST_CHECK_CLOSE( a_45z_cpy.angle(), 0.25 * M_PI, rel_tol );
  BOOST_CHECK_SMALL( norm_2( a_45z_cpy.axis() - vect< double, 3 >( 0.0, 0.0, 1.0 ) ), rel_tol );
  a_45z_cpy = r_45z;
  BOOST_CHECK_CLOSE( a_45z_cpy.angle(), 0.25 * M_PI, 10.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( a_45z_cpy.axis() - vect< double, 3 >( 0.0, 0.0, 1.0 ) ), rel_tol );
  a_45z_cpy = q_45z;
  BOOST_CHECK_CLOSE( a_45z_cpy.angle(), 0.25 * M_PI, 10.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( a_45z_cpy.axis() - vect< double, 3 >( 0.0, 0.0, 1.0 ) ), rel_tol );
  a_45z_cpy = e_45z;
  BOOST_CHECK_CLOSE( a_45z_cpy.angle(), 0.25 * M_PI, 10.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( a_45z_cpy.axis() - vect< double, 3 >( 0.0, 0.0, 1.0 ) ), rel_tol );

  BOOST_CHECK( ( a_ident == a_ident ) );
  BOOST_CHECK( ( a_ident != a_45z ) );

  BOOST_CHECK_CLOSE( trace( a_45z_cpy ), 1.0 + std::sqrt( 2.0 ), rel_tol );
  BOOST_CHECK_CLOSE( determinant( a_45z_cpy ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( elem_norm_2( a_45z_cpy.getSymPart() - r_45z.getSymPart() ), rel_tol );
  BOOST_CHECK_SMALL( elem_norm_2( a_45z_cpy.getSkewSymPart() - r_45z.getSkewSymPart() ), rel_tol );
  BOOST_CHECK_CLOSE( ( q_45z_cpy * invert( a_45z_cpy ) )[0], 1.0, rel_tol );
  BOOST_CHECK_CLOSE( ( q_45z_cpy * transpose( a_45z_cpy ) )[0], 1.0, rel_tol );

  trans_mat_3D< double > t_ident;
  BOOST_CHECK( is_diagonal( t_ident, rel_tol ) );
  BOOST_CHECK_CLOSE( elem_norm_2( t_ident ), std::sqrt( 4.0 ), rel_tol );
  double t_45z_array[]
    = {std::cos( 0.25 * M_PI ), std::sin( 0.25 * M_PI ), 0.0, 0.0, -std::sin( 0.25 * M_PI ), std::cos( 0.25 * M_PI ),
       0.0,                     0.0,                     0.0, 0.0, 1.0,                      0.0,
       0.0,                     0.0,                     0.0, 1.0};
  trans_mat_3D< double > t_45z( t_45z_array );
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z( 3, 3 ), 1.0, rel_tol );
  trans_mat_3D< double > t_45z_cpy( t_45z );
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_cpy ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 3, 3 ), 1.0, rel_tol );
  trans_mat_3D< double > t_45z_r( r_45z );
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_r ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_r( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_r( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_r( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_r( 3, 3 ), 1.0, rel_tol );
  trans_mat_3D< double > t_45z_q( q_45z );
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_q ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_q( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_q( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_q( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_q( 3, 3 ), 1.0, rel_tol );
  trans_mat_3D< double > t_45z_e( e_45z );
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_e ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_e( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_e( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_e( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_e( 3, 3 ), 1.0, rel_tol );
  trans_mat_3D< double > t_45z_a( a_45z );
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_a ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_a( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_a( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_a( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_a( 3, 3 ), 1.0, rel_tol );
  rot_mat_3D< double > t_45z_q_rot = t_45z_q.getRotMat();
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_q_rot ), std::sqrt( 3.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_q_rot( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_q_rot( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_q_rot( 2, 2 ), 1.0, rel_tol );
  mat< double, mat_structure::square > t_45z_q_mat = t_45z_q.getMat();
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_q_mat ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_q_mat( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_q_mat( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_q_mat( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_q_mat( 3, 3 ), 1.0, rel_tol );

  t_45z_cpy = t_45z;
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_cpy ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 3, 3 ), 1.0, rel_tol );
  t_45z_cpy = r_45z;
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_cpy ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 3, 3 ), 1.0, rel_tol );
  t_45z_cpy = q_45z;
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_cpy ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 3, 3 ), 1.0, rel_tol );
  t_45z_cpy = e_45z;
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_cpy ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 3, 3 ), 1.0, rel_tol );
  t_45z_cpy = a_45z;
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_cpy ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 3, 3 ), 1.0, rel_tol );

  t_45z_cpy = t_ident;
  t_45z_cpy *= t_45z;
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_cpy ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 3, 3 ), 1.0, rel_tol );
  t_45z_cpy = t_ident;
  t_45z_cpy *= r_45z;
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_cpy ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 3, 3 ), 1.0, rel_tol );
  t_45z_cpy = t_ident;
  t_45z_cpy *= q_45z;
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_cpy ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 3, 3 ), 1.0, rel_tol );
  t_45z_cpy = t_ident;
  t_45z_cpy *= e_45z;
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_cpy ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 3, 3 ), 1.0, rel_tol );
  t_45z_cpy = t_ident;
  t_45z_cpy *= a_45z;
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_cpy ), std::sqrt( 4.0 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_cpy( 3, 3 ), 1.0, rel_tol );


  BOOST_CHECK( ( t_ident == t_ident ) );
  BOOST_CHECK( ( t_ident != t_45z ) );

  trans_mat_3D< double > t_45z_123( r_45z, vect< double, 3 >( 1.0, 2.0, 3.0 ) );
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_123 ), std::sqrt( 18.0 ), 2.0 * rel_tol );
  BOOST_CHECK_CLOSE( t_45z_123( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_123( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_123( 0, 3 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_123( 1, 3 ), 2.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_123( 2, 3 ), 3.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_123( 3, 3 ), 1.0, rel_tol );

  BOOST_CHECK_CLOSE( trace( t_45z_123 ), 2.0 + std::sqrt( 2.0 ), rel_tol );
  BOOST_CHECK_CLOSE( determinant( t_45z_123 ), 1.0, rel_tol );
  mat< double, mat_structure::square > t_45z_mat_sym_skw( t_45z_123.getSymPart() + t_45z_123.getSkewSymPart() );
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_mat_sym_skw ), std::sqrt( 18.0 ), 2.0 * rel_tol );
  BOOST_CHECK_CLOSE( t_45z_mat_sym_skw( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_mat_sym_skw( 0, 1 ), -std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_mat_sym_skw( 0, 3 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_mat_sym_skw( 1, 3 ), 2.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_mat_sym_skw( 2, 3 ), 3.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_mat_sym_skw( 3, 3 ), 1.0, rel_tol );
  t_ident = t_45z_123 * invert( t_45z_123 );
  BOOST_CHECK_CLOSE( elem_norm_2( t_ident ), std::sqrt( 4.0 ), 2.0 * rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 0, 0 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( t_ident( 0, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 1, 1 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( t_ident( 0, 3 ), rel_tol );
  BOOST_CHECK_SMALL( t_ident( 1, 3 ), rel_tol );
  BOOST_CHECK_SMALL( t_ident( 2, 3 ), rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 3, 3 ), 1.0, rel_tol );
  t_ident = ( invert( t_45z_123 ).getMat() * t_45z_123 );
  BOOST_CHECK_CLOSE( elem_norm_2( t_ident ), std::sqrt( 4.0 ), 2.0 * rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 0, 0 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( t_ident( 0, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 1, 1 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( t_ident( 0, 3 ), rel_tol );
  BOOST_CHECK_SMALL( t_ident( 1, 3 ), rel_tol );
  BOOST_CHECK_SMALL( t_ident( 2, 3 ), rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 3, 3 ), 1.0, rel_tol );
  t_ident = ( invert( t_45z_123 ) * t_45z_123.getMat() );
  BOOST_CHECK_CLOSE( elem_norm_2( t_ident ), std::sqrt( 4.0 ), 2.0 * rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 0, 0 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( t_ident( 0, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 1, 1 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( t_ident( 0, 3 ), rel_tol );
  BOOST_CHECK_SMALL( t_ident( 1, 3 ), rel_tol );
  BOOST_CHECK_SMALL( t_ident( 2, 3 ), rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 3, 3 ), 1.0, rel_tol );

  mat< double, mat_structure::square > t_45z_123_t = transpose( t_45z_123 );
  BOOST_CHECK_CLOSE( elem_norm_2( t_45z_123_t ), std::sqrt( 18.0 ), 2.0 * rel_tol );
  BOOST_CHECK_CLOSE( t_45z_123_t( 0, 0 ), std::cos( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_123_t( 0, 1 ), std::sin( 0.25 * M_PI ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_123_t( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( t_45z_123_t( 0, 3 ), rel_tol );
  BOOST_CHECK_CLOSE( t_45z_123_t( 3, 1 ), 2.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_123_t( 3, 2 ), 3.0, rel_tol );
  BOOST_CHECK_CLOSE( t_45z_123_t( 3, 3 ), 1.0, rel_tol );

  t_ident = ( invert( t_45z_123 ) * r_45z ) * trans_mat_3D< double >( r_ident, -invert( t_45z_123 ).getTranslation() );
  BOOST_CHECK_CLOSE( elem_norm_2( t_ident ), std::sqrt( 4.0 ), 2.0 * rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 0, 0 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( t_ident( 0, 1 ), rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 1, 1 ), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 2, 2 ), 1.0, rel_tol );
  BOOST_CHECK_SMALL( t_ident( 0, 3 ), rel_tol );
  BOOST_CHECK_SMALL( t_ident( 1, 3 ), rel_tol );
  BOOST_CHECK_SMALL( t_ident( 2, 3 ), rel_tol );
  BOOST_CHECK_CLOSE( t_ident( 3, 3 ), 1.0, rel_tol );

  vect< double, 3 > v1_trans( t_45z_123 * v1 );
  BOOST_CHECK_CLOSE( v1_trans[0], 1.0, rel_tol );
  BOOST_CHECK_CLOSE( v1_trans[1], 2.0 + std::sqrt( 2.0 ), rel_tol );
  BOOST_CHECK_CLOSE( v1_trans[2], 5.0, rel_tol );
  v1_trans = t_45z_123.rotate( v1 );
  BOOST_CHECK_SMALL( v1_trans[0], rel_tol );
  BOOST_CHECK_CLOSE( v1_trans[1], std::sqrt( 2.0 ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( v1_trans[2], 2.0, rel_tol );
  vect< double, 4 > v3 = t_45z_123 * vect< double, 4 >( 1, 1, 2, 2 );

  BOOST_CHECK_CLOSE( v3[0], 2.0, rel_tol );
  BOOST_CHECK_CLOSE( v3[1], 4.0 + std::sqrt( 2.0 ), rel_tol );
  BOOST_CHECK_CLOSE( v3[2], 8.0, rel_tol );
  BOOST_CHECK_CLOSE( v3[3], 2.0, rel_tol );

  axis_angle< double > a_weird( 0.3241, vect< double, 3 >( 0.5, 0.5, sqrt( 0.5 ) ) );
  quaternion< double > q_weird;
  q_weird = a_weird;
  euler_angles_TB< double > e_weird;
  e_weird = q_weird;
  rot_mat_3D< double > r_weird;
  r_weird = q_weird;
  axis_angle< double > a_weird_q( q_weird );
  axis_angle< double > a_weird_e( e_weird );
  axis_angle< double > a_weird_r( r_weird );
  BOOST_CHECK_CLOSE( a_weird_q.angle(), a_weird.angle(), 100.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( a_weird_q.axis() - a_weird.axis() ), 10.0 * rel_tol );
  BOOST_CHECK_CLOSE( a_weird_e.angle(), a_weird.angle(), 500.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( a_weird_e.axis() - a_weird.axis() ), 20.0 * rel_tol );
  BOOST_CHECK_CLOSE( a_weird_r.angle(), a_weird.angle(), 100.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( a_weird_r.axis() - a_weird.axis() ), 20.0 * rel_tol );

  quaternion< double > q_res( q_45z * quaternion< double >( a_weird * a_45z ) );
  vect< double, 4 > v_a( q_res[0], q_res[1], q_res[2], q_res[3] );
  q_res = q_45z * ( q_weird * q_45z );
  vect< double, 4 > v_q( q_res[0], q_res[1], q_res[2], q_res[3] );
  q_res = quaternion< double >( q_45z * e_weird * r_45z );
  vect< double, 4 > v_e( q_res[0], q_res[1], q_res[2], q_res[3] );
  q_res = quaternion< double >( a_45z * r_weird * e_45z );
  vect< double, 4 > v_r( q_res[0], q_res[1], q_res[2], q_res[3] );
  BOOST_CHECK_SMALL( norm_2( v_a - v_q ), 2.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( v_a - v_e ), 2.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( v_a - v_r ), 2.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( v_q - v_e ), 2.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( v_q - v_r ), 2.0 * rel_tol );
  BOOST_CHECK_SMALL( norm_2( v_e - v_r ), 2.0 * rel_tol );
};

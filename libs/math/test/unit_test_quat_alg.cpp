
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


#include <ReaK/math/kinetostatics/quat_alg.hpp>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif


#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE quat_alg
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


BOOST_AUTO_TEST_CASE( quat_tests )
{
  using namespace ReaK;
  const double rel_tol = std::numeric_limits<double>::epsilon(); RK_UNUSED(rel_tol);
  
  quat<double> q_45z(std::cos(0.125 * M_PI), 0.0, 0.0, std::sin(0.125 * M_PI));
  BOOST_CHECK_CLOSE( q_45z[0], std::cos(0.125 * M_PI), rel_tol );
  BOOST_CHECK_SMALL( q_45z[1], rel_tol );
  BOOST_CHECK_SMALL( q_45z[2], rel_tol );
  BOOST_CHECK_CLOSE( q_45z[3], std::sin(0.125 * M_PI), rel_tol );
  quat<double> q_ident = q_45z * conj(q_45z);
  BOOST_CHECK_CLOSE( q_ident[0], 1.0, rel_tol );
  BOOST_CHECK_SMALL( q_ident[1], rel_tol );
  BOOST_CHECK_SMALL( q_ident[2], rel_tol );
  BOOST_CHECK_SMALL( q_ident[3], rel_tol );
  quat<double> q_zero = q_45z - q_45z;
  BOOST_CHECK_SMALL( q_zero[0], rel_tol );
  BOOST_CHECK_SMALL( q_zero[1], rel_tol );
  BOOST_CHECK_SMALL( q_zero[2], rel_tol );
  BOOST_CHECK_SMALL( q_zero[3], rel_tol );
  q_zero = q_45z + (-q_45z);
  BOOST_CHECK_SMALL( q_zero[0], rel_tol );
  BOOST_CHECK_SMALL( q_zero[1], rel_tol );
  BOOST_CHECK_SMALL( q_zero[2], rel_tol );
  BOOST_CHECK_SMALL( q_zero[3], rel_tol );
  quat<double> q_45z_cpy(q_45z);
  BOOST_CHECK_CLOSE( q_45z_cpy[0], std::cos(0.125 * M_PI), rel_tol );
  BOOST_CHECK_SMALL( q_45z_cpy[1], rel_tol );
  BOOST_CHECK_SMALL( q_45z_cpy[2], rel_tol );
  BOOST_CHECK_CLOSE( q_45z_cpy[3], std::sin(0.125 * M_PI), rel_tol );
  q_45z_cpy = q_45z;
  BOOST_CHECK_CLOSE( q_45z_cpy[0], std::cos(0.125 * M_PI), rel_tol );
  BOOST_CHECK_SMALL( q_45z_cpy[1], rel_tol );
  BOOST_CHECK_SMALL( q_45z_cpy[2], rel_tol );
  BOOST_CHECK_CLOSE( q_45z_cpy[3], std::sin(0.125 * M_PI), rel_tol );
  q_ident = q_45z * q_45z;
  BOOST_CHECK_CLOSE( q_ident[0], std::cos(0.25 * M_PI), 100.0 * rel_tol );
  BOOST_CHECK_SMALL( q_ident[1], rel_tol );
  BOOST_CHECK_SMALL( q_ident[2], rel_tol );
  BOOST_CHECK_CLOSE( q_ident[3], std::sin(0.25 * M_PI), 100.0 * rel_tol );
  q_ident *= conj(q_45z);
  BOOST_CHECK_CLOSE( q_ident[0], std::cos(0.125 * M_PI), 100.0 * rel_tol );
  BOOST_CHECK_SMALL( q_ident[1], rel_tol );
  BOOST_CHECK_SMALL( q_ident[2], rel_tol );
  BOOST_CHECK_CLOSE( q_ident[3], std::sin(0.125 * M_PI), 100.0 * rel_tol );
  q_ident *= invert(q_45z);
  BOOST_CHECK_CLOSE( q_ident[0], 1.0, rel_tol );
  BOOST_CHECK_SMALL( q_ident[1], rel_tol );
  BOOST_CHECK_SMALL( q_ident[2], rel_tol );
  BOOST_CHECK_SMALL( q_ident[3], rel_tol );
  BOOST_CHECK_CLOSE( norm_2_sqr(q_45z), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( norm_2(q_45z), 1.0, rel_tol );
  q_45z *= 2.0;
  q_45z = unit(q_45z);
  BOOST_CHECK_CLOSE( q_45z[0], std::cos(0.125 * M_PI), rel_tol );
  BOOST_CHECK_SMALL( q_45z[1], rel_tol );
  BOOST_CHECK_SMALL( q_45z[2], rel_tol );
  BOOST_CHECK_CLOSE( q_45z[3], std::sin(0.125 * M_PI), rel_tol );
  
  q_45z = sqrt( pow(q_45z, quat<double>(2.0)) );
  BOOST_CHECK_CLOSE( q_45z[0], std::cos(0.125 * M_PI), 2.0 * rel_tol );
  BOOST_CHECK_SMALL( q_45z[1], 2.0 * rel_tol );
  BOOST_CHECK_SMALL( q_45z[2], 2.0 * rel_tol );
  BOOST_CHECK_CLOSE( q_45z[3], std::sin(0.125 * M_PI), 2.0 * rel_tol );
  
  quat<double> temp = cos(q_45z) * cos(q_45z) + sin(q_45z) * sin(q_45z);
  BOOST_CHECK_CLOSE( temp[0], 1.0, 100.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[1], rel_tol );
  BOOST_CHECK_SMALL( temp[2], rel_tol );
  BOOST_CHECK_SMALL( temp[3], rel_tol );
  temp =  invert(cos(q_45z) * cos(q_45z)) - tan(q_45z) * tan(q_45z);
  BOOST_CHECK_CLOSE( temp[0], 1.0, 200.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[1], 10.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[2], 10.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[3], 10.0 * rel_tol );
  temp =  acos(cos(q_45z)) * invert(q_45z);
  BOOST_CHECK_CLOSE( temp[0], 1.0, 100.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[1], 10.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[2], 10.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[3], 10.0 * rel_tol );
  temp =  asin(sin(q_45z)) * invert(q_45z);
  BOOST_CHECK_CLOSE( temp[0], 1.0, 100.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[1], 10.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[2], 10.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[3], 10.0 * rel_tol );
  temp =  atan(tan(q_45z)) * invert(q_45z);
  BOOST_CHECK_CLOSE( temp[0], 1.0, 100.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[1], 10.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[2], 10.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[3], 10.0 * rel_tol );
  temp =  exp(q_45z) * exp(q_45z) - exp(q_45z + q_45z);
  BOOST_CHECK_SMALL( temp[0], 10.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[1], 10.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[2], 10.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[3], 10.0 * rel_tol );
  temp =  log(q_45z) + log(q_45z) - log(q_45z * q_45z);
  BOOST_CHECK_SMALL( temp[0], 10.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[1], 10.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[2], 10.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[3], 10.0 * rel_tol );
  
};
  
BOOST_AUTO_TEST_CASE( unit_quat_tests )
{
  using namespace ReaK;
  const double rel_tol = std::numeric_limits<double>::epsilon(); RK_UNUSED(rel_tol);
  
  unit_quat<double> q_45z(std::cos(0.125 * M_PI), 0.0, 0.0, std::sin(0.125 * M_PI));
  BOOST_CHECK_CLOSE( q_45z[0], std::cos(0.125 * M_PI), rel_tol );
  BOOST_CHECK_SMALL( q_45z[1], rel_tol );
  BOOST_CHECK_SMALL( q_45z[2], rel_tol );
  BOOST_CHECK_CLOSE( q_45z[3], std::sin(0.125 * M_PI), rel_tol );
  unit_quat<double> q_ident = q_45z * conj(q_45z);
  BOOST_CHECK_CLOSE( q_ident[0], 1.0, rel_tol );
  BOOST_CHECK_SMALL( q_ident[1], rel_tol );
  BOOST_CHECK_SMALL( q_ident[2], rel_tol );
  BOOST_CHECK_SMALL( q_ident[3], rel_tol );
  unit_quat<double> q_45z_cpy(q_45z);
  BOOST_CHECK_CLOSE( q_45z_cpy[0], std::cos(0.125 * M_PI), rel_tol );
  BOOST_CHECK_SMALL( q_45z_cpy[1], rel_tol );
  BOOST_CHECK_SMALL( q_45z_cpy[2], rel_tol );
  BOOST_CHECK_CLOSE( q_45z_cpy[3], std::sin(0.125 * M_PI), rel_tol );
  q_45z_cpy = q_45z;
  BOOST_CHECK_CLOSE( q_45z_cpy[0], std::cos(0.125 * M_PI), rel_tol );
  BOOST_CHECK_SMALL( q_45z_cpy[1], rel_tol );
  BOOST_CHECK_SMALL( q_45z_cpy[2], rel_tol );
  BOOST_CHECK_CLOSE( q_45z_cpy[3], std::sin(0.125 * M_PI), rel_tol );
  q_ident = q_45z * q_45z;
  BOOST_CHECK_CLOSE( q_ident[0], std::cos(0.25 * M_PI), 100.0 * rel_tol );
  BOOST_CHECK_SMALL( q_ident[1], rel_tol );
  BOOST_CHECK_SMALL( q_ident[2], rel_tol );
  BOOST_CHECK_CLOSE( q_ident[3], std::sin(0.25 * M_PI), 100.0 * rel_tol );
  q_ident *= conj(q_45z);
  BOOST_CHECK_CLOSE( q_ident[0], std::cos(0.125 * M_PI), 100.0 * rel_tol );
  BOOST_CHECK_SMALL( q_ident[1], rel_tol );
  BOOST_CHECK_SMALL( q_ident[2], rel_tol );
  BOOST_CHECK_CLOSE( q_ident[3], std::sin(0.125 * M_PI), 100.0 * rel_tol );
  q_ident *= invert(q_45z);
  BOOST_CHECK_CLOSE( q_ident[0], 1.0, rel_tol );
  BOOST_CHECK_SMALL( q_ident[1], rel_tol );
  BOOST_CHECK_SMALL( q_ident[2], rel_tol );
  BOOST_CHECK_SMALL( q_ident[3], rel_tol );
  BOOST_CHECK_CLOSE( norm_2_sqr(q_45z), 1.0, rel_tol );
  BOOST_CHECK_CLOSE( norm_2(q_45z), 1.0, rel_tol );
  
  vect<double,3> v_45z(0.0, 0.0, 0.125 * M_PI);
  quat<double> temp =  (exp(v_45z) * exp(v_45z)) * conj( exp(v_45z + v_45z) );
  BOOST_CHECK_CLOSE( temp[0], 1.0, 5.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[1], 5.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[2], 5.0 * rel_tol );
  BOOST_CHECK_SMALL( temp[3], 5.0 * rel_tol );
  vect<double,3> v_zero = log(q_45z) + log(q_45z) - log(q_45z * q_45z);
  BOOST_CHECK_SMALL( v_zero[0], 5.0 * rel_tol );
  BOOST_CHECK_SMALL( v_zero[1], 5.0 * rel_tol );
  BOOST_CHECK_SMALL( v_zero[2], 5.0 * rel_tol );
  
  q_45z = sqrt( q_45z * q_45z );
  BOOST_CHECK_CLOSE( q_45z[0], std::cos(0.125 * M_PI), 5.0 * rel_tol );
  BOOST_CHECK_SMALL( q_45z[1], 5.0 * rel_tol );
  BOOST_CHECK_SMALL( q_45z[2], 5.0 * rel_tol );
  BOOST_CHECK_CLOSE( q_45z[3], std::sin(0.125 * M_PI), 5.0 * rel_tol );
  
  q_45z = sqrt( pow(q_45z, 2.0) );
  BOOST_CHECK_CLOSE( q_45z[0], std::cos(0.125 * M_PI), 5.0 * rel_tol );
  BOOST_CHECK_SMALL( q_45z[1], 5.0 * rel_tol );
  BOOST_CHECK_SMALL( q_45z[2], 5.0 * rel_tol );
  BOOST_CHECK_CLOSE( q_45z[3], std::sin(0.125 * M_PI), 5.0 * rel_tol );
  
};
















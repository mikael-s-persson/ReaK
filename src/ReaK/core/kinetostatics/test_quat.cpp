
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


#include "quat_alg.hpp"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

int main() {
  ReaK::quat<double> q;
  using std::fabs;
  
  try {
    {
      ReaK::quat<double> q_45z(std::cos(0.125 * M_PI), 0.0, 0.0, std::sin(0.125 * M_PI));
      if( ( fabs( q_45z[0] - std::cos(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[3] - std::sin(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("q_45z does not have the correct elements!");
	return 1;
      };
      ReaK::quat<double> q_ident = q_45z * conj(q_45z);
      if( ( fabs( q_ident[0] - 1.0) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[3] ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("Composition of 45degZ-quaternion and its inverse is not the identity quaternion!");
        return 1;
      };
      ReaK::quat<double> q_zero = q_45z - q_45z;
      if( ( fabs( q_zero[0] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_zero[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_zero[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_zero[3] ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("q_zero does not have the correct elements!");
	return 1;
      };
      q_zero = q_45z + (-q_45z);
      if( ( fabs( q_zero[0] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_zero[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_zero[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_zero[3] ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("q_zero does not have the correct elements!");
	return 1;
      };
      ReaK::quat<double> q_45z_cpy(q_45z);
      if( ( fabs( q_45z_cpy[0] - std::cos(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z_cpy[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z_cpy[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z_cpy[3] - std::sin(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("q_45z copy does not have the correct elements!");
	return 1;
      };
      q_45z_cpy = q_45z;
      if( ( fabs( q_45z_cpy[0] - std::cos(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z_cpy[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z_cpy[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z_cpy[3] - std::sin(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("q_45z copy-assigned does not have the correct elements!");
	return 1;
      };
      q_ident = q_45z * q_45z;
      if( ( fabs( q_ident[0] - std::cos(0.25 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[3] - std::sin(0.25 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("q_45z * q_45z does not have the correct elements!");
	return 1;
      };
      q_ident *= conj(q_45z);
      if( ( fabs( q_ident[0] - std::cos(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[3] - std::sin(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("q_45z * q_45z * conj(q_45z) does not have the correct elements!");
	return 1;
      };
      q_ident *= invert(q_45z);
      if( ( fabs( q_ident[0] - 1.0) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[3] ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("q_45z * q_45z * invert(q_45z) * invert(q_45z) does not have the correct elements!");
	return 1;
      };
      if( ( fabs( norm_sqr(q_45z) - 1.0) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("norm_sqr(q_45z) does not have the correct value!");
	return 1;
      };
      if( ( fabs( norm(q_45z) - 1.0) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("norm(q_45z) does not have the correct value!");
	return 1;
      };
      q_45z *= 2.0;
      q_45z = unit(q_45z);
      if( ( fabs( q_45z[0] - std::cos(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[3] - std::sin(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("unit(q_45z) does not have the correct elements!");
	return 1;
      };
      
      q_45z = sqrt( pow(q_45z, ReaK::quat<double>(2.0)) );
      if( ( fabs( q_45z[0] - std::cos(0.125 * M_PI)) > 2.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[1] ) > 2.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[2] ) > 2.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[3] - std::sin(0.125 * M_PI)) > 2.0 * std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("sqrt(pow(q_45z,2.0)) does not have the correct elements!");
	return 1;
      };
      
      
      ReaK::quat<double> temp = cos(q_45z) * cos(q_45z) + sin(q_45z) * sin(q_45z);
      if( ( fabs( temp[0] - 1.0) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[3] ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("The identity: 'cos(q)^2 + sin(q)^2 = 1' does not hold!");
	return 1;
      };
      temp =  invert(cos(q_45z) * cos(q_45z)) - tan(q_45z) * tan(q_45z);
      if( ( fabs( temp[0] - 1.0) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[1] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[2] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[3] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("The identity: 'sec(q)^2 - tan(q)^2 = 1' does not hold!");
	return 1;
      };
      temp =  acos(cos(q_45z)) * invert(q_45z);
      if( ( fabs( temp[0] - 1.0) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[1] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[2] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[3] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("The identity: 'acos(cos(q)) = q' does not hold!");
	return 1;
      };
      temp =  asin(sin(q_45z)) * invert(q_45z);
      if( ( fabs( temp[0] - 1.0) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[1] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[2] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[3] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("The identity: 'asin(sin(q)) = q' does not hold!");
	return 1;
      };
      temp =  atan(tan(q_45z)) * invert(q_45z);
      if( ( fabs( temp[0] - 1.0) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[1] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[2] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[3] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("The identity: 'atan(tan(q)) = q' does not hold!");
	return 1;
      };
      temp =  exp(q_45z) * exp(q_45z) - exp(q_45z + q_45z);
      if( ( fabs( temp[0] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[1] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[2] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[3] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("The identity: 'exp(q) * exp(q) = exp(q + q)' does not hold!");
	return 1;
      };
      temp =  log(q_45z) + log(q_45z) - log(q_45z * q_45z);
      if( ( fabs( temp[0] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[1] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[2] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[3] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("The identity: 'log(q) + log(q) = log(q * q)' does not hold!");
	return 1;
      };
    };  

    {
      ReaK::unit_quat<double> q_45z(std::cos(0.125 * M_PI), 0.0, 0.0, std::sin(0.125 * M_PI));
      if( ( fabs( q_45z[0] - std::cos(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[3] - std::sin(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("uq_45z does not have the correct elements!");
	return 1;
      };
      ReaK::unit_quat<double> q_ident = q_45z * conj(q_45z);
      if( ( fabs( q_ident[0] - 1.0) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[3] ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("Composition of uq_45z and its inverse is not the identity quaternion!");
        return 1;
      };
      ReaK::unit_quat<double> q_45z_cpy(q_45z);
      if( ( fabs( q_45z_cpy[0] - std::cos(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z_cpy[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z_cpy[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z_cpy[3] - std::sin(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("uq_45z copy does not have the correct elements!");
	return 1;
      };
      q_45z_cpy = q_45z;
      if( ( fabs( q_45z_cpy[0] - std::cos(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z_cpy[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z_cpy[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z_cpy[3] - std::sin(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("uq_45z copy-assigned does not have the correct elements!");
	return 1;
      };
      q_ident = q_45z * q_45z;
      if( ( fabs( q_ident[0] - std::cos(0.25 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[3] - std::sin(0.25 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("uq_45z * uq_45z does not have the correct elements!");
	return 1;
      };
      q_ident *= conj(q_45z);
      if( ( fabs( q_ident[0] - std::cos(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[3] - std::sin(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("uq_45z * uq_45z * conj(uq_45z) does not have the correct elements!");
	return 1;
      };
      q_ident *= invert(q_45z);
      if( ( fabs( q_ident[0] - 1.0) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[3] ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("uq_45z * uq_45z * invert(uq_45z) * invert(uq_45z) does not have the correct elements!");
	return 1;
      };
      if( ( fabs( norm_sqr(q_45z) - 1.0) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("norm_sqr(uq_45z) does not have the correct value!");
	return 1;
      };
      if( ( fabs( norm(q_45z) - 1.0) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("norm(uq_45z) does not have the correct value!");
	return 1;
      };

      ReaK::vect<double,3> v_45z(0.0, 0.0, 0.125 * M_PI);
      ReaK::quat<double> temp =  (exp(v_45z) * exp(v_45z)) * conj( exp(v_45z + v_45z) );
      if( ( fabs( temp[0] - 1.0 ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[1] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[2] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( temp[3] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("The identity: 'exp(q) * exp(q) = exp(q + q)' does not hold!");
	return 1;
      };
      ReaK::vect<double,3> v_zero = log(q_45z) + log(q_45z) - log(q_45z * q_45z);
      if( ( fabs( v_zero[0] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( v_zero[1] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( v_zero[2] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("The identity: 'log(q) + log(q) = log(q * q)' does not hold!");
	return 1;
      };
      
      q_45z = sqrt( q_45z * q_45z );
      if( ( fabs( q_45z[0] - std::cos(0.125 * M_PI)) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[1] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[2] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[3] - std::sin(0.125 * M_PI)) > 5.0 * std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("sqrt( q_45z * q_45z ) does not have the correct elements!");
	return 1;
      };

      q_45z = sqrt( pow(q_45z, 2.0) );
      if( ( fabs( q_45z[0] - std::cos(0.125 * M_PI)) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[1] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[2] ) > 5.0 * std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[3] - std::sin(0.125 * M_PI)) > 5.0 * std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("sqrt(pow(q_45z,2.0)) does not have the correct elements!");
	return 1;
      };
     
    };  
  } catch(std::exception& e) {
    RK_ERROR("An exception has occurred during the math_gen test: '" << e.what() << "'");
    return 1;
  } catch(...) {
    RK_ERROR("An unexpected and unidentified exception has occurred during the math_gen test.");
    return 1;
  };
  
  return 0;
};
















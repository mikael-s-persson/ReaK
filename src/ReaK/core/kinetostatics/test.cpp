
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
#include <sstream>

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_norms.hpp"

#include "rotations.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979324
#define M_PI_2 1.57079632679489662
#endif

int main() {
  using std::fabs;
  using namespace ReaK;
  
  try {

    if(true){
      rot_mat_2D<double> r_ident;
      std::stringstream ss_r_ident; ss_r_ident << r_ident;
      if( ss_r_ident.str() != "(angle = 0)" ) {
	RK_ERROR("No-rotation2D matrix stream-output is incorrect!");
	return 1;
      };
      if( fabs(r_ident.getAngle()) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("No-rotation2D matrix does not have zero angle!");
	return 1;
      };
      if( ( fabs(r_ident(0,0) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_ident(0,1)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_ident(1,0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_ident(1,1) - 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("No-rotation2D matrix does not have identity elements!");
        return 1;
      };
      rot_mat_2D<double> r_45deg(0.25f * double(M_PI));
      if( fabs(r_45deg.getAngle() - 0.25f * double(M_PI)) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("45deg-rotation2D matrix does not have 45deg angle!");
	return 1;
      };
      if( ( fabs(r_45deg(0,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45deg(0,1) + std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45deg(1,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45deg(1,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("45deg-rotation2D matrix does not have correct elements!");
        return 1;
      };
      rot_mat_2D<double> r_45deg_cpy(r_45deg);
      if( fabs(r_45deg_cpy.getAngle() - 0.25f * double(M_PI)) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("45deg-rotation2D matrix copy does not have 45deg angle!");
	return 1;
      };
      mat<double,mat_structure::square> m_45deg = r_45deg.getMat();
      if( ( fabs(m_45deg(0,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m_45deg(0,1) + std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m_45deg(1,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m_45deg(1,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("45deg-rotation2D matrix does not produce a correct regular matrix!");
        return 1;
      };
      r_45deg_cpy.setAngle(0.0);
      if( fabs(r_45deg_cpy.getAngle()) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Rotation2D matrix set to zero does not have zero angle!");
	return 1;
      };
      r_45deg_cpy = r_45deg;
      rot_mat_2D<double> r_90deg = r_45deg * r_45deg_cpy;
      if( fabs(r_90deg.getAngle() - 0.5f * double(M_PI)) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Composition of two 45deg-rotation2D matrices does not have 90deg angle!");
	return 1;
      };
      r_90deg *= r_45deg;
      if( fabs(r_90deg.getAngle() - 0.75f * double(M_PI)) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Composition of three 45deg-rotation2D matrices does not have 135deg angle!");
	return 1;
      };
      vect<double,2> v1 = r_45deg * vect<double,2>(1.0,1.0);
      if( ( fabs(v1[0]) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(v1[1] - std::sqrt(2)) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Vector (1,1) rotated by 45deg does not produce the correct vector!");
        return 1;
      };
      v1 = vect<double,2>(1.0,1.0) * r_45deg;
      if( ( fabs(v1[1]) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(v1[0] - std::sqrt(2)) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Vector (1,1) transposed rotated by 45deg does not produce the correct vector!");
        return 1;
      };
      r_90deg *= invert(r_45deg);
      if( fabs(r_90deg.getAngle() - 0.5f * double(M_PI)) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Composition of a 135deg-rotation2D and a -45deg-rotation2D does not have 90deg angle!");
	return 1;
      };
      if( fabs(trace(r_45deg_cpy) - std::sqrt(2)) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Trace of a 45deg-rotation2D is not correct!");
	return 1;
      };
      if( fabs(determinant(r_45deg_cpy) - 1.0) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Determinant of a 45deg-rotation2D is not correct!");
	return 1;
      };
      const mat<double,mat_structure::symmetric> msym_45deg = r_45deg_cpy.getSymPart();
      if( ( fabs(msym_45deg(0,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(msym_45deg(0,1)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(msym_45deg(1,0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(msym_45deg(1,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Symmetric part of 45deg-rotation2D matrix is not correct!");
        return 1;
      };
      const mat<double,mat_structure::skew_symmetric> mskw_45deg = r_45deg_cpy.getSkewSymPart();
      if( ( fabs(mskw_45deg(0,0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mskw_45deg(0,1) + std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mskw_45deg(1,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mskw_45deg(1,1)) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Skew-Symmetric part of 45deg-rotation2D matrix is not correct!");
        return 1;
      };
      if( fabs((r_45deg * invert(r_45deg)).getAngle()) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Composition of a 45deg-rotation2D and an inverse 45deg-rotation2D does not have 0 angle!");
	return 1;
      };
      
      if(r_ident != r_ident) {
        RK_ERROR("r_ident is not equal to itself!");
	return 1;
      };
      if(r_45deg == r_ident) {
        RK_ERROR("r_45deg is equal to identity!");
	return 1;
      };

      mat<double,mat_structure::rectangular> m23(2,3);
      m23(0,0) = 1.0; m23(1,0) = 0.0;
      m23(0,1) = 0.0; m23(1,1) = 1.0;
      m23(0,2) = 1.0; m23(1,2) = 1.0;
      mat<double,mat_structure::rectangular> m23_r45(r_45deg * m23);
      if( ( fabs(m23_r45(0,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m23_r45(0,1) + std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m23_r45(1,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m23_r45(1,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m23_r45(0,2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m23_r45(1,2) - std::sqrt(2)) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("45deg-rotation2D matrix times a 2x3 matrix is not correct!");
        return 1;
      };
      mat<double,mat_structure::rectangular> m32(3,2);
      m32(0,0) = 1.0; m32(1,0) = 0.0; m32(2,0) = 1.0;
      m32(0,1) = 0.0; m32(1,1) = 1.0; m32(2,1) = 1.0;
      mat<double,mat_structure::rectangular> m32_r45(m32 * r_45deg);
      if( ( fabs(m32_r45(0,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m32_r45(0,1) + std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m32_r45(1,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m32_r45(1,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m32_r45(2,0) - std::sqrt(2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m32_r45(2,1)) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("A 3x2 matrix times a 45deg-rotation2D matrix is not correct!");
        return 1;
      };
      
      trans_mat_2D<double> t_ident;
      std::stringstream ss_t_ident; ss_t_ident << t_ident;
      if( ss_t_ident.str() != "(angle = 0; translation = (0; 0))" ) {
	RK_ERROR("The null-transformation2D stream-output is not correct!");
	return 1;
      };
      if( ( fabs(t_ident(0,0) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_ident(0,1)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_ident(0,2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_ident(1,0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_ident(1,1) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_ident(1,2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_ident(2,0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_ident(2,1)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_ident(2,2) - 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Null-transformation2D matrix does not have identity elements!");
        return 1;
      };
      trans_mat_2D<double> t_45deg_11(0.25f * double(M_PI),vect<double,2>(1.0,1.0));
      if( ( fabs(t_45deg_11(0,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_45deg_11(0,1) + std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_45deg_11(0,2) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_45deg_11(1,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_45deg_11(1,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_45deg_11(1,2) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_45deg_11(2,0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_45deg_11(2,1)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_45deg_11(2,2) - 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("(45deg,(1,1))-transformation2D does not have correct elements!");
        return 1;
      };
      mat<double,mat_structure::square> m_45deg_11 = t_45deg_11.getMat();
      if( ( fabs(m_45deg_11(0,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m_45deg_11(0,1) + std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m_45deg_11(0,2) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m_45deg_11(1,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m_45deg_11(1,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m_45deg_11(1,2) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m_45deg_11(2,0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m_45deg_11(2,1)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m_45deg_11(2,2) - 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("(45deg,(1,1))-transformation2D does not produce the correct regular matrix!");
        return 1;
      };
      trans_mat_2D<double> t_45deg_11_cpy(t_45deg_11);
      if( ( fabs(t_45deg_11_cpy.getAngle() - 0.25 * double(M_PI)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_45deg_11_cpy.getTranslation()[0] - 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("(45deg,(1,1))-transformation2D copy does not produce the correct angle or translation!");
	return 1;
      };
      t_45deg_11_cpy.setAngle(0.0);
      if( fabs(t_45deg_11_cpy.getAngle()) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("(0deg,(1,1))-transformation2D copy does not produce the correct angle!");
	return 1;
      };
      t_45deg_11_cpy.setRotMat(r_45deg);
      if( fabs(t_45deg_11_cpy.getAngle() - 0.25 * double(M_PI)) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("(45deg,(1,1))-transformation2D copy does not produce the correct angle!");
	return 1;
      };
      t_45deg_11_cpy.setTranslation(vect<double,2>(-1.0,-1.0));
      if( ( fabs(t_45deg_11_cpy.getTranslation()[0] + 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_45deg_11_cpy.getTranslation()[1] + 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("(45deg,(-1,-1))-transformation2D copy does not produce the correct translation!");
	return 1;
      };
      t_45deg_11_cpy.setTranslation(vect<double,2>(1.0,1.0));
      if( fabs(trace(t_45deg_11_cpy) - 1.0 - std::sqrt(2)) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Trace of (45deg,(1,1))-transformation2D is not correct!");
	return 1;
      };
      if( fabs(determinant(t_45deg_11_cpy) - 1.0) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Determinant of (45deg,(1,1))-transformation2D is not correct!");
	return 1;
      };
      const mat<double,mat_structure::symmetric> msym_45deg_11_cpy(t_45deg_11_cpy.getSymPart());
      if( ( fabs(msym_45deg_11_cpy(0,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(msym_45deg_11_cpy(0,1)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(msym_45deg_11_cpy(0,2) - 0.5) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(msym_45deg_11_cpy(1,0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(msym_45deg_11_cpy(1,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(msym_45deg_11_cpy(1,2) - 0.5) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(msym_45deg_11_cpy(2,0) - 0.5) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(msym_45deg_11_cpy(2,1) - 0.5) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(msym_45deg_11_cpy(2,2) - 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Symmetric part of (45deg,(1,1))-transformation2D is not correct!");
        return 1;
      };
      const mat<double,mat_structure::skew_symmetric> mskw_45deg_11_cpy(t_45deg_11_cpy.getSkewSymPart());
      if( ( fabs(mskw_45deg_11_cpy(0,0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mskw_45deg_11_cpy(0,1) + std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mskw_45deg_11_cpy(0,2) - 0.5) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mskw_45deg_11_cpy(1,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mskw_45deg_11_cpy(1,1)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mskw_45deg_11_cpy(1,2) - 0.5) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mskw_45deg_11_cpy(2,0) + 0.5) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mskw_45deg_11_cpy(2,1) + 0.5) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mskw_45deg_11_cpy(2,2)) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Skew-Symmetric part of (45deg,(1,1))-transformation2D is not correct! with value: " << mskw_45deg_11_cpy);
        return 1;
      };
      if( ( fabs(invert(t_45deg_11_cpy).getAngle() + 0.25 * double(M_PI)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(invert(t_45deg_11_cpy).getTranslation()[0] + std::sqrt(2)) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("(45deg,(1,1))-transformation2D copy does not produce the correct angle or translation!");
	return 1;
      };
      mat<double,mat_structure::square> mt_45deg_11_cpy(transpose(t_45deg_11_cpy));
      if( ( fabs(mt_45deg_11_cpy(0,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mt_45deg_11_cpy(0,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mt_45deg_11_cpy(0,2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mt_45deg_11_cpy(1,0) + std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mt_45deg_11_cpy(1,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mt_45deg_11_cpy(1,2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mt_45deg_11_cpy(2,0) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mt_45deg_11_cpy(2,1) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(mt_45deg_11_cpy(2,2) - 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("(45deg,(1,1))-transformation2D transposed does not produce the correct regular matrix!");
        return 1;
      };
      if( ( fabs((invert(t_45deg_11_cpy) * t_45deg_11_cpy).getAngle()) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs((invert(t_45deg_11_cpy) * t_45deg_11_cpy).getTranslation()[0]) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Composition of (45deg,(1,1))-transformation2D and its inverse does not produce a null-transformation2D!");
	return 1;
      };
      t_ident *= t_45deg_11_cpy;
      if( ( fabs(t_ident.getAngle() - 0.25 * double(M_PI)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_ident.getTranslation()[0] - 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Composition of identity and (45deg,(1,1))-transformation2D does not produce the correct transformation!");
	return 1;
      };
      t_ident *= invert(t_45deg_11_cpy);
      if( ( fabs(t_ident.getAngle()) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(t_ident.getTranslation()[0]) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Composition of identity, (45deg,(1,1))-transformation2D and its inverse does not produce the null-transformation!");
	return 1;
      };

      vect<double,2> v2 = t_45deg_11 * vect<double,2>(1.0,1.0);
      if( ( fabs(v2[0] - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(v2[1] - 1.0 - std::sqrt(2)) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Vector (1,1) transformed by (45deg,(1,1)) does not produce the correct vector!");
	return 1;
      };

      vect<double,3> v3 = t_45deg_11 * vect<double,3>(1.0,1.0,1.0);
      if( ( fabs(v3[0] - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(v3[1] - 1.0 - std::sqrt(2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(v3[2] - 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Vector (1,1) transformed by (45deg,(1,1)) does not produce the correct vector!");
	return 1;
      };
      v2 = t_45deg_11.rotate(vect<double,2>(1.0,1.0));
      if( ( fabs(v2[0]) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(v2[1] - std::sqrt(2)) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Vector (1,1) rotated by (45deg,(1,1)) does not produce the correct vector!");
	return 1;
      };

      mat<double,mat_structure::rectangular> m23_t(2,3);
      m23_t(0,0) = 1.0; m23_t(1,0) = 0.0;
      m23_t(0,1) = 0.0; m23_t(1,1) = 1.0;
      m23_t(0,2) = 1.0; m23_t(1,2) = 1.0;
      mat<double,mat_structure::rectangular> m23_r45_11(m23_t * t_45deg_11);
      if( ( fabs(m23_r45_11(0,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m23_r45_11(0,1) + std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m23_r45_11(0,2) - 2.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m23_r45_11(1,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m23_r45_11(1,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m23_r45_11(1,2) - 2.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("m23_t (rotated 45 degrees + (1,1)).transpose is not correct!");
        return 1;
      };
      mat<double,mat_structure::rectangular> m32_t(3,2);
      m32_t(0,0) = 1.0; m32_t(1,0) = 0.0; m32_t(2,0) = 1.0;
      m32_t(0,1) = 0.0; m32_t(1,1) = 1.0; m32_t(2,1) = 1.0;
      mat<double,mat_structure::rectangular> m32_r45_11(t_45deg_11 * m32_t);
      if( ( fabs(m32_r45_11(0,0) - std::sqrt(0.5) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m32_r45_11(0,1) + std::sqrt(0.5) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m32_r45_11(1,0) - std::sqrt(0.5) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m32_r45_11(1,1) - std::sqrt(0.5) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m32_r45_11(2,0) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(m32_r45_11(2,1) - 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("m32_t rotated 45 degrees + (1,1) is not correct!");
        return 1;
      };
      
      RK_NOTICE(2,"/!!!!!! TESTS OF ROTATION 2D PASSED !!!!!!/");
    };

    if(true){
    
      rot_mat_3D<double> r_ident;
      if( ( !is_diagonal(r_ident, std::numeric_limits<double>::epsilon()) ) ||
	  ( fabs( elem_norm_2(r_ident) - std::sqrt(3) ) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Null-rotation3D is not the identity matrix!");
        return 1;
      };
      double r_45z_a[] = {double(std::cos(0.25*M_PI)),double(std::sin(0.25*M_PI)),0.0,double(-std::sin(0.25*M_PI)),double(std::cos(0.25*M_PI)),0.0,0.0,0.0,1.0};
      rot_mat_3D<double> r_45z(r_45z_a);
      if( ( fabs(r_45z(0,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z(0,1) + std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z(0,2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z(1,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z(1,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z(1,2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z(2,0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z(2,1)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z(2,2) - 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("45degZ-rotation3D is not correct!");
        return 1;
      };
      rot_mat_3D<double> r_45z_cpy(r_45z);
      if( ( fabs(r_45z_cpy(0,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(0,1) + std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(0,2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(1,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(1,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(1,2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(2,0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(2,1)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(2,2) - 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("45degZ-rotation3D copy is not correct!");
        return 1;
      };
      r_45z_cpy = r_45z;
      if( ( fabs(r_45z_cpy(0,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(0,1) + std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(0,2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(1,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(1,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(1,2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(2,0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(2,1)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_45z_cpy(2,2) - 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("45degZ-rotation3D copy-assigned is not correct!");
        return 1;
      };
      r_ident *= r_45z;
      if( ( fabs(r_ident(0,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_ident(0,1) + std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_ident(0,2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_ident(1,0) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_ident(1,1) - std::sqrt(0.5)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_ident(1,2)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_ident(2,0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_ident(2,1)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(r_ident(2,2) - 1.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Null-rotation3D times 45degZ-rotation3D is not correct!");
        return 1;
      };
      r_ident *= invert(r_45z);
      if( ( !is_diagonal(r_ident, std::numeric_limits<double>::epsilon()) ) ||
	  ( fabs( elem_norm_2(r_ident) - std::sqrt(3) ) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Null-rotation3D (from composition of 45degZ-rotation3D and its inverse) is not the identity matrix!");
        return 1;
      };
      if( fabs(trace(r_45z_cpy) - std::sqrt(2) - 1.0) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Trace of the 45degZ-rotation3D is not correct!");
	return 1;
      };
      if( fabs(determinant(r_45z_cpy) - 1.0) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Determinant of the 45degZ-rotation3D is not correct!");
	return 1;
      };
      mat<double,mat_structure::symmetric> msym_45z(r_45z_cpy.getSymPart());
      if( ( !is_diagonal(msym_45z, std::numeric_limits<double>::epsilon()) ) ||
	  ( fabs( elem_norm_2(msym_45z) - std::sqrt(2) ) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Symmetric part of 45degZ-rotation3D is not correct!");
        return 1;
      };
      mat<double,mat_structure::skew_symmetric> mskw_45z(r_45z_cpy.getSkewSymPart());
      if( elem_norm_2(mskw_45z + transpose(mskw_45z)) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Skew-Symmetric part of 45degZ-rotation3D is not correct!");
        return 1;
      };
      if( elem_norm_2(r_45z_cpy.getMat() - r_45z_cpy) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("45degZ-rotation3D as a regular matrix is not the same!");
        return 1;
      };
      if( elem_norm_2(invert(r_45z_cpy).getMat() - transpose(r_45z_cpy)) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Transpose of 45degZ-rotation3D and the inverse are not the same!");
        return 1;
      };
      if( ( !is_diagonal((invert(r_45z_cpy) * r_45z_cpy), std::numeric_limits<double>::epsilon()) ) ||
	  ( fabs( elem_norm_2((invert(r_45z_cpy) * r_45z_cpy)) - std::sqrt(3) ) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Composition of 45degZ-rotation3D and its inverse is not the identity matrix!");
        return 1;
      };
    
      if( !(r_ident == r_ident) ) {
        RK_ERROR("Rotation identity is not equal to itself!");
	return 1;
      };
      if( !(r_ident != r_45z) ) {
        RK_ERROR("Rotation identity is equal to 45 deg rotation about z!");
	return 1;
      };
      vect<double,3> v1(1.0,1.0,2.0);
      if( norm( (r_45z * v1) - vect<double,3>(0.0,std::sqrt(2.0),2.0) ) > 1.15*std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("r_45z * v(1,1,2) is not correct!");
	return 1;
      };
      if( norm( (v1 * r_45z) - vect<double,3>(std::sqrt(2.0),0.0,2.0) ) > 1.15*std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("v(1,1,2) * r_45z is not correct!");
	return 1;
      };
    
      quaternion<double> q_45z(r_45z);
      if( ( fabs( q_45z[0] - std::cos(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_45z[3] - std::sin(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("q_45z does not have the correct elements!");
	return 1;
      };
      if( ( !is_diagonal((r_45z * invert(q_45z)), std::numeric_limits<double>::epsilon()) ) ||
	  ( fabs( elem_norm_2(r_45z * invert(q_45z)) - std::sqrt(3) ) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Composition of 45degZ-rotation3D and the inverse of 45degZ-quaternion is not the identity matrix!");
        return 1;
      };
      quaternion<double> q_ident;
      if( ( fabs( q_ident[0] - 1.0) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[3] ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("q_ident does not have the correct elements!");
	return 1;
      };
      quaternion<double> q_45z_cpy(q_45z);
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
      q_ident *= invert(q_45z);
      if( ( fabs( q_ident[0] - std::cos(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_ident[3] - std::sin(0.125 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("q_45z * q_45z * invert(q_45z) does not have the correct elements!");
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
      mat<double,mat_structure::square> rm_90z(q_45z * r_45z);
      if( ( fabs( rm_90z(0,0) ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( rm_90z(0,1) + 1.0) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( rm_90z(0,2) ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( rm_90z(1,0) - 1.0) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( rm_90z(1,1) ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( rm_90z(1,2) ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( rm_90z(2,0) ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( rm_90z(2,1) ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( rm_90z(2,2) - 1.0 ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("q_45z * r_45z does not have the correct elements!");
	return 1;
      };
    
      if( fabs( trace(q_45z_cpy) - 1.0 - std::sqrt(2) ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Trace of q_45z_cpy is not correct!");
	return 1;
      };
      if( fabs( determinant(q_45z_cpy) - 1.0 ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Determinant of q_45z_cpy is not correct!");
	return 1;
      };
      if( elem_norm_2(q_45z_cpy.getSymPart() - r_45z.getSymPart()) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Symmetric part of q_45z_cpy is not correct!");
	return 1;
      };
      if( elem_norm_2(q_45z_cpy.getSkewSymPart() - r_45z.getSkewSymPart()) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Skew-Symmetric part of q_45z_cpy is not correct!");
	return 1;
      };
      if( fabs( (q_45z_cpy * invert(q_45z_cpy))[0] - 1.0 ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("q_45z_cpy * invert(q_45z_cpy) is not correct!");
	return 1;
      };
      if( fabs( (q_45z_cpy * transpose(q_45z_cpy))[0] - 1.0 ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("q_45z_cpy * transpose(q_45z_cpy) is not correct!");
	return 1;
      };
      
      if( !(q_ident == q_ident) ) {
        RK_ERROR("q_ident is not equal to itself!");
	return 1;
      };
      if( !(q_ident != q_45z) ) {
        RK_ERROR("q_ident is equal to q_45z!");
        return 1;
      };
      if( norm( (q_45z * v1) - vect<double,3>(0.0,std::sqrt(2.0),2.0) ) > 1.5 * std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("q_45z * v(1,1,2) is not correct!");
	return 1;
      };

      euler_angles_TB<double> e_45z(0.25*M_PI,0.0,0.0);
      quaternion<double> q_e_90z(q_45z * e_45z);
      if( ( fabs( q_e_90z[0] - std::cos(0.25 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_e_90z[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_e_90z[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_e_90z[3] - std::sin(0.25 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("q_45z * e_45z does not have the correct elements!");
	return 1;
      };
      axis_angle<double> a_45z(0.25*M_PI,vect<double,3>(0.0,0.0,1.0));
      quaternion<double> q_a_90z(q_45z * a_45z);
      if( ( fabs( q_a_90z[0] - std::cos(0.25 * M_PI)) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_a_90z[1] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_a_90z[2] ) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs( q_a_90z[3] - std::sin(0.25 * M_PI)) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("q_45z * a_45z does not have the correct elements!");
	return 1;
      };
    
      euler_angles_TB<double> e_ident;
      if( ( fabs(e_ident.yaw()) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_ident.pitch()) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_ident.roll()) > std::numeric_limits<double>::epsilon() ) ) {
        RK_ERROR("e_ident does not have correct elements!");
        return 1;
      };
      euler_angles_TB<double> e_45z_cpy(e_45z);
      if( ( fabs(e_45z_cpy.yaw() - 0.25 * M_PI) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_45z_cpy.pitch()) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_45z_cpy.roll()) > std::numeric_limits<double>::epsilon() ) ) {
        RK_ERROR("e_45z copy does not have correct elements!");
        return 1;
      };
      euler_angles_TB<double> e_45z_r(r_45z);
      if( ( fabs(e_45z_cpy.yaw() - 0.25 * M_PI) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_45z_cpy.pitch()) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_45z_cpy.roll()) > std::numeric_limits<double>::epsilon() ) ) {
        RK_ERROR("e_45z created from r_45z does not have correct elements!");
        return 1;
      };
      euler_angles_TB<double> e_45z_q(q_45z);
      if( ( fabs(e_45z_cpy.yaw() - 0.25 * M_PI) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_45z_cpy.pitch()) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_45z_cpy.roll()) > std::numeric_limits<double>::epsilon() ) ) {
        RK_ERROR("e_45z created from q_45z does not have correct elements!");
        return 1;
      };
      if( elem_norm_2(e_45z_q.getRotMat().getMat() - r_45z.getMat()) > std::numeric_limits<double>::epsilon() ) {
        RK_ERROR("e_45z rotation matrix is not correct!");
        return 1;
      };
      if( elem_norm_2(e_45z_q.getMat() - r_45z.getMat()) > std::numeric_limits<double>::epsilon() ) {
        RK_ERROR("e_45z as matrix is not correct!");
        return 1;
      };
      e_45z_cpy = e_45z;
      if( ( fabs(e_45z_cpy.yaw() - 0.25 * M_PI) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_45z_cpy.pitch()) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_45z_cpy.roll()) > std::numeric_limits<double>::epsilon() ) ) {
        RK_ERROR("e_45z copy-assigned to e_45z does not have correct elements!");
        return 1;
      };
      e_45z_cpy = r_45z;
      if( ( fabs(e_45z_cpy.yaw() - 0.25 * M_PI) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_45z_cpy.pitch()) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_45z_cpy.roll()) > std::numeric_limits<double>::epsilon() ) ) {
        RK_ERROR("e_45z copy-assigned to r_45z does not have correct elements!");
        return 1;
      };
      e_45z_cpy = q_45z;
      if( ( fabs(e_45z_cpy.yaw() - 0.25 * M_PI) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_45z_cpy.pitch()) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_45z_cpy.roll()) > std::numeric_limits<double>::epsilon() ) ) {
        RK_ERROR("e_45z copy-assigned to q_45z does not have correct elements!");
        return 1;
      };
      e_45z_cpy = a_45z;
      if( ( fabs(e_45z_cpy.yaw() - 0.25 * M_PI) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_45z_cpy.pitch()) > std::numeric_limits<double>::epsilon() ) ||
          ( fabs(e_45z_cpy.roll()) > std::numeric_limits<double>::epsilon() ) ) {
        RK_ERROR("e_45z copy-assigned to a_45z does not have correct elements!");
        return 1;
      };

      if( !(e_ident == e_ident) ) {
        RK_ERROR("e_ident is not equal to itself!");
        return 1;
      };
      if( !(e_ident != e_45z) ) {
        RK_ERROR("e_ident is equal to e_45z!");
        return 1;
      };

      if( fabs( trace(e_45z_cpy) - 1.0 - std::sqrt(2) ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Trace of e_45z_cpy is not correct!");
	return 1;
      };
      if( fabs( determinant(e_45z_cpy) - 1.0 ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Determinant of e_45z_cpy is not correct!");
	return 1;
      };
      if( elem_norm_2(e_45z_cpy.getSymPart() - r_45z.getSymPart()) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Symmetric part of e_45z_cpy is not correct!");
	return 1;
      };
      if( elem_norm_2(e_45z_cpy.getSkewSymPart() - r_45z.getSkewSymPart()) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Skew-Symmetric part of e_45z_cpy is not correct!");
	return 1;
      };
      if( fabs( (q_45z_cpy * invert(e_45z_cpy))[0] - 1.0 ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("q_45z_cpy * invert(e_45z_cpy) is not correct!");
	return 1;
      };
      if( fabs( (q_45z_cpy * transpose(e_45z_cpy))[0] - 1.0 ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("q_45z_cpy * transpose(e_45z_cpy) is not correct!");
	return 1;
      };
      
      axis_angle<double> a_ident;
      if( fabs( a_ident.angle() ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("a_ident does not have the correct angle!");
	return 1;
      };
      axis_angle<double> a_45z_cpy(a_45z);
      if( ( fabs( a_45z_cpy.angle() - 0.25 * M_PI ) > std::numeric_limits<double>::epsilon() ) ||
          ( norm( a_45z_cpy.axis() - vect<double,3>(0.0,0.0,1.0) ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("a_45z copy does not have the correct angle-axis!");
	return 1;
      };
      axis_angle<double> a_45z_r(r_45z);
      if( ( fabs( a_45z_r.angle() - 0.25 * M_PI ) > std::numeric_limits<double>::epsilon() ) ||
          ( norm( a_45z_r.axis() - vect<double,3>(0.0,0.0,1.0) ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("a_45z created from r_45z does not have the correct angle-axis!");
	return 1;
      };
      axis_angle<double> a_45z_q(q_45z);
      if( ( fabs( a_45z_q.angle() - 0.25 * M_PI ) > std::numeric_limits<double>::epsilon() ) ||
          ( norm( a_45z_q.axis() - vect<double,3>(0.0,0.0,1.0) ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("a_45z created from q_45z does not have the correct angle-axis!");
	return 1;
      };
      axis_angle<double> a_45z_e(e_45z);
      if( ( fabs( a_45z_e.angle() - 0.25 * M_PI ) > std::numeric_limits<double>::epsilon() ) ||
          ( norm( a_45z_e.axis() - vect<double,3>(0.0,0.0,1.0) ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("a_45z created from e_45z does not have the correct angle-axis!");
	return 1;
      };
      if( elem_norm_2(a_45z_q.getRotMat().getMat() - r_45z.getMat()) > std::numeric_limits<double>::epsilon() ) {
        RK_ERROR("a_45z rotation matrix is not correct!");
        return 1;
      };
      if( elem_norm_2(a_45z_q.getMat() - r_45z.getMat()) > std::numeric_limits<double>::epsilon() ) {
        RK_ERROR("a_45z as matrix is not correct!");
        return 1;
      };
      a_45z_cpy = a_45z;
      if( ( fabs( a_45z_cpy.angle() - 0.25 * M_PI ) > std::numeric_limits<double>::epsilon() ) ||
          ( norm( a_45z_cpy.axis() - vect<double,3>(0.0,0.0,1.0) ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("a_45z copy-assigned does not have the correct angle-axis!");
	return 1;
      };
      a_45z_cpy = r_45z;
      if( ( fabs( a_45z_cpy.angle() - 0.25 * M_PI ) > std::numeric_limits<double>::epsilon() ) ||
          ( norm( a_45z_cpy.axis() - vect<double,3>(0.0,0.0,1.0) ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("a_45z copy-assigned from r_45z does not have the correct angle-axis!");
	return 1;
      };
      a_45z_cpy = q_45z;
      if( ( fabs( a_45z_cpy.angle() - 0.25 * M_PI ) > std::numeric_limits<double>::epsilon() ) ||
          ( norm( a_45z_cpy.axis() - vect<double,3>(0.0,0.0,1.0) ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("a_45z copy-assigned from q_45z does not have the correct angle-axis!");
	return 1;
      };
      a_45z_cpy = e_45z;
      if( ( fabs( a_45z_cpy.angle() - 0.25 * M_PI ) > std::numeric_limits<double>::epsilon() ) ||
          ( norm( a_45z_cpy.axis() - vect<double,3>(0.0,0.0,1.0) ) > std::numeric_limits<double>::epsilon() ) ){
	RK_ERROR("a_45z copy-assigned from e_45z does not have the correct angle-axis!");
	return 1;
      };
      
      if( !(a_ident == a_ident) ) {
        RK_ERROR("a_ident is not equal to itself!");
	return 1;
      };
      if( !(a_ident != a_45z) ) {
        RK_ERROR("a_ident is equal to a_45z!");
	return 1;
      };
      
      if( fabs( trace(a_45z_cpy) - 1.0 - std::sqrt(2) ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Trace of a_45z_cpy is not correct!");
	return 1;
      };
      if( fabs( determinant(a_45z_cpy) - 1.0 ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Determinant of a_45z_cpy is not correct!");
	return 1;
      };
      if( elem_norm_2(a_45z_cpy.getSymPart() - r_45z.getSymPart()) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Symmetric part of a_45z_cpy is not correct!");
	return 1;
      };
      if( elem_norm_2(a_45z_cpy.getSkewSymPart() - r_45z.getSkewSymPart()) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("Skew-Symmetric part of a_45z_cpy is not correct!");
	return 1;
      };
      if( fabs( (q_45z_cpy * invert(a_45z_cpy))[0] - 1.0 ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("q_45z_cpy * invert(a_45z_cpy) is not correct!");
	return 1;
      };
      if( fabs( (q_45z_cpy * transpose(a_45z_cpy))[0] - 1.0 ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("q_45z_cpy * transpose(a_45z_cpy) is not correct!");
	return 1;
      };


      trans_mat_3D<double> t_ident;
      if( ( !is_diagonal(t_ident, std::numeric_limits<double>::epsilon()) ) ||
	  ( fabs( elem_norm_2(t_ident) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Null-transformation3D is not the identity matrix!");
        return 1;
      };
      double t_45z_array[] = {std::cos(0.25*M_PI),std::sin(0.25*M_PI),0.0,0.0,-std::sin(0.25*M_PI),std::cos(0.25*M_PI),0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0};
      trans_mat_3D<double> t_45z(t_45z_array);
      if( ( fabs( elem_norm_2(t_45z) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from array is not correct!");
        return 1;
      };
      trans_mat_3D<double> t_45z_cpy(t_45z);
      if( ( fabs( elem_norm_2(t_45z_cpy) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_cpy(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_cpy(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from copy is not correct!");
        return 1;
      };
      trans_mat_3D<double> t_45z_r(r_45z);
      if( ( fabs( elem_norm_2(t_45z_r) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_r(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_r(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_r(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_r(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from rotation3D is not correct!");
        return 1;
      };
      trans_mat_3D<double> t_45z_q(q_45z);
      if( ( fabs( elem_norm_2(t_45z_q) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_q(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_q(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_q(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_q(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from quaternion is not correct!");
        return 1;
      };
      trans_mat_3D<double> t_45z_e(e_45z);
      if( ( fabs( elem_norm_2(t_45z_e) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_e(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_e(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_e(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_e(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from Euler-angles is not correct!");
        return 1;
      };
      trans_mat_3D<double> t_45z_a(a_45z);
      if( ( fabs( elem_norm_2(t_45z_a) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_a(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_a(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_a(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_a(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from axis-angle is not correct!");
        return 1;
      };
      rot_mat_3D<double> t_45z_q_rot = t_45z_q.getRotMat();
      if( ( fabs( elem_norm_2(t_45z_q_rot) - std::sqrt(3.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_q_rot(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_q_rot(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_q_rot(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("Rotation matrix from 45degz-transformation3D is not correct!");
        return 1;
      };
      mat<double,mat_structure::square> t_45z_q_mat = t_45z_q.getMat();
      if( ( fabs( elem_norm_2(t_45z_q_mat) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_q_mat(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_q_mat(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_q_mat(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_q_mat(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("Regular matrix from 45degz-transformation3D is not correct!");
        return 1;
      };
      
      t_45z_cpy = t_45z;
      if( ( fabs( elem_norm_2(t_45z_cpy) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_cpy(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_cpy(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from assignment is not correct!");
        return 1;
      };
      t_45z_cpy = r_45z;
      if( ( fabs( elem_norm_2(t_45z_cpy) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_cpy(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_cpy(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from assignment to rotation3D is not correct!");
        return 1;
      };
      t_45z_cpy = q_45z;
      if( ( fabs( elem_norm_2(t_45z_cpy) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_cpy(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_cpy(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from assignment to quaternion is not correct!");
        return 1;
      };
      t_45z_cpy = e_45z;
      if( ( fabs( elem_norm_2(t_45z_cpy) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_cpy(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_cpy(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from assignment to Euler-angles is not correct!");
        return 1;
      };
      t_45z_cpy = a_45z;
      if( ( fabs( elem_norm_2(t_45z_cpy) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_cpy(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_cpy(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from assignment to axis-angle is not correct!");
        return 1;
      };
      
      t_45z_cpy = t_ident;
      t_45z_cpy *= t_45z;
      if( ( fabs( elem_norm_2(t_45z_cpy) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_cpy(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_cpy(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from identity times 45degz-transformation3D is not correct!");
        return 1;
      };
      t_45z_cpy = t_ident;
      t_45z_cpy *= r_45z;
      if( ( fabs( elem_norm_2(t_45z_cpy) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_cpy(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_cpy(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from identity times 45degz-rotation3D is not correct!");
        return 1;
      };
      t_45z_cpy = t_ident;
      t_45z_cpy *= q_45z;
      if( ( fabs( elem_norm_2(t_45z_cpy) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_cpy(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_cpy(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from identity times 45degz-quaternion is not correct!");
        return 1;
      };
      t_45z_cpy = t_ident;
      t_45z_cpy *= e_45z;
      if( ( fabs( elem_norm_2(t_45z_cpy) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_cpy(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_cpy(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from identity times 45degz-euler-angles is not correct!");
        return 1;
      };
      t_45z_cpy = t_ident;
      t_45z_cpy *= a_45z;
      if( ( fabs( elem_norm_2(t_45z_cpy) - std::sqrt(4.0) ) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_cpy(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_cpy(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_cpy(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-transformation3D from identity times 45degz-axis-angle is not correct!");
        return 1;
      };
      
      
      if(!(t_ident == t_ident)) {
        RK_ERROR("t_ident is not equal to itself!");
	return 1;
      };
      if(!(t_ident != t_45z)) {
        RK_ERROR("t_ident is equal to t_45z!");
	return 1;
      };

      trans_mat_3D<double> t_45z_123(r_45z,vect<double,3>(1.0,2.0,3.0));
      if( ( fabs( elem_norm_2(t_45z_123) - std::sqrt(18.0) ) > 1.8*std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_123(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_123(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_123(0,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_123(1,3) - 2.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_123(2,3) - 3.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_123(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("45degz-(1,2,3)t-transformation3D is not correct!");
        return 1;
      };

      if( fabs( trace(t_45z_123) - 2.0 - std::sqrt(2.0) ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("The trace of the 45degz-(1,2,3)t-transformation3D is not correct!");
        return 1;
      };
      if( fabs( determinant(t_45z_123) - 1.0 ) > std::numeric_limits<double>::epsilon() ) {
	RK_ERROR("The determinant of the 45degz-(1,2,3)t-transformation3D is not correct!");
        return 1;
      };
      mat<double, mat_structure::square> t_45z_mat_sym_skw(t_45z_123.getSymPart() + t_45z_123.getSkewSymPart());
      if( ( fabs( elem_norm_2(t_45z_mat_sym_skw) - std::sqrt(18.0) ) > 1.8*std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_mat_sym_skw(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_mat_sym_skw(0,1) + std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_mat_sym_skw(0,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_mat_sym_skw(1,3) - 2.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_mat_sym_skw(2,3) - 3.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_mat_sym_skw(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("The skew-symmetric and symmetric parts of 45degz-(1,2,3)t-transformation3D don't add up to the correct matrix!");
        return 1;
      };
      t_ident = t_45z_123 * invert(t_45z_123);
      if( ( fabs( elem_norm_2(t_ident) - std::sqrt(4.0) ) > 1.8*std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_ident(0,0) - 1.0 ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_ident(0,1) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(1,1) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(0,3) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(1,3) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(2,3) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("The compositionof 45degz-(1,2,3)t-transformation3D with its inverse doesn't give the identity matrix!");
        return 1;
      };
      t_ident = (invert(t_45z_123).getMat() * t_45z_123);
      if( ( fabs( elem_norm_2(t_ident) - std::sqrt(4.0) ) > 1.8*std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_ident(0,0) - 1.0 ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_ident(0,1) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(1,1) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(0,3) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(1,3) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(2,3) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("The composition of 45degz-(1,2,3)t-transformation3D with its inverse as a regular matrix doesn't give the identity matrix!");
        return 1;
      };
      t_ident = (invert(t_45z_123) * t_45z_123.getMat());
      if( ( fabs( elem_norm_2(t_ident) - std::sqrt(4.0) ) > 1.8*std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_ident(0,0) - 1.0 ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_ident(0,1) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(1,1) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(0,3) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(1,3) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(2,3) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("The composition of 45degz-(1,2,3)t-transformation3D as a regular matrix with its inverse doesn't give the identity matrix!");
        return 1;
      };
      
      mat<double, mat_structure::square> t_45z_123_t = transpose(t_45z_123);
      if( ( fabs( elem_norm_2(t_45z_123_t) - std::sqrt(18.0) ) > 1.8*std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_45z_123_t(0,0) - std::cos(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_45z_123_t(0,1) - std::sin(0.25*M_PI) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_123_t(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_123_t(0,3) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_123_t(3,1) - 2.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_123_t(3,2) - 3.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_45z_123_t(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("The transpose of 45degz-(1,2,3)t-transformation3D doesn't give the correct matrix!");
        return 1;
      };
      
      t_ident = (invert(t_45z_123) * r_45z) * trans_mat_3D<double>(r_ident,-invert(t_45z_123).getTranslation());
      if( ( fabs( elem_norm_2(t_ident) - std::sqrt(4.0) ) > 1.8*std::numeric_limits<double>::epsilon() ) ||
	  ( fabs( t_ident(0,0) - 1.0 ) > std::numeric_limits< double >::epsilon() )  ||
	  ( fabs( t_ident(0,1) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(1,1) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(2,2) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(0,3) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(1,3) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(2,3) ) > std::numeric_limits< double >::epsilon() ) ||
	  ( fabs( t_ident(3,3) - 1.0 ) > std::numeric_limits< double >::epsilon() ) ) {
	RK_ERROR("The composition of inverse 45degz-(1,2,3)t-transformation3D with a 45degz-rotation3D and the inverse translation doesn't give the identity matrix!");
        return 1;
      };
      
      vect<double,3> v1_trans(t_45z_123 * v1);
      if( ( fabs(v1_trans[0] - 1.0) > std::numeric_limits<double>::epsilon() ) || 
	  ( fabs(v1_trans[1] - 2.0 - std::sqrt(2.0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(v1_trans[2] - 5.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("t_45z_123 * v(1,1,2) is not correct!");
        return 1;
      };
      v1_trans = t_45z_123.rotate(v1);
      if( ( fabs(v1_trans[0]) > std::numeric_limits<double>::epsilon() ) || 
	  ( fabs(v1_trans[1] - std::sqrt(2.0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(v1_trans[2] - 2.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("t_45z_123.rotate( v(1,1,2) ) is not correct!");
        return 1;
      };
      vect<double,4> v3 =  t_45z_123 * vect<double,4>(1,1,2,2);
      
      if( ( fabs(v3[0] - 2.0) > std::numeric_limits<double>::epsilon() ) || 
	  ( fabs(v3[1] - 4.0 - std::sqrt(2.0)) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(v3[2] - 8.0) > std::numeric_limits<double>::epsilon() ) ||
	  ( fabs(v3[3] - 2.0) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("t_45z_123 * v(1,1,2,2) is not correct!");
        return 1;
      };
      

      axis_angle<double> a_weird(0.3241,vect<double,3>(0.5,0.5,sqrt(0.5)));
      quaternion<double> q_weird;
      q_weird = a_weird;
      euler_angles_TB<double> e_weird;
      e_weird = q_weird;
      rot_mat_3D<double> r_weird;
      r_weird = q_weird;
      axis_angle<double> a_weird_q(q_weird);
      axis_angle<double> a_weird_e(e_weird);
      axis_angle<double> a_weird_r(r_weird);
      if( ( fabs(a_weird_q.angle() - a_weird.angle()) > std::numeric_limits<double>::epsilon() ) ||
	  ( norm(a_weird_q.axis() - a_weird.axis()) > std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Weird rotation does not translate correctly from axis-angle to quaternion and back!");
        return 1;
      };
      if( ( fabs(a_weird_e.angle() - a_weird.angle()) > 20.0*std::numeric_limits<double>::epsilon() ) ||
	  ( norm(a_weird_e.axis() - a_weird.axis()) > 20.0*std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Weird rotation does not translate correctly from axis-angle to Euler-angles and back!");
        return 1;
      };
      if( ( fabs(a_weird_r.angle() - a_weird.angle()) > 20.0*std::numeric_limits<double>::epsilon() ) ||
	  ( norm(a_weird_r.axis() - a_weird.axis()) > 20.0*std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Weird rotation does not translate correctly from axis-angle to rotation3D and back!");
        return 1;
      };
      
      quaternion<double> q_res(q_45z * quaternion<double>(a_weird * a_45z));
      vect<double,4> v_a(q_res[0],q_res[1],q_res[2],q_res[3]);
      q_res = q_45z * (q_weird * q_45z);
      vect<double,4> v_q(q_res[0],q_res[1],q_res[2],q_res[3]);
      q_res = quaternion<double>(q_45z * e_weird * r_45z);
      vect<double,4> v_e(q_res[0],q_res[1],q_res[2],q_res[3]);
      q_res = quaternion<double>(a_45z * r_weird * e_45z);
      vect<double,4> v_r(q_res[0],q_res[1],q_res[2],q_res[3]);
      if( ( norm( v_a - v_q ) > std::numeric_limits<double>::epsilon() ) ||
          ( norm( v_a - v_e ) > 1.25 * std::numeric_limits<double>::epsilon() ) ||
          ( norm( v_a - v_r ) > std::numeric_limits<double>::epsilon() ) ||
          ( norm( v_q - v_e ) > std::numeric_limits<double>::epsilon() ) ||
          ( norm( v_q - v_r ) > std::numeric_limits<double>::epsilon() ) ||
          ( norm( v_e - v_r ) > 1.25 * std::numeric_limits<double>::epsilon() ) ) {
	RK_ERROR("Inter-operability tests between axis-angle, euler-angles, quaternions and rotation3D have failed!");
        return 1;
      };

      RK_NOTICE(2,"/!!!!!! TESTS OF ROTATION 3D PASSED !!!!!!/");
    };
  
  } catch(std::exception& e) {
    RK_ERROR("An exception has occurred during the math_gen test: '" << e.what() << "'");
  } catch(...) {
    RK_ERROR("An unexpected and unidentified exception has occurred during the math_gen test.");
  };
  
  return 0;
};









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
#include "ReaK/math/kinetostatics/rotations.h"

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_matchers.h"
#include "ReaK/math/lin_alg/mat_norms.h"
#include "ReaK/math/lin_alg/vect_matchers.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <numbers>

namespace ReaK {
namespace {

using ::ReaK::testing::MatrixIsDiagonal;
using ::ReaK::testing::MatrixIsIdentity;
using ::ReaK::testing::MatrixIsNear;
using ::ReaK::testing::VectorIsNear;

const double rel_tol = 10.0 * std::numeric_limits<double>::epsilon();

auto R452D() {
  return rot_mat_2D<double>(0.25 * std::numbers::pi);
}

auto T452D11() {
  return trans_mat_2D<double>(0.25 * std::numbers::pi, vect<double, 2>(1.0, 1.0));
}

auto R45Z() {
  std::array<double, 9> r_45z_a = {double(std::cos(0.25 * std::numbers::pi)),
                                   double(std::sin(0.25 * std::numbers::pi)),
                                   0.0,
                                   double(-std::sin(0.25 * std::numbers::pi)),
                                   double(std::cos(0.25 * std::numbers::pi)),
                                   0.0,
                                   0.0,
                                   0.0,
                                   1.0};
  return rot_mat_3D<double>(r_45z_a.data());
}

auto Q45Z() {
  return quaternion<double>(R45Z());
}

auto E45Z() {
  return euler_angles_TB<double>(0.25 * std::numbers::pi, 0.0, 0.0);
}

auto A45Z() {
  return axis_angle<double>(0.25 * std::numbers::pi, vect<double, 3>(0.0, 0.0, 1.0));
}

auto T45Z() {
  std::array<double, 16> t_45z_array = {std::cos(0.25 * std::numbers::pi),
                                        std::sin(0.25 * std::numbers::pi),
                                        0.0,
                                        0.0,
                                        -std::sin(0.25 * std::numbers::pi),
                                        std::cos(0.25 * std::numbers::pi),
                                        0.0,
                                        0.0,
                                        0.0,
                                        0.0,
                                        1.0,
                                        0.0,
                                        0.0,
                                        0.0,
                                        0.0,
                                        1.0};
  return trans_mat_3D<double>(t_45z_array.data());
}

auto T45Z123() {
  return trans_mat_3D<double>(R45Z(), vect<double, 3>(1.0, 2.0, 3.0));
}

TEST(Rotation2d, Identity) {
  rot_mat_2D<double> r_ident;
  std::stringstream ss_r_ident;
  ss_r_ident << r_ident;
  EXPECT_EQ(ss_r_ident.str(), "(angle = 0)");
  EXPECT_NEAR(r_ident.getAngle(), 0.0, rel_tol);
  EXPECT_NEAR(r_ident(0, 0), 1.0, rel_tol);
  EXPECT_NEAR(r_ident(0, 1), 0.0, rel_tol);
  EXPECT_NEAR(r_ident(1, 0), 0.0, rel_tol);
  EXPECT_NEAR(r_ident(1, 1), 1.0, rel_tol);
  EXPECT_EQ(r_ident, r_ident);
}

TEST(Rotation2d, GettersAndSetters) {
  rot_mat_2D<double> r_45deg(R452D());
  EXPECT_NEAR(r_45deg.getAngle(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_NEAR(r_45deg(0, 0), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(r_45deg(0, 1), -std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(r_45deg(1, 0), std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(r_45deg(1, 1), std::sqrt(0.5), rel_tol);
  rot_mat_2D<double> r_45deg_cpy(r_45deg);
  EXPECT_NEAR(r_45deg_cpy.getAngle(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_NEAR(trace(r_45deg_cpy), std::sqrt(2.0), rel_tol);
  EXPECT_NEAR(determinant(r_45deg_cpy), 1.0, rel_tol);
  r_45deg_cpy.setAngle(0.0);
  EXPECT_NEAR(r_45deg_cpy.getAngle(), 0.0, rel_tol);
  EXPECT_NE(r_45deg, rot_mat_2D<double>());
}

TEST(Rotation2d, MatrixConversion) {
  mat<double, mat_structure::square> m_45deg = R452D().getMat();
  EXPECT_NEAR(m_45deg(0, 0), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(m_45deg(0, 1), -std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(m_45deg(1, 0), std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(m_45deg(1, 1), std::sqrt(0.5), rel_tol);
}

TEST(Rotation2d, Products) {
  rot_mat_2D<double> r_45deg(R452D());
  rot_mat_2D<double> r_90deg = r_45deg * r_45deg;
  EXPECT_NEAR(r_90deg.getAngle(), 0.5 * std::numbers::pi, 10.0 * rel_tol);
  r_90deg *= r_45deg;
  EXPECT_NEAR(r_90deg.getAngle(), 0.75 * std::numbers::pi, rel_tol);
  vect<double, 2> v1 = r_45deg * vect<double, 2>(1.0, 1.0);
  EXPECT_NEAR(v1[0], 0.0, rel_tol);
  EXPECT_NEAR(v1[1], std::sqrt(2.0), 10.0 * rel_tol);
  v1 = vect<double, 2>(1.0, 1.0) * r_45deg;
  EXPECT_NEAR(v1[1], 0.0, rel_tol);
  EXPECT_NEAR(v1[0], std::sqrt(2.0), 10.0 * rel_tol);
  r_90deg *= invert(r_45deg);
  EXPECT_NEAR(r_90deg.getAngle(), 0.5 * std::numbers::pi, 10.0 * rel_tol);
  EXPECT_NEAR((r_45deg * invert(r_45deg)).getAngle(), 0.0, rel_tol);
}

TEST(Rotation2d, SymMatrixConversion) {
  const mat<double, mat_structure::symmetric> msym_45deg = R452D().getSymPart();
  EXPECT_NEAR(msym_45deg(0, 0), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(msym_45deg(0, 1), 0.0, rel_tol);
  EXPECT_NEAR(msym_45deg(1, 0), 0.0, rel_tol);
  EXPECT_NEAR(msym_45deg(1, 1), std::sqrt(0.5), rel_tol);
}

TEST(Rotation2d, SkewSymMatrixConversion) {
  const mat<double, mat_structure::skew_symmetric> mskw_45deg =
      R452D().getSkewSymPart();
  EXPECT_NEAR(mskw_45deg(0, 0), 0.0, rel_tol);
  EXPECT_NEAR(mskw_45deg(0, 1), -std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(mskw_45deg(1, 0), std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(mskw_45deg(1, 1), 0.0, rel_tol);
}

TEST(Rotation2d, MatrixProducts) {
  mat<double, mat_structure::rectangular> m23(2, 3);
  m23(0, 0) = 1.0;
  m23(1, 0) = 0.0;
  m23(0, 1) = 0.0;
  m23(1, 1) = 1.0;
  m23(0, 2) = 1.0;
  m23(1, 2) = 1.0;
  mat<double, mat_structure::rectangular> m23_r45(R452D() * m23);
  EXPECT_NEAR(m23_r45(0, 0), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(m23_r45(0, 1), -std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(m23_r45(1, 0), std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(m23_r45(1, 1), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(m23_r45(0, 2), 0.0, rel_tol);
  EXPECT_NEAR(m23_r45(1, 2), std::sqrt(2.0), 10.0 * rel_tol);

  mat<double, mat_structure::rectangular> m32(3, 2);
  m32(0, 0) = 1.0;
  m32(1, 0) = 0.0;
  m32(2, 0) = 1.0;
  m32(0, 1) = 0.0;
  m32(1, 1) = 1.0;
  m32(2, 1) = 1.0;
  mat<double, mat_structure::rectangular> m32_r45(m32 * R452D());
  EXPECT_NEAR(m32_r45(0, 0), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(m32_r45(0, 1), -std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(m32_r45(1, 0), std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(m32_r45(1, 1), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(m32_r45(2, 1), 0.0, rel_tol);
  EXPECT_NEAR(m32_r45(2, 0), std::sqrt(2.0), 10.0 * rel_tol);
}

TEST(Transform2d, Identity) {
  trans_mat_2D<double> t_ident;
  std::stringstream ss_t_ident;
  ss_t_ident << t_ident;
  EXPECT_EQ(ss_t_ident.str(), "(angle = 0; translation = (0; 0))");
  EXPECT_THAT(t_ident, MatrixIsIdentity(rel_tol));
  EXPECT_EQ(t_ident, t_ident);
}

TEST(Transform2d, GettersAndSetters) {
  trans_mat_2D<double> t_45deg_11(T452D11());
  EXPECT_NEAR(t_45deg_11(0, 0), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(t_45deg_11(0, 1), -std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(t_45deg_11(0, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45deg_11(1, 0), std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(t_45deg_11(1, 1), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(t_45deg_11(1, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45deg_11(2, 0), 0.0, rel_tol);
  EXPECT_NEAR(t_45deg_11(2, 1), 0.0, rel_tol);
  EXPECT_NEAR(t_45deg_11(2, 2), 1.0, rel_tol);

  trans_mat_2D<double> t_45deg_11_cpy(t_45deg_11);
  EXPECT_NEAR(t_45deg_11_cpy.getAngle(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_NEAR(t_45deg_11_cpy.getTranslation()[0], 1.0, rel_tol);
  t_45deg_11_cpy.setAngle(0.0);
  EXPECT_NEAR(t_45deg_11_cpy.getAngle(), 0.0, rel_tol);
  t_45deg_11_cpy.setRotMat(R452D());
  EXPECT_NEAR(t_45deg_11_cpy.getAngle(), 0.25 * std::numbers::pi, rel_tol);
  t_45deg_11_cpy.setTranslation(vect<double, 2>(-1.0, -1.0));
  EXPECT_NEAR(t_45deg_11_cpy.getTranslation()[0], -1.0, rel_tol);
  EXPECT_NEAR(t_45deg_11_cpy.getTranslation()[1], -1.0, rel_tol);
  t_45deg_11_cpy.setTranslation(vect<double, 2>(1.0, 1.0));
  EXPECT_NEAR(trace(t_45deg_11_cpy), 1.0 + std::sqrt(2.0), rel_tol);
  EXPECT_NEAR(determinant(t_45deg_11_cpy), 1.0, rel_tol);
  EXPECT_NEAR(invert(t_45deg_11).getAngle(), -0.25 * std::numbers::pi, rel_tol);
  EXPECT_NEAR(invert(t_45deg_11).getTranslation()[0], -std::sqrt(2.0),
              10.0 * rel_tol);
}

TEST(Transform2d, MatrixConversion) {
  trans_mat_2D<double> t_45deg_11(T452D11());
  mat<double, mat_structure::square> m_45deg_11 = t_45deg_11.getMat();
  EXPECT_NEAR(m_45deg_11(0, 0), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(m_45deg_11(0, 1), -std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(m_45deg_11(0, 2), 1.0, rel_tol);
  EXPECT_NEAR(m_45deg_11(1, 0), std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(m_45deg_11(1, 1), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(m_45deg_11(1, 2), 1.0, rel_tol);
  EXPECT_NEAR(m_45deg_11(2, 0), 0.0, rel_tol);
  EXPECT_NEAR(m_45deg_11(2, 1), 0.0, rel_tol);
  EXPECT_NEAR(m_45deg_11(2, 2), 1.0, rel_tol);
  mat<double, mat_structure::square> mt_45deg_11(transpose(t_45deg_11));
  EXPECT_NEAR(mt_45deg_11(0, 0), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(mt_45deg_11(0, 1), std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(mt_45deg_11(0, 2), 0.0, rel_tol);
  EXPECT_NEAR(mt_45deg_11(1, 0), -std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(mt_45deg_11(1, 1), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(mt_45deg_11(1, 2), 0.0, rel_tol);
  EXPECT_NEAR(mt_45deg_11(2, 0), 1.0, rel_tol);
  EXPECT_NEAR(mt_45deg_11(2, 1), 1.0, rel_tol);
  EXPECT_NEAR(mt_45deg_11(2, 2), 1.0, rel_tol);
}

TEST(Transform2d, SymMatrixConversion) {
  const mat<double, mat_structure::symmetric> msym_45deg_11(
      T452D11().getSymPart());
  EXPECT_NEAR(msym_45deg_11(0, 0), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(msym_45deg_11(0, 1), 0.0, rel_tol);
  EXPECT_NEAR(msym_45deg_11(0, 2), 0.5, rel_tol);
  EXPECT_NEAR(msym_45deg_11(1, 0), 0.0, rel_tol);
  EXPECT_NEAR(msym_45deg_11(1, 1), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(msym_45deg_11(1, 2), 0.5, rel_tol);
  EXPECT_NEAR(msym_45deg_11(2, 0), 0.5, rel_tol);
  EXPECT_NEAR(msym_45deg_11(2, 1), 0.5, rel_tol);
  EXPECT_NEAR(msym_45deg_11(2, 2), 1.0, rel_tol);
}

TEST(Transform2d, SkewSymMatrixConversion) {
  const mat<double, mat_structure::skew_symmetric> mskw_45deg_11(
      T452D11().getSkewSymPart());
  EXPECT_NEAR(mskw_45deg_11(0, 0), 0.0, rel_tol);
  EXPECT_NEAR(mskw_45deg_11(0, 1), -std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(mskw_45deg_11(0, 2), 0.5, rel_tol);
  EXPECT_NEAR(mskw_45deg_11(1, 0), std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(mskw_45deg_11(1, 1), 0.0, rel_tol);
  EXPECT_NEAR(mskw_45deg_11(1, 2), 0.5, rel_tol);
  EXPECT_NEAR(mskw_45deg_11(2, 0), -0.5, rel_tol);
  EXPECT_NEAR(mskw_45deg_11(2, 1), -0.5, rel_tol);
  EXPECT_NEAR(mskw_45deg_11(2, 2), 0.0, rel_tol);
}

TEST(Transform2d, MatrixProducts) {
  trans_mat_2D<double> t_45deg_11(T452D11());
  EXPECT_NEAR((invert(t_45deg_11) * t_45deg_11).getAngle(), 0.0, rel_tol);
  EXPECT_NEAR((invert(t_45deg_11) * t_45deg_11).getTranslation()[0], 0.0,
              rel_tol);
  trans_mat_2D<double> t_ident;
  t_ident *= t_45deg_11;
  EXPECT_NEAR(t_ident.getAngle(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_NEAR(t_ident.getTranslation()[0], 1.0, rel_tol);
  t_ident *= invert(t_45deg_11);
  EXPECT_NEAR(t_ident.getAngle(), 0.0, rel_tol);
  EXPECT_NEAR(t_ident.getTranslation()[0], 0.0, rel_tol);

  vect<double, 2> v2 = t_45deg_11 * vect<double, 2>(1.0, 1.0);
  EXPECT_NEAR(v2[0], 1.0, rel_tol);
  EXPECT_NEAR(v2[1], 1.0 + std::sqrt(2.0), rel_tol);

  vect<double, 3> v3 = t_45deg_11 * vect<double, 3>(1.0, 1.0, 1.0);
  EXPECT_NEAR(v3[0], 1.0, rel_tol);
  EXPECT_NEAR(v3[1], 1.0 + std::sqrt(2.0), rel_tol);
  EXPECT_NEAR(v3[2], 1.0, rel_tol);
  v2 = t_45deg_11.rotate(vect<double, 2>(1.0, 1.0));
  EXPECT_NEAR(v2[0], 0.0, rel_tol);
  EXPECT_NEAR(v2[1], std::sqrt(2.0), 10.0 * rel_tol);

  mat<double, mat_structure::rectangular> m23_t(2, 3);
  m23_t(0, 0) = 1.0;
  m23_t(1, 0) = 0.0;
  m23_t(0, 1) = 0.0;
  m23_t(1, 1) = 1.0;
  m23_t(0, 2) = 1.0;
  m23_t(1, 2) = 1.0;
  mat<double, mat_structure::rectangular> m23_r45_11(m23_t * t_45deg_11);
  EXPECT_NEAR(m23_r45_11(0, 0), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(m23_r45_11(0, 1), -std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(m23_r45_11(0, 2), 2.0, rel_tol);
  EXPECT_NEAR(m23_r45_11(1, 0), std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(m23_r45_11(1, 1), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(m23_r45_11(1, 2), 2.0, rel_tol);

  mat<double, mat_structure::rectangular> m32_t(3, 2);
  m32_t(0, 0) = 1.0;
  m32_t(1, 0) = 0.0;
  m32_t(2, 0) = 1.0;
  m32_t(0, 1) = 0.0;
  m32_t(1, 1) = 1.0;
  m32_t(2, 1) = 1.0;
  mat<double, mat_structure::rectangular> m32_r45_11(t_45deg_11 * m32_t);
  EXPECT_NEAR(m32_r45_11(0, 0), std::sqrt(0.5) + 1.0, rel_tol);
  EXPECT_NEAR(m32_r45_11(0, 1), -std::sqrt(0.5) + 1.0, 100.0 * rel_tol);
  EXPECT_NEAR(m32_r45_11(1, 0), std::sqrt(0.5) + 1.0, rel_tol);
  EXPECT_NEAR(m32_r45_11(1, 1), std::sqrt(0.5) + 1.0, rel_tol);
  EXPECT_NEAR(m32_r45_11(2, 0), 1.0, rel_tol);
  EXPECT_NEAR(m32_r45_11(2, 1), 1.0, rel_tol);
}

TEST(Rotation3d, Identity) {
  rot_mat_3D<double> r_ident;
  EXPECT_TRUE(is_diagonal(r_ident, rel_tol));
  EXPECT_NEAR(elem_norm_2(r_ident), std::sqrt(3.0), rel_tol);
}

TEST(Rotation3d, Copy) {
  rot_mat_3D<double> r_45z = R45Z();
  EXPECT_NEAR(r_45z(0, 0), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(r_45z(0, 1), -std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(r_45z(0, 2), 0.0, rel_tol);
  EXPECT_NEAR(r_45z(1, 0), std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(r_45z(1, 1), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(r_45z(1, 2), 0.0, rel_tol);
  EXPECT_NEAR(r_45z(2, 0), 0.0, rel_tol);
  EXPECT_NEAR(r_45z(2, 1), 0.0, rel_tol);
  EXPECT_NEAR(r_45z(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(trace(r_45z), std::sqrt(2.0) + 1.0, rel_tol);
  EXPECT_NEAR(determinant(r_45z), 1.0, rel_tol);
  rot_mat_3D<double> r_45z_cpy(r_45z);
  EXPECT_NEAR(r_45z_cpy(0, 0), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(r_45z_cpy(0, 1), -std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(r_45z_cpy(0, 2), 0.0, rel_tol);
  EXPECT_NEAR(r_45z_cpy(1, 0), std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(r_45z_cpy(1, 1), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(r_45z_cpy(1, 2), 0.0, rel_tol);
  EXPECT_NEAR(r_45z_cpy(2, 0), 0.0, rel_tol);
  EXPECT_NEAR(r_45z_cpy(2, 1), 0.0, rel_tol);
  EXPECT_NEAR(r_45z_cpy(2, 2), 1.0, rel_tol);
  r_45z_cpy = r_45z;
  EXPECT_NEAR(r_45z_cpy(0, 0), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(r_45z_cpy(0, 1), -std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(r_45z_cpy(0, 2), 0.0, rel_tol);
  EXPECT_NEAR(r_45z_cpy(1, 0), std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(r_45z_cpy(1, 1), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(r_45z_cpy(1, 2), 0.0, rel_tol);
  EXPECT_NEAR(r_45z_cpy(2, 0), 0.0, rel_tol);
  EXPECT_NEAR(r_45z_cpy(2, 1), 0.0, rel_tol);
  EXPECT_NEAR(r_45z_cpy(2, 2), 1.0, rel_tol);
}

TEST(Rotation3d, Products) {
  rot_mat_3D<double> r_45z = R45Z();
  rot_mat_3D<double> r_ident;
  r_ident *= r_45z;
  EXPECT_NEAR(r_ident(0, 0), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(r_ident(0, 1), -std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(r_ident(0, 2), 0.0, rel_tol);
  EXPECT_NEAR(r_ident(1, 0), std::sqrt(0.5), 10.0 * rel_tol);
  EXPECT_NEAR(r_ident(1, 1), std::sqrt(0.5), rel_tol);
  EXPECT_NEAR(r_ident(1, 2), 0.0, rel_tol);
  EXPECT_NEAR(r_ident(2, 0), 0.0, rel_tol);
  EXPECT_NEAR(r_ident(2, 1), 0.0, rel_tol);
  EXPECT_NEAR(r_ident(2, 2), 1.0, rel_tol);
  r_ident *= invert(r_45z);
  EXPECT_THAT(r_ident, MatrixIsIdentity(rel_tol));
  vect<double, 3> v1(1.0, 1.0, 2.0);
  EXPECT_NEAR(norm_2((r_45z * v1) - vect<double, 3>(0.0, std::sqrt(2.0), 2.0)),
              0.0, 10.0 * rel_tol);
  EXPECT_NEAR(norm_2((v1 * r_45z) - vect<double, 3>(std::sqrt(2.0), 0.0, 2.0)),
              0.0, 10.0 * rel_tol);
}

TEST(Rotation3d, SymMatrixConversion) {
  rot_mat_3D<double> r_45z = R45Z();
  mat<double, mat_structure::symmetric> msym_45z(r_45z.getSymPart());
  EXPECT_THAT(msym_45z, MatrixIsDiagonal(rel_tol));
  EXPECT_NEAR(elem_norm_2(msym_45z), std::sqrt(2.0), rel_tol);
}

TEST(Rotation3d, SkewSymMatrixConversion) {
  rot_mat_3D<double> r_45z = R45Z();
  mat<double, mat_structure::skew_symmetric> mskw_45z(r_45z.getSkewSymPart());
  EXPECT_NEAR(elem_norm_2(mskw_45z + transpose(mskw_45z)), 0.0, rel_tol);
  EXPECT_NEAR(elem_norm_2(r_45z.getMat() - r_45z), 0.0, rel_tol);
  EXPECT_NEAR(elem_norm_2(invert(r_45z).getMat() - transpose(r_45z)), 0.0,
              rel_tol);
  EXPECT_THAT(invert(r_45z) * r_45z, MatrixIsIdentity(rel_tol));
}

TEST(Quaternion, Identity) {
  quaternion<double> q_ident;
  EXPECT_NEAR(q_ident[0], 1.0, rel_tol);
  EXPECT_NEAR(q_ident[1], 0.0, rel_tol);
  EXPECT_NEAR(q_ident[2], 0.0, rel_tol);
  EXPECT_NEAR(q_ident[3], 0.0, rel_tol);
}

TEST(Quaternion, Copy) {
  quaternion<double> q_45z = Q45Z();
  EXPECT_NEAR(q_45z[0], std::cos(0.125 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(q_45z[1], 0.0, rel_tol);
  EXPECT_NEAR(q_45z[2], 0.0, rel_tol);
  EXPECT_NEAR(q_45z[3], std::sin(0.125 * std::numbers::pi), 10.0 * rel_tol);
  quaternion<double> q_45z_cpy(q_45z);
  EXPECT_NEAR(q_45z_cpy[0], std::cos(0.125 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(q_45z_cpy[1], 0.0, rel_tol);
  EXPECT_NEAR(q_45z_cpy[2], 0.0, rel_tol);
  EXPECT_NEAR(q_45z_cpy[3], std::sin(0.125 * std::numbers::pi), 10.0 * rel_tol);
  q_45z_cpy = q_45z;
  EXPECT_NEAR(q_45z_cpy[0], std::cos(0.125 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(q_45z_cpy[1], 0.0, rel_tol);
  EXPECT_NEAR(q_45z_cpy[2], 0.0, rel_tol);
  EXPECT_NEAR(q_45z_cpy[3], std::sin(0.125 * std::numbers::pi), 10.0 * rel_tol);
}

TEST(Quaternion, Products) {
  quaternion<double> q_ident = Q45Z() * Q45Z();
  EXPECT_NEAR(q_ident[0], std::cos(0.25 * std::numbers::pi), 10.0 * rel_tol);
  EXPECT_NEAR(q_ident[1], 0.0, rel_tol);
  EXPECT_NEAR(q_ident[2], 0.0, rel_tol);
  EXPECT_NEAR(q_ident[3], std::sin(0.25 * std::numbers::pi), rel_tol);
  q_ident *= invert(Q45Z());
  EXPECT_NEAR(q_ident[0], std::cos(0.125 * std::numbers::pi), 10.0 * rel_tol);
  EXPECT_NEAR(q_ident[1], 0.0, rel_tol);
  EXPECT_NEAR(q_ident[2], 0.0, rel_tol);
  EXPECT_NEAR(q_ident[3], std::sin(0.125 * std::numbers::pi), 10.0 * rel_tol);
  q_ident *= invert(Q45Z());
  EXPECT_NEAR(q_ident[0], 1.0, 100.0 * rel_tol);
  EXPECT_NEAR(q_ident[1], 0.0, rel_tol);
  EXPECT_NEAR(q_ident[2], 0.0, rel_tol);
  EXPECT_NEAR(q_ident[3], 0.0, rel_tol);
  EXPECT_THAT(R45Z() * invert(Q45Z()), MatrixIsIdentity(rel_tol));
}

TEST(Quaternion, MatrixConversions) {
  quaternion<double> q_45z = Q45Z();
  mat<double, mat_structure::square> rm_90z(q_45z * R45Z());
  EXPECT_NEAR(rm_90z(0, 0), 0.0, rel_tol);
  EXPECT_NEAR(rm_90z(0, 1), -1.0, rel_tol);
  EXPECT_NEAR(rm_90z(0, 2), 0.0, rel_tol);
  EXPECT_NEAR(rm_90z(1, 0), 1.0, rel_tol);
  EXPECT_NEAR(rm_90z(1, 1), 0.0, rel_tol);
  EXPECT_NEAR(rm_90z(1, 2), 0.0, rel_tol);
  EXPECT_NEAR(rm_90z(2, 0), 0.0, rel_tol);
  EXPECT_NEAR(rm_90z(2, 1), 0.0, rel_tol);
  EXPECT_NEAR(rm_90z(2, 2), 1.0, rel_tol);

  EXPECT_NEAR(trace(q_45z), 1.0 + std::sqrt(2.0), rel_tol);
  EXPECT_NEAR(determinant(q_45z), 1.0, rel_tol);
  EXPECT_THAT(q_45z.getSymPart(), MatrixIsNear(R45Z().getSymPart(), rel_tol));
  EXPECT_THAT(q_45z.getSkewSymPart(),
              MatrixIsNear(R45Z().getSkewSymPart(), rel_tol));
  EXPECT_NEAR((q_45z * invert(q_45z))[0], 1.0, rel_tol);
  EXPECT_NEAR((q_45z * transpose(q_45z))[0], 1.0, rel_tol);
}

TEST(Quaternion, EulerAnglesTBConversions) {
  quaternion<double> q_e_90z(Q45Z() * E45Z());
  EXPECT_NEAR(q_e_90z[0], std::cos(0.25 * std::numbers::pi), 10.0 * rel_tol);
  EXPECT_NEAR(q_e_90z[1], 0.0, rel_tol);
  EXPECT_NEAR(q_e_90z[2], 0.0, rel_tol);
  EXPECT_NEAR(q_e_90z[3], std::sin(0.25 * std::numbers::pi), rel_tol);
}

TEST(Quaternion, AxisAngleConversions) {
  quaternion<double> q_a_90z(Q45Z() * A45Z());
  EXPECT_NEAR(q_a_90z[0], std::cos(0.25 * std::numbers::pi), 10.0 * rel_tol);
  EXPECT_NEAR(q_a_90z[1], 0.0, rel_tol);
  EXPECT_NEAR(q_a_90z[2], 0.0, rel_tol);
  EXPECT_NEAR(q_a_90z[3], std::sin(0.25 * std::numbers::pi), rel_tol);
}

TEST(Quaternion, VectorOperations) {
  vect<double, 3> v1(1.0, 1.0, 2.0);
  EXPECT_THAT(Q45Z() * v1,
              VectorIsNear(vect<double, 3>(0.0, std::sqrt(2.0), 2.0), rel_tol));
}

TEST(EulerAnglesTB, Identity) {
  euler_angles_TB<double> e_ident;
  EXPECT_NEAR(e_ident.yaw(), 0.0, rel_tol);
  EXPECT_NEAR(e_ident.pitch(), 0.0, rel_tol);
  EXPECT_NEAR(e_ident.roll(), 0.0, rel_tol);
}

TEST(EulerAnglesTB, Copy) {
  euler_angles_TB<double> e_45z_cpy(E45Z());
  EXPECT_NEAR(e_45z_cpy.yaw(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_NEAR(e_45z_cpy.pitch(), 0.0, rel_tol);
  EXPECT_NEAR(e_45z_cpy.roll(), 0.0, rel_tol);
  e_45z_cpy = E45Z();
  EXPECT_NEAR(e_45z_cpy.yaw(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_NEAR(e_45z_cpy.pitch(), 0.0, rel_tol);
  EXPECT_NEAR(e_45z_cpy.roll(), 0.0, rel_tol);
  EXPECT_NEAR(trace(e_45z_cpy), 1.0 + std::sqrt(2.0), rel_tol);
  EXPECT_NEAR(determinant(e_45z_cpy), 1.0, rel_tol);
}

TEST(EulerAnglesTB, Rotation3dConversion) {
  euler_angles_TB<double> e_45z(R45Z());
  EXPECT_NEAR(e_45z.yaw(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_NEAR(e_45z.pitch(), 0.0, rel_tol);
  EXPECT_NEAR(e_45z.roll(), 0.0, rel_tol);
  e_45z = R45Z();
  EXPECT_NEAR(e_45z.yaw(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_NEAR(e_45z.pitch(), 0.0, rel_tol);
  EXPECT_NEAR(e_45z.roll(), 0.0, rel_tol);
  EXPECT_THAT(e_45z.getRotMat().getMat(),
              MatrixIsNear(R45Z().getMat(), rel_tol));
  EXPECT_THAT(e_45z.getMat(), MatrixIsNear(R45Z().getMat(), rel_tol));
  EXPECT_THAT(e_45z.getSymPart(), MatrixIsNear(R45Z().getSymPart(), rel_tol));
  EXPECT_THAT(e_45z.getSkewSymPart(),
              MatrixIsNear(R45Z().getSkewSymPart(), rel_tol));
}

TEST(EulerAnglesTB, QuaternionConversion) {
  euler_angles_TB<double> e_45z(Q45Z());
  EXPECT_NEAR(e_45z.yaw(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_NEAR(e_45z.pitch(), 0.0, rel_tol);
  EXPECT_NEAR(e_45z.roll(), 0.0, rel_tol);
  e_45z = Q45Z();
  EXPECT_NEAR(e_45z.yaw(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_NEAR(e_45z.pitch(), 0.0, rel_tol);
  EXPECT_NEAR(e_45z.roll(), 0.0, rel_tol);
  EXPECT_NEAR((Q45Z() * invert(e_45z))[0], 1.0, rel_tol);
  EXPECT_NEAR((Q45Z() * transpose(e_45z))[0], 1.0, rel_tol);
}

TEST(EulerAnglesTB, AxisAngleConversion) {
  euler_angles_TB<double> e_45z(A45Z());
  EXPECT_NEAR(e_45z.yaw(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_NEAR(e_45z.pitch(), 0.0, rel_tol);
  EXPECT_NEAR(e_45z.roll(), 0.0, rel_tol);
  e_45z = A45Z();
  EXPECT_NEAR(e_45z.yaw(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_NEAR(e_45z.pitch(), 0.0, rel_tol);
  EXPECT_NEAR(e_45z.roll(), 0.0, rel_tol);
}

TEST(AxisAngle, Identity) {
  axis_angle<double> a_ident;
  EXPECT_NEAR(a_ident.angle(), 0.0, rel_tol);
}

TEST(AxisAngle, Copy) {
  axis_angle<double> a_45z_cpy(A45Z());
  EXPECT_NEAR(a_45z_cpy.angle(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_THAT(a_45z_cpy.axis(),
              VectorIsNear(vect<double, 3>(0.0, 0.0, 1.0), rel_tol));
  a_45z_cpy = A45Z();
  EXPECT_NEAR(a_45z_cpy.angle(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_THAT(a_45z_cpy.axis(),
              VectorIsNear(vect<double, 3>(0.0, 0.0, 1.0), rel_tol));
  EXPECT_NEAR(trace(a_45z_cpy), 1.0 + std::sqrt(2.0), rel_tol);
  EXPECT_NEAR(determinant(a_45z_cpy), 1.0, rel_tol);
}

TEST(AxisAngle, Rotation3dConversion) {
  axis_angle<double> a_45z(R45Z());
  EXPECT_NEAR(a_45z.angle(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_THAT(a_45z.axis(),
              VectorIsNear(vect<double, 3>(0.0, 0.0, 1.0), rel_tol));
  a_45z = R45Z();
  EXPECT_NEAR(a_45z.angle(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_THAT(a_45z.axis(),
              VectorIsNear(vect<double, 3>(0.0, 0.0, 1.0), rel_tol));
  EXPECT_THAT(a_45z.getRotMat().getMat(),
              MatrixIsNear(R45Z().getMat(), rel_tol));
  EXPECT_THAT(a_45z.getMat(), MatrixIsNear(R45Z().getMat(), rel_tol));
  EXPECT_THAT(a_45z.getSymPart(), MatrixIsNear(R45Z().getSymPart(), rel_tol));
  EXPECT_THAT(a_45z.getSkewSymPart(),
              MatrixIsNear(R45Z().getSkewSymPart(), rel_tol));
}

TEST(AxisAngle, QuaternionConversion) {
  axis_angle<double> a_45z(Q45Z());
  EXPECT_NEAR(a_45z.angle(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_THAT(a_45z.axis(),
              VectorIsNear(vect<double, 3>(0.0, 0.0, 1.0), rel_tol));
  a_45z = Q45Z();
  EXPECT_NEAR(a_45z.angle(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_THAT(a_45z.axis(),
              VectorIsNear(vect<double, 3>(0.0, 0.0, 1.0), rel_tol));
  EXPECT_NEAR((Q45Z() * invert(a_45z))[0], 1.0, rel_tol);
  EXPECT_NEAR((Q45Z() * transpose(a_45z))[0], 1.0, rel_tol);
}

TEST(AxisAngle, EulerAnglesTBConversion) {
  axis_angle<double> a_45z(E45Z());
  EXPECT_NEAR(a_45z.angle(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_THAT(a_45z.axis(),
              VectorIsNear(vect<double, 3>(0.0, 0.0, 1.0), rel_tol));
  a_45z = E45Z();
  EXPECT_NEAR(a_45z.angle(), 0.25 * std::numbers::pi, rel_tol);
  EXPECT_THAT(a_45z.axis(),
              VectorIsNear(vect<double, 3>(0.0, 0.0, 1.0), rel_tol));
}

TEST(Transform3d, Identity) {
  trans_mat_3D<double> t_ident;
  EXPECT_THAT(t_ident, MatrixIsIdentity(rel_tol));
}

TEST(Transform3d, Copy) {
  trans_mat_3D<double> t_45z(T45Z());
  EXPECT_NEAR(elem_norm_2(t_45z), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z(3, 3), 1.0, rel_tol);
  trans_mat_3D<double> t_45z_cpy(t_45z);
  EXPECT_NEAR(elem_norm_2(t_45z_cpy), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z_cpy(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_cpy(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_cpy(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z_cpy(3, 3), 1.0, rel_tol);
  t_45z_cpy = t_45z;
  EXPECT_NEAR(elem_norm_2(t_45z_cpy), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z_cpy(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_cpy(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_cpy(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z_cpy(3, 3), 1.0, rel_tol);
  t_45z_cpy = trans_mat_3D<double>();
  t_45z_cpy *= t_45z;
  EXPECT_NEAR(elem_norm_2(t_45z_cpy), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z_cpy(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_cpy(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_cpy(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z_cpy(3, 3), 1.0, rel_tol);
}

TEST(Transform3d, Rotation3dConversion) {
  trans_mat_3D<double> t_45z(R45Z());
  EXPECT_NEAR(elem_norm_2(t_45z), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z(3, 3), 1.0, rel_tol);
  t_45z = R45Z();
  EXPECT_NEAR(elem_norm_2(t_45z), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z(3, 3), 1.0, rel_tol);
  rot_mat_3D<double> t_45z_rot = t_45z.getRotMat();
  EXPECT_NEAR(elem_norm_2(t_45z_rot), std::sqrt(3.0), rel_tol);
  EXPECT_NEAR(t_45z_rot(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_rot(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_rot(2, 2), 1.0, rel_tol);
  mat<double, mat_structure::square> t_45z_mat = t_45z.getMat();
  EXPECT_NEAR(elem_norm_2(t_45z_mat), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z_mat(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_mat(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_mat(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z_mat(3, 3), 1.0, rel_tol);
  t_45z = trans_mat_3D<double>();
  t_45z *= R45Z();
  EXPECT_NEAR(elem_norm_2(t_45z), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z(3, 3), 1.0, rel_tol);
}

TEST(Transform3d, QuaternionConversion) {
  trans_mat_3D<double> t_45z(Q45Z());
  EXPECT_NEAR(elem_norm_2(t_45z), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z(3, 3), 1.0, rel_tol);
  t_45z = Q45Z();
  EXPECT_NEAR(elem_norm_2(t_45z), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z(3, 3), 1.0, rel_tol);
  t_45z = trans_mat_3D<double>();
  t_45z *= Q45Z();
  EXPECT_NEAR(elem_norm_2(t_45z), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z(3, 3), 1.0, rel_tol);
}

TEST(Transform3d, EulerAnglesTBConversion) {
  trans_mat_3D<double> t_45z(E45Z());
  EXPECT_NEAR(elem_norm_2(t_45z), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z(3, 3), 1.0, rel_tol);
  t_45z = E45Z();
  EXPECT_NEAR(elem_norm_2(t_45z), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z(3, 3), 1.0, rel_tol);
  t_45z = trans_mat_3D<double>();
  t_45z *= E45Z();
  EXPECT_NEAR(elem_norm_2(t_45z), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z(3, 3), 1.0, rel_tol);
}

TEST(Transform3d, AxisAngleConversion) {
  trans_mat_3D<double> t_45z(A45Z());
  EXPECT_NEAR(elem_norm_2(t_45z), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z(3, 3), 1.0, rel_tol);
  t_45z = A45Z();
  EXPECT_NEAR(elem_norm_2(t_45z), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z(3, 3), 1.0, rel_tol);
  t_45z = trans_mat_3D<double>();
  t_45z *= A45Z();
  EXPECT_NEAR(elem_norm_2(t_45z), std::sqrt(4.0), rel_tol);
  EXPECT_NEAR(t_45z(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z(3, 3), 1.0, rel_tol);
}

TEST(Transform3d, WithTranslation) {
  trans_mat_3D<double> t_45z_123(T45Z123());
  EXPECT_NEAR(elem_norm_2(t_45z_123), std::sqrt(18.0), 2.0 * rel_tol);
  EXPECT_NEAR(t_45z_123(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_123(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_123(0, 3), 1.0, rel_tol);
  EXPECT_NEAR(t_45z_123(1, 3), 2.0, rel_tol);
  EXPECT_NEAR(t_45z_123(2, 3), 3.0, rel_tol);
  EXPECT_NEAR(t_45z_123(3, 3), 1.0, rel_tol);
  EXPECT_NEAR(trace(t_45z_123), 2.0 + std::sqrt(2.0), rel_tol);
  EXPECT_NEAR(determinant(t_45z_123), 1.0, rel_tol);
  mat<double, mat_structure::square> t_45z_mat_sym_skw(
      t_45z_123.getSymPart() + t_45z_123.getSkewSymPart());
  EXPECT_NEAR(elem_norm_2(t_45z_mat_sym_skw), std::sqrt(18.0), 2.0 * rel_tol);
  EXPECT_NEAR(t_45z_mat_sym_skw(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_mat_sym_skw(0, 1), -std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_mat_sym_skw(0, 3), 1.0, rel_tol);
  EXPECT_NEAR(t_45z_mat_sym_skw(1, 3), 2.0, rel_tol);
  EXPECT_NEAR(t_45z_mat_sym_skw(2, 3), 3.0, rel_tol);
  EXPECT_NEAR(t_45z_mat_sym_skw(3, 3), 1.0, rel_tol);
  EXPECT_THAT(t_45z_123 * invert(t_45z_123), MatrixIsIdentity(rel_tol));
  EXPECT_THAT(invert(t_45z_123).getMat() * t_45z_123,
              MatrixIsIdentity(rel_tol));
  EXPECT_THAT(invert(t_45z_123) * t_45z_123.getMat(),
              MatrixIsIdentity(rel_tol));
  mat<double, mat_structure::square> t_45z_123_t = transpose(t_45z_123);
  EXPECT_NEAR(elem_norm_2(t_45z_123_t), std::sqrt(18.0), 2.0 * rel_tol);
  EXPECT_NEAR(t_45z_123_t(0, 0), std::cos(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_123_t(0, 1), std::sin(0.25 * std::numbers::pi), rel_tol);
  EXPECT_NEAR(t_45z_123_t(2, 2), 1.0, rel_tol);
  EXPECT_NEAR(t_45z_123_t(0, 3), 0.0, rel_tol);
  EXPECT_NEAR(t_45z_123_t(3, 1), 2.0, rel_tol);
  EXPECT_NEAR(t_45z_123_t(3, 2), 3.0, rel_tol);
  EXPECT_NEAR(t_45z_123_t(3, 3), 1.0, rel_tol);
  EXPECT_THAT((invert(t_45z_123) * R45Z()) *
                  trans_mat_3D<double>(rot_mat_3D<double>(),
                                       -invert(t_45z_123).getTranslation()),
              MatrixIsIdentity(rel_tol));
}

TEST(Transform3d, VectorOperations) {
  trans_mat_3D<double> t_45z_123(T45Z123());
  vect<double, 3> v1(1.0, 1.0, 2.0);
  vect<double, 3> v1_trans(t_45z_123 * v1);
  EXPECT_NEAR(v1_trans[0], 1.0, rel_tol);
  EXPECT_NEAR(v1_trans[1], 2.0 + std::sqrt(2.0), rel_tol);
  EXPECT_NEAR(v1_trans[2], 5.0, rel_tol);
  v1_trans = t_45z_123.rotate(v1);
  EXPECT_NEAR(v1_trans[0], 0.0, rel_tol);
  EXPECT_NEAR(v1_trans[1], std::sqrt(2.0), 10.0 * rel_tol);
  EXPECT_NEAR(v1_trans[2], 2.0, rel_tol);
  vect<double, 4> v3 = t_45z_123 * vect<double, 4>(1, 1, 2, 2);
  EXPECT_NEAR(v3[0], 2.0, rel_tol);
  EXPECT_NEAR(v3[1], 4.0 + std::sqrt(2.0), rel_tol);
  EXPECT_NEAR(v3[2], 8.0, rel_tol);
  EXPECT_NEAR(v3[3], 2.0, rel_tol);
}

TEST(Rotations3d, WeirdCases) {
  axis_angle<double> a_weird(0.3241, vect<double, 3>(0.5, 0.5, sqrt(0.5)));
  quaternion<double> q_weird;
  q_weird = a_weird;
  euler_angles_TB<double> e_weird;
  e_weird = q_weird;
  rot_mat_3D<double> r_weird;
  r_weird = q_weird;
  axis_angle<double> a_weird_q(q_weird);
  axis_angle<double> a_weird_e(e_weird);
  axis_angle<double> a_weird_r(r_weird);
  EXPECT_NEAR(a_weird_q.angle(), a_weird.angle(), rel_tol);
  EXPECT_THAT(a_weird_q.axis(), VectorIsNear(a_weird.axis(), rel_tol));
  EXPECT_NEAR(a_weird_e.angle(), a_weird.angle(), rel_tol);
  EXPECT_THAT(a_weird_e.axis(), VectorIsNear(a_weird.axis(), rel_tol));
  EXPECT_NEAR(a_weird_r.angle(), a_weird.angle(), rel_tol);
  EXPECT_THAT(a_weird_r.axis(), VectorIsNear(a_weird.axis(), rel_tol));

  quaternion<double> q_res(Q45Z() * quaternion<double>(a_weird * A45Z()));
  vect<double, 4> v_a(q_res[0], q_res[1], q_res[2], q_res[3]);
  q_res = Q45Z() * (q_weird * Q45Z());
  vect<double, 4> v_q(q_res[0], q_res[1], q_res[2], q_res[3]);
  q_res = quaternion<double>(Q45Z() * e_weird * R45Z());
  vect<double, 4> v_e(q_res[0], q_res[1], q_res[2], q_res[3]);
  q_res = quaternion<double>(A45Z() * r_weird * E45Z());
  vect<double, 4> v_r(q_res[0], q_res[1], q_res[2], q_res[3]);
  EXPECT_THAT(v_a, VectorIsNear(v_q, rel_tol));
  EXPECT_THAT(v_a, VectorIsNear(v_e, rel_tol));
  EXPECT_THAT(v_a, VectorIsNear(v_r, rel_tol));
  EXPECT_THAT(v_q, VectorIsNear(v_e, rel_tol));
  EXPECT_THAT(v_q, VectorIsNear(v_r, rel_tol));
  EXPECT_THAT(v_e, VectorIsNear(v_r, rel_tol));
}

}  // namespace
}  // namespace ReaK

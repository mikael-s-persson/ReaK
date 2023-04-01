
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

#include "ReaK/math/kinetostatics/quat_alg.h"
#include "ReaK/math/lin_alg/vect_matchers.h"
#include "gtest/gtest.h"

#include <cmath>

namespace ReaK {
namespace {

using ::ReaK::testing::VectorIsNear;

const double tolerance = 10.0 * std::numeric_limits<double>::epsilon();

TEST(QuatAlg, QuatTests) {
  quat<double> q_45z(std::cos(0.125 * M_PI), 0.0, 0.0, std::sin(0.125 * M_PI));
  EXPECT_THAT(q_45z, VectorIsNear(vect<double, 4>(std::cos(0.125 * M_PI), 0.0,
                                                  0.0, std::sin(0.125 * M_PI)),
                                  tolerance));
  quat<double> q_ident = q_45z * conj(q_45z);
  EXPECT_THAT(q_ident,
              VectorIsNear(vect<double, 4>(1.0, 0.0, 0.0, 0.0), tolerance));
  quat<double> q_zero = q_45z - q_45z;
  EXPECT_THAT(q_zero,
              VectorIsNear(vect<double, 4>(0.0, 0.0, 0.0, 0.0), tolerance));
  q_zero = q_45z + (-q_45z);
  EXPECT_THAT(q_zero,
              VectorIsNear(vect<double, 4>(0.0, 0.0, 0.0, 0.0), tolerance));
  quat<double> q_45z_cpy(q_45z);
  EXPECT_THAT(q_45z_cpy,
              VectorIsNear(vect<double, 4>(std::cos(0.125 * M_PI), 0.0, 0.0,
                                           std::sin(0.125 * M_PI)),
                           tolerance));
  q_45z_cpy = q_45z;
  EXPECT_THAT(q_45z_cpy,
              VectorIsNear(vect<double, 4>(std::cos(0.125 * M_PI), 0.0, 0.0,
                                           std::sin(0.125 * M_PI)),
                           tolerance));
  q_ident = q_45z * q_45z;
  EXPECT_THAT(q_ident, VectorIsNear(vect<double, 4>(std::cos(0.25 * M_PI), 0.0,
                                                    0.0, std::sin(0.25 * M_PI)),
                                    tolerance));
  q_ident *= conj(q_45z);
  EXPECT_THAT(q_ident,
              VectorIsNear(vect<double, 4>(std::cos(0.125 * M_PI), 0.0, 0.0,
                                           std::sin(0.125 * M_PI)),
                           tolerance));
  q_ident *= invert(q_45z);
  EXPECT_THAT(q_ident,
              VectorIsNear(vect<double, 4>(1.0, 0.0, 0.0, 0.0), tolerance));
  EXPECT_NEAR(norm_2_sqr(q_45z), 1.0, tolerance);
  EXPECT_NEAR(norm_2(q_45z), 1.0, tolerance);
  q_45z *= 2.0;
  q_45z = unit(q_45z);
  EXPECT_THAT(q_45z, VectorIsNear(vect<double, 4>(std::cos(0.125 * M_PI), 0.0,
                                                  0.0, std::sin(0.125 * M_PI)),
                                  tolerance));

  q_45z = sqrt(pow(q_45z, quat<double>(2.0)));
  EXPECT_THAT(q_45z, VectorIsNear(vect<double, 4>(std::cos(0.125 * M_PI), 0.0,
                                                  0.0, std::sin(0.125 * M_PI)),
                                  tolerance));

  quat<double> temp = cos(q_45z) * cos(q_45z) + sin(q_45z) * sin(q_45z);
  EXPECT_THAT(temp,
              VectorIsNear(vect<double, 4>(1.0, 0.0, 0.0, 0.0), tolerance));
  temp = invert(cos(q_45z) * cos(q_45z)) - tan(q_45z) * tan(q_45z);
  EXPECT_THAT(temp,
              VectorIsNear(vect<double, 4>(1.0, 0.0, 0.0, 0.0), tolerance));
  temp = acos(cos(q_45z)) * invert(q_45z);
  EXPECT_THAT(temp,
              VectorIsNear(vect<double, 4>(1.0, 0.0, 0.0, 0.0), tolerance));
  temp = asin(sin(q_45z)) * invert(q_45z);
  EXPECT_THAT(temp,
              VectorIsNear(vect<double, 4>(1.0, 0.0, 0.0, 0.0), tolerance));
  temp = atan(tan(q_45z)) * invert(q_45z);
  EXPECT_THAT(temp,
              VectorIsNear(vect<double, 4>(1.0, 0.0, 0.0, 0.0), tolerance));
  temp = exp(q_45z) * exp(q_45z) - exp(q_45z + q_45z);
  EXPECT_THAT(temp,
              VectorIsNear(vect<double, 4>(0.0, 0.0, 0.0, 0.0), tolerance));
  temp = log(q_45z) + log(q_45z) - log(q_45z * q_45z);
  EXPECT_THAT(temp,
              VectorIsNear(vect<double, 4>(0.0, 0.0, 0.0, 0.0), tolerance));
}

TEST(QuatAlg, UnitQuatTests) {
  unit_quat<double> q_45z(std::cos(0.125 * M_PI), 0.0, 0.0,
                          std::sin(0.125 * M_PI));
  EXPECT_THAT(q_45z, VectorIsNear(vect<double, 4>(std::cos(0.125 * M_PI), 0.0,
                                                  0.0, std::sin(0.125 * M_PI)),
                                  tolerance));
  unit_quat<double> q_ident = q_45z * conj(q_45z);
  EXPECT_THAT(q_ident,
              VectorIsNear(vect<double, 4>(1.0, 0.0, 0.0, 0.0), tolerance));
  unit_quat<double> q_45z_cpy(q_45z);
  EXPECT_THAT(q_45z_cpy,
              VectorIsNear(vect<double, 4>(std::cos(0.125 * M_PI), 0.0, 0.0,
                                           std::sin(0.125 * M_PI)),
                           tolerance));
  q_45z_cpy = q_45z;
  EXPECT_THAT(q_45z_cpy,
              VectorIsNear(vect<double, 4>(std::cos(0.125 * M_PI), 0.0, 0.0,
                                           std::sin(0.125 * M_PI)),
                           tolerance));
  q_ident = q_45z * q_45z;
  EXPECT_THAT(q_ident, VectorIsNear(vect<double, 4>(std::cos(0.25 * M_PI), 0.0,
                                                    0.0, std::sin(0.25 * M_PI)),
                                    tolerance));
  q_ident *= conj(q_45z);
  EXPECT_THAT(q_ident,
              VectorIsNear(vect<double, 4>(std::cos(0.125 * M_PI), 0.0, 0.0,
                                           std::sin(0.125 * M_PI)),
                           tolerance));
  q_ident *= invert(q_45z);
  EXPECT_THAT(q_ident,
              VectorIsNear(vect<double, 4>(1.0, 0.0, 0.0, 0.0), tolerance));
  EXPECT_NEAR(norm_2_sqr(q_45z), 1.0, tolerance);
  EXPECT_NEAR(norm_2(q_45z), 1.0, tolerance);

  vect<double, 3> v_45z(0.0, 0.0, 0.125 * M_PI);
  quat<double> temp = (exp(v_45z) * exp(v_45z)) * conj(exp(v_45z + v_45z));
  EXPECT_THAT(temp,
              VectorIsNear(vect<double, 4>(1.0, 0.0, 0.0, 0.0), tolerance));
  vect<double, 3> v_zero = log(q_45z) + log(q_45z) - log(q_45z * q_45z);
  EXPECT_THAT(v_zero, VectorIsNear(vect<double, 3>(0.0, 0.0, 0.0), tolerance));

  q_45z = sqrt(q_45z * q_45z);
  EXPECT_THAT(q_45z, VectorIsNear(vect<double, 4>(std::cos(0.125 * M_PI), 0.0,
                                                  0.0, std::sin(0.125 * M_PI)),
                                  tolerance));

  q_45z = sqrt(pow(q_45z, 2.0));
  EXPECT_THAT(q_45z, VectorIsNear(vect<double, 4>(std::cos(0.125 * M_PI), 0.0,
                                                  0.0, std::sin(0.125 * M_PI)),
                                  tolerance));
}

}  // namespace
}  // namespace ReaK

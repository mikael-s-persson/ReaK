
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

#include "ReaK/math/lin_alg/vect_alg.h"
#include "ReaK/core/base/defs.h"

#include <cstdio>
#include <fstream>
#include <iostream>

#include "gtest/gtest.h"

namespace ReaK {
namespace {

template <typename T>
class Vect3Test : public ::testing::Test {};

using Vect3TestTypes = ::testing::Types<vect<float, 3>, vect<double, 3>>;
TYPED_TEST_SUITE(Vect3Test, Vect3TestTypes);

TYPED_TEST(Vect3Test, Vect3Operators) {
  using Vector = TypeParam;
  using ValueType = vect_value_type_t<Vector>;
  using std::abs;

  Vector gravity_acc(0.0, -9.81, 0.0);
  EXPECT_EQ(gravity_acc.size(), 3);
  EXPECT_TRUE(gravity_acc.max_size());
  EXPECT_TRUE(gravity_acc.capacity());
  EXPECT_TRUE(!gravity_acc.empty());

  auto it = gravity_acc.begin();
  EXPECT_NEAR(*it, 0.0, std::numeric_limits<ValueType>::epsilon());
  EXPECT_EQ(gravity_acc.end() - it, 3);
  ++it;
  EXPECT_EQ(gravity_acc.end() - it, 2);

  const Vector& gravity_acc_ref = gravity_acc;
  auto cit = gravity_acc_ref.begin();
  EXPECT_NEAR(*cit, 0.0, std::numeric_limits<ValueType>::epsilon());
  EXPECT_EQ(gravity_acc.end() - cit, 3);
  ++cit;
  EXPECT_EQ(gravity_acc.end() - cit, 2);

  std::array<ValueType, 3> ones = {1.0, 1.0, 1.0};
  Vector ones_v(ones.data());
  EXPECT_NEAR(ones[1], ones_v[1], std::numeric_limits<ValueType>::epsilon());

  ValueType obj_mass(3.0);
  Vector gravity_force;
  gravity_force = (gravity_acc * obj_mass);
  EXPECT_NEAR(gravity_force[1], ValueType(-3.0 * 9.81),
              std::numeric_limits<ValueType>::epsilon());
  EXPECT_TRUE((gravity_force == obj_mass * gravity_acc));
  EXPECT_TRUE((gravity_force == obj_mass * gravity_acc));
  EXPECT_NEAR(norm_2_sqr(gravity_acc), ValueType(9.81 * 9.81),
              100.0 * std::numeric_limits<ValueType>::epsilon());
  EXPECT_NEAR(norm_2(gravity_acc), ValueType(9.81),
              10.0 * std::numeric_limits<ValueType>::epsilon());

  Vector gravity_dir(unit(gravity_acc));
  EXPECT_NEAR(norm_2(gravity_dir), ValueType(1.0),
              std::numeric_limits<ValueType>::epsilon());

  Vector displacement(1.0, 2.0, 3.0);
  ValueType gravity_potential = gravity_force * displacement;
  EXPECT_NEAR(gravity_potential, ValueType(-2.0 * 9.81 * 3.0),
              100.0 * std::numeric_limits<ValueType>::epsilon());

  Vector gravity_moment = displacement % gravity_force;
  EXPECT_NEAR((abs(gravity_moment[0] - 3.0 * 9.81 * 3.0) +
               abs(gravity_moment[2] + 3.0 * 9.81)),
              0.0, 100.0 * std::numeric_limits<ValueType>::epsilon());

  ValueType dist_sqr = norm_2_sqr(displacement);
  Vector displacement_inv(-displacement);
  EXPECT_NEAR(displacement_inv * displacement, -dist_sqr,
              std::numeric_limits<ValueType>::epsilon());

  Vector long_displacement(10.0, 20.0, 30.0);
  EXPECT_TRUE(colinear(displacement, long_displacement));
  EXPECT_TRUE((long_displacement > displacement));
  EXPECT_TRUE((displacement < long_displacement));
  EXPECT_TRUE((displacement != long_displacement));
  EXPECT_TRUE((displacement == displacement));

  long_displacement += displacement;
  EXPECT_NEAR(norm_2(long_displacement), 11.0 * norm_2(displacement),
              200.0 * std::numeric_limits<ValueType>::epsilon());
  long_displacement *= 2.0;
  EXPECT_NEAR(norm_2(long_displacement), 22.0 * norm_2(displacement),
              400.0 * std::numeric_limits<ValueType>::epsilon());
  long_displacement /= 2.0;
  EXPECT_NEAR(norm_2(long_displacement), 11.0 * norm_2(displacement),
              200.0 * std::numeric_limits<ValueType>::epsilon());
  long_displacement -= displacement;
  EXPECT_NEAR(norm_2(long_displacement), 10.0 * norm_2(displacement),
              200.0 * std::numeric_limits<ValueType>::epsilon());
}

template <typename T>
class VectNTest : public ::testing::Test {};

using VectNTestTypes = ::testing::Types<vect_n<float>, vect_n<double>>;
TYPED_TEST_SUITE(VectNTest, VectNTestTypes);

TYPED_TEST(VectNTest, VectNOperators) {
  using Vector = TypeParam;
  using ValueType = vect_value_type_t<Vector>;
  using std::abs;

  Vector gravity_acc(ValueType(0.0), ValueType(-9.81), ValueType(0.0));
  EXPECT_EQ(gravity_acc.size(), 3);
  EXPECT_TRUE(gravity_acc.max_size());
  EXPECT_TRUE(gravity_acc.capacity());
  EXPECT_TRUE(!gravity_acc.empty());

  auto it = gravity_acc.begin();
  EXPECT_NEAR(*it, 0.0, std::numeric_limits<ValueType>::epsilon());
  EXPECT_EQ(gravity_acc.end() - it, 3);
  ++it;
  EXPECT_EQ(gravity_acc.end() - it, 2);

  const Vector& gravity_acc_ref = gravity_acc;
  auto cit = gravity_acc_ref.begin();
  EXPECT_NEAR(*cit, 0.0, std::numeric_limits<ValueType>::epsilon());
  EXPECT_EQ(gravity_acc.end() - cit, 3);
  ++cit;
  EXPECT_EQ(gravity_acc.end() - cit, 2);

  ValueType obj_mass(3.0);
  Vector gravity_force;
  gravity_force = (gravity_acc * obj_mass);
  EXPECT_NEAR(gravity_force[1], ValueType(-3.0 * 9.81),
              std::numeric_limits<ValueType>::epsilon());
  EXPECT_TRUE((gravity_force == obj_mass * gravity_acc));
  EXPECT_TRUE((gravity_force == obj_mass * gravity_acc));
  EXPECT_NEAR(norm_2_sqr(gravity_acc), ValueType(9.81 * 9.81),
              100.0 * std::numeric_limits<ValueType>::epsilon());
  EXPECT_NEAR(norm_2(gravity_acc), ValueType(9.81),
              10.0 * std::numeric_limits<ValueType>::epsilon());

  Vector gravity_dir(unit(gravity_acc));
  EXPECT_NEAR(norm_2(gravity_dir), ValueType(1.0),
              std::numeric_limits<ValueType>::epsilon());

  Vector displacement(1.0, 2.0, 3.0);
  ValueType gravity_potential = gravity_force * displacement;
  EXPECT_NEAR(gravity_potential, ValueType(-2.0 * 9.81 * 3.0),
              100.0 * std::numeric_limits<ValueType>::epsilon());

  ValueType dist_sqr = norm_2_sqr(displacement);
  Vector displacement_inv(-displacement);
  EXPECT_NEAR(displacement_inv * displacement, -dist_sqr,
              std::numeric_limits<ValueType>::epsilon());

  Vector long_displacement(10.0, 20.0, 30.0);
  EXPECT_TRUE(colinear(displacement, long_displacement));
  EXPECT_TRUE((long_displacement > displacement));
  EXPECT_TRUE((displacement < long_displacement));
  EXPECT_TRUE((displacement != long_displacement));
  EXPECT_TRUE((displacement == displacement));

  long_displacement += displacement;
  EXPECT_NEAR(norm_2(long_displacement), 11.0 * norm_2(displacement),
              200.0 * std::numeric_limits<ValueType>::epsilon());
  long_displacement *= 2.0;
  EXPECT_NEAR(norm_2(long_displacement), 22.0 * norm_2(displacement),
              400.0 * std::numeric_limits<ValueType>::epsilon());
  long_displacement /= 2.0;
  EXPECT_NEAR(norm_2(long_displacement), 11.0 * norm_2(displacement),
              200.0 * std::numeric_limits<ValueType>::epsilon());
  long_displacement -= displacement;
  EXPECT_NEAR(norm_2(long_displacement), 10.0 * norm_2(displacement),
              200.0 * std::numeric_limits<ValueType>::epsilon());
}

TEST(Vect, VectConstructionTests) {
  EXPECT_NEAR(norm_2(ReaK::vect<double, 1>(1.0)), 1.0,
              std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(norm_2(ReaK::vect<double, 2>(0.0, 1.0)), 1.0,
              std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(norm_2(ReaK::vect<double, 3>(0.0, 0.0, 1.0)), 1.0,
              std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(norm_2(ReaK::vect<double, 4>(0.0, 0.0, 0.0, 1.0)), 1.0,
              std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(norm_2(ReaK::vect<double, 5>(0.0, 0.0, 0.0, 0.0, 1.0)), 1.0,
              std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(norm_2(ReaK::vect<double, 6>(0.0, 0.0, 0.0, 0.0, 0.0, 1.0)), 1.0,
              std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(norm_2(ReaK::vect<double, 7>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)),
              1.0, std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(
      norm_2(ReaK::vect<double, 8>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)),
      1.0, std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(norm_2(ReaK::vect<double, 9>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                           0.0, 1.0)),
              1.0, std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(norm_2(ReaK::vect<double, 10>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                            0.0, 0.0, 1.0)),
              1.0, std::numeric_limits<double>::epsilon());

  EXPECT_NEAR(norm_2(ReaK::vect_n<double>(0.0, 0.0, 1.0)), 1.0,
              std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(norm_2(ReaK::vect_n<double>(0.0, 0.0, 0.0, 1.0)), 1.0,
              std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(norm_2(ReaK::vect_n<double>(0.0, 0.0, 0.0, 0.0, 1.0)), 1.0,
              std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(norm_2(ReaK::vect_n<double>(0.0, 0.0, 0.0, 0.0, 0.0, 1.0)), 1.0,
              std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(norm_2(ReaK::vect_n<double>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)),
              1.0, std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(
      norm_2(ReaK::vect_n<double>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)), 1.0,
      std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(
      norm_2(ReaK::vect_n<double>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)),
      1.0, std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(norm_2(ReaK::vect_n<double>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.0, 1.0)),
              1.0, std::numeric_limits<double>::epsilon());
}

}  // namespace
}  // namespace ReaK

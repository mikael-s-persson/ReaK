
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

#include "ReaK/math/lin_alg/vect_matchers.h"
#include "ReaK/math/optimization/finite_diff_jacobians.h"
#include "ReaK/math/optimization/line_search.h"
#include "ReaK/math/optimization/mehrotra_method.h"
#include "ReaK/math/optimization/quadratic_programs.h"
#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
#include <sstream>

namespace ReaK::optim {
namespace {

using ::ReaK::testing::VectorIsNear;

const double desired_tolerance = 1e-8;
const double expected_tolerance = 1e-8;

class QPProblemsTest : public ::testing::Test {
 protected:
  QPProblemsTest();

  std::string GetProblemDescription(int i) {
    std::stringstream ss;
    ss << "Quadratic problem #" << i << " with:"
       << "\n  G = " << Gs[i] << "\n  c = " << cs[i] << "\n  A = " << As[i]
       << "\n  b = " << bs[i] << std::endl;
    return ss.str();
  }

  std::vector<mat<double, mat_structure::symmetric>> Gs;
  std::vector<mat<double, mat_structure::rectangular>> As;
  std::vector<vect_n<double>> bs;
  std::vector<vect_n<double>> cs;
  std::vector<vect_n<double>> xs;
  std::vector<vect_n<double>> x_sols;
};

QPProblemsTest::QPProblemsTest() {
  Gs.emplace_back(4.0, 1.5, 5.0, 8.0, 2.0, 12.0);
  cs.emplace_back(1.0, 1.0, 0.5);
  bs.emplace_back(1, 0.0);
  mat<double, mat_structure::rectangular> A(1, 3);
  A(0, 0) = -2.0;
  A(0, 1) = 1.0;
  A(0, 2) = 4.0;
  As.push_back(A);
  xs.emplace_back(1.0, 1.0, 1.0);
  x_sols.emplace_back(-0.10659045, -0.11634163, -0.02420982);

  Gs.emplace_back(6.0, 2.0, 1.0, 5.0, 2.0, 4.0);
  cs.emplace_back(-8.0, -3.0, -3.0);
  bs.emplace_back(vect<double, 2>(3.0, 0.0));
  A = mat<double, mat_structure::rectangular>(2, 3);
  A(0, 0) = 1.0;
  A(0, 1) = 0.0;
  A(0, 2) = 1.0;
  A(1, 0) = 0.0;
  A(1, 1) = 1.0;
  A(1, 2) = 1.0;
  As.push_back(A);
  xs.emplace_back(1.0, 1.0, 1.0);
  x_sols.emplace_back(2.0, -1.0, 1.0);
}

TEST_F(QPProblemsTest, NullSpaceQPMethod) {
  vect_n<double> x;
  for (std::size_t i = 0; i < Gs.size(); ++i) {
    x = xs[i];
    EXPECT_NO_THROW(optim::null_space_QP_method(As[i], bs[i], Gs[i], cs[i], x,
                                                desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(x_sols[i], expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(QPProblemsTest, ProjectedCGMethod) {
  vect_n<double> x;
  for (std::size_t i = 0; i < Gs.size(); ++i) {
    x = xs[i];
    EXPECT_NO_THROW(optim::projected_CG_method(As[i], bs[i], Gs[i], cs[i], x,
                                               100, desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(x_sols[i], expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(QPProblemsTest, MehrotraQPMethod) {
  vect_n<double> x;
  for (std::size_t i = 0; i < Gs.size(); ++i) {
    x = xs[i];
    EXPECT_NO_THROW(optim::mehrotra_QP_method(
        mat<double, mat_structure::rectangular>(0, x.size()), vect_n<double>(0),
        Gs[i], cs[i], As[i], bs[i], x, 100, desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(x_sols[i], expected_tolerance))
        << GetProblemDescription(i);
  }
}

}  // namespace
}  // namespace ReaK::optim

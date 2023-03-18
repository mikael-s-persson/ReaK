
/*
 *    Copyright 2023 Sven Mikael Persson
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

#include "ReaK/math/optimization/finite_diff_jacobians.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <functional>
#include <sstream>
#include <string>
#include <vector>

namespace ReaK::optim {
namespace {

using ::testing::DoubleNear;

const double expected_tolerance = 1e-3;

class FiniteDiffProblemsTest : public ::testing::Test {
 protected:
  FiniteDiffProblemsTest();

  std::string GetProblemDescription(int i) {
    std::stringstream ss;
    ss << "Function #" << i << " called: '" << func_names[i] << "'"
       << std::endl;
    return ss.str();
  }

  using FunctionPtr = std::function<double(double)>;
  std::vector<FunctionPtr> funcs;
  std::vector<double> func_sols;
  std::vector<double> func_starts;
  std::vector<std::string> func_names;
};

FiniteDiffProblemsTest::FiniteDiffProblemsTest() {
  func_names.emplace_back("Max Truss Section Stress");
  funcs.emplace_back([&](double x) {
    double sigma1 = 0.8165 / x;
    double sigma2 = 1.1154 / (1 - x);
    return (sigma1 > sigma2 ? sigma1 : sigma2);
  });
  func_sols.emplace_back(4.4616);
  func_starts.emplace_back(0.5);
}

TEST_F(FiniteDiffProblemsTest, ComputeJacobian2ptsForward) {
  ReaK::mat<double, ReaK::mat_structure::rectangular> J(1, 1);
  for (int i = 0; i < funcs.size(); ++i) {
    double x = func_starts[i];
    double y = funcs[i](x);
    J(0, 0) = 0.0;
    EXPECT_NO_THROW(compute_jacobian_2pts_forward(funcs[i], x, y, J))
        << GetProblemDescription(i);
    EXPECT_THAT(J(0, 0), DoubleNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(FiniteDiffProblemsTest, ComputeJacobian2ptsCentral) {
  ReaK::mat<double, ReaK::mat_structure::rectangular> J(1, 1);
  for (int i = 0; i < funcs.size(); ++i) {
    double x = func_starts[i];
    double y = funcs[i](x);
    J(0, 0) = 0.0;
    EXPECT_NO_THROW(compute_jacobian_2pts_central(funcs[i], x, y, J))
        << GetProblemDescription(i);
    EXPECT_THAT(J(0, 0), DoubleNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(FiniteDiffProblemsTest, ComputeJacobian5ptsCentral) {
  ReaK::mat<double, ReaK::mat_structure::rectangular> J(1, 1);
  for (int i = 0; i < funcs.size(); ++i) {
    double x = func_starts[i];
    double y = funcs[i](x);
    J(0, 0) = 0.0;
    EXPECT_NO_THROW(compute_jacobian_5pts_central(funcs[i], x, y, J))
        << GetProblemDescription(i);
    EXPECT_THAT(J(0, 0), DoubleNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
  }
}

}  // namespace
}  // namespace ReaK::optim

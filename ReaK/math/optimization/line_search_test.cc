
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

#include "ReaK/math/optimization/line_search.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <functional>
#include <sstream>
#include <string>
#include <vector>

namespace ReaK::optim {
namespace {

using ::testing::DoubleNear;

const double desired_tolerance = 1e-4;

class LineSearchProblemsTest : public ::testing::Test {
 protected:
  LineSearchProblemsTest();

  std::string GetProblemDescription(int i) {
    std::stringstream ss;
    ss << "Function #" << i << " called: '" << func_names[i] << "'"
       << std::endl;
    return ss.str();
  }

  int eval_count = 0;
  using FunctionPtr = std::function<double(double)>;
  std::vector<FunctionPtr> funcs;
  std::vector<double> func_sols;
  std::vector<double> func_lowers;
  std::vector<double> func_uppers;
  std::vector<std::string> func_names;
};

LineSearchProblemsTest::LineSearchProblemsTest() {
  func_names.emplace_back("Max Truss Section Stress");
  funcs.emplace_back([&](double x) {
    ++(this->eval_count);
    double sigma1 = 0.8165 / x;
    double sigma2 = 1.1154 / (1 - x);
    return (sigma1 > sigma2 ? sigma1 : sigma2);
  });
  func_sols.emplace_back(0.4227);
  func_lowers.emplace_back(0.30);
  func_uppers.emplace_back(0.48);
}

TEST_F(LineSearchProblemsTest, DichotomousSearch) {
  for (int i = 0; i < funcs.size(); ++i) {
    eval_count = 0;
    double l = func_lowers[i];
    double u = func_uppers[i];
    EXPECT_NO_THROW(dichotomous_search(funcs[i], l, u, desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT((l + u) * 0.5, DoubleNear(func_sols[i], desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_LT(eval_count, 32) << GetProblemDescription(i);
  }
}

TEST_F(LineSearchProblemsTest, GoldenSectionSearch) {
  for (int i = 0; i < funcs.size(); ++i) {
    eval_count = 0;
    double l = func_lowers[i];
    double u = func_uppers[i];
    EXPECT_NO_THROW(golden_section_search(funcs[i], l, u, desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT((l + u) * 0.5, DoubleNear(func_sols[i], desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_LT(eval_count, 32) << GetProblemDescription(i);
  }
}

TEST_F(LineSearchProblemsTest, FibonacciSearch) {
  for (int i = 0; i < funcs.size(); ++i) {
    eval_count = 0;
    double l = func_lowers[i];
    double u = func_uppers[i];
    EXPECT_NO_THROW(fibonacci_search(funcs[i], l, u, desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT((l + u) * 0.5, DoubleNear(func_sols[i], desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_LT(eval_count, 32) << GetProblemDescription(i);
  }
}

}  // namespace
}  // namespace ReaK::optim

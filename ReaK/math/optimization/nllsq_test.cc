
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

#include "ReaK/math/lin_alg/mat_svd_method.h"
#include "ReaK/math/lin_alg/vect_matchers.h"
#include "ReaK/math/optimization/finite_diff_jacobians.h"
#include "ReaK/math/optimization/gauss_newton_method.h"
#include "ReaK/math/optimization/jacobian_transpose_method.h"
#include "ReaK/math/optimization/levenberg_marquardt_method.h"
#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
#include <sstream>

namespace ReaK::optim {
namespace {

using ::ReaK::testing::VectorIsNear;

const double desired_tolerance = 1e-8;
const double desired_func_tolerance = 1e-14;
const double expected_tolerance = 1e-7;

class NLLSQProblemsTest : public ::testing::Test {
 protected:
  NLLSQProblemsTest();

  std::string GetProblemDescription(int i) {
    std::stringstream ss;
    ss << "Function #" << i << " called: '" << func_names[i] << "'"
       << std::endl;
    return ss.str();
  }

  using FunctionPtr = std::function<vect_n<double>(const vect_n<double>&)>;
  using JacFunctionPtr =
      std::function<void(mat<double, mat_structure::rectangular>&,
                         const vect_n<double>&, const vect_n<double>&)>;

  int eval_count = 0;
  int grad_count = 0;

  std::vector<FunctionPtr> funcs;
  std::vector<JacFunctionPtr> func_jacs;
  std::vector<vect_n<double>> func_sols;
  std::vector<vect_n<double>> func_starts;
  std::vector<vect_n<double>> func_lowers;
  std::vector<vect_n<double>> func_uppers;
  std::vector<std::string> func_names;
};

NLLSQProblemsTest::NLLSQProblemsTest() {
  func_names.emplace_back("Linear Function, full-rank");
  funcs.emplace_back([&](const vect_n<double>& x) {
    ++(this->eval_count);
    double sum = 0.0;
    for (double value : x) {
      sum += value;
    }
    vect_n<double> result(2 * x.size(), -1.0 - sum / double(x.size()));
    for (int i = 0; i < x.size(); ++i) {
      result[i] += x[i];
    }
    return result;
  });
  func_jacs.emplace_back([&](mat<double, mat_structure::rectangular>& J,
                             const vect_n<double>& x, const vect_n<double>& f) {
    ++(this->grad_count);
    J.set_col_count(x.size());
    J.set_row_count(f.size());
    for (int i = 0; i < f.size(); ++i) {
      for (int j = 0; j < x.size(); ++j) {
        J(i, j) = -2.0 / double(f.size());
      }
    }
    for (std::size_t i = 0; i < x.size(); ++i) {
      J(i, i) += 1.0;
    }
  });
  func_sols.emplace_back(10, -1.0);
  func_starts.emplace_back(10, 1.0);
  func_lowers.emplace_back(10, -4.0);
  func_uppers.emplace_back(10, 4.0);

  func_names.emplace_back("Linear Function, rank 1");
  funcs.emplace_back([&](const vect_n<double>& x) {
    ++(this->eval_count);
    double sum = 0.0;
    for (int i = 0; i < x.size(); ++i) {
      sum += double(i + 1) * x[i];
    }
    vect_n<double> result(x.size() * 2);
    for (int i = 0; i < result.size(); ++i) {
      result[i] = double(i + 1) * sum - 1.0;
    }
    return result;
  });
  func_jacs.emplace_back([&](mat<double, mat_structure::rectangular>& J,
                             const vect_n<double>& x, const vect_n<double>& f) {
    ++(this->grad_count);
    J.set_col_count(x.size());
    J.set_row_count(f.size());
    for (std::size_t i = 0; i < f.size(); ++i) {
      for (std::size_t j = 0; j < x.size(); ++j) {
        J(i, j) = double((i + 1) * (j + 1));
      }
    }
  });
  func_sols.emplace_back(10, 6.0 / double(41 * 11 * 10));
  func_starts.emplace_back(10, 1.0);
  func_lowers.emplace_back(10, -2.0);
  func_uppers.emplace_back(10, 2.0);

  func_names.emplace_back("Rosenbrock Function");
  funcs.emplace_back([&](const vect_n<double>& x) {
    ++(this->eval_count);
    vect_n<double> result(2);
    result[0] = 10.0 * (x[1] - x[0] * x[0]);
    result[1] = 1.0 - x[0];
    return result;
  });
  func_jacs.emplace_back([&](mat<double, mat_structure::rectangular>& J,
                             const vect_n<double>& x,
                             const vect_n<double>& /*unused*/) {
    ++(this->grad_count);
    J.set_col_count(2);
    J.set_row_count(2);
    J(0, 0) = -20.0 * x[0];
    J(0, 1) = 10.0;
    J(1, 0) = -1.0;
    J(1, 1) = 0.0;
  });
  func_sols.emplace_back(2, 1.0);
  func_starts.emplace_back(2, -0.5);
  func_lowers.emplace_back(2, -1.0);
  func_uppers.emplace_back(2, 2.0);

  func_names.emplace_back("Helical-valley Function");
  funcs.emplace_back([&](const vect_n<double>& x) {
    ++(this->eval_count);
    vect_n<double> result(3);
    double tmp = atan2(x[1], x[0]);
    result[0] = 10.0 * (x[2] - 10.0 * tmp);
    result[1] = 10.0 * (sqrt(x[0] * x[0] + x[1] * x[1]) - 1.0);
    result[2] = x[2];
    return result;
  });
  func_jacs.emplace_back([&](mat<double, mat_structure::rectangular>& J,
                             const vect_n<double>& x,
                             const vect_n<double>& /*unused*/) {
    ++(this->grad_count);
    J.set_col_count(3);
    J.set_row_count(3);

    J(0, 0) = 50.0 * x[1] / (M_PI * (x[0] * x[0] + x[1] * x[1]));
    J(0, 1) = -50.0 * x[0] / (M_PI * (x[0] * x[0] + x[1] * x[1]));
    J(0, 2) = 10.0;

    J(1, 0) = 10.0 * x[0] / (sqrt(x[0] * x[0] + x[1] * x[1]));
    J(1, 1) = 10.0 * x[1] / (sqrt(x[0] * x[0] + x[1] * x[1]));
    J(1, 2) = 0.0;

    J(2, 0) = 0.0;
    J(2, 1) = 0.0;
    J(2, 2) = 1.0;
  });
  func_sols.emplace_back(1.0, 0.0, 0.0);
  func_starts.emplace_back(-1.0, 0.0, 0.0);
  func_lowers.emplace_back(-2.0, -1.0, -1.0);
  func_uppers.emplace_back(2.0, 1.0, 1.0);

  func_names.emplace_back("Powell singular Function");
  funcs.emplace_back([&](const vect_n<double>& x) {
    ++(this->eval_count);
    vect_n<double> result(4);
    result[0] = x[0] + 10.0 * x[1];
    result[1] = sqrt(5.0) * (x[2] - x[3]);
    result[2] = (x[1] - 2.0 * x[2]) * (x[1] - 2.0 * x[2]);
    result[3] = sqrt(10.0) * (x[0] - x[3]) * (x[0] - x[3]);
    return result;
  });
  func_jacs.emplace_back([&](mat<double, mat_structure::rectangular>& J,
                             const vect_n<double>& x,
                             const vect_n<double>& /*unused*/) {
    ++(this->grad_count);
    J.set_col_count(4);
    J.set_row_count(4);
    J(0, 0) = 1.0;
    J(0, 1) = 10.0;
    J(0, 2) = 0.0;
    J(0, 3) = 0.0;

    J(1, 0) = 0.0;
    J(1, 1) = 0.0;
    J(1, 2) = sqrt(5.0);
    J(1, 3) = -sqrt(5.0);

    J(2, 0) = 0.0;
    J(2, 1) = 2.0 * (x[1] - 2.0 * x[2]);
    J(2, 2) = -4.0 * (x[1] - 2.0 * x[2]);
    J(2, 3) = 0.0;

    J(3, 0) = 2.0 * sqrt(10.0) * (x[0] - x[3]);
    J(3, 1) = 0.0;
    J(3, 2) = 0.0;
    J(3, 3) = -2.0 * sqrt(10.0) * (x[0] - x[3]);
  });
  func_sols.emplace_back(0.0, 0.0, 0.0, 0.0);
  func_starts.emplace_back(3.0, -1.0, 0.0, 1.0);
  func_lowers.emplace_back(-4.0, -2.0, -2.0, -2.0);
  func_uppers.emplace_back(4.0, 2.0, 2.0, 2.0);

  func_names.emplace_back("Freudenstein-Roth Function");
  funcs.emplace_back([&](const vect_n<double>& x) {
    ++(this->eval_count);
    vect_n<double> result(2);
    result[0] = -13.0 + x[0] + ((5.0 - x[1]) * x[1] - 2.0) * x[1];
    result[1] = -29.0 + x[0] + ((1.0 + x[1]) * x[1] - 14.0) * x[1];
    return result;
  });
  func_jacs.emplace_back([&](mat<double, mat_structure::rectangular>& J,
                             const vect_n<double>& x,
                             const vect_n<double>& /*unused*/) {
    ++(this->grad_count);
    J.set_col_count(2);
    J.set_row_count(2);
    J(0, 0) = 1.0;
    J(0, 1) = x[1] * (10.0 - 3.0 * x[1]) - 2.0;
    J(1, 0) = 1.0;
    J(1, 1) = x[1] * (2.0 + 3.0 * x[1]) - 14.0;
  });
  func_sols.emplace_back(5.0, 4.0);
  func_starts.emplace_back(0.5, -2.0);
  func_lowers.emplace_back(0.0, -5.0);
  func_uppers.emplace_back(8.0, 8.0);
}

TEST_F(NLLSQProblemsTest, GaussNewtonNLLSQ) {
  std::unordered_set<std::string> expected_failures = {
      "Linear Function, rank 1", "Helical-valley Function"};
  vect_n<double> x;
  vect_n<double> y;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    y = funcs[i](x);
    y *= 0.0;
    eval_count = 0;
    grad_count = 0;
    EXPECT_NO_THROW(optim::gauss_newton_nllsq(funcs[i], func_jacs[i], x, y, 200,
                                              desired_tolerance,
                                              desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(NLLSQProblemsTest, JacobianTransposeNLLSQ) {
  std::unordered_set<std::string> expected_failures = {
      "Linear Function, rank 1", "Rosenbrock Function",
      "Helical-valley Function", "Powell singular Function",
      "Freudenstein-Roth Function"};
  vect_n<double> x;
  vect_n<double> y;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    y = funcs[i](x);
    y *= 0.0;
    eval_count = 0;
    grad_count = 0;
    EXPECT_NO_THROW(optim::jacobian_transpose_nllsq(funcs[i], func_jacs[i], x,
                                                    y, 500, desired_tolerance,
                                                    desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(NLLSQProblemsTest, LevenbergMarquardtNLLSQ) {
  std::unordered_set<std::string> expected_failures = {
      "Linear Function, rank 1", "Helical-valley Function",
      "Powell singular Function", "Freudenstein-Roth Function"};
  vect_n<double> x;
  vect_n<double> y;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    y = funcs[i](x);
    y *= 0.0;
    eval_count = 0;
    grad_count = 0;
    EXPECT_NO_THROW(optim::levenberg_marquardt_nllsq(
        funcs[i], func_jacs[i], x, y, 200, 1e-4, desired_func_tolerance,
        desired_tolerance, desired_func_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(NLLSQProblemsTest, GaussNewtonNLLSQ_SVD) {
  std::unordered_set<std::string> expected_failures = {
      "Linear Function, rank 1", "Helical-valley Function"};
  vect_n<double> x;
  vect_n<double> y;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    y = funcs[i](x);
    y *= 0.0;
    eval_count = 0;
    grad_count = 0;
    EXPECT_NO_THROW(optim::make_gauss_newton_nllsq(funcs[i], func_jacs[i], y,
                                                   200, desired_tolerance,
                                                   desired_tolerance)
                        .set_lin_solver(SVD_linlsqsolver())(x))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(NLLSQProblemsTest, LevenbergMarquardtNLLSQ_SVD) {
  std::unordered_set<std::string> expected_failures = {
      "Linear Function, rank 1", "Helical-valley Function",
      "Powell singular Function", "Freudenstein-Roth Function"};
  vect_n<double> x;
  vect_n<double> y;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    y = funcs[i](x);
    y *= 0.0;
    eval_count = 0;
    grad_count = 0;
    EXPECT_NO_THROW(optim::make_levenberg_marquardt_nllsq(
                        funcs[i], func_jacs[i], y, 200, 1e-4,
                        desired_func_tolerance, desired_tolerance,
                        desired_func_tolerance)
                        .set_lin_solver(SVD_linlsqsolver())(x))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(NLLSQProblemsTest, GaussNewtonNLLSQ_Box) {
  std::unordered_set<std::string> expected_failures = {
      "Linear Function, rank 1", "Helical-valley Function"};
  vect_n<double> x;
  vect_n<double> y;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    y = funcs[i](x);
    y *= 0.0;
    eval_count = 0;
    grad_count = 0;
    auto lim_func = [l = func_lowers[i], u = func_uppers[i]](
                        const vect_n<double>& x, vect_n<double>& dx) {
      optim::box_limit_function(x, dx, l, u);
    };
    EXPECT_NO_THROW(optim::limited_gauss_newton_nllsq(
        funcs[i], func_jacs[i], x, y, 200, lim_func, desired_tolerance,
        desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(NLLSQProblemsTest, JacobianTransposeNLLSQ_Box) {
  std::unordered_set<std::string> expected_failures = {
      "Linear Function, rank 1", "Rosenbrock Function",
      "Helical-valley Function", "Powell singular Function",
      "Freudenstein-Roth Function"};
  vect_n<double> x;
  vect_n<double> y;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    y = funcs[i](x);
    y *= 0.0;
    eval_count = 0;
    grad_count = 0;
    auto lim_func = [l = func_lowers[i], u = func_uppers[i]](
                        const vect_n<double>& x, vect_n<double>& dx) {
      optim::box_limit_function(x, dx, l, u);
    };
    EXPECT_NO_THROW(optim::limited_jacobian_transpose_nllsq(
        funcs[i], func_jacs[i], x, y, 500, lim_func, desired_tolerance,
        desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(NLLSQProblemsTest, LevenbergMarquardtNLLSQ_Box) {
  std::unordered_set<std::string> expected_failures = {
      "Linear Function, rank 1", "Helical-valley Function",
      "Powell singular Function", "Freudenstein-Roth Function"};
  vect_n<double> x;
  vect_n<double> y;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    y = funcs[i](x);
    y *= 0.0;
    eval_count = 0;
    grad_count = 0;
    auto lim_func = [l = func_lowers[i], u = func_uppers[i]](
                        const vect_n<double>& x, vect_n<double>& dx) {
      optim::box_limit_function(x, dx, l, u);
    };
    EXPECT_NO_THROW(optim::limited_levenberg_marquardt_nllsq(
        funcs[i], func_jacs[i], x, y, lim_func, 200, 1e-4,
        desired_func_tolerance, desired_tolerance, desired_func_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
  }
}

}  // namespace
}  // namespace ReaK::optim

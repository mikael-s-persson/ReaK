
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

#include "ReaK/core/base/global_rng.h"
#include "ReaK/math/lin_alg/vect_matchers.h"
#include "ReaK/math/optimization/conjugate_gradient_methods.h"
#include "ReaK/math/optimization/nelder_mead_method.h"
#include "ReaK/math/optimization/quasi_newton_methods.h"
#include "ReaK/math/optimization/trust_region_search.h"
#include "gtest/gtest.h"

#include <cmath>
#include <iostream>
#include <sstream>

namespace ReaK::optim {
namespace {

using ::ReaK::testing::VectorIsNear;

const double desired_tolerance = 1e-8;
const double expected_tolerance = 1e-7;

class NLPUnconstrainedProblemsTest : public ::testing::Test {
 protected:
  NLPUnconstrainedProblemsTest();

  std::string GetProblemDescription(int i) {
    std::stringstream ss;
    ss << "Function #" << i << " called: '" << func_names[i] << "'"
       << std::endl;
    return ss.str();
  }

  using FunctionPtr = std::function<double(const vect<double, 2>&)>;
  using GradFunctionPtr =
      std::function<vect<double, 2>(const vect<double, 2>&)>;

  int eval_count = 0;
  int grad_count = 0;

  std::vector<FunctionPtr> funcs;
  std::vector<GradFunctionPtr> func_grads;
  std::vector<vect<double, 2>> func_sols;
  std::vector<vect<double, 2>> func_starts;
  std::vector<std::string> func_names;
};

NLPUnconstrainedProblemsTest::NLPUnconstrainedProblemsTest() {
  func_names.emplace_back("Easy Function");
  mat<double, mat_structure::symmetric, mat_alignment::column_major, 2, 2> Q(
      10.0, -2.0, 1.0);
  funcs.emplace_back([this, Q](const vect<double, 2>& x) {
    ++(this->eval_count);
    return x * (Q * x);
  });
  func_grads.emplace_back([this, Q](const vect<double, 2>& x) {
    ++(this->grad_count);
    return vect<double, 2>(2.0 * (Q * x));
  });
  func_sols.emplace_back(0.0, 0.0);
  func_starts.emplace_back(0.5, 1.0);

  func_names.emplace_back("Banana Function");
  funcs.emplace_back([this](const vect<double, 2>& x) {
    ++(this->eval_count);
    return (1.0 - x[0]) * (1.0 - x[0]) +
           100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
  });
  func_grads.emplace_back([this](const vect<double, 2>& x) {
    ++(this->grad_count);
    return vect<double, 2>{
        2.0 * (x[0] - 1.0) - 400.0 * x[0] * (x[1] - x[0] * x[0]),
        200.0 * (x[1] - x[0] * x[0])};
  });
  func_sols.emplace_back(1.0, 1.0);
  func_starts.emplace_back(0.5, 0.75);
}

// TODO BROKEN!
TEST_F(NLPUnconstrainedProblemsTest, DISABLED_NelderMeadMethod) {
  std::unordered_set<std::string> expected_failures = {};
  vect<double, 2> x;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    eval_count = 0;
    grad_count = 0;
    double initial_spread = 2.0;
    EXPECT_NO_THROW(nelder_mead_method(funcs[i], x, initial_spread,
                                       get_global_rng(), desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(func_grads[i](x),
                VectorIsNear(vect<double, 2>(0.0, 0.0), expected_tolerance))
        << GetProblemDescription(i);
  }
}

// TODO BROKEN!
TEST_F(NLPUnconstrainedProblemsTest, DISABLED_BFGSMethod) {
  std::unordered_set<std::string> expected_failures = {};
  vect<double, 2> x;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    eval_count = 0;
    grad_count = 0;
    EXPECT_NO_THROW(bfgs_method(funcs[i], func_grads[i], x, 100,
                                desired_tolerance, desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(func_grads[i](x),
                VectorIsNear(vect<double, 2>(0.0, 0.0), expected_tolerance))
        << GetProblemDescription(i);
  }
}

// TODO BROKEN!
TEST_F(NLPUnconstrainedProblemsTest, DISABLED_DFPMethod) {
  std::unordered_set<std::string> expected_failures = {};
  vect<double, 2> x;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    eval_count = 0;
    grad_count = 0;
    EXPECT_NO_THROW(dfp_method(funcs[i], func_grads[i], x, 100,
                               desired_tolerance, desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(func_grads[i](x),
                VectorIsNear(vect<double, 2>(0.0, 0.0), expected_tolerance))
        << GetProblemDescription(i);
  }
}

// TODO BROKEN!
TEST_F(NLPUnconstrainedProblemsTest, DISABLED_BroydenClassMethod) {
  std::unordered_set<std::string> expected_failures = {};
  vect<double, 2> x;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    eval_count = 0;
    grad_count = 0;
    EXPECT_NO_THROW(broyden_class_method(funcs[i], func_grads[i], x, 100, 0.5,
                                         desired_tolerance, desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(func_grads[i](x),
                VectorIsNear(vect<double, 2>(0.0, 0.0), expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(NLPUnconstrainedProblemsTest, FletcherReevesConjGradMethod) {
  std::unordered_set<std::string> expected_failures = {};
  vect<double, 2> x;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    eval_count = 0;
    grad_count = 0;
    EXPECT_NO_THROW(non_linear_conj_grad_method(
        funcs[i], func_grads[i], x, 400, fletcher_reeves_beta(),
        line_search_expand_and_zoom<double>(1e-4, 0.1), desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(func_grads[i](x),
                VectorIsNear(vect<double, 2>(0.0, 0.0), expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(NLPUnconstrainedProblemsTest, PolakRibiereConjGradMethod) {
  std::unordered_set<std::string> expected_failures = {};
  vect<double, 2> x;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    eval_count = 0;
    grad_count = 0;
    EXPECT_NO_THROW(non_linear_conj_grad_method(
        funcs[i], func_grads[i], x, 200, polak_ribiere_beta(),
        line_search_expand_and_zoom<double>(1e-4, 0.1), desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(func_grads[i](x),
                VectorIsNear(vect<double, 2>(0.0, 0.0), expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(NLPUnconstrainedProblemsTest, HestenesStiefelConjGradMethod) {
  std::unordered_set<std::string> expected_failures = {};
  vect<double, 2> x;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    eval_count = 0;
    grad_count = 0;
    EXPECT_NO_THROW(non_linear_conj_grad_method(
        funcs[i], func_grads[i], x, 300, hestenes_stiefel_beta(),
        line_search_expand_and_zoom<double>(1e-4, 0.1), desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(func_grads[i](x),
                VectorIsNear(vect<double, 2>(0.0, 0.0), expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(NLPUnconstrainedProblemsTest, DaiYuanConjGradMethod) {
  std::unordered_set<std::string> expected_failures = {};
  vect<double, 2> x;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    eval_count = 0;
    grad_count = 0;
    EXPECT_NO_THROW(non_linear_conj_grad_method(
        funcs[i], func_grads[i], x, 200, dai_yuan_beta(),
        line_search_expand_and_zoom<double>(1e-4, 0.1), desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(func_grads[i](x),
                VectorIsNear(vect<double, 2>(0.0, 0.0), expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(NLPUnconstrainedProblemsTest, HagerZhangConjGradMethod) {
  std::unordered_set<std::string> expected_failures = {};
  vect<double, 2> x;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    eval_count = 0;
    grad_count = 0;
    EXPECT_NO_THROW(non_linear_conj_grad_method(
        funcs[i], func_grads[i], x, 200, hager_zhang_beta(),
        line_search_expand_and_zoom<double>(1e-4, 0.1), desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(func_grads[i](x),
                VectorIsNear(vect<double, 2>(0.0, 0.0), expected_tolerance))
        << GetProblemDescription(i);
  }
}

TEST_F(NLPUnconstrainedProblemsTest, SR1TrustRegionQuasiNewtonMethod) {
  std::unordered_set<std::string> expected_failures = {};
  vect<double, 2> x;
  for (int i = 0; i < funcs.size(); ++i) {
    if (expected_failures.count(func_names[i]) == 1) {
      std::cout << "Skipping known failure: '" << func_names[i] << "'"
                << std::endl;
      continue;
    }
    x = func_starts[i];
    eval_count = 0;
    grad_count = 0;
    EXPECT_NO_THROW(quasi_newton_trust_region(
        funcs[i], func_grads[i], x, 0.5, 100, trust_region_solver_dogleg(),
        hessian_update_sr1(), no_limit_functor(), desired_tolerance,
        desired_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(x, VectorIsNear(func_sols[i], expected_tolerance))
        << GetProblemDescription(i);
    EXPECT_THAT(func_grads[i](x),
                VectorIsNear(vect<double, 2>(0.0, 0.0), expected_tolerance))
        << GetProblemDescription(i);
  }
}

}  // namespace
}  // namespace ReaK::optim

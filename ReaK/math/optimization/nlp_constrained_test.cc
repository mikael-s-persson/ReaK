
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
#include "ReaK/math/optimization/augmented_lagrangian_methods.h"
#include "ReaK/math/optimization/nl_interior_points_methods.h"
#include "ReaK/math/optimization/sequential_qp_methods.h"

#include <cmath>
#include <iostream>
#include <sstream>

namespace ReaK::optim {
namespace {

using ::ReaK::testing::VectorIsGreater;
using ::ReaK::testing::VectorIsNear;
using ::ReaK::testing::VectorIsZero;

const double desired_tolerance = 1e-8;
const double expected_tolerance = 1e-4;

// Simple macro for single use inside a scope to keep the exception error message.
// This is a substitute for EXPECT_NO_THROW to keep the error message but don't
// necessarily fail.
#define SCOPE_TRACE_EXCEPTION_CONCAT_IMPL(x, y) x##y
#define SCOPE_TRACE_EXCEPTION_CONCAT(x, y) \
  SCOPE_TRACE_EXCEPTION_CONCAT_IMPL(x, y)
#define SCOPE_TRACE_EXCEPTION(EXPR)                \
  std::string scope_trace_msg = "Succeeded!";      \
  try {                                            \
    EXPR;                                          \
  } catch (const std::exception& e) {              \
    scope_trace_msg = e.what();                    \
  }                                                \
  ::testing::ScopedTrace trace(__FILE__, __LINE__, \
                               "Optimizer finished with: " + scope_trace_msg);

class NLPConstrainedProblemsTest : public ::testing::Test {
 protected:
  NLPConstrainedProblemsTest();

  std::string GetProblemDescription(int i) {
    std::stringstream ss;
    vect_n<double> x = funcs_start[i];
    int M = funcs_g[i](x).size();
    int K = funcs_h[i](x).size();
    ss << "Function #" << i << " called: '" << funcs_name[i] << "' with "
       << x.size() << " variables, " << M << " equality constraints, and " << K
       << " inequality constraints." << std::endl;
    return ss.str();
  }

  void CheckGradientProjection(const vect_n<double>& x, int i,
                               double tolerance);

  bool IsVectorFinite(const vect_n<double>& x) {
    return std::all_of(x.begin(), x.end(),
                       [](double x) { return std::isfinite(x); });
  }

  using FunctionPtr = std::function<double(const vect_n<double>&)>;
  using GradFunctionPtr = std::function<vect_n<double>(const vect_n<double>&)>;
  using HessFunctionPtr =
      std::function<void(mat<double, mat_structure::symmetric>&,
                         const vect_n<double>&, double, const vect_n<double>&)>;
  using GFunctionPtr = std::function<vect_n<double>(const vect_n<double>&)>;
  using GJacFunctionPtr =
      std::function<void(mat<double, mat_structure::rectangular>&,
                         const vect_n<double>&, const vect_n<double>&)>;
  using HFunctionPtr = std::function<vect_n<double>(const vect_n<double>&)>;
  using HJacFunctionPtr =
      std::function<void(mat<double, mat_structure::rectangular>&,
                         const vect_n<double>&, const vect_n<double>&)>;

  int eval_count = 0;
  int grad_count = 0;

  std::vector<std::string> funcs_name;
  std::vector<FunctionPtr> funcs_f;
  std::vector<GradFunctionPtr> funcs_grad;
  std::vector<HessFunctionPtr> funcs_hess;
  std::vector<GFunctionPtr> funcs_g;
  std::vector<GJacFunctionPtr> funcs_g_jac;
  std::vector<HFunctionPtr> funcs_h;
  std::vector<HJacFunctionPtr> funcs_h_jac;
  std::vector<vect_n<double>> funcs_sol;
  std::vector<vect_n<double>> funcs_start;
  std::vector<vect_n<double>> funcs_lower;
  std::vector<vect_n<double>> funcs_upper;
};

void NLPConstrainedProblemsTest::CheckGradientProjection(
    const vect_n<double>& x, int i, double tolerance) {
  using std::abs;
  const int N = x.size();
  vect_n<double> f_grad = funcs_grad[i](x);
  vect_n<double> g_value = funcs_g[i](x);
  mat<double, mat_structure::rectangular> g_jac;
  funcs_g_jac[i](g_jac, x, g_value);
  vect_n<double> h_value = funcs_h[i](x);
  mat<double, mat_structure::rectangular> h_jac;
  funcs_h_jac[i](h_jac, x, h_value);
  mat<double, mat_structure::rectangular> h_jac_reduced;
  h_jac_reduced.set_col_count(N);
  for (int j = 0, k = 0; j < h_jac.get_row_count(); ++j) {
    if (abs(h_value[j]) < tolerance) {
      // Constraint j is active.
      h_jac_reduced.set_row_count(k + 1, true);
      sub(h_jac_reduced)(range(k, k + 1), range(0, N)) =
          sub(h_jac)(range(j, j + 1), range(0, N));
      ++k;
    }
  }
  vect_n<double> f_grad_remainder = f_grad;
  mat<double, mat_structure::rectangular> gh_jac_t(
      N, g_jac.get_row_count() + h_jac_reduced.get_row_count());
  sub(gh_jac_t)(range(0, N), range(0, g_jac.get_row_count())) =
      transpose_view(g_jac);
  sub(gh_jac_t)(range(0, N),
                range(g_jac.get_row_count(),
                      g_jac.get_row_count() + h_jac_reduced.get_row_count())) =
      transpose_view(h_jac_reduced);
  if (gh_jac_t.get_col_count() > 0) {
    mat<double, mat_structure::rectangular> gh_jac_t_pinv(
        gh_jac_t.get_col_count(), gh_jac_t.get_row_count());
    pseudoinvert_SVD(gh_jac_t, gh_jac_t_pinv, tolerance * 1e-3);
    f_grad_remainder -= gh_jac_t * (gh_jac_t_pinv * f_grad);
  }
  EXPECT_THAT(f_grad_remainder, VectorIsZero(10 * tolerance))
      << GetProblemDescription(i);
}

NLPConstrainedProblemsTest::NLPConstrainedProblemsTest() {
  using std::atan;
  using std::cos;
  using std::log;
  using std::pow;
  using std::sin;
  using std::sqrt;
  /**********************************************************************/
  funcs_name.emplace_back("Banana Function");
  // Objective function.
  funcs_f.emplace_back([this](const vect_n<double>& x) {
    ++(this->eval_count);
    return 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) +
           (1.0 - x[0]) * (1.0 - x[0]);
  });
  funcs_grad.emplace_back([this](const vect_n<double>& x) {
    ++(this->grad_count);
    vect_n<double> result(2);
    result[0] = -400.0 * x[0] * (x[1] - x[0] * x[0]) - 2.0 * (1.0 - x[0]);
    result[1] = 200.0 * (x[1] - x[0] * x[0]);
    return result;
  });
  funcs_hess.emplace_back([](mat<double, mat_structure::symmetric>& H,
                             const vect_n<double>& x, double f,
                             const vect_n<double>& x_grad) {
    H.set_col_count(2);
    H(0, 0) = 1200.0 * x[0] * x[0] - 400.0 * x[1] + 2.0;
    H(0, 1) = -400.0 * x[0];
    H(1, 1) = 200.0;
  });
  // Equality constraints.
  funcs_g.emplace_back(
      [](const vect_n<double>& /*unused*/) { return vect_n<double>(0); });
  funcs_g_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(0);
  });
  // Inequality constraints.
  funcs_h.emplace_back([](const vect_n<double>& x) {
    vect_n<double> result(1);
    result[0] = x[1] + 1.5;
    return result;
  });
  funcs_h_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(1);
    J(0, 0) = 0.0;
    J(0, 1) = 1.0;
  });
  // Variables.
  funcs_sol.emplace_back(1.0, 1.0);
  funcs_start.emplace_back(-2.0, 1.0);
  funcs_lower.emplace_back(-std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity());
  funcs_upper.emplace_back(std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity());

  /**********************************************************************/
  if (funcs_name.back() == "Banana Function") {
    funcs_name.emplace_back("Banana Function With Inequality Constraints");
    funcs_f.push_back(funcs_f.back());
    funcs_grad.push_back(funcs_grad.back());
    funcs_hess.push_back(funcs_hess.back());
    funcs_g.push_back(funcs_g.back());
    funcs_g_jac.push_back(funcs_g_jac.back());
    // Inequality constraints.
    funcs_h.emplace_back([](const vect_n<double>& x) {
      vect_n<double> result(1);
      result[0] = x[1] - 1.5;  // Make (1,1) solution infeasible.
      return result;
    });
    funcs_h_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                                const vect_n<double>& /*unused*/,
                                const vect_n<double>& /*unused*/) {
      J.set_col_count(2);
      J.set_row_count(1);
      J(0, 0) = 0.0;
      J(0, 1) = 1.0;
    });
    // Variables.
    funcs_sol.emplace_back(-1.2210274, 1.5);
    funcs_start.push_back(funcs_start.back());
    funcs_lower.push_back(funcs_lower.back());
    funcs_upper.push_back(funcs_upper.back());
  }

  /**********************************************************************/
  funcs_name.emplace_back("Function 3");
  // Objective function.
  funcs_f.emplace_back([this](const vect_n<double>& x) {
    ++(this->eval_count);
    return x[1] + 1.0e-5 * (x[1] - x[0]) * (x[1] - x[0]);
  });
  funcs_grad.emplace_back([this](const vect_n<double>& x) {
    ++(this->grad_count);
    vect_n<double> result(2);
    result[0] = -2.0e-5 * (x[1] - x[0]);
    result[1] = 1.0 - result[0];
    return result;
  });
  funcs_hess.emplace_back([](mat<double, mat_structure::symmetric>& H,
                             const vect_n<double>& x, double f,
                             const vect_n<double>& x_grad) {
    H.set_col_count(2);
    H(0, 0) = 2.0e-5;
    H(0, 1) = -2.0e-5;
    H(1, 1) = 2.0;
  });
  // Equality constraints.
  funcs_g.emplace_back(
      [](const vect_n<double>& /*unused*/) { return vect_n<double>(0); });
  funcs_g_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(0);
  });
  // Inequality constraints.
  funcs_h.emplace_back([](const vect_n<double>& x) {
    vect_n<double> result(1);
    result[0] = x[1];
    return result;
  });
  funcs_h_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(1);
    J(0, 0) = 0.0;
    J(0, 1) = 1.0;
  });
  // Variables.
  funcs_sol.emplace_back(0.0, 0.0);
  funcs_start.emplace_back(10.0, 1.0);
  funcs_lower.emplace_back(-std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity());
  funcs_upper.emplace_back(std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity());

  /**********************************************************************/
  funcs_name.emplace_back("Function 4");
  // Objective function.
  funcs_f.emplace_back([this](const vect_n<double>& x) {
    ++(this->eval_count);
    return (x[0] + 1.0) * (x[0] + 1.0) * (x[0] + 1.0) / 3.0 + x[1];
  });
  funcs_grad.emplace_back([this](const vect_n<double>& x) {
    ++(this->grad_count);
    vect_n<double> result(2);
    result[0] = (x[0] + 1.0) * (x[0] + 1.0);
    result[1] = 1.0;
    return result;
  });
  funcs_hess.emplace_back([](mat<double, mat_structure::symmetric>& H,
                             const vect_n<double>& x, double f,
                             const vect_n<double>& x_grad) {
    H.set_col_count(2);
    H(0, 0) = 2.0 * (x[0] + 1.0);
    H(0, 1) = 0.0;
    H(1, 1) = 0.0;
  });
  // Equality constraints.
  funcs_g.emplace_back(
      [](const vect_n<double>& /*unused*/) { return vect_n<double>(0); });
  funcs_g_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(0);
  });
  // Inequality constraints.
  funcs_h.emplace_back([](const vect_n<double>& x) {
    vect_n<double> result(2);
    result[0] = x[0] - 1.0;
    result[1] = x[1];
    return result;
  });
  funcs_h_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(2);
    J(0, 0) = 1.0;
    J(0, 1) = 0.0;
    J(1, 0) = 0.0;
    J(1, 1) = 1.0;
  });
  // Variables.
  funcs_sol.emplace_back(1.0, 0.0);
  funcs_start.emplace_back(1.125, 0.125);
  funcs_lower.emplace_back(-std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity());
  funcs_upper.emplace_back(std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity());

  /**********************************************************************/
  funcs_name.emplace_back("Function 5");
  // Objective function.
  funcs_f.emplace_back([this](const vect_n<double>& x) {
    ++(this->eval_count);
    return sin(x[0] + x[1]) + (x[0] - x[1]) * (x[0] - x[1]) - 1.5 * x[0] +
           2.5 * x[1] + 1.0;
  });
  funcs_grad.emplace_back([this](const vect_n<double>& x) {
    ++(this->grad_count);
    vect_n<double> result(2);
    double v1 = cos(x[0] + x[1]);
    double v2 = 2.0 * (x[0] - x[1]);
    result[0] = v1 + v2 - 1.5;
    result[1] = v1 - v2 + 2.5;
    return result;
  });
  funcs_hess.emplace_back([](mat<double, mat_structure::symmetric>& H,
                             const vect_n<double>& x, double f,
                             const vect_n<double>& x_grad) {
    H.set_col_count(2);
    double v1 = -sin(x[0] + x[1]);
    H(0, 0) = v1 + 2.0;
    H(0, 1) = v1 - 2.0;
    H(1, 1) = v1 + 2.0;
  });
  // Equality constraints.
  funcs_g.emplace_back(
      [](const vect_n<double>& /*unused*/) { return vect_n<double>(0); });
  funcs_g_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(0);
  });
  // Inequality constraints.
  funcs_h.emplace_back([](const vect_n<double>& x) {
    vect_n<double> result(4);
    result[0] = x[0] + 1.5;
    result[1] = x[1] + 3.0;
    result[2] = 4.0 - x[0];
    result[3] = 3.0 - x[1];
    return result;
  });
  funcs_h_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(4);
    J(0, 0) = 1.0;
    J(0, 1) = 0.0;
    J(1, 0) = 0.0;
    J(1, 1) = 1.0;
    J(2, 0) = -1.0;
    J(2, 1) = 0.0;
    J(3, 0) = 0.0;
    J(3, 1) = -1.0;
  });
  // Variables.
  funcs_sol.emplace_back(0.5 - 4.0 * atan(1.0) / 3.0,
                         0.5 - 4.0 * atan(1.0) / 3.0 - 1.0);
  funcs_start.emplace_back(0.0, 0.0);
  funcs_lower.emplace_back(-std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity());
  funcs_upper.emplace_back(std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity());

  /**********************************************************************/
  funcs_name.emplace_back("Function 6");
  // Objective function.
  funcs_f.emplace_back([this](const vect_n<double>& x) {
    ++(this->eval_count);
    return (1.0 - x[0]) * (1.0 - x[0]);
  });
  funcs_grad.emplace_back([this](const vect_n<double>& x) {
    ++(this->grad_count);
    vect_n<double> result(2);
    result[0] = -2.0 * (1.0 - x[0]);
    result[1] = 0.0;
    return result;
  });
  funcs_hess.emplace_back([](mat<double, mat_structure::symmetric>& H,
                             const vect_n<double>& x, double f,
                             const vect_n<double>& x_grad) {
    H.set_col_count(2);
    H(0, 0) = 2.0;
    H(0, 1) = 0.0;
    H(1, 1) = 0.0;
  });
  // Equality constraints.
  funcs_g.emplace_back([](const vect_n<double>& x) {
    vect_n<double> result(1);
    result[0] = 10.0 * (x[1] - x[0] * x[0]);
    return result;
  });
  funcs_g_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& x,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(1);
    J(0, 0) = -20.0 * x[0];
    J(0, 1) = 10.0;
  });
  // Inequality constraints.
  funcs_h.emplace_back(
      [](const vect_n<double>& x) { return vect_n<double>(0); });
  funcs_h_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(0);
  });
  // Variables.
  funcs_sol.emplace_back(1.0, 1.0);
  funcs_start.emplace_back(-1.2, 1.0);
  funcs_lower.emplace_back(-std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity());
  funcs_upper.emplace_back(std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity());

  /**********************************************************************/
  funcs_name.emplace_back("Function 7");
  // Objective function.
  funcs_f.emplace_back([this](const vect_n<double>& x) {
    ++(this->eval_count);
    return log(1.0 + x[0] * x[0]) - x[1];
  });
  funcs_grad.emplace_back([this](const vect_n<double>& x) {
    ++(this->grad_count);
    vect_n<double> result(2);
    result[0] = 2.0 * x[0] / (1.0 + x[0] * x[0]);
    result[1] = -1.0;
    return result;
  });
  funcs_hess.emplace_back([](mat<double, mat_structure::symmetric>& H,
                             const vect_n<double>& x, double f,
                             const vect_n<double>& x_grad) {
    H.set_col_count(2);
    double x2 = x[0] * x[0];
    H(0, 0) = 2.0 / (1.0 + x2) - 4.0 * x2 / ((1.0 + x2) * (1.0 + x2));
    H(0, 1) = 0.0;
    H(1, 1) = 0.0;
  });
  // Equality constraints.
  funcs_g.emplace_back([](const vect_n<double>& x) {
    vect_n<double> result(1);
    result[0] = (1.0 + x[0] * x[0]) * (1.0 + x[0] * x[0]) + x[1] * x[1] - 4.0;
    return result;
  });
  funcs_g_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& x,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(1);
    J(0, 0) = 4.0 * x[0] * (1.0 + x[0] * x[0]);
    J(0, 1) = 2.0 * x[1];
  });
  // Inequality constraints.
  funcs_h.emplace_back(
      [](const vect_n<double>& x) { return vect_n<double>(0); });
  funcs_h_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(0);
  });
  // Variables.
  funcs_sol.emplace_back(0.0, sqrt(3.0));
  funcs_start.emplace_back(2.0, 2.0);
  funcs_lower.emplace_back(-std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity());
  funcs_upper.emplace_back(std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity());

  /**********************************************************************/
  funcs_name.emplace_back("Function 10");
  // Objective function.
  funcs_f.emplace_back([this](const vect_n<double>& x) {
    ++(this->eval_count);
    return x[0] - x[1];
  });
  funcs_grad.emplace_back([this](const vect_n<double>& x) {
    ++(this->grad_count);
    vect_n<double> result(2);
    result[0] = 1.0;
    result[1] = -1.0;
    return result;
  });
  funcs_hess.emplace_back([](mat<double, mat_structure::symmetric>& H,
                             const vect_n<double>& x, double f,
                             const vect_n<double>& x_grad) {
    H.set_col_count(2);
    H(0, 0) = 0.0;
    H(0, 1) = 0.0;
    H(1, 1) = 0.0;
  });
  // Equality constraints.
  funcs_g.emplace_back(
      [](const vect_n<double>& /*unused*/) { return vect_n<double>(0); });
  funcs_g_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(0);
  });
  // Inequality constraints.
  funcs_h.emplace_back([](const vect_n<double>& x) {
    vect_n<double> result(1);
    result[0] = -3.0 * x[0] * x[0] + 2.0 * x[0] * x[1] - x[1] * x[1] + 1.0;
    return result;
  });
  funcs_h_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& x,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(1);
    J(0, 0) = -6.0 * x[0] + 2.0 * x[1];
    J(0, 1) = 2.0 * (x[0] - x[1]);
  });
  // Variables.
  funcs_sol.emplace_back(0.0, 1.0);
  funcs_start.emplace_back(-10.0, 10.0);
  funcs_lower.emplace_back(-std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity());
  funcs_upper.emplace_back(std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity());

  /**********************************************************************/
  funcs_name.emplace_back("Function 11");
  // Objective function.
  funcs_f.emplace_back([this](const vect_n<double>& x) {
    ++(this->eval_count);
    return (x[0] - 5.0) * (x[0] - 5.0) + x[1] * x[1] - 25.0;
  });
  funcs_grad.emplace_back([this](const vect_n<double>& x) {
    ++(this->grad_count);
    vect_n<double> result(2);
    result[0] = 2.0 * (x[0] - 5.0);
    result[1] = 2.0 * x[1];
    return result;
  });
  funcs_hess.emplace_back([](mat<double, mat_structure::symmetric>& H,
                             const vect_n<double>& x, double f,
                             const vect_n<double>& x_grad) {
    H.set_col_count(2);
    H(0, 0) = 2.0;
    H(0, 1) = 0.0;
    H(1, 1) = 2.0;
  });
  // Equality constraints.
  funcs_g.emplace_back(
      [](const vect_n<double>& /*unused*/) { return vect_n<double>(0); });
  funcs_g_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(0);
  });
  // Inequality constraints.
  funcs_h.emplace_back([](const vect_n<double>& x) {
    vect_n<double> result(1);
    result[0] = x[1] - x[0] * x[0];
    return result;
  });
  funcs_h_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& x,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(1);
    J(0, 0) = -2.0 * x[0];
    J(0, 1) = 1.0;
  });
  // Variables.
  auto p11_get_sol = []() {
    vect_n<double> result(2);
    double AEX = 7.5 * sqrt(6.0);
    double AW = pow(sqrt(AEX * AEX + 1.0) + AEX, 1.0 / 3.0);
    double QAW = AW * AW;
    result[0] = (AW - 1.0 / AW) / sqrt(6.0);
    result[1] = (QAW - 2.0 + 1.0 / QAW) / 6.0;
    return result;
  };
  funcs_sol.emplace_back(p11_get_sol());
  funcs_start.emplace_back(4.9, 0.1);
  funcs_lower.emplace_back(-std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity());
  funcs_upper.emplace_back(std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity());

  /**********************************************************************/
  funcs_name.emplace_back("Function 12");
  // Objective function.
  funcs_f.emplace_back([this](const vect_n<double>& x) {
    ++(this->eval_count);
    return 0.5 * x[0] * x[0] + x[1] * x[1] - x[0] * x[1] - 7.0 * x[0] -
           7.0 * x[1];
  });
  funcs_grad.emplace_back([this](const vect_n<double>& x) {
    ++(this->grad_count);
    vect_n<double> result(2);
    result[0] = x[0] - x[1] - 7.0;
    result[1] = 2.0 * x[1] - x[0] - 7.0;
    return result;
  });
  funcs_hess.emplace_back([](mat<double, mat_structure::symmetric>& H,
                             const vect_n<double>& x, double f,
                             const vect_n<double>& x_grad) {
    H.set_col_count(2);
    H(0, 0) = 1.0;
    H(0, 1) = -1.0;
    H(1, 1) = 2.0;
  });
  // Equality constraints.
  funcs_g.emplace_back(
      [](const vect_n<double>& /*unused*/) { return vect_n<double>(0); });
  funcs_g_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(0);
  });
  // Inequality constraints.
  funcs_h.emplace_back([](const vect_n<double>& x) {
    vect_n<double> result(1);
    result[0] = 25.0 - 4.0 * x[0] * x[0] - x[1] * x[1];
    return result;
  });
  funcs_h_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& x,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(1);
    J(0, 0) = -8.0 * x[0];
    J(0, 1) = -2.0 * x[1];
  });
  // Variables.
  funcs_sol.emplace_back(2.0, 3.0);
  funcs_start.emplace_back(0.0, 0.0);
  funcs_lower.emplace_back(-std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity());
  funcs_upper.emplace_back(std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity());

  /**********************************************************************/
  funcs_name.emplace_back("Function 14");
  // Objective function.
  funcs_f.emplace_back([this](const vect_n<double>& x) {
    ++(this->eval_count);
    return (x[0] - 2.0) * (x[0] - 2.0) + (x[1] - 1.0) * (x[1] - 1.0);
  });
  funcs_grad.emplace_back([this](const vect_n<double>& x) {
    ++(this->grad_count);
    vect_n<double> result(2);
    result[0] = 2.0 * (x[0] - 2.0);
    result[1] = 2.0 * (x[1] - 1.0);
    return result;
  });
  funcs_hess.emplace_back([](mat<double, mat_structure::symmetric>& H,
                             const vect_n<double>& x, double f,
                             const vect_n<double>& x_grad) {
    H.set_col_count(2);
    H(0, 0) = 2.0;
    H(0, 1) = 0.0;
    H(1, 1) = 2.0;
  });
  // Equality constraints.
  funcs_g.emplace_back([](const vect_n<double>& x) {
    vect_n<double> result(1);
    result[0] = x[0] - 2.0 * x[1] + 1.0;
    return result;
  });
  funcs_g_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(1);
    J(0, 0) = 1.0;
    J(0, 1) = -2.0;
  });
  // Inequality constraints.
  funcs_h.emplace_back([](const vect_n<double>& x) {
    vect_n<double> result(1);
    result[0] = 1.0 - x[0] * x[0] * 0.25 - x[1] * x[1];
    return result;
  });
  funcs_h_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& x,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(1);
    J(0, 0) = -0.5 * x[0];
    J(0, 1) = -2.0 * x[1];
  });
  // Variables.
  funcs_sol.emplace_back((sqrt(7.0) - 1.0) * 0.5, (sqrt(7.0) + 1.0) * 0.25);
  funcs_start.emplace_back(2.0, 2.0);
  funcs_lower.emplace_back(-std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity());
  funcs_upper.emplace_back(std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity());

  /**********************************************************************/
  funcs_name.emplace_back("Function 15");
  // Objective function.
  funcs_f.emplace_back([this](const vect_n<double>& x) {
    ++(this->eval_count);
    return (x[0] - 2.0) * (x[0] - 2.0) + (x[1] - 1.0) * (x[1] - 1.0);
  });
  funcs_grad.emplace_back([this](const vect_n<double>& x) {
    ++(this->grad_count);
    vect_n<double> result(2);
    result[0] = 2.0 * (x[0] - 2.0);
    result[1] = 2.0 * (x[1] - 1.0);
    return result;
  });
  funcs_hess.emplace_back([](mat<double, mat_structure::symmetric>& H,
                             const vect_n<double>& x, double f,
                             const vect_n<double>& x_grad) {
    H.set_col_count(2);
    H(0, 0) = 2.0;
    H(0, 1) = 0.0;
    H(1, 1) = 2.0;
  });
  // Equality constraints.
  funcs_g.emplace_back(
      [](const vect_n<double>& x) { return vect_n<double>(0); });
  funcs_g_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& /*unused*/,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(0);
  });
  // Inequality constraints.
  funcs_h.emplace_back([](const vect_n<double>& x) {
    vect_n<double> result(1);
    result[0] = 1.0 - x[0] * x[0] * 0.25 - x[1] * x[1];
    return result;
  });
  funcs_h_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& x,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(2);
    J.set_row_count(1);
    J(0, 0) = -0.5 * x[0];
    J(0, 1) = -2.0 * x[1];
  });
  // Variables.
  funcs_sol.emplace_back(1.66497, 0.55405);
  funcs_start.emplace_back(2.0, 2.0);
  funcs_lower.emplace_back(-std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity());
  funcs_upper.emplace_back(std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::infinity());

  /**********************************************************************/
  funcs_name.emplace_back("Robot 1");
  // Objective function.
  funcs_f.emplace_back([this](const vect_n<double>& x) {
    ++(this->eval_count);
    return x[4] * x[4] + x[5] * x[5] + x[6] * x[6] + x[7] * x[7] +
           0.001 * cos(x[2]) * cos(x[2]) + 0.001 * cos(x[3]) * cos(x[3]);
  });
  funcs_grad.emplace_back([this](const vect_n<double>& x) {
    ++(this->grad_count);
    vect_n<double> result(8);
    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = -2.0 * 0.001 * cos(x[2]) * sin(x[2]);
    result[3] = -2.0 * 0.001 * cos(x[3]) * sin(x[3]);
    result[4] = 2 * x[4];
    result[5] = 2 * x[5];
    result[6] = 2 * x[6];
    result[7] = 2 * x[7];
    return result;
  });
  funcs_hess.emplace_back([](mat<double, mat_structure::symmetric>& H,
                             const vect_n<double>& x, double f,
                             const vect_n<double>& x_grad) {
    H.set_col_count(8);
    H = mat<double, mat_structure::nil>(8, 8);
    H(0, 0) = 0.0;
    H(1, 1) = 0.0;
    H(2, 2) = 2.0 * 0.001 * sin(x[2]) * sin(x[2]) -
              2.0 * 0.001 * cos(x[2]) * cos(x[2]);
    H(3, 3) = 2.0 * 0.001 * sin(x[3]) * sin(x[3]) -
              2.0 * 0.001 * cos(x[3]) * cos(x[3]);
    H(4, 4) = 2.0;
    H(5, 5) = 2.0;
    H(6, 6) = 2.0;
    H(7, 7) = 2.0;
  });
  // Equality constraints.
  auto get_robot1_lengths = []() {
    vect_n<double> result(4);
    result[0] = 3.0;
    result[1] = 0.35;
    result[2] = 0.35;
    result[3] = 0.1;
    return result;
  };
  auto get_robot1_desired = []() {
    vect_n<double> result(6);
    result[0] = 2.0;
    result[1] = 0.5;
    result[2] = M_PI * 0.5;
    result[3] = 0.1;
    result[4] = 0.1;
    result[5] = 0.1;
    return result;
  };
  funcs_g.emplace_back([get_robot1_lengths,
                        get_robot1_desired](const vect_n<double>& x) {
    vect_n<double> result(6);
    vect_n<double> l = get_robot1_lengths();
    double s1 = l[1] * sin(x[1]);
    double s12 = l[2] * sin(x[1] + x[2]);
    double s123 = l[3] * sin(x[1] + x[2] + x[3]);
    double c1 = l[1] * cos(x[1]);
    double c12 = l[2] * cos(x[1] + x[2]);
    double c123 = l[3] * cos(x[1] + x[2] + x[3]);
    result = get_robot1_desired();
    result[0] -= l[0] * x[0] + c1 + c12 + c123;
    result[1] -= s1 + s12 + s123;
    result[2] -= x[1] + x[2] + x[3];
    result[3] -= l[0] * x[4] - s1 * x[5] - s12 * (x[5] + x[6]) -
                 s123 * (x[5] + x[6] + x[7]);
    result[4] -= c1 * x[5] + c12 * (x[5] + x[6]) + c123 * (x[5] + x[6] + x[7]);
    result[5] -= x[5] + x[6] + x[7];
    return result;
  });
  funcs_g_jac.emplace_back(
      [get_robot1_lengths](mat<double, mat_structure::rectangular>& J,
                           const vect_n<double>& x,
                           const vect_n<double>& /*unused*/) {
        J.set_col_count(8);
        J.set_row_count(6);
        vect_n<double> l = get_robot1_lengths();
        double s1 = l[1] * sin(x[1]);
        double s12 = l[2] * sin(x[1] + x[2]);
        double s123 = l[3] * sin(x[1] + x[2] + x[3]);
        double c1 = l[1] * cos(x[1]);
        double c12 = l[2] * cos(x[1] + x[2]);
        double c123 = l[3] * cos(x[1] + x[2] + x[3]);
        J(0, 0) = -l[0];
        J(0, 1) = s1 + s12 + s123;
        J(0, 2) = s12 + s123;
        J(0, 3) = s123;
        J(0, 4) = 0.0;
        J(0, 5) = 0.0;
        J(0, 6) = 0.0;
        J(0, 7) = 0.0;
        J(1, 0) = 0.0;
        J(1, 1) = -c1 - c12 - c123;
        J(1, 2) = -c12 - c123;
        J(1, 3) = -c123;
        J(1, 4) = 0.0;
        J(1, 5) = 0.0;
        J(1, 6) = 0.0;
        J(1, 7) = 0.0;
        J(2, 0) = 0.0;
        J(2, 1) = -1.0;
        J(2, 2) = -1.0;
        J(2, 3) = -1.0;
        J(2, 4) = 0.0;
        J(2, 5) = 0.0;
        J(2, 6) = 0.0;
        J(2, 7) = 0.0;

        s1 *= x[5];
        s12 *= (x[5] + x[6]);
        s123 *= (x[5] + x[6] + x[7]);
        c1 *= (x[5]);
        c12 *= (x[5] + x[6]);
        c123 *= (x[5] + x[6] + x[7]);
        J(3, 0) = 0.0;
        J(3, 1) = c1 + c12 + c123;
        J(3, 2) = c12 + c123;
        J(3, 3) = c123;
        J(3, 4) = J(0, 0);
        J(3, 5) = J(0, 1);
        J(3, 6) = J(0, 2);
        J(3, 7) = J(0, 3);
        J(4, 0) = 0.0;
        J(4, 1) = s1 + s12 + s123;
        J(4, 2) = s12 + s123;
        J(4, 3) = s123;
        J(4, 4) = J(1, 0);
        J(4, 5) = J(1, 1);
        J(4, 6) = J(1, 2);
        J(4, 7) = J(1, 3);
        J(5, 0) = 0.0;
        J(5, 1) = 0.0;
        J(5, 2) = 0.0;
        J(5, 3) = 0.0;
        J(5, 4) = J(2, 0);
        J(5, 5) = J(2, 1);
        J(5, 6) = J(2, 2);
        J(5, 7) = J(2, 3);
      });
  // Inequality constraints.
  funcs_h.emplace_back([](const vect_n<double>& x) {
    vect_n<double> result(2);
    result[0] = x[0];
    result[1] = 1.0 - x[0];
    return result;
  });
  funcs_h_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& x,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(8);
    J.set_row_count(2);
    J(0, 0) = 1.0;
    J(0, 1) = 0.0;
    J(0, 2) = 0.0;
    J(0, 3) = 0.0;
    J(0, 4) = 0.0;
    J(0, 5) = 0.0;
    J(0, 6) = 0.0;
    J(0, 7) = 0.0;
    J(1, 0) = -1.0;
    J(1, 1) = 0.0;
    J(1, 2) = 0.0;
    J(1, 3) = 0.0;
    J(1, 4) = 0.0;
    J(1, 5) = 0.0;
    J(1, 6) = 0.0;
    J(1, 7) = 0.0;
  });
  // Variables.
  // Multiple solutions, so only check gradient projection.
  funcs_sol.emplace_back(8, std::numeric_limits<double>::quiet_NaN());
  funcs_start.emplace_back(0.66, M_PI * 0.35, M_PI * 0.35, M_PI * 0.35, 0.03,
                           0.03, 0.03, 0.03);
  funcs_lower.emplace_back(8, -std::numeric_limits<double>::infinity());
  funcs_upper.emplace_back(8, std::numeric_limits<double>::infinity());

  /**********************************************************************/
  funcs_name.emplace_back("Robot 2");
  // Objective function.
  funcs_f.emplace_back([this](const vect_n<double>& x) {
    ++(this->eval_count);
    return cos(x[2]) * cos(x[2]) + cos(x[3]) * cos(x[3]);
  });
  funcs_grad.emplace_back([this](const vect_n<double>& x) {
    ++(this->grad_count);
    vect_n<double> result(4);
    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = -2.0 * cos(x[2]) * sin(x[2]);
    result[3] = -2.0 * cos(x[3]) * sin(x[3]);
    return result;
  });
  funcs_hess.emplace_back([](mat<double, mat_structure::symmetric>& H,
                             const vect_n<double>& x, double f,
                             const vect_n<double>& x_grad) {
    H.set_col_count(4);
    H = mat<double, mat_structure::nil>(4, 4);
    H(0, 0) = 0.0;
    H(1, 1) = 0.0;
    H(2, 2) = 2.0 * sin(x[2]) * sin(x[2]) - 2.0 * cos(x[2]) * cos(x[2]);
    H(3, 3) = 2.0 * sin(x[3]) * sin(x[3]) - 2.0 * cos(x[3]) * cos(x[3]);
  });
  // Equality constraints.
  auto get_robot2_lengths = []() {
    vect_n<double> result(4);
    result[0] = 3.0;
    result[1] = 0.35;
    result[2] = 0.35;
    result[3] = 0.1;
    return result;
  };
  auto get_robot2_desired = []() {
    vect_n<double> result(3);
    result[0] = 2.0;
    result[1] = 0.5;
    result[2] = M_PI * 0.5;
    return result;
  };
  funcs_g.emplace_back(
      [get_robot2_lengths, get_robot2_desired](const vect_n<double>& x) {
        vect_n<double> result(3);
        vect_n<double> l = get_robot2_lengths();
        double s1 = l[1] * sin(x[1]);
        double s12 = l[2] * sin(x[1] + x[2]);
        double s123 = l[3] * sin(x[1] + x[2] + x[3]);
        double c1 = l[1] * cos(x[1]);
        double c12 = l[2] * cos(x[1] + x[2]);
        double c123 = l[3] * cos(x[1] + x[2] + x[3]);
        result = get_robot2_desired();
        result[0] -= l[0] * x[0] + c1 + c12 + c123;
        result[1] -= s1 + s12 + s123;
        result[2] -= x[1] + x[2] + x[3];
        result[0] /= l[0];
        result[1] /= (l[1] + l[2] + l[3]);
        result[2] /= M_PI;
        return result;
      });
  funcs_g_jac.emplace_back(
      [get_robot2_lengths](mat<double, mat_structure::rectangular>& J,
                           const vect_n<double>& x,
                           const vect_n<double>& /*unused*/) {
        J.set_col_count(4);
        J.set_row_count(3);
        vect_n<double> l = get_robot2_lengths();
        double s1 = l[1] * sin(x[1]);
        double s12 = l[2] * sin(x[1] + x[2]);
        double s123 = l[3] * sin(x[1] + x[2] + x[3]);
        double c1 = l[1] * cos(x[1]);
        double c12 = l[2] * cos(x[1] + x[2]);
        double c123 = l[3] * cos(x[1] + x[2] + x[3]);
        J(0, 0) = -l[0];
        J(0, 1) = s1 + s12 + s123;
        J(0, 2) = s12 + s123;
        J(0, 3) = s123;
        J(1, 0) = 0.0;
        J(1, 1) = -c1 - c12 - c123;
        J(1, 2) = -c12 - c123;
        J(1, 3) = -c123;
        J(2, 0) = 0.0;
        J(2, 1) = -1.0;
        J(2, 2) = -1.0;
        J(2, 3) = -1.0;

        J(0, 0) /= l[0];
        J(0, 1) /= l[0];
        J(0, 2) /= l[0];
        J(0, 3) /= l[0];
        J(1, 0) /= l[1] + l[2] + l[3];
        J(1, 1) /= l[1] + l[2] + l[3];
        J(1, 2) /= l[1] + l[2] + l[3];
        J(1, 3) /= l[1] + l[2] + l[3];
        J(2, 0) /= M_PI;
        J(2, 1) /= M_PI;
        J(2, 2) /= M_PI;
        J(2, 3) /= M_PI;
      });
  // Inequality constraints.
  funcs_h.emplace_back([](const vect_n<double>& x) {
    vect_n<double> result(2);
    result[0] = x[0];
    result[1] = 1.0 - x[0];
    return result;
  });
  funcs_h_jac.emplace_back([](mat<double, mat_structure::rectangular>& J,
                              const vect_n<double>& x,
                              const vect_n<double>& /*unused*/) {
    J.set_col_count(4);
    J.set_row_count(2);
    J(0, 0) = 1.0;
    J(0, 1) = 0.0;
    J(0, 2) = 0.0;
    J(0, 3) = 0.0;
    J(1, 0) = -1.0;
    J(1, 1) = 0.0;
    J(1, 2) = 0.0;
    J(1, 3) = 0.0;
  });
  // Variables.
  // Multiple solutions, so only check gradient projection.
  funcs_sol.emplace_back(4, std::numeric_limits<double>::quiet_NaN());
  funcs_start.emplace_back(0.66, M_PI * 0.35, M_PI * 0.35, M_PI * 0.35);
  funcs_lower.emplace_back(4, -std::numeric_limits<double>::infinity());
  funcs_upper.emplace_back(4, std::numeric_limits<double>::infinity());
}

TEST_F(NLPConstrainedProblemsTest, TrustRegionNewtonMethodUnconstrained) {
  const double max_radius = 2.0;
  const int max_iter = 300;
  const double eta = 1e-3;
  std::unordered_set<std::string> expected_failures = {};
  vect_n<double> x;
  for (int i = 0; i < funcs_f.size(); ++i) {
    if (expected_failures.count(funcs_name[i]) == 1) {
      std::cout << "Skipping known failure: '" << funcs_name[i] << "'"
                << std::endl;
      continue;
    }
    x = funcs_start[i];
    if (!funcs_h[i](x).empty() || !funcs_g[i](x).empty()) {
      continue;  // Cannot handle constraints.
    }
    eval_count = 0;
    grad_count = 0;
    SCOPE_TRACE_EXCEPTION(make_newton_method_tr(
        funcs_f[i], funcs_grad[i], funcs_hess[i], max_radius, max_iter,
        desired_tolerance, desired_tolerance, eta)(x));
    if (IsVectorFinite(x)) {
      if (std::isfinite(funcs_sol[i][0])) {
        EXPECT_THAT(x, VectorIsNear(funcs_sol[i], expected_tolerance))
            << GetProblemDescription(i);
      }
      CheckGradientProjection(x, i, expected_tolerance);
      EXPECT_THAT(funcs_g[i](x), VectorIsZero(expected_tolerance))
          << GetProblemDescription(i);
      EXPECT_THAT(funcs_h[i](x),
                  VectorIsGreater(vect_n<double>(funcs_h[i](x).size(), 0.0),
                                  expected_tolerance))
          << GetProblemDescription(i);
    }
  }
}

TEST_F(NLPConstrainedProblemsTest, TrustRegionNewtonMethodConstrained) {
  const double max_radius = 2.0;
  const int max_iter = 1000;
  const double eta = 1e-3;
  // TODO This method fails at everything.
  std::unordered_set<std::string> expected_failures = {
      "Banana Function", "Banana Function With Inequality Constraints",
      "Function 3",      "Function 4",
      "Function 5",      "Function 6",
      "Function 7",      "Function 10",
      "Function 11",     "Function 12",
      "Function 14",     "Function 15",
      "Robot 1",         "Robot 2"};
  vect_n<double> x;
  for (int i = 0; i < funcs_f.size(); ++i) {
    if (expected_failures.count(funcs_name[i]) == 1) {
      std::cout << "Skipping known failure: '" << funcs_name[i] << "'"
                << std::endl;
      continue;
    }
    x = funcs_start[i];
    eval_count = 0;
    grad_count = 0;
    SCOPE_TRACE_EXCEPTION(
        make_constraint_newton_method_tr(
            funcs_f[i], funcs_grad[i], funcs_hess[i], max_radius, max_iter,
            desired_tolerance, desired_tolerance, eta)
            .set_limiter([l = funcs_lower[i], u = funcs_upper[i]](
                             const vect_n<double>& x, vect_n<double>& dx) {
              box_limit_function(x, dx, l, u);
            })
            .set_eq_constraints(funcs_g[i], funcs_g_jac[i])
            .set_ineq_constraints(funcs_h[i], funcs_h_jac[i])(x));
    if (IsVectorFinite(x)) {
      if (std::isfinite(funcs_sol[i][0])) {
        EXPECT_THAT(x, VectorIsNear(funcs_sol[i], expected_tolerance))
            << GetProblemDescription(i);
      }
      CheckGradientProjection(x, i, expected_tolerance);
      EXPECT_THAT(funcs_g[i](x), VectorIsZero(expected_tolerance))
          << GetProblemDescription(i);
      EXPECT_THAT(funcs_h[i](x),
                  VectorIsGreater(vect_n<double>(funcs_h[i](x).size(), 0.0),
                                  expected_tolerance))
          << GetProblemDescription(i);
    }
  }
}

TEST_F(NLPConstrainedProblemsTest,
       TrustRegionNewtonMethodConstrainedRegularized) {
  const double max_radius = 2.0;
  const int max_iter = 600;
  const double eta = 1e-3;
  const double regularization_value = 1e-8;
  std::unordered_set<std::string> expected_failures = {
      "Banana Function With Inequality Constraints",
      "Function 3",
      "Function 4",
      "Function 6",
      "Function 7",
      "Function 10",
      "Function 11",
      "Function 12",
      "Function 14",
      "Function 15",
      "Robot 1",
      "Robot 2"};
  vect_n<double> x;
  for (int i = 0; i < funcs_f.size(); ++i) {
    if (expected_failures.count(funcs_name[i]) == 1) {
      std::cout << "Skipping known failure: '" << funcs_name[i] << "'"
                << std::endl;
      continue;
    }
    x = funcs_start[i];
    eval_count = 0;
    grad_count = 0;
    SCOPE_TRACE_EXCEPTION(
        make_constraint_newton_method_tr(
            funcs_f[i], funcs_grad[i], funcs_hess[i], max_radius, max_iter,
            desired_tolerance, desired_tolerance, eta)
            .set_limiter([l = funcs_lower[i], u = funcs_upper[i]](
                             const vect_n<double>& x, vect_n<double>& dx) {
              box_limit_function(x, dx, l, u);
            })
            .set_eq_constraints(funcs_g[i], funcs_g_jac[i])
            .set_ineq_constraints(funcs_h[i], funcs_h_jac[i])
            .regularize(regularization_value)(x));
    if (IsVectorFinite(x)) {
      if (std::isfinite(funcs_sol[i][0])) {
        EXPECT_THAT(x, VectorIsNear(funcs_sol[i], expected_tolerance))
            << GetProblemDescription(i);
      }
      CheckGradientProjection(x, i, expected_tolerance);
      EXPECT_THAT(funcs_g[i](x), VectorIsZero(expected_tolerance))
          << GetProblemDescription(i);
      EXPECT_THAT(funcs_h[i](x),
                  VectorIsGreater(vect_n<double>(funcs_h[i](x).size(), 0.0),
                                  expected_tolerance))
          << GetProblemDescription(i);
    }
  }
}

TEST_F(NLPConstrainedProblemsTest, TrustRegionInteriorPointNewtonMethod) {
  const double max_radius = 1.0;
  const double mu = 0.1;
  const int max_iter = 300;
  const double eta = 1e-3;
  const double tau = 0.99;
  std::unordered_set<std::string> expected_failures = {
      "Banana Function", "Banana Function With Inequality Constraints",
      "Function 3",      "Function 4",
      "Function 5",      "Function 6",
      "Function 7",      "Function 10"};
  vect_n<double> x;
  for (int i = 0; i < funcs_f.size(); ++i) {
    if (expected_failures.count(funcs_name[i]) == 1) {
      std::cout << "Skipping known failure: '" << funcs_name[i] << "'"
                << std::endl;
      continue;
    }
    x = funcs_start[i];
    eval_count = 0;
    grad_count = 0;
    SCOPE_TRACE_EXCEPTION(
        make_nlip_newton_tr(funcs_f[i], funcs_grad[i], funcs_hess[i],
                            max_radius, mu, max_iter, desired_tolerance, eta,
                            tau)
            .set_limiter([l = funcs_lower[i], u = funcs_upper[i]](
                             const vect_n<double>& x, vect_n<double>& dx) {
              box_limit_function(x, dx, l, u);
            })
            .set_eq_constraints(funcs_g[i], funcs_g_jac[i])
            .set_ineq_constraints(funcs_h[i], funcs_h_jac[i])(x));
    if (IsVectorFinite(x)) {
      if (std::isfinite(funcs_sol[i][0])) {
        EXPECT_THAT(x, VectorIsNear(funcs_sol[i], expected_tolerance))
            << GetProblemDescription(i);
      }
      CheckGradientProjection(x, i, expected_tolerance);
      EXPECT_THAT(funcs_g[i](x), VectorIsZero(expected_tolerance))
          << GetProblemDescription(i);
      EXPECT_THAT(funcs_h[i](x),
                  VectorIsGreater(vect_n<double>(funcs_h[i](x).size(), 0.0),
                                  expected_tolerance))
          << GetProblemDescription(i);
    }
  }
}

TEST_F(NLPConstrainedProblemsTest, TrustRegionInteriorPointQuasiNewtonMethod) {
  const double max_radius = 1.0;
  const double mu = 0.1;
  const int max_iter = 300;
  const double eta = 1e-3;
  const double tau = 0.99;
  std::unordered_set<std::string> expected_failures = {
      "Banana Function", "Banana Function With Inequality Constraints",
      "Function 3",      "Function 4",
      "Function 5",      "Function 6",
      "Function 7"};
  vect_n<double> x;
  for (int i = 0; i < funcs_f.size(); ++i) {
    if (expected_failures.count(funcs_name[i]) == 1) {
      std::cout << "Skipping known failure: '" << funcs_name[i] << "'"
                << std::endl;
      continue;
    }
    x = funcs_start[i];
    eval_count = 0;
    grad_count = 0;
    SCOPE_TRACE_EXCEPTION(
        make_nlip_quasi_newton_tr(funcs_f[i], funcs_grad[i], max_radius, mu,
                                  max_iter, desired_tolerance, eta, tau)
            .set_limiter([l = funcs_lower[i], u = funcs_upper[i]](
                             const vect_n<double>& x, vect_n<double>& dx) {
              box_limit_function(x, dx, l, u);
            })
            .set_eq_constraints(funcs_g[i], funcs_g_jac[i])
            .set_ineq_constraints(funcs_h[i], funcs_h_jac[i])(x));
    if (IsVectorFinite(x)) {
      if (std::isfinite(funcs_sol[i][0])) {
        EXPECT_THAT(x, VectorIsNear(funcs_sol[i], expected_tolerance))
            << GetProblemDescription(i);
      }
      CheckGradientProjection(x, i, expected_tolerance);
      EXPECT_THAT(funcs_g[i](x), VectorIsZero(expected_tolerance))
          << GetProblemDescription(i);
      EXPECT_THAT(funcs_h[i](x),
                  VectorIsGreater(vect_n<double>(funcs_h[i](x).size(), 0.0),
                                  expected_tolerance))
          << GetProblemDescription(i);
    }
  }
}

TEST_F(NLPConstrainedProblemsTest, LineSearchInteriorPointNewtonMethod) {
  const double mu = 10.0;
  const int max_iter = 300;
  const double eta = 1e-1;
  const double tau = 0.95;
  std::unordered_set<std::string> expected_failures = {
      "Banana Function", "Banana Function With Inequality Constraints",
      "Function 3",      "Function 4",
      "Function 5",      "Function 10",
      "Function 11",     "Function 12",
      "Function 14",     "Function 15"};
  vect_n<double> x;
  for (int i = 0; i < funcs_f.size(); ++i) {
    if (expected_failures.count(funcs_name[i]) == 1) {
      std::cout << "Skipping known failure: '" << funcs_name[i] << "'"
                << std::endl;
      continue;
    }
    x = funcs_start[i];
    eval_count = 0;
    grad_count = 0;
    SCOPE_TRACE_EXCEPTION(
        make_nlip_newton_ls(funcs_f[i], funcs_grad[i], funcs_hess[i], mu,
                            max_iter, desired_tolerance, eta, tau)
            .set_eq_constraints(funcs_g[i], funcs_g_jac[i])
            .set_ineq_constraints(funcs_h[i], funcs_h_jac[i])(x));
    if (IsVectorFinite(x)) {
      if (std::isfinite(funcs_sol[i][0])) {
        EXPECT_THAT(x, VectorIsNear(funcs_sol[i], expected_tolerance))
            << GetProblemDescription(i);
      }
      CheckGradientProjection(x, i, expected_tolerance);
      EXPECT_THAT(funcs_g[i](x), VectorIsZero(expected_tolerance))
          << GetProblemDescription(i);
      EXPECT_THAT(funcs_h[i](x),
                  VectorIsGreater(vect_n<double>(funcs_h[i](x).size(), 0.0),
                                  expected_tolerance))
          << GetProblemDescription(i);
    }
  }
}

TEST_F(NLPConstrainedProblemsTest, LineSearchInteriorPointQuasiNewtonMethod) {
  const double mu = 10.0;
  const int max_iter = 100;
  const double eta = 1e-1;
  const double tau = 0.95;
  std::unordered_set<std::string> expected_failures = {
      "Banana Function", "Banana Function With Inequality Constraints",
      "Function 3",      "Function 4",
      "Function 5",      "Function 10",
      "Function 11",     "Function 12",
      "Function 14",     "Function 15",
      "Robot 1",         "Robot 2"};
  vect_n<double> x;
  for (int i = 0; i < funcs_f.size(); ++i) {
    if (expected_failures.count(funcs_name[i]) == 1) {
      std::cout << "Skipping known failure: '" << funcs_name[i] << "'"
                << std::endl;
      continue;
    }
    x = funcs_start[i];
    eval_count = 0;
    grad_count = 0;
    SCOPE_TRACE_EXCEPTION(
        make_nlip_quasi_newton_ls(funcs_f[i], funcs_grad[i], mu, max_iter,
                                  desired_tolerance, eta, tau)
            .set_eq_constraints(funcs_g[i], funcs_g_jac[i])
            .set_ineq_constraints(funcs_h[i], funcs_h_jac[i])(x));
    if (IsVectorFinite(x)) {
      if (std::isfinite(funcs_sol[i][0])) {
        EXPECT_THAT(x, VectorIsNear(funcs_sol[i], expected_tolerance))
            << GetProblemDescription(i);
      }
      CheckGradientProjection(x, i, expected_tolerance);
      EXPECT_THAT(funcs_g[i](x), VectorIsZero(expected_tolerance))
          << GetProblemDescription(i);
      EXPECT_THAT(funcs_h[i](x),
                  VectorIsGreater(vect_n<double>(funcs_h[i](x).size(), 0.0),
                                  expected_tolerance))
          << GetProblemDescription(i);
    }
  }
}

TEST_F(NLPConstrainedProblemsTest, ByrdOmojokunSQPNewtonMethod) {
  const double max_radius = 2.0;
  const int max_iter = 300;
  const double eta = 1e-3;
  const double rho = 0.8;
  std::unordered_set<std::string> expected_failures = {};
  vect_n<double> x;
  for (int i = 0; i < funcs_f.size(); ++i) {
    if (expected_failures.count(funcs_name[i]) == 1) {
      std::cout << "Skipping known failure: '" << funcs_name[i] << "'"
                << std::endl;
      continue;
    }
    x = funcs_start[i];
    if (!funcs_h[i](x).empty()) {
      continue;  // Cannot handle inequality constraints.
    }
    eval_count = 0;
    grad_count = 0;
    SCOPE_TRACE_EXCEPTION(
        make_bosqp_newton_tr(funcs_f[i], funcs_grad[i], funcs_hess[i],
                             max_radius, max_iter, desired_tolerance, eta, rho)
            .set_limiter([l = funcs_lower[i], u = funcs_upper[i]](
                             const vect_n<double>& x, vect_n<double>& dx) {
              box_limit_function(x, dx, l, u);
            })
            .set_eq_constraints(funcs_g[i], funcs_g_jac[i])(x));
    if (IsVectorFinite(x)) {
      if (std::isfinite(funcs_sol[i][0])) {
        EXPECT_THAT(x, VectorIsNear(funcs_sol[i], expected_tolerance))
            << GetProblemDescription(i);
      }
      CheckGradientProjection(x, i, expected_tolerance);
      EXPECT_THAT(funcs_g[i](x), VectorIsZero(expected_tolerance))
          << GetProblemDescription(i);
      EXPECT_THAT(funcs_h[i](x),
                  VectorIsGreater(vect_n<double>(funcs_h[i](x).size(), 0.0),
                                  expected_tolerance))
          << GetProblemDescription(i);
    }
  }
}

TEST_F(NLPConstrainedProblemsTest, ByrdOmojokunSQPQuasiNewtonMethod) {
  const double max_radius = 2.0;
  const int max_iter = 300;
  const double eta = 1e-3;
  const double rho = 0.8;
  std::unordered_set<std::string> expected_failures = {};
  vect_n<double> x;
  for (int i = 0; i < funcs_f.size(); ++i) {
    if (expected_failures.count(funcs_name[i]) == 1) {
      std::cout << "Skipping known failure: '" << funcs_name[i] << "'"
                << std::endl;
      continue;
    }
    x = funcs_start[i];
    if (!funcs_h[i](x).empty()) {
      continue;  // Cannot handle inequality constraints.
    }
    eval_count = 0;
    grad_count = 0;
    SCOPE_TRACE_EXCEPTION(
        make_bosqp_quasi_newton_tr(funcs_f[i], funcs_grad[i], max_radius,
                                   max_iter, desired_tolerance, eta, rho)
            .set_limiter([l = funcs_lower[i], u = funcs_upper[i]](
                             const vect_n<double>& x, vect_n<double>& dx) {
              box_limit_function(x, dx, l, u);
            })
            .set_eq_constraints(funcs_g[i], funcs_g_jac[i])(x));
    if (IsVectorFinite(x)) {
      if (std::isfinite(funcs_sol[i][0])) {
        EXPECT_THAT(x, VectorIsNear(funcs_sol[i], expected_tolerance))
            << GetProblemDescription(i);
      }
      CheckGradientProjection(x, i, expected_tolerance);
      EXPECT_THAT(funcs_g[i](x), VectorIsZero(expected_tolerance))
          << GetProblemDescription(i);
      EXPECT_THAT(funcs_h[i](x),
                  VectorIsGreater(vect_n<double>(funcs_h[i](x).size(), 0.0),
                                  expected_tolerance))
          << GetProblemDescription(i);
    }
  }
}

}  // namespace
}  // namespace ReaK::optim

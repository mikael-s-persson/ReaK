/**
 * \file conjugate_gradient_methods.h
 *
 * The following library is a collection of conjugate gradient methods.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
 */

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

#ifndef REAK_MATH_OPTIMIZATION_CONJUGATE_GRADIENT_METHODS_H_
#define REAK_MATH_OPTIMIZATION_CONJUGATE_GRADIENT_METHODS_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_num_exceptions.h"
#include "ReaK/math/lin_alg/vect_concepts.h"

#include <type_traits>

namespace ReaK::optim {

namespace detail {

template <typename T, typename Vector>
T fletcher_reeves_beta_impl(const Vector& dx, const Vector& dx_prev,
                            const Vector& /*unused*/) {
  return (dx * dx) / (dx_prev * dx_prev);
}

template <typename T, typename Vector>
T polak_ribiere_beta_impl(const Vector& dx, const Vector& dx_prev,
                          const Vector& /*unused*/) {
  return (dx * (dx - dx_prev)) / (dx_prev * dx_prev);
}

template <typename T, typename Vector>
T hestenes_stiefel_beta_impl(const Vector& dx, const Vector& dx_prev,
                             const Vector& /*unused*/) {
  Vector tmp = dx - dx_prev;
  return (dx * tmp) / (dx_prev * tmp);
}

template <typename T, typename Vector>
T dai_yuan_beta_impl(const Vector& dx, const Vector& dx_prev, const Vector& p) {
  return (dx * dx) / (p * (dx - dx_prev));
}

template <typename T, typename Vector>
T hager_zhang_beta_impl(const Vector& dx, const Vector& dx_prev,
                        const Vector& p) {
  Vector tmp = dx - dx_prev;
  T denom = tmp * p;
  return ((tmp - (T(2.0) * (tmp * tmp) / denom) * p) * dx) * (T(1.0) / denom);
}
}  // namespace detail

/**
 * This functor class can be used to compute the Fletcher-Reeves beta-value for the
 * conjugate gradient methods.
 * TEST PASSED
 */
struct fletcher_reeves_beta {

  /**
   * This computes the FR-beta value for the given delta-x and previous delta-x.
   * \param dx The current delta-x.
   * \param dx_prev The previous delta-x.
   * \param p The search direction.
   * \return The FR-beta value.
   */
  template <typename T>
  auto operator()(const T& dx, const T& dx_prev, const T& p) const {
    if constexpr (ReadableVector<T>) {
      return detail::fletcher_reeves_beta_impl<vect_value_type_t<T>>(
          dx, dx_prev, p);
    } else {
      return detail::fletcher_reeves_beta_impl<T>(dx, dx_prev, p);
    }
  }
};

/**
 * This functor class can be used to compute the Polak-Ribiere beta-value for the
 * conjugate gradient methods.
 * TEST PASSED
 */
struct polak_ribiere_beta {
  /**
   * This computes the PR-beta value for the given delta-x and previous delta-x.
   * \param dx The current delta-x.
   * \param dx_prev The previous delta-x.
   * \return The PR-beta value.
   */
  template <typename T>
  auto operator()(const T& dx, const T& dx_prev, const T& p) const {
    if constexpr (ReadableVector<T>) {
      return detail::polak_ribiere_beta_impl<vect_value_type_t<T>>(dx, dx_prev,
                                                                   p);
    } else {
      return detail::polak_ribiere_beta_impl<T>(dx, dx_prev, p);
    }
  }
};

/**
 * This functor class can be used to compute the Hestenes-Stiefel beta-value for the
 * conjugate gradient methods.
 * TEST PASSED
 */
struct hestenes_stiefel_beta {
  /**
   * This computes the HS-beta value for the given delta-x and previous delta-x.
   * \param dx The current delta-x.
   * \param dx_prev The previous delta-x.
   * \return The HS-beta value.
   */
  template <typename T>
  auto operator()(const T& dx, const T& dx_prev, const T& p) const {
    if constexpr (ReadableVector<T>) {
      return detail::hestenes_stiefel_beta_impl<vect_value_type_t<T>>(
          dx, dx_prev, p);
    } else {
      return detail::hestenes_stiefel_beta_impl<T>(dx, dx_prev, p);
    }
  }
};

/**
 * This functor class can be used to compute the Dai-Yuan beta-value for the
 * conjugate gradient methods.
 * TEST PASSED
 */
struct dai_yuan_beta {
  /**
   * This computes the DY-beta value for the given delta-x and previous delta-x.
   * \param dx The current delta-x.
   * \param dx_prev The previous delta-x.
   * \param p The search direction.
   * \return The DY-beta value.
   */
  template <typename T>
  auto operator()(const T& dx, const T& dx_prev, const T& p) const {
    if constexpr (ReadableVector<T>) {
      return detail::dai_yuan_beta_impl<vect_value_type_t<T>>(dx, dx_prev, p);
    } else {
      return detail::dai_yuan_beta_impl<T>(dx, dx_prev, p);
    }
  }
};

/**
 * This functor class can be used to compute the Hager-Zhang beta-value for the
 * conjugate gradient methods.
 * TEST PASSED
 */
struct hager_zhang_beta {
  /**
   * This overload computes the HZ-beta value for the given delta-x and previous delta-x
   * scalars.
   * \tparam T A scalar type.
   * \param dx The current delta-x value.
   * \param dx_prev The previous delta-x value.
   * \param p The search direction vector.
   * \return The HZ-beta value.
   */
  template <typename T>
  auto operator()(const T& dx, const T& dx_prev, const T& p) const {
    if constexpr (ReadableVector<T>) {
      return detail::hager_zhang_beta_impl<vect_value_type_t<T>>(dx, dx_prev,
                                                                 p);
    } else {
      return detail::hager_zhang_beta_impl<T>(dx, dx_prev, p);
    }
  }
};

/**
 * This function performs an iterative linear conjugate gradient method. This algorithm should find
 * the exact solution to the "Ax = b" problem (or minimizing "0.5 * x * A * x - x * b") in a finite
 * number of iterations. This algorithm only works under the assumption that the minimized function
 * is an exact quadratic, use the non-linear version if it is otherwise. Also, if the positive definite
 * matrix A is ill-conditioned, then it is preferrable to use the pre-conditionned version (see overload).
 * \test Must create a unit-test for this.
 * \param b The left-hand-side of the "Ax = b" equation.
 * \param A A positive definite symmetric matrix (should be well conditioned).
 * \param x The initial guess of the solution, and stores as output the final solution.
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance to be acheived by the residual error in the equation.
 */
template <WritableVector Vector, ReadableMatrix Matrix>
void linear_conj_grad_method(
    const Vector& b, const Matrix& A, Vector& x, unsigned int max_iter = 100,
    vect_value_type_t<Vector> abs_tol = vect_value_type_t<Vector>(1e-6)) {
  using ValueType = vect_value_type_t<Vector>;
  using std::sqrt;

  Vector r = b - A * x;
  ValueType rr = r * r;
  Vector p = r;
  Vector Ap = A * p;
  unsigned int k = 0;
  while (true) {
    ValueType alpha = rr / (p * Ap);
    ValueType beta = 1.0 / rr;
    x += alpha * p;
    if (++k > max_iter) {
      throw maximum_iteration(max_iter);
    }
    r -= alpha * Ap;
    rr = r * r;
    if (sqrt(rr) < abs_tol) {
      break;
    }
    beta *= rr;
    p *= beta;
    p += r;
    Ap = A * p;
  }
}

/**
 * This function performs an iterative linear conjugate gradient method with preconditionning. This
 * algorithm should find the exact solution to the "Ax = b" problem (or minimizing "0.5 * x * A * x - x * b")
 * in a finite number of iterations. This algorithm only works under the assumption that the minimized function
 * is an exact quadratic, use the non-linear version if it is otherwise. This overload version
 * is most appropriate if the positive definite matrix A is ill-conditioned. This algorithm uses a
 * matrix M_inv which pre-conditions the residual vector to improve the convergence, M_inv should be
 * chosen such that the condition number of the product "M_inv * A" is low (approaching 1.0).
 * \test Must create a unit-test for this.
 * \param b The left-hand-side of the "Ax = b" equation.
 * \param A A positive definite symmetric matrix (should be well conditioned).
 * \param M_inv The preconditionning positive definite symmetric matrix, should be chosen such that the condition number
 * of the product "M_inv * A" is low (approaching 1.0).
 * \param x The initial guess of the solution, and stores as output the final solution.
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance to be acheived by the residual error in the equation.
 */
template <WritableVector Vector, ReadableMatrix Matrix>
void linear_conj_grad_method(
    const Vector& b, const Matrix& A, const Matrix& M_inv, Vector& x,
    unsigned int max_iter = 100,
    vect_value_type_t<Vector> abs_tol = vect_value_type_t<Vector>(1e-6)) {
  using ValueType = vect_value_type_t<Vector>;
  using std::sqrt;

  Vector r = b - A * x;
  Vector z = M_inv * r;
  ValueType zr = z * r;
  Vector p = z;
  Vector Ap = A * p;
  unsigned int k = 0;
  while (true) {
    ValueType alpha = zr / (p * Ap);
    ValueType beta = 1.0 / zr;
    x += alpha * p;
    r -= alpha * Ap;
    z = M_inv * r;
    zr = z * r;
    if (++k > max_iter) {
      throw maximum_iteration(max_iter);
    }
    if (sqrt(zr) < abs_tol) {
      break;
    }
    beta *= zr;
    p *= beta;
    p += z;
    Ap = A * p;
  }
}

/**
 * This function performs a non-linear conjugate gradient method to find the minimum of a cost function.
 * This algorithm requires a function and its gradient, an initial guess, and methods to compute
 * the "beta" and "alpha" values. The beta calculator can be one of the following (or any similar functor):
 * fletcher_reeves_beta, polak_ribiere_beta, or hestenes_stiefel_beta. The functor needed to compute
 * the alpha value is a line-search functor, similar to line_search_backtracking (preferred),
 * bounded_line_srch_fibonacci, bounded_line_srch_gold_sect, or bounded_line_srch_dichotomous.
 * TEST PASSED
 * \tparam Function A functor type to represent the cost-function.
 * \tparam GradFunction A functor type to represent the gradient of the cost-function.
 * \tparam Vector The vector type of the independent variable on which the search is performed.
 * \tparam BetaCalculator The functor type for computing the beta value (see fletcher_reeves_beta, polak_ribiere_beta,
 * or hestenes_stiefel_beta).
 * \tparam LineSearcher The functor type to perform the line-search (see line_search_backtracking (preferred),
 * bounded_line_srch_fibonacci, bounded_line_srch_gold_sect, or bounded_line_srch_dichotomous).
 * \param f The cost-function to minimize.
 * \param df The gradient of the cost-function.
 * \param x The initial guess of the solution, and stores, as output, the obtained minimum.
 * \param max_iter The maximum number of iterations to perform.
 * \param get_beta The functor to use to compute the beta value.
 * \param get_alpha The functor to use to perform the line-search for the "alpha" value (compute the best update along
 * the search direction).
 * \param abs_tol The tolerance at which to stop the algorithm.
 */
template <typename Function, typename GradFunction, typename Vector,
          typename BetaCalculator, typename LineSearcher>
void non_linear_conj_grad_method(
    Function f, GradFunction df, Vector& x, unsigned int max_iter,
    BetaCalculator get_beta, LineSearcher get_alpha,
    vect_value_type_t<Vector> abs_tol = vect_value_type_t<Vector>(1e-6)) {
  using ValueType = vect_value_type_t<Vector>;
  using std::sqrt;

  Vector dx = -df(x);
  Vector Dx = dx;
  Vector dx_prev = dx;
  auto alpha = ValueType(1.0);
  auto beta = ValueType(1.0);
  auto adfp_prev = ValueType(-1.0);
  unsigned int k = 0;

  while (norm_2(Dx) > abs_tol) {
    ValueType alpha_0 = adfp_prev;
    adfp_prev = (Dx * dx);
    alpha_0 /= adfp_prev;
    alpha = get_alpha(f, df, ValueType(0.0), alpha_0, x, Dx, abs_tol);
    adfp_prev *= alpha;
    x += alpha * Dx;
    if (++k > max_iter) {
      throw maximum_iteration(max_iter);
    }
    dx_prev = dx;
    dx = -df(x);
    beta = get_beta(dx, dx_prev, Dx);
    if (beta < ValueType(0.0)) {
      beta = ValueType(0.0);
    }
    Dx *= beta;
    Dx += dx;
  }
}

}  // namespace ReaK::optim

#endif  // REAK_MATH_OPTIMIZATION_CONJUGATE_GRADIENT_METHODS_H_

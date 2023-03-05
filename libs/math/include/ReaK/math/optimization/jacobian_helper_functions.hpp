/**
 * \file jacobian_helper_functions.hpp
 *
 * The following library provides implementations of a number of helper functions related to
 * Jacobian matrices.
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

#ifndef REAK_JACOBIAN_HELPER_FUNCTIONS_HPP
#define REAK_JACOBIAN_HELPER_FUNCTIONS_HPP

#include "ReaK/core/base/defs.hpp"
#include "ReaK/math/lin_alg/mat_alg.hpp"
#include "ReaK/math/lin_alg/mat_svd_method.hpp"

#include "ReaK/math/optimization/optim_exceptions.hpp"

namespace ReaK::optim {

/**
 * Check the weight-jacobian of basis nonlinear functions in input_count variables
 * evaluated at current parameters, for consistency with the function itself.
 *
 * Based on fortran77 subroutine CHKDER by
 * Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
 * Argonne National Laboratory. MINPACK project. March 1980.
 *
 * The function does not perform reliably if cancellation or
 * rounding errors cause a severe loss of significance in the
 * evaluation of a function. therefore, none of the parameters
 * should be unusually small (in particular, zero) or any
 * other value which may cause loss of significance.
 * \test Must create a unit-test for this.
 *
 * \tparam Function The functor type that can represent the evaluation function.
 * \tparam Vector1 The input vector type.
 * \tparam Vector2 The output vector type.
 * \tparam Matrix The Jacobian matrix type.
 * \param f The functor object to evaluate the function at a given point.
 * \param x The point at which the jacobian was evaluated.
 * \param y The output vector of the function evaluated as the given point.
 * \param Jac The Jacobian matrix evaluated at the given point.
 * \param tol The numerical tolerance value below which a tolerance value can be considered zero.
 * \return An error vector which contains a measure of correctness of the gradients in the different
 *         output dimensions. A value of 1.0 at a given component indicates no significant loss of
 *         precision on the gradient on that coordinate. A value of 0.0 indicates that the gradient
 *         (or Jacobian) at this coordinate is incorrect. Intermediate values give a log-scale of the
 *         goodness of the gradient at the coordinate in question.
 */
template <typename Function, typename Vector1, typename Vector2,
          typename Matrix>
Vector2 check_jacobian_consistency(
    Function f, const Vector1& x, const Vector2& y, const Matrix& Jac,
    vect_value_type_t<Vector1> tol = vect_value_type_t<Vector1>(1e-6)) {
  using ValueType = vect_value_type_t<Vector1>;
  using std::abs;
  using std::log10;

  ValueType factor = 100.0;
  ValueType tol_sqr = tol * tol;

  int N = x.size();
  int M = y.size();

  /* compute pp */
  Vector1 x1 = x;
  for (int j = 0; j < N; ++j) {
    ValueType temp = tol * abs(x[j]);
    if (temp < tol) {
      temp = tol;
    }
    x1[j] += temp;
  }

  /* compute fvecp=func(pp) */
  Vector2 y1 = f(x1);

  ValueType epsf = factor * tol_sqr;
  ValueType epslog = log10(tol);

  Vector2 err = y;
  for (int i = 0; i < M; ++i) {
    err[i] = ValueType(0.0);
  }

  for (int j = 0; j < N; ++j) {
    ValueType temp = abs(x[j]);
    if (temp < tol) {
      temp = ValueType(1.0);
    }

    for (int i = 0; i < M; ++i) {
      err[i] += temp * Jac(i, j);
    }
  }

  for (int i = 0; i < M; ++i) {
    ValueType temp(1.0);
    if ((abs(y[i]) < tol_sqr) && (abs(y1[i]) < tol_sqr) &&
        (abs(y1[i] - y[i]) >= epsf * abs(y[i]))) {
      temp =
          tol * abs((y1[i] - y[i]) / tol - err[i]) / (abs(y[i]) + abs(y1[i]));
    }
    err[i] = ValueType(1.0);
    if ((temp > tol * tol) && (temp < tol)) {
      err[i] = (log10(temp) - epslog) / epslog;
    }
    if (temp >= tol) {
      err[i] = ValueType(0.0);
    }
  }

  return err;
}

/**
 * This function computes in C the covariance matrix corresponding to a least
 * squares fit. JtJ is the approximate Hessian at the solution (i.e. J^T*J, where
 * J is the jacobian at the solution), sumsq is the sum of squared residuals
 * (i.e. goodness of fit) at the solution, M is the number of parameters (variables)
 * and N the number of observations. JtJ can coincide with C.
 *
 * if JtJ is of full rank, C is computed as sumsq/(n-m)*(JtJ)^-1
 * otherwise C=sumsq/(n-r)*(JtJ)^+ where r is JtJ's rank and ^+ denotes
 * the pseudoinverse. The diagonal of C is made up from the estimates of
 * the variances of the estimated regression coefficients.
 * See the documentation of routine E04YCF from the NAG fortran lib
 *
 * The function returns the rank of JtJ if successful, 0 on error
 *
 * JtJ and C are MxM
 * \test Must create a unit-test for this.
 *
 * \tparam Matrix1 The matrix type which represents the approximate Hessian.
 * \tparam Matrix2 The matrix type which represents the resulting covariance matrix.
 * \param JtJ The approximate Hessian matrix.
 * \param C Stores as output the resulting covariance matrix.
 * \param sumsq The sum of squared residuals at the solution.
 * \param N The number of observations from which the approximate Hessian comes from.
 * \return The numerical rank of JtJ (a zero value is an error).
 *
 */
template <typename Matrix1, typename Matrix2>
mat_size_type_t<Matrix2> compute_covariance(const Matrix1& JtJ, Matrix2& C,
                                            mat_value_type_t<Matrix2> sumsq,
                                            mat_size_type_t<Matrix2> N) {
  using ValueType = mat_value_type_t<Matrix2>;

  mat<ValueType, mat_structure::square> U, V;
  mat<ValueType, mat_structure::diagonal> E;
  decompose_SVD(JtJ, U, E, V);
  int rnk = numrank_SVD(E);
  pseudoinvert_SVD(U, E, V, C);

  if (rnk) {
    C *= sumsq / ValueType(N - rnk);
  }

  return rnk;
}

}  // namespace ReaK::optim

#endif

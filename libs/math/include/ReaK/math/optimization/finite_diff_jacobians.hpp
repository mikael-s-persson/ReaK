/**
 * \file finite_diff_jacobians.hpp
 *
 * The following library provides implementations of finite difference methods for evaluating a Jacobian
 * for a function of multiple variables and multiple outputs (or single).
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

#ifndef REAK_FINITE_DIFF_JACOBIANS_HPP
#define REAK_FINITE_DIFF_JACOBIANS_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/lin_alg/mat_alg.hpp>

#include "optim_exceptions.hpp"

#include <type_traits>

namespace ReaK::optim {

namespace detail {

template <typename Function, typename Vector1, typename Vector2,
          typename Matrix>
void compute_jacobian_2pts_forward_impl(Function f, Vector1& x,
                                        const Vector2& y, Matrix& jac,
                                        vect_value_type_t<Vector1> delta) {
  using ValueType = vect_value_type_t<Vector1>;
  using std::abs;

  int N = x.size();
  int M = y.size();
  if (jac.get_row_count() != M) {
    jac.set_row_count(M);
  }
  if (jac.get_col_count() != N) {
    jac.set_col_count(N);
  }

  Vector2 y1 = y;

  for (int j = 0; j < N; ++j) {
    /* determine d=max(1E-04*|p[j]|, delta), see HZ */
    ValueType d = ValueType(1E-04) * abs(x[j]);
    if (d < delta) {
      d = delta;
    }

    ValueType tmp = x[j];
    x[j] += d;

    y1 = f(x);

    x[j] = tmp; /* restore */

    slice(jac)(range(0, M), j) = (y1 - y) * (1.0 / d);
  }
}

template <typename Function, typename Vector1, typename Vector2,
          typename Matrix>
void compute_jacobian_2pts_central_impl(Function f, Vector1& x,
                                        const Vector2& y, Matrix& jac,
                                        vect_value_type_t<Vector1> delta) {
  using ValueType = vect_value_type_t<Vector1>;
  using std::abs;

  int N = x.size();
  int M = y.size();
  if (jac.get_row_count() != M) {
    jac.set_row_count(M);
  }
  if (jac.get_col_count() != N) {
    jac.set_col_count(N);
  }

  Vector2 y_prev = y;
  Vector2 y_next = y;

  for (int j = 0; j < N; ++j) {
    /* determine d=max(1E-04*|p[j]|, delta), see HZ */
    ValueType d = ValueType(1E-04) * abs(x[j]);
    if (d < delta) {
      d = delta;
    }

    ValueType tmp = x[j];
    x[j] -= d;
    y_prev = f(x);

    x[j] = tmp + d;
    y_next = f(x);
    x[j] = tmp; /* restore */

    slice(jac)(range(0, M), j) = (y_next - y_prev) * (0.5 / d);
  }
}

template <typename Function, typename Vector1, typename Vector2,
          typename Matrix>
void compute_jacobian_5pts_central_impl(Function f, Vector1& x,
                                        const Vector2& y, Matrix& jac,
                                        vect_value_type_t<Vector1> delta) {
  using ValueType = vect_value_type_t<Vector1>;
  using std::abs;

  int N = x.size();
  int M = y.size();
  if (jac.get_row_count() != M) {
    jac.set_row_count(M);
  }
  if (jac.get_col_count() != N) {
    jac.set_col_count(N);
  }

  Vector2 y0 = y;
  Vector2 y1 = y;
  Vector2 y2 = y;
  Vector2 y3 = y;

  for (int i = 0; i < N; ++i) {
    /* determine d=max(1E-04*|p[j]|, delta), see HZ */
    ValueType d = ValueType(1E-04) * abs(x[i]);
    if (d < delta) {
      d = delta;
    }

    ValueType tmp = x[i];
    x[i] -= 2.0 * d;
    y0 = f(x);
    x[i] += d;
    y1 = f(x);
    x[i] += 2.0 * d;
    y2 = f(x);
    x[i] += d;
    y3 = f(x);

    x[i] = tmp;  // restore

    slice(jac)(range(0, M), i) =
        (y0 - 8.0 * (y1 - y2) - y3) * (1.0 / (12.0 * d));
  }
}
}  // namespace detail

template <typename Function, typename Vector1, typename Vector2,
          typename Matrix, typename Scalar = mat_value_type_t<Matrix>>
void compute_jacobian_2pts_forward(Function f, Vector1& x, const Vector2& y,
                                   Matrix& jac, Scalar delta = Scalar(1e-6)) {
  static_assert(is_writable_matrix_v<Matrix>);
  if constexpr (!is_readable_vector_v<Vector1>) {
    vect<Scalar, 1> x_tmp;
    x_tmp[0] = x;
    compute_jacobian_2pts_forward(
        [&f](const vect<Scalar, 1>& x) { return f(x[0]); }, x_tmp, y, jac,
        delta);
    x = x_tmp[0];
  } else if constexpr (!is_readable_vector_v<Vector2>) {
    vect<Scalar, 1> y_tmp;
    y_tmp[0] = y;
    compute_jacobian_2pts_forward(
        [&f](const Vector1& x) -> vect<Scalar, 1> {
          vect<Scalar, 1> result;
          result[0] = f(x);
          return result;
        },
        x, y_tmp, jac, delta);
  } else if constexpr (is_fully_writable_matrix_v<Matrix>) {
    detail::compute_jacobian_2pts_forward_impl(f, x, y, jac, delta);
  } else {
    mat<vect_value_type_t<Vector1>, mat_structure::rectangular> jac_tmp(
        y.size(), x.size());
    detail::compute_jacobian_2pts_forward_impl(f, x, y, jac_tmp, delta);
    jac = jac_tmp;
  }
}

template <typename Function, typename Vector1, typename Vector2,
          typename Matrix, typename Scalar = mat_value_type_t<Matrix>>
void compute_jacobian_2pts_central(Function f, Vector1& x, const Vector2& y,
                                   Matrix& jac, Scalar delta = Scalar(1e-6)) {
  static_assert(is_writable_matrix_v<Matrix>);
  if constexpr (!is_readable_vector_v<Vector1>) {
    vect<Scalar, 1> x_tmp;
    x_tmp[0] = x;
    compute_jacobian_2pts_central(
        [&f](const vect<Scalar, 1>& x) { return f(x[0]); }, x_tmp, y, jac,
        delta);
    x = x_tmp[0];
  } else if constexpr (!is_readable_vector_v<Vector2>) {
    vect<Scalar, 1> y_tmp;
    y_tmp[0] = y;
    compute_jacobian_2pts_central(
        [&f](const Vector1& x) -> vect<Scalar, 1> {
          vect<Scalar, 1> result;
          result[0] = f(x);
          return result;
        },
        x, y_tmp, jac, delta);
  } else if constexpr (is_fully_writable_matrix_v<Matrix>) {
    detail::compute_jacobian_2pts_central_impl(f, x, y, jac, delta);
  } else {
    mat<vect_value_type_t<Vector1>, mat_structure::rectangular> jac_tmp(
        y.size(), x.size());
    detail::compute_jacobian_2pts_central_impl(f, x, y, jac_tmp, delta);
    jac = jac_tmp;
  }
}

template <typename Function, typename Vector1, typename Vector2,
          typename Matrix, typename Scalar = mat_value_type_t<Matrix>>
void compute_jacobian_5pts_central(Function f, Vector1& x, const Vector2& y,
                                   Matrix& jac, Scalar delta = Scalar(1e-6)) {
  static_assert(is_writable_matrix_v<Matrix>);
  if constexpr (!is_readable_vector_v<Vector1>) {
    vect<Scalar, 1> x_tmp;
    x_tmp[0] = x;
    compute_jacobian_5pts_central(
        [&f](const vect<Scalar, 1>& x) { return f(x[0]); }, x_tmp, y, jac,
        delta);
    x = x_tmp[0];
  } else if constexpr (!is_readable_vector_v<Vector2>) {
    vect<Scalar, 1> y_tmp;
    y_tmp[0] = y;
    compute_jacobian_5pts_central(
        [&f](const Vector1& x) -> vect<Scalar, 1> {
          vect<Scalar, 1> result;
          result[0] = f(x);
          return result;
        },
        x, y_tmp, jac, delta);
  } else if constexpr (is_fully_writable_matrix_v<Matrix>) {
    detail::compute_jacobian_5pts_central_impl(f, x, y, jac, delta);
  } else {
    mat<vect_value_type_t<Vector1>, mat_structure::rectangular> jac_tmp(
        y.size(), x.size());
    detail::compute_jacobian_5pts_central_impl(f, x, y, jac_tmp, delta);
    jac = jac_tmp;
  }
}

}  // namespace ReaK::optim

#endif

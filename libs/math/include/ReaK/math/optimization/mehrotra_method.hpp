/**
 * \file mehrotra_method.hpp
 *
 * The following library is an implementation of the Mehrotra Method to solve a linear programming
 * problem. The algorithm follows the specification given by Nocedal's Numerical Optimization book.
 * The Mehrotra method is an interior-point predictor-corrector algorithm for solving standard
 * linear programming problems (linear optimization).
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

#ifndef REAK_MEHROTRA_METHOD_HPP
#define REAK_MEHROTRA_METHOD_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/math/lin_alg/mat_qr_decomp.hpp>

#include "optim_exceptions.hpp"
#include "quadratic_programs.hpp"

#include <vector>

namespace ReaK::optim {

/**
 * This function is an implementation of the Mehrotra Method to solve a linear programming
 * problem. The Mehrotra method is an interior-point predictor-corrector algorithm for solving standard
 * linear programming problems (linear optimization). It solves the following problem: \n
 * \n
 *           min c'x \n
 *               Ax = b \n
 *               x >= 0 \n
 * \n
 * The implementation was inspired from the algorithm described in the book:\n
 *   Nocedal, Numerical Optimization, 2nd Ed..
 * \test Must create a unit-test for this. So far, this method fails the tests.
 *
 * \tparam Matrix A general matrix type, should model the WritableMatrixConcept (and be fully-writable).
 * \tparam Vector1 A vector type, should model the WritableVectorConcept.
 * \tparam Vector2 A vector type, should model the WritableVectorConcept.
 * \param A The constraint matrix of dimension M*N.
 * \param b The b vector of dimension M.
 * \param c The cost vector of dimension N.
 * \param x Stores, as output, the optimal vector.
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance to which the solution is accepted (tolerance on the steps).
 *
 * \author Mikael Persson
 */
template <typename Matrix, typename Vector1, typename Vector2>
void mehrotra_method(
    const Matrix& A, const Vector1& b, const Vector2& c, Vector2& x,
    unsigned int max_iter = 100,
    vect_value_type_t<Vector1> abs_tol =
        std::numeric_limits<vect_value_type_t<Vector1>>::epsilon()) {
  using ValueType = vect_value_type_t<Vector1>;
  using std::abs;
  using std::swap;

  int N = c.size();
  int M = b.size();

  mat<ValueType, mat_structure::rectangular> A_tmp(transpose_view(A));
  mat<ValueType, mat_structure::rectangular> R(N, M);
  mat<ValueType, mat_structure::square> Q(N);
  decompose_QR(A_tmp, Q, R, abs_tol);
  mat<ValueType, mat_structure::rectangular> L =
      mat<ValueType, mat_structure::rectangular>(transpose_view(R));
  L.set_col_count(M, true);

  // first phase, obtain a starting feasible interior-point solution:
  Vector1 b_tmp = b;
  mat_vect_adaptor<Vector1> b_tmp_mat(b_tmp);
  ReaK::detail::backsub_Cholesky_impl(L, b_tmp_mat);
  x = b_tmp * A;

  Vector1 lambda = A * c;
  mat_vect_adaptor<Vector1> lambda_mat(lambda);
  ReaK::detail::backsub_Cholesky_impl(L, lambda_mat);

  Vector2 s = x - (lambda * A);

  ValueType del_x(0.0);
  ValueType del_s(0.0);
  for (int i = 0; i < N; ++i) {
    if (x[i] < del_x) {
      del_x = x[i];
    }
    if (s[i] < del_s) {
      del_s = s[i];
    }
  }
  del_x *= ValueType(-2.0 / 3.0);
  del_s *= ValueType(-2.0 / 3.0);
  for (int i = 0; i < N; ++i) {
    x[i] += del_x;
    s[i] += del_s;
  }
  del_x = ValueType(0.0);
  del_s = ValueType(0.0);
  ValueType sum(0.0);
  for (int i = 0; i < N; ++i) {
    sum += s[i] * x[i];
    del_x += s[i];
    del_s += x[i];
  }
  del_x = ValueType(0.5) * sum / del_x;
  del_s = ValueType(0.5) * sum / del_s;
  for (int i = 0; i < N; ++i) {
    x[i] += del_x;
    s[i] += del_s;
  }

  Vector1 r_b = b - A * x;
  Vector2 r_c = c - s - lambda * A;
  Vector2 x_s = x;
  for (int i = 0; i < N; ++i) {
    x_s /= s[i];
  }

  Vector2 dx = x;
  Vector1 dl = lambda;
  Vector2 ds = s;

  ValueType eta(0.9);
  ValueType mu = (x * s) / ValueType(N);
  unsigned int k = 0;

  // second phase, the predictor-corrector loop:
  do {
    dl = r_b;
    for (int i = 0; i < dx.size(); ++i) {
      dx[i] = r_c[i] * x_s[i] + x[i];
    }
    dl += A * dx;

    for (int i = 0; i < N; ++i) {
      ValueType tmp = sqrt(x_s[i]);
      for (int j = 0; j < M; ++j) {
        A_tmp(i, j) = A(j, i) * tmp;
      }
    }
    decompose_QR(A_tmp, Q, R, abs_tol);
    L = transpose_view(R);
    L.set_col_count(M, true);
    mat_vect_adaptor<Vector1> dl_mat(dl);
    ReaK::detail::backsub_Cholesky_impl(L, dl_mat);
    ds = r_c - (dl * A);
    for (int i = 0; i < N; ++i) {
      dx[i] = -x[i] - x_s[i] * ds[i];
    }

    ValueType a_p(1.0);
    ValueType a_d(1.0);
    for (int i = 0; i < N; ++i) {
      if (dx[i] < ValueType(0.0)) {
        ValueType tmp = -x[i] / dx[i];
        if (tmp < a_p) {
          a_p = tmp;
        }
      }
      if (ds[i] < ValueType(0.0)) {
        ValueType tmp = -s[i] / ds[i];
        if (tmp < a_d) {
          a_d = tmp;
        }
      }
    }
    ValueType mu_aff = ((x + a_p * dx) * (s + a_d * ds)) / ValueType(N);
    ValueType sigma = (mu_aff * mu_aff * mu_aff) / (mu * mu * mu);

    dl = r_b;
    for (int i = 0; i < N; ++i) {
      dx[i] = -x[i] - (dx[i] * ds[i] - sigma * mu) / s[i];
      ds[i] = r_c[i] * x_s[i];
    }
    dl += A * (ds - dx);
    ReaK::detail::backsub_Cholesky_impl(L, dl_mat);
    ds = r_c - (dl * A);
    for (int i = 0; i < N; ++i) {
      dx[i] -= x_s[i] * ds[i];
    }

    a_p = ValueType(1.0);
    a_d = ValueType(1.0);
    for (int i = 0; i < N; ++i) {
      if (dx[i] < ValueType(0.0)) {
        ValueType tmp = -eta * x[i] / dx[i];
        if (tmp < a_p) {
          a_p = tmp;
        }
      }
      if (ds[i] < ValueType(0.0)) {
        ValueType tmp = -eta * s[i] / ds[i];
        if (tmp < a_d) {
          a_d = tmp;
        }
      }
    }

    dx *= a_p;
    dl *= a_d;
    ds *= a_d;
    x += dx;
    lambda += dl;
    s += ds;
    r_b *= ValueType(1.0) - a_p;
    r_c *= ValueType(1.0) - a_d;
    mu = (x * s) / ValueType(N);
    x_s = x;
    for (int i = 0; i < x.size(); ++i) {
      x_s /= s[i];
    }

    eta += ValueType(0.25) * (ValueType(1.0) - eta);

    if (++k > max_iter) {
      throw maximum_iteration(max_iter);
    }
  } while ((dx * dx) > abs_tol);
}

/**
 * This function is an implementation of the Mehrotra Method to solve a quadratic programming
 * problem. The Mehrotra method is an interior-point predictor-corrector algorithm for solving standard
 * quadratic programming problems. It solves the following problem: \n
 * \n
 *           min 0.5 * x' G x + c'x \n
 *               Ax >= b \n
 *               Ex = d \n
 * \n
 * The implementation was inspired from the algorithm described in the book:\n
 *   Nocedal, Numerical Optimization, 2nd Ed..
 * \test Must create a unit-test for this. TEST PASSED for equality constraints only.
 *
 * \tparam Matrix1 A general matrix type, should model the WritableMatrixConcept (and be fully-writable).
 * \tparam Vector1 A vector type, should model the WritableVectorConcept.
 * \tparam Matrix2 A general matrix type, should model the WritableMatrixConcept (and be fully-writable).
 * \tparam Vector2 A vector type, should model the WritableVectorConcept.
 * \param A The constraint matrix of dimension M*N.
 * \param b The b vector of dimension M.
 * \param G The constraint matrix of dimension N*N.
 * \param c The cost vector of dimension N.
 * \param E The equality-constraint matrix of dimension K*N.
 * \param d The equality-constraint vector of dimension K.
 * \param x Stores, as output, the optimal vector.
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance to which the solution is accepted (tolerance on the steps).
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Vector1, typename Matrix2,
          typename Vector2, typename Matrix3, typename Vector3>
void mehrotra_QP_method(
    const Matrix1& A, const Vector1& b, const Matrix2& G, const Vector2& c,
    const Matrix3& E, const Vector3& d, Vector2& x, unsigned int max_iter,
    vect_value_type_t<Vector1> abs_tol =
        std::numeric_limits<vect_value_type_t<Vector1>>::epsilon()) {
  using ValueType = vect_value_type_t<Vector1>;
  using std::abs;
  using std::swap;

  if (A.get_row_count() == 0) {
    null_space_QP_method(E, d, G, c, x, abs_tol);
    return;
  }

  int N = c.size();
  int M = b.size();
  int K = d.size();

  mat<ValueType, mat_structure::rectangular> E_tmp(transpose_view(E));
  mat<ValueType, mat_structure::rectangular> E_R(N, K);
  mat<ValueType, mat_structure::square> E_Q(N);
  decompose_QR(E_tmp, E_Q, E_R, abs_tol);
  mat<ValueType, mat_structure::rectangular> E_L =
      mat<ValueType, mat_structure::rectangular>(transpose_view(E_R));
  E_L.set_col_count(K, true);
  ValueType E_L_trace(0.0);
  for (int i = 0; i < K; ++i) {
    E_L_trace += abs(E_L(i, i));
  }

  mat_const_sub_block<mat<ValueType, mat_structure::square>> E_Y(E_Q, N, K, 0,
                                                                 0);
  mat_const_sub_block<mat<ValueType, mat_structure::square>> E_Z(E_Q, N, N - K,
                                                                 0, K);
  Vector2 dx = x;
  mat_vect_adaptor<Vector2> dx_y(dx, K, 1, 0);
  mat_vect_adaptor<Vector2> dx_z(dx, N - K, 1, K);

  // first phase, obtain a starting feasible interior-point solution:
  Vector1 r_b = b;
  r_b -= A * x;
  mat_vect_adaptor<Vector1> r_b_mat(r_b);

  Vector2 r_c = -c;
  r_c -= G * x;
  mat_vect_adaptor<Vector2> r_c_mat(r_c);

  Vector3 r_e = d;
  mat_vect_adaptor<Vector3> r_e_mat(r_e);
  if (K > 0) {
    r_e -= E * x;
    dx_y = r_e_mat;
    ReaK::detail::forwardsub_L_impl(E_L, dx_y, abs_tol);

    mat<ValueType, mat_structure::rectangular> e_y = E_Y * dx_y;
    r_b_mat -= A * e_y;
    r_c_mat -= G * e_y;
  }

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      r_c[i] += A(j, i) * (ValueType(1.0) + r_b[j]);
    }
  }
  for (int i = 0; i < N - K; ++i) {
    dx_z(i, 0) = ValueType(0.0);
    for (int j = 0; j < N; ++j) {
      dx_z(i, 0) += E_Z(j, i) * r_c[j];
    }
  }
  mat<ValueType, mat_structure::symmetric> ZGZ(transpose_view(E_Z) * G * E_Z);

  mat<ValueType, mat_structure::rectangular> AZ(A * E_Z);
  mat<ValueType, mat_structure::symmetric> ZG_AAZ(transpose_view(AZ) * AZ +
                                                  ZGZ);
  mat<ValueType, mat_structure::square> LHS_L(N - K);

  decompose_Cholesky(ZG_AAZ, LHS_L, abs_tol);
  ReaK::detail::backsub_Cholesky_impl(LHS_L, dx_z);

  Vector1 dl = r_b;
  mat_vect_adaptor<Vector1> dl_mat(dl);
  dl_mat -= AZ * dx_z;

  Vector1 dy = -dl;
  for (int i = 0; i < M; ++i) {
    dy[i] -= ValueType(1.0);
  }

  Vector1 y = dy;
  Vector1 l = dl;
  for (int i = 0; i < M; ++i) {
    y[i] = abs(ValueType(1.0) + dy[i]);
    if (y[i] < ValueType(1.0)) {
      y[i] = ValueType(1.0);
    }
    l[i] = abs(ValueType(1.0) + dl[i]);
    if (l[i] < ValueType(1.0)) {
      l[i] = ValueType(1.0);
    }
  }

  ValueType mu = (y * l) / ValueType(M);
  ValueType eta(0.5);
  unsigned int k = 0;

  do {
    r_c = -c;
    r_c -= G * x;
    r_c += l * A;

    r_b = b;
    r_b -= A * x;

    r_e = d;
    if (K > 0) {
      r_e -= E * x;
      dx_y = r_e_mat;
      ReaK::detail::forwardsub_L_impl(E_L, dx_y, abs_tol);

      mat<ValueType, mat_structure::rectangular> e_y = E_Y * dx_y;
      r_b_mat -= A * e_y;
      r_c_mat -= G * e_y;
    }
    for (int i = 0; i < M; ++i) {
      r_b[i] *= l[i] / y[i];
    }

    dx_z = transpose_view(E_Z) * r_c_mat + transpose_view(AZ) * r_b_mat;

    mat<ValueType, mat_structure::rectangular> AZ_temp = AZ;
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N - K; ++j) {
        AZ_temp(i, j) *= l[i] / y[i];
      }
    }
    ZG_AAZ = transpose_view(AZ) * AZ_temp + ZGZ;
    decompose_Cholesky(ZG_AAZ, LHS_L, abs_tol);
    ReaK::detail::backsub_Cholesky_impl(LHS_L, dx_z);
    dl_mat = r_b_mat - AZ_temp * dx_z;
    for (int i = 0; i < M; ++i) {
      dy[i] = -y[i] * (dl[i] / l[i] + ValueType(1.0));
    }

    ValueType alpha = ValueType(0.995);
    for (int i = 0; i < M; ++i) {
      if (y[i] < -alpha * dy[i]) {
        alpha = -y[i] / dy[i];
      }
      if (l[i] < -alpha * dl[i]) {
        alpha = -l[i] / dl[i];
      }
    }
    ValueType mu_aff = (y + alpha * dy) * (l + alpha * dl) / ValueType(M);
    ValueType sigma = (mu_aff * mu_aff * mu_aff) / (mu * mu * mu);

    for (int i = 0; i < M; ++i) {
      r_b[i] += (sigma * mu - dl[i] * dy[i]) / y[i];
    }

    dx_z = transpose_view(E_Z) * r_c_mat + transpose_view(AZ) * r_b_mat;
    ReaK::detail::backsub_Cholesky_impl(LHS_L, dx_z);
    dl_mat = r_b_mat - AZ_temp * dx_z;
    for (int i = 0; i < M; ++i) {
      dy[i] = -y[i] * (dl[i] / l[i] + ValueType(1.0));
    }

    ValueType a_p = ValueType(0.995);
    ValueType a_d = ValueType(0.995);
    for (int i = 0; i < M; ++i) {
      if (a_p * dy[i] < -eta * y[i]) {
        a_p = -eta * y[i] / dy[i];
      }
      if (a_d * dl[i] < -eta * l[i]) {
        a_d = -eta * l[i] / dl[i];
      }
    }
    alpha = (a_p < a_d ? a_p : a_d);
    x += alpha * (E_Q * dx);
    l += alpha * dl;
    y += alpha * dy;

    eta += ValueType(0.25) * (ValueType(1.0) - eta);
    if (++k > max_iter) {
      throw maximum_iteration(max_iter);
    }
    if (alpha * norm_2(dx) < abs_tol) {
      break;
    }
  } while (true);
}

}  // namespace ReaK::optim

#endif

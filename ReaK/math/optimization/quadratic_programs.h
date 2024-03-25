/**
 * \file quadratic_programs.h
 *
 * The following library provides implementations of quadratic programming algorithms.
 * The algorithm follows the specification given by Nocedal's Numerical Optimization book.
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

#ifndef REAK_MATH_OPTIMIZATION_QUADRATIC_PROGRAMS_H_
#define REAK_MATH_OPTIMIZATION_QUADRATIC_PROGRAMS_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/mat_qr_decomp.h"
#include "ReaK/math/lin_alg/mat_views.h"

#include "ReaK/math/optimization/optim_exceptions.h"

#include <vector>

namespace ReaK::optim {

/**
 * This function is an implementation of the null-space direct
 * method for solving a quadratic optimization problem with equality constraints.
 * It solves the following problem: \n
 * \n
 *           min c'x + 0.5 * x' G x \n
 *               Ax = b \n
 * \n
 * The implementation was inspired from the algorithm described in the book:\n
 *   Nocedal, Numerical Optimization, 2nd Ed..
 * TEST PASSED
 *
 * \param A The constraint matrix of dimension M*N.
 * \param b The b vector of dimension M.
 * \param G The G matrix of dimension NxN (defines the quadratic function to minimize, should be positive
 *semi-definite).
 * \param c The cost vector of dimension N.
 * \param x Stores, as output, the optimal vector.
 * \param abs_tol The tolerance on the singularity of components of the matrices involved.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, WritableVector Vector1,
          ReadableMatrix Matrix2, WritableVector Vector2>
void null_space_QP_method(
    const Matrix1& A, const Vector1& b, const Matrix2& G, const Vector2& c,
    Vector2& x,
    vect_value_type_t<Vector1> abs_tol =
        std::numeric_limits<vect_value_type_t<Vector1>>::epsilon(),
    vect_value_type_t<Vector1> max_norm =
        std::numeric_limits<vect_value_type_t<Vector1>>::infinity(),
    Vector1* lambda = nullptr) {
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

  mat_const_sub_block<mat<ValueType, mat_structure::square>> Y(Q, N, M, 0, 0);
  mat_const_sub_block<mat<ValueType, mat_structure::square>> Z(Q, N, N - M, 0,
                                                               M);

  Vector2 x_tmp(x);

  // forward-sub with L to find p_y
  mat_vect_adaptor<Vector2> p_y(x_tmp, M, 1, 0);
  ValueType L_abstrace(0.0);
  for (int i = 0; i < M; ++i) {
    p_y(i, 0) = b[i];
    L_abstrace += abs(L(i, i));
  }
  ReaK::detail::forwardsub_L_impl(L, p_y, abs_tol);

  // solve for p_z:
  mat_vect_adaptor<Vector2> p_z(x_tmp, N - M, 1, M);
  mat<ValueType, mat_structure::rectangular> GY_py = G * Y * p_y;
  for (int i = 0; i < N - M; ++i) {
    p_z(i, 0) = ValueType(0.0);
    for (int j = 0; j < N; ++j) {
      p_z(i, 0) -= Z(j, i) * (c[j] + GY_py(j, 0));
    }
  }
  mat<ValueType, mat_structure::symmetric> ZGZ =
      mat<ValueType, mat_structure::symmetric>(transpose_view(Z) * G * Z);

  try {
    linsolve_Cholesky(ZGZ, p_z, abs_tol);
  } catch (singularity_error&) {
    for (int i = 0; i < N - M; ++i) {
      p_z(i, 0) = ValueType(0.0);
      for (int j = 0; j < N; ++j) {
        p_z(i, 0) -= Z(j, i) * (c[j] + GY_py(j, 0));
      }
    }
    ValueType zz(0.0);
    ValueType zBz(0.0);
    for (int i = 0; i < N - M; ++i) {
      zz += p_z(i, 0) * p_z(i, 0);
      for (int j = 0; j < N - M; ++j) {
        zBz += p_z(i, 0) * (ZGZ(i, j) * p_z(j, 0));
      }
    }
    zBz = zz / (abs(zBz) + abs_tol);
    if (zBz < ValueType(1.0)) {
      zz *= zBz * zBz;
      for (int i = 0; i < N - M; ++i) {
        p_z(i, 0) = p_z(i, 0) * zBz;
      }
    }
    if (zz > max_norm * max_norm) {
      zBz = max_norm / abs(zz);
      for (int i = 0; i < N - M; ++i) {
        p_z(i, 0) = p_z(i, 0) * zBz;
      }
    }
  }

  x = Q * x_tmp;

  if (lambda) {
    (*lambda) = (c + G * x) * Y;
    mat_vect_adaptor<Vector1> lambda_mat(*lambda);
    ReaK::detail::backsub_R_impl(transpose_view(L), lambda_mat, abs_tol);
  }
}

/**
 * This function is an implementation of the null-space direct
 * method for solving a quadratic optimization problem with equality constraints.
 * It solves the following problem: \n
 * \n
 *           min c'x + 0.5 * x' G x \n
 *               Ax = b \n
 * \n
 * The implementation was inspired from the algorithm described in the book:\n
 *   Nocedal, Numerical Optimization, 2nd Ed..
 * TEST PASSED
 *
 * \param A The constraint matrix of dimension M*N.
 * \param b The b vector of dimension M.
 * \param G The G matrix of dimension NxN (defines the quadratic function to minimize, should be positive
 *semi-definite).
 * \param c The cost vector of dimension N.
 * \param x Stores, as output, the optimal vector.
 * \param abs_tol The tolerance on the singularity of components of the matrices involved.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, WritableVector Vector1,
          ReadableMatrix Matrix2, WritableVector Vector2>
void null_space_RRQP_method(
    const Matrix1& A, const Vector1& b, const Matrix2& G, const Vector2& c,
    Vector2& x,
    vect_value_type_t<Vector1> abs_tol =
        std::numeric_limits<vect_value_type_t<Vector1>>::epsilon(),
    vect_value_type_t<Vector1> max_norm =
        std::numeric_limits<vect_value_type_t<Vector1>>::infinity(),
    Vector1* lambda = nullptr) {
  using ValueType = vect_value_type_t<Vector1>;
  using std::abs;
  using std::swap;

  int N = c.size();
  int M = b.size();

  mat<ValueType, mat_structure::rectangular> R(transpose_view(A));
  mat<ValueType, mat_structure::square> Q =
      mat<ValueType, mat_structure::square>(
          mat<ValueType, mat_structure::identity>(N));
  mat<ValueType, mat_structure::permutation> P(N);
  int K = ReaK::detail::decompose_RRQR_impl(R, &Q, P, abs_tol);
  mat<ValueType, mat_structure::rectangular> L =
      mat<ValueType, mat_structure::rectangular>(transpose_view(R));
  L.set_col_count(M, true);

  mat_const_sub_block<mat<ValueType, mat_structure::square>> Y(Q, N, K, 0, 0);
  mat_const_sub_block<mat<ValueType, mat_structure::square>> Z(Q, N, N - K, 0,
                                                               K);
  Vector2 x_tmp(c);

  // forward-sub with L to find p_y
  mat_vect_adaptor<Vector2> p_y(x_tmp, K, 1, 0);
  for (int i = 0; i < K; ++i) {
    p_y(i, 0) = b[i];
  }
  ReaK::detail::forwardsub_L_impl(sub(L)(range(0, K), range(0, K)), p_y,
                                  abs_tol);

  // solve for p_z:
  mat_vect_adaptor<Vector2> p_z(x_tmp, N - K, 1, K);
  mat<ValueType, mat_structure::rectangular> GY_py = G * Y * p_y;
  for (int i = 0; i < N - K; ++i) {
    p_z(i, 0) = ValueType(0.0);
    for (int j = 0; j < N; ++j) {
      p_z(i, 0) -= Z(j, i) * (c[j] + GY_py(j, 0));
    }
  }
  mat<ValueType, mat_structure::symmetric> ZGZ =
      mat<ValueType, mat_structure::symmetric>(transpose_view(Z) * G * Z);
  try {
    linsolve_Cholesky(ZGZ, p_z, abs_tol);
  } catch (singularity_error&) {
    for (int i = 0; i < N - K; ++i) {
      p_z(i, 0) = ValueType(0.0);
      for (int j = 0; j < N; ++j) {
        p_z(i, 0) -= Z(j, i) * (c[j] + GY_py(j, 0));
      }
    }
    ValueType zz(0.0);
    ValueType zBz(0.0);
    for (int i = 0; i < N - K; ++i) {
      zz += p_z(i, 0) * p_z(i, 0);
      for (int j = 0; j < N - K; ++j) {
        zBz += p_z(i, 0) * (ZGZ(i, j) * p_z(j, 0));
      }
    }
    zBz = zz / (abs(zBz) + abs_tol);
    if (zBz < ValueType(1.0)) {
      zz *= zBz * zBz;
      for (int i = 0; i < N - K; ++i) {
        p_z(i, 0) = p_z(i, 0) * zBz;
      }
    }
    if (zz > max_norm * max_norm) {
      zBz = max_norm / abs(zz);
      for (int i = 0; i < N - K; ++i) {
        p_z(i, 0) = p_z(i, 0) * zBz;
      }
    }
  }
  x = Q * x_tmp;
  x = P * x;
  if (lambda) {
    vect_n<ValueType> lam((c + G * x) * Y);
    mat_vect_adaptor<vect_n<ValueType>> lam_mat(lam);
    ReaK::detail::backsub_R_impl(
        sub(transpose_view(L))(range(0, K), range(0, K)), lam_mat, abs_tol);
    for (int i = 0; i < K; ++i) {
      (*lambda)[i] = lam[i];
    }
    for (int i = K; i < M; ++i) {
      (*lambda)[i] = 0.0;
    }
  }
}

/**
 * This function is an implementation of the projected Conjugate Gradient
 * method (with normal equations approach) for solving a quadratic optimization problem
 * with equality constraints. It solves the following problem: \n
 * \n
 *           min c'x + 0.5 * x' G x \n
 *               Ax = b \n
 * \n
 * The implementation was inspired from the algorithm described in the book:\n
 *   Nocedal, Numerical Optimization, 2nd Ed..
 * TEST PASSED
 *
 * \param A The constraint matrix of dimension M*N.
 * \param b The b vector of dimension M.
 * \param G The G matrix of dimension NxN (defines the quadratic function to minimize, should be positive
 *semi-definite).
 * \param c The cost vector of dimension N.
 * \param x Stores, as output, the optimal vector.
 * \param tol The relative tolerance on the singularity of components of the matrices involved.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, WritableVector Vector1,
          ReadableMatrix Matrix2, WritableVector Vector2>
void projected_CG_method(
    const Matrix1& A, const Vector1& b, const Matrix2& G, const Vector2& c,
    Vector2& x, unsigned int max_iter = 20,
    vect_value_type_t<Vector1> abs_tol =
        std::numeric_limits<vect_value_type_t<Vector1>>::epsilon(),
    Vector1* lambda = nullptr) {
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

  Vector1 b_tmp = b;
  mat_vect_adaptor<Vector1> b_tmp_mat(b_tmp);
  ReaK::detail::backsub_Cholesky_impl(L, b_tmp_mat);
  x = b_tmp * A;

  Vector2 r = G * x + c;
  Vector1 Ar = A * r;
  mat_vect_adaptor<Vector1> Ar_mat(Ar);
  ReaK::detail::backsub_Cholesky_impl(L, Ar_mat);
  Vector2 g = r;
  g -= Ar * A;
  Vector2 d = -g;
  Vector2 Gd = G * d;
  ValueType rg = r * g;

  unsigned int k = 0;
  while (abs(rg) > abs_tol) {
    ValueType alpha = rg / (d * Gd);
    x += alpha * d;
    r += alpha * Gd;
    Ar = A * r;
    ReaK::detail::backsub_Cholesky_impl(L, Ar_mat);
    g = r;
    g -= Ar * A;
    ValueType rg_p = r * g;
    ValueType beta = rg_p / rg;
    rg = rg_p;
    d *= beta;
    d -= g;
    Gd = G * d;
    if (++k > max_iter) {
      throw maximum_iteration(max_iter);
    }
  }

  if (lambda) {
    (*lambda) = (c + G * x) *
                mat_const_sub_block<mat<ValueType, mat_structure::square>>(
                    Q, N, M, 0, 0);
    mat_vect_adaptor<Vector1> lambda_mat(*lambda);
    ReaK::detail::backsub_R_impl(transpose_view(L), lambda_mat, abs_tol);
  }
}

}  // namespace ReaK::optim

#endif  // REAK_MATH_OPTIMIZATION_QUADRATIC_PROGRAMS_H_

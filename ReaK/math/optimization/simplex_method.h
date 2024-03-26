/**
 * \file simplex_method.h
 *
 * The following library is an implementation of the Simplex Method to solve a linear programming
 * problem. The algorithm follows that of Chvatal and Vasek 1983. This simplex method should not
 * be confused with the Nelder-Mead method which is also sometimes referred to as the simplex method
 * (see nelder_mead_method.hpp).
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

#ifndef REAK_MATH_OPTIMIZATION_SIMPLEX_METHOD_H_
#define REAK_MATH_OPTIMIZATION_SIMPLEX_METHOD_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_qr_decomp.h"

#include "ReaK/math/optimization/optim_exceptions.h"

#include <vector>

namespace ReaK::optim {

namespace detail {

template <typename Matrix, typename Vector, typename T>
void simplex_method_loop_impl(const Matrix& A, Matrix& B, const Vector& c,
                              Vector& c_B, Vector& x, const Vector& l,
                              const Vector& u, std::vector<int>& i_b,
                              std::vector<int>& i_n, T tol) {
  using ValueType = vect_value_type_t<Vector>;
  using std::swap;

  int N = i_n.size();
  int M = A.get_row_count();

  Matrix B_Q(M, M);
  Matrix B_R(M, M);
  decompose_QR(B, B_Q, B_R, tol);

  while (true) {
    // Step 1
    Vector y = c_B * B_Q;
    mat_vect_adaptor<Vector> y_mat(y);
    ReaK::detail::backsub_R_impl(B_R, y_mat, tol);

    // Step 2
    ValueType sum;
    int enter_var = N;
    for (int i = 0; i < N; ++i) {
      sum = y * slice(A)(range(int(0), M), i_n[i]);
      if (((sum < c[i_n[i]]) && (x[i_n[i]] < u[i_n[i]])) ||
          ((sum > c[i_n[i]]) && (x[i_n[i]] > l[i_n[i]]))) {
        enter_var = i;
        break;
      }
    }
    if (enter_var == N) {
      return;
    }
    // Step 3
    y = slice(A)(range(int(0), M), i_n[enter_var]);
    y = y * B_Q;
    ReaK::detail::backsub_R_impl(B_R, y_mat, tol);
    // Step 4
    if (sum < c[i_n[enter_var]]) {
      ValueType t_max = std::numeric_limits<ValueType>::infinity();
      int leave_var = 0;
      for (int i = 0; i < M; ++i) {
        if (y[i] > 0.0) {
          if (t_max * y[i] > x[i_b[i]] - l[i_b[i]]) {
            t_max = (x[i_b[i]] - l[i_b[i]]) / y[i];
            leave_var = i;
          }
        } else if (y[i] < 0.0) {
          if (t_max * y[i] < x[i_b[i]] - u[i_b[i]]) {
            t_max = (x[i_b[i]] - u[i_b[i]]) / y[i];
            leave_var = i;
          }
        }
      }
      if (t_max > u[i_n[enter_var]] - x[i_n[enter_var]]) {
        t_max = u[i_n[enter_var]] - x[i_n[enter_var]];
        x[i_n[enter_var]] = u[i_n[enter_var]];
        for (int i = 0; i < M; ++i) {
          x[i_b[i]] -= t_max * y[i];
        }
        continue;
      }
      if (t_max != std::numeric_limits<ValueType>::infinity()) {
        x[i_n[enter_var]] += t_max;
        for (int i = 0; i < M; ++i) {
          x[i_b[i]] -= t_max * y[i];
        }
        slice(B)(range(0, M), leave_var) =
            slice(A)(range(0, M), i_n[enter_var]);
        decompose_QR(B, B_Q, B_R, tol);
        c_B[leave_var] = c[i_n[enter_var]];
        swap(i_b[leave_var], i_n[enter_var]);
      } else {
        throw unbounded_problem(
            "Simplex method failed due to an unbounded search domain!");
      }
    } else if (sum > c[i_n[enter_var]]) {
      ValueType t_max = std::numeric_limits<ValueType>::infinity();
      int leave_var = 0;
      for (int i = 0; i < M; ++i) {
        if (y[i] > 0.0) {
          if (t_max * y[i] > u[i_b[i]] - x[i_b[i]]) {
            t_max = (u[i_b[i]] - x[i_b[i]]) / y[i];
            leave_var = i;
          }
        } else if (y[i] < 0.0) {
          if (t_max * y[i] < l[i_b[i]] - x[i_b[i]]) {
            t_max = (l[i_b[i]] - x[i_b[i]]) / y[i];
            leave_var = i;
          }
        }
      }
      if (t_max > x[i_n[enter_var]] - l[i_n[enter_var]]) {
        t_max = x[i_n[enter_var]] - l[i_n[enter_var]];
        x[i_n[enter_var]] = l[i_n[enter_var]];
        for (int i = 0; i < M; ++i) {
          x[i_b[i]] += t_max * y[i];
        }
        continue;
      }
      if (t_max != std::numeric_limits<ValueType>::infinity()) {
        x[i_n[enter_var]] -= t_max;
        for (int i = 0; i < M; ++i) {
          x[i_b[i]] += t_max * y[i];
        }
        slice(B)(range(0, M), leave_var) =
            slice(A)(range(0, M), i_n[enter_var]);
        decompose_QR(B, B_Q, B_R, tol);
        c_B[leave_var] = c[i_n[enter_var]];
        swap(i_b[leave_var], i_n[enter_var]);
      } else {
        throw unbounded_problem(
            "Simplex method failed due to an unbounded search domain!");
      }
    }
  }
}
}  // namespace detail

/**
 * This function is an implementation of the general two-phase revised simplex
 * method for bounded variables. It solves the following problem: \n
 * \n
 *           max c'x \n
 *               Ax = b \n
 *             l <= x <= u \n
 * \n
 * The implementation was inspired from the algorithm described in the book:\n
 *   Chvatal, Vasek, Linear Programming, W. H. Freeman and Company, 1983.
 * \test Must create a unit-test for this. So far, this method fails the tests.
 *
 * \param A The constraint matrix of dimension M*N.
 * \param b The b vector of dimension M.
 * \param c The cost vector of dimension N.
 * \param x0 The initial guess for the optimal vector, then stores, as output, the optimal vector.
 * \param l The lower-bound on the independent variables, if there is no lower bound for a variable, set to minus
 *infinity.
 * \param u The upper-bound on the independent variables, if there is no upper bound for a variable, set to infinity.
 *
 * \throw
 *
 * \author Mikael Persson
 */
template <FullyWritableMatrix Matrix, WritableVector Vector>
void simplex_method(
    const Matrix& A, const Vector& b, const Vector& c, Vector& x0,
    const Vector& l, const Vector& u,
    vect_value_type_t<Vector> tol =
        std::numeric_limits<vect_value_type_t<Vector>>::epsilon()) {
  using ValueType = vect_value_type_t<Vector>;
  using std::abs;
  using std::swap;

  int N = c.size();
  int M = b.size();

  Vector x(N + M);
  std::copy(x0.begin(), x0.end(), x.begin());

  Vector b_G = b;

  Vector c_G(N + M, 0.0);
  Vector c_B(M);

  Vector l_G(N + M);
  Vector u_G(N + M);
  std::copy(l.begin(), l.end(), l_G.begin());
  std::copy(u.begin(), u.end(), u_G.begin());

  Matrix A_G(M, N + M);
  sub(A_G)(range(0, M), range(0, N)) = A;
  sub(A_G)(range(0, M), range(N, N + M)) =
      mat<ValueType, mat_structure::identity>(M);

  Matrix B = Matrix(mat<ValueType, mat_structure::identity>(M));
  Matrix B_Q = Matrix(mat<ValueType, mat_structure::identity>(M));
  Matrix B_R = Matrix(mat<ValueType, mat_structure::identity>(M));

  std::vector<int> i_b(M);
  for (int i = 0; i < M; ++i) {
    i_b[i] = i + N;
  }
  std::vector<int> i_n(N);
  for (int i = 0; i < N; ++i) {
    i_n[i] = i;
  }

  for (int i = 0; i < N; ++i) {
    if (x[i] < l[i]) {
      x[i] = l[i];
    }
    if (x[i] > u[i]) {
      x[i] = u[i];
    }
  }

  bool FeasibleStart = true;
  for (int i = 0; i < M; ++i) {
    x[N + i] = b_G[i];
    for (int j = 0; j < N; ++j) {
      x[N + i] -= A(i, j) * x[j];
    }
    if (FeasibleStart) {
      FeasibleStart = (abs(x[N + i]) < tol);
    }
    if (x[N + i] >= 0.0) {
      l_G[N + i] = 0.0;
      u_G[N + i] = std::numeric_limits<ValueType>::infinity();
    } else {
      l_G[N + i] = -std::numeric_limits<ValueType>::infinity();
      u_G[N + i] = 0.0;
    }
    c_G[N + i] = -1.0;
    c_B[i] = -1.0;
  }

  if (!FeasibleStart) {
    // First-Phase
    detail::simplex_method_loop_impl(A_G, B, c_G, c_B, x, l_G, u_G, i_b, i_n,
                                     tol);

    // Did the first phase succeed?
    for (int i = 0; i < M; ++i) {
      if (abs(x[N + i]) > tol) {
        throw infeasible_problem(
            "Simplex method failed due to an empty search domain! No feasible "
            "solution exists!");
      }
    }
  }

  // Getting Rid of the Artificial Variables
  decompose_QR(B, B_Q, B_R, tol);
  for (int i = 0; i < M; ++i) {
    if (i_b[i] >= N) {
      Vector y;
      y = slice(A_G)(range(0, M), i_b[i]);
      y = y * B_Q;
      mat_vect_adaptor<Vector> y_mat(y);
      ReaK::detail::backsub_R_impl(B_R, y_mat, tol);
      for (int j = 0; j < N; ++j) {
        ValueType sum = y * slice(A_G)(range(0, M), i_n[j]);
        if ((sum != 0.0) && (i_n[j] < N)) {
          slice(B)(range(0, M), i) = slice(A_G)(range(0, M), i_n[j]);
          swap(i_b[i], i_n[j]);
          break;
        }
      }
      decompose_QR(B, B_Q, B_R, tol);
    }
  }
  // Getting Rid of the Redundant equations
  int RedundantCount = 0;
  for (int i = 0; i < M; ++i) {
    if (i_b[i] >= N) {
      // Must Delete Redundant Equation
      ++RedundantCount;
      --M;
      Matrix tempA(M, N);
      Matrix tempB(M, M);
      Vector tempb(M);

      tempB = ((sub(B)(range(0, i_b[i] - N), range(0, i)) &
                sub(B)(range(0, i_b[i] - N), range(i + 1, M + 1))) |
               (sub(B)(range(i_b[i] - N + 1, M + 1), range(0, i)) &
                sub(B)(range(i_b[i] - N + 1, M + 1), range(i + 1, M + 1))));

      tempA = ((sub(A_G)(range(0, i_b[i] - N), range(0, N))) |
               (sub(A_G)(range(i_b[i] - N + 1, M + 1), range(0, N))));

      std::copy(b_G.begin(), b_G.begin() + i_b[i] - N, tempb.begin());
      std::copy(b_G.begin() + i_b[i] - N + 1, b_G.end(),
                tempb.begin() + i_b[i] - N);
      swap(tempb, b_G);
      swap(tempA, A_G);
      swap(tempB, B);
    }
  }
  decompose_QR(B, B_Q, B_R, tol);

  // Prepare the variables for the second phase
  {
    x.resize(N);

    c_B.resize(M);
    for (int i = 0; i < M; ++i) {
      c_B[i] = c[i_b[i]];
    }

    std::vector<int> tempIB(M);
    std::vector<int> tempIN(N - M);
    {
      int j = 0;
      for (int i = 0; i < M + RedundantCount; ++i) {
        if (i_b[i] < N) {
          tempIB[j] = i_b[i];
          ++j;
        }
      }
      j = 0;
      for (int i = 0; i < N; ++i) {
        if (i_n[i] < N) {
          tempIN[j] = i_n[i];
          j++;
        }
      }
    }
    i_b.swap(tempIB);
    i_n.swap(tempIN);
  }

  detail::simplex_method_loop_impl(A_G, B, c, c_B, x, l, u, i_b, i_n, tol);

  x0 = x;
}

}  // namespace ReaK::optim

#endif  // REAK_MATH_OPTIMIZATION_SIMPLEX_METHOD_H_

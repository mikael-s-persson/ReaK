/**
 * \file mat_schur_decomp.h
 *
 * This library provides function templates to compute the Real-Schur Decomposition
 * and the Generalized Real-Schur Decomposition (or QZ decomposition). These algorithms
 * can be used as the main step of several algorithms including the QR algorithm and
 * generalized eigen-problems. Schur decompositions are also useful for several other
 * advanced matrix numerical methods (some which are planned for the future in this library).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2011
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

#ifndef REAK_MATH_LIN_ALG_MAT_SCHUR_DECOMP_H_
#define REAK_MATH_LIN_ALG_MAT_SCHUR_DECOMP_H_

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/mat_num_exceptions.h"

#include "ReaK/math/lin_alg/mat_hess_decomp.h"
#include "ReaK/math/lin_alg/mat_householder.h"

#include <type_traits>

namespace ReaK {

namespace detail {

/* TESTED and works. Golub & vanLoan Alg-8.3.2 */
template <typename Matrix1, typename Matrix2>
void symmetric_QR_step(Matrix1& T, Matrix2* Z) {
  using ValueType = mat_value_type_t<Matrix1>;
  using std::sqrt;

  int N = T.get_row_count();

  ValueType d = 0.5 * (T(N - 2, N - 2) - T(N - 1, N - 1));
  ValueType tsqr = T(N - 1, N - 2) * T(N - 1, N - 2);
  ValueType mu = T(N - 1, N - 1) -
                 tsqr / (d + (d > 0.0 ? 1.0 : -1.0) * sqrt(d * d + tsqr));
  ValueType x = T(0, 0) - mu;
  ValueType z = T(1, 0);
  for (int k = 0; k < N - 1; ++k) {
    givens_rot_matrix<ValueType> g(x, z,
                                   std::numeric_limits<ValueType>::epsilon());

    mat_sub_block<Matrix1> subT1(
        T, 2, (k < N - 2 ? (k > 0 ? 4 : 3) : (k > 0 ? 3 : 2)), k,
        (k > 0 ? k - 1 : 0));
    givens_rot_prod(g, subT1);  // G * A

    mat_sub_block<Matrix1> subT2(
        T, (k < N - 2 ? (k > 0 ? 4 : 3) : (k > 0 ? 3 : 2)), 2,
        (k > 0 ? k - 1 : 0), k);
    givens_rot_prod(subT2, transpose(g));  // A * G^T

    if (Z) {
      mat_sub_block<Matrix2> subZ(*Z, Z->get_row_count(), 2, 0, k);
      givens_rot_prod(subZ, transpose(g));  // Q_prev * P
    }

    if (k < N - 2) {
      x = T(k + 1, k);
      z = T(k + 2, k);
    }
  }
}

template <typename Matrix1, typename Matrix2>
void symmetric_QRalg_impl(Matrix1& T, Matrix2* Z,
                          mat_value_type_t<Matrix1> NumTol = 1E-8) {
  using ValueType = mat_value_type_t<Matrix1>;
  using std::abs;
  using std::sqrt;

  int N = T.get_row_count();

  ValueType absNumTol = 0.0;
  for (int i = 0; i < N; ++i) {
    absNumTol = abs(T(i, i));
  }
  absNumTol /= N;
  // Best relative tolerance is the square root of epsilon for the scalar type
  // since the QR step involves a division by the square-root off-diagonal term.
  absNumTol *=
      std::max(NumTol, sqrt(std::numeric_limits<ValueType>::epsilon()));

  detail::decompose_TriDiag_impl(T, Z, absNumTol);

  int q = N;
  int p = N;

  while (q > 0) {
    int i = N - 1;
    // find a trailing diagonal sub-matrix
    for (; i > 0; --i) {
      if (abs(T(i, i - 1)) < absNumTol) {
        T(i, i - 1) = ValueType(0.0);
      } else {
        break;
      }
    }
    q = i;
    if (i == 0) {  // break if it is entirely diagonal.
      break;
    }

    // find the middle, biggest unreduced tridiagonal matrix
    for (; i > 0; --i) {
      if (abs(T(i, i - 1)) < absNumTol) {
        T(i, i - 1) = ValueType(0.0);
        break;
      }
    }
    p = i;

    // set remaining sub-diagonals to zero if they are very small.
    for (i = p; i > 0; --i) {
      if (abs(T(i, i - 1)) < absNumTol) {
        T(i, i - 1) = ValueType(0.0);
      }
    }

    if (Z) {
      mat_sub_block<Matrix2> subZ(*Z, Z->get_row_count(), q - p + 1, 0, p);
      mat_sub_block<Matrix1> subT(T, q - p + 1, q - p + 1, p, p);
      symmetric_QR_step(subT, &subZ);
    } else {
      mat_sub_block<Matrix1> subT(T, q - p + 1, q - p + 1, p, p);
      symmetric_QR_step(subT, static_cast<Matrix2*>(nullptr));
    }
  }
}

template <typename Matrix1, typename Matrix2>
void francis_QR_step(Matrix1& H, Matrix2* Z, mat_size_type_t<Matrix1> Offset,
                     mat_value_type_t<Matrix1> NumTol) {
  using ValueType = mat_value_type_t<Matrix1>;

  householder_matrix<vect<ValueType, 3>> hhm;

  int N = H.get_row_count() - Offset;
  if (N < 3) {
    return;
  }

  ValueType s = H(Offset + N - 2, N - 2) + H(Offset + N - 1, N - 1);
  ValueType t = H(Offset + N - 2, N - 2) * H(Offset + N - 1, N - 1) -
                H(Offset + N - 2, N - 1) * H(Offset + N - 1, N - 2);
  vect<ValueType, 3> v(H(Offset, 0) * H(Offset, 0) +
                           H(Offset, 1) * H(Offset + 1, 0) - s * H(Offset, 0) +
                           t,
                       H(Offset + 1, 0) * (H(Offset, 0) + H(Offset + 1, 1) - s),
                       H(Offset + 1, 0) * H(Offset + 2, 1));
  for (int k = 0; k < N - 2; ++k) {
    hhm.set(v, NumTol);

    int q = (k == 0 ? 0 : k - 1);
    mat_sub_block<Matrix1> subH1(H, 3, H.get_col_count() - q, Offset + k, q);
    householder_prod(hhm, subH1);

    int r = (k == N - 3 ? N : k + 4);
    mat_sub_block<Matrix1> subH2(H, Offset + r, 3, 0, k);
    householder_prod(subH2, hhm);

    if (Z) {
      mat_sub_block<Matrix2> subZ(*Z, Z->get_row_count(), 3, 0, k);
      householder_prod(subZ, hhm);  // Q_prev * P
    }

    v[0] = H(Offset + k + 1, k);
    v[1] = H(Offset + k + 2, k);
    if (k < N - 3) {
      v[2] = H(Offset + k + 3, k);
    }
  }

  householder_matrix<vect<ValueType, 2>> hhm2(vect<ValueType, 2>(v[0], v[1]),
                                              NumTol);
  mat_sub_block<Matrix1> subH1(H, 2, H.get_col_count() - N + 3, Offset + N - 2,
                               N - 3);
  householder_prod(hhm2, subH1);
  mat_sub_block<Matrix1> subH2(H, Offset + N, 2, 0, N - 2);
  householder_prod(subH2, hhm2);

  if (Z) {
    mat_sub_block<Matrix2> subZ(*Z, Z->get_row_count(), 2, 0, N - 2);
    householder_prod(subZ, hhm2);  // Q_prev * P
  }
}

template <typename Matrix1, typename Matrix2>
void schur_decomp_impl(Matrix1& T, Matrix2* Q,
                       mat_value_type_t<Matrix1> NumTol = 1E-8) {
  using ValueType = mat_value_type_t<Matrix1>;
  using std::abs;

  int N = T.get_row_count();

  ValueType absNumTol = 0.0;
  for (int i = 0; i < N; ++i) {
    absNumTol = abs(T(i, i));
  }
  absNumTol *= NumTol / N;

  detail::decompose_Hess_impl(T, Q, absNumTol);

  int q = N;
  int p = N;

  while (q > 0) {
    bool last_off_diag_was_nil = true;
    int i = N - 1;

    // find a trailing quasi-upper-triangular sub-matrix
    for (; i > 0; --i) {
      if (abs(T(i, i - 1)) < absNumTol) {
        T(i, i - 1) = ValueType(0.0);
        q = i;
        last_off_diag_was_nil = true;
      } else {
        if (!last_off_diag_was_nil) {
          break;
        }
        last_off_diag_was_nil = false;
      }
    }
    if (i == 0) {  // break if it is entirely quasi-upper-triangular.
      break;
    }

    // find the middle, biggest unreduced upper-Hessenberg matrix
    p = 0;
    for (i = q - 1; i > 0; --i) {
      if (abs(T(i, i - 1)) < absNumTol) {
        T(i, i - 1) = ValueType(0.0);
        p = i;
        break;
      }
    }

    // set remaining sub-diagonals to zero if they are very small.
    for (i = p; i > 0; --i) {
      if (abs(T(i, i - 1)) < absNumTol) {
        T(i, i - 1) = ValueType(0.0);
      }
    }

    if (Q) {
      mat_sub_block<Matrix2> subQ(*Q, Q->get_row_count(),
                                  Q->get_col_count() - p, 0, p);
      mat_sub_block<Matrix1> subT(T, q, N - p, 0, p);
      francis_QR_step(subT, &subQ, p, absNumTol);
    } else {
      mat_sub_block<Matrix1> subT(T, q, N - p, 0, p);
      francis_QR_step(subT, static_cast<Matrix2*>(nullptr), p, absNumTol);
    }
  }
}

/*
 * This function deflates the shifted diagonal element (q,p) down to the lower right corner
 * of the upper-triangular matrix B while zeroing the lower right sub-diagonal element of the
 * upper-Hessenberg matrix A.
 * \param A The upper-Hessenberg matrix of the GEP (lhs), resulting in an upper-Hessenberg matrix
 *          with a zero at the lower right sub-diagonal term.
 * \param B The upper-triangular matrix of the GEP (rhs) with a zero value at the entry (q,p),
 *          resulting in an upper-triangular matrix with a zero at the lower right diagonal term.
 * \param Q Pointer to the matrix in which to accumulate the orthogonal pre-multiplied transformations.
 * \param Z Pointer to the matrix in which to accumulate the orthogonal post-multiplied transformations.
 * \param p The column-index of the diagonal element of B to be chased down to the lower-right corner.
 * \param q The row-index of the diagonal element of B to be chased down to the lower-right corner. Note
 *          that the reason for this parameter is that this function can be applied to the case where
 *          A and B include the upper-right part of a larger matrix and thus, the "diagonal" is actually
 *          shifted down. This was done because this deflation algorithm is called by other decomposition
 *          algorithms (e.g. QZ-algorithm) which operate on the sub-range of the diagonal but still want
 *          the orthogonal transformations to be applied to all the affected elements (upper-right part).
 */
template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4>
void deflate_hess_elem_down_impl(Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
                                 mat_size_type_t<Matrix1> p,
                                 mat_size_type_t<Matrix1> q) {
  using ValueType = mat_value_type_t<Matrix1>;
  int N = A.get_row_count();

  givens_rot_matrix<ValueType> G;

  // int t = N-2;

  for (; q < N - 1; ++p, ++q) {
    G.set(B(q, p + 1), B(q + 1, p + 1));

    mat_sub_block<Matrix1> subA1(
        A, 2, (p == 0 ? A.get_col_count() : A.get_col_count() - p + 1), q,
        (p == 0 ? p : p - 1));
    givens_rot_prod(G, subA1);  // G * A

    mat_sub_block<Matrix2> subB1(B, 2, B.get_col_count() - p - 1, q, p + 1);
    givens_rot_prod(G, subB1);  // G * B

    if (Q) {
      mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 2, 0, q);
      givens_rot_prod(subQ, transpose(G));  // Q_prev * G^T
    }

    if (p == 0) {
      continue;
    }

    G.set(-A(q + 1, p), A(q + 1, p - 1));
    G = transpose(G);

    mat_sub_block<Matrix2> subB2(B, q, 2, 0, p - 1);
    givens_rot_prod(subB2, G);  // B * G^T

    mat_sub_block<Matrix1> subA2(A, q + 2, 2, 0, p - 1);
    givens_rot_prod(subA2, G);  // A * G^T

    if (Z) {
      mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, p - 1);
      givens_rot_prod(subZ, G);  // Q_prev * G^T
    }
  }

  G.set(-A(N - 1, p), A(N - 1, p - 1));
  G = transpose(G);

  mat_sub_block<Matrix2> subB3(B, N - 1, 2, 0, p - 1);
  givens_rot_prod(subB3, G);  // B * G^T

  mat_sub_block<Matrix1> subA3(A, N, 2, 0, p - 1);
  givens_rot_prod(subA3, G);  // A * G^T

  if (Z) {
    mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, p - 1);
    givens_rot_prod(subZ, G);  // Q_prev * G^T
  }
}

template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4>
bool deflate_hess_all_down_impl(Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
                                mat_size_type_t<Matrix1> row_offset,
                                mat_value_type_t<Matrix1> NumTol) {
  using std::abs;
  int N = A.get_row_count();
  int M = N - row_offset;

  bool deflated_prev_time = true;
  for (int j = 0; j < M && deflated_prev_time; ++j) {
    deflated_prev_time = false;
    for (int i = M - j; i > 0; --i) {
      if (abs(B(row_offset + i - 1, i - 1)) < NumTol) {
        mat_sub_block<Matrix1> subA(A, N - j, A.get_col_count(), 0, 0);
        mat_sub_block<Matrix2> subB(B, N - j, A.get_col_count(), 0, 0);
        if (Q) {
          mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), N - j, 0, 0);
          if (Z) {
            mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), N - j, 0, 0);
            deflate_hess_elem_down_impl(subA, subB, &subQ, &subZ, i - 1,
                                        row_offset + i - 1);
          } else {
            deflate_hess_elem_down_impl(subA, subB, &subQ,
                                        static_cast<Matrix4*>(nullptr), i - 1,
                                        row_offset + i - 1);
          }
        } else {
          if (Z) {
            mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), N - j, 0, 0);
            deflate_hess_elem_down_impl(subA, subB,
                                        static_cast<Matrix3*>(nullptr), &subZ,
                                        i - 1, row_offset + i - 1);
          } else {
            deflate_hess_elem_down_impl(
                subA, subB, static_cast<Matrix3*>(nullptr),
                static_cast<Matrix4*>(nullptr), i - 1, row_offset + i - 1);
          }
        }
        deflated_prev_time = true;
        break;
      }
    }
    if (j == 0 && !deflated_prev_time) {
      return false;
    }
  }
  return true;
}

/*
 * This function deflates the shifted diagonal element (q,p) up to the upper left corner
 * of the upper-triangular matrix B while zeroing the upper left sub-diagonal element of the
 * upper-Hessenberg matrix A.
 * \param A The upper-Hessenberg matrix of the GEP (lhs), resulting in an upper-Hessenberg matrix
 *          with a zero at the lower right sub-diagonal term.
 * \param B The upper-triangular matrix of the GEP (rhs) with a zero value at the entry (q,p),
 *          resulting in an upper-triangular matrix with a zero at the lower right diagonal term.
 * \param Q Pointer to the matrix in which to accumulate the orthogonal pre-multiplied transformations.
 * \param Z Pointer to the matrix in which to accumulate the orthogonal post-multiplied transformations.
 * \param p The column-index of the diagonal element of B to be chased down to the lower-right corner.
 * \param q The row-index of the diagonal element of B to be chased down to the lower-right corner. Note
 *          that the reason for this parameter is that this function can be applied to the case where
 *          A and B include the upper-right part of a larger matrix and thus, the "diagonal" is actually
 *          shifted down. This was done because this deflation algorithm is called by other decomposition
 *          algorithms (e.g. QZ-algorithm) which operate on the sub-range of the diagonal but still want
 *          the orthogonal transformations to be applied to all the affected elements (upper-right part).
 */
template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4>
void deflate_hess_elem_up_impl(Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
                               mat_size_type_t<Matrix1> p,
                               mat_size_type_t<Matrix1> q) {
  using ValueType = mat_value_type_t<Matrix1>;
  int N = A.get_row_count();

  givens_rot_matrix<ValueType> G;

  // int t = N-2;

  for (; p > 0; --p, --q) {
    G.set(-B(q - 1, p), B(q - 1, p - 1));
    G = transpose(G);

    mat_sub_block<Matrix2> subB2(B, q, 2, 0, p - 1);
    givens_rot_prod(subB2, G);  // B * G^T

    mat_sub_block<Matrix1> subA2(A, (q == N - 1 ? N : q + 2), 2, 0, p - 1);
    givens_rot_prod(subA2, G);  // A * G^T

    if (Z) {
      mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, p - 1);
      givens_rot_prod(subZ, G);  // Q_prev * G^T
    }

    if (q == N - 1) {
      continue;
    }

    G.set(A(q, p - 1), A(q + 1, p - 1));

    mat_sub_block<Matrix1> subA1(A, 2, A.get_col_count() - p + 1, q, p - 1);
    givens_rot_prod(G, subA1);  // G * A

    mat_sub_block<Matrix2> subB1(B, 2, B.get_col_count() - p - 1, q, p + 1);
    givens_rot_prod(G, subB1);  // G * B

    if (Q) {
      mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 2, 0, q);
      givens_rot_prod(subQ, transpose(G));  // Q_prev * G^T
    }
  }

  G.set(A(q, p), A(q + 1, p));

  mat_sub_block<Matrix1> subA1(A, 2, A.get_col_count() - p, q, p);
  givens_rot_prod(G, subA1);  // G * A

  mat_sub_block<Matrix2> subB1(B, 2, B.get_col_count() - p - 1, q, p + 1);
  givens_rot_prod(G, subB1);  // G * B

  if (Q) {
    mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 2, 0, q);
    givens_rot_prod(subQ, transpose(G));  // Q_prev * G^T
  }
}

template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4>
bool deflate_hess_all_up_impl(Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
                              mat_size_type_t<Matrix1> row_offset,
                              mat_value_type_t<Matrix1> NumTol) {
  using std::abs;
  int N = A.get_row_count();
  int M = N - row_offset;

  bool deflated_prev_time = true;
  for (int j = 0; j < M && deflated_prev_time; ++j) {
    deflated_prev_time = false;
    for (int i = j; i < M; ++i) {
      if (abs(B(row_offset + i, i)) < NumTol) {
        mat_sub_block<Matrix1> subA(A, N, A.get_col_count() - j, 0, j);
        mat_sub_block<Matrix2> subB(B, N, A.get_col_count() - j, 0, j);
        if (Q) {
          mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), N, 0, 0);
          if (Z) {
            mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), N - j, 0, j);
            deflate_hess_elem_up_impl(subA, subB, &subQ, &subZ, i - j,
                                      row_offset + i);
          } else {
            deflate_hess_elem_up_impl(subA, subB, &subQ,
                                      static_cast<Matrix4*>(nullptr), i - j,
                                      row_offset + i);
          }
        } else {
          if (Z) {
            mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), N - j, 0, 0);
            deflate_hess_elem_up_impl(subA, subB,
                                      static_cast<Matrix3*>(nullptr), &subZ,
                                      i - j, row_offset + i);
          } else {
            deflate_hess_elem_up_impl(
                subA, subB, static_cast<Matrix3*>(nullptr),
                static_cast<Matrix4*>(nullptr), i - j, row_offset + i);
          }
        }
        deflated_prev_time = true;
        break;
      }
    }
    if (j == 0 && !deflated_prev_time) {
      return false;
    }
  }
  return true;
}

template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4>
void francis_QZ_step_impl(Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
                          mat_size_type_t<Matrix1> Offset,
                          mat_size_type_t<Matrix1> N,
                          vect<mat_value_type_t<Matrix1>, 3> v,
                          mat_value_type_t<Matrix1> NumTol) {
  using ValueType = mat_value_type_t<Matrix1>;

  householder_matrix<vect<ValueType, 3>> hhm;

  for (int k = 0; k < N - 2; ++k) {
    hhm.set(v, NumTol);

    mat_sub_block<Matrix1> subA1(A, 3, A.get_col_count(), Offset + k, 0);
    householder_prod(hhm, subA1);

    mat_sub_block<Matrix2> subB1(B, 3, B.get_col_count(), Offset + k, 0);
    householder_prod(hhm, subB1);

    if (Q) {
      mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 3, 0, k);
      householder_prod(subQ, hhm);  // Q * Q_k ^T
    }

    vect<ValueType, 3> v2(B(Offset + k + 2, k), B(Offset + k + 2, k + 1),
                          B(Offset + k + 2, k + 2));
    hhm.set_inv(v2, NumTol);

    mat_sub_block<Matrix1> subA2(A, A.get_row_count(), 3, 0, k);
    householder_prod(subA2, hhm);

    mat_sub_block<Matrix2> subB2(B, B.get_row_count(), 3, 0, k);
    householder_prod(subB2, hhm);

    if (Z) {
      mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 3, 0, k);
      householder_prod(subZ, hhm);  // Q_prev * P
    }

    householder_matrix<vect<ValueType, 2>> hhm2;
    hhm2.set_inv(
        vect<ValueType, 2>(B(Offset + k + 1, k), B(Offset + k + 1, k + 1)),
        NumTol);

    mat_sub_block<Matrix1> subA3(A, A.get_row_count(), 2, 0, k);
    householder_prod(subA3, hhm2);

    mat_sub_block<Matrix2> subB3(B, B.get_row_count(), 2, 0, k);
    householder_prod(subB3, hhm2);

    if (Z) {
      mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, k);
      householder_prod(subZ, hhm2);  // Q_prev * P
    }

    v[0] = A(Offset + k + 1, k);
    v[1] = A(Offset + k + 2, k);
    if (k < N - 3) {
      v[2] = A(Offset + k + 3, k);
    }
  }

  householder_matrix<vect<ValueType, 2>> hhm3(vect<ValueType, 2>(v[0], v[1]),
                                              NumTol);

  mat_sub_block<Matrix1> subA4(A, 2, A.get_col_count(), Offset + N - 2, 0);
  householder_prod(hhm3, subA4);

  mat_sub_block<Matrix2> subB4(B, 2, B.get_col_count(), Offset + N - 2, 0);
  householder_prod(hhm3, subB4);

  if (Q) {
    mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 2, 0, N - 2);
    householder_prod(subQ, hhm3);  // Q * Q_k ^T
  }

  hhm3.set_inv(
      vect<ValueType, 2>(B(Offset + N - 1, N - 2), B(Offset + N - 1, N - 1)),
      NumTol);

  mat_sub_block<Matrix1> subA5(A, A.get_row_count(), 2, 0, N - 2);
  householder_prod(subA5, hhm3);

  mat_sub_block<Matrix2> subB5(B, B.get_row_count(), 2, 0, N - 2);
  householder_prod(subB5, hhm3);

  if (Z) {
    mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, N - 2);
    householder_prod(subZ, hhm3);  // Q_prev * P
  }
}

template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4>
void francis_QZ_step(Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
                     mat_size_type_t<Matrix1> Offset,
                     mat_value_type_t<Matrix1> NumTol) {
  using ValueType = mat_value_type_t<Matrix1>;

  householder_matrix<vect<ValueType, 3>> hhm;

  int N = A.get_row_count() - Offset;
  if (N < 3) {
    return;
  }

  vect<ValueType, 3> v;

  v[0] = ((A(Offset + N - 2, N - 2) / B(Offset + N - 2, N - 2) -
           A(Offset, 0) / B(Offset, 0)) *
              (A(Offset + N - 1, N - 1) / B(Offset + N - 1, N - 1) -
               A(Offset, 0) / B(Offset, 0)) -
          (A(Offset + N - 2, N - 1) / B(Offset + N - 1, N - 1)) *
              (A(Offset + N - 1, N - 2) / B(Offset + N - 2, N - 2)) +
          (A(Offset + N - 1, N - 2) / B(Offset + N - 2, N - 2)) *
              (B(Offset + N - 2, N - 1) / B(Offset + N - 1, N - 1)) *
              (A(Offset, 0) / B(Offset, 0))) *
             (B(Offset, 0) / A(Offset + 1, 0)) +
         A(Offset, 1) / B(Offset + 1, 1) -
         (A(Offset, 0) / B(Offset, 0)) * (B(Offset, 1) / B(Offset + 1, 1));
  v[1] = (A(Offset + 1, 1) / B(Offset + 1, 1) - A(Offset, 0) / B(Offset, 0)) -
         (A(Offset + 1, 0) / B(Offset, 0)) * (B(Offset, 1) / B(Offset + 1, 1)) -
         (A(Offset + N - 2, N - 2) / B(Offset + N - 2, N - 2) -
          A(Offset, 0) / B(Offset, 0)) -
         (A(Offset + N - 1, N - 1) / B(Offset + N - 1, N - 1) -
          A(Offset, 0) / B(Offset, 0)) +
         (A(Offset + N - 1, N - 2) / B(Offset + N - 2, N - 2)) *
             (B(Offset + N - 2, N - 1) / B(Offset + N - 1, N - 1));
  v[2] = A(Offset + 2, 1) / B(Offset + 1, 1);

  francis_QZ_step_impl(A, B, Q, Z, Offset, N, v, NumTol);
}

template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4>
void gen_schur_decomp_impl(Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
                           mat_value_type_t<Matrix1> NumTol = 1E-8) {
  using ValueType = mat_value_type_t<Matrix1>;
  using std::abs;
  using std::sqrt;

  int N = A.get_row_count();
  if (N == 0) {
    return;
  }

  ValueType absNumTol = 0.0;
  for (int i = 0; i < N; ++i) {
    absNumTol += abs(A(i, i));
  }
  absNumTol *= NumTol / N;

  detail::reduce_HessTri_offset_impl(A, B, Q, Z, 0, absNumTol);

  deflate_hess_all_up_impl(A, B, Q, Z, 0, absNumTol);

  int q = N;
  int p = N;

  while (q > 0) {
    bool last_off_diag_was_nil = true;
    int i = N - 1;
    // find a trailing quasi-upper-triangular sub-matrix
    for (; i > 0; --i) {
      if (abs(A(i, i - 1)) < absNumTol) {
        A(i, i - 1) = ValueType(0.0);
        q = i;
        last_off_diag_was_nil = true;
      } else {
        if (!last_off_diag_was_nil) {
          break;
        }
        last_off_diag_was_nil = false;
      }
    }
    if (i == 0) {  // break if it is entirely quasi-upper-triangular.
      break;
    }

    // find the middle, biggest unreduced upper-Hessenberg matrix
    p = 0;  // in case the below loop never gets to the condition.
    for (i = q - 1; i > 0; --i) {
      if (abs(A(i, i - 1)) < absNumTol) {
        A(i, i - 1) = ValueType(0.0);
        p = i;
        break;
      }
    }

    // set remaining sub-diagonals to zero if they are very small.
    for (i = p; i > 0; --i) {
      if (abs(A(i, i - 1)) < absNumTol) {
        A(i, i - 1) = ValueType(0.0);
      }
    }

    mat_sub_block<Matrix1> subA(
        A, q, N - p, 0,
        p);  // the QZ step will affect the entire upper corner block (1:q,p:N)
    mat_sub_block<Matrix2> subB(B, q, N - p, 0, p);

    if (Q) {
      mat_sub_block<Matrix3> subQ(
          *Q, Q->get_row_count(), Q->get_col_count() - p, 0,
          p);  // Q_new will only change in the columns after p.
      if (Z) {
        mat_sub_block<Matrix3> subZ(
            *Z, Z->get_row_count(), Z->get_col_count() - p, 0,
            p);  // Z_new will only change in the columns after p.
        if (!deflate_hess_all_up_impl(subA, subB, &subQ, &subZ, p, absNumTol)) {
          francis_QZ_step(subA, subB, &subQ, &subZ, p, absNumTol);
        }
      } else {
        if (!deflate_hess_all_up_impl(subA, subB, &subQ,
                                      static_cast<Matrix4*>(nullptr), p,
                                      absNumTol)) {
          francis_QZ_step(subA, subB, &subQ, static_cast<Matrix4*>(nullptr), p,
                          absNumTol);
        }
      }
    } else {
      if (Z) {
        mat_sub_block<Matrix3> subZ(*Z, Z->get_row_count(),
                                    Z->get_col_count() - p, 0, p);
        if (!deflate_hess_all_up_impl(subA, subB,
                                      static_cast<Matrix3*>(nullptr), &subZ, p,
                                      absNumTol)) {
          francis_QZ_step(subA, subB, static_cast<Matrix3*>(nullptr), &subZ, p,
                          absNumTol);
        }
      } else {
        if (!deflate_hess_all_up_impl(
                subA, subB, static_cast<Matrix3*>(nullptr),
                static_cast<Matrix4*>(nullptr), p, absNumTol)) {
          francis_QZ_step(subA, subB, static_cast<Matrix3*>(nullptr),
                          static_cast<Matrix4*>(nullptr), p, absNumTol);
        }
      }
    }
  }

  // third step is to deal with the 2-2 blocks remaining and see if they can be reduced
  //  to 2 1-1 blocks (eigenvalues are real).
  // NOTE Not sure this is needed.
  for (q = 1; q < N; ++q) {
    if (abs(A(q, q - 1)) < absNumTol) {
      continue;
    }

    ValueType mu = A(q - 1, q - 1) / B(q - 1, q - 1);
    ValueType a_12 = A(q - 1, q) - mu * B(q - 1, q);
    ValueType a_22 = A(q, q) - mu * B(q, q);
    ValueType p_val =
        ValueType(0.5) * (a_22 / B(q, q) - (B(q - 1, q) * A(q, q - 1)) /
                                               (B(q - 1, q - 1) * B(q, q)));
    ValueType q_val = (A(q, q - 1) * a_12) / (B(q - 1, q - 1) * B(q, q));
    ValueType r = p_val * p_val + q_val;
    if (r < ValueType(0.0)) {
      continue;
    }

    ValueType l = mu + p_val + (p_val > 0 ? sqrt(r) : -sqrt(r));

    givens_rot_matrix<ValueType> G;

    {
      G.set(l * B(q, q) - A(q, q), A(q, q - 1));
      G = transpose(G);

      mat_sub_block<Matrix2> subB(B, q + 1, 2, 0, q - 1);
      givens_rot_prod(subB, G);  // B * G^T

      mat_sub_block<Matrix1> subA(A, q + 1, 2, 0, q - 1);
      givens_rot_prod(subA, G);  // A * G^T

      if (Z) {
        mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, q - 1);
        givens_rot_prod(subZ, G);  // Q_prev * G^T
      }
    }

    {
      G.set(B(q - 1, q - 1), B(q, q - 1));

      mat_sub_block<Matrix1> subA2(A, 2, N - q + 1, q - 1, q - 1);
      givens_rot_prod(G, subA2);  // G * A

      mat_sub_block<Matrix2> subB2(B, 2, N - q + 1, q - 1, q - 1);
      givens_rot_prod(G, subB2);  // G * B

      if (Q) {
        mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 2, 0, q - 1);
        givens_rot_prod(subQ, transpose(G));  // Q_prev * G^T
      }
    }
  }
}

}  // namespace detail

/**
 * Solves for the eigen-values of a matrix, using the symmetric QR algorithm (Golub and vanLoan Alg.-8.3.3).
 *
 * \param A square, symmetric matrix.
 * \param Q holds as output, the unitary square matrix Q.
 * \param D holds as output, the diagonal matrix D in A = Q D Q^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, FullyWritableMatrix Matrix2,
          WritableMatrix Matrix3>
void eigensolve_SymQR(const Matrix1& A, Matrix2& Q, Matrix3& D,
                      mat_value_type_t<Matrix1> NumTol = 1E-8) {
  if (A.get_row_count() != A.get_col_count()) {
    throw std::range_error(
        "Symmetric QR algorithm is only possible on a square (symmetric) "
        "matrix!");
  }

  Q = mat<mat_value_type_t<Matrix2>, mat_structure::identity>(
      A.get_row_count());
  if constexpr (FullyWritableMatrix<Matrix3>) {
    D = A;
    detail::symmetric_QRalg_impl(D, &Q, NumTol);
  } else {
    mat<mat_value_type_t<Matrix1>, mat_structure::square> D_tmp(A);
    detail::symmetric_QRalg_impl(D_tmp, &Q, NumTol);
    D = D_tmp;
  }
}

/**
 * Solves for the eigen-values of a matrix, using the symmetric QR algorithm (Golub and vanLoan Alg.-8.3.3).
 *
 * \param A square, symmetric matrix.
 * \param D holds as output, the diagonal matrix D in A = Q D Q^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, WritableMatrix Matrix2>
void eigensolve_SymQR(const Matrix1& A, Matrix2& D,
                      mat_value_type_t<Matrix1> NumTol = 1E-8) {
  if (A.get_row_count() != A.get_col_count()) {
    throw std::range_error(
        "Symmetric QR algorithm is only possible on a square (symmetric) "
        "matrix!");
  }

  if constexpr (FullyWritableMatrix<Matrix2>) {
    D = A;
    detail::symmetric_QRalg_impl(
        D,
        static_cast<mat<mat_value_type_t<Matrix1>, mat_structure::square>*>(
            nullptr),
        NumTol);
  } else {
    mat<mat_value_type_t<Matrix1>, mat_structure::square> D_tmp(A);
    detail::symmetric_QRalg_impl(
        D_tmp,
        static_cast<mat<mat_value_type_t<Matrix1>, mat_structure::square>*>(
            nullptr),
        NumTol);
    D = D_tmp;
  }
}

/**
 * Inverses a matrix, using the symmetric QR algorithm for eigenvalues (Golub and vanLoan Alg.-8.3.3).
 * Note that this function will output the pseudo-inverse if there are any zero eigenvalues.
 *
 * \param A square, symmetric matrix.
 * \param A_inv holds as output, the (pseudo-)inverse of matrix A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, WritableMatrix Matrix2>
void pseudoinvert_SymQR(const Matrix1& A, Matrix2& A_inv,
                        mat_value_type_t<Matrix1> NumTol = 1E-8) {
  using ValueType = mat_value_type_t<Matrix1>;
  using std::abs;
  if (A.get_row_count() != A.get_col_count()) {
    throw std::range_error(
        "Symmetric QR algorithm is only possible on a square (symmetric) "
        "matrix!");
  }

  mat<ValueType, mat_structure::square> Q(
      mat<ValueType, mat_structure::identity>(A.get_row_count()));
  mat<ValueType, mat_structure::square> D_tmp(A);
  detail::symmetric_QRalg_impl(D_tmp, &Q, NumTol);
  mat<ValueType, mat_structure::diagonal> D(A.get_col_count());
  for (int i = 0; i < A.get_col_count(); ++i) {
    if (abs(D_tmp(i, i)) > NumTol) {
      D(i, i) = ValueType(1.0) / D_tmp(i, i);
    } else {
      D(i, i) = ValueType(0.0);
    }
  }
  D_tmp = Q * D;
  A_inv = D_tmp * transpose_view(Q);
}

/**
 * Inverses a matrix, using the symmetric QR algorithm for eigenvalues (Golub and vanLoan Alg.-8.3.3).
 * Note that this function will output the pseudo-inverse if there are any zero eigenvalues.
 *
 * \param A square, symmetric matrix.
 * \param B holds, as input, the RHS of the linear system of equations,
 *          and, as output, the least-square solution of the system.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, FullyWritableMatrix Matrix2>
void linsolve_SymQR(const Matrix1& A, Matrix2& B,
                    mat_value_type_t<Matrix1> NumTol = 1E-8) {
  using ValueType = mat_value_type_t<Matrix1>;
  using std::abs;
  if (A.get_row_count() != A.get_col_count()) {
    throw std::range_error(
        "Symmetric QR algorithm is only possible on a square (symmetric) "
        "matrix!");
  }
  if (A.get_row_count() != B.get_row_count()) {
    throw std::range_error(
        "The linear system of equations has mismatched dimensions!");
  }

  mat<ValueType, mat_structure::square> Q(
      mat<ValueType, mat_structure::identity>(A.get_row_count()));
  mat<ValueType, mat_structure::square> D_tmp(A);
  detail::symmetric_QRalg_impl(D_tmp, &Q, NumTol);
  B = transpose_view(Q) * B;
  for (int i = 0; i < A.get_col_count(); ++i) {
    ValueType fact(0.0);
    if (abs(D_tmp(i, i)) > NumTol) {
      fact = ValueType(1.0) / D_tmp(i, i);
    }
    for (int j = 0; j < B.get_col_count(); ++j) {
      B(i, j) *= fact;
    }
  }
  B = Q * B;
}

/**
 * Functor to wrap a call to a symmetric QR-algorithm-based linear system solver.
 */
struct SymQR_linsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B,
                  mat_value_type_t<Matrix1> NumTol = 1E-8) {
    X = B;
    linsolve_SymQR(A, X, NumTol);
  }
};

/**
 * Solves for the eigen-values of a matrix, using the symmetric QR algorithm (Golub and vanLoan Alg.-8.3.3).
 *
 * \param A square, symmetric matrix.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1>
mat_value_type<Matrix1> condition_number_SymQR(
    const Matrix1& A, mat_value_type_t<Matrix1> NumTol = 1E-8) {
  using ValueType = mat_value_type_t<Matrix1>;
  using std::abs;
  if (A.get_row_count() != A.get_col_count()) {
    throw std::range_error(
        "Symmetric QR algorithm is only possible on a square (symmetric) "
        "matrix!");
  }

  mat<mat_value_type_t<Matrix1>, mat_structure::square> D_tmp(A);
  detail::symmetric_QRalg_impl(
      D_tmp,
      static_cast<mat<mat_value_type_t<Matrix1>, mat_structure::square>*>(
          nullptr),
      NumTol);
  ValueType l_max(0.0);
  ValueType l_min(std::numeric_limits<ValueType>::infinity());
  for (int i = 0; i < A.get_col_count(); ++i) {
    if (abs(D_tmp(i, i)) > l_max) {
      l_max = abs(D_tmp(i, i));
    }
    if (abs(D_tmp(i, i)) < l_min) {
      l_min = abs(D_tmp(i, i));
    }
  }
  if (l_min < 1e-10 * l_max) {
    return std::numeric_limits<ValueType>::infinity();
  }
  return l_max / l_min;
}

/**
 * Performs the Real Schur decomposition on a matrix, using the Francis QR-step method.
 *
 * \param A square matrix with row-count == column-count.
 * \param Q holds as output, the unitary square matrix Q.
 * \param T holds as output, the quasi-upper-triangular matrix T in A = Q T Q^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, FullyWritableMatrix Matrix2,
          FullyWritableMatrix Matrix3>
void decompose_RealSchur(const Matrix1& A, Matrix2& Q, Matrix3& T,
                         mat_value_type_t<Matrix1> NumTol = 1E-8) {
  if (A.get_row_count() != A.get_col_count()) {
    throw std::range_error(
        "Real Schur decomposition is only possible on a square matrix!");
  }

  Q = mat<mat_value_type_t<Matrix2>, mat_structure::identity>(
      A.get_row_count());
  T = A;
  detail::schur_decomp_impl(T, &Q, NumTol);
}

/**
 * Performs the Real Schur decomposition on a matrix, using the Francis QR-step method.
 *
 * \param A square matrix with row-count == column-count.
 * \param T holds as output, the quasi-upper-triangular matrix T in A = Q T Q^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, FullyWritableMatrix Matrix2>
void decompose_RealSchur(const Matrix1& A, Matrix2& T,
                         mat_value_type_t<Matrix1> NumTol = 1E-8) {
  if (A.get_row_count() != A.get_col_count()) {
    throw std::range_error(
        "Real Schur decomposition is only possible on a square matrix!");
  }

  T = A;
  detail::schur_decomp_impl(T, static_cast<Matrix2*>(nullptr), NumTol);
}

/**
 * Performs the Generalized Real Schur decomposition on a matrix, using the QZ-step method.
 *
 * \param A square matrix with row-count == column-count.
 * \param B square matrix with row-count == column-count.
 * \param Q holds as output, the unitary square matrix Q.
 * \param Z holds as output, the unitary square matrix Z.
 * \param T holds as output, the quasi-upper-triangular matrix T in A = Q T Z^T.
 * \param R holds as output, the upper-triangular matrix R in B = Q R Z^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2,
          FullyWritableMatrix Matrix3, FullyWritableMatrix Matrix4,
          FullyWritableMatrix Matrix5, FullyWritableMatrix Matrix6>
void decompose_GenRealSchur(const Matrix1& A, const Matrix2& B, Matrix3& Q,
                            Matrix4& Z, Matrix5& T, Matrix6& R,
                            mat_value_type_t<Matrix1> NumTol = 1E-8) {
  if ((A.get_row_count() != A.get_col_count()) &&
      (B.get_row_count() != B.get_col_count())) {
    throw std::range_error(
        "Generalized Real Schur decomposition is only possible on a square "
        "matrices!");
  }

  Q = mat<mat_value_type_t<Matrix3>, mat_structure::identity>(
      A.get_row_count());
  Z = mat<mat_value_type_t<Matrix4>, mat_structure::identity>(
      A.get_row_count());
  T = A;
  R = B;
  detail::gen_schur_decomp_impl(T, R, &Q, &Z, NumTol);
}

/**
 * Performs the Generalized Real Schur decomposition on a matrix, using the QZ-step method.
 *
 * \param A square matrix with row-count == column-count.
 * \param B square matrix with row-count == column-count.
 * \param T holds as output, the quasi-upper-triangular matrix T in A = Q T Z^T.
 * \param R holds as output, the upper-triangular matrix R in B = Q R Z^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2,
          FullyWritableMatrix Matrix3, FullyWritableMatrix Matrix4>
void decompose_GenRealSchur(const Matrix1& A, const Matrix2& B, Matrix3& T,
                            Matrix4& R,
                            mat_value_type_t<Matrix1> NumTol = 1E-8) {
  if ((A.get_row_count() != A.get_col_count()) &&
      (B.get_row_count() != B.get_col_count())) {
    throw std::range_error(
        "Generalized Real Schur decomposition is only possible on a square "
        "matrices!");
  }

  T = A;
  R = B;
  detail::gen_schur_decomp_impl(T, R, static_cast<Matrix3*>(nullptr),
                                static_cast<Matrix4*>(nullptr), NumTol);
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_SCHUR_DECOMP_H_

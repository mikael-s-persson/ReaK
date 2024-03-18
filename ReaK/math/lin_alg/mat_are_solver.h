/**
 * \file mat_are_solver.h
 *
 * This library provides function templates to solve Algebraic Riccati Equations (AREs) of
 * different kinds.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_MATH_LIN_ALG_MAT_ARE_SOLVER_H_
#define REAK_MATH_LIN_ALG_MAT_ARE_SOLVER_H_

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/mat_num_exceptions.h"

#include "ReaK/math/lin_alg/mat_hess_decomp.h"
#include "ReaK/math/lin_alg/mat_householder.h"
#include "ReaK/math/lin_alg/mat_schur_decomp.h"

#include "ReaK/math/lin_alg/mat_ctrl_decomp.h"

#include "ReaK/math/lin_alg/mat_balance.h"
#include "ReaK/math/lin_alg/mat_norms.h"

#include <type_traits>

namespace ReaK {

namespace detail {

template <typename Matrix1, typename Matrix2>
mat_value_type_t<Matrix1> get_norm_gen_eigens_impl(const Matrix1& A,
                                                   const Matrix2& B) {
  using ValueType = mat_value_type_t<Matrix1>;
  using SizeType = mat_size_type_t<Matrix1>;
  using std::abs;
  using std::sqrt;
  SizeType N = A.get_row_count();

  if (N == 1) {
    ValueType l = abs(A(0, 0));
    ValueType tmp = abs(B(0, 0));
    if (tmp < std::numeric_limits<ValueType>::epsilon() * l) {
      return std::numeric_limits<ValueType>::infinity();
    }
    return l / tmp;
  }
  ValueType l = abs(A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1));
  ValueType tmp = abs(B(0, 0) * B(1, 1) - B(1, 0) * B(0, 1));
  if (tmp < std::numeric_limits<ValueType>::epsilon() * l) {
    return std::numeric_limits<ValueType>::infinity();
  }
  return sqrt(l / tmp);
}

template <typename Matrix1, typename Matrix2>
mat_value_type_t<Matrix1> get_real_val_gen_eigens_impl(const Matrix1& A,
                                                       const Matrix2& B) {
  using ValueType = mat_value_type_t<Matrix1>;
  using SizeType = mat_size_type_t<Matrix1>;
  using std::abs;
  using std::sqrt;
  SizeType N = A.get_row_count();

  if (N == 1) {
    ValueType l = A(0, 0);
    ValueType tmp = B(0, 0);
    if (abs(tmp) < std::numeric_limits<ValueType>::epsilon() * abs(l)) {
      return (l < 0 ? -std::numeric_limits<ValueType>::infinity()
                    : std::numeric_limits<ValueType>::infinity());
    }
    return l / tmp;
  }
  ValueType l = abs(A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1));
  ValueType tmp = abs(B(0, 0) * B(1, 1) - B(1, 0) * B(0, 1));
  if (tmp < std::numeric_limits<ValueType>::epsilon() * l) {
    return std::numeric_limits<ValueType>::infinity();
  }
  ValueType mu = A(0, 0) / B(0, 0);
  ValueType a_22 = A(1, 1) - mu * B(1, 1);
  ValueType p = ValueType(0.5) *
                (a_22 / B(1, 1) - (B(0, 1) * A(1, 0)) / (B(0, 0) * B(1, 1)));
  return mu + p;
}

struct lesser_norm_eigen_first {

  template <typename Matrix1, typename Matrix2, typename Matrix3,
            typename Matrix4>
  int operator()(const Matrix1& A1, const Matrix2& B1, const Matrix3& A2,
                 const Matrix4& B2) {
    using ValueType = mat_value_type_t<Matrix1>;
    ValueType l1 = get_norm_gen_eigens_impl(A1, B1);
    ValueType l2 = get_norm_gen_eigens_impl(A2, B2);
    if (l1 < l2) {
      return 1;
    }
    if (l1 == l2) {
      return 0;
    }
    return -1;
  }
};

struct greater_norm_eigen_first {

  template <typename Matrix1, typename Matrix2, typename Matrix3,
            typename Matrix4>
  int operator()(const Matrix1& A1, const Matrix2& B1, const Matrix3& A2,
                 const Matrix4& B2) {
    using ValueType = mat_value_type_t<Matrix1>;
    ValueType l1 = get_norm_gen_eigens_impl(A1, B1);
    ValueType l2 = get_norm_gen_eigens_impl(A2, B2);
    if (l1 > l2) {
      return 1;
    }
    if (l1 == l2) {
      return 0;
    }
    return -1;
  }
};

struct lesser_real_val_eigen_first {

  template <typename Matrix1, typename Matrix2, typename Matrix3,
            typename Matrix4>
  int operator()(const Matrix1& A1, const Matrix2& B1, const Matrix3& A2,
                 const Matrix4& B2) {
    using ValueType = mat_value_type_t<Matrix1>;
    ValueType l1 = get_real_val_gen_eigens_impl(A1, B1);
    ValueType l2 = get_real_val_gen_eigens_impl(A2, B2);
    if (l1 < l2) {
      return 1;
    }
    if (l1 == l2) {
      return 0;
    }
    return -1;
  }
};

struct greater_real_val_eigen_first {

  template <typename Matrix1, typename Matrix2, typename Matrix3,
            typename Matrix4>
  int operator()(const Matrix1& A1, const Matrix2& B1, const Matrix3& A2,
                 const Matrix4& B2) {
    using ValueType = mat_value_type_t<Matrix1>;
    ValueType l1 = get_real_val_gen_eigens_impl(A1, B1);
    ValueType l2 = get_real_val_gen_eigens_impl(A2, B2);
    if (l1 > l2) {
      return 1;
    }
    if (l1 == l2) {
      return 0;
    }
    return -1;
  }
};

struct neg_real_val_eigen_first {

  template <typename Matrix1, typename Matrix2, typename Matrix3,
            typename Matrix4>
  int operator()(const Matrix1& A1, const Matrix2& B1, const Matrix3& A2,
                 const Matrix4& B2) {
    using ValueType = mat_value_type_t<Matrix1>;
    ValueType l1 = get_real_val_gen_eigens_impl(A1, B1);
    ValueType l2 = get_real_val_gen_eigens_impl(A2, B2);
    if (l1 < ValueType(0.0) && l2 >= ValueType(0.0)) {
      return 1;
    }
    if (l2 < ValueType(0.0) && l1 >= ValueType(0.0)) {
      return -1;
    }
    return 0;
  }
};

struct pos_real_val_eigen_first {

  template <typename Matrix1, typename Matrix2, typename Matrix3,
            typename Matrix4>
  int operator()(const Matrix1& A1, const Matrix2& B1, const Matrix3& A2,
                 const Matrix4& B2) {
    using ValueType = mat_value_type_t<Matrix1>;
    ValueType l1 = get_real_val_gen_eigens_impl(A1, B1);
    ValueType l2 = get_real_val_gen_eigens_impl(A2, B2);
    if (l1 > ValueType(0.0) && l2 <= ValueType(0.0)) {
      return 1;
    }
    if (l2 > ValueType(0.0) && l1 <= ValueType(0.0)) {
      return -1;
    }
    return 0;
  }
};

struct in_unit_circle_eigen_first {

  template <typename Matrix1, typename Matrix2, typename Matrix3,
            typename Matrix4>
  int operator()(const Matrix1& A1, const Matrix2& B1, const Matrix3& A2,
                 const Matrix4& B2) {
    using ValueType = mat_value_type_t<Matrix1>;
    ValueType l1 = get_norm_gen_eigens_impl(A1, B1);
    ValueType l2 = get_norm_gen_eigens_impl(A2, B2);
    if (l1 < ValueType(1.0) && l2 >= ValueType(1.0)) {
      return 1;
    }
    if (l2 < ValueType(1.0) && l1 >= ValueType(1.0)) {
      return -1;
    }
    return 0;
  }
};

struct out_unit_circle_eigen_first {

  template <typename Matrix1, typename Matrix2, typename Matrix3,
            typename Matrix4>
  int operator()(const Matrix1& A1, const Matrix2& B1, const Matrix3& A2,
                 const Matrix4& B2) {
    using ValueType = mat_value_type_t<Matrix1>;
    ValueType l1 = get_norm_gen_eigens_impl(A1, B1);
    ValueType l2 = get_norm_gen_eigens_impl(A2, B2);
    if (l1 < ValueType(1.0) && l2 >= ValueType(1.0)) {
      return -1;
    }
    if (l2 < ValueType(1.0) && l1 >= ValueType(1.0)) {
      return 1;
    }
    return 0;
  }
};

/// This function performs the swapping of two schur blocks in the real schur pencil (A,B).
/// This is the Case I in Van Dooren (1981), where the blocks to be swapped both have dimension 1.
template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4>
void swap_schur_blocks11_impl(
    Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
    mat_size_type_t<Matrix1> p,  // upper-left diagonal element.
    mat_size_type_t<Matrix1> row_offset) {
  using ValueType = mat_value_type_t<Matrix1>;
  using SizeType = mat_size_type_t<Matrix1>;
  using std::abs;

  givens_rot_matrix<ValueType> G;

  SizeType q = row_offset + p;

  bool reduce_A = abs(B(q + 1, p + 1)) < abs(A(q + 1, p + 1));

  ValueType x1 = A(row_offset + p + 1, p + 1) * B(row_offset + p, p) -
                 B(row_offset + p + 1, p + 1) * A(row_offset + p, p);
  ValueType x2 = A(row_offset + p + 1, p + 1) * B(row_offset + p, p + 1) -
                 B(row_offset + p + 1, p + 1) * A(row_offset + p, p + 1);

  G.set(-x2, x1);
  G = transpose(G);

  mat_sub_block<Matrix2> subB1(B, q + 2, 2, 0, p);
  givens_rot_prod(subB1, G);  // B * G^T

  mat_sub_block<Matrix1> subA1(A, q + 2, 2, 0, p);
  givens_rot_prod(subA1, G);  // A * G^T

  if (Z) {
    mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, p);
    givens_rot_prod(subZ, G);  // Q_prev * G^T
  }

  if (reduce_A) {
    G.set(A(q, p), A(q + 1, p));
  } else {
    G.set(B(q, p), B(q + 1, p));
  }

  mat_sub_block<Matrix1> subA2(A, 2, A.get_col_count() - p, q, p);
  givens_rot_prod(G, subA2);  // G * A

  mat_sub_block<Matrix2> subB2(B, 2, B.get_col_count() - p, q, p);
  givens_rot_prod(G, subB2);  // G * B

  if (Q) {
    mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 2, 0, q - row_offset);
    givens_rot_prod(subQ, transpose(G));  // Q_prev * G^T
  }
}

/// This function performs the swapping of two schur blocks in the real schur pencil (A,B).
/// This is the Case II in Van Dooren (1981), where the two blocks to be swapped both have
/// dimension 2 and 1, in that order.
template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4>
void swap_schur_blocks21_impl(
    Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
    mat_size_type_t<Matrix1> p,  // upper-left diagonal element.
    mat_size_type_t<Matrix1> row_offset) {
  using ValueType = mat_value_type_t<Matrix1>;
  using SizeType = mat_size_type_t<Matrix1>;
  using std::abs;

  givens_rot_matrix<ValueType> G;

  SizeType q = row_offset + p;

  bool reduce_A = abs(B(q + 2, p + 2)) < abs(A(q + 2, p + 2));

  ValueType a33 = A(row_offset + p + 2, p + 2);
  ValueType b33 = B(row_offset + p + 2, p + 2);

  {
    ValueType x11 = a33 * B(row_offset + p, p) - b33 * A(row_offset + p, p);
    ValueType x21 = -b33 * A(row_offset + p + 1, p);

    G.set(x11, x21);

    mat_sub_block<Matrix1> subA(A, 2, A.get_col_count() - p, q, p);
    givens_rot_prod(G, subA);  // G * A

    mat_sub_block<Matrix2> subB(B, 2, B.get_col_count() - p, q, p);
    givens_rot_prod(G, subB);  // G * B

    if (Q) {
      mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 2, 0, q - row_offset);
      givens_rot_prod(subQ, transpose(G));  // Q_prev * G^T
    }
  }

  // Annihilate x1 in R' * H
  {
    ValueType x1 =
        a33 * B(row_offset + p + 1, p + 1) - b33 * A(row_offset + p + 1, p + 1);
    ValueType x2 =
        a33 * B(row_offset + p + 1, p + 2) - b33 * A(row_offset + p + 1, p + 2);

    G.set(-x2, x1);
    G = transpose(G);

    mat_sub_block<Matrix2> subB(B, q + 3, 2, 0, p + 1);
    givens_rot_prod(subB, G);  // B * G^T

    mat_sub_block<Matrix1> subA(A, q + 3, 2, 0, p + 1);
    givens_rot_prod(subA, G);  // A * G^T

    if (Z) {
      mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, p + 1);
      givens_rot_prod(subZ, G);  // Q_prev * G^T
    }
  }

  {
    G.set(B(q + 1, p + 1), B(q + 2, p + 1));

    mat_sub_block<Matrix1> subA(A, 2, A.get_col_count() - p - 1, q + 1, p + 1);
    givens_rot_prod(G, subA);  // G * A

    mat_sub_block<Matrix2> subB(B, 2, B.get_col_count() - p - 1, q + 1, p + 1);
    givens_rot_prod(G, subB);  // G * B

    if (Q) {
      mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 2, 0,
                                  q + 1 - row_offset);
      givens_rot_prod(subQ, transpose(G));  // Q_prev * G^T
    }
  }

  {
    // Annihilate x_2 (is x1 here) in R' * H
    ValueType x1 = a33 * B(row_offset + p, p) - b33 * A(row_offset + p, p);
    ValueType x2 =
        a33 * B(row_offset + p, p + 1) - b33 * A(row_offset + p, p + 1);

    G.set(-x2, x1);
    G = transpose(G);

    mat_sub_block<Matrix2> subB(B, q + 2, 2, 0, p);
    givens_rot_prod(subB, G);  // B * G^T

    mat_sub_block<Matrix1> subA(A, q + 2, 2, 0, p);
    givens_rot_prod(subA, G);  // A * G^T

    if (Z) {
      mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, p);
      givens_rot_prod(subZ, G);  // Q_prev * G^T
    }
  }

  {
    if (reduce_A) {
      G.set(A(q, p), A(q + 1, p));
    } else {
      G.set(B(q, p), B(q + 1, p));
    }

    mat_sub_block<Matrix1> subA(A, 2, A.get_col_count() - p, q, p);
    givens_rot_prod(G, subA);  // G * A

    mat_sub_block<Matrix2> subB(B, 2, B.get_col_count() - p, q, p);
    givens_rot_prod(G, subB);  // G * B

    if (Q) {
      mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 2, 0, q - row_offset);
      givens_rot_prod(subQ, transpose(G));  // Q_prev * G^T
    }
  }
}

/// This function performs the swapping of two schur blocks in the real schur pencil (A,B).
/// This is the Case II in Van Dooren (1981), where the two blocks to be swapped both have
/// dimension 1 and 2, in that order.
template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4>
void swap_schur_blocks12_impl(
    Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
    mat_size_type_t<Matrix1> p,  // upper-left diagonal element.
    mat_size_type_t<Matrix1> row_offset) {
  using ValueType = mat_value_type_t<Matrix1>;
  using SizeType = mat_size_type_t<Matrix1>;
  using std::abs;

  givens_rot_matrix<ValueType> G;

  SizeType q = row_offset + p;

  bool reduce_A = abs(B(q + 2, p + 2)) < abs(A(q + 2, p + 2));

  ValueType a11 = A(q, p);
  ValueType b11 = B(q, p);

  {
    ValueType x33 = a11 * B(q + 2, p + 2) - b11 * A(q + 2, p + 2);
    ValueType x32 = -b11 * A(q + 2, p + 1);

    G.set(-x33, x32);
    G = transpose(G);

    mat_sub_block<Matrix2> subB(B, q + 3, 2, 0, p + 1);
    givens_rot_prod(subB, G);  // B * G^T

    mat_sub_block<Matrix1> subA(A, q + 3, 2, 0, p + 1);
    givens_rot_prod(subA, G);  // A * G^T

    if (Z) {
      mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, p + 1);
      givens_rot_prod(subZ, G);  // Q_prev * G^T
    }
  }

  {
    // Annihilate x1 in R' * H  [ x1; x2 ]
    ValueType x1 = a11 * B(q, p + 1) - b11 * A(q, p + 1);
    ValueType x2 = a11 * B(q + 1, p + 1) - b11 * A(q + 1, p + 1);

    G.set(x1, x2);

    mat_sub_block<Matrix1> subA(A, 2, A.get_col_count() - p, q, p);
    givens_rot_prod(G, subA);  // G * A

    mat_sub_block<Matrix2> subB(B, 2, B.get_col_count() - p, q, p);
    givens_rot_prod(G, subB);  // G * B

    if (Q) {
      mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 2, 0, q - row_offset);
      givens_rot_prod(subQ, transpose(G));  // Q_prev * G^T
    }
  }

  {
    G.set(-B(q + 1, p + 1), B(q + 1, p));
    G = transpose(G);

    mat_sub_block<Matrix2> subB(B, q + 2, 2, 0, p);
    givens_rot_prod(subB, G);  // B * G^T

    mat_sub_block<Matrix1> subA(A, q + 2, 2, 0, p);
    givens_rot_prod(subA, G);  // A * G^T

    if (Z) {
      mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, p);
      givens_rot_prod(subZ, G);  // Q_prev * G^T
    }
  }

  {
    // Annihilate x2 in R' * H  [ x1; x2 ]
    ValueType x1 = a11 * B(q + 1, p + 2) - b11 * A(q + 1, p + 2);
    ValueType x2 = a11 * B(q + 2, p + 2) - b11 * A(q + 2, p + 2);

    G.set(x1, x2);

    mat_sub_block<Matrix1> subA(A, 2, A.get_col_count() - p - 1, q + 1, p + 1);
    givens_rot_prod(G, subA);  // G * A

    mat_sub_block<Matrix2> subB(B, 2, B.get_col_count() - p - 1, q + 1, p + 1);
    givens_rot_prod(G, subB);  // G * B

    if (Q) {
      mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 2, 0,
                                  q + 1 - row_offset);
      givens_rot_prod(subQ, transpose(G));  // Q_prev * G^T
    }
  }

  {
    if (reduce_A) {
      G.set(-A(q + 2, p + 2), A(q + 2, p + 1));
    } else {
      G.set(-B(q + 2, p + 2), B(q + 2, p + 1));
    }
    G = transpose(G);

    mat_sub_block<Matrix2> subB(B, q + 3, 2, 0, p + 1);
    givens_rot_prod(subB, G);  // B * G^T

    mat_sub_block<Matrix1> subA(A, q + 3, 2, 0, p + 1);
    givens_rot_prod(subA, G);  // A * G^T

    if (Z) {
      mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, p + 1);
      givens_rot_prod(subZ, G);  // Q_prev * G^T
    }
  }
}

/// This function performs the swapping of two schur blocks in the real schur pencil (A,B).
/// This is the Case II in Van Dooren (1981), where the two blocks to be swapped both have
/// dimension 2 and 2, in that order.
/// \note This function uses the "direct swapping" method based on Kressner's work.
/// \note The implicit QZ-step method from Van Dooren does not seem to work at all.
template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4>
void swap_schur_blocks22_impl(
    Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
    mat_size_type_t<Matrix1> p,  // upper-left diagonal element.
    mat_size_type_t<Matrix1> row_offset, mat_value_type_t<Matrix1> absNumTol) {
  using ValueType = mat_value_type_t<Matrix1>;
  using SizeType = mat_size_type_t<Matrix1>;
  using std::abs;

  SizeType q = row_offset + p;

  mat<ValueType, mat_structure::rectangular> X(8, 1);
  {
    mat<ValueType, mat_structure::square> E(8, ValueType(0.0));
    sub(E)(range(0, 2), range(0, 2)) = sub(A)(range(q, q + 2), range(p, p + 2));
    sub(E)(range(2, 4), range(2, 4)) = sub(A)(range(q, q + 2), range(p, p + 2));
    sub(E)(range(4, 6), range(0, 2)) = sub(B)(range(q, q + 2), range(p, p + 2));
    sub(E)(range(6, 8), range(2, 4)) = sub(B)(range(q, q + 2), range(p, p + 2));
    E(0, 4) = -A(q + 2, p + 2);
    E(0, 6) = -A(q + 3, p + 2);
    E(1, 5) = -A(q + 2, p + 2);
    E(1, 7) = -A(q + 3, p + 2);
    E(2, 4) = -A(q + 2, p + 3);
    E(2, 6) = -A(q + 3, p + 3);
    E(3, 5) = -A(q + 2, p + 3);
    E(3, 7) = -A(q + 3, p + 3);
    E(4, 4) = -B(q + 2, p + 2);
    E(4, 6) = -B(q + 3, p + 2);
    E(5, 5) = -B(q + 2, p + 2);
    E(5, 7) = -B(q + 3, p + 2);
    E(6, 4) = -B(q + 2, p + 3);
    E(6, 6) = -B(q + 3, p + 3);
    E(7, 5) = -B(q + 2, p + 3);
    E(7, 7) = -B(q + 3, p + 3);

    mat<ValueType, mat_structure::rectangular> F(8, 1);
    F(0, 0) = A(q, p + 2);
    F(1, 0) = A(q + 1, p + 2);
    F(2, 0) = A(q, p + 3);
    F(3, 0) = A(q + 1, p + 3);
    F(4, 0) = B(q, p + 2);
    F(5, 0) = B(q + 1, p + 2);
    F(6, 0) = B(q, p + 3);
    F(7, 0) = B(q + 1, p + 3);

    linlsq_QR_impl(E, X, F, absNumTol);
  }

  mat<ValueType, mat_structure::rectangular> Y_lhs(4, 2);
  Y_lhs(0, 0) = -X(4, 0);
  Y_lhs(1, 0) = -X(5, 0);
  Y_lhs(0, 1) = -X(6, 0);
  Y_lhs(1, 1) = -X(7, 0);
  sub(Y_lhs)(range(2, 4), range(0, 2)) =
      mat<ValueType, mat_structure::identity>(2);
  // sub(Y_lhs)(range(2,4),range(0,2)) *= ValueType(0.1);

  mat<ValueType, mat_structure::square> V =
      mat<ValueType, mat_structure::square>(
          mat<ValueType, mat_structure::identity>(4));
  decompose_QR_impl(Y_lhs, &V, absNumTol);

  mat<ValueType, mat_structure::rectangular> X_lhs(2, 4);
  X_lhs(0, 2) = X(0, 0);
  X_lhs(1, 2) = X(1, 0);
  X_lhs(0, 3) = X(2, 0);
  X_lhs(1, 3) = X(3, 0);
  sub(X_lhs)(range(0, 2), range(0, 2)) =
      mat<ValueType, mat_structure::identity>(2);
  // sub(X_lhs)(range(0,2),range(0,2)) *= ValueType(0.1);

  mat<ValueType, mat_structure::square> W =
      mat<ValueType, mat_structure::square>(
          mat<ValueType, mat_structure::identity>(4));
  decompose_RQ_impl(X_lhs, &W, absNumTol);

  // Now, multiply V^T on (A,B,[Q])
  sub(A)(range(q, q + 4), range(p, A.get_col_count())) =
      transpose_view(V) * sub(A)(range(q, q + 4), range(p, A.get_col_count()));
  sub(B)(range(q, q + 4), range(p, B.get_col_count())) =
      transpose_view(V) * sub(B)(range(q, q + 4), range(p, B.get_col_count()));

  if (Q) {
    sub (*Q)(range(0, Q->get_row_count()),
             range(q - row_offset, q - row_offset + 4)) *= V;
  }

  // Now, multiply W on (A,B,[Z])
  sub(A)(range(0, q + 4), range(p, p + 4)) *= W;
  sub(B)(range(0, q + 4), range(p, p + 4)) *= W;

  if (Z) {
    sub (*Z)(range(0, Z->get_row_count()), range(p, p + 4)) *= W;
  }

  // Finally, reduce B back to triangular form (which may have been lost on the swap.

  givens_rot_matrix<ValueType> G;

  // Z34
  {
    G.set(-B(q + 3, p + 3), B(q + 3, p + 2));
    G = transpose(G);

    mat_sub_block<Matrix2> subB(B, q + 4, 2, 0, p + 2);
    givens_rot_prod(subB, G);  // B * G^T

    mat_sub_block<Matrix1> subA(A, q + 4, 2, 0, p + 2);
    givens_rot_prod(subA, G);  // A * G^T

    if (Z) {
      mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, p + 2);
      givens_rot_prod(subZ, G);  // Q_prev * G^T
    }
  }

  // Z12
  {
    G.set(-B(q + 1, p + 1), B(q + 1, p));
    G = transpose(G);

    mat_sub_block<Matrix2> subB(B, q + 2, 2, 0, p);
    givens_rot_prod(subB, G);  // B * G^T

    mat_sub_block<Matrix1> subA(A, q + 2, 2, 0, p);
    givens_rot_prod(subA, G);  // A * G^T

    if (Z) {
      mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, p);
      givens_rot_prod(subZ, G);  // Q_prev * G^T
    }
  }
}

template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4, typename CompareFunc>
void partition_schur_pencil_impl(Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
                                 CompareFunc compare,
                                 mat_value_type_t<Matrix1> NumTol) {
  using ValueType = mat_value_type_t<Matrix1>;
  using SizeType = mat_size_type_t<Matrix1>;
  using std::abs;
  SizeType N = A.get_row_count();

  ValueType absNumTol = 0.0;
  for (SizeType i = 0; i < N; ++i) {
    absNumTol += abs(A(i, i));
  }
  absNumTol *= NumTol / N;

  // This algorithm is basically an insertion sort based on the compare functor and
  //  a mix of 1-1 and 2-2 blocks to be swapped along the diagonal to sort the eigen-values.
  SizeType q = 0;
  while (q < N - 1) {
    SizeType p = ++q;
    if ((q < N - 1) && (abs(A(q + 1, q)) > absNumTol)) {
      ++q;
    }
    bool is_next_block_by2 = false;
    if ((p < N - 1) && (abs(A(p + 1, p)) > absNumTol)) {
      is_next_block_by2 = true;
    }
    while (true) {
      bool is_prev_block_by2 = false;
      if ((p > 1) && (abs(A(p - 1, p - 2)) > absNumTol)) {
        is_prev_block_by2 = true;
      }
      bool swap_needed = false;
      if (is_next_block_by2 && is_prev_block_by2) {
        swap_needed = (-1 == compare(sub(A)(range(p - 2, p), range(p - 2, p)),
                                     sub(B)(range(p - 2, p), range(p - 2, p)),
                                     sub(A)(range(p, p + 2), range(p, p + 2)),
                                     sub(B)(range(p, p + 2), range(p, p + 2))));
        if (swap_needed) {
          mat_sub_block<Matrix1> subA(A, p + 2, N - p + 2, 0, p - 2);
          mat_sub_block<Matrix2> subB(B, p + 2, N - p + 2, 0, p - 2);
          if (Q) {
            mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 4, 0, p - 2);
            if (Z) {
              mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 4, 0, p - 2);
              swap_schur_blocks22_impl(subA, subB, &subQ, &subZ, 0, p - 2,
                                       absNumTol);
            } else {
              swap_schur_blocks22_impl(subA, subB, &subQ,
                                       static_cast<Matrix4*>(nullptr), 0, p - 2,
                                       absNumTol);
            }
          } else {
            if (Z) {
              mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 4, 0, p - 2);
              swap_schur_blocks22_impl(subA, subB,
                                       static_cast<Matrix3*>(nullptr), &subZ, 0,
                                       p - 2, absNumTol);
            } else {
              swap_schur_blocks22_impl(
                  subA, subB, static_cast<Matrix3*>(nullptr),
                  static_cast<Matrix4*>(nullptr), 0, p - 2, absNumTol);
            }
          }
          p -= 2;
        }
      } else if (is_next_block_by2 && !is_prev_block_by2) {
        swap_needed = (-1 == compare(sub(A)(range(p - 1, p), range(p - 1, p)),
                                     sub(B)(range(p - 1, p), range(p - 1, p)),
                                     sub(A)(range(p, p + 2), range(p, p + 2)),
                                     sub(B)(range(p, p + 2), range(p, p + 2))));
        if (swap_needed) {
          mat_sub_block<Matrix1> subA(A, p + 2, N - p + 1, 0, p - 1);
          mat_sub_block<Matrix2> subB(B, p + 2, N - p + 1, 0, p - 1);
          if (Q) {
            mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 3, 0, p - 1);
            if (Z) {
              mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 3, 0, p - 1);
              swap_schur_blocks12_impl(subA, subB, &subQ, &subZ, 0, p - 1);
            } else {
              swap_schur_blocks12_impl(
                  subA, subB, &subQ, static_cast<Matrix4*>(nullptr), 0, p - 1);
            }
          } else {
            if (Z) {
              mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 3, 0, p - 1);
              swap_schur_blocks12_impl(
                  subA, subB, static_cast<Matrix3*>(nullptr), &subZ, 0, p - 1);
            } else {
              swap_schur_blocks12_impl(
                  subA, subB, static_cast<Matrix3*>(nullptr),
                  static_cast<Matrix4*>(nullptr), 0, p - 1);
            }
          }
          --p;
        }
      } else if (!is_next_block_by2 && is_prev_block_by2) {
        swap_needed = (-1 == compare(sub(A)(range(p - 2, p), range(p - 2, p)),
                                     sub(B)(range(p - 2, p), range(p - 2, p)),
                                     sub(A)(range(p, p + 1), range(p, p + 1)),
                                     sub(B)(range(p, p + 1), range(p, p + 1))));
        if (swap_needed) {
          mat_sub_block<Matrix1> subA(A, p + 1, N - p + 2, 0, p - 2);
          mat_sub_block<Matrix2> subB(B, p + 1, N - p + 2, 0, p - 2);
          if (Q) {
            mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 3, 0, p - 2);
            if (Z) {
              mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 3, 0, p - 2);
              swap_schur_blocks21_impl(subA, subB, &subQ, &subZ, 0, p - 2);
            } else {
              swap_schur_blocks21_impl(
                  subA, subB, &subQ, static_cast<Matrix4*>(nullptr), 0, p - 2);
            }
          } else {
            if (Z) {
              mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 3, 0, p - 2);
              swap_schur_blocks21_impl(
                  subA, subB, static_cast<Matrix3*>(nullptr), &subZ, 0, p - 2);
            } else {
              swap_schur_blocks21_impl(
                  subA, subB, static_cast<Matrix3*>(nullptr),
                  static_cast<Matrix4*>(nullptr), 0, p - 2);
            }
          }
          p -= 2;
        }
      } else if (!is_next_block_by2 && !is_prev_block_by2) {
        swap_needed = (-1 == compare(sub(A)(range(p - 1, p), range(p - 1, p)),
                                     sub(B)(range(p - 1, p), range(p - 1, p)),
                                     sub(A)(range(p, p + 1), range(p, p + 1)),
                                     sub(B)(range(p, p + 1), range(p, p + 1))));
        if (swap_needed) {
          mat_sub_block<Matrix1> subA(A, p + 1, N - p + 1, 0, p - 1);
          mat_sub_block<Matrix2> subB(B, p + 1, N - p + 1, 0, p - 1);
          if (Q) {
            mat_sub_block<Matrix3> subQ(*Q, Q->get_row_count(), 2, 0, p - 1);
            if (Z) {
              mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, p - 1);
              swap_schur_blocks11_impl(subA, subB, &subQ, &subZ, 0, p - 1);
            } else {
              swap_schur_blocks11_impl(
                  subA, subB, &subQ, static_cast<Matrix4*>(nullptr), 0, p - 1);
            }
          } else {
            if (Z) {
              mat_sub_block<Matrix4> subZ(*Z, Z->get_row_count(), 2, 0, p - 1);
              swap_schur_blocks11_impl(
                  subA, subB, static_cast<Matrix3*>(nullptr), &subZ, 0, p - 1);
            } else {
              swap_schur_blocks11_impl(
                  subA, subB, static_cast<Matrix3*>(nullptr),
                  static_cast<Matrix4*>(nullptr), 0, p - 1);
            }
          }
          --p;
        }
      }

      if (!swap_needed || p == 0) {
        break;
      }
    }
  }
}

}  // namespace detail

/// Solves the Continuous-time Algebraic Riccati Equation (for infinite horizon LQR).
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981).
/// This method first reduces the augmented (2n+m x 2n+m) pencil to a (2n x 2n) pencil
/// using a QR decomposition on the last (mxm) block-column (which has infinite eigenvalues).
/// Then, it performs a generalized real Schur decomposition of the pencil. Finally, it
/// reorders the eigenvalues in the pencil such that stable (within unit-circle) eigenvalues
/// percolate to the upper (nxn) pencil, which allows the extraction of the eigenvectors
/// spanning the stable subspace, which are, in turn, used to compute the unique solution P.
/// \n
/// $Q + A^T P + P A - P B R^{-1} B^T P = 0$
/// \n
/// The initial pencil is: lambda * (I 0 0; 0 I 0; 0 0 0) - (A 0 B; -Q -A^T 0; 0 B^T R)
/// \n
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A fully-writable (square) matrix type.
/// \param A square (n x n) matrix which represents state-to-state-derivative linear map.
/// \param B rectangular (n x m) matrix which represents input-to-state-derivative linear map.
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-semi-definite matrix which represents quadratic input penalty.
/// \param P holds as output, the nonnegative definite solution to Q + A^T P + P A - P B R^-1 B^T P = 0.
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the CARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableMatrix Matrix3,
          ReadableMatrix Matrix4, FullyWritableMatrix Matrix5>
void solve_care_problem(const Matrix1& A, const Matrix2& B, const Matrix3& Q,
                        const Matrix4& R, Matrix5& P,
                        mat_value_type_t<Matrix1> NumTol = 1E-8,
                        bool UseBalancing = false) {
  if ((A.get_row_count() != A.get_col_count()) ||
      (B.get_row_count() != A.get_row_count()) ||
      (Q.get_row_count() != Q.get_col_count()) ||
      (R.get_row_count() != R.get_col_count()) ||
      (B.get_col_count() != R.get_col_count())) {
    throw std::range_error(
        "The dimensions of the CARE system matrices do not match! Should be "
        "A(n x n), B(n x m), "
        "Q(n x n), and R(m x m).");
  }

  using ValueType = mat_value_type_t<Matrix1>;
  using SizeType = mat_size_type_t<Matrix1>;
  SizeType N = A.get_row_count();
  SizeType M = R.get_row_count();
  if ((N == 0) || (M == 0)) {
    return;
  }

  mat<ValueType, mat_structure::rectangular> R_tmp(N * 2 + M, M);
  sub(R_tmp)(range(0, M), range(0, M)) = R;
  sub(R_tmp)(range(M, M + N), range(0, M)) = B;
  mat<ValueType, mat_structure::square> Q_tmp =
      mat<ValueType, mat_structure::square>(N * 2 + M);
  sub(Q_tmp)(range(0, 2 * N), range(M, M + 2 * N)) =
      mat<ValueType, mat_structure::identity>(2 * N);
  sub(Q_tmp)(range(2 * N, 2 * N + M), range(0, M)) =
      mat<ValueType, mat_structure::identity>(M);

  detail::decompose_QR_impl(R_tmp, &Q_tmp, NumTol);
  Q_tmp = transpose(Q_tmp);

  mat<ValueType, mat_structure::rectangular> B_aug(2 * N, 2 * N);
  B_aug = sub(Q_tmp)(range(M, M + 2 * N), range(0, 2 * N));

  mat<ValueType, mat_structure::rectangular> A_aug(2 * N, 2 * N);
  sub(A_aug)(range(0, 2 * N), range(0, N)) =
      sub(Q_tmp)(range(M, M + 2 * N), range(0, N)) * A -
      sub(Q_tmp)(range(M, M + 2 * N), range(N, 2 * N)) * Q;
  sub(A_aug)(range(0, 2 * N), range(N, 2 * N)) =
      sub(Q_tmp)(range(M, M + 2 * N), range(2 * N, 2 * N + M)) *
          transpose_view(B) -
      sub(Q_tmp)(range(M, M + 2 * N), range(N, 2 * N)) * transpose_view(A);

  mat<ValueType, mat_structure::square> Q_aug =
      mat<ValueType, mat_structure::square>(
          mat<ValueType, mat_structure::identity>(2 * N));
  mat<ValueType, mat_structure::square> Z_aug =
      mat<ValueType, mat_structure::square>(
          mat<ValueType, mat_structure::identity>(2 * N));

  bool should_interchange = false;
  if (norm_1(A_aug) > norm_1(B_aug)) {
    should_interchange = true;
  }

  vect_n<int> Dl_aug(2 * N);
  vect_n<int> Dr_aug(2 * N);
  if (UseBalancing) {
    balance_pencil(A_aug, B_aug, Dl_aug, Dr_aug);
  }

  if (should_interchange) {
    detail::gen_schur_decomp_impl(B_aug, A_aug, &Q_aug, &Z_aug, NumTol);
  } else {
    detail::gen_schur_decomp_impl(A_aug, B_aug, &Q_aug, &Z_aug, NumTol);
  }

  if (should_interchange) {
    detail::partition_schur_pencil_impl(B_aug, A_aug, &Q_aug, &Z_aug,
                                        detail::neg_real_val_eigen_first(),
                                        NumTol);
  } else {
    detail::partition_schur_pencil_impl(A_aug, B_aug, &Q_aug, &Z_aug,
                                        detail::neg_real_val_eigen_first(),
                                        NumTol);
  }

  P.set_row_count(N);
  P.set_col_count(N);
  mat_sub_block<mat<ValueType, mat_structure::square>> subZ11(Z_aug, N, N, 0,
                                                              0);
  mat_sub_block<mat<ValueType, mat_structure::square>> subZ21(Z_aug, N, N, N,
                                                              0);
  try {
    linlsq_QR(transpose_view(subZ11), P, transpose_view(subZ21), NumTol);
  } catch (singularity_error& e) {
    RK_UNUSED(e);
    throw singularity_error(
        "The Continuous-time Algebraic Riccati Equation (CARE) cannot be "
        "solved! Usually "
        "indicates that the system is not stabilizable.");
  }

  if (UseBalancing) {
    apply_left_bal_inv_exp(sub(Dr_aug)[range(0, N)], P);
    apply_right_bal_exp(P, sub(Dr_aug)[range(N, 2 * N)]);
  }

  P += transpose(P);
  P *= ValueType(0.5);
}

/// Solves the Infinite-horizon Continuous-time Linear Quadratic Regulator (LQR) problem.
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981)
/// and implemented in the function solve_care_problem.
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A fully-writable matrix type.
/// \tparam Matrix6 A fully-writable (square) matrix type.
/// \param A square (n x n) matrix which represents state-to-state-derivative linear map.
/// \param B rectangular (n x m) matrix which represents input-to-state-derivative linear map.
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-definite matrix which represents quadratic input penalty.
/// \param K holds as output, the (mxn) LQR-optimal gain matrix for u = - K * (x_cur - x_ref).
/// \param P holds as output, the (nxn) nonnegative definite solution to P = F^T P F - F^T P G ( R + G^T P G )^{-1} G^T P
/// F + Q.
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the CARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableMatrix Matrix3,
          ReadableMatrix Matrix4, FullyWritableMatrix Matrix5, FullyWritableMatrix Matrix6>
void solve_IHCT_LQR(const Matrix1& A, const Matrix2& B, const Matrix3& Q,
                    const Matrix4& R, Matrix5& K, Matrix6& P,
                    mat_value_type_t<Matrix1> NumTol = 1E-8,
                    bool UseBalancing = false) {
  using ValueType = mat_value_type_t<Matrix1>;
  solve_care_problem(A, B, Q, R, P, NumTol, UseBalancing);

  mat<ValueType, mat_structure::rectangular> M_tmp(B.get_col_count(),
                                                   A.get_col_count());
  M_tmp = transpose_view(B) * P;
  linlsq_QR(R, K, M_tmp);
}

/// Solves the Infinite-horizon Continuous-time Linear Quadratic Regulator (LQR) problem.
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981)
/// and implemented in the function solve_care_problem.
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A fully-writable matrix type.
/// \param A square (n x n) matrix which represents state-to-state-derivative linear map.
/// \param B rectangular (n x m) matrix which represents input-to-state-derivative linear map.
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-definite matrix which represents quadratic input penalty.
/// \param K holds as output, the (mxn) LQR-optimal gain matrix for u = - K * (x_cur - x_ref).
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the CARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableMatrix Matrix3,
          ReadableMatrix Matrix4, WritableMatrix Matrix5>
void solve_IHCT_LQR(const Matrix1& A, const Matrix2& B, const Matrix3& Q,
                    const Matrix4& R, Matrix5& K,
                    mat_value_type_t<Matrix1> NumTol = 1E-8,
                    bool UseBalancing = false) {
  using ValueType = mat_value_type_t<Matrix1>;
  mat<ValueType, mat_structure::square> P =
      mat<ValueType, mat_structure::square>(A.get_row_count());
  solve_IHCT_LQR(A, B, Q, R, K, P, NumTol, UseBalancing);
}

/// Solves the Infinite-horizon Continuous-time Linear Quadratic Regulator (LQR) problem
/// with a controllability reduction prior to solving the Algebraic Riccati Equation (CARE).
/// This function should be used if the system is (or could be) uncontrollable in some
/// directions of motion (modes).
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981)
/// and implemented in the function solve_care_problem.
///
/// \note If the system is not completely controllable, the resulting controller will in effect ignore
/// the uncontrollable state errors, both for the computation of the control inputs and for the
/// cost-to-go metric (quadratic weighting matrix P). This means that both matrices K and P will
/// have a rank equal to that returned by this function (number of controllable states). This
/// implies that the possible rank-deficiency of the matrices must be considered when doing further
/// operations with these matrices (e.g., inversions or factorizations).
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A fully-writable matrix type.
/// \tparam Matrix6 A fully-writable (square) matrix type.
/// \tparam Matrix7 A fully-writable (square) matrix type.
/// \param A square (n x n) matrix which represents state-to-state-derivative linear map.
/// \param B rectangular (n x m) matrix which represents input-to-state-derivative linear map.
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-definite matrix which represents quadratic input penalty.
/// \param K holds as output, the (mxn) LQR-optimal gain matrix for u = - K * (x_cur - x_ref).
/// \param P holds as output, the (nxn) nonnegative definite solution to P = F^T P F - F^T P G ( R + G^T P G )^{-1} G^T P
///          F + Q.
/// \param Qr holds as output, the (nxn) orthogonal transformation that maps the state vector into an
///           equivalent state vector whose first r elements are controllable and N-r elements remaining are not.
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
/// \return The numerical rank of the system, i.e., the number of controllable states.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the CARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableMatrix Matrix3,
          ReadableMatrix Matrix4, FullyWritableMatrix Matrix5, FullyWritableMatrix Matrix6,
          FullyWritableMatrix Matrix7>
std::size_t solve_IHCT_LQR_with_reduction(
    const Matrix1& A, const Matrix2& B, const Matrix3& Q, const Matrix4& R,
    Matrix5& K, Matrix6& P, Matrix7& Qr,
    mat_value_type_t<Matrix1> NumTol = 1E-8, bool UseBalancing = false) {
  if ((A.get_row_count() != A.get_col_count()) ||
      (B.get_row_count() != A.get_row_count()) ||
      (Q.get_row_count() != Q.get_col_count()) ||
      (R.get_row_count() != R.get_col_count()) ||
      (B.get_col_count() != R.get_col_count())) {
    throw std::range_error(
        "Infinite-horizon Continuous-time LQR with reduction: The dimensions "
        "of the system "
        "matrices do not match! Should be A(n x n), B(n x m), Q(n x n), and "
        "R(m x m).");
  }

  using ValueType = mat_value_type_t<Matrix1>;
  using SizeType = mat_size_type_t<Matrix1>;

  SizeType N = A.get_row_count();
  SizeType M = B.get_col_count();

  mat<ValueType, mat_structure::rectangular> Br(B);
  mat<ValueType, mat_structure::rectangular> Ar(A);
  Qr = mat<ValueType, mat_structure::identity>(N);
  mat<ValueType, mat_structure::square> Zr(
      (mat<ValueType, mat_structure::identity>(M)));

  SizeType r = ctrl_reduction(Ar, Br, Qr, Zr, NumTol);

  // all states are controllable. Use the normal function (avoid unnecessary operations):
  if (r == N) {
    solve_IHCT_LQR(A, B, Q, R, K, P, NumTol, UseBalancing);
    return r;
  }

  // create the transformed quadratic penalty matrices:
  mat<ValueType, mat_structure::rectangular> Qrr(
      sub(Qr)(range(0, N), range(0, r)));
  mat<ValueType, mat_structure::square> R_reduced(transpose_view(Zr) * R * Zr);
  mat<ValueType, mat_structure::square> Q_reduced(transpose_view(Qrr) * Q *
                                                  Qrr);

  // solve the IHCT LQR problem over the controllable states:
  mat<ValueType, mat_structure::rectangular> P_reduced(r, r);
  mat<ValueType, mat_structure::rectangular> K_reduced(M, r);
  solve_IHCT_LQR(sub(Ar)(range(0, r), range(0, r)),
                 sub(Br)(range(0, r), range(0, M)), Q_reduced, R_reduced,
                 K_reduced, P_reduced, NumTol, UseBalancing);

  // compute the resulting cost-matrix and gain-matrix (note: they are positive-semi-definite and rank deficient,
  // respectively).
  P = Qrr * P_reduced * transpose_view(Qrr);
  K = Zr * K_reduced * transpose_view(Qrr);

  return r;
}

/// Solves the Infinite-horizon Continuous-time Linear Quadratic Regulator (LQR) problem
/// with a controllability reduction prior to solving the Algebraic Riccati Equation (CARE).
/// This function should be used if the system is (or could be) uncontrollable in some
/// directions of motion (modes).
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981)
/// and implemented in the function solve_care_problem.
///
/// \note If the system is not completely controllable, the resulting controller will in effect ignore
/// the uncontrollable state errors, both for the computation of the control inputs and for the
/// cost-to-go metric (quadratic weighting matrix P). This means that both matrices K and P will
/// have a rank equal to that returned by this function (number of controllable states). This
/// implies that the possible rank-deficiency of the matrices must be considered when doing further
/// operations with these matrices (e.g., inversions or factorizations).
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A fully-writable matrix type.
/// \tparam Matrix6 A fully-writable (square) matrix type.
/// \param A square (n x n) matrix which represents state-to-state-derivative linear map.
/// \param B rectangular (n x m) matrix which represents input-to-state-derivative linear map.
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-definite matrix which represents quadratic input penalty.
/// \param K holds as output, the (mxn) LQR-optimal gain matrix for u = - K * (x_cur - x_ref).
/// \param P holds as output, the (nxn) nonnegative definite solution to P = F^T P F - F^T P G ( R + G^T P G )^{-1} G^T P
///F + Q.
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
/// \return The numerical rank of the system, i.e., the number of controllable states.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the CARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableMatrix Matrix3,
          ReadableMatrix Matrix4, WritableMatrix Matrix5, WritableMatrix Matrix6>
std::size_t solve_IHCT_LQR_with_reduction(
    const Matrix1& A, const Matrix2& B, const Matrix3& Q, const Matrix4& R,
    Matrix5& K, Matrix6& P, mat_value_type_t<Matrix1> NumTol = 1E-8,
    bool UseBalancing = false) {
  using ValueType = mat_value_type_t<Matrix1>;
  mat<ValueType, mat_structure::square> Qr =
      mat<ValueType, mat_structure::square>(A.get_row_count());
  return solve_IHCT_LQR_with_reduction(A, B, Q, R, K, P, Qr, NumTol,
                                       UseBalancing);
}

/// Solves the Infinite-horizon Continuous-time Linear Quadratic Regulator (LQR) problem
/// with a controllability reduction prior to solving the Algebraic Riccati Equation (CARE).
/// This function should be used if the system is (or could be) uncontrollable in some
/// directions of motion (modes).
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981)
/// and implemented in the function solve_care_problem.
///
/// \note If the system is not completely controllable, the resulting controller will in effect ignore
/// the uncontrollable state errors, both for the computation of the control inputs and for the
/// cost-to-go metric (quadratic weighting matrix P). This means that both matrices K and P will
/// have a rank equal to that returned by this function (number of controllable states). This
/// implies that the possible rank-deficiency of the matrices must be considered when doing further
/// operations with these matrices (e.g., inversions or factorizations).
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A fully-writable matrix type.
/// \param A square (n x n) matrix which represents state-to-state-derivative linear map.
/// \param B rectangular (n x m) matrix which represents input-to-state-derivative linear map.
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-definite matrix which represents quadratic input penalty.
/// \param K holds as output, the (mxn) LQR-optimal gain matrix for u = - K * (x_cur - x_ref).
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
/// \return The numerical rank of the system, i.e., the number of controllable states.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the CARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableMatrix Matrix3,
          ReadableMatrix Matrix4, WritableMatrix Matrix5>
std::size_t solve_IHCT_LQR_with_reduction(
    const Matrix1& A, const Matrix2& B, const Matrix3& Q, const Matrix4& R,
    Matrix5& K, mat_value_type_t<Matrix1> NumTol = 1E-8,
    bool UseBalancing = false) {
  using ValueType = mat_value_type_t<Matrix1>;
  mat<ValueType, mat_structure::square> P =
      mat<ValueType, mat_structure::square>(A.get_row_count());
  return solve_IHCT_LQR_with_reduction(A, B, Q, R, K, P, NumTol, UseBalancing);
}

/// Solves the Infinite-horizon Continuous-time Affine Quadratic Regulator (AQR) problem.
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981)
/// and implemented in the function solve_care_problem.
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A fully-writable matrix type.
/// \tparam Matrix6 A fully-writable (square) matrix type.
/// \tparam Vector1 A readable vector type.
/// \tparam Vector2 A writable vector type.
/// \param A square (n x n) matrix which represents state-to-state-derivative linear map.
/// \param B rectangular (n x m) matrix which represents input-to-state-derivative linear map.
/// \param c a readable vector (n) which represents the constant term of the state-derivative expression.
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-definite matrix which represents quadratic input penalty.
/// \param K holds as output, the (mxn) LQR-optimal gain matrix for u = - K * (x_cur - x_ref).
/// \param P holds as output, the (nxn) nonnegative definite solution to P = F^T P F - F^T P G ( R + G^T P G )^{-1} G^T P
///          F + Q.
/// \param u_bias holds as output, the constant input term to apply to the system in addition to the feedback term (K *
///               (x_cur - x_ref)).
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the CARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableVector Vector1,
          ReadableMatrix Matrix3, ReadableMatrix Matrix4, FullyWritableMatrix Matrix5,
          FullyWritableMatrix Matrix6, WritableVector Vector2>
void solve_IHCT_AQR(const Matrix1& A, const Matrix2& B, const Vector1& c,
                    const Matrix3& Q, const Matrix4& R, Matrix5& K, Matrix6& P,
                    Vector2& u_bias, mat_value_type_t<Matrix1> NumTol = 1E-8,
                    bool UseBalancing = false) {
  using ValueType = mat_value_type_t<Matrix1>;

  solve_care_problem(A, B, Q, R, P, NumTol, UseBalancing);

  mat<ValueType, mat_structure::rectangular> M_tmp(B.get_col_count(),
                                                   A.get_col_count());
  M_tmp = transpose_view(B) * P;
  linlsq_QR(R, K, M_tmp);

  mat<ValueType, mat_structure::square> KB_A(
      transpose_view(K) * transpose_view(B) - transpose_view(A));
  mat_const_vect_adaptor<Vector1> c_v_m(c);
  mat<ValueType, mat_structure::rectangular> eta(c_v_m);
  linlsq_QR(KB_A, eta, c_v_m);
  mat_vect_adaptor<Vector2> u_bias_m(u_bias);
  linlsq_QR(R, u_bias_m, transpose_view(B) * eta);
}

/// Solves the Infinite-horizon Continuous-time Affine Quadratic Regulator (AQR) problem.
/// This function should be used if the system is (or could be) uncontrollable in some
/// directions of motion (modes).
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981)
/// and implemented in the function solve_care_problem.
///
/// \note If the system is not completely controllable, the resulting controller will in effect ignore
/// the uncontrollable state errors, both for the computation of the control inputs and for the
/// cost-to-go metric (quadratic weighting matrix P). This means that both matrices K and P will
/// have a rank equal to that returned by this function (number of controllable states). This
/// implies that the possible rank-deficiency of the matrices must be considered when doing further
/// operations with these matrices (e.g., inversions or factorizations).
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A fully-writable matrix type.
/// \tparam Matrix6 A fully-writable (square) matrix type.
/// \tparam Vector1 A readable vector type.
/// \tparam Vector2 A writable vector type.
/// \param A square (n x n) matrix which represents state-to-state-derivative linear map.
/// \param B rectangular (n x m) matrix which represents input-to-state-derivative linear map.
/// \param c a readable vector (n) which represents the constant term of the state-derivative expression.
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-definite matrix which represents quadratic input penalty.
/// \param K holds as output, the (mxn) LQR-optimal gain matrix for u = - K * (x_cur - x_ref).
/// \param P holds as output, the (nxn) nonnegative definite solution to P = F^T P F - F^T P G ( R + G^T P G )^{-1} G^T P
///          F + Q.
/// \param u_bias holds as output, the constant input term to apply to the system in addition to the feedback term (K *
///               (x_cur - x_ref)).
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the CARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableVector Vector1,
          ReadableMatrix Matrix3, ReadableMatrix Matrix4, FullyWritableMatrix Matrix5,
          FullyWritableMatrix Matrix6, WritableVector Vector2>
std::size_t solve_IHCT_AQR_with_reduction(
    const Matrix1& A, const Matrix2& B, const Vector1& c, const Matrix3& Q,
    const Matrix4& R, Matrix5& K, Matrix6& P, Vector2& u_bias,
    mat_value_type_t<Matrix1> NumTol = 1E-8, bool UseBalancing = false) {
  if ((A.get_row_count() != A.get_col_count()) ||
      (B.get_row_count() != A.get_row_count()) ||
      (Q.get_row_count() != Q.get_col_count()) ||
      (R.get_row_count() != R.get_col_count()) ||
      (B.get_col_count() != R.get_col_count())) {
    throw std::range_error(
        "Infinite-horizon Continuous-time AQR with reduction: The dimensions "
        "of the system "
        "matrices do not match! Should be A(n x n), B(n x m), Q(n x n), and "
        "R(m x m).");
  }

  using ValueType = mat_value_type_t<Matrix1>;
  using SizeType = mat_size_type_t<Matrix1>;

  SizeType N = A.get_row_count();
  SizeType M = B.get_col_count();
  mat<ValueType, mat_structure::rectangular> Br(B);
  mat<ValueType, mat_structure::rectangular> Ar(A);
  mat<ValueType, mat_structure::square> Qr(
      (mat<ValueType, mat_structure::identity>(N)));
  mat<ValueType, mat_structure::square> Zr(
      (mat<ValueType, mat_structure::identity>(M)));

  SizeType r = ctrl_reduction(Ar, Br, Qr, Zr, NumTol);

  // all states are controllable. Use the normal function (avoid unnecessary operations):
  if (r == N) {
    solve_IHCT_AQR(A, B, c, Q, R, K, P, u_bias, NumTol, UseBalancing);
    return r;
  }

  // create the transformed quadratic penalty matrices:
  mat<ValueType, mat_structure::rectangular> Qrr(
      sub(Qr)(range(0, N), range(0, r)));
  mat<ValueType, mat_structure::square> R_reduced(transpose_view(Zr) * R * Zr);
  mat<ValueType, mat_structure::square> Q_reduced(transpose_view(Qrr) * Q *
                                                  Qrr);

  // solve the IHCT AQR problem over the controllable states:
  mat<ValueType, mat_structure::rectangular> P_reduced(r, r);
  mat<ValueType, mat_structure::rectangular> K_reduced(M, r);
  solve_IHCT_AQR(sub(Ar)(range(0, r), range(0, r)),
                 sub(Br)(range(0, r), range(0, M)), transpose_view(Qrr) * c,
                 Q_reduced, R_reduced, K_reduced, P_reduced, u_bias, NumTol,
                 UseBalancing);

  // compute the resulting cost-matrix and gain-matrix (note: they are positive-semi-definite and rank deficient,
  // respectively).
  P = Qrr * P_reduced * transpose_view(Qrr);
  K = Zr * K_reduced * transpose_view(Qrr);
  u_bias = Zr * u_bias;

  return r;
}

/// Solves the Infinite-horizon Continuous-time Linear Quadratic Gaussian control (LQG) problem.
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981)
/// and implemented in the function solve_dare_problem.
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A readable matrix type.
/// \tparam Matrix7 A readable matrix type.
/// \tparam Matrix8 A fully-writable matrix type.
/// \tparam Matrix9 A fully-writable (square) matrix type.
/// \tparam Matrix10 A fully-writable matrix type.
/// \tparam Matrix11 A fully-writable (square) matrix type.
/// \param A square (n x n) matrix which represents state-to-state-derivative linear map.
/// \param B rectangular (n x m) matrix which represents input-to-state-derivative linear map.
/// \param C rectangular (l x n) matrix which represents state-to-output linear map.
/// \param V square (n x n) positive-semi-definite matrix which represents the covariance of the state disturbances.
/// \param W square (l x l) positive-definite matrix which represents the covariance of the measurement noise (additive).
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-definite matrix which represents quadratic input penalty.
/// \param K holds as output, the (nxl) Kalman gain matrix for x_posterior = x_prior + K * (y - x_prior).
/// \param P holds as output, the (nxn) nonnegative definite solution to P = F P F^T - F P H^T ( W + H P H^T )^{-1} H P
///          F^T + V.
/// \param L holds as output, the (mxn) LQR-optimal gain matrix for u = - K * (x_cur - x_ref).
/// \param S holds as output, the (nxn) nonnegative definite solution to S = F^T S F - F^T S G ( R + G^T S G )^{-1} G^T S
///          F + Q.
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the CARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableMatrix Matrix3,
          ReadableMatrix Matrix4, ReadableMatrix Matrix5, ReadableMatrix Matrix6,
          ReadableMatrix Matrix7, FullyWritableMatrix Matrix8, FullyWritableMatrix Matrix9,
          FullyWritableMatrix Matrix10, FullyWritableMatrix Matrix11>
void solve_IHCT_LQG(const Matrix1& A, const Matrix2& B, const Matrix3& C,
                    const Matrix4& V, const Matrix5& W, const Matrix6& Q,
                    const Matrix7& R, Matrix8& K, Matrix9& P, Matrix10& L,
                    Matrix11& S, mat_value_type_t<Matrix1> NumTol = 1E-8,
                    bool UseBalancing = false) {
  using ValueType = mat_value_type_t<Matrix1>;

  solve_care_problem(transpose_view(A), transpose_view(C), V, W, P, NumTol);

  mat<ValueType, mat_structure::rectangular> M1_tmp(C.get_row_count(),
                                                    A.get_col_count());
  M1_tmp = B * P;
  mat<ValueType, mat_structure::rectangular> Msol_tmp;
  linlsq_QR(W, Msol_tmp, M1_tmp);
  K = transpose_view(Msol_tmp);

  solve_care_problem(A, B, Q, R, S, NumTol, UseBalancing);

  mat<ValueType, mat_structure::rectangular> M2_tmp(B.get_col_count(),
                                                    A.get_col_count());
  M2_tmp = transpose_view(B) * S;
  linlsq_QR(R, L, M2_tmp);
}

/// Solves the Infinite-horizon Continuous-time Linear Quadratic Gaussian control (LQG) problem.
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981)
/// and implemented in the function solve_dare_problem.
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A readable matrix type.
/// \tparam Matrix7 A readable matrix type.
/// \tparam Matrix8 A fully-writable matrix type.
/// \tparam Matrix9 A fully-writable matrix type.
/// \param A square (n x n) matrix which represents state-to-state-derivative linear map.
/// \param B rectangular (n x m) matrix which represents input-to-state-derivative linear map.
/// \param C rectangular (l x n) matrix which represents state-to-output linear map.
/// \param V square (n x n) positive-semi-definite matrix which represents the covariance of the state disturbances.
/// \param W square (l x l) positive-definite matrix which represents the covariance of the measurement noise (additive).
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-definite matrix which represents quadratic input penalty.
/// \param K holds as output, the (nxl) Kalman gain matrix for x_posterior = x_prior + K * (y - x_prior).
/// \param L holds as output, the (mxn) LQR-optimal gain matrix for u = - K * (x_cur - x_ref).
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the CARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4, typename Matrix5, typename Matrix6,
          typename Matrix7, typename Matrix8, typename Matrix9>
void solve_IHCT_LQG(const Matrix1& A, const Matrix2& B, const Matrix3& C,
                    const Matrix4& V, const Matrix5& W, const Matrix6& Q,
                    const Matrix7& R, Matrix8& K, Matrix9& L,
                    mat_value_type_t<Matrix1> NumTol = 1E-8,
                    bool UseBalancing = false) {
  using ValueType = mat_value_type_t<Matrix1>;
  mat<ValueType, mat_structure::square> P =
      mat<ValueType, mat_structure::square>(A.get_row_count());
  mat<ValueType, mat_structure::square> S =
      mat<ValueType, mat_structure::square>(A.get_row_count());
  solve_IHCT_LQG(A, B, C, V, W, Q, R, K, P, L, S, NumTol, UseBalancing);
}

/// Solves the Discrete-time Algebraic Riccati Equation (for infinite horizon LQR).
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981).
/// This method first reduces the augmented (2n+m x 2n+m) pencil to a (2n x 2n) pencil
/// using a QR decomposition on the last (mxm) block-column (which has infinite eigenvalues).
/// Then, it performs a generalized real Schur decomposition of the pencil. Finally, it
/// reorders the eigenvalues in the pencil such that stable (within unit-circle) eigenvalues
/// percolate to the upper (nxn) pencil, which allows the extraction of the eigenvectors
/// spanning the stable subspace, which are, in turn, used to compute the unique solution P.
/// \n
/// $P = F^T P F - F^T P G ( R + G^T P G )^{-1} G^T P F + Q$
/// \n
/// The initial pencil is: lambda * (I 0 0; 0 F^T 0; 0 G^T 0) - (F 0 -G; -Q I 0; 0 0 R)
/// \n
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A fully-writable (square) matrix type.
/// \param F square (n x n) matrix which represents state-to-next-state linear map.
/// \param G rectangular (n x m) matrix which represents input-to-next-state linear map.
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-semi-definite matrix which represents quadratic input penalty.
/// \param P holds as output, the nonnegative definite solution to P = F^T P F - F^T P G ( R + G^T P G )^{-1} G^T P F + Q.
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the DARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableMatrix Matrix3,
          ReadableMatrix Matrix4, FullyWritableMatrix Matrix5>
void solve_dare_problem(const Matrix1& F, const Matrix2& G, const Matrix3& Q,
                        const Matrix4& R, Matrix5& P,
                        mat_value_type_t<Matrix1> NumTol = 1E-8,
                        bool UseBalancing = false) {
  if ((F.get_row_count() != F.get_col_count()) ||
      (G.get_row_count() != F.get_row_count()) ||
      (Q.get_row_count() != Q.get_col_count()) ||
      (R.get_row_count() != R.get_col_count()) ||
      (G.get_col_count() != R.get_col_count())) {
    throw std::range_error(
        "The dimensions of the DARE system matrices do not match! Should be "
        "F(n x n), G(n x m), "
        "Q(n x n), and R(m x m).");
  }

  using ValueType = mat_value_type_t<Matrix1>;
  using SizeType = mat_size_type_t<Matrix1>;
  SizeType N = F.get_row_count();
  SizeType M = R.get_row_count();

  mat<ValueType, mat_structure::rectangular> R_tmp(N * 2 + M, M);
  sub(R_tmp)(range(0, M), range(0, M)) = R;
  sub(R_tmp)(range(M, M + N), range(0, M)) = -G;
  mat<ValueType, mat_structure::square> Q_tmp =
      mat<ValueType, mat_structure::square>(N * 2 + M);
  sub(Q_tmp)(range(0, 2 * N), range(M, M + 2 * N)) =
      mat<ValueType, mat_structure::identity>(2 * N);
  sub(Q_tmp)(range(2 * N, 2 * N + M), range(0, M)) =
      mat<ValueType, mat_structure::identity>(M);

  detail::decompose_QR_impl(R_tmp, &Q_tmp, NumTol);
  Q_tmp = transpose(Q_tmp);

  mat<ValueType, mat_structure::rectangular> B_aug(2 * N, 2 * N);
  sub(B_aug)(range(0, 2 * N), range(0, N)) =
      sub(Q_tmp)(range(M, M + 2 * N), range(0, N));
  sub(B_aug)(range(0, 2 * N), range(N, 2 * N)) =
      sub(Q_tmp)(range(M, M + 2 * N), range(N, 2 * N)) * transpose_view(F) +
      sub(Q_tmp)(range(M, M + 2 * N), range(2 * N, 2 * N + M)) *
          transpose_view(G);

  mat<ValueType, mat_structure::rectangular> A_aug(2 * N, 2 * N);
  sub(A_aug)(range(0, 2 * N), range(0, N)) =
      sub(Q_tmp)(range(M, M + 2 * N), range(0, N)) * F -
      sub(Q_tmp)(range(M, M + 2 * N), range(N, 2 * N)) * Q;
  sub(A_aug)(range(0, 2 * N), range(N, 2 * N)) =
      sub(Q_tmp)(range(M, M + 2 * N), range(N, 2 * N));

  mat<ValueType, mat_structure::square> Q_aug =
      mat<ValueType, mat_structure::square>(
          mat<ValueType, mat_structure::identity>(2 * N));
  mat<ValueType, mat_structure::square> Z_aug =
      mat<ValueType, mat_structure::square>(
          mat<ValueType, mat_structure::identity>(2 * N));

  bool should_interchange = false;
  if (norm_1(A_aug) > norm_1(B_aug)) {
    should_interchange = true;
  }

  vect_n<int> Dl_aug(2 * N);
  vect_n<int> Dr_aug(2 * N);
  if (UseBalancing) {
    balance_pencil(A_aug, B_aug, Dl_aug, Dr_aug);
  }

  if (should_interchange) {
    detail::gen_schur_decomp_impl(B_aug, A_aug, &Q_aug, &Z_aug, NumTol);
  } else {
    detail::gen_schur_decomp_impl(A_aug, B_aug, &Q_aug, &Z_aug, NumTol);
  }

  if (should_interchange) {
    detail::partition_schur_pencil_impl(B_aug, A_aug, &Q_aug, &Z_aug,
                                        detail::out_unit_circle_eigen_first(),
                                        NumTol);
  } else {
    detail::partition_schur_pencil_impl(A_aug, B_aug, &Q_aug, &Z_aug,
                                        detail::in_unit_circle_eigen_first(),
                                        NumTol);
  }

  P.set_row_count(N);
  P.set_col_count(N);
  mat_sub_block<mat<ValueType, mat_structure::square>> subZ11(Z_aug, N, N, 0,
                                                              0);
  mat_sub_block<mat<ValueType, mat_structure::square>> subZ21(Z_aug, N, N, N,
                                                              0);

  try {
    linlsq_QR(transpose_view(subZ11), P, transpose_view(subZ21), NumTol);
  } catch (singularity_error& e) {
    RK_UNUSED(e);
    throw singularity_error(
        "The Discrete-time Algebraic Riccati Equation (DARE) cannot be solved! "
        "Usually indicates "
        "that the system is not stabilizable.");
  }

  if (UseBalancing) {
    apply_left_bal_inv_exp(sub(Dr_aug)[range(0, N)], P);
    apply_right_bal_exp(P, sub(Dr_aug)[range(N, 2 * N)]);
  }

  P += transpose(P);
  P *= ValueType(0.5);
}

/// Solves the Infinite-horizon Discrete-time Linear Quadratic Regulator (LQR) problem.
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981)
/// and implemented in the function solve_dare_problem.
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A fully-writable matrix type.
/// \tparam Matrix6 A fully-writable (square) matrix type.
/// \param F square (n x n) matrix which represents state-to-next-state linear map.
/// \param G rectangular (n x m) matrix which represents input-to-next-state linear map.
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-semi-definite matrix which represents quadratic input penalty.
/// \param K holds as output, the (mxn) LQR-optimal gain matrix for u = - K * (x_cur - x_ref).
/// \param P holds as output, the (nxn) nonnegative definite solution to P = F^T P F - F^T P G ( R + G^T P G )^{-1} G^T P F + Q.
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the DARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableMatrix Matrix3,
          ReadableMatrix Matrix4, FullyWritableMatrix Matrix5, FullyWritableMatrix Matrix6>
void solve_IHDT_LQR(const Matrix1& F, const Matrix2& G, const Matrix3& Q,
                    const Matrix4& R, Matrix5& K, Matrix6& P,
                    mat_value_type_t<Matrix1> NumTol = 1E-8,
                    bool UseBalancing = false) {
  solve_dare_problem(F, G, Q, R, P, NumTol, UseBalancing);

  mat<double, mat_structure::rectangular> M_tmp(R);
  M_tmp += transpose_view(G) * P * G;
  mat<double, mat_structure::rectangular> M2_tmp(G.get_col_count(),
                                                 F.get_col_count());
  M2_tmp = transpose_view(G) * P * F;
  linlsq_QR(M_tmp, K, M2_tmp);
}

/// Solves the Infinite-horizon Discrete-time Linear Quadratic Regulator (LQR) problem.
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981)
/// and implemented in the function solve_dare_problem.
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A fully-writable matrix type.
/// \param F square (n x n) matrix which represents state-to-next-state linear map.
/// \param G rectangular (n x m) matrix which represents input-to-next-state linear map.
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-semi-definite matrix which represents quadratic input penalty.
/// \param K holds as output, the (mxn) LQR-optimal gain matrix for u = - K * (x_cur - x_ref).
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the DARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableMatrix Matrix3,
          ReadableMatrix Matrix4, WritableMatrix Matrix5>
void solve_IHDT_LQR(const Matrix1& F, const Matrix2& G, const Matrix3& Q,
                    const Matrix4& R, Matrix5& K,
                    mat_value_type_t<Matrix1> NumTol = 1E-8,
                    bool UseBalancing = false) {
  using ValueType = mat_value_type_t<Matrix1>;
  mat<ValueType, mat_structure::square> P =
      mat<ValueType, mat_structure::square>(F.get_row_count());
  solve_IHDT_LQR(F, G, Q, R, K, P, NumTol, UseBalancing);
}

/// Solves the Infinite-horizon Discrete-time Linear Quadratic Gaussian control (LQG) problem.
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981)
/// and implemented in the function solve_dare_problem.
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A readable matrix type.
/// \tparam Matrix7 A readable matrix type.
/// \tparam Matrix8 A fully-writable matrix type.
/// \tparam Matrix9 A fully-writable (square) matrix type.
/// \tparam Matrix10 A fully-writable matrix type.
/// \tparam Matrix11 A fully-writable (square) matrix type.
/// \param F square (n x n) matrix which represents state-to-next-state linear map.
/// \param G rectangular (n x m) matrix which represents input-to-next-state linear map.
/// \param H rectangular (l x n) matrix which represents state-to-output linear map.
/// \param V square (n x n) positive-semi-definite matrix which represents the covariance of the state disturbances.
/// \param W square (l x l) positive-semi-definite matrix which represents the covariance of the measurement noise
///          (additive).
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-semi-definite matrix which represents quadratic input penalty.
/// \param K holds as output, the (nxl) Kalman gain matrix for x_posterior = x_prior + K * (y - x_prior).
/// \param P holds as output, the (nxn) nonnegative definite solution to P = F P F^T - F P H^T ( W + H P H^T )^{-1} H P
///          F^T + V.
/// \param L holds as output, the (mxn) LQR-optimal gain matrix for u = - K * (x_cur - x_ref).
/// \param S holds as output, the (nxn) nonnegative definite solution to S = F^T S F - F^T S G ( R + G^T S G )^{-1} G^T S
///          F + Q.
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the DARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableMatrix Matrix3,
          ReadableMatrix Matrix4, ReadableMatrix Matrix5, ReadableMatrix Matrix6,
          ReadableMatrix Matrix7, FullyWritableMatrix Matrix8, FullyWritableMatrix Matrix9,
          FullyWritableMatrix Matrix10, FullyWritableMatrix Matrix11>
void solve_IHDT_LQG(const Matrix1& F, const Matrix2& G, const Matrix3& H,
                    const Matrix4& V, const Matrix5& W, const Matrix6& Q,
                    const Matrix7& R, Matrix8& K, Matrix9& P, Matrix10& L,
                    Matrix11& S, mat_value_type_t<Matrix1> NumTol = 1E-8,
                    bool UseBalancing = false) {
  solve_dare_problem(transpose_view(F), transpose_view(H), V, W, P, NumTol,
                     UseBalancing);

  mat<double, mat_structure::rectangular> M_tmp = W;
  M_tmp += H * P * transpose_view(H);
  mat<double, mat_structure::rectangular> M2_tmp(H.get_row_count(),
                                                 F.get_col_count());
  M2_tmp = H * P * transpose_view(F);
  mat<double, mat_structure::rectangular> Msol_tmp;
  linlsq_QR(M_tmp, Msol_tmp, M2_tmp);
  K = transpose_view(Msol_tmp);

  solve_dare_problem(F, G, Q, R, S, NumTol, UseBalancing);

  mat<double, mat_structure::rectangular> M3_tmp = R;
  M3_tmp += transpose_view(G) * S * G;
  mat<double, mat_structure::rectangular> M4_tmp(G.get_col_count(),
                                                 F.get_col_count());
  M4_tmp = transpose_view(G) * S * F;
  linlsq_QR(M3_tmp, L, M4_tmp);
}

/// Solves the Infinite-horizon Discrete-time Linear Quadratic Gaussian control (LQG) problem.
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981)
/// and implemented in the function solve_dare_problem.
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A readable matrix type.
/// \tparam Matrix7 A readable matrix type.
/// \tparam Matrix8 A fully-writable matrix type.
/// \tparam Matrix9 A fully-writable matrix type.
/// \param F square (n x n) matrix which represents state-to-next-state linear map.
/// \param G rectangular (n x m) matrix which represents input-to-next-state linear map.
/// \param H rectangular (l x n) matrix which represents state-to-output linear map.
/// \param V square (n x n) positive-semi-definite matrix which represents the covariance of the state disturbances.
/// \param W square (l x l) positive-semi-definite matrix which represents the covariance of the measurement noise
///          (additive).
/// \param Q square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param R square (m x m) positive-semi-definite matrix which represents quadratic input penalty.
/// \param K holds as output, the (nxl) Kalman gain matrix for x_posterior = x_prior + K * (y - x_prior).
/// \param L holds as output, the (mxn) LQR-optimal gain matrix for u = - K * (x_cur - x_ref).
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the DARE problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableMatrix Matrix3,
          ReadableMatrix Matrix4, ReadableMatrix Matrix5, ReadableMatrix Matrix6,
          ReadableMatrix Matrix7, WritableMatrix Matrix8, WritableMatrix Matrix9>
void solve_IHDT_LQG(const Matrix1& F, const Matrix2& G, const Matrix3& H,
                    const Matrix4& V, const Matrix5& W, const Matrix6& Q,
                    const Matrix7& R, Matrix8& K, Matrix9& L,
                    mat_value_type_t<Matrix1> NumTol = 1E-8,
                    bool UseBalancing = false) {
  using ValueType = mat_value_type_t<Matrix1>;
  mat<ValueType, mat_structure::square> P =
      mat<ValueType, mat_structure::square>(F.get_row_count());
  mat<ValueType, mat_structure::square> S =
      mat<ValueType, mat_structure::square>(F.get_row_count());
  solve_IHDT_LQG(F, G, H, V, W, Q, R, K, P, L, S, NumTol, UseBalancing);
}

/// Solves the Continuous-time Spectral Factorisation of a system.
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981).
/// This method first reduces the augmented (2n+m x 2n+m) pencil to a (2n x 2n) pencil
/// using a QR decomposition on the last (mxm) block-column (which has infinite eigenvalues).
/// Then, it performs a generalized real Schur decomposition of the pencil. Finally, it
/// reorders the eigenvalues in the pencil such that stable (within unit-circle) eigenvalues
/// percolate to the upper (nxn) pencil, which allows the extraction of the eigenvectors
/// spanning the stable subspace, which are, in turn, used to compute the unique solution P.
/// \n
/// $B (D + D^T)^{-1} B^T + P (A - B (D + D^T)^{-1} C)^T + (A - B (D + D^T)^{-1} C) P + P C^T (D + D^T)^{-1} C P = 0$
/// \n
/// The initial pencil is: lambda * (I 0 0; 0 I 0; 0 0 0) - (A 0 B; 0 -A^T C^T; C -B^T (D + D^T))
/// \n
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A fully-writable matrix type.
/// \param A square (n x n) matrix which represents state-to-state-derivative linear map.
/// \param B rectangular (n x m) matrix which represents input-to-state-derivative linear map.
/// \param C rectangular (m x n) matrix which represents state-to-output linear map.
/// \param D square (m x m) matrix which represents input-to-output linear map.
/// \param P holds as output, the nonnegative definite solution.
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the CTSF problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableMatrix Matrix3,
          ReadableMatrix Matrix4, FullyWritableMatrix Matrix5>
void solve_ctsf_problem(const Matrix1& A, const Matrix2& B, const Matrix3& C,
                        const Matrix4& D, Matrix5& P,
                        mat_value_type_t<Matrix1> NumTol = 1E-8,
                        bool UseBalancing = false) {
  if ((A.get_row_count() != A.get_col_count()) ||
      (B.get_row_count() != A.get_row_count()) ||
      (C.get_col_count() != A.get_col_count()) ||
      (C.get_row_count() != D.get_col_count()) ||
      (D.get_row_count() != D.get_col_count()) ||
      (B.get_col_count() != D.get_col_count())) {
    throw std::range_error(
        "The dimensions of the CTSF system matrices do not match! Should be "
        "A(n x n), B(n x m), "
        "C(m x n), and D(m x m).");
  }

  using ValueType = mat_value_type_t<Matrix1>;
  using SizeType = mat_size_type_t<Matrix1>;
  SizeType N = A.get_row_count();
  SizeType M = D.get_row_count();

  mat<ValueType, mat_structure::rectangular> R_tmp(N * 2 + M, M);
  sub(R_tmp)(range(0, M), range(0, M)) = D + transpose_view(D);
  sub(R_tmp)(range(M, M + N), range(0, M)) = B;
  sub(R_tmp)(range(M + N, M + 2 * N), range(0, M)) = transpose_view(C);
  mat<ValueType, mat_structure::square> Q_tmp =
      mat<ValueType, mat_structure::square>(N * 2 + M);
  sub(Q_tmp)(range(0, 2 * N), range(M, M + 2 * N)) =
      mat<ValueType, mat_structure::identity>(2 * N);
  sub(Q_tmp)(range(2 * N, 2 * N + M), range(0, M)) =
      mat<ValueType, mat_structure::identity>(M);

  detail::decompose_QR_impl(R_tmp, &Q_tmp, NumTol);
  Q_tmp = transpose(Q_tmp);

  mat<ValueType, mat_structure::rectangular> B_aug(2 * N, 2 * N);
  B_aug = sub(Q_tmp)(range(M, M + 2 * N), range(0, 2 * N));

  mat<ValueType, mat_structure::rectangular> A_aug(2 * N, 2 * N);
  sub(A_aug)(range(0, 2 * N), range(0, N)) =
      sub(Q_tmp)(range(M, M + 2 * N), range(0, N)) * A +
      sub(Q_tmp)(range(M, M + 2 * N), range(2 * N, 2 * N + M)) * C;
  sub(A_aug)(range(0, 2 * N), range(N, 2 * N)) =
      -sub(Q_tmp)(range(M, M + 2 * N), range(2 * N, 2 * N + M)) *
          transpose_view(B) -
      sub(Q_tmp)(range(M, M + 2 * N), range(N, 2 * N)) * transpose_view(A);

  mat<ValueType, mat_structure::square> Q_aug =
      mat<ValueType, mat_structure::square>(
          mat<ValueType, mat_structure::identity>(2 * N));
  mat<ValueType, mat_structure::square> Z_aug =
      mat<ValueType, mat_structure::square>(
          mat<ValueType, mat_structure::identity>(2 * N));

  bool should_interchange = false;
  if (norm_1(A_aug) > norm_1(B_aug)) {
    should_interchange = true;
  }

  vect_n<int> Dl_aug(2 * N);
  vect_n<int> Dr_aug(2 * N);
  if (UseBalancing) {
    balance_pencil(A_aug, B_aug, Dl_aug, Dr_aug);
  }

  if (should_interchange) {
    detail::gen_schur_decomp_impl(B_aug, A_aug, &Q_aug, &Z_aug, NumTol);
  } else {
    detail::gen_schur_decomp_impl(A_aug, B_aug, &Q_aug, &Z_aug, NumTol);
  }

  if (should_interchange) {
    detail::partition_schur_pencil_impl(B_aug, A_aug, &Q_aug, &Z_aug,
                                        detail::neg_real_val_eigen_first(),
                                        NumTol);
  } else {
    detail::partition_schur_pencil_impl(A_aug, B_aug, &Q_aug, &Z_aug,
                                        detail::neg_real_val_eigen_first(),
                                        NumTol);
  }

  P.set_row_count(N);
  P.set_col_count(N);
  mat_sub_block<mat<ValueType, mat_structure::square>> subZ11(Z_aug, N, N, 0,
                                                              0);
  mat_sub_block<mat<ValueType, mat_structure::square>> subZ21(Z_aug, N, N, N,
                                                              0);

  try {
    linlsq_QR(transpose_view(subZ11), P, transpose_view(subZ21), NumTol);
  } catch (singularity_error& e) {
    throw singularity_error(
        "The Continuous-time Spectral Factorisation (CTSF) cannot be solved! "
        "Usually indicates "
        "that the system is not stabilizable.");
  }

  if (UseBalancing) {
    apply_left_bal_inv_exp(sub(Dr_aug)[range(0, N)], P);
    apply_right_bal_exp(P, sub(Dr_aug)[range(N, 2 * N)]);
  }

  P += transpose(P);
  P *= ValueType(0.5);
}

/// Solves the Discrete-time Spectral Factorisation of a system.
/// This implementation uses the QZ-algorithm approach as described in Van Dooren (1981).
/// This method first reduces the augmented (2n+m x 2n+m) pencil to a (2n x 2n) pencil
/// using a QR decomposition on the last (mxm) block-column (which has infinite eigenvalues).
/// Then, it performs a generalized real Schur decomposition of the pencil. Finally, it
/// reorders the eigenvalues in the pencil such that stable (within unit-circle) eigenvalues
/// percolate to the upper (nxn) pencil, which allows the extraction of the eigenvectors
/// spanning the stable subspace, which are, in turn, used to compute the unique solution P.
/// \n
/// $P = F P F^T + (G - F P H^T) ( J + J^T - H P H^T )^{-1} (G^T - H P F^T)$
/// \n
/// The initial pencil is: lambda * (I 0 0; 0 F^T 0; 0 G^T 0) - (F 0 -G; -Q I 0; 0 0 R)
/// \n
///
/// \tparam Matrix1 A readable matrix type.
/// \tparam Matrix2 A readable matrix type.
/// \tparam Matrix3 A readable matrix type.
/// \tparam Matrix4 A readable matrix type.
/// \tparam Matrix5 A fully-writable matrix type.
/// \param F square (n x n) matrix which represents state-to-next-state linear map.
/// \param G rectangular (n x m) matrix which represents input-to-next-state linear map.
/// \param H square (n x n) positive-definite matrix which represents quadratic state-error penalty.
/// \param J square (m x m) positive-semi-definite matrix which represents quadratic input penalty.
/// \param P holds as output, the nonnegative definite solution to P = F^T P F - F^T P G ( R + G^T P G )^{-1} G^T P F + Q.
/// \param NumTol tolerance for considering a value to be zero in avoiding divisions
///               by zero and singularities.
/// \param UseBalancing specifies whether balancing should be applied to the problem before performing
///                     the Schur decomposition. This can help for ill-conditioned systems to increase
///                     the final accuracy (accumulated round-off errors). However, for certain systems,
///                     balancing will worsen the results, so, use with caution. By default, no balancing
///                     is performed.
///
/// \throws std::range_error if the matrix dimensions are not consistent.
/// \throws singularity_error if the DTSF problem cannot be solved, usually because the system is not stabilizable.
///
/// \author Mikael Persson
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2, ReadableMatrix Matrix3,
          ReadableMatrix Matrix4, FullyWritableMatrix Matrix5>
void solve_dtsf_problem(const Matrix1& F, const Matrix2& G, const Matrix3& H,
                        const Matrix4& J, Matrix5& P,
                        mat_value_type_t<Matrix1> NumTol = 1E-8,
                        bool UseBalancing = false) {
  if ((F.get_row_count() != F.get_col_count()) ||
      (G.get_row_count() != F.get_row_count()) ||
      (H.get_col_count() != F.get_col_count()) ||
      (H.get_row_count() != J.get_col_count()) ||
      (J.get_row_count() != J.get_col_count()) ||
      (G.get_col_count() != J.get_col_count())) {
    throw std::range_error(
        "The dimensions of the DTSF system matrices do not match! Should be "
        "F(n x n), G(n x m), "
        "H(m x n), and J(m x m).");
  }

  using ValueType = mat_value_type_t<Matrix1>;
  using SizeType = mat_size_type_t<Matrix1>;
  SizeType N = F.get_row_count();
  SizeType M = J.get_row_count();

  mat<ValueType, mat_structure::rectangular> R_tmp(N * 2 + M, M);
  sub(R_tmp)(range(0, M), range(0, M)) = -(J + transpose_view(J));
  sub(R_tmp)(range(M, M + N), range(0, M)) = -G;
  sub(R_tmp)(range(M + N, M + 2 * N), range(0, M)) = -transpose_view(H);
  mat<ValueType, mat_structure::square> Q_tmp =
      mat<ValueType, mat_structure::square>(N * 2 + M);
  sub(Q_tmp)(range(0, 2 * N), range(M, M + 2 * N)) =
      mat<ValueType, mat_structure::identity>(2 * N);
  sub(Q_tmp)(range(2 * N, 2 * N + M), range(0, M)) =
      mat<ValueType, mat_structure::identity>(M);

  detail::decompose_QR_impl(R_tmp, &Q_tmp, NumTol);
  Q_tmp = transpose(Q_tmp);

  mat<ValueType, mat_structure::rectangular> B_aug(2 * N, 2 * N);
  sub(B_aug)(range(0, 2 * N), range(0, N)) =
      sub(Q_tmp)(range(M, M + 2 * N), range(0, N));
  sub(B_aug)(range(0, 2 * N), range(N, 2 * N)) =
      sub(Q_tmp)(range(M, M + 2 * N), range(N, 2 * N)) * transpose_view(F) +
      sub(Q_tmp)(range(M, M + 2 * N), range(2 * N, 2 * N + M)) *
          transpose_view(G);

  mat<ValueType, mat_structure::rectangular> A_aug(2 * N, 2 * N);
  sub(A_aug)(range(0, 2 * N), range(0, N)) =
      sub(Q_tmp)(range(M, M + 2 * N), range(0, N)) * F +
      sub(Q_tmp)(range(M, M + 2 * N), range(2 * N, 2 * N + M)) * H;
  sub(A_aug)(range(0, 2 * N), range(N, 2 * N)) =
      sub(Q_tmp)(range(M, M + 2 * N), range(N, 2 * N));

  mat<ValueType, mat_structure::square> Q_aug =
      mat<ValueType, mat_structure::square>(
          mat<ValueType, mat_structure::identity>(2 * N));
  mat<ValueType, mat_structure::square> Z_aug =
      mat<ValueType, mat_structure::square>(
          mat<ValueType, mat_structure::identity>(2 * N));

  bool should_interchange = false;
  if (norm_1(A_aug) > norm_1(B_aug)) {
    should_interchange = true;
  }

  vect_n<int> Dl_aug(2 * N);
  vect_n<int> Dr_aug(2 * N);
  if (UseBalancing) {
    balance_pencil(A_aug, B_aug, Dl_aug, Dr_aug);
  }

  if (should_interchange) {
    detail::gen_schur_decomp_impl(B_aug, A_aug, &Q_aug, &Z_aug, NumTol);
  } else {
    detail::gen_schur_decomp_impl(A_aug, B_aug, &Q_aug, &Z_aug, NumTol);
  }

  if (should_interchange) {
    detail::partition_schur_pencil_impl(B_aug, A_aug, &Q_aug, &Z_aug,
                                        detail::out_unit_circle_eigen_first(),
                                        NumTol);
  } else {
    detail::partition_schur_pencil_impl(A_aug, B_aug, &Q_aug, &Z_aug,
                                        detail::in_unit_circle_eigen_first(),
                                        NumTol);
  }

  P.set_row_count(N);
  P.set_col_count(N);
  mat_sub_block<mat<ValueType, mat_structure::square>> subZ11(Z_aug, N, N, 0,
                                                              0);
  mat_sub_block<mat<ValueType, mat_structure::square>> subZ21(Z_aug, N, N, N,
                                                              0);

  try {
    linlsq_QR(transpose_view(subZ11), P, transpose_view(subZ21), NumTol);
  } catch (singularity_error& e) {
    throw singularity_error(
        "The Discrete-time Spectral Factorisation (CTSF) cannot be solved! "
        "Usually indicates that "
        "the system is not stabilizable.");
  }

  if (UseBalancing) {
    apply_left_bal_inv_exp(sub(Dr_aug)[range(0, N)], P);
    apply_right_bal_exp(P, sub(Dr_aug)[range(N, 2 * N)]);
  }

  P += transpose(P);
  P *= ValueType(0.5);
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_ARE_SOLVER_H_

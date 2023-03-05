/**
 * \file mat_balance.hpp
 *
 * This library provides the matrix balancing function template which can be used
 * to scale a square matrix such that its rows and columns all have a comparable
 * order of magnitude. Applying this method before some matrix decomposition algorithms
 * increases the numerical stability and convergence rates (for iterative methods).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2011
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

#ifndef REAK_MAT_BALANCE_HPP
#define REAK_MAT_BALANCE_HPP

#include "ReaK/math/lin_alg/mat_alg.hpp"
#include "ReaK/math/lin_alg/mat_num_exceptions.hpp"

#include <type_traits>

namespace ReaK {

/**
 * Performs matrix balancing of a square matrix such that the sum of absolute values of row elements and
 * column elements match in order of magnitude. This algorithm produces a diagonal matrix that can
 * scale the matrix, as A_balanced = D^-1 A D
 *
 * \param A square matrix with row-count == column-count, and stores, as output, the balanced matrix.
 * \param D holds as output, the diagonal matrix D which balances A as a vector of base-2 exponents.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 *
 * Taken from Golub & Van Loan, "Matrix Computations" (3rd ed).
 */
template <typename Matrix1>
void balance(Matrix1& A, vect_n<int>& D) {
  static_assert(is_fully_writable_matrix_v<Matrix1>);
  if (A.get_row_count() < A.get_col_count()) {
    throw std::range_error(
        "Matrix balancing is only possible on a square matrix!");
  }

  using std::abs;
  using ValueType = mat_value_type_t<Matrix1>;
  using std::frexp;
  using std::ldexp;

  int N = A.get_row_count();
  D = vect_n<int>(N, 0);

  bool keep_going = true;

  while (keep_going) {

    keep_going = false;

    for (int i = 0; i < N; ++i) {
      ValueType row_mag = ValueType();
      ValueType col_mag = ValueType();

      for (int j = 0; j < N; ++j) {
        if (j != i) {
          col_mag += abs(A(j, i));
          row_mag += abs(A(i, j));
        }
      }

      if ((col_mag < std::numeric_limits<ValueType>::epsilon()) ||
          (row_mag < std::numeric_limits<ValueType>::epsilon())) {
        continue;
      }

      ValueType g = ldexp(row_mag, -1);
      int f = 0;
      ValueType s = col_mag + row_mag;

      while (col_mag < g) {
        ++f;
        col_mag = ldexp(col_mag, 2);
      }

      g = ldexp(row_mag, 1);
      while (col_mag > g) {
        --f;
        col_mag = ldexp(col_mag, -2);
      }

      if ((row_mag + col_mag) < ValueType(0.95) * ldexp(s, f)) {
        keep_going = true;
        for (int j = 0; j < N; ++j) {
          A(j, i) = ldexp(A(j, i), f);
        }
        for (int j = 0; j < N; ++j) {
          A(i, j) = ldexp(A(i, j), -f);
        }
        D[i] += f;
      }
    }
  }
}

/**
 * Performs matrix balancing of a matrix pencil such that the sum of absolute values of row elements and
 * column elements match in order of magnitude. This algorithm produces a left and right diagonal matrix
 * that can scale the matrix, as (A,B)_balanced = Dl^-1 (A,B) Dr
 *
 * \param A square matrix with row-count == column-count, and stores, as output, the balanced matrix.
 * \param B square matrix with row-count == column-count, and stores, as output, the balanced matrix.
 * \param Dl holds as output, the left-diagonal matrix Dl.
 * \param Dr holds as output, the left-diagonal matrix Dr.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 *
 * Taken from Lemonnier and Van Dooren 2006.
 */
template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4>
std::enable_if_t<is_writable_matrix_v<Matrix3> && is_writable_matrix_v<Matrix4>>
balance_pencil(Matrix1& A, Matrix2& B, Matrix3& Dl, Matrix4& Dr) {
  static_assert(is_fully_writable_matrix_v<Matrix1>);
  static_assert(is_fully_writable_matrix_v<Matrix2>);
  if ((A.get_row_count() != A.get_col_count()) ||
      (B.get_row_count() != A.get_row_count()) ||
      (B.get_row_count() != B.get_col_count())) {
    throw std::range_error(
        "Matrix pencil balancing is only possible on square matrices!");
  }

  using ValueType = mat_value_type_t<Matrix1>;
  using std::frexp;
  using std::ldexp;

  int N = A.get_row_count();

  Dl = mat<mat_value_type_t<Matrix3>, mat_structure::identity>(N);
  Dr = mat<mat_value_type_t<Matrix4>, mat_structure::identity>(N);

  mat<ValueType, mat_structure::square> M =
      mat<ValueType, mat_structure::square>(N);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      M(i, j) = A(i, j) * A(i, j) + B(i, j) * B(i, j);
    }
  }

  while (true) {
    int e_min = 0;
    int e_max = 0;

    for (int i = 0; i < N; ++i) {
      ValueType sum = ValueType(0.0);
      for (int j = 0; j < N; ++j) {
        sum += M(i, j);
      }
      int e = 0;
      frexp(sum, &e);
      e = -(e / 2);  // using integer arithmetic
      for (int j = 0; j < N; ++j) {
        A(i, j) = ldexp(A(i, j), e);
        B(i, j) = ldexp(B(i, j), e);
        M(i, j) = ldexp(M(i, j), 2 * e);
      }
      Dl(i, i) = ldexp(Dl(i, i), -e);
      if (e > e_max) {
        e_max = e;
      }
      if (e < e_min) {
        e_min = e;
      }
    }

    for (int i = 0; i < N; ++i) {
      ValueType sum = ValueType(0.0);
      for (int j = 0; j < N; ++j) {
        sum += M(j, i);
      }
      int e = 0;
      frexp(sum, &e);
      e = -(e / 2);  // using integer arithmetic
      for (int j = 0; j < N; ++j) {
        A(j, i) = ldexp(A(j, i), e);
        B(j, i) = ldexp(B(j, i), e);
        M(j, i) = ldexp(M(j, i), 2 * e);
      }
      Dr(i, i) = ldexp(Dr(i, i), e);
      if (e > e_max) {
        e_max = e;
      }
      if (e < e_min) {
        e_min = e;
      }
    }

    if (e_max <= e_min + 2) {
      break;
    }
  }
}

/**
 * Performs matrix balancing of a matrix pencil such that the sum of absolute values of row elements and
 * column elements match in order of magnitude. This algorithm produces a left and right vector
 * of base-2 exponents that scale the matrix, as (A,B)_balanced = Dl^-1 (A,B) Dr
 *
 * \param A square matrix with row-count == column-count, and stores, as output, the balanced matrix.
 * \param B square matrix with row-count == column-count, and stores, as output, the balanced matrix.
 * \param Dl holds as output, the left-diagonal matrix Dl as a vector of base-2 exponents.
 * \param Dr holds as output, the right-diagonal matrix Dr as a vector of base-2 exponents.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 *
 * Taken from Lemonnier and Van Dooren 2006.
 */
template <typename Matrix1, typename Matrix2>
void balance_pencil(Matrix1& A, Matrix2& B, vect_n<int>& Dl, vect_n<int>& Dr) {
  static_assert(is_fully_writable_matrix_v<Matrix1>);
  static_assert(is_fully_writable_matrix_v<Matrix2>);
  if ((A.get_row_count() != A.get_col_count()) ||
      (B.get_row_count() != A.get_row_count()) ||
      (B.get_row_count() != B.get_col_count())) {
    throw std::range_error(
        "Matrix pencil balancing is only possible on square matrices!");
  }

  using ValueType = mat_value_type_t<Matrix1>;
  using std::frexp;
  using std::ldexp;

  int N = A.get_row_count();

  Dl = vect_n<int>(N, 0);
  Dr = vect_n<int>(N, 0);

  mat<ValueType, mat_structure::square> M =
      mat<ValueType, mat_structure::square>(N);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      M(i, j) = A(i, j) * A(i, j) + B(i, j) * B(i, j);
    }
  }

  while (true) {
    int e_min = 0;
    int e_max = 0;

    for (int i = 0; i < N; ++i) {
      auto sum = ValueType(0.0);
      for (int j = 0; j < N; ++j) {
        sum += M(i, j);
      }
      int e = 0;
      frexp(sum, &e);
      e = -(e / 2);  // using integer arithmetic
      for (int j = 0; j < N; ++j) {
        A(i, j) = ldexp(A(i, j), e);
        B(i, j) = ldexp(B(i, j), e);
        M(i, j) = ldexp(M(i, j), 2 * e);
      }
      Dl[i] -= e;
      if (e > e_max) {
        e_max = e;
      }
      if (e < e_min) {
        e_min = e;
      }
    }

    for (int i = 0; i < N; ++i) {
      auto sum = ValueType(0.0);
      for (int j = 0; j < N; ++j) {
        sum += M(j, i);
      }
      int e = 0;
      frexp(sum, &e);
      e = -(e / 2);  // using integer arithmetic
      for (int j = 0; j < N; ++j) {
        A(j, i) = ldexp(A(j, i), e);
        B(j, i) = ldexp(B(j, i), e);
        M(j, i) = ldexp(M(j, i), 2 * e);
      }
      Dr[i] += e;
      if (e > e_max) {
        e_max = e;
      }
      if (e < e_min) {
        e_min = e;
      }
    }

    if (e_max <= e_min + 2) {
      break;
    }
  }
}

/**
 * Performs matrix balancing of a matrix pencil such that the sum of absolute values of row elements and
 * column elements match in order of magnitude. This algorithm produces a left and right diagonal matrix
 * that can scale the matrix, as (A,B)_balanced = Dl^-1 (A,B) Dr
 *
 * \param A square matrix with row-count == column-count, and stores, as output, the balanced matrix.
 * \param B square matrix with row-count == column-count, and stores, as output, the balanced matrix.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 *
 * Taken from Lemonnier and Van Dooren 2006.
 */
template <typename Matrix1, typename Matrix2>
void balance_pencil(Matrix1& A, Matrix2& B) {
  vect_n<int> Dl(A.get_row_count());
  vect_n<int> Dr(A.get_row_count());
  balance_pencil(A, B, Dl, Dr);
}

/**
 * This function applies a matrix balance factor from the left of a matrix. The
 * matrix balance factor is represented as a vector of integer base-2 exponents
 * to be applied with exact arithmetic.
 *
 * \tparam Matrix1 A fully-writable matrix type.
 * \tparam Vector1 An integer vector type.
 * \param D The vector of integer base-2 exponents to be applied to the rows of A.
 * \param A The matrix to which the balance factor is applied.
 *
 * \throw std::range_error If the dimensions of A and D do not match.
 */
template <typename Matrix1, typename Vector1>
void apply_left_bal_exp(const Vector1& D, Matrix1& A) {
  static_assert(is_fully_writable_matrix_v<Matrix1>);
  static_assert(is_readable_vector_v<Vector1>);
  using std::ldexp;
  if (A.get_row_count() != D.size()) {
    throw std::range_error(
        "Matrix balancing factor does not match the dimension of the matrix!");
  }

  for (int i = 0; i < A.get_row_count(); ++i) {
    for (int j = 0; j < A.get_col_count(); ++j) {
      A(i, j) = ldexp(A(i, j), D[i]);
    }
  }
}

/**
 * This function applies a matrix balance factor from the right of a matrix. The
 * matrix balance factor is represented as a vector of integer base-2 exponents
 * to be applied with exact arithmetic.
 *
 * \tparam Matrix1 A fully-writable matrix type.
 * \tparam Vector1 An integer vector type.
 * \param D The vector of integer base-2 exponents to be applied to the columns of A.
 * \param A The matrix to which the balance factor is applied.
 *
 * \throw std::range_error If the dimensions of A and D do not match.
 */
template <typename Matrix1, typename Vector1>
void apply_right_bal_exp(Matrix1& A, const Vector1& D) {
  static_assert(is_fully_writable_matrix_v<Matrix1>);
  static_assert(is_readable_vector_v<Vector1>);
  using std::ldexp;
  if (A.get_col_count() != D.size()) {
    throw std::range_error(
        "Matrix balancing factor does not match the dimension of the matrix!");
  }

  for (int j = 0; j < A.get_col_count(); ++j) {
    for (int i = 0; i < A.get_row_count(); ++i) {
      A(i, j) = ldexp(A(i, j), D[j]);
    }
  }
}

/**
 * This function applies an inverse matrix balance factor from the left of a matrix. The
 * matrix balance factor is represented as a vector of integer base-2 exponents
 * to be applied with exact arithmetic, the exponents are inverted.
 *
 * \tparam Matrix1 A fully-writable matrix type.
 * \tparam Vector1 An integer vector type.
 * \param D The vector of integer base-2 exponents to be applied, inverted, to the rows of A.
 * \param A The matrix to which the balance factor is applied.
 *
 * \throw std::range_error If the dimensions of A and D do not match.
 */
template <typename Matrix1, typename Vector1>
void apply_left_bal_inv_exp(const Vector1& D, Matrix1& A) {
  static_assert(is_fully_writable_matrix_v<Matrix1>);
  static_assert(is_readable_vector_v<Vector1>);
  using std::ldexp;
  if (A.get_row_count() != D.size()) {
    throw std::range_error(
        "Matrix balancing factor does not match the dimension of the matrix!");
  }

  for (int i = 0; i < A.get_row_count(); ++i) {
    for (int j = 0; j < A.get_col_count(); ++j) {
      A(i, j) = ldexp(A(i, j), -D[i]);
    }
  }
}

/**
 * This function applies an inverse matrix balance factor from the right of a matrix. The
 * matrix balance factor is represented as a vector of integer base-2 exponents
 * to be applied with exact arithmetic, the exponents are inverted.
 *
 * \tparam Matrix1 A fully-writable matrix type.
 * \tparam Vector1 An integer vector type.
 * \param D The vector of integer base-2 exponents to be applied, inverted, to the columns of A.
 * \param A The matrix to which the balance factor is applied.
 *
 * \throw std::range_error If the dimensions of A and D do not match.
 */
template <typename Matrix1, typename Vector1>
void apply_right_bal_inv_exp(Matrix1& A, const Vector1& D) {
  static_assert(is_fully_writable_matrix_v<Matrix1>);
  static_assert(is_readable_vector_v<Vector1>);
  using std::ldexp;
  if (A.get_col_count() != D.size()) {
    throw std::range_error(
        "Matrix balancing factor does not match the dimension of the matrix!");
  }

  for (int j = 0; j < A.get_col_count(); ++j) {
    for (int i = 0; i < A.get_row_count(); ++i) {
      A(i, j) = ldexp(A(i, j), -D[j]);
    }
  }
}

}  // namespace ReaK

#endif

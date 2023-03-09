/**
 * \file mat_gaussian_elim.h
 *
 * This library provides a number of functions related to performing a Gaussian elimination on a
 * matrix, e.g., to invert a matrix and to solve a linear system. Most implementations provided
 * are based on the PLU decomposition (Crout's method). PLU decomposition is pretty efficient and
 * is generally preferred if there is good reasons to believe that the matrix involved will always
 * be well-conditioned. If a matrix cannot be guaranteed to be well-conditioned, algorithms such as
 * QR-decomposition or methods specific to symmetric matrices are preferred (Cholesky), or even SVD.
 *
 * According to performance tests, PLU methods are as good as Cholesky methods in terms of speed.
 * And they are both the best for well-conditioned matrices. For ill-conditioned matrices, QR-decomposition
 * methods are only a little slower then PLU (about 20% slower, same time-complexity) but provide better
 * numerical stability. The Jacobi methods are significantly slower, but this implementation is in need
 * of a revision for performance enhancement. And, of course, SVD is also very slow (slightly faster than
 * Jacobi) but it is based on a LAPACK implementation that is very poorly written, and it has not been
 * updated since.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2011
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

#ifndef REAK_MATH_LIN_ALG_MAT_GAUSSIAN_ELIM_H_
#define REAK_MATH_LIN_ALG_MAT_GAUSSIAN_ELIM_H_

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_num_exceptions.h"

#include "ReaK/core/base/shared_object.h"

namespace ReaK {

/*************************************************************************
                        Gaussian Elimination
*************************************************************************/

/**
 * Inverts a matrix using the Gauss-Jordan elimination on the identity matrix.
 * \note that PLU decomposition or any other method is faster for matrix sizes of more than 20x20.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A fully-writable matrix type.
 * \param A A well-conditioned, square (Size x Size), real, full-rank matrix to be inverted.
 * \param A_inv The matrix which stores, as output, the inverse of A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero.
 *
 * \throws singularity_error if the matrix A is numerically singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not a square matrix.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
void invert_gaussian(const Matrix1& A, Matrix2& A_inv,
                     mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix1>);
  static_assert(is_writable_matrix_v<Matrix2>);

  if constexpr (!is_fully_writable_matrix_v<Matrix2>) {
    mat<mat_value_type_t<Matrix2>, mat_structure::rectangular> A_inv_tmp(
        A_inv.get_row_count(), A_inv.get_col_count());
    invert_gaussian(A, A_inv_tmp, NumTol);
    A_inv = A_inv_tmp;
  } else {
    using ValueType = mat_value_type_t<Matrix1>;
    using std::abs;

    if (A.get_col_count() != A.get_row_count()) {
      throw std::range_error("Inversion impossible! Matrix A is not square!");
    }

    int Size = A.get_col_count();
    Matrix2 tmp(A);
    A_inv = mat<ValueType, mat_structure::identity>(Size);

    for (int i = 0; i < Size; ++i) {
      if (abs(tmp(i, i)) < NumTol) {
        for (int j = i + 1; j <= Size; ++j) {
          if (j == Size) {
            throw singularity_error("M");
          }
          if (abs(tmp(j, i)) > NumTol) {
            for (int k = i; k < Size; ++k) {
              tmp(i, k) += tmp(j, k);
            }
            for (int k = 0; k < Size; ++k) {
              A_inv(i, k) += A_inv(j, k);
            }
            break;
          }
        }
      }

      ValueType s = tmp(i, i);
      for (int k = i; k < Size; ++k) {
        tmp(i, k) /= s;
      }
      for (int k = 0; k < Size; ++k) {
        A_inv(i, k) /= s;
      }

      for (int j = 0; j < Size; ++j) {
        if (i != j) {
          s = tmp(j, i);
          for (int k = i; k < Size; ++k) {
            tmp(j, k) -= s * tmp(i, k);
          }
          for (int k = 0; k < Size; ++k) {
            A_inv(j, k) -= s * A_inv(i, k);
          }
        }
      }
    }
  }
}

/*************************************************************************
             Permuted Lower- / Upper-triangular Decomposition
*************************************************************************/

namespace detail {

template <typename Matrix1, typename Matrix2, typename IndexVector>
void linsolve_PLU_impl(Matrix1& A, Matrix2& b, IndexVector& P,
                       mat_value_type_t<Matrix1> NumTol) {
  using std::abs;
  using std::swap;
  using ValueType = mat_value_type_t<Matrix1>;

  mat<ValueType, mat_structure::rectangular, mat_alignment::column_major> s(
      b.get_row_count(), b.get_col_count());
  int An = A.get_row_count();
  int bn = b.get_col_count();

  for (int i = 0; i < An; ++i) {
    P[i] = i;
    s(i, 0) = 0.0;
    for (int j = 0; j < An; ++j) {
      if (s(i, 0) < abs(A(i, j))) {
        s(i, 0) = abs(A(i, j));
      }
    }
  }

  for (int k = 0; k < An; ++k) {

    for (int i = k; i < An; ++i) {
      for (int j = 0; j < k; ++j) {
        A(i, k) -= A(i, j) * A(j, k);
      }
    }

    int temp_i = k;
    for (int i = k + 1; i < An; ++i) {
      if (abs(A(i, k) / s(i, 0)) > abs(A(temp_i, k) / s(temp_i, 0))) {
        temp_i = i;
      }
    }

    if (k != temp_i) {
      for (int i = 0; i < An; ++i) {
        swap(A(k, i), A(temp_i, i));
      }
      swap(s(k, 0), s(temp_i, 0));
      swap(P[k], P[temp_i]);
    }

    for (int j = k + 1; j < An; ++j) {
      for (int i = 0; i < k; ++i) {
        A(k, j) -= A(k, i) * A(i, j);
      }
      if (abs(A(k, k)) < NumTol) {
        throw singularity_error("A");
      }
      A(k, j) /= A(k, k);
    }
  }

  // Back-substitution
  for (int k = 0; k < An; ++k) {
    for (int l = 0; l < bn; ++l) {
      s(P[k], l) = b(P[k], l);
    }

    for (int l = 0; l < bn; ++l) {
      for (int j = 0; j < k; ++j) {
        s(k, l) -= A(k, j) * s(j, l);
      }
      s(k, l) /= A(k, k);
    }
  }

  b = s;

  for (int k = An; k > 0; --k) {
    for (int l = 0; l < bn; ++l) {
      for (int j = k; j < An; ++j) {
        b(k - 1, l) -= A(k - 1, j) * b(j, l);
      }
    }
  }
}
}  // namespace detail

/**
 * Solves the linear problem AX = B using PLU decomposition as defined by Crout`s method.
 * \note To solve a linear system involving vectors for X and B, use the Matrix-Vector Adaptors
 *(mat_vector_adaptor.hpp).
 *
 * \tparam Matrix1 A writable matrix type.
 * \tparam Matrix2 A writable matrix type.
 * \tparam IndexVector A writable vector type.
 * \param A well-conditioned, square (Size x Size), real, full-rank matrix which multiplies x.
 *          As output, A stores the LU decomposition, permutated by P.
 * \param b stores, as input, the RHS of the linear system of equation and stores, as output,
 *          the solution matrix X (Size x B_ColCount).
 * \param P vector of Size unsigned integer elements holding, as output, the permutations done
 *          the rows of matrix A during the decomposition to LU.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero.
 *
 * \throws singularity_error if the matrix A is numerically singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename IndexVector>
void linsolve_PLU(Matrix1& A, Matrix2& b, IndexVector& P,
                  mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_writable_matrix_v<Matrix1>);
  static_assert(is_writable_matrix_v<Matrix2> || is_writable_vector_v<Matrix2>);
  static_assert(is_writable_vector_v<IndexVector>);

  P.resize(A.get_col_count());

  auto Atmp = [&]() {
    if constexpr (is_fully_writable_matrix_v<Matrix1>) {
      return std::unique_ptr<Matrix1, null_deleter>(&A, null_deleter());
    } else {
      return std::make_unique<
          mat<mat_value_type_t<Matrix1>, mat_structure::square>>(A);
    }
  }();
  auto btmp = [&]() {
    if constexpr (is_fully_writable_matrix_v<Matrix2>) {
      return std::unique_ptr<Matrix2, null_deleter>(&b, null_deleter());
    } else if constexpr (is_writable_vector_v<Matrix2>) {
      return std::make_unique<
          mat_vect_adaptor<Matrix2, mat_alignment::column_major>>(b);
    } else {
      return std::make_unique<
          mat<mat_value_type_t<Matrix2>, mat_structure::rectangular>>(b);
    }
  }();

  if (Atmp->get_col_count() != Atmp->get_row_count()) {
    throw std::range_error(
        "PLU decomposition impossible! Matrix A is not square!");
  }
  if (btmp->get_row_count() != Atmp->get_col_count()) {
    throw std::range_error(
        "PLU decomposition impossible! Matrix b must have same row count as "
        "A!");
  }

  detail::linsolve_PLU_impl(*Atmp, *btmp, P, NumTol);

  if constexpr (!is_fully_writable_matrix_v<Matrix1>) {
    A = *Atmp;
  }
  if constexpr (!is_fully_writable_matrix_v<Matrix2>) {
    b = *btmp;
  }
}

/**
 * Solves the linear problem AX = B using PLU decomposition as defined by Crout`s method.
 * \note To solve a linear system involving vectors for X and B, use the Matrix-Vector Adaptors
 *(mat_vector_adaptor.hpp).
 *
 * \tparam Matrix1 A writable matrix type.
 * \tparam Matrix2 A writable matrix type.
 * \param A well-conditioned, square (Size x Size), real, full-rank matrix which multiplies x.
 *          As output, A stores the LU decomposition, permutated by P.
 * \param b stores, as input, the RHS of the linear system of equation and stores, as output,
 *          the solution matrix X (Size x B_ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero.
 *
 * \throws singularity_error if the matrix A is numerically singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename IndexVector>
void linsolve_PLU(Matrix1& A, Matrix2& b,
                  mat_value_type_t<Matrix1> NumTol = 1E-8) {
  vect_n<int> P;
  linsolve_PLU(A, b, P, NumTol);
}

/**
 * Functor to wrap a call to a PLU decomposition-based linear system solver.
 */
struct PLU_linsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B,
                  mat_value_type_t<Matrix1> NumTol = 1E-8) {
    using ValueType = mat_value_type_t<Matrix1>;
    mat<ValueType, mat_structure::square> A_tmp(A);
    X = B;
    linsolve_PLU(A_tmp, X, NumTol);
  }
};

/**
 * Inverts a matrix using PLU decomposition as defined by Crout`s method.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A writable matrix type.
 * \param A well-conditioned, square (Size x Size), real, full-rank matrix to be inverted.
 * \param A_inv The matrix which stores, as output, the inverse of A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero.
 *
 * \throws singularity_error if the matrix A is numerically singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not a square matrix.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
void invert_PLU(Matrix1 A, Matrix2& A_inv,
                mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_writable_matrix_v<Matrix1>);
  static_assert(is_writable_matrix_v<Matrix2>);
  using ValueType = mat_value_type_t<Matrix1>;

  A_inv = mat<ValueType, mat_structure::identity>(A.get_col_count());
  vect_n<int> P(A.get_col_count());
  linsolve_PLU(A, A_inv, P, NumTol);
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_GAUSSIAN_ELIM_H_

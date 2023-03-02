/**
 * \file mat_cholesky.hpp
 *
 * This library provides the Cholesky decomposition function template which can be used
 * to find the Cholesky factors of a symmetric matrix. The Cholesky decomposition will take
 * a symmetric, positive-definite matrix and find a lower-triangular matrix such that the
 * product of this factor and its transpose produce the original matrix. This method is
 * especially useful in contexts where well-conditioned positive-definite matrix are involved.
 * Following the decomposition, the obtained lower-triangular matrix factor can be used to
 * efficiently compute the inverse of the matrix or solve linear systems of equations represented
 * by the matrix and some right-hand-side. This library provides all those methods, overloaded
 * according to the types of matrices given to it (to do in-place algorithms when possible).
 *
 * According to performance tests, PLU methods are as good as Cholesky methods in terms of speed.
 * And they are both the best for well-conditioned matrices. For ill-conditioned matrices, QR-decomposition
 * methods are only a little slower then PLU/Cholesky (about 20% slower, same time-complexity) but provide better
 * numerical stability. The Jacobi methods are significantly slower, but this implementation is in need
 * of a revision for performance enhancement. And, of course, SVD is also very slow (slightly faster than
 * Jacobi) but it is based on a LAPACK implementation that is very poorly written, and it has not been
 * updated since.
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

#ifndef REAK_MAT_CHOLESKY_HPP
#define REAK_MAT_CHOLESKY_HPP

#include "mat_alg.hpp"
#include "mat_num_exceptions.hpp"

#include <type_traits>

namespace ReaK {

/*************************************************************************
                        Cholesky Decomposition
*************************************************************************/

namespace detail {

template <typename Matrix1, typename Matrix2>
void decompose_Cholesky_impl(const Matrix1& A, Matrix2& L,
                             mat_value_type_t<Matrix1> NumTol) {
  using std::sqrt;
  int N = A.get_row_count();
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < i; ++j) {
      L(i, j) = A(i, j);
      for (int k = 0; k < j; ++k) {
        L(i, j) -= L(i, k) * L(j, k);
      }
      L(i, j) /= L(j, j);
    }
    L(i, i) = A(i, i);
    for (int k = 0; k < i; ++k) {
      L(i, i) -= L(i, k) * L(i, k);
    }
    if (L(i, i) < NumTol) {
      throw singularity_error("A");
    }
    L(i, i) = sqrt(L(i, i));
  }
}

template <typename Matrix1>
void decompose_BandCholesky_impl(Matrix1& A, mat_size_type_t<Matrix1> p,
                                 mat_value_type_t<Matrix1> NumTol) {
  using std::sqrt;
  int N = A.get_row_count();
  for (int i = 0; i < N; ++i) {
    for (int j = (i > p ? i - p : 0); j < i; ++j) {
      int k = j + p;
      if (k >= N) {
        k = N - 1;
      }
      for (int l = i; l <= k; ++l) {
        A(l, i) -= A(i, j) * A(l, j);
      }
    }
    int k = i + p;
    if (k >= N) {
      k = N - 1;
    }
    if (A(i, i) < NumTol) {
      throw singularity_error("A");
    }
    A(i, i) = sqrt(A(i, i));
    for (int l = i + 1; l <= k; ++l) {
      A(l, i) /= A(i, i);
    }
  }
}

template <typename Matrix1>
void decompose_TriDiagLDL_impl(Matrix1& A, mat_value_type_t<Matrix1> NumTol) {
  using ValueType = mat_value_type_t<Matrix1>;
  using std::abs;
  int N = A.get_row_count();
  if (N == 0) {
    return;
  }
  if (abs(A(0, 0)) < NumTol) {
    throw singularity_error("A");
  }
  if (N > 1) {
    A(1, 0) /= A(0, 0);
  }
  for (int i = 1; i < N; ++i) {
    ValueType v_prev = A(i, i - 1) * A(i - 1, i - 1);
    ValueType v = A(i, i) - A(i, i - 1) * v_prev;
    if (abs(v) < NumTol) {
      throw singularity_error("A");
    }
    A(i, i) = v;
    if (i < N - 1) {
      A(i + 1, i) = (A(i + 1, i) - A(i + 1, i - 1) * v_prev) / v;
    }
  }
}

template <typename Matrix1>
void decompose_LDL_impl(Matrix1& A, mat_value_type_t<Matrix1> NumTol) {
  using ValueType = mat_value_type_t<Matrix1>;
  using std::abs;
  int N = A.get_row_count();
  vect_n<ValueType> v(N);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < i; ++j) {
      v[j] = A(i, j) * A(j, j);
    }
    v[i] = A(i, i);
    for (int j = 0; j < i; ++j) {
      v[i] -= A(i, j) * v[j];
    }
    A(i, i) = v[i];
    if (abs(v[i]) < NumTol) {
      throw singularity_error("A");
    }
    for (int j = i + 1; j < N; ++j) {
      for (int k = 0; k < i; ++k) {
        A(j, i) -= A(j, k) * v[k];
      }
      A(j, i) /= v[i];
    }
  }
}

template <typename Matrix1, typename Matrix2>
void backsub_Cholesky_impl(const Matrix1& L, Matrix2& B) {
  int N = L.get_row_count();
  int M = B.get_col_count();
  for (int j = 0; j < M; ++j) {  // for every column of B
    // Start solving L * Y = B
    for (int i = 0; i < N; ++i) {  // for every row of L
      for (int k = 0; k < i;
           ++k) {  // for every element of row i in L before the diagonal.
        B(i, j) -=
            L(i, k) *
            B(k, j);  // subtract to B(i,j) the product of L(i,k) * Y(k,j)
      }
      B(i, j) /=
          L(i, i);  // do Y(i,j) = (B(i,j) - sum_k(L(i,k) * Y(k,j))) / L(i,i)
    }
    // Then solve L.transpose() * X = Y
    for (int i = N; i > 0;
         --i) {  // for every row of L.transpose(), starting from the last
      for (
          int k = N; k > i;
          --k) {  // for every element of row i in L.tranpose(), after the diagonal.
        B(i - 1, j) -=
            L(k - 1, i - 1) *
            B(k - 1, j);  // subtract to B(i,j) the product of L(k,i) * Y(k,j)
      }
      B(i - 1, j) /= L(i - 1, i - 1);
    }
  }
}

template <typename Matrix1, typename Matrix2>
void backsub_BandCholesky_impl(const Matrix1& L, Matrix2& B,
                               mat_size_type_t<Matrix1> p) {
  int N = L.get_row_count();
  int M = B.get_col_count();
  for (int j = 0; j < M; ++j) {  // for every column of B
    // Start solving L * Y = B
    for (int i = 0; i < N; ++i) {  // for every row of L
      for (int k = (i > p ? i - p : 0); k < i;
           ++k) {  // for every element of row i in L before the diagonal.
        B(i, j) -=
            L(i, k) *
            B(k, j);  // subtract to B(i,j) the product of L(i,k) * Y(k,j)
      }
      B(i, j) /=
          L(i, i);  // do Y(i,j) = (B(i,j) - sum_k(L(i,k) * Y(k,j))) / L(i,i)
    }
    // Then solve transpose(L) * X = Y
    for (int i = N;
         i > 0;) {  // for every row of L.transpose(), starting from the last
      int l = --i + p;
      if (l >= N) {
        l = N - 1;
      }
      for (
          int k = l; k > i;
          --k) {  // for every element of row i in L.tranpose(), after the diagonal.
        B(i, j) -=
            L(k, i) *
            B(k, j);  // subtract to B(i,j) the product of L(k,i) * Y(k,j)
      }
      B(i, j) /= L(i, i);
    }
  }
}

template <typename Matrix1, typename Matrix2>
void backsub_LDL_impl(const Matrix1& L, Matrix2& B) {
  int N = L.get_row_count();
  int M = B.get_col_count();
  for (int j = 0; j < M; ++j) {  // for every column of B
    // Start solving L * Y = B
    for (int i = 0; i < N; ++i) {  // for every row of L
      for (int k = 0; k < i;
           ++k) {  // for every element of row i in L before the diagonal.
        B(i, j) -=
            L(i, k) *
            B(k, j);  // subtract to B(i,j) the product of L(i,k) * Y(k,j)
      }
    }
    // Solve D * X = Y;
    for (int i = 0; i < N; ++i) {
      B(i, j) /= L(i, i);
    }
    // Then solve transpose(L) * W = X
    for (int i = N;
         i > 0;) {  // for every row of L.transpose(), starting from the last
      --i;
      for (
          int k = N - 1; k > i;
          --k) {  // for every element of row i in L.tranpose(), after the diagonal.
        B(i, j) -=
            L(k, i) *
            B(k, j);  // subtract to B(i,j) the product of L(k,i) * Y(k,j)
      }
    }
  }
}

template <typename Matrix1, typename Matrix2>
void backsub_TriDiagLDL_impl(const Matrix1& L, Matrix2& B) {
  int N = L.get_row_count();
  int M = B.get_col_count();

  for (int j = 0; j < M; ++j) {  // for every column of B
    // Start solving L * Y = B
    for (int i = 1; i < N; ++i) {  // for every row of L
      // for every element of row i in L before the diagonal.
      B(i, j) -=
          L(i, i - 1) *
          B(i - 1, j);  // subtract to B(i,j) the product of L(i,k) * Y(k,j)
    }
    // Solve D * X = Y;
    for (int i = 0; i < N; ++i) {
      B(i, j) /= L(i, i);
    }
    // Then solve transpose(L) * W = X
    for (int i = N - 1; i > 0;
         --i) {  // for every row of L.transpose(), starting from the last
      // for every element of row i in L.tranpose(), after the diagonal.
      B(i - 1, j) -=
          L(i, i - 1) *
          B(i, j);  // subtract to B(i,j) the product of L(k,i) * Y(k,j)
    }
  }
}

}  // namespace detail

/**
 * Performs the Cholesky decomposition of A (positive-definite symmetric matrix) (Cholesky-Crout Algorithm).
 * returns a Lower-triangular matrix such that A = L * transpose(L).
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param L stores, as output, the lower-triangular matrix in A = L * transpose(L).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient) or not positive-definite.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility to ensure it's correct.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
void decompose_Cholesky(const Matrix1& A, Matrix2& L,
                        mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix1>);
  static_assert(is_writable_matrix_v<Matrix2>);
  if constexpr (is_diagonal_matrix_v<Matrix1>) {
    using std::sqrt;
    L = A;
    for (int i = 0; i < A.get_row_count(); ++i) {
      L(i, i) = sqrt(L(i, i));
    }
  } else if constexpr (is_fully_writable_matrix_v<Matrix2>) {
    L.set_col_count(A.get_col_count());
    detail::decompose_Cholesky_impl(A, L, NumTol);
  } else {
    using ValueType = mat_value_type_t<Matrix1>;
    mat<ValueType, mat_structure::square> L_tmp(A.get_row_count(),
                                                ValueType(0));
    detail::decompose_Cholesky_impl(A, L_tmp, NumTol);
    L = L_tmp;
  }
}

/**
 * Performs the Cholesky decomposition of A (positive-definite symmetric matrix) (Cholesky-Crout Algorithm),
 * and uses the result to compute the determinant of A.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 * \return the determinant of A, or 0 if it was judged to be singular or not positive-definite.
 *
 * \note the symmetry of the matrix A is not checked and thus it is
 *       the caller's responsibility to ensure it's correct.
 *
 * \author Mikael Persson
 */
template <typename Matrix>
mat_value_type_t<Matrix> determinant_Cholesky(
    const Matrix& A, mat_value_type_t<Matrix> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  mat<ValueType, mat_structure::square> L(A.get_row_count(), ValueType(0));
  try {
    decompose_Cholesky(A, L, NumTol);
  } catch (singularity_error& e) {
    return ValueType(0);
  }
  ValueType prod = ValueType(1);
  for (int i = 0; i < L.get_col_count(); ++i) {
    prod *= L(i, i);
  }
  return prod * prod;
}

/**
 * Solves the linear problem AX = B using Cholesky decomposition.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param b stores, as input, the RHS of the linear system of equations and stores, as output,
 *          the solution matrix X (Size x B_ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility.
 * \note if you wish to apply this method with a vector on the RHS, then use ReaK::mat_vect_adaptor (and related
 *classes).
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
void linsolve_Cholesky(const Matrix1& A, Matrix2& b,
                       mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix1>);
  static_assert(is_writable_matrix_v<Matrix2>);
  if (b.get_row_count() != A.get_col_count()) {
    throw std::range_error(
        "For linear equation solution, matrix b must have same row count as "
        "A!");
  }
  using ValueType = mat_value_type_t<Matrix1>;

  if constexpr (mat_traits<Matrix1>::structure == mat_structure::identity) {
    // nothing to do.
  } else if constexpr (is_diagonal_matrix_v<Matrix1>) {
    if constexpr (is_diagonal_matrix_v<Matrix2>) {
      for (int i = 0; i < A.get_row_count(); ++i) {
        if (A(i, i) < NumTol) {
          throw singularity_error("A");
        }
        b(i, i) /= A(i, i);
      }
    } else if constexpr (is_fully_writable_matrix_v<Matrix2>) {
      for (int i = 0; i < A.get_row_count(); ++i) {
        if (A(i, i) < NumTol) {
          throw singularity_error("A");
        }
        for (int j = 0; j < b.get_col_count(); ++j) {
          b(i, j) /= A(i, i);
        }
      }
    } else {
      mat<ValueType, mat_structure::rectangular> b_tmp(b);
      for (int i = 0; i < A.get_row_count(); ++i) {
        if (A(i, i) < NumTol) {
          throw singularity_error("A");
        }
        for (int j = 0; j < b.get_col_count(); ++j) {
          b_tmp(i, j) /= A(i, i);
        }
      }
      b = b_tmp;
    }

  } else if constexpr (is_fully_writable_matrix_v<Matrix2>) {
    mat<ValueType, mat_structure::square> L(A.get_row_count(), ValueType(0));
    detail::decompose_Cholesky_impl(A, L, NumTol);
    detail::backsub_Cholesky_impl(L, b);
  } else {
    mat<ValueType, mat_structure::square> L(A.get_row_count(), ValueType(0));
    detail::decompose_Cholesky_impl(A, L, NumTol);
    mat<mat_value_type_t<Matrix2>, mat_structure::rectangular> b_tmp(b);
    detail::backsub_Cholesky_impl(L, b_tmp);
    b = b_tmp;
  }
}

/**
 * Functor to wrap a call to a Cholesky-decomposition-based linear system solver.
 */
struct Cholesky_linsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B,
                  mat_value_type_t<Matrix1> NumTol = 1E-8) {
    X = B;
    linsolve_Cholesky(A, X, NumTol);
  }
};

/**
 * Inverts a matrix using Cholesky decomposition.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be inverted.
 * \param A_inv stores, as output, the inverse of matrix A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
void invert_Cholesky(const Matrix1& A, Matrix2& A_inv,
                     mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix1>);
  static_assert(is_writable_matrix_v<Matrix2>);
  using ValueType = mat_value_type_t<Matrix1>;

  mat<ValueType, mat_structure::square> L;
  decompose_Cholesky(A, L, NumTol);
  mat<ValueType, mat_structure::square> result(
      mat<ValueType, mat_structure::identity>(A.get_col_count()));
  detail::backsub_Cholesky_impl(L, result);
  A_inv = result;
}

/**
 * Performs the Band-Cholesky decomposition of A (positive-definite band-symmetric matrix).
 * returns a Lower-triangular matrix such that A = L * transpose(L).
 *
 * \param A real, positive-definite, band-symmetric, square, full-rank matrix to be decomposed.
 * \param L stores, as output, the lower-triangular matrix in A = L * transpose(L).
 * \param bandwidth the band-width of the diagonal band of the matrix (i.e., a tridiagonal matrix has bandwidth = 1).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient) or not positive-definite.
 *
 * \note the positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility to ensure it's correct.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
void decompose_BandCholesky(const Matrix1& A, Matrix2& L,
                            mat_size_type_t<Matrix1> bandwidth = 1,
                            mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix1>);
  static_assert(is_square_matrix_v<Matrix1>);
  static_assert(is_fully_writable_matrix_v<Matrix2>);
  using ValueType = mat_value_type_t<Matrix1>;
  L = A;
  detail::decompose_BandCholesky_impl(L, bandwidth, NumTol);
  for (int i = 1; i < L.get_col_count(); ++i) {
    for (int j = 0; j < i; ++j) {
      L(j, i) = ValueType(0.0);
    }
  }
}

/**
 * Performs the Band-Cholesky decomposition of A (positive-definite symmetric matrix),
 * and uses the result to compute the determinant of A.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 * \return the determinant of A, or 0 if it was judged to be singular or not positive-definite.
 *
 * \note the symmetry of the matrix A is not checked and thus it is
 *       the caller's responsibility to ensure it's correct.
 *
 * \author Mikael Persson
 */
template <typename Matrix>
mat_value_type_t<Matrix> determinant_BandCholesky(
    const Matrix& A, mat_size_type_t<Matrix> bandwidth,
    mat_value_type_t<Matrix> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  mat<ValueType, mat_structure::square> L(A);
  try {
    detail::decompose_BandCholesky_impl(L, bandwidth, NumTol);
  } catch (singularity_error& e) {
    return ValueType(0);
  }
  ValueType prod = ValueType(1);
  for (int i = 0; i < L.get_col_count(); ++i) {
    prod *= L(i, i);
  }
  return prod * prod;
}

/**
 * Solves the linear problem AX = B using Band-Cholesky decomposition.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param b stores, as input, the RHS of the linear system of equations and stores, as output,
 *          the solution matrix X (Size x B_ColCount).
 * \param bandwidth the band-width of the diagonal band of the matrix (i.e., a tridiagonal matrix has bandwidth = 1).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility.
 * \note if you wish to apply this method with a vector on the RHS, then use ReaK::mat_vect_adaptor (and related
 *classes).
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
void linsolve_BandCholesky(const Matrix1& A, Matrix2& b,
                           mat_size_type_t<Matrix1> bandwidth,
                           mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix1>);
  static_assert(is_square_matrix_v<Matrix1>);
  static_assert(is_fully_writable_matrix_v<Matrix2>);
  if (b.get_row_count() != A.get_col_count()) {
    throw std::range_error(
        "For linear equation solution, matrix b must have same row count as "
        "A!");
  }

  using ValueType = mat_value_type_t<Matrix1>;
  mat<ValueType, mat_structure::square> L(A);
  detail::decompose_BandCholesky_impl(L, bandwidth, NumTol);
  detail::backsub_BandCholesky_impl(L, b, bandwidth);
}

/**
 * Functor to wrap a call to a Cholesky-decomposition-based linear system solver.
 */
struct BandCholesky_linsolver {
  std::size_t band_width;
  explicit BandCholesky_linsolver(std::size_t aBandWidth = 1)
      : band_width(aBandWidth) {}

  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B,
                  mat_value_type_t<Matrix1> NumTol = 1E-8) {
    X = B;
    linsolve_BandCholesky(A, X, band_width, NumTol);
  }
};

/**
 * Inverts a matrix using Band-Cholesky decomposition.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be inverted.
 * \param A_inv stores, as output, the inverse of matrix A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
void invert_BandCholesky(const Matrix1& A, Matrix2& A_inv,
                         mat_size_type_t<Matrix1> bandwidth,
                         mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix1>);
  static_assert(is_writable_matrix_v<Matrix2>);
  using ValueType = mat_value_type_t<Matrix1>;
  mat<ValueType, mat_structure::square> result(
      mat<ValueType, mat_structure::identity>(A.get_col_count()));
  linsolve_BandCholesky(A, result, bandwidth, NumTol);
  A_inv = result;
}

/**
 * Performs the LDL decomposition of A (symmetric matrix).
 * returns a lower-triangular and diagonal matrix such that A = L * D * transpose(L).
 *
 * \param A real, band-symmetric, square, full-rank matrix to be decomposed.
 * \param L stores, as output, the lower-triangular matrix in A = L * D * transpose(L).
 * \param D stores, as output, the diagonal matrix in A = L * D * transpose(L).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient) or not positive-definite.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3>
void decompose_LDL(const Matrix1& A, Matrix2& L, Matrix3& D,
                   mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix1>);
  static_assert(is_square_matrix_v<Matrix1>);
  static_assert(is_fully_writable_matrix_v<Matrix2>);
  static_assert(is_writable_matrix_v<Matrix3>);
  using ValueType = mat_value_type_t<Matrix1>;
  L = A;
  detail::decompose_LDL_impl(L, NumTol);
  for (int i = 0; i < L.get_col_count(); ++i) {
    D(i, i) = L(i, i);
    L(i, i) = ValueType(1.0);
    for (int j = 0; j < i; ++j) {
      L(j, i) = ValueType(0.0);
    }
  }
}

/**
 * Performs the LDL decomposition of A (symmetric matrix),
 * and uses the result to compute the determinant of A.
 *
 * \param A real, symmetric, square, full-rank matrix to be decomposed.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 * \return the determinant of A, or 0 if it was judged to be singular or not positive-definite.
 *
 * \author Mikael Persson
 */
template <typename Matrix>
mat_value_type_t<Matrix> determinant_LDL(
    const Matrix& A, mat_value_type_t<Matrix> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  mat<ValueType, mat_structure::square> L(A);
  try {
    detail::decompose_LDL_impl(L, NumTol);
  } catch (singularity_error& e) {
    return ValueType(0);
  }
  ValueType prod = ValueType(1.0);
  for (int i = 0; i < L.get_col_count(); ++i) {
    prod *= L(i, i);
  }
  return prod;  // note: no squaring here since the diagonals are not square-roots (as in Cholesky).
}

/**
 * Solves the linear problem AX = B using LDL decomposition.
 *
 * \param A real, symmetric, square, full-rank matrix to be decomposed.
 * \param b stores, as input, the RHS of the linear system of equations and stores, as output,
 *          the solution matrix X (Size x B_ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \note if you wish to apply this method with a vector on the RHS, then use ReaK::mat_vect_adaptor (and related
 *classes).
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
void linsolve_LDL(const Matrix1& A, Matrix2& b,
                  mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix1>);
  static_assert(is_square_matrix_v<Matrix1>);
  static_assert(is_fully_writable_matrix_v<Matrix2>);
  if (b.get_row_count() != A.get_col_count()) {
    throw std::range_error(
        "For linear equation solution, matrix b must have same row count as "
        "A!");
  }

  using ValueType = mat_value_type_t<Matrix1>;
  mat<ValueType, mat_structure::square> L(A);
  detail::decompose_LDL_impl(L, NumTol);
  detail::backsub_LDL_impl(L, b);
}

/**
 * Functor to wrap a call to a Cholesky-decomposition-based linear system solver.
 */
struct LDL_linsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B,
                  mat_value_type_t<Matrix1> NumTol = 1E-8) {
    X = B;
    linsolve_LDL(A, X, NumTol);
  }
};

/**
 * Inverts a matrix using LDL decomposition.
 *
 * \param A real, symmetric, square, full-rank matrix to be inverted.
 * \param A_inv stores, as output, the inverse of matrix A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
void invert_LDL(const Matrix1& A, Matrix2& A_inv,
                mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix1>);
  static_assert(is_writable_matrix_v<Matrix2>);
  using ValueType = mat_value_type_t<Matrix1>;
  mat<ValueType, mat_structure::square> result(
      mat<ValueType, mat_structure::identity>(A.get_col_count()));
  linsolve_LDL(A, result, NumTol);
  A_inv = result;
}

/**
 * Performs the Tridiagonal-LDL decomposition of A (symmetric matrix).
 * returns a lower-triangular and diagonal matrix such that A = L * D * transpose(L).
 * \note This method is currently not working (bug).
 * \param A real, tridiagonal, square, full-rank matrix to be decomposed.
 * \param L stores, as output, the lower-triangular matrix in A = L * D * transpose(L).
 * \param D stores, as output, the diagonal matrix in A = L * D * transpose(L).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient) or not positive-definite.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3>
void decompose_TriDiagLDL(const Matrix1& A, Matrix2& L, Matrix3& D,
                          mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix1>);
  static_assert(is_square_matrix_v<Matrix1>);
  static_assert(is_fully_writable_matrix_v<Matrix2>);
  static_assert(is_writable_matrix_v<Matrix3>);
  using ValueType = mat_value_type_t<Matrix1>;
  L = A;
  detail::decompose_TriDiagLDL_impl(L, NumTol);
  for (int i = 0; i < L.get_col_count(); ++i) {
    D(i, i) = L(i, i);
    L(i, i) = ValueType(1.0);
    for (int j = 0; j < i; ++j) {
      L(j, i) = ValueType(0.0);
    }
  }
}

/**
 * Performs the Tridiagonal-LDL decomposition of A (symmetric matrix),
 * and uses the result to compute the determinant of A.
 *
 * \param A real, tridiagonal, square, full-rank matrix to be decomposed.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 * \return the determinant of A, or 0 if it was judged to be singular or not positive-definite.
 *
 * \author Mikael Persson
 */
template <typename Matrix>
mat_value_type_t<Matrix> determinant_TriDiagLDL(
    const Matrix& A, mat_value_type_t<Matrix> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  mat<ValueType, mat_structure::square> L(A);
  try {
    detail::decompose_TriDiagLDL_impl(L, NumTol);
  } catch (singularity_error& e) {
    return ValueType(0);
  }
  ValueType prod = ValueType(1.0);
  for (int i = 0; i < L.get_col_count(); ++i) {
    prod *= L(i, i);
  }
  return prod;  // note: no squaring here since the diagonals are not square-roots (as in Cholesky).
}

/**
 * Solves the linear problem AX = B using Tridiagonal-LDL decomposition.
 *
 * \param A real, tridiagonal, square, full-rank matrix to be decomposed.
 * \param b stores, as input, the RHS of the linear system of equations and stores, as output,
 *          the solution matrix X (Size x B_ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \note if you wish to apply this method with a vector on the RHS, then use ReaK::mat_vect_adaptor (and related
 *classes).
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
void linsolve_TriDiagLDL(const Matrix1& A, Matrix2& b,
                         mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix1>);
  static_assert(is_square_matrix_v<Matrix1>);
  static_assert(is_fully_writable_matrix_v<Matrix2>);
  if (b.get_row_count() != A.get_col_count()) {
    throw std::range_error(
        "For linear equation solution, matrix b must have same row count as "
        "A!");
  }

  using ValueType = mat_value_type_t<Matrix1>;
  mat<ValueType, mat_structure::square> L(A);
  detail::decompose_TriDiagLDL_impl(L, NumTol);
  detail::backsub_TriDiagLDL_impl(L, b);
}

/**
 * Functor to wrap a call to a Tridiagonal-LDL decomposition-based linear system solver.
 */
struct TriDiagLDL_linsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B,
                  mat_value_type_t<Matrix1> NumTol = 1E-8) {
    X = B;
    linsolve_TriDiagLDL(A, X, NumTol);
  }
};

/**
 * Inverts a matrix using Tridiagonal-LDL decomposition.
 *
 * \param A real, tridiagonal, square, full-rank matrix to be inverted.
 * \param A_inv stores, as output, the inverse of matrix A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
void invert_TriDiagLDL(const Matrix1& A, Matrix2& A_inv,
                       mat_value_type_t<Matrix1> NumTol = 1E-8) {
  static_assert(is_readable_matrix_v<Matrix1>);
  static_assert(is_writable_matrix_v<Matrix2>);
  using ValueType = mat_value_type_t<Matrix1>;
  mat<ValueType, mat_structure::square> result(
      mat<ValueType, mat_structure::identity>(A.get_col_count()));
  linsolve_TriDiagLDL(A, result, NumTol);
  A_inv = result;
}

}  // namespace ReaK

#endif

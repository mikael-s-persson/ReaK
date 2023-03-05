/**
 * \file mat_alg_general_hpp
 *
 * This library implements the general versions of many meta-functions (templates),
 * functions, and operators. These are meant to be used when no more-specialized
 * implementations exist for the matrix types involved.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date april 2011 (originally february 2010)
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

#ifndef REAK_MAT_OPERATORS_HPP
#define REAK_MAT_OPERATORS_HPP

#include "ReaK/core/base/defs.hpp"

//#include "ReaK/math/lin_alg/vect_alg.hpp"
#include "ReaK/math/lin_alg/mat_concepts.hpp"
#include "ReaK/math/lin_alg/mat_traits.hpp"
#include "ReaK/math/lin_alg/vect_concepts.hpp"

#include "ReaK/math/lin_alg/mat_op_results.hpp"

#include <type_traits>

#include <iomanip>

namespace ReaK {

/**
 * Prints a matrix to a standard output stream (<<) as "((a11; a12); (a21; a22))". \test PASSED
 */
template <typename Matrix>
std::enable_if_t<is_readable_matrix_v<Matrix>, std::ostream&> operator<<(
    std::ostream& out_stream, const Matrix& M) {
  out_stream << "(\n";
  if ((M.get_row_count() != 0) && (M.get_col_count() != 0)) {
    for (int i = 0; i < M.get_row_count(); ++i) {
      out_stream << "(" << std::setw(16) << M(i, 0);
      for (int j = 1; j < M.get_col_count(); ++j) {
        out_stream << "; " << std::setw(16) << M(i, j);
      }
      out_stream << ")";
      if (i != M.get_row_count() - 1) {
        out_stream << ";\n";
      } else {
        out_stream << "\n";
      }
    }
  }
  return (out_stream << ")");
}

/*******************************************************************************
                         Multiplication Operators
*******************************************************************************/

namespace detail {

template <typename Matrix1, typename Matrix2, typename ResultMatrix>
void dense_mat_multiply_impl(const Matrix1& M1, const Matrix2& M2,
                             ResultMatrix& MR) {
  using ValueType = mat_value_type_t<ResultMatrix>;
  for (int i = 0; i < M1.get_row_count(); ++i) {
    for (int jj = 0; jj < M2.get_col_count(); ++jj) {
      MR(i, jj) = ValueType(0.0);
      for (int j = 0; j < M1.get_col_count(); ++j) {
        MR(i, jj) += M1(i, j) * M2(j, jj);
      }
    }
  }
}

template <typename Matrix1, typename MatrixDiag, typename ResultMatrix>
void dense_diag_mat_multiply_impl(const Matrix1& M1, const MatrixDiag& M2,
                                  ResultMatrix& MR) {
  for (int i = 0; i < M1.get_row_count(); ++i) {
    for (int j = 0; j < M1.get_col_count(); ++j) {
      MR(i, j) = M1(i, j) * M2(j, j);
    }
  }
}

template <typename MatrixDiag, typename Matrix2, typename ResultMatrix>
void diag_dense_mat_multiply_impl(const MatrixDiag& M1, const Matrix2& M2,
                                  ResultMatrix& MR) {
  for (int i = 0; i < M2.get_row_count(); ++i) {
    for (int j = 0; j < M2.get_col_count(); ++j) {
      MR(i, j) = M1(i, i) * M2(i, j);
    }
  }
}

template <typename MatrixDiag1, typename MatrixDiag2, typename ResultMatrix>
void diag_diag_mat_multiply_impl(const MatrixDiag1& M1, const MatrixDiag2& M2,
                                 ResultMatrix& MR) {
  for (int i = 0; i < M1.get_row_count(); ++i) {
    MR(i, i) = M1(i, i) * M2(i, i);
  }
}

/* Multiplies two lower-triangular matrices and saves the result in-place into M2. */
template <typename MatrixLower1, typename MatrixLower2>
void inplace_lower_multiply_impl(const MatrixLower1& M1, MatrixLower2& M2) {
  using ValueType = mat_value_type_t<MatrixLower2>;

  for (int i = M1.get_row_count(); i > 0;) {
    --i;
    for (int j = 0; j <= i; ++j) {
      ValueType sum = ValueType(0.0);
      for (int k = j; k <= i; ++k) {
        sum += M1(i, k) * M2(k, j);
      }
      M2(i, j) = sum;
    }
  }
}

/* Multiplies two lower-triangular matrices and saves the result in-place into M2, and fills the upper part of M2 with
 * zeros. */
template <typename MatrixLower1, typename MatrixLower2>
void inplace_lower_multiply_with_fill_impl(const MatrixLower1& M1,
                                           MatrixLower2& M2) {
  using ValueType = mat_value_type_t<MatrixLower2>;

  for (int i = M1.get_row_count(); i > 0;) {
    --i;
    for (int j = 0; j <= i; ++j) {
      auto sum = ValueType(0.0);
      for (int k = j; k <= i; ++k) {
        sum += M1(i, k) * M2(k, j);
      }
      M2(i, j) = sum;
    }
    for (int j = i + 1; j < M2.get_col_count(); ++j) {
      M2(i, j) = ValueType(0.0);
    }
  }
}

}  // namespace detail

/**
 * General multiplication operator for any type of matrices. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M1 first matrix (first operand).
 * \param M2 second matrix (second operand).
 * \return General column-major matrix.
 * \throw std::range_error if the two matrix dimensions do not fit together.
 * \test PASSED
 */
template <typename Matrix1, typename Matrix2>
std::enable_if_t<is_readable_matrix_v<Matrix1> && is_readable_matrix_v<Matrix2>,
                 mat_product_result_t<Matrix1, Matrix2>>
operator*(const Matrix1& M1, const Matrix2& M2) {
  using result_type = mat_product_result_t<Matrix1, Matrix2>;
  if (M1.get_col_count() != M2.get_row_count()) {
    throw std::range_error("Matrix dimension mismatch.");
  }

  constexpr auto priority1 = mat_product_priority_v<Matrix1>;
  constexpr auto priority2 = mat_product_priority_v<Matrix2>;
  if constexpr ((priority1 == detail::product_priority_v<mat_structure::nil>) ||
                (priority2 == detail::product_priority_v<mat_structure::nil>)) {
    // nil and any
    return result_type{M1.get_row_count(), M2.get_col_count()};
  } else if constexpr ((priority2 ==
                        detail::product_priority_v<mat_structure::identity>)) {
    // any and identity
    return result_type{M1};
  } else if constexpr ((priority1 ==
                        detail::product_priority_v<mat_structure::identity>)) {
    // identity and any
    return result_type{M2};
  } else if constexpr ((priority2 ==
                        detail::product_priority_v<mat_structure::scalar>)) {
    // any and scalar
    result_type result(M1);
    if (M2.get_row_count() > 0) {
      result *= M2(0, 0);
    }
    return result;
  } else if constexpr ((priority1 ==
                        detail::product_priority_v<mat_structure::scalar>)) {
    // scalar and any
    result_type result(M2);
    if (M1.get_row_count() > 0) {
      result *= M1(0, 0);
    }
    return result;
  } else if constexpr (
      (priority1 == detail::product_priority_v<mat_structure::diagonal>)&&(
          priority2 == detail::product_priority_v<mat_structure::diagonal>)) {
    // diag and diag
    result_type result(M2);
    detail::diag_diag_mat_multiply_impl(M1, M2, result);
    return result;
  } else if constexpr (priority1 ==
                       detail::product_priority_v<mat_structure::diagonal>) {
    // diag and dense
    result_type result(M2);
    detail::diag_dense_mat_multiply_impl(M1, M2, result);
    return result;
  } else if constexpr (priority2 ==
                       detail::product_priority_v<mat_structure::diagonal>) {
    // dense and diag
    result_type result(M1);
    detail::dense_diag_mat_multiply_impl(M1, M2, result);
    return result;
  } else if constexpr (
      ((priority1 == detail::product_priority_v<mat_structure::symmetric>)&&(
          priority2 == detail::product_priority_v<mat_structure::symmetric>)) ||
      ((priority1 ==
        detail::product_priority_v<
            mat_structure::skew_symmetric>)&&(priority2 ==
                                              detail::product_priority_v<
                                                  mat_structure::
                                                      skew_symmetric>))) {
    // sym/skew and sym/skew
    return M1.multiply_with_same_mat(M2);
  } else if constexpr (priority1 == detail::product_priority_v<
                                        mat_structure::symmetric> ||
                       priority1 == detail::product_priority_v<
                                        mat_structure::skew_symmetric>) {
    // sym/skew and dense
    return M1.multiply_this_and_dense_mat(M2);
  } else if constexpr (priority2 == detail::product_priority_v<
                                        mat_structure::symmetric> ||
                       priority2 == detail::product_priority_v<
                                        mat_structure::skew_symmetric>) {
    // dense and sym/skew
    return M2.multiply_dense_and_this_mat(M1);
  } else if constexpr (is_square_matrix_v<Matrix2>) {
    // square/rect and square
    result_type result(M1);
    detail::dense_mat_multiply_impl(M1, M2, result);
    return result;
  } else if constexpr (is_square_matrix_v<Matrix1>) {
    // square and rect
    result_type result(M2);
    detail::dense_mat_multiply_impl(M1, M2, result);
    return result;
  } else {
    // rect and rect
    static_assert(is_resizable_matrix_v<result_type>);
    result_type result(M1);
    result.set_col_count(M2.get_col_count());
    detail::dense_mat_multiply_impl(M1, M2, result);
    return result;
  }
}

/**
 * Matrix multiplication operator with a scalar. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M The matrix (first operand).
 * \param S The scalar (second operand).
 * \return Column-major matrix equal to M * S.
 * \test PASSED
 */
template <typename Matrix, typename Scalar>
std::enable_if_t<is_writable_matrix_v<Matrix> &&
                     !is_readable_matrix_v<Scalar> &&
                     !is_readable_vector_v<Scalar>,
                 Matrix>
operator*(Matrix M, const Scalar& S) {
  M *= S;
  return M;
}

/**
 * Matrix multiplication operator with a scalar. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param S The scalar (first operand).
 * \param M The matrix (second operand).
 * \return Column-major matrix equal to S * M.
 * \test PASSED
 */
template <typename Matrix, typename Scalar>
std::enable_if_t<is_writable_matrix_v<Matrix> &&
                     !is_readable_matrix_v<Scalar> &&
                     !is_readable_vector_v<Scalar>,
                 Matrix>
operator*(const Scalar& S, Matrix M) {
  M *= S;
  return M;
}

/**
 * Matrix multiplication operator with a column vector. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M The matrix (first operand).
 * \param V The column vector (second operand).
 * \return Column vector equal to M * V.
 * \throw std::range_error if the matrix column count does not correspond to the vector dimension.
 * \test PASSED
 */
template <typename Matrix, typename Vector>
std::enable_if_t<is_readable_matrix_v<Matrix> && is_readable_vector_v<Vector>,
                 vect_copy_t<Vector>>
operator*(const Matrix& M, const Vector& V) {
  using result_type = vect_copy_t<Vector>;
  if (V.size() != M.get_col_count()) {
    throw std::range_error("Matrix dimension mismatch.");
  }
  using ValueType = mat_value_type_t<Matrix>;
  result_type result;
  if constexpr (is_resizable_vector_v<result_type>) {
    result.resize(M.get_row_count());
  } else {
    if (V.size() != M.get_row_count()) {
      throw std::range_error("Matrix dimension mismatch.");
    }
  }

  constexpr auto priority = mat_product_priority_v<Matrix>;
  if constexpr (priority == detail::product_priority_v<mat_structure::nil>) {
    // Leave result as zero vector.
  } else if constexpr (priority ==
                       detail::product_priority_v<mat_structure::identity>) {
    result = V;
  } else if constexpr (priority ==
                       detail::product_priority_v<mat_structure::scalar>) {
    if (V.size() > 0) {
      result = M(0, 0) * V;
    }
  } else if constexpr (
      priority == detail::product_priority_v<mat_structure::symmetric> ||
      priority == detail::product_priority_v<mat_structure::skew_symmetric> ||
      priority == detail::product_priority_v<mat_structure::diagonal>) {
    M.multiply_with_vector_rhs(V, result);
  } else {
    for (int i = 0; i < M.get_row_count(); ++i) {
      result[i] = ValueType(0.0);
      for (int j = 0; j < M.get_col_count(); ++j) {
        result[i] += M(i, j) * V[j];
      }
    }
  }
  return result;
}

/**
 * Matrix multiplication operator with a row vector. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param V The row vector (first operand).
 * \param M The matrix (second operand).
 * \return Row vector equal to V * M.
 * \throw std::range_error if the matrix row count does not correspond to the vector dimension.
 * \test PASSED
 */
template <typename Matrix, typename Vector>
std::enable_if_t<is_readable_matrix_v<Matrix> && is_readable_vector_v<Vector>,
                 vect_copy_t<Vector>>
operator*(const Vector& V, const Matrix& M) {
  using result_type = vect_copy_t<Vector>;
  if (V.size() != M.get_row_count()) {
    throw std::range_error("Matrix dimension mismatch.");
  }
  using ValueType = mat_value_type_t<Matrix>;
  result_type result;
  if constexpr (is_resizable_vector_v<result_type>) {
    result.resize(M.get_col_count());
  } else {
    if (V.size() != M.get_col_count()) {
      throw std::range_error("Matrix dimension mismatch.");
    }
  }

  constexpr auto priority = mat_product_priority_v<Matrix>;
  if constexpr (priority == detail::product_priority_v<mat_structure::nil>) {
    // Leave result as zero vector.
  } else if constexpr (priority ==
                       detail::product_priority_v<mat_structure::identity>) {
    result = V;
  } else if constexpr (priority ==
                       detail::product_priority_v<mat_structure::scalar>) {
    if (V.size() > 0) {
      result = M(0, 0) * V;
    }
  } else if constexpr (
      priority == detail::product_priority_v<mat_structure::symmetric> ||
      priority == detail::product_priority_v<mat_structure::skew_symmetric> ||
      priority == detail::product_priority_v<mat_structure::diagonal>) {
    M.multiply_with_vector_lhs(V, result);
  } else {
    for (int j = 0; j < M.get_col_count(); ++j) {
      result[j] = ValueType(0.0);
      for (int i = 0; i < M.get_row_count(); ++i) {
        result[j] += M(i, j) * V[i];
      }
    }
  }
  return result;
}

/*******************************************************************************
                         Addition / Subtraction Operators
*******************************************************************************/

/**
 * General (least-specialized) transpose function for any type of matrix.
 * \param M The matrix to be transposed.
 * \return The transpose of M.
 */
template <typename Matrix>
std::enable_if_t<is_readable_matrix_v<Matrix>,
                 mat_addition_result_t<Matrix, Matrix>>
transpose(const Matrix& M) {
  using result_type = mat_addition_result_t<Matrix, Matrix>;
  result_type result(M);
  for (int j = 0; j < result.get_col_count(); ++j) {
    for (int i = 0; i < result.get_row_count(); ++i) {
      result(i, j) = M(j, i);
    }
  }
  return result;
}

/**
 * General (least-specialized) transpose function for any type of matrix.
 * \param M The matrix to be transposed.
 * \return The transpose of M.
 */
template <typename Matrix>
std::enable_if_t<is_readable_matrix_v<Matrix>,
                 mat_addition_result_t<Matrix, Matrix>>
transpose_move(const Matrix& M) {
  return transpose(M);
}

/**
 * General (least-specialized) unary-negation operator for any type of matrices.
 * This is a default operator that will be called if no better special-purpose overload exists.
 * \param M first matrix (first operand).
 * \return General column-major matrix.
 * \throw std::range_error if the two matrix dimensions do not fit together.
 * \test PASSED
 */
template <typename Matrix>
std::enable_if_t<is_readable_matrix_v<Matrix>,
                 mat_addition_result_t<Matrix, Matrix>>
operator-(const Matrix& M) {
  using result_type = mat_addition_result_t<Matrix, Matrix>;
  result_type result(M);
  for (int j = 0; j < M.get_col_count(); ++j) {
    for (int i = 0; i < M.get_row_count(); ++i) {
      // this needs to be taking the original matrix (avoid cross-referenced elements).
      result(i, j) = -M(i, j);
    }
  }
  return result;
}

/**
 * General (least-specialized) addition operator for any type of matrices. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M1 first matrix (first operand).
 * \param M2 second matrix (second operand).
 * \return General column-major matrix.
 * \throw std::range_error if the two matrix dimensions do not fit together.
 * \test PASSED
 */
template <typename Matrix1, typename Matrix2>
std::enable_if_t<is_readable_matrix_v<Matrix1> && is_readable_matrix_v<Matrix2>,
                 mat_addition_result_t<Matrix1, Matrix2>>
operator+(const Matrix1& M1, const Matrix2& M2) {
  if ((M1.get_row_count() != M2.get_row_count()) ||
      (M1.get_col_count() != M2.get_col_count())) {
    throw std::range_error("Matrix dimension mismatch.");
  }
  using result_type = mat_addition_result_t<Matrix1, Matrix2>;

  constexpr auto priority1 = mat_addition_priority_v<Matrix1>;
  constexpr auto priority2 = mat_addition_priority_v<Matrix2>;
  if constexpr (priority1 == detail::addition_priority_v<mat_structure::nil>) {
    // nil + any
    return result_type(M2);
  } else if constexpr (priority2 ==
                       detail::addition_priority_v<mat_structure::nil>) {
    // any + nil
    return result_type(M1);
  } else if constexpr (
      (priority1 == detail::addition_priority_v<mat_structure::identity> ||
       priority1 ==
           detail::addition_priority_v<
               mat_structure::scalar>)&&(priority2 ==
                                             detail::addition_priority_v<
                                                 mat_structure::identity> ||
                                         priority2 ==
                                             detail::addition_priority_v<
                                                 mat_structure::scalar>)) {
    // identity/scalar + identity/scalar
    result_type result(M1.get_row_count());
    if (M1.get_row_count() > 0) {
      result(0, 0) = M1(0, 0) + M2(0, 0);
    }
    return result;
  } else if constexpr (
      (priority1 == detail::addition_priority_v<mat_structure::identity> ||
       priority1 ==
           detail::addition_priority_v<
               mat_structure::scalar>)&&(priority2 <=
                                         detail::addition_priority_v<
                                             mat_structure::diagonal>)) {
    // identity/scalar + dense
    result_type result(M2);
    for (int i = 0; i < result.get_row_count(); ++i) {
      result(i, i) += M1(0, 0);
    }
    return result;
  } else if constexpr (
      (priority1 <= detail::addition_priority_v<mat_structure::diagonal>)&&(
          priority2 == detail::addition_priority_v<mat_structure::identity> ||
          priority2 == detail::addition_priority_v<mat_structure::scalar>)) {
    // dense + identity/scalar
    result_type result(M1);
    for (int i = 0; i < result.get_row_count(); ++i) {
      result(i, i) += M2(0, 0);
    }
    return result;
  } else if constexpr (
      (priority1 == detail::addition_priority_v<mat_structure::diagonal>)&&(
          priority2 == detail::addition_priority_v<mat_structure::diagonal>)) {
    // diag + diag
    result_type result(M1);
    result += result_type(M2);
    return result;
  } else if constexpr (
      (priority1 == detail::addition_priority_v<mat_structure::diagonal>)&&(
          priority2 < detail::addition_priority_v<mat_structure::diagonal>)) {
    // diag + dense
    result_type result(M2);
    for (int i = 0; i < result.get_row_count(); ++i) {
      result(i, i) += M1(i, i);
    }
    return result;
  } else if constexpr (
      (priority1 < detail::addition_priority_v<mat_structure::diagonal>)&&(
          priority2 == detail::addition_priority_v<mat_structure::diagonal>)) {
    // dense + diag
    result_type result(M1);
    for (int i = 0; i < result.get_row_count(); ++i) {
      result(i, i) += M2(i, i);
    }
    return result;
  } else {
    // general case
    result_type result(M1);
    for (int j = 0; j < M1.get_col_count(); ++j) {
      for (int i = 0; i < M1.get_row_count(); ++i) {
        result(i, j) += M2(i, j);
      }
    }
    return result;
  }
}

/**
 * General substraction operator for any type of matrices. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M1 first matrix (first operand).
 * \param M2 second matrix (second operand).
 * \return General column-major matrix.
 * \throw std::range_error if the two matrix dimensions do not fit together.
 * \test PASSED
 */
template <typename Matrix1, typename Matrix2>
std::enable_if_t<is_readable_matrix_v<Matrix1> && is_readable_matrix_v<Matrix2>,
                 mat_addition_result_t<Matrix1, Matrix2>>
operator-(const Matrix1& M1, const Matrix2& M2) {
  if ((M1.get_row_count() != M2.get_row_count()) ||
      (M1.get_col_count() != M2.get_col_count())) {
    throw std::range_error("Matrix dimension mismatch.");
  }
  using result_type = mat_addition_result_t<Matrix1, Matrix2>;

  constexpr auto priority1 = mat_addition_priority_v<Matrix1>;
  constexpr auto priority2 = mat_addition_priority_v<Matrix2>;
  if constexpr (priority1 == detail::addition_priority_v<mat_structure::nil>) {
    // nil - any
    return result_type(-M2);
  } else if constexpr (priority2 ==
                       detail::addition_priority_v<mat_structure::nil>) {
    // any - nil
    return result_type(M1);
  } else if constexpr (
      (priority1 == detail::addition_priority_v<mat_structure::identity> ||
       priority1 ==
           detail::addition_priority_v<
               mat_structure::scalar>)&&(priority2 ==
                                             detail::addition_priority_v<
                                                 mat_structure::identity> ||
                                         priority2 ==
                                             detail::addition_priority_v<
                                                 mat_structure::scalar>)) {
    // identity/scalar - identity/scalar
    result_type result(M1.get_row_count());
    if (M1.get_row_count() > 0) {
      result(0, 0) = M1(0, 0) - M2(0, 0);
    }
    return result;
  } else if constexpr (
      (priority1 == detail::addition_priority_v<mat_structure::identity> ||
       priority1 ==
           detail::addition_priority_v<
               mat_structure::scalar>)&&(priority2 <=
                                         detail::addition_priority_v<
                                             mat_structure::diagonal>)) {
    // identity/scalar - dense
    result_type result(-M2);
    for (int i = 0; i < result.get_row_count(); ++i) {
      result(i, i) += M1(0, 0);
    }
    return result;
  } else if constexpr (
      (priority1 <= detail::addition_priority_v<mat_structure::diagonal>)&&(
          priority2 == detail::addition_priority_v<mat_structure::identity> ||
          priority2 == detail::addition_priority_v<mat_structure::scalar>)) {
    // dense - identity/scalar
    result_type result(M1);
    for (int i = 0; i < result.get_row_count(); ++i) {
      result(i, i) -= M2(0, 0);
    }
    return result;
  } else if constexpr (
      (priority1 == detail::addition_priority_v<mat_structure::diagonal>)&&(
          priority2 == detail::addition_priority_v<mat_structure::diagonal>)) {
    // diag - diag
    result_type result(M1);
    result -= result_type(M2);
    return result;
  } else if constexpr (
      (priority1 == detail::addition_priority_v<mat_structure::diagonal>)&&(
          priority2 < detail::addition_priority_v<mat_structure::diagonal>)) {
    // diag - dense
    result_type result(-M2);
    for (int i = 0; i < result.get_row_count(); ++i) {
      result(i, i) += M1(i, i);
    }
    return result;
  } else if constexpr (
      (priority1 < detail::addition_priority_v<mat_structure::diagonal>)&&(
          priority2 == detail::addition_priority_v<mat_structure::diagonal>)) {
    // dense - diag
    result_type result(M1);
    for (int i = 0; i < result.get_row_count(); ++i) {
      result(i, i) -= M2(i, i);
    }
    return result;
  } else {
    // general case
    result_type result(M1);
    for (int j = 0; j < M1.get_col_count(); ++j) {
      for (int i = 0; i < M1.get_row_count(); ++i) {
        result(i, j) -= M2(i, j);
      }
    }
    return result;
  }
}

/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

/**
 * Equality Comparison operator for general matrices, component-wise.
 * \test PASSED
 */
template <typename Matrix1, typename Matrix2>
std::enable_if_t<is_readable_matrix_v<Matrix1> && is_readable_matrix_v<Matrix2>,
                 bool>
operator==(const Matrix1& M1, const Matrix2& M2) {
  if ((M1.get_row_count() != M2.get_row_count()) ||
      (M1.get_col_count() != M2.get_col_count())) {
    return false;
  }
  for (int j = 0; j < M1.get_col_count(); ++j) {
    for (int i = 0; i < M1.get_row_count(); ++i) {
      if (M1(i, j) != M2(i, j)) {
        return false;
      }
    }
  }
  return true;
}

/**
 * Inequality Comparison operator for general matrices, component-wise.
 * \test PASSED
 */
template <typename Matrix1, typename Matrix2>
std::enable_if_t<is_readable_matrix_v<Matrix1> && is_readable_matrix_v<Matrix2>,
                 bool>
operator!=(const Matrix1& M1, const Matrix2& M2) {
  return !(M1 == M2);
}

}  // namespace ReaK

#endif

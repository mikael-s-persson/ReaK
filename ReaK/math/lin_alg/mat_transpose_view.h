/**
 * \file mat_transpose_view.h
 *
 * This library provides class templates to create transposed matrix views.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date December 2011
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

#ifndef REAK_MATH_LIN_ALG_MAT_TRANSPOSE_VIEW_H_
#define REAK_MATH_LIN_ALG_MAT_TRANSPOSE_VIEW_H_

#include "ReaK/math/lin_alg/mat_alg_general.h"

#include <type_traits>

namespace ReaK {

template <ReadableMatrix Matrix>
class mat_const_transpose_view {
 public:
  using self = mat_const_transpose_view<Matrix>;

  using value_type = mat_value_type_t<Matrix>;

  using reference = typename mat_traits<Matrix>::reference;
  using const_reference = typename mat_traits<Matrix>::const_reference;
  using pointer = typename mat_traits<Matrix>::pointer;
  using const_pointer = typename mat_traits<Matrix>::const_pointer;

  using col_iterator = typename mat_traits<Matrix>::col_iterator;
  using const_col_iterator = typename mat_traits<Matrix>::const_col_iterator;
  using row_iterator = typename mat_traits<Matrix>::row_iterator;
  using const_row_iterator = typename mat_traits<Matrix>::const_row_iterator;

  using size_type = mat_size_type_t<Matrix>;
  using difference_type = typename mat_traits<Matrix>::difference_type;

  static constexpr unsigned int static_row_count =
      mat_traits<Matrix>::static_col_count;
  static constexpr unsigned int static_col_count =
      mat_traits<Matrix>::static_row_count;
  static constexpr mat_alignment::tag alignment = mat_traits<Matrix>::alignment;
  static constexpr mat_structure::tag structure = mat_traits<Matrix>::structure;

 private:
  const Matrix* m;

  self& operator=(const self&);
  explicit mat_const_transpose_view(Matrix&&);

 public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_const_transpose_view(const Matrix& aM) : m(&aM) {}

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()(int i, int j) const { return (*m)(j, i); }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  int get_row_count() const noexcept { return m->get_col_count(); }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  int get_col_count() const noexcept { return m->get_row_count(); }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept {
    return {m->get_col_count(), m->get_row_count()};
  }

  /**
   * General negation operator for any type of matrices. This is a default operator
   * that will be called if no better special-purpose overload exists.
   * \return General column-major matrix.
   * \test PASSED
   */
  mat<value_type, mat_structure::rectangular> operator-() const {
    mat<value_type, mat_structure::rectangular> result(*this);
    for (int j = 0; j < result.get_col_count(); ++j) {
      for (int i = 0; i < result.get_row_count(); ++i) {
        result(i, j) = -result(i, j);
      }
    }
    return result;
  }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend mat<value_type, mat_structure::rectangular> transpose(const self& M) {
    return mat<value_type, mat_structure::rectangular>(*(M.m));
  }
};

template <typename Matrix>
static constexpr bool is_square_matrix_v<mat_const_transpose_view<Matrix>> = is_square_matrix_v<Matrix>;

template <typename Matrix>
static constexpr bool is_symmetric_matrix_v<mat_const_transpose_view<Matrix>> = is_symmetric_matrix_v<Matrix>;

template <typename Matrix>
static constexpr bool is_diagonal_matrix_v<mat_const_transpose_view<Matrix>> = is_diagonal_matrix_v<Matrix>;

template <WritableMatrix Matrix>
class mat_transpose_view {
 public:
  using self = mat_transpose_view<Matrix>;

  using value_type = mat_value_type_t<Matrix>;

  using reference = typename mat_traits<Matrix>::reference;
  using const_reference = typename mat_traits<Matrix>::const_reference;
  using pointer = typename mat_traits<Matrix>::pointer;
  using const_pointer = typename mat_traits<Matrix>::const_pointer;

  using col_iterator = typename mat_traits<Matrix>::col_iterator;
  using const_col_iterator = typename mat_traits<Matrix>::const_col_iterator;
  using row_iterator = typename mat_traits<Matrix>::row_iterator;
  using const_row_iterator = typename mat_traits<Matrix>::const_row_iterator;

  using size_type = mat_size_type_t<Matrix>;
  using difference_type = typename mat_traits<Matrix>::difference_type;

  static constexpr unsigned int static_row_count =
      mat_traits<Matrix>::static_col_count;
  static constexpr unsigned int static_col_count =
      mat_traits<Matrix>::static_row_count;
  static constexpr mat_alignment::tag alignment = mat_traits<Matrix>::alignment;
  static constexpr mat_structure::tag structure = mat_traits<Matrix>::structure;

 private:
  Matrix* m;

 public:
  /**
   * Constructs a tranposed view of a matrix.
   */
  explicit mat_transpose_view(Matrix& aM) : m(&aM) {}

  /**
   * Standard assignment operator.
   */
  template <ReadableMatrix Matrix2>
  self& operator=(const Matrix2& rhs) {
    *m = mat_const_transpose_view<Matrix2>(rhs);
    return *this;
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-write access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  reference operator()(int i, int j) { return (*m)(j, i); }
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()(int i, int j) const { return (*m)(j, i); }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * TEST PASSED
   */
  int get_row_count() const noexcept { return m->get_col_count(); }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * TEST PASSED
   */
  int get_col_count() const noexcept { return m->get_row_count(); }

  /**
   * Sets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * TEST PASSED
   */
  void set_row_count(int aRowCount) { m->set_col_count(aRowCount); }

  /**
   * Sets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * TEST PASSED
   */
  void get_col_count(int aColCount) { m->set_row_count(aColCount); }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept {
    return {m->get_col_count(), m->get_row_count()};
  }

  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator+=(const Matrix2& M) {
    (*m) += mat_const_transpose_view<Matrix2>(M);
    return *this;
  }

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator-=(const Matrix2& M) {
    (*m) -= mat_const_transpose_view<Matrix2>(M);
    return *this;
  }

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=(const value_type& S) {
    (*m) *= S;
    return *this;
  }

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator*=(const Matrix2& M) {
    if constexpr (!ReadableMatrix<Matrix2>) {
      return *this *= value_type(M);
    } else {
      *this = *this * M;
      return *this;
    }
  }

  /** WORKS FOR ALL
   * General negation operator for any type of matrices. This is a default operator
   * that will be called if no better special-purpose overload exists.
   * \return General column-major matrix.
   * \test PASSED
   */
  mat<value_type, mat_structure::rectangular> operator-() const {
    mat<value_type, mat_structure::rectangular> result(*this);
    for (int j = 0; j < result.get_col_count(); ++j) {
      for (int i = 0; i < result.get_row_count(); ++i) {
        result(i, j) = -result(i, j);
      }
    }
    return result;
  }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend mat<value_type, mat_structure::rectangular> transpose(const self& M) {
    return mat<value_type, mat_structure::rectangular>(*(M.m));
  }
};

template <typename Matrix>
static constexpr bool is_fully_writable_matrix_v<mat_transpose_view<Matrix>> = is_fully_writable_matrix_v<Matrix>;

template <typename Matrix>
static constexpr bool is_square_matrix_v<mat_transpose_view<Matrix>> = is_square_matrix_v<Matrix>;

template <typename Matrix>
static constexpr bool is_symmetric_matrix_v<mat_transpose_view<Matrix>> = is_symmetric_matrix_v<Matrix>;

template <typename Matrix>
static constexpr bool is_diagonal_matrix_v<mat_transpose_view<Matrix>> = is_diagonal_matrix_v<Matrix>;

template <WritableMatrix Matrix>
auto transpose_view(Matrix& M) {
  return mat_transpose_view<std::decay_t<Matrix>>(M);
}

template <ReadableMatrix Matrix>
auto transpose_view(const Matrix& M) {
  return mat_const_transpose_view<std::decay_t<Matrix>>(M);
}

template <ReadableMatrix Matrix>
auto transpose_view(Matrix&& M) {
  return mat_const_transpose_view<std::decay_t<Matrix>>(M);
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_TRANSPOSE_VIEW_H_

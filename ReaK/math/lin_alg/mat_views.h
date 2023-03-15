/**
 * \file mat_views.h
 *
 * This library provides a number of class templates to create matrix views. A matrix
 * view simply means that a sub-block of a matrix is used as if it was a matrix in its
 * own right. This can be very useful to set sub-blocks to other values or to use a
 * sub-block in a matrix expression (e.g. applying the operation on the entire matrix
 * is not practical or efficient).
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

#ifndef REAK_MATH_LIN_ALG_MAT_VIEWS_H_
#define REAK_MATH_LIN_ALG_MAT_VIEWS_H_

#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/mat_traits.h"
#include "ReaK/math/lin_alg/vect_concepts.h"
#include "ReaK/math/lin_alg/vect_views.h"

#include <type_traits>

namespace ReaK {

/**
 * This class template
 *
 *
 */
template <typename Matrix>
class mat_copy_sub_block {
 public:
  using self = mat_copy_sub_block<Matrix>;

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

  static constexpr unsigned int static_row_count = 0;
  static constexpr unsigned int static_col_count = 0;
  static constexpr mat_alignment::tag alignment = mat_traits<Matrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::rectangular;

 private:
  Matrix m;
  int rowOffset;
  int colOffset;
  int rowCount;
  int colCount;

 public:
  /**
   * Default constructor.
   */
  mat_copy_sub_block()
      : m(), rowOffset(0), colOffset(0), rowCount(0), colCount(0) {}

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_copy_sub_block(const Matrix& aM)
      : m(aM),
        rowOffset(0),
        colOffset(0),
        rowCount(aM.get_row_count()),
        colCount(aM.get_col_count()) {}

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aRowCount The number of rows for the sub-block.
   * \param aColCount The number of columns for the sub-block.
   * \param aRowOffset The row-offset from the start of the matrix.
   * \param aColOffset The column-offset from the start of the matrix.
   */
  mat_copy_sub_block(const Matrix& aM, int aRowCount, int aColCount,
                     int aRowOffset = 0, int aColOffset = 0)
      : m(aM),
        rowOffset(aRowOffset),
        colOffset(aColOffset),
        rowCount(aRowCount),
        colCount(aColCount) {}

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_copy_sub_block(Matrix&& aM)
      : m(std::move(aM)), rowOffset(0), colOffset(0), rowCount(0), colCount(0) {
    rowCount = m.get_row_count();
    colCount = m.get_col_count();
  }

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aRowCount The number of rows for the sub-block.
   * \param aColCount The number of columns for the sub-block.
   * \param aRowOffset The row-offset from the start of the matrix.
   * \param aColOffset The column-offset from the start of the matrix.
   */
  mat_copy_sub_block(Matrix&& aM, int aRowCount, int aColCount,
                     int aRowOffset = 0, int aColOffset = 0)
      : m(std::move(aM)),
        rowOffset(aRowOffset),
        colOffset(aColOffset),
        rowCount(aRowCount),
        colCount(aColCount) {}

  /**
   * Standard copy-constructor.
   */
  mat_copy_sub_block(const self& aObj) = default;

  /**
   * Standard move-constructor.
   */
  mat_copy_sub_block(self&& aObj) noexcept = default;

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.rowOffset, rhs.rowOffset);
    swap(lhs.colOffset, rhs.colOffset);
    swap(lhs.rowCount, rhs.rowCount);
    swap(lhs.colCount, rhs.colCount);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(self&& rhs) noexcept = default;
  self& operator=(const self& rhs) = default;

  /**
   * Standard assignment operator.
   */
  template <typename Matrix2>
  self& operator=(const Matrix2& rhs) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((rhs.get_row_count() != rowCount) ||
        (rhs.get_col_count() != colCount)) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    for (int j = 0; j < colCount; ++j) {
      for (int i = 0; i < rowCount; ++i) {
        m(rowOffset + i, colOffset + j) = rhs(i, j);
      }
    }
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
  reference operator()(int i, int j) { return m(rowOffset + i, colOffset + j); }
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()(int i, int j) const {
    return m(rowOffset + i, colOffset + j);
  }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  int get_row_count() const noexcept { return rowCount; }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  int get_col_count() const noexcept { return colCount; }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept { return {rowCount, colCount}; }

  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator+=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != colCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int j = 0; j < colCount; ++j) {
      for (int i = 0; i < rowCount; ++i) {
        m(rowOffset + i, colOffset + j) += M(i, j);
      }
    }
    return *this;
  }

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator-=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != colCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int j = 0; j < colCount; ++j) {
      for (int i = 0; i < rowCount; ++i) {
        m(rowOffset + i, colOffset + j) -= M(i, j);
      }
    }
    return *this;
  }

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=(const value_type& S) {
    for (int j = 0; j < colCount; ++j) {
      for (int i = 0; i < rowCount; ++i) {
        m(rowOffset + i, colOffset + j) *= S;
      }
    }
    return *this;
  }

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator*=(const Matrix2& M) {
    if constexpr (!is_readable_matrix_v<Matrix2>) {
      return *this *= value_type(M);
    } else {
      if ((M.get_col_count() != colCount) || (M.get_row_count() != colCount)) {
        throw std::range_error("Matrix Dimension Mismatch.");
      }
      *this = *this * M;
      return *this;
    }
  }
};

template <typename Matrix>
struct is_readable_matrix<mat_copy_sub_block<Matrix>> {
  static constexpr bool value = is_readable_matrix_v<Matrix>;
  using type = is_readable_matrix<Matrix>;
};

template <typename Matrix>
struct is_writable_matrix<mat_copy_sub_block<Matrix>> {
  static constexpr bool value = is_fully_writable_matrix_v<Matrix>;
  using type = is_fully_writable_matrix<Matrix>;
};

template <typename Matrix>
struct is_fully_writable_matrix<mat_copy_sub_block<Matrix>> {
  static constexpr bool value = is_fully_writable_matrix_v<Matrix>;
  using type = is_fully_writable_matrix<Matrix>;
};

template <typename Matrix>
struct is_row_resizable_matrix<mat_copy_sub_block<Matrix>> {
  static constexpr bool value = false;
  using type = is_row_resizable_matrix<mat_copy_sub_block<Matrix>>;
};

template <typename Matrix>
struct is_col_resizable_matrix<mat_copy_sub_block<Matrix>> {
  static constexpr bool value = false;
  using type = is_col_resizable_matrix<mat_copy_sub_block<Matrix>>;
};

template <typename Matrix>
class mat_sub_block {
 public:
  using self = mat_sub_block<Matrix>;

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

  static constexpr unsigned int static_row_count = 0;
  static constexpr unsigned int static_col_count = 0;
  static constexpr mat_alignment::tag alignment = mat_traits<Matrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::rectangular;

 private:
  Matrix* m;
  int rowOffset;
  int colOffset;
  int rowCount;
  int colCount;

 public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_sub_block(Matrix& aM)
      : m(&aM),
        rowOffset(0),
        colOffset(0),
        rowCount(aM.get_row_count()),
        colCount(aM.get_col_count()) {}

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aRowCount The number of rows for the sub-block.
   * \param aColCount The number of columns for the sub-block.
   * \param aRowOffset The row-offset from the start of the matrix.
   * \param aColOffset The column-offset from the start of the matrix.
   */
  mat_sub_block(Matrix& aM, int aRowCount, int aColCount, int aRowOffset = 0,
                int aColOffset = 0)
      : m(&aM),
        rowOffset(aRowOffset),
        colOffset(aColOffset),
        rowCount(aRowCount),
        colCount(aColCount) {}

  /**
   * Standard copy-constructor.
   */
  mat_sub_block(const self& aObj) = default;

  /**
   * Standard move-constructor.
   */
  mat_sub_block(self&& aObj) noexcept = default;

  /**
   * Standard swap function (shallow).
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.rowOffset, rhs.rowOffset);
    swap(lhs.colOffset, rhs.colOffset);
    swap(lhs.rowCount, rhs.rowCount);
    swap(lhs.colCount, rhs.colCount);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(const self& rhs) {
    if ((rhs.get_row_count() != rowCount) ||
        (rhs.get_col_count() != colCount)) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    for (int j = 0; j < colCount; ++j) {
      for (int i = 0; i < rowCount; ++i) {
        (*m)(rowOffset + i, colOffset + j) = rhs(i, j);
      }
    }
    return *this;
  }

  /**
   * Standard assignment operator.
   */
  template <typename Matrix2>
  self& operator=(const Matrix2& rhs) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((rhs.get_row_count() != rowCount) ||
        (rhs.get_col_count() != colCount)) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    for (int j = 0; j < colCount; ++j) {
      for (int i = 0; i < rowCount; ++i) {
        (*m)(rowOffset + i, colOffset + j) = rhs(i, j);
      }
    }
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
  reference operator()(int i, int j) {
    return (*m)(rowOffset + i, colOffset + j);
  }
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()(int i, int j) const {
    return (*m)(rowOffset + i, colOffset + j);
  }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  int get_row_count() const noexcept { return rowCount; }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  int get_col_count() const noexcept { return colCount; }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept { return {rowCount, colCount}; }

  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator+=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != colCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int j = 0; j < colCount; ++j) {
      for (int i = 0; i < rowCount; ++i) {
        (*m)(rowOffset + i, colOffset + j) += M(i, j);
      }
    }
    return *this;
  }

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator-=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != colCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int j = 0; j < colCount; ++j) {
      for (int i = 0; i < rowCount; ++i) {
        (*m)(rowOffset + i, colOffset + j) -= M(i, j);
      }
    }
    return *this;
  }

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=(const value_type& S) {
    for (int j = 0; j < colCount; ++j) {
      for (int i = 0; i < rowCount; ++i) {
        (*m)(rowOffset + i, colOffset + j) *= S;
      }
    }
    return *this;
  }

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator*=(const Matrix2& M) {
    if constexpr (!is_readable_matrix_v<Matrix2>) {
      return *this *= value_type(M);
    } else {
      if ((M.get_col_count() != colCount) || (M.get_row_count() != colCount)) {
        throw std::range_error("Matrix Dimension Mismatch.");
      }
      *this = *this * M;
      return *this;
    }
  }
};

template <typename Matrix>
struct is_readable_matrix<mat_sub_block<Matrix>> {
  static constexpr bool value = is_readable_matrix_v<Matrix>;
  using type = is_readable_matrix<Matrix>;
};

template <typename Matrix>
struct is_writable_matrix<mat_sub_block<Matrix>> {
  static constexpr bool value = is_fully_writable_matrix_v<Matrix>;
  using type = is_fully_writable_matrix<Matrix>;
};

template <typename Matrix>
struct is_fully_writable_matrix<mat_sub_block<Matrix>> {
  static constexpr bool value = is_fully_writable_matrix_v<Matrix>;
  using type = is_fully_writable_matrix<Matrix>;
};

template <typename Matrix>
struct is_row_resizable_matrix<mat_sub_block<Matrix>> {
  static constexpr bool value = false;
  using type = is_row_resizable_matrix<mat_sub_block<Matrix>>;
};

template <typename Matrix>
struct is_col_resizable_matrix<mat_sub_block<Matrix>> {
  static constexpr bool value = false;
  using type = is_col_resizable_matrix<mat_sub_block<Matrix>>;
};

template <typename Matrix>
class mat_const_sub_block {
 public:
  using self = mat_const_sub_block<Matrix>;

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

  static constexpr unsigned int static_row_count = 0;
  static constexpr unsigned int static_col_count = 0;
  static constexpr mat_alignment::tag alignment = mat_traits<Matrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::rectangular;

 private:
  const Matrix* m;
  int rowOffset;
  int colOffset;
  int rowCount;
  int colCount;

 public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_const_sub_block(const Matrix& aM)
      : m(&aM),
        rowOffset(0),
        colOffset(0),
        rowCount(aM.get_row_count()),
        colCount(aM.get_col_count()) {}

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aRowCount The number of rows for the sub-block.
   * \param aColCount The number of columns for the sub-block.
   * \param aRowOffset The row-offset from the start of the matrix.
   * \param aColOffset The column-offset from the start of the matrix.
   */
  mat_const_sub_block(const Matrix& aM, int aRowCount, int aColCount,
                      int aRowOffset = 0, int aColOffset = 0)
      : m(&aM),
        rowOffset(aRowOffset),
        colOffset(aColOffset),
        rowCount(aRowCount),
        colCount(aColCount) {}

  /**
   * Standard copy-constructor.
   */
  mat_const_sub_block(const self& aObj) = default;

  /**
   * Standard move-constructor.
   */
  mat_const_sub_block(self&& aObj) noexcept = default;

  self& operator=(const self&) = delete;
  explicit mat_const_sub_block(Matrix&&) = delete;
  mat_const_sub_block(Matrix&&, int, int, int aRowOffset = 0,
                      int aColOffset = 0) = delete;

  /**
   * Standard swap function (shallow).
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.rowOffset, rhs.rowOffset);
    swap(lhs.colOffset, rhs.colOffset);
    swap(lhs.rowCount, rhs.rowCount);
    swap(lhs.colCount, rhs.colCount);
  }

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
  value_type operator()(int i, int j) const {
    return (*m)(rowOffset + i, colOffset + j);
  }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  int get_row_count() const noexcept { return rowCount; }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  int get_col_count() const noexcept { return colCount; }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept { return {rowCount, colCount}; }
};

template <typename Matrix>
struct is_readable_matrix<mat_const_sub_block<Matrix>> {
  static constexpr bool value = is_readable_matrix_v<Matrix>;
  using type = is_readable_matrix<Matrix>;
};

template <typename Matrix>
struct is_writable_matrix<mat_const_sub_block<Matrix>> {
  static constexpr bool value = false;
  using type = is_writable_matrix<mat_const_sub_block<Matrix>>;
};

template <typename Matrix>
struct is_fully_writable_matrix<mat_const_sub_block<Matrix>> {
  static constexpr bool value = false;
  using type = is_fully_writable_matrix<mat_const_sub_block<Matrix>>;
};

template <typename Matrix>
struct is_row_resizable_matrix<mat_const_sub_block<Matrix>> {
  static constexpr bool value = false;
  using type = is_row_resizable_matrix<mat_const_sub_block<Matrix>>;
};

template <typename Matrix>
struct is_col_resizable_matrix<mat_const_sub_block<Matrix>> {
  static constexpr bool value = false;
  using type = is_col_resizable_matrix<mat_const_sub_block<Matrix>>;
};

template <typename Matrix>
struct mat_copy_sub_block_factory {
  Matrix m;
  explicit mat_copy_sub_block_factory(Matrix&& aM) : m(std::move(aM)) {}
  mat_copy_sub_block<Matrix> operator()(const std::pair<int, int>& rows,
                                        const std::pair<int, int>& cols) {
    return mat_copy_sub_block<Matrix>(std::move(m), rows.second - rows.first,
                                      cols.second - cols.first, rows.first,
                                      cols.first);
  }
};

template <typename Matrix>
struct mat_sub_block_factory {
  Matrix& m;
  explicit mat_sub_block_factory(Matrix& aM) : m(aM) {}
  mat_sub_block<Matrix> operator()(const std::pair<int, int>& rows,
                                   const std::pair<int, int>& cols) {
    return mat_sub_block<Matrix>(m, rows.second - rows.first,
                                 cols.second - cols.first, rows.first,
                                 cols.first);
  }
};

template <typename Matrix>
struct mat_const_sub_block_factory {
  const Matrix& m;
  explicit mat_const_sub_block_factory(const Matrix& aM) : m(aM) {}
  mat_const_sub_block<Matrix> operator()(const std::pair<int, int>& rows,
                                         const std::pair<int, int>& cols) {
    return mat_const_sub_block<Matrix>(m, rows.second - rows.first,
                                       cols.second - cols.first, rows.first,
                                       cols.first);
  }
};

template <typename Matrix>
std::enable_if_t<is_readable_matrix_v<Matrix>, mat_sub_block_factory<Matrix>>
sub(Matrix& M) {
  return mat_sub_block_factory<Matrix>(M);
}

template <typename Matrix>
std::enable_if_t<is_readable_matrix_v<Matrix>,
                 mat_const_sub_block_factory<Matrix>>
sub(const Matrix& M) {
  return mat_const_sub_block_factory<Matrix>(M);
}

template <typename Matrix>
std::enable_if_t<is_readable_matrix_v<Matrix>,
                 mat_copy_sub_block_factory<Matrix>>
sub_copy(const Matrix& M) {
  return mat_copy_sub_block_factory<Matrix>(M);
}

template <typename Matrix>
std::enable_if_t<is_readable_matrix_v<Matrix>,
                 mat_copy_sub_block_factory<Matrix>>
sub(Matrix&& M) {
  return mat_copy_sub_block_factory<Matrix>(std::move(M));
}

template <typename Matrix,
          mat_structure::tag Structure = mat_traits<Matrix>::structure>
class mat_copy_sub_sym_block {
  static_assert(sizeof(Matrix) == 0, "invalid matrix type");
};

template <typename Matrix,
          mat_structure::tag Structure = mat_traits<Matrix>::structure>
class mat_sub_sym_block {
  static_assert(sizeof(Matrix) == 0, "invalid matrix type");
};

template <typename Matrix,
          mat_structure::tag Structure = mat_traits<Matrix>::structure>
class mat_const_sub_sym_block {
  static_assert(sizeof(Matrix) == 0, "invalid matrix type");
};

template <typename Matrix>
class mat_copy_sub_sym_block<Matrix, mat_structure::symmetric> {
 public:
  using self = mat_sub_sym_block<Matrix, mat_structure::symmetric>;

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

  static constexpr unsigned int static_row_count = 0;
  static constexpr unsigned int static_col_count = 0;
  static constexpr mat_alignment::tag alignment = mat_traits<Matrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::symmetric;

 private:
  Matrix m;
  int offset;
  int rowCount;

 public:
  /**
   * Default constructor.
   */
  mat_copy_sub_sym_block() : m(), offset(0), rowCount(0) {}

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_copy_sub_sym_block(const Matrix& aM)
      : m(aM), offset(0), rowCount(aM.get_row_count()) {}

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_copy_sub_sym_block(const Matrix& aM, int aSize, int aOffset = 0)
      : m(aM), offset(aOffset), rowCount(aSize) {}

  /**
   * Standard copy-constructor.
   */
  explicit mat_copy_sub_sym_block(const self& aObj) = default;

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_copy_sub_sym_block(Matrix&& aM)
      : m(std::move(aM)), offset(0), rowCount(0) {
    rowCount = m.get_row_count();
  }

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_copy_sub_sym_block(Matrix&& aM, int aSize, int aOffset = 0)
      : m(std::move(aM)), offset(aOffset), rowCount(aSize) {}

  /**
   * Standard move-constructor.
   */
  explicit mat_copy_sub_sym_block(self&& aObj) noexcept = default;

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.offset, rhs.offset);
    swap(lhs.rowCount, rhs.rowCount);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(self&& rhs) noexcept = default;
  self& operator=(const self& rhs) = default;

  /**
   * Standard assignment operator.
   */
  template <typename Matrix2>
  self& operator=(const Matrix2& rhs) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((rhs.get_row_count() != rowCount) ||
        (rhs.get_col_count() != rowCount)) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    for (int j = 0; j < rowCount; ++j) {
      for (int i = j; i < rowCount; ++i) {
        m(offset + i, offset + j) = (rhs(i, j) + rhs(j, i)) * value_type(0.5);
      }
    }
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
  reference operator()(int i, int j) { return m(offset + i, offset + j); }
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()(int i, int j) const {
    return m(offset + i, offset + j);
  }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  int get_row_count() const noexcept { return rowCount; }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  int get_col_count() const noexcept { return rowCount; }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept { return {rowCount, rowCount}; }

  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator+=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int j = 0; j < rowCount; ++j) {
      for (int i = j; i < rowCount; ++i) {
        m(offset + i, offset + j) += (M(i, j) + M(j, i)) * value_type(0.5);
      }
    }
    return *this;
  }

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator-=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int j = 0; j < rowCount; ++j) {
      for (int i = j; i < rowCount; ++i) {
        m(offset + i, offset + j) -= (M(i, j) + M(j, i)) * value_type(0.5);
      }
    }
    return *this;
  }

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=(const value_type& S) {
    for (int j = 0; j < rowCount; ++j) {
      for (int i = j; i < rowCount; ++i) {
        m(offset + i, offset + j) *= S;
      }
    }
    return *this;
  }

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator*=(const Matrix2& M) {
    if constexpr (!is_readable_matrix_v<Matrix2>) {
      return *this *= value_type(M);
    } else {
      if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
        throw std::range_error("Matrix Dimension Mismatch.");
      }
      *this = *this * M;
      return *this;
    }
  }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose(const self& M) { return M; }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose_move(const self& M) { return M; }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend self&& transpose(self&& M) { return std::move(M); }

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace(const self& M) {
    value_type result(0.0);
    for (int i = 0; i < M.rowCount; ++i) {
      result += M(i, i);
    }
    return result;
  }
};

template <typename Matrix>
class mat_sub_sym_block<Matrix, mat_structure::symmetric> {
 public:
  using self = mat_sub_sym_block<Matrix, mat_structure::symmetric>;

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

  static constexpr unsigned int static_row_count = 0;
  static constexpr unsigned int static_col_count = 0;
  static constexpr mat_alignment::tag alignment = mat_traits<Matrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::symmetric;

 private:
  Matrix* m;
  int offset;
  int rowCount;

 public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_sub_sym_block(Matrix& aM)
      : m(&aM), offset(0), rowCount(aM.get_row_count()) {}

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_sub_sym_block(Matrix& aM, int aSize, int aOffset = 0)
      : m(&aM), offset(aOffset), rowCount(aSize) {}
  /**
   * Standard copy-constructor.
   */
  mat_sub_sym_block(const self& aObj) = default;

  /**
   * Standard move-constructor.
   */
  mat_sub_sym_block(self&& aObj) noexcept = default;

  /**
   * Standard swap function (shallow).
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.offset, rhs.offset);
    swap(lhs.rowCount, rhs.rowCount);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(const self& rhs) {
    if ((rhs.get_row_count() != rowCount) ||
        (rhs.get_col_count() != rowCount)) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    for (int j = 0; j < rowCount; ++j) {
      for (int i = j; i < rowCount; ++i) {
        (*m)(offset + i, offset + j) =
            (rhs(i, j) + rhs(j, i)) * value_type(0.5);
      }
    }
    return *this;
  }

  /**
   * Standard assignment operator.
   */
  template <typename Matrix2>
  self& operator=(const Matrix2& rhs) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((rhs.get_row_count() != rowCount) ||
        (rhs.get_col_count() != rowCount)) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    for (int j = 0; j < rowCount; ++j) {
      for (int i = j; i < rowCount; ++i) {
        (*m)(offset + i, offset + j) =
            (rhs(i, j) + rhs(j, i)) * value_type(0.5);
      }
    }
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
  reference operator()(int i, int j) { return (*m)(offset + i, offset + j); }
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()(int i, int j) const {
    return (*m)(offset + i, offset + j);
  }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  int get_row_count() const noexcept { return rowCount; }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  int get_col_count() const noexcept { return rowCount; }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept { return {rowCount, rowCount}; }

  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator+=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int j = 0; j < rowCount; ++j) {
      for (int i = j; i < rowCount; ++i) {
        (*m)(offset + i, offset + j) += (M(i, j) + M(j, i)) * value_type(0.5);
      }
    }
    return *this;
  }

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator-=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int j = 0; j < rowCount; ++j) {
      for (int i = j; i < rowCount; ++i) {
        (*m)(offset + i, offset + j) -= (M(i, j) + M(j, i)) * value_type(0.5);
      }
    }
    return *this;
  }

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=(const value_type& S) {
    for (int j = 0; j < rowCount; ++j) {
      for (int i = j; i < rowCount; ++i) {
        (*m)(offset + i, offset + j) *= S;
      }
    }
    return *this;
  }

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator*=(const Matrix2& M) {
    if constexpr (!is_readable_matrix_v<Matrix2>) {
      return *this *= value_type(M);
    } else {
      if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
        throw std::range_error("Matrix Dimension Mismatch.");
      }
      *this = *this * M;
      return *this;
    }
  }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose(const self& M) { return M; }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose_move(const self& M) { return M; }

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace(const self& M) {
    value_type result(0.0);
    for (int i = 0; i < M.rowCount; ++i) {
      result += M(i, i);
    }
    return result;
  }
};

template <typename Matrix>
class mat_const_sub_sym_block<Matrix, mat_structure::symmetric> {
 public:
  using self = mat_const_sub_sym_block<Matrix, mat_structure::symmetric>;

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

  static constexpr unsigned int static_row_count = 0;
  static constexpr unsigned int static_col_count = 0;
  static constexpr mat_alignment::tag alignment = mat_traits<Matrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::symmetric;

 private:
  const Matrix* m;
  int offset;
  int rowCount;

 public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_const_sub_sym_block(const Matrix& aM)
      : m(&aM), offset(0), rowCount(aM.get_row_count()) {}

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_const_sub_sym_block(const Matrix& aM, int aSize, int aOffset = 0)
      : m(&aM), offset(aOffset), rowCount(aSize) {}

  /**
   * Standard copy-constructor.
   */
  mat_const_sub_sym_block(const self& aObj) = default;

  /**
   * Standard move-constructor.
   */
  mat_const_sub_sym_block(self&& aObj) noexcept = default;

  self& operator=(const self&) = delete;
  self& operator=(self&&) = delete;

  explicit mat_const_sub_sym_block(Matrix&&) = delete;
  mat_const_sub_sym_block(Matrix&&, int, int aOffset = 0) = delete;

  /**
   * Standard swap function (shallow).
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.offset, rhs.offset);
    swap(lhs.rowCount, rhs.rowCount);
  }

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
  value_type operator()(int i, int j) const {
    return (*m)(offset + i, offset + j);
  }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  int get_row_count() const noexcept { return rowCount; }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  int get_col_count() const noexcept { return rowCount; }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept { return {rowCount, rowCount}; }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose(const self& M) { return M; }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose_move(const self& M) { return M; }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend self&& transpose(self&& M) { return std::move(M); }

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace(const self& M) {
    value_type result(0.0);
    for (int i = 0; i < M.rowCount; ++i) {
      result += M(i, i);
    }
    return result;
  }
};

template <typename Matrix>
class mat_copy_sub_sym_block<Matrix, mat_structure::skew_symmetric> {
 public:
  using self = mat_sub_sym_block<Matrix, mat_structure::skew_symmetric>;

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

  static constexpr unsigned int static_row_count = 0;
  static constexpr unsigned int static_col_count = 0;
  static constexpr mat_alignment::tag alignment = mat_traits<Matrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::skew_symmetric;

 private:
  Matrix m;
  int offset;
  int rowCount;

 public:
  /**
   * Default constructor.
   */
  mat_copy_sub_sym_block() : m(), offset(0), rowCount(0) {}

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_copy_sub_sym_block(const Matrix& aM)
      : m(aM), offset(0), rowCount(aM.get_row_count()) {}

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_copy_sub_sym_block(const Matrix& aM, int aSize, int aOffset = 0)
      : m(aM), offset(aOffset), rowCount(aSize) {}

  /**
   * Standard copy-constructor.
   */
  explicit mat_copy_sub_sym_block(const self& aObj) = default;

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_copy_sub_sym_block(Matrix&& aM)
      : m(std::move(aM)), offset(0), rowCount(0) {
    rowCount = m.get_row_count();
  }

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_copy_sub_sym_block(Matrix&& aM, int aSize, int aOffset = 0)
      : m(std::move(aM)), offset(aOffset), rowCount(aSize) {}

  /**
   * Standard move-constructor.
   */
  explicit mat_copy_sub_sym_block(self&& aObj) noexcept = default;

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.offset, rhs.offset);
    swap(lhs.rowCount, rhs.rowCount);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(self&& rhs) noexcept = default;
  self& operator=(const self& rhs) = default;

  /**
   * Standard assignment operator.
   */
  template <typename Matrix2>
  self& operator=(const Matrix2& rhs) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((rhs.get_row_count() != rowCount) ||
        (rhs.get_col_count() != rowCount)) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    for (int j = 0; j < rowCount; ++j) {
      for (int i = j; i < rowCount; ++i) {
        m(offset + i, offset + j) = (rhs(i, j) + rhs(j, i)) * value_type(0.5);
      }
    }
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
  reference operator()(int i, int j) { return m(offset + i, offset + j); }
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()(int i, int j) const {
    return m(offset + i, offset + j);
  }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  int get_row_count() const noexcept { return rowCount; }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  int get_col_count() const noexcept { return rowCount; }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept { return {rowCount, rowCount}; }

  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator+=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int j = 0; j < rowCount; ++j) {
      for (int i = j; i < rowCount; ++i) {
        m(offset + i, offset + j) += (M(i, j) + M(j, i)) * value_type(0.5);
      }
    }
    return *this;
  }

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator-=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int j = 0; j < rowCount; ++j) {
      for (int i = j; i < rowCount; ++i) {
        m(offset + i, offset + j) -= (M(i, j) + M(j, i)) * value_type(0.5);
      }
    }
    return *this;
  }

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=(const value_type& S) {
    for (int j = 0; j < rowCount; ++j) {
      for (int i = j; i < rowCount; ++i) {
        m(offset + i, offset + j) *= S;
      }
    }
    return *this;
  }

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator*=(const Matrix2& M) {
    if constexpr (!is_readable_matrix_v<Matrix2>) {
      return *this *= value_type(M);
    } else {
      if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
        throw std::range_error("Matrix Dimension Mismatch.");
      }
      *this = *this * M;
      return *this;
    }
  }

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace(const self& M) {
    value_type result(0.0);
    for (int i = 0; i < M.rowCount; ++i) {
      result += M(i, i);
    }
    return result;
  }
};

template <typename Matrix>
class mat_sub_sym_block<Matrix, mat_structure::skew_symmetric> {
 public:
  using self = mat_sub_sym_block<Matrix, mat_structure::skew_symmetric>;

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

  static constexpr unsigned int static_row_count = 0;
  static constexpr unsigned int static_col_count = 0;
  static constexpr mat_alignment::tag alignment = mat_traits<Matrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::skew_symmetric;

 private:
  Matrix* m;
  int offset;
  int rowCount;

 public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_sub_sym_block(Matrix& aM)
      : m(&aM), offset(0), rowCount(aM.get_row_count()) {}

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_sub_sym_block(Matrix& aM, int aSize, int aOffset = 0)
      : m(&aM), offset(aOffset), rowCount(aSize) {}

  /**
   * Standard copy-constructor.
   */
  mat_sub_sym_block(const self& aObj) = default;

  /**
   * Standard move-constructor.
   */
  mat_sub_sym_block(self&& aObj) noexcept = default;

  /**
   * Standard swap function (shallow).
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.offset, rhs.offset);
    swap(lhs.rowCount, rhs.rowCount);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(const self& rhs) {
    if ((rhs.get_row_count() != rowCount) ||
        (rhs.get_col_count() != rowCount)) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    for (int j = 1; j < rowCount; ++j) {
      for (int i = 0; i < j; ++i) {
        (*m)(offset + i, offset + j) =
            (rhs(i, j) - rhs(j, i)) * value_type(0.5);
      }
    }
    return *this;
  }

  /**
   * Standard assignment operator.
   */
  template <typename Matrix2>
  self& operator=(const Matrix2& rhs) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((rhs.get_row_count() != rowCount) ||
        (rhs.get_col_count() != rowCount)) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    for (int j = 1; j < rowCount; ++j) {
      for (int i = 0; i < j; ++i) {
        (*m)(offset + i, offset + j) =
            (rhs(i, j) - rhs(j, i)) * value_type(0.5);
      }
    }
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
  reference operator()(int i, int j) { return (*m)(offset + i, offset + j); }
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()(int i, int j) const {
    return (*m)(offset + i, offset + j);
  }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  int get_row_count() const noexcept { return rowCount; }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  int get_col_count() const noexcept { return rowCount; }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept { return {rowCount, rowCount}; }

  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator+=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int j = 0; j < rowCount; ++j) {
      for (int i = 0; i < j; ++i) {
        (*m)(offset + i, offset + j) += (M(i, j) - M(j, i)) * value_type(0.5);
      }
    }
    return *this;
  }

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator-=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int j = 0; j < rowCount; ++j) {
      for (int i = 0; i < j; ++i) {
        (*m)(offset + i, offset + j) -= (M(i, j) - M(j, i)) * value_type(0.5);
      }
    }
    return *this;
  }

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=(const value_type& S) {
    for (int j = 0; j < rowCount; ++j) {
      for (int i = 0; i < j; ++i) {
        (*m)(offset + i, offset + j) *= S;
      }
    }
    return *this;
  }

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator*=(const Matrix2& M) {
    if constexpr (!is_readable_matrix_v<Matrix2>) {
      return *this *= value_type(M);
    } else {
      if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
        throw std::range_error("Matrix Dimension Mismatch.");
      }
      *this = *this * M;
      return *this;
    }
  }

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace(const self& M) { return value_type(0.0); }
};

template <typename Matrix>
class mat_const_sub_sym_block<Matrix, mat_structure::skew_symmetric> {
 public:
  using self = mat_const_sub_sym_block<Matrix, mat_structure::skew_symmetric>;

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

  static constexpr unsigned int static_row_count = 0;
  static constexpr unsigned int static_col_count = 0;
  static constexpr mat_alignment::tag alignment = mat_traits<Matrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::skew_symmetric;

 private:
  const Matrix* m;
  int offset;
  int rowCount;

 public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_const_sub_sym_block(const Matrix& aM)
      : m(&aM), offset(0), rowCount(aM.get_row_count()) {}

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_const_sub_sym_block(const Matrix& aM, int aSize, int aOffset = 0)
      : m(&aM), offset(aOffset), rowCount(aSize) {}

  /**
   * Standard copy-constructor.
   */
  mat_const_sub_sym_block(const self& aObj) = default;

  /**
   * Standard move-constructor.
   */
  mat_const_sub_sym_block(self&& aObj) noexcept = default;

  self& operator=(const self&) = delete;

  explicit mat_const_sub_sym_block(Matrix&& aM) = delete;
  mat_const_sub_sym_block(Matrix&& aM, int aSize, int aOffset = 0) = delete;

  /**
   * Standard swap function (shallow).
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.offset, rhs.offset);
    swap(lhs.rowCount, rhs.rowCount);
  }

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
  value_type operator()(int i, int j) const {
    return (*m)(offset + i, offset + j);
  }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  int get_row_count() const noexcept { return rowCount; }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  int get_col_count() const noexcept { return rowCount; }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept { return {rowCount, rowCount}; }

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace(const self& M) { return value_type(0.0); }
};

template <typename Matrix>
class mat_copy_sub_sym_block<Matrix, mat_structure::diagonal> {
 public:
  using self = mat_copy_sub_sym_block<Matrix, mat_structure::diagonal>;

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

  static constexpr unsigned int static_row_count = 0;
  static constexpr unsigned int static_col_count = 0;
  static constexpr mat_alignment::tag alignment = mat_traits<Matrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::diagonal;

 private:
  Matrix m;
  int offset;
  int rowCount;

 public:
  /**
   * Default constructor.
   */
  mat_copy_sub_sym_block() : m(), offset(0), rowCount(0) {
    rowCount = m.get_row_count();
  }

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_copy_sub_sym_block(const Matrix& aM)
      : m(aM), offset(0), rowCount(aM.get_row_count()) {}

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_copy_sub_sym_block(const Matrix& aM, int aSize, int aOffset = 0)
      : m(aM), offset(aOffset), rowCount(aSize) {}

  /**
   * Standard copy-constructor.
   */
  mat_copy_sub_sym_block(const self& aObj) = default;

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_copy_sub_sym_block(Matrix&& aM)
      : m(std::move(aM)), offset(0), rowCount(0) {
    rowCount = m.get_row_count();
  }

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_copy_sub_sym_block(Matrix&& aM, int aSize, int aOffset = 0)
      : m(std::move(aM)), offset(aOffset), rowCount(aSize) {}

  /**
   * Standard move-constructor.
   */
  mat_copy_sub_sym_block(self&& aObj) noexcept = default;

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.offset, rhs.offset);
    swap(lhs.rowCount, rhs.rowCount);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(self&& rhs) noexcept = default;
  self& operator=(const self& rhs) = default;

  /**
   * Standard assignment operator.
   */
  template <typename Matrix2>
  self& operator=(const Matrix2& rhs) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((rhs.get_row_count() != rowCount) ||
        (rhs.get_col_count() != rowCount)) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    for (int i = 0; i < rowCount; ++i) {
      m(offset + i, offset + i) = rhs(i, i);
    }
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
  reference operator()(int i, int j) { return m(offset + i, offset + j); }
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()(int i, int j) const {
    return m(offset + i, offset + j);
  }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  int get_row_count() const noexcept { return rowCount; }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  int get_col_count() const noexcept { return rowCount; }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept { return {rowCount, rowCount}; }

  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator+=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int i = 0; i < rowCount; ++i) {
      m(offset + i, offset + i) += M(i, i);
    }
    return *this;
  }

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator-=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int i = 0; i < rowCount; ++i) {
      m(offset + i, offset + i) -= M(i, i);
    }
    return *this;
  }

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=(const value_type& S) {
    for (int i = 0; i < rowCount; ++i) {
      m(offset + i, offset + i) *= S;
    }
    return *this;
  }

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator*=(const Matrix2& M) {
    if constexpr (!is_readable_matrix_v<Matrix2>) {
      return *this *= value_type(M);
    } else {
      if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
        throw std::range_error("Matrix Dimension Mismatch.");
      }
      *this = *this * M;
      return *this;
    }
  }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose(const self& M) { return M; }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose_move(const self& M) { return M; }

  friend self&& transpose(self&& M) { return std::move(M); }

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace(const self& M) {
    value_type result(0.0);
    for (int i = 0; i < M.rowCount; ++i) {
      result += M(i, i);
    }
    return result;
  }
};

template <typename Matrix>
class mat_sub_sym_block<Matrix, mat_structure::diagonal> {
 public:
  using self = mat_sub_sym_block<Matrix, mat_structure::diagonal>;

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

  static constexpr unsigned int static_row_count = 0;
  static constexpr unsigned int static_col_count = 0;
  static constexpr mat_alignment::tag alignment = mat_traits<Matrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::diagonal;

 private:
  Matrix* m;
  int offset;
  int rowCount;

 public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_sub_sym_block(Matrix& aM)
      : m(&aM), offset(0), rowCount(aM.get_row_count()) {}

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_sub_sym_block(Matrix& aM, int aSize, int aOffset = 0)
      : m(&aM), offset(aOffset), rowCount(aSize) {}
  /**
   * Standard copy-constructor.
   */
  mat_sub_sym_block(const self& aObj) = default;

  /**
   * Standard move-constructor.
   */
  mat_sub_sym_block(self&& aObj) noexcept = default;

  /**
   * Standard swap function (shallow).
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.offset, rhs.offset);
    swap(lhs.rowCount, rhs.rowCount);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(const self& rhs) {
    if ((rhs.get_row_count() != rowCount) ||
        (rhs.get_col_count() != rowCount)) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    for (int i = 0; i < rowCount; ++i) {
      (*m)(offset + i, offset + i) = rhs(i, i);
    }
    return *this;
  }

  /**
   * Standard assignment operator.
   */
  template <typename Matrix2>
  self& operator=(const Matrix2& rhs) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((rhs.get_row_count() != rowCount) ||
        (rhs.get_col_count() != rowCount)) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    for (int i = 0; i < rowCount; ++i) {
      (*m)(offset + i, offset + i) = rhs(i, i);
    }
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
  reference operator()(int i, int j) { return (*m)(offset + i, offset + j); }
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()(int i, int j) const {
    return (*m)(offset + i, offset + j);
  }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  int get_row_count() const noexcept { return rowCount; }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  int get_col_count() const noexcept { return rowCount; }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept { return {rowCount, rowCount}; }

  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator+=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int i = 0; i < rowCount; ++i) {
      (*m)(offset + i, offset + i) += M(i, i);
    }
    return *this;
  }

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator-=(const Matrix2& M) {
    static_assert(is_readable_matrix_v<Matrix2>);
    if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int i = 0; i < rowCount; ++i) {
      (*m)(offset + i, offset + i) -= M(i, i);
    }
    return *this;
  }

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=(const value_type& S) {
    for (int i = 0; i < rowCount; ++i) {
      (*m)(offset + i, offset + i) *= S;
    }
    return *this;
  }

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template <typename Matrix2>
  self& operator*=(const Matrix2& M) {
    if constexpr (!is_readable_matrix_v<Matrix2>) {
      return *this *= value_type(M);
    } else {
      if ((M.get_col_count() != rowCount) || (M.get_row_count() != rowCount)) {
        throw std::range_error("Matrix Dimension Mismatch.");
      }
      *this = *this * M;
      return *this;
    }
  }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose(const self& M) { return M; }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose_move(const self& M) { return M; }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend self&& transpose(self&& M) { return std::move(M); }

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace(const self& M) {
    value_type result(0.0);
    for (int i = 0; i < M.rowCount; ++i) {
      result += M(i, i);
    }
    return result;
  }
};

template <typename Matrix>
class mat_const_sub_sym_block<Matrix, mat_structure::diagonal> {
 public:
  using self = mat_const_sub_sym_block<Matrix, mat_structure::diagonal>;

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

  static constexpr unsigned int static_row_count = 0;
  static constexpr unsigned int static_col_count = 0;
  static constexpr mat_alignment::tag alignment = mat_traits<Matrix>::alignment;
  static constexpr mat_structure::tag structure = mat_structure::diagonal;

 private:
  const Matrix* m;
  int offset;
  int rowCount;

 public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_const_sub_sym_block(const Matrix& aM)
      : m(&aM), offset(0), rowCount(aM.get_row_count()) {}

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_const_sub_sym_block(const Matrix& aM, int aSize, int aOffset = 0)
      : m(&aM), offset(aOffset), rowCount(aSize) {}

  /**
   * Standard copy-constructor.
   */
  mat_const_sub_sym_block(const self& aObj) = default;

  /**
   * Standard move-constructor.
   */
  mat_const_sub_sym_block(self&& aObj) noexcept = default;

  self& operator=(const self&) = delete;

  explicit mat_const_sub_sym_block(Matrix&&) = delete;
  mat_const_sub_sym_block(Matrix&&, int, int aOffset = 0) = delete;

  /**
   * Standard swap function (shallow).
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.offset, rhs.offset);
    swap(lhs.rowCount, rhs.rowCount);
  }

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
  value_type operator()(int i, int j) const {
    return (*m)(offset + i, offset + j);
  }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  int get_row_count() const noexcept { return rowCount; }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  int get_col_count() const noexcept { return rowCount; }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept { return {rowCount, rowCount}; }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose(const self& M) { return M; }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose_move(const self& M) { return M; }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend self&& transpose(self&& M) { return std::move(M); }

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace(const self& M) {
    value_type result(0.0);
    for (int i = 0; i < M.rowCount; ++i) {
      result += M(i, i);
    }
    return result;
  }
};

template <typename Matrix, mat_structure::tag Structure>
struct is_readable_matrix<mat_copy_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = is_readable_matrix_v<Matrix>;
  using type = is_readable_matrix<Matrix>;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_writable_matrix<mat_copy_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = is_writable_matrix_v<Matrix>;
  using type = is_writable_matrix<Matrix>;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_fully_writable_matrix<mat_copy_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = false;
  using type =
      is_fully_writable_matrix<mat_copy_sub_sym_block<Matrix, Structure>>;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_row_resizable_matrix<mat_copy_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = false;
  using type =
      is_row_resizable_matrix<mat_copy_sub_sym_block<Matrix, Structure>>;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_col_resizable_matrix<mat_copy_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = false;
  using type =
      is_col_resizable_matrix<mat_copy_sub_sym_block<Matrix, Structure>>;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_readable_matrix<mat_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = is_readable_matrix_v<Matrix>;
  using type = is_readable_matrix<Matrix>;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_writable_matrix<mat_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = is_writable_matrix_v<Matrix>;
  using type = is_writable_matrix<Matrix>;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_fully_writable_matrix<mat_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = false;
  using type = is_fully_writable_matrix<mat_sub_sym_block<Matrix, Structure>>;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_row_resizable_matrix<mat_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = false;
  using type = is_row_resizable_matrix<mat_sub_sym_block<Matrix, Structure>>;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_col_resizable_matrix<mat_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = false;
  using type = is_col_resizable_matrix<mat_sub_sym_block<Matrix, Structure>>;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_readable_matrix<mat_const_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = is_readable_matrix_v<Matrix>;
  using type = is_readable_matrix<Matrix>;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_writable_matrix<mat_const_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = false;
  using type = is_writable_matrix<mat_const_sub_sym_block<Matrix, Structure>>;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_fully_writable_matrix<mat_const_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = false;
  using type =
      is_fully_writable_matrix<mat_const_sub_sym_block<Matrix, Structure>>;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_row_resizable_matrix<mat_const_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = false;
  using type =
      is_row_resizable_matrix<mat_const_sub_sym_block<Matrix, Structure>>;
};

template <typename Matrix, mat_structure::tag Structure>
struct is_col_resizable_matrix<mat_const_sub_sym_block<Matrix, Structure>> {
  static constexpr bool value = false;
  using type =
      is_col_resizable_matrix<mat_const_sub_sym_block<Matrix, Structure>>;
};

template <typename Matrix>
struct mat_copy_sub_sym_block_factory {
  Matrix m;
  explicit mat_copy_sub_sym_block_factory(Matrix&& aM) : m(std::move(aM)) {}
  mat_copy_sub_sym_block<Matrix> operator()(const std::pair<int, int>& rows) {
    return mat_copy_sub_sym_block<Matrix>(std::move(m),
                                          rows.second - rows.first, rows.first);
  }
};

template <typename Matrix>
struct mat_sub_sym_block_factory {
  Matrix& m;
  explicit mat_sub_sym_block_factory(Matrix& aM) : m(aM) {}
  mat_sub_sym_block<Matrix> operator()(const std::pair<int, int>& rows) {
    return mat_sub_sym_block<Matrix>(m, rows.second - rows.first, rows.first);
  }
};

template <typename Matrix>
struct mat_const_sub_sym_block_factory {
  const Matrix& m;
  explicit mat_const_sub_sym_block_factory(const Matrix& aM) : m(aM) {}
  mat_const_sub_sym_block<Matrix> operator()(const std::pair<int, int>& rows) {
    return mat_const_sub_sym_block<Matrix>(m, rows.second - rows.first,
                                           rows.first);
  }
};

template <typename Matrix>
std::enable_if_t<is_readable_matrix_v<Matrix>,
                 mat_sub_sym_block_factory<Matrix>>
sub_sym(Matrix& M) {
  return mat_sub_sym_block_factory<Matrix>(M);
}

template <typename Matrix>
std::enable_if_t<is_readable_matrix_v<Matrix>,
                 mat_const_sub_sym_block_factory<Matrix>>
sub_sym(const Matrix& M) {
  return mat_const_sub_sym_block_factory<Matrix>(M);
}

template <typename Matrix>
std::enable_if_t<is_readable_matrix_v<Matrix>,
                 mat_copy_sub_sym_block_factory<Matrix>>
sub_sym_copy(const Matrix& M) {
  return mat_copy_sub_sym_block_factory<Matrix>(M);
}

template <typename Matrix>
std::enable_if_t<is_readable_matrix_v<Matrix>,
                 mat_copy_sub_sym_block_factory<Matrix>>
sub_sym(Matrix&& M) {
  return mat_copy_sub_sym_block_factory<Matrix>(std::move(M));
}
}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_VIEWS_H_

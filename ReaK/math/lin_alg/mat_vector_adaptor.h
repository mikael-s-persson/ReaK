/**
 * \file mat_vector_adaptor.h
 *
 * This library provides a number of classes which can be used to make a vector appear
 * as a matrix, with different alignments. These adaptors are useful when one has a vector
 * but wants to use a matrix method (i.e. solving a linear system of equations).
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

#ifndef REAK_MATH_LIN_ALG_MAT_VECTOR_ADAPTOR_H_
#define REAK_MATH_LIN_ALG_MAT_VECTOR_ADAPTOR_H_

#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/mat_slices.h"
#include "ReaK/math/lin_alg/mat_traits.h"
#include "ReaK/math/lin_alg/mat_views.h"
#include "ReaK/math/lin_alg/vect_concepts.h"
#include "ReaK/math/lin_alg/vect_views.h"

#include <type_traits>

namespace ReaK {

/**
 * This class template is the general template for a matrix-vector adaptor.
 * \tparam Vector A vector type (at least a readable vector type).
 * \tparam Alignment An enum describing the alignment of the matrix.
 */
template <WritableVector Vector,
          mat_alignment::tag Alignment = mat_alignment::column_major>
class mat_vect_adaptor {};

/**
 * This class template is a template for a matrix-vector adaptor with column-major alignment.
 * \tparam Vector A vector type (at least a readable vector type).
 */
template <WritableVector Vector>
class mat_vect_adaptor<Vector, mat_alignment::column_major> {
 public:
  using self = mat_vect_adaptor<Vector, mat_alignment::column_major>;

  using value_type = vect_value_type_t<Vector>;

  using reference = typename vect_traits<Vector>::reference;
  using const_reference = typename vect_traits<Vector>::const_reference;
  using pointer = typename vect_traits<Vector>::pointer;
  using const_pointer = typename vect_traits<Vector>::const_pointer;

  using col_iterator = typename vect_traits<Vector>::iterator;
  using const_col_iterator = typename vect_traits<Vector>::const_iterator;
  using row_iterator = typename vect_traits<Vector>::iterator;
  using const_row_iterator = typename vect_traits<Vector>::const_iterator;

  using size_type = typename vect_traits<Vector>::size_type;
  using difference_type = typename vect_traits<Vector>::difference_type;

  static constexpr unsigned int static_row_count =
      vect_traits<Vector>::dimensions;
  static constexpr unsigned int static_col_count = 1;
  static constexpr mat_alignment::tag alignment = mat_alignment::column_major;
  static constexpr mat_structure::tag structure = mat_structure::rectangular;

 private:
  Vector* v;           ///< Holds the reference to the vector.
  int offset;          ///< Holds the offset to start from in the vector.
  int rowCount;        ///< Holds the number of rows in the matrix.
  int colCount;        ///< Holds the number of columns in the matrix.
 public:
  /**
   * Constructs the adaptor with a given vector, taking the entire vector as the unique
   * column of a matrix.
   */
  explicit mat_vect_adaptor(Vector& aV)
      : v(&aV), offset(0), rowCount(aV.size()), colCount(1) {}

  /**
   * Constructs the adaptor with a given vector and parameters.
   * \param aV The vector to adapt to a matrix.
   * \param aRowCount The row-count of the resulting matrix.
   * \param aColCount The column-count of the resulting matrix.
   * \param aOffset The index into the vector from which to start the matrix elements.
   */
  mat_vect_adaptor(Vector& aV, int aRowCount, int aColCount, int aOffset = 0)
      : v(&aV), offset(aOffset), rowCount(aRowCount), colCount(aColCount) {}
  /**
   * Standard copy-constructor (shallow).
   */
  mat_vect_adaptor(const self& aObj)
      : v(aObj.v),
        offset(aObj.offset),
        rowCount(aObj.rowCount),
        colCount(aObj.colCount) {}

  /**
   * Standard swap function (shallow).
   */
  friend void swap(self& lhs, self& rhs) {
    using std::swap;
    swap(lhs.v, rhs.v);
    swap(lhs.offset, rhs.offset);
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
    int it = offset;
    for (int j = 0; j < colCount; ++j) {
      for (int i = 0; i < rowCount; ++i, ++it) {
        (*v)[it] = rhs(i, j);
      }
    }
    return *this;
  }

  /**
   * Standard assignment operator.
   */
  template <ReadableMatrix Matrix>
  self& operator=(const Matrix& rhs) {
    if ((rhs.get_row_count() != rowCount) ||
        (rhs.get_col_count() != colCount)) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    int it = offset;
    for (int j = 0; j < colCount; ++j) {
      for (int i = 0; i < rowCount; ++i, ++it) {
        (*v)[it] = rhs(i, j);
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
  reference operator()(int i, int j) { return (*v)[offset + j * rowCount + i]; }
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  const_reference operator()(int i, int j) const {
    return (*v)[offset + j * rowCount + i];
  }

  /**
   * Sub-matrix operator, accessor for read/write.
   * \test PASSED
   */
  mat_sub_block<self> operator()(const std::pair<int, int>& r,
                                 const std::pair<int, int>& c) {
    return sub(*this)(r, c);
  }

  /**
   * Sub-matrix operator, accessor for read only.
   * \test PASSED
   */
  mat_const_sub_block<self> operator()(const std::pair<int, int>& r,
                                       const std::pair<int, int>& c) const {
    return sub(*this)(r, c);
  }

  /**
   * Sub-matrix operator, accessor for read/write.
   * \test PASSED
   */
  mat_col_slice<self> operator()(int r, const std::pair<int, int>& c) {
    return slice(*this)(r, c);
  }

  /**
   * Sub-matrix operator, accessor for read only.
   * \test PASSED
   */
  mat_const_col_slice<self> operator()(int r,
                                       const std::pair<int, int>& c) const {
    return slice(*this)(r, c);
  }

  /**
   * Sub-matrix operator, accessor for read/write.
   * \test PASSED
   */
  mat_row_slice<self> operator()(const std::pair<int, int>& r, int c) {
    return slice(*this)(r, c);
  }

  /**
   * Sub-matrix operator, accessor for read only.
   * \test PASSED
   */
  mat_const_row_slice<self> operator()(const std::pair<int, int>& r,
                                       int c) const {
    return slice(*this)(r, c);
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
   * Sets the row-count (number of rows) of the matrix (however, it actually does nothing).
   */
  void set_row_count(int /*unused*/) noexcept {}
  /**
   * Sets the column-count (number of columns) of the matrix (however, it actually does nothing).
   */
  void set_col_count(int /*unused*/) noexcept {}

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept { return {rowCount, colCount}; }

  /** COL-MAJOR ONLY
*Add-and-store operator with standard semantics.
*\test PASSED
*/
  template <ReadableMatrix Matrix>
  self& operator+=(const Matrix& M) {
    if ((M.get_col_count() != colCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int it = offset;
    for (int j = 0; j < colCount; ++j) {
      for (int i = 0; i < rowCount; ++i, ++it) {
        (*v)[it] += M(i, j);
      }
    }
    return *this;
  }

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template <ReadableMatrix Matrix>
  self& operator-=(const Matrix& M) {
    if ((M.get_col_count() != colCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int it = offset;
    for (int j = 0; j < colCount; ++j) {
      for (int i = 0; i < rowCount; ++i, ++it) {
        (*v)[it] -= M(i, j);
      }
    }
    return *this;
  }

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=(const value_type& S) {
    for (int it = offset; it < offset + rowCount * colCount; ++it) {
      (*v)[it] *= S;
    }
    return *this;
  }

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template <ReadableMatrix Matrix>
  self& operator*=(const Matrix& M) {
    if ((M.get_col_count() != colCount) || (M.get_row_count() != colCount)) {
      throw std::range_error("Matrix Dimension Mismatch.");
    }
    *this = *this * M;
    return *this;
  }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend mat_vect_adaptor<Vector, mat_alignment::row_major> transpose(
      const self& M) {
    return mat_vect_adaptor<Vector, mat_alignment::row_major>(
        *M.v, M.offset, M.colCount, M.rowCount);
  }
};

template <WritableVector Vector>
class mat_vect_adaptor<Vector, mat_alignment::row_major> {
 public:
  using self = mat_vect_adaptor<Vector, mat_alignment::row_major>;

  using value_type = vect_value_type_t<Vector>;

  using reference = typename vect_traits<Vector>::reference;
  using const_reference = typename vect_traits<Vector>::const_reference;
  using pointer = typename vect_traits<Vector>::pointer;
  using const_pointer = typename vect_traits<Vector>::const_pointer;

  using col_iterator = typename vect_traits<Vector>::iterator;
  using const_col_iterator = typename vect_traits<Vector>::const_iterator;
  using row_iterator = typename vect_traits<Vector>::iterator;
  using const_row_iterator = typename vect_traits<Vector>::const_iterator;

  using size_type = typename vect_traits<Vector>::size_type;
  using difference_type = typename vect_traits<Vector>::difference_type;

  static constexpr unsigned int static_row_count = 1;
  static constexpr unsigned int static_col_count =
      vect_traits<Vector>::dimensions;
  static constexpr mat_alignment::tag alignment = mat_alignment::row_major;
  static constexpr mat_structure::tag structure = mat_structure::rectangular;

 private:
  Vector* v;           ///< Holds the reference to the vector.
  int offset;          ///< Holds the offset to start from in the vector.
  int rowCount;        ///< Holds the number of rows in the matrix.
  int colCount;        ///< Holds the number of columns in the matrix.
 public:
  /**
   * Constructs the adaptor with a given vector, taking the entire vector as the unique
   * column of a matrix.
   */
  explicit mat_vect_adaptor(Vector& aV)
      : v(&aV), offset(0), rowCount(1), colCount(aV.size()) {}

  /**
   * Constructs the adaptor with a given vector and parameters.
   * \param aV The vector to adapt to a matrix.
   * \param aRowCount The row-count of the resulting matrix.
   * \param aColCount The column-count of the resulting matrix.
   * \param aOffset The index into the vector from which to start the matrix elements.
   */
  mat_vect_adaptor(Vector& aV, int aRowCount, int aColCount, int aOffset = 0)
      : v(&aV), offset(aOffset), rowCount(aRowCount), colCount(aColCount) {}
  /**
   * Standard copy-constructor (shallow).
   */
  mat_vect_adaptor(const self& aObj)
      : v(aObj.v),
        offset(aObj.offset),
        rowCount(aObj.rowCount),
        colCount(aObj.colCount) {}

  /**
   * Standard swap function (shallow).
   */
  friend void swap(self& lhs, self& rhs) {
    using std::swap;
    swap(lhs.v, rhs.v);
    swap(lhs.offset, rhs.offset);
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
    int it = offset;
    for (int i = 0; i < rowCount; ++i) {
      for (int j = 0; j < colCount; ++j, ++it) {
        (*v)[it] = rhs(i, j);
      }
    }
    return *this;
  }

  /**
   * Standard assignment operator.
   */
  template <ReadableMatrix Matrix>
  self& operator=(const Matrix& rhs) {
    if ((rhs.get_row_count() != rowCount) ||
        (rhs.get_col_count() != colCount)) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
    int it = offset;
    for (int i = 0; i < rowCount; ++i) {
      for (int j = 0; j < colCount; ++j, ++it) {
        (*v)[it] = rhs(i, j);
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
  reference operator()(int i, int j) { return (*v)[offset + i * colCount + j]; }
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  const_reference operator()(int i, int j) const {
    return (*v)[offset + i * colCount + j];
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
   * Sets the row-count (number of rows) of the matrix (however, it actually does nothing).
   */
  void set_row_count(int /*unused*/) noexcept {}
  /**
   * Sets the column-count (number of columns) of the matrix (however, it actually does nothing).
   */
  void set_col_count(int /*unused*/) noexcept {}

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<int, int> size() const noexcept {
    return std::make_pair(rowCount, colCount);
  }

  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template <ReadableMatrix Matrix>
  self& operator+=(const Matrix& M) {
    if ((M.get_col_count() != colCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int it = offset;
    for (int i = 0; i < rowCount; ++i) {
      for (int j = 0; j < colCount; ++j, ++it) {
        (*v)[it] += M(i, j);
      }
    }
    return *this;
  }

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template <ReadableMatrix Matrix>
  self& operator-=(const Matrix& M) {
    if ((M.get_col_count() != colCount) || (M.get_row_count() != rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int it = offset;
    for (int i = 0; i < rowCount; ++i) {
      for (int j = 0; j < colCount; ++j, ++it) {
        (*v)[it] -= M(i, j);
      }
    }
    return *this;
  }

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=(const value_type& S) {
    for (int it = offset; it < offset + rowCount * colCount; ++it) {
      (*v)[it] *= S;
    }
    return *this;
  }

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template <ReadableMatrix Matrix>
  self& operator*=(const Matrix& M) {
    if ((M.get_col_count() != colCount) || (M.get_row_count() != colCount)) {
      throw std::range_error("Matrix Dimension Mismatch.");
    }
    *this = *this * M;
    return *this;
  }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend mat_vect_adaptor<Vector, mat_alignment::column_major> transpose(
      const self& M) {
    return mat_vect_adaptor<Vector, mat_alignment::column_major>(
        *M.v, M.offset, M.colCount, M.rowCount);
  }
};

template <typename Vector, mat_alignment::tag Alignment>
static constexpr bool is_fully_writable_matrix_v<mat_vect_adaptor<Vector, Alignment>> = WritableVector<Vector>;

template <ReadableVector Vector,
          mat_alignment::tag Alignment = mat_alignment::column_major>
class mat_const_vect_adaptor {};

template <ReadableVector Vector>
class mat_const_vect_adaptor<Vector, mat_alignment::column_major> {
 public:
  using self = mat_const_vect_adaptor<Vector, mat_alignment::column_major>;

  using value_type = vect_value_type_t<Vector>;

  using reference = typename vect_traits<Vector>::reference;
  using const_reference = typename vect_traits<Vector>::const_reference;
  using pointer = typename vect_traits<Vector>::pointer;
  using const_pointer = typename vect_traits<Vector>::const_pointer;

  using col_iterator = typename vect_traits<Vector>::iterator;
  using const_col_iterator = typename vect_traits<Vector>::const_iterator;
  using row_iterator = typename vect_traits<Vector>::iterator;
  using const_row_iterator = typename vect_traits<Vector>::const_iterator;

  using size_type = typename vect_traits<Vector>::size_type;
  using difference_type = typename vect_traits<Vector>::difference_type;

  static constexpr unsigned int static_row_count =
      vect_traits<Vector>::dimensions;
  static constexpr unsigned int static_col_count = 1;
  static constexpr mat_alignment::tag alignment = mat_alignment::column_major;
  static constexpr mat_structure::tag structure = mat_structure::rectangular;

 private:
  const Vector* v;     ///< Holds the reference to the vector.
  int offset;          ///< Holds the offset to start from in the vector.
  int rowCount;        ///< Holds the number of rows in the matrix.
  int colCount;        ///< Holds the number of columns in the matrix.

  self& operator=(const self&);  // non-assignable.

  explicit mat_const_vect_adaptor(Vector&&);
  mat_const_vect_adaptor(Vector&&, int, int, int aOffset = 0);

 public:
  /**
   * Constructs the adaptor with a given vector, taking the entire vector as the unique
   * column of a matrix.
   */
  explicit mat_const_vect_adaptor(const Vector& aV)
      : v(&aV), offset(0), rowCount(aV.size()), colCount(1) {}

  /**
   * Constructs the adaptor with a given vector and parameters.
   * \param aV The vector to adapt to a matrix.
   * \param aRowCount The row-count of the resulting matrix.
   * \param aColCount The column-count of the resulting matrix.
   * \param aOffset The index into the vector from which to start the matrix elements.
   */
  mat_const_vect_adaptor(const Vector& aV, int aRowCount, int aColCount,
                         int aOffset = 0)
      : v(&aV), offset(aOffset), rowCount(aRowCount), colCount(aColCount) {}
  /**
   * Standard copy-constructor (shallow).
   */
  mat_const_vect_adaptor(const self& aObj)
      : v(aObj.v),
        offset(aObj.offset),
        rowCount(aObj.rowCount),
        colCount(aObj.colCount) {}

  /**
   * Standard swap function (shallow).
   */
  friend void swap(self& lhs, self& rhs) {
    using std::swap;
    swap(lhs.v, rhs.v);
    swap(lhs.offset, rhs.offset);
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
  const_reference operator()(int i, int j) const {
    return (*v)[offset + j * rowCount + i];
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

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend mat_const_vect_adaptor<Vector, mat_alignment::row_major> transpose(
      const self& M) {
    return mat_const_vect_adaptor<Vector, mat_alignment::row_major>(
        *M.v, M.offset, M.colCount, M.rowCount);
  }
};

template <ReadableVector Vector>
class mat_const_vect_adaptor<Vector, mat_alignment::row_major> {
 public:
  using self = mat_const_vect_adaptor<Vector, mat_alignment::row_major>;

  using value_type = vect_value_type_t<Vector>;

  using reference = typename vect_traits<Vector>::reference;
  using const_reference = typename vect_traits<Vector>::const_reference;
  using pointer = typename vect_traits<Vector>::pointer;
  using const_pointer = typename vect_traits<Vector>::const_pointer;

  using col_iterator = typename vect_traits<Vector>::iterator;
  using const_col_iterator = typename vect_traits<Vector>::const_iterator;
  using row_iterator = typename vect_traits<Vector>::iterator;
  using const_row_iterator = typename vect_traits<Vector>::const_iterator;

  using size_type = typename vect_traits<Vector>::size_type;
  using difference_type = typename vect_traits<Vector>::difference_type;

  static constexpr unsigned int static_row_count = 1;
  static constexpr unsigned int static_col_count =
      vect_traits<Vector>::dimensions;
  static constexpr mat_alignment::tag alignment = mat_alignment::row_major;
  static constexpr mat_structure::tag structure = mat_structure::rectangular;

 private:
  const Vector* v;
  int offset;
  int rowCount;
  int colCount;

  self& operator=(const self&);  // non-assignable.

  explicit mat_const_vect_adaptor(Vector&&);
  mat_const_vect_adaptor(Vector&&, int, int, int aOffset = 0);

 public:
  /**
   * Constructs the adaptor with a given vector, taking the entire vector as the unique
   * column of a matrix.
   */
  explicit mat_const_vect_adaptor(const Vector& aV)
      : v(&aV), offset(0), rowCount(1), colCount(aV.size()) {}

  /**
   * Constructs the adaptor with a given vector and parameters.
   * \param aV The vector to adapt to a matrix.
   * \param aRowCount The row-count of the resulting matrix.
   * \param aColCount The column-count of the resulting matrix.
   * \param aOffset The index into the vector from which to start the matrix elements.
   */
  mat_const_vect_adaptor(const Vector& aV, int aRowCount, int aColCount,
                         int aOffset = 0)
      : v(&aV), offset(aOffset), rowCount(aRowCount), colCount(aColCount) {}
  /**
   * Standard copy-constructor (shallow).
   */
  mat_const_vect_adaptor(const self& aObj)
      : v(aObj.v),
        offset(aObj.offset),
        rowCount(aObj.rowCount),
        colCount(aObj.colCount) {}

  /**
   * Standard swap function (shallow).
   */
  friend void swap(self& lhs, self& rhs) {
    using std::swap;
    swap(lhs.v, rhs.v);
    swap(lhs.offset, rhs.offset);
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
  const_reference operator()(int i, int j) const {
    return (*v)[offset + i * colCount + j];
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

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend mat_const_vect_adaptor<Vector, mat_alignment::column_major> transpose(
      const self& M) {
    return mat_const_vect_adaptor<Vector, mat_alignment::column_major>(
        *M.v, M.offset, M.colCount, M.rowCount);
  }
};

template <WritableVector Vector>
struct mat_vect_adaptor_factory {
  Vector& v;
  explicit mat_vect_adaptor_factory(Vector& aV) : v(aV) {}
  mat_vect_adaptor<Vector, mat_alignment::row_major> operator()(
      int rowCount, const std::pair<int, int>& cols) {
    return mat_vect_adaptor<Vector, mat_alignment::row_major>(
        v, rowCount, cols.second - cols.first, cols.first);
  }
  mat_vect_adaptor<Vector, mat_alignment::column_major> operator()(
      const std::pair<int, int>& rows, int colCount) {
    return mat_vect_adaptor<Vector, mat_alignment::column_major>(
        v, rows.second - rows.first, colCount, rows.first);
  }
};

template <ReadableVector Vector>
struct mat_const_vect_adaptor_factory {
  const Vector& v;
  explicit mat_const_vect_adaptor_factory(const Vector& aV) : v(aV) {}
  mat_const_vect_adaptor<Vector, mat_alignment::row_major> operator()(
      int rowCount, const std::pair<int, int>& cols) {
    return mat_const_vect_adaptor<Vector, mat_alignment::row_major>(
        v, rowCount, cols.second - cols.first, cols.first);
  }
  mat_const_vect_adaptor<Vector, mat_alignment::column_major> operator()(
      const std::pair<int, int>& rows, int colCount) {
    return mat_const_vect_adaptor<Vector, mat_alignment::column_major>(
        v, rows.second - rows.first, colCount, rows.first);
  }
};

/**
 * This function template can be used to generate a matrix-vector adaptor. It constructs
 * an object of a factory functor which can generate either a column-major or row-major matrix
 * depending on the call syntax. The matrix is constructed via a range function as so:
 *
 *  make_mat(V)(range(3,7),2)  generates a column-major matrix starting from index 3, with 4 rows (index 3 to 6) and 2
 *columns.
 *
 *  make_mat(V)(2,range(3,7))  generates a row-major matrix starting from index 3, with 4 columns (index 3 to 6) and 2
 *rows.
 *
 * \tparam Vector A readable vector type.
 */
template <WritableVector Vector>
mat_vect_adaptor_factory<Vector> make_mat(Vector& V) {
  return mat_vect_adaptor_factory<Vector>(V);
}

/**
 * This function template can be used to generate a matrix-vector adaptor. It constructs
 * an object of a factory functor which can generate either a column-major or row-major matrix
 * depending on the call syntax. The matrix is constructed via a range function as so:
 *
 *  make_mat(V)(range(3,6),2)  generates a column-major matrix starting from index 3, with 4 rows (index 3 to 6) and 2
 *columns.
 *
 *  make_mat(V)(2,range(3,6))  generates a row-major matrix starting from index 3, with 4 columns (index 3 to 6) and 2
 *rows.
 *
 * \tparam Vector A readable vector type.
 */
template <ReadableVector Vector>
mat_const_vect_adaptor_factory<Vector> make_mat(const Vector& V) {
  return mat_const_vect_adaptor_factory<Vector>(V);
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_VECTOR_ADAPTOR_H_

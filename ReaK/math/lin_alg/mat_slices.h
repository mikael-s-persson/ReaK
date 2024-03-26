/**
 * \file mat_slices.h
 *
 * This library provides a number of classes and functions to take a slice of a matrix.
 * A slice is simply a reduction of the order (i.e. a matrix is a 2nd order tensor,
 * and a vector is a 1st order tensor). So, this library allows one to extract a row
 * or column from a matrix and have this row or column appear as a vector. In other words,
 * this library provides adaptors from matrices to vectors.
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

#ifndef REAK_MATH_LIN_ALG_MAT_SLICES_H_
#define REAK_MATH_LIN_ALG_MAT_SLICES_H_

#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/mat_traits.h"
#include "ReaK/math/lin_alg/vect_concepts.h"
#include "ReaK/math/lin_alg/vect_views.h"

#include "ReaK/math/lin_alg/vect_alg.h"

#include <type_traits>

namespace ReaK {

/// This class template can be used to view a column of a matrix as a vector (i.e. a row-slice). It
/// takes the matrix object by reference (internally by pointer).
///
/// Models: ReadableVector and WritableVector (if the matrix is writable)
/// \tparam Matrix A matrix type.
template <typename Matrix>
class mat_row_slice {
 public:
  using self = mat_row_slice<Matrix>;

  using value_type = mat_value_type_t<Matrix>;

  using reference = typename mat_traits<Matrix>::reference;
  using const_reference = typename mat_traits<Matrix>::const_reference;
  using pointer = typename mat_traits<Matrix>::pointer;
  using const_pointer = typename mat_traits<Matrix>::const_pointer;

  using iterator = typename mat_traits<Matrix>::row_iterator;
  using const_iterator = typename mat_traits<Matrix>::const_row_iterator;

  using size_type = mat_size_type_t<Matrix>;
  using difference_type = typename mat_traits<Matrix>::difference_type;

  static constexpr std::size_t dimensions = 0;

 private:
  Matrix* m;     ///< Holds the reference to the matrix object.
  int offset;    ///< Holds the offset from the start of the column.
  int colIndex;  ///< Holds the index of column of the slice.
  int count;     ///< Holds the number of elements of the column to take.
 public:
  /// Constructs the row-slice from a matrix M, taking the entire first column.
  /// \param aM The matrix from which to take the slice.
  explicit mat_row_slice(Matrix& aM)
      : m(&aM), offset(0), colIndex(0), count(aM.get_row_count()) {}

  /// Constructs the row-slice from a matrix M, taking the aColIndex column
  /// from aOffset with aCount elements.
  /// \param aM The matrix from which to take the slice.
  /// \param aColIndex The column to use for the slice.
  /// \param aOffset The offset into the column used for the slice.
  /// \param aCount The number of elements to take in the slice.
  mat_row_slice(Matrix& aM, int aColIndex, int aOffset, int aCount)
      : m(&aM), offset(aOffset), colIndex(aColIndex), count(aCount) {}

  /// Standard copy-constructor (shallow).
  mat_row_slice(const self& rhs) = default;

  /// Standard swap function (shallow)
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.offset, rhs.offset);
    swap(lhs.colIndex, rhs.colIndex);
    swap(lhs.count, rhs.count);
  }

  /// Standard assignment operator for any readable vector type.
  template <ReadableVector Vector>
  self& operator=(const Vector& v) {
    if (v.size() != count) {
      throw std::range_error("Vector dimensions mismatch.");
    }
    for (int i = 0; i < count; ++i) {
      (*m)(offset + i, colIndex) = v[i];
    }
    return *this;
  }

  self& operator=(const self& v) {
    if (v.size() != count) {
      throw std::range_error("Vector dimensions mismatch.");
    }
    for (int i = 0; i < count; ++i) {
      (*m)(offset + i, colIndex) = v[i];
    }
    return *this;
  }

  /// Standard assignment operator for any readable vector type.
  template <ReadableVector Vector>
  self& operator+=(const Vector& v) {
    if (v.size() != count) {
      throw std::range_error("Vector dimensions mismatch.");
    }
    for (int i = 0; i < count; ++i) {
      (*m)(offset + i, colIndex) += v[i];
    }
    return *this;
  }

  /// Standard assignment operator for any readable vector type.
  template <ReadableVector Vector>
  self& operator-=(const Vector& v) {
    if (v.size() != count) {
      throw std::range_error("Vector dimensions mismatch.");
    }
    for (int i = 0; i < count; ++i) {
      (*m)(offset + i, colIndex) -= v[i];
    }
    return *this;
  }

  /// Returns the size of the vector.
  int size() const { return count; }
  /// Returns the max-size of the vector.
  int max_size() const { return count; }
  /// Returns the capacity of the vector.
  int capacity() const { return count; }
  /// Resizes the vector.
  /// \param sz The new size for the vector.
  /// \param c The value to fill any additional elements to the vector.
  void resize(int sz, const_reference c = value_type()) const {
    RK_UNUSED(sz);
    RK_UNUSED(c);
  }
  /// Checks whether the vector is empty.
  bool empty() const { return false; }
  /// Reserve a certain amount of capacity for future additions.
  /// \param sz The new capacity for the vector.
  void reserve(int sz) const { RK_UNUSED(sz); }

  /// Returns an iterator to the first element of the vector.
  iterator begin() { return m->first_row(m->first_col() + colIndex) + offset; }
  /// Returns an const-iterator to the first element of the vector.
  const_iterator begin() const {
    return m->first_row(m->first_col() + colIndex) + offset;
  }
  /// Returns an iterator to the one-passed-last element of the vector.
  iterator end() {
    return m->first_row(m->first_col() + colIndex) + offset + count;
  }
  /// Returns an const-iterator to the one-passed-last element of the vector.
  const_iterator end() const {
    return m->first_row(m->first_col() + colIndex) + offset + count;
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Array indexing operator, accessor for read/write.
  reference operator[](int i) { return (*m)(offset + i, colIndex); }

  /// Array indexing operator, accessor for read only.
  const_reference operator[](int i) const { return (*m)(offset + i, colIndex); }

  /// Array indexing operator, accessor for read/write.
  reference operator()(int i) { return (*m)(offset + i, colIndex); }

  /// Array indexing operator, accessor for read only.
  const_reference operator()(int i) const { return (*m)(offset + i, colIndex); }
};

template <typename Matrix>
struct vect_copy<mat_row_slice<Matrix>> {
  using type = vect_n<vect_value_type_t<mat_row_slice<Matrix>>>;
};

/// This class template can be used to view a column of a matrix as a const-vector (i.e. a row-slice). It
/// takes the matrix object by const-reference (internally by pointer).
///
/// Models: ReadableVector
/// \tparam Matrix A matrix type.
template <typename Matrix>
class mat_const_row_slice {
 public:
  using self = mat_const_row_slice<Matrix>;

  using value_type = mat_value_type_t<Matrix>;

  using reference = typename mat_traits<Matrix>::reference;
  using const_reference = typename mat_traits<Matrix>::const_reference;
  using pointer = typename mat_traits<Matrix>::pointer;
  using const_pointer = typename mat_traits<Matrix>::const_pointer;

  using iterator = typename mat_traits<Matrix>::row_iterator;
  using const_iterator = typename mat_traits<Matrix>::const_row_iterator;

  using size_type = mat_size_type_t<Matrix>;
  using difference_type = typename mat_traits<Matrix>::difference_type;

  static constexpr std::size_t dimensions = 0;

 private:
  const Matrix* m;  ///< Holds the reference to the matrix object.
  int offset;       ///< Holds the offset from the start of the column.
  int colIndex;     ///< Holds the index of column of the slice.
  int count;        ///< Holds the number of elements of the column to take.

  explicit mat_const_row_slice(Matrix&&);
  mat_const_row_slice(Matrix&&, int, int, int);

 public:
  /// Constructs the row-slice from a matrix M, taking the entire first column.
  /// \param aM The matrix from which to take the slice.
  explicit mat_const_row_slice(const Matrix& aM)
      : m(&aM), offset(0), colIndex(0), count(aM.get_row_count()) {}

  /// Constructs the row-slice from a matrix M, taking the aColIndex column
  /// from aOffset with aCount elements.
  /// \param aM The matrix from which to take the slice.
  /// \param aColIndex The column to use for the slice.
  /// \param aOffset The offset into the column used for the slice.
  /// \param aCount The number of elements to take in the slice.
  mat_const_row_slice(const Matrix& aM, int aColIndex, int aOffset, int aCount)
      : m(&aM), offset(aOffset), colIndex(aColIndex), count(aCount) {}

  /// Standard copy-constructor (shallow).
  mat_const_row_slice(const self& rhs) = default;

  self& operator=(const self& rhs) = delete;

  /// Standard swap function (shallow)
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.offset, rhs.offset);
    swap(lhs.colIndex, rhs.colIndex);
    swap(lhs.count, rhs.count);
  }

  /// Returns the size of the vector.
  int size() const { return count; }
  /// Returns the max-size of the vector.
  int max_size() const { return count; }
  /// Returns the capacity of the vector.
  int capacity() const { return count; }
  /// Resizes the vector.
  /// \param sz The new size for the vector.
  /// \param c The value to fill any additional elements to the vector.
  void resize(int sz, const_reference c = value_type()) const {}
  /// Checks whether the vector is empty.
  bool empty() const { return false; }
  /// Reserve a certain amount of capacity for future additions.
  /// \param sz The new capacity for the vector.
  void reserve(int sz) const {}

  /// Returns an const-iterator to the first element of the vector.
  const_iterator begin() const {
    return m->first_row(m->first_col() + colIndex) + offset;
  }
  /// Returns an const-iterator to the one-passed-last element of the vector.
  const_iterator end() const {
    return m->first_row(m->first_col() + colIndex) + offset + count;
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Array indexing operator, accessor for read only.
  const_reference operator[](int i) const { return (*m)(offset + i, colIndex); }

  /// Array indexing operator, accessor for read only.
  const_reference operator()(int i) const { return (*m)(offset + i, colIndex); }
};

template <typename Matrix>
struct vect_copy<mat_const_row_slice<Matrix>> {
  using type = vect_n<vect_value_type_t<mat_const_row_slice<Matrix>>>;
};

/// This class template can be used to view a row of a matrix as a vector (i.e. a column-slice). It
/// takes the matrix object by reference (internally by pointer).
///
/// Models: ReadableVector and WritableVector (if the matrix is writable)
/// \tparam Matrix A matrix type.
template <typename Matrix>
class mat_col_slice {
 public:
  using self = mat_col_slice<Matrix>;

  using value_type = mat_value_type_t<Matrix>;

  using reference = typename mat_traits<Matrix>::reference;
  using const_reference = typename mat_traits<Matrix>::const_reference;
  using pointer = typename mat_traits<Matrix>::pointer;
  using const_pointer = typename mat_traits<Matrix>::const_pointer;

  using iterator = typename mat_traits<Matrix>::col_iterator;
  using const_iterator = typename mat_traits<Matrix>::const_col_iterator;

  using size_type = mat_size_type_t<Matrix>;
  using difference_type = typename mat_traits<Matrix>::difference_type;

  static constexpr std::size_t dimensions = 0;

 private:
  Matrix* m;     ///< Holds the reference to the matrix object.
  int offset;    ///< Holds the offset from the start of the row.
  int rowIndex;  ///< Holds the index of row of the slice.
  int count;     ///< Holds the number of elements of the row to take.
 public:
  /// Constructs the column-slice from a matrix M, taking the entire first row.
  /// \param aM The matrix from which to take the slice.
  explicit mat_col_slice(Matrix& aM)
      : m(&aM), offset(0), rowIndex(0), count(aM.get_col_count()) {}

  /// Constructs the column-slice from a matrix M, taking the aRowIndex row
  /// from aOffset with aCount elements.
  /// \param aM The matrix from which to take the slice.
  /// \param aRowIndex The row to use for the slice.
  /// \param aOffset The offset into the row used for the slice.
  /// \param aCount The number of elements to take in the slice.
  mat_col_slice(Matrix& aM, int aRowIndex, int aOffset, int aCount)
      : m(&aM), offset(aOffset), rowIndex(aRowIndex), count(aCount) {}

  /// Standard copy-constructor (shallow).
  mat_col_slice(const self& rhs) = default;

  /// Standard swap function (shallow)
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.offset, rhs.offset);
    swap(lhs.rowIndex, rhs.rowIndex);
    swap(lhs.count, rhs.count);
  }

  /// Standard assignment operator for any readable vector type.
  template <ReadableVector Vector>
  self& operator=(const Vector& v) {
    if (v.size() != count) {
      throw std::range_error("Vector dimensions mismatch.");
    }
    for (int i = 0; i < count; ++i) {
      (*m)(rowIndex, offset + i) = v[i];
    }
    return *this;
  }

  self& operator=(const self& v) {
    if (v.size() != count) {
      throw std::range_error("Vector dimensions mismatch.");
    }
    for (int i = 0; i < count; ++i) {
      (*m)(rowIndex, offset + i) = v[i];
    }
    return *this;
  }

  /// Standard assignment operator for any readable vector type.
  /// \tparam Vector A readable vector type.
  template <ReadableVector Vector>
  self& operator+=(const Vector& v) {
    if (v.size() != count) {
      throw std::range_error("Vector dimensions mismatch.");
    }
    for (int i = 0; i < count; ++i) {
      (*m)(rowIndex, offset + i) += v[i];
    }
    return *this;
  }

  /// Standard assignment operator for any readable vector type.
  /// \tparam Vector A readable vector type.
  template <ReadableVector Vector>
  self& operator-=(const Vector& v) {
    if (v.size() != count) {
      throw std::range_error("Vector dimensions mismatch.");
    }
    for (int i = 0; i < count; ++i) {
      (*m)(rowIndex, offset + i) -= v[i];
    }
    return *this;
  }

  /// Returns the size of the vector.
  int size() const { return count; }
  /// Returns the max-size of the vector.
  int max_size() const { return count; }
  /// Returns the capacity of the vector.
  int capacity() const { return count; }
  /// Resizes the vector.
  /// \param sz The new size for the vector.
  /// \param c The value to fill any additional elements to the vector.
  void resize(int sz, const_reference c = value_type()) const {}
  /// Checks whether the vector is empty.
  bool empty() const { return false; }
  /// Reserve a certain amount of capacity for future additions.
  /// \param sz The new capacity for the vector.
  void reserve(int sz) const {}

  /// Returns an iterator to the first element of the vector.
  iterator begin() { return m->first_col(m->first_row() + rowIndex) + offset; }
  /// Returns an const-iterator to the first element of the vector.
  const_iterator begin() const {
    return m->first_col(m->first_row() + rowIndex) + offset;
  }
  /// Returns an iterator to the one-passed-last element of the vector.
  iterator end() {
    return m->first_col(m->first_row() + rowIndex) + offset + count;
  }
  /// Returns an const-iterator to the one-passed-last element of the vector.
  const_iterator end() const {
    return m->first_col(m->first_row() + rowIndex) + offset + count;
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Array indexing operator, accessor for read/write.
  reference operator[](int i) { return (*m)(rowIndex, offset + i); }

  /// Array indexing operator, accessor for read only.
  const_reference operator[](int i) const { return (*m)(rowIndex, offset + i); }

  /// Array indexing operator, accessor for read/write.
  reference operator()(int i) { return (*m)(rowIndex, offset + i); }

  /// Array indexing operator, accessor for read only.
  const_reference operator()(int i) const { return (*m)(rowIndex, offset + i); }
};

template <typename Matrix>
struct vect_copy<mat_col_slice<Matrix>> {
  using type = vect_n<vect_value_type_t<mat_col_slice<Matrix>>>;
};

/// This class template can be used to view a row of a matrix as a vector (i.e. a column-slice). It
/// takes the matrix object by const-reference (internally by const-pointer).
///
/// Models: ReadableVector
/// \tparam Matrix A matrix type.
template <typename Matrix>
class mat_const_col_slice {
 public:
  using self = mat_const_col_slice<Matrix>;

  using value_type = mat_value_type_t<Matrix>;

  using reference = typename mat_traits<Matrix>::reference;
  using const_reference = typename mat_traits<Matrix>::const_reference;
  using pointer = typename mat_traits<Matrix>::pointer;
  using const_pointer = typename mat_traits<Matrix>::const_pointer;

  using iterator = typename mat_traits<Matrix>::col_iterator;
  using const_iterator = typename mat_traits<Matrix>::const_col_iterator;

  using size_type = mat_size_type_t<Matrix>;
  using difference_type = typename mat_traits<Matrix>::difference_type;

  static constexpr std::size_t dimensions = 0;

 private:
  const Matrix* m;  ///< Holds the reference to the matrix object.
  int offset;       ///< Holds the offset from the start of the row.
  int rowIndex;     ///< Holds the index of row of the slice.
  int count;        ///< Holds the number of elements of the row to take.

  explicit mat_const_col_slice(Matrix&&);
  mat_const_col_slice(Matrix&&, int, int, int);

 public:
  /// Constructs the column-slice from a matrix M, taking the entire first row.
  /// \param aM The matrix from which to take the slice.
  explicit mat_const_col_slice(const Matrix& aM)
      : m(&aM), offset(0), rowIndex(0), count(aM.get_col_count()) {}

  /// Constructs the column-slice from a matrix M, taking the aRowIndex row
  /// from aOffset with aCount elements.
  /// \param aM The matrix from which to take the slice.
  /// \param aRowIndex The row to use for the slice.
  /// \param aOffset The offset into the row used for the slice.
  /// \param aCount The number of elements to take in the slice.
  mat_const_col_slice(const Matrix& aM, int aRowIndex, int aOffset, int aCount)
      : m(&aM), offset(aOffset), rowIndex(aRowIndex), count(aCount) {}

  /// Standard copy-constructor (shallow).
  mat_const_col_slice(const self& rhs) = default;

  self& operator=(const self& rhs) = delete;

  /// Standard swap function (shallow)
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.m, rhs.m);
    swap(lhs.offset, rhs.offset);
    swap(lhs.rowIndex, rhs.rowIndex);
    swap(lhs.count, rhs.count);
  }

  /// Returns the size of the vector.
  int size() const { return count; }
  /// Returns the max-size of the vector.
  int max_size() const { return count; }
  /// Returns the capacity of the vector.
  int capacity() const { return count; }
  /// Resizes the vector.
  /// \param sz The new size for the vector.
  /// \param c The value to fill any additional elements to the vector.
  void resize(int sz, const_reference c = value_type()) const {}
  /// Checks whether the vector is empty.
  bool empty() const { return false; }
  /// Reserve a certain amount of capacity for future additions.
  /// \param sz The new capacity for the vector.
  void reserve(int sz) const {}

  /// Returns an const-iterator to the first element of the vector.
  const_iterator begin() const {
    return m->first_col(m->first_row() + rowIndex) + offset;
  }
  /// Returns an const-iterator to the one-passed-last element of the vector.
  const_iterator end() const {
    return m->first_col(m->first_row() + rowIndex) + offset + count;
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Array indexing operator, accessor for read only.
  const_reference operator[](int i) const { return (*m)(rowIndex, offset + i); }

  /// Array indexing operator, accessor for read only.
  const_reference operator()(int i) const { return (*m)(rowIndex, offset + i); }
};

template <ReadableMatrix Matrix>
struct vect_copy<mat_const_col_slice<Matrix>> {
  using type = vect_n<vect_value_type_t<mat_const_col_slice<Matrix>>>;
};

template <WritableMatrix Matrix>
struct mat_slice_factory {
  Matrix& m;
  explicit mat_slice_factory(Matrix& aM) : m(aM) {}
  mat_col_slice<Matrix> operator()(int row, const std::pair<int, int>& cols) {
    return mat_col_slice<Matrix>(m, row, cols.first, cols.second - cols.first);
  }
  mat_row_slice<Matrix> operator()(const std::pair<int, int>& rows, int col) {
    return mat_row_slice<Matrix>(m, col, rows.first, rows.second - rows.first);
  }
};

template <ReadableMatrix Matrix>
struct mat_const_slice_factory {
  const Matrix& m;
  explicit mat_const_slice_factory(const Matrix& aM) : m(aM) {}
  mat_const_col_slice<Matrix> operator()(int row,
                                         const std::pair<int, int>& cols) {
    return mat_const_col_slice<Matrix>(m, row, cols.first,
                                       cols.second - cols.first);
  }
  mat_const_row_slice<Matrix> operator()(const std::pair<int, int>& rows,
                                         int col) {
    return mat_const_row_slice<Matrix>(m, col, rows.first,
                                       rows.second - rows.first);
  }
};

/// This function template prepares a matrix to be sliced via an intermediate factor class.
/// Once the factor class has been created with this function, one can use the range() function
/// to define to row or column range to use in the slice. Given matrix M, one can create
/// a slice as follow: slice(M)(range(0,M.get_row_count()-1), 0) which creates a vector view
/// on the entire first column.
/// \tparam Matrix A readable matrix type.
/// \param M the matrix to slice.
/// \return A factory class to create a slice from a row or column range and a row or column index.
template <WritableMatrix Matrix>
mat_slice_factory<Matrix> slice(Matrix& M) {
  return mat_slice_factory<Matrix>(M);
}

/// This function template prepares a matrix to be sliced via an intermediate factor class.
/// Once the factor class has been created with this function, one can use the range() function
/// to define to row or column range to use in the slice. Given matrix M, one can create
/// a slice as follow: slice(M)(range(0,M.get_row_count()-1), 0) which creates a vector view
/// on the entire first column.
/// \tparam Matrix A readable matrix type.
/// \param M the matrix to slice.
/// \return A factory class to create a slice from a row or column range and a row or column index.
template <ReadableMatrix Matrix>
mat_const_slice_factory<Matrix> slice(const Matrix& M) {
  return mat_const_slice_factory<Matrix>(M);
}
}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_SLICES_H_

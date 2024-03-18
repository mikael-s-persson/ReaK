/**
 * \file mat_alg_square.h
 *
 * This library implements the specialization of the mat<> template for a
 * general square matrix (dynamic dimension) of both column-major and
 * row-major alignment. This matrix type fulfills all the general matrix
 * concepts (Readable, Writable, Fully-Writable, and Resizable).
 *
 * This library also implements transposition of matrices via alignment
 * switching (switching from column-major to row-major, or vice versa). This
 * is very efficient and can even avoid copies completely (via transpose(std::move())) on an
 * optimizing compiler.
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

#ifndef REAK_MATH_LIN_ALG_MAT_ALG_SQUARE_H_
#define REAK_MATH_LIN_ALG_MAT_ALG_SQUARE_H_

#include "ReaK/math/lin_alg/mat_alg_general.h"

#include <type_traits>
#include <utility>

namespace ReaK {

/// This class template specialization implements a matrix with square structure
/// and column-major alignment. This class is serializable and registered to the ReaK::rtti
/// system. This matrix type is dynamically resizable.
///
/// Models: ReadableMatrix, WritableMatrix, FullyWritableMatrix,
/// and ResizableMatrix.
///
/// \tparam T Arithmetic type of the elements of the matrix.
template <typename T, unsigned int RowCount>
class mat<T, mat_structure::square, mat_alignment::column_major, RowCount,
          RowCount> : public serializable {
 public:
  using self = mat<T, mat_structure::square, mat_alignment::column_major,
                   RowCount, RowCount>;
  static constexpr bool is_dynamic_size = (RowCount == 0);

  using value_type = T;
  using container_type =
      std::conditional_t<is_dynamic_size, std::vector<value_type>,
                         std::array<value_type, RowCount * RowCount>>;

  using reference = typename container_type::reference;
  using const_reference = typename container_type::const_reference;
  using pointer = typename container_type::pointer;
  using const_pointer = typename container_type::const_pointer;

  using col_iterator = stride_iterator<typename container_type::iterator>;
  using const_col_iterator =
      stride_iterator<typename container_type::const_iterator>;
  using row_iterator = typename container_type::iterator;
  using const_row_iterator = typename container_type::const_iterator;

  using size_type = int;
  using difference_type = int;

  static constexpr unsigned int static_row_count = RowCount;
  static constexpr unsigned int static_col_count = RowCount;
  static constexpr mat_alignment::tag alignment = mat_alignment::column_major;
  static constexpr mat_structure::tag structure = mat_structure::square;

 private:
  struct dynamic_data {
    container_type q = {};
    int rowCount = 0;
  };
  struct static_data {
    container_type q = {};
  };
  /// Array which holds all the values of the matrix (dimension: rowCount x rowCount).
  std::conditional_t<is_dynamic_size, dynamic_data, static_data> data;

 public:
  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /// Default constructor: sets all to zero.
  mat() : data() {}

  /// Constructor for a sized matrix.
  explicit mat(int aRowCount, const value_type& aFill = value_type()) {
    if constexpr (is_dynamic_size) {
      data.q.resize(aRowCount * aRowCount, aFill);
      data.rowCount = aRowCount;
    } else {
      if (aRowCount != RowCount) {
        throw std::range_error("Row count mismatch!");
      }
      std::fill(data.q.begin(), data.q.end(), aFill);
    }
  }

  mat(int aRowCount, int aColCount, const value_type& aFill = value_type())
      : mat(aRowCount, aFill) {
    if (aRowCount != aColCount) {
      throw std::range_error("Matrix dimension are not square.");
    }
  }

  /// Constructor for an identity matrix.
  mat(int aRowCount, bool aIdentity) : mat(aRowCount) {
    if (aIdentity) {
      int minN = get_row_count() * (get_row_count() + 1);
      for (int i = 0; i < minN; i += get_row_count() + 1) {
        data.q[i] = 1;
      }
    }
  }

  mat(const self& rhs) = default;
  mat(self&& rhs) = default;
  self& operator=(const self& rhs) = default;
  self& operator=(self&& rhs) = default;
  ~mat() override = default;

  /// Explicit constructor from a any type of matrix.
  template <ReadableMatrix Matrix>
  explicit mat(const Matrix& M)
      : mat(M.get_row_count()) {
    if (M.get_col_count() != M.get_row_count()) {
      throw std::range_error("Matrix is not square!");
    }
    auto it = data.q.begin();
    for (int j = 0; j < get_row_count(); ++j) {
      for (int i = 0; i < get_row_count(); ++i, ++it) {
        *it = M(i, j);
      }
    }
  }

  /// Constructor from a vector of column major values.
  mat(container_type Q, int aRowCount) {
    data.q = std::move(Q);
    if constexpr (is_dynamic_size) {
      data.rowCount = aRowCount;
    }
  }

  /// Constructs a 2x2 matrix from four elements.
  mat(const_reference a11, const_reference a12, const_reference a21,
      const_reference a22)
      : mat(2) {
    data.q[0] = a11;
    data.q[1] = a21;
    data.q[2] = a12;
    data.q[3] = a22;
  }

  /// Constructs a 3x3 matrix from nine elements.
  mat(const_reference a11, const_reference a12, const_reference a13,
      const_reference a21, const_reference a22, const_reference a23,
      const_reference a31, const_reference a32, const_reference a33)
      : mat(3) {
    data.q[0] = a11;
    data.q[1] = a21;
    data.q[2] = a31;
    data.q[3] = a12;
    data.q[4] = a22;
    data.q[5] = a32;
    data.q[6] = a13;
    data.q[7] = a23;
    data.q[8] = a33;
  }

  /// Constructs a 4x4 matrix from sixteen elements.
  mat(const_reference a11, const_reference a12, const_reference a13,
      const_reference a14, const_reference a21, const_reference a22,
      const_reference a23, const_reference a24, const_reference a31,
      const_reference a32, const_reference a33, const_reference a34,
      const_reference a41, const_reference a42, const_reference a43,
      const_reference a44)
      : mat(4) {
    data.q[0] = a11;
    data.q[1] = a21;
    data.q[2] = a31;
    data.q[3] = a41;
    data.q[4] = a12;
    data.q[5] = a22;
    data.q[6] = a32;
    data.q[7] = a42;
    data.q[8] = a13;
    data.q[9] = a23;
    data.q[10] = a33;
    data.q[11] = a43;
    data.q[12] = a14;
    data.q[13] = a24;
    data.q[14] = a34;
    data.q[15] = a44;
  }

  /// The standard swap function (works with ADL).
  friend void swap(self& m1, self& m2) noexcept {
    using std::swap;
    swap(m1.data.q, m2.data.q);
    if constexpr (is_dynamic_size) {
      swap(m1.data.rowCount, m2.data.rowCount);
    }
  }

  /// A swap function to swap the matrix with a container of values to fill the matrix.
  /// \param m1 The matrix to swap with the container.
  /// \param q2 The container that will be swapped with m1's internal container.
  /// \param rowCount2 The row-count corresponding to q2's data.
  friend void swap(self& m1, container_type& q2, int& rowCount2) noexcept {
    using std::swap;
    swap(m1.data.q, q2);
    if constexpr (is_dynamic_size) {
      swap(m1.data.rowCount, rowCount2);
    }
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Matrix indexing accessor for read-write access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  reference operator()(int i, int j) { return data.q[j * get_row_count() + i]; }
  /// Matrix indexing accessor for read-only access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  const_reference operator()(int i, int j) const {
    return data.q[j * get_row_count() + i];
  }

  /// Sub-matrix operator, accessor for read/write.
  mat_sub_block<self> operator()(const std::pair<int, int>& r,
                                 const std::pair<int, int>& c) {
    return sub(*this)(r, c);
  }

  /// Sub-matrix operator, accessor for read only.
  mat_const_sub_block<self> operator()(const std::pair<int, int>& r,
                                       const std::pair<int, int>& c) const {
    return sub(*this)(r, c);
  }

  /// Sub-matrix operator, accessor for read/write.
  mat_col_slice<self> operator()(int r, const std::pair<int, int>& c) {
    return slice(*this)(r, c);
  }

  /// Sub-matrix operator, accessor for read only.
  mat_const_col_slice<self> operator()(int r,
                                       const std::pair<int, int>& c) const {
    return slice(*this)(r, c);
  }

  /// Sub-matrix operator, accessor for read/write.
  mat_row_slice<self> operator()(const std::pair<int, int>& r, int c) {
    return slice(*this)(r, c);
  }

  /// Sub-matrix operator, accessor for read only.
  mat_const_row_slice<self> operator()(const std::pair<int, int>& r,
                                       int c) const {
    return slice(*this)(r, c);
  }

  /// Gets the row-count (number of rows) of the matrix.
  /// \return number of rows of the matrix.
  int get_row_count() const noexcept {
    if constexpr (is_dynamic_size) {
      return data.rowCount;
    } else {
      return RowCount;
    }
  }
  /// Gets the column-count (number of columns) of the matrix.
  /// \return number of columns of the matrix.
  int get_col_count() const noexcept { return get_row_count(); }

  /// Gets the row-count and column-count of the matrix, as a std::pair of values.
  /// \return the row-count and column-count of the matrix, as a std::pair of values.
  std::pair<int, int> size() const noexcept {
    return {get_row_count(), get_row_count()};
  }
  /// Sets the row-count and column-count of the matrix, via a std::pair of dimension values.
  /// \param sz new dimensions for the matrix.
  void resize(const std::pair<int, int>& sz) { set_row_count(sz.first, true); }

  /// Sets the row-count (number of rows) of the matrix.
  /// \param aRowCount new number of rows for the matrix.
  /// \param aPreserveData If true, the resizing will preserve all the data it can.
  void set_row_count(int aRowCount, bool aPreserveData = false) {
    if constexpr (!is_dynamic_size) {
      if (aRowCount != RowCount) {
        throw std::range_error("Row count mismatch!");
      }
    } else {
      if (aPreserveData) {
        if (aRowCount > get_row_count()) {
          for (int i = get_row_count(); i != 0; --i) {
            data.q.insert(data.q.begin() + i * get_row_count(),
                          aRowCount - get_row_count(), value_type(0));
          }
        } else if (aRowCount < get_row_count()) {
          for (int i = get_row_count(); i != 0; --i) {
            data.q.erase(data.q.begin() + (i - 1) * get_row_count() + aRowCount,
                         data.q.begin() + i * get_row_count());
          }
        } else {
          return;
        }
      }
      data.q.resize(aRowCount * aRowCount, value_type(0));
      data.rowCount = aRowCount;
    }
  }

  /// Sets the column-count (number of columns) of the matrix.
  /// \param aColCount new number of columns for the matrix.
  /// \param aPreserveData If true, the resizing will preserve all the data it can.
  void set_col_count(int aColCount, bool aPreserveData = false) {
    set_row_count(aColCount, aPreserveData);
  }

  row_iterator first_row() { return data.q.begin(); }
  const_row_iterator first_row() const { return data.q.begin(); }
  row_iterator last_row() { return data.q.begin() + get_row_count(); }
  const_row_iterator last_row() const {
    return data.q.begin() + get_row_count();
  }
  row_iterator first_row(col_iterator cit) {
    int diff = cit.base() - data.q.begin();
    return data.q.begin() + ((diff / get_row_count()) * get_row_count());
  }
  const_row_iterator first_row(const_col_iterator cit) const {
    int diff = cit.base() - data.q.begin();
    return data.q.begin() + ((diff / get_row_count()) * get_row_count());
  }
  row_iterator last_row(col_iterator cit) {
    int diff = cit.base() - data.q.begin();
    return data.q.begin() + ((diff / get_row_count() + 1) * get_row_count());
  }
  const_row_iterator last_row(const_col_iterator cit) const {
    int diff = cit.base() - data.q.begin();
    return data.q.begin() + ((diff / get_row_count() + 1) * get_row_count());
  }
  std::pair<row_iterator, row_iterator> rows() {
    return {first_row(), last_row()};
  }
  std::pair<const_row_iterator, const_row_iterator> rows() const {
    return {first_row(), last_row()};
  }
  std::pair<row_iterator, row_iterator> rows(col_iterator cit) {
    return {first_row(cit), last_row(cit)};
  }
  std::pair<const_row_iterator, const_row_iterator> rows(
      const_col_iterator cit) const {
    return {first_row(cit), last_row(cit)};
  }

  col_iterator first_col() { return {data.q.begin(), get_row_count()}; }
  const_col_iterator first_col() const {
    return {data.q.begin(), get_row_count()};
  }
  col_iterator last_col() {
    return {data.q.begin() + get_row_count() * get_row_count(),
            get_row_count()};
  }
  const_col_iterator last_col() const {
    return {data.q.begin() + get_row_count() * get_row_count(),
            get_row_count()};
  }
  col_iterator first_col(row_iterator rit) {
    return {data.q.begin() + ((rit - data.q.begin()) % get_row_count()),
            get_row_count()};
  }
  const_col_iterator first_col(const_row_iterator rit) const {
    return {data.q.begin() + ((rit - data.q.begin()) % get_row_count()),
            get_row_count()};
  }
  col_iterator last_col(row_iterator rit) {
    return {data.q.begin() + ((rit - data.q.begin()) % get_row_count()) +
                get_row_count() * get_row_count(),
            get_row_count()};
  }
  const_col_iterator last_col(const_row_iterator rit) const {
    return {data.q.begin() + ((rit - data.q.begin()) % get_row_count()) +
                get_row_count() * get_row_count(),
            get_row_count()};
  }
  std::pair<col_iterator, col_iterator> cols() {
    return {first_col(), last_col()};
  }
  std::pair<const_col_iterator, const_col_iterator> cols() const {
    return {first_col(), last_col()};
  }
  std::pair<col_iterator, col_iterator> cols(row_iterator rit) {
    return {first_col(rit), last_col(rit)};
  }
  std::pair<const_col_iterator, const_col_iterator> cols(
      const_row_iterator rit) const {
    return {first_col(rit), last_col(rit)};
  }

  /*******************************************************************************
                           Assignment Operators
  *******************************************************************************/

  /// Standard Assignment operator with standard semantics.
  /// Strong exception-safety.
  template <typename Matrix>
  self& operator=(const Matrix& M) {
    self tmp(M);
    swap(*this, tmp);
    return *this;
  }

  /// Add-and-store operator with standard semantics.
  template <ReadableMatrix Matrix>
  self& operator+=(const Matrix& M) {
    if ((M.get_col_count() != get_row_count()) ||
        (M.get_row_count() != get_row_count())) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    auto it = data.q.begin();
    for (int j = 0; j < get_row_count(); ++j) {
      for (int i = 0; i < get_row_count(); ++i, ++it) {
        *it += M(i, j);
      }
    }
    return *this;
  }

  /// Sub-and-store operator with standard semantics.
  template <ReadableMatrix Matrix>
  self& operator-=(const Matrix& M) {
    if ((M.get_col_count() != get_row_count()) ||
        (M.get_row_count() != get_row_count())) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    auto it = data.q.begin();
    for (int j = 0; j < get_row_count(); ++j) {
      for (int i = 0; i < get_row_count(); ++i, ++it) {
        *it -= M(i, j);
      }
    }
    return *this;
  }

  /// Multiply-and-store operator with standard semantics.
  template <typename RhsType>
  self& operator*=(const RhsType& rhs) {
    if constexpr (ReadableMatrix<RhsType>) {
      self result(*this * rhs);
      swap(*this, result);
      return *this;
    } else {
      for (auto& x : data.q) {
        x *= rhs;
      }
      return *this;
    }
  }

  /// General negation operator for any type of matrices. This is a default operator
  /// that will be called if no better special-purpose overload exists.
  /// \return General column-major matrix.
  self operator-() const {
    self result(*this);
    auto itr = result.data.q.begin();
    for (auto it = data.q.begin(); it != data.q.end(); ++it, ++itr) {
      *itr = -(*it);
    }
    return result;
  }

  /// Transposes the matrix M by a simple copy with a change of alignment.
  /// \param M The matrix to be transposed.
  /// \return The transpose of M.
  friend auto transpose(const self& M) {
    return mat<T, mat_structure::square, mat_alignment::row_major, RowCount,
               RowCount>(M.data.q, M.get_row_count());
  }

  /// Transposes the matrix M by simply moving the data of M into a matrix of different alignment.
  /// \param M The matrix to be transposed.
  /// \return The transpose of M.
  friend auto transpose(self&& M) {
    using std::swap;
    mat<T, mat_structure::square, mat_alignment::row_major, RowCount, RowCount>
        result;
    if constexpr (is_dynamic_size) {
      swap(result, M.data.q, M.data.rowCount);
    } else {
      int rowCount = RowCount;
      swap(result, M.data.q, rowCount);
    }
    return result;
  }

  /// Returns the trace of matrix M.
  /// \param M A matrix.
  /// \return the trace of matrix M.
  friend value_type trace(const self& M) {
    auto sum = value_type(0);
    for (int i = 0; i < M.data.q.size(); i += M.get_row_count() + 1) {
      sum += M.data.q[i];
    }
    return sum;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& std::pair<std::string, const container_type&>("q", data.q);
    if constexpr (is_dynamic_size) {
      A& std::pair<std::string, int>("rowCount", data.rowCount);
    }
  }
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& std::pair<std::string, container_type&>("q", data.q);
    if constexpr (is_dynamic_size) {
      A& std::pair<std::string, int&>("rowCount", data.rowCount);
    }
  }

  RK_RTTI_REGISTER_CLASS_1BASE(self, 1, serializable)
};

/// This class template specialization implements a matrix with square structure
/// and row-major alignment. This class is serializable and registered to the ReaK::rtti
/// system. This matrix type is dynamically resizable.
///
/// Models: ReadableMatrix, WritableMatrix, FullyWritableMatrix,
/// and ResizableMatrix.
///
/// \tparam T Arithmetic type of the elements of the matrix.
template <typename T, unsigned int RowCount>
class mat<T, mat_structure::square, mat_alignment::row_major, RowCount,
          RowCount> : public serializable {
 public:
  using self = mat<T, mat_structure::square, mat_alignment::row_major, RowCount,
                   RowCount>;
  static constexpr bool is_dynamic_size = (RowCount == 0);

  using value_type = T;
  using container_type =
      std::conditional_t<is_dynamic_size, std::vector<value_type>,
                         std::array<value_type, RowCount * RowCount>>;

  using reference = typename container_type::reference;
  using const_reference = typename container_type::const_reference;
  using pointer = typename container_type::pointer;
  using const_pointer = typename container_type::const_pointer;

  using row_iterator = stride_iterator<typename container_type::iterator>;
  using const_row_iterator =
      stride_iterator<typename container_type::const_iterator>;
  using col_iterator = typename container_type::iterator;
  using const_col_iterator = typename container_type::const_iterator;

  using size_type = int;
  using difference_type = int;

  static constexpr unsigned int static_row_count = RowCount;
  static constexpr unsigned int static_col_count = RowCount;
  static constexpr mat_alignment::tag alignment = mat_alignment::row_major;
  static constexpr mat_structure::tag structure = mat_structure::square;

 private:
  struct dynamic_data {
    container_type q = {};
    int rowCount = 0;
  };
  struct static_data {
    container_type q = {};
  };
  /// Array which holds all the values of the matrix (dimension: rowCount x rowCount).
  std::conditional_t<is_dynamic_size, dynamic_data, static_data> data;

 public:
  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /// Default constructor: sets all to zero.
  mat() : data() {}

  /// Constructor for a sized matrix.
  explicit mat(int aRowCount, const value_type& aFill = value_type(0.0)) {
    if constexpr (is_dynamic_size) {
      data.q.resize(aRowCount * aRowCount, aFill);
      data.rowCount = aRowCount;
    } else {
      if (aRowCount != RowCount) {
        throw std::range_error("Row count mismatch!");
      }
      std::fill(data.q.begin(), data.q.end(), aFill);
    }
  }

  mat(int aRowCount, int aColCount, const value_type& aFill = value_type(0.0))
      : mat(aRowCount, aFill) {
    if (aRowCount != aColCount) {
      throw std::range_error("Matrix dimension are not square.");
    }
  }

  /// Constructor for an identity matrix.
  mat(int aRowCount, bool aIdentity) : mat(aRowCount) {
    if (aIdentity) {
      int minN = get_row_count() * (get_row_count() + 1);
      for (int i = 0; i < minN; i += get_row_count() + 1) {
        data.q[i] = 1;
      }
    }
  }

  /// Default Copy/Move Constructors.
  mat(const self& rhs) = default;
  mat(self&& rhs) = default;
  self& operator=(const self& rhs) = default;
  self& operator=(self&& rhs) = default;
  ~mat() override = default;

  /// Explicit constructor from a any type of matrix.
  template <ReadableMatrix Matrix>
  explicit mat(const Matrix& M)
      : mat(M.get_row_count()) {
    if (M.get_col_count() != M.get_row_count()) {
      throw std::range_error("Matrix is not square!");
    }
    auto it = data.q.begin();
    for (int i = 0; i < get_row_count(); ++i) {
      for (int j = 0; j < get_row_count(); ++j, ++it) {
        *it = M(i, j);
      }
    }
  }

  /// Constructor from a vector of column major values.
  mat(container_type Q, int aRowCount) {
    data.q = std::move(Q);
    if constexpr (is_dynamic_size) {
      data.rowCount = aRowCount;
    }
  }

  /// Constructs a 2x2 matrix from four elements.
  mat(const value_type& a11, const value_type& a12, const value_type& a21,
      const value_type& a22)
      : mat(2) {
    data.q[0] = a11;
    data.q[1] = a12;
    data.q[2] = a21;
    data.q[3] = a22;
  }

  /// Constructs a 3x3 matrix from nine elements.
  mat(const value_type& a11, const value_type& a12, const value_type& a13,
      const value_type& a21, const value_type& a22, const value_type& a23,
      const value_type& a31, const value_type& a32, const value_type& a33)
      : mat(3) {
    data.q[0] = a11;
    data.q[1] = a12;
    data.q[2] = a13;
    data.q[3] = a21;
    data.q[4] = a22;
    data.q[5] = a23;
    data.q[6] = a31;
    data.q[7] = a32;
    data.q[8] = a33;
  }

  /// Constructs a 4x4 matrix from sixteen elements.
  mat(const value_type& a11, const value_type& a12, const value_type& a13,
      const value_type& a14, const value_type& a21, const value_type& a22,
      const value_type& a23, const value_type& a24, const value_type& a31,
      const value_type& a32, const value_type& a33, const value_type& a34,
      const value_type& a41, const value_type& a42, const value_type& a43,
      const value_type& a44)
      : mat(4) {
    data.q[0] = a11;
    data.q[1] = a12;
    data.q[2] = a13;
    data.q[3] = a14;
    data.q[4] = a21;
    data.q[5] = a22;
    data.q[6] = a23;
    data.q[7] = a24;
    data.q[8] = a31;
    data.q[9] = a32;
    data.q[10] = a33;
    data.q[11] = a34;
    data.q[12] = a41;
    data.q[13] = a42;
    data.q[14] = a43;
    data.q[15] = a44;
  }

  /// The standard swap function (works with ADL).
  friend void swap(self& m1, self& m2) noexcept {
    using std::swap;
    swap(m1.data.q, m2.data.q);
    if constexpr (is_dynamic_size) {
      swap(m1.data.rowCount, m2.data.rowCount);
    }
  }

  /// A swap function to swap the matrix with a container of values to fill the matrix.
  /// \param m1 The matrix to swap with the container.
  /// \param q2 The container that will be swapped with m1's internal container.
  /// \param rowCount2 The row-count corresponding to q2's data.
  friend void swap(self& m1, container_type& q2, int& rowCount2) noexcept {
    using std::swap;
    swap(m1.data.q, q2);
    if constexpr (is_dynamic_size) {
      swap(m1.data.rowCount, rowCount2);
    }
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Matrix indexing accessor for read-write access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  reference operator()(int i, int j) { return data.q[i * get_row_count() + j]; }
  /// Matrix indexing accessor for read-only access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  const_reference operator()(int i, int j) const {
    return data.q[i * get_row_count() + j];
  }

  /// Sub-matrix operator, accessor for read/write.
  mat_sub_block<self> operator()(const std::pair<int, int>& r,
                                 const std::pair<int, int>& c) {
    return sub(*this)(r, c);
  }

  /// Sub-matrix operator, accessor for read only.
  mat_const_sub_block<self> operator()(const std::pair<int, int>& r,
                                       const std::pair<int, int>& c) const {
    return sub(*this)(r, c);
  }

  /// Sub-matrix operator, accessor for read/write.
  mat_col_slice<self> operator()(int r, const std::pair<int, int>& c) {
    return slice(*this)(r, c);
  }

  /// Sub-matrix operator, accessor for read only.
  mat_const_col_slice<self> operator()(int r,
                                       const std::pair<int, int>& c) const {
    return slice(*this)(r, c);
  }

  /// Sub-matrix operator, accessor for read/write.
  mat_row_slice<self> operator()(const std::pair<int, int>& r, int c) {
    return slice(*this)(r, c);
  }

  /// Sub-matrix operator, accessor for read only.
  mat_const_row_slice<self> operator()(const std::pair<int, int>& r,
                                       int c) const {
    return slice(*this)(r, c);
  }

  /// Gets the row-count (number of rows) of the matrix.
  /// \return number of rows of the matrix.
  int get_row_count() const noexcept {
    if constexpr (is_dynamic_size) {
      return data.rowCount;
    } else {
      return RowCount;
    }
  }
  /// Gets the column-count (number of columns) of the matrix.
  /// \return number of columns of the matrix.
  int get_col_count() const noexcept { return get_row_count(); }

  /// Gets the row-count and column-count of the matrix, as a std::pair of values.
  /// \return the row-count and column-count of the matrix, as a std::pair of values.
  std::pair<int, int> size() const noexcept {
    return {get_row_count(), get_row_count()};
  }
  /// Sets the row-count and column-count of the matrix, via a std::pair of dimension values.
  /// \param sz new dimensions for the matrix.
  void resize(const std::pair<int, int>& sz) { set_col_count(sz.first, true); }

  /// Sets the column-count (number of columns) of the matrix.
  /// \param aColCount new number of columns for the matrix.
  /// \param aPreserveData If true, the resizing will preserve all the data it can.
  void set_col_count(int aColCount, bool aPreserveData = false) {
    if constexpr (!is_dynamic_size) {
      if (aColCount != RowCount) {
        throw std::range_error("Row/col count mismatch!");
      }
    } else {
      if (aPreserveData) {
        if (aColCount > data.rowCount) {
          for (int i = data.rowCount; i != 0; --i) {
            data.q.insert(data.q.begin() + i * data.rowCount,
                          aColCount - data.rowCount, value_type(0.0));
          }
        } else if (aColCount < data.rowCount) {
          for (int i = data.rowCount; i != 0; --i) {
            data.q.erase(data.q.begin() + (i - 1) * data.rowCount + aColCount,
                         data.q.begin() + i * data.rowCount);
          }
        } else {
          return;
        }
      }
      data.q.resize(aColCount * aColCount, value_type(0));
      data.rowCount = aColCount;
    }
  }

  /// Sets the row-count (number of rows) of the matrix.
  /// \param aRowCount new number of rows for the matrix.
  /// \param aPreserveData If true, the resizing will preserve all the data it can.
  void set_row_count(int aRowCount, bool aPreserveData = false) {
    set_col_count(aRowCount, aPreserveData);
  }

  row_iterator first_row() { return {data.q.begin(), get_row_count()}; }
  const_row_iterator first_row() const {
    return {data.q.begin(), get_row_count()};
  }
  row_iterator last_row() {
    return {data.q.begin() + get_row_count() * get_row_count(),
            get_row_count()};
  }
  const_row_iterator last_row() const {
    return {data.q.begin() + get_row_count() * get_row_count(),
            get_row_count()};
  }
  row_iterator first_row(col_iterator cit) {
    return {data.q.begin() + ((cit - data.q.begin()) % get_row_count()),
            get_row_count()};
  }
  const_row_iterator first_row(const_col_iterator cit) const {
    return {data.q.begin() + ((cit - data.q.begin()) % get_row_count()),
            get_row_count()};
  }
  row_iterator last_row(col_iterator cit) {
    return {data.q.begin() + ((cit - data.q.begin()) % get_row_count()) +
                get_row_count() * get_row_count(),
            get_row_count()};
  }
  const_row_iterator last_row(const_col_iterator cit) const {
    return {data.q.begin() + ((cit - data.q.begin()) % get_row_count()) +
                get_row_count() * get_row_count(),
            get_row_count()};
  }
  std::pair<row_iterator, row_iterator> rows() {
    return {first_row(), last_row()};
  }
  std::pair<const_row_iterator, const_row_iterator> rows() const {
    return {first_row(), last_row()};
  }
  std::pair<row_iterator, row_iterator> rows(col_iterator cit) {
    return {first_row(cit), last_row(cit)};
  }
  std::pair<const_row_iterator, const_row_iterator> rows(
      const_col_iterator cit) const {
    return {first_row(cit), last_row(cit)};
  }

  col_iterator first_col() { return {data.q.begin()}; }
  const_col_iterator first_col() const { return {data.q.begin()}; }
  col_iterator last_col() {
    return {data.q.begin() + get_row_count() * get_row_count()};
  }
  const_col_iterator last_col() const {
    return {data.q.begin() + get_row_count() * get_row_count()};
  }
  col_iterator first_col(row_iterator rit) {
    int diff = rit.base() - data.q.begin();
    return data.q.begin() + ((diff / get_row_count()) * get_row_count());
  }
  const_col_iterator first_col(const_row_iterator rit) const {
    int diff = rit.base() - data.q.begin();
    return data.q.begin() + ((diff / get_row_count()) * get_row_count());
  }
  col_iterator last_col(row_iterator rit) {
    int diff = rit.base() - data.q.begin();
    return data.q.begin() + ((diff / get_row_count() + 1) * get_row_count());
  }
  const_col_iterator last_col(const_row_iterator rit) const {
    int diff = rit.base() - data.q.begin();
    return data.q.begin() + ((diff / get_row_count() + 1) * get_row_count());
  }
  std::pair<col_iterator, col_iterator> cols() {
    return {first_col(), last_col()};
  }
  std::pair<const_col_iterator, const_col_iterator> cols() const {
    return {first_col(), last_col()};
  }
  std::pair<col_iterator, col_iterator> cols(row_iterator rit) {
    return {first_col(rit), last_col(rit)};
  }
  std::pair<const_col_iterator, const_col_iterator> cols(
      const_row_iterator rit) const {
    return {first_col(rit), last_col(rit)};
  }

  /*******************************************************************************
                           Assignment Operators
  *******************************************************************************/

  /// Standard Assignment operator with standard semantics.
  /// Strong exception safety.
  template <typename Matrix>
  self& operator=(const Matrix& M) {
    self tmp(M);
    swap(*this, tmp);
    return *this;
  }

  /// Add-and-store operator with standard semantics.
  template <ReadableMatrix Matrix>
  self& operator+=(const Matrix& M) {
    if ((M.get_col_count() != get_row_count()) ||
        (M.get_row_count() != get_row_count())) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    auto it = data.q.begin();
    for (int i = 0; i < get_row_count(); ++i) {
      for (int j = 0; j < get_row_count(); ++j, ++it) {
        *it += M(i, j);
      }
    }
    return *this;
  }

  /// Sub-and-store operator with standard semantics.
  template <ReadableMatrix Matrix>
  self& operator-=(const Matrix& M) {
    if ((M.get_col_count() != get_row_count()) ||
        (M.get_row_count() != get_row_count())) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    auto it = data.q.begin();
    for (int i = 0; i < get_row_count(); ++i) {
      for (int j = 0; j < get_row_count(); ++j, ++it) {
        *it -= M(i, j);
      }
    }
    return *this;
  }

  /// Scalar-multiply-and-store operator with standard semantics.
  template <typename RhsType>
  self& operator*=(const RhsType& rhs) {
    if constexpr (ReadableMatrix<RhsType>) {
      self result(*this * rhs);
      swap(*this, result);
      return *this;
    } else {
      for (auto& x : data.q) {
        x *= rhs;
      }
      return *this;
    }
  }

  /// General negation operator for any type of matrices. This is a default operator
  /// that will be called if no better special-purpose overload exists.
  /// \return General column-major matrix.
  self operator-() const {
    self result(*this);
    auto itr = result.data.q.begin();
    for (auto it = data.q.begin(); it != data.q.end(); ++it, ++itr) {
      *itr = -(*it);
    }
    return result;
  }

  /// Transposes the matrix M by a simple copy with a change of alignment.
  /// \param M The matrix to be transposed.
  /// \return The transpose of M.
  friend auto transpose(const self& M) {
    return mat<T, mat_structure::square, mat_alignment::column_major, RowCount,
               RowCount>(M.data.q, M.get_row_count());
  }

  /// Transposes the matrix M by simply moving the data of M into a matrix of different alignment.
  /// \param M The matrix to be transposed.
  /// \return The transpose of M.
  friend auto transpose(self&& M) {
    using std::swap;
    mat<T, mat_structure::square, mat_alignment::column_major, RowCount,
        RowCount>
        result;
    if constexpr (is_dynamic_size) {
      swap(result, M.data.q, M.data.rowCount);
    } else {
      int rowCount = RowCount;
      swap(result, M.data.q, rowCount);
    }
    return result;
  }

  /// Returns the trace of matrix M.
  /// \param M A diagonal matrix.
  /// \return the trace of matrix M.
  friend value_type trace(const self& M) {
    value_type sum = value_type(0);
    for (int i = 0; i < M.data.q.size(); i += M.get_row_count() + 1) {
      sum += M.data.q[i];
    }
    return sum;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& std::pair<std::string, const container_type&>("q", data.q);
    if constexpr (is_dynamic_size) {
      A& std::pair<std::string, int>("rowCount", data.rowCount);
    }
  }
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& std::pair<std::string, container_type&>("q", data.q);
    if constexpr (is_dynamic_size) {
      A& std::pair<std::string, int&>("rowCount", data.rowCount);
    }
  }

  RK_RTTI_REGISTER_CLASS_1BASE(self, 1, serializable)
};

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_ALG_SQUARE_H_

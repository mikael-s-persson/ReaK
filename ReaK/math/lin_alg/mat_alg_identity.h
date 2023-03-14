/**
 * \file mat_alg_identity.h
 *
 * This library declares matrix specializations for representing and manipulating identity matrices.
 * This library implements many overloaded operators that turn out to be more efficiently implemented
 * if specialized for the identity matrix case. All those overloads are automatically selected through
 * Sfinae switches, and the identity matrix class is simply a partial specialization of the "ReaK::mat"
 * class template, so, the burden on the user is minimal.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date april 2011
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

#ifndef REAK_MATH_LIN_ALG_MAT_ALG_IDENTITY_H_
#define REAK_MATH_LIN_ALG_MAT_ALG_IDENTITY_H_

#include "ReaK/math/lin_alg/mat_alg_general.h"
#include "ReaK/math/lin_alg/mat_alg_nil.h"

namespace ReaK {

/// This class implements a place-holder or interface-implementation to represent
/// a identity matrix (all entries zero except diagonal is all unity). This is useful to build for example a
/// block-matrix with some identity-matrix blocks, and, of course, requires minimal storage space.
///
/// Models: ReadableMatrixConcept and ResizableMatrixConcept.
///
/// \tparam T Arithmetic type of the elements of the matrix.
/// \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or
/// mat_alignment::column_major (default).
template <typename T, mat_alignment::tag Alignment, unsigned int RowCount>
class mat<T, mat_structure::identity, Alignment, RowCount, RowCount>
    : public serializable {
 public:
  using self = mat<T, mat_structure::identity, Alignment, RowCount, RowCount>;
  static constexpr bool is_dynamic_size = (RowCount == 0);

  using value_type = T;
  using container_type = void;

  using reference = void;
  using const_reference = T;
  using pointer = void;
  using const_pointer = void;

  using row_iterator = void;
  using const_row_iterator = void;
  using col_iterator = void;
  using const_col_iterator = void;

  using size_type = int;
  using difference_type = int;

  static constexpr unsigned int static_row_count = RowCount;
  static constexpr unsigned int static_col_count = RowCount;
  static constexpr mat_alignment::tag alignment = Alignment;
  static constexpr mat_structure::tag structure = mat_structure::identity;

 private:
  struct dynamic_data {
    int rowCount = 0;
  };
  struct static_data {};
  /// Hold run-time size rowCount.
  std::conditional_t<is_dynamic_size, dynamic_data, static_data> data;
 public:
  /// Default constructor. Sets dimensions to zero.
  mat() = default;
  /// Constructs an identity matrix to the given dimensions.
  explicit mat(int aRowCount) {
    if constexpr (is_dynamic_size) {
      data.rowCount = aRowCount;
    } else {
      if (aRowCount != RowCount) {
        throw std::range_error("Row count mismatch!");
      }
    }
  }

  mat(const self& rhs) = default;
  mat(self&& rhs) = default;
  self& operator=(const self& rhs) = default;
  self& operator=(self&& rhs) = default;
  ~mat() override = default;

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Matrix indexing accessor for read-only access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  const_reference operator()(int i, int j) const {
    if (i == j) {
      return value_type(1);
    }
    return value_type(0);
  }

  /// Sub-matrix operator, accessor for read only.
  mat_const_sub_block<self> operator()(const std::pair<int, int>& r,
                                       const std::pair<int, int>& c) const {
    return sub(*this)(r, c);
  }

  /// Sub-matrix operator, accessor for read only.
  mat_const_col_slice<self> operator()(int r,
                                       const std::pair<int, int>& c) const {
    return slice(*this)(r, c);
  }

  /// Sub-matrix operator, accessor for read only.
  mat_const_row_slice<self> operator()(const std::pair<int, int>& r,
                                       int c) const {
    return slice(*this)(r, c);
  }

  /// Gets the row-count (number of rows) of the matrix.
  /// \return number of rows of the matrix.
  int get_row_count() const {
    if constexpr (is_dynamic_size) {
      return data.rowCount;
    } else {
      return RowCount;
    }
  }

  /// Sets the row-count (number of rows) of the matrix.
  /// \param aRowCount new number of rows for the matrix.
  /// \param aPreserveData If true, the resizing will preserve all the data it can.
  void set_row_count(int aRowCount, bool /*unused*/ = false) {
    if constexpr (is_dynamic_size) {
      data.rowCount = aRowCount;
    } else {
      if (aRowCount != RowCount) {
        throw std::range_error("Row count mismatch!");
      }
    }
  }

  /// Gets the column-count (number of columns) of the matrix.
  /// \return number of columns of the matrix.
  int get_col_count() const { return get_row_count(); }

  /// Sets the column-count (number of columns) of the matrix.
  /// \param aColCount new number of columns for the matrix.
  /// \param aPreserveData If true, the resizing will preserve all the data it can.
  void set_col_count(int aColCount, bool /*unused*/ = false) {
    set_row_count(aColCount);
  }

  /// Negate the matrix.
  /// \return This matrix, by constant reference.
  mat<value_type, mat_structure::scalar> operator-() const {
    return mat<value_type, mat_structure::scalar, Alignment, RowCount, RowCount>(get_row_count(), value_type(-1));
  }

  /// Transpose the matrix.
  /// \param rhs the matrix to be transposed.
  /// \return The rhs matrix, by value.
  friend self transpose(const self& rhs) { return rhs; }

  /// Transpose and move the matrix.
  /// \param rhs the matrix to be transposed and moved (emptied).
  /// \return The rhs matrix, by value.
  friend self transpose_move(self& rhs) {
    self result;
    swap(result, rhs);
    return result;
  }

  /// Returns the trace of a matrix.
  /// \param M A matrix.
  /// \return the trace of matrix M.
  friend value_type trace(const self& M) {
    return static_cast<value_type>(M.get_row_count());
  }

  /// Appends a matrix to another.
  /// \param lhs the matrix to which 'rhs' will be appended to.
  /// \param rhs the matrix to append to 'lhs'.
  friend void append_block_diag(self& lhs, const self& rhs) {
    lhs.set_row_count(lhs.get_row_count() + rhs.get_row_count());
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    if constexpr (is_dynamic_size) {
      A& RK_SERIAL_SAVE_WITH_ALIAS("rowCount", data.rowCount);
    }
  }
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    if constexpr (is_dynamic_size) {
      A& RK_SERIAL_LOAD_WITH_ALIAS("rowCount", data.rowCount);
    }
  }

  RK_RTTI_REGISTER_CLASS_1BASE(self, 1, serializable)
};

template <typename T>
struct mat_identity {
  using type = mat<T, mat_structure::identity>;
};

template <typename T>
using mat_identity_t = typename mat_identity<T>::type;

template <typename T>
mat<T, mat_structure::identity> mat_ident(int aRowCount) {
  return mat<T, mat_structure::identity>(aRowCount);
}

/// Scalar multiplication, always results in a scalar matrix.
/// \param M some matrix.
/// \param S some scalar.
/// \return Scalar matrix.
template <typename T, mat_alignment::tag Alignment, unsigned int RowCount>
std::enable_if_t<!is_readable_vector_v<T> && !is_readable_matrix_v<T>,
                 mat<T, mat_structure::scalar, Alignment, RowCount, RowCount>>
operator*(
    const mat<T, mat_structure::identity, Alignment, RowCount, RowCount>& M,
    const T& S) {
  return mat<T, mat_structure::scalar, Alignment, RowCount, RowCount>(
      M.get_row_count(), S);
}

/// Scalar multiplication, always results in a scalar matrix.
/// \param S some scalar.
/// \param M a null-matrix.
/// \return Scalar matrix.
template <typename T, mat_alignment::tag Alignment, unsigned int RowCount>
std::enable_if_t<!is_readable_vector_v<T> && !is_readable_matrix_v<T>,
                 mat<T, mat_structure::scalar, Alignment, RowCount, RowCount>>
operator*(
    const T& S,
    const mat<T, mat_structure::identity, Alignment, RowCount, RowCount>& M) {
  return mat<T, mat_structure::scalar, Alignment, RowCount, RowCount>(
      M.get_row_count(), S);
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_ALG_IDENTITY_H_

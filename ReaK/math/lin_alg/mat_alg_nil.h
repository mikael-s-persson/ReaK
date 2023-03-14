/**
 * \file mat_alg_nil.h
 *
 * This library declares matrix specializations for representing and manipulating nil matrices.
 * This library implements many overloaded operators that turn out to be more efficiently implemented
 * if specialized for the nil matrix case. All those overloads are automatically selected through
 * Sfinae switches, and the nil matrix class is simply a partial specialization of the "ReaK::mat"
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

#ifndef REAK_MATH_LIN_ALG_MAT_ALG_NIL_H_
#define REAK_MATH_LIN_ALG_MAT_ALG_NIL_H_

#include "ReaK/math/lin_alg/mat_alg_general.h"

#include <type_traits>

namespace ReaK {

/// This class implements a place-holder or interface-implementation to represent
/// a nil matrix (all entries zero). This is useful to build for example a
/// block-matrix with some zero-matrix blocks, and, of course, the storage is minimal.
///
/// Models: ReadableMatrixConcept and ResizableMatrixConcept.
///
/// \tparam T Arithmetic type of the elements of the matrix.
/// \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or
/// mat_alignment::column_major (default).
template <typename T, mat_alignment::tag Alignment, unsigned int RowCount,
          unsigned int ColCount>
class mat<T, mat_structure::nil, Alignment, RowCount, ColCount>
    : public serializable {
 public:
  using self = mat<T, mat_structure::nil, Alignment, RowCount, ColCount>;

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
  static constexpr unsigned int static_col_count = ColCount;
  static constexpr mat_alignment::tag alignment = Alignment;
  static constexpr mat_structure::tag structure = mat_structure::nil;

 private:
  struct row_dynamic {
    int rowCount = 0;
  };
  struct row_static {};
  /// Hold run-time size rowCount.
  std::conditional_t<(RowCount == 0), row_dynamic, row_static> row_data;
  struct col_dynamic {
    int colCount = 0;
  };
  struct col_static {};
  /// Hold run-time size colCount.
  std::conditional_t<(ColCount == 0), col_dynamic, col_static> col_data;
 public:
  /// Default constructor. Sets dimensions to zero.
  mat() = default;
  /// Constructs a null matrix to the given dimensions.
  mat(int aRowCount, int aColCount) {
    if constexpr (RowCount == 0) {
      row_data.rowCount = aRowCount;
    } else {
      if (aRowCount != RowCount) {
        throw std::range_error("Row count mismatch!");
      }
    }
    if constexpr (ColCount == 0) {
      col_data.colCount = aColCount;
    } else {
      if (aColCount != ColCount) {
        throw std::range_error("Col count mismatch!");
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
  const_reference operator()(int i, int j) const { return T(0.0); }

  /// Sub-matrix operator, accessor for read only.
  mat_const_sub_block<self> operator()(
      const std::pair<int, int>& r,
      const std::pair<int, int>& c) const {
    return sub(*this)(r, c);
  }

  /// Sub-matrix operator, accessor for read only.
  mat_const_col_slice<self> operator()(
      int r, const std::pair<int, int>& c) const {
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
    if constexpr (RowCount == 0) {
      return row_data.rowCount;
    } else {
      return RowCount;
    }
  }

  /// Sets the row-count (number of rows) of the matrix.
  /// \param aRowCount new number of rows for the matrix.
  /// \param aPreserveData If true, the resizing will preserve all the data it can.
  void set_row_count(int aRowCount, bool /*unused*/ = false) {
    if constexpr (RowCount == 0) {
      row_data.rowCount = aRowCount;
    } else {
      if (aRowCount != RowCount) {
        throw std::range_error("Row count mismatch!");
      }
    }
  }

  /// Gets the column-count (number of columns) of the matrix.
  /// \return number of columns of the matrix.
  int get_col_count() const {
    if constexpr (ColCount == 0) {
      return col_data.colCount;
    } else {
      return ColCount;
    }
  }

  /// Sets the column-count (number of columns) of the matrix.
  /// \param aColCount new number of columns for the matrix.
  /// \param aPreserveData If true, the resizing will preserve all the data it can.
  void set_col_count(int aColCount, bool /*unused*/ = false) {
    if constexpr (ColCount == 0) {
      col_data.colCount = aColCount;
    } else {
      if (aColCount != ColCount) {
        throw std::range_error("Col count mismatch!");
      }
    }
  }

  /// Negate the matrix, has no effect of course.
  /// \return This matrix, by constant reference.
  const self& operator-() const { return *this; }

  /// Transposes the matrix M.
  /// \param rhs The nil matrix to be transposed.
  /// \return The transpose of rhs.
  friend auto transpose(const self& rhs) {
    return mat<value_type, mat_structure::nil, Alignment, ColCount, RowCount>(rhs.get_col_count(), rhs.get_row_count());
  }

  /// Transposes the matrix M.
  /// \param rhs The nil matrix to be transposed.
  /// \return The transpose of rhs.
  friend auto transpose_move(self& rhs) {
    return transpose(rhs);
  }

  /// Returns the trace of the matrix.
  /// \return the trace of the matrix.
  friend value_type trace(const self& /*unused*/) { return value_type(0); }

  /// Appends the matrix 'rhs' to the end of the matrix 'lhs', which are both nil matrices.
  /// \param lhs The nil matrix to which to append the other.
  /// \param rhs The nil matrix to be appended to 'lhs'.
  template <unsigned int SubRowCount, unsigned int SubColCount>
  friend void append_block_diag(self& lhs, const mat<value_type, mat_structure::nil, Alignment, SubRowCount, SubColCount>& rhs) {
    lhs.set_col_count(lhs.get_col_count() + rhs.get_col_count());
    lhs.set_row_count(lhs.get_row_count() + rhs.get_row_count());
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    if constexpr (RowCount == 0) {
      A & RK_SERIAL_SAVE_WITH_ALIAS("rowCount", row_data.rowCount);
    }
    if constexpr (ColCount == 0) {
      A & RK_SERIAL_SAVE_WITH_ALIAS("colCount", col_data.colCount);
    }
  }
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    if constexpr (RowCount == 0) {
      A & RK_SERIAL_LOAD_WITH_ALIAS("rowCount", row_data.rowCount);
    }
    if constexpr (ColCount == 0) {
      A & RK_SERIAL_LOAD_WITH_ALIAS("colCount", col_data.colCount);
    }
  }

  RK_RTTI_REGISTER_CLASS_1BASE(self, 1, serializable)
};

template <typename T>
struct mat_null {
  using type = mat<T, mat_structure::nil>;
};

template <typename T>
using mat_null_t = typename mat_null<T>::type;

template <typename T>
mat<T, mat_structure::nil> mat_nil(int aRowCount, int aColCount) {
  return mat<T, mat_structure::nil>(aRowCount, aColCount);
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_ALG_NIL_H_

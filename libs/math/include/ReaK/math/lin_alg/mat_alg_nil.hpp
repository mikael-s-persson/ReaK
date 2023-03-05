/**
 * \file mat_alg_nil.hpp
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

#ifndef REAK_MAT_ALG_NIL_HPP
#define REAK_MAT_ALG_NIL_HPP

#include "ReaK/math/lin_alg/mat_alg_general.hpp"

#include <type_traits>

namespace ReaK {

/**
 * This class implements a place-holder or interface-implementation to represent
 * a nil matrix (all entries zero). This is useful to build for example a
 * block-matrix with some zero-matrix blocks, and, of course, the storage is minimal.
 *
 * Models: ReadableMatrixConcept and ResizableMatrixConcept.
 *
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or
 *mat_alignment::column_major (default).
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
class mat<T, mat_structure::nil, Alignment, Allocator> : public serializable {
 public:
  using self = mat<T, mat_structure::nil, Alignment, Allocator>;
  using allocator_type = void;

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

  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  static constexpr std::size_t static_row_count = 0;
  static constexpr std::size_t static_col_count = 0;
  static constexpr mat_alignment::tag alignment = Alignment;
  static constexpr mat_structure::tag structure = mat_structure::nil;

 private:
  size_type rowCount;  ///< Row Count.
  size_type colCount;  ///< Column Count.
 public:
  /**
   * Default constructor. Sets dimensions to zero.
   */
  mat() : rowCount(0), colCount(0) {}
  /**
   * Constructs a null matrix to the given dimensions.
   */
  mat(size_type aRowCount, size_type aColCount)
      : rowCount(aRowCount), colCount(aColCount) {}

  mat(const self& rhs) : rowCount(rhs.rowCount), colCount(rhs.colCount) {}
  /**
   * Default destructor.
   */
  ~mat() override = default;

  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
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
  const_reference operator()(size_type i, size_type j) const { return T(0.0); }

  /**
   * Sub-matrix operator, accessor for read only.
   * \test PASSED
   */
  mat_const_sub_block<self> operator()(
      const std::pair<size_type, size_type>& r,
      const std::pair<size_type, size_type>& c) const {
    return sub(*this)(r, c);
  }

  /**
   * Sub-matrix operator, accessor for read only.
   * \test PASSED
   */
  mat_const_col_slice<self> operator()(
      size_type r, const std::pair<size_type, size_type>& c) const {
    return slice(*this)(r, c);
  }

  /**
   * Sub-matrix operator, accessor for read only.
   * \test PASSED
   */
  mat_const_row_slice<self> operator()(const std::pair<size_type, size_type>& r,
                                       size_type c) const {
    return slice(*this)(r, c);
  }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const { return rowCount; }

  /**
   * Sets the row-count (number of rows) of the matrix.
   * \param aRowCount new number of rows for the matrix.
   * \param aPreserveData If true, the resizing will preserve all the data it can.
   * \test PASSED
   */
  void set_row_count(size_type aRowCount, bool aPreserveData = false) {
    RK_UNUSED(aPreserveData);
    rowCount = aRowCount;
  }

  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const { return colCount; }

  /**
   * Sets the column-count (number of columns) of the matrix.
   * \param aColCount new number of columns for the matrix.
   * \param aPreserveData If true, the resizing will preserve all the data it can.
   * \test PASSED
   */
  void set_col_count(size_type aColCount, bool aPreserveData = false) {
    RK_UNUSED(aPreserveData);
    colCount = aColCount;
  }

  /**
   * Negate the matrix, has no effect of course.
   * \return This matrix, by constant reference.
   * \test PASSED
   */
  const self& operator-() const { return *this; }

  /**
   * Transposes the matrix M.
   * \param rhs The nil matrix to be transposed.
   * \return The transpose of rhs.
   */
  friend self transpose(self rhs) {
    using std::swap;
    swap(rhs.colCount, rhs.rowCount);
    return rhs;
  }

  /**
   * Transposes the matrix M.
   * \param rhs The nil matrix to be transposed.
   * \return The transpose of rhs.
   */
  friend self transpose_move(self& rhs) {
    self result(rhs.colCount, rhs.rowCount);
    return result;
  }

  /**
   * Returns the trace of the matrix.
   * \return the trace of the matrix.
   */
  friend value_type trace(const self& /*unused*/) { return value_type(0); }

  /**
   * Appends the matrix 'rhs' to the end of the matrix 'lhs', which are both nil matrices.
   * \param lhs The nil matrix to which to append the other.
   * \param rhs The nil matrix to be appended to 'lhs'.
   */
  friend void append_block_diag(self& lhs, const self& rhs) {
    lhs.colCount += rhs.colCount;
    lhs.rowCount += rhs.rowCount;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(rowCount) & RK_SERIAL_SAVE_WITH_NAME(colCount);
  }
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(rowCount) & RK_SERIAL_LOAD_WITH_NAME(colCount);
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

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_readable_matrix<mat<T, mat_structure::nil, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type =
      is_readable_matrix<mat<T, mat_structure::nil, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_writable_matrix<mat<T, mat_structure::nil, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = false;
  using type =
      is_writable_matrix<mat<T, mat_structure::nil, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_fully_writable_matrix<
    mat<T, mat_structure::nil, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = false;
  using type =
      is_writable_matrix<mat<T, mat_structure::nil, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_resizable_matrix<mat<T, mat_structure::nil, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type =
      is_resizable_matrix<mat<T, mat_structure::nil, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct has_allocator_matrix<mat<T, mat_structure::nil, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = false;
  using type =
      has_allocator_matrix<mat<T, mat_structure::nil, Alignment, Allocator>>;
};

}  // namespace ReaK

#endif

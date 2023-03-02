/**
 * \file mat_alg_identity.hpp
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

#ifndef REAK_MAT_ALG_IDENTITY_HPP
#define REAK_MAT_ALG_IDENTITY_HPP

#include "mat_alg_general.hpp"
#include "mat_alg_nil.hpp"

namespace ReaK {

/**
 * This class implements a place-holder or interface-implementation to represent
 * a identity matrix (all entries zero except diagonal is all unity). This is useful to build for example a
 * block-matrix with some identity-matrix blocks, and, of course, requires minimal storage space.
 *
 * Models: ReadableMatrixConcept and ResizableMatrixConcept.
 *
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or
 *mat_alignment::column_major (default).
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
class mat<T, mat_structure::identity, Alignment, Allocator>
    : public serializable {
 public:
  using self = mat<T, mat_structure::identity, Alignment, Allocator>;
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
  static constexpr mat_structure::tag structure = mat_structure::identity;

 private:
  size_type rowCount;  ///< Row Count.
 public:
  /**
   * Default constructor. Sets dimensions to zero.
   */
  mat() : rowCount(0) {}
  /**
   * Constructs an identity matrix to the given dimensions.
   */
  explicit mat(size_type aRowCount) : rowCount(aRowCount) {}

  mat(const self& rhs) : rowCount(rhs.rowCount) {}
  /**
   * Default destructor.
   */
  ~mat() override = default;

  /**
   * Standard swap function (works with ADL).
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
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
  const_reference operator()(int i, int j) const {
    if (i == j) {
      return value_type(1);
    }
    return value_type(0);
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
   * Sub-matrix operator, accessor for read only.
   * \test PASSED
   */
  mat_const_col_slice<self> operator()(int r,
                                       const std::pair<int, int>& c) const {
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
  size_type get_col_count() const { return rowCount; }

  /**
   * Sets the column-count (number of columns) of the matrix.
   * \param aColCount new number of columns for the matrix.
   * \param aPreserveData If true, the resizing will preserve all the data it can.
   * \test PASSED
   */
  void set_col_count(size_type aColCount, bool aPreserveData = false) {
    RK_UNUSED(aPreserveData);
    rowCount = aColCount;
  }

  /**
   * Negate the matrix.
   * \return This matrix, by constant reference.
   * \test PASSED
   */
  mat<value_type, mat_structure::scalar> operator-() const {
    return mat<value_type, mat_structure::scalar>(rowCount, value_type(-1));
  }

  /**
   * Transpose the matrix.
   * \param rhs the matrix to be transposed.
   * \return The rhs matrix, by value.
   * \test PASSED
   */
  friend self transpose(const self& rhs) { return rhs; }

  /**
   * Transpose and move the matrix.
   * \param rhs the matrix to be transposed and moved (emptied).
   * \return The rhs matrix, by value.
   * \test PASSED
   */
  friend self transpose_move(self& rhs) {
    self result;
    swap(result, rhs);
    return result;
  }

  /**
   * Returns the trace of a matrix.
   * \param M A matrix.
   * \return the trace of matrix M.
   * \test PASSED
   */
  friend value_type trace(const self& M) {
    return static_cast<value_type>(M.rowCount);
  }

  /**
   * Appends a matrix to another.
   * \param lhs the matrix to which 'rhs' will be appended to.
   * \param rhs the matrix to append to 'lhs'.
   * \test PASSED
   */
  friend void append_block_diag(self& lhs, const self& rhs) {
    lhs.rowCount += rhs.rowCount;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(rowCount);
  }
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(rowCount);
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

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_readable_matrix<
    mat<T, mat_structure::identity, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type =
      is_readable_matrix<mat<T, mat_structure::identity, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_writable_matrix<
    mat<T, mat_structure::identity, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = false;
  using type =
      is_writable_matrix<mat<T, mat_structure::identity, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_resizable_matrix<
    mat<T, mat_structure::identity, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type = is_resizable_matrix<
      mat<T, mat_structure::identity, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct has_allocator_matrix<
    mat<T, mat_structure::identity, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = false;
  using type = has_allocator_matrix<
      mat<T, mat_structure::identity, Alignment, Allocator>>;
};

/**
 * Scalar multiplication, always results in a scalar matrix.
 * \param M some matrix.
 * \param S some scalar.
 * \return Scalar matrix.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
std::enable_if_t<!is_readable_vector_v<T> && !is_readable_matrix_v<T>,
                 mat<T, mat_structure::scalar, Alignment, Allocator>>
operator*(const mat<T, mat_structure::identity, Alignment, Allocator>& M,
          const T& S) {
  return mat<T, mat_structure::scalar, Alignment, Allocator>(M.get_row_count(),
                                                             S);
}

/**
 * Scalar multiplication, always results in a scalar matrix.
 * \param S some scalar.
 * \param M a null-matrix.
 * \return Scalar matrix.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
std::enable_if_t<!is_readable_vector_v<T> && !is_readable_matrix_v<T>,
                 mat<T, mat_structure::scalar, Alignment, Allocator>>
operator*(const T& S,
          const mat<T, mat_structure::identity, Alignment, Allocator>& M) {
  return mat<T, mat_structure::scalar, Alignment, Allocator>(M.get_row_count(),
                                                             S);
}

}  // namespace ReaK

#endif

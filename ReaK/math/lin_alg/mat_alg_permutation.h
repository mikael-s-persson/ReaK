/**
 * \file mat_alg_permutation.h
 *
 * This library declares matrix specializations for representing and manipulating permutation matrices.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
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

#ifndef REAK_MATH_LIN_ALG_MAT_ALG_PERMUTATION_H_
#define REAK_MATH_LIN_ALG_MAT_ALG_PERMUTATION_H_

#include "ReaK/math/lin_alg/mat_alg_general.h"
#include "ReaK/math/lin_alg/mat_alg_identity.h"

#include <type_traits>

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
class mat<T, mat_structure::permutation, Alignment, Allocator>
    : public serializable {
 public:
  using self = mat<T, mat_structure::permutation, Alignment, Allocator>;
  using allocator_type = Allocator;

  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using container_type = std::vector<size_type>;

  using reference = void;
  using const_reference = T;
  using pointer = void;
  using const_pointer = void;

  using row_iterator = void;
  using const_row_iterator = void;
  using col_iterator = void;
  using const_col_iterator = void;

  static constexpr std::size_t static_row_count = 0;
  static constexpr std::size_t static_col_count = 0;
  static constexpr mat_alignment::tag alignment = Alignment;
  static constexpr mat_structure::tag structure = mat_structure::permutation;

 private:
  container_type idx;
  size_type rowCount;

 public:
  /**
   * Default constructor. Sets dimensions to zero.
   */
  mat() : rowCount(0) {}
  /**
   * Constructs an identity matrix to the given dimensions.
   */
  explicit mat(size_type aRowCount) : idx(aRowCount), rowCount(aRowCount) {
    for (size_type i = 0; i < rowCount; ++i) {
      idx[i] = i;
    }
  }

  /**
   * Standard swap function (works with ADL).
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    lhs.idx.swap(rhs.idx);
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
    if (idx[i] == j) {
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
    idx.resize(rowCount);
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
    idx.resize(rowCount);
  }

  /**
   * Get the allocator object of the underlying container.
   * \return The allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return idx.get_allocator(); }

  /**
   * Get the index of the row that should be in-place of the given row index.
   * Note that for column swap this mapping can simply be applied in reverse, that is,
   * the returned index is the index of the original column that should appear at the
   * given destination column index.
   * \param i The index of the original row (or destination column).
   * \return the index of the destination row (or original column).
   */
  int operator[](int i) const { return idx[i]; }

  /**
   * Append a row swap to this permutation matrix.
   * Essentially equivalent to permute(i,j) * (this).
   * \param i The first row involved in the swap.
   * \param j The second row involved in the swap.
   */
  void add_row_swap(int i, int j) {
    if (i == j) {
      return;
    }
    using std::swap;
    swap(idx[i], idx[j]);
  }

  /**
   * Append a column swap to this permutation matrix.
   * Essentially equivalent to (this) * permute(i,j).
   * Note also that it is generally more efficient to add row swaps (for example, if
   * you want to accumulated many column swaps, it is more efficient to accumulated them
   * as row swaps and then invert (or transpose) the final permutation matrix).
   * \param i The first column involved in the swap.
   * \param j The second column involved in the swap.
   */
  void add_column_swap(int i, int j) {
    if (i == j) {
      return;
    }
    using std::swap;
    int p_i = 0;
    int p_j = 0;
    for (int k = 0; k < rowCount; ++k) {
      if (idx[k] == i) {
        p_i = k;
      }
      if (idx[k] == j) {
        p_j = k;
      }
    }
    swap(idx[p_i], idx[p_j]);
  }

  /**
   * Negate the matrix.
   * \return This matrix, by constant reference.
   * \test PASSED
   */
  mat<value_type, mat_structure::scalar> operator-() const {
    return mat<value_type, mat_structure::scalar>(rowCount, value_type(-1));
  }

  friend self operator*(const self& M1, const self& M2) {
    if (M1.get_row_count() != M2.get_row_count()) {
      throw std::range_error("Matrix dimensions mismatch!");
    }
    self result(M1.get_row_count());
    for (int i = 0; i < M1.get_row_count(); ++i) {
      result.idx[i] = M2.idx[M1.idx[i]];
    }
    return result;
  }

  /**
   * Transpose the matrix.
   * \param rhs the matrix to be transposed.
   * \return The transpose matrix, by value.
   * \test PASSED
   */
  friend self transpose(const self& rhs) {
    self result(rhs.rowCount);
    for (int i = 0; i < rhs.rowCount; ++i) {
      result.idx[rhs.idx[i]] = i;
    }
    return result;
  }

  /**
   * Transpose and move the matrix.
   * \param rhs the matrix to be transposed and moved (emptied).
   * \return The transpose matrix, by value.
   * \test PASSED
   */
  friend self transpose_move(const self& rhs) { return transpose(rhs); }

  /**
   * Invert the matrix.
   * \param rhs the matrix to be inverted.
   * \return The inverse matrix, by value.
   * \test PASSED
   */
  friend self invert(const self& rhs) { return transpose(rhs); }

  /**
   * Returns the trace of a matrix.
   * \param M A matrix.
   * \return the trace of matrix M.
   * \test PASSED
   */
  friend value_type trace(const self& M) {
    value_type result(0);
    for (int i = 0; i < M.rowCount; ++i) {
      if (M.idx[i] == i) {
        result += value_type(1);
      }
    }
    return result;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(idx) &
        std::pair<std::string, size_type>("rowCount", rowCount);
  }
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(idx) &
        std::pair<std::string, size_type&>("rowCount", rowCount);
  }

  RK_RTTI_REGISTER_CLASS_1BASE(self, 1, serializable)
};

template <typename T>
struct mat_permutation {
  using type = mat<T, mat_structure::permutation>;
};

template <typename T>
using mat_permutation_t = typename mat_permutation<T>::type;

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_readable_matrix<
    mat<T, mat_structure::permutation, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type = is_readable_matrix<
      mat<T, mat_structure::permutation, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_writable_matrix<
    mat<T, mat_structure::permutation, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = false;
  using type = is_writable_matrix<
      mat<T, mat_structure::permutation, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_resizable_matrix<
    mat<T, mat_structure::permutation, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type = is_resizable_matrix<
      mat<T, mat_structure::permutation, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct has_allocator_matrix<
    mat<T, mat_structure::permutation, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type = has_allocator_matrix<
      mat<T, mat_structure::permutation, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_square_matrix<
    mat<T, mat_structure::permutation, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type = is_square_matrix<
      mat<T, mat_structure::permutation, Alignment, Allocator>>;
};

/**
 * Column-vector multiplication, returns a permutated vector.
 * \param M some permutation matrix.
 * \param V some vector.
 * \return A permutated vector.
 * \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
 */
template <typename T, typename Vector, mat_alignment::tag Alignment,
          typename Allocator>
std::enable_if_t<is_readable_vector_v<Vector>, vect_copy_t<Vector>> operator*(
    const mat<T, mat_structure::permutation, Alignment, Allocator>& M,
    const Vector& V) {
  if (V.size() != M.get_col_count()) {
    throw std::range_error("Matrix dimension mismatch.");
  }

  vect_copy_t<Vector> result(V);
  for (int i = 0; i < result.size(); ++i) {
    result[i] = V[M[i]];
  }
  return result;
}

/**
 * Row-vector multiplication with permutation matrix, returns a permutated vector.
 * \param V some row-vector.
 * \param M some permutation matrix.
 * \return A permutated vector.
 * \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
 */
template <typename T, typename Vector, mat_alignment::tag Alignment,
          typename Allocator>
std::enable_if_t<is_readable_vector_v<Vector>, vect_copy_t<Vector>> operator*(
    const Vector& V,
    const mat<T, mat_structure::permutation, Alignment, Allocator>& M) {
  if (V.size() != M.get_row_count()) {
    throw std::range_error("Matrix dimension mismatch.");
  }

  vect_copy_t<Vector> result(V);
  for (int i = 0; i < result.size(); ++i) {
    result[M[i]] = V[i];
  }
  return result;
}

/**
 * Scalar multiplication of a permutation matrix.
 * \param M some permutation matrix.
 * \param S some scalar.
 * \return A square matrix.
 * \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
std::enable_if_t<!is_readable_vector_v<T> && !is_readable_matrix_v<T>,
                 mat<T, mat_structure::square, Alignment, Allocator>>
operator*(const mat<T, mat_structure::permutation, Alignment, Allocator>& M,
          const T& S) {
  mat<T, mat_structure::square, Alignment, Allocator> result(M.get_row_count());
  for (int i = 0; i < M.get_row_count(); ++i) {
    result(i, M[i]) = S;
  }
  return result;
}

/**
 * Scalar multiplication of a permutation matrix.
 * \param S some scalar.
 * \param M some permutation matrix.
 * \return A square matrix.
 * \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
std::enable_if_t<!is_readable_vector_v<T> && !is_readable_matrix_v<T>,
                 mat<T, mat_structure::square, Alignment, Allocator>>
operator*(const T& S,
          const mat<T, mat_structure::permutation, Alignment, Allocator>& M) {
  mat<T, mat_structure::square, Alignment, Allocator> result(M.get_row_count());
  for (int i = 0; i < M.get_row_count(); ++i) {
    result(i, M[i]) = S;
  }
  return result;
}

/**
 * Matrix multiplication with permutation matrix (column permutation).
 * \param M1 some matrix.
 * \param M2 a permutation-matrix.
 * \return A permuted matrix.
 * \throw std::range_error if matrices' dimensions are not proper for multiplication.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator,
          typename Matrix>
std::enable_if_t<!is_readable_matrix_v<Matrix> &&
                     (mat_product_priority_v<Matrix> <
                      mat_product_priority_v<mat<T, mat_structure::permutation,
                                                 Alignment, Allocator>>),
                 mat<T, mat_structure::rectangular, Alignment, Allocator>>
operator*(const Matrix& M1,
          const mat<T, mat_structure::permutation, Alignment, Allocator>& M2) {
  if (M1.get_col_count() != M2.get_row_count()) {
    throw std::range_error("Matrix dimension mismatch.");
  }
  mat<T, mat_structure::rectangular, Alignment, Allocator> result(
      M1.get_row_count(), M1.get_col_count());
  for (int i = 0; i < result.get_col_count(); ++i) {
    for (int j = 0; j < result.get_row_count(); ++j) {
      result(j, M2[i]) = M1(j, i);
    }
  }
  return result;
}

/**
 * Matrix multiplication with permutation matrix (row permutation).
 * \param M1 a permutation-matrix.
 * \param M2 some matrix.
 * \return A permuted matrix.
 * \throw std::range_error if matrices' dimensions are not proper for multiplication.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator,
          typename Matrix>
std::enable_if_t<is_readable_matrix_v<Matrix> &&
                     (mat_product_priority_v<Matrix> <
                      mat_product_priority_v<mat<T, mat_structure::permutation,
                                                 Alignment, Allocator>>),
                 mat<T, mat_structure::rectangular, Alignment, Allocator>>
operator*(const mat<T, mat_structure::permutation, Alignment, Allocator>& M1,
          const Matrix& M2) {
  if (M1.get_col_count() != M2.get_row_count()) {
    throw std::range_error("Matrix dimension mismatch.");
  }
  mat<T, mat_structure::rectangular, Alignment, Allocator> result(
      M2.get_row_count(), M2.get_col_count());
  for (int j = 0; j < result.get_row_count(); ++j) {
    for (int i = 0; i < result.get_col_count(); ++i) {
      result(j, i) = M2(M1[j], i);
    }
  }
  return result;
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_ALG_PERMUTATION_H_

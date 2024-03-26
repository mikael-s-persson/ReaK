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
#include "ReaK/math/lin_alg/mat_op_results.h"

#include <type_traits>

namespace ReaK {

/**
 * This class implements a place-holder or interface-implementation to represent
 * a identity matrix (all entries zero except diagonal is all unity). This is useful to build for example a
 * block-matrix with some identity-matrix blocks, and, of course, requires minimal storage space.
 *
 * Models: ReadableMatrix and ResizableMatrix.
 *
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or
 *mat_alignment::column_major (default).
 */
template <typename T, mat_alignment::tag Alignment, unsigned int RowCount>
class mat<T, mat_structure::permutation, Alignment, RowCount, RowCount>
    : public serializable {
 public:
  using self =
      mat<T, mat_structure::permutation, Alignment, RowCount, RowCount>;
  static constexpr bool is_dynamic_size = (RowCount == 0);

  using value_type = T;
  using size_type = int;
  using difference_type = int;
  using container_type = std::conditional_t<is_dynamic_size, std::vector<int>,
                                            std::array<int, RowCount>>;

  using reference = void;
  using const_reference = T;
  using pointer = void;
  using const_pointer = void;

  using row_iterator = void;
  using const_row_iterator = void;
  using col_iterator = void;
  using const_col_iterator = void;

  static constexpr unsigned int static_row_count = RowCount;
  static constexpr unsigned int static_col_count = RowCount;
  static constexpr mat_alignment::tag alignment = Alignment;
  static constexpr mat_structure::tag structure = mat_structure::permutation;

 private:
  struct dynamic_data {
    container_type idx = {};
    int rowCount = 0;
  };
  struct static_data {
    container_type idx = {};
  };
  /// Array which holds all the permuted indices of the matrix (dimension: rowCount).
  std::conditional_t<is_dynamic_size, dynamic_data, static_data> data;

 public:
  /// Default constructor. Sets dimensions to zero.
  mat() = default;
  mat(const self& rhs) = default;
  mat(self&& rhs) = default;
  self& operator=(const self& rhs) = default;
  self& operator=(self&& rhs) = default;
  ~mat() override = default;
  /// Constructs an identity matrix to the given dimensions.
  explicit mat(int aRowCount) {
    if constexpr (is_dynamic_size) {
      data.idx.resize(aRowCount);
      data.rowCount = aRowCount;
    } else {
      if (aRowCount != RowCount) {
        throw std::range_error("Row count mismatch!");
      }
    }
    for (int i = 0; i < aRowCount; ++i) {
      data.idx[i] = i;
    }
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Matrix indexing accessor for read-only access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  const_reference operator()(int i, int j) const {
    if (data.idx[i] == j) {
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
      data.idx.resize(aRowCount);
      if (aRowCount < data.rowCount) {
        // Only choice is to reset to identity.
        for (int i = 0; i < aRowCount; ++i) {
          data.idx[i] = i;
        }
      } else if (aRowCount > data.rowCount) {
        // Set new elements to identity.
        for (int i = data.rowCount; i < aRowCount; ++i) {
          data.idx[i] = i;
        }
      }
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

  /// Get the index of the row that should be in-place of the given row index.
  /// Note that for column swap this mapping can simply be applied in reverse, that is,
  /// the returned index is the index of the original column that should appear at the
  /// given destination column index.
  /// \param i The index of the original row (or destination column).
  /// \return the index of the destination row (or original column).
  int operator[](int i) const { return data.idx[i]; }

  /// Append a row swap to this permutation matrix.
  /// Essentially equivalent to permute(i,j) * (this).
  /// \param i The first row involved in the swap.
  /// \param j The second row involved in the swap.
  void add_row_swap(int i, int j) {
    if (i == j) {
      return;
    }
    using std::swap;
    swap(data.idx[i], data.idx[j]);
  }

  /// Append a column swap to this permutation matrix.
  /// Essentially equivalent to (this) * permute(i,j).
  /// Note also that it is generally more efficient to add row swaps (for example, if
  /// you want to accumulated many column swaps, it is more efficient to accumulated them
  /// as row swaps and then invert (or transpose) the final permutation matrix).
  /// \param i The first column involved in the swap.
  /// \param j The second column involved in the swap.
  void add_column_swap(int i, int j) {
    if (i == j) {
      return;
    }
    using std::swap;
    int p_i = 0;
    int p_j = 0;
    for (int k = 0; k < get_row_count(); ++k) {
      if (data.idx[k] == i) {
        p_i = k;
      }
      if (data.idx[k] == j) {
        p_j = k;
      }
    }
    swap(data.idx[p_i], data.idx[p_j]);
  }

  /// Negate the matrix. Loses the permutation structure.
  auto operator-() const {
    return -mat<value_type, mat_structure::square, Alignment, RowCount,
                RowCount>(*this);
  }

  friend self operator*(const self& M1, const self& M2) {
    if constexpr (is_dynamic_size) {
      if (M1.get_row_count() != M2.get_row_count()) {
        throw std::range_error("Matrix dimensions mismatch!");
      }
    }
    self result(M1.get_row_count());
    for (int i = 0; i < M1.get_row_count(); ++i) {
      result.data.idx[i] = M2.data.idx[M1.data.idx[i]];
    }
    return result;
  }

  /// Transpose the matrix.
  /// \param rhs the matrix to be transposed.
  /// \return The transpose matrix, by value.
  friend self transpose(const self& rhs) {
    self result(rhs.get_row_count());
    for (int i = 0; i < rhs.get_row_count(); ++i) {
      result.data.idx[rhs.data.idx[i]] = i;
    }
    return result;
  }

  /// Invert the matrix.
  /// \param rhs the matrix to be inverted.
  /// \return The inverse matrix, by value.
  friend self invert(const self& rhs) { return transpose(rhs); }

  /// Returns the trace of a matrix.
  /// \param M A matrix.
  /// \return the trace of matrix M.
  friend value_type trace(const self& M) {
    value_type result(0);
    for (int i = 0; i < M.rowCount; ++i) {
      if (M.data.idx[i] == i) {
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
    A& RK_SERIAL_SAVE_WITH_ALIAS("idx", data.idx);
    if constexpr (is_dynamic_size) {
      A& RK_SERIAL_SAVE_WITH_ALIAS("rowCount", data.rowCount);
    }
  }
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_ALIAS("idx", data.idx);
    if constexpr (is_dynamic_size) {
      A& RK_SERIAL_LOAD_WITH_ALIAS("rowCount", data.rowCount);
    }
  }

  RK_RTTI_REGISTER_CLASS_1BASE(self, 1, serializable)
};

template <typename T>
struct mat_permutation {
  using type = mat<T, mat_structure::permutation>;
};

template <typename T>
using mat_permutation_t = typename mat_permutation<T>::type;

/// Column-vector multiplication, returns a permutated vector.
/// \param M some permutation matrix.
/// \param V some vector.
/// \return A permutated vector.
/// \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
template <typename T, ReadableVector Vector, mat_alignment::tag Alignment,
          unsigned int RowCount>
vect_copy_t<Vector> operator*(
    const mat<T, mat_structure::permutation, Alignment, RowCount, RowCount>& M,
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

/// Row-vector multiplication with permutation matrix, returns a permutated vector.
/// \param V some row-vector.
/// \param M some permutation matrix.
/// \return A permutated vector.
/// \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
template <typename T, ReadableVector Vector, mat_alignment::tag Alignment,
          unsigned int RowCount>
vect_copy_t<Vector> operator*(const Vector& V,
                              const mat<T, mat_structure::permutation,
                                        Alignment, RowCount, RowCount>& M) {
  if (V.size() != M.get_row_count()) {
    throw std::range_error("Matrix dimension mismatch.");
  }

  vect_copy_t<Vector> result(V);
  for (int i = 0; i < result.size(); ++i) {
    result[M[i]] = V[i];
  }
  return result;
}

/// Scalar multiplication of a permutation matrix.
/// \param M some permutation matrix.
/// \param S some scalar.
/// \return A square matrix.
/// \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
template <typename T, mat_alignment::tag Alignment, unsigned int RowCount>
std::enable_if_t<!ReadableVector<T> && !ReadableMatrix<T>,
                 mat<T, mat_structure::square, Alignment, RowCount, RowCount>>
operator*(
    const mat<T, mat_structure::permutation, Alignment, RowCount, RowCount>& M,
    const T& S) {
  mat<T, mat_structure::square, Alignment, RowCount, RowCount> result(
      M.get_row_count());
  for (int i = 0; i < M.get_row_count(); ++i) {
    result(i, M[i]) = S;
  }
  return result;
}

/// Scalar multiplication of a permutation matrix.
/// \param S some scalar.
/// \param M some permutation matrix.
/// \return A square matrix.
/// \throw std::range_error if matrix and vector dimensions are not proper for multiplication.
template <typename T, mat_alignment::tag Alignment, unsigned int RowCount>
std::enable_if_t<!ReadableVector<T> && !ReadableMatrix<T>,
                 mat<T, mat_structure::square, Alignment, RowCount, RowCount>>
operator*(const T& S, const mat<T, mat_structure::permutation, Alignment,
                                RowCount, RowCount>& M) {
  mat<T, mat_structure::square, Alignment, RowCount, RowCount> result(
      M.get_row_count());
  for (int i = 0; i < M.get_row_count(); ++i) {
    result(i, M[i]) = S;
  }
  return result;
}

/// Matrix multiplication with permutation matrix (column permutation).
/// \param M1 some matrix.
/// \param M2 a permutation-matrix.
/// \return A permuted matrix.
/// \throw std::range_error if matrices' dimensions are not proper for multiplication.
template <typename T, mat_alignment::tag Alignment, unsigned int RowCount,
          ReadableMatrix Matrix>
std::enable_if_t<
    (mat_product_priority_v<Matrix> <
     mat_product_priority_v<
         mat<T, mat_structure::permutation, Alignment, RowCount, RowCount>>),
    mat_product_result_t<Matrix, mat<T, mat_structure::permutation, Alignment,
                                     RowCount, RowCount>>>
operator*(const Matrix& M1, const mat<T, mat_structure::permutation, Alignment,
                                      RowCount, RowCount>& M2) {
  if (M1.get_col_count() != M2.get_row_count()) {
    throw std::range_error("Matrix dimension mismatch.");
  }
  mat_product_result_t<
      Matrix, mat<T, mat_structure::permutation, Alignment, RowCount, RowCount>>
      result(M1.get_row_count(), M1.get_col_count());
  for (int i = 0; i < result.get_col_count(); ++i) {
    for (int j = 0; j < result.get_row_count(); ++j) {
      result(j, M2[i]) = M1(j, i);
    }
  }
  return result;
}

/// Matrix multiplication with permutation matrix (row permutation).
/// \param M1 a permutation-matrix.
/// \param M2 some matrix.
/// \return A permuted matrix.
/// \throw std::range_error if matrices' dimensions are not proper for multiplication.
template <typename T, mat_alignment::tag Alignment, unsigned int RowCount,
          ReadableMatrix Matrix>
std::enable_if_t<(mat_product_priority_v<Matrix> <
                  mat_product_priority_v<mat<T, mat_structure::permutation,
                                             Alignment, RowCount, RowCount>>),
                 mat_product_result_t<mat<T, mat_structure::permutation,
                                          Alignment, RowCount, RowCount>,
                                      Matrix>>
operator*(
    const mat<T, mat_structure::permutation, Alignment, RowCount, RowCount>& M1,
    const Matrix& M2) {
  if (M1.get_col_count() != M2.get_row_count()) {
    throw std::range_error("Matrix dimension mismatch.");
  }
  mat_product_result_t<
      mat<T, mat_structure::permutation, Alignment, RowCount, RowCount>, Matrix>
      result(M2.get_row_count(), M2.get_col_count());
  for (int j = 0; j < result.get_row_count(); ++j) {
    for (int i = 0; i < result.get_col_count(); ++i) {
      result(j, i) = M2(M1[j], i);
    }
  }
  return result;
}

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_ALG_PERMUTATION_H_

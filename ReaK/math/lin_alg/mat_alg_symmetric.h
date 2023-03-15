/**
 * \file mat_alg_symmetric.h
 *
 * This library implements the specialization of the mat<> template for a
 * general symmetric matrix (dynamic dimension). This matrix type fulfills the matrix
 * concepts of Readable, Writable, and Resizable, but not FullyWritable.
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

#ifndef REAK_MATH_LIN_ALG_MAT_ALG_SYMMETRIC_H_
#define REAK_MATH_LIN_ALG_MAT_ALG_SYMMETRIC_H_

#include "ReaK/math/lin_alg/mat_alg_general.h"

#include <type_traits>

namespace ReaK {

/// This class holds a symmetric matrix. This class will hold only the upper-triangular part
/// since the lower part is assumed to be equal to the upper one.
///
/// Models: ReadableMatrixConcept, WritableMatrixConcept, and ResizableMatrixConcept.
///
/// \tparam T Arithmetic type of the elements of the matrix.
/// \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or
/// mat_alignment::column_major (default).
template <typename T, mat_alignment::tag Alignment, unsigned int RowCount>
class mat<T, mat_structure::symmetric, Alignment, RowCount, RowCount>
    : public serializable {
 public:
  using self = mat<T, mat_structure::symmetric, Alignment, RowCount, RowCount>;
  static constexpr bool is_dynamic_size = (RowCount == 0);

  using value_type = T;
  using container_type =
      std::conditional_t<is_dynamic_size, std::vector<value_type>,
                         std::array<value_type, RowCount * RowCount>>;

  using reference = typename container_type::reference;
  using const_reference = typename container_type::const_reference;
  using pointer = typename container_type::pointer;
  using const_pointer = typename container_type::const_pointer;

  using col_iterator = void;
  using const_col_iterator = void;
  using row_iterator = void;
  using const_row_iterator = void;

  using size_type = int;
  using difference_type = int;

  static constexpr unsigned int static_row_count = RowCount;
  static constexpr unsigned int static_col_count = RowCount;
  static constexpr mat_alignment::tag alignment = Alignment;
  static constexpr mat_structure::tag structure = mat_structure::symmetric;

  template <typename OtherT, mat_structure::tag OtherStructure,
            mat_alignment::tag OtherAlignment, unsigned int OtherRowCount,
            unsigned int OtherColCount>
  friend class mat;

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

  static int mat_triangular_size(int Size) {
    return (Size * (Size - 1)) / 2 + Size;
  }

 public:
  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/
  /// Default constructor: sets all to zero.
  mat() : data() {}

  /// Constructor for a sized matrix.
  explicit mat(int aRowCount, const value_type& aFill = value_type(0.0)) {
    if constexpr (is_dynamic_size) {
      data.q.resize(mat_triangular_size(aRowCount), aFill);
      data.rowCount = aRowCount;
    } else {
      if (aRowCount != RowCount) {
        throw std::range_error("Row count mismatch!");
      }
      std::fill(data.q.begin(), data.q.end(), aFill);
    }
  }

  /// Constructor for an identity matrix.
  mat(int aRowCount, bool aIdentity) : mat(aRowCount) {
    if (aIdentity) {
      int k = 0;
      for (int i = 0; i < get_row_count(); k += ++i) {
        data.q[k + i] = 1.0;
      }
    }
  }

  /// Default Copy/Move Constructors.
  mat(const self& rhs) = default;
  mat(self&& rhs) = default;
  self& operator=(const self& rhs) = default;
  self& operator=(self&& rhs) = default;
  ~mat() override = default;

  /// The standard swap function (works with ADL).
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.data.q, rhs.data.q);
    if constexpr (is_dynamic_size) {
      swap(lhs.data.rowCount, rhs.data.rowCount);
    }
  }

  /// Explicit constructor from any type of matrix. The "(M + M.transpose) / 2" is applied to guarantee symmetry.
  template <typename Matrix>
  explicit mat(const Matrix& M,
               std::enable_if_t<is_readable_matrix_v<Matrix> &&
                                    !is_symmetric_matrix_v<Matrix> &&
                                    !std::is_same_v<Matrix, self>,
                                void*>
                   dummy = nullptr)
      : mat(std::max(M.get_row_count(), M.get_col_count())) {
    int k = 0;
    int i = 0;
    const int min_size = std::min(M.get_row_count(), M.get_col_count());
    for (; i < min_size; k += ++i) {
      for (int j = 0; j < i; ++j) {
        data.q[k + j] = value_type(0.5) * (M(j, i) + M(i, j));
      }
      data.q[k + i] = M(i, i);
    }
    if (M.get_row_count() > M.get_col_count()) {
      for (; i < get_row_count(); k += ++i) {
        for (int j = 0; j < min_size; ++j) {
          data.q[k + j] = value_type(0.5) * M(i, j);
        }
      }
    } else {
      for (; i < get_row_count(); k += ++i) {
        for (int j = 0; j < min_size; ++j) {
          data.q[k + j] = value_type(0.5) * M(j, i);
        }
      }
    }
  }

  /// Explicit constructor from any type of matrix. The "(M + M.transpose) / 2" is applied to guarantee symmetry.
  template <typename Matrix>
  explicit mat(const Matrix& M,
               std::enable_if_t<is_readable_matrix_v<Matrix> &&
                                    is_symmetric_matrix_v<Matrix> &&
                                    !std::is_same_v<Matrix, self>,
                                void*>
                   dummy = nullptr)
      : mat(M.get_row_count()) {
    int k = 0;
    int i = 0;
    for (; i < get_row_count(); k += ++i) {
      for (int j = 0; j < i; ++j) {
        data.q[k + j] = M(i, j);
      }
      data.q[k + i] = M(i, i);
    }
  }

  /// Constructs a 2x2 symmetric matrix from three elements.
  mat(const_reference a11, const_reference a12, const_reference a22) : mat(2) {
    data.q[0] = a11;
    data.q[1] = a12;
    data.q[2] = a22;
  }

  /// Constructs a 3x3 symmetric matrix from six elements.
  mat(const_reference a11, const_reference a12, const_reference a13,
      const_reference a22, const_reference a23, const_reference a33)
      : mat(3) {
    data.q[0] = a11;
    data.q[1] = a12;
    data.q[2] = a22;
    data.q[3] = a13;
    data.q[4] = a23;
    data.q[5] = a33;
  }

  /// Constructs a 4x4 symmetric matrix from ten elements.
  mat(const_reference a11, const_reference a12, const_reference a13,
      const_reference a14, const_reference a22, const_reference a23,
      const_reference a24, const_reference a33, const_reference a34,
      const_reference a44)
      : mat(4) {
    data.q[0] = a11;
    data.q[1] = a12;
    data.q[2] = a22;
    data.q[3] = a13;
    data.q[4] = a23;
    data.q[5] = a33;
    data.q[6] = a14;
    data.q[7] = a24;
    data.q[8] = a34;
    data.q[9] = a44;
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Matrix indexing accessor for read-write access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  reference operator()(int i, int j) {
    if (i > j) {
      return data.q[mat_triangular_size(i) + j];
    }
    return data.q[mat_triangular_size(j) + i];
  }

  /// Matrix indexing accessor for read-only access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  const_reference operator()(int i, int j) const {
    if (i > j) {
      return data.q[mat_triangular_size(i) + j];
    }
    return data.q[mat_triangular_size(j) + i];
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
      data.q.resize(mat_triangular_size(aRowCount), value_type(0.0));
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

  /// Gets the row-count and column-count of the matrix, as a std::pair of values.
  /// \return the row-count and column-count of the matrix, as a std::pair of values.
  std::pair<int, int> size() const noexcept {
    return std::make_pair(get_row_count(), get_row_count());
  }
  /// Sets the row-count and column-count of the matrix, via a std::pair of dimension values.
  /// \param sz new dimensions for the matrix.
  void resize(const std::pair<int, int>& sz) { set_row_count(sz.first); }

  /*******************************************************************************
                           Assignment Operators
  *******************************************************************************/

  /// Standard Assignment operator with a matrix of any type. The "(M + M.transpose) / 2" formula is applied to guarantee
  /// symmetry.
  template <typename Matrix>
  self& operator=(const Matrix& M) {
    self tmp(M);
    swap(*this, tmp);
    return *this;
  }

  /// Add-and-store operator with standard semantics.
  template <mat_alignment::tag Align2>
  self& operator+=(const mat<value_type, mat_structure::symmetric, Align2,
                             RowCount, RowCount>& M) {
    if (M.get_row_count() != get_row_count()) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int k = 0;
    for (int i = 0; i < get_row_count(); k += ++i) {
      for (int j = 0; j <= i; ++j) {
        data.q[k + j] += M(i, j);
      }
    }
    return *this;
  }

  /// Add-and-store operator with standard semantics.
  template <mat_alignment::tag Align2>
  self& operator+=(const mat<value_type, mat_structure::diagonal, Align2,
                             RowCount, RowCount>& M) {
    if (M.get_row_count() != get_row_count()) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int k = 0;
    for (int i = 0; i < get_row_count(); k += ++i) {
      data.q[k + i] += M(i, i);
    }
    return *this;
  }

  /// Sub-and-store operator with standard semantics.
  template <mat_alignment::tag Align2>
  self& operator-=(const mat<value_type, mat_structure::symmetric, Align2,
                             RowCount, RowCount>& M) {
    if (M.get_row_count() != get_row_count()) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int k = 0;
    for (int i = 0; i < get_row_count(); k += ++i) {
      for (int j = 0; j <= i; ++j) {
        data.q[k + j] -= M(i, j);
      }
    }
    return *this;
  }

  /// Add-and-store operator with standard semantics.
  template <mat_alignment::tag Align2>
  self& operator-=(const mat<value_type, mat_structure::diagonal, Align2,
                             RowCount, RowCount>& M) {
    if (M.get_row_count() != get_row_count()) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int k = 0;
    for (int i = 0; i < get_row_count(); k += ++i) {
      data.q[k + i] -= M(i, i);
    }
    return *this;
  }

  /// Scalar-multiply-and-store operator with standard semantics.
  self& operator*=(const T& S) {
    for (auto& v : data.q) {
      v *= S;
    }
    return *this;
  }

  /// Negation operator for any type of matrices. This is a default operator
  /// that will be called if no better special-purpose overload exists.
  /// \return Symmetric matrix.
  self operator-() const {
    self result(*this);
    auto itr = result.data.q.begin();
    for (auto it = data.q.begin(); it != data.q.end(); ++it, ++itr) {
      *itr = -(*it);
    }
    return result;
  }

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /// Add two matrices.
  /// \param M1 the first matrix (first operand).
  /// \param M2 the other matrix (second operand).
  /// \return symmetric matrix sum of M1 and M2.
  /// \throw std::range_error if the two matrix dimensions don't match.
  friend self operator+(const self& M1, const self& M2) {
    if constexpr (is_dynamic_size) {
      if (M1.get_row_count() != M2.get_row_count()) {
        throw std::range_error("Matrix dimension mismatch.");
      }
    }
    self result(M1.get_row_count());
    int k = 0;
    for (int i = 0; i < M1.get_row_count(); k += ++i) {
      for (int j = 0; j <= i; ++j) {
        result.data.q[k + j] = M1.data.q[k + j] + M2.data.q[k + j];
      }
    }
    return result;
  }

  /// Sub two matrices.
  /// \param M1 the first matrix (first operand).
  /// \param M2 the other matrix (second operand).
  /// \return symmetric matrix difference of M1 and M2.
  /// \throw std::range_error if the two matrix dimensions don't match.
  friend self operator-(const self& M1, const self& M2) {
    if constexpr (is_dynamic_size) {
      if (M1.get_row_count() != M2.get_row_count()) {
        throw std::range_error("Matrix dimension mismatch.");
      }
    }
    self result(M1.get_row_count());
    int k = 0;
    for (int i = 0; i < M1.get_row_count(); k += ++i) {
      for (int j = 0; j <= i; ++j) {
        result.data.q[k + j] = M1.data.q[k + j] - M2.data.q[k + j];
      }
    }
    return result;
  }

  /// General Matrix multiplication.
  /// \param M the other matrix (second operand).
  /// \return general matrix multiplication result of this and M.
  /// \throw std::range_error if the two matrix dimensions don't match.
  template <typename Matrix>
  auto multiply_this_and_dense_mat(const Matrix& M2) const {
    mat_product_result_t<self, Matrix> result{M2.get_row_count(),
                                              M2.get_col_count()};
    for (int k = 0, i = 0; i < get_row_count(); k += ++i) {
      for (int l = 0; l < M2.get_col_count(); ++l) {
        for (int j = 0; j < i; ++j) {
          result(j, l) += data.q[k + j] * M2(i, l);
          result(i, l) += data.q[k + j] * M2(j, l);
        }
        result(i, l) += data.q[k + i] * M2(i, l);
      }
    }
    return result;  // NRVO
  }

  /// General Matrix multiplication.
  /// \param M the other matrix (second operand).
  /// \return general matrix multiplication result of this and M.
  /// \throw std::range_error if the two matrix dimensions don't match.
  template <typename Matrix>
  auto multiply_dense_and_this_mat(const Matrix& M1) const {
    mat_product_result_t<Matrix, self> result{M1.get_row_count(),
                                              M1.get_col_count()};
    for (int k = 0, i = 0; i < M1.get_col_count(); k += ++i) {
      for (int l = 0; l < M1.get_row_count(); ++l) {
        for (int j = 0; j < i; ++j) {
          result(l, i) += data.q[k + j] * M1(l, j);
          result(l, j) += data.q[k + j] * M1(l, i);
        }
        result(l, i) += data.q[k + i] * M1(l, i);
      }
    }
    return result;  // NRVO
  }

  /// Symmetric Matrix multiplication.
  /// \param M1 the first symmetric matrix (first operand).
  /// \param M2 the other symmetric matrix (second operand).
  /// \return square matrix, result of M1 times M2.
  /// \throw std::range_error if the two matrix dimensions don't match.
  auto multiply_with_same_mat(const self& M2) const {
    mat_product_result_t<self, self> result{get_row_count()};
    int k = 0;
    int i = 0;
    for (; i < get_row_count(); k += ++i) {
      int h = 0;
      int l = 0;
      for (; l <= i; h += ++l) {
        int m = 0;
        int j = 0;
        for (; j < l; m += ++j) {
          result(j, l) += data.q[k + j] * M2.data.q[k + l];
          result(i, l) += data.q[k + j] * M2.data.q[h + j];
        }
        for (; j < i; m += ++j) {
          result(j, l) += data.q[k + j] * M2.data.q[k + l];
          result(i, l) += data.q[k + j] * M2.data.q[m + l];
        }
        result(i, l) += data.q[k + i] * M2.data.q[k + l];
      }
      for (; l < M2.get_row_count(); h += ++l) {
        for (int j = 0; j < i; ++j) {
          result(j, l) += data.q[k + j] * M2.data.q[h + i];
          result(i, l) += data.q[k + j] * M2.data.q[h + j];
        }
        result(i, l) += data.q[k + i] * M2.data.q[h + i];
      }
    }
    return result;
  }

  /// Multiplication by a column-vector.
  /// \param M the matrix.
  /// \param V the column vector.
  /// \return column-vector equal to this * V.
  /// \throw std::range_error if this matrix and the vector dimensions don't match.
  template <typename Vector, typename Result>
  void multiply_with_vector_rhs(const Vector& V, Result& R) const {
    int k = 0;
    int i = 0;
    for (; i < R.size(); k += ++i) {
      for (int j = 0; j < i; ++j) {
        R[i] += data.q[k + j] * V[j];
        R[j] += data.q[k + j] * V[i];
      }
      R[i] += data.q[k + i] * V[i];
    }
  }
  template <typename Vector, typename Result>
  void multiply_with_vector_lhs(const Vector& V, Result& R) const {
    multiply_with_vector_rhs(V, R);
  }

  /*******************************************************************************
                           Special Methods
  *******************************************************************************/
  /// Extracts a sub-matrix from this matrix.
  /// \param M A symmetric matrix.
  /// \param aRowOffset Number of rows before the start of the sub-matrix rows.
  /// \param aColOffset Number of columns before the start of the sub-matrix columns.
  /// \param aRowCountOut Number of rows of the sub-matrix.
  /// \param aColCountOut Number of columns of the sub-matrix.
  /// \return The sub-matrix contained in this matrix.
  /// \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
  friend mat<value_type, mat_structure::rectangular, Alignment> get_block(
      const self& M, int aRowOffset, int aColOffset, int aRowCountOut,
      int aColCountOut) {
    if ((aRowOffset + aRowCountOut > M.get_row_count()) ||
        (aColOffset + aColCountOut > M.get_row_count())) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    mat<value_type, mat_structure::rectangular, Alignment> result(aRowCountOut,
                                                                  aColCountOut);
    int k = mat_triangular_size(aColOffset);
    for (int j = 0; j < aColCountOut; k += (++j + aColOffset)) {
      int h = mat_triangular_size(aRowOffset);
      int i = 0;
      for (; ((i < aRowCountOut) && (i + aRowOffset <= j + aColOffset));
           h += (++i + aRowOffset)) {
        result(i, j) = M.data.q[k + i + aRowOffset];
      }
      for (; i < aRowCountOut; h += (++i + aRowOffset)) {
        result(i, j) = M.data.q[h + j + aColOffset];
      }
    }
    return result;
  }

  /// Extracts a sub-matrix from this matrix.
  /// \param M A symmetric matrix.
  /// \param aRowOffset Number of rows before the start of the sub-matrix rows.
  /// \param aColOffset Number of columns before the start of the sub-matrix columns.
  /// \param aSizeOut Number of rows and columns of the sub-matrix.
  /// \return The sub-matrix contained in this matrix.
  /// \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
  friend mat<value_type, mat_structure::square, Alignment> get_block(
      const self& M, int aRowOffset, int aColOffset, int aSizeOut) {
    if ((aRowOffset + aSizeOut > M.get_row_count()) ||
        (aColOffset + aSizeOut > M.get_row_count())) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    mat<value_type, mat_structure::square, Alignment> result(aSizeOut);
    int k = mat_triangular_size(aColOffset);
    for (int j = 0; j < aSizeOut; k += (++j + aColOffset)) {
      int h = mat_triangular_size(aRowOffset);
      int i = 0;
      for (; ((i < aSizeOut) && (i + aRowOffset <= j + aColOffset));
           h += (++i + aRowOffset)) {
        result(i, j) = M.data.q[k + i + aRowOffset];
      }
      for (; i < aSizeOut; h += (++i + aRowOffset)) {
        result(i, j) = M.data.q[h + j + aColOffset];
      }
    }
    return result;
  }

  /// Extracts a symmetric sub-matrix from this matrix.
  /// \param M A symmetric matrix.
  /// \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
  /// \param aSizeOut Number of rows/columns of the sub-matrix.
  /// \return The symmetric sub-matrix contained in this matrix.
  /// \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
  friend mat<value_type, mat_structure::symmetric, Alignment> get_block(
      const self& M, int aDiagOffset, int aSizeOut) {
    if (aDiagOffset + aSizeOut > M.get_row_count()) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    mat<value_type, mat_structure::symmetric, Alignment> result(aSizeOut);
    int k = mat_triangular_size(aDiagOffset);
    int k_out = 0;
    for (int i = 0; i < aSizeOut; k_out += ++i, k += (i + aDiagOffset)) {
      for (int j = 0; j <= i; ++j) {
        result.data.q[k_out + j] = M.data.q[k + j + aDiagOffset];
      }
    }
    return result;
  }

  /// Sets the sub-part of this matrix to a symmetric sub-matrix M.
  /// \param M A symmetric sub-matrix that will be written in the sub-part of this matrix.
  /// \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
  /// \return This matrix, by reference.
  /// \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
  template <unsigned int SubRowCount>
  friend self& set_block(self& M,
                         const mat<value_type, mat_structure::symmetric,
                                   Alignment, SubRowCount, SubRowCount>& subM,
                         int aDiagOffset) {
    if (aDiagOffset + subM.get_row_count() > M.get_row_count()) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int k = mat_triangular_size(aDiagOffset);
    int k_in = 0;
    for (int i = 0; i < subM.get_row_count();
         k_in += ++i, k += (i + aDiagOffset)) {
      for (int j = 0; j <= i; ++j) {
        M.data.q[k + j + aDiagOffset] = subM.data.q[k_in + j];
      }
    }
    return M;
  }

  /// Appends the matrix 'rhs' to the end of the matrix 'lhs', which are both symmetric matrices.
  /// \param lhs The symmetric matrix to which to append the other.
  /// \param rhs The symmetric matrix to be appended to 'lhs'.
  template <unsigned int SubRowCount>
  friend void append_block_diag(
      self& lhs, const mat<value_type, mat_structure::symmetric, Alignment,
                           SubRowCount, SubRowCount>& rhs) {
    static_assert(is_dynamic_size);
    int oldCount = lhs.get_row_count();
    lhs.set_row_count(oldCount + rhs.get_row_count());
    set_block(lhs, rhs, oldCount);
  }

  /// Transposes the matrix M (which has no effect since M is symmetric, simply copies it).
  /// \param M The symmetric matrix to be transposed.
  /// \return The transpose of M.
  friend self transpose(const self& M) { return M; }

  /// Transposes the matrix M in a potentially destructive way (move-semantics, pre-C++0x).
  /// \param M The symmetric matrix to be transposed and moved.
  /// \return The transpose of M.
  friend self transpose_move(self& M) {
    self result;
    swap(result, M);
    return result;
  }

  /// Returns the trace of matrix M.
  /// \param M A symmetric matrix.
  /// \return the trace of matrix M.
  friend value_type trace(const self& M) {
    value_type sum = value_type(0);
    int k = 0;
    for (int i = 0; i < M.get_row_count(); k += ++i) {
      sum += M.data.q[k + i];
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

#endif  // REAK_MATH_LIN_ALG_MAT_ALG_SYMMETRIC_H_

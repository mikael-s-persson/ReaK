/**
 * \file mat_alg_skew_symmetric.h
 *
 * This library implements the specialization of the mat<> template for a
 * skew-symmetric matrix (dynamic dimension). This matrix type fulfills the matrix
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

#ifndef REAK_MATH_LIN_ALG_MAT_ALG_SKEW_SYMMETRIC_H_
#define REAK_MATH_LIN_ALG_MAT_ALG_SKEW_SYMMETRIC_H_

#include "ReaK/math/lin_alg/mat_alg_general.h"
#include "ReaK/math/lin_alg/mat_op_results.h"

#include <type_traits>

namespace ReaK {

/// This class holds a skew-symmetric matrix. This class will hold only the strict upper-triangular part
/// since the lower part is assumed to be equal to the negative of the upper one.
///
/// Models: ReadableMatrix, WritableMatrix, and ResizableMatrix.
///
/// \tparam T Arithmetic type of the elements of the matrix.
/// \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or
/// mat_alignment::column_major (default).
template <typename T, mat_alignment::tag Alignment, unsigned int RowCount>
class mat<T, mat_structure::skew_symmetric, Alignment, RowCount, RowCount>
    : public serializable {
 public:
  using self =
      mat<T, mat_structure::skew_symmetric, Alignment, RowCount, RowCount>;
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
  static constexpr mat_structure::tag structure = mat_structure::skew_symmetric;

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
  mat(int aRowCount, const value_type& aFill) {
    if constexpr (is_dynamic_size) {
      data.q.resize(mat_triangular_size(aRowCount - 1), aFill);
      data.rowCount = aRowCount;
    } else {
      if (aRowCount != RowCount) {
        throw std::range_error("Row count mismatch!");
      }
      std::fill(data.q.begin(), data.q.end(), aFill);
    }
  }
  explicit mat(int aRowCount) : mat(aRowCount, value_type(0.0)) {}

  /// Default Copy/Move Constructors.
  mat(const self& rhs) = default;
  mat(self&& rhs) = default;
  self& operator=(const self& rhs) = default;
  self& operator=(self&& rhs) = default;
  ~mat() override = default;

  /// Explicit constructor from any type of matrix. The "(M - M.transpose) / 2" is applied to guarantee skew-symmetry.
  template <ReadableMatrix Matrix>
  explicit mat(const Matrix& M)
      : mat(std::max(M.get_row_count(), M.get_col_count())) {
    int k = 0;
    int i = 1;
    const int min_size = std::min(M.get_row_count(), M.get_col_count());
    for (; i < min_size; k += i++) {
      for (int j = 0; j < i; ++j) {
        data.q[k + j] = value_type(0.5) * (M(j, i) - M(i, j));
      }
    }
    if (M.get_row_count() > M.get_col_count()) {
      for (; i < get_row_count(); k += i++) {
        for (int j = 0; j < min_size; ++j) {
          data.q[k + j] = value_type(-0.5) * M(i, j);
        }
      }
    } else {
      for (; i < get_row_count(); k += i++) {
        for (int j = 0; j < min_size; ++j) {
          data.q[k + j] = value_type(0.5) * M(j, i);
        }
      }
    }
  }

  /// Constructor from a symmetric matrix, i.e., takes the skew-symmetric part of a symmetric matrix, which is null.
  template <mat_alignment::tag Align2>
  explicit mat(const mat<value_type, mat_structure::symmetric, Align2, RowCount,
                         RowCount>& M)
      : mat(M.get_row_count()) {}

  /// Constructs a 2x2 skew-symmetric matrix from one element.
  explicit mat(const_reference a12) : mat(2) { data.q[0] = a12; }

  /// Constructs a 3x3 skew-symmetric matrix from 3 elements.
  mat(const_reference a12, const_reference a13, const_reference a23) : mat(3) {
    data.q[0] = a12;
    data.q[1] = a13;
    data.q[2] = a23;
  }

  /// Constructs a 4x4 skew-symmetric matrix from six elements.
  mat(const_reference a12, const_reference a13, const_reference a14,
      const_reference a23, const_reference a24, const_reference a34)
      : mat(4) {
    data.q[0] = a12;
    data.q[1] = a13;
    data.q[2] = a23;
    data.q[3] = a14;
    data.q[4] = a24;
    data.q[5] = a34;
  }

  /// Explicit constructor of a skew-symmetric matrix from a 3D vector (cross-product matrix).
  /// \throw std::range_error if the size of the vector is not 3.
  explicit mat(const vect_n<T>& V) : mat(3) {
    if (V.size() != 3) {
      throw std::range_error(
          "To construct a skew-matrix from a vector, that vector must have "
          "dimension 3");
    }
    data.q[0] = -V[2];
    data.q[1] = V[1];
    data.q[2] = -V[0];
  }

  /// Explicit constructor of a skew-symmetric matrix from a 3D vector (cross-product matrix).
  explicit mat(const vect<T, 3>& V) : mat(3) {
    data.q[0] = -V[2];
    data.q[1] = V[1];
    data.q[2] = -V[0];
  }

  /// The standard swap function (works with ADL).
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.data.q, rhs.data.q);
    if constexpr (is_dynamic_size) {
      swap(lhs.data.rowCount, rhs.data.rowCount);
    }
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Matrix indexing accessor for read-write access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  /// \throw std::range_error if the element accessed cannot be written to (diagonal elements).
  reference operator()(int i, int j) {
    if (i > j) {
      return data.q[mat_triangular_size(i - 1) + j];
    }
    if (i < j) {
      return data.q[mat_triangular_size(j - 1) + i];
    }
    throw std::range_error(
        "Cannot set the elements of the diagonal of a skew-symmetric "
        "matrix!");
  }

  /// Matrix indexing accessor for read-only access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  value_type operator()(int i, int j) const {
    if (i > j) {
      return -data.q[mat_triangular_size(i - 1) + j];
    }
    if (i < j) {
      return data.q[mat_triangular_size(j - 1) + i];
    }
    return value_type(0.0);
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
      data.q.resize(mat_triangular_size(aRowCount - 1), value_type(0));
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
  void set_col_count(unsigned int aColCount, bool /*unused*/ = false) {
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

  /// Standard Assignment operator with a symmetric matrix.
  template <ReadableMatrix Matrix>
  self& operator=(const Matrix& M) {
    self tmp(M);
    swap(*this, tmp);
    return *this;
  }

  /// Add-and-store operator with standard semantics.
  /// \param M the other matrix to be added to this.
  /// \return this matrix by reference.
  /// \throw std::range_error if the matrix dimensions don't match.
  self& operator+=(const self& M) {
    if constexpr (is_dynamic_size) {
      if (M.get_row_count() != get_row_count()) {
        throw std::range_error("Matrix dimension mismatch.");
      }
    }
    auto it = data.q.begin();
    for (auto cit = M.data.q.begin(); it != data.q.end(); ++it, ++cit) {
      *it += *cit;
    }
    return *this;
  }

  /// Sub-and-store operator with standard semantics.
  /// \param M the other matrix to be substracted from this.
  /// \return this matrix by reference.
  /// \throw std::range_error if the matrix dimensions don't match.
  self& operator-=(const self& M) {
    if constexpr (is_dynamic_size) {
      if (M.get_row_count() != get_row_count()) {
        throw std::range_error("Matrix dimension mismatch.");
      }
    }
    auto it = data.q.begin();
    for (auto cit = M.data.q.begin(); it != data.q.end(); ++it, ++cit) {
      *it -= *cit;
    }
    return *this;
  }

  /// Scalar-multiply-and-store operator with standard semantics.
  /// \param S the scalar to be multiplied to this.
  /// \return this matrix by reference.
  self& operator*=(const T& S) {
    for (auto& x : data.q) {
      x *= S;
    }
    return *this;
  }

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /// Addition operator with standard semantics.
  /// \param M the other matrix to be added to this.
  /// \return the matrix sum.
  /// \throw std::range_error if the matrix dimensions don't match.
  friend self operator+(self M1, const self& M2) {
    M1 += M2;
    return M1;
  }

  /// Negation operator with standard semantics.
  /// \return the negative of this matrix.
  self operator-() const {
    self result(get_row_count());
    auto it = result.data.q.begin();
    for (auto cit = data.q.begin(); cit != data.q.end(); ++it, ++cit) {
      *it -= *cit;
    }
    return result;
  }

  /// Substraction operator with standard semantics.
  /// \param M the other matrix to be substracted from this.
  /// \return the matrix difference.
  /// \throw std::range_error if the matrix dimensions don't match.
  friend self operator-(self M1, const self& M2) {
    M1 -= M2;
    return M1;
  }

  /// Multiplication operator with standard semantics.
  /// \param M1 the first matrix (the skew-symmetric one).
  /// \param M2 the other matrix.
  /// \return the matrix multiplication result, this * M.
  /// \throw std::range_error if the matrix dimensions don't match.
  template <ReadableMatrix Matrix>
  auto multiply_this_and_dense_mat(const Matrix& M2) const {
    mat_product_result_t<self, Matrix> result{M2.get_row_count(),
                                              M2.get_col_count()};
    int k = 0;
    int i = 1;
    for (; i < get_row_count(); k += i++) {
      for (int l = 0; l < M2.get_col_count(); ++l) {
        for (int j = 0; j < i; ++j) {
          result(j, l) += data.q[k + j] * M2(i, l);
          result(i, l) -= data.q[k + j] * M2(j, l);
        }
      }
    }
    return result;  // NRVO
  }

  /// Multiplication operator with standard semantics.
  /// \param M1 the first matrix.
  /// \param M2 the other matrix (the skew-symmetric one).
  /// \return the matrix multiplication result.
  /// \throw std::range_error if the matrix dimensions don't match.
  template <ReadableMatrix Matrix>
  auto multiply_dense_and_this_mat(const Matrix& M1) const {
    mat_product_result_t<Matrix, self> result{M1.get_row_count(),
                                              M1.get_col_count()};
    int k = 0;
    int i = 1;
    for (; i < get_row_count(); k += i++) {
      for (int l = 0; l < M1.get_row_count(); ++l) {
        for (int j = 0; j < i; ++j) {
          result(l, j) -= data.q[k + j] * M1(l, i);
          result(l, i) += data.q[k + j] * M1(l, j);
        }
      }
    }
    return result;
  }

  /// Multiplication operator with a skew-symmetric matrix.
  /// \param M the other skew-symmetric matrix to be multiplied by this.
  /// \return the matrix multiplication result, this * M.
  /// \throw std::range_error if the matrix dimensions don't match.
  auto multiply_with_same_mat(const self& M2) const {
    mat_product_result_t<self, self> result{get_row_count()};
    int k = 0;
    int i = 1;
    for (; i < get_row_count(); k += i++) {
      int h = 0;
      int l = 0;
      for (; l <= i; h += l++) {
        int m = 0;
        int j = 0;
        for (; j < l; m += j++) {
          result(j, l) -= data.q[k + j] * M2.data.q[k + l];
          result(i, l) -= data.q[k + j] * M2.data.q[h + j];
        }
        result(j, l) -= data.q[k + j] * M2.data.q[k + l];
        for (m += j++; j < i; m += j++) {
          result(j, l) -= data.q[k + j] * M2.data.q[k + l];
          result(i, l) += data.q[k + j] * M2.data.q[m + l];
        }
      }
      for (; l < M2.get_row_count(); h += l++) {
        int j = 0;
        for (; j < i; j++) {
          result(j, l) += data.q[k + j] * M2.data.q[h + i];
          result(i, l) -= data.q[k + j] * M2.data.q[h + j];
        }
      }
    }
    return result;
  }

  /// Multiplication by a column-vector.
  /// \param M the matrix.
  /// \param V the column vector.
  /// \return column-vector equal to this * V.
  /// \throw std::range_error if this matrix and the vector dimensions don't match.
  template <ReadableVector Vector, WritableVector Result>
  void multiply_with_vector_rhs(const Vector& V, Result& R) const {
    int k = 0;
    int i = 1;
    for (; i < R.size(); k += i++) {
      for (int j = 0; j < i; ++j) {
        R[j] += data.q[k + j] * V[i];
        R[i] -= data.q[k + j] * V[j];
      }
    }
  }

  /// Multiplication by a row-vector.
  /// \param M the matrix.
  /// \param V the row vector.
  /// \return column-vector equal to this * V.
  /// \throw std::range_error if this matrix and the vector dimensions don't match.
  template <ReadableVector Vector, WritableVector Result>
  void multiply_with_vector_lhs(const Vector& V, Result& R) const {
    int k = 0;
    int i = 1;
    for (; i < R.size(); k += i++) {
      for (int j = 0; j < i; ++j) {
        R[j] -= data.q[k + j] * V[i];
        R[i] += data.q[k + j] * V[j];
      }
    }
  }

  /*******************************************************************************
                           Special Methods
  *******************************************************************************/

  /// Extracts a skew-symmetric sub-matrix from this matrix.
  /// \param M The skew-symmetric matrix from which the sub-block is obtained.
  /// \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
  /// \param aSizeOut Number of rows/columns of the sub-matrix.
  /// \return The skew-symmetric sub-matrix contained in this matrix.
  /// \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
  friend mat<value_type, mat_structure::skew_symmetric, Alignment> get_block(
      const self& M, int aDiagOffset, int aSizeOut) {
    if (aDiagOffset + aSizeOut > M.get_row_count()) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    mat<value_type, mat_structure::skew_symmetric, Alignment> result(aSizeOut);
    int k = mat_triangular_size(aDiagOffset);
    int k_out = 0;
    for (int i = 1; i < aSizeOut; k += (i + aDiagOffset), k_out += i++) {
      for (int j = 0; j < i; ++j) {
        result.data.q[k_out + j] = M.data.q[k + j + aDiagOffset];
      }
    }
    return result;
  }

  /// Sets the sub-part of this matrix to a skew-symmetric sub-matrix M.
  /// \param M A skew-symmetric sub-matrix that will be written in the sub-part of this matrix.
  /// \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
  /// \return This matrix, by reference.
  /// \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
  template <unsigned int SubRowCount>
  friend self& set_block(self& M,
                         const mat<value_type, mat_structure::skew_symmetric,
                                   Alignment, SubRowCount, SubRowCount>& subM,
                         int aDiagOffset) {
    if (aDiagOffset + subM.get_row_count() > M.get_row_count()) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int k = mat_triangular_size(aDiagOffset);
    int k_in = 0;
    for (int i = 1; i < subM.get_row_count();
         k += (i + aDiagOffset), k_in += i++) {
      for (int j = 0; j <= i; ++j) {
        M.data.q[k + j + aDiagOffset] = subM.data.q[k_in + j];
      }
    }
    return M;
  }

  /// Appends the matrix 'rhs' to the end of the matrix 'lhs'.
  /// \param lhs The matrix to which to append the other.
  /// \param rhs The matrix to be appended to 'lhs'.
  template <unsigned int SubRowCount>
  friend void append_block_diag(
      self& lhs, const mat<value_type, mat_structure::skew_symmetric, Alignment,
                           SubRowCount, SubRowCount>& rhs) {
    static_assert(is_dynamic_size);
    int oldCount = lhs.get_row_count();
    lhs.set_row_count(oldCount + rhs.get_row_count(), true);
    set_block(lhs, rhs, oldCount);
  }

  /// Transposes the matrix M.
  /// \param M The matrix to be transposed.
  /// \return The transpose of M.
  friend self transpose(const self& M) { return -M; }

  /// Transposes and moves the matrix M.
  /// \param M The matrix to be transposed and moved.
  /// \return The transpose of M.
  friend self transpose(self&& M) {
    self result(std::move(M));
    for (auto& x : result.data.q) {
      x = -x;
    }
    return result;
  }

  /// Returns the trace of matrix M.
  /// \param M A diagonal matrix.
  /// \return the trace of matrix M.
  friend value_type trace(const self& M) { return value_type(0); }

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

#define RK_CREATE_SUBSKEWMATRIX_TRANSPOSE_OPERATORS(SUBMATRIX)            \
  template <typename Matrix>                                              \
  mat<mat_value_type_t<Matrix>, mat_structure::skew_symmetric> transpose( \
      const SUBMATRIX<Matrix, mat_structure::skew_symmetric>& M) {        \
    return -M;                                                            \
  }

RK_CREATE_SUBSKEWMATRIX_TRANSPOSE_OPERATORS(mat_copy_sub_sym_block)
RK_CREATE_SUBSKEWMATRIX_TRANSPOSE_OPERATORS(mat_sub_sym_block)
RK_CREATE_SUBSKEWMATRIX_TRANSPOSE_OPERATORS(mat_const_sub_sym_block)

#undef RK_CREATE_SUBSKEWMATRIX_TRANSPOSE_OPERATORS

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_ALG_SKEW_SYMMETRIC_H_

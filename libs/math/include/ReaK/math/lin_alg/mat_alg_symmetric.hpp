/**
 * \file mat_alg_symmetric.hpp
 *
 * This library implements the specialization of the mat<> template for a
 * general symmetric matrix (dynamic dimension). This matrix type fulfills the matrix
 * concepts of Readable, Writable, Resizable and DynAlloc, but not FullyWritable.
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

#ifndef REAK_MAT_ALG_SYMMETRIC_HPP
#define REAK_MAT_ALG_SYMMETRIC_HPP

#include "ReaK/math/lin_alg/mat_alg_general.hpp"

#include <type_traits>

namespace ReaK {

//....not good, the following: (ill-defined)

template <mat_alignment::tag Alignment>
struct mat_indexer<mat_structure::symmetric, Alignment> {
  int rowCount;
  explicit mat_indexer<mat_structure::symmetric, Alignment>(int aRowCount)
      : rowCount(aRowCount) {}
  int mat_triangular_size(int Size) { return (Size * (Size - 1)) / 2 + Size; }
  int operator()(int i, int j) const {
    if (i > j) {
      return mat_triangular_size(i) + j;
    }
    return mat_triangular_size(j) + i;
  }
};

/**
 * This class holds a symmetric matrix. This class will hold only the upper-triangular part
 * since the lower part is assumed to be equal to the upper one.
 *
 * Models: ReadableMatrixConcept, WritableMatrixConcept, ResizableMatrixConcept, and DynAllocMatrixConcept.
 *
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or
 *mat_alignment::column_major (default).
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T, mat_alignment::tag Alignment, typename Allocator>
class mat<T, mat_structure::symmetric, Alignment, Allocator>
    : public serializable {
 public:
  using self = mat<T, mat_structure::symmetric, Alignment, Allocator>;
  using allocator_type = Allocator;

  using value_type = T;
  using container_type = std::vector<value_type, allocator_type>;

  using reference = typename container_type::reference;
  using const_reference = typename container_type::const_reference;
  using pointer = typename container_type::pointer;
  using const_pointer = typename container_type::const_pointer;

  using col_iterator = void;
  using const_col_iterator = void;
  using row_iterator = void;
  using const_row_iterator = void;

  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  static constexpr std::size_t static_row_count = 0;
  static constexpr std::size_t static_col_count = 0;
  static constexpr mat_alignment::tag alignment = Alignment;
  static constexpr mat_structure::tag structure = mat_structure::symmetric;

 private:
  /// Holds the array of scalar entries.
  container_type q;
  /// Holds the dimension, both row and column count are equal to size.
  size_type rowCount;

  static size_type mat_triangular_size(size_type Size) {
    return (Size * (Size - 1)) / 2 + Size;
  }

 public:
  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/
  /**
   * Default constructor: sets all to zero.
   * \test PASSED
   */
  explicit mat(const allocator_type& aAlloc = allocator_type())
      : q(0, value_type(), aAlloc), rowCount(0) {}

  /**
   * Constructor for a sized matrix.
   * \test PASSED
   */
  explicit mat(size_type aRowCount, T aFill = 0,
               const allocator_type& aAlloc = allocator_type())
      : q(mat_triangular_size(aRowCount), aFill, aAlloc), rowCount(aRowCount) {}

  /**
   * Constructor for an identity matrix.
   * \test PASSED
   */
  mat(size_type aRowCount, bool aIdentity,
      const allocator_type& aAlloc = allocator_type())
      : q(mat_triangular_size(aRowCount), T(0.0), aAlloc), rowCount(aRowCount) {
    if (aIdentity) {
      int k = 0;
      for (int i = 0; i < rowCount; k += ++i) {
        q[k + i] = 1.0;
      }
    }
  }

  /**
   * Standard Copy Constructor with standard semantics.
   * \test PASSED
   */
  mat(const self& M) : q(M.q), rowCount(M.rowCount) {}

  /**
   * Standard Copy Constructor with standard semantics.
   * \test PASSED
   */
  mat(self&& M) noexcept : q(std::move(M.q)), rowCount(std::move(M.rowCount)) {}

  /**
   * The standard swap function (works with ADL).
   */
  friend void swap(self& M1, self& M2) noexcept {
    using std::swap;
    swap(M1.q, M2.q);
    swap(M1.rowCount, M2.rowCount);
  }

  /**
   * Explicit constructor from any type of matrix. The "(M + M.transpose) / 2" is applied to guarantee symmetry.
   * \test PASSED
   */
  template <typename Matrix>
  explicit mat(const Matrix& M,
               std::enable_if_t<is_readable_matrix_v<Matrix> &&
                                    !is_symmetric_matrix_v<Matrix> &&
                                    !std::is_same_v<Matrix, self>,
                                void*>
                   dummy = nullptr)
      : q(mat_triangular_size((M.get_row_count() > M.get_col_count()
                                   ? M.get_row_count()
                                   : M.get_col_count())),
          T(0.0)),
        rowCount((M.get_row_count() > M.get_col_count() ? M.get_row_count()
                                                        : M.get_col_count())) {
    int k = 0;
    int i = 0;
    int min_size = (M.get_row_count() > M.get_col_count() ? M.get_col_count()
                                                          : M.get_row_count());
    for (; i < min_size; k += ++i) {
      for (int j = 0; j < i; ++j) {
        q[k + j] = value_type(0.5) * (M(j, i) + M(i, j));
      }
      q[k + i] = M(i, i);
    }
    if (M.get_row_count() > M.get_col_count()) {
      for (; i < rowCount; k += ++i) {
        for (int j = 0; j < min_size; ++j) {
          q[k + j] = value_type(0.5) * M(i, j);
        }
      }
    } else {
      for (; i < rowCount; k += ++i) {
        for (int j = 0; j < min_size; ++j) {
          q[k + j] = value_type(0.5) * M(j, i);
        }
      }
    }
  }

  /**
   * Explicit constructor from any type of matrix. The "(M + M.transpose) / 2" is applied to guarantee symmetry.
   * \test PASSED
   */
  template <typename Matrix>
  explicit mat(const Matrix& M,
               std::enable_if_t<is_readable_matrix_v<Matrix> &&
                                    is_symmetric_matrix_v<Matrix> &&
                                    !std::is_same_v<Matrix, self>,
                                void*>
                   dummy = nullptr)
      : q(mat_triangular_size(M.get_row_count()), T(0.0)),
        rowCount(M.get_row_count()) {
    int k = 0;
    int i = 0;
    for (; i < rowCount; k += ++i) {
      for (int j = 0; j < i; ++j) {
        q[k + j] = M(i, j);
      }
      q[k + i] = M(i, i);
    }
  }

  /**
   * Destructor.
   * \test PASSED
   */
  ~mat() override = default;

  /**
   * Constructs a 2x2 symmetric matrix from three elements.
   * \test PASSED
   */
  mat(const_reference a11, const_reference a12, const_reference a22)
      : q(3), rowCount(2) {
    q[0] = a11;
    q[1] = a12;
    q[2] = a22;
  }

  /**
   * Constructs a 3x3 symmetric matrix from six elements.
   * \test PASSED
   */
  mat(const_reference a11, const_reference a12, const_reference a13,
      const_reference a22, const_reference a23, const_reference a33)
      : q(6), rowCount(3) {
    q[0] = a11;
    q[1] = a12;
    q[2] = a22;
    q[3] = a13;
    q[4] = a23;
    q[5] = a33;
  }

  /**
   * Constructs a 4x4 symmetric matrix from ten elements.
   * \test PASSED
   */
  mat(const_reference a11, const_reference a12, const_reference a13,
      const_reference a14, const_reference a22, const_reference a23,
      const_reference a24, const_reference a33, const_reference a34,
      const_reference a44)
      : q(10), rowCount(4) {
    q[0] = a11;
    q[1] = a12;
    q[2] = a22;
    q[3] = a13;
    q[4] = a23;
    q[5] = a33;
    q[6] = a14;
    q[7] = a24;
    q[8] = a34;
    q[9] = a44;
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-write access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  reference operator()(int i, int j) {
    if (i > j) {
      return q[mat_triangular_size(i) + j];
    }
    return q[mat_triangular_size(j) + i];
  }

  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  const_reference operator()(int i, int j) const {
    if (i > j) {
      return q[mat_triangular_size(i) + j];
    }
    return q[mat_triangular_size(j) + i];
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
    q.resize(mat_triangular_size(aRowCount), T(0.0));
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
    q.resize(mat_triangular_size(aColCount), T(0.0));
    rowCount = aColCount;
  }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<size_type, size_type> size() const noexcept {
    return std::make_pair(rowCount, rowCount);
  }
  /**
   * Sets the row-count and column-count of the matrix, via a std::pair of dimension values.
   * \param sz new dimensions for the matrix.
   * \test PASSED
   */
  void resize(const std::pair<size_type, size_type>& sz) {
    set_row_count(sz.first, true);
  }

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return q.get_allocator(); }

  /*******************************************************************************
                           Assignment Operators
  *******************************************************************************/

  /**
   * Standard Assignment operator with a symmetric matrix.
   * \test PASSED
   */
  self& operator=(self M) {
    swap(*this, M);
    return *this;
  }

  /**
   * Standard Assignment operator with a matrix of any type. The "(M + M.transpose) / 2" formula is applied to guarantee
   * symmetry.
   * \test PASSED
   */
  template <typename Matrix>
  self& operator=(const Matrix& M) {
    self tmp(M);
    swap(*this, tmp);
    return *this;
  }

  /**
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template <mat_alignment::tag Align2, typename Alloc2>
  self& operator+=(const mat<T, mat_structure::symmetric, Align2, Alloc2>& M) {
    if (M.get_row_count() != rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int k = 0;
    for (int i = 0; i < rowCount; k += ++i) {
      for (int j = 0; j <= i; ++j) {
        q[k + j] += M(i, j);
      }
    }
    return *this;
  }

  /**
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template <mat_alignment::tag Align2, typename Alloc2>
  self& operator+=(const mat<T, mat_structure::diagonal, Align2, Alloc2>& M) {
    if (M.get_row_count() != rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int k = 0;
    for (int i = 0; i < rowCount; k += ++i) {
      q[k + i] += M(i, i);
    }
    return *this;
  }

  /**
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template <mat_alignment::tag Align2, typename Alloc2>
  self& operator-=(const mat<T, mat_structure::symmetric, Align2, Alloc2>& M) {
    if (M.get_row_count() != rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int k = 0;
    for (int i = 0; i < rowCount; k += ++i) {
      for (int j = 0; j <= i; ++j) {
        q[k + j] -= M(i, j);
      }
    }
    return *this;
  }

  /**
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template <mat_alignment::tag Align2, typename Alloc2>
  self& operator-=(const mat<T, mat_structure::diagonal, Align2, Alloc2>& M) {
    if (M.get_row_count() != rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int k = 0;
    for (int i = 0; i < rowCount; k += ++i) {
      q[k + i] -= M(i, i);
    }
    return *this;
  }

  /**
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=(const T& S) {
    for (auto& v : q) {
      v *= S;
    }
    return *this;
  }

  /**
   * Negation operator for any type of matrices. This is a default operator
   * that will be called if no better special-purpose overload exists.
   * \return Symmetric matrix.
   * \test PASSED
   */
  self operator-() const {
    self result(*this);
    auto itr = result.q.begin();
    for (auto it = q.begin(); it != q.end(); ++it, ++itr) {
      *itr = -(*it);
    }
    return result;
  }

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /**
   * Add two matrices.
   * \param M1 the first matrix (first operand).
   * \param M2 the other matrix (second operand).
   * \return symmetric matrix sum of M1 and M2.
   * \throw std::range_error if the two matrix dimensions don't match.
   * \test PASSED
   */
  friend self operator+(const self& M1, const self& M2) {
    if (M1.rowCount != M2.rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    self result(M1.rowCount);
    int k = 0;
    for (int i = 0; i < M1.rowCount; k += ++i) {
      for (int j = 0; j <= i; ++j) {
        result.q[k + j] = M1.q[k + j] + M2.q[k + j];
      }
    }
    return result;
  }

  /**
   * Sub two matrices.
   * \param M1 the first matrix (first operand).
   * \param M2 the other matrix (second operand).
   * \return symmetric matrix difference of M1 and M2.
   * \throw std::range_error if the two matrix dimensions don't match.
   * \test PASSED
   */
  friend self operator-(const self& M1, const self& M2) {
    if (M1.rowCount != M2.rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    self result(M1.rowCount);
    int k = 0;
    for (int i = 0; i < M1.rowCount; k += ++i) {
      for (int j = 0; j <= i; ++j) {
        result.q[k + j] = M1.q[k + j] - M2.q[k + j];
      }
    }
    return result;
  }

  /**
   * General Matrix multiplication.
   * \param M the other matrix (second operand).
   * \return general matrix multiplication result of this and M.
   * \throw std::range_error if the two matrix dimensions don't match.
   */
  template <typename Matrix>
  auto multiply_this_and_dense_mat(const Matrix& M2) const {
    using ValueType = mat_value_type_t<Matrix>;
    mat_product_result_t<self, Matrix> result{M2.get_row_count(),
                                              M2.get_col_count(), ValueType(0),
                                              M2.get_allocator()};
    for (int k = 0, i = 0; i < rowCount; k += ++i) {
      for (int l = 0; l < M2.get_col_count(); ++l) {
        for (int j = 0; j < i; ++j) {
          result(j, l) += q[k + j] * M2(i, l);
          result(i, l) += q[k + j] * M2(j, l);
        }
        result(i, l) += q[k + i] * M2(i, l);
      }
    }
    return result;  // NRVO
  }

  /**
   * General Matrix multiplication.
   * \param M the other matrix (second operand).
   * \return general matrix multiplication result of this and M.
   * \throw std::range_error if the two matrix dimensions don't match.
   */
  template <typename Matrix>
  auto multiply_dense_and_this_mat(const Matrix& M1) const {
    using ValueType = mat_value_type_t<Matrix>;
    mat_product_result_t<Matrix, self> result{M1.get_row_count(),
                                              M1.get_col_count(), ValueType(0),
                                              M1.get_allocator()};
    for (int k = 0, i = 0; i < M1.get_col_count(); k += ++i) {
      for (int l = 0; l < M1.get_row_count(); ++l) {
        for (int j = 0; j < i; ++j) {
          result(l, i) += q[k + j] * M1(l, j);
          result(l, j) += q[k + j] * M1(l, i);
        }
        result(l, i) += q[k + i] * M1(l, i);
      }
    }
    return result;  // NRVO
  }

  /**
   * Symmetric Matrix multiplication.
   * \param M1 the first symmetric matrix (first operand).
   * \param M2 the other symmetric matrix (second operand).
   * \return square matrix, result of M1 times M2.
   * \throw std::range_error if the two matrix dimensions don't match.
   */
  auto multiply_with_same_mat(const self& M2) const {
    mat_product_result_t<self, self> result{rowCount, T(0), get_allocator()};
    int k = 0;
    int i = 0;
    for (; i < rowCount; k += ++i) {
      int h = 0;
      int l = 0;
      for (; l <= i; h += ++l) {
        int m = 0;
        int j = 0;
        for (; j < l; m += ++j) {
          result(j, l) += q[k + j] * M2.q[k + l];
          result(i, l) += q[k + j] * M2.q[h + j];
        }
        for (; j < i; m += ++j) {
          result(j, l) += q[k + j] * M2.q[k + l];
          result(i, l) += q[k + j] * M2.q[m + l];
        }
        result(i, l) += q[k + i] * M2.q[k + l];
      }
      for (; l < M2.rowCount; h += ++l) {
        for (int j = 0; j < i; ++j) {
          result(j, l) += q[k + j] * M2.q[h + i];
          result(i, l) += q[k + j] * M2.q[h + j];
        }
        result(i, l) += q[k + i] * M2.q[h + i];
      }
    }
    return result;
  }

  /**
   * Multiplication by a column-vector.
   * \param M the matrix.
   * \param V the column vector.
   * \return column-vector equal to this * V.
   * \throw std::range_error if this matrix and the vector dimensions don't match.
   */
  template <typename Vector, typename Result>
  void multiply_with_vector_rhs(const Vector& V, Result& R) const {
    int k = 0;
    int i = 0;
    for (; i < R.size(); k += ++i) {
      for (int j = 0; j < i; ++j) {
        R[i] += q[k + j] * V[j];
        R[j] += q[k + j] * V[i];
      }
      R[i] += q[k + i] * V[i];
    }
  }
  template <typename Vector, typename Result>
  void multiply_with_vector_lhs(const Vector& V, Result& R) const {
    multiply_with_vector_rhs(V, R);
  }

  /*******************************************************************************
                           Special Methods
  *******************************************************************************/
  /**
   * Extracts a sub-matrix from this matrix.
   * \param M A symmetric matrix.
   * \param aRowOffset Number of rows before the start of the sub-matrix rows.
   * \param aColOffset Number of columns before the start of the sub-matrix columns.
   * \param aRowCountOut Number of rows of the sub-matrix.
   * \param aColCountOut Number of columns of the sub-matrix.
   * \return The sub-matrix contained in this matrix.
   * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
   */
  friend mat<T, mat_structure::rectangular, Alignment, Allocator> get_block(
      const self& M, int aRowOffset, int aColOffset, int aRowCountOut,
      int aColCountOut) {
    if ((aRowOffset + aRowCountOut > M.rowCount) ||
        (aColOffset + aColCountOut > M.rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    mat<T, mat_structure::rectangular, Alignment, Allocator> result(
        aRowCountOut, aColCountOut, T(0), M.get_allocator());
    int k = mat_triangular_size(aColOffset);
    for (int j = 0; j < aColCountOut; k += (++j + aColOffset)) {
      int h = mat_triangular_size(aRowOffset);
      int i = 0;
      for (; ((i < aRowCountOut) && (i + aRowOffset <= j + aColOffset));
           h += (++i + aRowOffset)) {
        result(i, j) = M.q[k + i + aRowOffset];
      }
      for (; i < aRowCountOut; h += (++i + aRowOffset)) {
        result(i, j) = M.q[h + j + aColOffset];
      }
    }
    return result;
  }

  /**
   * Extracts a sub-matrix from this matrix.
   * \param M A symmetric matrix.
   * \param aRowOffset Number of rows before the start of the sub-matrix rows.
   * \param aColOffset Number of columns before the start of the sub-matrix columns.
   * \param aSizeOut Number of rows and columns of the sub-matrix.
   * \return The sub-matrix contained in this matrix.
   * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
   */
  friend mat<T, mat_structure::square, Alignment, Allocator> get_block(
      const self& M, int aRowOffset, int aColOffset, int aSizeOut) {
    if ((aRowOffset + aSizeOut > M.rowCount) ||
        (aColOffset + aSizeOut > M.rowCount)) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    mat<T, mat_structure::square, Alignment, Allocator> result(
        aSizeOut, T(0), M.get_allocator());
    int k = mat_triangular_size(aColOffset);
    for (int j = 0; j < aSizeOut; k += (++j + aColOffset)) {
      int h = mat_triangular_size(aRowOffset);
      int i = 0;
      for (; ((i < aSizeOut) && (i + aRowOffset <= j + aColOffset));
           h += (++i + aRowOffset)) {
        result(i, j) = M.q[k + i + aRowOffset];
      }
      for (; i < aSizeOut; h += (++i + aRowOffset)) {
        result(i, j) = M.q[h + j + aColOffset];
      }
    }
    return result;
  }

  /**
   * Extracts a symmetric sub-matrix from this matrix.
   * \param M A symmetric matrix.
   * \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
   * \param aSizeOut Number of rows/columns of the sub-matrix.
   * \return The symmetric sub-matrix contained in this matrix.
   * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
   */
  friend self get_block(const self& M, int aDiagOffset, int aSizeOut) {
    if (aDiagOffset + aSizeOut > M.rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    self result(aSizeOut);
    int k = mat_triangular_size(aDiagOffset);
    int k_out = 0;
    for (int i = 0; i < aSizeOut; k_out += ++i, k += (i + aDiagOffset)) {
      for (int j = 0; j <= i; ++j) {
        result.q[k_out + j] = M.q[k + j + aDiagOffset];
      }
    }
    return result;
  }

  /** Sets the sub-part of this matrix to a symmetric sub-matrix M.
   * \param M A symmetric sub-matrix that will be written in the sub-part of this matrix.
   * \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
   * \return This matrix, by reference.
   * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
   */
  friend self& set_block(self& M, const self& subM, int aDiagOffset) {
    if (aDiagOffset + subM.rowCount > M.rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    int k = mat_triangular_size(aDiagOffset);
    int k_in = 0;
    for (int i = 0; i < subM.rowCount; k_in += ++i, k += (i + aDiagOffset)) {
      for (int j = 0; j <= i; ++j) {
        M.q[k + j + aDiagOffset] = subM.q[k_in + j];
      }
    }
    return M;
  }

  /**
   * Appends the matrix 'rhs' to the end of the matrix 'lhs', which are both symmetric matrices.
   * \param lhs The symmetric matrix to which to append the other.
   * \param rhs The symmetric matrix to be appended to 'lhs'.
   */
  friend void append_block_diag(self& lhs, const self& rhs) {
    size_type oldCount = lhs.get_col_count();
    lhs.set_col_count(oldCount + rhs.get_col_count(), true);
    set_block(lhs, rhs, oldCount);
  }

  /**
   * Transposes the matrix M (which has no effect since M is symmetric, simply copies it).
   * \param M The symmetric matrix to be transposed.
   * \return The transpose of M.
   */
  friend self transpose(const self& M) { return M; }

  /**
   * Transposes the matrix M in a potentially destructive way (move-semantics, pre-C++0x).
   * \param M The symmetric matrix to be transposed and moved.
   * \return The transpose of M.
   */
  friend self transpose_move(self& M) {
    self result;
    swap(result, M);
    return result;
  }

  /**
   * Returns the trace of matrix M.
   * \param M A symmetric matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace(const self& M) {
    value_type sum = value_type(0);
    int k = 0;
    for (int i = 0; i < M.rowCount; k += ++i) {
      sum += M.q[k + i];
    }
    return sum;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& std::pair<std::string, const std::vector<T>&>("q", q) &
        std::pair<std::string, std::size_t>("rowCount", rowCount);
  }
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& std::pair<std::string, std::vector<T>&>("q", q) &
        std::pair<std::string, std::size_t&>("rowCount", rowCount);
  }

  RK_RTTI_REGISTER_CLASS_1BASE(self, 1, serializable)
};

#if 0
extern template class mat< double, mat_structure::symmetric >;
extern template class mat< float, mat_structure::symmetric >;

extern template mat< double, mat_structure::symmetric >& mat< double, mat_structure::symmetric >::
  operator+=( const mat< double, mat_structure::symmetric >& M );
extern template mat< double, mat_structure::symmetric >& mat< double, mat_structure::symmetric >::
  operator+=( const mat< double, mat_structure::diagonal >& M );
extern template mat< double, mat_structure::symmetric >& mat< double, mat_structure::symmetric >::
  operator-=( const mat< double, mat_structure::symmetric >& M );
extern template mat< double, mat_structure::symmetric >& mat< double, mat_structure::symmetric >::
  operator-=( const mat< double, mat_structure::diagonal >& M );

extern template mat< float, mat_structure::symmetric >& mat< float, mat_structure::symmetric >::
  operator+=( const mat< float, mat_structure::symmetric >& M );
extern template mat< float, mat_structure::symmetric >& mat< float, mat_structure::symmetric >::
  operator+=( const mat< float, mat_structure::diagonal >& M );
extern template mat< float, mat_structure::symmetric >& mat< float, mat_structure::symmetric >::
  operator-=( const mat< float, mat_structure::symmetric >& M );
extern template mat< float, mat_structure::symmetric >& mat< float, mat_structure::symmetric >::
  operator-=( const mat< float, mat_structure::diagonal >& M );
#endif
}  // namespace ReaK

#endif

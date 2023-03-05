/**
 * \file mat_alg_diagonal.hpp
 *
 * This library declares matrix specializations for representing and manipulating diagonal matrices.
 * This library implements many overloaded operators that turn out to be more efficiently implemented
 * if specialized for the diagonal matrix case. All those overloads are automatically selected through
 * Sfinae switches, and the diagonal matrix class is simply a partial specialization of the "ReaK::mat"
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

#ifndef REAK_MAT_ALG_DIAGONAL_HPP
#define REAK_MAT_ALG_DIAGONAL_HPP

#include "ReaK/math/lin_alg/mat_alg_general.hpp"

namespace ReaK {

template <mat_alignment::tag Alignment>
struct mat_indexer<mat_structure::diagonal, Alignment> {
  int rowCount;
  explicit mat_indexer<mat_structure::diagonal, Alignment>(int aRowCount)
      : rowCount(aRowCount) {}
  int operator()(int i, int j) const { return i; }
};

/**
 * This class holds a diagonal matrix. This class will hold only the diagonal.
 *
 * Models: ReadableMatrixConcept, WritableMatrixConcept, ResizableMatrixConcept, and DynAllocMatrixConcept.
 *
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or
 *mat_alignment::column_major (default).
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <class T, mat_alignment::tag Alignment, typename Allocator>
class mat<T, mat_structure::diagonal, Alignment, Allocator>
    : public serializable {
 public:
  using self = mat<T, mat_structure::diagonal, Alignment, Allocator>;
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
  using difference_type = typename container_type::difference_type;

  static constexpr std::size_t static_row_count = 0;
  static constexpr std::size_t static_col_count = 0;
  static constexpr mat_alignment::tag alignment = Alignment;
  static constexpr mat_structure::tag structure = mat_structure::diagonal;

 private:
  /// Holds the array of scalar entries.
  container_type q;
  /// Holds the dimension, both row and column count are equal to size.
  size_type rowCount;

 public:
  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /**
   * Constructor for a sized matrix.
   */
  explicit mat(size_type aRowCount, value_type aFill = value_type(0),
               const allocator_type& aAlloc = allocator_type())
      : q(aRowCount, aFill, aAlloc), rowCount(aRowCount) {}

  /**
   * Default constructor: sets all to zero.
   */
  explicit mat(const allocator_type& aAlloc = allocator_type())
      : mat(0, value_type(0), aAlloc) {}

  /**
   * Constructor for an identity matrix.
   */
  mat(size_type aRowCount, bool aIdentity,
      const allocator_type& aAlloc = allocator_type())
      : mat(aRowCount, (aIdentity ? value_type(1) : value_type(0)), aAlloc) {}

  /**
   * Standard Copy Constructor with standard semantics.
   */
  mat(const self& M) = default;

  /**
   * Standard Move Constructor with standard semantics.
   */
  mat(self&& M) noexcept = default;

  /**
   * Constructor from a vector of size n.
   */
  template <typename Vector>
  explicit mat(
      const Vector& V, const allocator_type& aAlloc = allocator_type(),
      std::enable_if_t<
          is_readable_vector_v<Vector> && !std::is_same_v<Vector, self>, void*>
          dummy = nullptr)
      : q(V.begin(), V.end(), aAlloc), rowCount(V.size()) {}

  /**
   * Constructor from a general matrix, copying only the diagonal part.
   */
  template <typename Matrix>
  explicit mat(
      const Matrix& M, const allocator_type& aAlloc = allocator_type(),
      std::enable_if_t<
          is_readable_matrix_v<Matrix> && !std::is_same_v<Matrix, self>, void*>
          dummy = nullptr)
      : q((M.get_row_count() < M.get_col_count() ? M.get_row_count()
                                                 : M.get_col_count()),
          T(0.0), aAlloc),
        rowCount((M.get_row_count() < M.get_col_count() ? M.get_row_count()
                                                        : M.get_col_count())) {
    for (int i = 0; i < rowCount; ++i) {
      q[i] = M(i, i);
    }
  }

  /**
   * Destructor.
   * \test PASSED
   */
  ~mat() override = default;

  /**
   * Swap friend-function that allows ADL and efficient swapping of two matrices.
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.q, rhs.q);
    swap(lhs.rowCount, rhs.rowCount);
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-write access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \throw std::range_error if the row and column index are not the same (cannot write off-diagonal terms).
   * \test PASSED
   */
  reference operator()(int i, int j) {
    if (i == j) {
      return q[i];
    }
    throw std::range_error(
        "Cannot write to the off-diagonal terms of a diagonal matrix!");
  }

  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()(int i, int j) const {
    if (i == j) {
      return q[i];
    }
    return T(0.0);
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
    q.resize(aRowCount, value_type(0.0));
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
    q.resize(aColCount, value_type(0.0));
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
   * Standard Assignment operator with a diagonal matrix.
   */
  self& operator=(self M) {
    swap(*this, M);
    return *this;
  }

  /**
   * Standard Assignment operator with a general matrix. Copying only the diagonal part of M.
   */
  template <typename Matrix>
  self& operator=(const Matrix& M) {
    self tmp(M);
    swap(*this, tmp);
    return *this;
  }

  /**
   * Add-and-store operator with standard semantics.
   * \param M the other matrix to be added to this.
   * \return this matrix by reference.
   * \throw std::range_error if the matrix dimensions don't match.
   */
  self& operator+=(const self& M) {
    if (M.rowCount != rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int i = 0; i < rowCount; ++i) {
      q[i] += M.q[i];
    }
    return *this;
  }

  /**
   * Sub-and-store operator with standard semantics.
   * \param M the other matrix to be substracted from this.
   * \return this matrix by reference.
   * \throw std::range_error if the matrix dimensions don't match.
   */
  self& operator-=(const self& M) {
    if (M.rowCount != rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int i = 0; i < rowCount; ++i) {
      q[i] -= M.q[i];
    }
    return *this;
  }

  /**
   * Scalar-multiply-and-store operator with standard semantics.
   * \param S the scalar to be multiplied to this.
   * \return this matrix by reference.
   */
  self& operator*=(const T& S) {
    for (auto& v : q) {
      v *= S;
    }
    return *this;
  }

  /**
   * Matrix-multiply-and-store operator with a diagonal matrix.
   * \param M the other matrix to be multiplied with this.
   * \return this matrix by reference.
   * \throw std::range_error if the matrix dimensions don't match.
   */
  self& operator*=(const self& M) {
    if (rowCount != M.rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int i = 0; i < rowCount; ++i) {
      q[i] *= M.q[i];
    }
    return *this;
  }

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /**
   * Addition operator with standard semantics.
   * \param M the other matrix to be added to this.
   * \return the matrix sum of this and M.
   * \throw std::range_error if the matrix dimensions don't match.
   */
  friend self operator+(const self& M1, const self& M2) {
    if (M1.rowCount != M2.rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    self result(M1.rowCount);
    for (int i = 0; i < M1.rowCount; ++i) {
      result.q[i] = M1.q[i] + M2.q[i];
    }
    return result;
  }

  /**
   * Negation operator with standard semantics.
   * \return the negative of this matrix sum.
   */
  self operator-() const {
    self result(rowCount);
    for (int i = 0; i < rowCount; ++i) {
      result.q[i] = -q[i];
    }
    return result;
  }

  /**
   * Substraction operator with standard semantics.
   * \param M the other matrix to be substracted from this.
   * \return the matrix difference of this and M.
   * \throw std::range_error if the matrix dimensions don't match.
   */
  friend self operator-(const self& M1, const self& M2) {
    if (M1.rowCount != M2.rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    self result(M1.rowCount);
    for (int i = 0; i < M1.rowCount; ++i) {
      result.q[i] = M1.q[i] - M2.q[i];
    }
    return result;
  }

  /**
   * Matrix multiplication operator with a column vector.
   * \param M a diagonal matrix.
   * \param V the vector to be multiplied with this.
   * \return the matrix-vector product of this and V.
   * \throw std::range_error if the matrix-vector dimensions don't match.
   */
  template <typename Vector, typename Result>
  void multiply_with_vector_rhs(const Vector& V, Result& R) const {
    for (int i = 0; i < R.size(); ++i) {
      R[i] = V[i] * q[i];
    }
  }

  /**
   * Matrix multiplication operator with a row vector.
   * \param M a diagonal matrix.
   * \param V the vector to be multiplied with M.
   * \return the matrix-vector product of this and V.
   * \throw std::range_error if the matrix-vector dimensions don't match.
   */
  template <typename Vector, typename Result>
  void multiply_with_vector_lhs(const Vector& V, Result& R) const {
    multiply_with_vector_rhs(V, R);
  }

  /*******************************************************************************
                           Special Methods
  *******************************************************************************/

  /**
   * Extracts a diagonal sub-matrix from this matrix.
   * \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
   * \param aSizeOut Number of rows/columns of the sub-matrix.
   * \return The diagonal sub-matrix contained in this matrix.
   * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
   */
  friend self get_block(const self& M, int aDiagOffset, int aSizeOut) {
    if (aDiagOffset + aSizeOut > M.rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    self result(aSizeOut, value_type(0), M.get_allocator());
    for (int i = 0; i < aSizeOut; ++i) {
      result.q[i] = M.q[i + aDiagOffset];
    }
    return result;
  }

  /** Sets the sub-block of a matrix to a diagnoal sub-matrix M.
   * \param M A diagonal matrix to which the sub-block will be set.
   * \param subM A diagonal sub-matrix that will be written in the sub-part of this matrix.
   * \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
   * \return This matrix, by reference.
   * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
   */
  friend self& set_block(self& M, const self& subM, int aDiagOffset) {
    if (aDiagOffset + subM.rowCount > M.rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int i = 0; i < subM.rowCount; ++i) {
      M.q[i + aDiagOffset] = subM.q[i];
    }
    return M;
  }

  /**
   * Appends the matrix 'rhs' to the end of the matrix 'lhs', which are both diagonal matrices.
   * \param lhs The diagonal matrix to which to append the other.
   * \param rhs The diagonal matrix to be appended to 'lhs'.
   */
  friend void append_block_diag(self& lhs, const self& rhs) {
    size_type oldCount = lhs.get_col_count();
    lhs.set_col_count(oldCount + rhs.get_col_count(), true);
    set_block(lhs, rhs, oldCount);
  }

  /**
   * Transposes the matrix M (which has no effect since M is diagonal, simply copies it).
   * \param M The diagonal matrix to be transposed.
   * \return The transpose of M.
   */
  friend self transpose(const self& M) { return M; }

  /**
   * Transposes the matrix M in a potentially destructive way (move-semantics, pre-C++0x).
   * \param M The diagonal matrix to be transposed and moved.
   * \return The transpose of M.
   */
  friend self transpose_move(self& M) {
    self result;
    swap(result, M);
    return result;
  }

  /**
   * Transposes the matrix M in a potentially destructive way (move-semantics, C++0x).
   * \param M The diagonal matrix to be transposed and moved.
   * \return The transpose of M.
   */
  friend self transpose(self&& M) { return self(std::move(M)); }

  /**
   * Returns the trace of matrix M.
   * \param M A diagonal matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace(const self& M) {
    value_type sum = value_type(0);
    for (auto& v : M.q) {
      sum += v;
    }
    return sum;
  }

  /**
   * Inverts this matrix.
   */
  void invert() {
    for (auto& v : q) {
      v = value_type(1.0) / v;
    }
  }

  /**
   * Inverts the matrix M.
   * \param M The diagonal matrix to be inverted.
   * \return The inverse of M.
   */
  friend self invert(self M) {
    M.invert();
    return M;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(q) & RK_SERIAL_SAVE_WITH_NAME(rowCount);
  }
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(q) & RK_SERIAL_LOAD_WITH_NAME(rowCount);
  }

  RK_RTTI_REGISTER_CLASS_1BASE(self, 1, serializable)
};

#if 0

extern template class mat< double, mat_structure::diagonal >;
extern template class mat< float, mat_structure::diagonal >;

#endif
};  // namespace ReaK

#endif

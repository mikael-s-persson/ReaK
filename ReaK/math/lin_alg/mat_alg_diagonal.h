/**
 * \file mat_alg_diagonal.h
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

#ifndef REAK_MATH_LIN_ALG_MAT_ALG_DIAGONAL_H_
#define REAK_MATH_LIN_ALG_MAT_ALG_DIAGONAL_H_

#include "ReaK/math/lin_alg/mat_alg_general.h"
#include "ReaK/math/lin_alg/vect_alg.h"

namespace ReaK {

/// This class holds a diagonal matrix. This class will hold only the diagonal.
///
/// Models: ReadableMatrixConcept, WritableMatrixConcept, and ResizableMatrixConcept.
///
/// \tparam T Arithmetic type of the elements of the matrix.
/// \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or
/// mat_alignment::column_major (default).
template <class T, mat_alignment::tag Alignment, unsigned int RowCount>
class mat<T, mat_structure::diagonal, Alignment, RowCount, RowCount>
    : public serializable {
 public:
  using self = mat<T, mat_structure::diagonal, Alignment, RowCount, RowCount>;
  static constexpr bool is_dynamic_size = (RowCount == 0);

  using value_type = T;
  using container_type = vect<value_type, RowCount>;

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
  static constexpr mat_structure::tag structure = mat_structure::diagonal;

  template <typename OtherT, mat_structure::tag OtherStructure,
            mat_alignment::tag OtherAlignment, unsigned int OtherRowCount,
            unsigned int OtherColCount>
  friend class mat;

 private:
  /// Holds the vector of scalar entries.
  container_type q;

 public:
  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /// Constructor for a sized matrix.
  explicit mat(int aRowCount, value_type aFill = value_type(0)) {
    if constexpr (is_dynamic_size) {
      q.resize(aRowCount, aFill);
    } else {
      if (aRowCount != RowCount) {
        throw std::range_error("Row count mismatch!");
      }
      std::fill(q.begin(), q.end(), aFill);
    }
  }

  /// Default constructor: sets all to zero.
  mat() : mat(0) {}

  /// Constructor for an identity matrix.
  mat(int aRowCount, bool aIdentity)
      : mat(aRowCount, (aIdentity ? value_type(1) : value_type(0))) {}

  /// Standard Copy Constructor with standard semantics.
  mat(const self& M) = default;
  mat(self&& M) noexcept = default;
  self& operator=(const self& rhs) = default;
  self& operator=(self&& rhs) = default;

  /// Constructor from a vector of size n.
  template <typename Vector>
  explicit mat(
      const Vector& V,
      std::enable_if_t<
          is_readable_vector_v<Vector> && !std::is_same_v<Vector, self>, void*>
          dummy = nullptr)
      : q(V) {}

  /// Constructor from a general matrix, copying only the diagonal part.
  template <typename Matrix>
  explicit mat(
      const Matrix& M,
      std::enable_if_t<
          is_readable_matrix_v<Matrix> && !std::is_same_v<Matrix, self>, void*>
          dummy = nullptr)
      : mat(std::min(M.get_row_count(), M.get_col_count())) {
    for (int i = 0; i < q.size(); ++i) {
      q[i] = M(i, i);
    }
  }

  ~mat() override = default;

  /// Swap friend-function that allows ADL and efficient swapping of two matrices.
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.q, rhs.q);
  }

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Matrix indexing accessor for read-write access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  /// \throw std::range_error if the row and column index are not the same (cannot write off-diagonal terms).
  reference operator()(int i, int j) {
    if (i == j) {
      return q[i];
    }
    throw std::range_error(
        "Cannot write to the off-diagonal terms of a diagonal matrix!");
  }

  /// Matrix indexing accessor for read-only access.
  /// \param i Row index.
  /// \param j Column index.
  /// \return the element at the given position.
  value_type operator()(int i, int j) const {
    if (i == j) {
      return q[i];
    }
    return T(0.0);
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
  int get_row_count() const { return q.size(); }

  /// Sets the row-count (number of rows) of the matrix.
  /// \param aRowCount new number of rows for the matrix.
  /// \param aPreserveData If true, the resizing will preserve all the data it can.
  void set_row_count(int aRowCount, bool /*unused*/ = false) {
    q.resize(aRowCount, value_type(0.0));
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

  /// Standard Assignment operator with a general matrix. Copying only the diagonal part of M.
  template <typename Matrix>
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
    q += M.q;
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
    q -= M.q;
    return *this;
  }

  /// Scalar-multiply-and-store operator with standard semantics.
  /// \param S the scalar to be multiplied to this.
  /// \return this matrix by reference.
  self& operator*=(const T& S) {
    q *= S;
    return *this;
  }

  /// Matrix-multiply-and-store operator with a diagonal matrix.
  /// \param M the other matrix to be multiplied with this.
  /// \return this matrix by reference.
  /// \throw std::range_error if the matrix dimensions don't match.
  self& operator*=(const self& M) {
    if constexpr (is_dynamic_size) {
      if (get_row_count() != M.get_row_count()) {
        throw std::range_error("Matrix dimension mismatch.");
      }
    }
    for (int i = 0; i < get_row_count(); ++i) {
      q[i] *= M.q[i];
    }
    return *this;
  }

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /// Addition operator with standard semantics.
  /// \param M the other matrix to be added to this.
  /// \return the matrix sum of this and M.
  /// \throw std::range_error if the matrix dimensions don't match.
  friend self operator+(const self& M1, const self& M2) {
    if constexpr (is_dynamic_size) {
      if (M1.get_row_count() != M2.get_row_count()) {
        throw std::range_error("Matrix dimension mismatch.");
      }
    }
    return self(M1.q + M2.q);
  }

  /// Negation operator with standard semantics.
  /// \return the negative of this matrix sum.
  self operator-() const { return self(-q); }

  /// Substraction operator with standard semantics.
  /// \param M the other matrix to be substracted from this.
  /// \return the matrix difference of this and M.
  /// \throw std::range_error if the matrix dimensions don't match.
  friend self operator-(const self& M1, const self& M2) {
    if constexpr (is_dynamic_size) {
      if (M1.get_row_count() != M2.get_row_count()) {
        throw std::range_error("Matrix dimension mismatch.");
      }
    }
    return self(M1.q - M2.q);
  }

  /// Matrix multiplication operator with a column vector.
  /// \param M a diagonal matrix.
  /// \param V the vector to be multiplied with this.
  /// \return the matrix-vector product of this and V.
  /// \throw std::range_error if the matrix-vector dimensions don't match.
  template <typename Vector, typename Result>
  void multiply_with_vector_rhs(const Vector& V, Result& R) const {
    R = elem_product(q, V);
  }

  /// Matrix multiplication operator with a row vector.
  /// \param M a diagonal matrix.
  /// \param V the vector to be multiplied with M.
  /// \return the matrix-vector product of this and V.
  /// \throw std::range_error if the matrix-vector dimensions don't match.
  template <typename Vector, typename Result>
  void multiply_with_vector_lhs(const Vector& V, Result& R) const {
    R = elem_product(q, V);
  }

  /*******************************************************************************
                           Special Methods
  *******************************************************************************/

  /// Extracts a diagonal sub-matrix from this matrix.
  /// \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
  /// \param aSizeOut Number of rows/columns of the sub-matrix.
  /// \return The diagonal sub-matrix contained in this matrix.
  /// \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
  friend mat<value_type, mat_structure::diagonal> get_block(const self& M,
                                                            int aDiagOffset,
                                                            int aSizeOut) {
    if (aDiagOffset + aSizeOut > M.get_row_count()) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    mat<value_type, mat_structure::diagonal> result(aSizeOut);
    for (int i = 0; i < aSizeOut; ++i) {
      result.q[i] = M.q[i + aDiagOffset];
    }
    return result;
  }

  /// Sets the sub-block of a matrix to a diagnoal sub-matrix M.
  /// \param M A diagonal matrix to which the sub-block will be set.
  /// \param subM A diagonal sub-matrix that will be written in the sub-part of this matrix.
  /// \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
  /// \return This matrix, by reference.
  /// \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
  template <unsigned int SubRowCount>
  friend self& set_block(self& M,
                         const mat<value_type, mat_structure::diagonal,
                                   Alignment, SubRowCount, SubRowCount>& subM,
                         int aDiagOffset) {
    if (aDiagOffset + subM.get_row_count() > M.get_row_count()) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    for (int i = 0; i < subM.get_row_count(); ++i) {
      M.q[i + aDiagOffset] = subM.q[i];
    }
    return M;
  }

  /// Appends the matrix 'rhs' to the end of the matrix 'lhs', which are both diagonal matrices.
  /// \param lhs The diagonal matrix to which to append the other.
  /// \param rhs The diagonal matrix to be appended to 'lhs'.
  friend void append_block_diag(self& lhs, const self& rhs) {
    int oldCount = lhs.get_row_count();
    lhs.set_row_count(oldCount + rhs.get_row_count(), true);
    set_block(lhs, rhs, oldCount);
  }

  /// Transposes the matrix M (which has no effect since M is diagonal, simply copies it).
  friend self transpose(const self& M) { return M; }

  /// Transposes the matrix M in a potentially destructive way (move-semantics).
  friend self transpose_move(self& M) {
    self result;
    swap(result, M);
    return result;
  }

  /// Transposes the matrix M in a potentially destructive way (move-semantics).
  friend self transpose(self&& M) { return self(std::move(M)); }

  /// Returns the trace of matrix M.
  friend value_type trace(const self& M) {
    value_type sum = value_type(0);
    for (auto& v : M.q) {
      sum += v;
    }
    return sum;
  }

  /// Inverts this matrix.
  void invert() {
    for (auto& v : q) {
      v = value_type(1.0) / v;
    }
  }

  /// Inverts the matrix M.
  friend self invert(self M) {
    M.invert();
    return M;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(q);
  }
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(q);
  }

  RK_RTTI_REGISTER_CLASS_1BASE(self, 1, serializable)
};

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_ALG_DIAGONAL_H_

/**
 * \file mat_alg_scalar.hpp
 *
 * This library declares matrix specializations for representing and manipulating scalar matrices.
 * This library implements many overloaded operators that turn out to be more efficiently implemented
 * if specialized for the scalar matrix case. All those overloads are automatically selected through
 * Sfinae switches, and the scalar matrix class is simply a partial specialization of the "ReaK::mat"
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

#ifndef REAK_MAT_ALG_SCALAR_HPP
#define REAK_MAT_ALG_SCALAR_HPP

#include "mat_alg_general.hpp"

#include <type_traits>

namespace ReaK {

/**
 * This class holds a scalar matrix. This class will hold only the scalar value and the dimension.
 *
 * Models: ReadableMatrixConcept, WritableMatrixConcept, and ResizableMatrixConcept.
 *
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or
 *mat_alignment::column_major (default).
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <class T, mat_alignment::tag Alignment, typename Allocator>
class mat<T, mat_structure::scalar, Alignment, Allocator>
    : public serializable {
 public:
  using self = mat<T, mat_structure::scalar, Alignment, Allocator>;
  using allocator_type = void;

  using value_type = T;

  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;

  using col_iterator = void;
  using const_col_iterator = void;
  using row_iterator = void;
  using const_row_iterator = void;

  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  static constexpr std::size_t static_row_count = 0;
  static constexpr std::size_t static_col_count = 0;
  static constexpr mat_alignment::tag alignment = Alignment;
  static constexpr mat_structure::tag structure = mat_structure::scalar;

 private:
  /// Holds the scalar entry.
  value_type q;
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
               const Allocator& aAlloc = Allocator())
      : q(aFill), rowCount(aRowCount) {}

  /**
   * Default constructor: sets all to zero.
   */
  explicit mat(const Allocator& aAlloc = Allocator())
      : mat(0, value_type(0), aAlloc) {}

  /**
   * Constructor for an identity matrix.
   */
  mat(size_type aRowCount, bool aIdentity,
      const Allocator& aAlloc = Allocator())
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
   * Constructor from a general matrix, conserving only the trace value.
   */
  template <typename Matrix>
  explicit mat(
      const Matrix& M, const Allocator& aAlloc = Allocator(),
      std::enable_if_t<
          is_readable_matrix_v<Matrix> && !std::is_same_v<Matrix, self>, void*>
          dummy = nullptr)
      : q(0.0),
        rowCount((M.get_row_count() < M.get_col_count() ? M.get_row_count()
                                                        : M.get_col_count())) {
    q = trace(M) / value_type(rowCount);
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
      return q;
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
      return q;
    }
    return value_type(0.0);
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
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<size_type, size_type> size() const noexcept {
    return {rowCount, rowCount};
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
  allocator_type get_allocator() const {}

  /*******************************************************************************
                           Assignment Operators
  *******************************************************************************/

  /**
   * Standard Assignment operator with a scalar matrix.
   */
  self& operator=(self M) {
    swap(*this, M);
    return *this;
  }

  /**
   * Standard Assignment operator with a general matrix. Copying only the scalar part of M.
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
    q += M.q;
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
    q -= M.q;
    return *this;
  }

  /**
   * Scalar-multiply-and-store operator with standard semantics.
   * \param S the scalar to be multiplied to this.
   * \return this matrix by reference.
   */
  self& operator*=(const T& S) {
    q *= S;
    return *this;
  }

  /**
   * Matrix-multiply-and-store operator with a scalar matrix.
   * \param M the other matrix to be multiplied with this.
   * \return this matrix by reference.
   * \throw std::range_error if the matrix dimensions don't match.
   */
  self& operator*=(const self& M) {
    if (rowCount != M.rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    q *= M.q;
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
    return self(M1.rowCount, M1.q + M2.q);
  }

  /**
   * Negation operator with standard semantics.
   * \return the negative of this matrix sum.
   */
  self operator-() const { return self(rowCount, -q); }

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
    return self(M1.rowCount, M1.q - M2.q);
  }

  /*******************************************************************************
                           Special Methods
  *******************************************************************************/

  /**
   * Extracts a scalar sub-matrix from this matrix.
   * \param aDiagOffset Number of rows/columns before the start of the sub-matrix rows/columns.
   * \param aSizeOut Number of rows/columns of the sub-matrix.
   * \return The scalar sub-matrix contained in this matrix.
   * \throw std::range_error If the sub-matrix's dimensions and position does not fit within this matrix.
   */
  friend self get_block(const self& M, int aDiagOffset, int aSizeOut) {
    if (aDiagOffset + aSizeOut > M.rowCount) {
      throw std::range_error("Matrix dimension mismatch.");
    }
    return self(aSizeOut, value_type(M.q));
  }

  /**
   * Transposes the matrix M (which has no effect since M is scalar, simply copies it).
   * \param M The scalar matrix to be transposed.
   * \return The transpose of M.
   */
  friend self transpose(const self& M) { return M; }

  /**
   * Transposes the matrix M in a potentially destructive way (move-semantics, pre-C++0x).
   * \param M The scalar matrix to be transposed and moved.
   * \return The transpose of M.
   */
  friend self transpose_move(self& M) {
    self result;
    swap(result, M);
    return result;
  }

  /**
   * Transposes the matrix M in a potentially destructive way (move-semantics, C++0x).
   * \param M The scalar matrix to be transposed and moved.
   * \return The transpose of M.
   */
  friend self transpose(self&& M) { return self(std::move(M)); }

  /**
   * Returns the trace of matrix M.
   * \param M A scalar matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace(const self& M) { return M.rowCount * M.q; }

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

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_readable_matrix<mat<T, mat_structure::scalar, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type =
      is_readable_matrix<mat<T, mat_structure::scalar, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_writable_matrix<mat<T, mat_structure::scalar, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type =
      is_writable_matrix<mat<T, mat_structure::scalar, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_fully_writable_matrix<
    mat<T, mat_structure::scalar, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = false;
  using type =
      is_writable_matrix<mat<T, mat_structure::scalar, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct is_resizable_matrix<
    mat<T, mat_structure::scalar, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type =
      is_resizable_matrix<mat<T, mat_structure::scalar, Alignment, Allocator>>;
};

template <typename T, mat_alignment::tag Alignment, typename Allocator>
struct has_allocator_matrix<
    mat<T, mat_structure::scalar, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = false;
  using type =
      has_allocator_matrix<mat<T, mat_structure::scalar, Alignment, Allocator>>;
};

}  // namespace ReaK

#endif

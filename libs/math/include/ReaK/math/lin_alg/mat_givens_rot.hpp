/**
 * \file mat_givens_rot.hpp
 *
 * This library provides the means to perform Givens rotations on matrices, efficiently.
 * This library provides a class template that can be constructed from the relevant elements
 * to be rotated by the Givens rotation, and it creates a representation of a Givens rotation
 * that can then be used to multiply a matrix (or sub-matrix) efficiently (with minimal calculations).
 *
 * Givens rotations are useful in a number of matrix numerical methods. It is a fundamental construct.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2011
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

#ifndef REAK_MAT_GIVENS_ROT_HPP
#define REAK_MAT_GIVENS_ROT_HPP

#include "mat_alg.hpp"
#include "mat_num_exceptions.hpp"

#include <type_traits>

namespace ReaK {

// Forward declarations.
template <typename T>
class givens_rot_matrix;

template <typename Matrix, typename T>
void givens_rot_prod(const givens_rot_matrix<T>& G, Matrix& A, int j = 0,
                     int k = 1);

template <typename Matrix, typename T>
void givens_rot_prod(Matrix& A, const givens_rot_matrix<T>& G, int j = 0,
                     int k = 1);

/**
 * This class represents a Givens rotation as a 2x2 rotation matrix. It can be constructed
 * by providing the two components of a column-vector (a,b) such that the resulting Givens
 * rotation G has the effect of zero-ing the second element: G * (a, b) = (r, 0)
 * (where (x,y) are a column-vectors). To apply a Givens rotation efficiently, it is preferred
 * to use the givens_rot_prod function templates, although this Givens rotation matrix models
 * a readable matrix concept and can thus be involved in other matrix arithmetic (although rarely
 * useful). In order to apply the Givens rotation to a matrix with larger dimensions, which is most
 * often the case, one can use the classes provided in mat_views.hpp and mat_composite_adaptor.hpp to
 * create such sub-matrices and composites of them to extract only the couple of rows or columns affected
 * by the rotation (the resulting sub-matrix (or matrix-view) should have the appropriate dimensions, either
 * row-count of 2 or column-count of 2, depending on the order of rotation). This class also provides
 * friend functions to perform the transposition of the Givens rotation (and effectively, its inverse).
 *
 * Models: ReadableMatrixConcept.
 *
 * \tparam T The value-type of the elements of the matrix.
 */
template <typename T>
class givens_rot_matrix {
 public:
  using self = givens_rot_matrix<T>;
  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using allocator_type = void;

  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;
  using iterator = T*;
  using const_iterator = const T*;

  using col_iterator = void;
  using const_col_iterator = void;
  using row_iterator = void;
  using const_row_iterator = void;

  static constexpr std::size_t static_row_count = 0;
  static constexpr std::size_t static_col_count = 0;
  static constexpr mat_alignment::tag alignment = mat_alignment::column_major;
  static constexpr mat_structure::tag structure = mat_structure::orthogonal;

  value_type c;  ///< Holds the cosine of the angle of rotation.
  value_type s;  ///< Holds the sine of the angle of rotation.

 private:
  void calculate_givens(const value_type& NumTol) {
    using std::abs;
    using std::sqrt;
    if (abs(s) < NumTol * abs(c)) {
      c = value_type(1.0);
      s = value_type(0.0);
      return;
    }

    if (abs(c) < abs(s)) {
      value_type tau = -c / s;
      s = value_type(1.0) / sqrt(value_type(1.0) + tau * tau);
      c = s * tau;
    } else {
      value_type tau = -s / c;
      c = value_type(1.0) / sqrt(value_type(1.0) + tau * tau);
      s = c * tau;
    }
  }

 public:
  /**
   * Set the Givens rotation by providing the two components of a column-vector (aA,aB) such that the
   * resulting Givens rotation G has the effect of zero-ing the second element: G * (a, b) = (r, 0)
   * (where (x,y) are a column-vectors).
   * \param aA The first component of the vector.
   * \param aB The second component of the vector.
   * \param NumTol The numerical tolerance used to assume a value to be zero (avoid division by zero).
   */
  void set(const value_type& aA, const value_type& aB,
           const value_type& NumTol = value_type(1E-8)) {
    c = aA;
    s = aB;
    calculate_givens(NumTol);
  }

  /**
   * Default constructor.
   */
  givens_rot_matrix() : c(1.0), s(0.0) {}

  /**
   * Constructs the Givens rotation by providing the two components of a column-vector (aA,aB) such that the
   * resulting Givens rotation G has the effect of zero-ing the second element: G * (a, b) = (r, 0)
   * (where (x,y) are a column-vectors).
   * \param aA The first component of the vector.
   * \param aB The second component of the vector.
   * \param NumTol The numerical tolerance used to assume a value to be zero (avoid division by zero).
   */
  explicit givens_rot_matrix(const value_type& aA, const value_type& aB,
                             const value_type& NumTol = value_type(1E-8))
      : c(aA), s(aB) {
    calculate_givens(NumTol);
  }

  // Default copy-constructor and assignment operators are correct.

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const { return 2; }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const { return 2; }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<size_type, size_type> size() const noexcept {
    return std::make_pair(2, 2);
  }

  /**
   * Returns the allocator object of the underlying container, which is none at all in this case.
   * \return the allocator object of the underlying container, which is none at all in this case.
   */
  allocator_type get_allocator() const {}

  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()(size_type i, size_type j) const {
    return (i == j ? c : (i < j ? -s : s));
  }

  template <typename Matrix>
  friend Matrix operator*(const Matrix& A, const self& G) {
    static_assert(is_fully_writable_matrix_v<Matrix>);
    Matrix result(A);
    givens_rot_prod(result, G);
    return result;
  }

  template <typename Matrix>
  friend Matrix operator*(const self& G, const Matrix& A) {
    static_assert(is_fully_writable_matrix_v<Matrix>);
    Matrix result(A);
    givens_rot_prod(G, result);
    return result;
  }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend self transpose(self M) {
    M.s = -M.s;
    return M;
  }
  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend self transpose_move(self M) {
    M.s = -M.s;
    return M;
  }
  /**
   * Returns the trace of the matrix M.
   * \param M The matrix.
   * \return The trace of M.
   */
  friend value_type trace(const self& M) { return M.c + M.c; }
};

template <typename T>
struct is_readable_matrix<givens_rot_matrix<T>> {
  static constexpr bool value = true;
  using type = is_readable_matrix<givens_rot_matrix<T>>;
};

template <typename T>
struct is_writable_matrix<givens_rot_matrix<T>> {
  static constexpr bool value = false;
  using type = is_writable_matrix<givens_rot_matrix<T>>;
};

template <typename T>
struct is_resizable_matrix<givens_rot_matrix<T>> {
  static constexpr bool value = false;
  using type = is_resizable_matrix<givens_rot_matrix<T>>;
};

template <typename T>
struct has_allocator_matrix<givens_rot_matrix<T>> {
  static constexpr bool value = false;
  using type = has_allocator_matrix<givens_rot_matrix<T>>;
};

template <typename T>
struct mat_product_priority<givens_rot_matrix<T>> {
  static constexpr std::size_t value =
      detail::product_priority<mat_structure::diagonal>::value + 1;
};

template <typename T>
struct mat_addition_priority<givens_rot_matrix<T>> {
  static constexpr std::size_t value =
      detail::addition_priority<mat_structure::square>::value;
};

template <typename T>
struct is_square_matrix<givens_rot_matrix<T>> {
  static constexpr bool value = true;
  using type = is_square_matrix<givens_rot_matrix<T>>;
};

template <typename T>
struct is_symmetric_matrix<givens_rot_matrix<T>> {
  static constexpr bool value = false;
  using type = is_symmetric_matrix<givens_rot_matrix<T>>;
};

template <typename T>
struct is_diagonal_matrix<givens_rot_matrix<T>> {
  static constexpr bool value = false;
  using type = is_diagonal_matrix<givens_rot_matrix<T>>;
};

/**
 * This function template allows for efficient post-multiplication of a matrix with a
 * Givens rotation matrix. This is generally more efficient then to perform a generic
 * matrix multiplication (and it is done in-place).
 * \tparam Matrix A fully-writable matrix type.
 * \tparam T A value-type which is compatible with the value-type of the Matrix type (for arithmetic).
 * \param A The matrix to be multiplied by the Givens rotation, stores, as output, the resulting matrix.
 * \param G The Givens rotation which will post-multiply A.
 * \param j The first column of A affected by the rotation (default 0).
 * \param k The second column of A affector by the rotation (default 1).
 */
template <typename Matrix, typename T>
void givens_rot_prod(Matrix& A, const givens_rot_matrix<T>& G, int j, int k) {
  static_assert(is_fully_writable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  using std::abs;

  if (abs(G.s) < std::numeric_limits<ValueType>::epsilon()) {
    return;
  }
  for (int i = 0; i < A.get_row_count(); ++i) {
    ValueType t0 = A(i, j);
    ValueType t1 = A(i, k);
    A(i, j) = G.c * t0 + G.s * t1;
    A(i, k) = G.c * t1 - G.s * t0;
  }
}

/**
 * This function template allows for efficient pre-multiplication of a matrix with a
 * Givens rotation matrix. This is generally more efficient then to perform a generic
 * matrix multiplication (and it is done in-place).
 * \tparam Matrix A fully-writable matrix type.
 * \tparam T A value-type which is compatible with the value-type of the Matrix type (for arithmetic).
 * \param G The Givens rotation which will pre-multiply A.
 * \param A The matrix to be multiplied by the Givens rotation, stores, as output, the resulting matrix.
 * \param j The first row of A affected by the rotation (default 0).
 * \param k The second row of A affector by the rotation (default 1).
 */
template <typename Matrix, typename T>
void givens_rot_prod(const givens_rot_matrix<T>& G, Matrix& A, int j, int k) {
  static_assert(is_fully_writable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  using std::abs;

  if (abs(G.s) < std::numeric_limits<ValueType>::epsilon()) {
    return;
  }
  for (int i = 0; i < A.get_col_count(); ++i) {
    ValueType t0 = A(j, i);
    ValueType t1 = A(k, i);
    A(j, i) = G.c * t0 - G.s * t1;
    A(k, i) = G.s * t0 + G.c * t1;
  }
}
}  // namespace ReaK

#endif

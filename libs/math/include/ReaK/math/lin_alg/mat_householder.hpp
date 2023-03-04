/**
 * \file mat_householder.hpp
 *
 * This library provides the means to perform Householder reflections on matrices, efficiently.
 * This library provides a class template that can be constructed from the relevant elements
 * to be reflected by the Householder reflection, and it creates a representation of a Householder reflection
 * that can then be used to multiply a matrix (or sub-matrix) efficiently (with minimal calculations).
 *
 * Householder reflection are useful in a number of matrix numerical methods. It is a fundamental construct.
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

#ifndef REAK_MAT_HOUSEHOLDER_HPP
#define REAK_MAT_HOUSEHOLDER_HPP

#include "mat_alg.hpp"
#include "mat_num_exceptions.hpp"

#include <type_traits>

namespace ReaK {

// Forward declarations.
template <typename Vector>
class householder_matrix;

template <typename Matrix, typename Vector>
void householder_prod(const householder_matrix<Vector>& P, Matrix& A);

template <typename Matrix, typename Vector>
void householder_prod(Matrix& A, const householder_matrix<Vector>& P);

/**
 * This class represents a Householder reflection as a NxN matrix. It can be constructed
 * by providing the column-vector u such that the resulting Householder reflection H has
 * the effect of zero-ing all but the first element: H * u = (a, 0) (where u is a column-vector
 * and a is some scalar value). To apply a Householder reflection efficiently, it is preferred
 * to use the householder_prod function templates, although this Householder reflection matrix models
 * a readable matrix concept and can thus be involved in other matrix arithmetic (although rarely
 * useful). In order to apply the Householder reflection to a matrix with larger dimensions, which
 * is most often the case, one can use the classes provided in mat_views.hpp and
 * mat_composite_adaptor.hpp to create such sub-matrices and composites of them to extract only
 * the rows or columns affected by the reflection (the resulting sub-matrix (or matrix-view)
 * should have the appropriate dimensions). This class also provides friend functions to
 * perform the transposition of the Householder reflection (and effectively, its inverse).
 *
 * Models: ReadableMatrixConcept.
 *
 * \tparam Vector The vector type to store the householder vector.
 */
template <typename Vector>
class householder_matrix {
 public:
  using self = householder_matrix<Vector>;
  using value_type = vect_value_type_t<Vector>;
  using size_type = typename vect_traits<Vector>::size_type;
  using difference_type = typename vect_traits<Vector>::difference_type;
  using allocator_type = typename vect_traits<Vector>::allocator_type;

  using pointer = typename vect_traits<Vector>::pointer;
  using const_pointer = typename vect_traits<Vector>::const_pointer;
  using reference = typename vect_traits<Vector>::reference;
  using const_reference = typename vect_traits<Vector>::const_reference;
  using iterator = typename vect_traits<Vector>::iterator;
  using const_iterator = typename vect_traits<Vector>::const_iterator;

  using col_iterator = void;
  using const_col_iterator = void;
  using row_iterator = void;
  using const_row_iterator = void;

  static constexpr std::size_t static_row_count = 0;
  static constexpr std::size_t static_col_count = 0;
  static constexpr mat_alignment::tag alignment = mat_alignment::column_major;
  static constexpr mat_structure::tag structure = mat_structure::symmetric;

  value_type beta;  ///< Holds the householder coefficient.
  Vector v;         ///< Holds the householder vector.

 private:
  void calculate_hhv(const value_type& NumTol) {
    if (v.size() == 0) {
      return;
    }

    using std::sqrt;

    auto sigma = value_type(0.0);
    for (int i = 1; i < v.size(); ++i) {
      sigma += v[i] * v[i];
    }
    if (sigma < NumTol * NumTol) {
      v[0] = value_type(1.0);
      beta = value_type(0.0);
      return;
    }
    value_type mu = sqrt(sigma + v[0] * v[0]);
    if (v[0] < NumTol) {
      v[0] -= mu;
    } else {
      v[0] = -sigma / (v[0] + mu);
    }
    beta = 2.0 / (sigma + v[0] * v[0]);
  }

  void calculate_inv_hhv(const value_type& NumTol) {
    if (v.size() == 0) {
      return;
    }

    using std::sqrt;

    auto sigma = value_type(0.0);
    for (int i = 0; i < v.size() - 1; ++i) {
      sigma += v[i] * v[i];
    }
    if (sigma < NumTol * NumTol) {
      v[v.size() - 1] = value_type(1.0);
      beta = value_type(0.0);
      return;
    }
    value_type mu = sqrt(sigma + v[v.size() - 1] * v[v.size() - 1]);
    if (v[v.size() - 1] < NumTol) {
      v[v.size() - 1] -= mu;
    } else {
      v[v.size() - 1] = -sigma / (v[v.size() - 1] + mu);
    }
    beta = 2.0 / (sigma + v[v.size() - 1] * v[v.size() - 1]);
  }

 public:
  /**
   * Set the Householder reflection by providing a column-vector u such that the
   * resulting Householder reflection H has the effect of zero-ing the second element: H * u = (a, 0).
   * \tparam Vector2 A readable vector type.
   * \param aE The vector from which to calculate the Householder reflection (vector u).
   * \param NumTol The numerical tolerance used to assume a value to be zero (avoid division by zero).
   */
  template <typename Vector2>
  void set(const Vector2& aE, const value_type& NumTol = value_type(1E-8)) {
    static_assert(is_readable_vector_v<Vector2>);
    beta = value_type(0.0);
    v = aE;
    calculate_hhv(NumTol);
  }
  /**
   * Set the Householder reflection by providing a column-vector u such that the
   * resulting Householder reflection H has the effect of zero-ing the first element: H * u = (0, a).
   * \note This constructs an "opposite" Householder reflection (can be useful for row elimintations).
   * \tparam Vector2 A readable vector type.
   * \param aE The vector from which to calculate the Householder reflection (vector u).
   * \param NumTol The numerical tolerance used to assume a value to be zero (avoid division by zero).
   */
  template <typename Vector2>
  void set_inv(const Vector2& aE, const value_type& NumTol = value_type(1E-8)) {
    static_assert(is_readable_vector_v<Vector2>);
    beta = value_type(0.0);
    v = aE;
    calculate_inv_hhv(NumTol);
  }
  /**
   * Default constructor.
   */
  explicit householder_matrix(const value_type& aBeta,
                              const Vector& aV = Vector())
      : beta(aBeta), v(aV) {}

  householder_matrix() : householder_matrix(value_type(0.0)) {}

  /**
   * Constructs the Householder reflection by providing a column-vector u such that the
   * resulting Householder reflection H has the effect of zero-ing the second element: H * u = (a, 0).
   * \tparam Vector2 A readable vector type.
   * \param aE The vector from which to calculate the Householder reflection (vector u).
   * \param NumTol The numerical tolerance used to assume a value to be zero (avoid division by zero).
   */
  template <typename Vector2>
  explicit householder_matrix(const Vector2& aE,
                              const value_type& NumTol = value_type(1E-8))
      : beta(0.0), v(aE) {
    static_assert(is_readable_vector_v<Vector2>);
    calculate_hhv(NumTol);
  }

  /**
   * Standard copy-constructor.
   */
  householder_matrix(const self& rhs) : beta(rhs.beta), v(rhs.v) {}

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(lhs.beta, rhs.beta);
    swap(lhs.v, rhs.v);
  }

  /**
   * Standard assignment operator (copy-and-swap and move-and-swap).
   */
  self& operator=(self rhs) {
    swap(*this, rhs);
    return *this;
  }

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const { return v.size(); }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const { return v.size(); }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<size_type, size_type> size() const noexcept {
    return std::make_pair(v.size(), v.size());
  }

  /**
   * Returns the allocator object of the underlying container, which is none at all in this case.
   * \return the allocator object of the underlying container, which is none at all in this case.
   */
  allocator_type get_allocator() const { return v.get_allocator(); }

  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()(int i, int j) const {
    return (i == j ? value_type(1.0) : value_type(0.0)) - beta * v[i] * v[j];
  }

  template <typename Matrix>
  friend Matrix operator*(const Matrix& A, const self& P) {
    static_assert(is_fully_writable_matrix_v<Matrix>);
    Matrix result(A);
    householder_prod(result, P);
    return result;
  }

  template <typename Matrix>
  friend Matrix operator*(const self& P, const Matrix& A) {
    static_assert(is_fully_writable_matrix_v<Matrix>);
    Matrix result(A);
    householder_prod(P, result);
    return result;
  }

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose(const self& M) { return M; }
  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend self transpose_move(self& M) {
    self result;
    swap(M, result);
    return result;
  }

  /**
   * Returns the trace of the matrix M.
   * \param M The matrix.
   * \return The trace of M.
   */
  friend value_type trace(const self& M) {
    return value_type(M.v.size()) - M.beta * (M.v * M.v);
  }
};

template <typename Vector>
struct is_readable_matrix<householder_matrix<Vector>> {
  static constexpr bool value = true;
  using type = is_readable_matrix<householder_matrix<Vector>>;
};

template <typename Vector>
struct is_writable_matrix<householder_matrix<Vector>> {
  static constexpr bool value = false;
  using type = is_writable_matrix<householder_matrix<Vector>>;
};

template <typename Vector>
struct is_resizable_matrix<householder_matrix<Vector>> {
  static constexpr bool value = false;
  using type = is_resizable_matrix<householder_matrix<Vector>>;
};

template <typename Vector>
struct has_allocator_matrix<householder_matrix<Vector>> {
  static constexpr bool value = has_allocator_vector<Vector>::value;
  using type = has_allocator_matrix<householder_matrix<Vector>>;
};

template <typename Vector>
struct mat_product_priority<householder_matrix<Vector>> {
  static constexpr std::size_t value =
      detail::product_priority<mat_structure::diagonal>::value + 1;
};

template <typename Vector>
struct mat_addition_priority<householder_matrix<Vector>> {
  static constexpr std::size_t value =
      detail::addition_priority<mat_structure::square>::value;
};

template <typename Vector>
struct is_square_matrix<householder_matrix<Vector>> {
  static constexpr bool value = true;
  using type = is_square_matrix<householder_matrix<Vector>>;
};

template <typename Vector>
struct is_symmetric_matrix<householder_matrix<Vector>> {
  static constexpr bool value = true;
  using type = is_symmetric_matrix<householder_matrix<Vector>>;
};

template <typename Vector>
struct is_diagonal_matrix<householder_matrix<Vector>> {
  static constexpr bool value = false;
  using type = is_diagonal_matrix<householder_matrix<Vector>>;
};

/**
 * This function template allows for efficient post-multiplication of a matrix with a
 * Householder reflection matrix. This is generally more efficient then to perform a generic
 * matrix multiplication (and it is done in-place).
 * \tparam Matrix A fully-writable matrix type.
 * \tparam Vector A vector-type which is compatible with the value-type of the Matrix type (for arithmetic).
 * \param A The matrix to be multiplied by the Householder reflection, stores, as output, the resulting matrix.
 * \param P The Householder reflection which will post-multiply A.
 */
template <typename Matrix, typename Vector>
void householder_prod(Matrix& A, const householder_matrix<Vector>& P) {
  static_assert(is_fully_writable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  using std::abs;

  if (abs(P.beta) < std::numeric_limits<ValueType>::epsilon()) {
    return;
  }
  for (int i = 0; i < A.get_row_count(); ++i) {
    auto temp = ValueType(0);
    for (int j = 0; j < P.v.size(); ++j) {
      temp += P.beta * A(i, j) * P.v[j];
    }
    for (int j = 0; j < P.v.size(); ++j) {
      A(i, j) -= temp * P.v[j];
    }
  }
}

/**
 * This function template allows for efficient pre-multiplication of a matrix with a
 * Householder reflection matrix. This is generally more efficient then to perform a generic
 * matrix multiplication (and it is done in-place).
 * \tparam Matrix A fully-writable matrix type.
 * \tparam Vector A vector-type which is compatible with the value-type of the Matrix type (for arithmetic).
 * \param A The matrix to be multiplied by the Householder reflection, stores, as output, the resulting matrix.
 * \param P The Householder reflection which will pre-multiply A.
 */
template <typename Matrix, typename Vector>
void householder_prod(const householder_matrix<Vector>& P, Matrix& A) {
  static_assert(is_fully_writable_matrix_v<Matrix>);
  using ValueType = mat_value_type_t<Matrix>;
  using std::abs;

  if (abs(P.beta) < std::numeric_limits<ValueType>::epsilon()) {
    return;
  }
  for (int i = 0; i < A.get_col_count(); ++i) {
    auto temp = ValueType(0);
    for (int j = 0; j < P.v.size(); ++j) {
      temp += P.beta * A(j, i) * P.v[j];
    }
    for (int j = 0; j < P.v.size(); ++j) {
      A(j, i) -= temp * P.v[j];
    }
  }
}
}  // namespace ReaK

#endif

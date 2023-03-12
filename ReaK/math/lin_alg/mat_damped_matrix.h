/**
 * \file mat_damped_matrix.h
 *
 * This library provides an adaptor class that represents the addition of a diagonal matrix and
 * square matrix, i.e. a so-called damped matrix.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
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

#ifndef REAK_MATH_LIN_ALG_MAT_DAMPED_MATRIX_H_
#define REAK_MATH_LIN_ALG_MAT_DAMPED_MATRIX_H_

#include "ReaK/math/lin_alg/mat_alg_general.h"

namespace ReaK {

/**
 * This class template forms the addition of two matrices, one diagonal and one square.
 *
 * Models: ReadableMatrixConcept.
 *
 * \tparam SquareMatrix Matrix type for the left matrix.
 * \tparam DiagMatrix Matrix type for the right matrix.
 */
template <typename SquareMatrix, typename DiagMatrix>
class mat_damped_matrix {
 public:
  using self = mat_damped_matrix<SquareMatrix, DiagMatrix>;

  using value_type = mat_value_type_t<SquareMatrix>;

  using reference = typename mat_traits<SquareMatrix>::reference;
  using const_reference = typename mat_traits<SquareMatrix>::const_reference;
  using pointer = typename mat_traits<SquareMatrix>::pointer;
  using const_pointer = typename mat_traits<SquareMatrix>::const_pointer;

  using col_iterator = typename mat_traits<SquareMatrix>::col_iterator;
  using const_col_iterator =
      typename mat_traits<SquareMatrix>::const_col_iterator;
  using row_iterator = typename mat_traits<SquareMatrix>::row_iterator;
  using const_row_iterator =
      typename mat_traits<SquareMatrix>::const_row_iterator;

  using size_type = mat_size_type_t<SquareMatrix>;
  using difference_type = typename mat_traits<SquareMatrix>::difference_type;

  static constexpr unsigned int static_row_count =
      MatStaticSizeIfExpectedEqual(mat_traits<SquareMatrix>::static_row_count,
                                   mat_traits<DiagMatrix>::static_row_count);
  static constexpr unsigned int static_col_count = static_row_count;
  static constexpr mat_alignment::tag alignment =
      mat_traits<SquareMatrix>::alignment;
  static constexpr mat_structure::tag structure =
      mat_traits<SquareMatrix>::structure;

 private:
  const SquareMatrix* m_sqr;  ///< Holds the left part of the matrix.
  const DiagMatrix* m_diag;   ///< Holds the right part of the matrix.
 public:
  /**
   * Parametrized constructor.
   * \param aML Matrix to fill the left part of the matrix.
   * \param aMR Matrix to fill the right part of the matrix.
   */
  mat_damped_matrix(const SquareMatrix& aMSqr, const DiagMatrix& aMDiag)
      : m_sqr(&aMSqr), m_diag(&aMDiag) {
    if (m_sqr->get_row_count() != m_diag->get_row_count()) {
      throw std::range_error("Matrix dimensions mismatch.");
    }
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
  value_type operator()(int i, int j) const {
    if (i == j) {
      return (*m_sqr)(i, i) + (*m_diag)(i, i);
    }
    return (*m_sqr)(i, j);
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
  size_type get_row_count() const noexcept { return m_sqr->get_row_count(); }
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const noexcept { return m_sqr->get_col_count(); }

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair<size_type, size_type> size() const noexcept {
    return {m_sqr->get_row_count(), m_sqr->get_col_count()};
  }
};

template <typename SquareMatrix, typename DiagMatrix>
struct is_readable_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>> {
  static constexpr bool value =
      is_readable_matrix_v<SquareMatrix> && is_readable_matrix_v<DiagMatrix>;
  using type = is_readable_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>>;
};

template <typename SquareMatrix, typename DiagMatrix>
struct is_writable_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>> {
  static constexpr bool value = false;
  using type = is_writable_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>>;
};

template <typename SquareMatrix, typename DiagMatrix>
struct is_fully_writable_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>> {
  static constexpr bool value = false;
  using type =
      is_fully_writable_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>>;
};

template <typename SquareMatrix, typename DiagMatrix>
struct is_row_resizable_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>> {
  static constexpr bool value = false;
  using type =
      is_row_resizable_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>>;
};

template <typename SquareMatrix, typename DiagMatrix>
struct is_col_resizable_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>> {
  static constexpr bool value = false;
  using type =
      is_col_resizable_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>>;
};

template <typename SquareMatrix, typename DiagMatrix>
mat_damped_matrix<SquareMatrix, DiagMatrix> make_damped_matrix(
    const SquareMatrix& aMSqr, const DiagMatrix& aMDiag) {
  return mat_damped_matrix<SquareMatrix, DiagMatrix>(aMSqr, aMDiag);
};

template <typename SquareMatrix, typename DiagMatrix>
struct is_square_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>> {
  using value_type = bool;
  static constexpr bool value = is_square_matrix_v<SquareMatrix>;
  using type = is_square_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>>;
};

template <typename SquareMatrix, typename DiagMatrix>
struct is_symmetric_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>> {
  using value_type = bool;
  static constexpr bool value = is_symmetric_matrix_v<SquareMatrix>;
  using type = is_symmetric_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>>;
};

template <typename SquareMatrix, typename DiagMatrix>
struct is_diagonal_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>> {
  using value_type = bool;
  static constexpr bool value = is_diagonal_matrix_v<SquareMatrix>;
  using type = is_diagonal_matrix<mat_damped_matrix<SquareMatrix, DiagMatrix>>;
};

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_DAMPED_MATRIX_H_

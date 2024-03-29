/**
 * \file mat_concepts.h
 *
 * This header declares the various concepts to which a Matrix class can be expected
 * to fulfill. These concepts include ReadableMatrix, WritableMatrix,
 * and ResizableMatrix. All these concepts are also
 * paired with meta-functions that can evaluate whether a Matrix class fulfill the
 * concept or not, and return a compile-time constant bool (on the model of
 * std::bool_constant class). Note that these meta-functions cannot really check the
 * concepts directly (this is impossible in current and future C++ standard versions, might
 * eventually be part of the standard, but not in the forseeable future). These meta-functions
 * are constant meta-functions that always return false, the implementer of a given matrix
 * class should also define specializations of these meta-functions with their appropriate
 * return values.
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

#ifndef REAK_MATH_LIN_ALG_MAT_CONCEPTS_H_
#define REAK_MATH_LIN_ALG_MAT_CONCEPTS_H_

#include "ReaK/math/lin_alg/mat_traits.h"

#include <concepts>
#include <iterator>
#include <type_traits>

namespace ReaK {

/**
 * This concept will fail to be instantiated if the Matrix class does not model
 * the readable matrix concept, meaning that its (i,j) operator can return a readable
 * value, and that its row and column count can be obtained (by calling get_row_count() and get_col_count()).
 *
 * Required expressions for Matrix m:
 *  e = m(i,j)   read access to the elements of m by a row and column index.
 *  s = m.get_col_count()  can obtain the number of columns of the matrix.
 *  s = m.get_row_count()  can obtain the number of rows of the matrix.
 */
template <typename Matrix>
concept ReadableMatrix = requires(const Matrix& m) {
  { m(0, 0) } -> std::convertible_to<mat_value_type_t<Matrix>>;
  { m.get_col_count() } -> std::integral<>;
  { m.get_row_count() } -> std::integral<>;
};

// Legacy
template <typename Matrix>
static constexpr bool is_readable_matrix_v = ReadableMatrix<Matrix>;

/**
 * This concept will fail to be instantiated if the Matrix class does not model
 * the writable matrix concept, meaning that its (i,j) operator can return a writable
 * value, and that it fulfills the ReadableMatrix.
 *
 * Required expressions for Matrix m in addition to that of ReadableMatrix:
 *  m(i,j) = e;   write access to the elements of m by a row and column index.
 */
template <typename Matrix>
concept WritableMatrix = ReadableMatrix<Matrix>&& requires(Matrix& m) {
  { m(0, 0) } -> std::assignable_from<mat_value_type_t<Matrix>>;
};

// Legacy
template <typename Matrix>
static constexpr bool is_writable_matrix_v = WritableMatrix<Matrix>;

/**
 * This meta-function evaluates whether a Matrix class fulfills the WritableMatrix and
 * can be considered as "fully writable" meaning that all the (i,j) values are independent and writable,
 * however, it does not attempt to instantiate the template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results
 * in a false value, and the implementer of a matrix class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new matrix class.
 */
template <typename Matrix>
static constexpr bool is_fully_writable_matrix_v = false;

template <typename Matrix>
concept FullyWritableMatrix =
    WritableMatrix<Matrix>&& is_fully_writable_matrix_v<Matrix>;

/**
 * This concept will fail to be instantiated if the Matrix class does not model
 * the row-resizable matrix concept, meaning that its row counts can be
 * set to some values which results in the matrix having at least that row
 * count (the functions are set_row_count()).
 *
 * Required expressions for Matrix m:
 *  m.set_row_count(s)  can set the number of rows of the matrix.
 */
template <typename Matrix>
concept RowResizableMatrix = ReadableMatrix<Matrix>&& requires(Matrix& m,
                                                               int sz) {
  m.set_row_count(sz);
};

// Legacy
template <typename Matrix>
static constexpr bool is_row_resizable_matrix_v = RowResizableMatrix<Matrix>;

/**
 * This concept will fail to be instantiated if the Matrix class does not model
 * the column-resizable matrix concept, meaning that its column counts can be
 * set to some values which results in the matrix having at least that
 * column count (the functions are set_col_count()).
 *
 * Required expressions for Matrix m:
 *  m.set_col_count(s)  can set the number of columns of the matrix.
 */
template <typename Matrix>
concept ColResizableMatrix = ReadableMatrix<Matrix>&& requires(Matrix& m,
                                                               int sz) {
  m.set_col_count(sz);
};

// Legacy
template <typename Matrix>
static constexpr bool is_col_resizable_matrix_v = ColResizableMatrix<Matrix>;

/**
 * This meta-function evaluates whether a Matrix class is a square matrix. The implementer of
 * a matrix class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new matrix class.
 */
template <typename Matrix>
static constexpr bool is_square_matrix_v = false;

template <typename Matrix>
concept SquareMatrix = ReadableMatrix<Matrix>&& is_square_matrix_v<Matrix>;

/**
 * This meta-function evaluates whether a Matrix class is a symmetric matrix. The implementer of
 * a matrix class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new matrix class.
 */
template <typename Matrix>
static constexpr bool is_symmetric_matrix_v = false;

template <typename Matrix>
concept SymmetricMatrix = SquareMatrix<Matrix>&& is_symmetric_matrix_v<Matrix>;

/**
 * This meta-function evaluates whether a Matrix class is a diagonal matrix. The implementer of
 * a matrix class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new matrix class.
 */
template <typename Matrix>
static constexpr bool is_diagonal_matrix_v = false;

template <typename Matrix>
concept DiagonalMatrix = SymmetricMatrix<Matrix>&& is_diagonal_matrix_v<Matrix>;

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_CONCEPTS_H_

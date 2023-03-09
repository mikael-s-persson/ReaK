/**
 * \file mat_concepts.h
 *
 * This header declares the various concepts to which a Matrix class can be expected
 * to fulfill. These concepts include ReadableMatrixConcept, WritableMatrixConcept,
 * ResizableMatrixConcept, and DynAllocMatrixConcept. All these concepts are also
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
#include "boost/concept/requires.hpp"
#include "boost/concept_check.hpp"

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
struct ReadableMatrixConcept {
  Matrix m;

  mat_size_type_t<Matrix> s;
  mat_value_type_t<Matrix> e;

  using self = ReadableMatrixConcept<Matrix>;

  BOOST_CONCEPT_USAGE(ReadableMatrixConcept) {
    e = m(0, 0);  // can be indexed and given an rvalue
    s = m.get_col_count();
    s = m.get_row_count();
  }
};

/**
 * This meta-function evaluates whether a Matrix class fulfills the ReadableMatrixConcept,
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results
 * in a false value, and the implementer of a matrix class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new matrix class.
 */
template <typename Matrix>
struct is_readable_matrix {
  using value_type = bool;
  static constexpr bool value = false;
  using type = is_readable_matrix<Matrix>;
};

template <typename Matrix>
static constexpr bool is_readable_matrix_v = is_readable_matrix<Matrix>::value;

/**
 * This concept will fail to be instantiated if the Matrix class does not model
 * the writable matrix concept, meaning that its (i,j) operator can return a writable
 * value, and that it fulfills the ReadableMatrixConcept.
 *
 * Required expressions for Matrix m in addition to that of ReadableMatrixConcept:
 *  m(i,j) = e;   write access to the elements of m by a row and column index.
 */
template <typename Matrix>
struct WritableMatrixConcept
    : ReadableMatrixConcept<Matrix> {  // must also be readable.
  Matrix m;

  mat_value_type_t<Matrix> r;

  BOOST_CONCEPT_USAGE(WritableMatrixConcept) {
    m(0, 0) = r;  // can be indexed and given an lvalue
  }
};

/**
 * This meta-function evaluates whether a Matrix class fulfills the WritableMatrixConcept,
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results
 * in a false value, and the implementer of a matrix class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new matrix class.
 */
template <typename Matrix>
struct is_writable_matrix {
  using value_type = bool;
  static constexpr bool value = false;
  using type = is_writable_matrix<Matrix>;
};

template <typename Matrix>
static constexpr bool is_writable_matrix_v = is_writable_matrix<Matrix>::value;

/**
 * This meta-function evaluates whether a Matrix class fulfills the WritableMatrixConcept and
 * can be considered as "fully writable" meaning that all the (i,j) values are independent and writable,
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results
 * in a false value, and the implementer of a matrix class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new matrix class.
 */
template <typename Matrix>
struct is_fully_writable_matrix {
  using value_type = bool;
  static constexpr bool value = false;
  using type = is_fully_writable_matrix<Matrix>;
};

template <typename Matrix>
static constexpr bool is_fully_writable_matrix_v =
    is_fully_writable_matrix<Matrix>::value;

/**
 * This concept will fail to be instantiated if the Matrix class does not model
 * the resizable matrix concept, meaning that its row and column counts can be
 * set to some values which results in the matrix having at least that row or
 * column count (the functions are set_row_count() and set_col_count()).
 *
 * Required expressions for Matrix m:
 *  m.set_col_count(s)  can set the number of columns of the matrix.
 *  m.set_row_count(s)  can set the number of rows of the matrix.
 */
template <typename Matrix>
struct ResizableMatrixConcept {
  Matrix m;

  mat_size_type_t<Matrix> sz;

  BOOST_CONCEPT_USAGE(ResizableMatrixConcept) {
    m.set_row_count(sz);
    m.set_col_count(sz);
  }
};

/**
 * This meta-function evaluates whether a Matrix class fulfills the ResizableMatrixConcept,
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results
 * in a false value, and the implementer of a matrix class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new matrix class.
 */
template <typename Matrix>
struct is_resizable_matrix {
  using value_type = bool;
  static constexpr bool value = false;
  using type = is_resizable_matrix<Matrix>;
};

template <typename Matrix>
static constexpr bool is_resizable_matrix_v =
    is_resizable_matrix<Matrix>::value;

/**
 * This concept will fail to be instantiated if the Matrix class does not model
 * the dynamically allocated matrix concept, meaning that it stores its elements in
 * dynamically allocated memory. The only requirement to fulfill this concept is
 * to be able to obtain the allocate object associated to a matrix (here "allocator" is
 * used with the exact same meaning as "STL allocators").
 *
 * Required expressions for Matrix m:
 *  al = m.get_allocator()  can obtain the allocator object of the matrix.
 */
template <typename Matrix>
struct DynAllocMatrixConcept : ResizableMatrixConcept<Matrix> {
  Matrix m;

  typename mat_traits<Matrix>::allocator_type al;

  BOOST_CONCEPT_USAGE(DynAllocMatrixConcept) { al = m.get_allocator(); }
};

/**
 * This meta-function evaluates whether a Matrix class fulfills the DynAllocMatrixConcept,
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results
 * in a false value, and the implementer of a matrix class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new matrix class.
 */
template <typename Matrix>
struct has_allocator_matrix {
  using value_type = bool;
  static constexpr bool value = false;
  using type = has_allocator_matrix<Matrix>;
};

template <typename Matrix>
static constexpr bool has_allocator_matrix_v =
    has_allocator_matrix<Matrix>::value;

/**
 * This meta-function evaluates whether a Matrix class is a square matrix. The implementer of
 * a matrix class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new matrix class.
 */
template <typename Matrix>
struct is_square_matrix {
  using value_type = bool;
  static constexpr bool value = false;
  using type = is_square_matrix<Matrix>;
};

template <typename Matrix>
static constexpr bool is_square_matrix_v = is_square_matrix<Matrix>::value;

/**
 * This meta-function evaluates whether a Matrix class is a symmetric matrix. The implementer of
 * a matrix class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new matrix class.
 */
template <typename Matrix>
struct is_symmetric_matrix {
  using value_type = bool;
  static constexpr bool value = false;
  using type = is_symmetric_matrix<Matrix>;
};

template <typename Matrix>
static constexpr bool is_symmetric_matrix_v =
    is_symmetric_matrix<Matrix>::value;

/**
 * This meta-function evaluates whether a Matrix class is a diagonal matrix. The implementer of
 * a matrix class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new matrix class.
 */
template <typename Matrix>
struct is_diagonal_matrix {
  using value_type = bool;
  static constexpr bool value = false;
  using type = is_diagonal_matrix<Matrix>;
};

template <typename Matrix>
static constexpr bool is_diagonal_matrix_v = is_diagonal_matrix<Matrix>::value;

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_CONCEPTS_H_

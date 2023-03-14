/**
 * \file mat_traits.h
 *
 * This header declares a number of standard matrix type traits. These include
 * matrix alignment tags (column-major and row-major), matrix structure tags
 * (rectangular, square, symmetric, diagonal, skew-symmetric, etc.), some ReaK::rtti
 * template specializations for these tags, the main mat_traits<> template which
 * defines the nested typedefs related to a matrix class (and required by an implementation
 * of a matrix class), and, finally, a series of meta-functions to compute the product-priority
 * additive-priority of matrix classes based on their structural tags.
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

#ifndef REAK_MATH_LIN_ALG_MAT_TRAITS_H_
#define REAK_MATH_LIN_ALG_MAT_TRAITS_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/core/rtti/rtti.h"

#include <type_traits>
#include "boost/concept_check.hpp"

namespace ReaK {

namespace mat_alignment {

/**
 * These tags are used to mark a matrix as having a column-major or row-major alignment
 * in memory (note that array-of-array is not an option, due to the obvious performance problems
 * with such a memory model for matrix representations).
 */
enum tag { column_major = 1, row_major = 2 };

}  // namespace mat_alignment

namespace mat_structure {

/**
 * These tags are used to identify a matrix class by its underlying structure. The most
 * general structure is, of course, the rectangular matrix (any size, any shape, and densely
 * populated). Then, the frequently used tags are square, symmetric, skew_symmetric, and
 * diagonal. Additionally, the nil and identity tags signify that a matrix type is constraint
 * to never be anything other than nil or identity (obviously, these types will be read-only
 * and do not require storage other than the size information).
 */
enum tag {
  rectangular = 1,
  square = 2,
  symmetric = 3,
  skew_symmetric = 4,
  diagonal = 5,
  upper_triangular = 6,
  lower_triangular = 7,
  orthogonal = 8,
  tridiagonal = 9,
  nil = 10,
  identity = 11,
  scalar = 12,
  permutation = 13
};

}  // namespace mat_structure

/// Gets the static matrix size (either row or col) from two sizes that are expected to be equal
/// at run-time. In other words, if either is non-zero (statically sized matrix dimension), it
/// returns that value, expecting that the dynamically sized matrix will have the same size at
/// run-time, thus it is safe to have a statically sized result.
inline constexpr unsigned int MatStaticSizeIfExpectedEqual(unsigned int sz1,
                                                           unsigned int sz2) {
  return ((sz1 == 0 && sz2 == 0) ? 0 : std::max(sz1, sz2));
}

/// Gets the static matrix size (either row or col) from two sizes that will be concatenated.
/// In other words, if either is zero (dynamically sized matrix dimension), it
/// returns zero since no static overall size could be resolved in that case.
inline constexpr unsigned int MatStaticSizeIfConcat(unsigned int sz1,
                                                    unsigned int sz2) {
  return ((sz1 == 0 || sz2 == 0) ? 0 : (sz1 + sz2));
}

/*
 * The following are ReaK::rtti class templates which are used to associate
 * type information to the matrix tags (mat_alignment::tag and mat_structure::tag).
 *
 * Please ignore, unless interested in the inner-workings of the ReaK::rtti system.
 */
namespace rtti {

template <>
struct get_type_id<
    std::integral_constant<mat_alignment::tag, mat_alignment::column_major>> {
  static constexpr unsigned int ID = 1;
  static constexpr auto type_name = std::string_view{"column_major"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <>
struct get_type_id<
    std::integral_constant<mat_alignment::tag, mat_alignment::row_major>> {
  static constexpr unsigned int ID = 2;
  static constexpr auto type_name = std::string_view{"row_major"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <mat_alignment::tag U, typename Tail>
struct get_type_info<std::integral_constant<mat_alignment::tag, U>, Tail> {
  using type = type_id<std::integral_constant<mat_alignment::tag, U>,
                       typename Tail::type>;
  static constexpr auto type_name = ct_concat_v<
      get_type_id<std::integral_constant<mat_alignment::tag, U>>::type_name,
      get_type_name_tail<Tail>::value>;
};

template <>
struct get_type_id<
    std::integral_constant<mat_structure::tag, mat_structure::rectangular>> {
  static constexpr unsigned int ID = 1;
  static constexpr auto type_name = std::string_view{"rectangular"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <>
struct get_type_id<
    std::integral_constant<mat_structure::tag, mat_structure::square>> {
  static constexpr unsigned int ID = 2;
  static constexpr auto type_name = std::string_view{"square"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <>
struct get_type_id<
    std::integral_constant<mat_structure::tag, mat_structure::symmetric>> {
  static constexpr unsigned int ID = 3;
  static constexpr auto type_name = std::string_view{"symmetric"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <>
struct get_type_id<
    std::integral_constant<mat_structure::tag, mat_structure::skew_symmetric>> {
  static constexpr unsigned int ID = 4;
  static constexpr auto type_name = std::string_view{"skew_symmetric"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <>
struct get_type_id<
    std::integral_constant<mat_structure::tag, mat_structure::diagonal>> {
  static constexpr unsigned int ID = 5;
  static constexpr auto type_name = std::string_view{"diagonal"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <>
struct get_type_id<std::integral_constant<mat_structure::tag,
                                          mat_structure::upper_triangular>> {
  static constexpr unsigned int ID = 6;
  static constexpr auto type_name = std::string_view{"upper_triangular"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <>
struct get_type_id<std::integral_constant<mat_structure::tag,
                                          mat_structure::lower_triangular>> {
  static constexpr unsigned int ID = 7;
  static constexpr auto type_name = std::string_view{"lower_triangular"};
  static construct_ptr CreatePtr() { return nullptr; }
};

template <>
struct get_type_id<
    std::integral_constant<mat_structure::tag, mat_structure::orthogonal>> {
  static constexpr unsigned int ID = 8;
  static constexpr auto type_name = std::string_view{"orthogonal"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <>
struct get_type_id<
    std::integral_constant<mat_structure::tag, mat_structure::tridiagonal>> {
  static constexpr unsigned int ID = 9;
  static constexpr auto type_name = std::string_view{"tridiagonal"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <>
struct get_type_id<
    std::integral_constant<mat_structure::tag, mat_structure::nil>> {
  static constexpr unsigned int ID = 10;
  static constexpr auto type_name = std::string_view{"nil"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <>
struct get_type_id<
    std::integral_constant<mat_structure::tag, mat_structure::identity>> {
  static constexpr unsigned int ID = 11;
  static constexpr auto type_name = std::string_view{"identity"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <>
struct get_type_id<
    std::integral_constant<mat_structure::tag, mat_structure::scalar>> {
  static constexpr unsigned int ID = 12;
  static constexpr auto type_name = std::string_view{"scalar"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <>
struct get_type_id<
    std::integral_constant<mat_structure::tag, mat_structure::permutation>> {
  static constexpr unsigned int ID = 13;
  static constexpr auto type_name = std::string_view{"permutation"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

template <mat_structure::tag U, typename Tail>
struct get_type_info<std::integral_constant<mat_structure::tag, U>, Tail> {
  using type = type_id<std::integral_constant<mat_structure::tag, U>,
                       typename Tail::type>;
  static constexpr auto type_name = ct_concat_v<
      get_type_id<std::integral_constant<mat_structure::tag, U>>::type_name,
      get_type_name_tail<Tail>::value>;
};

}  // namespace rtti

/**
 * This type-traits definition describes the nested typedefs that are expected
 * from an implementation of a matrix class. They are mostly inspired from
 * STL-containers' traits. The simplest why to implement a new matrix class
 * that conforms with these mat_traits requirements is to provide all those
 * nested typedefs and static constants. The other alternative, which is non-intrusive,
 * is to define a specialization for mat_traits for the new matrix class, providing
 * all the public members as seen in this general trait template.
 */
template <typename Matrix>
struct mat_traits {
  /// The type of the elements of the matrix.
  using value_type = typename Matrix::value_type;
  /// The type for a reference to an element of the matrix
  using reference = typename Matrix::reference;
  /// The type for a const reference to an element of the matrix.
  using const_reference = typename Matrix::const_reference;
  /// The type for a pointer to an element of the matrix.
  using pointer = typename Matrix::pointer;
  /// The type for a const pointer to an element of the matrix.
  using const_pointer = typename Matrix::const_pointer;

  /// The type of the row-iterator for the matrix (a row-iterator goes from one row to another (on the same column)).
  using row_iterator = typename Matrix::row_iterator;
  /// The type of the const row-iterator for the matrix (a row-iterator goes from one row to another (on the same
  /// column)).
  using const_row_iterator = typename Matrix::const_row_iterator;
  /// The type of the column-iterator for the matrix (a column-iterator goes from one column to another (on the same
  /// row)).
  using col_iterator = typename Matrix::col_iterator;
  /// The type of the const column-iterator for the matrix (a column-iterator goes from one column to another (on the
  /// same row)).
  using const_col_iterator = typename Matrix::const_col_iterator;

  /// The type of the size descriptors (or index descriptors) of the matrix class.
  using size_type = typename Matrix::size_type;
  /// The type of the difference between two size (or index) descriptors of the matrix class.
  using difference_type = typename Matrix::difference_type;

  /// The static row count. Should be 0 if the row count is dynamic.
  static constexpr unsigned int static_row_count = Matrix::static_row_count;
  /// The static column count. Should be 0 if the column count is dynamic.
  static constexpr unsigned int static_col_count = Matrix::static_col_count;
  /// The alignment tag.
  static constexpr mat_alignment::tag alignment = Matrix::alignment;
  /// The structure tag.
  static constexpr mat_structure::tag structure = Matrix::structure;
};

template <typename Matrix, typename = void>
struct mat_value_type {
  using type = double;
};

// Only get value type without relying on mat_traits.
template <typename Matrix>
struct mat_value_type<Matrix,
                      std::void_t<decltype(std::declval<Matrix>()(0, 0))>> {
  using type =
      std::decay_t<decltype(std::declval<std::add_const_t<Matrix>>()(0, 0))>;
};

template <typename Matrix>
using mat_value_type_t = typename mat_value_type<Matrix>::type;

template <typename Matrix, typename = void>
struct mat_size_type {
  using type = int;
};

// Only get size type without relying on mat_traits.
template <typename Matrix>
struct mat_size_type<
    Matrix, std::void_t<decltype(std::declval<Matrix>().get_col_count())>> {
  using type = std::decay_t<
      decltype(std::declval<std::add_const_t<Matrix>>().get_col_count())>;
};

template <typename Matrix>
using mat_size_type_t = typename mat_size_type<Matrix>::type;

/*
 * This detail section includes several template specializations which are used
 * to define the priority of certain matrix structures when it comes to doing
 * additions (or subtractions) and multiplications. These meta-functions are
 * used to perform Sfinae-based switching between overloaded versions of operators
 * and functions based on the types. For example, when multiplying a symmetric
 * matrix with a diagonal matrix, surely, the general matrix multiplication operator
 * should be turned off by sfinae, additionally, the diagonal matrix multiplication
 * operator should be selected over the symmetric matrix multiplication operator. These
 * meta-functions allow this kind of situations to be detected and the correct overload
 * to be activated (by turning off the others (via Sfinae) which have lower priority).
 *
 * Please ignore this section unless interested in the inner-workings of Sfinae-based
 * switching of operator overload selection for the matrix classes.
 */
namespace detail {

template <mat_structure::tag Structure>
struct product_priority {
  using value_type = int;
  static constexpr int value = 0;
  using type = product_priority<Structure>;
};

template <mat_structure::tag Structure>
static constexpr int product_priority_v =
    product_priority<Structure>::value;

template <>
struct product_priority<mat_structure::rectangular> {
  using value_type = int;
  static constexpr int value = 1;
  using type = product_priority<mat_structure::rectangular>;
};

template <>
struct product_priority<mat_structure::square> {
  using value_type = int;
  static constexpr int value = 2;
  using type = product_priority<mat_structure::square>;
};

template <>
struct product_priority<mat_structure::symmetric> {
  using value_type = int;
  static constexpr int value = 10;
  using type = product_priority<mat_structure::symmetric>;
};

template <>
struct product_priority<mat_structure::skew_symmetric> {
  using value_type = int;
  static constexpr int value = 11;
  using type = product_priority<mat_structure::skew_symmetric>;
};

template <>
struct product_priority<mat_structure::diagonal> {
  using value_type = int;
  static constexpr int value = 40;
  using type = product_priority<mat_structure::diagonal>;
};

template <>
struct product_priority<mat_structure::scalar> {
  using value_type = int;
  static constexpr int value = 41;
  using type = product_priority<mat_structure::scalar>;
};

template <>
struct product_priority<mat_structure::upper_triangular> {
  using value_type = int;
  static constexpr int value = 20;
  using type = product_priority<mat_structure::upper_triangular>;
};

template <>
struct product_priority<mat_structure::lower_triangular> {
  using value_type = int;
  static constexpr int value = 21;
  using type = product_priority<mat_structure::lower_triangular>;
};

template <>
struct product_priority<mat_structure::orthogonal> {
  using value_type = int;
  static constexpr int value = 3;
  using type = product_priority<mat_structure::orthogonal>;
};

template <>
struct product_priority<mat_structure::tridiagonal> {
  using value_type = int;
  static constexpr int value = 30;
  using type = product_priority<mat_structure::tridiagonal>;
};

template <>
struct product_priority<mat_structure::nil> {
  using value_type = int;
  static constexpr int value = 50;
  using type = product_priority<mat_structure::nil>;
};

template <>
struct product_priority<mat_structure::identity> {
  using value_type = int;
  static constexpr int value = 49;
  using type = product_priority<mat_structure::identity>;
};

template <>
struct product_priority<mat_structure::permutation> {
  using value_type = int;
  static constexpr int value = 48;
  using type = product_priority<mat_structure::permutation>;
};

template <mat_structure::tag Structure>
struct addition_priority {
  using value_type = int;
  static constexpr int value = 0;
  using type = addition_priority<Structure>;
};

template <mat_structure::tag Structure>
static constexpr int addition_priority_v =
    addition_priority<Structure>::value;

template <>
struct addition_priority<mat_structure::rectangular> {
  using value_type = int;
  static constexpr int value = 1;
  using type = addition_priority<mat_structure::rectangular>;
};

template <>
struct addition_priority<mat_structure::square> {
  using value_type = int;
  static constexpr int value = 2;
  using type = addition_priority<mat_structure::square>;
};

template <>
struct addition_priority<mat_structure::symmetric> {
  using value_type = int;
  static constexpr int value = 10;
  using type = addition_priority<mat_structure::symmetric>;
};

template <>
struct addition_priority<mat_structure::skew_symmetric> {
  using value_type = int;
  static constexpr int value = 11;
  using type = addition_priority<mat_structure::skew_symmetric>;
};

template <>
struct addition_priority<mat_structure::diagonal> {
  using value_type = int;
  static constexpr int value = 35;
  using type = addition_priority<mat_structure::diagonal>;
};

template <>
struct addition_priority<mat_structure::scalar> {
  using value_type = int;
  static constexpr int value = 40;
  using type = addition_priority<mat_structure::scalar>;
};

template <>
struct addition_priority<mat_structure::identity> {
  using value_type = int;
  static constexpr int value = 40;
  using type = addition_priority<mat_structure::identity>;
};

template <>
struct addition_priority<mat_structure::permutation> {
  using value_type = int;
  static constexpr int value = 3;
  using type = addition_priority<mat_structure::permutation>;
};

template <>
struct addition_priority<mat_structure::upper_triangular> {
  using value_type = int;
  static constexpr int value = 20;
  using type = addition_priority<mat_structure::upper_triangular>;
};

template <>
struct addition_priority<mat_structure::lower_triangular> {
  using value_type = int;
  static constexpr int value = 21;
  using type = addition_priority<mat_structure::lower_triangular>;
};

template <>
struct addition_priority<mat_structure::orthogonal> {
  using value_type = int;
  static constexpr int value = 3;
  using type = addition_priority<mat_structure::orthogonal>;
};

template <>
struct addition_priority<mat_structure::tridiagonal> {
  using value_type = int;
  static constexpr int value = 30;
  using type = addition_priority<mat_structure::tridiagonal>;
};

template <>
struct addition_priority<mat_structure::nil> {
  using value_type = int;
  static constexpr int value = 50;
  using type = addition_priority<mat_structure::nil>;
};

}  // namespace detail

// Specialized in mat_alg_general.
template <typename Matrix>
struct mat_product_priority {
  using value_type = int;
  static constexpr int value = 0;
  using type = mat_product_priority<Matrix>;
};

template <typename Matrix>
static constexpr int mat_product_priority_v =
    mat_product_priority<Matrix>::value;

// Specialized in mat_alg_general.
template <typename Matrix>
struct mat_addition_priority {
  using value_type = int;
  static constexpr int value = 0;
  using type = mat_addition_priority<Matrix>;
};

template <typename Matrix>
static constexpr int mat_addition_priority_v =
    mat_addition_priority<Matrix>::value;

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_TRAITS_H_

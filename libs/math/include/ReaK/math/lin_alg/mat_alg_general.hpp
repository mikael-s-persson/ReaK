/**
 * \file mat_alg_general_hpp
 *
 * This library implements the general versions of many meta-functions (templates),
 * functions, and operators. These are meant to be used when no more-specialized
 * implementations exist for the matrix types involved.
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

#ifndef REAK_MAT_ALG_GENERAL_HPP
#define REAK_MAT_ALG_GENERAL_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/serializable.hpp>
#include <ReaK/core/rtti/so_register_type.hpp>

#include "mat_composite_adaptor.hpp"
#include "mat_concepts.hpp"
#include "mat_slices.hpp"
#include "mat_traits.hpp"
#include "mat_vector_adaptor.hpp"
#include "mat_views.hpp"
#include "stride_iterator.hpp"
#include "vect_alg.hpp"
#include "vect_concepts.hpp"

#include <type_traits>

#include <boost/concept_check.hpp>

namespace ReaK {

/**
 * This class is the general class template for all matrix classes in the ReaK linear algebra
 * libraries. The general template itself should never be used and will cause a compilation
 * error if it is, all useful matrix class templates are, in fact, partial specializations of
 * this general class template.
 *
 * Models: ReadableMatrixConcept.
 *
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Structure Enum which defines the structure of the matrix, see mat_structure::tag.
 * \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or
 *mat_alignment::column_major (default).
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T, mat_structure::tag Structure = mat_structure::rectangular,
          mat_alignment::tag Alignment = mat_alignment::column_major,
          typename Allocator = std::allocator<T>>
class mat {
  char this_specialization_is_not_available_or_possible[0];  // NOLINT
};

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, typename Allocator>
struct is_readable_matrix<mat<T, Structure, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type = is_readable_matrix<mat<T, Structure, Alignment, Allocator>>;
};

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, typename Allocator>
struct is_writable_matrix<mat<T, Structure, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type = is_writable_matrix<mat<T, Structure, Alignment, Allocator>>;
};

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, typename Allocator>
struct is_resizable_matrix<mat<T, Structure, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type = is_resizable_matrix<mat<T, Structure, Alignment, Allocator>>;
};

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, typename Allocator>
struct has_allocator_matrix<mat<T, Structure, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type = has_allocator_matrix<mat<T, Structure, Alignment, Allocator>>;
};

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, typename Allocator>
struct mat_product_priority<mat<T, Structure, Alignment, Allocator>> {
  using value_type = std::size_t;
  static constexpr std::size_t value =
      detail::product_priority<Structure>::value;
  using type = detail::product_priority<Structure>;
};

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, typename Allocator>
struct mat_addition_priority<mat<T, Structure, Alignment, Allocator>> {
  using value_type = std::size_t;
  static constexpr std::size_t value =
      detail::addition_priority<Structure>::value;
  using type = detail::addition_priority<Structure>;
};

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, typename Allocator>
struct is_square_matrix<mat<T, Structure, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = ((Structure != mat_structure::rectangular &&
                                  (Structure != mat_structure::nil)));
  using type = is_square_matrix<mat<T, Structure, Alignment, Allocator>>;
};

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, typename Allocator>
struct is_symmetric_matrix<mat<T, Structure, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = ((Structure == mat_structure::symmetric) ||
                                 (Structure == mat_structure::diagonal) ||
                                 (Structure == mat_structure::tridiagonal) ||
                                 (Structure == mat_structure::identity));
  using type = is_symmetric_matrix<mat<T, Structure, Alignment, Allocator>>;
};

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, typename Allocator>
struct is_diagonal_matrix<mat<T, Structure, Alignment, Allocator>> {
  using value_type = bool;
  static constexpr bool value = ((Structure == mat_structure::diagonal ||
                                  (Structure == mat_structure::identity)));
  using type = is_diagonal_matrix<mat<T, Structure, Alignment, Allocator>>;
};

template <typename T, mat_structure::tag Structure = mat_structure::rectangular,
          unsigned int RowCount = 1, unsigned int ColCount = RowCount,
          mat_alignment::tag Alignment = mat_alignment::column_major>
class mat_fix {
  char this_specialization_is_not_available_or_possible[0];  // NOLINT
};

template <typename T, mat_structure::tag Structure, unsigned int RowCount,
          unsigned int ColCount, mat_alignment::tag Alignment>
struct is_readable_matrix<
    mat_fix<T, Structure, RowCount, ColCount, Alignment>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type =
      is_readable_matrix<mat_fix<T, Structure, RowCount, ColCount, Alignment>>;
};

template <typename T, mat_structure::tag Structure, unsigned int RowCount,
          unsigned int ColCount, mat_alignment::tag Alignment>
struct is_writable_matrix<
    mat_fix<T, Structure, RowCount, ColCount, Alignment>> {
  using value_type = bool;
  static constexpr bool value = true;
  using type =
      is_writable_matrix<mat_fix<T, Structure, RowCount, ColCount, Alignment>>;
};

template <typename T, mat_structure::tag Structure, unsigned int RowCount,
          unsigned int ColCount, mat_alignment::tag Alignment>
struct is_resizable_matrix<
    mat_fix<T, Structure, RowCount, ColCount, Alignment>> {
  using value_type = bool;
  static constexpr bool value = false;
  using type =
      is_resizable_matrix<mat_fix<T, Structure, RowCount, ColCount, Alignment>>;
};

template <typename T, mat_structure::tag Structure, unsigned int RowCount,
          unsigned int ColCount, mat_alignment::tag Alignment>
struct has_allocator_matrix<
    mat_fix<T, Structure, RowCount, ColCount, Alignment>> {
  using value_type = bool;
  static constexpr bool value = false;
  using type = has_allocator_matrix<
      mat_fix<T, Structure, RowCount, ColCount, Alignment>>;
};

template <typename T, mat_structure::tag Structure, unsigned int RowCount,
          unsigned int ColCount, mat_alignment::tag Alignment>
struct mat_product_priority<
    mat_fix<T, Structure, RowCount, ColCount, Alignment>> {
  using value_type = std::size_t;
  static constexpr std::size_t value =
      detail::product_priority<Structure>::value;
  using type = detail::product_priority<Structure>;
};

template <typename T, mat_structure::tag Structure, unsigned int RowCount,
          unsigned int ColCount, mat_alignment::tag Alignment>
struct mat_addition_priority<
    mat_fix<T, Structure, RowCount, ColCount, Alignment>> {
  using value_type = std::size_t;
  static constexpr std::size_t value =
      detail::addition_priority<Structure>::value;
  using type = detail::addition_priority<Structure>;
};

template <typename T, mat_structure::tag Structure, unsigned int RowCount,
          unsigned int ColCount, mat_alignment::tag Alignment>
struct is_square_matrix<mat_fix<T, Structure, RowCount, ColCount, Alignment>> {
  using value_type = bool;
  static constexpr bool value = (RowCount == ColCount);
  using type =
      is_square_matrix<mat_fix<T, Structure, RowCount, ColCount, Alignment>>;
};

template <typename T, mat_structure::tag Structure, unsigned int RowCount,
          unsigned int ColCount, mat_alignment::tag Alignment>
struct is_symmetric_matrix<
    mat_fix<T, Structure, RowCount, ColCount, Alignment>> {
  using value_type = bool;
  static constexpr bool value = ((Structure == mat_structure::symmetric) ||
                                 (Structure == mat_structure::diagonal) ||
                                 (Structure == mat_structure::tridiagonal) ||
                                 (Structure == mat_structure::identity));
  using type =
      is_symmetric_matrix<mat_fix<T, Structure, RowCount, ColCount, Alignment>>;
};

template <typename T, mat_structure::tag Structure, unsigned int RowCount,
          unsigned int ColCount, mat_alignment::tag Alignment>
struct is_diagonal_matrix<
    mat_fix<T, Structure, RowCount, ColCount, Alignment>> {
  using value_type = bool;
  static constexpr bool value = ((Structure == mat_structure::diagonal ||
                                  (Structure == mat_structure::identity)));
  using type =
      is_diagonal_matrix<mat_fix<T, Structure, RowCount, ColCount, Alignment>>;
};

template <mat_structure::tag Structure, mat_alignment::tag Alignment>
struct mat_indexer {};

namespace rtti {

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, typename Allocator>
struct get_type_id<mat<T, Structure, Alignment, Allocator>> {
  static constexpr unsigned int ID = 0x00000012;
  static constexpr auto type_name = std::string_view{"mat"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }

  using save_type = const serializable&;
  using load_type = serializable&;
};

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, typename Allocator, typename Tail>
struct get_type_info<mat<T, Structure, Alignment, Allocator>, Tail> {
  using type =
      type_id<mat<T, Structure, Alignment, Allocator>,
              typename get_type_info_seq<
                  T, std::integral_constant<mat_structure::tag, Structure>,
                  std::integral_constant<mat_alignment::tag, Alignment>>::
                  template with_tail<Tail>::type::type>;
  static constexpr auto type_name = ct_concat_v<
      get_type_id<mat<T, Structure, Alignment, Allocator>>::type_name,
      lsl_left_bracket,
      get_type_info_seq<
          T, std::integral_constant<mat_structure::tag, Structure>,
          std::integral_constant<mat_alignment::tag, Alignment>>::type_name,
      lsl_right_bracket, get_type_name_tail<Tail>::value>;
};

template <typename T, mat_structure::tag Structure, unsigned int RowCount,
          unsigned int ColCount, mat_alignment::tag Alignment>
struct get_type_id<mat_fix<T, Structure, RowCount, ColCount, Alignment>> {
  static constexpr unsigned int ID = 0x00000013;
  static constexpr auto type_name = std::string_view{"mat_fix"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }

  using save_type = const serializable&;
  using load_type = serializable&;
};

template <typename T, mat_structure::tag Structure, unsigned int RowCount,
          unsigned int ColCount, mat_alignment::tag Alignment, typename Tail>
struct get_type_info<mat_fix<T, Structure, RowCount, ColCount, Alignment>,
                     Tail> {
  using type =
      type_id<mat_fix<T, Structure, RowCount, ColCount, Alignment>,
              typename get_type_info_seq<
                  T, std::integral_constant<mat_structure::tag, Structure>,
                  std::integral_constant<unsigned int, RowCount>,
                  std::integral_constant<unsigned int, ColCount>,
                  std::integral_constant<mat_alignment::tag, Alignment>>::
                  template with_tail<Tail>::type::type>;
  static constexpr auto type_name = ct_concat_v<
      get_type_id<
          mat_fix<T, Structure, RowCount, ColCount, Alignment>>::type_name,
      lsl_left_bracket,
      get_type_info_seq<
          T, std::integral_constant<mat_structure::tag, Structure>,
          std::integral_constant<unsigned int, RowCount>,
          std::integral_constant<unsigned int, ColCount>,
          std::integral_constant<mat_alignment::tag, Alignment>>::type_name,
      lsl_right_bracket, get_type_name_tail<Tail>::value>;
};

}  // namespace rtti

}  // namespace ReaK

#endif

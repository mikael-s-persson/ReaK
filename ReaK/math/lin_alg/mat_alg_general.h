/**
 * \file mat_alg_general.h
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

#ifndef REAK_MATH_LIN_ALG_MAT_ALG_GENERAL_H_
#define REAK_MATH_LIN_ALG_MAT_ALG_GENERAL_H_

#include "ReaK/core/rtti/so_register_type.h"
#include "ReaK/core/serialization/serializable.h"
#include "ReaK/math/lin_alg/mat_composite_adaptor.h"
#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/mat_slices.h"
#include "ReaK/math/lin_alg/mat_traits.h"
#include "ReaK/math/lin_alg/mat_vector_adaptor.h"
#include "ReaK/math/lin_alg/mat_views.h"
#include "ReaK/math/lin_alg/stride_iterator.h"
#include "ReaK/math/lin_alg/vect_alg.h"
#include "ReaK/math/lin_alg/vect_concepts.h"

#include <type_traits>

namespace ReaK {

/**
 * This class is the general class template for all matrix classes in the ReaK linear algebra
 * libraries. The general template itself should never be used and will cause a compilation
 * error if it is, all useful matrix class templates are, in fact, partial specializations of
 * this general class template.
 *
 * \tparam T Arithmetic type of the elements of the matrix.
 * \tparam Structure Enum which defines the structure of the matrix, see mat_structure::tag.
 * \tparam Alignment Enum which defines the memory alignment of the matrix. Either mat_alignment::row_major or
 *mat_alignment::column_major (default).
 * \tparam RowCount Compile-time row count (on stack), or 0 for dynamically sized (on heap).
 * \tparam ColCount Compile-time column count (on stack), or 0 for dynamically sized (on heap).
 */
template <typename T, mat_structure::tag Structure = mat_structure::rectangular,
          mat_alignment::tag Alignment = mat_alignment::column_major,
          unsigned int RowCount = 0, unsigned int ColCount = 0>
class mat {
  char this_specialization_is_not_available_or_possible[0];  // NOLINT
};

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, unsigned int RowCount,
          unsigned int ColCount>
static constexpr bool is_fully_writable_matrix_v<
    mat<T, Structure, Alignment, RowCount, ColCount>> =
    (Structure == mat_structure::rectangular) ||
    (Structure == mat_structure::square);

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, unsigned int RowCount,
          unsigned int ColCount>
struct mat_product_priority<mat<T, Structure, Alignment, RowCount, ColCount>> {
  using value_type = int;
  static constexpr int value = detail::product_priority<Structure>::value;
  using type = detail::product_priority<Structure>;
};

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, unsigned int RowCount,
          unsigned int ColCount>
struct mat_addition_priority<mat<T, Structure, Alignment, RowCount, ColCount>> {
  using value_type = int;
  static constexpr int value = detail::addition_priority<Structure>::value;
  using type = detail::addition_priority<Structure>;
};

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, unsigned int RowCount,
          unsigned int ColCount>
static constexpr bool
    is_square_matrix_v<mat<T, Structure, Alignment, RowCount, ColCount>> =
        ((Structure != mat_structure::rectangular &&
          (Structure != mat_structure::nil))) ||
        ((RowCount != 0) && (ColCount != 0) && (RowCount == ColCount));

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, unsigned int RowCount,
          unsigned int ColCount>
static constexpr bool
    is_symmetric_matrix_v<mat<T, Structure, Alignment, RowCount, ColCount>> =
        ((Structure == mat_structure::symmetric) ||
         (Structure == mat_structure::diagonal) ||
         (Structure == mat_structure::tridiagonal) ||
         (Structure == mat_structure::identity));

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, unsigned int RowCount,
          unsigned int ColCount>
static constexpr bool
    is_diagonal_matrix_v<mat<T, Structure, Alignment, RowCount, ColCount>> =
        ((Structure == mat_structure::diagonal) ||
         (Structure == mat_structure::identity));

namespace rtti {

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, unsigned int RowCount,
          unsigned int ColCount>
struct get_type_id<mat<T, Structure, Alignment, RowCount, ColCount>> {
  static constexpr std::uint32_t id = 0x00000012;
  static constexpr auto type_name = std::string_view{"mat"};
  static construct_ptr create_ptr() noexcept { return nullptr; }

  using save_type = const serializable&;
  using load_type = serializable&;
};

template <typename T, mat_structure::tag Structure,
          mat_alignment::tag Alignment, unsigned int RowCount,
          unsigned int ColCount, typename Tail>
struct get_type_info<mat<T, Structure, Alignment, RowCount, ColCount>, Tail> {
  using structure_ic = std::integral_constant<mat_structure::tag, Structure>;
  using alignment_ic = std::integral_constant<mat_alignment::tag, Alignment>;
  using row_count_ic = std::integral_constant<unsigned int, RowCount>;
  using col_count_ic = std::integral_constant<unsigned int, ColCount>;
  using arg_type_info_seq = std::conditional_t<
      (RowCount == 0) && (ColCount == 0),
      so_type_details::get_type_info_seq<T, structure_ic, alignment_ic>,
      so_type_details::get_type_info_seq<T, structure_ic, alignment_ic,
                                         row_count_ic, col_count_ic>>;
  using type = so_type_details::type_id<
      mat<T, Structure, Alignment, RowCount, ColCount>,
      typename arg_type_info_seq::template with_tail<Tail>::type::type>;
  static constexpr auto type_name = ct_concat_v<
      get_type_id<mat<T, Structure, Alignment, RowCount, ColCount>>::type_name,
      lsl_left_bracket, arg_type_info_seq::type_name, lsl_right_bracket,
      get_type_name_tail<Tail>::value>;
};

}  // namespace rtti

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_ALG_GENERAL_H_

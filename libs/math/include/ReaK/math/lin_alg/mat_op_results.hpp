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

#ifndef REAK_MAT_OP_RESULTS_HPP
#define REAK_MAT_OP_RESULTS_HPP

#include <ReaK/core/base/defs.hpp>

#include "mat_alg_general.hpp"
#include "mat_concepts.hpp"
#include "mat_traits.hpp"

#include <type_traits>

#include <boost/concept_check.hpp>

namespace ReaK {

namespace detail {

template <mat_structure::tag Structure1, mat_structure::tag Structure2>
struct product_result_structure {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::square;
  using type = product_result_structure<Structure1, Structure2>;
};

template <mat_structure::tag Structure2>
struct product_result_structure<mat_structure::rectangular, Structure2> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::rectangular;
  using type = product_result_structure<mat_structure::rectangular, Structure2>;
};

template <mat_structure::tag Structure1>
struct product_result_structure<Structure1, mat_structure::rectangular> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::rectangular;
  using type = product_result_structure<Structure1, mat_structure::rectangular>;
};

template <>
struct product_result_structure<mat_structure::rectangular,
                                mat_structure::rectangular> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::rectangular;
  using type = product_result_structure<mat_structure::rectangular,
                                        mat_structure::rectangular>;
};

template <mat_structure::tag Structure2>
struct product_result_structure<mat_structure::nil, Structure2> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::nil;
  using type = product_result_structure<mat_structure::nil, Structure2>;
};

template <mat_structure::tag Structure1>
struct product_result_structure<Structure1, mat_structure::nil> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::nil;
  using type = product_result_structure<Structure1, mat_structure::nil>;
};

template <>
struct product_result_structure<mat_structure::nil,
                                mat_structure::rectangular> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::nil;
  using type =
      product_result_structure<mat_structure::nil, mat_structure::rectangular>;
};

template <>
struct product_result_structure<mat_structure::rectangular,
                                mat_structure::nil> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::nil;
  using type =
      product_result_structure<mat_structure::rectangular, mat_structure::nil>;
};

template <>
struct product_result_structure<mat_structure::nil, mat_structure::identity> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::nil;
  using type =
      product_result_structure<mat_structure::nil, mat_structure::rectangular>;
};

template <>
struct product_result_structure<mat_structure::identity, mat_structure::nil> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::nil;
  using type =
      product_result_structure<mat_structure::rectangular, mat_structure::nil>;
};

template <>
struct product_result_structure<mat_structure::nil, mat_structure::nil> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::nil;
  using type = product_result_structure<mat_structure::nil, mat_structure::nil>;
};

template <mat_structure::tag Structure2>
struct product_result_structure<mat_structure::identity, Structure2> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = Structure2;
  using type = product_result_structure<mat_structure::identity, Structure2>;
};

template <mat_structure::tag Structure1>
struct product_result_structure<Structure1, mat_structure::identity> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = Structure1;
  using type = product_result_structure<Structure1, mat_structure::identity>;
};

template <>
struct product_result_structure<mat_structure::identity,
                                mat_structure::rectangular> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::rectangular;
  using type = product_result_structure<mat_structure::identity,
                                        mat_structure::rectangular>;
};

template <>
struct product_result_structure<mat_structure::rectangular,
                                mat_structure::identity> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::rectangular;
  using type = product_result_structure<mat_structure::rectangular,
                                        mat_structure::identity>;
};

template <>
struct product_result_structure<mat_structure::identity,
                                mat_structure::identity> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::identity;
  using type = product_result_structure<mat_structure::identity,
                                        mat_structure::identity>;
};

template <bool AreTheseMatrices, typename ResultValueType, typename Matrix1,
          typename Matrix2>
struct mat_product_result_impl {
  using type =
      mat<ResultValueType,
          detail::product_result_structure<
              mat_traits<Matrix1>::structure,
              mat_traits<Matrix2>::structure>::type::value,
          mat_traits<Matrix1>::alignment,
          typename std::allocator_traits<std::conditional_t<
              has_allocator_matrix_v<Matrix1>,
              typename mat_traits<Matrix1>::allocator_type,
              std::conditional_t<has_allocator_matrix<Matrix2>::value,
                                 typename mat_traits<Matrix2>::allocator_type,
                                 std::allocator<ResultValueType>>>>::
              template rebind_alloc<ResultValueType>>;
};

template <typename ResultValueType, typename Matrix1, typename Matrix2>
struct mat_product_result_impl<false, ResultValueType, Matrix1, Matrix2> {
  using type = mat<ResultValueType, mat_structure::rectangular>;
};

template <mat_structure::tag Structure1, mat_structure::tag Structure2>
struct addition_result_structure {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::square;
  using type = addition_result_structure<Structure1, Structure2>;
};

template <>
struct addition_result_structure<mat_structure::rectangular,
                                 mat_structure::rectangular> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::rectangular;
  using type = addition_result_structure<mat_structure::rectangular,
                                         mat_structure::rectangular>;
};

template <>
struct addition_result_structure<mat_structure::nil, mat_structure::nil> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::nil;
  using type =
      addition_result_structure<mat_structure::nil, mat_structure::nil>;
};

template <>
struct addition_result_structure<mat_structure::rectangular,
                                 mat_structure::nil> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::rectangular;
  using type =
      addition_result_structure<mat_structure::rectangular, mat_structure::nil>;
};

template <>
struct addition_result_structure<mat_structure::nil,
                                 mat_structure::rectangular> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::rectangular;
  using type =
      addition_result_structure<mat_structure::nil, mat_structure::rectangular>;
};

template <>
struct addition_result_structure<mat_structure::symmetric,
                                 mat_structure::symmetric> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::symmetric;
  using type = addition_result_structure<mat_structure::symmetric,
                                         mat_structure::symmetric>;
};

template <>
struct addition_result_structure<mat_structure::diagonal,
                                 mat_structure::symmetric> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::symmetric;
  using type = addition_result_structure<mat_structure::diagonal,
                                         mat_structure::symmetric>;
};

template <>
struct addition_result_structure<mat_structure::symmetric,
                                 mat_structure::diagonal> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::symmetric;
  using type = addition_result_structure<mat_structure::symmetric,
                                         mat_structure::diagonal>;
};

template <>
struct addition_result_structure<mat_structure::identity,
                                 mat_structure::symmetric> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::symmetric;
  using type = addition_result_structure<mat_structure::identity,
                                         mat_structure::symmetric>;
};

template <>
struct addition_result_structure<mat_structure::symmetric,
                                 mat_structure::identity> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::symmetric;
  using type = addition_result_structure<mat_structure::symmetric,
                                         mat_structure::identity>;
};

template <>
struct addition_result_structure<mat_structure::nil, mat_structure::symmetric> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::symmetric;
  using type =
      addition_result_structure<mat_structure::nil, mat_structure::symmetric>;
};

template <>
struct addition_result_structure<mat_structure::symmetric, mat_structure::nil> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::symmetric;
  using type =
      addition_result_structure<mat_structure::symmetric, mat_structure::nil>;
};

template <>
struct addition_result_structure<mat_structure::diagonal,
                                 mat_structure::diagonal> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::diagonal;
  using type = addition_result_structure<mat_structure::diagonal,
                                         mat_structure::diagonal>;
};

template <>
struct addition_result_structure<mat_structure::identity,
                                 mat_structure::diagonal> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::diagonal;
  using type = addition_result_structure<mat_structure::identity,
                                         mat_structure::diagonal>;
};

template <>
struct addition_result_structure<mat_structure::diagonal,
                                 mat_structure::identity> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::diagonal;
  using type = addition_result_structure<mat_structure::diagonal,
                                         mat_structure::identity>;
};

template <>
struct addition_result_structure<mat_structure::nil, mat_structure::identity> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::scalar;
  using type =
      addition_result_structure<mat_structure::nil, mat_structure::identity>;
};

template <>
struct addition_result_structure<mat_structure::identity, mat_structure::nil> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::scalar;
  using type =
      addition_result_structure<mat_structure::identity, mat_structure::nil>;
};

template <>
struct addition_result_structure<mat_structure::identity,
                                 mat_structure::identity> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::scalar;
  using type = addition_result_structure<mat_structure::identity,
                                         mat_structure::identity>;
};

template <>
struct addition_result_structure<mat_structure::scalar,
                                 mat_structure::identity> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::scalar;
  using type =
      addition_result_structure<mat_structure::scalar, mat_structure::identity>;
};

template <>
struct addition_result_structure<mat_structure::identity,
                                 mat_structure::scalar> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::scalar;
  using type =
      addition_result_structure<mat_structure::identity, mat_structure::scalar>;
};

template <>
struct addition_result_structure<mat_structure::scalar, mat_structure::scalar> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::scalar;
  using type =
      addition_result_structure<mat_structure::scalar, mat_structure::scalar>;
};

template <mat_structure::tag Structure2>
struct addition_result_structure<mat_structure::nil, Structure2> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = Structure2;
  using type = addition_result_structure<mat_structure::nil, Structure2>;
};

template <mat_structure::tag Structure1>
struct addition_result_structure<Structure1, mat_structure::nil> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = Structure1;
  using type = addition_result_structure<Structure1, mat_structure::nil>;
};

template <>
struct addition_result_structure<mat_structure::skew_symmetric,
                                 mat_structure::skew_symmetric> {
  using value_type = mat_structure::tag;
  static constexpr mat_structure::tag value = mat_structure::skew_symmetric;
  using type = addition_result_structure<mat_structure::skew_symmetric,
                                         mat_structure::skew_symmetric>;
};

template <bool AreTheseMatrices, typename ResultValueType, typename Matrix1,
          typename Matrix2>
struct mat_addition_result_impl {
  using type =
      mat<ResultValueType,
          detail::addition_result_structure<
              mat_traits<Matrix1>::structure,
              mat_traits<Matrix2>::structure>::type::value,
          mat_traits<Matrix1>::alignment,
          typename std::allocator_traits<std::conditional_t<
              has_allocator_matrix<Matrix1>::value,
              typename mat_traits<Matrix1>::allocator_type,
              std::conditional_t<has_allocator_matrix<Matrix2>::value,
                                 typename mat_traits<Matrix2>::allocator_type,
                                 std::allocator<ResultValueType>>>>::
              template rebind_alloc<ResultValueType>>;
};

template <typename ResultValueType, typename Matrix1, typename Matrix2>
struct mat_addition_result_impl<false, ResultValueType, Matrix1, Matrix2> {
  using type = mat<ResultValueType, mat_structure::rectangular>;
};
};  // namespace detail

template <typename Matrix1, typename Matrix2>
struct mat_product_result {
  using M1ValueType = typename std::conditional<
      is_readable_matrix_v<Matrix1>, mat_traits<Matrix1>,
      mat_traits<mat<double, mat_structure::rectangular>>>::type::value_type;
  using M2ValueType = typename std::conditional<
      is_readable_matrix_v<Matrix2>, mat_traits<Matrix2>,
      mat_traits<mat<double, mat_structure::rectangular>>>::type::value_type;
  using ResultValueType =
      decltype(std::declval<M1ValueType>() * std::declval<M2ValueType>() +
               std::declval<M1ValueType>() * std::declval<M2ValueType>());
  using type = typename detail::mat_product_result_impl<
      is_readable_matrix_v<Matrix1> && is_readable_matrix_v<Matrix2>,
      ResultValueType, Matrix1, Matrix2>::type;
};

template <typename Matrix1, typename Matrix2>
using mat_product_result_t =
    typename mat_product_result<Matrix1, Matrix2>::type;

template <typename Matrix1, typename Matrix2>
struct mat_addition_result {
  using M1ValueType = typename std::conditional<
      is_readable_matrix_v<Matrix1>, mat_traits<Matrix1>,
      mat_traits<mat<double, mat_structure::rectangular>>>::type::value_type;
  using M2ValueType = typename std::conditional<
      is_readable_matrix_v<Matrix2>, mat_traits<Matrix2>,
      mat_traits<mat<double, mat_structure::rectangular>>>::type::value_type;
  using ResultValueType =
      decltype(std::declval<M1ValueType>() + std::declval<M2ValueType>());
  using type = typename detail::mat_addition_result_impl<
      is_readable_matrix_v<Matrix1> && is_readable_matrix_v<Matrix2>,
      ResultValueType, Matrix1, Matrix2>::type;
};

template <typename Matrix1, typename Matrix2>
using mat_addition_result_t =
    typename mat_addition_result<Matrix1, Matrix2>::type;

};  // namespace ReaK

#endif

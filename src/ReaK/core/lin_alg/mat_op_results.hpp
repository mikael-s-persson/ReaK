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

#include "mat_concepts.hpp"
#include "mat_traits.hpp"
#include "mat_alg_general.hpp"

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/concept_check.hpp>

#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/comparison.hpp>


namespace ReaK {


namespace detail {

template <mat_structure::tag Structure1, mat_structure::tag Structure2>
struct product_result_structure {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::square);
  typedef product_result_structure<Structure1,Structure2> type;
};

template <mat_structure::tag Structure2>
struct product_result_structure<mat_structure::rectangular,Structure2> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::rectangular);
  typedef product_result_structure<mat_structure::rectangular,Structure2> type;
};

template <mat_structure::tag Structure1>
struct product_result_structure<Structure1,mat_structure::rectangular> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::rectangular);
  typedef product_result_structure<Structure1,mat_structure::rectangular> type;
};

template <>
struct product_result_structure<mat_structure::rectangular,mat_structure::rectangular> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::rectangular);
  typedef product_result_structure<mat_structure::rectangular,mat_structure::rectangular> type;
};

template <mat_structure::tag Structure2>
struct product_result_structure<mat_structure::nil,Structure2> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::nil);
  typedef product_result_structure<mat_structure::nil,Structure2> type;
};

template <mat_structure::tag Structure1>
struct product_result_structure<Structure1,mat_structure::nil> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::nil);
  typedef product_result_structure<Structure1,mat_structure::nil> type;
};

template <>
struct product_result_structure<mat_structure::nil,mat_structure::rectangular> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::nil);
  typedef product_result_structure<mat_structure::nil,mat_structure::rectangular> type;
};

template <>
struct product_result_structure<mat_structure::rectangular,mat_structure::nil> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::nil);
  typedef product_result_structure<mat_structure::rectangular,mat_structure::nil> type;
};

template <>
struct product_result_structure<mat_structure::nil,mat_structure::identity> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::nil);
  typedef product_result_structure<mat_structure::nil,mat_structure::rectangular> type;
};

template <>
struct product_result_structure<mat_structure::identity,mat_structure::nil> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::nil);
  typedef product_result_structure<mat_structure::rectangular,mat_structure::nil> type;
};

template <>
struct product_result_structure<mat_structure::nil,mat_structure::nil> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::nil);
  typedef product_result_structure<mat_structure::nil,mat_structure::nil> type;
};

template <mat_structure::tag Structure2>
struct product_result_structure<mat_structure::identity,Structure2> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = Structure2);
  typedef product_result_structure<mat_structure::identity,Structure2> type;
};

template <mat_structure::tag Structure1>
struct product_result_structure<Structure1,mat_structure::identity> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = Structure1);
  typedef product_result_structure<Structure1,mat_structure::identity> type;
};

template <>
struct product_result_structure<mat_structure::identity,mat_structure::rectangular> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::rectangular);
  typedef product_result_structure<mat_structure::identity,mat_structure::rectangular> type;
};

template <>
struct product_result_structure<mat_structure::rectangular,mat_structure::identity> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::rectangular);
  typedef product_result_structure<mat_structure::rectangular,mat_structure::identity> type;
};

template <>
struct product_result_structure<mat_structure::identity,mat_structure::identity> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::identity);
  typedef product_result_structure<mat_structure::identity,mat_structure::identity> type;
};


template <bool AreTheseMatrices, typename ResultValueType, typename Matrix1, typename Matrix2>
struct mat_product_result_impl {
  typedef mat< 
    ResultValueType,
    detail::product_result_structure< 
      mat_traits<Matrix1>::structure, 
          mat_traits<Matrix2>::structure 
        >::type::value,
    mat_traits<Matrix1>::alignment,
        typename boost::mpl::if_< 
          has_allocator_matrix< Matrix1 >,
          typename mat_traits< Matrix1 >::allocator_type,
          typename boost::mpl::if_< 
            has_allocator_matrix< Matrix2 >,
            typename mat_traits< Matrix2 >::allocator_type,
            std::allocator<ResultValueType> 
          >::type 
        >::type::template rebind<ResultValueType>::other 
  > type;
};

template <typename ResultValueType, typename Matrix1, typename Matrix2>
struct mat_product_result_impl<false, ResultValueType, Matrix1, Matrix2> {
  typedef mat<ResultValueType, mat_structure::rectangular> type;
};



template <mat_structure::tag Structure1, mat_structure::tag Structure2>
struct addition_result_structure {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::square);
  typedef addition_result_structure<Structure1,Structure2> type;
};

template <>
struct addition_result_structure<mat_structure::rectangular, mat_structure::rectangular> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::rectangular);
  typedef addition_result_structure<mat_structure::rectangular,mat_structure::rectangular> type;
};

template <>
struct addition_result_structure<mat_structure::nil, mat_structure::nil> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::nil);
  typedef addition_result_structure<mat_structure::nil,mat_structure::nil> type;
};

template <>
struct addition_result_structure<mat_structure::rectangular, mat_structure::nil> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::rectangular);
  typedef addition_result_structure<mat_structure::rectangular,mat_structure::nil> type;
};

template <>
struct addition_result_structure<mat_structure::nil, mat_structure::rectangular> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::rectangular);
  typedef addition_result_structure<mat_structure::nil,mat_structure::rectangular> type;
};

template <>
struct addition_result_structure<mat_structure::symmetric, mat_structure::symmetric> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::symmetric);
  typedef addition_result_structure<mat_structure::symmetric,mat_structure::symmetric> type;
};

template <>
struct addition_result_structure<mat_structure::diagonal, mat_structure::symmetric> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::symmetric);
  typedef addition_result_structure<mat_structure::diagonal,mat_structure::symmetric> type;
};

template <>
struct addition_result_structure<mat_structure::symmetric, mat_structure::diagonal> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::symmetric);
  typedef addition_result_structure<mat_structure::symmetric,mat_structure::diagonal> type;
};

template <>
struct addition_result_structure<mat_structure::identity, mat_structure::symmetric> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::symmetric);
  typedef addition_result_structure<mat_structure::identity,mat_structure::symmetric> type;
};

template <>
struct addition_result_structure<mat_structure::symmetric, mat_structure::identity> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::symmetric);
  typedef addition_result_structure<mat_structure::symmetric,mat_structure::identity> type;
};

template <>
struct addition_result_structure<mat_structure::nil, mat_structure::symmetric> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::symmetric);
  typedef addition_result_structure<mat_structure::nil,mat_structure::symmetric> type;
};

template <>
struct addition_result_structure<mat_structure::symmetric, mat_structure::nil> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::symmetric);
  typedef addition_result_structure<mat_structure::symmetric,mat_structure::nil> type;
};

template <>
struct addition_result_structure<mat_structure::diagonal, mat_structure::diagonal> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::diagonal);
  typedef addition_result_structure<mat_structure::diagonal,mat_structure::diagonal> type;
};

template <>
struct addition_result_structure<mat_structure::identity, mat_structure::diagonal> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::diagonal);
  typedef addition_result_structure<mat_structure::identity,mat_structure::diagonal> type;
};

template <>
struct addition_result_structure<mat_structure::diagonal, mat_structure::identity> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::diagonal);
  typedef addition_result_structure<mat_structure::diagonal,mat_structure::identity> type;
};

template <>
struct addition_result_structure<mat_structure::nil, mat_structure::identity> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::scalar);
  typedef addition_result_structure<mat_structure::nil,mat_structure::identity> type;
};

template <>
struct addition_result_structure<mat_structure::identity, mat_structure::nil> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::scalar);
  typedef addition_result_structure<mat_structure::identity,mat_structure::nil> type;
};

template <>
struct addition_result_structure<mat_structure::identity, mat_structure::identity> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::scalar);
  typedef addition_result_structure<mat_structure::identity,mat_structure::identity> type;
};

template <>
struct addition_result_structure<mat_structure::scalar, mat_structure::identity> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::scalar);
  typedef addition_result_structure<mat_structure::scalar,mat_structure::identity> type;
};

template <>
struct addition_result_structure<mat_structure::identity, mat_structure::scalar> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::scalar);
  typedef addition_result_structure<mat_structure::identity,mat_structure::scalar> type;
};

template <>
struct addition_result_structure<mat_structure::scalar, mat_structure::scalar> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::scalar);
  typedef addition_result_structure<mat_structure::scalar,mat_structure::scalar> type;
};

template <mat_structure::tag Structure2>
struct addition_result_structure<mat_structure::nil, Structure2> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = Structure2);
  typedef addition_result_structure<mat_structure::nil,Structure2> type;
};

template <mat_structure::tag Structure1>
struct addition_result_structure<Structure1, mat_structure::nil> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = Structure1);
  typedef addition_result_structure<Structure1,mat_structure::nil> type;
};

template <>
struct addition_result_structure<mat_structure::skew_symmetric, mat_structure::skew_symmetric> {
  typedef boost::mpl::integral_c_tag tag;
  typedef mat_structure::tag value_type;
  BOOST_STATIC_CONSTANT(mat_structure::tag, value = mat_structure::skew_symmetric);
  typedef addition_result_structure<mat_structure::skew_symmetric,mat_structure::skew_symmetric> type;
};




template <bool AreTheseMatrices, typename ResultValueType, typename Matrix1, typename Matrix2>
struct mat_addition_result_impl {
  typedef mat< 
    ResultValueType,
    detail::addition_result_structure< 
      mat_traits<Matrix1>::structure, 
      mat_traits<Matrix2>::structure 
        >::type::value,
        mat_traits<Matrix1>::alignment,
        typename boost::mpl::if_< 
          has_allocator_matrix< Matrix1 >,
          typename mat_traits< Matrix1 >::allocator_type,
          typename boost::mpl::if_< 
            has_allocator_matrix< Matrix2 >,
            typename mat_traits< Matrix2 >::allocator_type,
            std::allocator<ResultValueType> 
          >::type 
    >::type::template rebind<ResultValueType>::other
  > type;
};

template <typename ResultValueType, typename Matrix1, typename Matrix2>
struct mat_addition_result_impl<false, ResultValueType, Matrix1, Matrix2> {
  typedef mat<ResultValueType, mat_structure::rectangular> type;
};



};



template <typename Matrix1, typename Matrix2>
struct mat_product_result {
  typedef typename boost::mpl::if_<
        is_readable_matrix<Matrix1>,
        mat_traits<Matrix1>,
        mat_traits< mat<double, mat_structure::rectangular> > >::type::value_type M1ValueType;
  typedef typename boost::mpl::if_<
        is_readable_matrix<Matrix2>,
        mat_traits<Matrix2>,
        mat_traits< mat<double, mat_structure::rectangular> > >::type::value_type M2ValueType;
#ifndef BOOST_NO_CXX11_DECLTYPE
  typedef decltype( M1ValueType() * M2ValueType() + M1ValueType() * M2ValueType() ) ResultValueType;
#else
  typedef M1ValueType ResultValueType;
#endif
  typedef typename detail::mat_product_result_impl<
    boost::mpl::and_< is_readable_matrix<Matrix1>, is_readable_matrix<Matrix2> >::type::value,
        ResultValueType, Matrix1, Matrix2>::type type;
};


template <typename Matrix1, typename Matrix2>
struct mat_addition_result {
  typedef typename boost::mpl::if_<
        is_readable_matrix<Matrix1>,
        mat_traits<Matrix1>,
        mat_traits< mat<double, mat_structure::rectangular> > >::type::value_type M1ValueType;
  typedef typename boost::mpl::if_<
        is_readable_matrix<Matrix2>,
        mat_traits<Matrix2>,
        mat_traits< mat<double, mat_structure::rectangular> > >::type::value_type M2ValueType;
#ifndef BOOST_NO_CXX11_DECLTYPE
  typedef decltype( M1ValueType() + M2ValueType() ) ResultValueType;
#else
  typedef M1ValueType ResultValueType;
#endif
  typedef typename detail::mat_addition_result_impl<
    boost::mpl::and_< is_readable_matrix<Matrix1>, is_readable_matrix<Matrix2> >::type::value,
        ResultValueType, Matrix1, Matrix2>::type type;
};





};

#endif














/**
 * \file tensor_mat_adaptor.hpp
 * 
 * This header declares the necessary type-traits, meta-functions and function overloads to make a 
 * matrix class (as in ReadableMatrixConcept) usable as a second-order tensor.
 * 
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2012
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

#ifndef REAK_TENSOR_MAT_ADAPTOR_HPP
#define REAK_TENSOR_MAT_ADAPTOR_HPP

#include <ReaK/core/base/defs.hpp>

#include "tensor_traits.hpp"
#include "tensor_concepts.hpp"

#include <ReaK/core/rtti/rtti.hpp>

namespace ReaK {
  
  
namespace detail {
  
  
  template <mat_structure::tag Structure>
  struct translate_to_tensor_structure {
    BOOST_STATIC_CONSTANT(tensor_structure::tag, value = tensor_structure::rectangular);
  };
  
  template <>
  struct translate_to_tensor_structure<mat_structure::square> {
    BOOST_STATIC_CONSTANT(tensor_structure::tag, value = tensor_structure::square);
  };
  
  template <>
  struct translate_to_tensor_structure<mat_structure::symmetric> {
    BOOST_STATIC_CONSTANT(tensor_structure::tag, value = tensor_structure::square);
  };
  
  template <>
  struct translate_to_tensor_structure<mat_structure::skew_symmetric> {
    BOOST_STATIC_CONSTANT(tensor_structure::tag, value = tensor_structure::square);
  };
  
  template <>
  struct translate_to_tensor_structure<mat_structure::diagonal> {
    BOOST_STATIC_CONSTANT(tensor_structure::tag, value = tensor_structure::diagonal);
  };
  
  template <>
  struct translate_to_tensor_structure<mat_structure::upper_triangular> {
    BOOST_STATIC_CONSTANT(tensor_structure::tag, value = tensor_structure::square);
  };
  
  template <>
  struct translate_to_tensor_structure<mat_structure::lower_triangular> {
    BOOST_STATIC_CONSTANT(tensor_structure::tag, value = tensor_structure::square);
  };
  
  template <>
  struct translate_to_tensor_structure<mat_structure::orthogonal> {
    BOOST_STATIC_CONSTANT(tensor_structure::tag, value = tensor_structure::square);
  };
  
  template <>
  struct translate_to_tensor_structure<mat_structure::tridiagonal> {
    BOOST_STATIC_CONSTANT(tensor_structure::tag, value = tensor_structure::square);
  };
  
  template <>
  struct translate_to_tensor_structure<mat_structure::nil> {
    BOOST_STATIC_CONSTANT(tensor_structure::tag, value = tensor_structure::nil);
  };
  
  template <>
  struct translate_to_tensor_structure<mat_structure::identity> {
    BOOST_STATIC_CONSTANT(tensor_structure::tag, value = tensor_structure::identity);
  };
  
  template <>
  struct translate_to_tensor_structure<mat_structure::scalar> {
    BOOST_STATIC_CONSTANT(tensor_structure::tag, value = tensor_structure::scalar);
  };
  
  template <>
  struct translate_to_tensor_structure<mat_structure::permutation> {
    BOOST_STATIC_CONSTANT(tensor_structure::tag, value = tensor_structure::permutation);
  };
  
  
  
  template <mat_alignment::tag Alignment>
  struct translate_to_tensor_alignment {
    BOOST_STATIC_CONSTANT(tensor_alignment::tag, value = tensor_alignment::hi_dim_major);
  };
  
  template <>
  struct translate_to_tensor_alignment<mat_alignment::row_major> {
    BOOST_STATIC_CONSTANT(tensor_alignment::tag, value = tensor_alignment::low_dim_major);
  };
  
};


/**
 * This tensor type-traits definition is an adaptor which allows mat classes to act as 
 * order-one tensors.
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator>
struct tensor_traits< mat<T,Structure,Alignment,Allocator> > {
  typedef mat<T,Structure,Alignment,Allocator> mat_type;
  /// The type of the elements of the tensor.
  typedef typename mat_traits< mat_type >::value_type value_type;
  /// The type for a reference to an element of the tensor
  typedef typename mat_traits< mat_type >::reference reference;
  /// The type for a const reference to an element of the tensor.
  typedef typename mat_traits< mat_type >::const_reference const_reference;
  /// The type for a pointer to an element of the tensor.
  typedef typename mat_traits< mat_type >::pointer pointer;
  /// The type for a const pointer to an element of the tensor.
  typedef typename mat_traits< mat_type >::const_pointer const_pointer;
  /// The type of the allocator for the tensor (can be void if the tensor does not have an allocator).
  typedef typename mat_traits< mat_type >::allocator_type allocator_type;
  
  /// The type of the size descriptors (or index descriptors) of the tensor class.
  typedef typename mat_traits< mat_type >::size_type size_type;
  /// The type of the difference between two size (or index) descriptors of the tensor class.
  typedef typename mat_traits< mat_type >::difference_type difference_type;
  
  /// The alignment tag.
  BOOST_STATIC_CONSTANT(tensor_alignment::tag, alignment = detail::translate_to_tensor_alignment<Alignment>::value);
  /// The structure tag.
  BOOST_STATIC_CONSTANT(tensor_structure::tag, structure = detail::translate_to_tensor_structure<Structure>::value);
  /// The tensor order.
  BOOST_STATIC_CONSTANT(std::size_t, order = 2);
};

/**
 * This tensor-dimension type-traits definition is an adaptor which allows mat classes to act as 
 * order-one tensors.
 */
template <std::size_t Dim, typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator>
struct tensor_dim_traits<Dim, mat<T,Structure,Alignment,Allocator> > {
  
  /// The type of the iterator for the tensor and the given dimension.
  typedef void iterator;
  /// The type of the const iterator for the tensor and the given dimension.
  typedef void const_iterator;
  
  /// The static count in the given dimension. Should be 0 if the count is dynamic.
  BOOST_STATIC_CONSTANT(std::size_t, static_count = 1);
  
};

/**
 * This tensor-dimension type-traits definition is an adaptor which allows vect_n classes to act as 
 * order-one tensors.
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator>
struct tensor_dim_traits<0, mat<T,Structure,Alignment,Allocator> > {
  
  /// The type of the iterator for the tensor and the given dimension.
  typedef typename mat_traits< mat<T,Structure,Alignment,Allocator> >::row_iterator iterator;
  /// The type of the const iterator for the tensor and the given dimension.
  typedef typename mat_traits< mat<T,Structure,Alignment,Allocator> >::const_row_iterator const_iterator;
  
  /// The static count in the given dimension. Should be 0 if the count is dynamic.
  BOOST_STATIC_CONSTANT(std::size_t, static_count = mat_traits< mat<T,Structure,Alignment,Allocator> >::static_row_count);
  
};

/**
 * This tensor-dimension type-traits definition is an adaptor which allows vect_n classes to act as 
 * order-one tensors.
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator>
struct tensor_dim_traits<1, mat<T,Structure,Alignment,Allocator> > {
  
  /// The type of the iterator for the tensor and the given dimension.
  typedef typename mat_traits< mat<T,Structure,Alignment,Allocator> >::col_iterator iterator;
  /// The type of the const iterator for the tensor and the given dimension.
  typedef typename mat_traits< mat<T,Structure,Alignment,Allocator> >::const_col_iterator const_iterator;
  
  /// The static count in the given dimension. Should be 0 if the count is dynamic.
  BOOST_STATIC_CONSTANT(std::size_t, static_count = mat_traits< mat<T,Structure,Alignment,Allocator> >::static_col_count);
  
};



template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator>
struct is_readable_tensor< mat<T,Structure,Alignment,Allocator> > :
  is_readable_matrix< mat<T,Structure,Alignment,Allocator> > { };

template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator>
struct is_writable_tensor< mat<T,Structure,Alignment,Allocator> > :
  is_writable_matrix< mat<T,Structure,Alignment,Allocator> > { };

template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator>
struct is_fully_writable_tensor< mat<T,Structure,Alignment,Allocator> > 
  is_fully_writable_matrix< mat<T,Structure,Alignment,Allocator> > { };

template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator>
struct is_resizable_tensor< mat<T,Structure,Alignment,Allocator> > : 
  is_resizable_matrix< mat<T,Structure,Alignment,Allocator> > { };

template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator>
struct has_allocator_tensor< mat<T,Structure,Alignment,Allocator> > : 
  has_allocator_matrix< mat<T,Structure,Alignment,Allocator> > { };

  
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator>
struct is_diagonal_tensor< mat<T,Structure,Alignment,Allocator> > : 
  is_diagonal_matrix< mat<T,Structure,Alignment,Allocator> > { };
  
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator>
struct is_square_tensor< mat<T,Structure,Alignment,Allocator> > : 
  is_square_matrix< mat<T,Structure,Alignment,Allocator> > { };



template <std::size_t Dim, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_matrix<Matrix>,
    boost::mpl::equal_to< 
      boost::mpl::size_t<Dim>,
      boost::mpl::size_t<0>
    >
  >,
mat_traits< Matrix > >::size_type size(const Matrix& m) {
  return m.get_row_count();
};

template <std::size_t Dim, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_matrix<Matrix>,
    boost::mpl::equal_to< 
      boost::mpl::size_t<Dim>,
      boost::mpl::size_t<1>
    >
  >,
mat_traits< Matrix > >::size_type size(const Matrix& m) {
  return m.get_col_count();
};

template <std::size_t Dim, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_matrix<Matrix>,
    boost::mpl::greater< 
      boost::mpl::size_t<Dim>,
      boost::mpl::size_t<1>
    >
  >,
mat_traits< Matrix > >::type::size_type size(const Matrix& v) {
  return 1;
};


template <std::size_t Dim, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_resizable_matrix<Matrix>,
    boost::mpl::equal_to< 
      boost::mpl::size_t<Dim>,
      boost::mpl::size_t<0>
    >
  >,
void >::type resize(Matrix& m, std::size_t sz) {
  m.set_row_count(sz);
};

template <std::size_t Dim, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_resizable_matrix<Matrix>,
    boost::mpl::equal_to< 
      boost::mpl::size_t<Dim>,
      boost::mpl::size_t<1>
    >
  >,
void >::type resize(Matrix& m, std::size_t sz) {
  m.set_col_count(sz);
};

template <std::size_t Dim, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_resizable_matrix<Matrix>,
    boost::mpl::greater< 
      boost::mpl::size_t<Dim>,
      boost::mpl::size_t<1>
    >
  >,
void >::type resize(Matrix&, std::size_t) { };



};




#endif










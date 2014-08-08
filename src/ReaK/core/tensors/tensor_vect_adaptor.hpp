/**
 * \file tensor_vect_adaptor.hpp
 * 
 * This header declares the necessary type-traits, meta-functions and function overloads to make a 
 * vector class (as in ReadableVectorConcept) usable as a first-order tensor.
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

#ifndef REAK_TENSOR_VECT_ADAPTOR_HPP
#define REAK_TENSOR_VECT_ADAPTOR_HPP

#include <ReaK/core/base/defs.hpp>

#include "tensor_traits.hpp"
#include "tensor_concepts.hpp"

#include <ReaK/core/lin_alg/vect_alg.hpp>

#include <ReaK/core/rtti/rtti.hpp>

namespace ReaK {
  


/**
 * This tensor type-traits definition is an adaptor which allows vect_n classes to act as 
 * order-one tensors.
 */
template <typename T, typename Allocator>
struct tensor_traits< vect_n<T,Allocator> > {
  /// The type of the elements of the tensor.
  typedef typename vect_traits< vect_n<T,Allocator> >::value_type value_type;
  /// The type for a reference to an element of the tensor
  typedef typename vect_traits< vect_n<T,Allocator> >::reference reference;
  /// The type for a const reference to an element of the tensor.
  typedef typename vect_traits< vect_n<T,Allocator> >::const_reference const_reference;
  /// The type for a pointer to an element of the tensor.
  typedef typename vect_traits< vect_n<T,Allocator> >::pointer pointer;
  /// The type for a const pointer to an element of the tensor.
  typedef typename vect_traits< vect_n<T,Allocator> >::const_pointer const_pointer;
  /// The type of the allocator for the tensor (can be void if the tensor does not have an allocator).
  typedef typename vect_traits< vect_n<T,Allocator> >::allocator_type allocator_type;
  
  /// The type of the size descriptors (or index descriptors) of the tensor class.
  typedef typename vect_traits< vect_n<T,Allocator> >::size_type size_type;
  /// The type of the difference between two size (or index) descriptors of the tensor class.
  typedef typename vect_traits< vect_n<T,Allocator> >::difference_type difference_type;
  
  /// The alignment tag.
  BOOST_STATIC_CONSTANT(tensor_alignment::tag, alignment = tensor_alignment::hi_dim_major);
  /// The structure tag.
  BOOST_STATIC_CONSTANT(tensor_structure::tag, structure = tensor_structure::rectangular);
  /// The tensor order.
  BOOST_STATIC_CONSTANT(std::size_t, order = 1);
};

/**
 * This tensor-dimension type-traits definition is an adaptor which allows vect_n classes to act as 
 * order-one tensors.
 */
template <std::size_t Dim, typename T, typename Allocator>
struct tensor_dim_traits<Dim, vect_n<T,Allocator> > {
  
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
template <typename T, typename Allocator>
struct tensor_dim_traits<0, vect_n<T,Allocator> > {
  
  /// The type of the iterator for the tensor and the given dimension.
  typedef typename vect_traits< vect_n<T,Allocator> >::iterator iterator;
  /// The type of the const iterator for the tensor and the given dimension.
  typedef typename vect_traits< vect_n<T,Allocator> >::const_iterator const_iterator;
  
  /// The static count in the given dimension. Should be 0 if the count is dynamic.
  BOOST_STATIC_CONSTANT(std::size_t, static_count = vect_traits< vect_n<T,Allocator> >::dimensions);
  
};


template <typename T, typename Allocator>
struct is_readable_tensor< vect_n<T,Allocator> > :
  is_readable_vector< vect_n<T,Allocator> > { };

template <typename T, typename Allocator>
struct is_writable_tensor< vect_n<T,Allocator> > :
  is_writable_vector< vect_n<T,Allocator> > { };

template <typename T, typename Allocator>
struct is_fully_writable_tensor< vect_n<T,Allocator> > 
  is_writable_vector< vect_n<T,Allocator> > { };

template <typename T, typename Allocator>
struct is_resizable_tensor< vect_n<T,Allocator> > : 
  is_resizable_vector< vect_n<T,Allocator> > { };

template <typename T, typename Allocator>
struct has_allocator_tensor< vect_n<T,Allocator> > : 
  has_allocator_vector< vect_n<T,Allocator> > { };





/**
 * This tensor type-traits definition is an adaptor which allows vect_n classes to act as 
 * order-one tensors.
 */
template <typename T, typename Allocator>
struct tensor_traits< std::vector<T,Allocator> > {
  /// The type of the elements of the tensor.
  typedef typename vect_traits< std::vector<T,Allocator> >::value_type value_type;
  /// The type for a reference to an element of the tensor
  typedef typename vect_traits< std::vector<T,Allocator> >::reference reference;
  /// The type for a const reference to an element of the tensor.
  typedef typename vect_traits< std::vector<T,Allocator> >::const_reference const_reference;
  /// The type for a pointer to an element of the tensor.
  typedef typename vect_traits< std::vector<T,Allocator> >::pointer pointer;
  /// The type for a const pointer to an element of the tensor.
  typedef typename vect_traits< std::vector<T,Allocator> >::const_pointer const_pointer;
  /// The type of the allocator for the tensor (can be void if the tensor does not have an allocator).
  typedef typename vect_traits< std::vector<T,Allocator> >::allocator_type allocator_type;
  
  /// The type of the size descriptors (or index descriptors) of the tensor class.
  typedef typename vect_traits< std::vector<T,Allocator> >::size_type size_type;
  /// The type of the difference between two size (or index) descriptors of the tensor class.
  typedef typename vect_traits< std::vector<T,Allocator> >::difference_type difference_type;
  
  /// The alignment tag.
  BOOST_STATIC_CONSTANT(tensor_alignment::tag, alignment = tensor_alignment::hi_dim_major);
  /// The structure tag.
  BOOST_STATIC_CONSTANT(tensor_structure::tag, structure = tensor_structure::rectangular);
  /// The tensor order.
  BOOST_STATIC_CONSTANT(std::size_t, order = 1);
};

/**
 * This tensor-dimension type-traits definition is an adaptor which allows vect_n classes to act as 
 * order-one tensors.
 */
template <std::size_t Dim, typename T, typename Allocator>
struct tensor_dim_traits<Dim, std::vector<T,Allocator> > {
  
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
template <typename T, typename Allocator>
struct tensor_dim_traits<0, std::vector<T,Allocator> > {
  
  /// The type of the iterator for the tensor and the given dimension.
  typedef typename vect_traits< std::vector<T,Allocator> >::iterator iterator;
  /// The type of the const iterator for the tensor and the given dimension.
  typedef typename vect_traits< std::vector<T,Allocator> >::const_iterator const_iterator;
  
  /// The static count in the given dimension. Should be 0 if the count is dynamic.
  BOOST_STATIC_CONSTANT(std::size_t, static_count = vect_traits< std::vector<T,Allocator> >::dimensions);
  
};



template <typename T, typename Allocator>
struct is_readable_tensor< std::vector<T,Allocator> > :
  is_readable_vector< std::vector<T,Allocator> > { };

template <typename T, typename Allocator>
struct is_writable_tensor< std::vector<T,Allocator> > :
  is_writable_vector< std::vector<T,Allocator> > { };

template <typename T, typename Allocator>
struct is_fully_writable_tensor< std::vector<T,Allocator> > 
  is_writable_vector< std::vector<T,Allocator> > { };

template <typename T, typename Allocator>
struct is_resizable_tensor< std::vector<T,Allocator> > : 
  is_resizable_vector< std::vector<T,Allocator> > { };

template <typename T, typename Allocator>
struct has_allocator_tensor< std::vector<T,Allocator> > : 
  has_allocator_vector< std::vector<T,Allocator> > { };





/**
 * This tensor type-traits definition is an adaptor which allows vect_n classes to act as 
 * order-one tensors.
 */
template <typename T, unsigned int N>
struct tensor_traits< vect<T,N> > {
  /// The type of the elements of the tensor.
  typedef typename vect_traits< vect<T,N> >::value_type value_type;
  /// The type for a reference to an element of the tensor
  typedef typename vect_traits< vect<T,N> >::reference reference;
  /// The type for a const reference to an element of the tensor.
  typedef typename vect_traits< vect<T,N> >::const_reference const_reference;
  /// The type for a pointer to an element of the tensor.
  typedef typename vect_traits< vect<T,N> >::pointer pointer;
  /// The type for a const pointer to an element of the tensor.
  typedef typename vect_traits< vect<T,N> >::const_pointer const_pointer;
  /// The type of the allocator for the tensor (can be void if the tensor does not have an allocator).
  typedef typename vect_traits< vect<T,N> >::allocator_type allocator_type;
  
  /// The type of the size descriptors (or index descriptors) of the tensor class.
  typedef typename vect_traits< vect<T,N> >::size_type size_type;
  /// The type of the difference between two size (or index) descriptors of the tensor class.
  typedef typename vect_traits< vect<T,N> >::difference_type difference_type;
  
  /// The alignment tag.
  BOOST_STATIC_CONSTANT(tensor_alignment::tag, alignment = tensor_alignment::hi_dim_major);
  /// The structure tag.
  BOOST_STATIC_CONSTANT(tensor_structure::tag, structure = tensor_structure::rectangular);
  /// The tensor order.
  BOOST_STATIC_CONSTANT(std::size_t, order = 1);
};

/**
 * This tensor-dimension type-traits definition is an adaptor which allows vect_n classes to act as 
 * order-one tensors.
 */
template <std::size_t Dim, typename T, unsigned int N>
struct tensor_dim_traits<Dim, vect<T,N> > {
  
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
template <typename T, unsigned int N>
struct tensor_dim_traits<0, vect<T,N> > {
  
  /// The type of the iterator for the tensor and the given dimension.
  typedef typename vect_traits< vect<T,N> >::iterator iterator;
  /// The type of the const iterator for the tensor and the given dimension.
  typedef typename vect_traits< vect<T,N> >::const_iterator const_iterator;
  
  /// The static count in the given dimension. Should be 0 if the count is dynamic.
  BOOST_STATIC_CONSTANT(std::size_t, static_count = vect_traits< vect_n<T,Allocator> >::dimensions);
  
};



template <typename T, unsigned int N>
struct is_readable_tensor< vect<T,N> > :
  is_readable_vector< vect<T,N> > { };

template <typename T, unsigned int N>
struct is_writable_tensor< vect<T,N> > :
  is_writable_vector< vect<T,N> > { };

template <typename T, unsigned int N>
struct is_fully_writable_tensor< vect<T,N> > 
  is_writable_vector< vect<T,N> > { };

template <typename T, unsigned int N>
struct is_resizable_tensor< vect<T,N> > : 
  is_resizable_vector< vect<T,N> > { };

template <typename T, unsigned int N>
struct has_allocator_tensor< vect<T,N> > : 
  has_allocator_vector< vect<T,N> > { };

  


template <std::size_t Dim, typename Vector>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_vector<Vector>,
    boost::mpl::equal_to< 
      boost::mpl::size_t<Dim>,
      boost::mpl::size_t<0>
    >
  >,
vect_traits< Vector > >::size_type size(const Vector& v) {
  return v.size();
};

template <std::size_t Dim, typename Vector>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_vector<Vector>,
    boost::mpl::not_equal_to< 
      boost::mpl::size_t<Dim>,
      boost::mpl::size_t<0>
    >
  >,
vect_traits< Vector > >::type::size_type size(const Vector& v) {
  return 1;
};


template <std::size_t Dim, typename Vector>
  boost::mpl::and_<
    is_resizable_vector<Vector>,
    boost::mpl::equal_to< 
      boost::mpl::size_t<Dim>,
      boost::mpl::size_t<0>
    >
  >,
void >::type resize(Vector& v, std::size_t sz) {
  v.resize(sz);
};

template <std::size_t Dim, typename Vector>
  boost::mpl::and_<
    is_resizable_vector<Vector>,
    boost::mpl::not_equal_to< 
      boost::mpl::size_t<Dim>,
      boost::mpl::size_t<0>
    >
  >,
void >::type resize(Vector&, std::size_t) { };



};




#endif










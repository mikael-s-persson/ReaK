/**
 * \file tensor_traits.hpp
 * 
 * This header declares a number of standard tensor type traits. These include 
 * tensor alignment tags (low-dim-major and hi-dim-major), tensor structure tags 
 * (rectangular, square, diagonal, scalar, nil, identity, etc.), some ReaK::rtti
 * template specializations for these tags, the main tensor_traits<> template which
 * defines the nested typedefs related to a tensor class (and required by an implementation
 * of a tensor class), and, finally, a series of meta-functions to compute the product-priority
 * additive-priority of tensor classes based on their structural tags. 
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

#ifndef REAK_TENSOR_TRAITS_HPP
#define REAK_TENSOR_TRAITS_HPP

#include <ReaK/core/base/defs.hpp>

#include <ReaK/core/rtti/rtti.hpp>

namespace ReaK {
  
  
namespace tensor_alignment {
  
  /**
   * These tags are used to mark a tensor as having a low-dim-major or hi-dim-major alignment
   * in memory. For a vector (column-vector), which is a first-order tensor, both alignments are 
   * the same. For a matrix (second-order tensor), the hi_dim_major corresponds to column-major, and 
   * low_dim_major corresponds to row-major (because rows are the first index (the lower dimension), they
   * are stored first (i.e., major)).
   */
  enum tag {
    hi_dim_major = 1,
    low_dim_major = 2
  };
};

namespace tensor_structure {
  
  /**
   * These tags are used to identify a tensor class by its underlying structure. The most 
   * general structure is, of course, the rectangular tensor (any size, any shape, and densely 
   * populated). Then, the frequently used tags are square, scalar and 
   * diagonal. Additionally, the nil and identity tags signify that a tensor type is constraint
   * to never be anything other than nil or identity (obviously, these types will be read-only
   * and do not require storage other than the size information).
   */
  enum tag {
    rectangular = 1,
    square = 2,
    diagonal = 3, // ?
    nil = 4,
    identity = 5, // ?
    scalar = 6, // ?
    permutation = 7 // ?
  };
};
  
/*
 * The following are ReaK::rtti class templates which are used to associate 
 * type information to the matrix tags (mat_alignment::tag and mat_structure::tag).
 * 
 * Please ignore, unless interested in the inner-workings of the ReaK::rtti system.
 */
namespace rtti {
  
template <>
struct get_type_id< boost::mpl::integral_c<tensor_alignment::tag, tensor_alignment::hi_dim_major > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 1);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("hi_dim_major");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "hi_dim_major"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<tensor_alignment::tag, tensor_alignment::low_dim_major > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 2);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("low_dim_major");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "low_dim_major"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};

template <mat_alignment::tag U, typename Tail>
struct get_type_info< boost::mpl::integral_c<tensor_alignment::tag, U >, Tail > {
  typedef type_id<boost::mpl::integral_c<tensor_alignment::tag, U >, typename Tail::type> type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< boost::mpl::integral_c<tensor_alignment::tag, U > >::type_name
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< boost::mpl::integral_c<tensor_alignment::tag, U > >::type_name();
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};

  
template <>
struct get_type_id< boost::mpl::integral_c<tensor_structure::tag, tensor_structure::rectangular > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 1);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("rectangular");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "rectangular"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<tensor_structure::tag, tensor_structure::square > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 2);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("square");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "square"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<tensor_structure::tag, tensor_structure::diagonal > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 3);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("diagonal");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "diagonal"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<tensor_structure::tag, tensor_structure::nil > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 4);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("nil");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "nil"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<tensor_structure::tag, tensor_structure::identity > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 5);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("identity");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "identity"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<tensor_structure::tag, tensor_structure::scalar > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 6);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("scalar");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "scalar"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};

template <>
struct get_type_id< boost::mpl::integral_c<tensor_structure::tag, tensor_structure::permutation > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 7);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("permutation");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "permutation"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
};


template <mat_structure::tag U, typename Tail>
struct get_type_info< boost::mpl::integral_c<tensor_structure::tag, U >, Tail > {
  typedef type_id<boost::mpl::integral_c<tensor_structure::tag, U >, typename Tail::type> type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< boost::mpl::integral_c<tensor_structure::tag, U > >::type_name
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< boost::mpl::integral_c<tensor_structure::tag, U > >::type_name();
    result += get_type_name_tail<Tail>::value();
    return result; //NRVO
  };
#endif
};

  
};


/**
 * This type-traits definition describes the nested typedefs that are expected 
 * from an implementation of a tensor class. They are mostly inspired from 
 * STL-containers' traits. The simplest why to implement a new tensor class 
 * that conforms with these tensor_traits requirements is to provide all those 
 * nested typedefs and static constants. The other alternative, which is non-intrusive,
 * is to define a specialization for tensor_traits for the new tensor class, providing
 * all the public members as seen in this general trait template.
 */
template <typename Tensor>
struct tensor_traits {
  /// The type of the elements of the tensor.
  typedef typename Tensor::value_type value_type;
  /// The type for a reference to an element of the tensor
  typedef typename Tensor::reference reference;
  /// The type for a const reference to an element of the tensor.
  typedef typename Tensor::const_reference const_reference;
  /// The type for a pointer to an element of the tensor.
  typedef typename Tensor::pointer pointer;
  /// The type for a const pointer to an element of the tensor.
  typedef typename Tensor::const_pointer const_pointer;
  /// The type of the allocator for the tensor (can be void if the tensor does not have an allocator).
  typedef typename Tensor::allocator_type allocator_type;
  
  /// The type of the size descriptors (or index descriptors) of the tensor class.
  typedef typename Tensor::size_type size_type;
  /// The type of the difference between two size (or index) descriptors of the tensor class.
  typedef typename Tensor::difference_type difference_type;
  
  /// The alignment tag.
  BOOST_STATIC_CONSTANT(tensor_alignment::tag, alignment = Tensor::alignment);
  /// The structure tag.
  BOOST_STATIC_CONSTANT(tensor_structure::tag, structure = Tensor::structure);
    /// The tensor order.
  BOOST_STATIC_CONSTANT(std::size_t, order = Tensor::order);
  
};


/**
 * This type-traits definition describes the nested typedefs that are expected 
 * from an implementation of a tensor class and a given query dimension.
 */
template <std::size_t Dim, typename Tensor>
struct tensor_dim_traits {
  
  /// The type of the iterator for the tensor and the given dimension.
  typedef typename Tensor::template dim<Dim>::iterator iterator;
  /// The type of the const iterator for the tensor and the given dimension.
  typedef typename Tensor::template dim<Dim>::const_iterator const_iterator;
  
  /// The static count in the given dimension. Should be 0 if the count is dynamic.
  BOOST_STATIC_CONSTANT(std::size_t, static_count = Tensor::template dim<Dim>::static_count);
  
};





};




#endif










/**
 * \file tensor_alg_general_hpp
 * 
 * This library implements the general versions of many meta-functions (templates), 
 * functions, and operators. These are meant to be used when no more-specialized 
 * implementations exist for the tensor types involved.
 * 
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_TENSOR_ALG_GENERAL_HPP
#define REAK_TENSOR_ALG_GENERAL_HPP

#include <ReaK/core/base/defs.hpp>
#include "tensor_concepts.hpp"
#include "tensor_traits.hpp"
#include <ReaK/core/lin_alg/stride_iterator.hpp>

#include <ReaK/core/base/serializable.hpp>
#include <ReaK/core/rtti/so_register_type.hpp>

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/concept_check.hpp>

#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/comparison.hpp>


namespace ReaK {


/**
 * This class is the general class template for all tensor classes in the ReaK tensor
 * libraries. The general template itself should never be used and will cause a compilation 
 * error if it is, all useful tensor class templates are, in fact, partial specializations of
 * this general class template.
 * 
 * Models: ReadableTensorConcept.
 * 
 * \tparam T          Arithmetic type of the elements of the tensor.
 * \tparam Order      The order of the tensor.
 * \tparam Structure  Enum which defines the structure of the tensor, see tensor_structure::tag.
 * \tparam Alignment  Enum which defines the memory alignment of the tensor. Either tensor_alignment::hi_dim_major or tensor_alignment::low_dim_major (default).
 * \tparam Allocator  Standard allocator class (as in the STL), the default is std::allocator<T>.
 */  
template <typename T, 
          unsigned int Order,
          tensor_structure::tag Structure = tensor_structure::rectangular,
          tensor_alignment::tag Alignment = tensor_alignment::hi_dim_major,
          typename Allocator = std::allocator<T> >
class tensor { char this_specialization_is_not_available_or_possible[0]; };

template <typename T, 
          unsigned int Order,
          tensor_structure::tag Structure,
          tensor_alignment::tag Alignment,
          typename Allocator>
struct is_readable_tensor< tensor<T,Order,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_tensor< tensor<T,Order,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          unsigned int Order,
          tensor_structure::tag Structure,
          tensor_alignment::tag Alignment,
          typename Allocator>
struct is_writable_tensor< tensor<T,Order,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_writable_tensor< tensor<T,Order,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          unsigned int Order,
          tensor_structure::tag Structure,
          tensor_alignment::tag Alignment,
          typename Allocator>
struct is_fully_writable_tensor< tensor<T,Order,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_fully_writable_tensor< tensor<T,Order,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          unsigned int Order,
          tensor_structure::tag Structure,
          tensor_alignment::tag Alignment,
          typename Allocator>
struct is_resizable_tensor< tensor<T,Order,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_resizable_tensor< tensor<T,Order,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          unsigned int Order,
          tensor_structure::tag Structure,
          tensor_alignment::tag Alignment,
          typename Allocator>
struct has_allocator_tensor< tensor<T,Order,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef has_allocator_tensor< tensor<T,Order,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          unsigned int Order,
          tensor_structure::tag Structure,
          tensor_alignment::tag Alignment,
          typename Allocator>
struct is_square_tensor< tensor<T,Order,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = ((Structure != tensor_structure::rectangular) && (Structure != tensor_structure::nil)));
  typedef is_square_tensor< tensor<T,Order,Structure,Alignment,Allocator> > type;
};

template <typename T, 
          unsigned int Order,
          tensor_structure::tag Structure,
          tensor_alignment::tag Alignment,
          typename Allocator>
struct is_diagonal_tensor< tensor<T,Order,Structure,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = ((Structure == tensor_structure::diagonal) || (Structure == tensor_structure::identity)));
  typedef is_diagonal_tensor< tensor<T,Order,Structure,Alignment,Allocator> > type;
};




namespace rtti {

template <typename T,
          unsigned int Order,
          tensor_structure::tag Structure, 
          tensor_alignment::tag Alignment,
          typename Allocator>
struct get_type_id< tensor<T,Order,Structure,Alignment,Allocator> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000030);
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA("tensor");
#else
  static const char* type_name() BOOST_NOEXCEPT { return "tensor"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return NULL; };
  
  typedef const serializable& save_type;
  typedef serializable& load_type;
};

template <typename T, 
          unsigned int Order,
          tensor_structure::tag Structure, 
          tensor_alignment::tag Alignment, 
          typename Allocator, 
          typename Tail>
struct get_type_info< tensor<T,Order,Structure,Alignment,Allocator>, Tail > {
  typedef type_id< tensor<T,Order,Structure,Alignment,Allocator>, 
    typename get_type_info_seq<T,
      boost::mpl::integral_c<unsigned int,Order>,
      boost::mpl::integral_c<tensor_structure::tag,Structure>,
      boost::mpl::integral_c<tensor_alignment::tag,Alignment> >::template with_tail<Tail>::type::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< tensor<T,Order,Structure,Alignment,Allocator> >::type_name
    + lsl_left_bracket + get_type_info_seq<T,
      boost::mpl::integral_c<unsigned int,Order>,
      boost::mpl::integral_c<tensor_structure::tag,Structure>,
      boost::mpl::integral_c<tensor_alignment::tag,Alignment> >::type_name + lsl_right_bracket
    + get_type_name_tail<Tail>::value;
#else
  static std::string type_name() { 
    std::string result = get_type_id< mat<T,Structure,Alignment,Allocator> >::type_name();
    result += "<";
    result += get_type_info_seq<T,
      boost::mpl::integral_c<unsigned int,Order>,
      boost::mpl::integral_c<tensor_structure::tag,Structure>,
      boost::mpl::integral_c<tensor_alignment::tag,Alignment> >::type_name();
    result += ">";
    result += get_type_name_tail<Tail>::value(); 
    return result; //NRVO
  };
#endif
};


};







};


#endif











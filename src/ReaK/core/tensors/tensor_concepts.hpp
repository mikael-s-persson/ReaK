/**
 * \file tensor_concepts.hpp
 * 
 * This header declares the various concepts to which a Tensor class can be expected 
 * to fulfill. These concepts include ReadableTensorConcept, WritableTensorConcept,
 * ResizableTensorConcept, and DynAllocTensorConcept. All these concepts are also 
 * paired with meta-functions that can evaluate whether a Tensor class fulfill the 
 * concept or not, and return a compile-time constant bool (on the model of 
 * boost::mpl::bool_ class). Note that these meta-functions cannot really check the 
 * concepts directly (this is impossible in current and future C++ standard versions, might
 * eventually be part of the standard, but not in the forseeable future). These meta-functions
 * are constant meta-functions that always return false, the implementer of a given tensor
 * class should also define specializations of these meta-functions with their appropriate 
 * return values.
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

#ifndef REAK_TENSOR_CONCEPTS_HPP
#define REAK_TENSOR_CONCEPTS_HPP

#include "tensor_traits.hpp"
#include <boost/concept_check.hpp>

#include <boost/type_traits.hpp>


namespace ReaK {

/**
 * This concept will fail to be instantiated if the Tensor class does not model 
 * the readable tensor concept.
 * 
 * Required expressions for Tensor t:
 *  e = t(i, ..., i_n-1)   read access to the elements of t with indices.
 *  s = size< 0 >(t)       can obtain the number of elements along dimension 0 of the tensor.
 *  ...
 *  s = size< n-1 >(t)     can obtain the number of elements along dimension n-1 of the tensor.
 */
template <std::size_t Order, typename Tensor>
struct ReadableTensorConcept {
  BOOST_CONCEPT_USAGE(ReadableMatrixConcept) {
    char high_order_tensors_are_not_supported[0];
  };
};

template <typename Tensor>
struct ReadableTensorConcept<0,Tensor> { 
  
  BOOST_CONCEPT_USAGE(ReadableTensorConcept) {
  };
};

template <typename Tensor>
struct ReadableTensorConcept<1,Tensor> {
  Tensor t;
  
  typename tensor_traits<Tensor>::size_type s;
  typename tensor_traits<Tensor>::value_type e;
    
  BOOST_CONCEPT_USAGE(ReadableTensorConcept) {
    e = t(0);
    s = size<0>(t);
  };
  
};

template <typename Tensor>
struct ReadableTensorConcept<2,Tensor> {
  Tensor t;
  
  typename tensor_traits<Tensor>::size_type s;
  typename tensor_traits<Tensor>::value_type e;
    
  BOOST_CONCEPT_USAGE(ReadableTensorConcept) {
    e = t(0,0);
    s = size<0>(t);
    s = size<1>(t);
  };
  
};

template <typename Tensor>
struct ReadableTensorConcept<3,Tensor> {
  Tensor t;
  
  typename tensor_traits<Tensor>::size_type s;
  typename tensor_traits<Tensor>::value_type e;
    
  BOOST_CONCEPT_USAGE(ReadableTensorConcept) {
    e = t(0,0,0);
    s = size<0>(t);
    s = size<1>(t);
    s = size<2>(t);
  };
  
};

template <typename Tensor>
struct ReadableTensorConcept<4,Tensor> {
  Tensor t;
  
  typename tensor_traits<Tensor>::size_type s;
  typename tensor_traits<Tensor>::value_type e;
    
  BOOST_CONCEPT_USAGE(ReadableTensorConcept) {
    e = t(0,0,0,0);
    s = size<0>(t);
    s = size<1>(t);
    s = size<2>(t);
    s = size<3>(t);
  };
  
};



/**
 * This meta-function evaluates whether a Tensor class fulfills the ReadableTensorConcept, 
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results 
 * in a false value, and the implementer of a tensor class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new tensor class.
 */
template <typename Tensor>
struct is_readable_tensor {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_readable_tensor<Tensor> type;
};

/**
 * This concept will fail to be instantiated if the Tensor class does not model 
 * the writable tensor concept, it must also fulfill the ReadableTensorConcept.
 * 
 * Required expressions for Tensor t in addition to that of ReadableTensorConcept:
 *  m(i, ..., i_n-1) = e;   write access to the elements of the given indices.
 */
template <std::size_t Order, typename Tensor>
struct WritableTensorConcept { 
  BOOST_CONCEPT_USAGE(WritableTensorConcept) {
    char high_order_tensors_are_not_supported[0];
  };
};

template <typename Tensor>
struct WritableTensorConcept<0,Tensor> : 
    ReadableTensorConcept<0,Tensor> { //must also be readable.
  Tensor t;
  
  typename tensor_traits<Tensor>::value_type r;
  
  BOOST_CONCEPT_USAGE(WritableTensorConcept) {
  };
  
};

template <typename Tensor>
struct WritableTensorConcept<1,Tensor> : 
    ReadableTensorConcept<1,Tensor> { //must also be readable.
  Tensor t;
  
  typename tensor_traits<Tensor>::value_type r;
  
  BOOST_CONCEPT_USAGE(WritableTensorConcept) {
    m(0) = r; //can be indexed and given an lvalue
  };
  
};

template <typename Tensor>
struct WritableTensorConcept<2,Tensor> : 
    ReadableTensorConcept<2,Tensor> { //must also be readable.
  Tensor t;
  
  typename tensor_traits<Tensor>::value_type r;
  
  BOOST_CONCEPT_USAGE(WritableTensorConcept) {
    m(0,0) = r; //can be indexed and given an lvalue
  };
  
};

template <typename Tensor>
struct WritableTensorConcept<3,Tensor> : 
    ReadableTensorConcept<3,Tensor> { //must also be readable.
  Tensor t;
  
  typename tensor_traits<Tensor>::value_type r;
  
  BOOST_CONCEPT_USAGE(WritableTensorConcept) {
    m(0,0,0) = r; //can be indexed and given an lvalue
  };
  
};

template <typename Tensor>
struct WritableTensorConcept<4,Tensor> : 
    ReadableTensorConcept<4,Tensor> { //must also be readable.
  Tensor t;
  
  typename tensor_traits<Tensor>::value_type r;
  
  BOOST_CONCEPT_USAGE(WritableTensorConcept) {
    m(0,0,0,0) = r; //can be indexed and given an lvalue
  };
  
};



/**
 * This meta-function evaluates whether a Tensor class fulfills the WritableTensorConcept, 
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results 
 * in a false value, and the implementer of a tensor class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new tensor class.
 */
template <typename Tensor>
struct is_writable_tensor {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_tensor<Tensor> type;
};


/**
 * This meta-function evaluates whether a Tensor class fulfills the WritableTensorConcept and
 * can be considered as "fully writable" meaning that all the elements are independent and writable, 
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results 
 * in a false value, and the implementer of a tensor class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new tensor class.
 */
template <typename Tensor>
struct is_fully_writable_tensor {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_fully_writable_tensor<Tensor> type;
};


/**
 * This concept will fail to be instantiated if the Tensor class does not model 
 * the resizable tensor concept, meaning that its sizes can be 
 * set to some values which results in the tensor having at least that size in the 
 * given dimension.
 * 
 * Required expressions for Tensor t:
 *  size< 0 >(t) = s;    can set the size of the dimension 0 of the tensor.
 *  ...
 *  size< n-1 >(t) = s;  can set the size of the dimension n-1 of the tensor.
 */
template <std::size_t Order, typename Tensor>
struct ResizableTensorConcept { 
  Tensor t;
  
  typename tensor_traits<Tensor>::size_type sz;
  
  BOOST_CONCEPT_USAGE(ResizableTensorConcept) {
    char high_order_tensors_are_not_supported[0];
  };
  
};

template <typename Tensor>
struct ResizableTensorConcept<0,Tensor> { 
  Tensor t;
  
  typename tensor_traits<Tensor>::size_type sz;
  
  BOOST_CONCEPT_USAGE(ResizableTensorConcept) {
  };
  
};

template <typename Tensor>
struct ResizableTensorConcept<1,Tensor> { 
  Tensor t;
  
  typename tensor_traits<Tensor>::size_type sz;
  
  BOOST_CONCEPT_USAGE(ResizableTensorConcept) {
    resize<0>(t,sz);
  };
  
};

template <typename Tensor>
struct ResizableTensorConcept<2,Tensor> { 
  Tensor t;
  
  typename tensor_traits<Tensor>::size_type sz;
  
  BOOST_CONCEPT_USAGE(ResizableTensorConcept) {
    resize<0>(t,sz);
    resize<1>(t,sz);
  };
  
};

template <typename Tensor>
struct ResizableTensorConcept<3,Tensor> { 
  Tensor t;
  
  typename tensor_traits<Tensor>::size_type sz;
  
  BOOST_CONCEPT_USAGE(ResizableTensorConcept) {
    resize<0>(t,sz);
    resize<1>(t,sz);
    resize<2>(t,sz);
  };
  
};

template <typename Tensor>
struct ResizableTensorConcept<4,Tensor> { 
  Tensor t;
  
  typename tensor_traits<Tensor>::size_type sz;
  
  BOOST_CONCEPT_USAGE(ResizableTensorConcept) {
    resize<0>(t,sz);
    resize<1>(t,sz);
    resize<2>(t,sz);
    resize<3>(t,sz);
  };
  
};


/**
 * This meta-function evaluates whether a Tensor class fulfills the ResizableTensorConcept, 
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results 
 * in a false value, and the implementer of a tensor class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new tensor class.
 */
template <typename Tensor>
struct is_resizable_tensor {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_tensor<Tensor> type;
};

/**
 * This concept will fail to be instantiated if the Tensor class does not model 
 * the dynamically allocated tensor concept, meaning that it stores its elements in 
 * dynamically allocated memory. The only requirement to fulfill this concept is 
 * to be able to obtain the allocate object associated to a tensor (here "allocator" is
 * used with the exact same meaning as "STL allocators").
 * 
 * Required expressions for Tensor t:
 *  al = t.get_allocator()  can obtain the allocator object of the tensor.
 */
template <std::size_t Order, typename Tensor>
struct DynAllocTensorConcept : ResizableTensorConcept<Order,Tensor> { 
  Tensor t;
  
  typename tensor_traits<Tensor>::allocator_type al;
  
  BOOST_CONCEPT_USAGE(DynAllocTensorConcept) {
    al = t.get_allocator();
  };
  
};


/**
 * This meta-function evaluates whether a Tensor class fulfills the DynAllocTensorConcept, 
 * however, it does not attempt to instantiate the Concept template (because no technique can
 * be used to catch the failed instantiation properly), instead, the default version results 
 * in a false value, and the implementer of a tensor class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new tensor class.
 */
template <typename Tensor>
struct has_allocator_tensor {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef has_allocator_tensor<Tensor> type;
};





/**
 * This meta-function evaluates whether a Tensor class is a square tensor. The implementer of 
 * a tensor class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new tensor class.
 */
template <typename Tensor>
struct is_square_tensor {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_square_tensor<Tensor> type;
};


/**
 * This meta-function evaluates whether a Tensor class is a diagonal tensor. The implementer of 
 * a tensor class is required to provide a specialization
 * if he wants this meta-function to evaluate to true for that new tensor class.
 */
template <typename Tensor>
struct is_diagonal_tensor {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_diagonal_tensor<Tensor> type;
};







};



#endif









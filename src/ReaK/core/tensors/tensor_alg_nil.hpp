/**
 * \file tensor_alg_nil.hpp
 * 
 * This library implements the specialization of the tensor<> template for a 
 * general nil tensor (all zero values). This tensor type fulfills the tensor 
 * concepts (Readable and Resizable).
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

#ifndef REAK_TENSOR_ALG_NIL_HPP
#define REAK_TENSOR_ALG_NIL_HPP

#include "tensor_alg_general.hpp"

namespace ReaK {
  

template <typename T, 
          unsigned int Order,
	  tensor_alignment::tag Alignment,
	  typename Allocator>
struct is_writable_tensor< tensor<T,Order,tensor_structure::nil,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_tensor< tensor<T,Order,tensor_structure::nil,Alignment,Allocator> > type;
};

template <typename T, 
          unsigned int Order,
	  tensor_alignment::tag Alignment,
	  typename Allocator>
struct is_fully_writable_tensor< tensor<T,Order,tensor_structure::nil,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_fully_writable_tensor< tensor<T,Order,tensor_structure::nil,Alignment,Allocator> > type;
};

template <typename T, 
          unsigned int Order,
	  tensor_alignment::tag Alignment,
	  typename Allocator>
struct has_allocator_tensor< tensor<T,Order,tensor_structure::nil,Alignment,Allocator> > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef has_allocator_tensor< tensor<T,Order,tensor_structure::nil,Alignment,Allocator> > type;
};


/**
 * This class template specialization implements a tensor with zero elements. This class is 
 * serializable and registered to the ReaK::rtti system. This tensor type is resizable.
 * 
 * Models: ReadableTensorConcept, and ResizableTensorConcept.
 * 
 * \tparam T Arithmetic type of the elements of the tensor.
 * \tparam Alignment Alignment of the elements stored in this tensor class.
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T,
	  tensor_alignment::tag Alignment,
	  typename Allocator>
class tensor<T,3,tensor_structure::nil,Alignment,Allocator> : public serialization::serializable {
  public:    
    
    typedef tensor<T,3,tensor_structure::nil,Alignment,Allocator> self;
    typedef void allocator_type;
    
    typedef T value_type;
    typedef void container_type;
    
    typedef void reference;
    typedef T const_reference;
    typedef void pointer;
    typedef void const_pointer;
    
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, order = 3);
    BOOST_STATIC_CONSTANT(tensor_alignment::tag, alignment = Alignment);
    BOOST_STATIC_CONSTANT(tensor_structure::tag, structure = tensor_structure::nil);
    
    template <unsigned int Dim>
    struct dim {
      typedef void iterator;
      typedef void const_iterator;
      BOOST_STATIC_CONSTANT(std::size_t, static_count = 0);
    };
    
  
  private:
    size_type counts[3]; ///< Counts.
  public:  
    
/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default constructor. Sets dimensions to zero.
     */
    tensor() 
#ifdef RK_ENABLE_CXX0X_FEATURES
	     : counts{0,0,0} { };
#else
	     { counts[0] = 0; counts[1] = 0; counts[2] = 0; };
#endif
    /**
     * Constructs a null matrix to the given dimensions.
     */
    tensor(size_type aCount0, size_type aCount1, size_type aCount2) 
#ifdef RK_ENABLE_CXX0X_FEATURES
	   : counts{aCount0, aCount1, aCount2} { };
#else
	   { counts[0] = aCount0; counts[1] = aCount1; counts[2] = aCount2; };
#endif
    
    tensor(const self& rhs) { std::copy(&M.counts[0],&M.counts[0] + 3,&counts[0]); };
    /**
     * Default destructor.
     */
    ~tensor() { };
    
    /**
     * The standard swap function (works with ADL).
     */
    friend void swap(self& lhs,self& rhs) throw() {
      using std::swap;
      swap(m1.counts[0],m2.counts[0]);
      swap(m1.counts[1],m2.counts[1]);
      swap(m1.counts[2],m2.counts[2]);
    };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Tensor indexing accessor for read-only access.
     * \param i Index-0.
     * \param j Index-1.
     * \param k Index-2.
     * \return the element at the given position.
     */
    const_reference operator()(size_type i, size_type j, size_type k) const { RK_UNUSED(i); RK_UNUSED(j); RK_UNUSED(k); return value_type(0); };

    /**
     * Gets the count for the given dimension of the tensor.
     * \param m The tensor along of which the size is sought.
     * \return number of elements along the given dimension of the tensor.
     * \tparam Dim The dimension along which the size is sought.
     */
    template <unsigned int Dim>
    friend 
    typename boost::enable_if<
      boost::mpl::less< boost::mpl::size_t< Dim >, boost::mpl::size_t< 3 > >,
    size_type >::type size(const self& m) {
      return m.counts[Dim];
    };
    
    /**
     * Gets the count for the given dimension of the tensor.
     * \param m The tensor along of which the size is sought.
     * \return number of elements along the given dimension of the tensor.
     * \tparam Dim The dimension along which the size is sought.
     */
    template <unsigned int Dim>
    friend 
    typename boost::disable_if<
      boost::mpl::less< boost::mpl::size_t< Dim >, boost::mpl::size_t< 3 > >,
    size_type >::type size(const self& m) { RK_UNUSED(m);
      return 1;
    };
    

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/
    
    /** 
     * Scalar-multiply-and-store operator with standard semantics.
     */
    self& operator *=(const value_type& S) {
      return *this;
    };
    
    /** 
     * General negation operator for any type of tensors. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General hi-dim-major tensor.
     */
    self operator -() const {
      return *this;
    };
    
    
    /**
     * Transposes the tensor M by a simple copy with a change of alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend self transpose(const self& M) {
      return M;
    };
    
    /**
     * Transposes the tensor M by simply moving the data of M into a tensor of different alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend self transpose_move(self& M) {
      return M;
    };

#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Transposes the matrix M by simply moving the data of M into a matrix of different alignment.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend self transpose(self&& M) {
      return M;
    };    
#endif
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & std::pair<std::string, unsigned int>("counts[0]",counts[0])
	& std::pair<std::string, unsigned int>("counts[1]",counts[1])
	& std::pair<std::string, unsigned int>("counts[2]",counts[2]);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      unsigned int tmp;
      A & std::pair<std::string, unsigned int&>("counts[0]",tmp);
      counts[0] = tmp;
      A & std::pair<std::string, unsigned int&>("counts[1]",tmp);
      counts[1] = tmp;
      A & std::pair<std::string, unsigned int&>("counts[2]",tmp);
      counts[2] = tmp;
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serialization::serializable)
  
    
  
};






/**
 * This class template specialization implements a tensor with zero elements. This class is 
 * serializable and registered to the ReaK::rtti system. This tensor type is resizable.
 * 
 * Models: ReadableTensorConcept, and ResizableTensorConcept.
 * 
 * \tparam T Arithmetic type of the elements of the tensor.
 * \tparam Alignment Alignment of the elements stored in this tensor class.
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T,
	  tensor_alignment::tag Alignment,
	  typename Allocator>
class tensor<T,4,tensor_structure::nil,Alignment,Allocator> : public serialization::serializable {
  public:    
    
    typedef tensor<T,4,tensor_structure::nil,Alignment,Allocator> self;
    typedef void allocator_type;
    
    typedef T value_type;
    typedef void container_type;
    
    typedef void reference;
    typedef T const_reference;
    typedef void pointer;
    typedef void const_pointer;
    
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, order = 4);
    BOOST_STATIC_CONSTANT(tensor_alignment::tag, alignment = Alignment);
    BOOST_STATIC_CONSTANT(tensor_structure::tag, structure = tensor_structure::nil);
    
    template <unsigned int Dim>
    struct dim {
      typedef void iterator;
      typedef void const_iterator;
      BOOST_STATIC_CONSTANT(std::size_t, static_count = 0);
    };
    
  
  private:
    size_type counts[4]; ///< Counts.
  public:  
    
/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default constructor. Sets dimensions to zero.
     */
    tensor() 
#ifdef RK_ENABLE_CXX0X_FEATURES
	     : counts{0,0,0,0} { };
#else
	     { counts[0] = 0; counts[1] = 0; counts[2] = 0; counts[3] = 0; };
#endif
    /**
     * Constructs a null matrix to the given dimensions.
     */
    tensor(size_type aCount0, size_type aCount1, size_type aCount2, size_type aCount3) 
#ifdef RK_ENABLE_CXX0X_FEATURES
	   : counts{aCount0, aCount1, aCount2, aCount3} { };
#else
	   { counts[0] = aCount0; counts[1] = aCount1; counts[2] = aCount2; counts[3] = aCount3; };
#endif
    
    tensor(const self& rhs) { std::copy(&M.counts[0],&M.counts[0] + 4,&counts[0]); };
    /**
     * Default destructor.
     */
    ~tensor() { };
    
    /**
     * The standard swap function (works with ADL).
     */
    friend void swap(self& lhs,self& rhs) throw() {
      using std::swap;
      swap(m1.counts[0],m2.counts[0]);
      swap(m1.counts[1],m2.counts[1]);
      swap(m1.counts[2],m2.counts[2]);
      swap(m1.counts[3],m2.counts[3]);
    };
    
/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Tensor indexing accessor for read-only access.
     * \param i Index-0.
     * \param j Index-1.
     * \param k Index-2.
     * \param l Index-3.
     * \return the element at the given position.
     */
    const_reference operator()(size_type i, size_type j, size_type k, size_type l) const { RK_UNUSED(i); RK_UNUSED(j); RK_UNUSED(k); RK_UNUSED(l); return value_type(0); };

    /**
     * Gets the count for the given dimension of the tensor.
     * \param m The tensor along of which the size is sought.
     * \return number of elements along the given dimension of the tensor.
     * \tparam Dim The dimension along which the size is sought.
     */
    template <unsigned int Dim>
    friend 
    typename boost::enable_if<
      boost::mpl::less< boost::mpl::size_t< Dim >, boost::mpl::size_t< 4 > >,
    size_type >::type size(const self& m) {
      return m.counts[Dim];
    };
    
    /**
     * Gets the count for the given dimension of the tensor.
     * \param m The tensor along of which the size is sought.
     * \return number of elements along the given dimension of the tensor.
     * \tparam Dim The dimension along which the size is sought.
     */
    template <unsigned int Dim>
    friend 
    typename boost::disable_if<
      boost::mpl::less< boost::mpl::size_t< Dim >, boost::mpl::size_t< 4 > >,
    size_type >::type size(const self& m) { RK_UNUSED(m);
      return 1;
    };
    

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/
    
    /** 
     * Scalar-multiply-and-store operator with standard semantics.
     */
    self& operator *=(const value_type& S) {
      return *this;
    };
    
    /** 
     * General negation operator for any type of tensors. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General hi-dim-major tensor.
     */
    self operator -() const {
      return *this;
    };
    
    
    /**
     * Transposes the tensor M by a simple copy with a change of alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend self transpose(const self& M) {
      return M;
    };
    
    /**
     * Transposes the tensor M by simply moving the data of M into a tensor of different alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend self transpose_move(self& M) {
      return M;
    };

#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Transposes the matrix M by simply moving the data of M into a matrix of different alignment.
     * \param M The matrix to be transposed.
     * \return The transpose of M.
     */
    friend self transpose(self&& M) {
      return M;
    };    
#endif
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & std::pair<std::string, unsigned int>("counts[0]",counts[0])
	& std::pair<std::string, unsigned int>("counts[1]",counts[1])
	& std::pair<std::string, unsigned int>("counts[2]",counts[2])
	& std::pair<std::string, unsigned int>("counts[3]",counts[3]);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      unsigned int tmp;
      A & std::pair<std::string, unsigned int&>("counts[0]",tmp);
      counts[0] = tmp;
      A & std::pair<std::string, unsigned int&>("counts[1]",tmp);
      counts[1] = tmp;
      A & std::pair<std::string, unsigned int&>("counts[2]",tmp);
      counts[2] = tmp;
      A & std::pair<std::string, unsigned int&>("counts[3]",tmp);
      counts[3] = tmp;
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serialization::serializable)
  
};



  
};

#endif








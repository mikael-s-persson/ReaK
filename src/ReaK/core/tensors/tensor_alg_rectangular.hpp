/**
 * \file tensor_alg_rectangular.hpp
 * 
 * This library implements the specialization of the tensor<> template for a 
 * general rectangular tensor (dense, dynamic dimensions) of both hi-dim-major and 
 * low-dim-major alignment. This tensor type fulfills all the general tensor 
 * concepts (Readable, Writable, Fully-Writable, Resizable and DynAlloc).
 * 
 * This library also implements transposition of tensors via alignment 
 * switching (switching from hi-dim-major to low-dim-major, or vice versa). This 
 * is very efficient and can even avoid copies completely (via transpose_move()) on an
 * optimizing compiler.
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

#ifndef REAK_TENSOR_ALG_RECTANGULAR_HPP
#define REAK_TENSOR_ALG_RECTANGULAR_HPP

#include "tensor_alg_general.hpp"

namespace ReaK {
  



/**
 * This class template specialization implements a tensor with rectangular structure
 * and hi-dim-major alignment. This class is serializable and registered to the ReaK::rtti
 * system. This tensor type is dynamically resizable.
 * 
 * Models: ReadableTensorConcept, WritableTensorConcept, FullyWritableTensorConcept, 
 * ResizableTensorConcept, and DynAllocTensorConcept.
 * 
 * \tparam T Arithmetic type of the elements of the tensor.
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T,
          typename Allocator>
class tensor<T,3,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator> : public serialization::serializable {
  public:    
    
    typedef tensor<T,3,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator> self;
    typedef Allocator allocator_type;
    
    typedef T value_type;
    typedef std::vector<value_type,allocator_type> container_type;
    
    typedef typename container_type::reference reference;
    typedef typename container_type::const_reference const_reference;
    typedef typename container_type::pointer pointer;
    typedef typename container_type::const_pointer const_pointer;

    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, order = 3);
    BOOST_STATIC_CONSTANT(tensor_alignment::tag, alignment = tensor_alignment::hi_dim_major);
    BOOST_STATIC_CONSTANT(tensor_structure::tag, structure = tensor_structure::rectangular);
    
    template <unsigned int Dim>
    struct dim {
      typedef void iterator;
      typedef void const_iterator;
      BOOST_STATIC_CONSTANT(std::size_t, static_count = 0);
    };
    
  
  private:
    std::vector<value_type,allocator_type> q; ///< Array which holds all the values of the tensor (dimension: counts[0] x counts[1] x counts[2]).
    size_type counts[3]; ///< Counts.
  public:  
    
/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default constructor: sets all to zero.
     */
    tensor(const allocator_type& aAlloc = allocator_type()) :
             q(0,value_type(),aAlloc)
#ifndef BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX
             , counts{0,0,0} { };
#else
             { counts[0] = 0; counts[1] = 0; counts[2] = 0; };
#endif

    /**
     * Constructor for a sized tensor.
     */
    tensor(size_type aCount0, size_type aCount1, size_type aCount2, 
           const value_type& aFill = value_type(), const allocator_type& aAlloc = allocator_type()) :
           q(aCount0 * aCount1 * aCount2,aFill,aAlloc)
#ifndef BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX
           , counts{aCount0, aCount1, aCount2} { };
#else
           { counts[0] = aCount0; counts[1] = aCount1; counts[2] = aCount2; };
#endif

    /**
     * Standard Copy Constructor with standard semantics.
     */
    tensor(const self& M) :
             q(M.q) { std::copy(&M.counts[0],&M.counts[0] + 3,&counts[0]); };
         
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Standard Copy Constructor with standard semantics.
     */
    tensor(self&& M) :
             q(std::move(M.q)) { std::move(&M.counts[0],&M.counts[0] + 3,&counts[0]); };
#endif

    /**
     * Explicit constructor from a any type of tensor.
     */
    template <typename Tensor>
    explicit tensor(const Tensor& M, typename boost::enable_if< 
                                                 boost::mpl::and_<
                                                   is_readable_tensor<Tensor>, 
                                                   boost::mpl::not_< boost::is_same<Tensor,self> >
                                                 >, void* >::type dummy = NULL) :
             q(size<0>(M) * size<1>(M) * size<2>(M),T(0.0)) {
      counts[0] = size<0>(M); counts[1] = size<1>(M); counts[2] = size<2>(M);
      typename container_type::iterator it = q.begin();
      for(size_type i = 0; i < counts[2]; ++i)
        for(size_type j = 0; j < counts[1]; ++j)
          for(size_type k = 0; k < counts[0]; ++k, ++it)
            *it = M(k,j,i);
    };

    /**
     * Constructor from a vector of hi-dim-major values.
     */
    tensor(const container_type& Q, size_type aCount0, size_type aCount1, size_type aCount2) :
             q(Q)
#ifndef BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX
           , counts{aCount0, aCount1, aCount2} { };
#else
           { counts[0] = aCount0; counts[1] = aCount1; counts[2] = aCount2; };
#endif

    /**
     * Destructor.
     */
    ~tensor() { };
    
    /**
     * The standard swap function (works with ADL).
     */
    friend void swap(self& m1, self& m2) throw() {
      using std::swap;
      swap(m1.q,m2.q);
      swap(m1.counts[0],m2.counts[0]);
      swap(m1.counts[1],m2.counts[1]);
      swap(m1.counts[2],m2.counts[2]);
    };
    
    /**
     * A swap function to swap the tensor with a container of values to fill the tensor.
     * \param m1 The tensor to swap with the container.
     * \param q2 The container that will be swapped with m1's internal container.
     * \param Count02 The count0 corresponding to q2's data.
     * \param Count12 The count1 corresponding to q2's data.
     * \param Count22 The count2 corresponding to q2's data.
     */
    friend void swap(self& m1, container_type& q2, size_type& Count02, size_type& Count12, size_type& Count22) throw() {
      using std::swap;
      swap(m1.q,q2);
      swap(m1.counts[0],Count02);
      swap(m1.counts[1],Count12);
      swap(m1.counts[2],Count22);
    };
    
    /**
     * Standard copy-assignment operator (and move-assignment operator, for C++0x). Uses the copy-and-swap (and
     * move-and-swap) idiom.
     */
    self& operator=(self rhs) {
      swap(*this, rhs);
      return *this;
    };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Tensor indexing accessor for read-write access.
     * \param i Index-0.
     * \param j Index-1.
     * \param k Index-2.
     * \return the element at the given position.
     */
    reference operator()(size_type i, size_type j, size_type k) { return q[k * counts[0] * counts[1] + j * counts[0] + i]; };
    /**
     * Tensor indexing accessor for read-only access.
     * \param i Index-0.
     * \param j Index-1.
     * \param k Index-2.
     * \return the element at the given position.
     */
    const_reference operator()(size_type i, size_type j, size_type k) const { return q[k * counts[0] * counts[1] + j * counts[0] + i]; };

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
    
    /**
     * Sets the count for the given dimension of the tensor.
     * \param number of elements along the given dimension of the tensor.
     * \tparam Dim The dimension along which the resize is to apply.
     */
    template <unsigned int Dim>
    friend 
    typename boost::enable_if<
      boost::mpl::less< boost::mpl::size_t< Dim >, boost::mpl::size_t< 3 > >,
    void >::type resize(self& m, size_type sz) {
      size_type new_counts[3] = {m.counts[0], m.counts[1], m.counts[2]};
      new_counts[Dim] = sz;
      self new_m(new_counts[0],new_counts[1],new_counts[2]);
      for(size_type k = 0; k < std::min(new_m.counts[2], m.counts[2]); ++k)
        for(size_type j = 0; j < std::min(new_m.counts[1], m.counts[1]); ++j)
          for(size_type i = 0; i < std::min(new_m.counts[0], m.counts[0]); ++i)
            new_m.q[(k * new_m.counts[1] + j) * new_m.counts[0] + i] = m.q[(k * m.counts[1] + j) * m.counts[0] + i]; 
      swap(new_m,m);
    };
    
    /**
     * Sets the count for the given dimension of the tensor.
     * \param number of elements along the given dimension of the tensor.
     * \tparam Dim The dimension along which the resize is to apply.
     */
    template <unsigned int Dim>
    friend 
    typename boost::disable_if<
      boost::mpl::less< boost::mpl::size_t< Dim >, boost::mpl::size_t< 3 > >,
    void >::type resize(self& m, size_type sz) { RK_UNUSED(m); RK_UNUSED(sz); };
    
    
    /**
     * Returns the allocator object of the underlying container.
     * \return the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return q.get_allocator(); };
    
    

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /** 
     * Standard Assignment operator with standard semantics.
     * Strong exception-safety.
     */
    template <typename Tensor>
    self& operator =(const Tensor& M) {
      self tmp(M);
      swap(*this,tmp);
      return *this;
    };

    /** 
     * Add-and-store operator with standard semantics.
     */
    template <typename Tensor>
    typename boost::enable_if< 
      is_readable_tensor<Tensor>,
    self& >::type operator +=(const Tensor& M) {
      if((size<0>(M) != counts[0]) || 
         (size<1>(M) != counts[1]) || 
         (size<2>(M) != counts[2]))
        throw std::range_error("Tensor dimensions mismatch.");
      typename container_type::iterator it = q.begin();
      for(size_type k = 0; k < counts[2]; ++k)
        for(size_type j = 0; j < counts[1]; ++j)
          for(size_type i = 0; i < counts[0]; ++i, ++it)
            *it += M(i,j,k);
      return *this;
    };

    /** 
     * Sub-and-store operator with standard semantics.
     */
    template <typename Tensor>
    typename boost::enable_if< 
      is_readable_tensor<Tensor>,
    self& >::type operator -=(const Tensor& M) {
      if((size<0>(M) != counts[0]) || 
         (size<1>(M) != counts[1]) || 
         (size<2>(M) != counts[2]))
        throw std::range_error("Tensor dimensions mismatch.");
      typename container_type::iterator it = q.begin();
      for(size_type k = 0; k < counts[2]; ++k)
        for(size_type j = 0; j < counts[1]; ++j)
          for(size_type i = 0; i < counts[0]; ++i, ++it)
            *it -= M(i,j,k);
      return *this;
    };

    /** 
     * Scalar-multiply-and-store operator with standard semantics.
     */
    self& operator *=(const value_type& S) {
      for(typename container_type::iterator it = q.begin();it != q.end();++it)
        *it *= S;
      return *this;
    };
    
    /** 
     * General negation operator for any type of tensors. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General hi-dim-major tensor.
     */
    self operator -() const {
      self result(*this);
      typename container_type::iterator itr = result.q.begin();
      for(typename container_type::const_iterator it = q.begin(); it != q.end(); ++it, ++itr)
        *itr = -(*it);
      return result;
    };
    
    
    
    /**
     * Transposes the tensor M by a simple copy with a change of alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend tensor<T,3,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator> transpose(const self& M) {
      return tensor<T,3,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator>(M.q,M.counts[2],M.counts[1],M.counts[0]);
    };
    
    /**
     * Transposes the tensor M by simply moving the data of M into a tensor of different alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend tensor<T,3,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator> transpose_move(self& M) {
      tensor<T,3,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator> result;
      using std::swap;
      swap(result,M.q,M.counts[2],M.counts[1],M.counts[0]);
      return result;
    };

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Transposes the tensor M by simply moving the data of M into a tensor of different alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend tensor<T,3,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator> transpose(self&& M) {
      tensor<T,3,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator> result;
      using std::swap;
      swap(result,M.q,M.counts[2],M.counts[1],M.counts[0]);
      return result;
    };    
#endif
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & std::pair<std::string, const std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int>("counts[0]",counts[0])
        & std::pair<std::string, unsigned int>("counts[1]",counts[1])
        & std::pair<std::string, unsigned int>("counts[2]",counts[2]);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      unsigned int tmp;
      A & std::pair<std::string, std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int&>("counts[0]",tmp);
      counts[0] = tmp;
      A & std::pair<std::string, unsigned int&>("counts[1]",tmp);
      counts[1] = tmp;
      A & std::pair<std::string, unsigned int&>("counts[2]",tmp);
      counts[2] = tmp;
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serialization::serializable)
  
    
  
};


/**
 * This class template specialization implements a tensor with rectangular structure
 * and low-dim-major alignment. This class is serializable and registered to the ReaK::rtti
 * system. This tensor type is dynamically resizable.
 * 
 * Models: ReadableTensorConcept, WritableTensorConcept, FullyWritableTensorConcept, 
 * ResizableTensorConcept, and DynAllocTensorConcept.
 * 
 * \tparam T Arithmetic type of the elements of the tensor.
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T,
          typename Allocator>
class tensor<T,3,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator> : public serialization::serializable {
  public:    
    
    typedef tensor<T,3,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator> self;
    typedef Allocator allocator_type;
    
    typedef T value_type;
    typedef std::vector<value_type,allocator_type> container_type;
    
    typedef typename container_type::reference reference;
    typedef typename container_type::const_reference const_reference;
    typedef typename container_type::pointer pointer;
    typedef typename container_type::const_pointer const_pointer;
  
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, order = 3);
    BOOST_STATIC_CONSTANT(tensor_alignment::tag, alignment = tensor_alignment::low_dim_major);
    BOOST_STATIC_CONSTANT(tensor_structure::tag, structure = tensor_structure::rectangular);
    
    
    template <unsigned int Dim>
    struct dim {
      typedef void iterator;
      typedef void const_iterator;
      BOOST_STATIC_CONSTANT(std::size_t, static_count = 0);
    };
    
  
  private:
    std::vector<value_type,allocator_type> q; ///< Array which holds all the values of the tensor (dimension: counts[0] x counts[1] x counts[2]).
    size_type counts[3]; ///< Counts.
  public:  
    
/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default constructor: sets all to zero.
     */
    tensor(const allocator_type& aAlloc = allocator_type()) :
             q(0,value_type(),aAlloc)
#ifndef BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX
             , counts{0,0,0} { };
#else
             { counts[0] = 0; counts[1] = 0; counts[2] = 0; };
#endif

    /**
     * Constructor for a sized tensor.
     */
    tensor(size_type aCount0, size_type aCount1, size_type aCount2, 
           const value_type& aFill = value_type(), const allocator_type& aAlloc = allocator_type()) :
           q(aCount0 * aCount1 * aCount2,aFill,aAlloc)
#ifndef BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX
           , counts{aCount0, aCount1, aCount2} { };
#else
           { counts[0] = aCount0; counts[1] = aCount1; counts[2] = aCount2; };
#endif

    /**
     * Standard Copy Constructor with standard semantics.
     */
    tensor(const self& M) :
             q(M.q) { std::copy(&M.counts[0],&M.counts[0] + 3,&counts[0]); };
         
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Standard Copy Constructor with standard semantics.
     */
    tensor(self&& M) :
             q(std::move(M.q)) { std::move(&M.counts[0],&M.counts[0] + 3,&counts[0]); };
#endif
        
    /**
     * Explicit constructor from a any type of tensor.
     */
    template <typename Tensor>
    explicit tensor(const Tensor& M, typename boost::enable_if< 
                                                 boost::mpl::and_<
                                                   is_readable_tensor<Tensor>, 
                                                   boost::mpl::not_< boost::is_same<Tensor,self> >
                                                 >, void* >::type dummy = NULL) :
             q(size<0>(M) * size<1>(M) * size<2>(M),T(0.0)) {
      counts[0] = size<0>(M); counts[1] = size<1>(M); counts[2] = size<2>(M);
      typename container_type::iterator it = q.begin();
      for(size_type i = 0; i < counts[0]; ++i)
        for(size_type j = 0; j < counts[1]; ++j)
          for(size_type k = 0; k < counts[2]; ++k, ++it)
            *it = M(i,j,k);
    };

    /**
     * Constructor from a vector of low-dim-major values.
     */
    tensor(const container_type& Q, size_type aCount0, size_type aCount1, size_type aCount2) :
             q(Q)
#ifndef BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX
           , counts{aCount0, aCount1, aCount2} { };
#else
           { counts[0] = aCount0; counts[1] = aCount1; counts[2] = aCount2; };
#endif

    /**
     * Destructor.
     */
    ~tensor() { };
    
    /**
     * The standard swap function (works with ADL).
     */
    friend void swap(self& m1, self& m2) throw() {
      using std::swap;
      swap(m1.q,m2.q);
      swap(m1.counts[0],m2.counts[0]);
      swap(m1.counts[1],m2.counts[1]);
      swap(m1.counts[2],m2.counts[2]);
    };
    
    /**
     * A swap function to swap the tensor with a container of values to fill the tensor.
     * \param m1 The tensor to swap with the container.
     * \param q2 The container that will be swapped with m1's internal container.
     * \param Count02 The count0 corresponding to q2's data.
     * \param Count12 The count1 corresponding to q2's data.
     * \param Count22 The count2 corresponding to q2's data.
     */
    friend void swap(self& m1, container_type& q2, size_type& Count02, size_type& Count12, size_type& Count22) throw() {
      using std::swap;
      swap(m1.q,q2);
      swap(m1.counts[0],Count02);
      swap(m1.counts[1],Count12);
      swap(m1.counts[2],Count22);
    };
    
    /**
     * Standard copy-assignment operator (and move-assignment operator, for C++0x). Uses the copy-and-swap (and
     * move-and-swap) idiom.
     */
    self& operator=(self rhs) {
      swap(*this, rhs);
      return *this;
    };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Tensor indexing accessor for read-write access.
     * \param i Index-0.
     * \param j Index-1.
     * \param k Index-2.
     * \return the element at the given position.
     */
    reference operator()(size_type i, size_type j, size_type k) { return q[(i * counts[1] + j) * counts[2] + k]; };
    /**
     * Tensor indexing accessor for read-only access.
     * \param i Index-0.
     * \param j Index-1.
     * \param k Index-2.
     * \return the element at the given position.
     */
    const_reference operator()(size_type i, size_type j, size_type k) const { return q[(i * counts[1] + j) * counts[2] + k]; };

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
    
    /**
     * Sets the count for the given dimension of the tensor.
     * \param number of elements along the given dimension of the tensor.
     * \tparam Dim The dimension along which the resize is to apply.
     */
    template <unsigned int Dim>
    friend 
    typename boost::enable_if<
      boost::mpl::less< boost::mpl::size_t< Dim >, boost::mpl::size_t< 3 > >,
    void >::type resize(self& m, size_type sz) {
      size_type new_counts[3] = {m.counts[0], m.counts[1], m.counts[2]};
      new_counts[Dim] = sz;
      self new_m(new_counts[0],new_counts[1],new_counts[2]);
      for(size_type i = 0; i < std::min(new_m.counts[0], m.counts[0]); ++i)
        for(size_type j = 0; j < std::min(new_m.counts[1], m.counts[1]); ++j)
          for(size_type k = 0; k < std::min(new_m.counts[2], m.counts[2]); ++k)
            new_m.q[(i * new_m.counts[1] + j) * new_m.counts[2] + k] = m.q[(i * m.counts[1] + j) * m.counts[2] + k]; 
      swap(new_m,m);
    };
    
    /**
     * Sets the count for the given dimension of the tensor.
     * \param number of elements along the given dimension of the tensor.
     * \tparam Dim The dimension along which the resize is to apply.
     */
    template <unsigned int Dim>
    friend 
    typename boost::disable_if<
      boost::mpl::less< boost::mpl::size_t< Dim >, boost::mpl::size_t< 3 > >,
    void >::type resize(self& m, size_type sz) { RK_UNUSED(m); RK_UNUSED(sz); };
    

    /**
     * Returns the allocator object of the underlying container.
     * \return the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return q.get_allocator(); };

    
/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /** 
     * Standard Assignment operator with standard semantics.
     * Strong exception-safety.
     */
    template <typename Tensor>
    self& operator =(const Tensor& M) {
      self tmp(M);
      swap(*this,tmp);
      return *this;
    };

    /** 
     * Add-and-store operator with standard semantics.
     */
    template <typename Tensor>
    typename boost::enable_if< 
      is_readable_tensor<Tensor>,
    self& >::type operator +=(const Tensor& M) {
      if((size<0>(M) != counts[0]) || 
         (size<1>(M) != counts[1]) || 
         (size<2>(M) != counts[2]))
        throw std::range_error("Tensor dimensions mismatch.");
      typename container_type::iterator it = q.begin();
      for(size_type i = 0; i < counts[0]; ++i)
        for(size_type j = 0; j < counts[1]; ++j)
          for(size_type k = 0; k < counts[2]; ++k, ++it)
            *it += M(i,j,k);
      return *this;
    };

    /** 
     * Sub-and-store operator with standard semantics.
     */
    template <typename Tensor>
    typename boost::enable_if< 
      is_readable_tensor<Tensor>,
    self& >::type operator -=(const Tensor& M) {
      if((size<0>(M) != counts[0]) || 
         (size<1>(M) != counts[1]) || 
         (size<2>(M) != counts[2]))
        throw std::range_error("Tensor dimensions mismatch.");
      typename container_type::iterator it = q.begin();
      for(size_type i = 0; i < counts[0]; ++i)
        for(size_type j = 0; j < counts[1]; ++j)
          for(size_type k = 0; k < counts[2]; ++k, ++it)
            *it -= M(i,j,k);
      return *this;
    };

    /** 
     * Scalar-multiply-and-store operator with standard semantics.
     */
    self& operator *=(const value_type& S) {
      for(typename container_type::iterator it = q.begin();it != q.end();++it)
        *it *= S;
      return *this;
    };
    
    /** 
     * General negation operator for any type of tensors. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General hi-dim-major tensor.
     */
    self operator -() const {
      self result(*this);
      typename container_type::iterator itr = result.q.begin();
      for(typename container_type::const_iterator it = q.begin(); it != q.end(); ++it, ++itr)
        *itr = -(*it);
      return result;
    };
    
    /**
     * Transposes the tensor M by a simple copy with a change of alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend tensor<T,3,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator> transpose(const self& M) {
      return tensor<T,3,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator>(M.q,M.counts[2],M.counts[1],M.counts[0]);
    };
    
    /**
     * Transposes the tensor M by simply moving the data of M into a tensor of different alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend tensor<T,3,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator> transpose_move(self& M) {
      tensor<T,3,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator> result;
      using std::swap;
      swap(result,M.q,M.counts[2],M.counts[1],M.counts[0]);
      return result;
    };

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Transposes the tensor M by simply moving the data of M into a tensor of different alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend tensor<T,3,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator> transpose(self&& M) {
      tensor<T,3,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator> result;
      using std::swap;
      swap(result,M.q,M.counts[2],M.counts[1],M.counts[0]);
      return result;
    };    
#endif

    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & std::pair<std::string, const std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int>("counts[0]",counts[0])
        & std::pair<std::string, unsigned int>("counts[1]",counts[1])
        & std::pair<std::string, unsigned int>("counts[2]",counts[2]);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      unsigned int tmp;
      A & std::pair<std::string, std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int&>("counts[0]",tmp);
      counts[0] = tmp;
      A & std::pair<std::string, unsigned int&>("counts[1]",tmp);
      counts[1] = tmp;
      A & std::pair<std::string, unsigned int&>("counts[2]",tmp);
      counts[2] = tmp;
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serialization::serializable)
  
};











/**
 * This class template specialization implements a tensor with rectangular structure
 * and hi-dim-major alignment. This class is serializable and registered to the ReaK::rtti
 * system. This tensor type is dynamically resizable.
 * 
 * Models: ReadableTensorConcept, WritableTensorConcept, FullyWritableTensorConcept, 
 * ResizableTensorConcept, and DynAllocTensorConcept.
 * 
 * \tparam T Arithmetic type of the elements of the tensor.
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T,
          typename Allocator>
class tensor<T,4,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator> : public serialization::serializable {
  public:    
    
    typedef tensor<T,4,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator> self;
    typedef Allocator allocator_type;
    
    typedef T value_type;
    typedef std::vector<value_type,allocator_type> container_type;
    
    typedef typename container_type::reference reference;
    typedef typename container_type::const_reference const_reference;
    typedef typename container_type::pointer pointer;
    typedef typename container_type::const_pointer const_pointer;

    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, order = 4);
    BOOST_STATIC_CONSTANT(tensor_alignment::tag, alignment = tensor_alignment::hi_dim_major);
    BOOST_STATIC_CONSTANT(tensor_structure::tag, structure = tensor_structure::rectangular);
    
    template <unsigned int Dim>
    struct dim {
      typedef void iterator;
      typedef void const_iterator;
      BOOST_STATIC_CONSTANT(std::size_t, static_count = 0);
    };
    
  
  private:
    std::vector<value_type,allocator_type> q; ///< Array which holds all the values of the tensor (dimension: counts[0] x counts[1] x counts[2]).
    size_type counts[4]; ///< Counts.
  public:  
    
/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default constructor: sets all to zero.
     */
    tensor(const allocator_type& aAlloc = allocator_type()) :
             q(0,value_type(),aAlloc)
#ifndef BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX
             , counts{0,0,0,0} { };
#else
             { counts[0] = 0; counts[1] = 0; counts[2] = 0; counts[3] = 0; };
#endif

    /**
     * Constructor for a sized tensor.
     */
    tensor(size_type aCount0, size_type aCount1, size_type aCount2, size_type aCount3, 
           const value_type& aFill = value_type(), const allocator_type& aAlloc = allocator_type()) :
           q(aCount0 * aCount1 * aCount2 * aCount3,aFill,aAlloc)
#ifndef BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX
           , counts{aCount0, aCount1, aCount2, aCount3} { };
#else
           { counts[0] = aCount0; counts[1] = aCount1; counts[2] = aCount2; counts[3] = aCount3; };
#endif

    /**
     * Standard Copy Constructor with standard semantics.
     */
    tensor(const self& M) :
             q(M.q) { std::copy(&M.counts[0],&M.counts[0] + 4,&counts[0]); };
         
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Standard Copy Constructor with standard semantics.
     */
    tensor(self&& M) :
             q(std::move(M.q)) { std::move(&M.counts[0],&M.counts[0] + 4,&counts[0]); };
#endif

    /**
     * Explicit constructor from a any type of tensor.
     */
    template <typename Tensor>
    explicit tensor(const Tensor& M, typename boost::enable_if< 
                                                 boost::mpl::and_<
                                                   is_readable_tensor<Tensor>, 
                                                   boost::mpl::not_< boost::is_same<Tensor,self> >
                                                 >, void* >::type dummy = NULL) :
             q(size<0>(M) * size<1>(M) * size<2>(M) * size<3>(M),T(0.0)) {
      counts[0] = size<0>(M); counts[1] = size<1>(M); counts[2] = size<2>(M); counts[3] = size<3>(M);
      typename container_type::iterator it = q.begin();
      for(size_type i = 0; i < counts[3]; ++i)
        for(size_type j = 0; j < counts[2]; ++j)
          for(size_type k = 0; k < counts[1]; ++k)
            for(size_type l = 0; l < counts[0]; ++l, ++it)
              *it = M(l,k,j,i);
    };

    /**
     * Constructor from a vector of hi-dim-major values.
     */
    tensor(const container_type& Q, size_type aCount0, size_type aCount1, size_type aCount2, size_type aCount3) :
             q(Q)
#ifndef BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX
           , counts{aCount0, aCount1, aCount2, aCount3} { };
#else
           { counts[0] = aCount0; counts[1] = aCount1; counts[2] = aCount2; counts[3] = aCount3; };
#endif

    /**
     * Destructor.
     */
    ~tensor() { };
    
    /**
     * The standard swap function (works with ADL).
     */
    friend void swap(self& m1, self& m2) throw() {
      using std::swap;
      swap(m1.q,m2.q);
      swap(m1.counts[0],m2.counts[0]);
      swap(m1.counts[1],m2.counts[1]);
      swap(m1.counts[2],m2.counts[2]);
      swap(m1.counts[3],m2.counts[3]);
    };
    
    /**
     * A swap function to swap the tensor with a container of values to fill the tensor.
     * \param m1 The tensor to swap with the container.
     * \param q2 The container that will be swapped with m1's internal container.
     * \param Count02 The count0 corresponding to q2's data.
     * \param Count12 The count1 corresponding to q2's data.
     * \param Count22 The count2 corresponding to q2's data.
     * \param Count32 The count2 corresponding to q2's data.
     */
    friend void swap(self& m1, container_type& q2, size_type& Count02, size_type& Count12, size_type& Count22, size_type& Count32) throw() {
      using std::swap;
      swap(m1.q,q2);
      swap(m1.counts[0],Count02);
      swap(m1.counts[1],Count12);
      swap(m1.counts[2],Count22);
      swap(m1.counts[3],Count32);
    };
    
    /**
     * Standard copy-assignment operator (and move-assignment operator, for C++0x). Uses the copy-and-swap (and
     * move-and-swap) idiom.
     */
    self& operator=(self rhs) {
      swap(*this, rhs);
      return *this;
    };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Tensor indexing accessor for read-write access.
     * \param i Index-0.
     * \param j Index-1.
     * \param k Index-2.
     * \param l Index-3.
     * \return the element at the given position.
     */
    reference operator()(size_type i, size_type j, size_type k, size_type l) { return q[((l * counts[2] + k) * counts[1] + j) * counts[0] + i]; };
    /**
     * Tensor indexing accessor for read-only access.
     * \param i Index-0.
     * \param j Index-1.
     * \param k Index-2.
     * \param l Index-3.
     * \return the element at the given position.
     */
    const_reference operator()(size_type i, size_type j, size_type k, size_type l) const { return q[((l * counts[2] + k) * counts[1] + j) * counts[0] + i]; };

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
    
    /**
     * Sets the count for the given dimension of the tensor.
     * \param number of elements along the given dimension of the tensor.
     * \tparam Dim The dimension along which the resize is to apply.
     */
    template <unsigned int Dim>
    friend 
    typename boost::enable_if<
      boost::mpl::less< boost::mpl::size_t< Dim >, boost::mpl::size_t< 4 > >,
    void >::type resize(self& m, size_type sz) {
      size_type new_counts[4] = {m.counts[0], m.counts[1], m.counts[2], m.counts[3]};
      new_counts[Dim] = sz;
      self new_m(new_counts[0],new_counts[1],new_counts[2],new_counts[3]);
      for(size_type l = 0; l < std::min(new_m.counts[3], m.counts[3]); ++l)
        for(size_type k = 0; k < std::min(new_m.counts[2], m.counts[2]); ++k)
          for(size_type j = 0; j < std::min(new_m.counts[1], m.counts[1]); ++j)
            for(size_type i = 0; i < std::min(new_m.counts[0], m.counts[0]); ++i)
              new_m.q[((l * new_m.counts[2] + k) * new_m.counts[1] + j) * new_m.counts[0] + i] = m.q[((l * m.counts[2] + k) * m.counts[1] + j) * m.counts[0] + i]; 
      swap(new_m,m);
    };
    
    /**
     * Sets the count for the given dimension of the tensor.
     * \param number of elements along the given dimension of the tensor.
     * \tparam Dim The dimension along which the resize is to apply.
     */
    template <unsigned int Dim>
    friend 
    typename boost::disable_if<
      boost::mpl::less< boost::mpl::size_t< Dim >, boost::mpl::size_t< 4 > >,
    void >::type resize(self& m, size_type sz) { RK_UNUSED(m); RK_UNUSED(sz); };
    
    
    /**
     * Returns the allocator object of the underlying container.
     * \return the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return q.get_allocator(); };
    
    

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /** 
     * Standard Assignment operator with standard semantics.
     * Strong exception-safety.
     */
    template <typename Tensor>
    self& operator =(const Tensor& M) {
      self tmp(M);
      swap(*this,tmp);
      return *this;
    };

    /** 
     * Add-and-store operator with standard semantics.
     */
    template <typename Tensor>
    typename boost::enable_if< 
      is_readable_tensor<Tensor>,
    self& >::type operator +=(const Tensor& M) {
      if((size<0>(M) != counts[0]) || 
         (size<1>(M) != counts[1]) || 
         (size<2>(M) != counts[2]) || 
         (size<3>(M) != counts[3]))
        throw std::range_error("Tensor dimensions mismatch.");
      typename container_type::iterator it = q.begin();
      for(size_type l = 0; l < counts[3]; ++l)
        for(size_type k = 0; k < counts[2]; ++k)
          for(size_type j = 0; j < counts[1]; ++j)
            for(size_type i = 0; i < counts[0]; ++i, ++it)
              *it += M(i,j,k,l);
      return *this;
    };

    /** 
     * Sub-and-store operator with standard semantics.
     */
    template <typename Tensor>
    typename boost::enable_if< 
      is_readable_tensor<Tensor>,
    self& >::type operator -=(const Tensor& M) {
      if((size<0>(M) != counts[0]) || 
         (size<1>(M) != counts[1]) || 
         (size<2>(M) != counts[2]) || 
         (size<2>(M) != counts[3]))
        throw std::range_error("Tensor dimensions mismatch.");
      typename container_type::iterator it = q.begin();
      for(size_type l = 0; l < counts[3]; ++l)
        for(size_type k = 0; k < counts[2]; ++k)
          for(size_type j = 0; j < counts[1]; ++j)
            for(size_type i = 0; i < counts[0]; ++i, ++it)
              *it -= M(i,j,k,l);
      return *this;
    };

    /** 
     * Scalar-multiply-and-store operator with standard semantics.
     */
    self& operator *=(const value_type& S) {
      for(typename container_type::iterator it = q.begin();it != q.end();++it)
        *it *= S;
      return *this;
    };
    
    /** 
     * General negation operator for any type of tensors. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General hi-dim-major tensor.
     */
    self operator -() const {
      self result(*this);
      typename container_type::iterator itr = result.q.begin();
      for(typename container_type::const_iterator it = q.begin(); it != q.end(); ++it, ++itr)
        *itr = -(*it);
      return result;
    };
    
    
    
    /**
     * Transposes the tensor M by a simple copy with a change of alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend tensor<T,4,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator> transpose(const self& M) {
      return tensor<T,4,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator>(M.q,M.counts[3],M.counts[2],M.counts[1],M.counts[0]);
    };
    
    /**
     * Transposes the tensor M by simply moving the data of M into a tensor of different alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend tensor<T,4,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator> transpose_move(self& M) {
      tensor<T,4,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator> result;
      using std::swap;
      swap(result,M.q,M.counts[3],M.counts[2],M.counts[1],M.counts[0]);
      return result;
    };

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Transposes the tensor M by simply moving the data of M into a tensor of different alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend tensor<T,4,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator> transpose(self&& M) {
      tensor<T,4,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator> result;
      using std::swap;
      swap(result,M.q,M.counts[3],M.counts[2],M.counts[1],M.counts[0]);
      return result;
    };    
#endif
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & std::pair<std::string, const std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int>("counts[0]",counts[0])
        & std::pair<std::string, unsigned int>("counts[1]",counts[1])
        & std::pair<std::string, unsigned int>("counts[2]",counts[2])
        & std::pair<std::string, unsigned int>("counts[3]",counts[3]);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      unsigned int tmp;
      A & std::pair<std::string, std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int&>("counts[0]",tmp);
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


/**
 * This class template specialization implements a tensor with rectangular structure
 * and low-dim-major alignment. This class is serializable and registered to the ReaK::rtti
 * system. This tensor type is dynamically resizable.
 * 
 * Models: ReadableTensorConcept, WritableTensorConcept, FullyWritableTensorConcept, 
 * ResizableTensorConcept, and DynAllocTensorConcept.
 * 
 * \tparam T Arithmetic type of the elements of the tensor.
 * \tparam Allocator Standard allocator class (as in the STL), the default is std::allocator<T>.
 */
template <typename T,
          typename Allocator>
class tensor<T,4,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator> : public serialization::serializable {
  public:    
    
    typedef tensor<T,4,tensor_structure::rectangular,tensor_alignment::low_dim_major,Allocator> self;
    typedef Allocator allocator_type;
    
    typedef T value_type;
    typedef std::vector<value_type,allocator_type> container_type;
    
    typedef typename container_type::reference reference;
    typedef typename container_type::const_reference const_reference;
    typedef typename container_type::pointer pointer;
    typedef typename container_type::const_pointer const_pointer;
  
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, order = 4);
    BOOST_STATIC_CONSTANT(tensor_alignment::tag, alignment = tensor_alignment::low_dim_major);
    BOOST_STATIC_CONSTANT(tensor_structure::tag, structure = tensor_structure::rectangular);
    
    
    template <unsigned int Dim>
    struct dim {
      typedef void iterator;
      typedef void const_iterator;
      BOOST_STATIC_CONSTANT(std::size_t, static_count = 0);
    };
    
  
  private:
    std::vector<value_type,allocator_type> q; ///< Array which holds all the values of the tensor (dimension: counts[0] x counts[1] x counts[2]).
    size_type counts[4]; ///< Counts.
  public:  
    
/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/

    /**
     * Default constructor: sets all to zero.
     */
    tensor(const allocator_type& aAlloc = allocator_type()) :
             q(0,value_type(),aAlloc)
#ifndef BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX
             , counts{0,0,0,0} { };
#else
             { counts[0] = 0; counts[1] = 0; counts[2] = 0; counts[3] = 0; };
#endif

    /**
     * Constructor for a sized tensor.
     */
    tensor(size_type aCount0, size_type aCount1, size_type aCount2, size_type aCount3, 
           const value_type& aFill = value_type(), const allocator_type& aAlloc = allocator_type()) :
           q(aCount0 * aCount1 * aCount2 * aCount3,aFill,aAlloc)
#ifndef BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX
           , counts{aCount0, aCount1, aCount2, aCount3} { };
#else
           { counts[0] = aCount0; counts[1] = aCount1; counts[2] = aCount2; counts[3] = aCount3; };
#endif

    /**
     * Standard Copy Constructor with standard semantics.
     */
    tensor(const self& M) :
             q(M.q) { std::copy(&M.counts[0],&M.counts[0] + 4,&counts[0]); };
         
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Standard Copy Constructor with standard semantics.
     */
    tensor(self&& M) :
             q(std::move(M.q)) { std::move(&M.counts[0],&M.counts[0] + 4,&counts[0]); };
#endif
        
    /**
     * Explicit constructor from a any type of tensor.
     */
    template <typename Tensor>
    explicit tensor(const Tensor& M, typename boost::enable_if< 
                                                 boost::mpl::and_<
                                                   is_readable_tensor<Tensor>, 
                                                   boost::mpl::not_< boost::is_same<Tensor,self> >
                                                 >, void* >::type dummy = NULL) :
             q(size<0>(M) * size<1>(M) * size<2>(M) * size<3>(M),T(0.0)) {
      counts[0] = size<0>(M); counts[1] = size<1>(M); counts[2] = size<2>(M); counts[3] = size<3>(M);
      typename container_type::iterator it = q.begin();
      for(size_type i = 0; i < counts[0]; ++i)
        for(size_type j = 0; j < counts[1]; ++j)
          for(size_type k = 0; k < counts[2]; ++k)
            for(size_type l = 0; l < counts[3]; ++l, ++it)
              *it = M(i,j,k,l);
    };

    /**
     * Constructor from a vector of low-dim-major values.
     */
    tensor(const container_type& Q, size_type aCount0, size_type aCount1, size_type aCount2, size_type aCount3) :
             q(Q)
#ifndef BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX
           , counts{aCount0, aCount1, aCount2, aCount3} { };
#else
           { counts[0] = aCount0; counts[1] = aCount1; counts[2] = aCount2; counts[3] = aCount3; };
#endif

    /**
     * Destructor.
     */
    ~tensor() { };
    
    /**
     * The standard swap function (works with ADL).
     */
    friend void swap(self& m1, self& m2) throw() {
      using std::swap;
      swap(m1.q,m2.q);
      swap(m1.counts[0],m2.counts[0]);
      swap(m1.counts[1],m2.counts[1]);
      swap(m1.counts[2],m2.counts[2]);
      swap(m1.counts[3],m2.counts[3]);
    };
    
    /**
     * A swap function to swap the tensor with a container of values to fill the tensor.
     * \param m1 The tensor to swap with the container.
     * \param q2 The container that will be swapped with m1's internal container.
     * \param Count02 The count0 corresponding to q2's data.
     * \param Count12 The count1 corresponding to q2's data.
     * \param Count22 The count2 corresponding to q2's data.
     * \param Count32 The count3 corresponding to q2's data.
     */
    friend void swap(self& m1, container_type& q2, size_type& Count02, size_type& Count12, size_type& Count22, size_type& Count32) throw() {
      using std::swap;
      swap(m1.q,q2);
      swap(m1.counts[0],Count02);
      swap(m1.counts[1],Count12);
      swap(m1.counts[2],Count22);
      swap(m1.counts[3],Count32);
    };
    
    /**
     * Standard copy-assignment operator (and move-assignment operator, for C++0x). Uses the copy-and-swap (and
     * move-and-swap) idiom.
     */
    self& operator=(self rhs) {
      swap(*this, rhs);
      return *this;
    };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Tensor indexing accessor for read-write access.
     * \param i Index-0.
     * \param j Index-1.
     * \param k Index-2.
     * \param l Index-3.
     * \return the element at the given position.
     */
    reference operator()(size_type i, size_type j, size_type k, size_type l) { return q[((i * counts[1] + j) * counts[2] + k) * counts[3] + l]; };
    /**
     * Tensor indexing accessor for read-only access.
     * \param i Index-0.
     * \param j Index-1.
     * \param k Index-2.
     * \param l Index-3.
     * \return the element at the given position.
     */
    const_reference operator()(size_type i, size_type j, size_type k, size_type l) const { return q[((i * counts[1] + j) * counts[2] + k) * counts[3] + l]; };

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
    
    /**
     * Sets the count for the given dimension of the tensor.
     * \param number of elements along the given dimension of the tensor.
     * \tparam Dim The dimension along which the resize is to apply.
     */
    template <unsigned int Dim>
    friend 
    typename boost::enable_if<
      boost::mpl::less< boost::mpl::size_t< Dim >, boost::mpl::size_t< 4 > >,
    void >::type resize(self& m, size_type sz) {
      size_type new_counts[3] = {m.counts[0], m.counts[1], m.counts[2], m.counts[3]};
      new_counts[Dim] = sz;
      self new_m(new_counts[0],new_counts[1],new_counts[2],new_counts[3]);
      for(size_type i = 0; i < std::min(new_m.counts[0], m.counts[0]); ++i)
        for(size_type j = 0; j < std::min(new_m.counts[1], m.counts[1]); ++j)
          for(size_type k = 0; k < std::min(new_m.counts[2], m.counts[2]); ++k)
            for(size_type l = 0; l < std::min(new_m.counts[3], m.counts[3]); ++l)
              new_m.q[((i * new_m.counts[1] + j) * new_m.counts[2] + k) * new_m.counts[3] + l] = m.q[((i * m.counts[1] + j) * m.counts[2] + k) * m.counts[3] + l]; 
      swap(new_m,m);
    };
    
    /**
     * Sets the count for the given dimension of the tensor.
     * \param number of elements along the given dimension of the tensor.
     * \tparam Dim The dimension along which the resize is to apply.
     */
    template <unsigned int Dim>
    friend 
    typename boost::disable_if<
      boost::mpl::less< boost::mpl::size_t< Dim >, boost::mpl::size_t< 4 > >,
    void >::type resize(self& m, size_type sz) { RK_UNUSED(m); RK_UNUSED(sz); };
    

    /**
     * Returns the allocator object of the underlying container.
     * \return the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return q.get_allocator(); };

    
/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /** 
     * Standard Assignment operator with standard semantics.
     * Strong exception-safety.
     */
    template <typename Tensor>
    self& operator =(const Tensor& M) {
      self tmp(M);
      swap(*this,tmp);
      return *this;
    };

    /** 
     * Add-and-store operator with standard semantics.
     */
    template <typename Tensor>
    typename boost::enable_if< 
      is_readable_tensor<Tensor>,
    self& >::type operator +=(const Tensor& M) {
      if((size<0>(M) != counts[0]) || 
         (size<1>(M) != counts[1]) || 
         (size<2>(M) != counts[2]) || 
         (size<3>(M) != counts[3]))
        throw std::range_error("Tensor dimensions mismatch.");
      typename container_type::iterator it = q.begin();
      for(size_type i = 0; i < counts[0]; ++i)
        for(size_type j = 0; j < counts[1]; ++j)
          for(size_type k = 0; k < counts[2]; ++k)
            for(size_type l = 0; l < counts[3]; ++l, ++it)
              *it += M(i,j,k,l);
      return *this;
    };

    /** 
     * Sub-and-store operator with standard semantics.
     */
    template <typename Tensor>
    typename boost::enable_if< 
      is_readable_tensor<Tensor>,
    self& >::type operator -=(const Tensor& M) {
      if((size<0>(M) != counts[0]) || 
         (size<1>(M) != counts[1]) || 
         (size<2>(M) != counts[2]) || 
         (size<3>(M) != counts[3]))
        throw std::range_error("Tensor dimensions mismatch.");
      typename container_type::iterator it = q.begin();
      for(size_type i = 0; i < counts[0]; ++i)
        for(size_type j = 0; j < counts[1]; ++j)
          for(size_type k = 0; k < counts[2]; ++k)
            for(size_type l = 0; l < counts[3]; ++l, ++it)
              *it -= M(i,j,k,l);
      return *this;
    };

    /** 
     * Scalar-multiply-and-store operator with standard semantics.
     */
    self& operator *=(const value_type& S) {
      for(typename container_type::iterator it = q.begin();it != q.end();++it)
        *it *= S;
      return *this;
    };
    
    /** 
     * General negation operator for any type of tensors. This is a default operator
     * that will be called if no better special-purpose overload exists.
     * \return General hi-dim-major tensor.
     */
    self operator -() const {
      self result(*this);
      typename container_type::iterator itr = result.q.begin();
      for(typename container_type::const_iterator it = q.begin(); it != q.end(); ++it, ++itr)
        *itr = -(*it);
      return result;
    };
    
    /**
     * Transposes the tensor M by a simple copy with a change of alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend tensor<T,4,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator> transpose(const self& M) {
      return tensor<T,4,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator>(M.q,M.counts[3],M.counts[2],M.counts[1],M.counts[0]);
    };
    
    /**
     * Transposes the tensor M by simply moving the data of M into a tensor of different alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend tensor<T,4,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator> transpose_move(self& M) {
      tensor<T,4,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator> result;
      using std::swap;
      swap(result,M.q,M.counts[3],M.counts[2],M.counts[1],M.counts[0]);
      return result;
    };

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Transposes the tensor M by simply moving the data of M into a tensor of different alignment.
     * \param M The tensor to be transposed.
     * \return The transpose of M.
     */
    friend tensor<T,4,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator> transpose(self&& M) {
      tensor<T,4,tensor_structure::rectangular,tensor_alignment::hi_dim_major,Allocator> result;
      using std::swap;
      swap(result,M.q,M.counts[3],M.counts[2],M.counts[1],M.counts[0]);
      return result;
    };    
#endif

    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & std::pair<std::string, const std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int>("counts[0]",counts[0])
        & std::pair<std::string, unsigned int>("counts[1]",counts[1])
        & std::pair<std::string, unsigned int>("counts[2]",counts[2])
        & std::pair<std::string, unsigned int>("counts[3]",counts[3]);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      unsigned int tmp;
      A & std::pair<std::string, std::vector<T>&>("q",q)
        & std::pair<std::string, unsigned int&>("counts[0]",tmp);
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








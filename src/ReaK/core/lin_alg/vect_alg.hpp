/**
 * \file vect_alg.hpp
 * 
 * This library provides classes to represent fixed and variable sized vectors. These 
 * vectors can be used in linear algebra expressions, with standard vector algebra semantics.
 * Generally, vectors are considered to be column-vectors (but pre-multiplication to a matrix
 * turns assumes them to a row-vector for the duration of the operation). The interface 
 * of the vectors are also generally compatible with STL containers (and the underlying 
 * container for the variable-size vector is a std::vector class template).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2011
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

#ifndef REAK_VECT_ALG_HPP
#define REAK_VECT_ALG_HPP


#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <limits>

#include "base/defs.hpp"
#include "base/serializable.hpp"
#include "rtti/so_register_type.hpp"
#include "rtti/typed_primitives.hpp"

#include "vect_concepts.hpp"

#include <boost/config.hpp>
#include <boost/static_assert.hpp>

//#include "rk_math_traits.hpp"

namespace ReaK {

/**
 * This class implements a fixed-size templated vector class which holds components of primitive type.
 */
template <class T, unsigned int Size>
class vect : public serialization::serializable {
  public:
    typedef vect<T,Size> self;
    
    typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef std::allocator<T> allocator_type; //just in case it is cast to variable-length vector.
  
    typedef pointer iterator;
    typedef const_pointer const_iterator;
  
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = Size);
    
    T q[Size]; /**< Components of the vector.*/

    //typedef array_order<1> order;
    
    /**
     * Returns the size of the vector.
     */
    size_type size() const { return Size; };
    /**
     * Returns the max-size of the vector.
     */
    size_type max_size() const { return Size; };
    /**
     * Returns the capacity of the vector.
     */
    size_type capacity() const { return Size; };
    /**
     * Resizes the vector.
     */
    void resize(size_type sz, T c = T()) const { };
    /**
     * Checks if the vector is empty.
     */
    bool empty() const { return false; };
    /**
     * Reserve a capacity for the vector.
     */
    void reserve(size_type sz) const { };
    
    /**
     * Returns an iterator to the first element of the vector.
     */
    iterator begin() { return q; };
    /**
     * Returns a const-iterator to the first element of the vector.
     */
    const_iterator begin() const { return q; };
    /**
     * Returns an iterator to the one-past-last element of the vector.
     */
    iterator end() { return q + Size; };
    /**
     * Returns a const-iterator to the one-past-last element of the vector.
     */
    const_iterator end() const { return q + Size; };

/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/
    /**
     * Default constructor: sets all to zero.
     * \test PASSED
     */
    vect() {
      for(size_type i = 0;i < Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor from an array of values of type T.
     * \test PASSED
     */
    vect(const_pointer Q) {
      for(size_type i=0;i < Size;++i)
	q[i] = Q[i];
      return;
    };

    /**
     * Standard Copy Constructor with standard semantics.
     * \test PASSED
     */
    vect(const self& V) {
      for(size_type i=0;i < Size;++i)
	q[i] = V.q[i];
      return;
    };

    /**
     * Destructor.
     * \test PASSED
     */
    ~vect() { };

    /**
     * Constructor for 1 values.
     * \test PASSED
     */
    explicit vect(const_reference Q1) {
      BOOST_STATIC_ASSERT(Size >= 1);
      q[0] = Q1;
      for(size_type i=1;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 2 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2) {
      BOOST_STATIC_ASSERT(Size >= 2);
      q[0] = Q1;
      q[1] = Q2;
      for(size_type i=2;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 3 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3) {
      BOOST_STATIC_ASSERT(Size >= 3);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      for(size_type i=3;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 4 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4) {
      BOOST_STATIC_ASSERT(Size >= 4);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      for(size_type i=4;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 5 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,
	 const_reference Q4,const_reference Q5) {
      BOOST_STATIC_ASSERT(Size >= 5);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      for(size_type i=5;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 6 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,
	 const_reference Q4,const_reference Q5,const_reference Q6) {
      BOOST_STATIC_ASSERT(Size >= 6);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      for(size_type i=6;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 7 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,
	 const_reference Q5,const_reference Q6,const_reference Q7) {
      BOOST_STATIC_ASSERT(Size >= 7);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      for(size_type i=7;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 8 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,
	 const_reference Q5,const_reference Q6,const_reference Q7,const_reference Q8) {
      BOOST_STATIC_ASSERT(Size >= 8);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      for(size_type i=8;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 9 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	 const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9) {
      BOOST_STATIC_ASSERT(Size >= 9);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      for(size_type i=9;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 10 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	 const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10) {
      BOOST_STATIC_ASSERT(Size >= 10);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      for(size_type i=10;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 11 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	 const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	 const_reference Q11) {
      BOOST_STATIC_ASSERT(Size >= 11);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      for(size_type i=11;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 12 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	 const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	 const_reference Q11, const_reference Q12) {
      BOOST_STATIC_ASSERT(Size >= 12);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      for(size_type i=12;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 13 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	 const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	 const_reference Q11, const_reference Q12, const_reference Q13) {
      BOOST_STATIC_ASSERT(Size >= 13);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      for(size_type i=13;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 14 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	 const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	 const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14) {
      BOOST_STATIC_ASSERT(Size >= 14);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      q[13] = Q14;
      for(size_type i=14;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 15 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	 const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	 const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15) {
      BOOST_STATIC_ASSERT(Size >= 15);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      q[13] = Q14;
      q[14] = Q15;
      for(size_type i=15;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 16 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	 const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	 const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
	 const_reference Q16) {
      BOOST_STATIC_ASSERT(Size >= 16);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      q[13] = Q14;
      q[14] = Q15;
      q[15] = Q16;
      for(size_type i=16;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 17 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	 const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	 const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
	 const_reference Q16, const_reference Q17) {
      BOOST_STATIC_ASSERT(Size >= 17);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      q[13] = Q14;
      q[14] = Q15;
      q[15] = Q16;
      q[16] = Q17;
      for(size_type i=17;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 18 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	 const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	 const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
	 const_reference Q16, const_reference Q17, const_reference Q18) {
      BOOST_STATIC_ASSERT(Size >= 18);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      q[13] = Q14;
      q[14] = Q15;
      q[15] = Q16;
      q[16] = Q17;
      q[17] = Q18;
      for(size_type i=18;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 19 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	 const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	 const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
	 const_reference Q16, const_reference Q17, const_reference Q18, const_reference Q19) {
      BOOST_STATIC_ASSERT(Size >= 19);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      q[13] = Q14;
      q[14] = Q15;
      q[15] = Q16;
      q[16] = Q17;
      q[17] = Q18;
      q[18] = Q19;
      for(size_type i=19;i<Size;++i)
	q[i] = value_type();
      return;
    };

    /**
     * Constructor for 20 values.
     * \test PASSED
     */
    vect(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	 const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	 const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
	 const_reference Q16, const_reference Q17, const_reference Q18, const_reference Q19, const_reference Q20) {
      BOOST_STATIC_ASSERT(Size >= 20);
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      q[13] = Q14;
      q[14] = Q15;
      q[15] = Q16;
      q[16] = Q17;
      q[17] = Q18;
      q[18] = Q19;
      q[19] = Q20;
      for(size_type i=20;i<Size;++i)
	q[i] = value_type();
      return;
    };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Array indexing operator, accessor for read/write.
     * \test PASSED
     */
    reference operator [](size_type i) {
      return q[i];
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator [](size_type i) const {
      return q[i];
    };

    /**
     * Array indexing operator, accessor for read/write.
     * \test PASSED
     */
    reference operator ()(size_type i) {
      return q[i];
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator ()(size_type i) const {
      return q[i];
    };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /**
     * Standard assignment operator.
     * \test PASSED
     */
    self& operator =(const self& V) {
      for(size_type i=0;i<Size;++i)
	q[i] = V[i];
      return *this;
    };
    
    /**
     * Standard assignment operator.
     * \test PASSED
     */
    template <typename Vector>
    typename boost::enable_if_c< is_readable_vector<Vector>::value &&
                                 !boost::is_same<Vector,self>::value,
    self& >::type operator =(const Vector& V) {
      if(Size != V.size())
        throw std::range_error("Vector size mismatch.");
      for(size_type i=0; i < Size; ++i)
	q[i] = V[i];
      return *this;
    };

    /**
     * Standard add-and-store operator.
     * \test PASSED
     */
    self& operator +=(const self& V) {
      for(size_type i=0;i< Size;++i)
	q[i] += V[i];
      return *this;
    };

    /**
     * Standard sub-and-store operator.
     * \test PASSED
     */
    self& operator -=(const self& V) {
      for(size_type i=0;i<Size;++i)
	q[i] -= V[i];
      return *this;
    };

    /**
     * Scalar multiply-and-store operator for gain.
     * \test PASSED
     */
    template <typename U>
    self& operator *=(const U& S) {
      for(size_type i=0;i<Size;++i)
	q[i] = value_type(q[i] * S);
      return *this;
    };

    /**
     * Scalar divide-and-store operator for gain.
     * \test PASSED
     */
    template <typename U>
    self& operator /=(const U& S) {
      for(size_type i=0;i<Size;++i)
	q[i] = value_type(q[i] / S);
      return *this;
    };

/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      for(size_type i=0;i<Size;++i)
        A & std::pair<std::string, typename ReaK::rtti::get_type_id<T>::save_type >("q",q[i]);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      for(size_type i=0;i<Size;++i)
        A & std::pair<std::string, typename ReaK::rtti::get_type_id<T>::load_type >("q",q[i]);
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serialization::serializable)

};

namespace rtti {

template <typename T, unsigned int Size>
struct get_type_id< vect<T,Size> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000011);
  static std::string type_name() { return "vect"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const serialization::serializable& save_type;
  typedef serialization::serializable& load_type;
};

template <typename T, unsigned int Size, typename Tail>
struct get_type_info< vect<T,Size>, Tail > {
  typedef detail::type_id< vect<T,Size> , typename get_type_info<T, 
                                                   get_type_info<boost::mpl::integral_c<unsigned int,Size> , Tail> >::type> type;
  static std::string type_name() { return get_type_id< vect<T,Size> >::type_name() + "<" + get_type_id<T>::type_name() + "," + get_type_id< boost::mpl::integral_c<unsigned int,Size> >::type_name() + ">" + "," + Tail::type_name(); };
};

};


template <typename T, unsigned int Size>
struct is_readable_vector< vect<T,Size> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_vector< vect<T,Size> > type;
};

template <typename T, unsigned int Size>
struct is_writable_vector< vect<T,Size> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_writable_vector< vect<T,Size> > type;
};

template <typename T, unsigned int Size>
struct is_resizable_vector< vect<T,Size> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_vector< vect<T,Size> > type;
};


template <typename T, unsigned int Size>
struct has_allocator_vector< vect<T,Size> > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef has_allocator_vector< vect<T,Size> > type;
};







/*******************************************************************************
                         Basic Constructors
*******************************************************************************/

template <typename T>
vect<T,1> make_vect(const T& Q1) {
  return vect<T,1>(Q1);
};

template <typename T>
vect<T,2> make_vect(const T& Q1,const T& Q2) {
  return vect<T,2>(Q1,Q2);
};

template <typename T>
vect<T,3> make_vect(const T& Q1,const T& Q2,const T& Q3) {
  return vect<T,3>(Q1,Q2,Q3);
};

template <typename T>
vect<T,4> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4) {
  return vect<T,4>(Q1,Q2,Q3,Q4);
};

template <typename T>
vect<T,5> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5) {
  return vect<T,5>(Q1,Q2,Q3,Q4,Q5);
};

template <typename T>
vect<T,6> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6) {
  return vect<T,6>(Q1,Q2,Q3,Q4,Q5,Q6);
};

template <typename T>
vect<T,7> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7) {
  return vect<T,7>(Q1,Q2,Q3,Q4,Q5,Q6,Q7);
};

template <typename T>
vect<T,8> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8) {
  return vect<T,8>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8);
};

template <typename T>
vect<T,9> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8,const T& Q9) {
  return vect<T,9>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9);
};

template <typename T>
vect<T,10> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8,const T& Q9,const T& Q10) {
  return vect<T,10>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10);
};

template <typename T>
vect<T,11> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8,const T& Q9,const T& Q10,const T& Q11) {
  return vect<T,11>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11);
};

template <typename T>
vect<T,12> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8,const T& Q9,const T& Q10,const T& Q11,const T& Q12) {
  return vect<T,12>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12);
};

template <typename T>
vect<T,13> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8,const T& Q9,const T& Q10,const T& Q11,const T& Q12,const T& Q13) {
  return vect<T,13>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13);
};

template <typename T>
vect<T,14> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8,const T& Q9,const T& Q10,const T& Q11,const T& Q12,const T& Q13,const T& Q14) {
  return vect<T,14>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14);
};

template <typename T>
vect<T,15> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8,const T& Q9,const T& Q10,const T& Q11,const T& Q12,const T& Q13,const T& Q14,const T& Q15) {
  return vect<T,15>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15);
};

template <typename T>
vect<T,16> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8,const T& Q9,const T& Q10,const T& Q11,const T& Q12,const T& Q13,const T& Q14,const T& Q15,const T& Q16) {
  return vect<T,16>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16);
};

template <typename T>
vect<T,17> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8,const T& Q9,const T& Q10,const T& Q11,const T& Q12,const T& Q13,const T& Q14,const T& Q15,const T& Q16,const T& Q17) {
  return vect<T,17>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16,Q17);
};

template <typename T>
vect<T,18> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8,const T& Q9,const T& Q10,const T& Q11,const T& Q12,const T& Q13,const T& Q14,const T& Q15,const T& Q16,const T& Q17,const T& Q18) {
  return vect<T,18>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16,Q17,Q18);
};

template <typename T>
vect<T,19> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8,const T& Q9,const T& Q10,const T& Q11,const T& Q12,const T& Q13,const T& Q14,const T& Q15,const T& Q16,const T& Q17,const T& Q18,const T& Q19) {
  return vect<T,19>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16,Q17,Q18,Q19);
};

template <typename T>
vect<T,20> make_vect(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8,const T& Q9,const T& Q10,const T& Q11,const T& Q12,const T& Q13,const T& Q14,const T& Q15,const T& Q16,const T& Q17,const T& Q18,const T& Q19,const T& Q20) {
  return vect<T,20>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16,Q17,Q18,Q19,Q20);
};




/*******************************************************************************
                         Basic Functions
*******************************************************************************/


/**
 * Square magnitude of the vector.
 * \test PASSED
 */
template <typename T, unsigned int Size>
T norm_sqr(const vect<T,Size>& v) {
  T sum(0);
  for(unsigned int i=0;i<Size;++i)
    sum += v.q[i]*v.q[i];
  return sum;
};

/**
 * Magnitude of the vector.
 * \test PASSED
 */
template <typename T, unsigned int Size>
T norm(const vect<T,Size>& v) {
  return std::sqrt( norm_sqr(v) );
};


/**
 * Unit vector in the same direction.
 * \test PASSED
 */
template <typename T, unsigned int Size>
vect<T,Size> unit(const vect<T,Size>& v) {
  T m = norm(v);
  T result[Size];
  for(unsigned int i=0; i<Size; ++i)
    result[i] = v.q[i] / m;
  return vect<T,Size>(result);
};

/**
 * Checks if two vectors are colinear.
 * \test PASSED
 */
template <typename T, unsigned int Size>
bool colinear(const vect<T,Size>& v1, const vect<T,Size>& v2) {
  T tmp_mag1 = norm(v1);
  T tmp_mag2 = norm(v2);
  T tmp_comb = norm(v1 + v2);
  return (((tmp_mag1 + tmp_mag2) * (T(1.0) - std::numeric_limits<T>::epsilon()) <= tmp_comb) || 
          (std::fabs(tmp_mag1 - tmp_mag2) * (T(1.0) + std::numeric_limits<T>::epsilon()) >= tmp_comb));
};


/*******************************************************************************
                         Basic Operators
*******************************************************************************/

/**
 * Add two vectors.
 * \test PASSED
 */
template <typename T, unsigned int Size>
vect<T,Size> operator +(const vect<T,Size>& v1, const vect<T,Size>& v2) {
  T result[Size];
  for(unsigned int i=0;i<Size;++i)
    result[i] = v1[i] + v2[i];
  return vect<T,Size>(result);
};

/**
 * Invert the vector.
 * \test PASSED
 */
template <typename T, unsigned int Size>
vect<T,Size> operator -(const vect<T,Size>& v) {
  T result[Size];
  for(unsigned int i=0;i<Size;++i)
    result[i] = -v[i];
  return vect<T,Size>(result);
};

/**
 * Sub two vectors.
 * \test PASSED
 */
template <typename T, unsigned int Size>
vect<T,Size> operator -(const vect<T,Size>& v1, const vect<T,Size>& v2) {
  T result[Size];
  for(unsigned int i=0;i<Size;++i)
    result[i] = v1[i] - v2[i];
  return vect<T,Size>(result);
};

/**
 * Dot Product.
 * \test PASSED
 */
template <typename T, unsigned int Size>
T operator *(const vect<T,Size>& v1, const vect<T,Size>& v2) {
  T result(0);
  for(unsigned int i=0;i<Size;++i)
    result += v1[i] * v2[i];
  return result;
};

/**
 * Scalar-vector product.
 * \test PASSED
 */
template <typename T, unsigned int Size, typename U>
vect<T,Size> operator *(const vect<T,Size>& V, const U& S) {
  T result[Size];
  for(unsigned int i=0;i<Size;++i)
    result[i] = V[i] * S;
  return vect<T,Size>(result);
};

/**
 * Scalar-vector product.
 * \test PASSED
 */
template <typename T, unsigned int Size, typename U>
vect<T,Size> operator *(const U& S,const vect<T,Size>& V) {
  T result[Size];
  for(unsigned int i=0;i<Size;++i)
    result[i] = V[i] * S;
  return vect<T,Size>(result);
};

/**
 * Scalar-vector division.
 * \test PASSED
 */
template <typename T, unsigned int Size, typename U>
vect<T,Size> operator /(const vect<T,Size>& V, const U& S) {
  T result[Size];
  for(unsigned int i=0;i<Size;++i)
    result[i] = V[i] / S;
  return vect<T,Size>(result);
};

/**
 * Sub two vectors. For functional interfaces.
 * \test PASSED
 */
template <typename T, unsigned int Size>
vect<T,Size> diff(const vect<T,Size>& v1, const vect<T,Size>& v2) {
  return v1 - v2;
};

/**
 * Add two vectors. For functional interfaces.
 * \test PASSED
 */
template <typename T, unsigned int Size>
vect<T,Size> add(const vect<T,Size>& v1, const vect<T,Size>& v2) {
  return v1 + v2;
};


/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

/**
 * Equality Comparison operator, component-wise.
 * \test PASSED
 */
template <typename T, unsigned int Size>
bool operator ==(const vect<T,Size>& v1, const vect<T,Size>& v2) {
  for(unsigned int i=0;i<Size;++i)
    if(v1[i] != v2[i])
      return false;
  return true;
};

/**
 * Inequality Comparison operator, component-wise.
 * \test PASSED
 */
template <typename T, unsigned int Size>
bool operator !=(const vect<T,Size>& v1, const vect<T,Size>& v2) {
  for(unsigned int i=0;i<Size;++i)
    if(v1[i] != v2[i])
      return true;
  return false;
};

/**
 * Greater-than Comparison operator, Euclidean norm.
 * \test PASSED
 */
template <typename T, unsigned int Size>
bool operator >(const vect<T,Size>& v1, const vect<T,Size>& v2) {
  return (norm_sqr(v1) > norm_sqr(v2));
};

/**
 * Smaller-than Comparison operator, Euclidean norm.
 * \test PASSED
 */
template <typename T, unsigned int Size>
bool operator <(const vect<T,Size>& v1, const vect<T,Size>& v2) {
  return (norm_sqr(v1) < norm_sqr(v2));
};

/**
 * Greater-or-equal Comparison operator, Euclidean norm.
 * \test PASSED
 */
template <typename T, unsigned int Size>
bool operator >=(const vect<T,Size>& v1, const vect<T,Size>& v2) {
  return (norm_sqr(v1) >= norm_sqr(v2));
};

/**
 * Smaller-or-equal Comparison operator, Euclidean norm.
 * \test PASSED
 */
template <typename T, unsigned int Size>
bool operator <=(const vect<T,Size>& v1, const vect<T,Size>& v2) {
  return (norm_sqr(v1) <= norm_sqr(v2));
};
    


/*******************************************************************************
                  Streaming, Serializing and RTTI Functionalities
*******************************************************************************/

/**
 * Prints a fixed-size vector to an output stream as "(v1; v2; v3; ..; vN)".
 * \test PASSED
 */
template <typename T,unsigned int Size>
std::ostream& operator <<(std::ostream& out_stream,const vect<T,Size>& V) {
  out_stream << "(";
  if(Size > 0)
    out_stream << V[0];
  for(unsigned int i=1;i<Size;++i)
    out_stream << "; " << V[i];
  return out_stream << ")";
};





/*******************************************************************************
                         Special Vector Products / Operators
*******************************************************************************/


/**
 * 2D Cross-Product.
 * \test PASSED
 */
template <class T>
T operator %(const vect<T,2>& V1,const vect<T,2>& V2) {
  return V1[0] * V2[1] - V1[1] * V2[0];
};

/**
 * 2D Cross-Product.
 * \test PASSED
 */
template <class T>
vect<T,2> operator %(const T& S, const vect<T,2>& V) {
  vect<T,2> result;
  result[0] = -V[1]*S;
  result[1] =  V[0]*S;
  return result;
};

/**
 * 2D Cross-Product.
 * \test PASSED
 */
template <class T>
vect<T,2> operator %(const vect<T,2>& V,const T& S) {
  vect<T,2> result;
  result[0] =  V[1]*S;
  result[1] = -V[0]*S;
  return result;
};

/**
 * 3D Cross-Product.
 * \test PASSED
 */
template <class T>
vect<T,3> operator %(const vect<T,3>& V1 , const vect<T,3>& V2) {
  vect<T,3> result;
  result[0] = V1[1]*V2[2] - V1[2]*V2[1];
  result[1] = V1[2]*V2[0] - V1[0]*V2[2];
  result[2] = V1[0]*V2[1] - V1[1]*V2[0];
  return result;
};







/**
 * This class implements a variable-size templated vector class which holds components of dimensional quantities.
 */
template <typename T, typename Allocator = std::allocator<T> >
class vect_n : public serialization::serializable {
  public:
    
    typedef vect_n<T,Allocator> self;
    
    typedef T value_type;
    typedef typename std::vector<T,Allocator>::reference reference;
    typedef typename std::vector<T,Allocator>::const_reference const_reference;
    typedef typename std::vector<T,Allocator>::pointer pointer;
    typedef typename std::vector<T,Allocator>::const_pointer const_pointer;
    typedef Allocator allocator_type;
  
    typedef typename std::vector<T,Allocator>::iterator iterator;
    typedef typename std::vector<T,Allocator>::const_iterator const_iterator;
  
    typedef typename std::vector<T,Allocator>::size_type size_type;
    typedef typename std::vector<T,Allocator>::difference_type difference_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);

    std::vector<T,Allocator> q; /**< Components of the vector. */
    
    /**
     * Returns the size of the vector.
     */
    size_type size() const { return q.size(); };
    /**
     * Returns the max-size of the vector.
     */
    size_type max_size() const { return q.max_size(); };
    /**
     * Returns the capacity of the vector.
     */
    size_type capacity() const { return q.capacity(); };
    /**
     * Resizes the vector.
     */
    void resize(size_type sz, T c = T()) { q.resize(sz,c); };
    /**
     * Checks if the vector is empty.
     */
    bool empty() const { return q.empty(); };
    /**
     * Reserve a capacity for the vector.
     */
    void reserve(size_type sz) const { q.reserve(sz); };
    
    /**
     * Returns an iterator to the first element of the vector.
     */
    iterator begin() { return q.begin(); };
    /**
     * Returns a const-iterator to the first element of the vector.
     */
    const_iterator begin() const { return q.begin(); };
    /**
     * Returns an iterator to the one-past-last element of the vector.
     */
    iterator end() { return q.end(); };
    /**
     * Returns a const-iterator to the one-past-last element of the vector.
     */
    const_iterator end() const { return q.end(); };

/*******************************************************************************
                         Constructors / Destructors
*******************************************************************************/
    /**
     * Default constructor.
     * \test PASSED
     */
    vect_n(const allocator_type& aAlloc = allocator_type()) : q(aAlloc) { };

    /**
     * Constructor for a given size or dimension of vector.
     * \test PASSED
     */
    explicit vect_n(size_type aSize,const_reference Fill = value_type(),const allocator_type& aAlloc = allocator_type()) : q(aSize,Fill,aAlloc) { };

    /**
     * Constructor from an array of values of type "value_type".
     * \test PASSED
     */
    template <typename OtherAllocator>
    explicit vect_n(const std::vector<value_type,OtherAllocator>& Q,const allocator_type& aAlloc = allocator_type()) : q(Q.begin(),Q.end(),aAlloc) { };
    
    /**
     * Constructor from a forward iterator of values of type "value_type".
     * \test PASSED
     */
    template <typename InputIter>
    explicit vect_n(InputIter first,InputIter last,const allocator_type& aAlloc = allocator_type()) : q(first,last,aAlloc) { };

    /**
     * Standard Copy Constructor with standard semantics.
     * \test PASSED
     */
    vect_n(const self& V) : q(V.begin(),V.end(),V.get_allocator()) {};
    
    /**
     * Allocator-agnostic Copy Constructor with standard semantics.
     * \test PASSED
     */
    template <typename OtherAllocator>
    vect_n(const vect_n<value_type,OtherAllocator>& V,const allocator_type& aAlloc = allocator_type()) : q(V.begin(),V.end(),aAlloc) {};

    /**
     * Constructor from a fixed-length vector.
     */
    template <unsigned int Size>
    vect_n(const vect<value_type,Size>& V,const allocator_type& aAlloc = allocator_type()) : q(V.begin(),V.end(),aAlloc) { };

    /**
     * Destructor.
     * \test PASSED
     */
    ~vect_n() { };

    /**
     * Constructor for 3 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3) : q(3) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      return;
    };

    /**
     * Constructor for 4 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4) : q(4) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      return;
    };

    /**
     * Constructor for 5 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,
	   const_reference Q4,const_reference Q5) : q(5) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      return;
    };

    /**
     * Constructor for 6 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,
	   const_reference Q4,const_reference Q5,const_reference Q6) : q(6) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      return;
    };

    /**
     * Constructor for 7 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,
	   const_reference Q5,const_reference Q6,const_reference Q7) : q(7) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      return;
    };

    /**
     * Constructor for 8 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,
	   const_reference Q5,const_reference Q6,const_reference Q7,const_reference Q8) : q(8) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      return;
    };

    /**
     * Constructor for 9 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	   const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9) : q(9) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      return;
    };

    /**
     * Constructor for 10 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	   const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10) : q(10) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      return;
    };

    /**
     * Constructor for 11 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	   const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	   const_reference Q11) : q(11) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      return;
    };

    /**
     * Constructor for 12 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	   const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	   const_reference Q11,const_reference Q12) : q(12) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      return;
    };

    /**
     * Constructor for 13 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	   const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	   const_reference Q11,const_reference Q12,const_reference Q13) : q(13) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      return;
    };

    /**
     * Constructor for 14 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	   const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	   const_reference Q11,const_reference Q12,const_reference Q13,const_reference Q14) : q(14) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      q[13] = Q14;
      return;
    };

    /**
     * Constructor for 15 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	   const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	   const_reference Q11,const_reference Q12,const_reference Q13,const_reference Q14,const_reference Q15) : q(15) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      q[13] = Q14;
      q[14] = Q15;
      return;
    };

    /**
     * Constructor for 16 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	   const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	   const_reference Q11,const_reference Q12,const_reference Q13,const_reference Q14,const_reference Q15,
	   const_reference Q16) : q(16) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      q[13] = Q14;
      q[14] = Q15;
      q[15] = Q16;
      return;
    };

    /**
     * Constructor for 17 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	   const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	   const_reference Q11,const_reference Q12,const_reference Q13,const_reference Q14,const_reference Q15,
	   const_reference Q16,const_reference Q17) : q(17) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      q[13] = Q14;
      q[14] = Q15;
      q[15] = Q16;
      q[16] = Q17;
      return;
    };

    /**
     * Constructor for 18 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	   const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	   const_reference Q11,const_reference Q12,const_reference Q13,const_reference Q14,const_reference Q15,
	   const_reference Q16,const_reference Q17,const_reference Q18) : q(18) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      q[13] = Q14;
      q[14] = Q15;
      q[15] = Q16;
      q[16] = Q17;
      q[17] = Q18;
      return;
    };

    /**
     * Constructor for 19 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	   const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	   const_reference Q11,const_reference Q12,const_reference Q13,const_reference Q14,const_reference Q15,
	   const_reference Q16,const_reference Q17,const_reference Q18,const_reference Q19) : q(19) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      q[13] = Q14;
      q[14] = Q15;
      q[15] = Q16;
      q[16] = Q17;
      q[17] = Q18;
      q[18] = Q19;
      return;
    };

    /**
     * Constructor for 20 values.
     * \test PASSED
     */
    vect_n(const_reference Q1,const_reference Q2,const_reference Q3,const_reference Q4,const_reference Q5,
	   const_reference Q6,const_reference Q7,const_reference Q8,const_reference Q9,const_reference Q10,
	   const_reference Q11,const_reference Q12,const_reference Q13,const_reference Q14,const_reference Q15,
	   const_reference Q16,const_reference Q17,const_reference Q18,const_reference Q19,const_reference Q20) : q(20) {
      q[0] = Q1;
      q[1] = Q2;
      q[2] = Q3;
      q[3] = Q4;
      q[4] = Q5;
      q[5] = Q6;
      q[6] = Q7;
      q[7] = Q8;
      q[8] = Q9;
      q[9] = Q10;
      q[10] = Q11;
      q[11] = Q12;
      q[12] = Q13;
      q[13] = Q14;
      q[14] = Q15;
      q[15] = Q16;
      q[16] = Q17;
      q[17] = Q18;
      q[18] = Q19;
      q[19] = Q20;
      return;
    };

/*******************************************************************************
                         Accessors and Methods
*******************************************************************************/

    /**
     * Array indexing operator, accessor for read/write. <
     * \test PASSED
     */
    reference operator [](size_type i) {
      return q[i];
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator [](size_type i) const {
      return q[i];
    };
    
    /**
     * Array indexing operator, accessor for read/write. <
     * \test PASSED
     */
    reference operator ()(size_type i) {
      return q[i];
    };

    /**
     * Array indexing operator, accessor for read only.
     * \test PASSED
     */
    const_reference operator ()(size_type i) const {
      return q[i];
    };
    
    /**
     * Returns the allocator object of the underlying container.
     */
    allocator_type get_allocator() const { return q.get_allocator(); };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

    /**
     * Standard assignment operator.
     * \test PASSED
     */
    self& operator =(const self& V) {
      q.assign(V.q.begin(),V.q.end());
      return *this;
    };
    
    /**
     * Standard assignment operator.
     * \test PASSED
     */
    template <typename Vector>
    typename boost::enable_if_c< is_readable_vector<Vector>::value &&
                                 !boost::is_same<Vector,self>::value,
    self& >::type operator =(const Vector& V) {
      q.resize(V.size());
      for(size_type i=0; i < q.size(); ++i)
	q[i] = V[i];
      return *this;
    };

    /**
     * Standard add-and-store operator.
     * \test PASSED
     */
    self& operator +=(const self& V) {
      if(q.size() != V.q.size())
        throw std::range_error("Vector size mismatch.");
      for(size_type i=0;i<q.size();++i)
	q[i] += V[i];
      return *this;
    };

    /**
     * Standard sub-and-store operator.
     * \test PASSED
     */
    self& operator -=(const self& V) {
      if(q.size() != V.q.size())
        throw std::range_error("Vector size mismatch.");
      for(size_type i=0;i<q.size();++i)
	q[i] -= V[i];
      return *this;
    };

    /**
     * Scalar multiply-and-store operator for gain.
     * \test PASSED
     */
    self& operator *=(const_reference S) {
      for(size_type i=0;i<q.size();++i)
	q[i] *= S;
      return *this;
    };

    /**
     * Scalar divide-and-store operator for gain.
     * \test PASSED
     */
    self& operator /=(const_reference S) {
      for(size_type i=0;i<q.size();++i)
	q[i] /= S;
      return *this;
    };

/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & std::pair<std::string, const std::vector<T>&>("q",q);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & std::pair<std::string, std::vector<T>&>("q",q);
    };
    
    RK_RTTI_REGISTER_CLASS_1BASE(self,1,serialization::serializable)

};


namespace rtti {

template <typename T,typename Allocator>
struct get_type_id< vect_n<T,Allocator> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x00000010);
  static std::string type_name() { return "vect_n"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const serialization::serializable& save_type;
  typedef serialization::serializable& load_type;
};

template <typename T, typename Allocator, typename Tail>
struct get_type_info< vect_n<T,Allocator>, Tail > {
  typedef detail::type_id< vect_n<T,Allocator> , typename get_type_info<T, Tail>::type > type;
  static std::string type_name() { return get_type_id< vect_n<T,Allocator> >::type_name() + "<" + get_type_id<T>::type_name() + ">" + "," + Tail::type_name(); };
};

};



template <typename T,typename Allocator>
struct is_readable_vector< vect_n<T,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_vector< vect_n<T,Allocator> > type;
};

template <typename T,typename Allocator>
struct is_writable_vector< vect_n<T,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_writable_vector< vect_n<T,Allocator> > type;
};

template <typename T,typename Allocator>
struct is_resizable_vector< vect_n<T,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_resizable_vector< vect_n<T,Allocator> > type;
};

template <typename T,typename Allocator>
struct has_allocator_vector< vect_n<T,Allocator> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef has_allocator_vector< vect_n<T,Allocator> > type;
};




/*******************************************************************************
                         Basic Constructors
*******************************************************************************/

template <typename T>
vect_n<T> make_vect_n(const T& Q1,const T& Q2,const T& Q3) {
  return vect_n<T>(Q1,Q2,Q3);
};

template <typename T>
vect_n<T> make_vect_n(const T& Q1,const T& Q2,const T& Q3,const T& Q4) {
  return vect_n<T>(Q1,Q2,Q3,Q4);
};

template <typename T>
vect_n<T> make_vect_n(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5) {
  return vect_n<T>(Q1,Q2,Q3,Q4,Q5);
};

template <typename T>
vect_n<T> make_vect_n(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6) {
  return vect_n<T>(Q1,Q2,Q3,Q4,Q5,Q6);
};

template <typename T>
vect_n<T> make_vect_n(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7) {
  return vect_n<T>(Q1,Q2,Q3,Q4,Q5,Q6,Q7);
};

template <typename T>
vect_n<T> make_vect_n(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8) {
  return vect_n<T>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8);
};

template <typename T>
vect_n<T> make_vect_n(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8,const T& Q9) {
  return vect_n<T>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9);
};

template <typename T>
vect_n<T> make_vect_n(const T& Q1,const T& Q2,const T& Q3,const T& Q4,const T& Q5,const T& Q6,const T& Q7,const T& Q8,const T& Q9,const T& Q10) {
  return vect_n<T>(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10);
};


/*******************************************************************************
                         Basic Functions
*******************************************************************************/

/**
 * Magnitude of the vector.
 * \test PASSED
 */
template <typename T, typename Allocator>
T norm(const vect_n<T,Allocator>& v) {
  return std::sqrt( norm_sqr(v) );
};

/**
 * Square magnitude of the vector.
 * \test PASSED
 */
template <typename T, typename Allocator>
T norm_sqr(const vect_n<T,Allocator>& v) {
  T sum(0);
  for(unsigned int i=0;i<v.size();++i)
    sum += v[i] * v[i];
  return sum;
};

/**
 * Unit vector in the same direction.
 * \test PASSED
 */
template <typename T, typename Allocator>
vect_n<T,Allocator> unit(const vect_n<T,Allocator>& v) {
  T m = norm(v);
  return v / m;
};

/**
 * Checks if two vectors are colinear.
 * \test PASSED
 */
template <typename T, typename Allocator>
bool colinear(const vect_n<T,Allocator>& v1, const vect_n<T,Allocator>& v2) {
  T tmp_mag2 = norm(v2);
  T tmp_mag1 = norm(v1);
  T tmp_comb = norm(v1 + v2);
  return (((tmp_mag1 + tmp_mag2) * (T(1.0) - std::numeric_limits<T>::epsilon()) < tmp_comb) || 
          (std::fabs(tmp_mag1 - tmp_mag2) * (T(1.0) + std::numeric_limits<T>::epsilon()) > tmp_comb));
};





/*******************************************************************************
                         Basic Operators
*******************************************************************************/

/**
 * Add two vectors.
 * \test PASSED
 */
template <typename T, typename Allocator>
vect_n<T,Allocator> operator +(const vect_n<T,Allocator>& v1, const vect_n<T,Allocator>& v2) {
  if(v1.size() != v2.size())
    throw std::range_error("Vector size mismatch.");
  std::vector<T,Allocator> result(v1.size());
  for(unsigned int i=0;i<v1.size();++i)
    result[i] = v1[i] + v2[i];
  return vect_n<T,Allocator>(result);
};

    /**
     * Invert the vector.
     * \test PASSED
     */
template <typename T, typename Allocator>
vect_n<T,Allocator> operator -(const vect_n<T,Allocator>& v) {
  std::vector<T,Allocator> result(v.size());
  for(unsigned int i=0;i<v.size();++i)
    result[i] = -v[i];
  return vect_n<T,Allocator>(result);
};

    /**
     * Sub two vectors.
     * \test PASSED
     */
template <typename T, typename Allocator>
vect_n<T,Allocator> operator -(const vect_n<T,Allocator>& v1, const vect_n<T,Allocator>& v2) {
  if(v1.size() != v2.size())
    throw std::range_error("Vector size mismatch.");
  std::vector<T,Allocator> result(v1.size());
  for(unsigned int i=0;i<v1.size();++i)
    result[i] = v1[i] - v2[i];
  return vect_n<T,Allocator>(result);
};

    /**
     * Dot Product.
     * \test PASSED
     */
template <typename T, typename Allocator>
T operator *(const vect_n<T,Allocator>& v1, const vect_n<T,Allocator>& v2) {
  if(v1.size() != v2.size())
    throw std::range_error("Vector size mismatch.");
  T result(0);
  for(unsigned int i=0;i<v1.size();++i)
    result += v1[i] * v2[i];
  return result;
};

    /**
     * Scalar-vector product.
     * \test PASSED
     */
template <typename T, typename Allocator>
vect_n<T,Allocator> operator *(const vect_n<T,Allocator>& v, const T& S) {
  std::vector<T,Allocator> result(v.size());
  for(unsigned int i=0;i<v.size();++i)
    result[i] = v[i] * S;
  return vect_n<T,Allocator>(result);
};

    /**
     * Scalar-vector product.
     * \test PASSED
     */
template <typename T, typename Allocator>
vect_n<T,Allocator> operator *(const T& S, const vect_n<T,Allocator>& v) {
  std::vector<T,Allocator> result(v.size());
  for(unsigned int i=0;i<v.size();++i)
    result[i] = v[i] * S;
  return vect_n<T,Allocator>(result);
};

    /**
     * Scalar-vector division.
     * \test PASSED
     */
template <typename T, typename Allocator>
vect_n<T,Allocator> operator /(const vect_n<T,Allocator>& v, const T& S) {
  std::vector<T,Allocator> result(v.size());
  for(unsigned int i=0;i<v.size();++i)
    result[i] = v[i] / S;
  return vect_n<T,Allocator>(result);
};


/**
 * Sub two vectors. For functional interfaces.
 * \test PASSED
 */
template <typename T, typename Allocator>
vect_n<T,Allocator> diff(const vect_n<T,Allocator>& v1, const vect_n<T,Allocator>& v2) {
  return v1 - v2;
};

/**
 * Add two vectors. For functional interfaces.
 * \test PASSED
 */
template <typename T, typename Allocator>
vect_n<T,Allocator> add(const vect_n<T,Allocator>& v1, const vect_n<T,Allocator>& v2) {
  return v1 + v2;
};

/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

    /**
     * Equality Comparison operator, component-wise.
     * \test PASSED
     */
template <typename T, typename Allocator>
bool operator ==(const vect_n<T,Allocator>& v1, const vect_n<T,Allocator>& v2) {
  if(v1.size() != v2.size())
    return false;
  for(unsigned int i=0;i<v1.size();++i)
    if(v1[i] != v2[i])
      return false;
  return true;
};

    /**
     * Inequality Comparison operator, component-wise.
     * \test PASSED
     */
template <typename T, typename Allocator>
bool operator !=(const vect_n<T,Allocator>& v1, const vect_n<T,Allocator>& v2) {
  if(v1.size() != v2.size())
    return true;
  for(unsigned int i=0;i<v1.size();++i)
    if(v1[i] != v2[i])
      return true;
  return false;
};

    /**
     * Greater-than Comparison operator, Euclidean norm.
     * \test PASSED
     */
template <typename T, typename Allocator>
bool operator >(const vect_n<T,Allocator>& v1, const vect_n<T,Allocator>& v2) {
  return (norm_sqr(v1) > norm_sqr(v2));
};

    /**
     * Smaller-than Comparison operator, Euclidean norm.
     * \test PASSED
     */
template <typename T, typename Allocator>
bool operator <(const vect_n<T,Allocator>& v1, const vect_n<T,Allocator>& v2) {
  return (norm_sqr(v1) < norm_sqr(v2));
};

    /**
     * Greater-or-equal Comparison operator, Euclidean norm.
     * \test PASSED
     */
template <typename T, typename Allocator>
bool operator >=(const vect_n<T,Allocator>& v1, const vect_n<T,Allocator>& v2) {
  return (norm_sqr(v1) >= norm_sqr(v2));
};

    /**
     * Smaller-or-equal Comparison operator, Euclidean norm.
     * \test PASSED
     */
template <typename T, typename Allocator>
bool operator <=(const vect_n<T,Allocator>& v1, const vect_n<T,Allocator>& v2) {
  return (norm_sqr(v1) <= norm_sqr(v2));
};




/**
 * Prints a variable-size vector to an output stream as "(v1; v2; v3; ..; vN)".
 * \test PASSED
 */
template <typename T, typename Allocator>
std::ostream& operator <<(std::ostream& out_stream,const vect_n<T,Allocator>& V) {
  out_stream << "(";
  if(V.q.size() > 0)
    out_stream << V[0];
  for(unsigned int i=1;i<V.size();++i)
    out_stream << "; " << V[i];
  return out_stream << ")";
};



};




#endif

















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
#include <algorithm>

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/rtti/so_register_type.hpp>
#include <ReaK/core/rtti/typed_primitives.hpp>
#include <ReaK/core/serialization/archiver.hpp>

#include "vect_concepts.hpp"
#include "vect_views.hpp"
#include "vect_index_iterator.hpp"
#include "mat_concepts.hpp"

#include <boost/static_assert.hpp>
#include <boost/mpl/and.hpp>
#include <boost/utility/enable_if.hpp>

namespace ReaK {


/**
 * This class implements a fixed-index vector component of primitive type.
 * This class is mainly meant to be used with the "vect" class template.
 */
template < typename T, unsigned int Index >
class vect_component {
public:
  typedef vect_component< T, Index > self;

  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;
  typedef std::size_t size_type;

  value_type q;

  explicit vect_component( const_reference Q = 0 ) BOOST_NOEXCEPT : q( Q ){};

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  value_type operator[]( size_type i ) const BOOST_NOEXCEPT { return ( i == Index ? q : value_type() ); };

  /**
   * Call operator, accessor for read only.
   * \test PASSED
   */
  value_type operator()( size_type i ) const BOOST_NOEXCEPT { return ( i == Index ? q : value_type() ); };


  friend self operator+( const self& lhs, const self& rhs ) BOOST_NOEXCEPT { return self( lhs.q + rhs.q ); };

  friend self operator-( const self& lhs, const self& rhs ) BOOST_NOEXCEPT { return self( lhs.q - rhs.q ); };

  friend self operator-( const self& lhs ) BOOST_NOEXCEPT { return self( -lhs.q ); };

  self& operator+=( const self& rhs ) BOOST_NOEXCEPT {
    q += rhs.q;
    return *this;
  };

  self& operator-=( const self& rhs ) BOOST_NOEXCEPT {
    q -= rhs.q;
    return *this;
  };

  friend self operator*(const self& lhs, const_reference rhs)BOOST_NOEXCEPT { return self( lhs.q * rhs ); };

  friend self operator*(const_reference lhs, const self& rhs)BOOST_NOEXCEPT { return self( lhs * rhs.q ); };

  self& operator*=( const_reference rhs ) BOOST_NOEXCEPT {
    q *= rhs;
    return *this;
  };
};

static const vect_component< double, 0 > vect_i = vect_component< double, 0 >( 1.0 );
static const vect_component< double, 1 > vect_j = vect_component< double, 1 >( 1.0 );
static const vect_component< double, 2 > vect_k = vect_component< double, 2 >( 1.0 );


template < typename T, typename Allocator >
class vect_n; // forward-declaration.


/**
 * This class implements a fixed-size templated vector class which holds components of primitive type.
 */
template < typename T, unsigned int Size >
class vect {
public:
  typedef vect< T, Size > self;

  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef std::allocator< T > allocator_type; // just in case it is cast to variable-length vector.

  typedef pointer iterator;
  typedef const_pointer const_iterator;

  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, dimensions = Size );

  T q[Size]; /**< Components of the vector.*/

  /**
   * Returns the size of the vector.
   */
  size_type size() const BOOST_NOEXCEPT { return Size; };
  /**
   * Returns the max-size of the vector.
   */
  size_type max_size() const BOOST_NOEXCEPT { return Size; };
  /**
   * Returns the capacity of the vector.
   */
  size_type capacity() const BOOST_NOEXCEPT { return Size; };
  /**
   * Resizes the vector.
   */
  void resize( size_type sz, T c = T() ) const BOOST_NOEXCEPT{};
  /**
   * Checks if the vector is empty.
   */
  bool empty() const BOOST_NOEXCEPT { return false; };
  /**
   * Reserve a capacity for the vector.
   */
  void reserve( size_type sz ) const BOOST_NOEXCEPT{};

  /**
   * Returns an iterator to the first element of the vector.
   */
  iterator begin() BOOST_NOEXCEPT { return q; };
  /**
   * Returns a const-iterator to the first element of the vector.
   */
  const_iterator begin() const BOOST_NOEXCEPT { return q; };
  /**
   * Returns an iterator to the one-past-last element of the vector.
   */
  iterator end() BOOST_NOEXCEPT { return q + Size; };
  /**
   * Returns a const-iterator to the one-past-last element of the vector.
   */
  const_iterator end() const BOOST_NOEXCEPT { return q + Size; };

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/
  /**
   * Default constructor: sets all to zero.
   * \test PASSED
   */
  vect() BOOST_NOEXCEPT {
    for( pointer p_i = q, p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  template < typename U, typename Allocator >
  vect( const vect_n< U, Allocator >& V ) BOOST_NOEXCEPT;

  /**
   * Constructor from an array of values of type T.
   * \test PASSED
   */
  vect( const_pointer Q ) BOOST_NOEXCEPT {
    for( pointer p_i = q, p_end = q + Size; p_i < p_end; ++p_i, ++Q )
      *p_i = *Q;
    return;
  };

#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
  /**
   * Standard Copy Constructor with standard semantics.
   * \test PASSED
   */
  vect( const self& V ) BOOST_NOEXCEPT {
    const_pointer Q = V.q;
    for( pointer p_i = q, p_end = q + Size; p_i < p_end; ++p_i, ++Q )
      *p_i = *Q;
    return;
  };
#else
  vect( const self& V ) BOOST_NOEXCEPT = default;
#endif

  template < typename U, unsigned int OtherSize >
  explicit vect( const vect< U, OtherSize >& rhs ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= OtherSize );
    const_pointer Q = rhs.q;
    pointer p_i = q;
    for( pointer p_end = q + OtherSize; p_i < p_end; ++p_i, ++Q )
      *p_i = *Q;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  template < typename U, unsigned int Index >
  explicit vect( const vect_component< U, Index >& rhs ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size > Index );
    for( pointer p_i = q, p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    q[Index] = rhs.q;
    return;
  };


#ifdef BOOST_NO_CXX11_VARIADIC_TEMPLATES

  /**
   * Constructor for 1 values.
   * \test PASSED
   */
  explicit vect( const_reference Q1 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 1 );
    pointer p_i = q;
    *p_i++ = Q1;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 2 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 2 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 3 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 3 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 4 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 4 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 5 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4,
        const_reference Q5 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 5 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 6 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 6 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 7 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6, const_reference Q7 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 7 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    *p_i++ = Q7;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 8 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6, const_reference Q7, const_reference Q8 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 8 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    *p_i++ = Q7;
    *p_i++ = Q8;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 9 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 9 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    *p_i++ = Q7;
    *p_i++ = Q8;
    *p_i++ = Q9;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 10 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9,
        const_reference Q10 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 10 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    *p_i++ = Q7;
    *p_i++ = Q8;
    *p_i++ = Q9;
    *p_i++ = Q10;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 11 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
        const_reference Q11 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 11 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    *p_i++ = Q7;
    *p_i++ = Q8;
    *p_i++ = Q9;
    *p_i++ = Q10;
    *p_i++ = Q11;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 12 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
        const_reference Q11, const_reference Q12 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 12 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    *p_i++ = Q7;
    *p_i++ = Q8;
    *p_i++ = Q9;
    *p_i++ = Q10;
    *p_i++ = Q11;
    *p_i++ = Q12;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 13 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
        const_reference Q11, const_reference Q12, const_reference Q13 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 13 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    *p_i++ = Q7;
    *p_i++ = Q8;
    *p_i++ = Q9;
    *p_i++ = Q10;
    *p_i++ = Q11;
    *p_i++ = Q12;
    *p_i++ = Q13;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 14 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
        const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 14 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    *p_i++ = Q7;
    *p_i++ = Q8;
    *p_i++ = Q9;
    *p_i++ = Q10;
    *p_i++ = Q11;
    *p_i++ = Q12;
    *p_i++ = Q13;
    *p_i++ = Q14;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 15 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
        const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14,
        const_reference Q15 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 15 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    *p_i++ = Q7;
    *p_i++ = Q8;
    *p_i++ = Q9;
    *p_i++ = Q10;
    *p_i++ = Q11;
    *p_i++ = Q12;
    *p_i++ = Q13;
    *p_i++ = Q14;
    *p_i++ = Q15;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 16 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
        const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
        const_reference Q16 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 16 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    *p_i++ = Q7;
    *p_i++ = Q8;
    *p_i++ = Q9;
    *p_i++ = Q10;
    *p_i++ = Q11;
    *p_i++ = Q12;
    *p_i++ = Q13;
    *p_i++ = Q14;
    *p_i++ = Q15;
    *p_i++ = Q16;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 17 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
        const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
        const_reference Q16, const_reference Q17 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 17 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    *p_i++ = Q7;
    *p_i++ = Q8;
    *p_i++ = Q9;
    *p_i++ = Q10;
    *p_i++ = Q11;
    *p_i++ = Q12;
    *p_i++ = Q13;
    *p_i++ = Q14;
    *p_i++ = Q15;
    *p_i++ = Q16;
    *p_i++ = Q17;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 18 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
        const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
        const_reference Q16, const_reference Q17, const_reference Q18 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 18 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    *p_i++ = Q7;
    *p_i++ = Q8;
    *p_i++ = Q9;
    *p_i++ = Q10;
    *p_i++ = Q11;
    *p_i++ = Q12;
    *p_i++ = Q13;
    *p_i++ = Q14;
    *p_i++ = Q15;
    *p_i++ = Q16;
    *p_i++ = Q17;
    *p_i++ = Q18;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 19 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
        const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
        const_reference Q16, const_reference Q17, const_reference Q18, const_reference Q19 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 19 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    *p_i++ = Q7;
    *p_i++ = Q8;
    *p_i++ = Q9;
    *p_i++ = Q10;
    *p_i++ = Q11;
    *p_i++ = Q12;
    *p_i++ = Q13;
    *p_i++ = Q14;
    *p_i++ = Q15;
    *p_i++ = Q16;
    *p_i++ = Q17;
    *p_i++ = Q18;
    *p_i++ = Q19;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

  /**
   * Constructor for 20 values.
   * \test PASSED
   */
  vect( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
        const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
        const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
        const_reference Q16, const_reference Q17, const_reference Q18, const_reference Q19,
        const_reference Q20 ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size >= 20 );
    pointer p_i = q;
    *p_i++ = Q1;
    *p_i++ = Q2;
    *p_i++ = Q3;
    *p_i++ = Q4;
    *p_i++ = Q5;
    *p_i++ = Q6;
    *p_i++ = Q7;
    *p_i++ = Q8;
    *p_i++ = Q9;
    *p_i++ = Q10;
    *p_i++ = Q11;
    *p_i++ = Q12;
    *p_i++ = Q13;
    *p_i++ = Q14;
    *p_i++ = Q15;
    *p_i++ = Q16;
    *p_i++ = Q17;
    *p_i++ = Q18;
    *p_i++ = Q19;
    *p_i++ = Q20;
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

#else

private:
  template < typename Arg1 >
  static void set_value_impl( pointer& pval, const Arg1& a1 ) BOOST_NOEXCEPT {
    *pval++ = value_type( a1 );
  };

  template < typename Arg1, typename... Args >
  static void set_value_impl( pointer& pval, const Arg1& a1, const Args&... tail ) BOOST_NOEXCEPT {
    *pval++ = value_type( a1 );
    set_value_impl( pval, tail... );
  };

public:
  /**
   * Constructor for Size values.
   * \test PASSED
   */
  template < typename... Args >
  vect( const value_type& a1, Args... args ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size > sizeof...( Args ) );
    pointer p_i = q;
    set_value_impl( p_i, a1, args... );
    for( pointer p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    return;
  };

#endif

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Array indexing operator, accessor for read/write.
   * \test PASSED
   */
  reference operator[](size_type i)BOOST_NOEXCEPT { return q[i]; };

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator[]( size_type i ) const BOOST_NOEXCEPT { return q[i]; };

  /**
   * Sub-vector operator, accessor for read/write.
   * \test PASSED
   */
  vect_ref_view< self > operator[]( const std::pair< size_type, size_type >& r ) { return sub( *this )[r]; };

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_const_ref_view< self > operator[]( const std::pair< size_type, size_type >& r ) const {
    return sub( *this )[r];
  };

  /**
   * Array indexing operator, accessor for read/write.
   * \test PASSED
   */
  reference operator()( size_type i ) BOOST_NOEXCEPT { return q[i]; };

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator()( size_type i ) const BOOST_NOEXCEPT { return q[i]; };

/*******************************************************************************
                         Assignment Operators
*******************************************************************************/

#ifdef BOOST_NO_DEFAULTED_FUNCTIONS
  /**
   * Standard assignment operator.
   * \test PASSED
   */
  self& operator=( const self& V ) BOOST_NOEXCEPT {
    const_pointer Q = V.q;
    for( pointer p_i = q, p_end = q + Size; p_i < p_end; ++p_i, ++Q )
      *p_i = *Q;
    return *this;
  };
#else
  self& operator=( const self& V ) BOOST_NOEXCEPT = default;
#endif

  /**
   * Standard assignment operator.
   * \test PASSED
   */
  template < typename Vector >
  typename boost::enable_if< boost::mpl::and_< is_readable_vector< Vector >,
                                               boost::mpl::not_< boost::is_same< Vector, self > > >,
                             self& >::type
    operator=( const Vector& V ) {
    if( Size != V.size() )
      throw std::range_error( "Vector size mismatch." );
    for( size_type i = 0; i < Size; ++i )
      q[i] = V[i];
    return *this;
  };

  template < typename U, unsigned int Index >
  self& operator=( const vect_component< U, Index >& V ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size > Index );
    for( pointer p_i = q, p_end = q + Size; p_i < p_end; ++p_i )
      *p_i = value_type();
    q[Index] = V.q;
  };

  template < typename U, unsigned int Index >
  friend self operator+( self lhs, const vect_component< U, Index >& rhs ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size > Index );
    lhs.q[Index] += rhs.q;
    return lhs;
  };

  template < typename U, unsigned int Index >
  self& operator+=( const vect_component< U, Index >& rhs ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size > Index );
    q[Index] += rhs.q;
    return *this;
  };

  template < typename U, unsigned int Index >
  friend self operator+( const vect_component< U, Index >& lhs, self rhs ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size > Index );
    rhs.q[Index] += lhs.q;
    return rhs;
  };

  template < typename U, unsigned int Index >
  friend self operator-( self lhs, const vect_component< U, Index >& rhs ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size > Index );
    lhs.q[Index] -= rhs.q;
    return lhs;
  };

  template < typename U, unsigned int Index >
  self& operator-=( const vect_component< U, Index >& rhs ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size > Index );
    q[Index] -= rhs.q;
    return *this;
  };

  template < typename U, unsigned int Index >
  friend self operator-( const vect_component< U, Index >& lhs, const self& rhs ) BOOST_NOEXCEPT {
    BOOST_STATIC_ASSERT( Size > Index );
    self result = -rhs;
    result.q[Index] += lhs.q;
    return result;
  };
};

namespace rtti {

template < typename T, unsigned int Size >
struct get_type_id< vect< T, Size > > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000011 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "vect" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "vect"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef const vect< T, Size >& save_type;
  typedef vect< T, Size >& load_type;
};

template < typename T, unsigned int Size, typename Tail >
struct get_type_info< vect< T, Size >, Tail > {
  typedef type_id< vect< T, Size >,
                   typename get_type_info< T, get_type_info< boost::mpl::integral_c< unsigned int, Size >,
                                                             Tail > >::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< vect< T, Size > >::type_name + lsl_left_bracket
                                          + get_type_id< T >::type_name + lsl_comma
                                          + get_type_id< boost::mpl::integral_c< unsigned int, Size > >::type_name
                                          + lsl_right_bracket + get_type_name_tail< Tail >::value;
#else
  static std::string type_name() {
    std::string result = get_type_id< vect< T, Size > >::type_name();
    result += "<";
    result += get_type_id< T >::type_name();
    result += ",";
    result += get_type_id< boost::mpl::integral_c< unsigned int, Size > >::type_name();
    result += ">";
    result += get_type_name_tail< Tail >::value();
    return result; // NRVO
  };
#endif
};
};

namespace serialization {
template < typename T, unsigned int Size >
iarchive& operator>>( iarchive& in, vect< T, Size >& v ) {
  for( unsigned int i = 0; i != Size; ++i )
    in >> v[i];
  return in;
};

template < typename T, unsigned int Size >
iarchive& operator&( iarchive& in, const std::pair< std::string, vect< T, Size >& >& v ) {
  for( unsigned int i = 0; i != Size; ++i ) {
    std::stringstream s_stream;
    s_stream << v.first << "_q[" << i << "]";
    in& RK_SERIAL_LOAD_WITH_ALIAS( s_stream.str(), v.second[i] );
  };
  return in;
};

template < typename T, unsigned int Size >
oarchive& operator<<( oarchive& out, const vect< T, Size >& v ) {
  for( unsigned int i = 0; i != Size; ++i )
    out << v[i];
  return out;
};

template < typename T, unsigned int Size >
oarchive& operator&( oarchive& out, const std::pair< std::string, const vect< T, Size >& >& v ) {
  for( unsigned int i = 0; i != Size; ++i ) {
    std::stringstream s_stream;
    s_stream << v.first << "_q[" << i << "]";
    out& RK_SERIAL_SAVE_WITH_ALIAS( s_stream.str(), v.second[i] );
  };
  return out;
};
};


template < typename T, unsigned int Size >
struct is_readable_vector< vect< T, Size > > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_vector< vect< T, Size > > type;
};

template < typename T, unsigned int Size >
struct is_writable_vector< vect< T, Size > > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_writable_vector< vect< T, Size > > type;
};

template < typename T, unsigned int Size >
struct is_resizable_vector< vect< T, Size > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_vector< vect< T, Size > > type;
};


template < typename T, unsigned int Size >
struct has_allocator_vector< vect< T, Size > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef has_allocator_vector< vect< T, Size > > type;
};


/*******************************************************************************
                         Basic Constructors
*******************************************************************************/

template < typename T >
vect< T, 1 > make_vect( const T& Q1 ) BOOST_NOEXCEPT {
  return vect< T, 1 >( Q1 );
};

template < typename T >
vect< T, 2 > make_vect( const T& Q1, const T& Q2 ) BOOST_NOEXCEPT {
  return vect< T, 2 >( Q1, Q2 );
};

template < typename T >
vect< T, 3 > make_vect( const T& Q1, const T& Q2, const T& Q3 ) BOOST_NOEXCEPT {
  return vect< T, 3 >( Q1, Q2, Q3 );
};

template < typename T >
vect< T, 4 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4 ) BOOST_NOEXCEPT {
  return vect< T, 4 >( Q1, Q2, Q3, Q4 );
};

template < typename T >
vect< T, 5 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5 ) BOOST_NOEXCEPT {
  return vect< T, 5 >( Q1, Q2, Q3, Q4, Q5 );
};

template < typename T >
vect< T, 6 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6 ) BOOST_NOEXCEPT {
  return vect< T, 6 >( Q1, Q2, Q3, Q4, Q5, Q6 );
};

template < typename T >
vect< T, 7 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6,
                        const T& Q7 ) BOOST_NOEXCEPT {
  return vect< T, 7 >( Q1, Q2, Q3, Q4, Q5, Q6, Q7 );
};

template < typename T >
vect< T, 8 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                        const T& Q8 ) BOOST_NOEXCEPT {
  return vect< T, 8 >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8 );
};

template < typename T >
vect< T, 9 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                        const T& Q8, const T& Q9 ) BOOST_NOEXCEPT {
  return vect< T, 9 >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9 );
};

template < typename T >
vect< T, 10 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                         const T& Q8, const T& Q9, const T& Q10 ) BOOST_NOEXCEPT {
  return vect< T, 10 >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10 );
};

template < typename T >
vect< T, 11 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                         const T& Q8, const T& Q9, const T& Q10, const T& Q11 ) BOOST_NOEXCEPT {
  return vect< T, 11 >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11 );
};

template < typename T >
vect< T, 12 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                         const T& Q8, const T& Q9, const T& Q10, const T& Q11, const T& Q12 ) BOOST_NOEXCEPT {
  return vect< T, 12 >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12 );
};

template < typename T >
vect< T, 13 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                         const T& Q8, const T& Q9, const T& Q10, const T& Q11, const T& Q12,
                         const T& Q13 ) BOOST_NOEXCEPT {
  return vect< T, 13 >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13 );
};

template < typename T >
vect< T, 14 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                         const T& Q8, const T& Q9, const T& Q10, const T& Q11, const T& Q12, const T& Q13,
                         const T& Q14 ) BOOST_NOEXCEPT {
  return vect< T, 14 >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13, Q14 );
};

template < typename T >
vect< T, 15 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                         const T& Q8, const T& Q9, const T& Q10, const T& Q11, const T& Q12, const T& Q13, const T& Q14,
                         const T& Q15 ) BOOST_NOEXCEPT {
  return vect< T, 15 >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13, Q14, Q15 );
};

template < typename T >
vect< T, 16 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                         const T& Q8, const T& Q9, const T& Q10, const T& Q11, const T& Q12, const T& Q13, const T& Q14,
                         const T& Q15, const T& Q16 ) BOOST_NOEXCEPT {
  return vect< T, 16 >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13, Q14, Q15, Q16 );
};

template < typename T >
vect< T, 17 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                         const T& Q8, const T& Q9, const T& Q10, const T& Q11, const T& Q12, const T& Q13, const T& Q14,
                         const T& Q15, const T& Q16, const T& Q17 ) BOOST_NOEXCEPT {
  return vect< T, 17 >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13, Q14, Q15, Q16, Q17 );
};

template < typename T >
vect< T, 18 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                         const T& Q8, const T& Q9, const T& Q10, const T& Q11, const T& Q12, const T& Q13, const T& Q14,
                         const T& Q15, const T& Q16, const T& Q17, const T& Q18 ) BOOST_NOEXCEPT {
  return vect< T, 18 >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13, Q14, Q15, Q16, Q17, Q18 );
};

template < typename T >
vect< T, 19 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                         const T& Q8, const T& Q9, const T& Q10, const T& Q11, const T& Q12, const T& Q13, const T& Q14,
                         const T& Q15, const T& Q16, const T& Q17, const T& Q18, const T& Q19 ) BOOST_NOEXCEPT {
  return vect< T, 19 >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13, Q14, Q15, Q16, Q17, Q18, Q19 );
};

template < typename T >
vect< T, 20 > make_vect( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                         const T& Q8, const T& Q9, const T& Q10, const T& Q11, const T& Q12, const T& Q13, const T& Q14,
                         const T& Q15, const T& Q16, const T& Q17, const T& Q18, const T& Q19,
                         const T& Q20 ) BOOST_NOEXCEPT {
  return vect< T, 20 >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13, Q14, Q15, Q16, Q17, Q18, Q19, Q20 );
};


/*******************************************************************************
                         Basic Operators
*******************************************************************************/


/**
 * Sub two vectors. For functional interfaces.
 * \test PASSED
 */
template < typename T, unsigned int Size >
vect< T, Size > diff( const vect< T, Size >& v1, const vect< T, Size >& v2 ) BOOST_NOEXCEPT {
  return v1 - v2;
};

/**
 * Add two vectors. For functional interfaces.
 * \test PASSED
 */
template < typename T, unsigned int Size >
vect< T, Size > add( const vect< T, Size >& v1, const vect< T, Size >& v2 ) BOOST_NOEXCEPT {
  return v1 + v2;
};


/*******************************************************************************
                         Special Vector Products / Operators
*******************************************************************************/


/**
 * 2D Cross-Product.
 * \test PASSED
 */
template < class T >
T operator%( const vect< T, 2 >& V1, const vect< T, 2 >& V2 ) BOOST_NOEXCEPT {
  return V1[0] * V2[1] - V1[1] * V2[0];
};

template < class T >
T operator%( const vect_component< T, 0 >& V1, const vect_component< T, 0 >& V2 ) BOOST_NOEXCEPT {
  return T( 0.0 );
};

template < class T >
T operator%( const vect_component< T, 1 >& V1, const vect_component< T, 1 >& V2 ) BOOST_NOEXCEPT {
  return T( 0.0 );
};

template < class T >
vect_component< T, 2 > operator%( const vect_component< T, 0 >& V1, const vect_component< T, 1 >& V2 ) BOOST_NOEXCEPT {
  return vect_component< T, 2 >( V1.q * V2.q );
};

template < class T >
vect_component< T, 2 > operator%( const vect_component< T, 1 >& V1, const vect_component< T, 0 >& V2 ) BOOST_NOEXCEPT {
  return vect_component< T, 2 >( -V1.q * V2.q );
};

/**
 * 2D Cross-Product.
 * \test PASSED
 */
template < class T >
vect< T, 2 > operator%( const T& S, const vect< T, 2 >& V ) BOOST_NOEXCEPT {
  vect< T, 2 > result;
  result[0] = -V[1] * S;
  result[1] = V[0] * S;
  return result;
};

template < class T >
vect_component< T, 1 > operator%( const T& S, const vect_component< T, 0 >& V ) BOOST_NOEXCEPT {
  return vect_component< T, 1 >( V.q * S );
};

template < class T >
vect_component< T, 0 > operator%( const T& S, const vect_component< T, 1 >& V ) BOOST_NOEXCEPT {
  return vect_component< T, 0 >( -V.q * S );
};

/**
 * 2D Cross-Product.
 * \test PASSED
 */
template < class T >
vect< T, 2 > operator%( const vect< T, 2 >& V, const T& S ) BOOST_NOEXCEPT {
  vect< T, 2 > result;
  result[0] = V[1] * S;
  result[1] = -V[0] * S;
  return result;
};

template < class T >
vect_component< T, 1 > operator%( const vect_component< T, 0 >& V, const T& S ) BOOST_NOEXCEPT {
  return vect_component< T, 1 >( -V.q * S );
};

template < class T >
vect_component< T, 0 > operator%( const vect_component< T, 1 >& V, const T& S ) BOOST_NOEXCEPT {
  return vect_component< T, 0 >( V.q * S );
};

/**
 * 3D Cross-Product.
 * \test PASSED
 */
template < class T >
vect< T, 3 > operator%( const vect< T, 3 >& V1, const vect< T, 3 >& V2 ) BOOST_NOEXCEPT {
  vect< T, 3 > result;
  result[0] = V1[1] * V2[2] - V1[2] * V2[1];
  result[1] = V1[2] * V2[0] - V1[0] * V2[2];
  result[2] = V1[0] * V2[1] - V1[1] * V2[0];
  return result;
};

template < class T >
vect_component< T, 1 > operator%( const vect_component< T, 0 >& V1, const vect_component< T, 2 >& V2 ) BOOST_NOEXCEPT {
  return vect_component< T, 1 >( -V1.q * V2.q );
};

template < class T >
vect_component< T, 1 > operator%( const vect_component< T, 2 >& V1, const vect_component< T, 0 >& V2 ) BOOST_NOEXCEPT {
  return vect_component< T, 1 >( V1.q * V2.q );
};

template < class T >
vect_component< T, 0 > operator%( const vect_component< T, 1 >& V1, const vect_component< T, 2 >& V2 ) BOOST_NOEXCEPT {
  return vect_component< T, 2 >( V1.q * V2.q );
};

template < class T >
vect_component< T, 0 > operator%( const vect_component< T, 2 >& V1, const vect_component< T, 1 >& V2 ) BOOST_NOEXCEPT {
  return vect_component< T, 2 >( -V1.q * V2.q );
};

template < class T >
T operator%( const vect_component< T, 2 >& V1, const vect_component< T, 2 >& V2 ) BOOST_NOEXCEPT {
  return T( 0.0 );
};

template < class T >
vect< T, 3 > operator%( const vect< T, 3 >& V1, const vect_component< T, 0 >& V2 ) BOOST_NOEXCEPT {
  vect< T, 3 > result;
  result[0] = T( 0.0 );
  result[1] = V1[2] * V2.q;
  result[2] = -V1[1] * V2.q;
  return result;
};

template < class T >
vect< T, 3 > operator%( const vect< T, 3 >& V1, const vect_component< T, 1 >& V2 ) BOOST_NOEXCEPT {
  vect< T, 3 > result;
  result[0] = -V1[2] * V2.q;
  result[1] = T( 0.0 );
  result[2] = V1[0] * V2.q;
  return result;
};

template < class T >
vect< T, 3 > operator%( const vect< T, 3 >& V1, const vect_component< T, 2 >& V2 ) BOOST_NOEXCEPT {
  vect< T, 3 > result;
  result[0] = V1[1] * V2.q;
  result[1] = -V1[0] * V2.q;
  result[2] = T( 0.0 );
  return result;
};

template < class T >
vect< T, 3 > operator%( const vect_component< T, 0 >& V1, const vect< T, 3 >& V2 ) BOOST_NOEXCEPT {
  vect< T, 3 > result;
  result[0] = T( 0.0 );
  result[1] = -V2[2] * V1.q;
  result[2] = V2[1] * V1.q;
  return result;
};

template < class T >
vect< T, 3 > operator%( const vect_component< T, 1 >& V1, const vect< T, 3 >& V2 ) BOOST_NOEXCEPT {
  vect< T, 3 > result;
  result[0] = V2[2] * V1.q;
  result[1] = T( 0.0 );
  result[2] = -V2[0] * V1.q;
  return result;
};

template < class T >
vect< T, 3 > operator%( const vect_component< T, 2 >& V1, const vect< T, 3 >& V2 ) BOOST_NOEXCEPT {
  vect< T, 3 > result;
  result[0] = -V2[1] * V1.q;
  result[1] = V2[0] * V1.q;
  result[2] = T( 0.0 );
  return result;
};


/**
 * This class implements a variable-size templated vector class which holds components of dimensional quantities.
 */
template < typename T, typename Allocator = std::allocator< T > >
class vect_n {
public:
  typedef vect_n< T, Allocator > self;

  typedef T value_type;
  typedef typename std::vector< T, Allocator >::reference reference;
  typedef typename std::vector< T, Allocator >::const_reference const_reference;
  typedef typename std::vector< T, Allocator >::pointer pointer;
  typedef typename std::vector< T, Allocator >::const_pointer const_pointer;
  typedef Allocator allocator_type;

  typedef typename std::vector< T, Allocator >::iterator iterator;
  typedef typename std::vector< T, Allocator >::const_iterator const_iterator;

  typedef typename std::vector< T, Allocator >::size_type size_type;
  typedef typename std::vector< T, Allocator >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, dimensions = 0 );

  std::vector< T, Allocator > q; /**< Components of the vector. */

  /**
   * Returns the size of the vector.
   */
  size_type size() const BOOST_NOEXCEPT { return q.size(); };
  /**
   * Returns the max-size of the vector.
   */
  size_type max_size() const BOOST_NOEXCEPT { return q.max_size(); };
  /**
   * Returns the capacity of the vector.
   */
  size_type capacity() const BOOST_NOEXCEPT { return q.capacity(); };
  /**
   * Resizes the vector.
   */
  void resize( size_type sz, T c = T() ) { q.resize( sz, c ); };
  /**
   * Checks if the vector is empty.
   */
  bool empty() const BOOST_NOEXCEPT { return q.empty(); };
  /**
   * Reserve a capacity for the vector.
   */
  void reserve( size_type sz ) { q.reserve( sz ); };

  /**
   * Returns an iterator to the first element of the vector.
   */
  iterator begin() BOOST_NOEXCEPT { return q.begin(); };
  /**
   * Returns a const-iterator to the first element of the vector.
   */
  const_iterator begin() const BOOST_NOEXCEPT { return q.begin(); };
  /**
   * Returns an iterator to the one-past-last element of the vector.
   */
  iterator end() BOOST_NOEXCEPT { return q.end(); };
  /**
   * Returns a const-iterator to the one-past-last element of the vector.
   */
  const_iterator end() const BOOST_NOEXCEPT { return q.end(); };

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/
  /**
   * Default constructor.
   * \test PASSED
   */
  vect_n( const allocator_type& aAlloc = allocator_type() ) : q( aAlloc ){};

  /**
   * Constructor for a given size or dimension of vector.
   * \test PASSED
   */
  explicit vect_n( size_type aSize, const_reference Fill = value_type(),
                   const allocator_type& aAlloc = allocator_type() )
      : q( aSize, Fill, aAlloc ){};

  /**
   * Constructor from an array of values of type "value_type".
   * \test PASSED
   */
  template < typename OtherAllocator >
  explicit vect_n( const std::vector< value_type, OtherAllocator >& Q, const allocator_type& aAlloc = allocator_type() )
      : q( Q.begin(), Q.end(), aAlloc ){};

  /**
   * Constructor from a forward iterator of values of type "value_type".
   * \test PASSED
   */
  template < typename InputIter >
  explicit vect_n( InputIter first, InputIter last, const allocator_type& aAlloc = allocator_type() )
      : q( first, last, aAlloc ){};

  /**
   * Standard Copy Constructor with standard semantics.
   * \test PASSED
   */
  vect_n( const self& V ) : q( V.begin(), V.end(), V.get_allocator() ){};

  /**
   * Allocator-agnostic Copy Constructor with standard semantics.
   * \test PASSED
   */
  template < typename OtherAllocator >
  vect_n( const vect_n< value_type, OtherAllocator >& V, const allocator_type& aAlloc = allocator_type() )
      : q( V.begin(), V.end(), aAlloc ){};

  /**
   * Constructor from a fixed-length vector.
   */
  template < unsigned int Size >
  vect_n( const vect< value_type, Size >& V, const allocator_type& aAlloc = allocator_type() )
      : q( V.begin(), V.end(), aAlloc ){};

#ifdef BOOST_NO_CXX11_VARIADIC_TEMPLATES

  /**
   * Constructor for 3 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3 ) : q( 3 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    return;
  };

  /**
   * Constructor for 4 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4 ) : q( 4 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    return;
  };

  /**
   * Constructor for 5 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5 )
      : q( 5 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    return;
  };

  /**
   * Constructor for 6 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6 )
      : q( 6 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    return;
  };

  /**
   * Constructor for 7 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6, const_reference Q7 )
      : q( 7 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    *it++ = Q7;
    return;
  };

  /**
   * Constructor for 8 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6, const_reference Q7, const_reference Q8 )
      : q( 8 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    *it++ = Q7;
    *it++ = Q8;
    return;
  };

  /**
   * Constructor for 9 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9 )
      : q( 9 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    *it++ = Q7;
    *it++ = Q8;
    *it++ = Q9;
    return;
  };

  /**
   * Constructor for 10 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10 )
      : q( 10 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    *it++ = Q7;
    *it++ = Q8;
    *it++ = Q9;
    *it++ = Q10;
    return;
  };

  /**
   * Constructor for 11 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
          const_reference Q11 )
      : q( 11 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    *it++ = Q7;
    *it++ = Q8;
    *it++ = Q9;
    *it++ = Q10;
    *it++ = Q11;
    return;
  };

  /**
   * Constructor for 12 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
          const_reference Q11, const_reference Q12 )
      : q( 12 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    *it++ = Q7;
    *it++ = Q8;
    *it++ = Q9;
    *it++ = Q10;
    *it++ = Q11;
    *it++ = Q12;
    return;
  };

  /**
   * Constructor for 13 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
          const_reference Q11, const_reference Q12, const_reference Q13 )
      : q( 13 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    *it++ = Q7;
    *it++ = Q8;
    *it++ = Q9;
    *it++ = Q10;
    *it++ = Q11;
    *it++ = Q12;
    *it++ = Q13;
    return;
  };

  /**
   * Constructor for 14 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
          const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14 )
      : q( 14 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    *it++ = Q7;
    *it++ = Q8;
    *it++ = Q9;
    *it++ = Q10;
    *it++ = Q11;
    *it++ = Q12;
    *it++ = Q13;
    *it++ = Q14;
    return;
  };

  /**
   * Constructor for 15 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
          const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15 )
      : q( 15 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    *it++ = Q7;
    *it++ = Q8;
    *it++ = Q9;
    *it++ = Q10;
    *it++ = Q11;
    *it++ = Q12;
    *it++ = Q13;
    *it++ = Q14;
    *it++ = Q15;
    return;
  };

  /**
   * Constructor for 16 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
          const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
          const_reference Q16 )
      : q( 16 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    *it++ = Q7;
    *it++ = Q8;
    *it++ = Q9;
    *it++ = Q10;
    *it++ = Q11;
    *it++ = Q12;
    *it++ = Q13;
    *it++ = Q14;
    *it++ = Q15;
    *it++ = Q16;
    return;
  };

  /**
   * Constructor for 17 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
          const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
          const_reference Q16, const_reference Q17 )
      : q( 17 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    *it++ = Q7;
    *it++ = Q8;
    *it++ = Q9;
    *it++ = Q10;
    *it++ = Q11;
    *it++ = Q12;
    *it++ = Q13;
    *it++ = Q14;
    *it++ = Q15;
    *it++ = Q16;
    *it++ = Q17;
    return;
  };

  /**
   * Constructor for 18 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
          const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
          const_reference Q16, const_reference Q17, const_reference Q18 )
      : q( 18 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    *it++ = Q7;
    *it++ = Q8;
    *it++ = Q9;
    *it++ = Q10;
    *it++ = Q11;
    *it++ = Q12;
    *it++ = Q13;
    *it++ = Q14;
    *it++ = Q15;
    *it++ = Q16;
    *it++ = Q17;
    *it++ = Q18;
    return;
  };

  /**
   * Constructor for 19 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
          const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
          const_reference Q16, const_reference Q17, const_reference Q18, const_reference Q19 )
      : q( 19 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    *it++ = Q7;
    *it++ = Q8;
    *it++ = Q9;
    *it++ = Q10;
    *it++ = Q11;
    *it++ = Q12;
    *it++ = Q13;
    *it++ = Q14;
    *it++ = Q15;
    *it++ = Q16;
    *it++ = Q17;
    *it++ = Q18;
    *it++ = Q19;
    return;
  };

  /**
   * Constructor for 20 values.
   * \test PASSED
   */
  vect_n( const_reference Q1, const_reference Q2, const_reference Q3, const_reference Q4, const_reference Q5,
          const_reference Q6, const_reference Q7, const_reference Q8, const_reference Q9, const_reference Q10,
          const_reference Q11, const_reference Q12, const_reference Q13, const_reference Q14, const_reference Q15,
          const_reference Q16, const_reference Q17, const_reference Q18, const_reference Q19, const_reference Q20 )
      : q( 20 ) {
    iterator it = q.begin();
    *it++ = Q1;
    *it++ = Q2;
    *it++ = Q3;
    *it++ = Q4;
    *it++ = Q5;
    *it++ = Q6;
    *it++ = Q7;
    *it++ = Q8;
    *it++ = Q9;
    *it++ = Q10;
    *it++ = Q11;
    *it++ = Q12;
    *it++ = Q13;
    *it++ = Q14;
    *it++ = Q15;
    *it++ = Q16;
    *it++ = Q17;
    *it++ = Q18;
    *it++ = Q19;
    *it++ = Q20;
    return;
  };

#else

private:
  static void set_value_impl( iterator& pval, const_reference a1 ) BOOST_NOEXCEPT { *pval++ = a1; };

  template < typename... Args >
  static void set_value_impl( iterator& pval, const_reference a1, const Args&... tail ) BOOST_NOEXCEPT {
    *pval++ = a1;
    set_value_impl( pval, tail... );
  };

public:
  /**
   * Constructor for Size values.
   * \test PASSED
   */
  template < typename... Args >
  vect_n( const_reference a1, const_reference a2, const_reference a3, Args... args ) BOOST_NOEXCEPT
    : q( 3 + sizeof...( Args ) ) {
    iterator p_i = q.begin();
    set_value_impl( p_i, a1, a2, a3, args... );
    return;
  };

#endif

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Array indexing operator, accessor for read/write. <
   * \test PASSED
   */
  reference operator[](size_type i)BOOST_NOEXCEPT { return q[i]; };

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator[]( size_type i ) const BOOST_NOEXCEPT { return q[i]; };

  /**
   * Sub-vector operator, accessor for read/write.
   * \test PASSED
   */
  vect_ref_view< self > operator[](const std::pair< size_type, size_type >& r)BOOST_NOEXCEPT {
    return sub( *this )[r];
  };

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_const_ref_view< self > operator[]( const std::pair< size_type, size_type >& r ) const BOOST_NOEXCEPT {
    return sub( *this )[r];
  };

  /**
   * Array indexing operator, accessor for read/write. <
   * \test PASSED
   */
  reference operator()( size_type i ) BOOST_NOEXCEPT { return q[i]; };

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator()( size_type i ) const BOOST_NOEXCEPT { return q[i]; };

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
  self& operator=( const self& V ) {
    q.assign( V.q.begin(), V.q.end() );
    return *this;
  };

  /**
   * Standard assignment operator.
   * \test PASSED
   */
  template < typename Vector >
  typename boost::enable_if< boost::mpl::and_< is_readable_vector< Vector >,
                                               boost::mpl::not_< boost::is_same< Vector, self > > >,
                             self& >::type
    operator=( const Vector& V ) {
    q.resize( V.size() );
    for( size_type i = 0; i < q.size(); ++i )
      q[i] = V[i];
    return *this;
  };
};


namespace rtti {

template < typename T, typename Allocator >
struct get_type_id< vect_n< T, Allocator > > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000010 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "vect_n" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "vect_n"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef const vect_n< T, Allocator >& save_type;
  typedef vect_n< T, Allocator >& load_type;
};

template < typename T, typename Allocator, typename Tail >
struct get_type_info< vect_n< T, Allocator >, Tail > {
  typedef type_id< vect_n< T, Allocator >, typename get_type_info< T, Tail >::type > type;
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = get_type_id< vect_n< T, Allocator > >::type_name + lsl_left_bracket
                                          + get_type_id< T >::type_name + lsl_right_bracket
                                          + get_type_name_tail< Tail >::value;
#else
  static std::string type_name() {
    std::string result = get_type_id< vect_n< T, Allocator > >::type_name();
    result += "<";
    result += get_type_id< T >::type_name();
    result += ">";
    result += get_type_name_tail< Tail >::value();
    return result; // NRVO
  };
#endif
};
};

namespace serialization {
template < typename T, typename Allocator >
iarchive& operator>>( iarchive& in, vect_n< T, Allocator >& v ) {
  in >> v.q;
  return in;
};

template < typename T, typename Allocator >
iarchive& operator&( iarchive& in, const std::pair< std::string, vect_n< T, Allocator >& >& v ) {
  in& RK_SERIAL_LOAD_WITH_ALIAS( v.first, v.second.q );
  return in;
};

template < typename T, typename Allocator >
oarchive& operator<<( oarchive& out, const vect_n< T, Allocator >& v ) {
  out << v.q;
  return out;
};

template < typename T, typename Allocator >
oarchive& operator&( oarchive& out, const std::pair< std::string, const vect_n< T, Allocator >& >& v ) {
  out& RK_SERIAL_SAVE_WITH_ALIAS( v.first, v.second.q );
  return out;
};
};


template < typename T, typename Allocator >
struct is_readable_vector< vect_n< T, Allocator > > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_vector< vect_n< T, Allocator > > type;
};

template < typename T, typename Allocator >
struct is_writable_vector< vect_n< T, Allocator > > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_writable_vector< vect_n< T, Allocator > > type;
};

template < typename T, typename Allocator >
struct is_resizable_vector< vect_n< T, Allocator > > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_resizable_vector< vect_n< T, Allocator > > type;
};

template < typename T, typename Allocator >
struct has_allocator_vector< vect_n< T, Allocator > > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef has_allocator_vector< vect_n< T, Allocator > > type;
};


template < typename Vector >
struct vect_copy< vect_ref_view< Vector > > {
  typedef vect_n< typename vect_traits< Vector >::value_type > type;
};

template < typename Vector >
struct vect_copy< vect_const_ref_view< Vector > > {
  typedef vect_n< typename vect_traits< Vector >::value_type > type;
};


/*******************************************************************************
                         Basic Constructors
*******************************************************************************/

template < typename T >
vect_n< T > make_vect_n( const T& Q1, const T& Q2, const T& Q3 ) {
  return vect_n< T >( Q1, Q2, Q3 );
};

template < typename T >
vect_n< T > make_vect_n( const T& Q1, const T& Q2, const T& Q3, const T& Q4 ) {
  return vect_n< T >( Q1, Q2, Q3, Q4 );
};

template < typename T >
vect_n< T > make_vect_n( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5 ) {
  return vect_n< T >( Q1, Q2, Q3, Q4, Q5 );
};

template < typename T >
vect_n< T > make_vect_n( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6 ) {
  return vect_n< T >( Q1, Q2, Q3, Q4, Q5, Q6 );
};

template < typename T >
vect_n< T > make_vect_n( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7 ) {
  return vect_n< T >( Q1, Q2, Q3, Q4, Q5, Q6, Q7 );
};

template < typename T >
vect_n< T > make_vect_n( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                         const T& Q8 ) {
  return vect_n< T >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8 );
};

template < typename T >
vect_n< T > make_vect_n( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                         const T& Q8, const T& Q9 ) {
  return vect_n< T >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9 );
};

template < typename T >
vect_n< T > make_vect_n( const T& Q1, const T& Q2, const T& Q3, const T& Q4, const T& Q5, const T& Q6, const T& Q7,
                         const T& Q8, const T& Q9, const T& Q10 ) {
  return vect_n< T >( Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10 );
};


template < typename T, unsigned int Size >
template < typename U, typename Allocator >
vect< T, Size >::vect( const vect_n< U, Allocator >& V ) BOOST_NOEXCEPT {
  if( Size > V.size() ) {
    std::copy( V.begin(), V.end(), q );
    std::fill( q + V.size(), q + Size, T() );
  } else {
    std::copy( V.begin(), V.begin() + Size, q );
  };
};


/**
 * This class implements a variable-size templated vector class in which all vector-elements
 * have the same value.
 */
template < typename T, std::size_t Size = 0 >
class vect_scalar {
public:
  typedef vect_scalar< T > self;

  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef std::allocator< T > allocator_type;

  typedef void iterator;
  typedef vect_index_const_iter< self > const_iterator;

  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, dimensions = Size );

private:
  value_type q; /**< Components of the vector. */
public:
  /**
   * Returns the size of the vector.
   */
  size_type size() const BOOST_NOEXCEPT { return Size; };
  /**
   * Returns the max-size of the vector.
   */
  size_type max_size() const BOOST_NOEXCEPT { return Size; };
  /**
   * Returns the capacity of the vector.
   */
  size_type capacity() const BOOST_NOEXCEPT { return Size; };
  /**
   * Resizes the vector.
   */
  void resize( size_type sz, T c = T() ) BOOST_NOEXCEPT{};
  /**
   * Checks if the vector is empty.
   */
  bool empty() const BOOST_NOEXCEPT { return true; };
  /**
   * Reserve a capacity for the vector.
   */
  void reserve( size_type sz ) const BOOST_NOEXCEPT{};

  /**
   * Returns a const-iterator to the first element of the vector.
   */
  const_iterator begin() const BOOST_NOEXCEPT { return const_iterator( *this, 0 ); };
  /**
   * Returns a const-iterator to the one-past-last element of the vector.
   */
  const_iterator end() const BOOST_NOEXCEPT { return const_iterator( *this, Size ); };

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /**
   * Constructor for a given value for the vector.
   * \test PASSED
   */
  explicit vect_scalar( const_reference aFill = value_type() ) BOOST_NOEXCEPT : q( aFill ){};

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator[]( size_type i ) const BOOST_NOEXCEPT { return q; };

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_const_ref_view< self > operator[]( const std::pair< size_type, size_type >& r ) const BOOST_NOEXCEPT {
    return sub( *this )[r];
  };

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator()( size_type i ) const BOOST_NOEXCEPT { return q; };
};


template < typename T, std::size_t Size >
struct is_readable_vector< vect_scalar< T, Size > > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_vector< vect_scalar< T, Size > > type;
};

template < typename T, std::size_t Size >
struct is_writable_vector< vect_scalar< T, Size > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_vector< vect_scalar< T, Size > > type;
};

template < typename T, std::size_t Size >
struct is_resizable_vector< vect_scalar< T, Size > > {
  BOOST_STATIC_CONSTANT( bool, value = ( Size == 0 ) );
  typedef is_resizable_vector< vect_scalar< T, Size > > type;
};

template < typename T, std::size_t Size >
struct has_allocator_vector< vect_scalar< T, Size > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef has_allocator_vector< vect_scalar< T, Size > > type;
};


template < typename T, std::size_t Size >
struct vect_copy< vect_scalar< T, Size > > {
  typedef vect_n< T > type;
};


/**
 * This class implements a variable-size templated vector class in which all vector-elements
 * have the same value.
 */
template < typename T >
class vect_scalar< T, 0 > {
public:
  typedef vect_scalar< T, 0 > self;

  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef std::allocator< T > allocator_type;

  typedef void iterator;
  typedef vect_index_const_iter< self > const_iterator;

  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, dimensions = 0 );

private:
  value_type q; /**< Components of the vector. */
  size_type count;

public:
  /**
   * Returns the size of the vector.
   */
  size_type size() const BOOST_NOEXCEPT { return count; };
  /**
   * Returns the max-size of the vector.
   */
  size_type max_size() const BOOST_NOEXCEPT { return std::numeric_limits< size_type >::max(); };
  /**
   * Returns the capacity of the vector.
   */
  size_type capacity() const BOOST_NOEXCEPT { return std::numeric_limits< size_type >::max(); };
  /**
   * Resizes the vector.
   */
  void resize( size_type sz, T c = T() ) BOOST_NOEXCEPT {
    count = sz;
    q = c;
  };
  /**
   * Checks if the vector is empty.
   */
  bool empty() const BOOST_NOEXCEPT { return ( count == 0 ); };
  /**
   * Reserve a capacity for the vector.
   */
  void reserve( size_type sz ) const BOOST_NOEXCEPT{};

  /**
   * Returns a const-iterator to the first element of the vector.
   */
  const_iterator begin() const BOOST_NOEXCEPT { return const_iterator( *this, 0 ); };
  /**
   * Returns a const-iterator to the one-past-last element of the vector.
   */
  const_iterator end() const BOOST_NOEXCEPT { return const_iterator( *this, count ); };

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/
  /**
   * Default constructor.
   * \test PASSED
   */
  vect_scalar() BOOST_NOEXCEPT : q(), count( 0 ){};

  /**
   * Constructor for a given size or dimension of vector.
   * \test PASSED
   */
  explicit vect_scalar( size_type aSize, const_reference aFill = value_type() ) BOOST_NOEXCEPT : q( aFill ),
                                                                                                 count( aSize ){};

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator[]( size_type i ) const BOOST_NOEXCEPT { return q; };

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_const_ref_view< self > operator[]( const std::pair< size_type, size_type >& r ) const BOOST_NOEXCEPT {
    return sub( *this )[r];
  };

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator()( size_type i ) const BOOST_NOEXCEPT { return q; };
};


/*******************************************************************************
                         Basic Functions
*******************************************************************************/


/**
 * Square magnitude of the vector.
 */
template < typename Vector >
auto norm_2_sqr( const Vector& v )
  BOOST_NOEXCEPT -> typename std::decay< decltype( v[v.size()] * v[v.size()] ) >::type {
  typedef typename vect_traits< Vector >::value_type ValueType;
  typedef typename vect_traits< Vector >::size_type SizeType;
  ValueType sum( 0.0 );
  for( SizeType i = 0; i < v.size(); ++i )
    sum += v[i] * v[i];
  return sum;
};

/**
 * Magnitude of the vector.
 */
template < typename Vector >
auto norm_2( const Vector& v ) BOOST_NOEXCEPT -> typename std::decay< decltype( v[v.size()] ) >::type {
  using std::sqrt;
  return sqrt( norm_2_sqr( v ) );
};

/**
 * Infinite norm of the vector.
 */
template < typename Vector >
auto norm_inf( const Vector& v ) BOOST_NOEXCEPT -> typename std::decay< decltype( v[v.size()] ) >::type {
  typedef typename vect_traits< Vector >::value_type ValueType;
  typedef typename vect_traits< Vector >::size_type SizeType;
  using std::fabs;
  ValueType result( 0.0 );
  for( SizeType i = 0; i < v.size(); ++i )
    if( result < fabs( v[i] ) )
      result = fabs( v[i] );
  return result;
};

/**
 * Square magnitude of the vector.
 */
template < typename Vector >
auto norm_1( const Vector& v ) BOOST_NOEXCEPT -> typename std::decay< decltype( v[v.size()] ) >::type {
  typedef typename vect_traits< Vector >::value_type ValueType;
  typedef typename vect_traits< Vector >::size_type SizeType;
  using std::fabs;
  ValueType sum( 0.0 );
  for( SizeType i = 0; i < v.size(); ++i )
    sum += fabs( v[i] );
  return sum;
};

/**
 * Unit vector in the same direction.
 */
template < typename Vector >
auto unit( const Vector& v ) -> typename std::decay< decltype( v / v[v.size()] ) >::type {
  typename std::decay< decltype( v / v[v.size()] ) >::type result;
  result = v;
  result /= norm_2( result );
  return result;
};

/**
 * Checks if two vectors are colinear.
 * \test PASSED
 */
template < typename Vector1, typename Vector2 >
bool colinear( const Vector1& v1, const Vector2& v2 ) BOOST_NOEXCEPT {
  typedef typename vect_traits< Vector1 >::value_type ValueType;
  using std::fabs;
  ValueType tmp_mag2 = norm_2( v2 );
  ValueType tmp_mag1 = norm_2( v1 );
  ValueType tmp_comb = norm_2( v1 + v2 );
  return ( ( ( tmp_mag1 + tmp_mag2 )
             * ( ValueType( 1.0 ) - ValueType( 10.0 ) * std::numeric_limits< ValueType >::epsilon() ) < tmp_comb )
           || ( fabs( tmp_mag1 - tmp_mag2 ) * ( ValueType( 1.0 ) + std::numeric_limits< ValueType >::epsilon() )
                > tmp_comb ) );
};


/*******************************************************************************
                         Basic Operators
*******************************************************************************/


/**
 * Standard add-and-store operator.
 * \test PASSED
 */
template < typename Vector1, typename Vector2 >
typename boost::enable_if< boost::mpl::and_< is_writable_vector< Vector1 >, is_readable_vector< Vector2 > >,
                           Vector1& >::type
  operator+=( Vector1& v1, const Vector2& v2 ) {
  if( v1.size() != v2.size() )
    throw std::range_error( "Vector size mismatch." );
  typedef typename vect_traits< Vector1 >::size_type SizeType;
  for( SizeType i = 0; i < v1.size(); ++i )
    v1[i] = v1[i] + v2[i];
  return v1;
};

/**
 * Standard sub-and-store operator.
 * \test PASSED
 */
template < typename Vector1, typename Vector2 >
typename boost::enable_if< boost::mpl::and_< is_writable_vector< Vector1 >, is_readable_vector< Vector2 > >,
                           Vector1& >::type
  operator-=( Vector1& v1, const Vector2& v2 ) {
  if( v1.size() != v2.size() )
    throw std::range_error( "Vector size mismatch." );
  typedef typename vect_traits< Vector1 >::size_type SizeType;
  for( SizeType i = 0; i < v1.size(); ++i )
    v1[i] = v1[i] - v2[i];
  return v1;
};

/**
 * Scalar multiply-and-store operator for gain.
 * \test PASSED
 */
template < typename T, typename Vector >
typename boost::enable_if< boost::mpl::and_< is_writable_vector< Vector >,
                                             boost::mpl::not_< is_readable_vector< T > > >,
                           Vector& >::type
  operator*=( Vector& v, const T& S ) BOOST_NOEXCEPT {
  typedef typename vect_traits< Vector >::size_type SizeType;
  typedef typename vect_traits< Vector >::value_type ValueType;
  for( SizeType i = 0; i < v.size(); ++i )
    v[i] *= ValueType( S );
  return v;
};

/**
 * Scalar divide-and-store operator for gain.
 * \test PASSED
 */
template < typename T, typename Vector >
typename boost::enable_if< boost::mpl::and_< is_writable_vector< Vector >,
                                             boost::mpl::not_< is_readable_vector< T > > >,
                           Vector& >::type
  operator/=( Vector& v, const T& S ) BOOST_NOEXCEPT {
  typedef typename vect_traits< Vector >::size_type SizeType;
  typedef typename vect_traits< Vector >::value_type ValueType;
  for( SizeType i = 0; i < v.size(); ++i )
    v[i] /= ValueType( S );
  return v;
};


/**
 * Add two vectors.
 * \test PASSED
 */
template < typename Vector1, typename Vector2 >
typename boost::enable_if< boost::mpl::and_< is_writable_vector< Vector1 >, is_readable_vector< Vector2 > >,
                           typename vect_copy< Vector1 >::type >::type
  operator+( const Vector1& v1, const Vector2& v2 ) {
  if( v1.size() != v2.size() )
    throw std::range_error( "Vector size mismatch." );
  typename vect_copy< Vector1 >::type result;
  result = v1;
  result += v2;
  return result;
};

/**
 * Invert the vector.
 * \test PASSED
 */
template < typename Vector >
typename boost::enable_if< is_readable_vector< Vector >, typename vect_copy< Vector >::type >::type
  operator-( const Vector& v ) {
  typedef typename vect_traits< Vector >::size_type SizeType;
  typename vect_copy< Vector >::type result;
  result = v;
  for( SizeType i = 0; i < v.size(); ++i )
    result[i] = -result[i];
  return result;
};

/**
 * Sub two vectors.
 * \test PASSED
 */
template < typename Vector1, typename Vector2 >
typename boost::enable_if< boost::mpl::and_< is_readable_vector< Vector1 >, is_readable_vector< Vector2 > >,
                           typename vect_copy< Vector1 >::type >::type
  operator-( const Vector1& v1, const Vector2& v2 ) {
  if( v1.size() != v2.size() )
    throw std::range_error( "Vector size mismatch." );
  typename vect_copy< Vector1 >::type result;
  result = v1;
  result -= v2;
  return result;
};

/**
 * Dot Product.
 * \test PASSED
 */
template < typename Vector1, typename Vector2 >
auto operator*( const Vector1& v1, const Vector2& v2 )
  -> typename boost::enable_if< boost::mpl::and_< is_readable_vector< Vector1 >, is_readable_vector< Vector2 > >,
                                decltype( v1[v1.size()] * v2[v2.size()] ) >::type {
  if( v1.size() != v2.size() )
    throw std::range_error( "Vector size mismatch." );
  typedef typename vect_traits< Vector1 >::size_type SizeType;
  typename vect_traits< Vector1 >::value_type result( 0 );
  for( SizeType i = 0; i < v1.size(); ++i )
    result += v1[i] * v2[i];
  return result;
};

/**
 * Scalar-vector product.
 * \test PASSED
 */
template < typename T, typename Vector >
typename boost::enable_if< boost::mpl::and_< is_readable_vector< Vector >, boost::mpl::not_< is_readable_vector< T > >,
                                             boost::mpl::not_< is_readable_matrix< T > > >,
                           typename vect_copy< Vector >::type >::type
  operator*( const Vector& v, const T& S ) {
  typename vect_copy< Vector >::type result;
  result = v;
  result *= S;
  return result;
};

/**
 * Scalar-vector product.
 * \test PASSED
 */
template < typename T, typename Vector >
typename boost::enable_if< boost::mpl::and_< is_readable_vector< Vector >, boost::mpl::not_< is_readable_vector< T > >,
                                             boost::mpl::not_< is_readable_matrix< T > > >,
                           typename vect_copy< Vector >::type >::type
  operator*( const T& S, const Vector& v ) {
  typename vect_copy< Vector >::type result;
  result = v;
  result *= S;
  return result;
};

/**
 * Scalar-vector division.
 * \test PASSED
 */
template < typename T, typename Vector >
typename boost::enable_if< boost::mpl::and_< is_readable_vector< Vector >, boost::mpl::not_< is_readable_vector< T > >,
                                             boost::mpl::not_< is_readable_matrix< T > > >,
                           typename vect_copy< Vector >::type >::type
  operator/( const Vector& v, const T& S ) {
  typename vect_copy< Vector >::type result;
  result = v;
  result /= S;
  return result;
};


/**
 * Element-wise product of two vectors.
 * \test PASSED
 */
template < typename Vector1, typename Vector2 >
typename boost::enable_if< boost::mpl::and_< is_writable_vector< Vector1 >, is_readable_vector< Vector2 > >,
                           typename vect_copy< Vector1 >::type >::type
  elem_product( const Vector1& v1, const Vector2& v2 ) {
  if( v1.size() != v2.size() )
    throw std::range_error( "Vector size mismatch." );
  typedef typename vect_traits< Vector1 >::size_type SizeType;
  typename vect_copy< Vector1 >::type result;
  result = v1;
  for( SizeType i = 0; i < result.size(); ++i )
    result[i] *= v2[i];
  return result;
};


/**
 * Sub two vectors. For functional interfaces.
 * \test PASSED
 */
template < typename T, typename Allocator >
vect_n< T, Allocator > diff( const vect_n< T, Allocator >& v1, const vect_n< T, Allocator >& v2 ) {
  return v1 - v2;
};

/**
 * Add two vectors. For functional interfaces.
 * \test PASSED
 */
template < typename T, typename Allocator >
vect_n< T, Allocator > add( const vect_n< T, Allocator >& v1, const vect_n< T, Allocator >& v2 ) {
  return v1 + v2;
};

/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

/**
 * Equality Comparison operator, component-wise.
 * \test PASSED
 */
template < typename Vector1, typename Vector2 >
typename boost::enable_if< boost::mpl::and_< is_readable_vector< Vector1 >, is_readable_vector< Vector2 > >,
                           bool >::type
  operator==( const Vector1& v1, const Vector2& v2 ) BOOST_NOEXCEPT {
  typedef typename vect_traits< Vector1 >::size_type SizeType;
  if( v1.size() != v2.size() )
    return false;
  for( SizeType i = 0; i < v1.size(); ++i )
    if( v1[i] != v2[i] )
      return false;
  return true;
};

/**
 * Inequality Comparison operator, component-wise.
 * \test PASSED
 */
template < typename Vector1, typename Vector2 >
typename boost::enable_if< boost::mpl::and_< is_readable_vector< Vector1 >, is_readable_vector< Vector2 > >,
                           bool >::type
  operator!=( const Vector1& v1, const Vector2& v2 ) BOOST_NOEXCEPT {
  typedef typename vect_traits< Vector1 >::size_type SizeType;
  if( v1.size() != v2.size() )
    return true;
  for( SizeType i = 0; i < v1.size(); ++i )
    if( v1[i] != v2[i] )
      return true;
  return false;
};

/**
 * Greater-than Comparison operator, Euclidean norm.
 * \test PASSED
 */
template < typename Vector1, typename Vector2 >
typename boost::enable_if< boost::mpl::and_< is_readable_vector< Vector1 >, is_readable_vector< Vector2 > >,
                           bool >::type
  operator>( const Vector1& v1, const Vector2& v2 ) BOOST_NOEXCEPT {
  return ( norm_2_sqr( v1 ) > norm_2_sqr( v2 ) );
};

/**
 * Smaller-than Comparison operator, Euclidean norm.
 * \test PASSED
 */
template < typename Vector1, typename Vector2 >
typename boost::enable_if< boost::mpl::and_< is_readable_vector< Vector1 >, is_readable_vector< Vector2 > >,
                           bool >::type
  operator<( const Vector1& v1, const Vector2& v2 ) BOOST_NOEXCEPT {
  return ( norm_2_sqr( v1 ) < norm_2_sqr( v2 ) );
};

/**
 * Greater-or-equal Comparison operator, Euclidean norm.
 * \test PASSED
 */
template < typename Vector1, typename Vector2 >
typename boost::enable_if< boost::mpl::and_< is_readable_vector< Vector1 >, is_readable_vector< Vector2 > >,
                           bool >::type
  operator>=( const Vector1& v1, const Vector2& v2 ) BOOST_NOEXCEPT {
  return ( norm_2_sqr( v1 ) >= norm_2_sqr( v2 ) );
};

/**
 * Smaller-or-equal Comparison operator, Euclidean norm.
 * \test PASSED
 */
template < typename Vector1, typename Vector2 >
typename boost::enable_if< boost::mpl::and_< is_readable_vector< Vector1 >, is_readable_vector< Vector2 > >,
                           bool >::type
  operator<=( const Vector1& v1, const Vector2& v2 ) BOOST_NOEXCEPT {
  return ( norm_2_sqr( v1 ) <= norm_2_sqr( v2 ) );
};


/**
 * Prints a variable-size vector to an output stream as "(v1; v2; v3; ..; vN)".
 * \test PASSED
 */
template < typename Vector >
typename boost::enable_if< is_readable_vector< Vector >, std::ostream& >::type operator<<( std::ostream& out_stream,
                                                                                           const Vector& V ) {
  typedef typename vect_traits< Vector >::size_type SizeType;
  out_stream << "(";
  if( V.size() > 0 )
    out_stream << V[0];
  for( SizeType i = 1; i < V.size(); ++i )
    out_stream << "; " << V[i];
  return out_stream << ")";
};

/**
 * Reads a variable-size vector to an input stream as "(v1; v2; v3; ..; vN)".
 * \test PASSED
 */
template < typename T >
std::istream& operator>>( std::istream& in_stream, vect_n< T >& V ) {
  typedef typename vect_traits< vect_n< T > >::size_type SizeType;
  std::string tmp_str;
  std::getline( in_stream, tmp_str, '(' ); // skip to opening bracket.
  std::getline( in_stream, tmp_str, ')' ); // read to closing bracket.
  SizeType sz = std::count( tmp_str.begin(), tmp_str.end(), ';' ) + 1;
  std::stringstream ss( tmp_str );
  V.resize( sz );
  std::string tmp2;
  for( SizeType i = 0; ss >> V[i]; ++i )
    std::getline( ss, tmp2, ';' );
  return in_stream;
};

/**
 * Reads a variable-size vector to an input stream as "(v1; v2; v3; ..; vN)".
 * \test PASSED
 */
template < typename T, unsigned int Size >
std::istream& operator>>( std::istream& in_stream, vect< T, Size >& V ) {
  typedef typename vect_traits< vect< T, Size > >::size_type SizeType;
  std::string tmp_str;
  std::getline( in_stream, tmp_str, '(' ); // skip to opening bracket.
  std::getline( in_stream, tmp_str, ')' ); // read to closing bracket.
  std::stringstream ss( tmp_str );
  std::string tmp2;
  for( SizeType i = 0; ( i < Size ) && ( ss >> V[i] ); ++i )
    std::getline( ss, tmp2, ';' );
  return in_stream;
};


#if 0
// #ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

extern template class vect<double,2>;
extern template class vect<double,3>;
extern template class vect<double,4>;
extern template class vect<double,6>;
extern template class vect_n<double>;

extern template class vect<float,2>;
extern template class vect<float,3>;
extern template class vect<float,4>;
extern template class vect<float,6>;
extern template class vect_n<float>;


extern template vect<double,2> diff(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
extern template vect<double,3> diff(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
extern template vect<double,4> diff(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
extern template vect<double,6> diff(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
extern template vect_n<double> diff(const vect_n<double>& v1, const vect_n<double>& v2);

extern template vect<double,2> add(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
extern template vect<double,3> add(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
extern template vect<double,4> add(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
extern template vect<double,6> add(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
extern template vect_n<double> add(const vect_n<double>& v1, const vect_n<double>& v2);

extern template double operator %(const vect<double,2>& V1,const vect<double,2>& V2) BOOST_NOEXCEPT;
extern template vect<double,2> operator %(const double& S, const vect<double,2>& V) BOOST_NOEXCEPT;
extern template vect<double,2> operator %(const vect<double,2>& V,const double& S) BOOST_NOEXCEPT;
extern template vect<double,3> operator %(const vect<double,3>& V1 , const vect<double,3>& V2) BOOST_NOEXCEPT;


extern template vect<double,2>& operator +=(vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
extern template vect<double,3>& operator +=(vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
extern template vect<double,4>& operator +=(vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
extern template vect<double,6>& operator +=(vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
extern template vect_n<double>& operator +=(vect_n<double>& v1, const vect_n<double>& v2);

extern template vect<double,2>& operator -=(vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
extern template vect<double,3>& operator -=(vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
extern template vect<double,4>& operator -=(vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
extern template vect<double,6>& operator -=(vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
extern template vect_n<double>& operator -=(vect_n<double>& v1, const vect_n<double>& v2);

extern template vect<double,2>& operator *=(vect<double,2>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect<double,3>& operator *=(vect<double,3>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect<double,4>& operator *=(vect<double,4>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect<double,6>& operator *=(vect<double,6>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect_n<double>& operator *=(vect_n<double>& v1, const double& v2) BOOST_NOEXCEPT;

extern template vect<double,2>& operator /=(vect<double,2>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect<double,3>& operator /=(vect<double,3>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect<double,4>& operator /=(vect<double,4>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect<double,6>& operator /=(vect<double,6>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect_n<double>& operator /=(vect_n<double>& v1, const double& v2) BOOST_NOEXCEPT;

extern template vect<double,2> operator +(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
extern template vect<double,3> operator +(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
extern template vect<double,4> operator +(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
extern template vect<double,6> operator +(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
extern template vect_n<double> operator +(const vect_n<double>& v1, const vect_n<double>& v2);

extern template vect<double,2> operator -(const vect<double,2>& v1) BOOST_NOEXCEPT;
extern template vect<double,3> operator -(const vect<double,3>& v1) BOOST_NOEXCEPT;
extern template vect<double,4> operator -(const vect<double,4>& v1) BOOST_NOEXCEPT;
extern template vect<double,6> operator -(const vect<double,6>& v1) BOOST_NOEXCEPT;
extern template vect_n<double> operator -(const vect_n<double>& v1);

extern template vect<double,2> operator -(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
extern template vect<double,3> operator -(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
extern template vect<double,4> operator -(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
extern template vect<double,6> operator -(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
extern template vect_n<double> operator -(const vect_n<double>& v1, const vect_n<double>& v2);

extern template double operator *(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
extern template double operator *(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
extern template double operator *(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
extern template double operator *(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
extern template double operator *(const vect_n<double>& v1, const vect_n<double>& v2);

extern template vect<double,2> operator *(const vect<double,2>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect<double,3> operator *(const vect<double,3>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect<double,4> operator *(const vect<double,4>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect<double,6> operator *(const vect<double,6>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect_n<double> operator *(const vect_n<double>& v1, const double& v2);

extern template vect<double,2> operator *(const double& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
extern template vect<double,3> operator *(const double& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
extern template vect<double,4> operator *(const double& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
extern template vect<double,6> operator *(const double& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
extern template vect_n<double> operator *(const double& v1, const vect_n<double>& v2);

extern template vect<double,2> operator /(const vect<double,2>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect<double,3> operator /(const vect<double,3>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect<double,4> operator /(const vect<double,4>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect<double,6> operator /(const vect<double,6>& v1, const double& v2) BOOST_NOEXCEPT;
extern template vect_n<double> operator /(const vect_n<double>& v1, const double& v2);


extern template bool operator ==(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
extern template bool operator ==(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
extern template bool operator ==(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
extern template bool operator ==(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
extern template bool operator ==(const vect_n<double>& v1, const vect_n<double>& v2) BOOST_NOEXCEPT;

extern template bool operator !=(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
extern template bool operator !=(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
extern template bool operator !=(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
extern template bool operator !=(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
extern template bool operator !=(const vect_n<double>& v1, const vect_n<double>& v2) BOOST_NOEXCEPT;

extern template bool operator <=(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
extern template bool operator <=(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
extern template bool operator <=(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
extern template bool operator <=(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
extern template bool operator <=(const vect_n<double>& v1, const vect_n<double>& v2) BOOST_NOEXCEPT;

extern template bool operator >=(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
extern template bool operator >=(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
extern template bool operator >=(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
extern template bool operator >=(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
extern template bool operator >=(const vect_n<double>& v1, const vect_n<double>& v2) BOOST_NOEXCEPT;

extern template bool operator <(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
extern template bool operator <(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
extern template bool operator <(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
extern template bool operator <(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
extern template bool operator <(const vect_n<double>& v1, const vect_n<double>& v2) BOOST_NOEXCEPT;

extern template bool operator >(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
extern template bool operator >(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
extern template bool operator >(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
extern template bool operator >(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
extern template bool operator >(const vect_n<double>& v1, const vect_n<double>& v2) BOOST_NOEXCEPT;


extern template std::ostream& operator <<(std::ostream& out_stream, const vect<double,2>& V);
extern template std::ostream& operator <<(std::ostream& out_stream, const vect<double,3>& V);
extern template std::ostream& operator <<(std::ostream& out_stream, const vect<double,4>& V);
extern template std::ostream& operator <<(std::ostream& out_stream, const vect<double,6>& V);
extern template std::ostream& operator <<(std::ostream& out_stream, const vect_n<double>& V);



extern template vect<float,2> diff(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
extern template vect<float,3> diff(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
extern template vect<float,4> diff(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
extern template vect<float,6> diff(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
extern template vect_n<float> diff(const vect_n<float>& v1, const vect_n<float>& v2);

extern template vect<float,2> add(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
extern template vect<float,3> add(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
extern template vect<float,4> add(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
extern template vect<float,6> add(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
extern template vect_n<float> add(const vect_n<float>& v1, const vect_n<float>& v2);

extern template float operator %(const vect<float,2>& V1,const vect<float,2>& V2) BOOST_NOEXCEPT;
extern template vect<float,2> operator %(const float& S, const vect<float,2>& V) BOOST_NOEXCEPT;
extern template vect<float,2> operator %(const vect<float,2>& V,const float& S) BOOST_NOEXCEPT;
extern template vect<float,3> operator %(const vect<float,3>& V1 , const vect<float,3>& V2) BOOST_NOEXCEPT;


extern template vect<float,2>& operator +=(vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
extern template vect<float,3>& operator +=(vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
extern template vect<float,4>& operator +=(vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
extern template vect<float,6>& operator +=(vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
extern template vect_n<float>& operator +=(vect_n<float>& v1, const vect_n<float>& v2);

extern template vect<float,2>& operator -=(vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
extern template vect<float,3>& operator -=(vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
extern template vect<float,4>& operator -=(vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
extern template vect<float,6>& operator -=(vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
extern template vect_n<float>& operator -=(vect_n<float>& v1, const vect_n<float>& v2);

extern template vect<float,2>& operator *=(vect<float,2>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect<float,3>& operator *=(vect<float,3>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect<float,4>& operator *=(vect<float,4>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect<float,6>& operator *=(vect<float,6>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect_n<float>& operator *=(vect_n<float>& v1, const float& v2) BOOST_NOEXCEPT;

extern template vect<float,2>& operator /=(vect<float,2>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect<float,3>& operator /=(vect<float,3>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect<float,4>& operator /=(vect<float,4>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect<float,6>& operator /=(vect<float,6>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect_n<float>& operator /=(vect_n<float>& v1, const float& v2) BOOST_NOEXCEPT;

extern template vect<float,2> operator +(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
extern template vect<float,3> operator +(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
extern template vect<float,4> operator +(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
extern template vect<float,6> operator +(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
extern template vect_n<float> operator +(const vect_n<float>& v1, const vect_n<float>& v2);

extern template vect<float,2> operator -(const vect<float,2>& v1) BOOST_NOEXCEPT;
extern template vect<float,3> operator -(const vect<float,3>& v1) BOOST_NOEXCEPT;
extern template vect<float,4> operator -(const vect<float,4>& v1) BOOST_NOEXCEPT;
extern template vect<float,6> operator -(const vect<float,6>& v1) BOOST_NOEXCEPT;
extern template vect_n<float> operator -(const vect_n<float>& v1);

extern template vect<float,2> operator -(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
extern template vect<float,3> operator -(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
extern template vect<float,4> operator -(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
extern template vect<float,6> operator -(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
extern template vect_n<float> operator -(const vect_n<float>& v1, const vect_n<float>& v2);

extern template float operator *(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
extern template float operator *(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
extern template float operator *(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
extern template float operator *(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
extern template float operator *(const vect_n<float>& v1, const vect_n<float>& v2);

extern template vect<float,2> operator *(const vect<float,2>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect<float,3> operator *(const vect<float,3>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect<float,4> operator *(const vect<float,4>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect<float,6> operator *(const vect<float,6>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect_n<float> operator *(const vect_n<float>& v1, const float& v2);

extern template vect<float,2> operator *(const float& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
extern template vect<float,3> operator *(const float& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
extern template vect<float,4> operator *(const float& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
extern template vect<float,6> operator *(const float& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
extern template vect_n<float> operator *(const float& v1, const vect_n<float>& v2);

extern template vect<float,2> operator /(const vect<float,2>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect<float,3> operator /(const vect<float,3>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect<float,4> operator /(const vect<float,4>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect<float,6> operator /(const vect<float,6>& v1, const float& v2) BOOST_NOEXCEPT;
extern template vect_n<float> operator /(const vect_n<float>& v1, const float& v2);


extern template bool operator ==(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
extern template bool operator ==(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
extern template bool operator ==(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
extern template bool operator ==(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
extern template bool operator ==(const vect_n<float>& v1, const vect_n<float>& v2) BOOST_NOEXCEPT;

extern template bool operator !=(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
extern template bool operator !=(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
extern template bool operator !=(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
extern template bool operator !=(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
extern template bool operator !=(const vect_n<float>& v1, const vect_n<float>& v2) BOOST_NOEXCEPT;

extern template bool operator <=(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
extern template bool operator <=(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
extern template bool operator <=(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
extern template bool operator <=(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
extern template bool operator <=(const vect_n<float>& v1, const vect_n<float>& v2) BOOST_NOEXCEPT;

extern template bool operator >=(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
extern template bool operator >=(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
extern template bool operator >=(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
extern template bool operator >=(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
extern template bool operator >=(const vect_n<float>& v1, const vect_n<float>& v2) BOOST_NOEXCEPT;

extern template bool operator <(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
extern template bool operator <(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
extern template bool operator <(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
extern template bool operator <(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
extern template bool operator <(const vect_n<float>& v1, const vect_n<float>& v2) BOOST_NOEXCEPT;

extern template bool operator >(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
extern template bool operator >(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
extern template bool operator >(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
extern template bool operator >(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
extern template bool operator >(const vect_n<float>& v1, const vect_n<float>& v2) BOOST_NOEXCEPT;


extern template std::ostream& operator <<(std::ostream& out_stream, const vect<float,2>& V);
extern template std::ostream& operator <<(std::ostream& out_stream, const vect<float,3>& V);
extern template std::ostream& operator <<(std::ostream& out_stream, const vect<float,4>& V);
extern template std::ostream& operator <<(std::ostream& out_stream, const vect<float,6>& V);
extern template std::ostream& operator <<(std::ostream& out_stream, const vect_n<float>& V);


#endif
};


#endif

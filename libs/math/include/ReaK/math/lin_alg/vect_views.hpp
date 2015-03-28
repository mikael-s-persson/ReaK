/**
 * \file vect_views.hpp
 *
 * This library provides a number of class templates to create vector views. A vector
 * view simply means that a sub-part of a vector is used as if it was a vector in its
 * own right. This can be very useful to set sub-parts to other values or to use a
 * sub-part in a vector expression (e.g. applying the operation on the entire vector
 * is not practical or efficient).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date December 2011
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

#ifndef REAK_VECT_VIEWS_HPP
#define REAK_VECT_VIEWS_HPP

#include "vect_concepts.hpp"
#include "vect_traits.hpp"
#include <boost/static_assert.hpp>
#include <boost/utility/enable_if.hpp>

#include <stdexcept>

namespace ReaK {

/**
 * This function can be used to generate a range of indices for the matrices.
 * A range is simple a pair of first and last indices (notice that the last index is
 * included in the range).
 */
inline std::pair< std::size_t, std::size_t > range( std::size_t aFirst, std::size_t aLast ) {
  return std::pair< std::size_t, std::size_t >( aFirst, aLast );
};


/**
 * This class template constructs a sub-vector which represents part of the vector.
 * This class takes a const reference to the given vector.
 * \tparam Vector A readable vector type.
 */
template < typename Vector >
class vect_const_ref_view {
public:
  typedef vect_const_ref_view< Vector > self;
  typedef typename vect_traits< Vector >::allocator_type allocator_type;

  typedef typename vect_traits< Vector >::value_type value_type;

  typedef typename vect_traits< Vector >::reference reference;
  typedef typename vect_traits< Vector >::const_reference const_reference;
  typedef typename vect_traits< Vector >::pointer pointer;
  typedef typename vect_traits< Vector >::const_pointer const_pointer;

  typedef typename vect_traits< Vector >::iterator iterator;
  typedef typename vect_traits< Vector >::const_iterator const_iterator;

  typedef typename vect_traits< Vector >::size_type size_type;
  typedef typename vect_traits< Vector >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, dimensions = vect_traits< Vector >::dimensions );

  BOOST_CONCEPT_ASSERT( ( ReadableVectorConcept< Vector > ) );

private:
  const Vector* v;
  size_type offset;
  size_type count;

  self& operator=( const self& );
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
  explicit vect_const_ref_view( Vector&& );

  vect_const_ref_view( Vector&&, size_type, size_type aOffset = 0 );
#endif

public:
  /**
   * Constructs the sub-vector which represents the entire vector.
   * \param aV The vector from which the sub-part is taken.
   */
  explicit vect_const_ref_view( const Vector& aV ) BOOST_NOEXCEPT : v( &aV ), offset( 0 ), count( aV.size() ){};

  /**
   * Constructs the sub-vector which represents part of the vector.
   * \param aV The vector from which the sub-part is taken.
   * \param aCount The number of elements for the sub-vector.
   * \param aOffset The offset from the start of the vector.
   */
  vect_const_ref_view( const Vector& aV, size_type aCount, size_type aOffset = 0 ) BOOST_NOEXCEPT : v( &aV ),
                                                                                                    offset( aOffset ),
                                                                                                    count( aCount ){};


  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Vector indexing accessor for read-only access.
   * \param i Index.
   * \return the element at the given position.
   * TEST PASSED
   */
  value_type operator[]( size_type i ) const BOOST_NOEXCEPT { return ( *v )[offset + i]; };

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_const_ref_view< self > operator[]( const std::pair< size_type, size_type >& r ) const BOOST_NOEXCEPT {
    return vect_const_ref_view< self >( *this, r.second - r.first, r.first );
  };

  /**
   * Vector indexing operator, accessor for read only.
   * TEST PASSED
   */
  value_type operator()( size_type i ) const BOOST_NOEXCEPT { return ( *v )[offset + i]; };

  /**
   * Gets the size of the vector.
   * \return number of elements of the vector.
   * TEST PASSED
   */
  size_type size() const BOOST_NOEXCEPT { return count; };
  /**
   * Returns the max-size of the vector.
   */
  size_type max_size() const BOOST_NOEXCEPT { return count; };
  /**
   * Returns the capacity of the vector.
   */
  size_type capacity() const BOOST_NOEXCEPT { return count; };
  /**
   * Checks if the vector is empty.
   */
  bool empty() const BOOST_NOEXCEPT { return ( count == 0 ); };

  /**
   * Returns a const-iterator to the first element of the vector.
   */
  const_iterator begin() const BOOST_NOEXCEPT { return v->begin() + offset; };
  /**
   * Returns a const-iterator to the one-past-last element of the vector.
   */
  const_iterator end() const BOOST_NOEXCEPT { return v->begin() + offset + count; };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return v->get_allocator(); };
};


template < typename Vector >
struct is_readable_vector< vect_const_ref_view< Vector > > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_vector< Vector >::value );
  typedef is_readable_vector< Vector > type;
};

template < typename Vector >
struct is_writable_vector< vect_const_ref_view< Vector > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_vector< vect_const_ref_view< Vector > > type;
};

template < typename Vector >
struct is_resizable_vector< vect_const_ref_view< Vector > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_vector< vect_const_ref_view< Vector > > type;
};

template < typename Vector >
struct has_allocator_vector< vect_const_ref_view< Vector > > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_vector< Vector >::value );
  typedef has_allocator_vector< Vector > type;
};


/**
 * This class template constructs a sub-vector which represents part of the vector.
 * This class takes a reference to the given vector.
 * \tparam Vector A readable vector type.
 */
template < typename Vector >
class vect_ref_view {
public:
  typedef vect_ref_view< Vector > self;
  typedef typename vect_traits< Vector >::allocator_type allocator_type;

  typedef typename vect_traits< Vector >::value_type value_type;

  typedef typename vect_traits< Vector >::reference reference;
  typedef typename vect_traits< Vector >::const_reference const_reference;
  typedef typename vect_traits< Vector >::pointer pointer;
  typedef typename vect_traits< Vector >::const_pointer const_pointer;

  typedef typename vect_traits< Vector >::iterator iterator;
  typedef typename vect_traits< Vector >::const_iterator const_iterator;

  typedef typename vect_traits< Vector >::size_type size_type;
  typedef typename vect_traits< Vector >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, dimensions = vect_traits< Vector >::dimensions );

  BOOST_CONCEPT_ASSERT( ( ReadableVectorConcept< Vector > ) );

private:
  Vector* v;
  size_type offset;
  size_type count;

public:
  /**
   * Constructs the sub-vector which represents the entire vector.
   * \param aV The vector from which the sub-part is taken.
   */
  explicit vect_ref_view( Vector& aV ) BOOST_NOEXCEPT : v( &aV ), offset( 0 ), count( aV.size() ){};

  /**
   * Constructs the sub-vector which represents part of the vector.
   * \param aV The vector from which the sub-part is taken.
   * \param aCount The number of elements for the sub-part.
   * \param aOffset The offset from the start of the vector.
   */
  vect_ref_view( Vector& aV, size_type aCount, size_type aOffset = 0 ) BOOST_NOEXCEPT : v( &aV ),
                                                                                        offset( aOffset ),
                                                                                        count( aCount ){};

  /**
   * Standard assignment operator.
   */
  template < typename Vector2 >
  typename boost::enable_if< is_readable_vector< Vector2 >, self& >::type operator=( const Vector2& rhs ) {
    if( rhs.size() != count )
      throw std::range_error( "Vector dimensions mismatch." );
    for( size_type i = 0; i < count; ++i )
      ( *v )[offset + i] = rhs[i];
    return *this;
  };

  /**
   * Standard assignment operator.
   */
  template < typename Vector2 >
  typename boost::enable_if< is_readable_vector< Vector2 >, self& >::type operator+=( const Vector2& rhs ) {
    if( rhs.size() != count )
      throw std::range_error( "Vector dimensions mismatch." );
    for( size_type i = 0; i < count; ++i )
      ( *v )[offset + i] += rhs[i];
    return *this;
  };

  /**
   * Standard assignment operator.
   */
  template < typename Vector2 >
  typename boost::enable_if< is_readable_vector< Vector2 >, self& >::type operator-=( const Vector2& rhs ) {
    if( rhs.size() != count )
      throw std::range_error( "Vector dimensions mismatch." );
    for( size_type i = 0; i < count; ++i )
      ( *v )[offset + i] -= rhs[i];
    return *this;
  };

  /**
   * Standard assignment operator.
   */
  template < typename Scalar >
  self& operator*=( const Scalar& rhs ) BOOST_NOEXCEPT {
    for( size_type i = 0; i < count; ++i )
      ( *v )[offset + i] *= rhs;
    return *this;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Vector indexing accessor for read-write access.
   * \param i Index.
   * \return the element at the given position.
   * TEST PASSED
   */
  reference operator[](size_type i)BOOST_NOEXCEPT { return ( *v )[offset + i]; };
  /**
   * Vector indexing accessor for read-only access.
   * \param i Index.
   * \return the element at the given position.
   * TEST PASSED
   */
  value_type operator[]( size_type i ) const BOOST_NOEXCEPT { return ( *v )[offset + i]; };

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_ref_view< self > operator[](const std::pair< size_type, size_type >& r)BOOST_NOEXCEPT {
    return vect_ref_view< self >( *this, r.second - r.first, r.first );
  };

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_const_ref_view< self > operator[]( const std::pair< size_type, size_type >& r ) const BOOST_NOEXCEPT {
    return vect_const_ref_view< self >( *this, r.second - r.first, r.first );
  };

  /**
   * Vector indexing operator, accessor for read/write.
   * TEST PASSED
   */
  reference operator()( size_type i ) BOOST_NOEXCEPT { return ( *v )[offset + i]; };

  /**
   * Vector indexing operator, accessor for read only.
   * TEST PASSED
   */
  value_type operator()( size_type i ) const BOOST_NOEXCEPT { return ( *v )[offset + i]; };

  /**
   * Gets the size of the vector.
   * \return number of elements of the vector.
   * TEST PASSED
   */
  size_type size() const BOOST_NOEXCEPT { return count; };
  /**
   * Returns the max-size of the vector.
   */
  size_type max_size() const BOOST_NOEXCEPT { return count; };
  /**
   * Returns the capacity of the vector.
   */
  size_type capacity() const BOOST_NOEXCEPT { return count; };
  /**
   * Resizes the vector.
   */
  void resize( size_type sz, value_type c = value_type() ) const BOOST_NOEXCEPT{};
  /**
   * Checks if the vector is empty.
   */
  bool empty() const BOOST_NOEXCEPT { return ( count == 0 ); };
  /**
   * Reserve a capacity for the vector.
   */
  void reserve( size_type sz ) const BOOST_NOEXCEPT{};

  /**
   * Returns an iterator to the first element of the vector.
   */
  iterator begin() BOOST_NOEXCEPT { return v->begin() + offset; };
  /**
   * Returns a const-iterator to the first element of the vector.
   */
  const_iterator begin() const BOOST_NOEXCEPT { return v->begin() + offset; };
  /**
   * Returns an iterator to the one-past-last element of the vector.
   */
  iterator end() BOOST_NOEXCEPT { return v->begin() + offset + count; };
  /**
   * Returns a const-iterator to the one-past-last element of the vector.
   */
  const_iterator end() const BOOST_NOEXCEPT { return v->begin() + offset + count; };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return v->get_allocator(); };
};


template < typename Vector >
struct is_readable_vector< vect_ref_view< Vector > > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_vector< Vector >::value );
  typedef is_readable_vector< Vector > type;
};

template < typename Vector >
struct is_writable_vector< vect_ref_view< Vector > > {
  BOOST_STATIC_CONSTANT( bool, value = is_writable_vector< Vector >::value );
  typedef is_writable_vector< Vector > type;
};

template < typename Vector >
struct is_resizable_vector< vect_ref_view< Vector > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_vector< vect_ref_view< Vector > > type;
};

template < typename Vector >
struct has_allocator_vector< vect_ref_view< Vector > > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_vector< Vector >::value );
  typedef has_allocator_vector< Vector > type;
};


/**
 * This class template constructs a sub-vector which represents part of the vector.
 * This class makes a copy of the given vector (it is mainly meant to harmonize syntax
 * when rvalue vectors are involved, requires C++11).
 * \tparam Vector A readable vector type.
 */
template < typename Vector >
class vect_copy_view {
public:
  typedef vect_copy_view< Vector > self;
  typedef typename vect_traits< Vector >::allocator_type allocator_type;

  typedef typename vect_traits< Vector >::value_type value_type;

  typedef typename vect_traits< Vector >::reference reference;
  typedef typename vect_traits< Vector >::const_reference const_reference;
  typedef typename vect_traits< Vector >::pointer pointer;
  typedef typename vect_traits< Vector >::const_pointer const_pointer;

  typedef typename vect_traits< Vector >::iterator iterator;
  typedef typename vect_traits< Vector >::const_iterator const_iterator;

  typedef typename vect_traits< Vector >::size_type size_type;
  typedef typename vect_traits< Vector >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, dimensions = vect_traits< Vector >::dimensions );

  BOOST_CONCEPT_ASSERT( ( ReadableVectorConcept< Vector > ) );

private:
  Vector v;
  size_type offset;
  size_type count;

public:
  /**
   * Default constructor.
   */
  vect_copy_view() : v(), offset( 0 ), count( 0 ){};

  /**
   * Constructs the sub-vector which represents the entire vector.
   */
  explicit vect_copy_view( const Vector& aV ) : v( aV ), offset( 0 ), count( aV.size() ){};

  /**
   * Constructs the sub-vector which represents part of the vector.
   * \param aV The vector from which the sub-block is taken.
   * \param aCount The number of elements for the sub-vector.
   * \param aOffset The offset from the start of the vector.
   */
  vect_copy_view( const Vector& aV, size_type aCount, size_type aOffset = 0 )
      : v( aV ), offset( aOffset ), count( aCount ){};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit vect_copy_view( Vector&& aV ) : v( std::move( aV ) ), offset( 0 ), count( 0 ) { count = v.size(); };

  /**
   * Constructs the sub-vector which represents part of the vector.
   * \param aV The vector from which the sub-block is taken.
   * \param aCount The number of elements for the sub-vector.
   * \param aOffset The offset from the start of the vector.
   */
  vect_copy_view( Vector&& aV, size_type aCount, size_type aOffset = 0 )
      : v( std::move( aV ) ), offset( aOffset ), count( aCount ){};
#endif

  /**
   * Standard swap function.
   */
  friend void swap( self& lhs, self& rhs ) BOOST_NOEXCEPT {
    using std::swap;
    swap( lhs.v, rhs.v );
    swap( lhs.offset, rhs.offset );
    swap( lhs.count, rhs.count );
    return;
  };

  /**
   * Standard assignment operator.
   */
  self& operator=( self rhs ) BOOST_NOEXCEPT {
    swap( *this, rhs );
    return *this;
  };

  /**
   * Standard assignment operator.
   */
  template < typename Vector2 >
  typename boost::enable_if_c< is_readable_vector< Vector2 >::value && !boost::is_same< Vector2, self >::value,
                               self& >::type
    operator=( const Vector2& rhs ) {
    if( rhs.size() != count )
      throw std::range_error( "Vector dimensions mismatch." );
    for( size_type i = 0; i < count; ++i )
      v[offset + i] = rhs[i];
    return *this;
  };

  /**
   * Standard assignment operator.
   */
  template < typename Vector2 >
  typename boost::enable_if_c< is_readable_vector< Vector2 >::value, self& >::type operator+=( const Vector2& rhs ) {
    if( rhs.size() != count )
      throw std::range_error( "Vector dimensions mismatch." );
    for( size_type i = 0; i < count; ++i )
      v[offset + i] += rhs[i];
    return *this;
  };

  /**
   * Standard assignment operator.
   */
  template < typename Vector2 >
  typename boost::enable_if_c< is_readable_vector< Vector2 >::value, self& >::type operator-=( const Vector2& rhs ) {
    if( rhs.size() != count )
      throw std::range_error( "Vector dimensions mismatch." );
    for( size_type i = 0; i < count; ++i )
      v[offset + i] -= rhs[i];
    return *this;
  };

  /**
   * Standard assignment operator.
   */
  template < typename Scalar >
  self& operator*=( const Scalar& rhs ) BOOST_NOEXCEPT {
    for( size_type i = 0; i < count; ++i )
      v[offset + i] *= rhs;
    return *this;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Vector indexing accessor for read-write access.
   * \param i Index.
   * \return the element at the given position.
   * TEST PASSED
   */
  reference operator[](size_type i)BOOST_NOEXCEPT { return v[offset + i]; };
  /**
   * Vector indexing accessor for read-only access.
   * \param i Index.
   * \return the element at the given position.
   * TEST PASSED
   */
  value_type operator[]( size_type i ) const BOOST_NOEXCEPT { return v[offset + i]; };

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_ref_view< self > operator[](const std::pair< size_type, size_type >& r)BOOST_NOEXCEPT {
    return vect_ref_view< self >( *this, r.second - r.first, r.first );
  };

  /**
   * Sub-vector operator, accessor for read only.
   * \test PASSED
   */
  vect_const_ref_view< self > operator[]( const std::pair< size_type, size_type >& r ) const BOOST_NOEXCEPT {
    return vect_const_ref_view< self >( *this, r.second - r.first, r.first );
  };

  /**
   * Vector indexing operator, accessor for read/write.
   * TEST PASSED
   */
  reference operator()( size_type i ) BOOST_NOEXCEPT { return v[offset + i]; };

  /**
   * Vector indexing operator, accessor for read only.
   * TEST PASSED
   */
  value_type operator()( size_type i ) const BOOST_NOEXCEPT { return v[offset + i]; };

  /**
   * Gets the size of the vector.
   * \return number of elements of the vector.
   * TEST PASSED
   */
  size_type size() const BOOST_NOEXCEPT { return count; };
  /**
   * Returns the max-size of the vector.
   */
  size_type max_size() const BOOST_NOEXCEPT { return count; };
  /**
   * Returns the capacity of the vector.
   */
  size_type capacity() const BOOST_NOEXCEPT { return count; };
  /**
   * Resizes the vector.
   */
  void resize( size_type sz, value_type c = value_type() ) const BOOST_NOEXCEPT{};
  /**
   * Checks if the vector is empty.
   */
  bool empty() const BOOST_NOEXCEPT { return ( count == 0 ); };
  /**
   * Reserve a capacity for the vector.
   */
  void reserve( size_type sz ) const BOOST_NOEXCEPT{};

  /**
   * Returns an iterator to the first element of the vector.
   */
  iterator begin() BOOST_NOEXCEPT { return v.begin() + offset; };
  /**
   * Returns a const-iterator to the first element of the vector.
   */
  const_iterator begin() const BOOST_NOEXCEPT { return v.begin() + offset; };
  /**
   * Returns an iterator to the one-past-last element of the vector.
   */
  iterator end() BOOST_NOEXCEPT { return v.begin() + offset + count; };
  /**
   * Returns a const-iterator to the one-past-last element of the vector.
   */
  const_iterator end() const BOOST_NOEXCEPT { return v.begin() + offset + count; };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return v.get_allocator(); };
};


template < typename Vector >
struct is_readable_vector< vect_copy_view< Vector > > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_vector< Vector >::value );
  typedef is_readable_vector< Vector > type;
};

template < typename Vector >
struct is_writable_vector< vect_copy_view< Vector > > {
  BOOST_STATIC_CONSTANT( bool, value = is_writable_vector< Vector >::value );
  typedef is_writable_vector< Vector > type;
};

template < typename Vector >
struct is_resizable_vector< vect_copy_view< Vector > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_vector< vect_copy_view< Vector > > type;
};

template < typename Vector >
struct has_allocator_vector< vect_copy_view< Vector > > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_vector< Vector >::value );
  typedef has_allocator_vector< Vector > type;
};


template < typename Vector >
struct vect_copy_view_factory {
  typedef typename vect_traits< Vector >::size_type size_type;

  Vector v;
  vect_copy_view_factory( const Vector& aV ) : v( aV ){};
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
  vect_copy_view< Vector > operator[]( const std::pair< size_type, size_type >& indices ) {
    return vect_copy_view< Vector >( v, indices.second - indices.first, indices.first );
  };
#else
  vect_copy_view_factory( Vector&& aV ) : v( std::move( aV ) ){};
  vect_copy_view< Vector > operator[]( const std::pair< size_type, size_type >& indices ) {
    return vect_copy_view< Vector >( std::move( v ), indices.second - indices.first, indices.first );
  };
#endif
};

template < typename Vector >
struct vect_ref_view_factory {
  typedef typename vect_traits< Vector >::size_type size_type;

  Vector& v;
  vect_ref_view_factory( Vector& aV ) BOOST_NOEXCEPT : v( aV ){};
  vect_ref_view< Vector > operator[](const std::pair< size_type, size_type >& indices)BOOST_NOEXCEPT {
    return vect_ref_view< Vector >( v, indices.second - indices.first, indices.first );
  };
};

template < typename Vector >
struct vect_const_ref_view_factory {
  typedef typename vect_traits< Vector >::size_type size_type;

  const Vector& v;
  vect_const_ref_view_factory( const Vector& aV ) BOOST_NOEXCEPT : v( aV ){};
  vect_const_ref_view< Vector > operator[](const std::pair< size_type, size_type >& indices)BOOST_NOEXCEPT {
    return vect_const_ref_view< Vector >( v, indices.second - indices.first, indices.first );
  };
};


template < typename Vector >
typename boost::enable_if< is_readable_vector< Vector >, vect_ref_view_factory< Vector > >::type
  sub( Vector& V ) BOOST_NOEXCEPT {
  return vect_ref_view_factory< Vector >( V );
};

template < typename Vector >
typename boost::enable_if< is_readable_vector< Vector >, vect_const_ref_view_factory< Vector > >::type
  sub( const Vector& V ) BOOST_NOEXCEPT {
  return vect_const_ref_view_factory< Vector >( V );
};

template < typename Vector >
typename boost::enable_if< is_readable_vector< Vector >, vect_copy_view_factory< Vector > >::type
  sub_copy( const Vector& V ) {
  return vect_copy_view_factory< Vector >( V );
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
template < typename Vector >
typename boost::enable_if< is_readable_vector< Vector >, vect_copy_view_factory< Vector > >::type sub( Vector&& V ) {
  return vect_copy_view_factory< Vector >( std::move( V ) );
};
#endif
};

#endif

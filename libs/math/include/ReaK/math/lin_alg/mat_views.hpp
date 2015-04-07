/**
 * \file mat_views.hpp
 *
 * This library provides a number of class templates to create matrix views. A matrix
 * view simply means that a sub-block of a matrix is used as if it was a matrix in its
 * own right. This can be very useful to set sub-blocks to other values or to use a
 * sub-block in a matrix expression (e.g. applying the operation on the entire matrix
 * is not practical or efficient).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2011
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

#ifndef REAK_MAT_VIEWS_HPP
#define REAK_MAT_VIEWS_HPP

#include "mat_concepts.hpp"
#include "mat_traits.hpp"
#include "vect_concepts.hpp"
#include "vect_views.hpp"

#include <boost/static_assert.hpp>

namespace ReaK {


/**
 * This class template
 *
 *
 */
template < typename Matrix >
class mat_copy_sub_block {
public:
  typedef mat_copy_sub_block< Matrix > self;
  typedef typename mat_traits< Matrix >::allocator_type allocator_type;

  typedef typename mat_traits< Matrix >::value_type value_type;

  typedef typename mat_traits< Matrix >::reference reference;
  typedef typename mat_traits< Matrix >::const_reference const_reference;
  typedef typename mat_traits< Matrix >::pointer pointer;
  typedef typename mat_traits< Matrix >::const_pointer const_pointer;

  typedef typename mat_traits< Matrix >::col_iterator col_iterator;
  typedef typename mat_traits< Matrix >::const_col_iterator const_col_iterator;
  typedef typename mat_traits< Matrix >::row_iterator row_iterator;
  typedef typename mat_traits< Matrix >::const_row_iterator const_row_iterator;

  typedef typename mat_traits< Matrix >::size_type size_type;
  typedef typename mat_traits< Matrix >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 0 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_traits< Matrix >::alignment );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_structure::rectangular );

private:
  Matrix m;
  size_type rowOffset;
  size_type colOffset;
  size_type rowCount;
  size_type colCount;

public:
  /**
   * Default constructor.
   */
  mat_copy_sub_block() : m(), rowOffset( 0 ), colOffset( 0 ), rowCount( 0 ), colCount( 0 ){};

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_copy_sub_block( const Matrix& aM )
      : m( aM ), rowOffset( 0 ), colOffset( 0 ), rowCount( aM.get_row_count() ), colCount( aM.get_col_count() ){};

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aRowCount The number of rows for the sub-block.
   * \param aColCount The number of columns for the sub-block.
   * \param aRowOffset The row-offset from the start of the matrix.
   * \param aColOffset The column-offset from the start of the matrix.
   */
  mat_copy_sub_block( const Matrix& aM, size_type aRowCount, size_type aColCount, size_type aRowOffset = 0,
                      size_type aColOffset = 0 )
      : m( aM ), rowOffset( aRowOffset ), colOffset( aColOffset ), rowCount( aRowCount ), colCount( aColCount ){};

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_copy_sub_block( Matrix&& aM )
      : m( std::move( aM ) ), rowOffset( 0 ), colOffset( 0 ), rowCount( 0 ), colCount( 0 ) {
    rowCount = m.get_row_count();
    colCount = m.get_col_count();
  };

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aRowCount The number of rows for the sub-block.
   * \param aColCount The number of columns for the sub-block.
   * \param aRowOffset The row-offset from the start of the matrix.
   * \param aColOffset The column-offset from the start of the matrix.
   */
  mat_copy_sub_block( Matrix&& aM, size_type aRowCount, size_type aColCount, size_type aRowOffset = 0,
                      size_type aColOffset = 0 )
      : m( std::move( aM ) ), rowOffset( aRowOffset ), colOffset( aColOffset ), rowCount( aRowCount ),
        colCount( aColCount ){};

  /**
   * Standard copy-constructor.
   */
  mat_copy_sub_block( const self& aObj )
      : m( aObj.m ), rowOffset( aObj.rowOffset ), colOffset( aObj.colOffset ), rowCount( aObj.rowCount ),
        colCount( aObj.colCount ){};

  /**
   * Standard move-constructor.
   */
  mat_copy_sub_block( self&& aObj )
      : m( std::move( aObj.m ) ), rowOffset( aObj.rowOffset ), colOffset( aObj.colOffset ), rowCount( aObj.rowCount ),
        colCount( aObj.colCount ){};

  /**
   * Standard swap function.
   */
  friend void swap( self& lhs, self& rhs ) throw() {
    using std::swap;
    swap( lhs.m, rhs.m );
    swap( lhs.rowOffset, rhs.rowOffset );
    swap( lhs.colOffset, rhs.colOffset );
    swap( lhs.rowCount, rhs.rowCount );
    swap( lhs.colCount, rhs.colCount );
    return;
  };

  /**
   * Standard assignment operator.
   */
  self& operator=( self rhs ) {
    swap( *this, rhs );
    return *this;
  };

  /**
   * Standard assignment operator.
   */
  template < typename Matrix2 >
  typename boost::enable_if< boost::mpl::and_< is_readable_matrix< Matrix2 >,
                                               boost::mpl::not_< boost::is_same< Matrix2, self > > >,
                             self& >::type
    operator=( const Matrix2& rhs ) {
    if( ( rhs.get_row_count() != rowCount ) || ( rhs.get_col_count() != colCount ) )
      throw std::range_error( "Matrix dimensions mismatch." );
    for( size_type j = 0; j < colCount; ++j )
      for( size_type i = 0; i < rowCount; ++i )
        m( rowOffset + i, colOffset + j ) = rhs( i, j );
    return *this;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-write access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  reference operator()( size_type i, size_type j ) { return m( rowOffset + i, colOffset + j ); };
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()( size_type i, size_type j ) const { return m( rowOffset + i, colOffset + j ); };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const throw() { return rowCount; };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const throw() { return colCount; };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() { return std::make_pair( rowCount, colCount ); };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return m.get_allocator(); };


  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator+=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != colCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type j = 0; j < colCount; ++j )
      for( size_type i = 0; i < rowCount; ++i )
        m( rowOffset + i, colOffset + j ) += M( i, j );
    return *this;
  };

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator-=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != colCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type j = 0; j < colCount; ++j )
      for( size_type i = 0; i < rowCount; ++i )
        m( rowOffset + i, colOffset + j ) -= M( i, j );
    return *this;
  };

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=( const value_type& S ) {
    for( size_type j = 0; j < colCount; ++j )
      for( size_type i = 0; i < rowCount; ++i )
        m( rowOffset + i, colOffset + j ) *= S;
    return *this;
  };

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator*=( const Matrix2& M ) {
    if( ( M.get_col_count() != colCount ) || ( M.get_row_count() != colCount ) )
      throw std::range_error( "Matrix Dimension Mismatch." );
    *this = *this * M;
    return *this;
  };
};


template < typename Matrix >
struct is_readable_matrix< mat_copy_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix< Matrix >::value );
  typedef is_readable_matrix< Matrix > type;
};

template < typename Matrix >
struct is_writable_matrix< mat_copy_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix< Matrix >::value );
  typedef is_fully_writable_matrix< Matrix > type;
};

template < typename Matrix >
struct is_fully_writable_matrix< mat_copy_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix< Matrix >::value );
  typedef is_fully_writable_matrix< Matrix > type;
};

template < typename Matrix >
struct is_resizable_matrix< mat_copy_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_copy_sub_block< Matrix > > type;
};

template < typename Matrix >
struct has_allocator_matrix< mat_copy_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix< Matrix >::value );
  typedef has_allocator_matrix< Matrix > type;
};


template < typename Matrix >
class mat_sub_block {
public:
  typedef mat_sub_block< Matrix > self;
  typedef typename mat_traits< Matrix >::allocator_type allocator_type;

  typedef typename mat_traits< Matrix >::value_type value_type;

  typedef typename mat_traits< Matrix >::reference reference;
  typedef typename mat_traits< Matrix >::const_reference const_reference;
  typedef typename mat_traits< Matrix >::pointer pointer;
  typedef typename mat_traits< Matrix >::const_pointer const_pointer;

  typedef typename mat_traits< Matrix >::col_iterator col_iterator;
  typedef typename mat_traits< Matrix >::const_col_iterator const_col_iterator;
  typedef typename mat_traits< Matrix >::row_iterator row_iterator;
  typedef typename mat_traits< Matrix >::const_row_iterator const_row_iterator;

  typedef typename mat_traits< Matrix >::size_type size_type;
  typedef typename mat_traits< Matrix >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 0 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_traits< Matrix >::alignment );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_structure::rectangular );

private:
  Matrix* m;
  size_type rowOffset;
  size_type colOffset;
  size_type rowCount;
  size_type colCount;

public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_sub_block( Matrix& aM )
      : m( &aM ), rowOffset( 0 ), colOffset( 0 ), rowCount( aM.get_row_count() ), colCount( aM.get_col_count() ){};

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aRowCount The number of rows for the sub-block.
   * \param aColCount The number of columns for the sub-block.
   * \param aRowOffset The row-offset from the start of the matrix.
   * \param aColOffset The column-offset from the start of the matrix.
   */
  mat_sub_block( Matrix& aM, size_type aRowCount, size_type aColCount, size_type aRowOffset = 0,
                 size_type aColOffset = 0 )
      : m( &aM ), rowOffset( aRowOffset ), colOffset( aColOffset ), rowCount( aRowCount ), colCount( aColCount ){};

  /**
   * Standard copy-constructor.
   */
  mat_sub_block( const self& aObj )
      : m( aObj.m ), rowOffset( aObj.rowOffset ), colOffset( aObj.colOffset ), rowCount( aObj.rowCount ),
        colCount( aObj.colCount ){};

  /**
   * Standard move-constructor.
   */
  mat_sub_block( self&& aObj )
      : m( aObj.m ), rowOffset( aObj.rowOffset ), colOffset( aObj.colOffset ), rowCount( aObj.rowCount ),
        colCount( aObj.colCount ){};

  /**
   * Standard swap function (shallow).
   */
  friend void swap( self& lhs, self& rhs ) throw() {
    using std::swap;
    swap( lhs.m, rhs.m );
    swap( lhs.rowOffset, rhs.rowOffset );
    swap( lhs.colOffset, rhs.colOffset );
    swap( lhs.rowCount, rhs.rowCount );
    swap( lhs.colCount, rhs.colCount );
    return;
  };

  /**
   * Standard assignment operator.
   */
  self& operator=( const self& rhs ) {
    if( ( rhs.get_row_count() != rowCount ) || ( rhs.get_col_count() != colCount ) )
      throw std::range_error( "Matrix dimensions mismatch." );
    for( size_type j = 0; j < colCount; ++j )
      for( size_type i = 0; i < rowCount; ++i )
        ( *m )( rowOffset + i, colOffset + j ) = rhs( i, j );
    return *this;
  };

  /**
   * Standard assignment operator.
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator=( const Matrix2& rhs ) {
    if( ( rhs.get_row_count() != rowCount ) || ( rhs.get_col_count() != colCount ) )
      throw std::range_error( "Matrix dimensions mismatch." );
    for( size_type j = 0; j < colCount; ++j )
      for( size_type i = 0; i < rowCount; ++i )
        ( *m )( rowOffset + i, colOffset + j ) = rhs( i, j );
    return *this;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-write access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  reference operator()( size_type i, size_type j ) { return ( *m )( rowOffset + i, colOffset + j ); };
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()( size_type i, size_type j ) const { return ( *m )( rowOffset + i, colOffset + j ); };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const throw() { return rowCount; };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const throw() { return colCount; };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() { return std::make_pair( rowCount, colCount ); };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return m->get_allocator(); };


  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator+=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != colCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type j = 0; j < colCount; ++j )
      for( size_type i = 0; i < rowCount; ++i )
        ( *m )( rowOffset + i, colOffset + j ) += M( i, j );
    return *this;
  };

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator-=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != colCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type j = 0; j < colCount; ++j )
      for( size_type i = 0; i < rowCount; ++i )
        ( *m )( rowOffset + i, colOffset + j ) -= M( i, j );
    return *this;
  };

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=( const value_type& S ) {
    for( size_type j = 0; j < colCount; ++j )
      for( size_type i = 0; i < rowCount; ++i )
        ( *m )( rowOffset + i, colOffset + j ) *= S;
    return *this;
  };

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator*=( const Matrix2& M ) {
    if( ( M.get_col_count() != colCount ) || ( M.get_row_count() != colCount ) )
      throw std::range_error( "Matrix Dimension Mismatch." );
    *this = *this * M;
    return *this;
  };
};


template < typename Matrix >
struct is_readable_matrix< mat_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix< Matrix >::value );
  typedef is_readable_matrix< Matrix > type;
};

template < typename Matrix >
struct is_writable_matrix< mat_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix< Matrix >::value );
  typedef is_fully_writable_matrix< Matrix > type;
};

template < typename Matrix >
struct is_fully_writable_matrix< mat_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix< Matrix >::value );
  typedef is_fully_writable_matrix< Matrix > type;
};

template < typename Matrix >
struct is_resizable_matrix< mat_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_sub_block< Matrix > > type;
};

template < typename Matrix >
struct has_allocator_matrix< mat_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix< Matrix >::value );
  typedef has_allocator_matrix< Matrix > type;
};


template < typename Matrix >
class mat_const_sub_block {
public:
  typedef mat_const_sub_block< Matrix > self;
  typedef typename mat_traits< Matrix >::allocator_type allocator_type;

  typedef typename mat_traits< Matrix >::value_type value_type;

  typedef typename mat_traits< Matrix >::reference reference;
  typedef typename mat_traits< Matrix >::const_reference const_reference;
  typedef typename mat_traits< Matrix >::pointer pointer;
  typedef typename mat_traits< Matrix >::const_pointer const_pointer;

  typedef typename mat_traits< Matrix >::col_iterator col_iterator;
  typedef typename mat_traits< Matrix >::const_col_iterator const_col_iterator;
  typedef typename mat_traits< Matrix >::row_iterator row_iterator;
  typedef typename mat_traits< Matrix >::const_row_iterator const_row_iterator;

  typedef typename mat_traits< Matrix >::size_type size_type;
  typedef typename mat_traits< Matrix >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 0 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_traits< Matrix >::alignment );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_structure::rectangular );

private:
  const Matrix* m;
  size_type rowOffset;
  size_type colOffset;
  size_type rowCount;
  size_type colCount;

  self& operator=( const self& );
  explicit mat_const_sub_block( Matrix&& );
  mat_const_sub_block( Matrix&&, size_type, size_type, size_type aRowOffset = 0, size_type aColOffset = 0 );

public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_const_sub_block( const Matrix& aM )
      : m( &aM ), rowOffset( 0 ), colOffset( 0 ), rowCount( aM.get_row_count() ), colCount( aM.get_col_count() ){};

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aRowCount The number of rows for the sub-block.
   * \param aColCount The number of columns for the sub-block.
   * \param aRowOffset The row-offset from the start of the matrix.
   * \param aColOffset The column-offset from the start of the matrix.
   */
  mat_const_sub_block( const Matrix& aM, size_type aRowCount, size_type aColCount, size_type aRowOffset = 0,
                       size_type aColOffset = 0 )
      : m( &aM ), rowOffset( aRowOffset ), colOffset( aColOffset ), rowCount( aRowCount ), colCount( aColCount ){};

  /**
   * Standard copy-constructor.
   */
  mat_const_sub_block( const self& aObj )
      : m( aObj.m ), rowOffset( aObj.rowOffset ), colOffset( aObj.colOffset ), rowCount( aObj.rowCount ),
        colCount( aObj.colCount ){};

  /**
   * Standard move-constructor.
   */
  mat_const_sub_block( self&& aObj )
      : m( aObj.m ), rowOffset( aObj.rowOffset ), colOffset( aObj.colOffset ), rowCount( aObj.rowCount ),
        colCount( aObj.colCount ){};

  /**
   * Standard swap function (shallow).
   */
  friend void swap( self& lhs, self& rhs ) throw() {
    using std::swap;
    swap( lhs.m, rhs.m );
    swap( lhs.rowOffset, rhs.rowOffset );
    swap( lhs.colOffset, rhs.colOffset );
    swap( lhs.rowCount, rhs.rowCount );
    swap( lhs.colCount, rhs.colCount );
    return;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()( size_type i, size_type j ) const { return ( *m )( rowOffset + i, colOffset + j ); };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const throw() { return rowCount; };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const throw() { return colCount; };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() { return std::make_pair( rowCount, colCount ); };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return ( *m ).get_allocator(); };
};


template < typename Matrix >
struct is_readable_matrix< mat_const_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix< Matrix >::value );
  typedef is_readable_matrix< Matrix > type;
};

template < typename Matrix >
struct is_writable_matrix< mat_const_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< mat_const_sub_block< Matrix > > type;
};

template < typename Matrix >
struct is_fully_writable_matrix< mat_const_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_fully_writable_matrix< mat_const_sub_block< Matrix > > type;
};

template < typename Matrix >
struct is_resizable_matrix< mat_const_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_const_sub_block< Matrix > > type;
};

template < typename Matrix >
struct has_allocator_matrix< mat_const_sub_block< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix< Matrix >::value );
  typedef has_allocator_matrix< Matrix > type;
};


template < typename Matrix >
struct mat_copy_sub_block_factory {
  typedef typename mat_traits< Matrix >::size_type size_type;

  Matrix m;
  mat_copy_sub_block_factory( Matrix&& aM ) : m( std::move( aM ) ){};
  mat_copy_sub_block< Matrix > operator()( const std::pair< size_type, size_type >& rows,
                                           const std::pair< size_type, size_type >& cols ) {
    return mat_copy_sub_block< Matrix >( std::move( m ), rows.second - rows.first, cols.second - cols.first, rows.first,
                                         cols.first );
  };
};

template < typename Matrix >
struct mat_sub_block_factory {
  typedef typename mat_traits< Matrix >::size_type size_type;

  Matrix& m;
  mat_sub_block_factory( Matrix& aM ) : m( aM ){};
  mat_sub_block< Matrix > operator()( const std::pair< size_type, size_type >& rows,
                                      const std::pair< size_type, size_type >& cols ) {
    return mat_sub_block< Matrix >( m, rows.second - rows.first, cols.second - cols.first, rows.first, cols.first );
  };
};

template < typename Matrix >
struct mat_const_sub_block_factory {
  typedef typename mat_traits< Matrix >::size_type size_type;

  const Matrix& m;
  mat_const_sub_block_factory( const Matrix& aM ) : m( aM ){};
  mat_const_sub_block< Matrix > operator()( const std::pair< size_type, size_type >& rows,
                                            const std::pair< size_type, size_type >& cols ) {
    return mat_const_sub_block< Matrix >( m, rows.second - rows.first, cols.second - cols.first, rows.first,
                                          cols.first );
  };
};


template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, mat_sub_block_factory< Matrix > >::type sub( Matrix& M ) {
  return mat_sub_block_factory< Matrix >( M );
};

template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, mat_const_sub_block_factory< Matrix > >::type
  sub( const Matrix& M ) {
  return mat_const_sub_block_factory< Matrix >( M );
};

template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, mat_copy_sub_block_factory< Matrix > >::type
  sub_copy( const Matrix& M ) {
  return mat_copy_sub_block_factory< Matrix >( M );
};

template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, mat_copy_sub_block_factory< Matrix > >::type
  sub( Matrix&& M ) {
  return mat_copy_sub_block_factory< Matrix >( std::move( M ) );
};


template < typename Matrix, mat_structure::tag Structure = mat_traits< Matrix >::structure >
class mat_copy_sub_sym_block {
  char invalid_matrix_type[0];
};

template < typename Matrix, mat_structure::tag Structure = mat_traits< Matrix >::structure >
class mat_sub_sym_block {
  char invalid_matrix_type[0];
};

template < typename Matrix, mat_structure::tag Structure = mat_traits< Matrix >::structure >
class mat_const_sub_sym_block {
  char invalid_matrix_type[0];
};


template < typename Matrix >
class mat_copy_sub_sym_block< Matrix, mat_structure::symmetric > {
public:
  typedef mat_sub_sym_block< Matrix, mat_structure::symmetric > self;
  typedef typename mat_traits< Matrix >::allocator_type allocator_type;

  typedef typename mat_traits< Matrix >::value_type value_type;

  typedef typename mat_traits< Matrix >::reference reference;
  typedef typename mat_traits< Matrix >::const_reference const_reference;
  typedef typename mat_traits< Matrix >::pointer pointer;
  typedef typename mat_traits< Matrix >::const_pointer const_pointer;

  typedef typename mat_traits< Matrix >::col_iterator col_iterator;
  typedef typename mat_traits< Matrix >::const_col_iterator const_col_iterator;
  typedef typename mat_traits< Matrix >::row_iterator row_iterator;
  typedef typename mat_traits< Matrix >::const_row_iterator const_row_iterator;

  typedef typename mat_traits< Matrix >::size_type size_type;
  typedef typename mat_traits< Matrix >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 0 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_traits< Matrix >::alignment );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_structure::symmetric );

private:
  Matrix m;
  size_type offset;
  size_type rowCount;

public:
  /**
   * Default constructor.
   */
  mat_copy_sub_sym_block() : m(), offset( 0 ), rowCount( 0 ){};

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_copy_sub_sym_block( const Matrix& aM ) : m( aM ), offset( 0 ), rowCount( aM.get_row_count() ){};

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_copy_sub_sym_block( const Matrix& aM, size_type aSize, size_type aOffset = 0 )
      : m( aM ), offset( aOffset ), rowCount( aSize ){};

  /**
   * Standard copy-constructor.
   */
  mat_copy_sub_sym_block( const self& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_copy_sub_sym_block( Matrix&& aM ) : m( std::move( aM ) ), offset( 0 ), rowCount( 0 ) {
    rowCount = m.get_row_count();
  };

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_copy_sub_sym_block( Matrix&& aM, size_type aSize, size_type aOffset = 0 )
      : m( std::move( aM ) ), offset( aOffset ), rowCount( aSize ){};

  /**
   * Standard move-constructor.
   */
  mat_copy_sub_sym_block( self&& aObj ) : m( std::move( aObj.m ) ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard swap function.
   */
  friend void swap( self& lhs, self& rhs ) throw() {
    using std::swap;
    swap( lhs.m, rhs.m );
    swap( lhs.offset, rhs.offset );
    swap( lhs.rowCount, rhs.rowCount );
    return;
  };

  /**
   * Standard assignment operator.
   */
  self& operator=( self rhs ) {
    swap( *this, rhs );
    return *this;
  };


  /**
   * Standard assignment operator.
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator=( const Matrix2& rhs ) {
    if( ( rhs.get_row_count() != rowCount ) || ( rhs.get_col_count() != rowCount ) )
      throw std::range_error( "Matrix dimensions mismatch." );
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = j; i < rowCount; ++i )
        m( offset + i, offset + j ) = ( rhs( i, j ) + rhs( j, i ) ) * value_type( 0.5 );
    return *this;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-write access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  reference operator()( size_type i, size_type j ) { return m( offset + i, offset + j ); };
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()( size_type i, size_type j ) const { return m( offset + i, offset + j ); };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const throw() { return rowCount; };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const throw() { return rowCount; };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() { return std::make_pair( rowCount, rowCount ); };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return m.get_allocator(); };


  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator+=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = j; i < rowCount; ++i )
        m( offset + i, offset + j ) += ( M( i, j ) + M( j, i ) ) * value_type( 0.5 );
    return *this;
  };

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator-=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = j; i < rowCount; ++i )
        m( offset + i, offset + j ) -= ( M( i, j ) + M( j, i ) ) * value_type( 0.5 );
    return *this;
  };

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=( const value_type& S ) {
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = j; i < rowCount; ++i )
        m( offset + i, offset + j ) *= S;
    return *this;
  };

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator*=( const Matrix2& M ) {
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix Dimension Mismatch." );
    *this = *this * M;
    return *this;
  };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose( const self& M ) { return M; };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose_move( const self& M ) { return M; };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend self&& transpose( self&& M ) { return std::move( M ); };

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace( const self& M ) {
    value_type result( 0.0 );
    for( size_type i = 0; i < M.rowCount; ++i )
      result += M( i, i );
    return result;
  };
};


template < typename Matrix >
class mat_sub_sym_block< Matrix, mat_structure::symmetric > {
public:
  typedef mat_sub_sym_block< Matrix, mat_structure::symmetric > self;
  typedef typename mat_traits< Matrix >::allocator_type allocator_type;

  typedef typename mat_traits< Matrix >::value_type value_type;

  typedef typename mat_traits< Matrix >::reference reference;
  typedef typename mat_traits< Matrix >::const_reference const_reference;
  typedef typename mat_traits< Matrix >::pointer pointer;
  typedef typename mat_traits< Matrix >::const_pointer const_pointer;

  typedef typename mat_traits< Matrix >::col_iterator col_iterator;
  typedef typename mat_traits< Matrix >::const_col_iterator const_col_iterator;
  typedef typename mat_traits< Matrix >::row_iterator row_iterator;
  typedef typename mat_traits< Matrix >::const_row_iterator const_row_iterator;

  typedef typename mat_traits< Matrix >::size_type size_type;
  typedef typename mat_traits< Matrix >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 0 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_traits< Matrix >::alignment );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_structure::symmetric );

private:
  Matrix* m;
  size_type offset;
  size_type rowCount;

public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_sub_sym_block( Matrix& aM ) : m( &aM ), offset( 0 ), rowCount( aM.get_row_count() ){};

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_sub_sym_block( Matrix& aM, size_type aSize, size_type aOffset = 0 )
      : m( &aM ), offset( aOffset ), rowCount( aSize ){};
  /**
   * Standard copy-constructor.
   */
  mat_sub_sym_block( const self& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard move-constructor.
   */
  mat_sub_sym_block( self&& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard swap function (shallow).
   */
  friend void swap( self& lhs, self& rhs ) throw() {
    using std::swap;
    swap( lhs.m, rhs.m );
    swap( lhs.offset, rhs.offset );
    swap( lhs.rowCount, rhs.rowCount );
    return;
  };

  /**
   * Standard assignment operator.
   */
  self& operator=( const self& rhs ) {
    if( ( rhs.get_row_count() != rowCount ) || ( rhs.get_col_count() != rowCount ) )
      throw std::range_error( "Matrix dimensions mismatch." );
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = j; i < rowCount; ++i )
        ( *m )( offset + i, offset + j ) = ( rhs( i, j ) + rhs( j, i ) ) * value_type( 0.5 );
    return *this;
  };

  /**
   * Standard assignment operator.
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator=( const Matrix2& rhs ) {
    if( ( rhs.get_row_count() != rowCount ) || ( rhs.get_col_count() != rowCount ) )
      throw std::range_error( "Matrix dimensions mismatch." );
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = j; i < rowCount; ++i )
        ( *m )( offset + i, offset + j ) = ( rhs( i, j ) + rhs( j, i ) ) * value_type( 0.5 );
    return *this;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-write access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  reference operator()( size_type i, size_type j ) { return ( *m )( offset + i, offset + j ); };
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()( size_type i, size_type j ) const { return ( *m )( offset + i, offset + j ); };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const throw() { return rowCount; };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const throw() { return rowCount; };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() { return std::make_pair( rowCount, rowCount ); };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return m->get_allocator(); };


  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator+=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = j; i < rowCount; ++i )
        ( *m )( offset + i, offset + j ) += ( M( i, j ) + M( j, i ) ) * value_type( 0.5 );
    return *this;
  };

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator-=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = j; i < rowCount; ++i )
        ( *m )( offset + i, offset + j ) -= ( M( i, j ) + M( j, i ) ) * value_type( 0.5 );
    return *this;
  };

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=( const value_type& S ) {
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = j; i < rowCount; ++i )
        ( *m )( offset + i, offset + j ) *= S;
    return *this;
  };

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator*=( const Matrix2& M ) {
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix Dimension Mismatch." );
    *this = *this * M;
    return *this;
  };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose( const self& M ) { return M; };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose_move( const self& M ) { return M; };

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace( const self& M ) {
    value_type result( 0.0 );
    for( size_type i = 0; i < M.rowCount; ++i )
      result += M( i, i );
    return result;
  };
};


template < typename Matrix >
class mat_const_sub_sym_block< Matrix, mat_structure::symmetric > {
public:
  typedef mat_const_sub_sym_block< Matrix, mat_structure::symmetric > self;
  typedef typename mat_traits< Matrix >::allocator_type allocator_type;

  typedef typename mat_traits< Matrix >::value_type value_type;

  typedef typename mat_traits< Matrix >::reference reference;
  typedef typename mat_traits< Matrix >::const_reference const_reference;
  typedef typename mat_traits< Matrix >::pointer pointer;
  typedef typename mat_traits< Matrix >::const_pointer const_pointer;

  typedef typename mat_traits< Matrix >::col_iterator col_iterator;
  typedef typename mat_traits< Matrix >::const_col_iterator const_col_iterator;
  typedef typename mat_traits< Matrix >::row_iterator row_iterator;
  typedef typename mat_traits< Matrix >::const_row_iterator const_row_iterator;

  typedef typename mat_traits< Matrix >::size_type size_type;
  typedef typename mat_traits< Matrix >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 0 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_traits< Matrix >::alignment );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_structure::symmetric );

private:
  const Matrix* m;
  size_type offset;
  size_type rowCount;

  self& operator=( const self& );

  mat_const_sub_sym_block( Matrix&& );
  mat_const_sub_sym_block( Matrix&&, size_type, size_type aOffset = 0 );

public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  mat_const_sub_sym_block( const Matrix& aM ) : m( &aM ), offset( 0 ), rowCount( aM.get_row_count() ){};

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_const_sub_sym_block( const Matrix& aM, size_type aSize, size_type aOffset = 0 )
      : m( &aM ), offset( aOffset ), rowCount( aSize ){};

  /**
   * Standard copy-constructor.
   */
  mat_const_sub_sym_block( const self& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard move-constructor.
   */
  mat_const_sub_sym_block( self&& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard swap function (shallow).
   */
  friend void swap( self& lhs, self& rhs ) throw() {
    using std::swap;
    swap( lhs.m, rhs.m );
    swap( lhs.offset, rhs.offset );
    swap( lhs.rowCount, rhs.rowCount );
    return;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()( size_type i, size_type j ) const { return ( *m )( offset + i, offset + j ); };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const throw() { return rowCount; };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const throw() { return rowCount; };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() { return std::make_pair( rowCount, rowCount ); };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return m->get_allocator(); };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose( const self& M ) { return M; };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose_move( const self& M ) { return M; };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend self&& transpose( self&& M ) { return std::move( M ); };

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace( const self& M ) {
    value_type result( 0.0 );
    for( size_type i = 0; i < M.rowCount; ++i )
      result += M( i, i );
    return result;
  };
};


template < typename Matrix >
class mat_copy_sub_sym_block< Matrix, mat_structure::skew_symmetric > {
public:
  typedef mat_sub_sym_block< Matrix, mat_structure::skew_symmetric > self;
  typedef typename mat_traits< Matrix >::allocator_type allocator_type;

  typedef typename mat_traits< Matrix >::value_type value_type;

  typedef typename mat_traits< Matrix >::reference reference;
  typedef typename mat_traits< Matrix >::const_reference const_reference;
  typedef typename mat_traits< Matrix >::pointer pointer;
  typedef typename mat_traits< Matrix >::const_pointer const_pointer;

  typedef typename mat_traits< Matrix >::col_iterator col_iterator;
  typedef typename mat_traits< Matrix >::const_col_iterator const_col_iterator;
  typedef typename mat_traits< Matrix >::row_iterator row_iterator;
  typedef typename mat_traits< Matrix >::const_row_iterator const_row_iterator;

  typedef typename mat_traits< Matrix >::size_type size_type;
  typedef typename mat_traits< Matrix >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 0 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_traits< Matrix >::alignment );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_structure::skew_symmetric );

private:
  Matrix m;
  size_type offset;
  size_type rowCount;

public:
  /**
   * Default constructor.
   */
  mat_copy_sub_sym_block() : m(), offset( 0 ), rowCount( 0 ){};

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_copy_sub_sym_block( const Matrix& aM ) : m( aM ), offset( 0 ), rowCount( aM.get_row_count() ){};

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_copy_sub_sym_block( const Matrix& aM, size_type aSize, size_type aOffset = 0 )
      : m( aM ), offset( aOffset ), rowCount( aSize ){};

  /**
   * Standard copy-constructor.
   */
  mat_copy_sub_sym_block( const self& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_copy_sub_sym_block( Matrix&& aM ) : m( std::move( aM ) ), offset( 0 ), rowCount( 0 ) {
    rowCount = m.get_row_count();
  };

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_copy_sub_sym_block( Matrix&& aM, size_type aSize, size_type aOffset = 0 )
      : m( std::move( aM ) ), offset( aOffset ), rowCount( aSize ){};

  /**
   * Standard move-constructor.
   */
  mat_copy_sub_sym_block( self&& aObj ) : m( std::move( aObj.m ) ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard swap function.
   */
  friend void swap( self& lhs, self& rhs ) throw() {
    using std::swap;
    swap( lhs.m, rhs.m );
    swap( lhs.offset, rhs.offset );
    swap( lhs.rowCount, rhs.rowCount );
    return;
  };

  /**
   * Standard assignment operator.
   */
  self& operator=( self rhs ) {
    swap( *this, rhs );
    return *this;
  };


  /**
   * Standard assignment operator.
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator=( const Matrix2& rhs ) {
    if( ( rhs.get_row_count() != rowCount ) || ( rhs.get_col_count() != rowCount ) )
      throw std::range_error( "Matrix dimensions mismatch." );
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = j; i < rowCount; ++i )
        m( offset + i, offset + j ) = ( rhs( i, j ) + rhs( j, i ) ) * value_type( 0.5 );
    return *this;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-write access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  reference operator()( size_type i, size_type j ) { return m( offset + i, offset + j ); };
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()( size_type i, size_type j ) const { return m( offset + i, offset + j ); };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const throw() { return rowCount; };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const throw() { return rowCount; };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() { return std::make_pair( rowCount, rowCount ); };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return m.get_allocator(); };


  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator+=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = j; i < rowCount; ++i )
        m( offset + i, offset + j ) += ( M( i, j ) + M( j, i ) ) * value_type( 0.5 );
    return *this;
  };

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator-=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = j; i < rowCount; ++i )
        m( offset + i, offset + j ) -= ( M( i, j ) + M( j, i ) ) * value_type( 0.5 );
    return *this;
  };

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=( const value_type& S ) {
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = j; i < rowCount; ++i )
        m( offset + i, offset + j ) *= S;
    return *this;
  };

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator*=( const Matrix2& M ) {
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix Dimension Mismatch." );
    *this = *this * M;
    return *this;
  };

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace( const self& M ) {
    value_type result( 0.0 );
    for( size_type i = 0; i < M.rowCount; ++i )
      result += M( i, i );
    return result;
  };
};


template < typename Matrix >
class mat_sub_sym_block< Matrix, mat_structure::skew_symmetric > {
public:
  typedef mat_sub_sym_block< Matrix, mat_structure::skew_symmetric > self;
  typedef typename mat_traits< Matrix >::allocator_type allocator_type;

  typedef typename mat_traits< Matrix >::value_type value_type;

  typedef typename mat_traits< Matrix >::reference reference;
  typedef typename mat_traits< Matrix >::const_reference const_reference;
  typedef typename mat_traits< Matrix >::pointer pointer;
  typedef typename mat_traits< Matrix >::const_pointer const_pointer;

  typedef typename mat_traits< Matrix >::col_iterator col_iterator;
  typedef typename mat_traits< Matrix >::const_col_iterator const_col_iterator;
  typedef typename mat_traits< Matrix >::row_iterator row_iterator;
  typedef typename mat_traits< Matrix >::const_row_iterator const_row_iterator;

  typedef typename mat_traits< Matrix >::size_type size_type;
  typedef typename mat_traits< Matrix >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 0 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_traits< Matrix >::alignment );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_structure::skew_symmetric );

private:
  Matrix* m;
  size_type offset;
  size_type rowCount;

public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  mat_sub_sym_block( Matrix& aM ) : m( &aM ), offset( 0 ), rowCount( aM.get_row_count() ){};

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_sub_sym_block( Matrix& aM, size_type aSize, size_type aOffset = 0 )
      : m( &aM ), offset( aOffset ), rowCount( aSize ){};

  /**
   * Standard copy-constructor.
   */
  mat_sub_sym_block( const self& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard move-constructor.
   */
  mat_sub_sym_block( self&& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard swap function (shallow).
   */
  friend void swap( self& lhs, self& rhs ) throw() {
    using std::swap;
    swap( lhs.m, rhs.m );
    swap( lhs.offset, rhs.offset );
    swap( lhs.rowCount, rhs.rowCount );
    return;
  };

  /**
   * Standard assignment operator.
   */
  self& operator=( const self& rhs ) {
    if( ( rhs.get_row_count() != rowCount ) || ( rhs.get_col_count() != rowCount ) )
      throw std::range_error( "Matrix dimensions mismatch." );
    for( size_type j = 1; j < rowCount; ++j )
      for( size_type i = 0; i < j; ++i )
        ( *m )( offset + i, offset + j ) = ( rhs( i, j ) - rhs( j, i ) ) * value_type( 0.5 );
    return *this;
  };

  /**
   * Standard assignment operator.
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator=( const Matrix2& rhs ) {
    if( ( rhs.get_row_count() != rowCount ) || ( rhs.get_col_count() != rowCount ) )
      throw std::range_error( "Matrix dimensions mismatch." );
    for( size_type j = 1; j < rowCount; ++j )
      for( size_type i = 0; i < j; ++i )
        ( *m )( offset + i, offset + j ) = ( rhs( i, j ) - rhs( j, i ) ) * value_type( 0.5 );
    return *this;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-write access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  reference operator()( size_type i, size_type j ) { return ( *m )( offset + i, offset + j ); };
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()( size_type i, size_type j ) const { return ( *m )( offset + i, offset + j ); };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const throw() { return rowCount; };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const throw() { return rowCount; };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() { return std::make_pair( rowCount, rowCount ); };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return m->get_allocator(); };


  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator+=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = 0; i < j; ++i )
        ( *m )( offset + i, offset + j ) += ( M( i, j ) - M( j, i ) ) * value_type( 0.5 );
    return *this;
  };

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator-=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = 0; i < j; ++i )
        ( *m )( offset + i, offset + j ) -= ( M( i, j ) - M( j, i ) ) * value_type( 0.5 );
    return *this;
  };

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=( const value_type& S ) {
    for( size_type j = 0; j < rowCount; ++j )
      for( size_type i = 0; i < j; ++i )
        ( *m )( offset + i, offset + j ) *= S;
    return *this;
  };

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator*=( const Matrix2& M ) {
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix Dimension Mismatch." );
    *this = *this * M;
    return *this;
  };

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace( const self& M ) { return value_type( 0.0 ); };
};


template < typename Matrix >
class mat_const_sub_sym_block< Matrix, mat_structure::skew_symmetric > {
public:
  typedef mat_const_sub_sym_block< Matrix, mat_structure::skew_symmetric > self;
  typedef typename mat_traits< Matrix >::allocator_type allocator_type;

  typedef typename mat_traits< Matrix >::value_type value_type;

  typedef typename mat_traits< Matrix >::reference reference;
  typedef typename mat_traits< Matrix >::const_reference const_reference;
  typedef typename mat_traits< Matrix >::pointer pointer;
  typedef typename mat_traits< Matrix >::const_pointer const_pointer;

  typedef typename mat_traits< Matrix >::col_iterator col_iterator;
  typedef typename mat_traits< Matrix >::const_col_iterator const_col_iterator;
  typedef typename mat_traits< Matrix >::row_iterator row_iterator;
  typedef typename mat_traits< Matrix >::const_row_iterator const_row_iterator;

  typedef typename mat_traits< Matrix >::size_type size_type;
  typedef typename mat_traits< Matrix >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 0 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_traits< Matrix >::alignment );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_structure::skew_symmetric );

private:
  const Matrix* m;
  size_type offset;
  size_type rowCount;

  self& operator=( const self& );

  mat_const_sub_sym_block( Matrix&& aM );
  mat_const_sub_sym_block( Matrix&& aM, size_type aSize, size_type aOffset = 0 );

public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  mat_const_sub_sym_block( const Matrix& aM ) : m( &aM ), offset( 0 ), rowCount( aM.get_row_count() ){};

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_const_sub_sym_block( const Matrix& aM, size_type aSize, size_type aOffset = 0 )
      : m( &aM ), offset( aOffset ), rowCount( aSize ){};

  /**
   * Standard copy-constructor.
   */
  mat_const_sub_sym_block( const self& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard move-constructor.
   */
  mat_const_sub_sym_block( self&& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard swap function (shallow).
   */
  friend void swap( self& lhs, self& rhs ) throw() {
    using std::swap;
    swap( lhs.m, rhs.m );
    swap( lhs.offset, rhs.offset );
    swap( lhs.rowCount, rhs.rowCount );
    return;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()( size_type i, size_type j ) const { return ( *m )( offset + i, offset + j ); };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const throw() { return rowCount; };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const throw() { return rowCount; };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() { return std::make_pair( rowCount, rowCount ); };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return m->get_allocator(); };

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace( const self& M ) { return value_type( 0.0 ); };
};


template < typename Matrix >
class mat_copy_sub_sym_block< Matrix, mat_structure::diagonal > {
public:
  typedef mat_copy_sub_sym_block< Matrix, mat_structure::diagonal > self;
  typedef typename mat_traits< Matrix >::allocator_type allocator_type;

  typedef typename mat_traits< Matrix >::value_type value_type;

  typedef typename mat_traits< Matrix >::reference reference;
  typedef typename mat_traits< Matrix >::const_reference const_reference;
  typedef typename mat_traits< Matrix >::pointer pointer;
  typedef typename mat_traits< Matrix >::const_pointer const_pointer;

  typedef typename mat_traits< Matrix >::col_iterator col_iterator;
  typedef typename mat_traits< Matrix >::const_col_iterator const_col_iterator;
  typedef typename mat_traits< Matrix >::row_iterator row_iterator;
  typedef typename mat_traits< Matrix >::const_row_iterator const_row_iterator;

  typedef typename mat_traits< Matrix >::size_type size_type;
  typedef typename mat_traits< Matrix >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 0 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_traits< Matrix >::alignment );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_structure::diagonal );

private:
  Matrix m;
  size_type offset;
  size_type rowCount;

public:
  /**
   * Default constructor.
   */
  mat_copy_sub_sym_block() : m(), offset( 0 ), rowCount( 0 ) { rowCount = m.get_row_count(); };

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  mat_copy_sub_sym_block( const Matrix& aM ) : m( aM ), offset( 0 ), rowCount( aM.get_row_count() ){};

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_copy_sub_sym_block( const Matrix& aM, size_type aSize, size_type aOffset = 0 )
      : m( aM ), offset( aOffset ), rowCount( aSize ){};

  /**
   * Standard copy-constructor.
   */
  mat_copy_sub_sym_block( const self& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  mat_copy_sub_sym_block( Matrix&& aM ) : m( std::move( aM ) ), offset( 0 ), rowCount( 0 ) {
    rowCount = m.get_row_count();
  };

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_copy_sub_sym_block( Matrix&& aM, size_type aSize, size_type aOffset = 0 )
      : m( std::move( aM ) ), offset( aOffset ), rowCount( aSize ){};

  /**
   * Standard move-constructor.
   */
  mat_copy_sub_sym_block( self&& aObj ) : m( std::move( aObj.m ) ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard swap function.
   */
  friend void swap( self& lhs, self& rhs ) throw() {
    using std::swap;
    swap( lhs.m, rhs.m );
    swap( lhs.offset, rhs.offset );
    swap( lhs.rowCount, rhs.rowCount );
    return;
  };

  /**
   * Standard assignment operator.
   */
  self& operator=( self rhs ) {
    swap( *this, rhs );
    return *this;
  };

  /**
   * Standard assignment operator.
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator=( const Matrix2& rhs ) {
    if( ( rhs.get_row_count() != rowCount ) || ( rhs.get_col_count() != rowCount ) )
      throw std::range_error( "Matrix dimensions mismatch." );
    for( size_type i = 0; i < rowCount; ++i )
      m( offset + i, offset + i ) = rhs( i, i );
    return *this;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-write access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  reference operator()( size_type i, size_type j ) { return m( offset + i, offset + j ); };
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()( size_type i, size_type j ) const { return m( offset + i, offset + j ); };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const throw() { return rowCount; };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const throw() { return rowCount; };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() { return std::make_pair( rowCount, rowCount ); };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return m.get_allocator(); };


  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator+=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type i = 0; i < rowCount; ++i )
      m( offset + i, offset + i ) += M( i, i );
    return *this;
  };

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator-=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type i = 0; i < rowCount; ++i )
      m( offset + i, offset + i ) -= M( i, i );
    return *this;
  };

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=( const value_type& S ) {
    for( size_type i = 0; i < rowCount; ++i )
      m( offset + i, offset + i ) *= S;
    return *this;
  };

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator*=( const Matrix2& M ) {
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix Dimension Mismatch." );
    *this = *this * M;
    return *this;
  };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose( const self& M ) { return M; };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose_move( const self& M ) { return M; };

  friend self&& transpose( self&& M ) { return std::move( M ); };

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace( const self& M ) {
    value_type result( 0.0 );
    for( size_type i = 0; i < M.rowCount; ++i )
      result += M( i, i );
    return result;
  };
};


template < typename Matrix >
class mat_sub_sym_block< Matrix, mat_structure::diagonal > {
public:
  typedef mat_sub_sym_block< Matrix, mat_structure::diagonal > self;
  typedef typename mat_traits< Matrix >::allocator_type allocator_type;

  typedef typename mat_traits< Matrix >::value_type value_type;

  typedef typename mat_traits< Matrix >::reference reference;
  typedef typename mat_traits< Matrix >::const_reference const_reference;
  typedef typename mat_traits< Matrix >::pointer pointer;
  typedef typename mat_traits< Matrix >::const_pointer const_pointer;

  typedef typename mat_traits< Matrix >::col_iterator col_iterator;
  typedef typename mat_traits< Matrix >::const_col_iterator const_col_iterator;
  typedef typename mat_traits< Matrix >::row_iterator row_iterator;
  typedef typename mat_traits< Matrix >::const_row_iterator const_row_iterator;

  typedef typename mat_traits< Matrix >::size_type size_type;
  typedef typename mat_traits< Matrix >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 0 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_traits< Matrix >::alignment );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_structure::diagonal );

private:
  Matrix* m;
  size_type offset;
  size_type rowCount;

public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  mat_sub_sym_block( Matrix& aM ) : m( &aM ), offset( 0 ), rowCount( aM.get_row_count() ){};

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_sub_sym_block( Matrix& aM, size_type aSize, size_type aOffset = 0 )
      : m( &aM ), offset( aOffset ), rowCount( aSize ){};
  /**
   * Standard copy-constructor.
   */
  mat_sub_sym_block( const self& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard move-constructor.
   */
  mat_sub_sym_block( self&& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard swap function (shallow).
   */
  friend void swap( self& lhs, self& rhs ) throw() {
    using std::swap;
    swap( lhs.m, rhs.m );
    swap( lhs.offset, rhs.offset );
    swap( lhs.rowCount, rhs.rowCount );
    return;
  };

  /**
   * Standard assignment operator.
   */
  self& operator=( const self& rhs ) {
    if( ( rhs.get_row_count() != rowCount ) || ( rhs.get_col_count() != rowCount ) )
      throw std::range_error( "Matrix dimensions mismatch." );
    for( size_type i = 0; i < rowCount; ++i )
      ( *m )( offset + i, offset + i ) = rhs( i, i );
    return *this;
  };

  /**
   * Standard assignment operator.
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator=( const Matrix2& rhs ) {
    if( ( rhs.get_row_count() != rowCount ) || ( rhs.get_col_count() != rowCount ) )
      throw std::range_error( "Matrix dimensions mismatch." );
    for( size_type i = 0; i < rowCount; ++i )
      ( *m )( offset + i, offset + i ) = rhs( i, i );
    return *this;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-write access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  reference operator()( size_type i, size_type j ) { return ( *m )( offset + i, offset + j ); };
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()( size_type i, size_type j ) const { return ( *m )( offset + i, offset + j ); };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const throw() { return rowCount; };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const throw() { return rowCount; };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() { return std::make_pair( rowCount, rowCount ); };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return m->get_allocator(); };


  /** COL-MAJOR ONLY
   * Add-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator+=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type i = 0; i < rowCount; ++i )
      ( *m )( offset + i, offset + i ) += M( i, i );
    return *this;
  };

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator-=( const Matrix2& M ) {
    BOOST_CONCEPT_ASSERT( (ReadableMatrixConcept< Matrix2 >));
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix dimension mismatch." );
    for( size_type i = 0; i < rowCount; ++i )
      ( *m )( offset + i, offset + i ) -= M( i, i );
    return *this;
  };

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=( const value_type& S ) {
    for( size_type i = 0; i < rowCount; ++i )
      ( *m )( offset + i, offset + i ) *= S;
    return *this;
  };

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator*=( const Matrix2& M ) {
    if( ( M.get_col_count() != rowCount ) || ( M.get_row_count() != rowCount ) )
      throw std::range_error( "Matrix Dimension Mismatch." );
    *this = *this * M;
    return *this;
  };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose( const self& M ) { return M; };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose_move( const self& M ) { return M; };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend self&& transpose( self&& M ) { return std::move( M ); };

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace( const self& M ) {
    value_type result( 0.0 );
    for( size_type i = 0; i < M.rowCount; ++i )
      result += M( i, i );
    return result;
  };
};


template < typename Matrix >
class mat_const_sub_sym_block< Matrix, mat_structure::diagonal > {
public:
  typedef mat_const_sub_sym_block< Matrix, mat_structure::diagonal > self;
  typedef typename mat_traits< Matrix >::allocator_type allocator_type;

  typedef typename mat_traits< Matrix >::value_type value_type;

  typedef typename mat_traits< Matrix >::reference reference;
  typedef typename mat_traits< Matrix >::const_reference const_reference;
  typedef typename mat_traits< Matrix >::pointer pointer;
  typedef typename mat_traits< Matrix >::const_pointer const_pointer;

  typedef typename mat_traits< Matrix >::col_iterator col_iterator;
  typedef typename mat_traits< Matrix >::const_col_iterator const_col_iterator;
  typedef typename mat_traits< Matrix >::row_iterator row_iterator;
  typedef typename mat_traits< Matrix >::const_row_iterator const_row_iterator;

  typedef typename mat_traits< Matrix >::size_type size_type;
  typedef typename mat_traits< Matrix >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 0 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_traits< Matrix >::alignment );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_structure::diagonal );

private:
  const Matrix* m;
  size_type offset;
  size_type rowCount;

  self& operator=( const self& );

  mat_const_sub_sym_block( Matrix&& );
  mat_const_sub_sym_block( Matrix&&, size_type, size_type aOffset = 0 );

public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  mat_const_sub_sym_block( const Matrix& aM ) : m( &aM ), offset( 0 ), rowCount( aM.get_row_count() ){};

  /**
   * Constructs the sub-matrix which represents part of the matrix.
   * \param aM The matrix from which the sub-block is taken.
   * \param aSize The number of rows for the sub-block.
   * \param aOffset The row-offset from the start of the matrix.
   */
  mat_const_sub_sym_block( const Matrix& aM, size_type aSize, size_type aOffset = 0 )
      : m( &aM ), offset( aOffset ), rowCount( aSize ){};

  /**
   * Standard copy-constructor.
   */
  mat_const_sub_sym_block( const self& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard move-constructor.
   */
  mat_const_sub_sym_block( self&& aObj ) : m( aObj.m ), offset( aObj.offset ), rowCount( aObj.rowCount ){};

  /**
   * Standard swap function (shallow).
   */
  friend void swap( self& lhs, self& rhs ) throw() {
    using std::swap;
    swap( lhs.m, rhs.m );
    swap( lhs.offset, rhs.offset );
    swap( lhs.rowCount, rhs.rowCount );
    return;
  };

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()( size_type i, size_type j ) const { return ( *m )( offset + i, offset + j ); };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const throw() { return rowCount; };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const throw() { return rowCount; };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() { return std::make_pair( rowCount, rowCount ); };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return m->get_allocator(); };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose( const self& M ) { return M; };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend const self& transpose_move( const self& M ) { return M; };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend self&& transpose( self&& M ) { return std::move( M ); };

  /**
   * Returns the trace of matrix M.
   * \param M A matrix.
   * \return the trace of matrix M.
   */
  friend value_type trace( const self& M ) {
    value_type result( 0.0 );
    for( size_type i = 0; i < M.rowCount; ++i )
      result += M( i, i );
    return result;
  };
};


template < typename Matrix, mat_structure::tag Structure >
struct is_readable_matrix< mat_copy_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix< Matrix >::value );
  typedef is_readable_matrix< Matrix > type;
};

template < typename Matrix, mat_structure::tag Structure >
struct is_writable_matrix< mat_copy_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = is_writable_matrix< Matrix >::value );
  typedef is_writable_matrix< Matrix > type;
};

template < typename Matrix, mat_structure::tag Structure >
struct is_fully_writable_matrix< mat_copy_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_fully_writable_matrix< mat_copy_sub_sym_block< Matrix, Structure > > type;
};

template < typename Matrix, mat_structure::tag Structure >
struct is_resizable_matrix< mat_copy_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_copy_sub_sym_block< Matrix, Structure > > type;
};

template < typename Matrix, mat_structure::tag Structure >
struct has_allocator_matrix< mat_copy_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix< Matrix >::value );
  typedef has_allocator_matrix< Matrix > type;
};


template < typename Matrix, mat_structure::tag Structure >
struct is_readable_matrix< mat_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix< Matrix >::value );
  typedef is_readable_matrix< Matrix > type;
};

template < typename Matrix, mat_structure::tag Structure >
struct is_writable_matrix< mat_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = is_writable_matrix< Matrix >::value );
  typedef is_writable_matrix< Matrix > type;
};

template < typename Matrix, mat_structure::tag Structure >
struct is_fully_writable_matrix< mat_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_fully_writable_matrix< mat_sub_sym_block< Matrix, Structure > > type;
};

template < typename Matrix, mat_structure::tag Structure >
struct is_resizable_matrix< mat_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_sub_sym_block< Matrix, Structure > > type;
};

template < typename Matrix, mat_structure::tag Structure >
struct has_allocator_matrix< mat_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix< Matrix >::value );
  typedef has_allocator_matrix< Matrix > type;
};


template < typename Matrix, mat_structure::tag Structure >
struct is_readable_matrix< mat_const_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix< Matrix >::value );
  typedef is_readable_matrix< Matrix > type;
};

template < typename Matrix, mat_structure::tag Structure >
struct is_writable_matrix< mat_const_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< mat_const_sub_sym_block< Matrix, Structure > > type;
};

template < typename Matrix, mat_structure::tag Structure >
struct is_fully_writable_matrix< mat_const_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_fully_writable_matrix< mat_const_sub_sym_block< Matrix, Structure > > type;
};

template < typename Matrix, mat_structure::tag Structure >
struct is_resizable_matrix< mat_const_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_const_sub_sym_block< Matrix, Structure > > type;
};

template < typename Matrix, mat_structure::tag Structure >
struct has_allocator_matrix< mat_const_sub_sym_block< Matrix, Structure > > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix< Matrix >::value );
  typedef has_allocator_matrix< Matrix > type;
};


template < typename Matrix >
struct mat_copy_sub_sym_block_factory {
  typedef typename mat_traits< Matrix >::size_type size_type;

  Matrix m;
  mat_copy_sub_sym_block_factory( Matrix&& aM ) : m( std::move( aM ) ){};
  mat_copy_sub_sym_block< Matrix > operator()( const std::pair< size_type, size_type >& rows ) {
    return mat_copy_sub_sym_block< Matrix >( std::move( m ), rows.second - rows.first, rows.first );
  };
};

template < typename Matrix >
struct mat_sub_sym_block_factory {
  typedef typename mat_traits< Matrix >::size_type size_type;

  Matrix& m;
  mat_sub_sym_block_factory( Matrix& aM ) : m( aM ){};
  mat_sub_sym_block< Matrix > operator()( const std::pair< size_type, size_type >& rows ) {
    return mat_sub_sym_block< Matrix >( m, rows.second - rows.first, rows.first );
  };
};

template < typename Matrix >
struct mat_const_sub_sym_block_factory {
  typedef typename mat_traits< Matrix >::size_type size_type;

  const Matrix& m;
  mat_const_sub_sym_block_factory( const Matrix& aM ) : m( aM ){};
  mat_const_sub_sym_block< Matrix > operator()( const std::pair< size_type, size_type >& rows ) {
    return mat_const_sub_sym_block< Matrix >( m, rows.second - rows.first, rows.first );
  };
};


template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, mat_sub_sym_block_factory< Matrix > >::type
  sub_sym( Matrix& M ) {
  return mat_sub_sym_block_factory< Matrix >( M );
};

template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, mat_const_sub_sym_block_factory< Matrix > >::type
  sub_sym( const Matrix& M ) {
  return mat_const_sub_sym_block_factory< Matrix >( M );
};

template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, mat_copy_sub_sym_block_factory< Matrix > >::type
  sub_sym_copy( const Matrix& M ) {
  return mat_copy_sub_sym_block_factory< Matrix >( M );
};

template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, mat_copy_sub_sym_block_factory< Matrix > >::type
  sub_sym( Matrix&& M ) {
  return mat_copy_sub_sym_block_factory< Matrix >( std::move( M ) );
};
};

#endif

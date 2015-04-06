/**
 * \file mat_transpose_view.hpp
 *
 * This library provides class templates to create transposed matrix views.
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

#ifndef REAK_MAT_TRANSPOSE_VIEW_HPP
#define REAK_MAT_TRANSPOSE_VIEW_HPP

#include "mat_alg_general.hpp"
#include <boost/static_assert.hpp>

namespace ReaK {


template < typename Matrix >
class mat_const_transpose_view {
public:
  typedef mat_const_transpose_view< Matrix > self;
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
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_traits< Matrix >::structure );

private:
  const Matrix* m;

  self& operator=( const self& );
  explicit mat_const_transpose_view( Matrix&& );

public:
  /**
   * Constructs the sub-matrix which represents the entire matrix.
   */
  explicit mat_const_transpose_view( const Matrix& aM ) : m( &aM ){};

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
  value_type operator()( size_type i, size_type j ) const { return ( *m )( j, i ); };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const throw() { return m->get_col_count(); };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const throw() { return m->get_row_count(); };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() {
    return std::make_pair( m->get_col_count(), m->get_row_count() );
  };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return ( *m ).get_allocator(); };

  /**
   * General negation operator for any type of matrices. This is a default operator
   * that will be called if no better special-purpose overload exists.
   * \return General column-major matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::rectangular > operator-() const {
    mat< value_type, mat_structure::rectangular > result( *this );
    for( size_type j = 0; j < result.get_col_count(); ++j )
      for( size_type i = 0; i < result.get_row_count(); ++i )
        result( i, j ) = -result( i, j );
    return result;
  };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend mat< value_type, mat_structure::rectangular > transpose( const self& M ) {
    return mat< value_type, mat_structure::rectangular >( *( M.m ) );
  };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend mat< value_type, mat_structure::rectangular > transpose_move( const self& M ) {
    return mat< value_type, mat_structure::rectangular >( *( M.m ) );
  };
};


template < typename Matrix >
struct is_readable_matrix< mat_const_transpose_view< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix< Matrix >::value );
  typedef is_readable_matrix< Matrix > type;
};

template < typename Matrix >
struct is_writable_matrix< mat_const_transpose_view< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< mat_const_transpose_view< Matrix > > type;
};

template < typename Matrix >
struct is_fully_writable_matrix< mat_const_transpose_view< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_fully_writable_matrix< mat_const_transpose_view< Matrix > > type;
};

template < typename Matrix >
struct is_resizable_matrix< mat_const_transpose_view< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_const_transpose_view< Matrix > > type;
};

template < typename Matrix >
struct has_allocator_matrix< mat_const_transpose_view< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix< Matrix >::value );
  typedef has_allocator_matrix< Matrix > type;
};

template < typename Matrix >
struct is_square_matrix< mat_const_transpose_view< Matrix > > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = is_square_matrix< Matrix >::value );
  typedef is_square_matrix< Matrix > type;
};

template < typename Matrix >
struct is_symmetric_matrix< mat_const_transpose_view< Matrix > > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = is_symmetric_matrix< Matrix >::value );
  typedef is_symmetric_matrix< Matrix > type;
};

template < typename Matrix >
struct is_diagonal_matrix< mat_const_transpose_view< Matrix > > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = is_diagonal_matrix< Matrix >::value );
  typedef is_diagonal_matrix< Matrix > type;
};


template < typename Matrix >
class mat_transpose_view {
public:
  typedef mat_transpose_view< Matrix > self;
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
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_traits< Matrix >::structure );

private:
  Matrix* m;

public:
  /**
   * Constructs a tranposed view of a matrix.
   */
  explicit mat_transpose_view( Matrix& aM ) : m( &aM ){};

  /**
   * Standard assignment operator.
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator=( const Matrix2& rhs ) {
    *m = mat_const_transpose_view< Matrix2 >( rhs );
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
  reference operator()( size_type i, size_type j ) { return ( *m )( j, i ); };
  /**
   * Matrix indexing accessor for read-only access.
   * \param i Row index.
   * \param j Column index.
   * \return the element at the given position.
   * \test PASSED
   */
  value_type operator()( size_type i, size_type j ) const { return ( *m )( j, i ); };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * TEST PASSED
   */
  size_type get_row_count() const throw() { return m->get_col_count(); };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * TEST PASSED
   */
  size_type get_col_count() const throw() { return m->get_row_count(); };

  /**
   * Sets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * TEST PASSED
   */
  void set_row_count( size_type aRowCount ) { m->set_col_count( aRowCount ); };

  /**
   * Sets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * TEST PASSED
   */
  void get_col_count( size_type aColCount ) { m->set_row_count( aColCount ); };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() {
    return std::make_pair( m->get_col_count(), m->get_row_count() );
  };

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
    ( *m ) += mat_const_transpose_view< Matrix2 >( M );
    return *this;
  };

  /** COL-MAJOR ONLY
   * Sub-and-store operator with standard semantics.
   * \test PASSED
   */
  template < typename Matrix2 >
  self& operator-=( const Matrix2& M ) {
    ( *m ) -= mat_const_transpose_view< Matrix2 >( M );
    return *this;
  };

  /** WORKS FOR ALL
   * Scalar-multiply-and-store operator with standard semantics.
   * \test PASSED
   */
  self& operator*=( const value_type& S ) {
    ( *m ) *= S;
    return *this;
  };

  /** WORKS FOR ALL
   * General Matrix multiplication.
   * \test PASSED
   */
  template < typename Matrix2 >
  typename boost::enable_if< is_readable_matrix< Matrix2 >, self& >::type operator*=( const Matrix2& M ) {
    *this = *this * M;
    return *this;
  };

  /** WORKS FOR ALL
   * General negation operator for any type of matrices. This is a default operator
   * that will be called if no better special-purpose overload exists.
   * \return General column-major matrix.
   * \test PASSED
   */
  mat< value_type, mat_structure::rectangular > operator-() const {
    mat< value_type, mat_structure::rectangular > result( *this );
    for( size_type j = 0; j < result.get_col_count(); ++j )
      for( size_type i = 0; i < result.get_row_count(); ++i )
        result( i, j ) = -result( i, j );
    return result;
  };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend mat< value_type, mat_structure::rectangular > transpose( const self& M ) {
    return mat< value_type, mat_structure::rectangular >( *( M.m ) );
  };

  /**
   * Transposes the matrix M.
   * \param M The matrix to be transposed.
   * \return The transpose of M.
   */
  friend mat< value_type, mat_structure::rectangular > transpose_move( const self& M ) {
    return mat< value_type, mat_structure::rectangular >( *( M.m ) );
  };
};


template < typename Matrix >
struct is_readable_matrix< mat_transpose_view< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = is_readable_matrix< Matrix >::value );
  typedef is_readable_matrix< Matrix > type;
};

template < typename Matrix >
struct is_writable_matrix< mat_transpose_view< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix< Matrix >::value );
  typedef is_fully_writable_matrix< Matrix > type;
};

template < typename Matrix >
struct is_fully_writable_matrix< mat_transpose_view< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = is_fully_writable_matrix< Matrix >::value );
  typedef is_fully_writable_matrix< Matrix > type;
};

template < typename Matrix >
struct is_resizable_matrix< mat_transpose_view< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = is_resizable_matrix< Matrix >::value );
  typedef is_resizable_matrix< Matrix > type;
};

template < typename Matrix >
struct has_allocator_matrix< mat_transpose_view< Matrix > > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix< Matrix >::value );
  typedef has_allocator_matrix< Matrix > type;
};


template < typename Matrix >
struct is_square_matrix< mat_transpose_view< Matrix > > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = is_square_matrix< Matrix >::value );
  typedef is_square_matrix< Matrix > type;
};

template < typename Matrix >
struct is_symmetric_matrix< mat_transpose_view< Matrix > > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = is_symmetric_matrix< Matrix >::value );
  typedef is_symmetric_matrix< Matrix > type;
};

template < typename Matrix >
struct is_diagonal_matrix< mat_transpose_view< Matrix > > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = is_diagonal_matrix< Matrix >::value );
  typedef is_diagonal_matrix< Matrix > type;
};


template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, mat_transpose_view< Matrix > >::type
  transpose_view( Matrix& M ) {
  return mat_transpose_view< Matrix >( M );
};

template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, mat_const_transpose_view< Matrix > >::type
  transpose_view( const Matrix& M ) {
  return mat_const_transpose_view< Matrix >( M );
};

template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, mat_const_transpose_view< Matrix > >::type
  transpose_view( Matrix&& M ) {
  return mat_const_transpose_view< Matrix >( std::move( M ) );
};
};


#endif

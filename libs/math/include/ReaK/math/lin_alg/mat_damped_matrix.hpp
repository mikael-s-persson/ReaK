/**
 * \file mat_damped_matrix.hpp
 *
 * This library provides an adaptor class that represents the addition of a diagonal matrix and
 * square matrix, i.e. a so-called damped matrix.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
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

#ifndef REAK_MAT_DAMPED_MATRIX_HPP
#define REAK_MAT_DAMPED_MATRIX_HPP


#include "mat_alg_general.hpp"

#include <boost/static_assert.hpp>

namespace ReaK {

/**
 * This class template forms the addition of two matrices, one diagonal and one square.
 *
 * Models: ReadableMatrixConcept.
 *
 * \tparam SquareMatrix Matrix type for the left matrix.
 * \tparam DiagMatrix Matrix type for the right matrix.
 */
template < typename SquareMatrix, typename DiagMatrix >
class mat_damped_matrix {
public:
  typedef mat_damped_matrix< SquareMatrix, DiagMatrix > self;
  typedef typename mat_traits< SquareMatrix >::allocator_type allocator_type;

  typedef typename mat_traits< SquareMatrix >::value_type value_type;

  typedef typename mat_traits< SquareMatrix >::reference reference;
  typedef typename mat_traits< SquareMatrix >::const_reference const_reference;
  typedef typename mat_traits< SquareMatrix >::pointer pointer;
  typedef typename mat_traits< SquareMatrix >::const_pointer const_pointer;

  typedef typename mat_traits< SquareMatrix >::col_iterator col_iterator;
  typedef typename mat_traits< SquareMatrix >::const_col_iterator const_col_iterator;
  typedef typename mat_traits< SquareMatrix >::row_iterator row_iterator;
  typedef typename mat_traits< SquareMatrix >::const_row_iterator const_row_iterator;

  typedef typename mat_traits< SquareMatrix >::size_type size_type;
  typedef typename mat_traits< SquareMatrix >::difference_type difference_type;

  BOOST_STATIC_CONSTANT( std::size_t, static_row_count = 0 );
  BOOST_STATIC_CONSTANT( std::size_t, static_col_count = 0 );
  BOOST_STATIC_CONSTANT( mat_alignment::tag, alignment = mat_traits< SquareMatrix >::alignment );
  BOOST_STATIC_CONSTANT( mat_structure::tag, structure = mat_traits< SquareMatrix >::structure );

private:
  const SquareMatrix* m_sqr; ///< Holds the left part of the matrix.
  const DiagMatrix* m_diag;  ///< Holds the right part of the matrix.
public:
  /**
   * Parametrized constructor.
   * \param aML Matrix to fill the left part of the matrix.
   * \param aMR Matrix to fill the right part of the matrix.
   */
  mat_damped_matrix( const SquareMatrix& aMSqr, const DiagMatrix& aMDiag ) : m_sqr( &aMSqr ), m_diag( &aMDiag ) {
    if( m_sqr->get_row_count() != m_diag->get_row_count() )
      throw std::range_error( "Matrix dimensions mismatch." );
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
  value_type operator()( size_type i, size_type j ) const {
    if( i == j )
      return ( *m_sqr )( i, i ) + ( *m_diag )( i, i );
    else
      return ( *m_sqr )( i, j );
  };

  /**
   * Sub-matrix operator, accessor for read only.
   * \test PASSED
   */
  mat_const_sub_block< self > operator()( const std::pair< size_type, size_type >& r,
                                          const std::pair< size_type, size_type >& c ) const {
    return sub ( *this )( r, c );
  };

  /**
   * Sub-matrix operator, accessor for read only.
   * \test PASSED
   */
  mat_const_col_slice< self > operator()( size_type r, const std::pair< size_type, size_type >& c ) const {
    return slice ( *this )( r, c );
  };

  /**
   * Sub-matrix operator, accessor for read only.
   * \test PASSED
   */
  mat_const_row_slice< self > operator()( const std::pair< size_type, size_type >& r, size_type c ) const {
    return slice ( *this )( r, c );
  };

  /**
   * Gets the row-count (number of rows) of the matrix.
   * \return number of rows of the matrix.
   * \test PASSED
   */
  size_type get_row_count() const throw() { return m_sqr->get_row_count(); };
  /**
   * Gets the column-count (number of columns) of the matrix.
   * \return number of columns of the matrix.
   * \test PASSED
   */
  size_type get_col_count() const throw() { return m_sqr->get_col_count(); };

  /**
   * Gets the row-count and column-count of the matrix, as a std::pair of values.
   * \return the row-count and column-count of the matrix, as a std::pair of values.
   * \test PASSED
   */
  std::pair< size_type, size_type > size() const throw() {
    return std::make_pair( m_sqr->get_row_count(), m_sqr->get_col_count() );
  };

  /**
   * Returns the allocator object of the underlying container.
   * \return the allocator object of the underlying container.
   */
  allocator_type get_allocator() const { return m_sqr->get_allocator(); };
};


template < typename SquareMatrix, typename DiagMatrix >
struct is_readable_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > {
  BOOST_STATIC_CONSTANT( bool,
                         value = is_readable_matrix< SquareMatrix >::value && is_readable_matrix< DiagMatrix >::value );
  typedef is_readable_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > type;
};

template < typename SquareMatrix, typename DiagMatrix >
struct is_writable_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_writable_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > type;
};

template < typename SquareMatrix, typename DiagMatrix >
struct is_fully_writable_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_fully_writable_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > type;
};

template < typename SquareMatrix, typename DiagMatrix >
struct is_resizable_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > {
  BOOST_STATIC_CONSTANT( bool, value = false );
  typedef is_resizable_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > type;
};

template < typename SquareMatrix, typename DiagMatrix >
struct has_allocator_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > {
  BOOST_STATIC_CONSTANT( bool, value = has_allocator_matrix< SquareMatrix >::value );
  typedef has_allocator_matrix< SquareMatrix > type;
};


template < typename SquareMatrix, typename DiagMatrix >
mat_damped_matrix< SquareMatrix, DiagMatrix > make_damped_matrix( const SquareMatrix& aMSqr,
                                                                  const DiagMatrix& aMDiag ) {
  return mat_damped_matrix< SquareMatrix, DiagMatrix >( aMSqr, aMDiag );
};


template < typename SquareMatrix, typename DiagMatrix >
struct is_square_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = is_square_matrix< SquareMatrix >::value );
  typedef is_square_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > type;
};


template < typename SquareMatrix, typename DiagMatrix >
struct is_symmetric_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = is_symmetric_matrix< SquareMatrix >::value );
  typedef is_symmetric_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > type;
};


template < typename SquareMatrix, typename DiagMatrix >
struct is_diagonal_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > {
  typedef boost::mpl::integral_c_tag tag;
  typedef bool value_type;
  BOOST_STATIC_CONSTANT( bool, value = is_diagonal_matrix< SquareMatrix >::value );
  typedef is_diagonal_matrix< mat_damped_matrix< SquareMatrix, DiagMatrix > > type;
};
};

#endif

/**
 * \file mat_comparisons.hpp
 *
 * This library provides several functions to compute the various matrix comparisons.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date September 2011
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

#ifndef REAK_MAT_COMPARISONS_HPP
#define REAK_MAT_COMPARISONS_HPP

#include "mat_traits.hpp"
#include "mat_concepts.hpp"

#include <boost/utility/enable_if.hpp>

namespace ReaK {


/**
 * This function template computes the element-wise comparison of two matrices.
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A readable matrix type.
 * \param M1 A matrix for which the 1-norm is sought.
 * \param M2 A matrix for which the 1-norm is sought.
 * \param NumTol The numerical tolerance to consider a value to be zero.
 * \return true iff both matrices are the same, within the given tolerance.
 */
template < typename Matrix1, typename Matrix2 >
typename boost::enable_if_c< is_readable_matrix< Matrix1 >::value && is_readable_matrix< Matrix2 >::value, bool >::type
  is_equal_mat( const Matrix1& M1, const Matrix2& M2, typename mat_traits< Matrix1 >::value_type NumTol
                                                      = typename mat_traits< Matrix1 >::value_type( 1E-8 ) ) {
  if( ( M1.get_row_count() != M2.get_row_count() ) || ( M1.get_col_count() != M2.get_col_count() ) )
    return false;

  typedef typename mat_traits< Matrix1 >::size_type SizeType;
  using std::fabs;

  for( SizeType i = 0; i < M1.get_row_count(); ++i )
    for( SizeType j = 0; j < M1.get_col_count(); ++j )
      if( fabs( M1( i, j ) - M2( i, j ) ) > NumTol )
        return false;

  return true;
};


/**
 * Verifies that the matrix A is diagonal up to a tolerance.
 * \tparam Matrix A readable matrix type.
 * \param A A matrix to verify for being diagonal.
 * \param NumTol The numerical tolerance to consider a value to be zero.
 * \return true iff the matrix is diagonal, within the given tolerance.
 */
template < class Matrix >
typename boost::enable_if_c< is_readable_matrix< Matrix >::value, bool >::type
  is_diagonal( const Matrix& A,
               typename mat_traits< Matrix >::value_type NumTol = typename mat_traits< Matrix >::value_type( 1E-8 ) ) {
  if( A.get_row_count() != A.get_col_count() )
    return false;

  typedef typename mat_traits< Matrix >::size_type SizeType;
  using std::fabs;

  for( SizeType i = 0; i < A.get_row_count(); i++ )
    for( SizeType j = 0; j < i; j++ )
      if( ( fabs( A( j, i ) ) > NumTol ) || ( fabs( A( i, j ) ) > NumTol ) )
        return false;

  return true;
};

/**
 * Verifies that the matrix A is symmetric up to a tolerance.
 * \tparam Matrix A readable matrix type.
 * \param A A matrix to verify for being symmetric.
 * \param NumTol The numerical tolerance to consider a value to be zero.
 * \return true iff the matrix is symmetric, within the given tolerance.
 */
template < class Matrix >
typename boost::enable_if_c< is_readable_matrix< Matrix >::value, bool >::type
  is_symmetric( const Matrix& A,
                typename mat_traits< Matrix >::value_type NumTol = typename mat_traits< Matrix >::value_type( 1E-8 ) ) {
  if( A.get_row_count() != A.get_col_count() )
    return false;

  typedef typename mat_traits< Matrix >::size_type SizeType;
  using std::fabs;

  for( SizeType i = 0; i < A.get_row_count(); i++ )
    for( SizeType j = 0; j < i; j++ )
      if( fabs( A( j, i ) - A( i, j ) ) > NumTol )
        return false;

  return true;
};

/**
 * Verifies that the matrix A is the identity up to a tolerance.
 * \tparam Matrix A readable matrix type.
 * \param A A matrix to verify for being the identity.
 * \param NumTol The numerical tolerance to consider a value to be zero.
 * \return true iff the matrix is the identity, within the given tolerance.
 */
template < class Matrix >
typename boost::enable_if_c< is_readable_matrix< Matrix >::value, bool >::type
  is_identity_mat( const Matrix& A, typename mat_traits< Matrix >::value_type NumTol
                                    = typename mat_traits< Matrix >::value_type( 1E-8 ) ) {
  if( A.get_row_count() != A.get_col_count() )
    return false;

  typedef typename mat_traits< Matrix >::value_type ValueType;
  typedef typename mat_traits< Matrix >::size_type SizeType;
  using std::fabs;

  for( SizeType i = 0; i < A.get_row_count(); i++ )
    for( SizeType j = 0; j < i; j++ )
      if( ( fabs( A( j, i ) ) > NumTol ) || ( fabs( A( i, j ) ) > NumTol ) )
        return false;

  for( SizeType i = 0; i < A.get_row_count(); i++ )
    if( fabs( A( i, i ) - ValueType( 1.0 ) ) > NumTol )
      return false;

  return true;
};

/**
 * Verifies that the matrix A is null up to a tolerance.
 * \tparam Matrix A readable matrix type.
 * \param A A matrix to verify for being null.
 * \param NumTol The numerical tolerance to consider a value to be zero.
 * \return true iff the matrix is null, within the given tolerance.
 */
template < class Matrix >
typename boost::enable_if_c< is_readable_matrix< Matrix >::value, bool >::type
  is_null_mat( const Matrix& A,
               typename mat_traits< Matrix >::value_type NumTol = typename mat_traits< Matrix >::value_type( 1E-8 ) ) {
  typedef typename mat_traits< Matrix >::size_type SizeType;
  using std::fabs;

  for( SizeType i = 0; i < A.get_row_count(); ++i )
    for( SizeType j = 0; j < A.get_col_count(); ++j )
      if( fabs( A( i, j ) ) > NumTol )
        return false;

  return true;
};


/**
 * Verifies that the matrix A is upper-triangular up to a tolerance.
 * \tparam Matrix A readable matrix type.
 * \param A A matrix to verify for being upper-triangular.
 * \param NumTol The numerical tolerance to consider a value to be zero.
 * \return true iff the matrix is upper-triangular, within the given tolerance.
 */
template < class Matrix >
typename boost::enable_if_c< is_readable_matrix< Matrix >::value, bool >::type
  is_upper_triangular( const Matrix& A, typename mat_traits< Matrix >::value_type NumTol
                                        = typename mat_traits< Matrix >::value_type( 1E-8 ) ) {

  typedef typename mat_traits< Matrix >::size_type SizeType;
  using std::fabs;
  SizeType N = A.get_row_count();
  if( N > A.get_col_count() )
    N = A.get_col_count();

  for( SizeType i = 0; i + 1 < N; i++ )
    for( SizeType j = i + 1; j < A.get_row_count(); j++ )
      if( fabs( A( j, i ) ) > NumTol )
        return false;

  return true;
};

/**
 * Verifies that the matrix A is lower-triangular up to a tolerance.
 * \tparam Matrix A readable matrix type.
 * \param A A matrix to verify for being lower-triangular.
 * \param NumTol The numerical tolerance to consider a value to be zero.
 * \return true iff the matrix is lower-triangular, within the given tolerance.
 */
template < class Matrix >
typename boost::enable_if_c< is_readable_matrix< Matrix >::value, bool >::type
  is_lower_triangular( const Matrix& A, typename mat_traits< Matrix >::value_type NumTol
                                        = typename mat_traits< Matrix >::value_type( 1E-8 ) ) {

  typedef typename mat_traits< Matrix >::size_type SizeType;
  using std::fabs;

  for( SizeType i = 1; i < A.get_col_count(); i++ )
    for( SizeType j = 0; ( j < i ) && ( j < A.get_row_count() ); j++ )
      if( fabs( A( j, i ) ) > NumTol )
        return false;

  return true;
};


/**
 * Verifies that the matrix A is upper-hessenberg up to a tolerance.
 * \tparam Matrix A readable matrix type.
 * \param A A matrix to verify for being upper-hessenberg.
 * \param NumTol The numerical tolerance to consider a value to be zero.
 * \return true iff the matrix is upper-hessenberg, within the given tolerance.
 */
template < class Matrix >
typename boost::enable_if_c< is_readable_matrix< Matrix >::value, bool >::type
  is_upper_hessenberg( const Matrix& A, typename mat_traits< Matrix >::value_type NumTol
                                        = typename mat_traits< Matrix >::value_type( 1E-8 ) ) {

  typedef typename mat_traits< Matrix >::size_type SizeType;
  using std::fabs;
  SizeType N = A.get_row_count();
  if( N > A.get_col_count() )
    N = A.get_col_count();

  for( SizeType i = 0; i + 2 < N; i++ )
    for( SizeType j = i + 2; j < A.get_row_count(); j++ )
      if( fabs( A( j, i ) ) > NumTol )
        return false;

  return true;
};

/**
 * Verifies that the matrix A is lower-hessenberg up to a tolerance.
 * \tparam Matrix A readable matrix type.
 * \param A A matrix to verify for being lower-hessenberg.
 * \param NumTol The numerical tolerance to consider a value to be zero.
 * \return true iff the matrix is lower-hessenberg, within the given tolerance.
 */
template < class Matrix >
typename boost::enable_if_c< is_readable_matrix< Matrix >::value, bool >::type
  is_lower_hessenberg( const Matrix& A, typename mat_traits< Matrix >::value_type NumTol
                                        = typename mat_traits< Matrix >::value_type( 1E-8 ) ) {

  typedef typename mat_traits< Matrix >::size_type SizeType;
  using std::fabs;

  for( SizeType i = 2; i < A.get_col_count(); i++ )
    for( SizeType j = 0; ( j + 1 < i ) && ( j < A.get_row_count() ); j++ )
      if( fabs( A( j, i ) ) > NumTol )
        return false;

  return true;
};


/**
 * Verifies that the matrix A is tri-diagonal up to a tolerance.
 * \tparam Matrix A readable matrix type.
 * \param A A matrix to verify for being tri-diagonal.
 * \param NumTol The numerical tolerance to consider a value to be zero.
 * \return true iff the matrix is tri-diagonal, within the given tolerance.
 */
template < class Matrix >
typename boost::enable_if_c< is_readable_matrix< Matrix >::value, bool >::type
  is_tri_diagonal( const Matrix& A, typename mat_traits< Matrix >::value_type NumTol
                                    = typename mat_traits< Matrix >::value_type( 1E-8 ) ) {
  if( A.get_row_count() != A.get_col_count() )
    return false;

  typedef typename mat_traits< Matrix >::size_type SizeType;
  using std::fabs;

  for( SizeType i = 2; i < A.get_col_count(); i++ )
    for( SizeType j = 0; ( j + 1 < i ); j++ )
      if( ( fabs( A( j, i ) ) > NumTol ) || ( fabs( A( i, j ) ) > NumTol ) )
        return false;

  return true;
};
};

#endif

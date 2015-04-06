/**
 * \file mat_norms.hpp
 *
 * This library provides several functions to compute the various matrix norms.
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

#ifndef REAK_MAT_NORMS_HPP
#define REAK_MAT_NORMS_HPP

#include "mat_traits.hpp"
#include "mat_concepts.hpp"

namespace ReaK {


/**
 * This function template computes the 1-norm of a matrix.
 * \tparam Matrix A readable matrix type.
 * \param M A matrix for which the 1-norm is sought.
 * \return the 1-norm of matrix M.
 */
template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, typename Matrix::value_type >::type
  norm_1( const Matrix& M ) {
  typedef typename mat_traits< Matrix >::value_type ValueType;
  typedef typename mat_traits< Matrix >::size_type SizeType;
  using std::fabs;

  ValueType max = ValueType();
  for( SizeType j = 0; j < M.get_col_count(); ++j ) {
    ValueType sum = ValueType();
    for( SizeType i = 0; i < M.get_row_count(); ++i )
      sum += fabs( M( i, j ) );
    if( sum > max )
      max = sum;
  };
  return max;
};


/**
 * This function template computes the infinity-norm of a matrix.
 * \tparam Matrix A readable matrix type.
 * \param M A matrix for which the infinity-norm is sought.
 * \return the infinity-norm of matrix M.
 */
template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, typename Matrix::value_type >::type
  norm_inf( const Matrix& M ) {
  typedef typename mat_traits< Matrix >::value_type ValueType;
  typedef typename mat_traits< Matrix >::size_type SizeType;
  using std::fabs;

  ValueType max = ValueType();
  for( SizeType i = 0; i < M.get_row_count(); ++i ) {
    ValueType sum = ValueType();
    for( SizeType j = 0; j < M.get_col_count(); ++j )
      sum += fabs( M( i, j ) );
    if( sum > max )
      max = sum;
  };
  return max;
};


/**
 * This function template computes the element-wise 2-norm of a matrix.
 * \tparam Matrix A readable matrix type.
 * \param M A matrix for which the element-wise 2-norm is sought.
 * \return the element-wise 2-norm of matrix M.
 */
template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, typename Matrix::value_type >::type
  elem_norm_2( const Matrix& M ) {
  typedef typename mat_traits< Matrix >::value_type ValueType;
  typedef typename mat_traits< Matrix >::size_type SizeType;
  using std::sqrt;

  ValueType sum = ValueType();
  for( SizeType i = 0; i < M.get_row_count(); ++i )
    for( SizeType j = 0; j < M.get_col_count(); ++j )
      sum += M( i, j ) * M( i, j );
  return sqrt( sum );
};


/**
 * This function template computes the Frobenius-norm of a matrix.
 * \tparam Matrix A readable matrix type.
 * \param M A matrix for which the Frobenius-norm is sought.
 * \return the Frobenius-norm of matrix M.
 */
template < typename Matrix >
typename Matrix::value_type frobenius_norm( const Matrix& M ) {
  return elem_norm_2( M );
};


/**
 * This function template computes the element-wise infinity-norm of a matrix.
 * \tparam Matrix A readable matrix type.
 * \param M A matrix for which the element-wise infinity-norm is sought.
 * \return the element-wise infinity-norm of matrix M.
 */
template < typename Matrix >
typename boost::enable_if< is_readable_matrix< Matrix >, typename Matrix::value_type >::type
  elem_norm_max( const Matrix& M ) {
  typedef typename mat_traits< Matrix >::value_type ValueType;
  typedef typename mat_traits< Matrix >::size_type SizeType;
  using std::fabs;

  ValueType max = ValueType();
  for( SizeType i = 0; i < M.get_row_count(); ++i )
    for( SizeType j = 0; j < M.get_col_count(); ++j )
      if( max < fabs( M( i, j ) ) )
        max = fabs( M( i, j ) );
  return max;
};


#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

extern template double norm_1( const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M );
extern template double norm_1( const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M );
extern template double norm_1( const mat< double, mat_structure::square, mat_alignment::column_major >& M );
extern template double norm_1( const mat< double, mat_structure::square, mat_alignment::row_major >& M );
extern template double norm_1( const mat< double, mat_structure::symmetric >& M );
extern template double norm_1( const mat< double, mat_structure::skew_symmetric >& M );
extern template double norm_1( const mat< double, mat_structure::diagonal >& M );
extern template double norm_1( const mat< double, mat_structure::scalar >& M );
extern template double norm_1( const mat< double, mat_structure::identity >& M );
extern template double norm_1( const mat< double, mat_structure::nil >& M );

extern template double norm_inf( const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M );
extern template double norm_inf( const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M );
extern template double norm_inf( const mat< double, mat_structure::square, mat_alignment::column_major >& M );
extern template double norm_inf( const mat< double, mat_structure::square, mat_alignment::row_major >& M );
extern template double norm_inf( const mat< double, mat_structure::symmetric >& M );
extern template double norm_inf( const mat< double, mat_structure::skew_symmetric >& M );
extern template double norm_inf( const mat< double, mat_structure::diagonal >& M );
extern template double norm_inf( const mat< double, mat_structure::scalar >& M );
extern template double norm_inf( const mat< double, mat_structure::identity >& M );
extern template double norm_inf( const mat< double, mat_structure::nil >& M );

extern template double elem_norm_2( const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M );
extern template double elem_norm_2( const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M );
extern template double elem_norm_2( const mat< double, mat_structure::square, mat_alignment::column_major >& M );
extern template double elem_norm_2( const mat< double, mat_structure::square, mat_alignment::row_major >& M );
extern template double elem_norm_2( const mat< double, mat_structure::symmetric >& M );
extern template double elem_norm_2( const mat< double, mat_structure::skew_symmetric >& M );
extern template double elem_norm_2( const mat< double, mat_structure::diagonal >& M );
extern template double elem_norm_2( const mat< double, mat_structure::scalar >& M );
extern template double elem_norm_2( const mat< double, mat_structure::identity >& M );
extern template double elem_norm_2( const mat< double, mat_structure::nil >& M );

extern template double elem_norm_max( const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M );
extern template double elem_norm_max( const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M );
extern template double elem_norm_max( const mat< double, mat_structure::square, mat_alignment::column_major >& M );
extern template double elem_norm_max( const mat< double, mat_structure::square, mat_alignment::row_major >& M );
extern template double elem_norm_max( const mat< double, mat_structure::symmetric >& M );
extern template double elem_norm_max( const mat< double, mat_structure::skew_symmetric >& M );
extern template double elem_norm_max( const mat< double, mat_structure::diagonal >& M );
extern template double elem_norm_max( const mat< double, mat_structure::scalar >& M );
extern template double elem_norm_max( const mat< double, mat_structure::identity >& M );
extern template double elem_norm_max( const mat< double, mat_structure::nil >& M );


extern template float norm_1( const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M );
extern template float norm_1( const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M );
extern template float norm_1( const mat< float, mat_structure::square, mat_alignment::column_major >& M );
extern template float norm_1( const mat< float, mat_structure::square, mat_alignment::row_major >& M );
extern template float norm_1( const mat< float, mat_structure::symmetric >& M );
extern template float norm_1( const mat< float, mat_structure::skew_symmetric >& M );
extern template float norm_1( const mat< float, mat_structure::diagonal >& M );
extern template float norm_1( const mat< float, mat_structure::scalar >& M );
extern template float norm_1( const mat< float, mat_structure::identity >& M );
extern template float norm_1( const mat< float, mat_structure::nil >& M );

extern template float norm_inf( const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M );
extern template float norm_inf( const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M );
extern template float norm_inf( const mat< float, mat_structure::square, mat_alignment::column_major >& M );
extern template float norm_inf( const mat< float, mat_structure::square, mat_alignment::row_major >& M );
extern template float norm_inf( const mat< float, mat_structure::symmetric >& M );
extern template float norm_inf( const mat< float, mat_structure::skew_symmetric >& M );
extern template float norm_inf( const mat< float, mat_structure::diagonal >& M );
extern template float norm_inf( const mat< float, mat_structure::scalar >& M );
extern template float norm_inf( const mat< float, mat_structure::identity >& M );
extern template float norm_inf( const mat< float, mat_structure::nil >& M );

extern template float elem_norm_2( const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M );
extern template float elem_norm_2( const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M );
extern template float elem_norm_2( const mat< float, mat_structure::square, mat_alignment::column_major >& M );
extern template float elem_norm_2( const mat< float, mat_structure::square, mat_alignment::row_major >& M );
extern template float elem_norm_2( const mat< float, mat_structure::symmetric >& M );
extern template float elem_norm_2( const mat< float, mat_structure::skew_symmetric >& M );
extern template float elem_norm_2( const mat< float, mat_structure::diagonal >& M );
extern template float elem_norm_2( const mat< float, mat_structure::scalar >& M );
extern template float elem_norm_2( const mat< float, mat_structure::identity >& M );
extern template float elem_norm_2( const mat< float, mat_structure::nil >& M );

extern template float elem_norm_max( const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M );
extern template float elem_norm_max( const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M );
extern template float elem_norm_max( const mat< float, mat_structure::square, mat_alignment::column_major >& M );
extern template float elem_norm_max( const mat< float, mat_structure::square, mat_alignment::row_major >& M );
extern template float elem_norm_max( const mat< float, mat_structure::symmetric >& M );
extern template float elem_norm_max( const mat< float, mat_structure::skew_symmetric >& M );
extern template float elem_norm_max( const mat< float, mat_structure::diagonal >& M );
extern template float elem_norm_max( const mat< float, mat_structure::scalar >& M );
extern template float elem_norm_max( const mat< float, mat_structure::identity >& M );
extern template float elem_norm_max( const mat< float, mat_structure::nil >& M );


#endif
};

#endif

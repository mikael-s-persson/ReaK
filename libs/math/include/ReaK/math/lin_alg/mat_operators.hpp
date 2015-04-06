/**
 * \file mat_alg_general_hpp
 *
 * This library implements the general versions of many meta-functions (templates),
 * functions, and operators. These are meant to be used when no more-specialized
 * implementations exist for the matrix types involved.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date april 2011 (originally february 2010)
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

#ifndef REAK_MAT_OPERATORS_HPP
#define REAK_MAT_OPERATORS_HPP

#include <ReaK/core/base/defs.hpp>

//#include "vect_alg.hpp"
#include "vect_concepts.hpp"
#include "mat_concepts.hpp"
#include "mat_traits.hpp"

#include "mat_op_results.hpp"

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/concept_check.hpp>

#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/comparison.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/not_equal_to.hpp>

#include <iomanip>

namespace ReaK {


/**
 * Prints a matrix to a standard output stream (<<) as "((a11; a12); (a21; a22))". \test PASSED
 */
template < typename Matrix >
typename boost::enable_if_c< is_readable_matrix< Matrix >::value, std::ostream& >::type
  operator<<( std::ostream& out_stream, const Matrix& M ) {
  out_stream << "(\n";
  if( ( M.get_row_count() != 0 ) && ( M.get_col_count() != 0 ) ) {
    for( unsigned int i = 0; i < M.get_row_count(); ++i ) {
      out_stream << "(" << std::setw( 16 ) << M( i, 0 );
      for( unsigned int j = 1; j < M.get_col_count(); ++j ) {
        out_stream << "; " << std::setw( 16 ) << M( i, j );
      };
      out_stream << ")";
      if( i != M.get_row_count() - 1 )
        out_stream << ";\n";
      else
        out_stream << "\n";
    };
  };
  return ( out_stream << ")" );
};


/*******************************************************************************
                         Multiplication Operators
*******************************************************************************/


namespace detail {


template < typename Matrix1, typename Matrix2, typename ResultMatrix >
void dense_mat_multiply_impl( const Matrix1& M1, const Matrix2& M2, ResultMatrix& MR ) {
  typedef typename mat_traits< ResultMatrix >::value_type ValueType;
  typedef typename mat_traits< ResultMatrix >::size_type SizeType;
  for( SizeType i = 0; i < M1.get_row_count(); ++i ) {
    for( SizeType jj = 0; jj < M2.get_col_count(); ++jj ) {
      MR( i, jj ) = ValueType( 0.0 );
      for( SizeType j = 0; j < M1.get_col_count(); ++j )
        MR( i, jj ) += M1( i, j ) * M2( j, jj );
    };
  };
};

template < typename Matrix1, typename MatrixDiag, typename ResultMatrix >
void dense_diag_mat_multiply_impl( const Matrix1& M1, const MatrixDiag& M2, ResultMatrix& MR ) {
  typedef typename mat_traits< ResultMatrix >::size_type SizeType;
  for( SizeType i = 0; i < M1.get_row_count(); ++i )
    for( SizeType j = 0; j < M1.get_col_count(); ++j )
      MR( i, j ) = M1( i, j ) * M2( j, j );
};

template < typename MatrixDiag, typename Matrix2, typename ResultMatrix >
void diag_dense_mat_multiply_impl( const MatrixDiag& M1, const Matrix2& M2, ResultMatrix& MR ) {
  typedef typename mat_traits< ResultMatrix >::size_type SizeType;
  for( SizeType i = 0; i < M2.get_row_count(); ++i )
    for( SizeType j = 0; j < M2.get_col_count(); ++j )
      MR( i, j ) = M1( i, i ) * M2( i, j );
};

template < typename MatrixDiag1, typename MatrixDiag2, typename ResultMatrix >
void diag_diag_mat_multiply_impl( const MatrixDiag1& M1, const MatrixDiag2& M2, ResultMatrix& MR ) {
  typedef typename mat_traits< ResultMatrix >::size_type SizeType;
  for( SizeType i = 0; i < M1.get_row_count(); ++i )
    MR( i, i ) = M1( i, i ) * M2( i, i );
};


/* Multiplies two lower-triangular matrices and saves the result in-place into M2. */
template < typename MatrixLower1, typename MatrixLower2 >
void inplace_lower_multiply_impl( const MatrixLower1& M1, MatrixLower2& M2 ) {
  typedef typename mat_traits< MatrixLower2 >::size_type SizeType;
  typedef typename mat_traits< MatrixLower2 >::value_type ValueType;

  for( SizeType i = M1.get_row_count(); i > 0; ) {
    --i;
    for( SizeType j = 0; j <= i; ++j ) {
      ValueType sum = ValueType( 0.0 );
      for( SizeType k = j; k <= i; ++k )
        sum += M1( i, k ) * M2( k, j );
      M2( i, j ) = sum;
    };
  };
};

/* Multiplies two lower-triangular matrices and saves the result in-place into M2, and fills the upper part of M2 with
 * zeros. */
template < typename MatrixLower1, typename MatrixLower2 >
void inplace_lower_multiply_with_fill_impl( const MatrixLower1& M1, MatrixLower2& M2 ) {
  typedef typename mat_traits< MatrixLower2 >::size_type SizeType;
  typedef typename mat_traits< MatrixLower2 >::value_type ValueType;

  for( SizeType i = M1.get_row_count(); i > 0; ) {
    --i;
    for( SizeType j = 0; j <= i; ++j ) {
      ValueType sum = ValueType( 0.0 );
      for( SizeType k = j; k <= i; ++k )
        sum += M1( i, k ) * M2( k, j );
      M2( i, j ) = sum;
    };
    for( SizeType j = i + 1; j < M2.get_col_count(); ++j )
      M2( i, j ) = ValueType( 0.0 );
  };
};
};


/**
 * General multiplication operator for any type of matrices. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M1 first matrix (first operand).
 * \param M2 second matrix (second operand).
 * \return General column-major matrix.
 * \throw std::range_error if the two matrix dimensions do not fit together.
 * \test PASSED
 */
template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                                             boost::mpl::and_< boost::mpl::not_< is_square_matrix< Matrix1 > >,
                                                               boost::mpl::not_< is_square_matrix< Matrix2 > > >,
                                             boost::mpl::less< mat_product_priority< Matrix1 >,
                                                               detail::product_priority< mat_structure::diagonal > >,
                                             boost::mpl::less< mat_product_priority< Matrix2 >,
                                                               detail::product_priority< mat_structure::diagonal > > >,
                           mat_product_result< Matrix1, Matrix2 > >::type
  operator*( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_product_result< Matrix1, Matrix2 >::type result_type;
  BOOST_STATIC_ASSERT( ( is_fully_writable_matrix< result_type >::value ) );
  BOOST_STATIC_ASSERT( ( is_resizable_matrix< result_type >::value ) );
  if( M1.get_col_count() != M2.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );

  result_type result( M1 );
  result.set_col_count( M2.get_col_count() );
  detail::dense_mat_multiply_impl( M1, M2, result );
  return result;
};

// rect and square matrix.
template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                                             boost::mpl::and_< boost::mpl::not_< is_square_matrix< Matrix1 > >,
                                                               is_square_matrix< Matrix2 > >,
                                             boost::mpl::less< mat_product_priority< Matrix1 >,
                                                               detail::product_priority< mat_structure::diagonal > >,
                                             boost::mpl::less< mat_product_priority< Matrix2 >,
                                                               detail::product_priority< mat_structure::diagonal > > >,
                           mat_product_result< Matrix1, Matrix2 > >::type
  operator*( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_product_result< Matrix1, Matrix2 >::type result_type;
  BOOST_STATIC_ASSERT( ( is_fully_writable_matrix< result_type >::value ) );
  if( M1.get_col_count() != M2.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );

  result_type result( M1 );
  detail::dense_mat_multiply_impl( M1, M2, result );
  return result;
};

// square and rect matrix.
template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                                             boost::mpl::and_< is_square_matrix< Matrix1 >,
                                                               boost::mpl::not_< is_square_matrix< Matrix2 > > >,
                                             boost::mpl::less< mat_product_priority< Matrix1 >,
                                                               detail::product_priority< mat_structure::diagonal > >,
                                             boost::mpl::less< mat_product_priority< Matrix2 >,
                                                               detail::product_priority< mat_structure::diagonal > > >,
                           mat_product_result< Matrix1, Matrix2 > >::type
  operator*( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_product_result< Matrix1, Matrix2 >::type result_type;
  BOOST_STATIC_ASSERT( ( is_fully_writable_matrix< result_type >::value ) );
  if( M1.get_col_count() != M2.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );

  result_type result( M2 );
  detail::dense_mat_multiply_impl( M1, M2, result );
  return result;
};

// square and square matrix.
template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::and_< is_square_matrix< Matrix1 >, is_square_matrix< Matrix2 > >,
                               boost::mpl::less< mat_product_priority< Matrix1 >,
                                                 detail::product_priority< mat_structure::diagonal > >,
                               boost::mpl::less< mat_product_priority< Matrix2 >,
                                                 detail::product_priority< mat_structure::diagonal > > >,
             mat_product_result< Matrix1, Matrix2 > >::type
  operator*( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_product_result< Matrix1, Matrix2 >::type result_type;
  BOOST_STATIC_ASSERT( ( is_fully_writable_matrix< result_type >::value ) );
  if( M1.get_col_count() != M2.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );

  result_type result( M1 );
  detail::dense_mat_multiply_impl( M1, M2, result );
  return result;
};

// dense and diagonal matrix.
template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::less< mat_product_priority< Matrix1 >,
                                                 detail::product_priority< mat_structure::diagonal > >,
                               boost::mpl::equal_to< mat_product_priority< Matrix2 >,
                                                     detail::product_priority< mat_structure::diagonal > > >,
             mat_product_result< Matrix1, Matrix2 > >::type
  operator*( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_product_result< Matrix1, Matrix2 >::type result_type;
  BOOST_STATIC_ASSERT( ( is_fully_writable_matrix< result_type >::value ) );
  if( M1.get_col_count() != M2.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );

  result_type result( M1 );
  detail::dense_diag_mat_multiply_impl( M1, M2, result );
  return result;
};

// diagonal and dense matrix.
template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::equal_to< mat_product_priority< Matrix1 >,
                                                     detail::product_priority< mat_structure::diagonal > >,
                               boost::mpl::less< mat_product_priority< Matrix2 >,
                                                 detail::product_priority< mat_structure::diagonal > > >,
             mat_product_result< Matrix1, Matrix2 > >::type
  operator*( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_product_result< Matrix1, Matrix2 >::type result_type;
  BOOST_STATIC_ASSERT( ( is_fully_writable_matrix< result_type >::value ) );
  if( M1.get_col_count() != M2.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );

  result_type result( M2 );
  detail::diag_dense_mat_multiply_impl( M1, M2, result );
  return result;
};

// diagonal and diagonal matrix.
template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::equal_to< mat_product_priority< Matrix1 >,
                                                     detail::product_priority< mat_structure::diagonal > >,
                               boost::mpl::equal_to< mat_product_priority< Matrix2 >,
                                                     detail::product_priority< mat_structure::diagonal > > >,
             mat_product_result< Matrix1, Matrix2 > >::type
  operator*( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_product_result< Matrix1, Matrix2 >::type result_type;
  BOOST_STATIC_ASSERT( ( is_fully_writable_matrix< result_type >::value ) );
  if( M1.get_col_count() != M2.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );

  result_type result( M2 );
  detail::diag_diag_mat_multiply_impl( M1, M2, result );
  return result;
};


// nil matrix.
template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::
               and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                     boost::mpl::or_< boost::mpl::equal_to< mat_product_priority< Matrix1 >,
                                                            detail::product_priority< mat_structure::nil > >,
                                      boost::mpl::equal_to< mat_product_priority< Matrix2 >,
                                                            detail::product_priority< mat_structure::nil > > > >,
             mat_product_result< Matrix1, Matrix2 > >::type
  operator*( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_product_result< Matrix1, Matrix2 >::type result_type;
  BOOST_STATIC_ASSERT( ( is_resizable_matrix< result_type >::value ) );
  if( M1.get_col_count() != M2.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );

  result_type result;
  result.set_row_count( M1.get_row_count() );
  result.set_col_count( M2.get_col_count() );
  return result;
};


// scalar and normal matrix.
template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                                             boost::mpl::equal_to< mat_product_priority< Matrix1 >,
                                                                   detail::product_priority< mat_structure::scalar > >,
                                             boost::mpl::less< mat_product_priority< Matrix2 >,
                                                               detail::product_priority< mat_structure::scalar > > >,
                           mat_product_result< Matrix1, Matrix2 > >::type
  operator*( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_product_result< Matrix1, Matrix2 >::type result_type;
  if( M1.get_col_count() != M2.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );

  result_type result( M2 );
  if( M1.get_row_count() )
    result *= M1( 0, 0 );
  return result;
};


// normal and scalar matrix.
template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::less_equal< mat_product_priority< Matrix1 >,
                                                       detail::product_priority< mat_structure::scalar > >,
                               boost::mpl::equal_to< mat_product_priority< Matrix2 >,
                                                     detail::product_priority< mat_structure::scalar > > >,
             mat_product_result< Matrix1, Matrix2 > >::type
  operator*( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_product_result< Matrix1, Matrix2 >::type result_type;
  if( M1.get_col_count() != M2.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );

  result_type result( M1 );
  if( M1.get_row_count() )
    result *= M2( 0, 0 );
  return result;
};


// identity and normal matrix.
template < typename Matrix1, typename Matrix2 >
typename boost::
  lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::equal_to< mat_product_priority< Matrix1 >,
                                                     detail::product_priority< mat_structure::identity > >,
                               boost::mpl::less< mat_product_priority< Matrix2 >,
                                                 detail::product_priority< mat_structure::identity > > >,
             mat_product_result< Matrix1, Matrix2 > >::type
  operator*( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_product_result< Matrix1, Matrix2 >::type result_type;
  if( M1.get_col_count() != M2.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );

  return result_type( M2 );
};


// normal and identity matrix.
template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::less_equal< mat_product_priority< Matrix1 >,
                                                       detail::product_priority< mat_structure::identity > >,
                               boost::mpl::equal_to< mat_product_priority< Matrix2 >,
                                                     detail::product_priority< mat_structure::identity > > >,
             mat_product_result< Matrix1, Matrix2 > >::type
  operator*( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_product_result< Matrix1, Matrix2 >::type result_type;
  if( M1.get_col_count() != M2.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );

  return result_type( M1 );
};


/**
 * Matrix multiplication operator with a scalar. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M The matrix (first operand).
 * \param S The scalar (second operand).
 * \return Column-major matrix equal to M * S.
 * \test PASSED
 */
template < typename Matrix, typename Scalar >
typename boost::enable_if< boost::mpl::and_< is_writable_matrix< Matrix >,
                                             boost::mpl::not_< is_readable_matrix< Scalar > >,
                                             boost::mpl::not_< is_readable_vector< Scalar > > >,
                           Matrix >::type
  operator*( Matrix M, const Scalar& S ) {
  M *= S;
  return M;
};

/**
 * Matrix multiplication operator with a scalar. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param S The scalar (first operand).
 * \param M The matrix (second operand).
 * \return Column-major matrix equal to S * M.
 * \test PASSED
 */
template < typename Matrix, typename Scalar >
typename boost::enable_if< boost::mpl::and_< is_writable_matrix< Matrix >,
                                             boost::mpl::not_< is_readable_matrix< Scalar > >,
                                             boost::mpl::not_< is_readable_vector< Scalar > > >,
                           Matrix >::type
  operator*( const Scalar& S, Matrix M ) {
  M *= S;
  return M;
};


/*******************************************************************************
                         Addition / Subtraction Operators
*******************************************************************************/


/**
 * General (least-specialized) unary-negation operator for any type of matrices.
 * This is a default operator that will be called if no better special-purpose overload exists.
 * \param M first matrix (first operand).
 * \return General column-major matrix.
 * \throw std::range_error if the two matrix dimensions do not fit together.
 * \test PASSED
 */
template < typename Matrix >
typename boost::lazy_enable_if< is_readable_matrix< Matrix >, mat_addition_result< Matrix, Matrix > >::type
  operator-( const Matrix& M ) {
  typedef typename mat_addition_result< Matrix, Matrix >::type result_type;
  typedef typename mat_traits< result_type >::size_type SizeType;
  result_type result( M );
  for( SizeType j = 0; j < M.get_col_count(); ++j )
    for( SizeType i = 0; i < M.get_row_count(); ++i )
      result( i, j ) = -M( i, j ); // this needs to be taking the original matrix (avoid cross-referenced elements).
  return result;
};


/**
 * General (least-specialized) transpose function for any type of matrix.
 * \param M The matrix to be transposed.
 * \return The transpose of M.
 */
template < typename Matrix >
typename boost::lazy_enable_if< is_readable_matrix< Matrix >, mat_addition_result< Matrix, Matrix > >::type
  transpose( const Matrix& M ) {
  typedef typename mat_addition_result< Matrix, Matrix >::type result_type;
  typedef typename mat_traits< result_type >::size_type SizeType;
  result_type result( M );
  for( SizeType j = 0; j < result.get_col_count(); ++j )
    for( SizeType i = 0; i < result.get_row_count(); ++i )
      result( i, j ) = M( j, i );
  return result;
};

/**
 * General (least-specialized) transpose function for any type of matrix.
 * \param M The matrix to be transposed.
 * \return The transpose of M.
 */
template < typename Matrix >
typename boost::lazy_enable_if< is_readable_matrix< Matrix >, mat_addition_result< Matrix, Matrix > >::type
  transpose_move( const Matrix& M ) {
  return transpose( M );
};


/**
 * General (least-specialized) addition operator for any type of matrices. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M1 first matrix (first operand).
 * \param M2 second matrix (second operand).
 * \return General column-major matrix.
 * \throw std::range_error if the two matrix dimensions do not fit together.
 * \test PASSED
 */
template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::less< mat_addition_priority< Matrix1 >,
                                                 detail::addition_priority< mat_structure::upper_triangular > >,
                               boost::mpl::less< mat_addition_priority< Matrix2 >,
                                                 detail::addition_priority< mat_structure::upper_triangular > > >,
             mat_addition_result< Matrix1, Matrix2 > >::type
  operator+( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_addition_result< Matrix1, Matrix2 >::type result_type;
  if( ( M1.get_row_count() != M2.get_row_count() ) || ( M1.get_col_count() != M2.get_col_count() ) )
    throw std::range_error( "Matrix dimension mismatch." );
  typedef typename mat_traits< result_type >::size_type SizeType;
  result_type result( M1 );
  for( SizeType j = 0; j < M1.get_col_count(); ++j )
    for( SizeType i = 0; i < M1.get_row_count(); ++i )
      result( i, j ) += M2( i, j );
  return result;
};


/**
 * General substraction operator for any type of matrices. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M1 first matrix (first operand).
 * \param M2 second matrix (second operand).
 * \return General column-major matrix.
 * \throw std::range_error if the two matrix dimensions do not fit together.
 * \test PASSED
 */
template < typename Matrix1, typename Matrix2 >
typename boost::
  lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::less< mat_addition_priority< Matrix1 >,
                                                 detail::addition_priority< mat_structure::upper_triangular > >,
                               boost::mpl::less< mat_addition_priority< Matrix2 >,
                                                 detail::addition_priority< mat_structure::upper_triangular > > >,
             mat_addition_result< Matrix1, Matrix2 > >::type
  operator-( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_addition_result< Matrix1, Matrix2 >::type result_type;
  if( ( M1.get_row_count() != M2.get_row_count() ) || ( M1.get_col_count() != M2.get_col_count() ) )
    throw std::range_error( "Matrix dimension mismatch." );
  typedef typename mat_traits< result_type >::size_type SizeType;
  result_type result( M1 );
  for( SizeType j = 0; j < M1.get_col_count(); ++j )
    for( SizeType i = 0; i < M1.get_row_count(); ++i )
      result( i, j ) -= M2( i, j );
  return result;
};


template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::equal_to< mat_addition_priority< Matrix1 >,
                                                     detail::addition_priority< mat_structure::diagonal > >,
                               boost::mpl::equal_to< mat_addition_priority< Matrix2 >,
                                                     detail::addition_priority< mat_structure::diagonal > > >,
             mat_addition_result< Matrix1, Matrix2 > >::type
  operator+( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_addition_result< Matrix1, Matrix2 >::type result_type;
  if( M1.get_row_count() != M2.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );
  result_type result( M1 );
  result += result_type( M2 );
  return result;
};


template < typename Matrix1, typename Matrix2 >
typename boost::
  lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::equal_to< mat_addition_priority< Matrix1 >,
                                                     detail::addition_priority< mat_structure::diagonal > >,
                               boost::mpl::equal_to< mat_addition_priority< Matrix2 >,
                                                     detail::addition_priority< mat_structure::diagonal > > >,
             mat_addition_result< Matrix1, Matrix2 > >::type
  operator-( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_addition_result< Matrix1, Matrix2 >::type result_type;
  if( M1.get_row_count() != M2.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );
  result_type result( M1 );
  result -= result_type( M2 );
  return result;
};


template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::equal_to< mat_addition_priority< Matrix1 >,
                                                     detail::addition_priority< mat_structure::diagonal > >,
                               boost::mpl::less< mat_addition_priority< Matrix2 >,
                                                 detail::addition_priority< mat_structure::diagonal > > >,
             mat_addition_result< Matrix1, Matrix2 > >::type
  operator+( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_addition_result< Matrix1, Matrix2 >::type result_type;
  if( ( M1.get_row_count() != M2.get_row_count() ) || ( M1.get_col_count() != M2.get_col_count() ) )
    throw std::range_error( "Matrix dimension mismatch." );
  result_type result( M2 );
  typedef typename mat_traits< result_type >::size_type SizeType;
  for( SizeType i = 0; i < result.get_row_count(); ++i )
    result( i, i ) += M1( i, i );
  return result;
};


template < typename Matrix1, typename Matrix2 >
typename boost::
  lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::equal_to< mat_addition_priority< Matrix1 >,
                                                     detail::addition_priority< mat_structure::diagonal > >,
                               boost::mpl::less< mat_addition_priority< Matrix2 >,
                                                 detail::addition_priority< mat_structure::diagonal > > >,
             mat_addition_result< Matrix1, Matrix2 > >::type
  operator-( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_addition_result< Matrix1, Matrix2 >::type result_type;
  if( ( M1.get_row_count() != M2.get_row_count() ) || ( M1.get_col_count() != M2.get_col_count() ) )
    throw std::range_error( "Matrix dimension mismatch." );
  result_type result( -M2 );
  typedef typename mat_traits< result_type >::size_type SizeType;
  for( SizeType i = 0; i < result.get_row_count(); ++i )
    result( i, i ) += M1( i, i );
  return result;
};


template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::less< mat_addition_priority< Matrix1 >,
                                                 detail::addition_priority< mat_structure::diagonal > >,
                               boost::mpl::equal_to< mat_addition_priority< Matrix2 >,
                                                     detail::addition_priority< mat_structure::diagonal > > >,
             mat_addition_result< Matrix1, Matrix2 > >::type
  operator+( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_addition_result< Matrix1, Matrix2 >::type result_type;
  if( ( M1.get_row_count() != M2.get_row_count() ) || ( M1.get_col_count() != M2.get_col_count() ) )
    throw std::range_error( "Matrix dimension mismatch." );
  result_type result( M1 );
  typedef typename mat_traits< result_type >::size_type SizeType;
  for( SizeType i = 0; i < result.get_row_count(); ++i )
    result( i, i ) += M2( i, i );
  return result;
};


template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::less< mat_addition_priority< Matrix1 >,
                                                 detail::addition_priority< mat_structure::diagonal > >,
                               boost::mpl::equal_to< mat_addition_priority< Matrix2 >,
                                                     detail::addition_priority< mat_structure::diagonal > > >,
             mat_addition_result< Matrix1, Matrix2 > >::type
  operator-( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_addition_result< Matrix1, Matrix2 >::type result_type;
  if( ( M1.get_row_count() != M2.get_row_count() ) || ( M1.get_col_count() != M2.get_col_count() ) )
    throw std::range_error( "Matrix dimension mismatch." );
  result_type result( M1 );
  typedef typename mat_traits< result_type >::size_type SizeType;
  for( SizeType i = 0; i < result.get_row_count(); ++i )
    result( i, i ) -= M2( i, i );
  return result;
};


template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::equal_to< mat_addition_priority< Matrix1 >,
                                                     detail::addition_priority< mat_structure::nil > >,
                               boost::mpl::not_equal_to< mat_addition_priority< Matrix2 >,
                                                         detail::addition_priority< mat_structure::nil > > >,
             mat_addition_result< Matrix1, Matrix2 > >::type
  operator+( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_addition_result< Matrix1, Matrix2 >::type result_type;
  if( ( M1.get_row_count() != M2.get_row_count() ) || ( M1.get_col_count() != M2.get_col_count() ) )
    throw std::range_error( "Matrix dimension mismatch." );
  return result_type( M2 );
};


template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                               boost::mpl::equal_to< mat_addition_priority< Matrix1 >,
                                                     detail::addition_priority< mat_structure::nil > >,
                               boost::mpl::not_equal_to< mat_addition_priority< Matrix2 >,
                                                         detail::addition_priority< mat_structure::nil > > >,
             mat_addition_result< Matrix1, Matrix2 > >::type
  operator-( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_addition_result< Matrix1, Matrix2 >::type result_type;
  if( ( M1.get_row_count() != M2.get_row_count() ) || ( M1.get_col_count() != M2.get_col_count() ) )
    throw std::range_error( "Matrix dimension mismatch." );
  return result_type( -M2 );
};


template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                                             boost::mpl::equal_to< mat_addition_priority< Matrix2 >,
                                                                   detail::addition_priority< mat_structure::nil > > >,
                           mat_addition_result< Matrix1, Matrix2 > >::type
  operator+( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_addition_result< Matrix1, Matrix2 >::type result_type;
  if( ( M1.get_row_count() != M2.get_row_count() ) || ( M1.get_col_count() != M2.get_col_count() ) )
    throw std::range_error( "Matrix dimension mismatch." );
  return result_type( M1 );
};


template < typename Matrix1, typename Matrix2 >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix1 >, is_readable_matrix< Matrix2 >,
                                             boost::mpl::equal_to< mat_addition_priority< Matrix2 >,
                                                                   detail::addition_priority< mat_structure::nil > > >,
                           mat_addition_result< Matrix1, Matrix2 > >::type
  operator-( const Matrix1& M1, const Matrix2& M2 ) {
  typedef typename mat_addition_result< Matrix1, Matrix2 >::type result_type;
  if( ( M1.get_row_count() != M2.get_row_count() ) || ( M1.get_col_count() != M2.get_col_count() ) )
    throw std::range_error( "Matrix dimension mismatch." );
  return result_type( M1 );
};


#if 0
/**
 * General multiplication operator for any type of matrices. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M1 first matrix (first operand).
 * \param M2 second matrix (second operand).
 * \return General column-major matrix.
 * \throw std::range_error if the two matrix dimensions do not fit together.
 * \test PASSED
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix2>::value &&
                             boost::mpl::less< boost::mpl::integral_c< std::size_t, mat_product_priority< mat<T,Structure,Alignment,Allocator> >::value>,
                                               boost::mpl::integral_c< std::size_t, detail::product_priority<mat_structure::diagonal>::value> >::value &&
                             boost::mpl::less< boost::mpl::integral_c< std::size_t, mat_product_priority<Matrix2>::value>,
                                               boost::mpl::integral_c< std::size_t, detail::product_priority<mat_structure::diagonal>::value> >::value,
 mat< T, mat_structure::rectangular, Alignment, Allocator > >::type
  operator *(const Matrix2& M1, const mat<T,Structure,Alignment,Allocator>& M2) {
    typedef mat< T, mat_structure::rectangular, Alignment, Allocator > result_type;
    if(M1.get_col_count() != M2.get_row_count())
      throw std::range_error("Matrix dimension mismatch.");
    typedef typename mat_traits<result_type>::value_type ValueType;
    typedef typename mat_traits<result_type>::size_type SizeType;
    result_type result(M1.get_row_count(),M2.get_col_count(),T(0),M2.get_allocator());
    for(SizeType i=0;i<M1.get_row_count();++i) {
      for(SizeType jj=0;jj<M2.get_col_count();++jj) {
        for(SizeType j=0;j<M1.get_col_count();++j) {
          result(i,jj) += M1(i,j) * M2(j,jj);
        };
      };
    };
    return result;
  };
#endif
/**
 * Matrix multiplication operator with a column vector. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param M The matrix (first operand).
 * \param V The column vector (second operand).
 * \return Column vector equal to M * V.
 * \throw std::range_error if the matrix column count does not correspond to the vector dimension.
 * \test PASSED
 */
template < typename Matrix, typename Vector >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix >, is_readable_vector< Vector > >,
                           vect_copy< Vector > >::type
  operator*( const Matrix& M, const Vector& V ) {
  typedef typename vect_copy< Vector >::type result_type;
  if( V.size() != M.get_col_count() )
    throw std::range_error( "Matrix dimension mismatch." );
  typedef typename mat_traits< Matrix >::value_type ValueType;
  typedef typename mat_traits< Matrix >::size_type SizeType;
  result_type result( M.get_row_count() );
  for( SizeType i = 0; i < M.get_row_count(); ++i ) {
    result[i] = ValueType( 0.0 );
    for( SizeType j = 0; j < M.get_col_count(); ++j )
      result[i] += M( i, j ) * V[j];
  };
  return result;
};

/**
 * Matrix multiplication operator with a row vector. This is a default operator
 * that will be called if no better special-purpose overload exists.
 * \param V The row vector (first operand).
 * \param M The matrix (second operand).
 * \return Row vector equal to V * M.
 * \throw std::range_error if the matrix row count does not correspond to the vector dimension.
 * \test PASSED
 */
template < typename Matrix, typename Vector >
typename boost::lazy_enable_if< boost::mpl::and_< is_readable_matrix< Matrix >, is_readable_vector< Vector > >,
                           vect_copy< Vector > >::type
  operator*( const Vector& V, const Matrix& M ) {
  typedef typename vect_copy< Vector >::type result_type;
  if( V.size() != M.get_row_count() )
    throw std::range_error( "Matrix dimension mismatch." );
  typedef typename mat_traits< Matrix >::value_type ValueType;
  typedef typename mat_traits< Matrix >::size_type SizeType;
  result_type result( M.get_col_count() );
  for( SizeType j = 0; j < M.get_col_count(); ++j ) {
    result[j] = ValueType( 0.0 );
    for( SizeType i = 0; i < M.get_row_count(); ++i ) {
      result[j] += M( i, j ) * V[i];
    };
  };
  return result;
};


#if 0
/**
 * Multiplication by a column-vector (fixed-size).
 * \param M the matrix (square)
 * \param V the column vector.
 * \return column-vector equal to M * V.
 * \throw std::range_error if this matrix and the vector dimensions don't match.
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, unsigned int Size>
vect<T,Size> operator *(const mat<T,Structure,Alignment,Allocator>& M,const vect<T,Size>& V) {
  if((Size != M.get_col_count()) || (Size != M.get_row_count()))
    throw std::range_error("Matrix dimension mismatch.");
  typedef typename mat_traits< mat<T,Structure,Alignment,Allocator> >::value_type ValueType;
  typedef typename mat_traits< mat<T,Structure,Alignment,Allocator> >::size_type SizeType;
  vect<T,Size> result;
  for(SizeType i=0;i<Size;++i) {
    result[i] = 0;
    for(SizeType j=0;j<Size;++j)
      result[i] += M(i,j) * V[j];
  };
  return result;
};

/**
 * Multiplication by a row-vector (fixed-size).
 * \param M the matrix (square)
 * \param V the column vector.
 * \return row-vector equal to V * M.
 * \throw std::range_error if this matrix and the vector dimensions don't match.
 */
template <typename T, mat_structure::tag Structure, mat_alignment::tag Alignment, typename Allocator, unsigned int Size>
vect<T,Size> operator *(const vect<T,Size>& V,const mat<T,Structure,Alignment,Allocator>& M) {
  if((Size != M.get_col_count()) || (Size != M.get_row_count()))
    throw std::range_error("Matrix dimension mismatch.");
  typedef typename mat_traits< mat<T,Structure,Alignment,Allocator> >::value_type ValueType;
  typedef typename mat_traits< mat<T,Structure,Alignment,Allocator> >::size_type SizeType;
  vect<T,Size> result;
  for(SizeType j=0;j<Size;++j) {
    result[j] = 0;
    for(SizeType i=0;i<Size;++i)
      result[j] += M(i,j) * V[i];
  };
  return result;
};
#endif


/*******************************************************************************
                         Comparison Operators
*******************************************************************************/

/**
 * Equality Comparison operator for general matrices, component-wise.
 * \test PASSED
 */
template < typename Matrix1, typename Matrix2 >
typename boost::enable_if_c< is_readable_matrix< Matrix1 >::value && is_readable_matrix< Matrix2 >::value, bool >::type
  operator==( const Matrix1& M1, const Matrix2& M2 ) {
  if( ( M1.get_row_count() != M2.get_row_count() ) || ( M1.get_col_count() != M2.get_col_count() ) )
    return false;
  typedef typename mat_traits< Matrix1 >::size_type SizeType;
  for( SizeType j = 0; j < M1.get_col_count(); ++j )
    for( SizeType i = 0; i < M1.get_row_count(); ++i )
      if( M1( i, j ) != M2( i, j ) )
        return false;
  return true;
};

/**
 * Inequality Comparison operator for general matrices, component-wise.
 * \test PASSED
 */
template < typename Matrix1, typename Matrix2 >
typename boost::enable_if_c< is_readable_matrix< Matrix1 >::value && is_readable_matrix< Matrix2 >::value, bool >::type
  operator!=( const Matrix1& M1, const Matrix2& M2 ) {
  if( ( M1.get_row_count() != M2.get_row_count() ) || ( M1.get_col_count() != M2.get_col_count() ) )
    return true;
  typedef typename mat_traits< Matrix1 >::size_type SizeType;
  for( SizeType j = 0; j < M1.get_col_count(); ++j )
    for( SizeType i = 0; i < M1.get_row_count(); ++i )
      if( M1( i, j ) != M2( i, j ) )
        return true;
  return false;
};


#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE


extern template std::ostream&
  operator<<( std::ostream& out_stream,
              const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M );
extern template std::ostream&
  operator<<( std::ostream& out_stream, const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M );
extern template std::ostream& operator<<( std::ostream& out_stream,
                                          const mat< double, mat_structure::square, mat_alignment::column_major >& M );
extern template std::ostream& operator<<( std::ostream& out_stream,
                                          const mat< double, mat_structure::square, mat_alignment::row_major >& M );
extern template std::ostream& operator<<( std::ostream& out_stream, const mat< double, mat_structure::symmetric >& M );
extern template std::ostream& operator<<( std::ostream& out_stream,
                                          const mat< double, mat_structure::skew_symmetric >& M );
extern template std::ostream& operator<<( std::ostream& out_stream, const mat< double, mat_structure::nil >& M );
extern template std::ostream& operator<<( std::ostream& out_stream, const mat< double, mat_structure::identity >& M );
extern template std::ostream& operator<<( std::ostream& out_stream, const mat< double, mat_structure::scalar >& M );
extern template std::ostream& operator<<( std::ostream& out_stream, const mat< double, mat_structure::diagonal >& M );
extern template std::ostream& operator<<( std::ostream& out_stream,
                                          const mat< double, mat_structure::permutation >& M );

extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< double, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< double, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< double, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< double, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M2 );


// rect and square matrix.
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< double, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< double, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< double, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::square, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< double, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::square, mat_alignment::row_major >& M2 );

// NOTE: These specializations appear to be ambiguous...
// extern template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::column_major>,
// mat<double,mat_structure::symmetric> >::type operator *<
// mat<double,mat_structure::rectangular,mat_alignment::column_major> >(const
// mat<double,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<double,mat_structure::symmetric>&
// M2);
// extern template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::row_major>,
// mat<double,mat_structure::symmetric> >::type operator *<
// mat<double,mat_structure::rectangular,mat_alignment::row_major> >(const
// mat<double,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<double,mat_structure::symmetric>& M2);

// extern template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::column_major>,
// mat<double,mat_structure::skew_symmetric> >::type operator *<
// mat<double,mat_structure::rectangular,mat_alignment::column_major> >(const
// mat<double,mat_structure::rectangular,mat_alignment::column_major>& M1, const
// mat<double,mat_structure::skew_symmetric>& M2);
// extern template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::row_major>,
// mat<double,mat_structure::skew_symmetric> >::type operator *<
// mat<double,mat_structure::rectangular,mat_alignment::row_major> >(const
// mat<double,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<double,mat_structure::skew_symmetric>&
// M2);


// square and rect matrix.
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::column_major >,
                                    mat< double, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::column_major >,
                                    mat< double, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::row_major >,
                                    mat< double, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::row_major >,
                                    mat< double, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M2 );


// square and square matrix.
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::column_major >,
                                    mat< double, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::column_major >,
                                    mat< double, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::square, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::row_major >,
                                    mat< double, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::row_major >,
                                    mat< double, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::square, mat_alignment::row_major >& M2 );

extern template mat_product_result< mat< double, mat_structure::symmetric >,
                                    mat< double, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::symmetric >& M1,
             const mat< double, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::symmetric >,
                                    mat< double, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::symmetric >& M1,
             const mat< double, mat_structure::square, mat_alignment::row_major >& M2 );

extern template mat_product_result< mat< double, mat_structure::skew_symmetric >,
                                    mat< double, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::skew_symmetric >& M1,
             const mat< double, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::skew_symmetric >,
                                    mat< double, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::skew_symmetric >& M1,
             const mat< double, mat_structure::square, mat_alignment::row_major >& M2 );

extern template mat_product_result< mat< double, mat_structure::symmetric >,
                                    mat< double, mat_structure::skew_symmetric > >::type
  operator*( const mat< double, mat_structure::symmetric >& M1,
             const mat< double, mat_structure::skew_symmetric >& M2 );
extern template mat_product_result< mat< double, mat_structure::skew_symmetric >,
                                    mat< double, mat_structure::symmetric > >::type
  operator*( const mat< double, mat_structure::skew_symmetric >& M1,
             const mat< double, mat_structure::symmetric >& M2 );

extern template mat_product_result< mat< double, mat_structure::symmetric >,
                                    mat< double, mat_structure::symmetric > >::type
  operator*( const mat< double, mat_structure::symmetric >& M1, const mat< double, mat_structure::symmetric >& M2 );
extern template mat_product_result< mat< double, mat_structure::skew_symmetric >,
                                    mat< double, mat_structure::skew_symmetric > >::type
  operator*( const mat< double, mat_structure::skew_symmetric >& M1,
             const mat< double, mat_structure::skew_symmetric >& M2 );


// dense and diagonal matrix.
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< double, mat_structure::diagonal > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::diagonal >& M2 );
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< double, mat_structure::diagonal > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::diagonal >& M2 );
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::column_major >,
                                    mat< double, mat_structure::diagonal > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::diagonal >& M2 );
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::row_major >,
                                    mat< double, mat_structure::diagonal > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::diagonal >& M2 );
extern template mat_product_result< mat< double, mat_structure::symmetric >,
                                    mat< double, mat_structure::diagonal > >::type
  operator*( const mat< double, mat_structure::symmetric >& M1, const mat< double, mat_structure::diagonal >& M2 );
extern template mat_product_result< mat< double, mat_structure::skew_symmetric >,
                                    mat< double, mat_structure::diagonal > >::type
  operator*( const mat< double, mat_structure::skew_symmetric >& M1, const mat< double, mat_structure::diagonal >& M2 );

// diagonal and dense matrix.
extern template mat_product_result< mat< double, mat_structure::diagonal >,
                                    mat< double, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::diagonal >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::diagonal >,
                                    mat< double, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::diagonal >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::diagonal >,
                                    mat< double, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::diagonal >& M1,
             const mat< double, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::diagonal >,
                                    mat< double, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::diagonal >& M1,
             const mat< double, mat_structure::square, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::diagonal >,
                                    mat< double, mat_structure::symmetric > >::type
  operator*( const mat< double, mat_structure::diagonal >& M1, const mat< double, mat_structure::symmetric >& M2 );
extern template mat_product_result< mat< double, mat_structure::diagonal >,
                                    mat< double, mat_structure::skew_symmetric > >::type
  operator*( const mat< double, mat_structure::diagonal >& M1, const mat< double, mat_structure::skew_symmetric >& M2 );

// diagonal and diagonal matrix.
extern template mat_product_result< mat< double, mat_structure::diagonal >,
                                    mat< double, mat_structure::diagonal > >::type
  operator*( const mat< double, mat_structure::diagonal >& M1, const mat< double, mat_structure::diagonal >& M2 );


// dense and nil matrix.
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< double, mat_structure::nil > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::nil >& M2 );
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< double, mat_structure::nil > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::nil >& M2 );
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::column_major >,
                                    mat< double, mat_structure::nil > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::nil >& M2 );
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::row_major >,
                                    mat< double, mat_structure::nil > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::nil >& M2 );
extern template mat_product_result< mat< double, mat_structure::symmetric >, mat< double, mat_structure::nil > >::type
  operator*( const mat< double, mat_structure::symmetric >& M1, const mat< double, mat_structure::nil >& M2 );
extern template mat_product_result< mat< double, mat_structure::skew_symmetric >,
                                    mat< double, mat_structure::nil > >::type
  operator*( const mat< double, mat_structure::skew_symmetric >& M1, const mat< double, mat_structure::nil >& M2 );

// nil and dense matrix.
extern template mat_product_result< mat< double, mat_structure::nil >,
                                    mat< double, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::nil >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::nil >,
                                    mat< double, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::nil >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::nil >,
                                    mat< double, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::nil >& M1,
             const mat< double, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::nil >,
                                    mat< double, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::nil >& M1,
             const mat< double, mat_structure::square, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::nil >, mat< double, mat_structure::symmetric > >::type
  operator*( const mat< double, mat_structure::nil >& M1, const mat< double, mat_structure::symmetric >& M2 );
extern template mat_product_result< mat< double, mat_structure::nil >,
                                    mat< double, mat_structure::skew_symmetric > >::type
  operator*( const mat< double, mat_structure::nil >& M1, const mat< double, mat_structure::skew_symmetric >& M2 );

// nil and nil matrix.
extern template mat_product_result< mat< double, mat_structure::nil >, mat< double, mat_structure::nil > >::type
  operator*( const mat< double, mat_structure::nil >& M1, const mat< double, mat_structure::nil >& M2 );


// dense and scalar matrix.
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< double, mat_structure::scalar > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::scalar >& M2 );
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< double, mat_structure::scalar > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::scalar >& M2 );
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::column_major >,
                                    mat< double, mat_structure::scalar > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::scalar >& M2 );
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::row_major >,
                                    mat< double, mat_structure::scalar > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::scalar >& M2 );
extern template mat_product_result< mat< double, mat_structure::symmetric >,
                                    mat< double, mat_structure::scalar > >::type
  operator*( const mat< double, mat_structure::symmetric >& M1, const mat< double, mat_structure::scalar >& M2 );
extern template mat_product_result< mat< double, mat_structure::skew_symmetric >,
                                    mat< double, mat_structure::scalar > >::type
  operator*( const mat< double, mat_structure::skew_symmetric >& M1, const mat< double, mat_structure::scalar >& M2 );

// scalar and dense matrix.
extern template mat_product_result< mat< double, mat_structure::scalar >,
                                    mat< double, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::scalar >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::scalar >,
                                    mat< double, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::scalar >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::scalar >,
                                    mat< double, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::scalar >& M1,
             const mat< double, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::scalar >,
                                    mat< double, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::scalar >& M1,
             const mat< double, mat_structure::square, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::scalar >,
                                    mat< double, mat_structure::symmetric > >::type
  operator*( const mat< double, mat_structure::scalar >& M1, const mat< double, mat_structure::symmetric >& M2 );
extern template mat_product_result< mat< double, mat_structure::scalar >,
                                    mat< double, mat_structure::skew_symmetric > >::type
  operator*( const mat< double, mat_structure::scalar >& M1, const mat< double, mat_structure::skew_symmetric >& M2 );

// scalar and scalar matrix.
extern template mat_product_result< mat< double, mat_structure::scalar >, mat< double, mat_structure::scalar > >::type
  operator*( const mat< double, mat_structure::scalar >& M1, const mat< double, mat_structure::scalar >& M2 );


// dense and identity matrix.
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< double, mat_structure::identity > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::identity >& M2 );
extern template mat_product_result< mat< double, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< double, mat_structure::identity > >::type
  operator*( const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::identity >& M2 );
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::column_major >,
                                    mat< double, mat_structure::identity > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< double, mat_structure::identity >& M2 );
extern template mat_product_result< mat< double, mat_structure::square, mat_alignment::row_major >,
                                    mat< double, mat_structure::identity > >::type
  operator*( const mat< double, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< double, mat_structure::identity >& M2 );
extern template mat_product_result< mat< double, mat_structure::symmetric >,
                                    mat< double, mat_structure::identity > >::type
  operator*( const mat< double, mat_structure::symmetric >& M1, const mat< double, mat_structure::identity >& M2 );
extern template mat_product_result< mat< double, mat_structure::skew_symmetric >,
                                    mat< double, mat_structure::identity > >::type
  operator*( const mat< double, mat_structure::skew_symmetric >& M1, const mat< double, mat_structure::identity >& M2 );

// identity and dense matrix.
extern template mat_product_result< mat< double, mat_structure::identity >,
                                    mat< double, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::identity >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::identity >,
                                    mat< double, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::identity >& M1,
             const mat< double, mat_structure::rectangular, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::identity >,
                                    mat< double, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< double, mat_structure::identity >& M1,
             const mat< double, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::identity >,
                                    mat< double, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< double, mat_structure::identity >& M1,
             const mat< double, mat_structure::square, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< double, mat_structure::identity >,
                                    mat< double, mat_structure::symmetric > >::type
  operator*( const mat< double, mat_structure::identity >& M1, const mat< double, mat_structure::symmetric >& M2 );
extern template mat_product_result< mat< double, mat_structure::identity >,
                                    mat< double, mat_structure::skew_symmetric > >::type
  operator*( const mat< double, mat_structure::identity >& M1, const mat< double, mat_structure::skew_symmetric >& M2 );

// identity and identity matrix.
extern template mat_product_result< mat< double, mat_structure::identity >,
                                    mat< double, mat_structure::identity > >::type
  operator*( const mat< double, mat_structure::identity >& M1, const mat< double, mat_structure::identity >& M2 );


extern template std::ostream&
  operator<<( std::ostream& out_stream,
              const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M );
extern template std::ostream& operator<<( std::ostream& out_stream,
                                          const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M );
extern template std::ostream& operator<<( std::ostream& out_stream,
                                          const mat< float, mat_structure::square, mat_alignment::column_major >& M );
extern template std::ostream& operator<<( std::ostream& out_stream,
                                          const mat< float, mat_structure::square, mat_alignment::row_major >& M );
extern template std::ostream& operator<<( std::ostream& out_stream, const mat< float, mat_structure::symmetric >& M );
extern template std::ostream& operator<<( std::ostream& out_stream,
                                          const mat< float, mat_structure::skew_symmetric >& M );
extern template std::ostream& operator<<( std::ostream& out_stream, const mat< float, mat_structure::nil >& M );
extern template std::ostream& operator<<( std::ostream& out_stream, const mat< float, mat_structure::identity >& M );
extern template std::ostream& operator<<( std::ostream& out_stream, const mat< float, mat_structure::scalar >& M );
extern template std::ostream& operator<<( std::ostream& out_stream, const mat< float, mat_structure::diagonal >& M );
extern template std::ostream& operator<<( std::ostream& out_stream, const mat< float, mat_structure::permutation >& M );

extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< float, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< float, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< float, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< float, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M2 );


// rect and square matrix.
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< float, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< float, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< float, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::square, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< float, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::square, mat_alignment::row_major >& M2 );

// extern template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::column_major>,
// mat<float,mat_structure::symmetric> >::type operator *(const
// mat<float,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<float,mat_structure::symmetric>&
// M2);
// extern template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::row_major>,
// mat<float,mat_structure::symmetric> >::type operator *(const
// mat<float,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<float,mat_structure::symmetric>& M2);

// extern template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::column_major>,
// mat<float,mat_structure::skew_symmetric> >::type operator *(const
// mat<float,mat_structure::rectangular,mat_alignment::column_major>& M1, const
// mat<float,mat_structure::skew_symmetric>& M2);
// extern template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::row_major>,
// mat<float,mat_structure::skew_symmetric> >::type operator *(const
// mat<float,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<float,mat_structure::skew_symmetric>&
// M2);


// square and rect matrix.
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::column_major >,
                                    mat< float, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::column_major >,
                                    mat< float, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::row_major >,
                                    mat< float, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::row_major >,
                                    mat< float, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M2 );


// square and square matrix.
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::column_major >,
                                    mat< float, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::column_major >,
                                    mat< float, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::square, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::row_major >,
                                    mat< float, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::row_major >,
                                    mat< float, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::square, mat_alignment::row_major >& M2 );

extern template mat_product_result< mat< float, mat_structure::symmetric >,
                                    mat< float, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::symmetric >& M1,
             const mat< float, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::symmetric >,
                                    mat< float, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::symmetric >& M1,
             const mat< float, mat_structure::square, mat_alignment::row_major >& M2 );

extern template mat_product_result< mat< float, mat_structure::skew_symmetric >,
                                    mat< float, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::skew_symmetric >& M1,
             const mat< float, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::skew_symmetric >,
                                    mat< float, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::skew_symmetric >& M1,
             const mat< float, mat_structure::square, mat_alignment::row_major >& M2 );

extern template mat_product_result< mat< float, mat_structure::symmetric >,
                                    mat< float, mat_structure::skew_symmetric > >::type
  operator*( const mat< float, mat_structure::symmetric >& M1, const mat< float, mat_structure::skew_symmetric >& M2 );
extern template mat_product_result< mat< float, mat_structure::skew_symmetric >,
                                    mat< float, mat_structure::symmetric > >::type
  operator*( const mat< float, mat_structure::skew_symmetric >& M1, const mat< float, mat_structure::symmetric >& M2 );

extern template mat_product_result< mat< float, mat_structure::symmetric >,
                                    mat< float, mat_structure::symmetric > >::type
  operator*( const mat< float, mat_structure::symmetric >& M1, const mat< float, mat_structure::symmetric >& M2 );
extern template mat_product_result< mat< float, mat_structure::skew_symmetric >,
                                    mat< float, mat_structure::skew_symmetric > >::type
  operator*( const mat< float, mat_structure::skew_symmetric >& M1,
             const mat< float, mat_structure::skew_symmetric >& M2 );


// dense and diagonal matrix.
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< float, mat_structure::diagonal > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::diagonal >& M2 );
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< float, mat_structure::diagonal > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::diagonal >& M2 );
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::column_major >,
                                    mat< float, mat_structure::diagonal > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::diagonal >& M2 );
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::row_major >,
                                    mat< float, mat_structure::diagonal > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::diagonal >& M2 );
extern template mat_product_result< mat< float, mat_structure::symmetric >,
                                    mat< float, mat_structure::diagonal > >::type
  operator*( const mat< float, mat_structure::symmetric >& M1, const mat< float, mat_structure::diagonal >& M2 );
extern template mat_product_result< mat< float, mat_structure::skew_symmetric >,
                                    mat< float, mat_structure::diagonal > >::type
  operator*( const mat< float, mat_structure::skew_symmetric >& M1, const mat< float, mat_structure::diagonal >& M2 );

// diagonal and dense matrix.
extern template mat_product_result< mat< float, mat_structure::diagonal >,
                                    mat< float, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::diagonal >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::diagonal >,
                                    mat< float, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::diagonal >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::diagonal >,
                                    mat< float, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::diagonal >& M1,
             const mat< float, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::diagonal >,
                                    mat< float, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::diagonal >& M1,
             const mat< float, mat_structure::square, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::diagonal >,
                                    mat< float, mat_structure::symmetric > >::type
  operator*( const mat< float, mat_structure::diagonal >& M1, const mat< float, mat_structure::symmetric >& M2 );
extern template mat_product_result< mat< float, mat_structure::diagonal >,
                                    mat< float, mat_structure::skew_symmetric > >::type
  operator*( const mat< float, mat_structure::diagonal >& M1, const mat< float, mat_structure::skew_symmetric >& M2 );

// diagonal and diagonal matrix.
extern template mat_product_result< mat< float, mat_structure::diagonal >, mat< float, mat_structure::diagonal > >::type
  operator*( const mat< float, mat_structure::diagonal >& M1, const mat< float, mat_structure::diagonal >& M2 );


// dense and nil matrix.
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< float, mat_structure::nil > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::nil >& M2 );
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< float, mat_structure::nil > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::nil >& M2 );
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::column_major >,
                                    mat< float, mat_structure::nil > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::nil >& M2 );
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::row_major >,
                                    mat< float, mat_structure::nil > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::nil >& M2 );
extern template mat_product_result< mat< float, mat_structure::symmetric >, mat< float, mat_structure::nil > >::type
  operator*( const mat< float, mat_structure::symmetric >& M1, const mat< float, mat_structure::nil >& M2 );
extern template mat_product_result< mat< float, mat_structure::skew_symmetric >,
                                    mat< float, mat_structure::nil > >::type
  operator*( const mat< float, mat_structure::skew_symmetric >& M1, const mat< float, mat_structure::nil >& M2 );

// nil and dense matrix.
extern template mat_product_result< mat< float, mat_structure::nil >,
                                    mat< float, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::nil >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::nil >,
                                    mat< float, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::nil >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::nil >,
                                    mat< float, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::nil >& M1,
             const mat< float, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::nil >,
                                    mat< float, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::nil >& M1,
             const mat< float, mat_structure::square, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::nil >, mat< float, mat_structure::symmetric > >::type
  operator*( const mat< float, mat_structure::nil >& M1, const mat< float, mat_structure::symmetric >& M2 );
extern template mat_product_result< mat< float, mat_structure::nil >,
                                    mat< float, mat_structure::skew_symmetric > >::type
  operator*( const mat< float, mat_structure::nil >& M1, const mat< float, mat_structure::skew_symmetric >& M2 );

// nil and nil matrix.
extern template mat_product_result< mat< float, mat_structure::nil >, mat< float, mat_structure::nil > >::type
  operator*( const mat< float, mat_structure::nil >& M1, const mat< float, mat_structure::nil >& M2 );


// dense and scalar matrix.
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< float, mat_structure::scalar > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::scalar >& M2 );
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< float, mat_structure::scalar > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::scalar >& M2 );
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::column_major >,
                                    mat< float, mat_structure::scalar > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::scalar >& M2 );
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::row_major >,
                                    mat< float, mat_structure::scalar > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::scalar >& M2 );
extern template mat_product_result< mat< float, mat_structure::symmetric >, mat< float, mat_structure::scalar > >::type
  operator*( const mat< float, mat_structure::symmetric >& M1, const mat< float, mat_structure::scalar >& M2 );
extern template mat_product_result< mat< float, mat_structure::skew_symmetric >,
                                    mat< float, mat_structure::scalar > >::type
  operator*( const mat< float, mat_structure::skew_symmetric >& M1, const mat< float, mat_structure::scalar >& M2 );

// scalar and dense matrix.
extern template mat_product_result< mat< float, mat_structure::scalar >,
                                    mat< float, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::scalar >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::scalar >,
                                    mat< float, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::scalar >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::scalar >,
                                    mat< float, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::scalar >& M1,
             const mat< float, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::scalar >,
                                    mat< float, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::scalar >& M1,
             const mat< float, mat_structure::square, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::scalar >, mat< float, mat_structure::symmetric > >::type
  operator*( const mat< float, mat_structure::scalar >& M1, const mat< float, mat_structure::symmetric >& M2 );
extern template mat_product_result< mat< float, mat_structure::scalar >,
                                    mat< float, mat_structure::skew_symmetric > >::type
  operator*( const mat< float, mat_structure::scalar >& M1, const mat< float, mat_structure::skew_symmetric >& M2 );

// scalar and scalar matrix.
extern template mat_product_result< mat< float, mat_structure::scalar >, mat< float, mat_structure::scalar > >::type
  operator*( const mat< float, mat_structure::scalar >& M1, const mat< float, mat_structure::scalar >& M2 );


// dense and identity matrix.
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::column_major >,
                                    mat< float, mat_structure::identity > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::identity >& M2 );
extern template mat_product_result< mat< float, mat_structure::rectangular, mat_alignment::row_major >,
                                    mat< float, mat_structure::identity > >::type
  operator*( const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::identity >& M2 );
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::column_major >,
                                    mat< float, mat_structure::identity > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::column_major >& M1,
             const mat< float, mat_structure::identity >& M2 );
extern template mat_product_result< mat< float, mat_structure::square, mat_alignment::row_major >,
                                    mat< float, mat_structure::identity > >::type
  operator*( const mat< float, mat_structure::square, mat_alignment::row_major >& M1,
             const mat< float, mat_structure::identity >& M2 );
extern template mat_product_result< mat< float, mat_structure::symmetric >,
                                    mat< float, mat_structure::identity > >::type
  operator*( const mat< float, mat_structure::symmetric >& M1, const mat< float, mat_structure::identity >& M2 );
extern template mat_product_result< mat< float, mat_structure::skew_symmetric >,
                                    mat< float, mat_structure::identity > >::type
  operator*( const mat< float, mat_structure::skew_symmetric >& M1, const mat< float, mat_structure::identity >& M2 );

// identity and dense matrix.
extern template mat_product_result< mat< float, mat_structure::identity >,
                                    mat< float, mat_structure::rectangular, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::identity >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::identity >,
                                    mat< float, mat_structure::rectangular, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::identity >& M1,
             const mat< float, mat_structure::rectangular, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::identity >,
                                    mat< float, mat_structure::square, mat_alignment::column_major > >::type
  operator*( const mat< float, mat_structure::identity >& M1,
             const mat< float, mat_structure::square, mat_alignment::column_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::identity >,
                                    mat< float, mat_structure::square, mat_alignment::row_major > >::type
  operator*( const mat< float, mat_structure::identity >& M1,
             const mat< float, mat_structure::square, mat_alignment::row_major >& M2 );
extern template mat_product_result< mat< float, mat_structure::identity >,
                                    mat< float, mat_structure::symmetric > >::type
  operator*( const mat< float, mat_structure::identity >& M1, const mat< float, mat_structure::symmetric >& M2 );
extern template mat_product_result< mat< float, mat_structure::identity >,
                                    mat< float, mat_structure::skew_symmetric > >::type
  operator*( const mat< float, mat_structure::identity >& M1, const mat< float, mat_structure::skew_symmetric >& M2 );

// identity and identity matrix.
extern template mat_product_result< mat< float, mat_structure::identity >, mat< float, mat_structure::identity > >::type
  operator*( const mat< float, mat_structure::identity >& M1, const mat< float, mat_structure::identity >& M2 );


#endif
};


#endif

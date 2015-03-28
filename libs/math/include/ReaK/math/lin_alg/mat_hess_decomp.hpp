/**
 * \file mat_hess_decomp.hpp
 *
 * This library provides methods to perform the Upper-Hessenberg decomposition of a matrix as
 * well as the Hessenberg-Triangular reduction of two matrices. These algorithms form the first
 * step in many algorithms such as the QR algorithm, the QZ algorithm, the Real-Schur decomposition,
 * etc. These algorithms were implemented as described in Golub and van Loan's classic book.
 *
 * This library simply provides various versions of the same algorithm (same underlying implementation,
 * with different interfaces). Versions differ based mostly on the type of matrices fed to the function
 * overloads and whether the orthogonal matrices that achieve the decompositions are required or not
 * (since explicitly forming those matrices can be expensive and should be avoided if it's not needed).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2011
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

#ifndef REAK_MAT_HESS_DECOMP_HPP
#define REAK_MAT_HESS_DECOMP_HPP

#include "mat_alg.hpp"
#include "mat_num_exceptions.hpp"

#include "mat_householder.hpp"
#include "mat_givens_rot.hpp"

#include "mat_qr_decomp.hpp"

namespace ReaK {


/*************************************************************************
                          Hessenberg Decompositions
*************************************************************************/

namespace detail {


template < typename Matrix1, typename Matrix2 >
void decompose_Hess_impl( Matrix1& H, Matrix2* Q, typename mat_traits< Matrix1 >::value_type absNumTol ) {

  typedef typename mat_traits< Matrix1 >::value_type ValueType;
  typedef typename mat_traits< Matrix1 >::size_type SizeType;
  SizeType N = H.get_row_count();
  householder_matrix< vect_n< ValueType > > hhm;

  for( SizeType i = 0; i + 2 < N; ++i ) {

    hhm.set( mat_row_slice< Matrix1 >( H, i, i + 1, N - i - 1 ), absNumTol );

    mat_sub_block< Matrix1 > subH1( H, N - i - 1, N - i, i + 1, i );
    householder_prod( hhm, subH1 ); // P * H

    mat_sub_block< Matrix1 > subH2( H, N, N - i - 1, 0, i + 1 );
    householder_prod( subH2, hhm ); // H * P

    if( Q ) {
      mat_sub_block< Matrix2 > subQ( *Q, N, N - i - 1, 0, i + 1 );
      householder_prod( subQ, hhm ); // Q_prev * P
    };
  };
};

/* tested and working. Golub and vanLoan Alg.-8.3.1 */
template < typename Matrix1, typename Matrix2 >
void decompose_TriDiag_impl( Matrix1& T, Matrix2* Q, typename mat_traits< Matrix1 >::value_type absNumTol ) {
  typedef typename mat_traits< Matrix1 >::value_type ValueType;
  typedef typename mat_traits< Matrix1 >::size_type SizeType;
  using std::sqrt;
  SizeType N = T.get_row_count();
  householder_matrix< vect_n< ValueType > > hhm;

  vect_n< ValueType > p( N );
  vect_n< ValueType > w( N );

  for( SizeType i = 0; i + 2 < N; ++i ) {

    hhm.set( mat_row_slice< Matrix1 >( T, i, i + 1, N - i - 1 ), absNumTol );

    mat_sub_block< Matrix1 > subT1( T, N - i, N - i - 1, i, i + 1 );
    householder_prod( subT1, hhm ); // Q_prev * P

    mat_sub_block< Matrix1 > subT2( T, N - i - 1, N - i, i + 1, i );
    householder_prod( hhm, subT2 ); // Q_prev * P

    if( Q ) {
      // mat_sub_block<Matrix2> subQ(*Q,N - i - 1,N,i + 1,0);
      // householder_prod(hhm,subQ); // Q_prev * P
      mat_sub_block< Matrix2 > subQ( *Q, N, N - i - 1, 0, i + 1 );
      householder_prod( subQ, hhm ); // Q_prev * P
    };
  };
};


template < typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4 >
void reduce_HessTri_offset_impl( Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
                                 typename mat_traits< Matrix1 >::size_type row_offset,
                                 typename mat_traits< Matrix1 >::value_type absNumTol ) {
  typedef typename mat_traits< Matrix1 >::value_type ValueType;
  typedef typename mat_traits< Matrix1 >::size_type SizeType;
  using std::fabs;

  SizeType N = A.get_row_count() - row_offset;

  givens_rot_matrix< ValueType > G;

  //   std::cout << "Hess-Tri: Before QR: A = " << A << std::endl;
  //   std::cout << "Hess-Tri: Before QR: B = " << B << std::endl;

  householder_matrix< vect_n< ValueType > > hhm;

  for( SizeType i = 0; i + 1 < N; ++i ) {

    hhm.set( mat_row_slice< Matrix2 >( B, i, row_offset + i, N - i ), absNumTol );

    mat_sub_block< Matrix2 > subB( B, N - i, B.get_col_count() - i, row_offset + i, i );
    householder_prod( hhm, subB ); // P * R

    mat_sub_block< Matrix1 > subA( A, N - i, A.get_col_count(), row_offset + i, 0 );
    householder_prod( hhm, subA ); // P * R

    if( Q ) {
      mat_sub_block< Matrix3 > subQ( *Q, Q->get_row_count(), N - i, 0, i );
      householder_prod( subQ, hhm ); // Q_prev * P
    };
  };

  //   std::cout << "Hess-Tri: After QR: A = " << A << std::endl;
  //   std::cout << "Hess-Tri: After QR: B = " << B << std::endl;

  for( SizeType j = 0; j + 2 < N; ++j ) {
    for( SizeType i = N - 1; i > j + 1; --i ) {
      if( fabs( A( i, j ) ) < absNumTol )
        continue;

      G.set( A( row_offset + i - 1, j ), A( row_offset + i, j ) );

      mat_sub_block< Matrix1 > subA1( A, 2, A.get_col_count() - j, row_offset + i - 1, j );
      givens_rot_prod( G, subA1 ); // G * A

      mat_sub_block< Matrix2 > subB1( B, 2, B.get_col_count() - i + 1, row_offset + i - 1, i - 1 );
      givens_rot_prod( G, subB1 ); // G * B

      if( Q ) {
        mat_sub_block< Matrix3 > subQ( *Q, Q->get_row_count(), 2, 0, i - 1 );
        givens_rot_prod( subQ, transpose( G ) ); // Q_prev * G^T
      };


      //       std::cout << "Hess-Tri: Step 1: A = " << A << std::endl;
      //       std::cout << "Hess-Tri: Step 1: B = " << B << std::endl;


      G.set( -B( row_offset + i, i ), B( row_offset + i, i - 1 ) );
      G = transpose( G );

      mat_sub_block< Matrix2 > subB2( B, i + 1, 2, row_offset, i - 1 );
      givens_rot_prod( subB2, G ); // B * G^T

      mat_sub_block< Matrix1 > subA2( A, N, 2, row_offset, i - 1 );
      givens_rot_prod( subA2, G ); // A * G^T

      if( Z ) {
        mat_sub_block< Matrix4 > subZ( *Z, Z->get_row_count(), 2, 0, i - 1 );
        givens_rot_prod( subZ, G ); // Q_prev * G^T
      };

      //       std::cout << "Hess-Tri: Step 2: A = " << A << std::endl;
      //       std::cout << "Hess-Tri: Step 2: B = " << B << std::endl;
    };
  };
};

template < typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4 >
void reduce_HessTri_impl( Matrix1& A, Matrix2& B, Matrix3* Q, Matrix4* Z,
                          typename mat_traits< Matrix1 >::value_type absNumTol ) {
  reduce_HessTri_offset_impl( A, B, Q, Z, 0, absNumTol );
};
};


/**
 * Performs the Upper-Hessenberg decomposition on a matrix, using the Householder method.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A fully-writable matrix type.
 * \tparam Matrix3 A fully-writable matrix type.
 * \param A square matrix with row-count == column-count.
 * \param Q holds as output, the unitary square matrix Q.
 * \param H holds as output, the upper-hessenberg matrix R in A = Q H Q^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template < typename Matrix1, typename Matrix2, typename Matrix3 >
typename boost::enable_if_c< is_readable_matrix< Matrix1 >::value && is_fully_writable_matrix< Matrix2 >::value
                             && is_fully_writable_matrix< Matrix3 >::value,
                             void >::type
  decompose_Hess( const Matrix1& A, Matrix2& Q, Matrix3& H, typename mat_traits< Matrix1 >::value_type NumTol = 1E-8 ) {
  if( A.get_row_count() != A.get_col_count() )
    throw std::range_error( "Upper-Hessenberg decomposition is only possible on a square matrix!" );

  Q = mat< typename mat_traits< Matrix3 >::value_type, mat_structure::identity >( A.get_row_count() );
  H = A;
  detail::decompose_Hess_impl( H, &Q, NumTol );
};


/**
 * Performs the Upper-Hessenberg decomposition on a matrix, using the Householder method.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A writable matrix type.
 * \tparam Matrix3 A writable matrix type.
 * \param A square matrix with row-count == column-count.
 * \param Q holds as output, the unitary square matrix Q.
 * \param H holds as output, the upper-hessenberg matrix R in A = Q H Q^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template < typename Matrix1, typename Matrix2, typename Matrix3 >
typename boost::enable_if_c< is_readable_matrix< Matrix1 >::value && is_writable_matrix< Matrix2 >::value
                             && is_writable_matrix< Matrix3 >::value && !is_fully_writable_matrix< Matrix2 >::value
                             && !is_fully_writable_matrix< Matrix3 >::value,
                             void >::type
  decompose_Hess( const Matrix1& A, Matrix2& Q, Matrix3& H, typename mat_traits< Matrix1 >::value_type NumTol = 1E-8 ) {
  if( A.get_row_count() != A.get_col_count() )
    throw std::range_error( "Upper-Hessenberg decomposition is only possible on a square matrix!" );

  mat< typename mat_traits< Matrix3 >::value_type, mat_structure::square > Qtmp(
    mat< typename mat_traits< Matrix2 >::value_type, mat_structure::identity >( A.get_row_count() ) );
  mat< typename mat_traits< Matrix2 >::value_type, mat_structure::square > Htmp( A );
  detail::decompose_Hess_impl( Htmp, &Qtmp, NumTol );
  H = Htmp;
  Q = Qtmp;
};


/**
 * Performs the Upper-Hessenberg decomposition on a matrix, using the Householder method.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A fully-writable matrix type.
 * \param A square matrix with row-count == column-count.
 * \param H holds as output, the upper-hessenberg matrix R in A = Q H Q^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template < typename Matrix1, typename Matrix2 >
typename boost::enable_if_c< is_readable_matrix< Matrix1 >::value && is_fully_writable_matrix< Matrix2 >::value,
                             void >::type
  decompose_Hess( const Matrix1& A, Matrix2& H, typename mat_traits< Matrix1 >::value_type NumTol = 1E-8 ) {
  if( A.get_row_count() != A.get_col_count() )
    throw std::range_error( "Upper-Hessenberg decomposition is only possible on a square matrix!" );

  H = A;
  detail::decompose_Hess_impl( H, static_cast< Matrix2* >( NULL ), NumTol );
};


/**
 * Performs the Upper-Hessenberg decomposition on a matrix, using the Householder method.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A writable matrix type.
 * \param A square matrix with row-count == column-count.
 * \param H holds as output, the upper-hessenberg matrix R in A = Q H Q^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template < typename Matrix1, typename Matrix2 >
typename boost::enable_if_c< is_readable_matrix< Matrix1 >::value && is_writable_matrix< Matrix2 >::value
                             && !is_fully_writable_matrix< Matrix2 >::value,
                             void >::type
  decompose_Hess( const Matrix1& A, Matrix2& H, typename mat_traits< Matrix1 >::value_type NumTol = 1E-8 ) {
  if( A.get_row_count() != A.get_col_count() )
    throw std::range_error( "Upper-Hessenberg decomposition is only possible on a square matrix!" );

  mat< typename mat_traits< Matrix2 >::value_type, mat_structure::square > Htmp( A );
  detail::decompose_Hess_impl( Htmp, static_cast< Matrix2* >( NULL ), NumTol );
  H = Htmp;
};


/**
 * Performs the Hessenberg-Triangular reduction on a matrices A and B, using the Givens rotations.
 * Given two square matrices, A and B, this algorithm produces matrices H and R, where H is
 * upper-Hessenberg and R is upper-triangular, and are "similar" to matrices A and B through
 * the transformation H = Q^T * A * Z and R = Q^T * B * Z, where both Q and Z are orthogonal.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A readable matrix type.
 * \tparam Matrix3 A fully-writable matrix type.
 * \tparam Matrix4 A fully-writable matrix type.
 * \tparam Matrix5 A fully-writable matrix type.
 * \tparam Matrix6 A fully-writable matrix type.
 * \param A square matrix with row-count == column-count.
 * \param B square matrix with row-count == column-count.
 * \param H holds as output, the upper-hessenberg matrix H in H = Q^T A Z.
 * \param R holds as output, the upper-triangular matrix R in R = Q^T B Z.
 * \param Q holds as output, the unitary square matrix Q.
 * \param Z holds as output, the unitary square matrix Z.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template < typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4, typename Matrix5, typename Matrix6 >
typename boost::enable_if_c< is_readable_matrix< Matrix1 >::value && is_readable_matrix< Matrix2 >::value
                             && is_fully_writable_matrix< Matrix3 >::value && is_fully_writable_matrix< Matrix4 >::value
                             && is_fully_writable_matrix< Matrix5 >::value
                             && is_fully_writable_matrix< Matrix6 >::value,
                             void >::type
  reduce_HessTri( const Matrix1& A, const Matrix2& B, Matrix3& H, Matrix4& R, Matrix5& Q, Matrix6& Z,
                  typename mat_traits< Matrix1 >::value_type NumTol = 1E-8 ) {
  if( ( A.get_row_count() != A.get_col_count() ) || ( B.get_row_count() != B.get_col_count() ) )
    throw std::range_error( "Hessenberg-Triangular reduction is only possible on square matrices!" );

  H = A;
  R = B;
  Q = mat< typename mat_traits< Matrix5 >::value_type, mat_structure::identity >( A.get_row_count() );
  Z = mat< typename mat_traits< Matrix6 >::value_type, mat_structure::identity >( A.get_row_count() );
  detail::reduce_HessTri_offset_impl( H, R, &Q, &Z, 0, NumTol );
};


/**
 * Performs the Hessenberg-Triangular reduction on a matrices A and B, using the Givens rotations.
 * Given two square matrices, A and B, this algorithm produces matrices H and R, where H is
 * upper-Hessenberg and R is upper-triangular, and are "similar" to matrices A and B through
 * the transformation H = Q^T * A * Z and R = Q^T * B * Z, where both Q and Z are orthogonal.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A readable matrix type.
 * \tparam Matrix3 A writable matrix type.
 * \tparam Matrix4 A writable matrix type.
 * \tparam Matrix5 A writable matrix type.
 * \tparam Matrix6 A writable matrix type.
 * \param A square matrix with row-count == column-count.
 * \param B square matrix with row-count == column-count.
 * \param H holds as output, the upper-hessenberg matrix H in H = Q^T A Z.
 * \param R holds as output, the upper-triangular matrix R in R = Q^T B Z.
 * \param Q holds as output, the unitary square matrix Q.
 * \param Z holds as output, the unitary square matrix Z.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template < typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4, typename Matrix5, typename Matrix6 >
typename boost::enable_if_c< is_readable_matrix< Matrix1 >::value && is_readable_matrix< Matrix2 >::value
                             && is_writable_matrix< Matrix3 >::value && is_writable_matrix< Matrix4 >::value
                             && is_writable_matrix< Matrix5 >::value && is_writable_matrix< Matrix6 >::value
                             && !is_fully_writable_matrix< Matrix3 >::value
                             && !is_fully_writable_matrix< Matrix4 >::value
                             && !is_fully_writable_matrix< Matrix5 >::value
                             && !is_fully_writable_matrix< Matrix6 >::value,
                             void >::type
  reduce_HessTri( const Matrix1& A, const Matrix2& B, Matrix3& H, Matrix4& R, Matrix5& Q, Matrix6& Z,
                  typename mat_traits< Matrix1 >::value_type NumTol = 1E-8 ) {
  if( ( A.get_row_count() != A.get_col_count() ) || ( B.get_row_count() != B.get_col_count() ) )
    throw std::range_error( "Hessenberg-Triangular reduction is only possible on square matrices!" );

  mat< typename mat_traits< Matrix3 >::value_type, mat_structure::square > Htmp( A );
  mat< typename mat_traits< Matrix4 >::value_type, mat_structure::square > Rtmp( B );
  mat< typename mat_traits< Matrix5 >::value_type, mat_structure::square > Qtmp(
    mat< typename mat_traits< Matrix5 >::value_type, mat_structure::identity >( A.get_row_count() ) );
  mat< typename mat_traits< Matrix6 >::value_type, mat_structure::square > Ztmp(
    mat< typename mat_traits< Matrix6 >::value_type, mat_structure::identity >( A.get_row_count() ) );
  detail::reduce_HessTri_offset_impl( Htmp, Rtmp, &Qtmp, &Ztmp, 0, NumTol );
  H = Htmp;
  R = Rtmp;
  Q = Qtmp;
  Z = Ztmp;
};


/**
 * Performs the Hessenberg-Triangular reduction on a matrices A and B, using the Givens rotations.
 * Given two square matrices, A and B, this algorithm produces matrices H and R, where H is
 * upper-Hessenberg and R is upper-triangular, and are "similar" to matrices A and B through
 * the transformation H = Q^T * A * Z and R = Q^T * B * Z, where both Q and Z are orthogonal.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A readable matrix type.
 * \tparam Matrix3 A fully-writable matrix type.
 * \tparam Matrix4 A fully-writable matrix type.
 * \param A square matrix with row-count == column-count.
 * \param B square matrix with row-count == column-count.
 * \param H holds as output, the upper-hessenberg matrix H in H = Q^T A Z.
 * \param R holds as output, the upper-triangular matrix R in R = Q^T B Z.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template < typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4 >
typename boost::enable_if_c< is_readable_matrix< Matrix1 >::value && is_readable_matrix< Matrix2 >::value
                             && is_fully_writable_matrix< Matrix3 >::value
                             && is_fully_writable_matrix< Matrix4 >::value,
                             void >::type
  reduce_HessTri( const Matrix1& A, const Matrix2& B, Matrix3& H, Matrix4& R,
                  typename mat_traits< Matrix1 >::value_type NumTol = 1E-8 ) {
  if( ( A.get_row_count() != A.get_col_count() ) || ( B.get_row_count() != B.get_col_count() ) )
    throw std::range_error( "Hessenberg-Triangular reduction is only possible on square matrices!" );

  H = A;
  R = B;
  detail::reduce_HessTri_offset_impl( H, R, static_cast< Matrix3* >( NULL ), static_cast< Matrix4* >( NULL ), 0,
                                      NumTol );
};


/**
 * Performs the Hessenberg-Triangular reduction on a matrices A and B, using the Givens rotations.
 * Given two square matrices, A and B, this algorithm produces matrices H and R, where H is
 * upper-Hessenberg and R is upper-triangular, and are "similar" to matrices A and B through
 * the transformation H = Q^T * A * Z and R = Q^T * B * Z, where both Q and Z are orthogonal.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A readable matrix type.
 * \tparam Matrix3 A writable matrix type.
 * \tparam Matrix4 A writable matrix type.
 * \param A square matrix with row-count == column-count.
 * \param B square matrix with row-count == column-count.
 * \param H holds as output, the upper-hessenberg matrix H in H = Q^T A Z.
 * \param R holds as output, the upper-triangular matrix R in R = Q^T B Z.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template < typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4 >
typename boost::enable_if_c< is_readable_matrix< Matrix1 >::value && is_readable_matrix< Matrix2 >::value
                             && is_writable_matrix< Matrix3 >::value && is_writable_matrix< Matrix4 >::value
                             && !is_fully_writable_matrix< Matrix3 >::value
                             && !is_fully_writable_matrix< Matrix4 >::value,
                             void >::type
  reduce_HessTri( const Matrix1& A, const Matrix2& B, Matrix3& H, Matrix4& R,
                  typename mat_traits< Matrix1 >::value_type NumTol = 1E-8 ) {
  if( ( A.get_row_count() != A.get_col_count() ) || ( B.get_row_count() != B.get_col_count() ) )
    throw std::range_error( "Hessenberg-Triangular reduction is only possible on square matrices!" );

  mat< typename mat_traits< Matrix3 >::value_type, mat_structure::square > Htmp( A );
  mat< typename mat_traits< Matrix4 >::value_type, mat_structure::square > Rtmp( B );
  detail::reduce_HessTri_offset_impl( Htmp, Rtmp, static_cast< Matrix3* >( NULL ), static_cast< Matrix4* >( NULL ), 0,
                                      NumTol );
  H = Htmp;
  R = Rtmp;
};


#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

extern template void decompose_Hess( const mat< double, mat_structure::square >& A,
                                     mat< double, mat_structure::square >& Q, mat< double, mat_structure::square >& H,
                                     double NumTol );
extern template void decompose_Hess( const mat< double, mat_structure::rectangular >& A,
                                     mat< double, mat_structure::rectangular >& Q,
                                     mat< double, mat_structure::rectangular >& H, double NumTol );
extern template void decompose_Hess( const mat< double, mat_structure::rectangular >& A,
                                     mat< double, mat_structure::square >& Q,
                                     mat< double, mat_structure::rectangular >& H, double NumTol );

extern template void decompose_Hess( const mat< double, mat_structure::square >& A,
                                     mat< double, mat_structure::square >& H, double NumTol );
extern template void decompose_Hess( const mat< double, mat_structure::rectangular >& A,
                                     mat< double, mat_structure::square >& H, double NumTol );
extern template void decompose_Hess( const mat< double, mat_structure::rectangular >& A,
                                     mat< double, mat_structure::rectangular >& H, double NumTol );

extern template void reduce_HessTri( const mat< double, mat_structure::square >& A,
                                     const mat< double, mat_structure::square >& B,
                                     mat< double, mat_structure::square >& H, mat< double, mat_structure::square >& R,
                                     mat< double, mat_structure::square >& Q, mat< double, mat_structure::square >& Z,
                                     double NumTol );
extern template void reduce_HessTri( const mat< double, mat_structure::rectangular >& A,
                                     const mat< double, mat_structure::rectangular >& B,
                                     mat< double, mat_structure::rectangular >& H,
                                     mat< double, mat_structure::rectangular >& R,
                                     mat< double, mat_structure::square >& Q, mat< double, mat_structure::square >& Z,
                                     double NumTol );
extern template void reduce_HessTri( const mat< double, mat_structure::rectangular >& A,
                                     const mat< double, mat_structure::rectangular >& B,
                                     mat< double, mat_structure::rectangular >& H,
                                     mat< double, mat_structure::rectangular >& R,
                                     mat< double, mat_structure::rectangular >& Q,
                                     mat< double, mat_structure::rectangular >& Z, double NumTol );

extern template void reduce_HessTri( const mat< double, mat_structure::square >& A,
                                     const mat< double, mat_structure::square >& B,
                                     mat< double, mat_structure::square >& H, mat< double, mat_structure::square >& R,
                                     double NumTol );
extern template void reduce_HessTri( const mat< double, mat_structure::rectangular >& A,
                                     const mat< double, mat_structure::rectangular >& B,
                                     mat< double, mat_structure::rectangular >& H,
                                     mat< double, mat_structure::rectangular >& R, double NumTol );


extern template void decompose_Hess( const mat< float, mat_structure::square >& A,
                                     mat< float, mat_structure::square >& Q, mat< float, mat_structure::square >& H,
                                     float NumTol );
extern template void decompose_Hess( const mat< float, mat_structure::rectangular >& A,
                                     mat< float, mat_structure::rectangular >& Q,
                                     mat< float, mat_structure::rectangular >& H, float NumTol );
extern template void decompose_Hess( const mat< float, mat_structure::rectangular >& A,
                                     mat< float, mat_structure::square >& Q,
                                     mat< float, mat_structure::rectangular >& H, float NumTol );

extern template void decompose_Hess( const mat< float, mat_structure::square >& A,
                                     mat< float, mat_structure::square >& H, float NumTol );
extern template void decompose_Hess( const mat< float, mat_structure::rectangular >& A,
                                     mat< float, mat_structure::square >& H, float NumTol );
extern template void decompose_Hess( const mat< float, mat_structure::rectangular >& A,
                                     mat< float, mat_structure::rectangular >& H, float NumTol );

extern template void reduce_HessTri( const mat< float, mat_structure::square >& A,
                                     const mat< float, mat_structure::square >& B,
                                     mat< float, mat_structure::square >& H, mat< float, mat_structure::square >& R,
                                     mat< float, mat_structure::square >& Q, mat< float, mat_structure::square >& Z,
                                     float NumTol );
extern template void
  reduce_HessTri( const mat< float, mat_structure::rectangular >& A, const mat< float, mat_structure::rectangular >& B,
                  mat< float, mat_structure::rectangular >& H, mat< float, mat_structure::rectangular >& R,
                  mat< float, mat_structure::square >& Q, mat< float, mat_structure::square >& Z, float NumTol );
extern template void reduce_HessTri( const mat< float, mat_structure::rectangular >& A,
                                     const mat< float, mat_structure::rectangular >& B,
                                     mat< float, mat_structure::rectangular >& H,
                                     mat< float, mat_structure::rectangular >& R,
                                     mat< float, mat_structure::rectangular >& Q,
                                     mat< float, mat_structure::rectangular >& Z, float NumTol );

extern template void reduce_HessTri( const mat< float, mat_structure::square >& A,
                                     const mat< float, mat_structure::square >& B,
                                     mat< float, mat_structure::square >& H, mat< float, mat_structure::square >& R,
                                     float NumTol );
extern template void reduce_HessTri( const mat< float, mat_structure::rectangular >& A,
                                     const mat< float, mat_structure::rectangular >& B,
                                     mat< float, mat_structure::rectangular >& H,
                                     mat< float, mat_structure::rectangular >& R, float NumTol );


#endif
};

#endif


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

#ifndef MAT_HESS_DECOMP_HPP
#define MAT_HESS_DECOMP_HPP

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




template <typename Matrix1, typename Matrix2>
void decompose_Hess_impl(Matrix1& H, Matrix2* Q, typename mat_traits<Matrix1>::value_type NumTol)
{
  
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  SizeType N = H.get_row_count();
  householder_matrix< vect_n<ValueType> > hhm;
  
  SizeType t = N-2;

  for(SizeType i=0;i<t;++i) {
    
    hhm.set(mat_row_slice<Matrix1>(H,i,i+1,N - i - 1),NumTol);
    
    mat_sub_block<Matrix1> subH1(H,N - i - 1,N - i,i+1,i);
    householder_prod(hhm,subH1); // P * H
    
    mat_sub_block<Matrix1> subH2(H,N,N - i - 1,0,i+1);
    householder_prod(subH2,hhm); // H * P
    
    if(Q) {
      mat_sub_block<Matrix2> subQ(*Q,N,N - i,0,i);
      householder_prod(subQ,hhm); // Q_prev * P
    };

  };
  
};
  



template <typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4>
void reduce_HessTri_impl(Matrix1& A, Matrix2& B, Matrix3* Q, Matrix3* Z, typename mat_traits<Matrix1>::value_type NumTol)
{
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  SizeType N = H.get_row_count();
  
  givens_rot_matrix<ValueType> G;
  
  mat<ValueType,mat_structure::square> Q_init(mat<ValueType,mat_structure::identity>(N));
  decompose_QR_impl(B,&Q_init,NumTol);
  
  if(Q)
    *Q *= Q_init;
  A = transpose(Q_init) * A;
  
  SizeType t = N-2;

  for(SizeType j=0;j<t;++j) {
    for(SizeType i = N-1; i>j+1; --i) {
      G.set(A(i-1,j),A(i,j));
    
      mat_sub_block<Matrix1> subA1(A,2,N - j,i-1,j);
      givens_rot_prod(G,subA1); // G * A
    
      mat_sub_block<Matrix2> subB1(B,2,N - i + 1,i-1,i-1);
      givens_rot_prod(G,subB1); // G * B
    
      if(Q) {
        mat_sub_block<Matrix3> subQ(*Q,N,2,0,i-1);
        givens_rot_prod(subQ,transpose(G)); // Q_prev * G^T
      };
      
      G.set(-B(i,i),B(i,i-1));
      G = transpose(G);
    
      mat_sub_block<Matrix2> subB2(B,i+1,2,0,i-1);
      givens_rot_prod(subB2,G); // B * G^T
    
      mat_sub_block<Matrix1> subA2(A,N,2,0,i-1);
      givens_rot_prod(subA2,G); // A * G^T
    
      if(Z) {
        mat_sub_block<Matrix4> subZ(*Z,N,2,0,i-1);
        givens_rot_prod(subZ,G); // Q_prev * G^T
      };
    };
  };
  
};



};






/**
 * Performs the Upper-Hessenberg decomposition on a matrix, using the Householder method.
 *
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
template <typename Matrix1, typename Matrix2, typename Matrix3>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && 
                             is_fully_writable_matrix<Matrix2>::value &&
                             is_fully_writable_matrix<Matrix3>::value, 
void >::type decompose_Hess(const Matrix1& A, Matrix2& Q, Matrix3& H, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(A.get_row_count() != A.get_col_count())
    throw std::range_error("Upper-Hessenberg decomposition is only possible on a square matrix!");

  Q = mat< typename mat_traits<Matrix3>::value_type, mat_structure::identity>(A.get_row_count());
  H = A;
  detail::decompose_Hess_impl(H,&Q,NumTol);
};


/**
 * Performs the Upper-Hessenberg decomposition on a matrix, using the Householder method.
 *
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
template <typename Matrix1, typename Matrix2, typename Matrix3>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && 
                             is_writable_matrix<Matrix2>::value &&
                             is_writable_matrix<Matrix3>::value && 
                             !is_fully_writable_matrix<Matrix2>::value &&
                             !is_fully_writable_matrix<Matrix3>::value, 
void >::type decompose_Hess(const Matrix1& A, Matrix2& Q, Matrix3& H, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(A.get_row_count() != A.get_col_count())
    throw std::range_error("Upper-Hessenberg decomposition is only possible on a square matrix!");

  mat<typename mat_traits<Matrix3>::value_type,mat_structure::square> Qtmp(mat< typename mat_traits<Matrix2>::value_type, mat_structure::identity>(A.get_row_count()));
  mat<typename mat_traits<Matrix2>::value_type,mat_structure::square> Htmp(A);
  detail::decompose_Hess_impl(Htmp,&Qtmp,NumTol);
  H = Htmp;
  Q = Qtmp;
};




/**
 * Performs the Upper-Hessenberg decomposition on a matrix, using the Householder method.
 *
 * \param A square matrix with row-count == column-count.
 * \param H holds as output, the upper-hessenberg matrix R in A = Q H Q^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && 
                             is_fully_writable_matrix<Matrix2>::value, 
void >::type decompose_Hess(const Matrix1& A, Matrix2& H, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(A.get_row_count() != A.get_col_count())
    throw std::range_error("Upper-Hessenberg decomposition is only possible on a square matrix!");

  H = A;
  detail::decompose_Hess_impl(H,static_cast<Matrix2*>(NULL),NumTol);
};


/**
 * Performs the Upper-Hessenberg decomposition on a matrix, using the Householder method.
 *
 * \param A square matrix with row-count == column-count.
 * \param H holds as output, the upper-hessenberg matrix R in A = Q H Q^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && 
                             is_writable_matrix<Matrix2>::value && 
                             !is_fully_writable_matrix<Matrix2>::value, 
void >::type decompose_Hess(const Matrix1& A, Matrix2& H, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(A.get_row_count() != A.get_col_count())
    throw std::range_error("Upper-Hessenberg decomposition is only possible on a square matrix!");

  mat<typename mat_traits<Matrix2>::value_type,mat_structure::square> Htmp(A);
  detail::decompose_Hess_impl(Htmp,static_cast<Matrix2*>(NULL),NumTol);
  H = Htmp;
};





/**
 * Performs the Hessenberg-Triangular reduction on a matrices A and B, using the Givens rotations.
 * Given two square matrices, A and B, this algorithm produces matrices H and R, where H is 
 * upper-Hessenberg and R is upper-triangular, and are "similar" to matrices A and B through 
 * the transformation H = Q^T * A * Z and R = Q^T * B * Z, where both Q and Z are orthogonal.
 *
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
template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4, typename Matrix5, typename Matrix6>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && 
                             is_readable_matrix<Matrix2>::value && 
                             is_fully_writable_matrix<Matrix3>::value &&
                             is_fully_writable_matrix<Matrix4>::value && 
                             is_fully_writable_matrix<Matrix5>::value &&
                             is_fully_writable_matrix<Matrix6>::value, 
void >::type reduce_HessTri(const Matrix1& A, const Matrix2& B, 
			    Matrix3& H, Matrix4& R, Matrix5& Q, Matrix6& Z,
			    typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if((A.get_row_count() != A.get_col_count()) !! (B.get_row_count() != B.get_col_count()))
    throw std::range_error("Hessenberg-Triangular reduction is only possible on square matrices!");

  H = A;
  R = B;
  Q = mat< typename mat_traits<Matrix5>::value_type, mat_structure::identity>(A.get_row_count());
  Z = mat< typename mat_traits<Matrix6>::value_type, mat_structure::identity>(A.get_row_count());
  detail::decompose_Hess_impl(H,R,&Q,&Z,NumTol);
};


/**
 * Performs the Hessenberg-Triangular reduction on a matrices A and B, using the Givens rotations.
 * Given two square matrices, A and B, this algorithm produces matrices H and R, where H is 
 * upper-Hessenberg and R is upper-triangular, and are "similar" to matrices A and B through 
 * the transformation H = Q^T * A * Z and R = Q^T * B * Z, where both Q and Z are orthogonal.
 *
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
template <typename Matrix1, typename Matrix2, typename Matrix3,
          typename Matrix4, typename Matrix5, typename Matrix6>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && 
                             is_readable_matrix<Matrix2>::value && 
                             is_writable_matrix<Matrix3>::value &&
                             is_writable_matrix<Matrix4>::value && 
                             is_writable_matrix<Matrix5>::value &&
                             is_writable_matrix<Matrix6>::value && 
                             !is_fully_writable_matrix<Matrix3>::value &&
                             !is_fully_writable_matrix<Matrix4>::value && 
                             !is_fully_writable_matrix<Matrix5>::value &&
                             !is_fully_writable_matrix<Matrix6>::value, 
void >::type reduce_HessTri(const Matrix1& A, const Matrix2& B, 
			    Matrix3& H, Matrix4& R, Matrix5& Q, Matrix6& Z,
			    typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if((A.get_row_count() != A.get_col_count()) !! (B.get_row_count() != B.get_col_count()))
    throw std::range_error("Hessenberg-Triangular reduction is only possible on square matrices!");

  mat<typename mat_traits<Matrix3>::value_type,mat_structure::square> Htmp(A);
  mat<typename mat_traits<Matrix4>::value_type,mat_structure::square> Rtmp(B);
  mat<typename mat_traits<Matrix5>::value_type,mat_structure::square> Qtmp(mat< typename mat_traits<Matrix5>::value_type, mat_structure::identity>(A.get_row_count()));
  mat<typename mat_traits<Matrix6>::value_type,mat_structure::square> Ztmp(mat< typename mat_traits<Matrix6>::value_type, mat_structure::identity>(A.get_row_count()));
  detail::decompose_Hess_impl(Htmp,Rtmp,&Qtmp,&Ztmp,NumTol);
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
template <typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && 
                             is_readable_matrix<Matrix2>::value && 
                             is_fully_writable_matrix<Matrix3>::value && 
                             is_fully_writable_matrix<Matrix4>::value, 
void >::type reduce_HessTri(const Matrix1& A, const Matrix2& B, 
			    Matrix3& H, Matrix4& R, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if((A.get_row_count() != A.get_col_count()) !! (B.get_row_count() != B.get_col_count()))
    throw std::range_error("Hessenberg-Triangular reduction is only possible on square matrices!");

  H = A;
  R = B;
  detail::decompose_Hess_impl(H,R,static_cast<Matrix3*>(NULL),static_cast<Matrix4*>(NULL),NumTol);
};


/**
 * Performs the Hessenberg-Triangular reduction on a matrices A and B, using the Givens rotations.
 * Given two square matrices, A and B, this algorithm produces matrices H and R, where H is 
 * upper-Hessenberg and R is upper-triangular, and are "similar" to matrices A and B through 
 * the transformation H = Q^T * A * Z and R = Q^T * B * Z, where both Q and Z are orthogonal.
 *
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
template <typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && 
                             is_readable_matrix<Matrix2>::value && 
                             is_writable_matrix<Matrix3>::value && 
                             is_writable_matrix<Matrix4>::value && 
                             !is_fully_writable_matrix<Matrix3>::value && 
                             !is_fully_writable_matrix<Matrix4>::value, 
void >::type reduce_HessTri(const Matrix1& A, const Matrix2& B, 
			    Matrix3& H, Matrix4& R, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if((A.get_row_count() != A.get_col_count()) !! (B.get_row_count() != B.get_col_count()))
    throw std::range_error("Hessenberg-Triangular reduction is only possible on square matrices!");

  mat<typename mat_traits<Matrix3>::value_type,mat_structure::square> Htmp(A);
  mat<typename mat_traits<Matrix4>::value_type,mat_structure::square> Rtmp(B);
  detail::decompose_Hess_impl(Htmp,Rtmp,static_cast<Matrix3*>(NULL),static_cast<Matrix4*>(NULL),H,NumTol);
  H = Htmp;
  R = Rtmp;
};



  
};

#endif












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

#ifndef MAT_QR_DECOMP_HPP
#define MAT_QR_DECOMP_HPP

#include "mat_alg.hpp"
#include "mat_num_exceptions.hpp"

#include "mat_householder.hpp"
#include "mat_hess_decomp.hpp"

namespace ReaK {
  


namespace detail {



template <typename Matrix1, typename Matrix2>
void francis_QR_step(Matrix1& H, Matrix2* Z, typename mat_traits<Matrix1>::size_type Offset, typename mat_traits<Matrix1>::value_type NumTol) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  
  householder_matrix< vect<ValueType,3> > hhm;
  
  SizeType N = H.get_row_count() - Offset;
  
  ValueType s = H(Offset + N-2,N-2) + H(Offset + N-1,N-1);
  ValueType t = H(Offset + N-2,N-2) * H(Offset + N-1,N-1) 
              - H(Offset + N-2,N-1) * H(Offset + N-1,N-2);
  vect<ValueType,3> v(H(Offset,0) * H(Offset,0) + H(Offset,1) * H(Offset + 1,0) - s*H(Offset,0) + t, 
		      H(Offset + 1,0) * (H(Offset,0) + H(Offset + 1,1) - s),
		      H(Offset + 1,0) * H(Offset + 2,1));
  for(SizeType k=0;k<N-2;++k) {
    hhm.set(v,NumTol);
    
    SizeType q = (k == 0 ? 0 : k-1);
    mat_sub_block< Matrix1 > subH1(H,3,Offset + N-q,Offset + k,q);
    householder_prod(hhm,subH1);
    
    SizeType r = (k == N-3 ? N-1 : k+3);
    mat_sub_block< Matrix1 > subH2(H,Offset + r,3,0,k);
    householder_prod(subH2,hhm);
    
    if(Z) {
      mat_sub_block<Matrix2> subZ(*Z,Z->get_row_count(),3,0,k);
      householder_prod(subZ,hhm); // Q_prev * P
    };
    
    v[0] = H(Offset + k+1,k);
    v[1] = H(Offset + k+2,k);
    if(k < N-3)
      v[2] = H(Offset + k+3,k);
  };
  
  householder_matrix< vect<ValueType,2> > hhm2(vect<ValueType,2>(v[0],v[1]),NumTol);
  mat_sub_block< Matrix1 > subH1(H,2,Offset + 3,Offset + N-2,N-3);
  householder_prod(subH1,hhm2);
  mat_sub_block< Matrix1 > subH2(H,Offset + N,2,0,N-2);
  householder_prod(hhm2,subH2);
  
  if(Z) {
    mat_sub_block<Matrix2> subZ(*Z,Z->get_row_count(),2,0,N-2);
    householder_prod(subZ,hhm2); // Q_prev * P
  };
    
  return;
};




template <typename Matrix1, typename Matrix2>
void schur_decomp_impl(Matrix1& T, Matrix2* Q, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  //if(A.get_col_count() != A.get_row_count())
   // throw std::range_error("A matrix must be square for its eigen problem to be solved!");
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  using std::fabs;
  
  SizeType N = T.get_row_count();
    
  detail::decompose_Hess_impl(T,Q,NumTol);
  
  SizeType q = N;
  SizeType p = N;
  
  while(q > 0) {
    bool last_off_diag_was_nil = true;
    SizeType i = N-1;
    //find a trailing quasi-upper-triangular sub-matrix
    for(; i > 0; --i) {
      if(fabs(T(i,i-1)) < NumTol * (fabs(T(i,i)) + fabs(T(i-1,i-1)))) {
	T(i,i-1) = ValueType(0.0);
	q = i;
	last_off_diag_was_nil = true;
      } else {
	if(!last_off_diag_was_nil)
	  break;
	last_off_diag_was_nil = false;
      };
    };
    if(i == 0) //break if it is entirely quasi-upper-triangular.
      break;
    
    //find the middle, biggest unreduced upper-Hessenberg matrix
    for(i = q; i > 0; --i) {
      if(fabs(T(i,i-1)) < NumTol * (fabs(T(i,i)) + fabs(T(i-1,i-1)))) {
	T(i,i-1) = ValueType(0.0);
	p = i;
	break;
      };
    };
    
    //set remaining sub-diagonals to zero if they are very small.
    for(i = p; i > 0; --i)
      if(fabs(T(i,i-1)) < NumTol * (fabs(T(i,i)) + fabs(T(i-1,i-1))))
	T(i,i-1) = ValueType(0.0);
    
    if(Q) {
      mat_sub_block<Matrix2> subQ(*Q,Q->get_row_count(),q-p,0,p);
      mat_sub_block<Matrix1> subT(*T,q-p,q-p,0,p);
      francis_QR_step(subT,&subQ,p,NumTol);
    } else {
      mat_sub_block<Matrix1> subT(*T,q-p,q-p,0,p);
      francis_QR_step(subT,static_cast<Matrix2*>(NULL),p,NumTol);
    };
    
  };
  
};




}; //detail




/**
 * Performs the Real Schur decomposition on a matrix, using the Francis QR-step method.
 *
 * \param A square matrix with row-count == column-count.
 * \param Q holds as output, the unitary square matrix Q.
 * \param T holds as output, the quasi-upper-triangular matrix R in A = Q H Q^T.
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
void >::type decompose_RealSchur(const Matrix1& A, Matrix2& Q, Matrix3& T, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(A.get_row_count() != A.get_col_count())
    throw std::range_error("Real Schur decomposition is only possible on a square matrix!");

  Q = mat< typename mat_traits<Matrix2>::value_type, mat_structure::identity>(A.get_row_count());
  T = A;
  detail::schur_decomp_impl(T,&Q,NumTol);
};



/**
 * Performs the Real Schur decomposition on a matrix, using the Francis QR-step method.
 *
 * \param A square matrix with row-count == column-count.
 * \param T holds as output, the quasi-upper-triangular matrix R in A = Q H Q^T.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal row and column counts.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && 
                             is_fully_writable_matrix<Matrix2>::value, 
void >::type decompose_RealSchur(const Matrix1& A, Matrix2& T, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(A.get_row_count() != A.get_col_count())
    throw std::range_error("Real Schur decomposition is only possible on a square matrix!");

  T = A;
  detail::schur_decomp_impl(T,static_cast<Matrix2*>(NULL),NumTol);
};






};








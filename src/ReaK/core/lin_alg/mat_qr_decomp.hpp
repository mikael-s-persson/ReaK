/**
 * \file mat_qr_decomp.hpp
 * 
 * This library provides a number of functions related to performing a QR-decomposition on a 
 * matrix, e.g., to invert a matrix, to pseudo-invert a matrix, to solve a linear system with 
 * least-square error and to find the determinant. Most implementations provided 
 * are based on the QR-decomposition (via Householder reflections). QR-decomposition is pretty efficient and 
 * is generally preferred if there is reasons to believe that the matrix involved is not always
 * well-conditioned. If a matrix cannot be guaranteed to be well-conditioned,
 * QR-decomposition is preferred to Gaussian 
 * elimination methods (exact methods, like PLU decomposition).
 * 
 * According to performance tests, PLU methods are as good as Cholesky methods in terms of speed.
 * And they are both the best for well-conditioned matrices. For ill-conditioned matrices, QR-decomposition
 * methods are only a little slower then PLU (about 20% slower, same time-complexity) but provide better
 * numerical stability. The Jacobi methods are significantly slower, but this implementation is in need 
 * of a revision for performance enhancement. And, of course, SVD is also very slow (slightly faster than 
 * Jacobi) but it is based on a LAPACK implementation that is very poorly written, and it has not been 
 * updated since.
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

#ifndef REAK_MAT_QR_DECOMP_HPP
#define REAK_MAT_QR_DECOMP_HPP

#include "mat_alg.hpp"
#include "mat_num_exceptions.hpp"

#include "mat_householder.hpp"

#include "mat_cholesky.hpp"

namespace ReaK {
  

  

/*************************************************************************
                Stabilized Gram-Schmidt Orthogonalization
*************************************************************************/

/**
 * Transforms the columns of A through the Stabilized Gram-Schmidt orthogonalization method.
 * It can normalize the columns or not.
 *
 * \param A rectangular matrix with row-count >= column-count which stores, as input, the
 *          a real full-rank matrix, and stores, as output, the orthonormal column vectors.
 * \param Normalize choose to have normal vectors or not. Note that both algorithm have the
 *                  same cost.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A does not have equal-or-more rows than columns.
 *
 * \author Mikael Persson
 */
template <typename Matrix>
typename boost::enable_if_c< is_fully_writable_matrix<Matrix>::value,
void >::type orthogonalize_StableGramSchmidt(Matrix& A, bool Normalize = false, typename mat_traits<Matrix>::value_type NumTol = 1E-8) {
  if(A.get_row_count() < A.get_col_count())
    throw std::range_error("Orthogonalization only possible on a matrix with row-count >= column-count!");
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::size_type SizeType;
  using std::fabs;
  using std::sqrt;
  
  SizeType N = A.get_row_count();
  SizeType M = A.get_col_count();
  ValueType u;
  vect_n<ValueType> v(Normalize ? M : 0);

  for(SizeType i=0;i<M;++i) {
    for(SizeType j=0;j<i;++j) {
      u = 0.0;
      for(SizeType k=0;k<N;++k) {
        u += A(k,i) * A(k,j);
      };
      if (Normalize) {
        for(SizeType k=0;k<N;k++)
          A(k,i) -= u * A(k,j);
      } else {
        for(SizeType k=0;k<N;k++)
          A(k,i) -= u * A(k,j) / v[j];
      };
    };
    if (Normalize) {
      u = 0.0;
      for(SizeType k=0;k<N;k++)
        u += A(k,i) * A(k,i);
      if(fabs(u) < NumTol)
        throw singularity_error("A");
      u = sqrt(u);
      for(SizeType k=0;k<N;k++)
        A(k,i) /= u;
    } else {
      v[i] = 0.0;
      for(SizeType k=0;k<N;++k)
        v[i] += A(k,i) * A(k,i);
      if(fabs(v[i]) < NumTol)
        throw singularity_error("A");
    };
  };
};




/*************************************************************************
                          QR Decomposition
*************************************************************************/

namespace detail {



template <typename Matrix1, typename Matrix2>
void decompose_QR_impl(Matrix1& A, Matrix2* Q, typename mat_traits<Matrix1>::value_type NumTol)
{
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  SizeType N = A.get_row_count();
  SizeType M = A.get_col_count();
  householder_matrix< vect_n<ValueType> > hhm;
  
  SizeType t = (N-1 > M ? M : N-1);

  for(SizeType i = 0; i < t; ++i) {
    
    hhm.set(mat_row_slice<Matrix1>(A,i,i,N - i),NumTol);
    
    mat_sub_block<Matrix1> subA(A,N - i,M - i,i,i);
    householder_prod(hhm,subA); // P * R
    
    if(Q) {
      mat_sub_block<Matrix2> subQ(*Q,Q->get_row_count(),N - i,0,i);
      householder_prod(subQ,hhm); // Q_prev * P
    };
  };
  
};



template <typename Matrix1, typename Matrix2>
void decompose_RQ_impl(Matrix1& A, Matrix2* Q, typename mat_traits<Matrix1>::value_type NumTol)
{
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  SizeType N = A.get_row_count();
  SizeType M = A.get_col_count();
  householder_matrix< vect_n<ValueType> > hhm;
  
  SizeType t = (N > M-1 ? M-1 : N);

  for(SizeType i = 0; i < t; ++i) {
    
    hhm.set_inv(mat_col_slice< Matrix1 >(A, N-i-1, 0, M-i), NumTol);
    
    mat_sub_block<Matrix1> subA(A,N - i,M - i,0,0);
    householder_prod(subA,hhm); // P * R
    
    if(Q) {
      mat_sub_block<Matrix2> subQ(*Q,Q->get_row_count(),M - i,0,0);
      householder_prod(subQ,hhm); // Q_prev * P
    };
  };
  
};



template <typename Matrix1, typename Matrix2>
void backsub_R_impl(const Matrix1& R, Matrix2& b, typename mat_traits<Matrix1>::value_type NumTol)
{
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  using std::fabs;
  SizeType N = R.get_row_count();
  SizeType M = R.get_col_count();
  N = (N > M ? M : N);

  //back-substitution
  for(int i=N-1;i>=0;--i) {
    for(SizeType j=0;j<b.get_col_count();++j) {
      ValueType sum = b(i,j);
      for(SizeType k=i+1;k<N;++k)
        sum -= b(k,j) * R(i,k);
      if(fabs(R(i,i)) < NumTol)
        throw singularity_error("R");
      b(i,j) = sum / R(i,i);
    };
  };
};  


template <typename Matrix1, typename Matrix2>
void forwardsub_L_impl(const Matrix1& L, Matrix2& B, typename mat_traits<Matrix1>::value_type NumTol) {
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  using std::fabs;
  SizeType N = L.get_row_count();
  SizeType M = B.get_col_count();
  for(SizeType j=0;j<M;++j) { //for every column of B
    //Start solving L * Y = B
    for(SizeType i=0;i<N;++i) { //for every row of L
      for(SizeType k=0;k<i;++k) //for every element of row i in L before the diagonal.
        B(i,j) -= L(i,k) * B(k,j); // subtract to B(i,j) the product of L(i,k) * Y(k,j)
      if(fabs(L(i,i)) < NumTol)
        throw singularity_error("L");
      B(i,j) /= L(i,i); // do Y(i,j) = (B(i,j) - sum_k(L(i,k) * Y(k,j))) / L(i,i)
    };
  };
};


/* This implementation is that of Golub and vanLoan, the QR with column pivoting. 
 * Produces P-Q-R as follows:  A P = Q R
   This algorithm is not guaranteed to always reveal the rank. */
template <typename Matrix1, typename Matrix2>
typename mat_traits<Matrix1>::size_type decompose_RRQR_impl(Matrix1& A, Matrix2* Q, mat<typename mat_traits<Matrix1>::value_type,mat_structure::permutation>& P, typename mat_traits<Matrix1>::value_type NumTol)
{
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  using std::swap;
  SizeType N = A.get_row_count();
  SizeType M = A.get_col_count();
  P.set_row_count(M);
  householder_matrix< vect_n<ValueType> > hhm;
  
  SizeType t = (N-1 > M ? M : N-1);
  
  vect_n<ValueType> c(M);
  for(SizeType i = 0; i < M; ++i)
    c[i] = slice(A)(range(0,N),i) * slice(A)(range(0,N),i);
  
  for(SizeType i=0;i<t;++i) {
    
    ValueType tau(0.0);
    SizeType k = i;
    for(SizeType j = i; j < M; ++j) {
      if(c[j] > tau) {
        tau = c[j];
        k = j;
      };
    };
    if(tau < NumTol)
      return i;
    if( k != i ) {
      P.add_column_swap(i,k);
      for(SizeType j = 0; j < N; ++j)
        swap(A(j,i),A(j,k));
      swap(c[i],c[k]);
    };
    
    hhm.set(mat_row_slice<Matrix1>(A,i,i,N - i),NumTol);
    
    mat_sub_block<Matrix1> subA(A,N - i,M - i,i,i);
    householder_prod(hhm,subA); // P * R
    
    if(Q) {
      mat_sub_block<Matrix2> subQ(*Q,Q->get_row_count(),N - i,0,i);
      householder_prod(subQ,hhm); // Q_prev * P
    };
    
    for(SizeType j = i+1; j < M; ++j)
      c[j] -= A(i,j) * A(i,j);
  };
  return t;
};


/* This implementation is that of Gu and Eisenstat (1994 tech. report), the Strong RRQR algorithm. 
   This algorithm is guaranteed to always reveal the rank but will incur higher cost than the simpler 
   column pivoting QR. */
template <typename Matrix1, typename Matrix2>
typename mat_traits<Matrix1>::size_type 
  decompose_StrongRRQR_impl(Matrix1& A, Matrix2* Q, 
                            mat<typename mat_traits<Matrix1>::value_type,mat_structure::permutation>& P, 
                            typename mat_traits<Matrix1>::value_type f, 
                            typename mat_traits<Matrix1>::value_type NumTol)
{
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  using std::swap; using std::sqrt;
  using std::fabs;
  SizeType N = A.get_row_count();
  SizeType M = A.get_col_count();
  P.set_row_count(M);
  householder_matrix< vect_n<ValueType> > hhm;
  
  SizeType t = (N-1 > M ? M : N-1);
  
  mat<ValueType,mat_structure::rectangular> AB(0,M);
  mat<ValueType,mat_structure::rectangular> AB_prev(0,M);
  
  vect_n<ValueType> c(M);
  for(SizeType i = 0; i < M; ++i)
    c[i] = slice(A)(range(0,N),i) * slice(A)(range(0,N),i);
  vect_n<ValueType> w(t,ValueType(0.0));
  
  for(SizeType i=0;i<t;++i) {
    
    ValueType tau(0.0);
    SizeType k = i;
    for(SizeType j = i; j < M; ++j) {
      if(c[j] > tau) {
        tau = c[j];
        k = j;
      };
    };
    
    if(tau < NumTol)
      return i;
    
    if( k != i ) {
      P.add_column_swap(i,k);
      for(SizeType j = 0; j < N; ++j)
        swap(A(j,i),A(j,k));
      swap(c[i],c[k]);
      for(SizeType j = 0; j < i; ++j)
        swap(AB(j,0),AB(j,k-i));
    };
    
    // perform the QR update:
    hhm.set(mat_row_slice<Matrix1>(A,i,i,N - i),NumTol);
    
    mat_sub_block<Matrix1> subA(A,N - i,M - i,i,i);
    householder_prod(hhm,subA); // P * R
    
    if(Q) {
      mat_sub_block<Matrix2> subQ(*Q,Q->get_row_count(),N - i,0,i);
      householder_prod(subQ,hhm); // Q_prev * P
    };
    
    //update c, w, AB.
    AB_prev.set_row_count(i+1);
    AB_prev.set_col_count(M-i-1);
    for(SizeType j = i+1; j < M; ++j) {
      for(SizeType l = 0; l < i; ++l)
        AB_prev(l,j-i-1) = AB(l,j-i) - AB(l,0) * A(i,j) / tau;
      c[j] = sqrt( c[j] * c[j] - A(i,j) * A(i,j) );
      AB_prev(i,j-i-1) = A(i,j) / tau;
    };
    for(SizeType l = 0; l < i; ++l)
      w[l] = sqrt(ValueType(1.0) / ( ValueType(1.0) / (w[l] * w[l]) + (AB(l,0) * AB(l,0)) / (tau * tau)));
    w[i] = tau;
    swap(AB_prev,AB);
    
    if( i == t-1 )
      break;
    
    while( true ) {
      //compute rho
      ValueType rho(0.0);
      SizeType r_i = 0;
      SizeType r_j = i+1;
      for(SizeType j = i+1; j < M; ++j) {
        for(SizeType l = 0; l <= i; ++l) {
          if(rho < fabs(AB(l,j-i-1))) {
            rho = fabs(AB(l,j-i-1));
            r_i = l; r_j = j;
          };
          if(rho < c[j] / w[l]) {
            rho = fabs(AB(l,j-i-1));
            r_i = l; r_j = j;
          };
        };
      };
      
      //break if rho is satisfactory.
      if(rho <= f)
        break;
      
      if( r_j != i+1 ) {
        //swap r_j and i+1:
        P.add_column_swap(r_j,i+1);
        for(SizeType j = 0; j < N; ++j)
          swap(A(j,i+1),A(j,r_j));
        swap(c[i+1],c[r_j]);
        for(SizeType j = 0; j < i; ++j)
          swap(AB(j,0),AB(j,r_j-i-1));
      };
      
      if( r_i != i ) {
        //swap i and r_i:
        P.add_column_swap(r_i,i);
        swap(w[r_i],w[i]);
        for(SizeType j = 0; j < M-i-1; ++j)
          swap(AB(r_i,j),AB(i,j));
        for(SizeType j = 0; j <= i; ++j)
          swap(A(j,r_i),A(j,i));
        
        //do a QR pass on the sub-matrix:
        if(Q) {
          mat_sub_block<Matrix2> subQ2(*Q,Q->get_row_count(),i-r_i+1,0,r_i);
          mat_sub_block<Matrix1> subA2(A,i-r_i+1,M-r_i,r_i,r_i);
          decompose_QR_impl(subA2,&subQ2,NumTol);
        } else {
          mat_sub_block<Matrix1> subA2(A,i-r_i+1,M-r_i,r_i,r_i);
          decompose_QR_impl(subA2,static_cast<Matrix2*>(NULL),NumTol);
        };
        
      };
      
      // perform the QR update:
      hhm.set(mat_row_slice<Matrix1>(A,i+1,i+1,N - i-1),NumTol);
    
      mat_sub_block<Matrix1> subA3(A,N - i-1,M - i-1,i+1,i+1);
      householder_prod(hhm,subA3); // P * R
    
      if(Q) {
        mat_sub_block<Matrix2> subQ3(*Q,Q->get_row_count(),N - i-1,0,i+1);
        householder_prod(subQ3,hhm); // Q_prev * P
      };
      
      for(SizeType j = i+2; j < M; ++j)
        c[j] = sqrt( c[j] * c[j] - A(i,j) * A(i,j) );
      
      //modify the diagonal block of columns i and i + 1.
      ValueType gamma = A(i,i);
      ValueType mu = A(i,i+1) / gamma;
      ValueType nu = A(i+1,i+1) / gamma;
      ValueType rho2 = mu * mu + nu * nu;
      vect_n<ValueType> c2(M-i-2);
      for(SizeType j = i+2; j < M; ++j)
        c2[j-i-2] = A(i+1,j);
      mat<ValueType,mat_structure::rectangular> u(i,1);
      for(SizeType j = 0; j < i; ++j)
        u(j,0) = A(j,i);
      backsub_R_impl(sub(A)(range(0,i),range(0,i)),u,NumTol);
      
      P.add_column_swap(i,i+1);
      for(SizeType j = 0; j <= i+1; ++j)
        swap(A(j,i),A(j,i+1));
      
      hhm.set(mat_row_slice<Matrix1>(A,i,i,2),NumTol);
    
      mat_sub_block<Matrix1> subA4(A,2,M - i,i,i);
      householder_prod(hhm,subA4); // P * R
    
      if(Q) {
        mat_sub_block<Matrix2> subQ4(*Q,Q->get_row_count(),2,0,i);
        householder_prod(subQ4,hhm); // Q_prev * P
      };
      
      for(SizeType j = 0; j < i; ++j)
        w[j] = sqrt( w[j] * w[j] + ((AB(j,0) + mu * u(j,0)) * (AB(j,0) + mu * u(j,0))) / (A(i,i) * A(i,i)) - (u(j,0) * u(j,0)) / (gamma * gamma));
      w[i] = A(i,i);
      
      c[i+1] = A(i+1,i+1);
      for(SizeType j = i+2; j < M; ++j)
        c[j] = sqrt( c[j] * c[j] + A(i+1,j) * A(i+1,j) - c2[j] * c2[j] );
      
      for(SizeType j = i+2; j < M; ++j) {
        for(SizeType l = 0; l < i; ++l)
          AB(l,j-i-1) += (nu * u(l,0) * A(i+1,j) - AB(l,0) * A(i,j) ) / A(i,i);
        AB(i,j-i-1) = A(i,j) / A(i,i);
      };
      for(SizeType j = 0; j < i; ++j)
        AB(j,0) = (nu * nu * u(j,0) - mu * AB(j,0)) / rho2;
      AB(i,0) = mu / rho2;
      
    };
  };
  return t;
};



template <typename Matrix1, typename Matrix2, typename Matrix3>
void linlsq_QR_impl(const Matrix1& A, Matrix2& x,const Matrix3& b, typename mat_traits<Matrix1>::value_type NumTol)
{
  using std::fabs;
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  SizeType N = A.get_row_count();
  SizeType M = A.get_col_count();
  mat<ValueType,mat_structure::rectangular> R = mat<ValueType,mat_structure::rectangular>(A);
  householder_matrix< vect_n<ValueType> > hhm;
  
  mat<ValueType, mat_structure::rectangular> b_store = mat<ValueType, mat_structure::rectangular>(b);

  SizeType t = (N-1 > M ? M : N-1);

  for(SizeType i=0;i<t;++i) {
    
    hhm.set(mat_row_slice< mat<ValueType,mat_structure::rectangular> >(R,i,i,N - i),NumTol);
    
    mat_sub_block< mat<ValueType,mat_structure::rectangular> > subR(R,N - i,M - i,i,i);
    householder_prod(hhm,subR); // P * R
    
    mat_sub_block< mat<ValueType, mat_structure::rectangular> > subb(b_store,b_store.get_row_count() - i,b_store.get_col_count(),i,0);
    householder_prod(hhm,subb); // P * b
    
  };

  //back-substitution
  x.set_row_count(M);
  x.set_col_count(b_store.get_col_count());
  for(int i=M-1;i>=0;--i) {
    for(SizeType j=0;j<b_store.get_col_count();++j) {
      ValueType sum = b_store(i,j);
      for(SizeType k=i+1;k<M;++k)
        sum -= x(k,j) * R(i,k);
      if(fabs(R(i,i)) < NumTol)
        throw singularity_error("R");
      x(i,j) = sum / R(i,i);
    };
  };

};


};

/**
 * Performs the QR decomposition on a matrix, using Householder reflections approach.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A fully-writable matrix type.
 * \tparam Matrix3 A fully-writable matrix type.
 * \param A rectangular matrix with row-count >= column-count, a real full-rank matrix.
 * \param Q holds as output, the orthogonal rectangular matrix Q.
 * \param R holds as output, the upper-triangular or right-triangular matrix R in A = QR.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal-or-more rows than columns.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && 
                             is_fully_writable_matrix<Matrix2>::value &&
                             is_fully_writable_matrix<Matrix3>::value, 
void >::type decompose_QR(const Matrix1& A, Matrix2& Q, Matrix3& R, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(A.get_row_count() < A.get_col_count())
    throw std::range_error("QR decomposition is only possible on a matrix with row-count >= column-count!");

  typedef typename mat_traits<Matrix1>::value_type ValueType;
  
  Q = mat<ValueType,mat_structure::identity>(A.get_row_count());
  R = A;

  detail::decompose_QR_impl(R,&Q,NumTol);
};

/**
 * Performs the QR decomposition on a matrix, using Householder reflections approach.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A fully-writable matrix type.
 * \tparam Matrix3 A writable matrix type.
 * \param A rectangular matrix with row-count >= column-count, a real full-rank matrix.
 * \param Q holds as output, the orthogonal rectangular matrix Q.
 * \param R holds as output, the upper-triangular or right-triangular matrix R in A = QR.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal-or-more rows than columns.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && 
                             is_fully_writable_matrix<Matrix2>::value &&
                             is_writable_matrix<Matrix3>::value &&
                             !is_fully_writable_matrix<Matrix3>::value &&
                             (mat_traits<Matrix3>::structure == mat_structure::upper_triangular), 
void >::type decompose_QR(const Matrix1& A, Matrix2& Q, Matrix3& R, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(A.get_row_count() < A.get_col_count())
    throw std::range_error("QR decomposition is only possible on a matrix with row-count >= column-count!");

  typedef typename mat_traits<Matrix1>::value_type ValueType;
  
  Q = mat<ValueType,mat_structure::identity>(A.get_row_count());
  mat<typename mat_traits<Matrix3>::value_type, mat_structure::rectangular> R_tmp(A);

  detail::decompose_QR_impl(R_tmp,&Q,NumTol);
  R = R_tmp;
};



/**
 * Computes the determinant via QR decomposition of a matrix, using Householder reflections approach.
 *
 * \tparam Matrix A readable matrix type.
 * \param A real square matrix.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 * \return determinant of A, if A is singular, then the determinant is zero but no exception is thrown.
 *
 * \throws std::range_error if the matrix A is not square.
 *
 * \author Mikael Persson
 */
template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
typename mat_traits<Matrix>::value_type >::type determinant_QR(const Matrix& A, typename mat_traits<Matrix>::value_type NumTol = 1E-8) {
  if(A.get_row_count() != A.get_col_count())
    throw std::range_error("Determinant is only defined for a square matrix!");
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::size_type SizeType;
  
  mat<ValueType,mat_structure::square> R(A);
  detail::decompose_QR_impl(R,static_cast<mat<ValueType,mat_structure::square>*>(NULL),NumTol);

  ValueType result(1.0);
  for(SizeType i=0;i<R.get_row_count();++i)
    result *= R(i,i);
  return result;
};


/**
 * Solves the linear least square problem (AX \approx B or X = min_X(||AX - B||)) via Householder reflections.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A fully-writable matrix type.
 * \tparam Matrix3 A readable matrix type.
 * \param A rectangular matrix with row-count >= column-count.
 * \param x stores the solution matrix as output (ColCount x ColCount2).
 * \param b stores the RHS of the linear system of equation (RowCount x ColCount2).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal-or-more rows than columns or if b's
 *                          row count does not match that of A or if x's row count does not match the
 *                          column count of A.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value &&
                             is_readable_matrix<Matrix3>::value,
void >::type linlsq_QR(const Matrix1& A, Matrix2& x, const Matrix3& b, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(A.get_row_count() < A.get_col_count())
    throw std::range_error("Linear Least-square solution is only possible on a matrix with row-count >= column-count!");
  if(A.get_row_count() != b.get_row_count())
    throw std::range_error("Linear Least-square solution is only possible if row count of b is equal to row count of A!");

  detail::linlsq_QR_impl(A,x,b,NumTol);
};


/**
 * Functor to wrap a call to a QR decomposition-based linear-least-square solver.
 */
struct QR_linlsqsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
    linlsq_QR(A,X,B,NumTol);
  };
};



/**
 * Solves the linear minimum-norm problem (AX = B with min_X(||X||)) via QR decomposition.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A fully-writable matrix type.
 * \tparam Matrix3 A readable matrix type.
 * \param A rectangular matrix with row-count <= column-count.
 * \param x stores the solution matrix as output (ColCount x ColCount2).
 * \param b stores the RHS of the linear system of equation (RowCount x ColCount2).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal-or-more rows than columns or if b's
 *                          row count does not match that of A or if x's row count does not match the
 *                          column count of A.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value &&
                             is_readable_matrix<Matrix3>::value,
void >::type minnorm_QR(const Matrix1& A, Matrix2& x, const Matrix3& b, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(A.get_row_count() > A.get_col_count())
    throw std::range_error("Linear Minimum-Norm solution is only possible on a matrix with row-count <= column-count!");
  if(A.get_row_count() != b.get_row_count())
    throw std::range_error("Linear Minimum-Norm solution is only possible if row count of b is equal to row count of A!");
  
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  
  mat<ValueType,mat_structure::rectangular> R(A);
  mat_transpose_view< mat<ValueType,mat_structure::rectangular> > R_t(R);
  mat<ValueType,mat_structure::square> Q = mat<ValueType,mat_structure::square>(mat<ValueType,mat_structure::identity>(A.get_col_count()));
  detail::decompose_QR_impl(R_t, &Q, NumTol);
  
  mat<ValueType,mat_structure::rectangular> b_tmp(b);
  detail::forwardsub_L_impl(R, b_tmp, NumTol);
  x = sub(Q)(range(0, A.get_col_count()),range(0, A.get_row_count())) * b_tmp;
  
};


/**
 * Functor to wrap a call to a QR decomposition-based linear minimum-norm solver.
 */
struct QR_minnormsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
    minnorm_QR(A,X,B,NumTol);
  };
};




/**
 * Performs back-substitution to solve R x = b, where R is an upper-triangular matrix.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A fully-writable matrix type.
 * \tparam Matrix3 A readable matrix type.
 * \param R is an upper-triangular matrix.
 * \param x stores the solution matrix as output, with same dimension as b (zero-padding is assumed in the difference.
 * \param b stores the RHS of the linear system of equation.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if b's row count does not match that of R or if x's row count does not match the
 *                          column count of R.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value &&
                             is_readable_matrix<Matrix3>::value,
void >::type backsub_R(const Matrix1& R, Matrix2& x, const Matrix3& b, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(R.get_row_count() > b.get_row_count())
    throw std::range_error("Back-substitution is only possible if row count of b is equal to row count of R!");
  
  x = b;
  detail::backsub_R_impl(R,x,NumTol);
};


/**
 * Performs back-substitution to solve R x = b, where R is an upper-triangular matrix.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A fully-writable matrix type.
 * \param R is an upper-triangular matrix.
 * \param x stores the solution matrix as output, with same dimension as b (zero-padding is assumed in the difference), 
            also stores b as input.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if b's row count does not match that of R or if x's row count does not match the
 *                          column count of R.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value,
void >::type backsub_R(const Matrix1& R, Matrix2& x, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(R.get_row_count() > x.get_row_count())
    throw std::range_error("Back-substitution is only possible if row count of b is equal to row count of R!");
  
  detail::backsub_R_impl(R,x,NumTol);
};


/**
 * Computes the pseudo-inverse of a matrix via Householder reflections (left/right Moore-Penrose Pseudo-Inverse).
 * A_pinv = (A^T A)^-1 A^T (if M >= N)
 * A_pinv = A^T (A A^T)^-1 (if M < N)
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A fully-writable matrix type.
 * \param A real rectangular matrix.
 * \param A_pinv real rectangular matrix which is the pseudo-inverse of A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal-or-more rows than columns.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value, 
void >::type pseudoinvert_QR(const Matrix1& A, Matrix2& A_pinv, typename mat_traits<Matrix1>::value_type NumTol = 1E-8)  {
  
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  
  if(A.get_row_count() < A.get_col_count()) {
    minnorm_QR(A, A_pinv, mat<ValueType,mat_structure::identity>(A.get_row_count()), NumTol);
  } else {
    linlsq_QR(A, A_pinv, mat<ValueType,mat_structure::identity>(A.get_row_count()), NumTol);
  };
  
};






/**
 * Solves the linear least square problem (AX \approx B or X = min_X(||AX - B||)) via Householder 
 * reflections and a column-pivot strategy to reveal the column-rank of A.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A fully-writable matrix type.
 * \tparam Matrix3 A readable matrix type.
 * \param A rectangular matrix with row-count >= column-count.
 * \param x stores the solution matrix as output (ColCount x ColCount2).
 * \param b stores the RHS of the linear system of equation (RowCount x ColCount2).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal-or-more rows than columns or if b's
 *                          row count does not match that of A or if x's row count does not match the
 *                          column count of A.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value &&
                             is_readable_matrix<Matrix3>::value,
void >::type linlsq_RRQR(const Matrix1& A, Matrix2& x, const Matrix3& b, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(A.get_row_count() < A.get_col_count())
    throw std::range_error("Linear Least-square solution is only possible on a matrix with row-count >= column-count!");
  if(A.get_row_count() != b.get_row_count())
    throw std::range_error("Linear Least-square solution is only possible if row count of b is equal to row count of A!");

  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  
  mat<ValueType,mat_structure::rectangular> R(A);
  mat<ValueType,mat_structure::square> Q = mat<ValueType,mat_structure::square>(mat<ValueType,mat_structure::identity>(A.get_row_count()));
  mat<ValueType,mat_structure::permutation> P(A.get_col_count());
  
  SizeType K = detail::decompose_RRQR_impl(R,&Q,P,NumTol);
  mat_sub_block< mat<ValueType,mat_structure::square> > subQ(Q, A.get_row_count(), A.get_col_count());
  mat<ValueType,mat_structure::rectangular> x_tmp(transpose_view(subQ) * b);
  
  if(K < A.get_col_count()) {
    // A is rank-deficient.
    mat<ValueType,mat_structure::rectangular> x_tmp2(R.get_col_count(), x_tmp.get_col_count());
    minnorm_QR(sub(R)(range(0,K),range(0,R.get_col_count())), x_tmp2, sub(x_tmp)(range(0,K),range(0,x_tmp.get_col_count())), NumTol);
    x_tmp = x_tmp2;
  } else {
    // A is full-rank.
    detail::backsub_R_impl(R,x_tmp,NumTol);
  };
  
  x = P * x_tmp;
};


/**
 * Functor to wrap a call to a QR decomposition-based linear-least-square solver.
 */
struct RRQR_linlsqsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
    linlsq_RRQR(A,X,B,NumTol);
  };
};





/**
 * Solves the linear minimum-norm problem (AX = B with min_X(||X||)) via RRQR decomposition.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A fully-writable matrix type.
 * \tparam Matrix3 A readable matrix type.
 * \param A rectangular matrix with row-count <= column-count.
 * \param x stores the solution matrix as output (ColCount x ColCount2).
 * \param b stores the RHS of the linear system of equation (RowCount x ColCount2).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal-or-more rows than columns or if b's
 *                          row count does not match that of A or if x's row count does not match the
 *                          column count of A.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value &&
                             is_readable_matrix<Matrix3>::value,
void >::type minnorm_RRQR(const Matrix1& A, Matrix2& x, const Matrix3& b, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(A.get_row_count() > A.get_col_count())
    throw std::range_error("Linear Minimum-norm solution is only possible on a matrix with row-count <= column-count!");
  if(A.get_row_count() != b.get_row_count())
    throw std::range_error("Linear Minimum-norm solution is only possible if row count of b is equal to row count of A!");

  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  
  mat<ValueType,mat_structure::rectangular> R(A);
  mat_transpose_view< mat<ValueType,mat_structure::rectangular> > R_t(R);
  mat<ValueType,mat_structure::square> Q = mat<ValueType,mat_structure::square>(mat<ValueType,mat_structure::identity>(A.get_col_count()));
  
  mat<ValueType,mat_structure::permutation> P(A.get_col_count());
  SizeType K = detail::decompose_RRQR_impl(R_t, &Q, P, NumTol);
  
  mat<ValueType,mat_structure::rectangular> b_tmp;
  b_tmp = transpose(P) * b;
  if(K < A.get_row_count()) {
    mat<ValueType,mat_structure::rectangular> x_tmp(K, b.get_col_count());
    detail::linlsq_QR_impl(sub(R)(range(0,A.get_row_count()),range(0,K)), x_tmp, b_tmp, NumTol);
    for(SizeType i = 0; i < K; ++i)
      for(SizeType j = 0; j < b_tmp.get_col_count(); ++j)
        b_tmp(i,j) = x_tmp(i,j);
    for(SizeType i = K; i < A.get_row_count(); ++i)
      for(SizeType j = 0; j < b_tmp.get_col_count(); ++j)
        b_tmp(i,j) = ValueType(0.0);
  } else {
    detail::forwardsub_L_impl(R, b_tmp, NumTol);
  };
  
  x = sub(Q)(range(0, A.get_col_count()),range(0, A.get_row_count())) * b_tmp;
};


/**
 * Functor to wrap a call to a RRQR decomposition-based linear minimum-norm solver.
 */
struct RRQR_minnormsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
    minnorm_RRQR(A,X,B,NumTol);
  };
};

/**
 * Computes the pseudo-inverse of a matrix via RRQR decomposition.
 * A_pinv = (A^T A)^-1 A^T (if M >= N)
 * A_pinv = A^T (A A^T)^-1 (if M < N)
 * 
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A fully-writable matrix type.
 * \param A real rectangular matrix with row-count.
 * \param A_pinv real rectangular matrix which is the pseudo-inverse of A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the matrix A does not have equal-or-more rows than columns.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value, 
void >::type pseudoinvert_RRQR(const Matrix1& A, Matrix2& A_pinv, typename mat_traits<Matrix1>::value_type NumTol = 1E-8)  {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  
  if(A.get_row_count() < A.get_col_count()) {
    minnorm_RRQR(A, A_pinv, mat<ValueType,mat_structure::identity>(A.get_row_count()), NumTol);
  } else {
    linlsq_RRQR(A, A_pinv, mat<ValueType,mat_structure::identity>(A.get_row_count()), NumTol);
  };
};







#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE


extern template void decompose_QR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& Q, mat<double,mat_structure::rectangular>& R, double NumTol);
extern template void decompose_QR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& Q, mat<double,mat_structure::rectangular>& R, double NumTol);
extern template void decompose_QR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& Q, mat<double,mat_structure::square>& R, double NumTol);

extern template double determinant_QR(const mat<double,mat_structure::rectangular>& A, double NumTol);
extern template double determinant_QR(const mat<double,mat_structure::square>& A, double NumTol);

extern template void linlsq_QR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
extern template void linlsq_QR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
extern template void linlsq_QR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& x, const mat<double,mat_structure::square>& b, double NumTol);

extern template void minnorm_QR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
extern template void minnorm_QR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
extern template void minnorm_QR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& x, const mat<double,mat_structure::square>& b, double NumTol);

extern template void pseudoinvert_QR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
extern template void pseudoinvert_QR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& A_pinv, double NumTol);
extern template void pseudoinvert_QR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
extern template void pseudoinvert_QR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& A_pinv, double NumTol);

extern template void linlsq_RRQR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
extern template void linlsq_RRQR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
extern template void linlsq_RRQR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& x, const mat<double,mat_structure::square>& b, double NumTol);

extern template void minnorm_RRQR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
extern template void minnorm_RRQR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
extern template void minnorm_RRQR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& x, const mat<double,mat_structure::square>& b, double NumTol);

extern template void pseudoinvert_RRQR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
extern template void pseudoinvert_RRQR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& A_pinv, double NumTol);
extern template void pseudoinvert_RRQR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
extern template void pseudoinvert_RRQR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& A_pinv, double NumTol);



extern template void decompose_QR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& Q, mat<float,mat_structure::rectangular>& R, float NumTol);
extern template void decompose_QR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& Q, mat<float,mat_structure::rectangular>& R, float NumTol);
extern template void decompose_QR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& Q, mat<float,mat_structure::square>& R, float NumTol);

extern template float determinant_QR(const mat<float,mat_structure::rectangular>& A, float NumTol);
extern template float determinant_QR(const mat<float,mat_structure::square>& A, float NumTol);

extern template void linlsq_QR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
extern template void linlsq_QR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
extern template void linlsq_QR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& x, const mat<float,mat_structure::square>& b, float NumTol);

extern template void minnorm_QR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
extern template void minnorm_QR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
extern template void minnorm_QR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& x, const mat<float,mat_structure::square>& b, float NumTol);

extern template void pseudoinvert_QR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
extern template void pseudoinvert_QR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& A_pinv, float NumTol);
extern template void pseudoinvert_QR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
extern template void pseudoinvert_QR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& A_pinv, float NumTol);

extern template void linlsq_RRQR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
extern template void linlsq_RRQR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
extern template void linlsq_RRQR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& x, const mat<float,mat_structure::square>& b, float NumTol);

extern template void minnorm_RRQR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
extern template void minnorm_RRQR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
extern template void minnorm_RRQR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& x, const mat<float,mat_structure::square>& b, float NumTol);

extern template void pseudoinvert_RRQR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
extern template void pseudoinvert_RRQR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& A_pinv, float NumTol);
extern template void pseudoinvert_RRQR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
extern template void pseudoinvert_RRQR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& A_pinv, float NumTol);


#endif




  
  
};

#endif












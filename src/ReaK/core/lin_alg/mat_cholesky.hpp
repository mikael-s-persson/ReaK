/**
 * \file mat_cholesky.hpp
 * 
 * This library provides the Cholesky decomposition function template which can be used 
 * to find the Cholesky factors of a symmetric matrix. The Cholesky decomposition will take 
 * a symmetric, positive-definite matrix and find a lower-triangular matrix such that the 
 * product of this factor and its transpose produce the original matrix. This method is 
 * especially useful in contexts where well-conditioned positive-definite matrix are involved.
 * Following the decomposition, the obtained lower-triangular matrix factor can be used to 
 * efficiently compute the inverse of the matrix or solve linear systems of equations represented
 * by the matrix and some right-hand-side. This library provides all those methods, overloaded 
 * according to the types of matrices given to it (to do in-place algorithms when possible).
 * 
 * According to performance tests, PLU methods are as good as Cholesky methods in terms of speed.
 * And they are both the best for well-conditioned matrices. For ill-conditioned matrices, QR-decomposition
 * methods are only a little slower then PLU/Cholesky (about 20% slower, same time-complexity) but provide better
 * numerical stability. The Jacobi methods are significantly slower, but this implementation is in need 
 * of a revision for performance enhancement. And, of course, SVD is also very slow (slightly faster than 
 * Jacobi) but it is based on a LAPACK implementation that is very poorly written, and it has not been 
 * updated since.
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

#ifndef REAK_MAT_CHOLESKY_HPP
#define REAK_MAT_CHOLESKY_HPP

#include "mat_alg.hpp"
#include "mat_num_exceptions.hpp"

namespace ReaK {


/*************************************************************************
                        Cholesky Decomposition
*************************************************************************/

namespace detail {

template <typename Matrix1, typename Matrix2>
void decompose_Cholesky_impl(const Matrix1& A, Matrix2& L, typename mat_traits<Matrix1>::value_type NumTol) 
{
  using std::sqrt;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  SizeType N = A.get_row_count();
  for(SizeType i=0;i<N;++i) {
    for(SizeType j=0;j<i;++j) {
      L(i,j) = A(i,j);
      for(SizeType k=0;k<j;++k)
        L(i,j) -= L(i,k) * L(j,k);
      L(i,j) /= L(j,j);
    };
    L(i,i) = A(i,i);
    for(SizeType k=0;k<i;++k) {
      L(i,i) -= L(i,k) * L(i,k);
    };
    if(L(i,i) < NumTol)
      throw singularity_error("A");
    L(i,i) = sqrt(L(i,i));
  };
};

template <typename Matrix1>
void decompose_BandCholesky_impl(Matrix1& A, typename mat_traits<Matrix1>::size_type p, typename mat_traits<Matrix1>::value_type NumTol) 
{
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  using std::sqrt;
  SizeType N = A.get_row_count();
  for(SizeType i = 0; i < N; ++i) {
    for(SizeType j = (i > p ? i - p : 0); j < i; ++j) {
      SizeType k = j + p;
      if(k >= N)
        k = N-1;
      for(SizeType l = i; l <= k; ++l)
        A(l,i) -= A(i,j) * A(l,j);
    };
    SizeType k = i + p;
    if(k >= N)
      k = N-1;
    if(A(i,i) < NumTol)
      throw singularity_error("A");
    A(i,i) = sqrt(A(i,i));
    for(SizeType l = i+1; l <= k; ++l)
      A(l,i) /= A(i,i);
  };
};

template <typename Matrix1>
void decompose_TriDiagLDL_impl(Matrix1& A, typename mat_traits<Matrix1>::value_type NumTol) 
{
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  using std::fabs;
  SizeType N = A.get_row_count();
  if(N == 0) return;
  if(fabs(A(0,0)) < NumTol)
    throw singularity_error("A");
  if(N > 1)
    A(1,0) /= A(0,0);
  for(SizeType i = 1; i < N; ++i) {
    ValueType v_prev = A(i, i-1) * A(i-1, i-1);
    ValueType v = A(i,i) - A(i,i-1) * v_prev;
    if(fabs(v) < NumTol)
      throw singularity_error("A");
    A(i,i) = v;
    if(i < N-1)
      A(i+1,i) = ( A(i+1,i) - A(i+1,i-1) * v_prev ) / v;
  };
};

template <typename Matrix1>
void decompose_LDL_impl(Matrix1& A, typename mat_traits<Matrix1>::value_type NumTol) 
{
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  using std::fabs;
  SizeType N = A.get_row_count();
  vect_n<ValueType> v(N);
  for(SizeType i = 0; i < N; ++i) {
    for(SizeType j = 0; j < i; ++j)
      v[j] = A(i,j) * A(j,j);
    v[i] = A(i,i);
    for(SizeType j = 0; j < i; ++j)
      v[i] -= A(i,j) * v[j];
    A(i,i) = v[i];
    if(fabs(v[i]) < NumTol)
      throw singularity_error("A");
    for(SizeType j = i+1; j < N; ++j) {
      for(SizeType k = 0; k < i; ++k)
        A(j,i) -= A(j,k) * v[k];
      A(j,i) /= v[i];
    };
  };
};


template <typename Matrix1, typename Matrix2>
void backsub_Cholesky_impl(const Matrix1& L, Matrix2& B) {
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  SizeType N = L.get_row_count();
  SizeType M = B.get_col_count();
  for(SizeType j=0;j<M;++j) { //for every column of B
    //Start solving L * Y = B
    for(SizeType i=0;i<N;++i) { //for every row of L
      for(SizeType k=0;k<i;++k) //for every element of row i in L before the diagonal.
        B(i,j) -= L(i,k) * B(k,j); // subtract to B(i,j) the product of L(i,k) * Y(k,j)
      B(i,j) /= L(i,i); // do Y(i,j) = (B(i,j) - sum_k(L(i,k) * Y(k,j))) / L(i,i)
    };
    // Then solve L.transpose() * X = Y
    for(int i=N-1;i>=0;--i) { //for every row of L.transpose(), starting from the last
      for(int k=N-1;k>i;--k) //for every element of row i in L.tranpose(), after the diagonal.
        B(i,j) -= L(k,i) * B(k,j); // subtract to B(i,j) the product of L(k,i) * Y(k,j)
      B(i,j) /= L(i,i);
    };
  };
};


template <typename Matrix1, typename Matrix2>
void backsub_BandCholesky_impl(const Matrix1& L, Matrix2& B, typename mat_traits<Matrix1>::size_type p) {
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  SizeType N = L.get_row_count();
  SizeType M = B.get_col_count();
  for(SizeType j = 0; j < M; ++j) { //for every column of B
    //Start solving L * Y = B
    for(SizeType i = 0; i < N; ++i) { //for every row of L
      for(SizeType k = (i > p ? i - p : 0); k < i; ++k) //for every element of row i in L before the diagonal.
        B(i,j) -= L(i,k) * B(k,j); // subtract to B(i,j) the product of L(i,k) * Y(k,j)
      B(i,j) /= L(i,i); // do Y(i,j) = (B(i,j) - sum_k(L(i,k) * Y(k,j))) / L(i,i)
    };
    // Then solve transpose(L) * X = Y
    for(SizeType i = N; i > 0; ) { //for every row of L.transpose(), starting from the last
      SizeType l = --i + p;
      if(l >= N)
        l = N-1;
      for(SizeType k = l; k > i; --k) //for every element of row i in L.tranpose(), after the diagonal.
        B(i,j) -= L(k,i) * B(k,j); // subtract to B(i,j) the product of L(k,i) * Y(k,j)
      B(i,j) /= L(i,i);
    };
  };
};


template <typename Matrix1, typename Matrix2>
void backsub_LDL_impl(const Matrix1& L, Matrix2& B) {
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  SizeType N = L.get_row_count();
  SizeType M = B.get_col_count();
  for(SizeType j = 0; j < M; ++j) { //for every column of B
    //Start solving L * Y = B
    for(SizeType i = 0; i < N; ++i) //for every row of L
      for(SizeType k = 0; k < i; ++k) //for every element of row i in L before the diagonal.
        B(i,j) -= L(i,k) * B(k,j); // subtract to B(i,j) the product of L(i,k) * Y(k,j)
    // Solve D * X = Y;
    for(SizeType i = 0; i < N; ++i)
      B(i,j) /= L(i,i);
    // Then solve transpose(L) * W = X
    for(SizeType i = N; i > 0; ) { //for every row of L.transpose(), starting from the last
      --i;
      for(SizeType k = N-1; k > i; --k) //for every element of row i in L.tranpose(), after the diagonal.
        B(i,j) -= L(k,i) * B(k,j); // subtract to B(i,j) the product of L(k,i) * Y(k,j)
    };
  };
};

template <typename Matrix1, typename Matrix2>
void backsub_TriDiagLDL_impl(const Matrix1& L, Matrix2& B) {
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  SizeType N = L.get_row_count();
  SizeType M = B.get_col_count();
  
  for(SizeType j = 0; j < M; ++j) { //for every column of B
    //Start solving L * Y = B
    for(SizeType i = 1; i < N; ++i) //for every row of L
                                    //for every element of row i in L before the diagonal.
      B(i,j) -= L(i,i-1) * B(i-1,j); // subtract to B(i,j) the product of L(i,k) * Y(k,j)
    // Solve D * X = Y;
    for(SizeType i = 0; i < N; ++i)
      B(i,j) /= L(i,i);
    // Then solve transpose(L) * W = X
    for(SizeType i = N-1; i > 0; --i) //for every row of L.transpose(), starting from the last
                                      //for every element of row i in L.tranpose(), after the diagonal.
      B(i-1,j) -= L(i,i-1) * B(i,j); // subtract to B(i,j) the product of L(k,i) * Y(k,j)
  };
};



#if 0
template <typename Matrix1, typename Matrix2, typename Vector1, typename Vector2>
void decompose_PLTLP_impl(Matrix1& A, 
                          Matrix2& L, 
                          mat< typename mat_traits<Matrix1>::value_type, mat_structure::permutation >& P,
                          Vector1& e, 
                          Vector2& d, 
                          typename mat_traits<Matrix1>::value_type NumTol) 
{
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  using std::fabs; using std::swap;
  SizeType N = A.get_row_count();
  L = mat<ValueType, mat_structure::identity>(N);
  P = mat<ValueType, mat_structure::permutation>(N);
  if(N < 3) {
    for(SizeType i = 0; i < N; ++i)
      d[i] = A(i,i);
    for(SizeType i = 0; i + 1 < N; ++i)
      e[i] = A(i+1,i);
  };
  for(SizeType i = 1; i < N; ++i)
    L(i,0) = 1.0;
  vect_n<ValueType> v(N);
  vect_n<ValueType> h(N);
  vect_n<ValueType> l(N);
  l[0] = 0.0; 
  for(SizeType j = 0; j < N; ++j) {
    // first, compute h.
    if(j == 0) {
      d[0] = h[0] = A(0,0);
    } else if(j == 1) {
      h[0] = e[0];
      d[1] = h[1] = A(1,1);
    } else {
      l[j] = 1.0;
      for(SizeType i = 1; i < j; ++i)
        l[i] = L(j,i);
      h[0] = l[1] * e[0];
      h[j] = A(j,j);
      for(SizeType k = 1; k < j; ++k) {
        h[k] = e[k-1] * l[k-1] + d[k] * l[k] + e[k] * l[k+1];
        h[j] -= l[k] * h[k];
      };
      d[j] = h[j] - e[j-1] * L(j,j-1);
    };
    
    if(j < N-1) {
      ValueType max_v = -1.0;
      SizeType q = j+1;
      for(SizeType i = j+1; i < N; ++i) {
        v[i] = A(i,j);
        for(SizeType k = 0; k <= j; ++k)
          v[i] -= L(i,k) * h[k];
        if(fabs(v[i]) > max_v) {
          max_v = fabs(v[i]);
          q = i;
        };
      };
      if(q != j+1) {
        // swap everything
        P.add_column_swap(j+1,q);
        swap(v[j+1],v[q]);
        for(SizeType i = 1; i <= j; ++i)
          swap(L(j+1,i),L(q,i));
        for(SizeType i = j+1; i < N; ++i) {
          swap(A(j+1,i),A(q,i));
          swap(A(i,j+1),A(i,q));
        };
      };
      e[j] = v[j+1];
    };
    
    if(j < N-2) {
      for(SizeType i = j + 2; i < N; ++i)
        L(i,j+1) = v[i];
      if(fabs(v[j+1]) > NumTol) {
        for(SizeType i = j + 2; i < N; ++i)
          L(i,j+1) /= v[j+1];
      };
    };
    
  };
};

template <typename Matrix1, typename Vector1, typename Vector2, typename Matrix2>
void backsub_PLTLP_impl(const Matrix1& L, 
                        const mat< typename mat_traits<Matrix1>::value_type, mat_structure::permutation >& P,
                        const Vector1& e_orig, 
                        const Vector2& d_orig, 
                        Matrix2& B,
                        typename mat_traits<Matrix1>::value_type NumTol) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  SizeType N = L.get_row_count();
  SizeType M = B.get_col_count();
  
  B = P * B;
  
  for(SizeType j = 0; j < M; ++j) { //for every column of B
    //Start solving L * Y = B
    for(SizeType i = 0; i < N; ++i) //for every row of L
      for(SizeType k = 0; k < i; ++k) //for every element of row i in L before the diagonal.
        B(i,j) -= L(i,k) * B(k,j); // subtract to B(i,j) the product of L(i,k) * Y(k,j)
  };
  
  vect_n<ValueType> d(N);
  vect_n<ValueType> e(N-1);
  d[0] = d_orig[0];
  if(fabs(d[0]) < NumTol)
    throw singularity_error("A");
  for(SizeType i = 1; i < N; ++i) {
    ValueType tmp = (d[i] = d_orig[i] - e_orig[i-1] * e_orig[i-1] / d[i-1]);
    if(fabs(tmp) < NumTol)
      throw singularity_error("A");
    e[i-1] = e_orig[i-1] / tmp;
  };
  for(SizeType i = 0; i < N-1; ++i) {
    e[i] = e_orig[i] / d[i];
    d[i+1] = d_orig[i+1] - e_orig[i] * e[i];
  };
  
  for(SizeType j = 0; j < M; ++j) { //for every column of B
    //Start solving L * Y = B
    for(SizeType i = 1; i < N; ++i) 
      B(i,j) -= e[i-1] * B(i-1,j);
    B(N-1,j) /= d[N-1];
    // Then solve transpose(L) * W = X
    for(SizeType i = N-1; i > 0; ) { //for every row of L.transpose(), starting from the last
      --i;
      B(i,j) = B(i,j) / d[i] - e[i] * B(i+1,j); // subtract to B(i,j) the product of L(k,i) * Y(k,j)
    };
  };
  
  for(SizeType j = 0; j < M; ++j) { //for every column of B
    // Then solve transpose(L) * W = X
    for(SizeType i = N; i > 0; ) { //for every row of L.transpose(), starting from the last
      --i;
      for(SizeType k = N-1; k > i; --k) //for every element of row i in L.tranpose(), after the diagonal.
        B(i,j) -= L(k,i) * B(k,j); // subtract to B(i,j) the product of L(k,i) * Y(k,j)
    };
  };
  
  B = B * transpose(P);
};
#endif

};


/**
 * Performs the Cholesky decomposition of A (positive-definite symmetric matrix) (Cholesky-Crout Algorithm).
 * returns a Lower-triangular matrix such that A = L * transpose(L).
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param L stores, as output, the lower-triangular matrix in A = L * transpose(L).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient) or not positive-definite.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility to ensure it's correct.
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             ((mat_traits<Matrix1>::structure == mat_structure::square) ||
                              (mat_traits<Matrix1>::structure == mat_structure::symmetric) ||
                              (mat_traits<Matrix1>::structure == mat_structure::tridiagonal)) &&
                             is_writable_matrix<Matrix2>::value &&
                             is_resizable_matrix<Matrix2>::value &&
                             (mat_traits<Matrix2>::structure != mat_structure::lower_triangular),
void >::type decompose_Cholesky(const Matrix1& A, Matrix2& L, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  L.set_col_count(A.get_col_count());
  detail::decompose_Cholesky_impl(A,L,NumTol);
};

/**
 * Performs the Cholesky decomposition of A (positive-definite symmetric matrix) (Cholesky-Crout Algorithm).
 * returns a Lower-triangular matrix such that A = L * L.transpose().
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param L stores, as output, the lower-triangular matrix in A = L * transpose(L).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient) or not positive-definite.
 *
 * \note the positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility to ensure it's correct.
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             ((mat_traits<Matrix1>::structure == mat_structure::square) ||
                              (mat_traits<Matrix1>::structure == mat_structure::symmetric) ||
                              (mat_traits<Matrix1>::structure == mat_structure::tridiagonal)) &&
                             is_writable_matrix<Matrix2>::value &&
                             (mat_traits<Matrix2>::structure == mat_structure::lower_triangular),
void >::type decompose_Cholesky(const Matrix1& A, Matrix2& L, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  mat<ValueType, mat_structure::square> L_tmp(A.get_row_count(),ValueType(0));
  detail::decompose_Cholesky_impl(A,L_tmp,NumTol);
  L = L_tmp;
};

/**
 * Performs the Cholesky decomposition of A (positive-definite symmetric matrix) (Cholesky-Crout Algorithm).
 * returns a Lower-triangular matrix such that A = L * L.transpose().
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param L stores, as output, the lower-triangular matrix in A = L * transpose(L).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient) or not positive-definite.
 *
 * \note the positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility to ensure it's correct.
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             ((mat_traits<Matrix1>::structure == mat_structure::diagonal) || 
                              (mat_traits<Matrix1>::structure == mat_structure::identity)) &&
                             is_writable_matrix<Matrix2>::value,
void >::type decompose_Cholesky(const Matrix1& A, Matrix2& L, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  using std::sqrt;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  L = A;
  for(SizeType i=0;i < A.get_row_count();++i) 
    L(i,i) = sqrt(L(i,i));
};


/**
 * Performs the Cholesky decomposition of A (positive-definite symmetric matrix) (Cholesky-Crout Algorithm),
 * and uses the result to compute the determinant of A.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 * \return the determinant of A, or 0 if it was judged to be singular or not positive-definite.
 *
 * \note the symmetry of the matrix A is not checked and thus it is
 *       the caller's responsibility to ensure it's correct.
 * 
 * \author Mikael Persson
 */
template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
typename mat_traits<Matrix>::value_type >::type determinant_Cholesky(const Matrix& A, typename mat_traits<Matrix>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::value_type SizeType;
  mat<ValueType, mat_structure::square> L(A.get_row_count(),ValueType(0));
  try {
    decompose_Cholesky(A,L,NumTol);
  } catch(singularity_error& e) {
    return ValueType(0);
  };
  ValueType prod = ValueType(1);
  for(SizeType i = 0; i < L.get_col_count(); ++i)
    prod *= L(i,i);
  return prod * prod;
};

/**
 * Solves the linear problem AX = B using Cholesky decomposition.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param b stores, as input, the RHS of the linear system of equations and stores, as output,
 *          the solution matrix X (Size x B_ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility.
 * \note if you wish to apply this method with a vector on the RHS, then use ReaK::mat_vect_adaptor (and related classes).
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             ((mat_traits<Matrix1>::structure == mat_structure::square) ||
                              (mat_traits<Matrix1>::structure == mat_structure::symmetric) ||
                              (mat_traits<Matrix1>::structure == mat_structure::tridiagonal)) &&
                             is_fully_writable_matrix<Matrix2>::value,
void >::type linsolve_Cholesky(const Matrix1& A, Matrix2& b, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(b.get_row_count() != A.get_col_count())
    throw std::range_error("For linear equation solution, matrix b must have same row count as A!");

  typedef typename mat_traits<Matrix1>::value_type ValueType;
  mat<ValueType,mat_structure::square> L(A.get_row_count(),ValueType(0));
  detail::decompose_Cholesky_impl(A,L,NumTol);
  detail::backsub_Cholesky_impl(L,b);
};

/**
 * Solves the linear problem AX = B using Cholesky decomposition.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param b stores, as input, the RHS of the linear system of equation and stores, as output,
 *          the solution matrix x (Size x B_ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility.
 * \note if you wish to apply this method with a vector on the RHS, then use ReaK::mat_vect_adaptor (and related classes).
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             ((mat_traits<Matrix1>::structure == mat_structure::square) ||
                              (mat_traits<Matrix1>::structure == mat_structure::symmetric) ||
                              (mat_traits<Matrix1>::structure == mat_structure::tridiagonal)) &&
                             is_writable_matrix<Matrix2>::value &&
                             !is_fully_writable_matrix<Matrix2>::value,
void >::type linsolve_Cholesky(const Matrix1& A, Matrix2& b, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(b.get_row_count() != A.get_col_count())
    throw std::range_error("For linear equation solution, matrix b must have same row count as A!");

  typedef typename mat_traits<Matrix1>::value_type ValueType;
  mat<ValueType,mat_structure::square> L(A.get_row_count(),ValueType(0));
  detail::decompose_Cholesky_impl(A,L,NumTol);
  mat<typename mat_traits<Matrix2>::value_type, mat_structure::rectangular> b_tmp(b);
  detail::backsub_Cholesky_impl(L,b_tmp);
  b = b_tmp;
};

/**
 * Solves the linear problem AX = B using Cholesky decomposition.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param b stores, as input, the RHS of the linear system of equation and stores, as output,
 *          the solution matrix x (Size x B_ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility.
 * \note if you wish to apply this method with a vector on the RHS, then use ReaK::mat_vect_adaptor (and related classes).
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             (mat_traits<Matrix1>::structure == mat_structure::diagonal) &&
                             is_fully_writable_matrix<Matrix2>::value,
void >::type linsolve_Cholesky(const Matrix1& A, Matrix2& b, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(b.get_row_count() != A.get_col_count())
    throw std::range_error("For linear equation solution, matrix b must have same row count as A!");

  typedef typename mat_traits<Matrix1>::size_type SizeType;
  for(SizeType i = 0; i < A.get_row_count(); ++i) {
    if(A(i,i) < NumTol)
      throw singularity_error("A");
    for(SizeType j = 0; j < b.get_col_count(); ++j) {
      b(i,j) /= A(i,i);
    };
  };
};

/**
 * Solves the linear problem AX = B using Cholesky decomposition.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param b stores, as input, the RHS of the linear system of equation and stores, as output,
 *          the solution matrix x (Size x B_ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility.
 * \note if you wish to apply this method with a vector on the RHS, then use ReaK::mat_vect_adaptor (and related classes).
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             (mat_traits<Matrix1>::structure == mat_structure::diagonal) &&
                             is_writable_matrix<Matrix2>::value &&
                             !is_fully_writable_matrix<Matrix2>::value &&
                             (mat_traits<Matrix2>::structure != mat_structure::diagonal),
void >::type linsolve_Cholesky(const Matrix1& A, Matrix2& b, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(b.get_row_count() != A.get_col_count())
    throw std::range_error("For linear equation solution, matrix b must have same row count as A!");

  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  mat<ValueType,mat_structure::rectangular> b_tmp(b);
  for(SizeType i = 0; i < A.get_row_count(); ++i) {
    if(A(i,i) < NumTol)
      throw singularity_error("A");
    for(SizeType j = 0; j < b.get_col_count(); ++j) {
      b_tmp(i,j) /= A(i,i);
    };
  };
  b = b_tmp;
};

/**
 * Solves the linear problem AX = B using Cholesky decomposition.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param b stores, as input, the RHS of the linear system of equation and stores, as output,
 *          the solution matrix x (Size x B_ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility.
 * \note if you wish to apply this method with a vector on the RHS, then use ReaK::mat_vect_adaptor (and related classes).
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             (mat_traits<Matrix1>::structure == mat_structure::diagonal) &&
                             is_writable_matrix<Matrix2>::value && 
                             (mat_traits<Matrix2>::structure == mat_structure::diagonal),
void >::type linsolve_Cholesky(const Matrix1& A, Matrix2& b, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(b.get_row_count() != A.get_col_count())
    throw std::range_error("For linear equation solution, matrix b must have same row count as A!");

  typedef typename mat_traits<Matrix1>::size_type SizeType;
  for(SizeType i = 0; i < A.get_row_count(); ++i) {
    if(A(i,i) < NumTol)
      throw singularity_error("A");
    b(i,i) /= A(i,i);
  };
};

/**
 * Solves the linear problem AX = B using Cholesky decomposition.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param b stores, as input, the RHS of the linear system of equation and stores, as output,
 *          the solution matrix x (Size x B_ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility.
 * \note if you wish to apply this method with a vector on the RHS, then use ReaK::mat_vect_adaptor (and related classes).
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             (mat_traits<Matrix1>::structure == mat_structure::identity),
void >::type linsolve_Cholesky(const Matrix1& A, Matrix2& b, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(b.get_row_count() != A.get_col_count())
    throw std::range_error("For linear equation solution, matrix b must have same row count as A!");
  //nothing to do.
};


/**
 * Functor to wrap a call to a Cholesky-decomposition-based linear system solver.
 */
struct Cholesky_linsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
    X = B;
    linsolve_Cholesky(A,X,NumTol);
  };
};






/**
 * Inverts a matrix using Cholesky decomposition.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be inverted.
 * \param A_inv stores, as output, the inverse of matrix A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             ((mat_traits<Matrix1>::structure == mat_structure::square) ||
                              (mat_traits<Matrix1>::structure == mat_structure::symmetric) ||
                              (mat_traits<Matrix1>::structure == mat_structure::tridiagonal) ||
                              (mat_traits<Matrix1>::structure == mat_structure::diagonal) ||
                              (mat_traits<Matrix1>::structure == mat_structure::identity)) &&
                             is_writable_matrix<Matrix2>::value,
void >::type invert_Cholesky(const Matrix1& A, Matrix2& A_inv, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  
  mat<ValueType,mat_structure::square> L;
  decompose_Cholesky(A,L,NumTol);
  mat<ValueType,mat_structure::square> result(mat<ValueType,mat_structure::identity>(A.get_col_count()));
  detail::backsub_Cholesky_impl(L,result);
  A_inv = result;
};





/**
 * Performs the Band-Cholesky decomposition of A (positive-definite band-symmetric matrix).
 * returns a Lower-triangular matrix such that A = L * transpose(L).
 *
 * \param A real, positive-definite, band-symmetric, square, full-rank matrix to be decomposed.
 * \param L stores, as output, the lower-triangular matrix in A = L * transpose(L).
 * \param bandwidth the band-width of the diagonal band of the matrix (i.e., a tridiagonal matrix has bandwidth = 1).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient) or not positive-definite.
 *
 * \note the positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility to ensure it's correct.
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_square_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value,
void >::type decompose_BandCholesky(const Matrix1& A, Matrix2& L, typename mat_traits<Matrix1>::size_type bandwidth = 1, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  L = A;
  detail::decompose_BandCholesky_impl(L,bandwidth,NumTol);
  for(SizeType i = 1; i < L.get_col_count(); ++i)
    for(SizeType j = 0; j < i; ++j)
      L(j,i) = ValueType(0.0);
};


/**
 * Performs the Band-Cholesky decomposition of A (positive-definite symmetric matrix),
 * and uses the result to compute the determinant of A.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 * \return the determinant of A, or 0 if it was judged to be singular or not positive-definite.
 *
 * \note the symmetry of the matrix A is not checked and thus it is
 *       the caller's responsibility to ensure it's correct.
 * 
 * \author Mikael Persson
 */
template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
typename mat_traits<Matrix>::value_type >::type determinant_BandCholesky(const Matrix& A, 
                                                                         typename mat_traits<Matrix>::size_type bandwidth,
                                                                         typename mat_traits<Matrix>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::value_type SizeType;
  mat<ValueType, mat_structure::square> L(A);
  try {
    detail::decompose_BandCholesky_impl(L,bandwidth,NumTol);
  } catch(singularity_error& e) {
    return ValueType(0);
  };
  ValueType prod = ValueType(1);
  for(SizeType i = 0; i < L.get_col_count(); ++i)
    prod *= L(i,i);
  return prod * prod;
};

/**
 * Solves the linear problem AX = B using Band-Cholesky decomposition.
 * 
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be decomposed.
 * \param b stores, as input, the RHS of the linear system of equations and stores, as output,
 *          the solution matrix X (Size x B_ColCount).
 * \param bandwidth the band-width of the diagonal band of the matrix (i.e., a tridiagonal matrix has bandwidth = 1).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility.
 * \note if you wish to apply this method with a vector on the RHS, then use ReaK::mat_vect_adaptor (and related classes).
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_square_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value,
void >::type linsolve_BandCholesky(const Matrix1& A, Matrix2& b, 
                                   typename mat_traits<Matrix1>::size_type bandwidth,
                                   typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(b.get_row_count() != A.get_col_count())
    throw std::range_error("For linear equation solution, matrix b must have same row count as A!");
  
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  mat<ValueType,mat_structure::square> L(A);
  detail::decompose_BandCholesky_impl(L,bandwidth,NumTol);
  detail::backsub_BandCholesky_impl(L,b,bandwidth);
};

/**
 * Functor to wrap a call to a Cholesky-decomposition-based linear system solver.
 */
struct BandCholesky_linsolver {
  std::size_t band_width;
  BandCholesky_linsolver(std::size_t aBandWidth = 1) : band_width(aBandWidth) { };
  
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
    X = B;
    linsolve_BandCholesky(A,X,band_width,NumTol);
  };
};


/**
 * Inverts a matrix using Band-Cholesky decomposition.
 *
 * \param A real, positive-definite, symmetric, square, full-rank matrix to be inverted.
 * \param A_inv stores, as output, the inverse of matrix A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square.
 *
 * \note the symmetry or positive-definitiveness of the matrix A is not checked and thus it is
 *       the caller's responsibility.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             ((mat_traits<Matrix1>::structure == mat_structure::square) ||
                              (mat_traits<Matrix1>::structure == mat_structure::symmetric) ||
                              (mat_traits<Matrix1>::structure == mat_structure::tridiagonal) ||
                              (mat_traits<Matrix1>::structure == mat_structure::diagonal) ||
                              (mat_traits<Matrix1>::structure == mat_structure::identity)) &&
                             is_writable_matrix<Matrix2>::value,
void >::type invert_BandCholesky(const Matrix1& A, Matrix2& A_inv, 
                                 typename mat_traits<Matrix1>::size_type bandwidth,
                                 typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  mat<ValueType,mat_structure::square> result(mat<ValueType,mat_structure::identity>(A.get_col_count()));
  linsolve_BandCholesky(A,result,bandwidth,NumTol);
  A_inv = result;
};





/**
 * Performs the LDL decomposition of A (symmetric matrix).
 * returns a lower-triangular and diagonal matrix such that A = L * D * transpose(L).
 * 
 * \param A real, band-symmetric, square, full-rank matrix to be decomposed.
 * \param L stores, as output, the lower-triangular matrix in A = L * D * transpose(L).
 * \param D stores, as output, the diagonal matrix in A = L * D * transpose(L).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient) or not positive-definite.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_square_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value &&
                             is_writable_matrix<Matrix3>::value,
void >::type decompose_LDL(const Matrix1& A, Matrix2& L, Matrix3& D, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  L = A;
  detail::decompose_LDL_impl(L,NumTol);
  for(SizeType i = 0; i < L.get_col_count(); ++i) {
    D(i,i) = L(i,i);
    L(i,i) = ValueType(1.0);
    for(SizeType j = 0; j < i; ++j)
      L(j,i) = ValueType(0.0);
  };
};


/**
 * Performs the LDL decomposition of A (symmetric matrix),
 * and uses the result to compute the determinant of A.
 *
 * \param A real, symmetric, square, full-rank matrix to be decomposed.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 * \return the determinant of A, or 0 if it was judged to be singular or not positive-definite.
 *
 * \author Mikael Persson
 */
template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
typename mat_traits<Matrix>::value_type >::type determinant_LDL(const Matrix& A, 
                                                                typename mat_traits<Matrix>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::value_type SizeType;
  mat<ValueType, mat_structure::square> L(A);
  try {
    detail::decompose_LDL_impl(L,NumTol);
  } catch(singularity_error& e) {
    return ValueType(0);
  };
  ValueType prod = ValueType(1.0);
  for(SizeType i = 0; i < L.get_col_count(); ++i)
    prod *= L(i,i);
  return prod; // note: no squaring here since the diagonals are not square-roots (as in Cholesky).
};

/**
 * Solves the linear problem AX = B using LDL decomposition.
 * 
 * \param A real, symmetric, square, full-rank matrix to be decomposed.
 * \param b stores, as input, the RHS of the linear system of equations and stores, as output,
 *          the solution matrix X (Size x B_ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 * 
 * \note if you wish to apply this method with a vector on the RHS, then use ReaK::mat_vect_adaptor (and related classes).
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_square_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value,
void >::type linsolve_LDL(const Matrix1& A, Matrix2& b, 
                          typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(b.get_row_count() != A.get_col_count())
    throw std::range_error("For linear equation solution, matrix b must have same row count as A!");
  
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  mat<ValueType,mat_structure::square> L(A);
  detail::decompose_LDL_impl(L,NumTol);
  detail::backsub_LDL_impl(L,b);
};

/**
 * Functor to wrap a call to a Cholesky-decomposition-based linear system solver.
 */
struct LDL_linsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
    X = B;
    linsolve_LDL(A,X,NumTol);
  };
};


/**
 * Inverts a matrix using LDL decomposition.
 *
 * \param A real, symmetric, square, full-rank matrix to be inverted.
 * \param A_inv stores, as output, the inverse of matrix A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square.
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             ((mat_traits<Matrix1>::structure == mat_structure::square) ||
                              (mat_traits<Matrix1>::structure == mat_structure::symmetric) ||
                              (mat_traits<Matrix1>::structure == mat_structure::tridiagonal) ||
                              (mat_traits<Matrix1>::structure == mat_structure::diagonal) ||
                              (mat_traits<Matrix1>::structure == mat_structure::identity)) &&
                             is_writable_matrix<Matrix2>::value,
void >::type invert_LDL(const Matrix1& A, Matrix2& A_inv, 
                        typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  mat<ValueType,mat_structure::square> result(mat<ValueType,mat_structure::identity>(A.get_col_count()));
  linsolve_LDL(A,result,NumTol);
  A_inv = result;
};







/**
 * Performs the Tridiagonal-LDL decomposition of A (symmetric matrix).
 * returns a lower-triangular and diagonal matrix such that A = L * D * transpose(L).
 * \note This method is currently not working (bug).
 * \param A real, tridiagonal, square, full-rank matrix to be decomposed.
 * \param L stores, as output, the lower-triangular matrix in A = L * D * transpose(L).
 * \param D stores, as output, the diagonal matrix in A = L * D * transpose(L).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient) or not positive-definite.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_square_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value &&
                             is_writable_matrix<Matrix3>::value,
void >::type decompose_TriDiagLDL(const Matrix1& A, Matrix2& L, Matrix3& D, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  L = A;
  detail::decompose_TriDiagLDL_impl(L,NumTol);
  for(SizeType i = 0; i < L.get_col_count(); ++i) {
    D(i,i) = L(i,i);
    L(i,i) = ValueType(1.0);
    for(SizeType j = 0; j < i; ++j)
      L(j,i) = ValueType(0.0);
  };
};


/**
 * Performs the Tridiagonal-LDL decomposition of A (symmetric matrix),
 * and uses the result to compute the determinant of A.
 *
 * \param A real, tridiagonal, square, full-rank matrix to be decomposed.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 * \return the determinant of A, or 0 if it was judged to be singular or not positive-definite.
 *
 * \author Mikael Persson
 */
template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
typename mat_traits<Matrix>::value_type >::type determinant_TriDiagLDL(const Matrix& A, 
                                                                       typename mat_traits<Matrix>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  typedef typename mat_traits<Matrix>::value_type SizeType;
  mat<ValueType, mat_structure::square> L(A);
  try {
    detail::decompose_TriDiagLDL_impl(L,NumTol);
  } catch(singularity_error& e) {
    return ValueType(0);
  };
  ValueType prod = ValueType(1.0);
  for(SizeType i = 0; i < L.get_col_count(); ++i)
    prod *= L(i,i);
  return prod; // note: no squaring here since the diagonals are not square-roots (as in Cholesky).
};

/**
 * Solves the linear problem AX = B using Tridiagonal-LDL decomposition.
 * 
 * \param A real, tridiagonal, square, full-rank matrix to be decomposed.
 * \param b stores, as input, the RHS of the linear system of equations and stores, as output,
 *          the solution matrix X (Size x B_ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 * 
 * \note if you wish to apply this method with a vector on the RHS, then use ReaK::mat_vect_adaptor (and related classes).
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_square_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value,
void >::type linsolve_TriDiagLDL(const Matrix1& A, Matrix2& b, 
                          typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(b.get_row_count() != A.get_col_count())
    throw std::range_error("For linear equation solution, matrix b must have same row count as A!");
  
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  mat<ValueType,mat_structure::square> L(A);
  detail::decompose_TriDiagLDL_impl(L,NumTol);
  detail::backsub_TriDiagLDL_impl(L,b);
};

/**
 * Functor to wrap a call to a Tridiagonal-LDL decomposition-based linear system solver.
 */
struct TriDiagLDL_linsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
    X = B;
    linsolve_TriDiagLDL(A,X,NumTol);
  };
};


/**
 * Inverts a matrix using Tridiagonal-LDL decomposition.
 *
 * \param A real, tridiagonal, square, full-rank matrix to be inverted.
 * \param A_inv stores, as output, the inverse of matrix A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square.
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             ((mat_traits<Matrix1>::structure == mat_structure::square) ||
                              (mat_traits<Matrix1>::structure == mat_structure::symmetric) ||
                              (mat_traits<Matrix1>::structure == mat_structure::tridiagonal) ||
                              (mat_traits<Matrix1>::structure == mat_structure::diagonal) ||
                              (mat_traits<Matrix1>::structure == mat_structure::identity)) &&
                             is_writable_matrix<Matrix2>::value,
void >::type invert_TriDiagLDL(const Matrix1& A, Matrix2& A_inv, 
                        typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  mat<ValueType,mat_structure::square> result(mat<ValueType,mat_structure::identity>(A.get_col_count()));
  linsolve_TriDiagLDL(A,result,NumTol);
  A_inv = result;
};







#if 0
/**
 * Performs the Permuted LTL decomposition of A (symmetric matrix).
 * returns a permutation matrix, a unit-lower-triangular matrix and 
 * tridiagonal matrix such that P * A * transpose(P) = L * T * transpose(L).
 * Golub and vanLoan Alg.-4.4.1
 * 
 * \note This method is currently not working (bug).
 * 
 * \param A real, tridiagonal, square, full-rank matrix to be decomposed.
 * \param L stores, as output, the lower-triangular matrix in P * A * transpose(P) = L * T * transpose(L).
 * \param T stores, as output, the tridiagonal matrix in P * A * transpose(P) = L * T * transpose(L).
 * \param P stores, as output, the permutation matrix in P * A * transpose(P) = L * T * transpose(L).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient) or not positive-definite.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_square_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value &&
                             is_writable_matrix<Matrix3>::value &&
                             is_writable_matrix<Matrix4>::value,
void >::type decompose_PLTLP(const Matrix1& A, Matrix2& L, Matrix3& T, Matrix4& P, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  mat<ValueType, mat_structure::square> A_tmp(A);
  mat<ValueType, mat_structure::permutation> P_tmp(A.get_col_count());
  vect_n<ValueType> e(A.get_col_count());
  vect_n<ValueType> d(A.get_col_count());
  detail::decompose_PLTLP_impl(A_tmp,L,P_tmp,e,d,NumTol);
  for(SizeType i = 0; i < L.get_col_count(); ++i) {
    T(i,i) = d[i];
    if(i > 0)
      T(i-1,i) = T(i,i-1) = e[i-1];
    L(i,i) = ValueType(1.0);
    for(SizeType j = 0; j < i; ++j)
      L(j,i) = ValueType(0.0);
  };
  P = P_tmp;
};


/**
 * Solves the linear problem AX = B using Permuted LTL decomposition.
 * 
 * \param A real, tridiagonal, square, full-rank matrix to be decomposed.
 * \param b stores, as input, the RHS of the linear system of equations and stores, as output,
 *          the solution matrix X (Size x B_ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 * 
 * \note if you wish to apply this method with a vector on the RHS, then use ReaK::mat_vect_adaptor (and related classes).
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_square_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value,
void >::type linsolve_PLTLP(const Matrix1& A, Matrix2& b, 
                          typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(b.get_row_count() != A.get_col_count())
    throw std::range_error("For linear equation solution, matrix b must have same row count as A!");
  
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  mat<ValueType, mat_structure::square> A_tmp(A);
  mat<ValueType,mat_structure::square> L(A);
  mat<ValueType, mat_structure::permutation> P_tmp(A.get_col_count());
  vect_n<ValueType> e(A.get_col_count());
  vect_n<ValueType> d(A.get_col_count());
  detail::decompose_PLTLP_impl(A_tmp, L, P_tmp, e, d, NumTol);
  detail::backsub_PLTLP_impl(L, P_tmp, e, d, b, NumTol);
};

/**
 * Functor to wrap a call to a Permuted LTL decomposition-based linear system solver.
 */
struct PLTLP_linsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
    X = B;
    linsolve_PLTLP(A,X,NumTol);
  };
};


/**
 * Inverts a matrix using Permuted LTL decomposition.
 *
 * \param A real, symmetric, square, full-rank matrix to be inverted.
 * \param A_inv stores, as output, the inverse of matrix A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws singularity_error if the matrix A is singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square.
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             ((mat_traits<Matrix1>::structure == mat_structure::square) ||
                              (mat_traits<Matrix1>::structure == mat_structure::symmetric) ||
                              (mat_traits<Matrix1>::structure == mat_structure::tridiagonal) ||
                              (mat_traits<Matrix1>::structure == mat_structure::diagonal) ||
                              (mat_traits<Matrix1>::structure == mat_structure::identity)) &&
                             is_writable_matrix<Matrix2>::value,
void >::type invert_PLTLP(const Matrix1& A, Matrix2& A_inv, 
                        typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  mat<ValueType,mat_structure::square> result(mat<ValueType,mat_structure::identity>(A.get_col_count()));
  linsolve_PLTLP(A,result,NumTol);
  A_inv = result;
};
#endif




#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

extern template void decompose_Cholesky(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& L, double NumTol);
extern template void decompose_Cholesky(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::square>& L, double NumTol);

extern template void linsolve_Cholesky(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& b, double NumTol);
extern template void linsolve_Cholesky(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::square>& b, double NumTol);

extern template void linsolve_Cholesky(const mat<double,mat_structure::square>& A, mat<double,mat_structure::symmetric>& b, double NumTol);
extern template void linsolve_Cholesky(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::symmetric>& b, double NumTol);

extern template void invert_Cholesky(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& A_inv, double NumTol);
extern template void invert_Cholesky(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::symmetric>& A_inv, double NumTol);


extern template void decompose_Cholesky(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& L, float NumTol);
extern template void decompose_Cholesky(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::square>& L, float NumTol);

extern template void linsolve_Cholesky(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& b, float NumTol);
extern template void linsolve_Cholesky(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::square>& b, float NumTol);

extern template void linsolve_Cholesky(const mat<float,mat_structure::square>& A, mat<float,mat_structure::symmetric>& b, float NumTol);
extern template void linsolve_Cholesky(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::symmetric>& b, float NumTol);

extern template void invert_Cholesky(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& A_inv, float NumTol);
extern template void invert_Cholesky(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::symmetric>& A_inv, float NumTol);


#endif





};

#endif




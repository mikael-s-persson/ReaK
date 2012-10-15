/**
 * \file mat_svd_method.hpp
 * 
 * This library provides a number of functions related to performing a Singular Value Decomposition (SVD)
 * on a matrix, e.g., to invert a matrix, to pseudo-invert a matrix, to solve a linear system with 
 * least-square error and to find the eigen-values of a symmetric matrix. SVD is not very efficient but 
 * very powerful. The SVD implementation used here is based on the implementation from CLARAty, 
 * developed by the Jet Propulsion Laboratory.
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

#ifndef REAK_MAT_SVD_METHOD_HPP
#define REAK_MAT_SVD_METHOD_HPP

#include "mat_alg.hpp"
#include "mat_num_exceptions.hpp"


namespace ReaK {





/*************************************************************************
                    Singular Value Decomposition (SVD)
*************************************************************************/

/**
 * Singular Value Decomposition.
 *
 * For an N-by-M matrix A with N >= M, the singular value decomposition is
 * an N-by-M orthogonal matrix U, an M-by-M diagonal matrix S, and
 * an M-by-M orthogonal matrix V so that A = U*S*V'.
 *
 * The singular values, sigma(k) = S(k, k), are ordered so that
 * sigma(0) >= sigma(1) >= ... >= sigma(M-1).
 *
 * The singular value decompostion always exists, so the constructor will
 * never fail.  The matrix condition number and the effective numerical
 * rank can be computed from this decomposition.
 *
 * NIST Disclaimer:
 * This software was developed at the National Institute of Standards and
 * Technology (NIST) by employees of the Federal Government in the course
 * of their official duties. Pursuant to title 17 Section 105 of the
 * United States Code this software is not subject to copyright
 * protection and is in the public domain. NIST assumes no responsibility
 * whatsoever for its use by other parties, and makes no guarantees,
 * expressed or implied, about its quality, reliability, or any other
 * characteristic.
 *
 * &copy; 2006, Jet Propulsion Laboratory, California Institute of Technology<br>
 *
 * \param A rectangular matrix (RowCount x ColCount).
 * \param U output unitary matrix (RowCount x min(RowCount,ColCount)).
 * \param E vector of singular values, sorted in decreasing order (size = min(RowCount,ColCount))
 * \param V output unitary matrix (ColCount x ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \author NIST
 */
template <typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value &&
                             is_writable_matrix<Matrix3>::value &&
                             (mat_traits<Matrix3>::structure == mat_structure::diagonal) &&
                             is_fully_writable_matrix<Matrix4>::value, 
void >::type decompose_SVD(const Matrix1& A, Matrix2& U, Matrix3& E, Matrix4& V, typename mat_traits<Matrix1>::value_type NumTol = 1E-15) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  using std::swap;
  using std::fabs; 
  using std::sqrt;

  mat<ValueType,mat_structure::rectangular> At(A);
  
  SizeType N = At.get_row_count();
  SizeType M = At.get_col_count();
  SizeType nu = M;
  if(nu > N) nu = N;
  if((N == 0) || (M == 0))
    throw std::range_error("SVD-Decomp: Matrix A has 0 rows or columns!");
  SizeType nct,nrt;
  //SizeType nct = min(N-1,M);
  //SizeType nrt = max(0,min(int(M-2),int(N)));
  if(N >= M+1)
    nct = M;
  else
    nct = N-1;
  if(M >= N+2)
    nrt = N;
  else {
    if(M >= 2)
      nrt = M-2;
    else
      nrt = 0;
  };
  SizeType max_iter = nct;
  if(max_iter < nrt) max_iter = nrt;

  if((U.get_row_count() != N) || (U.get_col_count() != nu)) {
    U.set_row_count(N);
    U.set_col_count(M);
    for(SizeType i = 0; i < N; ++i) {
      for(SizeType j = 0; j < M; ++j) {
	if(i == j)
	  U(i,i) = 1;
	else
	  U(i,j) = 0;
      };
    };
  };
  if((V.get_row_count() != M) || (V.get_col_count() != M))
    V = mat<ValueType,mat_structure::identity>(M);
  if((E.get_row_count() != nu) || (E.get_col_count() != nu))
    E = mat<ValueType,mat_structure::identity>(nu);
  vect_n<ValueType> e(M);
  vect_n<ValueType> work(N);


  // Reduce A to bidiagonal form, storing the diagonal elements
  // in E and the super-diagonal elements in e.
  for (SizeType k=0;k < max_iter;++k) {

    if (k < nct) {

      // Compute the transformation for the k-th column and
      // place the k-th diagonal in E(k).
      // Compute 2-norm of k-th column without under/overflow.
      E(k,k) = 0;
      for (SizeType i=k;i < N;++i)
        E(k,k) = sqrt(E(k,k)*E(k,k) + At(i,k)*At(i,k));

      if (fabs(E(k,k)) > NumTol) {
        if (At(k,k) < 0.0)
          E(k,k) = -E(k,k);
        for (SizeType i=k;i < N;++i)
          At(i,k) /= E(k,k);
        At(k,k) += 1.0;
      };
      E(k,k) = -E(k,k);
    };
    for (SizeType j=k+1;j < M;++j) {
      if ((k < nct) && (fabs(E(k,k)) > NumTol))  {

        // Apply the transformation.
        ValueType t(0.0);
        for (SizeType i=k;i < N;++i)
          t += At(i,k)*At(i,j);

        t = -t/At(k,k);

        for (SizeType i=k;i < N;++i)
	  At(i,j) += t*At(i,k);
      };

      // Place the k-th row of A into e for the
      // subsequent calculation of the row transformation.
      e[j] = At(k,j);
    };
    if (k < nct) {
      // Place the transformation in U for subsequent back
      // multiplication.
      for (SizeType i=k;i < N;++i)
        U(i,k) = At(i,k);
    };
    if (k < nrt) {

      // Compute the k-th row transformation and place the
      // k-th super-diagonal in e(k).
      // Compute 2-norm without under/overflow.
      e[k] = 0;
      for (SizeType i=k+1;i < M;++i)
        e[k] = sqrt(e[k]*e[k] + e[i]*e[i]);
      if (fabs(e[k]) > NumTol) {
        if (e[k+1] < 0.0)
          e[k] = -e[k];
        for (SizeType i=k+1;i < M;++i)
          e[i] /= e[k];
        e[k+1] += 1.0;
      };
      e[k] = -e[k];
      if ((k+1 < N) & (fabs(e[k]) > NumTol)) {

        // Apply the transformation.

        for (SizeType i=k+1;i < N;++i)
          work[i] = 0.0;
        for (SizeType j=k+1;j < M;++j)
          for (SizeType i=k+1;i < N;++i)
            work[i] += e[j]*At(i,j);
        for (SizeType j=k+1;j < M;++j) {
          ValueType t = -e[j]/e[k+1];
          for (SizeType i=k+1;i < N;++i)
            At(i,j) += t*work[i];
        };
      };

        // Place the transformation in _V for subsequent
        // back multiplication.

      for (SizeType i=k+1;i < M;++i)
        V(i,k) = e[i];
    };
  };

  // Set up the final bidiagonal matrix or order p.

  SizeType p = M;
  if(p > N+1) p = N+1;
  if (nct < M)
    E(nct,nct) = At(nct,nct);
  if (N < p)
    E(p-1,p-1) = 0.0;
  if (nrt+1 < p)
    e[nrt] = At(nrt,p-1);

  e[p-1] = 0.0;

  // Generate U.

  for (SizeType j=nct;j < nu;++j) {
    for (SizeType i=0;i<N;++i)
      U(i,j) = 0.0;
    U(j,j) = 1.0;
  };

  for (int k = nct-1; k >= 0;--k) {
    if (fabs(E(k,k)) > NumTol) {
      for (SizeType j=k+1;j < nu;++j) {
        ValueType t = 0;
        for (SizeType i=k;i < N;++i)
          t += U(i,k)*U(i,j);
        t = -t/U(k,k);
        for (SizeType i=k;i < N;++i)
          U(i,j) += t*U(i,k);
      };
      for (SizeType i=k;i < N;++i)
        U(i,k) = -U(i,k);
      U(k,k) = 1.0 + U(k,k);
      for (SizeType i=0;int(i) < k-1;++i)
        U(i,k) = 0.0;
    } else {
      for (SizeType i=0;i < N;++i)
        U(i,k) = 0.0;
      U(k,k) = 1.0;
    };
  };

  // Generate V.

  for (int k=M-1;k >= 0;--k) {
    if ((SizeType(k) < nrt) & (fabs(e[k]) > NumTol)) {
      for (SizeType j=k+1;j < nu;++j) {
        ValueType t = 0;
        for (SizeType i=k+1;i < M;++i)
          t += V(i,k)*V(i,j);
        t = -t/V(k+1,k);
        for (SizeType i=k+1;i < M;++i)
          V(i,j) += t*V(i,k);
      };
    };
    for (SizeType i=0;i < M;++i)
      V(i,k) = 0.0;
    V(k,k) = 1.0;
  };

  // Main iteration loop for the singular values.
  int pp = p-1;
  int iter = 0;
  while (p > 0) {
    int k=0;
    int kase=0;

    // Here is where a test for too many iterations would go.

    // This section of the program inspects for
    // negligible elements in the E and e arrays.  On
    // completion the variables kase and k are set as follows.

    // kase = 1     if E(p) and e(k-1) are negligible and k<p
    // kase = 2     if E(k) is negligible and k<p
    // kase = 3     if e(k-1) is negligible, k<p, and
    //              E(k), ..., E(p)
    //              are not negligible (qr step).
    // kase = 4     if e(p-1) is negligible (convergence).

    for (k=p-2;k >= -1;--k) {
      if (k == -1)
        break;
      if (fabs(e[k]) <= NumTol*(fabs(E(k,k)) + fabs(E(k+1,k+1)))) {
        e[k] = 0.0;
        break;
      };
    };
    if (k == int(p-2)) {
      kase = 4;
    } else {
      int ks;
      for (ks = p-1;ks >= k;--ks) {
        if (ks == k) {
          break;
        };
        ValueType t = (ks != int(p) ? fabs(e[ks]) : 0.0) + (ks != k+1 ? fabs(e[ks-1]) : 0.0);
        if (fabs(E(ks,ks)) <= NumTol*t)  {
          E(ks,ks) = 0.0;
          break;
        };
      };
      if (ks == k) {
        kase = 3;
      } else if (ks == int(p-1)) {
        kase = 1;
      } else {
        kase = 2;
        k = ks;
      };
    };
    k++;

    // Perform the task indicated by kase.

    switch (kase) {

      // Deflate negligible E(p).
    case 1: {
      ValueType f = e[p-2];
      e[p-2] = 0.0;
      for (int j = p-2;j >= k;--j) {
        ValueType t = sqrt(E(j,j)*E(j,j) + f*f);
        ValueType cs = E(j,j)/t;
        ValueType sn = f/t;
        E(j,j) = t;
        if (j != k) {
          f = -sn*e[j-1];
          e[j-1] = cs*e[j-1];
        };
        for (SizeType i=0;i < M;++i) {
          t = cs*V(i,j) + sn*V(i,p-1);
          V(i,p-1) = -sn*V(i,j) + cs*V(i,p-1);
          V(i,j) = t;
        };
      };
    };
    break;

      // Split at negligible E(k).

    case 2: {
      ValueType f = e[k-1];
      e[k-1] = 0.0;
      for (SizeType j=k;j < p;++j) {
        ValueType t = sqrt(E(j,j)*E(j,j) + f*f);
        ValueType cs = E(j,j)/t;
        ValueType sn = f/t;
        E(j,j) = t;
        f = -sn*e[j];
        e[j] = cs*e[j];

        for (SizeType i=0;i < N;++i) {
          t = cs*U(i,j) + sn*U(i,k-1);
          U(i,k-1) = -sn*U(i,j) + cs*U(i,k-1);
          U(i,j) = t;
        };
      };
    };
    break;

      // Perform one qr step.

    case 3: {

      // Calculate the shift.

      //ValueType scale = MAX(MAX(MAX(MAX(fabs(E_a[p-1]),fabs(E_a[p-2])),fabs(e[p-2])),fabs(E_a[k])),fabs(e[k]));
      ValueType scale = fabs(E(p-1,p-1));
      if(scale < fabs(E(p-2,p-2))) scale = fabs(E(p-2,p-2));
      if(scale < fabs(E(k,k))) scale = fabs(E(k,k));
      if(scale < fabs(e[p-2])) scale = fabs(e[p-2]);
      if(scale < fabs(e[k])) scale = fabs(e[k]);
      ValueType sp = E(p-1,p-1)/scale;
      ValueType spm1 = E(p-2,p-2)/scale;
      ValueType epm1 = e[p-2]/scale;
      ValueType sk = E(k,k)/scale;
      ValueType ek = e[k]/scale;
      ValueType b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
      ValueType c = (sp*epm1)*(sp*epm1);
      ValueType shift = 0.0;
      if ((b != 0.0) | (c != 0.0)) {
        shift = sqrt(b*b + c);
        if (b < 0.0)
          shift = -shift;
        shift = c/(b + shift);
      };
      ValueType f = (sk + sp)*(sk - sp) + shift;
      ValueType g = sk*ek;

      // Chase zeros.

      for (SizeType j=k;j < p-1;++j) {
        ValueType t = sqrt(f*f + g*g);
        ValueType cs = f/t;
        ValueType sn = g/t;
        if (j != SizeType(k))
          e[j-1] = t;
        f = cs*E(j,j) + sn*e[j];
        e[j] = cs*e[j] - sn*E(j,j);
        g = sn*E(j+1,j+1);
        E(j+1,j+1) = cs*E(j+1,j+1);

        for (SizeType i=0;i < M;++i) {
          t = cs*V(i,j) + sn*V(i,j+1);
          V(i,j+1) = -sn*V(i,j) + cs*V(i,j+1);
          V(i,j) = t;
        };
        t = sqrt(f*f + g*g);
        cs = f/t;
        sn = g/t;
        E(j,j) = t;
        f = cs*e[j] + sn*E(j+1,j+1);
        E(j+1,j+1) = -sn*e[j] + cs*E(j+1,j+1);
        g = sn*e[j+1];
        e[j+1] = cs*e[j+1];
        if (j < N-1) {
          for (SizeType i=0;i < N;++i) {
            t = cs*U(i,j) + sn*U(i,j+1);
            U(i,j+1) = -sn*U(i,j) + cs*U(i,j+1);
            U(i,j) = t;
          };
        };
      };
      e[p-2] = f;
      iter = iter + 1;
    };
      break;

      // Convergence.

    case 4: {

      // Make the singular values positive.

      if (E(k,k) <= 0.0) {
        E(k,k) = (E(k,k) < 0.0 ? -E(k,k) : 0.0);
        for (SizeType i=0;i <= SizeType(pp);++i)
          V(i,k) = -V(i,k);
      };

      // Order the singular values.

      while (k < pp) {
        if (E(k,k) >= E(k+1,k+1))
          break;
	swap(E(k,k),E(k+1,k+1));
        if (k < int(M-1))
          for (SizeType i = 0;i < M;++i)
	    swap(V(i,k),V(i,k+1));
        if (k < int(N-1))
          for (SizeType i=0;i < N;++i)
	    swap(U(i,k),U(i,k+1));
        k++;
      };
      iter = 0;
      p--;
    };
      break;
    };
  };

  return;
};

/**
 * Singular Value Decomposition.
 *
 * For an N-by-M matrix A with N >= M, the singular value decomposition is
 * an N-by-M orthogonal matrix U, an M-by-M diagonal matrix S, and
 * an M-by-M orthogonal matrix V so that A = U*S*V'.
 *
 * The singular values, sigma(k) = S(k, k), are ordered so that
 * sigma(0) >= sigma(1) >= ... >= sigma(M-1).
 *
 * The singular value decompostion always exists, so the constructor will
 * never fail.  The matrix condition number and the effective numerical
 * rank can be computed from this decomposition.
 *
 * NIST Disclaimer:
 * This software was developed at the National Institute of Standards and
 * Technology (NIST) by employees of the Federal Government in the course
 * of their official duties. Pursuant to title 17 Section 105 of the
 * United States Code this software is not subject to copyright
 * protection and is in the public domain. NIST assumes no responsibility
 * whatsoever for its use by other parties, and makes no guarantees,
 * expressed or implied, about its quality, reliability, or any other
 * characteristic.
 *
 * &copy; 2006, Jet Propulsion Laboratory, California Institute of Technology<br>
 *
 * \param A rectangular matrix (RowCount x ColCount).
 * \param U output unitary matrix (RowCount x min(RowCount,ColCount)).
 * \param E vector of singular values, sorted in decreasing order (size = min(RowCount,ColCount))
 * \param V output unitary matrix (ColCount x ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \author NIST
 */
template <typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_writable_matrix<Matrix2>::value &&
                             is_writable_matrix<Matrix3>::value &&
                             (mat_traits<Matrix3>::structure == mat_structure::diagonal) &&
                             is_writable_matrix<Matrix4>::value &&
                             (!is_fully_writable_matrix<Matrix2>::value ||
                              !is_fully_writable_matrix<Matrix4>::value), 
void >::type decompose_SVD(const Matrix1& A, Matrix2& U, Matrix3& E, Matrix4& V, typename mat_traits<Matrix1>::value_type NumTol = 1E-15) {
  mat<typename mat_traits<Matrix2>::value_type, mat_structure::rectangular> U_tmp(U);
  mat<typename mat_traits<Matrix4>::value_type, mat_structure::rectangular> V_tmp(V);
  decompose_SVD(A,U_tmp,E,V_tmp,NumTol);
  U = U_tmp;
  V = V_tmp;
};
  
  

/**
 * This function returns the two norm of a singular value decomposition, i.e. the
 * highest singular value.
 *
 * \param E vector of singular values, sorted in decreasing order.
 * \return the two-norm.
 *
 * \throws std::range_error if the vector E is empty.
 *
 * \author Mikael Persson
 */
template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
typename mat_traits<Matrix>::value_type >::type two_norm_SVD(const Matrix& E) {
  if(E.get_row_count() == 0)
    throw std::range_error("No singular values available for 2-norm evaluation!");
  return E(0,0);
};

/**
 * This function returns the two norm (induced norm) of a matrix via the computation of 
 * its singular value decomposition.
 * 
 * \param M rectangular matrix (RowCount x ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 * \return the two-norm.
 */
template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
typename Matrix::value_type >::type norm_2(const Matrix& A, typename Matrix::value_type NumTol = 1E-15) {
  typedef typename mat_traits<Matrix>::value_type ValueType;
  
  int nu = (A.get_row_count() > A.get_col_count() ? A.get_col_count() : A.get_row_count());
  mat<ValueType,mat_structure::rectangular> U(A.get_row_count(),nu);
  mat<ValueType,mat_structure::diagonal> E(nu);
  mat<ValueType,mat_structure::square> V(A.get_col_count());

  decompose_SVD(A,U,E,V,NumTol);
  
  return E(0,0); 
};

/**
 * This function returns the condition number of a singular value decomposition,
 * i.e. max(E)/min(E).
 *
 * \param E vector of singular values, sorted in decreasing order.
 * \return the condition number.
 *
 * \throws std::range_error if the vector E is empty.
 *
 * \author Mikael Persson
 */
template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
typename mat_traits<Matrix>::value_type >::type condition_number_SVD(const Matrix& E) {
  if(E.get_row_count() == 0)
    throw std::range_error("No singular values available for condition number evaluation!");
  return E(0,0) / E(E.get_row_count()-1,E.get_row_count()-1);
};

/**
 * This function computes the effective numerical rank of a singular value
 * decomposition.
 *
 * \param E vector of singular values, sorted in decreasing order.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 * \return the numerical rank.
 *
 * \author Mikael Persson
 */
template <typename Matrix>
typename boost::enable_if_c< is_readable_matrix<Matrix>::value,
typename mat_traits<Matrix>::size_type >::type numrank_SVD(const Matrix& E, typename mat_traits<Matrix>::value_type NumTol = 1E-8) {
  using std::fabs;
  typename mat_traits<Matrix>::size_type r = 0;
  for (typename mat_traits<Matrix>::size_type i=0;i<E.get_row_count();++i) {
    if (fabs(E(i,i)) > fabs(E(0,0))*NumTol)
      ++r;
  };
  return r;
};

/**
 * Computes the pseudo-inverse of a RowCount x ColCount matrix already SV-Decomposed.
 *
 * \param U input unitary matrix (RowCount x min(RowCount,ColCount)).
 * \param E vector of singular values, sorted in decreasing order.
 * \param V input unitary matrix (ColCount x ColCount).
 * \param A_pinv the pseudo-inverse (ColCount x RowCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the dimensions don't match.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_readable_matrix<Matrix2>::value &&
                             is_readable_matrix<Matrix3>::value &&
                             is_fully_writable_matrix<Matrix4>::value,
void >::type pseudoinvert_SVD(const Matrix1& U, const Matrix2& E, const Matrix3& V, Matrix4& A_pinv, typename mat_traits<Matrix2>::value_type NumTol = 1E-15) {
  if((U.get_col_count() != E.get_row_count()) || (E.get_row_count() > V.get_row_count()))
    throw std::range_error("Dimensions of the U E V matrices don't match a singular value decomposition!");
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  using std::fabs;
  
  A_pinv.set_row_count(V.get_col_count());
  A_pinv.set_col_count(U.get_row_count());
  for(SizeType i=0;i<V.get_col_count();++i) {
    for(SizeType j=0;j<U.get_row_count();++j) {
      A_pinv(i,j) = 0.0;
      for(SizeType k=0;k<E.get_row_count();++k) {
        if(fabs(E(k,k)) > fabs(E(0,0)) * NumTol)
          A_pinv(i,j) += V(i,k) * U(j,k) / E(k,k);
      };
    };
  };
};

/**
 * Computes the pseudo-inverse of a RowCount x ColCount matrix already SV-Decomposed.
 *
 * \param U input unitary matrix (RowCount x min(RowCount,ColCount)).
 * \param E vector of singular values, sorted in decreasing order.
 * \param V input unitary matrix (ColCount x ColCount).
 * \param A_pinv the pseudo-inverse (ColCount x RowCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the dimensions don't match.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename Matrix3, typename Matrix4>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_readable_matrix<Matrix2>::value &&
                             is_readable_matrix<Matrix3>::value &&
                             is_writable_matrix<Matrix4>::value &&
                             !is_fully_writable_matrix<Matrix4>::value,
void >::type pseudoinvert_SVD(const Matrix1& U, const Matrix2& E, const Matrix3& V, Matrix4& A_pinv, typename mat_traits<Matrix2>::value_type NumTol = 1E-15) {
  mat<typename mat_traits<Matrix1>::value_type, mat_structure::rectangular> A_pinv_tmp(V.get_col_count(),U.get_row_count());
  pseudoinvert_SVD(U,E,V,A_pinv_tmp,NumTol);
  A_pinv = A_pinv_tmp;
};



/**
 * Computes the pseudo-inverse of a NxM matrix A, by performing Singular Value Decomposition.
 *
 * \param A rectangular matrix (RowCount x ColCount).
 * \param A_pinv the pseudo-inverse (ColCount x RowCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_fully_writable_matrix<Matrix2>::value, 
void >::type pseudoinvert_SVD(const Matrix1& A, Matrix2& A_pinv, typename mat_traits<Matrix1>::value_type NumTol = 1E-15) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  using std::fabs;
  
  int nu = (A.get_row_count() > A.get_col_count() ? A.get_col_count() : A.get_row_count());
  mat<ValueType,mat_structure::rectangular> U(A.get_row_count(),nu);
  mat<ValueType,mat_structure::diagonal> E(nu);
  mat<ValueType,mat_structure::square> V(A.get_col_count());

  decompose_SVD(A,U,E,V,NumTol);

  A_pinv.set_row_count(V.get_col_count());
  A_pinv.set_col_count(U.get_row_count());

  for(SizeType i=0;i<V.get_col_count();++i) {
    for(SizeType j=0;j<U.get_row_count();++j) {
      A_pinv(i,j) = 0.0;
      for(SizeType k=0;k<E.get_row_count();++k) {
        if(fabs(E(k,k)) > fabs(E(0,0)) * NumTol)
          A_pinv(i,j) += V(i,k) * U(j,k) / E(k,k);
      };
    };
  };
};

/**
 * Computes the pseudo-inverse of a NxM matrix A, by performing Singular Value Decomposition.
 *
 * \param A rectangular matrix (RowCount x ColCount).
 * \param A_pinv the pseudo-inverse (ColCount x RowCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_writable_matrix<Matrix2>::value &&
                             !is_fully_writable_matrix<Matrix2>::value, 
void >::type pseudoinvert_SVD(const Matrix1& A, Matrix2& A_pinv, typename mat_traits<Matrix1>::value_type NumTol = 1E-15) {
  mat<typename mat_traits<Matrix1>::value_type, mat_structure::rectangular> A_pinv_tmp(A.get_col_count(),A.get_row_count());
  pseudoinvert_SVD(A,A_pinv_tmp,NumTol);
  A_pinv = A_pinv_tmp;
};
  
  
/**
 * Functor to wrap a call to a SVD-based linear-least-square solver.
 */
struct SVD_linlsqsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
    mat<typename mat_traits<Matrix1>::value_type, mat_structure::rectangular> A_pinv(A.get_col_count(),A.get_row_count());
    pseudoinvert_SVD(A,A_pinv,NumTol);
    X = A_pinv * B;
  };
};

  




#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

extern template void decompose_SVD(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& U, mat<double,mat_structure::diagonal>& E, mat<double,mat_structure::rectangular>& V, double NumTol);
extern template void decompose_SVD(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& U, mat<double,mat_structure::diagonal>& E, mat<double,mat_structure::square>& V, double NumTol);
extern template void decompose_SVD(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& U, mat<double,mat_structure::diagonal>& E, mat<double,mat_structure::square>& V, double NumTol);
extern template void decompose_SVD(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& U, mat<double,mat_structure::diagonal>& E, mat<double,mat_structure::square>& V, double NumTol);

extern template void pseudoinvert_SVD(const mat<double,mat_structure::rectangular>& U, const mat<double,mat_structure::diagonal>& E, const mat<double,mat_structure::rectangular>& V, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
extern template void pseudoinvert_SVD(const mat<double,mat_structure::rectangular>& U, const mat<double,mat_structure::diagonal>& E, const mat<double,mat_structure::square>& V, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
extern template void pseudoinvert_SVD(const mat<double,mat_structure::square>& U, const mat<double,mat_structure::diagonal>& E, const mat<double,mat_structure::square>& V, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
extern template void pseudoinvert_SVD(const mat<double,mat_structure::square>& U, const mat<double,mat_structure::diagonal>& E, const mat<double,mat_structure::square>& V, mat<double,mat_structure::square>& A_pinv, double NumTol);

extern template void pseudoinvert_SVD(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
extern template void pseudoinvert_SVD(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& A_pinv, double NumTol);
extern template void pseudoinvert_SVD(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
extern template void pseudoinvert_SVD(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& A_pinv, double NumTol);


extern template void decompose_SVD(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& U, mat<float,mat_structure::diagonal>& E, mat<float,mat_structure::rectangular>& V, float NumTol);
extern template void decompose_SVD(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& U, mat<float,mat_structure::diagonal>& E, mat<float,mat_structure::square>& V, float NumTol);
extern template void decompose_SVD(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& U, mat<float,mat_structure::diagonal>& E, mat<float,mat_structure::square>& V, float NumTol);
extern template void decompose_SVD(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& U, mat<float,mat_structure::diagonal>& E, mat<float,mat_structure::square>& V, float NumTol);

extern template void pseudoinvert_SVD(const mat<float,mat_structure::rectangular>& U, const mat<float,mat_structure::diagonal>& E, const mat<float,mat_structure::rectangular>& V, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
extern template void pseudoinvert_SVD(const mat<float,mat_structure::rectangular>& U, const mat<float,mat_structure::diagonal>& E, const mat<float,mat_structure::square>& V, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
extern template void pseudoinvert_SVD(const mat<float,mat_structure::square>& U, const mat<float,mat_structure::diagonal>& E, const mat<float,mat_structure::square>& V, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
extern template void pseudoinvert_SVD(const mat<float,mat_structure::square>& U, const mat<float,mat_structure::diagonal>& E, const mat<float,mat_structure::square>& V, mat<float,mat_structure::square>& A_pinv, float NumTol);

extern template void pseudoinvert_SVD(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
extern template void pseudoinvert_SVD(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& A_pinv, float NumTol);
extern template void pseudoinvert_SVD(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
extern template void pseudoinvert_SVD(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& A_pinv, float NumTol);
 


#endif





  

};

#endif






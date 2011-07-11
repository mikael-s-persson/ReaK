
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

  for(SizeType i=0;i<t;++i) {
    
    hhm.set(mat_row_slice<Matrix1>(A,i,i,N - i),NumTol);
    
    mat_sub_block<Matrix1> subA(A,N - i,M - i,i,i);
    householder_prod(hhm,subA); // P * R
    
    if(Q) {
      mat_sub_block<Matrix2> subQ(*Q,N,N - i,0,i);
      householder_prod(subQ,hhm); // Q_prev * P
    };
  };
  
};

  



template <typename Matrix1, typename Matrix2, typename Matrix3>
void linlsq_QR_impl(const Matrix1& A, Matrix2& x,const Matrix3& b, typename mat_traits<Matrix1>::value_type NumTol)
{
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
      x(i,j) = sum / R(i,i);
    };
  };

};


};

/**
 * Performs the QR decomposition on a matrix, using Householder reflections approach.
 *
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
 * Computes the pseudo-inverse of a matrix via Householder reflections.
 *
 * \param A real rectangular matrix with row-count >= column-count.
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
  if(A.get_row_count() < A.get_col_count())
    throw std::range_error("Pseudo-inverse with QR is only possible on a matrix with row-count >= column-count!");

  detail::linlsq_QR_impl(A,A_pinv,mat<typename mat_traits<Matrix1>::value_type,mat_structure::identity>(A.get_row_count()),NumTol);
};









#if 0


/**
 * Computes the eigen-values / -vectors via QR decomposition of a matrix, using Householder reflections approach.
 *
 * \param A real square matrix.
 * \param E holds the unsorted eigenvalue vector on its diagonal.
 * \param Q holds as output, the eigenvectors corresponding to the list of eigenvalues in the diagonal E.
 * \param aMaxIter defines the maximum number of iterations to perform.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws maximum_iteration if the method has reached maximum iteration aMaxIter without converging.
 * \throws std::range_error if the matrix A is not square.
 *
 * \author Mikael Persson
 */
template <typename Matrix1,typename Matrix2, typename Matrix3>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value &&
                             is_writable_matrix<Matrix2>::value &&
                             is_writable_matrix<Matrix3>::value ,
void >::type eigensolve_QR(const Matrix1& A, Matrix2& E, Matrix3& Q, unsigned int aMaxIter, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(A.get_col_count() != A.get_row_count())
    throw std::range_error("A matrix must be square for its eigen problem to be solved!");
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  using std::fabs;
  
  mat<ValueType,mat_structure::square> A_tmp(A);
  mat<ValueType,mat_structure::square> R( mat<ValueType,mat_structure::identity>(A.get_col_count()) );
  mat<ValueType,mat_structure::square> Qk( mat<ValueType,mat_structure::identity>(A.get_col_count()) );
  mat<ValueType,mat_structure::square> Q_tmp( mat<ValueType,mat_structure::identity>(A.get_col_count()) );
  
  SizeType N = A_tmp.get_row_count();
  
  RK_NOTICE(2,"A_tmp -1 is " << A_tmp);
  
  unsigned int c(0);
  do {
    bool is_reduced = false;
    for(SizeType i = N-1;i>0;--i) {
      if(fabs(A_tmp(i,i-1)) < NumTol * (fabs(A_tmp(i-1,i-1)) + fabs(A_tmp(i,i)))) {
	mat<ValueType,mat_structure::square> E_11(i);
	mat<ValueType,mat_structure::square> Q_11(i);
	if(i > 1)
	  eigensolve_QR(mat_const_sub_block< mat<ValueType,mat_structure::square> >(A_tmp,i,i,0,0),E_11,Q_11,aMaxIter,NumTol);
	else {
	  E_11(0,0) = A_tmp(0,0);
	  Q_11(0,0) = ValueType(1.0);
	};
	mat<ValueType,mat_structure::square> E_22(N-i);
	mat<ValueType,mat_structure::square> Q_22(N-i);
	if(N-i > 1)
	  eigensolve_QR(mat_const_sub_block< mat<ValueType,mat_structure::square> >(A_tmp,N-i,N-i,i,i),E_22,Q_22,aMaxIter,NumTol);
	else {
	  E_22(0,0) = A_tmp(N-1,N-1);
	  Q_22(0,0) = ValueType(1.0);
	};
	set_block(A_tmp,transpose(Q_11) * mat_const_sub_block< mat<ValueType,mat_structure::square> >(A_tmp,i,N-i,0,i) * Q_22, 0, i);
	set_block(A_tmp,E_11,0,0);
	set_block(A_tmp,E_22,i,i);
	
	mat<ValueType,mat_structure::nil> Q_12(i,N-i);
	mat<ValueType,mat_structure::nil> Q_21(N-i,i);
	Q_tmp *= mat_const_ref_horiz_cat< 
	           mat_const_ref_vert_cat< mat<ValueType,mat_structure::square>, mat<ValueType,mat_structure::nil> >,
		   mat_const_ref_vert_cat< mat<ValueType,mat_structure::nil>, mat<ValueType,mat_structure::square> > >(
		     mat_const_ref_vert_cat< mat<ValueType,mat_structure::square>, mat<ValueType,mat_structure::nil> >(Q_11,Q_21),
		     mat_const_ref_vert_cat< mat<ValueType,mat_structure::nil>, mat<ValueType,mat_structure::square> >(Q_12,Q_22)
		  );
	
	is_reduced = true;
	break;
      };
    };
    
    if(!is_reduced) {
      ValueType mu = A_tmp(N-1,N-1);
      detail::decompose_QR_impl(A_tmp - mat<double,mat_structure::diagonal>(N,mu),Qk,R,NumTol);
      A_tmp = R * Qk + mat<double,mat_structure::diagonal>(N,mu);
      Q_tmp *= Qk;
    };
    
    RK_NOTICE(2,"A_tmp " << c << " is " << A_tmp);

    if (++c == aMaxIter)
      throw maximum_iteration(aMaxIter);

  } while (!is_diagonal(A_tmp,NumTol));
  Q = Q_tmp;
  E = A_tmp;
};

#endif


  
  
};

#endif












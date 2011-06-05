
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



template <typename Matrix1, typename Matrix2, typename Matrix3>
void decompose_QR_impl(const Matrix1& A, Matrix2& Q, Matrix3& R, typename mat_traits<Matrix1>::value_type NumTol)
{
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  ValueType alpha;
  SizeType N = A.get_row_count();
  SizeType M = A.get_col_count();
  mat<ValueType,mat_structure::square> Qt; Qt = mat<ValueType,mat_structure::identity>(N);
  mat<ValueType,mat_structure::square> M_tmp; M_tmp = mat<ValueType,mat_structure::identity>(N);
  mat<ValueType,mat_structure::square> Q_store; Q_store = mat<ValueType,mat_structure::identity>(N);
  mat<ValueType,mat_structure::rectangular> R_store(A);
  vect_n<ValueType> u(N);
  ValueType u_norm;
  SizeType t = (N-1 > M ? M : N-1);

  for(SizeType i=0;i<t;++i) {
    //Compute alpha
    alpha = 0.0;
    for(SizeType j=i;j<N;++j)
      alpha += R_store(j,i) * R_store(j,i);
    alpha = -sqrt(alpha);

    //Compute u & ||u||^2
    u[i] = R_store(i,i) - alpha;
    u_norm = u[i] * u[i];
    for(SizeType j=i+1;j<N;++j) {
      u[j] = R_store(j,i);
      u_norm += u[j] * u[j];
    };

    if(u_norm > NumTol) {
      //Compute Qt_i+1
      for(SizeType j=0;j<i;++j) {
        for(SizeType k=j+1;k<N;++k) {
	  Qt(j,k) = (Qt(k,j) = 0.0);
        };
        Qt(j,j) = 1.0;
      };
      for(SizeType j=i;j<N;++j) {
        for(SizeType k=j+1;k<N;++k) {
  	  Qt(j,k) = (Qt(k,j) = -2.0 * u[j] * u[k] / u_norm);
        };
        Qt(j,j) = 1.0 - 2.0 * u[j] * u[j] / u_norm;
      };

      //Compute R_i+1
      for(SizeType j=i;j<N;++j) { //for every row of Qt below i
        for(SizeType k=i;k<M;++k) { //for every column in R_store right of i
	  M_tmp(j,k) = 0.0;
	  for(SizeType l=i;l<N;++l) { //for every element in Qt rows and every element in R_store columns
	    M_tmp(j,k) += Qt(j,l) * R_store(l,k);
	  };
        };
      };
      for(SizeType k=i;k<M;++k) { //for every column in R right of i
        for(SizeType j=i;j<N;++j) { //for every row of R below i
          R_store(j,k) = M_tmp(j,k);
        };
      };

      //Compute Q_i+1, i.e. Qtt = Qt
      for(SizeType j=i;j<N;++j) { //for every row of Qt below i
        for(SizeType k=0;k<N;++k) { //for every column in Q_store
  	  M_tmp(j,k) = 0.0;
	  for(SizeType l=i;l<N;++l) { //for every element in Qt rows and every element in Q_store columns
	    M_tmp(j,k) += Qt(j,l) * Q_store(l,k);
	  };
        };
      };
      for(SizeType k=0;k<N;++k) { //for every column in Q_store
        for(SizeType j=i;j<N;++j) { //for every row of Q_store below i
          Q_store(j,k) = M_tmp(j,k);
        };
      };
    };
  };

  //Finally, copy the results.
  for(SizeType i=0;i<M;++i) {
    for(SizeType j=0;j<N;++j)
      Q(j,i) = Q_store(i,j); //note that this is a transpose because Q_store is a compound of non-transposed Qt multiplications.
    for(SizeType j=0;j<M;++j)
      R(j,i) = R_store(j,i);
  };
 
  
};
  


template <typename Matrix1, typename Matrix2, typename Matrix3>
void linlsq_QR_impl(const Matrix1& A, Matrix2& x,const Matrix3& b, typename mat_traits<Matrix1>::value_type NumTol)
{
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  using std::sqrt;
  
  ValueType alpha;
  SizeType N = A.get_row_count();
  SizeType M = A.get_col_count();
  mat<ValueType, mat_structure::square> Q_t; Q_t = mat<ValueType,mat_structure::identity>(N);
  mat<ValueType, mat_structure::square> M_tmp; M_tmp = mat<ValueType,mat_structure::identity>(N);
  mat<ValueType, mat_structure::square> Q_store; Q_store = mat<ValueType,mat_structure::identity>(N);
  mat<ValueType, mat_structure::rectangular> R_store(A);
  mat<ValueType, mat_structure::rectangular> b_store(b);
  mat<ValueType, mat_structure::rectangular> b_tmp(b);
  vect_n<ValueType> u(N);
  ValueType u_norm;
  SizeType t = (N-1 > M ? M : N-1);

  for(SizeType i=0;i<t;++i) {
    //Compute alpha
    alpha = 0.0;
    for(SizeType j=i;j<N;++j)
      alpha += R_store(j,i) * R_store(j,i);
    alpha = -sqrt(alpha);

    //Compute u & ||u||^2
    u[i] = R_store(i,i) - alpha;
    u_norm = u[i] * u[i];
    for(SizeType j=i+1;j<N;++j) {
      u[j] = R_store(j,i);
      u_norm += u[j] * u[j];
    };

    if(u_norm > NumTol) {
      //Compute Qt_i+1
      for(SizeType j=0;j<i;++j) {
        for(SizeType k=j+1;k<N;++k)
          Q_t(j,k) = (Q_t(k,j) = 0.0);
        Q_t(j,j) = 1.0;
      };
      for(SizeType j=i;j<N;++j) {
        for(SizeType k=j+1;k<N;++k)
          Q_t(j,k) = (Q_t(k,j) = -2.0 * u[j] * u[k] / u_norm);
        Q_t(j,j) = 1.0 - 2.0 * u[j] * u[j] / u_norm;
      };

      //Compute R_i+1
      for(SizeType j=i;j<N;++j) { //for every row of Qt below i
        for(SizeType k=i;k<M;++k) { //for every column in R_store right of i
          M_tmp(j,k) = 0.0;
          for(SizeType l=i;l<N;++l) //for every element in Qt rows and every element in R_store columns
            M_tmp(j,k) += Q_t(j,l) * R_store(l,k);
        };
      };
      for(SizeType k=i;k<M;++k) { //for every column in R right of i
        for(SizeType j=i;j<N;++j) //for every row of R below i
          R_store(j,k) = M_tmp(j,k);
      };

      //Compute b_i+1
      for(SizeType j=i;j<N;++j) { // for every row of Qt below i
        for(SizeType k=0;k<b.get_col_count();++k) { //for every column of b_store
          b_tmp(j,k) = 0.0;
          for(SizeType l=i;l<N;++l) //for every non-zero element in Qt row j and every element in b_store column k
            b_tmp(j,k) += Q_t(j,l) * b_store(l,k);
        };
      };
      for(SizeType k=0;k<b.get_col_count();++k)
        for(SizeType l=i;l<N;++l)
          b_store(l,k) = b_tmp(l,k);
    };
  };

  //back-substitution
  x.set_row_count(M);
  x.set_col_count(b.get_col_count());
  for(int i=M-1;i>=0;--i) {
    for(SizeType j=0;j<b.get_col_count();++j) {
      ValueType sum = b_store(i,j);
      for(SizeType k=i+1;k<M;++k)
        sum -= x(k,j) * R_store(i,k);
      x(i,j) = sum / R_store(i,i);
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

  if((A.get_row_count() != Q.get_row_count()) || (A.get_col_count() != Q.get_col_count())) {
    Q.set_row_count(A.get_row_count());
    Q.set_col_count(A.get_col_count());
  };
  R.set_row_count(A.get_col_count());
  R.set_col_count(A.get_col_count());

  detail::decompose_QR_impl(A,Q,R,NumTol);
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

  if((A.get_row_count() != Q.get_row_count()) || (A.get_col_count() != Q.get_col_count())) {
    Q.set_row_count(A.get_row_count());
    Q.set_col_count(A.get_col_count());
  };
  mat<typename mat_traits<Matrix3>::value_type, mat_structure::square> R_tmp(A.get_col_count());

  detail::decompose_QR_impl(A,Q,R_tmp,NumTol);
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
  
  mat<ValueType,mat_structure::rectangular> Q(A.get_row_count(),A.get_row_count());
  mat<ValueType,mat_structure::rectangular> R(A.get_row_count(),A.get_row_count());
  detail::decompose_QR_impl(A,Q,R,NumTol);

  ValueType result(1.0);
  for(SizeType i=0;i<A.get_row_count();++i)
    result *= R(i,i);
  return result;
};



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
  
  mat<ValueType,mat_structure::square> A_tmp(A);
  mat<ValueType,mat_structure::square> R( mat<ValueType,mat_structure::identity>(A.get_col_count()) );
  mat<ValueType,mat_structure::square> Qk( mat<ValueType,mat_structure::identity>(A.get_col_count()) );
  mat<ValueType,mat_structure::square> Q_tmp( mat<ValueType,mat_structure::identity>(A.get_col_count()) );
  
  unsigned int c(0);
  do {
    detail::decompose_QR_impl(A_tmp,Qk,R,NumTol);
    A_tmp = R * Qk;
    Q_tmp *= Qk;

    if (++c == aMaxIter)
      throw maximum_iteration(aMaxIter);

  } while (!is_diagonal(A_tmp,NumTol));
  Q = Q_tmp;
  E = A_tmp;
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


  
  
};

#endif












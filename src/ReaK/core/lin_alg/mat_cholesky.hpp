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

  typedef typename mat_traits<Matrix1>::value_type ValueType;
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

  typedef typename mat_traits<Matrix1>::value_type ValueType;
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





#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

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




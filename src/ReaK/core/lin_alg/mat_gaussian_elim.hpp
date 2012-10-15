/**
 * \file mat_gaussian_elim.hpp
 * 
 * This library provides a number of functions related to performing a Gaussian elimination on a 
 * matrix, e.g., to invert a matrix and to solve a linear system. Most implementations provided 
 * are based on the PLU decomposition (Crout's method). PLU decomposition is pretty efficient and 
 * is generally preferred if there is good reasons to believe that the matrix involved will always
 * be well-conditioned. If a matrix cannot be guaranteed to be well-conditioned, algorithms such as 
 * QR-decomposition or methods specific to symmetric matrices are preferred (Cholesky), or even SVD.
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

#ifndef REAK_MAT_GAUSSIAN_ELIM_HPP
#define REAK_MAT_GAUSSIAN_ELIM_HPP

#include "mat_alg.hpp"
#include "mat_num_exceptions.hpp"

namespace ReaK {



/*************************************************************************
                        Gaussian Elimination
*************************************************************************/

/**
 * Inverts a matrix using the Gauss-Jordan elimination on the identity matrix.
 * \note that PLU decomposition or any other method is faster for matrix sizes of more than 20x20.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A fully-writable matrix type.
 * \param A A well-conditioned, square (Size x Size), real, full-rank matrix to be inverted.
 * \param A_inv The matrix which stores, as output, the inverse of A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero.
 *
 * \throws singularity_error if the matrix A is numerically singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not a square matrix.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && is_fully_writable_matrix<Matrix2>::value,
void >::type invert_gaussian(const Matrix1& A, Matrix2& A_inv, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  
  using std::fabs;
   
  if(A.get_col_count() != A.get_row_count())
    throw std::range_error("Inversion impossible! Matrix A is not square!");
  
  SizeType Size = A.get_col_count();
  Matrix2 tmp(A);
  A_inv = mat<ValueType,mat_structure::identity>(Size);

  for(SizeType i=0;i<Size;++i)
  {
    if(fabs(tmp(i,i)) < NumTol) {
      for(SizeType j=i+1;j<=Size;++j)
      {
        if(j == Size) {
          throw singularity_error("M");
        };
        if(fabs(tmp(j,i)) > NumTol) {
          for(SizeType k=i;k<Size;++k)
	    tmp(i,k) += tmp(j,k);
	  for(SizeType k=0;k<Size;++k)
            A_inv(i,k) += A_inv(j,k);
          break;
        };
      };
    };

    ValueType s = tmp(i,i);
    for(SizeType k=i;k<Size;++k)
      tmp(i,k) /= s;
    for(SizeType k=0;k<Size;++k)
      A_inv(i,k) /= s;

    for(SizeType j=0;j<Size;++j)
    {
      if (i != j) {
        s = tmp(j,i);
        for(SizeType k=i;k<Size;++k)
          tmp(j,k) -= s*tmp(i,k);
	for(SizeType k=0;k<Size;++k)
	  A_inv(j,k) -= s*A_inv(i,k);
      };
    };
  };
};


/**
 * Inverts a matrix using the Gauss-Jordan elimination on the identity matrix.
 * \note that PLU decomposition or any other method is faster for matrix sizes of more than 20x20.
 * 
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 Any writable matrix type.
 * \param A well-conditioned, square (Size x Size), real, full-rank matrix to be inverted.
 * \param A_inv The matrix which stores, as output, the inverse of A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero.
 *
 * \throws singularity_error if the matrix A is numerically singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not a square matrix.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_readable_matrix<Matrix1>::value && !is_fully_writable_matrix<Matrix2>::value,
void >::type invert_gaussian(const Matrix1& A, Matrix2& A_inv, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix2>::value_type ValueType2;
  mat<ValueType2,mat_structure::rectangular> A_inv_tmp(A_inv.get_row_count(),A_inv.get_col_count());
  invert_gaussian(A,A_inv_tmp,NumTol);
  A_inv = A_inv_tmp;
};





/*************************************************************************
             Permuted Lower- / Upper-triangular Decomposition
*************************************************************************/


namespace detail {
  
template <typename Matrix1, typename Matrix2, typename IndexVector>
void linsolve_PLU_impl(Matrix1& A, Matrix2& b, IndexVector& P, typename mat_traits<Matrix1>::value_type NumTol) {
  using std::swap;
  using std::fabs;
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  typedef typename mat_traits<Matrix1>::size_type SizeType;
  
  mat<ValueType,mat_structure::rectangular,mat_alignment::column_major> s(b.get_row_count(),b.get_col_count());
  SizeType An = A.get_row_count();
  SizeType bn = b.get_col_count();
  
  for(SizeType i=0;i<An;++i) {
    P[i] = i;
    s(i,0) = 0.0;
    for(SizeType j=0;j<An;++j)
      if(s(i,0) < fabs(A(i,j)))
	s(i,0) = fabs(A(i,j));
  };

  for(SizeType k=0;k<An;++k) {
    
    for(SizeType i=k;i<An;++i)
      for(SizeType j=0;j<k;++j)
        A(i,k) -= A(i,j)*A(j,k);

    SizeType temp_i=k;
    for(SizeType i=k+1;i<An;++i)
      if(fabs(A(i,k) / s(i,0)) > fabs(A(temp_i,k) / s(temp_i,0)))
	temp_i = i;

    if(k != temp_i) {
      for(SizeType i=0;i<An;++i)
	swap(A(k,i),A(temp_i,i));
      swap(s(k,0), s(temp_i,0));
      swap(P[k], P[temp_i]);
    };

    for(SizeType j=k+1;j<An;++j) {
      for(SizeType i=0;i<k;++i)
        A(k,j) -= A(k,i) * A(i,j);
      if (fabs(A(k,k)) < NumTol)
	throw singularity_error("A");
      A(k,j) /= A(k,k);
    };
  };

  // Back-substitution
  for(SizeType k=0;k<An;++k) {
    for(SizeType l=0;l<bn;++l)
      s(P[k],l) = b(P[k],l);

    for(SizeType l=0;l<bn;++l) {
      for(SizeType j=0;j<k;++j)
        s(k,l) -= A(k,j) * s(j,l);
      s(k,l) /= A(k,k);
    };
  };

  b = s;

  for(int k=An-1;k>=0;--k)
    for(SizeType l=0;l<bn;++l)
      for(SizeType j=k+1;j<An;++j)
        b(k,l) -= A(k,j) * b(j,l);
 
  return;
};


template <typename Matrix1, typename Matrix2, typename IndexVector>
typename boost::enable_if_c< is_fully_writable_matrix< Matrix1 >::value && 
                             is_fully_writable_matrix< Matrix2 >::value, 
void >::type linsolve_PLU_dispatch(Matrix1& A, Matrix2& b, IndexVector& P, typename mat_traits<Matrix1>::value_type NumTol) {
  linsolve_PLU_impl(A,b,P,NumTol);
};

template <typename Matrix1, typename Matrix2, typename IndexVector>
typename boost::enable_if_c< !is_fully_writable_matrix< Matrix1 >::value && 
                             is_fully_writable_matrix< Matrix2 >::value, 
void >::type linsolve_PLU_dispatch(Matrix1& A, Matrix2& b, IndexVector& P, typename mat_traits<Matrix1>::value_type NumTol) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  mat<ValueType,mat_structure::square> A_tmp(A);
  linsolve_PLU_impl(A_tmp,b,P,NumTol);
  A = A_tmp;
};

template <typename Matrix1, typename Matrix2, typename IndexVector>
typename boost::enable_if_c< !is_fully_writable_matrix< Matrix1 >::value && 
                             !is_fully_writable_matrix< Matrix2 >::value, 
void >::type linsolve_PLU_dispatch(Matrix1& A, Matrix2& b, IndexVector& P, typename mat_traits<Matrix1>::value_type NumTol) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  mat<ValueType,mat_structure::square> A_tmp(A);
  typedef typename mat_traits<Matrix2>::value_type ValueType2;
  mat<ValueType2,mat_structure::rectangular> b_tmp(b);
  linsolve_PLU_impl(A_tmp,b_tmp,P,NumTol);
  A = A_tmp;
  b = b_tmp;
};

template <typename Matrix1, typename Matrix2, typename IndexVector>
typename boost::enable_if_c< is_fully_writable_matrix< Matrix1 >::value && 
                             !is_fully_writable_matrix< Matrix2 >::value, 
void >::type linsolve_PLU_dispatch(Matrix1& A, Matrix2& b, IndexVector& P, typename mat_traits<Matrix1>::value_type NumTol) {
  typedef typename mat_traits<Matrix2>::value_type ValueType2;
  mat<ValueType2,mat_structure::rectangular> b_tmp(b);
  linsolve_PLU_impl(A,b_tmp,P,NumTol);
  b = b_tmp;
};

  
};




/**
 * Solves the linear problem AX = B using PLU decomposition as defined by Crout`s method.
 * \note To solve a linear system involving vectors for X and B, use the Matrix-Vector Adaptors (mat_vector_adaptor.hpp).
 *
 * \tparam Matrix1 A writable matrix type.
 * \tparam Matrix2 A writable matrix type.
 * \tparam IndexVector A writable vector type.
 * \param A well-conditioned, square (Size x Size), real, full-rank matrix which multiplies x.
 *          As output, A stores the LU decomposition, permutated by P.
 * \param b stores, as input, the RHS of the linear system of equation and stores, as output,
 *          the solution matrix X (Size x B_ColCount).
 * \param P vector of Size unsigned integer elements holding, as output, the permutations done
 *          the rows of matrix A during the decomposition to LU.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero.
 *
 * \throws singularity_error if the matrix A is numerically singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename IndexVector>
typename boost::enable_if_c< is_writable_matrix< Matrix1 >::value && 
                             is_writable_matrix< Matrix2 >::value &&
                             is_writable_vector< IndexVector >::value, 
void >::type linsolve_PLU(Matrix1& A, Matrix2& b, IndexVector& P, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  if(A.get_col_count() != A.get_row_count())
    throw std::range_error("PLU decomposition impossible! Matrix A is not square!");
  if(b.get_row_count() != A.get_col_count())
    throw std::range_error("PLU decomposition impossible! Matrix b must have same row count as A!");

  P.resize(A.get_col_count());
  detail::linsolve_PLU_dispatch(A,b,P,NumTol);
};


/**
 * Solves the linear problem AX = B using PLU decomposition as defined by Crout`s method.
 * \note To solve a linear system involving vectors for X and B, use the Matrix-Vector Adaptors (mat_vector_adaptor.hpp).
 * 
 * \tparam Matrix1 A writable matrix type.
 * \tparam Matrix2 A writable matrix type.
 * \param A well-conditioned, square (Size x Size), real, full-rank matrix which multiplies x.
 *          As output, A stores the LU decomposition, permutated by P.
 * \param b stores, as input, the RHS of the linear system of equation and stores, as output,
 *          the solution matrix X (Size x B_ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero.
 *
 * \throws singularity_error if the matrix A is numerically singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2, typename IndexVector>
typename boost::enable_if_c< is_writable_matrix< Matrix1 >::value && 
                             is_writable_matrix< Matrix2 >::value, 
void >::type linsolve_PLU(Matrix1& A, Matrix2& b, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  vect_n<unsigned int> P;
  linsolve_PLU(A,b,P,NumTol);
};

/**
 * Solves the linear problem AX = B using PLU decomposition as defined by Crout`s method.
 *
 * \tparam Matrix A writable matrix type.
 * \tparam Vector A writable vector type.
 * \tparam IndexVector A writable vector type.
 * \param A well-conditioned, square (Size x Size), real, full-rank matrix which multiplies x.
 *          As output, A stores the LU decomposition, permutated by P.
 * \param b stores, as input, the RHS of the linear system of equation and stores, as output,
 *          the solution vector X.
 * \param P vector of Size unsigned integer elements holding, as output, the permutations done
 *          the rows of matrix A during the decomposition to LU.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero.
 *
 * \throws singularity_error if the matrix A is numerically singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \author Mikael Persson
 */
template <typename Matrix, typename Vector, typename IndexVector>
typename boost::enable_if_c< is_writable_matrix< Matrix >::value && 
                             is_writable_vector< Vector >::value &&
                             is_writable_vector< IndexVector >::value, 
void >::type linsolve_PLU(Matrix& A, Vector& b, IndexVector& P, typename mat_traits<Matrix>::value_type NumTol = 1E-8) {
  if(A.get_col_count() != A.get_row_count())
    throw std::range_error("PLU decomposition impossible! Matrix A is not square!");
  if(b.size() != A.get_col_count())
    throw std::range_error("PLU decomposition impossible! Matrix b must have same row count as A!");

  P.resize(A.get_col_count());
  mat_vect_adaptor<Vector,mat_alignment::column_major> b_mat(b);
  detail::linsolve_PLU_dispatch(A,b_mat,P,NumTol);
};

/**
 * Solves the linear problem AX = B using PLU decomposition as defined by Crout`s method.
 *
 * \tparam Matrix A writable matrix type.
 * \tparam Vector A writable vector type.
 * \param A well-conditioned, square (Size x Size), real, full-rank matrix which multiplies x.
 *          As output, A stores the LU decomposition, permutated by P.
 * \param b stores, as input, the RHS of the linear system of equation and stores, as output,
 *          the solution vector X.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero.
 *
 * \throws singularity_error if the matrix A is numerically singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not square or if b's row count does not match that of A.
 *
 * \author Mikael Persson
 */
template <typename Matrix, typename Vector, typename IndexVector>
typename boost::enable_if_c< is_writable_matrix< Matrix >::value && 
                             is_writable_vector< Vector >::value, 
void >::type linsolve_PLU(Matrix& A, Vector& b, typename mat_traits<Matrix>::value_type NumTol = 1E-8) {
  vect_n<unsigned int> P;
  linsolve_PLU(A,b,P,NumTol);
};


/**
 * Functor to wrap a call to a PLU decomposition-based linear system solver.
 */
struct PLU_linsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
    typedef typename mat_traits<Matrix1>::value_type ValueType;
    mat<ValueType, mat_structure::square> A_tmp(A);
    X = B;
    linsolve_PLU(A_tmp,X,NumTol);
  };
};


/**
 * Inverts a matrix using PLU decomposition as defined by Crout`s method.
 *
 * \tparam Matrix1 A readable matrix type.
 * \tparam Matrix2 A writable matrix type.
 * \param A well-conditioned, square (Size x Size), real, full-rank matrix to be inverted.
 * \param A_inv The matrix which stores, as output, the inverse of A.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero.
 *
 * \throws singularity_error if the matrix A is numerically singular (or rank-deficient).
 * \throws std::range_error if the matrix A is not a square matrix.
 *
 * \author Mikael Persson
 */
template <typename Matrix1, typename Matrix2>
typename boost::enable_if_c< is_writable_matrix< Matrix1 >::value && 
                             is_writable_matrix< Matrix2 >::value, 
void >::type invert_PLU(Matrix1 A,Matrix2& A_inv, typename mat_traits<Matrix1>::value_type NumTol = 1E-8) {
  typedef typename mat_traits<Matrix1>::value_type ValueType;
  
  if(A.get_col_count() != A.get_row_count())
    throw std::range_error("PLU decomposition impossible! Matrix A is not square!");

  A_inv = mat<ValueType,mat_structure::identity>(A.get_col_count());
  vect_n<unsigned int> P(A.get_col_count());
  detail::linsolve_PLU_dispatch(A,A_inv,P,NumTol);
};










#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

extern template void invert_gaussian(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& A_inv, double NumTol);
extern template void invert_gaussian(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& A_inv, double NumTol);
extern template void invert_gaussian(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& A_inv, double NumTol);
extern template void invert_gaussian(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& A_inv, double NumTol);

extern template void linsolve_PLU(mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& b, vect_n<unsigned int>& P, double NumTol);
extern template void linsolve_PLU(mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& b, vect_n<unsigned int>& P, double NumTol);
extern template void linsolve_PLU(mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& b, vect_n<unsigned int>& P, double NumTol);
extern template void linsolve_PLU(mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& b, vect_n<unsigned int>& P, double NumTol);

extern template void linsolve_PLU(mat<double,mat_structure::rectangular>& A, vect_n<double>& b, vect_n<unsigned int>& P, double NumTol);
extern template void linsolve_PLU(mat<double,mat_structure::square>& A, vect_n<double>& b, vect_n<unsigned int>& P, double NumTol);

extern template void invert_PLU(mat<double,mat_structure::rectangular> A, mat<double,mat_structure::rectangular>& A_inv, double NumTol);
extern template void invert_PLU(mat<double,mat_structure::rectangular> A, mat<double,mat_structure::square>& A_inv, double NumTol);
extern template void invert_PLU(mat<double,mat_structure::square> A, mat<double,mat_structure::rectangular>& A_inv, double NumTol);
extern template void invert_PLU(mat<double,mat_structure::square> A, mat<double,mat_structure::square>& A_inv, double NumTol);


extern template void invert_gaussian(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& A_inv, float NumTol);
extern template void invert_gaussian(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& A_inv, float NumTol);
extern template void invert_gaussian(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& A_inv, float NumTol);
extern template void invert_gaussian(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& A_inv, float NumTol);

extern template void linsolve_PLU(mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& b, vect_n<unsigned int>& P, float NumTol);
extern template void linsolve_PLU(mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& b, vect_n<unsigned int>& P, float NumTol);
extern template void linsolve_PLU(mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& b, vect_n<unsigned int>& P, float NumTol);
extern template void linsolve_PLU(mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& b, vect_n<unsigned int>& P, float NumTol);

extern template void linsolve_PLU(mat<float,mat_structure::rectangular>& A, vect_n<float>& b, vect_n<unsigned int>& P, float NumTol);
extern template void linsolve_PLU(mat<float,mat_structure::square>& A, vect_n<float>& b, vect_n<unsigned int>& P, float NumTol);

extern template void invert_PLU(mat<float,mat_structure::rectangular> A, mat<float,mat_structure::rectangular>& A_inv, float NumTol);
extern template void invert_PLU(mat<float,mat_structure::rectangular> A, mat<float,mat_structure::square>& A_inv, float NumTol);
extern template void invert_PLU(mat<float,mat_structure::square> A, mat<float,mat_structure::rectangular>& A_inv, float NumTol);
extern template void invert_PLU(mat<float,mat_structure::square> A, mat<float,mat_structure::square>& A_inv, float NumTol);

#endif










};


#endif










/**
 * \file quadratic_programs.hpp
 *
 * The following library provides implementations of quadratic programming algorithms. 
 * The algorithm follows the specification given by Nocedal's Numerical Optimization book.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
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

#ifndef REAK_QUADRATIC_PROGRAMS_HPP
#define REAK_QUADRATIC_PROGRAMS_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/math/lin_alg/mat_qr_decomp.hpp>
#include <ReaK/math/lin_alg/mat_views.hpp>

#include "optim_exceptions.hpp"

#include <vector>

namespace ReaK {
  
  
namespace optim {
  


/**
 * This function is an implementation of the null-space direct 
 * method for solving a quadratic optimization problem with equality constraints. 
 * It solves the following problem: \n
 * \n
 *           min c'x + 0.5 * x' G x \n
 *               Ax = b \n
 * \n
 * The implementation was inspired from the algorithm described in the book:\n
 *   Nocedal, Numerical Optimization, 2nd Ed..
 * TEST PASSED
 * 
 * \tparam Matrix A general matrix type, should model the WritableMatrixConcept (and be fully-writable).
 * \tparam Vector1 A vector type, should model the WritableVectorConcept.
 * \tparam Vector2 A vector type, should model the WritableVectorConcept.
 * \param A The constraint matrix of dimension M*N.
 * \param b The b vector of dimension M.
 * \param G The G matrix of dimension NxN (defines the quadratic function to minimize, should be positive semi-definite).
 * \param c The cost vector of dimension N.
 * \param x Stores, as output, the optimal vector.
 * \param abs_tol The tolerance on the singularity of components of the matrices involved.
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Vector1, typename Matrix2, typename Vector2>
void null_space_QP_method(const Matrix1& A, const Vector1& b, 
                         const Matrix2&  G, const Vector2& c, Vector2& x,
                         typename vect_traits<Vector1>::value_type abs_tol = std::numeric_limits<typename vect_traits<Vector1>::value_type>::epsilon(),
                         typename vect_traits<Vector1>::value_type max_norm = std::numeric_limits<typename vect_traits<Vector1>::value_type>::infinity(),
                         Vector1* lambda = NULL) {
  typedef typename vect_traits<Vector1>::value_type ValueType;
  typedef typename vect_traits<Vector1>::size_type SizeType;
  using std::swap;
  using std::fabs;
  
  SizeType N = c.size();
  SizeType M = b.size();
  
  mat<ValueType,mat_structure::rectangular> A_tmp(transpose_view(A));
  mat<ValueType,mat_structure::rectangular> R(N,M);
  mat<ValueType,mat_structure::square> Q(N);
  decompose_QR(A_tmp,Q,R,abs_tol);
  mat<ValueType,mat_structure::rectangular> L = mat<ValueType,mat_structure::rectangular>(transpose_view(R));
  L.set_col_count(M,true);
  
  mat_const_sub_block< mat<ValueType,mat_structure::square> > Y(Q, N, M, 0, 0);
  mat_const_sub_block< mat<ValueType,mat_structure::square> > Z(Q, N, N - M, 0, M);
  //RK_NOTICE(2," Found the null-space basis to be Z = " << Z);
  
  Vector2 x_tmp(x);
  
  //forward-sub with L to find p_y
  mat_vect_adaptor< Vector2 > p_y(x_tmp,M,1,0);
  ValueType L_abstrace(0.0);
  for(SizeType i = 0; i < M; ++i) {
    p_y(i,0) = b[i];
    L_abstrace += fabs(L(i,i));
  };
  ReaK::detail::forwardsub_L_impl(L,p_y,abs_tol);
  
  //RK_NOTICE(2," Found the particular solution to Ax = b as x = " << (Y * p_y) << " Ax = " << A * (Y * p_y));
  
  //solve for p_z:
  mat_vect_adaptor< Vector2 > p_z(x_tmp,N - M,1,M);
  mat<ValueType,mat_structure::rectangular> GY_py = G * Y * p_y;
  for(SizeType i = 0; i < N-M; ++i) {
    p_z(i,0) = ValueType(0.0);
    for(SizeType j = 0; j < N; ++j)
      p_z(i,0) -= Z(j,i) * (c[j] + GY_py(j,0));
  };
  mat<ValueType, mat_structure::symmetric> ZGZ = mat<ValueType, mat_structure::symmetric>(transpose_view(Z) * G * Z);

  try {
    linsolve_Cholesky(ZGZ, p_z, abs_tol);
  } catch(singularity_error&) {
    for(SizeType i = 0; i < N-M; ++i) {
      p_z(i,0) = ValueType(0.0);
      for(SizeType j = 0; j < N; ++j)
        p_z(i,0) -= Z(j,i) * (c[j] + GY_py(j,0));
    };
    ValueType zz(0.0);
    ValueType zBz(0.0);
    for(SizeType i = 0; i < N-M; ++i) {
      zz += p_z(i,0) * p_z(i,0);
      for(SizeType j = 0; j < N-M; ++j)
        zBz += p_z(i,0) * (ZGZ(i,j) * p_z(j,0));
    };
    zBz = zz / (fabs(zBz) + abs_tol);
    if(zBz < ValueType(1.0)) {
      zz *= zBz * zBz;
      for(SizeType i = 0; i < N-M; ++i)
        p_z(i,0) = p_z(i,0) * zBz;
    };
    if( zz > max_norm * max_norm ) {
      zBz = max_norm / fabs(zz);
      for(SizeType i = 0; i < N-M; ++i)
        p_z(i,0) = p_z(i,0) * zBz;
    };
  };
  
  x = Q * x_tmp;
  
  if(lambda) {
    (*lambda) = (c + G * x) * Y;
    mat_vect_adaptor< Vector1 > lambda_mat(*lambda);
    ReaK::detail::backsub_R_impl(transpose_view(L),lambda_mat,abs_tol);
  };
  
};


/**
 * This function is an implementation of the null-space direct 
 * method for solving a quadratic optimization problem with equality constraints. 
 * It solves the following problem: \n
 * \n
 *           min c'x + 0.5 * x' G x \n
 *               Ax = b \n
 * \n
 * The implementation was inspired from the algorithm described in the book:\n
 *   Nocedal, Numerical Optimization, 2nd Ed..
 * TEST PASSED
 * 
 * \tparam Matrix A general matrix type, should model the WritableMatrixConcept (and be fully-writable).
 * \tparam Vector1 A vector type, should model the WritableVectorConcept.
 * \tparam Vector2 A vector type, should model the WritableVectorConcept.
 * \param A The constraint matrix of dimension M*N.
 * \param b The b vector of dimension M.
 * \param G The G matrix of dimension NxN (defines the quadratic function to minimize, should be positive semi-definite).
 * \param c The cost vector of dimension N.
 * \param x Stores, as output, the optimal vector.
 * \param abs_tol The tolerance on the singularity of components of the matrices involved.
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Vector1, typename Matrix2, typename Vector2>
void null_space_RRQP_method(const Matrix1& A, const Vector1& b, 
                            const Matrix2&  G, const Vector2& c, Vector2& x,
                            typename vect_traits<Vector1>::value_type abs_tol = std::numeric_limits<typename vect_traits<Vector1>::value_type>::epsilon(),
                            typename vect_traits<Vector1>::value_type max_norm = std::numeric_limits<typename vect_traits<Vector1>::value_type>::infinity(),
                            Vector1* lambda = NULL) {
  typedef typename vect_traits<Vector1>::value_type ValueType;
  typedef typename vect_traits<Vector1>::size_type SizeType;
  using std::swap;
  using std::fabs;
  
  SizeType N = c.size();
  SizeType M = b.size();
  
  mat<ValueType,mat_structure::rectangular> R(transpose_view(A));
  mat<ValueType,mat_structure::square> Q = mat<ValueType,mat_structure::square>(mat<ValueType,mat_structure::identity>(N));
  mat<ValueType,mat_structure::permutation> P(N);
  SizeType K = ReaK::detail::decompose_RRQR_impl(R,&Q,P,abs_tol);
  mat<ValueType,mat_structure::rectangular> L = mat<ValueType,mat_structure::rectangular>(transpose_view(R));
  L.set_col_count(M,true);
  
  mat_const_sub_block< mat<ValueType,mat_structure::square> > Y(Q, N, K, 0, 0);
  mat_const_sub_block< mat<ValueType,mat_structure::square> > Z(Q, N, N - K, 0, K);
  //RK_NOTICE(2," Found the null-space basis to be Z = " << Z);
  //RK_NOTICE(1," reached! with P.get_col_count() = " << P.get_col_count() << " K = " << K << " N = " << N);
  Vector2 x_tmp(c);
  
  //forward-sub with L to find p_y
  mat_vect_adaptor< Vector2 > p_y(x_tmp,K,1,0);
  for(SizeType i = 0; i < K; ++i)
    p_y(i,0) = b[i];
  //RK_NOTICE(1," reached!");
  ReaK::detail::forwardsub_L_impl(sub(L)(range(0,K),range(0,K)),p_y,abs_tol);
  //RK_NOTICE(1," reached!");
  //RK_NOTICE(2," Found the particular solution to Ax = b as x = " << (Y * p_y) << " Ax = " << A * (Y * p_y));
  
  //solve for p_z:
  mat_vect_adaptor< Vector2 > p_z(x_tmp,N - K,1,K);
  mat<ValueType,mat_structure::rectangular> GY_py = G * Y * p_y;
  for(SizeType i = 0; i < N-K; ++i) {
    p_z(i,0) = ValueType(0.0);
    for(SizeType j = 0; j < N; ++j)
      p_z(i,0) -= Z(j,i) * (c[j] + GY_py(j,0));
  };
  //RK_NOTICE(1," reached!");
  mat<ValueType, mat_structure::symmetric> ZGZ = mat<ValueType, mat_structure::symmetric>(transpose_view(Z) * G * Z);
  //RK_NOTICE(1," reached!");
  try {
    linsolve_Cholesky(ZGZ, p_z, abs_tol);
    //RK_NOTICE(1," reached!");
  } catch(singularity_error&) {
    //RK_NOTICE(1," reached!");
    for(SizeType i = 0; i < N-K; ++i) {
      p_z(i,0) = ValueType(0.0);
      for(SizeType j = 0; j < N; ++j)
        p_z(i,0) -= Z(j,i) * (c[j] + GY_py(j,0));
    };
    //RK_NOTICE(1," reached!");
    ValueType zz(0.0);
    ValueType zBz(0.0);
    for(SizeType i = 0; i < N-K; ++i) {
      zz += p_z(i,0) * p_z(i,0);
      for(SizeType j = 0; j < N-K; ++j)
        zBz += p_z(i,0) * (ZGZ(i,j) * p_z(j,0));
    };
    //RK_NOTICE(1," reached!");
    zBz = zz / (fabs(zBz) + abs_tol);
    if(zBz < ValueType(1.0)) {
      zz *= zBz * zBz;
      for(SizeType i = 0; i < N-K; ++i)
        p_z(i,0) = p_z(i,0) * zBz;
    };
    //RK_NOTICE(1," reached!");
    if( zz > max_norm * max_norm ) {
      zBz = max_norm / fabs(zz);
      for(SizeType i = 0; i < N-K; ++i)
        p_z(i,0) = p_z(i,0) * zBz;
    };
  };
  //RK_NOTICE(1," reached!");
  x = Q * x_tmp;
  //RK_NOTICE(1," reached! with size(x) = " << x.size() << " and P.get_col_count() = " << P.get_col_count() << " size(A) = (" << A.get_row_count() << ", " << A.get_col_count() << ")");
  //RK_NOTICE(1,"   with size(R) = (" << R.get_row_count() << ", " << R.get_col_count() << ") and Q.get_col_count() = " << Q.get_col_count());
  x = P * x;
  //RK_NOTICE(1," reached!");
  if(lambda) {
    vect_n<ValueType> lam( (c + G * x) * Y );
    mat_vect_adaptor< vect_n<ValueType> > lam_mat(lam);
    RK_NOTICE(1," reached!");
    ReaK::detail::backsub_R_impl(sub(transpose_view(L))(range(0,K),range(0,K)),lam_mat,abs_tol);
    RK_NOTICE(1," reached!");
    for(SizeType i = 0; i < K; ++i)
      (*lambda)[i] = lam[i];
    for(SizeType i = K; i < M; ++i)
      (*lambda)[i] = 0.0;
  };
  //RK_NOTICE(1," reached!");
};



/**
 * This function is an implementation of the projected Conjugate Gradient 
 * method (with normal equations approach) for solving a quadratic optimization problem 
 * with equality constraints. It solves the following problem: \n
 * \n
 *           min c'x + 0.5 * x' G x \n
 *               Ax = b \n
 * \n
 * The implementation was inspired from the algorithm described in the book:\n
 *   Nocedal, Numerical Optimization, 2nd Ed..
 * TEST PASSED
 * 
 * \tparam Matrix A general matrix type, should model the WritableMatrixConcept (and be fully-writable).
 * \tparam Vector1 A vector type, should model the WritableVectorConcept.
 * \tparam Vector2 A vector type, should model the WritableVectorConcept.
 * \param A The constraint matrix of dimension M*N.
 * \param b The b vector of dimension M.
 * \param G The G matrix of dimension NxN (defines the quadratic function to minimize, should be positive semi-definite).
 * \param c The cost vector of dimension N.
 * \param x Stores, as output, the optimal vector.
 * \param tol The relative tolerance on the singularity of components of the matrices involved.
 * 
 * \author Mikael Persson
 */
template <typename Matrix1, typename Vector1, typename Matrix2, typename Vector2>
void projected_CG_method(const Matrix1& A, const Vector1& b, 
                        const Matrix2&  G, const Vector2& c, Vector2& x, unsigned int max_iter = 20,
                        typename vect_traits<Vector1>::value_type abs_tol = std::numeric_limits<typename vect_traits<Vector1>::value_type>::epsilon(),
                        Vector1* lambda = NULL) {
  typedef typename vect_traits<Vector1>::value_type ValueType;
  typedef typename vect_traits<Vector1>::size_type SizeType;
  using std::swap;
  using std::fabs;
  
  SizeType N = c.size();
  SizeType M = b.size();
  
  mat<ValueType,mat_structure::rectangular> A_tmp(transpose_view(A));
  mat<ValueType,mat_structure::rectangular> R(N,M);
  mat<ValueType,mat_structure::square> Q(N);
  decompose_QR(A_tmp,Q,R,abs_tol);
  mat<ValueType,mat_structure::rectangular> L = mat<ValueType,mat_structure::rectangular>(transpose_view(R));
  L.set_col_count(M,true);
  
  Vector1 b_tmp = b;
  mat_vect_adaptor<Vector1> b_tmp_mat(b_tmp);
  ReaK::detail::backsub_Cholesky_impl(L,b_tmp_mat);
  x = b_tmp * A;
  
  Vector2 r = G * x + c;
  Vector1 Ar = A * r;
  mat_vect_adaptor<Vector1> Ar_mat(Ar);
  ReaK::detail::backsub_Cholesky_impl(L,Ar_mat);
  Vector2 g = r; g -= Ar * A;
  Vector2 d = -g;
  Vector2 Gd = G * d;
  ValueType rg = r * g;
  
  unsigned int k = 0;
  while(fabs(rg) > abs_tol) {
    ValueType alpha = rg / (d * Gd);
    x += alpha * d;
    r += alpha * Gd;
    Ar = A * r;
    ReaK::detail::backsub_Cholesky_impl(L,Ar_mat);
    g = r; g -= Ar * A;
    ValueType rg_p = r * g;
    ValueType beta = rg_p / rg;
    rg = rg_p;
    d *= beta; d -= g;
    Gd = G * d;
    if(++k > max_iter)
      throw maximum_iteration(max_iter);
  };
  
  if(lambda) {
    (*lambda) = (c + G *x) * mat_const_sub_block< mat<ValueType,mat_structure::square> >(Q, N, M, 0, 0);
    mat_vect_adaptor< Vector1 > lambda_mat(*lambda);
    ReaK::detail::backsub_R_impl(transpose_view(L),lambda_mat,abs_tol);
  };
  
};






};

};


#endif








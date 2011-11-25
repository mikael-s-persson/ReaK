/**
 * \file simplex_method.hpp
 *
 * The following library is an implementation of the Simplex Method to solve a linear programming 
 * problem. The algorithm follows that of Chvatal and Vasek 1983.
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

#ifndef REAK_SIMPLEX_METHOD_HPP
#define REAK_SIMPLEX_METHOD_HPP

#include "base/defs.hpp"

#include "optim_exceptions.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_qr_decomp.hpp"

#include <vector>

namespace ReaK {
  
  
namespace optim {
  
  
namespace detail {
  
  
  template <typename Matrix, typename Vector>
  void simplex_method_loop_impl(const Matrix& A, Matrix& B, const Vector& c, Vector& c_B, 
                                Vector& x, const Vector& l, const Vector& u,
				std::vector< typename vect_traits<Vector>::size_type >& i_b,
				std::vector< typename vect_traits<Vector>::size_type >& i_n) {
    typedef typename vect_traits<Vector>::value_type ValueType;
    typedef typename vect_traits<Vector>::size_type SizeType;
    using std::swap;
  
    SizeType N = i_n.size();
    SizeType M = A.get_row_count();
    
    Matrix B_Q(M,M);
    Matrix B_R(M,M);
    decompose_QR(B,B_Q,B_R);
  
    while(true) {
      // Step 1
      Vector y = c_B * B_Q;
      backsub_R(B_R,y);
    
      // Step 2
      ValueType sum;
      SizeType enter_var = N;
      for(SizeType i=0; i < N; ++i) {
	sum = y * slice(A)(range(0,M-1),i_n[i]);
        if (((sum < c[i_n[i]]) && (x[i_n[i]] < u[i_n[i]])) 
	    || ((sum > c[i_n[i]]) && (x[i_n[i]] > l[i_n[i]]))) {
          enter_var = i;
          break;
        };
      };
      if (enter_var == N)
        return;
      // Step 3
      y = slice(A)(range(0,M-1),i_n[enter_var]);
      y = y * B_Q;
      backsub_R(B_R,y);
      // Step 4
      if (sum < c[i_n[enter_var]]) {
        ValueType t_max = std::numeric_limits<ValueType>::infinity();
        SizeType leave_var = 0;
        for(SizeType i=0; i < M; ++i) {
          if(y[i] > 0.0) {
            if (t_max * y[i] > x[i_b[i]] - l[i_b[i]]) {
              t_max = (x[i_b[i]] - l[i_b[i]]) / y[i];
              leave_var = i;
            };
          } else if (y[i] < 0.0) {
            if (t_max * y[i] < x[i_b[i]] - u[i_b[i]]) {
              t_max = (x[i_b[i]] - u[i_b[i]]) / y[i];
              leave_var = i;
            };
          };
        };
        if (t_max > u[i_n[enter_var]] - x[i_n[enter_var]]) {
          t_max = u[i_n[enter_var]] - x[i_n[enter_var]];
          x[i_n[enter_var]] = u[i_n[enter_var]];
          for (SizeType i=0; i < M; ++i)
            x[i_b[i]] -= t_max * y[i];
          continue;
        } else if (t_max != std::numeric_limits<ValueType>::infinity()) {
          x[i_n[enter_var]] += t_max;
          for (SizeType i=0; i < M; ++i)
            x[i_b[i]] -= t_max * y[i];
	  slice(B)(range(0,M-1),leave_var) = slice(A)(range(0,M-1),i_n[enter_var]);
	  decompose_QR(B,B_Q,B_R);
          c_B[leave_var] = c[i_n[enter_var]];
	  swap(i_b[leave_var],i_n[enter_var]);
        } else
          throw unbounded_problem("Simplex method failed due to an unbounded search domain!");
      } else if (sum > c[i_n[enter_var]]) {
        ValueType t_max = std::numeric_limits<ValueType>::infinity();
        SizeType leave_var = 0;
        for(SizeType i=0; i < M; ++i) {
          if(y[i] > 0.0) {
            if (t_max * y[i] > u[i_b[i]] - x[i_b[i]]) {
              t_max = (u[i_b[i]] - x[i_b[i]]) / y[i];
              leave_var = i;
            };
          } else if (y[i] < 0.0) {
            if (t_max * y[i] < l[i_b[i]] - x[i_b[i]]) {
              t_max = (l[i_b[i]] - x[i_b[i]]) / y[i];
              leave_var = i;
            };
          };
        };
        if (t_max > x[i_n[enter_var]] - l[i_n[enter_var]]) {
          t_max = x[i_n[enter_var]] - l[i_n[enter_var]];
          x[i_n[enter_var]] = l[i_n[enter_var]];
          for (SizeType i=0; i < M; ++i)
            x[i_b[i]] += t_max * y[i];
          continue;
        } else if (t_max != std::numeric_limits<ValueType>::infinity()) {
          x[i_n[enter_var]] -= t_max;
          for (SizeType i=0; i < M; ++i)
            x[i_b[i]] += t_max * y[i];
	  slice(B)(range(0,M-1),leave_var) = slice(A)(range(0,M-1),i_n[enter_var]);
	  decompose_QR(B,B_Q,B_R);
          c_B[leave_var] = c[i_n[enter_var]];
	  swap(i_b[leave_var],i_n[enter_var]);
        } else
          throw unbounded_problem("Simplex method failed due to an unbounded search domain!");
      };

    };
  };
  
  
};


/**
 * This function is an implementation of the general two-phase revised simplex
 * method for bounded variables. It solves the following problem: \n
 * \n
 *           max c'x \n
 *               Ax = b \n
 *             l <= x <= u \n
 * \n
 * The implementation was inspired from the algorithm described in the book:\n
 *   Chvatal, Vasek, Linear Programming, W. H. Freeman and Company, 1983.
 * 
 * \tparam Matrix A general matrix type, should model the WritableMatrixConcept (and be fully-writable).
 * \tparam Vector A vector type, should model the WritableVectorConcept.
 * \param A The constraint matrix of dimension M*N.
 * \param b The b vector of dimension M.
 * \param c The cost vector of dimension N.
 * \param x0 The initial guess for the optimal vector, then stores, as output, the optimal vector.
 * \param l The lower-bound on the independent variables, if there is no lower bound for a variable, set to minus infinity.
 * \param u The upper-bound on the independent variables, if there is no upper bound for a variable, set to infinity.
 *
 * \throw 
 * 
 * \author Mikael Persson
 */
template <typename Matrix, typename Vector>
int simplex_method(const Matrix& A, const Vector& b, const Vector& c, 
                   Vector& x0, const Vector& l, const Vector& u,
		   typename vect_traits<Vector>::value_type tol = std::numeric_limits<typename vect_traits<Vector>::value_type>::epsilon()) {
  typedef typename vect_traits<Vector>::value_type ValueType;
  typedef typename vect_traits<Vector>::size_type SizeType;
  using std::swap;
  using std::fabs;
  
  SizeType N = c.size();
  SizeType M = b.size();
  
  Vector x(N + M);
  std::copy(x0.begin(), x0.end(), x.begin());
  
  Vector b_G = b;
  
  Vector c_G(N+M, 0.0);
  Vector c_B(M);
  
  Vector l_G(N+M);
  Vector u_G(N+M);
  std::copy(l.begin(), l.end(), l_G.begin());
  std::copy(u.begin(), u.end(), u_G.begin());
  
  Matrix A_G(M,N+M);
  sub(A_G)(range(0,M-1),range(0,N-1)) = A;
  sub(A_G)(range(0,M-1),range(N,N+M-1)) = mat<ValueType,mat_structure::identity>(M);
  
  Matrix B = mat<ValueType,mat_structure::identity>(M);
  Matrix B_Q = mat<ValueType,mat_structure::identity>(M);
  Matrix B_R = mat<ValueType,mat_structure::identity>(M);
  
  std::vector< SizeType > i_b(M);
  for(SizeType i=0; i < M; ++i)
    i_b[i] = i + N;
  std::vector< SizeType > i_n(N);
  for(SizeType i=0; i < N; ++i)
    i_n[i] = i;
  
  for(SizeType i=0; i < N; ++i) {
    if(x[i] < l[i])
      x[i] = l[i];
    if(x[i] > u[i])
      x[i] = u[i];
  };

  bool FeasibleStart = true;
  for(SizeType i=0; i < M; ++i) {
    x[N+i] = b_G[i];
    for(SizeType j=0; j < N; ++j)
      x[N+i] -= A(i,j) * x[j];
    FeasibleStart &&= (fabs(x[N+i]) < tol);
    if (x[N+i] >= 0.0) {
      l_G[N+i] = 0.0;
      u_G[N+i] = std::numeric_limits<ValueType>::infinity();
    } else {
      l_G[N+i] = -std::numeric_limits<ValueType>::infinity();
      u_G[N+i] = 0.0;
    };
    c_G[N+i] = -1.0;
    c_B[i] = -1.0;
  };
  
  if (!FeasibleStart) {
    // First-Phase
    detail::simplex_method_loop_impl(A_G,B,c_G,c_B,x,l_G,u_G,i_b,i_n);
    
    // Did the first phase succeed?
    for (SizeType i=0; i < M; ++i)
      if (fabs(x[N+i]) > tol)
        throw infeasible_problem("Simplex method failed due to an empty search domain! No feasible solution exists!");
  };
  
  // Getting Rid of the Artificial Variables
  decompose_QR(B,B_Q,B_R);
  for(SizeType i=0; i < M; ++i) {
    if(i_b[i] >= N) {
      Vector y = slice(A_G)(range(0,M-1),i_b[i]);
      y = y * B_Q;
      backsub_R(B_R,y);
      for(SizeType j=0; j < N; ++j) {
        T sum = y * slice(A_G)(range(0,M-1),i_n[j]);
        if ((sum != 0.0) && (i_n[j] < N)) {
	  slice(B)(range(0,M-1),i) = slice(A_G)(range(0,M-1),i_n[j]);
          swap(i_b[i],i_n[j]);
	  break;
        };
      };
      decompose_QR(B,B_Q,B_R);
    };
  };
  // Getting Rid of the Redundant equations
  SizeType RedundantCount = 0;
  for(SizeType i=0; i < M; ++i) {
    if(i_b[i] >= N) {
      // Must Delete Redundant Equation
      ++RedundantCount;
      --M;
      Matrix tempA(M,N);
      Matrix tempB(M,M);
      Vector tempb(M);
      
      tempB = sub(B)(range(0,i_b[i]-N-1),range(0,i-1)) & sub(B)(range(0,i_b[i]-N-1),range(i+1,M)) |
              sub(B)(range(i_b[i]-N+1,M),range(0,i-1)) & sub(B)(range(i_b[i]-N+1,M),range(i+1,M)) ;
      
      tempA = sub(A_G)(range(0,i_b[i]-N-1),range(0,N-1)) |
              sub(A_G)(range(i_b[i]-N+1,M),range(0,N-1));
      
      std::copy(b_G.begin(),b_G.begin() + i_b[i] - N, tempb.begin());
      std::copy(b_G.begin() + i_b[i] - N + 1, b_G.end(), tempb.begin() + i_b[i] - N);
      swap(tempb,b_G);
      swap(tempA,A_G);
      swap(tempB,B);
    };
  };
  decompose_QR(B,B_Q,B_R);

  // Prepare the variables for the second phase
  {
    x.resize(N);

    c_B.resize(M);
    for(SizeType i=0; i < M; ++i)
      c_B[i] = c[i_b[i]];

    std::vector<SizeType> tempIB(M);
    std::vector<SizeType> tempIN(N-M);
    {
      SizeType j = 0;
      for(SizeType i=0; i < M + RedundantCount; ++i) {
        if(i_b[i] < N) {
          tempIB[j] = i_b[i];
          ++j;
        };
      };
      j = 0;
      for(SizeType i=0; i < N; ++i) {
        if(i_n[i] < N) {
          tempIN[j] = i_n[i];
          j++;
        };
      };
    };
    i_b.swap(tempIB);
    i_n.swap(tempIN);
  };
  
  detail::simplex_method_loop_impl(A_G,B,c,c_B,x,l,u,i_b,i_n);
  
  x0 = x;
  
};







};

};


#endif








/**
 * \file finite_diff_jacobians.hpp
 *
 * The following library provides implementations of finite difference methods for evaluating a Jacobian
 * for a function of multiple variables and multiple outputs (or single).
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


namespace ReaK {
  
  
namespace optim {
  
  
namespace detail {
  
  

template <typename Function, typename Vector1, typename Vector2, typename Matrix>
void compute_jacobian_2pts_forward_impl(Function f, Vector1& x, const Vector2& y, Matrix& jac, typename vect_traits<Vector1>::value_type delta) {
  typedef typename vect_traits<Vector1>::value_type ValueType;
  typedef typename vect_traits<Vector1>::size_type SizeType;
  using std::fabs;
  
  SizeType N = x.size();
  SizeType M = y.size();
  if(jac.get_row_count() != M) 
    jac.set_row_count(M);
  if(jac.get_col_count() != N)
    jac.set_col_count(N);
  
  Vector2 y1 = y;

  for(SizeType j=0; j < N; ++j) {
    /* determine d=max(1E-04*|p[j]|, delta), see HZ */
    ValueType d = ValueType(1E-04) * fabs(x[j]);
    if(d < delta)
      d = delta;

    ValueType tmp = x[j];
    x[j] += d;

    y1 = f(x);

    x[j] = tmp; /* restore */

    slice(jac)(range(0,M-1),j) = (y1 - y) * (1.0 / d);
    
  };
};

template <typename Function, typename Vector1, typename Vector2, typename Matrix>
void compute_jacobian_2pts_central_impl(Function f, Vector1& x, const Vector2& y, Matrix& jac, typename vect_traits<Vector1>::value_type delta) {
  typedef typename vect_traits<Vector1>::value_type ValueType;
  typedef typename vect_traits<Vector1>::size_type SizeType;
  using std::fabs;
  
  SizeType N = x.size();
  SizeType M = y.size();
  if(jac.get_row_count() != M) 
    jac.set_row_count(M);
  if(jac.get_col_count() != N)
    jac.set_col_count(N);

  Vector2 y_prev = y;
  Vector2 y_next = y;

  for(SizeType j=0; j < N; ++j) {
    /* determine d=max(1E-04*|p[j]|, delta), see HZ */
    ValueType d= ValueType(1E-04) * fabs(x[j]);
    if( d < delta )
      d = delta;

    ValueType tmp = x[j];
    x[j] -= d;
    y_prev = f(x);
    
    x[j] = tmp + d;
    y_next = f(x);
    x[j]=tmp; /* restore */
    
    slice(jac)(range(0,M-1),j) = (y_next - y_prev) * (0.5 / d);
  };
};

template <typename Function, typename Vector1, typename Vector2, typename Matrix>
void compute_jacobian_5pts_central_impl(Function f, Vector1& x, const Vector2& y, Matrix& jac, typename vect_traits<Vector1>::value_type delta) {
  typedef typename vect_traits<Vector1>::value_type ValueType;
  typedef typename vect_traits<Vector1>::size_type SizeType;
  using std::fabs;
  
  SizeType N = x.size();
  SizeType M = y.size();
  if(jac.get_row_count() != M) 
    jac.set_row_count(M);
  if(jac.get_col_count() != N)
    jac.set_col_count(N);
  
  Vector2 y0 = y;
  Vector2 y1 = y;
  Vector2 y2 = y;
  Vector2 y3 = y;
  
  for(SizeType i=0; i < N; ++i) {
    /* determine d=max(1E-04*|p[j]|, delta), see HZ */
    ValueType d = ValueType(1E-04) * fabs(x[i]);
    if( d < delta )
      d = delta;
    
    ValueType tmp = x[i];
    x[i] -= 2.0 * d;
    y0 = f(x);
    x[i] += d;
    y1 = f(x);
    x[i] += 2.0 * d;
    y2 = f(x);
    x[i] += d;
    y3 = f(x);
    
    x[i] = tmp; // restore
    
    slice(jac)(range(0,M-1),i) = (y0 - 8.0 * (y1 + y2) - y3) * (1.0 / (12.0 * d));
  };
};


};








};

};

#endif









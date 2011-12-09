/**
 * \file broyden_method.hpp
 * 
 * This library provides a root-finder function that uses the Broyden's method (good or fast) 
 * for solving a system of simultaneous non-linear equations (i.e. finding the root of a N-dimensional
 * function).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
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

#ifndef REAK_BROYDEN_METHOD_HPP
#define REAK_BROYDEN_METHOD_HPP

#include <limits>
#include <cmath>

#include "lin_alg/mat_num_exceptions.hpp"
#include "lin_alg/vect_concepts.hpp"
#include "lin_alg/mat_alg_square.hpp"

namespace ReaK {


/*************************************************************************
                        Broyden's Root Finding Methods
*************************************************************************/


/**
 * This function template performs Broyden's good method for the root of a function. This assumes 
 * that the function is monotonic and has a unique root near the given bounds.
 * \tparam Vector A vector type of the independent and dependent value of the function.
 * \tparam RootedFunction A unary functor type.
 * \param x_prev The x(-1) of the search for the root.
 * \param x0 The x(0) of the search for the root, also stores as output, the root of the function.
 * \param f The functor of which the root is sought.
 * \param tol The tolerance on the absolute value of the root.
 * \param max_iter The maximum number of iterations until abandonment.
 * \throw maximum_iteration If the maximum number of iterations is reached before convergence.
 * \throw singularity_error If a stationary point is reached.
 */
template <typename Vector,
	  typename RootedFunction>
void broyden_good_method(const Vector& x_prev, Vector& x0, RootedFunction f, const T& tol = std::numeric_limits<T>::epsilon(), std::size_t max_iter = 50) 
{
  using std::fabs;
  typedef typename vect_traits<Vector>::value_type ValueType;
  typedef typename vect_traits<Vector>::size_type SizeType;
  
  Vector dx = x0 - x_prev;
  Vector y0 = f(x0);
  Vector dy = y0 - f(x_prev);
  std::size_t iter = 0;
  
  mat<ValueType, mat_structure::square> J_inv = mat<ValueType, mat_structure::identity>(dx.size());
  Vector Jdy = dy;
  ValueType denom = dx * dy;
  if( fabs(denom) < tol )
    throw singularity_error("Broyden's good method failed due to a stationary point!");
  Vector dxJ = dx;
  for(SizeType i = 0; i < dx.size(); ++i)
    for(SizeType j = 0; j < dx.size(); ++j)
      J_inv(i,j) += (dx[i] - dy[i]) * dx[j] / denom;
  
  while(true) {
    
    dx = -J_inv * y0;
    x0 += dx;
    
    if(norm_2(dx) < tol)
      return;
  
    if( ++iter > max_iter ) 
      throw maximum_iteration("Broyden's good method diverged, as detected by reaching the maximum iteration limit!");
  
    dy = f(x0) - y0;
    y0 += dy;
    
    Jdy = J_inv * dy;
    denom = dx * Jdy;
    if( fabs(denom) < tol )
      throw singularity_error("Broyden's good method failed due to a stationary point!");
    dxJ = dx * J_inv;
    
    for(SizeType i = 0; i < dx.size(); ++i)
      for(SizeType j = 0; j < dx.size(); ++j)
        J_inv(i,j) += (dx[i] - Jdy[i]) * dxJ[j] / denom;
    
  };
  
  return x0;
};



/**
 * This function template performs Broyden's fast method for the root of a function. This assumes 
 * that the function is monotonic and has a unique root near the given bounds.
 * \tparam Vector A vector type of the independent and dependent value of the function.
 * \tparam RootedFunction A unary functor type.
 * \param x_prev The x(-1) of the search for the root.
 * \param x0 The x(0) of the search for the root, also stores as output, the root of the function.
 * \param f The functor of which the root is sought.
 * \param tol The tolerance on the absolute value of the root.
 * \param max_iter The maximum number of iterations until abandonment.
 * \throw maximum_iteration If the maximum number of iterations is reached before convergence.
 * \throw singularity_error If a stationary point is reached.
 */
template <typename Vector,
	  typename RootedFunction>
void broyden_fast_method(const Vector& x_prev, Vector& x0, RootedFunction f, const T& tol = std::numeric_limits<T>::epsilon(), std::size_t max_iter = 50) 
{
  using std::fabs;
  typedef typename vect_traits<Vector>::value_type ValueType;
  typedef typename vect_traits<Vector>::size_type SizeType;
  
  Vector dx = x0 - x_prev;
  Vector y0 = f(x0);
  Vector dy = y0 - f(x_prev);
  std::size_t iter = 0;
  
  mat<ValueType, mat_structure::square> J_inv = mat<ValueType, mat_structure::identity>(dx.size());
  Vector Jdy = dy;
  ValueType denom = dy * dy;
  if( fabs(denom) < tol )
    throw singularity_error("Broyden's fast method failed due to a stationary point!");
  for(SizeType i = 0; i < dx.size(); ++i)
    for(SizeType j = 0; j < dx.size(); ++j)
      J_inv(i,j) += (dx[i] - dy[i]) * dy[j] / denom;
  
  while(true) {
    
    dx = -J_inv * y0;
    x0 += dx;
    
    if(norm_2(dx) < tol)
      return;
  
    if( ++iter > max_iter ) 
      throw maximum_iteration("Broyden's fast method diverged, as detected by reaching the maximum iteration limit!");
  
    dy = f(x0) - y0;
    y0 += dy;
    
    denom = dy * dy;
    if( fabs(denom) < tol )
      throw singularity_error("Broyden's fast method failed due to a stationary point!");
    Jdy = J_inv * dy;
    for(SizeType i = 0; i < dx.size(); ++i)
      for(SizeType j = 0; j < dx.size(); ++j)
        J_inv(i,j) += (dx[i] - Jdy[i]) * dy[j] / denom;
    
  };
  
  return x0;
};



};

#endif




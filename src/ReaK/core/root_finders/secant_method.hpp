/**
 * \file secant_method.hpp
 * 
 * This library provides a root-finder function that uses the secant method.
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

#ifndef REAK_SECANT_METHOD_HPP
#define REAK_SECANT_METHOD_HPP

#include <limits>
#include <cmath>

namespace ReaK {


/*************************************************************************
                        Secant Root Finding Method
*************************************************************************/


/**
 * This function template performs a secant method search for the root of a function. This assumes 
 * that the function is monotonic and has a unique root between the two given bounds.
 * \tparam T A scalar value type of the independent and dependent value of the function.
 * \tparam RootedFunction A unary functor type.
 * \param low_bound The lower bound of the search for the root.
 * \param hi_bound The upper bound of the search for the root.
 * \param f The functor of which the root is sought.
 * \param tol The tolerance on the absolute value of the root.
 * \return The independent variable value at which the root was found, or the bound value that was reached while seeking the root (and diverging).
 */
template <typename T,
	  typename RootedFunction>
T secant_method(const T& low_bound, const T& hi_bound, RootedFunction f, const T& tol = std::numeric_limits<T>::epsilon()) 
{
  using std::fabs;
  
  T x0 = low_bound;
  T x1 = hi_bound;
  T x0_value = f(x0);
  T x1_value = f(x1);
  
  if( x0_value * x1_value > 0.0 )
    return (fabs(x0_value) < fabs(x1_value) ? x0 : x1);
  
  while((fabs(x0 - x1) > tol) && (fabs(x0_value - x1_value) > tol)) {
    
    T x2 = x1 - x1_value * (x1 - x0) / (x1_value - x0_value);
    
    if(x2 > hi_bound)
      return hi_bound;
    if(x2 < low_bound)
      return low_bound;
    
    x0 = x1; x0_value = x1_value;
    x1 = x2; x1_value = f(x1);
    
  };
  
  return x1;
};



};

#endif




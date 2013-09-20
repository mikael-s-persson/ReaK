/**
 * \file special_secant_methods.hpp
 * 
 * This library provides root-finder functions that use the secant method (Ford-3) within
 * special algorithms to deal with multiple roots and predicates.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date September 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_SPECIAL_SECANT_METHODS_HPP
#define REAK_SPECIAL_SECANT_METHODS_HPP

#include <limits>
#include <cmath>
#include <algorithm>

namespace ReaK {


/**
 * This function template performs the Ford-3 algorithm for the root of a function. This assumes 
 * that the function is monotonic and has a unique root between the two given bounds. The Ford-3
 * method is a Illinois-type method with a refined choice of scaling factor that acheives super-linear
 * convergence.
 * \tparam T A scalar value type of the independent and dependent value of the function.
 * \tparam RootedFunction A unary functor type.
 * \param low_bound The lower bound of the search for the root, and store as output the lower-bound of the narrowed down interval in which the root exists.
 * \param hi_bound The upper bound of the search for the root, and store as output the upper-bound of the narrowed down interval in which the root exists.
 * \param f The functor of which the root is sought.
 * \param tol The tolerance on the size of the narrowed-down intervale in which the root exists.
 */
template <typename T,
          typename RootedFunction>
void ford3_method(T& low_bound, T& hi_bound, RootedFunction f, const T& tol = std::numeric_limits<T>::epsilon()) 
{
  using std::fabs;
  
  T x0_value = f(low_bound);
  T x1_value = f(hi_bound);
  
   if( x0_value * x1_value > 0.0 )
     return;
  
  int last_retained_bound = 0;
  T abs_tol = tol * fabs(hi_bound - low_bound);
  T abs_f_tol = tol * (fabs(x0_value) + fabs(x1_value));
  
  while(fabs(low_bound - hi_bound) > abs_tol) {
    
    T x2 = low_bound / (T(1.0) - x0_value / x1_value) + hi_bound / (T(1.0) - x1_value / x0_value);
    T x2_value = f(x2);
    
    if(x2_value * x1_value > 0.0) {
      hi_bound = x2;
      if(fabs(x2_value) < abs_f_tol) {
        low_bound = hi_bound;
        return;
      };
      if(last_retained_bound == -1) {
        x0_value *= 1.0 - x2_value / (x1_value * (1.0 - x2_value / x0_value));
      };
      x1_value = x2_value;
      last_retained_bound = -1;
    } else {
      low_bound = x2;
      if(fabs(x2_value) < abs_f_tol) {
        hi_bound = low_bound;
        return;
      };
      if(last_retained_bound == 1) {
        x1_value *= 1.0 - x2_value / (x0_value * (1.0 - x2_value / x1_value));
      };
      x0_value = x2_value;
      last_retained_bound = 1;
    };
    
  };
  
};


};

#endif




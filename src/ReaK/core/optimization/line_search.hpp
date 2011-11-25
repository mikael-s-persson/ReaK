/**
 * \file line_search.hpp
 *
 * The following library is a collection of line-search optimization algorithms. 
 * For unimodal cost functions (monotonic derivative, i.e. a single global optimum), 
 * three classic algorithms are found in this library: Dichotomous search, Golden-Section
 * search and Fibonacci search.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date September 2010
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

#ifndef REAK_LINE_SEARCH_HPP
#define REAK_LINE_SEARCH_HPP

#include "base/defs.hpp"

#include <vector>

namespace ReaK {
  
  
namespace optim {


namespace detail { 
 
  template <typename T, typename Function>
  T dichotomous_search_impl(Function f, T& low_bound, T& up_bound, T delta, T tol) {
    using std::fabs;
    while (true) {
      T midpoint = (low_bound + up_bound) * T(0.5);
      if(fabs(up_bound - low_bound) < tol)
        return f(midpoint);
      T f1 = f(midpoint - delta);
      T f2 = f(midpoint + delta);
      if( f1 > f2 ) {
	low_bound = midpoint - delta;
	delta *= T(0.5);
      } else {
	up_bound = midpoint + delta;
	delta *= T(0.5);
      };
    };
  };
  
  const double GoldenRatioPhi = 1.618033988;
  
  template <typename T, typename Function>
  T golden_section_search_impl(Function f, T& low_bound, T& up_bound, T tol) {
    using std::fabs;
    
    T mid_value = low_bound + (up_bound - low_bound) / GoldenRatioPhi;
    T mid_cost = f(mid_value);
    
    while (true) {
      if(fabs(low_bound - up_bound) < tol)
        return f((low_bound + up_bound) * T(0.5));
  
      T test_value = mid_value + (up_bound - mid_value) / GoldenRatioPhi;
      T test_cost = f(test_value);
      
      if(test_cost < mid_cost) {
	low_bound = mid_value; 
	mid_value = test_value; mid_cost = test_cost;
      } else {
	up_bound = low_bound;
	low_bound = test_value;
      };
    };
  };
  
  
  std::vector<int>::reverse_iterator get_fibonacci_iter_for_tolerance(double low_bound, double up_bound, double tol) {
    using std::fabs;
    static int first_fib[2] = {0,1};
    static std::vector<int> fib_seq(first_fib,first_fib + 2);
    while(fabs(up_bound - low_bound) > fib_seq.back() * tol)
      fib_seq.push_back(fib_seq.back() + (*(fib_seq.rbegin()+1)));
    return fib_seq.rbegin();
  };
  
  
  template <typename T, typename Function>
  T fibonacci_search_impl(Function f, T& low_bound, T& up_bound, T tol) {
    std::vector<int>::reverse_iterator fib_iter = get_fibonacci_iter_for_tolerance(low_bound, up_bound, tol);
    
    T Lstar2 = (up_bound - low_bound) * T(*(fib_iter+2)) / T(*fib_iter);
    T x1 = low_bound + Lstar2; T x1_value = f(x1);
    T x2 = up_bound  - Lstar2; T x2_value = f(x2);
    while(true) {
      ++fib_iter;
      if((*fib_iter) == 2)
        return x1_value;
      if(x1_value > x2_value) {
        low_bound = x1;
	x1 = x2; x1_value = x2_value;
	x2 = up_bound - T(*(fib_iter+2)) / T(*fib_iter) * (up_bound - low_bound); 
	x2_value = f(x2);
      } else {
        up_bound = x2;
	x2 = x1; x2_value = x1_value;
	x1 = low_bound + T(*(fib_iter+2)) / T(*fib_iter) * (up_bound - low_bound); 
	x1_value = f(x1);
      };
    };
  };
  
};
  
/**
 * This function performs a 1D optimum Dichotomous search on a unimodal cost function to find the value for 
 * the lowest cost (both value and cost are returned by this function as a pair). The search guarantees
 * an uncertainty interval of length tol and limits the search to within (low_bound .. up_bound).
 * \tparam T The value type for the scalar that represents both the independent variable and the cost values.
 * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
 * \param f The cost functor.
 * \param low_bound The lower bound for the search.
 * \param up_bound the upper bound for the search.
 * \param tol the uncertainty on the X value at which the optimization should stop.
 * \return the cost that is the optimum point found, within uncertainty interval (low_bound, up_bound).
 */
template <typename T, typename Function>
T dichotomous_search(Function f, T& low_bound, T& up_bound, T tol) {
  return detail::dichotomous_search_impl(f,low_bound,up_bound,T((up_bound - low_bound) * 0.1),tol);
};

/**
 * This function performs a 1D optimum Golden-Section search on a unimodal cost function to find the value for 
 * the lowest cost (both value and cost are returned by this function as a pair). The search guarantees
 * an uncertainty interval of length tol and limits the search to within (bound1 .. bound2) (the bounds 
 * do not need to be sorted).
 * \tparam T The value type for the scalar that represents both the independent variable and the cost values.
 * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
 * \param f the cost functor.
 * \param bound1 the first bound for the search.
 * \param bound2 the second bound for the search.
 * \param tol the uncertainty on the X value at which the optimization should stop.
 * \return the cost that is the optimum point found, within uncertainty interval (bound1, bound2).
 */
template <typename T, typename Function>
T golden_section_search(Function f, T& bound1, T& bound2, T tol) {
  return detail::golden_section_search_impl(f, bound1, bound2, tol);
};

/**
 * This function performs a 1D optimum Fibonacci search on a unimodal cost function aFunc to find the value for 
 * the lowest cost (both value and cost are returned by this function as a pair). The search guarantees
 * an uncertainty interval of length aToleranceX and limits the search to within (aLowBound .. aUpBound).
 * \tparam T The value type for the scalar that represents both the independent variable and the cost values.
 * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
 * \param f the cost functor.
 * \param low_bound the lower bound for the search.
 * \param up_bound the upper bound for the search.
 * \param tol the uncertainty on the X value at which the optimization should stop.
 * \return the cost that is the optimum point found, within uncertainty interval (low_bound, up_bound).
 */
template <typename T, typename Function>
T fibonacci_search(Function f, T& low_bound, T& up_bound, T tol) {
  using std::swap;
  if(low_bound > up_bound)
    swap(low_bound,up_bound);
  
  return detail::fibonacci_search_impl(f,low_bound,up_bound,tol);
};


};

};


#endif

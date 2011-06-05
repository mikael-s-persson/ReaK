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

#ifndef LINE_SEARCH_HPP
#define LINE_SEARCH_HPP

#include "function_types.hpp"

namespace ReaK {
  
  
namespace optim {


/**
 * This function performs a 1D optimum Dichotomous search on a unimodal cost function aFunc to find the value for 
 * the lowest cost (both value and cost are returned by this function as a pair). The search guarantees
 * an uncertainty interval of length aToleranceX and limits the search to within (aLowBound .. aUpBound).
 * \param aFunc the cost function-object (an object which implements the interface cost_function_1D).
 * \param aLowBound the lower bound for the search (must be lower than aUpBound).
 * \param aUpBound the upper bound for the search (must be greater than aLowBound).
 * \param aToleranceX the uncertainty on the X value at which the optimization should stop.
 * \return a value and cost pair that is the optimum point found, within uncertainty interval of length aToleranceX.
 */
template <class T>
value_cost_pair<T> DichotomousSearch(const cost_function_1D<T>& aFunc, T aLowBound, T aUpBound, T aToleranceX);

/**
 * This function performs a 1D optimum Golden-Section search on a unimodal cost function aFunc to find the value for 
 * the lowest cost (both value and cost are returned by this function as a pair). The search guarantees
 * an uncertainty interval of length aToleranceX and limits the search to within (aBound1 .. aBound2).
 * \param aFunc the cost function-object (an object which implements the interface cost_function_1D).
 * \param aBound1 the first bound for the search.
 * \param aBound2 the second bound for the search.
 * \param aToleranceX the uncertainty on the X value at which the optimization should stop.
 * \return a value and cost pair that is the optimum point found, within uncertainty interval of length aToleranceX.
 */
template <class T>
value_cost_pair<T> GoldenSectionSearch(const cost_function_1D<T>& aFunc, T aBound1, T aBound2, T aToleranceX);

/**
 * This function performs a 1D optimum Fibonacci search on a unimodal cost function aFunc to find the value for 
 * the lowest cost (both value and cost are returned by this function as a pair). The search guarantees
 * an uncertainty interval of length aToleranceX and limits the search to within (aLowBound .. aUpBound).
 * \param aFunc the cost function-object (an object which implements the interface cost_function_1D).
 * \param aLowBound the lower bound for the search (must be lower than aUpBound).
 * \param aUpBound the upper bound for the search (must be greater than aLowBound).
 * \param aToleranceX the uncertainty on the X value at which the optimization should stop.
 * \return a value and cost pair that is the optimum point found, within uncertainty interval of length aToleranceX.
 */
template <class T>
value_cost_pair<T> FibonacciSearch(const cost_function_1D<T>& aFunc, T aLowBound, T aUpBound, T aToleranceX);


};

};


#endif


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


#include "line_search.hpp"

#include "math/complex_math.hpp"
#include <vector>
#include <iostream>

namespace ReaK {

namespace optim {

  
  
template <class T>
value_cost_pair<T> DichotomousSearchImpl(const cost_function_1D<T>& aFunc, T aLowBound, T aUpBound, T aToleranceX, T aDelta) {
  using std::fabs;
  T midpoint = (aLowBound + aUpBound) * T(0.5);
  if(fabs(aUpBound - aLowBound) < aToleranceX)
    return value_cost_pair<T>(midpoint,aFunc.computeCost(midpoint));
  T f1 = aFunc.computeCost(midpoint - aDelta);
  T f2 = aFunc.computeCost(midpoint + aDelta);
  if( f1 > f2 )
    return DichotomousSearchImpl(aFunc,midpoint - aDelta,aUpBound,aToleranceX,aDelta * T(0.5));
  else
    return DichotomousSearchImpl(aFunc,aLowBound,midpoint + aDelta,aToleranceX,aDelta * T(0.5));
};
  
template <class T>
value_cost_pair<T> DichotomousSearch(const cost_function_1D<T>& aFunc, T aLowBound, T aUpBound, T aToleranceX) {
  return DichotomousSearchImpl(aFunc,aLowBound,aUpBound,aToleranceX,T((aUpBound - aLowBound) * 0.1));
};
  
template value_cost_pair<float> DichotomousSearch(const cost_function_1D<float>& aFunc, float aLowBound, float aUpBound, float aToleranceX);
template value_cost_pair<double> DichotomousSearch(const cost_function_1D<double>& aFunc, double aLowBound, double aUpBound, double aToleranceX);
  


const double GoldenRatioPhi = 1.618033988;
  
template <class T>
value_cost_pair<T> GoldenSectionSearchImpl(const cost_function_1D<T>& aFunc, const value_cost_pair<T>& A, const value_cost_pair<T>& C, const value_cost_pair<T>& B, T aToleranceX) {
  using std::fabs;
  if(fabs(A.value - B.value) < fabs(aToleranceX))
    return value_cost_pair<T>((A.value + B.value) / 2, aFunc.computeCost((A.value + B.value) / 2));
  
  value_cost_pair<T> D(C.value + (B.value - C.value) / GoldenRatioPhi, T(0));
  D.cost = aFunc.computeCost(D.value);
  
  if(D.cost < C.cost)
    return GoldenSectionSearchImpl(aFunc,C,D,B,aToleranceX);
  else
    return GoldenSectionSearchImpl(aFunc,D,C,A,aToleranceX);
};
  
template <class T>
value_cost_pair<T> GoldenSectionSearch(const cost_function_1D<T>& aFunc, T aBound1, T aBound2, T aToleranceX) {
  value_cost_pair<T> A(aBound1, T(0));
  value_cost_pair<T> B(aBound2, T(0));
  value_cost_pair<T> C(A.value + (B.value - A.value) / GoldenRatioPhi, T(0));
  C.cost = aFunc.computeCost(C.value);
  return GoldenSectionSearchImpl(aFunc, A, C, B, aToleranceX);
};

template value_cost_pair<float> GoldenSectionSearch(const cost_function_1D<float>& aFunc, float aBound1, float aBound2, float aToleranceX);
template value_cost_pair<double> GoldenSectionSearch(const cost_function_1D<double>& aFunc, double aBound1, double aBound2, double aToleranceX);



template <class T>
value_cost_pair<T> FibonacciSearchImpl(const cost_function_1D<T>& aFunc, T aLowBound, T aUpBound, 
                                       value_cost_pair<T> X1, value_cost_pair<T> X2, 
                                       std::vector<int>::reverse_iterator aFibIter) {
  if((*aFibIter) == 2)
    return X1; //termination condition when fib(n-j+2) == 2 (meaning X1 and X2 have converged to the same value).
  if(X1.cost > X2.cost) {
    T tmp = aUpBound - T(*(aFibIter+2)) / T(*aFibIter) * (aUpBound - X1.value);
    return FibonacciSearchImpl(aFunc,X1.value,aUpBound,X2,value_cost_pair<T>(tmp,aFunc.computeCost(tmp)),aFibIter+1);
  } else {
    T tmp = aLowBound + T(*(aFibIter+2)) / T(*aFibIter) * (X2.value - aLowBound);
    return FibonacciSearchImpl(aFunc,aLowBound,X2.value,value_cost_pair<T>(tmp,aFunc.computeCost(tmp)),X1,aFibIter+1);
  };
};

template <class T>
value_cost_pair<T> FibonacciSearch(const cost_function_1D<T>& aFunc, T aLowBound, T aUpBound, T aToleranceX) {
  using std::fabs;
  using std::swap;
  if(aLowBound > aUpBound)
    swap(aLowBound,aUpBound);
  
  std::vector<int> fib_seq;
  fib_seq.push_back(0);
  fib_seq.push_back(1);
  while(fabs(aUpBound - aLowBound) > fib_seq.back() * aToleranceX)
    fib_seq.push_back(fib_seq.back() + (*(fib_seq.rbegin()+1)));
  
  std::vector<int>::reverse_iterator FibIter = fib_seq.rbegin();
  T Lstar2 = (aUpBound - aLowBound) * (*(FibIter+2)) / (*FibIter);
  value_cost_pair<T> X1(aLowBound + Lstar2,aFunc.computeCost(aLowBound + Lstar2));
  value_cost_pair<T> X2(aUpBound - Lstar2,aFunc.computeCost(aUpBound - Lstar2));
  
  return FibonacciSearchImpl(aFunc,aLowBound,aUpBound,X1,X2,FibIter+1);
};

template value_cost_pair<float> FibonacciSearch(const cost_function_1D<float>& aFunc, float aBound1, float aBound2, float aToleranceX);
template value_cost_pair<double> FibonacciSearch(const cost_function_1D<double>& aFunc, double aBound1, double aBound2, double aToleranceX);




};

};




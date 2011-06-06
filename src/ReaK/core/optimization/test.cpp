
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

#include <cmath>
#include "line_search.hpp"

#include <iostream>

template <class T>
class MaxTrussSectionStress : public ReaK::optim::cost_function_1D<T> {
  public:
    mutable int evalCount;
    MaxTrussSectionStress() : evalCount(0) { };
    virtual T computeCost(const T& aParameter) const {
      T sigma1 = 0.8165 / aParameter;
      T sigma2 = 1.1154 / (1 - aParameter);
      evalCount++;
      return (sigma1 > sigma2 ? sigma1 : sigma2);
    };
  
};


int main() {

  std::cout << "The actual optimal parameter from analytical calculation is 0.42264." << std::endl;
  MaxTrussSectionStress<double> DichotomousCF;
  std::cout << "The Dichotomous Search has found: " 
            << ReaK::optim::DichotomousSearch(DichotomousCF, 0.30, 0.48, 0.0018).value;
  std::cout << " with " << DichotomousCF.evalCount << " cost function evaluations." << std::endl;
  MaxTrussSectionStress<double> GoldenSectionCF;
  std::cout << "The Golden Section Search has found: " 
            << ReaK::optim::GoldenSectionSearch(GoldenSectionCF, 0.30, 0.48, 0.0018).value;
  std::cout << " with " << GoldenSectionCF.evalCount << " cost function evaluations." << std::endl;
  MaxTrussSectionStress<double> FibonacciCF;
  std::cout << "The Fibonacci Search has found: " 
            << ReaK::optim::FibonacciSearch(FibonacciCF, 0.30, 0.48, 0.0018).value;
  std::cout << " with " << FibonacciCF.evalCount << " cost function evaluations." << std::endl;


};









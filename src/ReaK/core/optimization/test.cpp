
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

static int evalCount;

double max_truss_section_stress(double x) {
  double sigma1 = 0.8165 / x;
  double sigma2 = 1.1154 / (1 - x);
  evalCount++;
  return (sigma1 > sigma2 ? sigma1 : sigma2);  
};


int main() {

  std::cout << "The actual optimal parameter from analytical calculation is 0.42264." << std::endl;
  
  evalCount = 0;
  double l = 0.30; double u = 0.48;
  ReaK::optim::dichotomous_search(max_truss_section_stress, l, u, 0.0018);
  std::cout << "The Dichotomous Search has found: " 
            << ((l + u) * 0.5);
  std::cout << " with " << evalCount << " cost function evaluations." << std::endl;
  
  evalCount = 0;
  l = 0.30; u = 0.48;
  ReaK::optim::golden_section_search(max_truss_section_stress, l, u, 0.0018);
  std::cout << "The Golden Section Search has found: " << ((l + u) * 0.5);
  std::cout << " with " << evalCount << " cost function evaluations." << std::endl;
  
  evalCount = 0;
  l = 0.30; u = 0.48;
  ReaK::optim::fibonacci_search(max_truss_section_stress, l, u, 0.0018);
  std::cout << "The Fibonacci Search has found: " << ((l + u) * 0.5);
  std::cout << " with " << evalCount << " cost function evaluations." << std::endl;

  return 0;
};









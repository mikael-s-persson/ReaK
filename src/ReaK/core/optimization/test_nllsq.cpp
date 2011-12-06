
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

#include "finite_diff_jacobians.hpp"

#include "levenberg_marquardt_method.hpp"
#include "jacobian_transpose_method.hpp"
#include "gauss_newton_method.hpp"

#include <iostream>
#include <cmath>

#ifndef M_PI 
#define M_PI 3.14159265359
#endif

static int evalCount;
static int gradCount;

using namespace ReaK;

vect_n<double> p01_f(const vect_n<double>& x) {
  ++evalCount;
  std::size_t N = x.size();
  std::size_t M = N * 2;
  double sum = 0.0;
  for(std::size_t i = 0; i < N; ++i)
    sum += x[i];
  vect_n<double> result(M,-1.0 - 2.0 * sum / double(M));
  for(std::size_t i = 0; i < N; ++i)
    result[i] += x[i];
  return result;
};

void p01_j(mat<double,mat_structure::rectangular>& J, const vect_n<double>& x, const vect_n<double>& f) {
  ++gradCount;
  std::size_t M = f.size();
  std::size_t N = x.size();
  J.set_col_count(N);
  J.set_row_count(M);
  for(std::size_t i = 0; i < M; ++i)
    for(std::size_t j = 0; j < N; ++j)
      J(i,j) = -2.0 / double(M);
  for(std::size_t i = 0; i < N; ++i)
    J(i,i) += 1.0;
};

vect_n<double> p01_sol = vect_n<double>(10,-1.0);
vect_n<double> p01_start = vect_n<double>(10,1.0);
vect_n<double> p01_lower = vect_n<double>(10,-4.0);
vect_n<double> p01_upper = vect_n<double>(10,4.0);
std::string p01_name = "Linear Function, full-rank";

vect_n<double> p02_f(const vect_n<double>& x) {
  ++evalCount;
  std::size_t N = x.size();
  std::size_t M = N * 2;
  double sum = 0.0;
  for(std::size_t i = 0; i < N; ++i)
    sum += double(i+1) * x[i];
  vect_n<double> result(M);
  for(std::size_t i = 0; i < M; ++i)
    result[i] = double(i+1) * sum - 1.0;
  return result;
};

void p02_j(mat<double,mat_structure::rectangular>& J, const vect_n<double>& x, const vect_n<double>& f) {
  ++gradCount;
  std::size_t M = f.size();
  std::size_t N = x.size();
  J.set_col_count(N);
  J.set_row_count(M);
  for(std::size_t i = 0; i < M; ++i)
    for(std::size_t j = 0; j < N; ++j)
      J(i,j) = double((i+1) * (j+1));
};

vect_n<double> p02_sol = vect_n<double>(10,6.0 / double( 41 * 11 * 10 ) );
vect_n<double> p02_start = vect_n<double>(10,1.0);
vect_n<double> p02_lower = vect_n<double>(10,-2.0);
vect_n<double> p02_upper = vect_n<double>(10,2.0);
std::string p02_name = "Linear Function, rank 1";

vect_n<double> p04_f(const vect_n<double>& x) {
  ++evalCount;
  vect_n<double> result(2);
  result[0] = 10.0 * ( x[1] - x[0] * x[0] );
  result[1] = 1.0 - x[0];
  return result;
};

void p04_j(mat<double,mat_structure::rectangular>& J, const vect_n<double>& x, const vect_n<double>&) {
  ++gradCount;
  J.set_col_count(2);
  J.set_row_count(2);
  J(0,0) = -20.0 * x[0];
  J(0,1) = 10.0;
  J(1,0) = -1.0;
  J(1,1) = 0.0;
};

vect_n<double> p04_sol = vect_n<double>(2, 1.0);
vect_n<double> p04_start = vect_n<double>(2,-0.5);
vect_n<double> p04_lower = vect_n<double>(2,-1.0);
vect_n<double> p04_upper = vect_n<double>(2,2.0);
std::string p04_name = "Rosenbrock Function";

vect_n<double> p05_f(const vect_n<double>& x) {
  ++evalCount;
  vect_n<double> result(3);
  double tmp = atan2(x[1], x[0]);
  
  result[0] = 10.0 * ( x[2] - 10.0 * tmp );
  result[1] = 10.0 * ( sqrt( x[0] * x[0] + x[1] * x[1] ) - 1.0 );
  result[2] = x[2];
  return result;
};

void p05_j(mat<double,mat_structure::rectangular>& J, const vect_n<double>& x, const vect_n<double>&) {
  ++gradCount;
  J.set_col_count(3);
  J.set_row_count(3);
  
  J(0,0) = 50.0 * x[1] / ( M_PI * ( x[0] * x[0] + x[1] * x[1] ) );
  J(0,1) = -50.0 * x[0] / ( M_PI * ( x[0] * x[0] + x[1] * x[1] ) );
  J(0,2) = 10.0;
  
  J(1,0) = 10.0 * x[0] / ( sqrt( x[0] * x[0] + x[1] * x[1] ) );
  J(1,1) = 10.0 * x[1] / ( sqrt( x[0] * x[0] + x[1] * x[1] ) );
  J(1,2) = 0.0;
  
  J(2,0) = 0.0;
  J(2,1) = 0.0;
  J(2,2) = 1.0;
};

vect_n<double> p05_sol = vect_n<double>(1.0,0.0,0.0);
vect_n<double> p05_start = vect_n<double>(-1.0,0.0,0.0);
vect_n<double> p05_lower = vect_n<double>(-2.0,-1.0,-1.0);
vect_n<double> p05_upper = vect_n<double>(2.0,1.0,1.0);
std::string p05_name = "Helical-valley Function";

vect_n<double> p06_f(const vect_n<double>& x) {
  ++evalCount;
  vect_n<double> result(4);
  result[0] = x[0] + 10.0 * x[1];
  result[1] = sqrt( 5.0 ) * ( x[2] - x[3] );
  result[2] = ( x[1] - 2.0 * x[2] ) * ( x[1] - 2.0 * x[2] );
  result[3] = sqrt( 10.0 ) * ( x[0] - x[3] ) * ( x[0] - x[3] );
  return result;
};

void p06_j(mat<double,mat_structure::rectangular>& J, const vect_n<double>& x, const vect_n<double>&) {
  ++gradCount;
  J.set_col_count(4);
  J.set_row_count(4);
  J(0,0) = 1.0;
  J(0,1) = 10.0;
  J(0,2) = 0.0;
  J(0,3) = 0.0;
  
  J(1,0) = 0.0;
  J(1,1) = 0.0;
  J(1,2) = sqrt( 5.0 );
  J(1,3) = -sqrt( 5.0 );
  
  J(2,0) = 0.0;
  J(2,1) = 2.0 * (x[1] - 2.0 * x[2]);
  J(2,2) = -4.0 * (x[1] - 2.0 * x[2]);
  J(2,3) = 0.0;
  
  J(3,0) = 2.0 * sqrt( 10.0 ) * ( x[0] - x[3] );
  J(3,1) = 0.0;
  J(3,2) = 0.0;
  J(3,3) = -2.0 * sqrt( 10.0 ) * ( x[0] - x[3] );
};

vect_n<double> p06_sol = vect_n<double>(4, 0.0);
vect_n<double> p06_start = vect_n<double>(3.0,-1.0,0.0,1.0);
vect_n<double> p06_lower = vect_n<double>(-4.0,-2.0,-2.0,-2.0);
vect_n<double> p06_upper = vect_n<double>(4.0,2.0,2.0,2.0);
std::string p06_name = "Powell singular Function";

vect_n<double> p07_f(const vect_n<double>& x) {
  ++evalCount;
  vect_n<double> result(2);
  result[0] = -13.0 + x[0] + ( ( 5.0 - x[1] ) * x[1] - 2.0 ) * x[1];
  result[1] = -29.0 + x[0] + ( ( 1.0 + x[1] ) * x[1] - 14.0 ) * x[1];
  return result;
};

void p07_j(mat<double,mat_structure::rectangular>& J, const vect_n<double>& x, const vect_n<double>&) {
  ++gradCount;
  J.set_col_count(2);
  J.set_row_count(2);
  J(0,0) = 1.0;
  J(0,1) = x[1] * (10.0 - 3.0 * x[1]) - 2.0;
  J(1,0) = 1.0;
  J(1,1) = x[1] * (2.0 + 3.0 * x[1]) - 14.0;
};

vect_n<double> p07_sol = vect_n<double>(vect<double,2>(5.0, 4.0));
vect_n<double> p07_start = vect_n<double>(vect<double,2>(0.5, -2.0));
vect_n<double> p07_lower = vect_n<double>(vect<double,2>(0.0,-5.0));
vect_n<double> p07_upper = vect_n<double>(vect<double,2>(8.0,8.0));
std::string p07_name = "Freudenstein-Roth Function";




int main() {
  
  typedef vect_n<double> (*FunctionPtr)(const vect_n<double>&);
  typedef void (*JacFunctionPtr)(mat<double,mat_structure::rectangular>&,const vect_n<double>&,const vect_n<double>&);
  
  std::vector< FunctionPtr > funcs;
  std::vector< JacFunctionPtr > func_jacs;
  std::vector< vect_n<double> > func_sols;
  std::vector< vect_n<double> > func_starts;
  std::vector< vect_n<double> > func_lowers;
  std::vector< vect_n<double> > func_uppers;
  std::vector< std::string > func_names;
  
  funcs.push_back(&p01_f);
  func_jacs.push_back(&p01_j);
  func_sols.push_back(p01_sol);
  func_starts.push_back(p01_start);
  func_lowers.push_back(p01_lower);
  func_uppers.push_back(p01_upper);
  func_names.push_back(p01_name);
  
  funcs.push_back(&p02_f);
  func_jacs.push_back(&p02_j);
  func_sols.push_back(p02_sol);
  func_starts.push_back(p02_start);
  func_lowers.push_back(p02_lower);
  func_uppers.push_back(p02_upper);
  func_names.push_back(p02_name);
  
  funcs.push_back(&p04_f);
  func_jacs.push_back(&p04_j);
  func_sols.push_back(p04_sol);
  func_starts.push_back(p04_start);
  func_lowers.push_back(p04_lower);
  func_uppers.push_back(p04_upper);
  func_names.push_back(p04_name);
  
  funcs.push_back(&p05_f);
  func_jacs.push_back(&p05_j);
  func_sols.push_back(p05_sol);
  func_starts.push_back(p05_start);
  func_lowers.push_back(p05_lower);
  func_uppers.push_back(p05_upper);
  func_names.push_back(p05_name);
  
  funcs.push_back(&p06_f);
  func_jacs.push_back(&p06_j);
  func_sols.push_back(p06_sol);
  func_starts.push_back(p06_start);
  func_lowers.push_back(p06_lower);
  func_uppers.push_back(p06_upper);
  func_names.push_back(p06_name);
  
  funcs.push_back(&p07_f);
  func_jacs.push_back(&p07_j);
  func_sols.push_back(p07_sol);
  func_starts.push_back(p07_start);
  func_lowers.push_back(p07_lower);
  func_uppers.push_back(p07_upper);
  func_names.push_back(p07_name);
  
  vect_n<double> x;
  vect_n<double> y;
  for(std::size_t i = 0; i < funcs.size(); ++i) {
    std::cout << "****************************************************************************************" << std::endl;
    std::cout << "Searching for least-square of function #" << i << " called: '" << func_names[i] << "'" << std::endl;
    
    x = func_starts[i];
    y = funcs[i](x); y *= 0.0;
    evalCount = 0; gradCount = 0;
    try {
      optim::gauss_newton_nllsq(funcs[i],func_jacs[i],x,y,200,1e-8);
      std::cout << "  Gauss-Newton method gives:\n"
                << "    x = " << x << " with error = " << norm(x - func_sols[i]) << "\n"
	        << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << std::endl;
    } catch(std::exception& e) {
      std::cout << "  Gauss-Newton method failed with error: " << e.what() << std::endl;
    };
    
    x = func_starts[i];
    y = funcs[i](x); y *= 0.0;
    evalCount = 0; gradCount = 0;
    try {
      optim::jacobian_transpose_nllsq(funcs[i],func_jacs[i],x,y,500,1e-8);
      std::cout << "  Jacobian-Transpose method gives:\n"
                << "    x = " << x << " with error = " << norm(x - func_sols[i]) << "\n"
	        << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << std::endl;
    } catch(std::exception& e) {
      std::cout << "  Jacobian-Transpose method failed with error: " << e.what() << std::endl;
    };
    
    x = func_starts[i];
    y = funcs[i](x); y *= 0.0;
    evalCount = 0; gradCount = 0;
    try {
      optim::levenberg_marquardt_nllsq(funcs[i],func_jacs[i],x,y,100,1e-4,1e-14,1e-8,1e-14);
      std::cout << "  Levenberg-Marquardt method gives:\n"
                << "    x = " << x << " with error = " << norm(x - func_sols[i]) << "\n"
                << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << std::endl;
    } catch(std::exception& e) {
      std::cout << "  Levenberg-Marquardt method failed with error: " << e.what() << std::endl;
    };
    
    
    std::cout << "*********** with box-constraints **************************************************" << std::endl;
    
    x = func_starts[i];
    y = funcs[i](x); y *= 0.0;
    evalCount = 0; gradCount = 0;
    try {
      optim::limited_gauss_newton_nllsq(funcs[i],func_jacs[i],x,y,200,boost::bind(optim::box_limit_function< vect_n<double> >,_1,_2,func_lowers[i],func_uppers[i]),1e-8);
      std::cout << "  Box-limited Gauss-Newton method gives:\n"
                << "    x = " << x << " with error = " << norm(x - func_sols[i]) << "\n"
	        << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << std::endl;
    } catch(std::exception& e) {
      std::cout << "  Box-limited Gauss-Newton method failed with error: " << e.what() << std::endl;
    };
    
    x = func_starts[i];
    y = funcs[i](x); y *= 0.0;
    evalCount = 0; gradCount = 0;
    try {
      optim::limited_jacobian_transpose_nllsq(funcs[i],func_jacs[i],x,y,500,boost::bind(optim::box_limit_function< vect_n<double> >,_1,_2,func_lowers[i],func_uppers[i]),1e-8);
      std::cout << "  Box-limited Jacobian-Transpose method gives:\n"
                << "    x = " << x << " with error = " << norm(x - func_sols[i]) << "\n"
	        << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << std::endl;
    } catch(std::exception& e) {
      std::cout << "  Box-limited Jacobian-Transpose method failed with error: " << e.what() << std::endl;
    };
    
    x = func_starts[i];
    y = funcs[i](x); y *= 0.0;
    evalCount = 0; gradCount = 0;
    try {
      optim::limited_levenberg_marquardt_nllsq(funcs[i],func_jacs[i],x,y,boost::bind(optim::box_limit_function< vect_n<double> >,_1,_2,func_lowers[i],func_uppers[i]),100,1e-4,1e-14,1e-8,1e-14);
      std::cout << "  Box-limited Levenberg-Marquardt method gives:\n"
                << "    x = " << x << " with error = " << norm(x - func_sols[i]) << "\n"
                << "    eval-count = " << evalCount << " and grad-eval-count = " << gradCount << std::endl;
    } catch(std::exception& e) {
      std::cout << "  Box-limited Levenberg-Marquardt method failed with error: " << e.what() << std::endl;
    };
    
    
    std::cout << "****************************************************************************************" << std::endl;
  };
  
  
  return 0;
};









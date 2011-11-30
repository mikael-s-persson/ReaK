
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

#include "conjugate_gradient_methods.hpp"
#include "nelder_mead_method.hpp"
#include "bfgs_methods.hpp"
#include "trust_region_search.hpp"

#include "path_planning/global_rng.hpp"

#include <iostream>

static int evalCount;
static int gradCount;

double max_truss_section_stress(double x) {
  double sigma1 = 0.8165 / x;
  double sigma2 = 1.1154 / (1 - x);
  evalCount++;
  return (sigma1 > sigma2 ? sigma1 : sigma2);  
};


double banana_function(const ReaK::vect<double,2>& x) {
  evalCount++;
  return (1.0 - x[0]) * (1.0 - x[0]) + 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
};

ReaK::vect<double,2> banana_function_grad(const ReaK::vect<double,2>& x) {
  gradCount++;
  return ReaK::vect<double,2>(2.0 * (x[0] - 1.0) - 400.0 * x[0] * (x[1] - x[0] * x[0]),
                              200.0 * (x[1] - x[0] * x[0]));
};

double easy_function(const ReaK::vect<double,2>& x) {
  evalCount++;
  static ReaK::mat<double,ReaK::mat_structure::symmetric> Q(10.0, -2.0, 1.0);
  return x * (Q * x);
};

ReaK::vect<double,2> easy_function_grad(const ReaK::vect<double,2>& x) {
  gradCount++;
  static ReaK::mat<double,ReaK::mat_structure::symmetric> Q(10.0, -2.0, 1.0);
  return 2.0 * (Q * x);
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

  ReaK::mat<double,ReaK::mat_structure::rectangular> J(1,1);
  double x = 0.5;
    
  std::cout << "Jacobian at point x = 0.5 is:" << std::endl;
  ReaK::optim::compute_jacobian_2pts_forward(max_truss_section_stress,x,max_truss_section_stress(x),J);
  std::cout << "  2-pts Forward: " << J << std::endl;
  
  ReaK::optim::compute_jacobian_2pts_central(max_truss_section_stress,x,max_truss_section_stress(x),J);
  std::cout << "  2-pts Central: " << J << std::endl;
  
  ReaK::optim::compute_jacobian_5pts_central(max_truss_section_stress,x,max_truss_section_stress(x),J);
  std::cout << "  5-pts Central: " << J << std::endl;
  
  std::cout << std::endl << std::endl;
  std::cout << "Testing optimization methods on the Banana Function (optimum at (1,1), with value " << banana_function(ReaK::vect<double,2>(1.0,1.0)) << ")" << std::endl;
  
//   evalCount = 0;
  ReaK::vect<double,2> x_2D(0.5,1.0);
//   std::cout << "  Nelder-Mead method started at " << x_2D << std::endl;
  double initial_spread = 0.5;
//   for(unsigned int i = 0; i < 20; ++i) {
//     ReaK::optim::nelder_mead_method(banana_function, x_2D, initial_spread, ReaK::pp::get_global_rng(), 1e-7);
//     initial_spread *= 0.8;
//   };
//   std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations." << std::endl;
//   std::cout << "    function gives: " << banana_function(x_2D) << " gradient gives: " << banana_function_grad(x_2D) << std::endl;
//   
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  BFGS method started at " << x_2D << std::endl;
  ReaK::optim::bfgs_method(banana_function,banana_function_grad, x_2D, 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << banana_function(x_2D) << " gradient gives: " << banana_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  DFP method started at " << x_2D << std::endl;
  ReaK::optim::dfp_method(banana_function,banana_function_grad, x_2D, 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << banana_function(x_2D) << " gradient gives: " << banana_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  Broyden-class method started at " << x_2D << std::endl;
  ReaK::optim::broyden_class_method(banana_function,banana_function_grad, x_2D, 0.5, 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << banana_function(x_2D) << " gradient gives: " << banana_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  FR CG method started at " << x_2D << std::endl;
  ReaK::optim::non_linear_conj_grad_method(banana_function,banana_function_grad, x_2D, 
					   ReaK::optim::fletcher_reeves_beta(), 
					   ReaK::optim::line_search_expand_and_zoom<double>(1e-4,0.1), 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << banana_function(x_2D) << " gradient gives: " << banana_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  PR CG method started at " << x_2D << std::endl;
  ReaK::optim::non_linear_conj_grad_method(banana_function,banana_function_grad, x_2D, 
					   ReaK::optim::polak_ribiere_beta(), 
					   ReaK::optim::line_search_expand_and_zoom<double>(1e-4,0.1), 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << banana_function(x_2D) << " gradient gives: " << banana_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  HS CG method started at " << x_2D << std::endl;
  ReaK::optim::non_linear_conj_grad_method(banana_function,banana_function_grad, x_2D, 
					   ReaK::optim::hestenes_stiefel_beta(), 
					   ReaK::optim::line_search_expand_and_zoom<double>(1e-4,0.1), 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << banana_function(x_2D) << " gradient gives: " << banana_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  DY CG method started at " << x_2D << std::endl;
  ReaK::optim::non_linear_conj_grad_method(banana_function,banana_function_grad, x_2D, 
					   ReaK::optim::dai_yuan_beta(), 
					   ReaK::optim::line_search_expand_and_zoom<double>(1e-4,0.1), 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << banana_function(x_2D) << " gradient gives: " << banana_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  HZ CG method started at " << x_2D << std::endl;
  ReaK::optim::non_linear_conj_grad_method(banana_function,banana_function_grad, x_2D, 
					   ReaK::optim::hager_zhang_beta(), 
					   ReaK::optim::line_search_expand_and_zoom<double>(1e-4,0.1), 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << banana_function(x_2D) << " gradient gives: " << banana_function_grad(x_2D) << std::endl;
  
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  SR1 trust-region method started at " << x_2D << std::endl;
  ReaK::optim::quasi_newton_trust_region(banana_function,banana_function_grad, x_2D, 0.5, 
					 ReaK::optim::trust_region_solver_dogleg(),
					 ReaK::optim::hessian_update_sr1(), 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << banana_function(x_2D) << " gradient gives: " << banana_function_grad(x_2D) << std::endl;
  
  
  
  
  std::cout << std::endl << std::endl;
  std::cout << "Testing optimization methods on the Easy Function (optimum at (0,0), with value 0.0)" << std::endl;
  
  evalCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  Nelder-Mead method started at " << x_2D << std::endl;
  initial_spread = 0.5;
  for(unsigned int i = 0; i < 20; ++i) {
    ReaK::optim::nelder_mead_method(easy_function, x_2D, initial_spread, ReaK::pp::get_global_rng(), 1e-6);
    initial_spread *= 0.8;
  };
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations." << std::endl;
  std::cout << "    function gives: " << easy_function(x_2D) << " gradient gives: " << easy_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  BFGS method started at " << x_2D << std::endl;
  ReaK::optim::bfgs_method(easy_function,easy_function_grad, x_2D, 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << easy_function(x_2D) << " gradient gives: " << easy_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  DFP method started at " << x_2D << std::endl;
  ReaK::optim::dfp_method(easy_function,easy_function_grad, x_2D, 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << easy_function(x_2D) << " gradient gives: " << easy_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  Broyden-class method started at " << x_2D << std::endl;
  ReaK::optim::broyden_class_method(easy_function,easy_function_grad, x_2D, 0.5, 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << easy_function(x_2D) << " gradient gives: " << easy_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  FR CG method started at " << x_2D << std::endl;
  ReaK::optim::non_linear_conj_grad_method(easy_function,easy_function_grad, x_2D, 
					   ReaK::optim::fletcher_reeves_beta(), 
					   ReaK::optim::line_search_expand_and_zoom<double>(1e-4,0.1), 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << easy_function(x_2D) << " gradient gives: " << easy_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  PR CG method started at " << x_2D << std::endl;
  ReaK::optim::non_linear_conj_grad_method(easy_function,easy_function_grad, x_2D, 
					   ReaK::optim::polak_ribiere_beta(), 
					   ReaK::optim::line_search_expand_and_zoom<double>(1e-4,0.1), 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << easy_function(x_2D) << " gradient gives: " << easy_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  HS CG method started at " << x_2D << std::endl;
  ReaK::optim::non_linear_conj_grad_method(easy_function,easy_function_grad, x_2D, 
					   ReaK::optim::hestenes_stiefel_beta(), 
					   ReaK::optim::line_search_expand_and_zoom<double>(1e-4,0.1), 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << easy_function(x_2D) << " gradient gives: " << easy_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  DY CG method started at " << x_2D << std::endl;
  ReaK::optim::non_linear_conj_grad_method(easy_function,easy_function_grad, x_2D, 
					   ReaK::optim::dai_yuan_beta(), 
					   ReaK::optim::line_search_expand_and_zoom<double>(1e-4,0.1), 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << easy_function(x_2D) << " gradient gives: " << easy_function_grad(x_2D) << std::endl;
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  HZ CG method started at " << x_2D << std::endl;
  ReaK::optim::non_linear_conj_grad_method(easy_function,easy_function_grad, x_2D, 
					   ReaK::optim::hager_zhang_beta(), 
					   ReaK::optim::line_search_expand_and_zoom<double>(1e-4,0.1), 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << easy_function(x_2D) << " gradient gives: " << easy_function_grad(x_2D) << std::endl;
  
  
  evalCount = 0;
  gradCount = 0;
  x_2D = ReaK::vect<double,2>(0.5,1.0);
  std::cout << "  SR1 trust-region method started at " << x_2D << std::endl;
  ReaK::optim::quasi_newton_trust_region(easy_function,easy_function_grad, x_2D, 0.5, 
					 ReaK::optim::trust_region_solver_dogleg(),
					 ReaK::optim::hessian_update_sr1(), 1e-7);
  std::cout << "    found optimum: " << x_2D << " after " << evalCount << " function evaluations and " << gradCount << " gradient evaluations." << std::endl;
  std::cout << "    function gives: " << easy_function(x_2D) << " gradient gives: " << easy_function_grad(x_2D) << std::endl;
  
  
  
  return 0;
};









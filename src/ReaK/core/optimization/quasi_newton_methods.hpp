/**
 * \file quasi_newton_methods.hpp
 *
 * The following library provides methods to perform non-linear optimization based on Quasi-Newton methods.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
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

#ifndef REAK_QUASI_NEWTON_METHODS_HPP
#define REAK_QUASI_NEWTON_METHODS_HPP

#include "base/defs.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_num_exceptions.hpp"


namespace ReaK {
  
  
namespace optim {


/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search 
 * direction, an approximation of the Hessian (without restarts), and using a line-search approach 
 * (satisfying the strong Wolfe conditions). Note that this function is the underlying optimization 
 * loop used by most quasi-Newton methods (see bfgs_method, dfp_method, and broyden_class_method).
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam LineSearcher A functor type that can perform a line-search which satisfy the strong Wolfe conditions.
 * \tparam InvHessianUpdater A functor type that can update the inverse of a Hessian given the changes in the solution and the function gradient.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param get_alpha The functor that can perform a line-search which satisfy the strong Wolfe conditions.
 * \param update_inv_hessian The functor that can update the inverse of a Hessian given the changes in the solution and the function gradient.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename Vector, typename LineSearcher, typename InvHessianUpdater>
void quasi_newton_line_search(Function f, GradFunction df, Vector& x, LineSearcher get_alpha, InvHessianUpdater update_inv_hessian, typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  typedef typename vect_traits<Vector>::value_type ValueType;
  using std::sqrt; using std::fabs;
  
  ValueType abs_tol = sqrt(x * x) * tol;
  ValueType x_value = f(x);
  Vector x_grad = -df(x);
  abs_tol += sqrt(x_grad * x_grad) * tol;
  
  Vector p = x_grad;
  ValueType alpha = ValueType(1.0);
  ValueType pxg = p * x_grad;
  if((f(x + p) > x_value - ValueType(1e-4) * pxg) || (fabs(p * df(x + p)) > ValueType(0.9) * fabs(pxg)))
    alpha = get_alpha(f,df,ValueType(0.0),ValueType(2.0),x,p,tol);
  
  Vector s = alpha * p;
  Vector y = x_grad;
  x += s;
  x_value = f(x);
  x_grad = -df(x);
  y -= x_grad;
  mat<ValueType,mat_structure::square> H(mat<ValueType,mat_structure::scalar>(x.size(),(y * s) / (y * y)));
  
  while( norm(x_grad) > abs_tol ) {
    update_inv_hessian(H,s,y);
    p = H * x_grad;
    // check Wolfe for alpha 1.0
    alpha = ValueType(1.0);
    pxg = p * x_grad;
    if((f(x + p) > x_value - ValueType(1e-4) * pxg) || (fabs(p * df(x + p)) > ValueType(0.9) * fabs(pxg)))
      alpha = get_alpha(f,df,ValueType(0.0),ValueType(2.0),x,p,tol);
    s = alpha * p;
    x += s;
    y = x_grad;
    x_value = f(x);
    x_grad = -df(x);
    y -= x_grad;
  };
};







/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search 
 * direction, an approximation of the Hessian (without restarts), and using a trust-region approach. 
 * Note that this function is the underlying optimization loop used by most quasi-Newton 
 * methods (see bfgs_method, dfp_method, and broyden_class_method).
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam TrustRegionSolver A functor type that can find the solution to the approximate quadratic over the current trust-region.
 * \tparam HessianUpdater A functor type that can update the Hessian given the changes in the solution and the function gradient.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param solve_step The functor that can find the solution to the approximate quadratic over the current trust-region.
 * \param update_hessian The functor that can update the Hessian given the changes in the solution and the function gradient.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 * \param eta The tolerance on the ratio between actual reduction and predicted reduction in order to accept a given step.
 * \param r_tol The tolerance on the alignment of the step and the residual in order to trigger an update of the approximate Hessian matrix.
 */
template <typename Function, typename GradFunction, typename Vector, typename TrustRegionSolver, typename HessianUpdater>
void quasi_newton_trust_region(Function f, GradFunction df, Vector& x, typename vect_traits<Vector>::value_type max_radius, TrustRegionSolver solve_step, HessianUpdater update_hessian, typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6), typename vect_traits<Vector>::value_type eta = typename vect_traits<Vector>::value_type(1e-4), typename vect_traits<Vector>::value_type r_tol = typename vect_traits<Vector>::value_type(1e-6)) {
  typedef typename vect_traits<Vector>::value_type ValueType;
  using std::sqrt; using std::fabs;
  
  ValueType abs_tol = sqrt(x * x) * tol;
  ValueType radius = ValueType(0.5) * max_radius;

  ValueType x_value = f(x);
  Vector x_grad = df(x);
  abs_tol += sqrt(x_grad * x_grad) * tol;
  
  Vector p = -x_grad;
  ValueType norm_p;
  mat<ValueType,mat_structure::square> B(mat<ValueType,mat_structure::identity>(x.size()));
  solve_step(x_grad,B,p,norm_p,radius,tol);
  Vector xt = x; xt += p;
  ValueType xt_value = f(xt);
  Vector xt_grad = df(xt);
  ValueType aredux = x_value - xt_value;
  ValueType predux = -(x_grad * p + ValueType(0.5) * (p * (B * p)));
  
  Vector y = xt_grad; y -= x_grad;
  B = mat<ValueType,mat_structure::scalar>(x.size(),(y * p) / (p * p));
  Vector y_Bp = y; y_Bp -= B * p;
  
  while( true ) {
    if( fabs(p * y_Bp) >= r_tol * norm_p * norm(y_Bp) )
      update_hessian(B,p,y);
    
    ValueType ratio = aredux / predux;
    if( ratio > ValueType(0.75) ) {
      if(norm_p > ValueType(0.8) * radius) {
	radius *= ValueType(2.0);
	if(radius > max_radius)
	  radius = max_radius;
      };
    } else if( ratio < ValueType(0.1) ) {
      radius *= ValueType(0.5);
    };
    if( ratio > eta ) {
      x = xt;
      x_value = xt_value;
      x_grad = xt_grad;
    };
    if(norm(x_grad) < abs_tol)
      return;
    solve_step(x_grad,B,p,norm_p,radius,tol);
    xt = x; xt += p;
    xt_value = f(xt);
    xt_grad = df(xt);
    aredux = x_value - xt_value;
    predux = -(x_grad * p + ValueType(0.5) * (p * (B * p)));
    y = xt_grad; y -= x_grad;
    y_Bp = y; y_Bp -= B * p;
  };
};





};

};


#endif


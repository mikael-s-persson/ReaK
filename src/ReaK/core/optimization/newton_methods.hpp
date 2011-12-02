/**
 * \file newton_methods.hpp
 *
 * The following library provides methods to perform non-linear optimization based on Newton methods.
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

#ifndef REAK_NEWTON_METHODS_HPP
#define REAK_NEWTON_METHODS_HPP

#include "base/defs.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_num_exceptions.hpp"

#include "newton_search_directions.hpp"
#include "limit_functions.hpp"
#include "line_search.hpp"
#include "trust_region_search.hpp"

namespace ReaK {
  
  
namespace optim {
  
  
namespace detail {
  
  
  template <typename Function, typename GradFunction, typename HessianFunction, 
            typename Vector, typename LineSearcher, 
	    typename LimitFunction, typename NewtonDirectioner>
  void newton_method_ls_impl(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, LineSearcher get_alpha, LimitFunction impose_limits, NewtonDirectioner get_direction, typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
    typedef typename vect_traits<Vector>::value_type ValueType;
    using std::sqrt; using std::fabs;
    
    ValueType abs_tol = (x * x) * tol;
    ValueType x_value = f(x);
    Vector x_grad = df(x);
    ValueType x_grad_norm_sqr = x_grad * x_grad;
    ValueType abs_grad_tol = x_grad_norm_sqr * tol;
  
    Vector p = x_grad;
    mat<ValueType,mat_structure::symmetric> H(x.size());
    fill_hessian(H,x,x_value,x_grad);
  
    while( x_grad_norm_sqr > abs_grad_tol ) {
      get_direction(H,x_grad,p,tol);
      // check Wolfe for alpha 1.0
      ValueType alpha = ValueType(1.0);
      ValueType pxg = p * x_grad;
      if((f(x + p) > x_value + ValueType(1e-4) * pxg) || (fabs(p * df(x + p)) > ValueType(0.4) * fabs(pxg)))
        alpha = get_alpha(f,df,ValueType(0.0),ValueType(2.0),x,p,tol);
      p *= alpha;
      impose_limits(x,p);
      if((p * p) > abs_tol)
        return;
    
      x += p;
      x_value = f(x);
      x_grad = df(x);
      x_grad_norm_sqr = x_grad * x_grad;
      fill_hessian(H,x,x_value,x_grad);
    };
  };
  
  
  
  template <typename Function, typename GradFunction, typename HessianFunction, typename Vector, typename TrustRegionSolver, typename LimitFunction>
  void newton_method_tr_impl(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, typename vect_traits<Vector>::value_type max_radius, TrustRegionSolver solve_step, LimitFunction impose_limits, typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6), typename vect_traits<Vector>::value_type eta = typename vect_traits<Vector>::value_type(1e-4), typename vect_traits<Vector>::value_type r_tol = typename vect_traits<Vector>::value_type(1e-6)) {
    typedef typename vect_traits<Vector>::value_type ValueType;
    using std::sqrt; using std::fabs;
  
    ValueType abs_tol = (x * x) * tol;
    ValueType radius = ValueType(0.5) * max_radius;
    ValueType x_value = f(x);
    Vector x_grad = df(x);
    ValueType x_grad_norm_sqr = x_grad * x_grad;
    ValueType abs_grad_tol = x_grad_norm_sqr * tol;
  
    Vector p = x_grad;
    ValueType norm_p = std::numeric_limits<ValueType>::max();
    mat<ValueType,mat_structure::symmetric> H(x.size());
    fill_hessian(H,x,x_value,x_grad);
    
    Vector xt = x;
    
    while( (norm_p > abs_tol) && (x_grad_norm_sqr > abs_grad_tol) ) {
      solve_step(x_grad,H,p,norm_p,radius,tol);
      impose_limits(x,p);
      xt = x; xt += p;
      ValueType xt_value = f(xt);
      ValueType aredux = x_value - xt_value;
      ValueType predux = -(x_grad * p + ValueType(0.5) * (p * (B * p)));
    
      
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
      if( ratio > eta ) {  //the step is accepted.
        x = xt;
        x_value = xt_value;
        x_grad = df(x);
        x_grad_norm_sqr = x_grad * x_grad;
        fill_hessian(H,x,x_value,x_grad);
      };
    };
  };

  
  
};



/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a newton search 
 * direction and using a line-search approach (satisfying the strong Wolfe conditions).
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename Vector>
void newton_method_ls(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, 
		      typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  detail::newton_method_ls_impl(
    f,df,fill_hessian,x,
    line_search_expand_and_zoom<typename vect_traits<Vector>::value_type>(1e-4,0.4),
    no_limit_functor(),newton_direction,tol);
  
};



/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a regularized 
 * newton search direction and using a line-search approach (satisfying the strong Wolfe conditions).
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param tau The initial relative damping factor to regularize the Hessian matrix.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename Vector>
void reg_newton_method_ls(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, 
			  typename vect_traits<Vector>::value_type tau = typename vect_traits<Vector>::value_type(1e-3), 
			  typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  detail::newton_method_ls_impl(
    f,df,fill_hessian,x,
    line_search_expand_and_zoom<typename vect_traits<Vector>::value_type>(1e-4,0.4),
    no_limit_functor(),regularized_newton_directioner<typename vect_traits<Vector>::value_type>(tau),tol);
  
};



/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a newton search 
 * direction, a function to impose simple limits on the search domain, and a line-search 
 * approach (satisfying the strong Wolfe conditions).
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying independent vector-space.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename Vector, typename LimitFunction>
void limited_newton_method_ls(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, LimitFunction impose_limits, typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  detail::newton_method_ls_impl(
    f,df,fill_hessian,x,
    line_search_expand_and_zoom<typename vect_traits<Vector>::value_type>(1e-4,0.4),
    impose_limits,newton_direction,tol);
  
};



/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a regularized 
 * newton search direction, a function to impose simple limits on the search domain, and a line-search 
 * approach (satisfying the strong Wolfe conditions).
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying independent vector-space.
 * \param tau The initial relative damping factor to regularize the Hessian matrix.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename Vector, typename LimitFunction>
void limited_reg_newton_method_ls(Function f, GradFunction df, HessianFunction fill_hessian, 
				  Vector& x, LimitFunction impose_limits, 
				  typename vect_traits<Vector>::value_type tau = typename vect_traits<Vector>::value_type(1e-3), 
				  typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  detail::newton_method_ls_impl(
    f,df,fill_hessian,x,
    line_search_expand_and_zoom<typename vect_traits<Vector>::value_type>(1e-4,0.4),
    impose_limits,regularized_newton_directioner<typename vect_traits<Vector>::value_type>(tau),tol);
  
};




/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a newton search 
 * direction and using a trust-region approach (dogleg solver).
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename Vector>
void newton_method_tr(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, 
		      typename vect_traits<Vector>::value_type max_radius,
		      typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  detail::newton_method_tr_impl(
    f,df,fill_hessian,x,max_radius,
    trust_region_solver_dogleg(),
    no_limit_functor(),tol);
  
};

/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a 
 * regularized newton search direction and using a trust-region approach (dogleg solver).
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param tau The initial relative damping factor to regularize the Hessian matrix.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename Vector>
void reg_newton_method_tr(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, 
		          typename vect_traits<Vector>::value_type max_radius,
		          typename vect_traits<Vector>::value_type tau = typename vect_traits<Vector>::value_type(1e-3),
		          typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  detail::newton_method_tr_impl(
    f,df,fill_hessian,x,max_radius,
    trust_region_solver_dogleg_reg<typename vect_traits<Vector>::value_type>(tau),
    no_limit_functor(),tol);
  
};


/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a newton search 
 * direction, a function to impose simple limits on the search domain and using a trust-region 
 * approach (dogleg solver).
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying independent vector-space.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename Vector, typename LimitFunction>
void limited_newton_method_tr(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, 
		              typename vect_traits<Vector>::value_type max_radius, LimitFunction impose_limits, 
		              typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  detail::newton_method_tr_impl(
    f,df,fill_hessian,x,max_radius,
    trust_region_solver_dogleg(),
    impose_limits,tol);
  
};

/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a 
 * regularized newton search direction, a function to impose simple limits on the search domain
 * and using a trust-region approach (dogleg solver).
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying independent vector-space.
 * \param tau The initial relative damping factor to regularize the Hessian matrix.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename Vector, typename LimitFunction>
void limited_reg_newton_method_tr(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, 
		                  typename vect_traits<Vector>::value_type max_radius, LimitFunction impose_limits, 
		                  typename vect_traits<Vector>::value_type tau = typename vect_traits<Vector>::value_type(1e-3),
		                  typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  detail::newton_method_tr_impl(
    f,df,fill_hessian,x,max_radius,
    trust_region_solver_dogleg_reg<typename vect_traits<Vector>::value_type>(tau),
    impose_limits,tol);
  
};




};

};


#endif


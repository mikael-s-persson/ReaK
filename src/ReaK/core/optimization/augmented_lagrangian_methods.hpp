/**
 * \file augmented_lagrangian_methods.hpp
 *
 * The following library provides methods to perform constrained non-linear optimization using 
 * methods of augmented Lagrangians (bound-constrained lagrangians).
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date December 2011
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

#ifndef REAK_AUGMENTED_LAGRANGIAN_METHODS_HPP
#define REAK_AUGMENTED_LAGRANGIAN_METHODS_HPP

#include "base/defs.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_num_exceptions.hpp"

#include "newton_search_directions.hpp"
#include "limit_functions.hpp"
#include "line_search.hpp"
#include "trust_region_search.hpp"

#include "newton_methods.hpp"

namespace ReaK {
  
  
namespace optim {
  
  
namespace detail {
  
  
  template <typename Function, typename GradFunction, typename HessianFunction, 
            typename Vector, typename EqFunction, typename EqJacFunction, 
	    typename IneqFunction, typename IneqJacFunction,
	    typename TrustRegionSolver, typename LimitFunction>
  void bcl_newton_method_tr_impl(Function f, GradFunction df, HessianFunction fill_hessian,  
				 EqFunction g, EqJacFunction fill_g_jac,
				 IneqFunction h, IneqJacFunction fill_h_jac,
			         Vector& x, typename vect_traits<Vector>::value_type max_radius, unsigned int max_iter,
			         TrustRegionSolver solve_step, LimitFunction impose_limits, 
			         typename vect_traits<Vector>::value_type abs_tol = typename vect_traits<Vector>::value_type(1e-6), 
			         typename vect_traits<Vector>::value_type abs_grad_tol = typename vect_traits<Vector>::value_type(1e-6), 
				 typename vect_traits<Vector>::value_type kappa = typename vect_traits<Vector>::value_type(1e-4)) {
    typedef typename vect_traits<Vector>::value_type ValueType;
    typedef typename vect_traits<Vector>::size_type SizeType;
    using std::sqrt; using std::fabs;
    
    SizeType N = x.size();
    Vector h_value = h(x);
    Vector ht_value = h_value;
    Vector s = h_value;
    SizeType K = h_value.size();
    Vector g_value = g(x);
    Vector gt_value = g_value;
    SizeType M = g_value.size();
    
    if((M == 0) && (K == 0)) { //this means it is an unconstrained problem.
      newton_method_tr_impl(f,df,fill_hessian,x,max_radius,max_iter,solve_step,impose_limits,abs_tol,abs_grad_tol,kappa);
      return;
    };
      
    ValueType mu(10.0);
    
    mat<ValueType, mat_structure::rectangular> Jac_g(M,N);
    fill_g_jac(Jac_g,x,g_value);
    mat<ValueType, mat_structure::rectangular> Jac_h(K,N);
    fill_h_jac(Jac_h,x,h_value);
    mat<ValueType, mat_structure::rectangular> muJac_h(-mu * Jac_h);
  
    Vector l_g(M, ValueType(1.0)); //choose initial lambda vector to be 1.0.
    Vector l_h(K, ValueType(1.0));
  
    //ValueType c_norm_star = (norm_2(g_value) + norm_2(h_value)) / ValueType(K + M);
    ValueType eta_star = abs_tol;
    //ValueType norm_star = norm_2(x) / ValueType(N) + ValueType(1.0);
    ValueType omega_star = abs_tol;
  
    ValueType omega = 1.0 / mu;
    ValueType eta = pow(mu,ValueType(-0.1));
    
    ValueType radius = ValueType(0.5) * max_radius;
    ValueType x_value = f(x);
    Vector x_grad = df(x);
    
    ValueType L_value = x_value - l_g * g_value + ValueType(0.5) * mu * (g_value * g_value); 
                      // + 0.5 * mu * (h_value - s) * (h_value - s) - l_h * (h_value - s)
    Vector L_grad; L_grad.resize(N+K);
    L_grad[range(0,N-1)] = x_grad - l_g * Jac_g - l_h * Jac_h + mu * (g_value * Jac_g);
    L_grad[range(N,N+K-1)] = l_h; // - mu * (h_value - s)
    
    Vector p = L_grad;
    Vector p_x = x;
    Vector p_s = s;
    ValueType norm_p = std::numeric_limits<ValueType>::max();
    ValueType norm_p_total(0.0);
    mat<ValueType,mat_structure::symmetric> H(x.size());
    fill_hessian(H,x,x_value,x_grad);
    
    mat<ValueType,mat_structure::symmetric> JJ(transpose_view(Jac_g) * Jac_g + transpose_view(Jac_h) * Jac_h);
    JJ *= mu;
    mat<ValueType,mat_structure::symmetric> L_H = H; L_H += JJ;
    
    mat<ValueType,mat_structure::symmetric> La_H(N+K);
    sub(La_H)(range(0,N-1),range(0,N-1)) = H; sub(La_H)(range(0,N-1),range(0,N-1)) += JJ;
    sub(La_H)(range(0,N-1),range(N,N+K-1)) = transpose_view(muJac_h);
    sub(La_H)(range(N,N+K-1),range(N,N+K-1)) = mat<ValueType,mat_structure::scalar>(K,mu);
    
//     mat_vert_cat<
//       mat_const_ref_horiz_cat< mat<ValueType,mat_structure::symmetric>, mat_transpose_view< mat<ValueType,mat_structure::rectangular> > >,
//       mat_const_ref_horiz_cat< mat<ValueType,mat_structure::rectangular>, mat<ValueType,mat_structure::scalar> >
//     > La_H(
//       mat_const_ref_horiz_cat< mat<ValueType,mat_structure::symmetric>, mat_transpose_view< mat<ValueType,mat_structure::rectangular> > >(L_H,muJac_h_t),
//       mat_const_ref_horiz_cat< mat<ValueType,mat_structure::rectangular>, mat<ValueType,mat_structure::scalar> >(muJac_h,mu_mat)
//     );
      
    
    Vector xt = x;
    Vector st = s;
    
    ValueType c_norm = sqrt(g_value * g_value);
    ValueType c_norm_max(c_norm);
    
    unsigned int k = 0;
    
    while (true) {
    
      do {
        solve_step(L_grad,La_H,p,norm_p,radius,abs_tol);
	p_x = p[range(0,N-1)];
        p_s = p[range(N,N+K-1)];
	impose_limits(x,p_x);  //impose limits on the x-space
	for(SizeType i = 0; i < K; ++i) 
	  if(s[i] + p_s[i] < ValueType(0.0))
	    p_s[i] = -s[i];  //impose non-negativity on the s-space.
	p[range(0,N-1)] = p_x;
        p[range(N,N+K-1)] = p_s;
	xt = x; xt += p_x;
	norm_p = norm_2(p);
	st = s; st += p_s;
	gt_value = g(xt);
	ht_value = h(xt) - st;
	ValueType xt_value = f(xt);
	ValueType Lt_value = f(xt) - l_g * gt_value - l_h * ht_value + ValueType(0.5) * mu * (gt_value * gt_value + ht_value * ht_value);
        ValueType aredux = L_value - Lt_value;
        ValueType predux = -(L_grad * p + ValueType(0.5) * (p * (La_H * p)));
    
        ValueType ratio = aredux / predux;
        //RK_NOTICE(1," ratio = " << ratio << " kappa = " << kappa << " aredux = " << aredux << " norm_p = " << norm_p);
        if( ratio > ValueType(0.75) ) {
          if(norm_p > ValueType(0.8) * radius) {
            radius *= ValueType(2.0);
            if(radius > max_radius)
              radius = max_radius;
          };
        } else if( ratio < ValueType(0.1) ) {
          radius *= ValueType(0.5);
        };
	if( ratio > kappa ) {  //the step is accepted.
          x = xt;
	  norm_p_total += norm_p;
	  if(norm_p < omega * norm_p_total)
	    break;
          x_value = xt_value;
	  L_value = Lt_value;
	  x_grad = df(x);
	  s = st;
	  g_value = gt_value;
	  h_value = ht_value;
	  fill_h_jac(Jac_h,x,h_value + s);
	  muJac_h = -mu * Jac_h;
	  sub(La_H)(range(0,N-1),range(N,N+K-1)) = transpose_view(muJac_h);
	  fill_g_jac(Jac_g,x,g_value);
	  L_grad[range(0,N-1)] = x_grad - l_g * Jac_g - l_h * Jac_h + mu * g_value * Jac_g + mu * h_value * Jac_h;
	  L_grad[range(N,N+K-1)] = l_h - mu * h_value;
          fill_hessian(H,x,x_value,x_grad);
	  JJ = transpose_view(Jac_g) * Jac_g + transpose_view(Jac_h) * Jac_h; JJ *= mu;
	  sub(La_H)(range(0,N-1),range(0,N-1)) = H; sub(La_H)(range(0,N-1),range(0,N-1)) += JJ;
	  c_norm = sqrt(h_value * h_value + g_value * g_value);
	  if(c_norm_max < c_norm)
	    c_norm_max = c_norm;
	  //RK_NOTICE(1," L_value = " << L_value << " x = " << x << " c_norm = " << c_norm << " norm_p = " << norm_p);
        };
	if( ++k > max_iter )
	  throw maximum_iteration(max_iter);
      } while(true);
      
      //if( ++k > max_iter )
	//throw maximum_iteration(max_iter);
      
      if( c_norm < eta * c_norm_max ) {
      
        if((c_norm < eta_star) && (norm_p < omega_star)) //in theory, this test should involve the projected gradient, but I think that the step-size is sufficient.
	  return;
      
	l_g -= mu * g_value;
	l_h -= mu * h_value;
	L_value = x_value - l_g * g_value - l_h * h_value + ValueType(0.5) * mu * (g_value * g_value + h_value * h_value);
        eta *= pow(mu,ValueType(-0.9));
	if(eta < eta_star)
	  eta = eta_star;
        omega /= mu;
	if(omega < omega_star)
	  omega = omega_star;
      } else {
        mu *= ValueType(100.0);
	sub(La_H)(range(N,N+K-1),range(N,N+K-1)) = mat<ValueType,mat_structure::scalar>(K,mu);
	muJac_h = -mu * Jac_h;
	sub(La_H)(range(0,N-1),range(N,N+K-1)) = transpose_view(muJac_h);
        L_value = x_value - l_g * g_value - l_h * h_value + ValueType(0.5) * mu * (g_value * g_value + h_value * h_value);
        eta = pow(mu,ValueType(-0.1));
	if(eta * c_norm_max < eta_star)
	  eta = eta_star / c_norm_max;
        omega = 1.0 / mu;
	if(omega * norm_p_total < omega_star)
	  omega = omega_star / norm_p_total;
      };
      L_grad[range(0,N-1)] = x_grad - l_g * Jac_g - l_h * Jac_h + mu * g_value * Jac_g + mu * h_value * Jac_h;
      L_grad[range(N,N+K-1)] = l_h - mu * h_value;
      
    };
  };


  
  
  struct impose_positive_slacks {
    unsigned int K;
    impose_positive_slacks(unsigned int aK) : K(aK) { };
    template <typename Vector>
    void operator()(const Vector& x, Vector& dx) const {
      typedef typename vect_traits<Vector>::size_type SizeType;
      typedef typename vect_traits<Vector>::value_type ValueType;
      for(SizeType i = x.size() - K; i < x.size(); ++i) {
	if(x[i] + dx[i] < ValueType(0.0))
	  dx[i] = -x[i];
      };
    };
  };
  
  
};


/**
 * This functor class can be used in-place of a constraint functor to represent "no-constraints".
 */
struct no_constraint_functor {
  template <typename Vector>
  Vector operator()(const Vector&) const {
    return Vector(0);
  };
};
  
/**
 * This functor class can be used in-place of a constraint Jacobian functor to represent "no-constraints".
 */
struct no_constraint_jac_functor {
  template <typename Matrix, typename Vector>
  void operator()(Matrix&, const Vector&, const Vector&) const { };
};



/**
 * This functor is a factory class to construct a Newton-method optimizer routine that uses a 
 * trust-region approach and incorporates equality and inequality constraints via the augmented 
 * lagrangian method. Use make_newton_method_tr to
 * construct this without having to specify the template arguments explicitly.
 * This algorithm uses the Bound-constrained
 * formulation of the augmented lagrangian method, roughly as described in Nocedal's Numerical Optimization 
 * book (i.e. similar to LANCELOT). This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Hess(f)(x), Jac(g)(x), and Jac(h)(x).\n
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \tparam EqFunction The functor type of the equality constraints function (vector function).
 * \tparam EqJacFunction The functor type of the equality constraints jacobian function.
 * \tparam IneqFunction The functor type of the inequality constraints function (vector function).
 * \tparam IneqJacFunction The functor type of the inequality constraints jacobian function.
 * \tparam TrustRegionSolver A functor type that can solve for a solution step within a trust-region (see trust_region_solver_dogleg for an example).
 * \tparam LimitFunction A functor type that can impose limits on a proposed solution step (see no_limit_functor or box_limit_function for examples).
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename T,
          typename EqFunction = no_constraint_functor, typename EqJacFunction = no_constraint_jac_functor, 
	  typename IneqFunction = no_constraint_functor, typename IneqJacFunction = no_constraint_jac_functor,
          typename TrustRegionSolver = trust_region_solver_dogleg, 
	  typename LimitFunction = no_limit_functor>
struct constraint_newton_method_tr_factory {
  Function f;
  GradFunction df;
  HessianFunction fill_hessian;
  T max_radius;
  unsigned int max_iter;
  T abs_tol;
  T abs_grad_tol;
  T eta;
  EqFunction g;
  EqJacFunction fill_g_jac;
  IneqFunction h;
  IneqJacFunction fill_h_jac;
  TrustRegionSolver solve_step;
  LimitFunction impose_limits;
  
  /**
   * Parametrized constructor of the factory object.
   * \param aF The function to minimize.
   * \param aDf The gradient of the function to minimize.
   * \param aFillHessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized. 
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aG The function that computes the equality constraints.
   * \param aFillGJac The function that computes the jacobian matrix of the equality constaints.
   * \param aH The function that computes the inequality constraints.
   * \param aFillHJac The function that computes the jacobian matrix of the inequality constraints.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aGradTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aSolveStep The functor that can solve for the step to take within the trust-region.
   * \param aImposeLimits The functor that can impose simple limits on the search domain (e.g. box-constraints or non-negativity).
   */
  constraint_newton_method_tr_factory(Function aF, GradFunction aDf,
				      HessianFunction aFillHessian, T aMaxRadius, unsigned int aMaxIter,
				      EqFunction aG = EqFunction(), EqJacFunction aFillGJac = EqJacFunction(),
				      IneqFunction aH = IneqFunction(), IneqJacFunction aFillHJac = IneqJacFunction(),
				      T aTol = T(1e-6), T aGradTol = T(1e-6),
				      T aEta = T(1e-4),
				      TrustRegionSolver aSolveStep = TrustRegionSolver(),
				      LimitFunction aImposeLimits = LimitFunction()) :
				      f(aF), df(aDf), fill_hessian(aFillHessian),
				      max_radius(aMaxRadius), max_iter(aMaxIter), 
				      abs_tol(aTol), abs_grad_tol(aGradTol), eta(aEta),
				      g(aG), fill_g_jac(aFillGJac), h(aH), fill_h_jac(aFillHJac),
				      solve_step(aSolveStep),
				      impose_limits(aImposeLimits) { };
  /**
   * This function finds the minimum of a function, given its derivative and Hessian, 
   * using a newton search direction and using a trust-region approach.
   * \tparam Vector The vector type of the independent variable for the function.
   * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
   */
  template <typename Vector>
  void operator()(Vector& x) const {
    detail::bcl_newton_method_tr_impl(
      f, df, fill_hessian, g, fill_g_jac, 
      h, fill_h_jac, x, max_radius, max_iter,
      solve_step,
      impose_limits,abs_tol,abs_grad_tol,eta);
  };
    
  /**
   * This function remaps the factory to one which will use a regularized solver within the trust-region.
   * You should regularize the matrix only if there are reasons to expect the Hessian to be near-singular.
   * \param tau The initial relative damping factor to regularize the Hessian matrix.
   */
  constraint_newton_method_tr_factory<Function,GradFunction,HessianFunction,T,
                                      EqFunction, EqJacFunction, IneqFunction, IneqJacFunction,
                                      trust_region_solver_dogleg_reg<T>, LimitFunction>
    regularize(const T& tau) const {
    return constraint_newton_method_tr_factory<Function,GradFunction,HessianFunction,T,
                                               EqFunction, EqJacFunction, IneqFunction, IneqJacFunction,
                                               trust_region_solver_dogleg_reg<T>, LimitFunction>(f,df,fill_hessian,
					                                                         max_radius, max_iter,
												 g, fill_g_jac, h, fill_h_jac,
												 abs_tol, abs_grad_tol, eta,
										                 trust_region_solver_dogleg_reg<T>(tau),
										                 impose_limits);
  };
    
  /**
   * This function remaps the factory to one which will use the given solver within the trust-region.
   * \tparam NewTrustRegionSolver A new functor type that can solve for a solution step within a trust-region (see trust_region_solver_dogleg for an example).
   * \param new_solver The functor that can solve for the step to take within the trust-region.
   */
  template <typename NewTrustRegionSolver>
  constraint_newton_method_tr_factory<Function,GradFunction,HessianFunction,T,
                                      EqFunction, EqJacFunction, IneqFunction, IneqJacFunction,
                                      NewTrustRegionSolver, LimitFunction>
    set_tr_solver(NewTrustRegionSolver new_solver) const {
    return constraint_newton_method_tr_factory<Function,GradFunction,HessianFunction,T,
                                               EqFunction, EqJacFunction, IneqFunction, IneqJacFunction,
                                               NewTrustRegionSolver, LimitFunction>(f,df,fill_hessian,
					                                            max_radius, max_iter,
									            g, fill_g_jac, h, fill_h_jac,
									            abs_tol, abs_grad_tol, eta,
									            new_solver,impose_limits);
  };
    
  /**
   * This function remaps the factory to one which will use the given limit-function for the search domain.
   * Using a limit-function boils down to a gradient projection method, and thus, for more complex 
   * constraints (such as non-linear ones or mix of equality-inequality constraints), please use a 
   * constraint optimization method instead (see augmented_lagrangian_methods.hpp for example).
   * \tparam NewLimitFunction A new functor type that can impose limits on a proposed solution step (see no_limit_functor or box_limit_function for examples).
   * \param new_limits The functor that can impose simple limits on the search domain (i.e. using this boils down to a gradient projection method, for more complex constraints please use a constraint optimization method instead).
   */
  template <typename NewLimitFunction>
  constraint_newton_method_tr_factory<Function,GradFunction,HessianFunction,T,
                                      EqFunction, EqJacFunction, IneqFunction, IneqJacFunction,
                                      TrustRegionSolver, NewLimitFunction>
    set_limiter(NewLimitFunction new_limits) const {
    return constraint_newton_method_tr_factory<Function,GradFunction,HessianFunction,T,
                                               EqFunction, EqJacFunction, IneqFunction, IneqJacFunction,
                                               TrustRegionSolver, NewLimitFunction>(f,df,fill_hessian,
					                                            max_radius, max_iter,
									            g, fill_g_jac, h, fill_h_jac,
									            abs_tol, abs_grad_tol, eta,
									            solve_step,new_limits);
  };
    
  /**
   * This function remaps the factory to one which will use the given equality constraints for 
   * the search domain. The equality constraints must be formulated at g(x) = 0.
   * \tparam NewEqFunction The functor type of the equality constraints function (vector function).
   * \tparam NewEqJacFunction The functor type of the equality constraints jacobian function.
   * \param new_g The function that computes the equality constraints.
   * \param new_fill_g_jac The function that computes the jacobian matrix of the equality constaints.
   */
  template <typename NewEqFunction, typename NewEqJacFunction>
  constraint_newton_method_tr_factory<Function,GradFunction,HessianFunction,T,
                                      NewEqFunction, NewEqJacFunction, IneqFunction, IneqJacFunction,
                                      TrustRegionSolver, LimitFunction>
    set_eq_constraints(NewEqFunction new_g, NewEqJacFunction new_fill_g_jac) const {
    return constraint_newton_method_tr_factory<Function,GradFunction,HessianFunction,T,
                                               NewEqFunction, NewEqJacFunction, IneqFunction, IneqJacFunction,
                                               TrustRegionSolver, LimitFunction>(f,df,fill_hessian,
					                                         max_radius, max_iter,
									         new_g, new_fill_g_jac, h, fill_h_jac,
									         abs_tol, abs_grad_tol, eta,
									         solve_step,impose_limits);
  };
    
  /**
   * This function remaps the factory to one which will use the given inequality constraints for 
   * the search domain. The inequality constraints must be formulated at h(x) >= 0.
   * \tparam NewIneqFunction The functor type of the equality constraints function (vector function).
   * \tparam NewIneqJacFunction The functor type of the equality constraints jacobian function.
   * \param new_h The function that computes the inequality constraints.
   * \param new_fill_h_jac The function that computes the jacobian matrix of the inequality constaints.
   */
  template <typename NewIneqFunction, typename NewIneqJacFunction>
  constraint_newton_method_tr_factory<Function,GradFunction,HessianFunction,T,
                                      EqFunction, EqJacFunction, NewIneqFunction, NewIneqJacFunction,
                                      TrustRegionSolver, LimitFunction>
    set_ineq_constraints(NewIneqFunction new_h, NewIneqJacFunction new_fill_h_jac) const {
    return constraint_newton_method_tr_factory<Function,GradFunction,HessianFunction,T,
                                               EqFunction, EqJacFunction, NewIneqFunction, NewIneqJacFunction,
                                               TrustRegionSolver, LimitFunction>(f,df,fill_hessian,
					                                         max_radius, max_iter,
									         g, fill_g_jac, new_h, new_fill_h_jac,
									         abs_tol, abs_grad_tol, eta,
									         solve_step,impose_limits);
  };
    
};

/**
 * This function template creates a factory object to construct a Constraint Newton-method optimizer 
 * routine that uses a trust-region approach. To set the constraints, manipulate the returned factory
 * object (e.g. set_eq_constraints() etc.). 
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized. 
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance on the norm of the step size.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 * \param eta The tolerance on the decrease in order to accept a step in the trust region.
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename T>
constraint_newton_method_tr_factory<Function,GradFunction,HessianFunction,T> 
  make_constraint_newton_method_tr(Function f, GradFunction df,
                                   HessianFunction fill_hessian, T max_radius, unsigned int max_iter,
			           T abs_tol = T(1e-6), T abs_grad_tol = T(1e-6),
			           T eta = T(1e-4)) {
  return constraint_newton_method_tr_factory<Function,GradFunction,HessianFunction,T>(f,df,fill_hessian,max_radius,max_iter,no_constraint_functor(),no_constraint_jac_functor(),no_constraint_functor(),no_constraint_jac_functor(),abs_tol,abs_grad_tol,eta);
};





/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a newton search 
 * direction and using a trust-region approach (dogleg solver). This algorithm uses the Bound-constrained
 * formulation of the augmented lagrangian method, roughly as described in Nocedal's Numerical Optimization 
 * book (i.e. similar to LANCELOT). This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Hess(f)(x), Jac(g)(x), and Jac(h)(x).\n
 * See algorithms eq_cnstr_newton_method_tr and ineq_cnstr_newton_method_tr for versions that take only an
 * equality or inequality constraint, respectively (no_constraint_functor and no_constraint_jac_functor can
 * also be used in-place of the contraint functions if none are desired).
 * \test Must create a unit-test for this.
 * 
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function (should be a resizable vector, modeling the ResizableVectorConcept).
 * \tparam EqFunction The functor type of the equality constraints function (vector function).
 * \tparam EqJacFunction The functor type of the equality constraints jacobian function.
 * \tparam IneqFunction The functor type of the inequality constraints function (vector function).
 * \tparam IneqJacFunction The functor type of the inequality constraints jacobian function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param g The function that computes the equality constraints.
 * \param fill_g_jac The function that computes the jacobian matrix of the equality constaints.
 * \param h The function that computes the inequality constraints.
 * \param fill_h_jac The function that computes the jacobian matrix of the inequality constraints.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance on the norm of the step size.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename Vector, 
          typename EqFunction, typename EqJacFunction, 
	  typename IneqFunction, typename IneqJacFunction>
void constraint_newton_method_tr(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, 
				 EqFunction g, EqJacFunction fill_g_jac,
				 IneqFunction h, IneqJacFunction fill_h_jac,
		                 typename vect_traits<Vector>::value_type max_radius, unsigned int max_iter = 100,
		                 typename vect_traits<Vector>::value_type abs_tol = typename vect_traits<Vector>::value_type(1e-6),
		                 typename vect_traits<Vector>::value_type abs_grad_tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  detail::bcl_newton_method_tr_impl(
    f, df, fill_hessian, g, fill_g_jac, 
    h, fill_h_jac, x, max_radius, max_iter,
    trust_region_solver_dogleg(),
    no_limit_functor(),abs_tol,abs_grad_tol
  );
  
};

/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a newton search 
 * direction and using a trust-region approach (dogleg solver). This algorithm uses the Bound-constrained
 * formulation of the augmented lagrangian method, roughly as described in Nocedal's Numerical Optimization 
 * book (i.e. similar to LANCELOT). This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 * \n
 *  given grad(f)(x), Hess(f)(x), and Jac(g)(x).\n
 * \test Must create a unit-test for this.
 * 
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function (should be a resizable vector, modeling the ResizableVectorConcept).
 * \tparam EqFunction The functor type of the equality constraints function (vector function).
 * \tparam EqJacFunction The functor type of the equality constraints jacobian function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param g The function that computes the equality constraints.
 * \param fill_g_jac The function that computes the jacobian matrix of the equality constaints.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance on the norm of the step size.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename Vector, 
          typename EqFunction, typename EqJacFunction>
void eq_cnstr_newton_method_tr(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, 
			       EqFunction g, EqJacFunction fill_g_jac,
			       typename vect_traits<Vector>::value_type max_radius, unsigned int max_iter = 100,
		               typename vect_traits<Vector>::value_type abs_tol = typename vect_traits<Vector>::value_type(1e-6),
		               typename vect_traits<Vector>::value_type abs_grad_tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  detail::bcl_newton_method_tr_impl(
    f, df, fill_hessian, g, fill_g_jac, 
    no_constraint_functor(), no_constraint_jac_functor(), x, max_radius, max_iter,
    trust_region_solver_dogleg(),
    no_limit_functor(),abs_tol,abs_grad_tol
  );
  
};

/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a newton search 
 * direction and using a trust-region approach (dogleg solver). This algorithm uses the Bound-constrained
 * formulation of the augmented lagrangian method, roughly as described in Nocedal's Numerical Optimization 
 * book (i.e. similar to LANCELOT). This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Hess(f)(x), and Jac(h)(x).\n
 * \test Must create a unit-test for this.
 * 
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function (should be a resizable vector, modeling the ResizableVectorConcept).
 * \tparam IneqFunction The functor type of the inequality constraints function (vector function).
 * \tparam IneqJacFunction The functor type of the inequality constraints jacobian function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param h The function that computes the inequality constraints.
 * \param fill_h_jac The function that computes the jacobian matrix of the inequality constraints.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance on the norm of the step size.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename Vector, 
          typename IneqFunction, typename IneqJacFunction>
void ineq_cnstr_newton_method_tr(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, 
				 IneqFunction h, IneqJacFunction fill_h_jac,
		                 typename vect_traits<Vector>::value_type max_radius, unsigned int max_iter = 100,
		                 typename vect_traits<Vector>::value_type abs_tol = typename vect_traits<Vector>::value_type(1e-6),
		                 typename vect_traits<Vector>::value_type abs_grad_tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  detail::bcl_newton_method_tr_impl(
    f, df, fill_hessian, 
    no_constraint_functor(), no_constraint_jac_functor(),
    h, fill_h_jac, x, max_radius, max_iter,
    trust_region_solver_dogleg(),
    no_limit_functor(),abs_tol,abs_grad_tol
  );
  
};



};

};


#endif


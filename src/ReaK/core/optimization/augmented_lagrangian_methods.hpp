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
			         Vector& x, typename vect_traits<Vector>::value_type max_radius, 
			         TrustRegionSolver solve_step, LimitFunction impose_limits, 
			         typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6), 
				 typename vect_traits<Vector>::value_type kappa = typename vect_traits<Vector>::value_type(1e-4)) {
    typedef typename vect_traits<Vector>::value_type ValueType;
    typedef typename vect_traits<Vector>::size_type SizeType;
    using std::sqrt; using std::fabs;
    
    SizeType N = x.size();
    Vector h_value = h(x);
    Vector s = -h_value;
    SizeType K = h_value.size();
    Vector g_value = g(x);
    SizeType M = g_value.size();
    
    mat<T, mat_structure::rectangular> Jac_g(M,N);
    fill_g_jac(Jac_g,x,g_value);
    mat<T, mat_structure::rectangular> Jac_h(K,N);
    fill_h_jac(Jac_h,x,h_value);
  
    if((M == 0) && (K == 0)) { //this means it is an unconstrained problem.
      newton_method_tr_impl(f,df,fill_hessian,x,max_radius,solve_step,impose_limits,tol);
      return;
    };
    
    ValueType mu(10.0);
    
    Vector l_g(M, ValueType(1.0)); //choose initial lambda vector to be 1.0.
    Vector l_h(K, ValueType(1.0));
  
    ValueType c_norm_star = (norm(g_value) + norm(h_value)) / ValueType(K + M);
    ValueType eta_star = tol * c_norm_star;
    ValueType norm_star = norm(x) / ValueType(N) + ValueType(1.0);
    ValueType omega_star = tol * norm_star;
  
    ValueType omega = norm_star / mu;
    ValueType eta = pow(mu,ValueType(-0.1)) * c_norm_star;
    
    ValueType radius = ValueType(0.5) * max_radius;
    ValueType x_value = f(x);
    Vector x_grad = df(x);
    
    ValueType L_value = x_value - l_g * g_value + ValueType(0.5) * mu * (g_value * g_value); 
                      // + 0.5 * mu * (h_value - s) * (h_value - s) - l_h * (h_value - s)
    Vector L_grad_s = l_h; // - mu * (h_value - s)
    Vector rhs_grad = x_grad - l_g * Jac_g + mu * (g_value * Jac_g);
    
    Vector p = rhs_grad;
    Vector p_s = s;
    ValueType norm_p = std::numeric_limits<ValueType>::max();
    mat<ValueType,mat_structure::symmetric> H(x.size());
    fill_hessian(H,x,x_value,x_grad);
    
    mat<ValueType,mat_structure::symmetric> JJ_g(transpose(Jac_g) * Jac_g);
    JJ_g *= mu;
    mat<ValueType,mat_structure::symmetric> L_H = H; L_H += JJ_g;
    
    Vector xt = x;
    Vector st = s;
    
    while (true) {
    
      while(norm_p > omega) {
        solve_step(rhs_grad,L_H,p,norm_p,radius,tol);
        p_s = (ValueType(-1.0) / mu) * L_grad_s + Jac_h * p;
	impose_limits(x,p);  //impose limits on the x-space
	for(SizeType i = 0; i < K; ++i) 
	  if(s[i] + p_s[i] < ValueType(0.0))
	    p_s[i] = -s[i];  //impose non-negativity on the s-space.
        xt = x; xt += p;
	norm_p = norm(p);
	st = s; st += p_s;
	g_value = g(xt);
	h_value = h(xt) - st;
	ValueType xt_value = f(xt);
	ValueType Lt_value = f(xt) - l_g * g_value - l_h * h_value + ValueType(0.5) * mu * (g_value * g_value + h_value * h_value);
        ValueType aredux = L_value - Lt_value;
        ValueType predux = -(rhs_grad * p + ValueType(0.5) * (p * (L_H * p)) + L_grad_s * p_s + ValueType(2.5) * (L_grad_s * L_grad_s));
        
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
        if( ratio > kappa ) {  //the step is accepted.
          x = xt;
          x_value = xt_value;
	  L_value = Lt_value;
          x_grad = df(x);
          fill_hessian(H,x,x_value,x_grad);
	  s = st;
	  fill_h_jac(Jac_h,x,h_value + s);
	  fill_g_jac(Jac_g,x,g_value);
	  JJ_g = transpose(Jac_g) * Jac_g; JJ_g *= mu;
	  L_H = H; L_H += JJ_g;
	  L_grad_s = l_h - mu * h_value;
	  rhs_grad = x_grad - l_g * Jac_g + mu * g_value * Jac_g;
        };
      };
      
      
      ValueType c_norm = sqrt(h_value * h_value + g_value * g_value);
      if( c_norm < eta ) {
      
        if((c_norm < eta_star) && (norm_p < omega_star)) //in theory, this test should involve the projected gradient, but I think that the step-size is sufficient.
	  return;
      
	l_g -= mu * g_value;
	l_h -= mu * h_value;
        eta *= pow(mu,ValueType(-0.9));
        omega /= mu;
      } else {
        mu *= ValueType(100.0);
        eta = pow(mu,ValueType(-0.1)) * c_norm_star;
        omega = norm_star / mu;
      };
      
    };
  };


  
  
  template <typename Function, typename Vector, 
            typename EqFunction, typename IneqFunction, 
	    typename T>
  T compute_lagrangian_impl(Function f, const Vector& x, 
                            EqFunction g, IneqFunction h,
			    const Vector& l, T mu,
			    unsigned int N, unsigned int K) {
    typedef typename vect_traits<Vector>::size_type SizeType;
    
    Vector x_r(N);
    for(SizeType i = 0; i < N; ++i)
      x_r[i] = x[i];
    T result = f(x_r);
    Vector c_g = g(x_r);
    SizeType M = c_g.size();
    Vector c_h = h(x_r);
    for(SizeType i = 0; i < K; ++i)
      c_h -= x[i + N];
    for(SizeType i = 0; i < M; ++i)
      result += T(0.5) * mu * c_g[i] * c_g[i] - l[i] * c_g[i];
    for(SizeType i = 0; i < K; ++i)
      result += T(0.5) * mu * c_h[i] * c_h[i] - l[i + M] * c_h[i];
    return result;
  };
  
  
  template <typename GradFunction, typename Vector, 
            typename EqFunction, typename EqJacFunction, 
	    typename IneqFunction, typename IneqJacFunction,
	    typename T>
  Vector compute_lagrangian_grad_impl(GradFunction df, const Vector& x,  
                                      EqFunction g, EqJacFunction fill_g_jac,
				      IneqFunction h, IneqJacFunction fill_h_jac,
			              const Vector& l, T mu,
			              unsigned int N, unsigned int K) {
    typedef typename vect_traits<Vector>::size_type SizeType;
    
    Vector x_r(N);
    for(SizeType i = 0; i < N; ++i)
      x_r[i] = x[i];
    Vector g_r = df(x_r);
    Vector c_g = g(x_r);
    SizeType M = c_g.size();
    Vector c_h = h(x_r);
    
    mat<T, mat_structure::rectangular> Jac_g(M,N);
    fill_g_jac(Jac_g,x_r,c_g);
    mat<T, mat_structure::rectangular> Jac_h(K,N);
    fill_h_jac(Jac_h,x_r,c_h);
    for(SizeType i = 0; i < K; ++i)
      c_h -= x[i + N];
    
    for(SizeType j = 0; j < N; ++j) {
      for(SizeType i = 0; i < M; ++i)
        g_r[j] += mu * c_g[i] * Jac_g(i,j) - l[i] * Jac_g(i,j);
      for(SizeType i = 0; i < K; ++i)
        g_r[j] += mu * c_h[i] * Jac_h(i,j) - l[i + M] * Jac_h(i,j);
    };
    
    Vector result(x.size());
    for(SizeType i = 0; i < N; ++i)
      result[i] = g_r[i];
    for(SizeType i = 0; i < K; ++j)
      result[i+N] = l[i + M] - mu * c_h[i];
    
    return result;
  };
  
  
  template <typename Function, typename GradFunction, typename HessianFunction, 
            typename EqFunction, typename EqJacFunction,
	    typename IneqFunction, typename IneqJacFunction, typename T>
  struct lagrangian_hess_filler_impl {
    Function f;
    GradFunction df;
    HessianFunction fill_hessian;
    EqFunction g;
    EqJacFunction fill_g_jac;
    IneqFunction h;
    IneqJacFunction fill_h_jac;
    T mu;
    unsigned int N, K;
    
    lagrangian_hess_filler_impl(Function aF, GradFunction aDf, HessianFunction aFill_hessian,
				EqFunction aG, EqJacFunction aFill_g_jac,
				IneqFunction aH, IneqJacFunction aFill_h_jac,
			        T aMu, unsigned int aN, unsigned int aK) :
			        f(aF), df(aDf), fill_hessian(aFill_hessian), 
			        g(aG), fill_g_jac(aFill_g_jac),
			        h(aH), fill_h_jac(aFill_h_jac),
			        mu(aMu), N(aN), K(aK) { };
     
    template <typename HessianMatrix, typename Vector>
    void operator()(HessianMatrix& H, const Vector& x, const T& x_value, const Vector& x_grad) const {
      typedef typename vect_traits<Vector>::size_type SizeType;
      if(H.get_row_count() != N + K)
        H.set_row_count(N+K);
      if(H.get_col_count() != N + K)
        H.set_col_count(N+K);
    
      Vector x_r(N);
      for(SizeType i = 0; i < N; ++i)
        x_r[i] = x[i];
      fill_hessian(sub(H)(range(0,N-1),range(0,N-1)),x_r,f(x),df(x));
    
      Vector c_g = g(x_r);
      SizeType M = c_g.size();
      Vector c_h = h(x_r);
    
      mat<T, mat_structure::rectangular> Jac_g(M,N);
      fill_g_jac(Jac_g,x_r,c_g);
      mat<T, mat_structure::rectangular> Jac_h(K,N);
      fill_h_jac(Jac_h,x_r,c_h);
    
      mat<T, mat_structure::rectangular> Jac_h_t = transpose(Jac_h);
      sub(H)(range(0,N-1),range(0,N-1)) += mu * mat<T,mat_structure::symmetric>(transpose(Jac_g) * Jac_g + Jac_h_t * Jac_h);
      sub(H)(range(0,N-1),range(N,N+K-1)) = mu * Jac_h_t;
      sub(H)(range(N,N+K-1),range(0,N-1)) = mu * Jac_h;
      sub(H)(range(N,N+K-1),range(N,N+K-1)) = mat<T,mat_structure::scalar>(K,mu);
      
    };
    
  };
  
  template <typename Function, typename GradFunction, typename HessianFunction, 
            typename EqFunction, typename EqJacFunction,
	    typename IneqFunction, typename IneqJacFunction, typename T>
  lagrangian_hess_filler_impl<Function, GradFunction, HessianFunction, 
                              EqFunction, EqJacFunction,
	                      IneqFunction,IneqJacFunction,T> 
    make_lagrangian_hess_filler_impl(
      Function aF, GradFunction aDf, HessianFunction aFill_hessian,
      EqFunction aG, EqJacFunction aFill_g_jac,
      IneqFunction aH, IneqJacFunction aFill_h_jac,
      T aMu, unsigned int aN, unsigned int aK) {
    return lagrangian_hess_filler_impl<Function, GradFunction, HessianFunction, 
                                       EqFunction, EqJacFunction,
	                               IneqFunction,IneqJacFunction,T>(
			                 aF, aDf, aFill_hessian, 
			                 aG, aFill_g_jac,
			                 aH, aFill_h_jac,
			                 aMu, aN, aK);
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
  T tol;
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
   * \param aG The function that computes the equality constraints.
   * \param aFillGJac The function that computes the jacobian matrix of the equality constaints.
   * \param aH The function that computes the inequality constraints.
   * \param aFillHJac The function that computes the jacobian matrix of the inequality constraints.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aSolveStep The functor that can solve for the step to take within the trust-region.
   * \param aImposeLimits The functor that can impose simple limits on the search domain (e.g. box-constraints or non-negativity).
   */
  constraint_newton_method_tr_factory(Function aF, GradFunction aDf,
				      HessianFunction aFillHessian, T aMaxRadius,
				      EqFunction aG = EqFunction(), EqJacFunction aFillGJac = EqJacFunction(),
				      IneqFunction aH = IneqFunction(), IneqJacFunction aFillHJac = IneqJacFunction(),
				      T aTol = T(1e-6),
				      T aEta = T(1e-4),
				      TrustRegionSolver aSolveStep = TrustRegionSolver(),
				      LimitFunction aImposeLimits = LimitFunction()) :
				      f(aF), df(aDf), fill_hessian(aFillHessian),
				      max_radius(aMaxRadius), tol(aTol), eta(aEta),
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
      h, fill_h_jac, x, max_radius,
      solve_step,
      impose_limits,tol,eta);
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
					                                                         max_radius, 
												 g, fill_g_jac, h, fill_h_jac,
												 tol, eta,
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
					                                            max_radius, 
									            g, fill_g_jac, h, fill_h_jac,
									            tol, eta,
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
					                                            max_radius, 
									            g, fill_g_jac, h, fill_h_jac,
									            tol, eta,
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
					                                         max_radius, 
									         new_g, new_fill_g_jac, h, fill_h_jac,
									         tol, eta,
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
					                                         max_radius, 
									         g, fill_g_jac, new_h, new_fill_h_jac,
									         tol, eta,
									         solve_step,impose_limits);
  };
    
};

/**
 * This function template creates a factory object to construct a Constraint Newton-method optimizer 
 * routine that uses a trust-region approach. To set the constraints, manipulate the returned factory
 * object (e.g. set_eq_constraints() etc.). 
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized. 
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 * \param eta The tolerance on the decrease in order to accept a step in the trust region.
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename T>
constraint_newton_method_tr_factory<Function,GradFunction,HessianFunction,T> 
  make_constraint_newton_method_tr(Function f, GradFunction df,
                                   HessianFunction fill_hessian, T max_radius,
			           T tol = T(1e-6),
			           T eta = T(1e-4)) {
  return constraint_newton_method_tr_factory<Function,GradFunction,HessianFunction,T>(f,df,fill_hessian,max_radius,no_constraint_functor(),no_constraint_jac_functor(),no_constraint_functor(),no_constraint_jac_functor(),tol,eta);
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
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename Vector, 
          typename EqFunction, typename EqJacFunction, 
	  typename IneqFunction, typename IneqJacFunction>
void constraint_newton_method_tr(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, 
				 EqFunction g, EqJacFunction fill_g_jac,
				 IneqFunction h, IneqJacFunction fill_h_jac,
		                 typename vect_traits<Vector>::value_type max_radius,
		                 typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  detail::bcl_newton_method_tr_impl(
    f, df, fill_hessian, g, fill_g_jac, 
    h, fill_h_jac, x, max_radius,
    trust_region_solver_dogleg(),
    no_limit_functor(),tol
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
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename Vector, 
          typename EqFunction, typename EqJacFunction>
void eq_cnstr_newton_method_tr(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, 
			       EqFunction g, EqJacFunction fill_g_jac,
			       typename vect_traits<Vector>::value_type max_radius,
		               typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  detail::bcl_newton_method_tr_impl(
    f, df, fill_hessian, g, fill_g_jac, 
    detail::no_constraint_functor(), detail::no_constraint_jac_functor(), x, max_radius,
    trust_region_solver_dogleg(),
    no_limit_functor(),tol
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
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename Vector, 
          typename IneqFunction, typename IneqJacFunction>
void ineq_cnstr_newton_method_tr(Function f, GradFunction df, HessianFunction fill_hessian, Vector& x, 
				 IneqFunction h, IneqJacFunction fill_h_jac,
		                 typename vect_traits<Vector>::value_type max_radius,
		                 typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  detail::bcl_newton_method_tr_impl(
    f, df, fill_hessian, 
    detail::no_constraint_functor(), detail::no_constraint_jac_functor(),
    h, fill_h_jac, x, max_radius,
    trust_region_solver_dogleg(),
    no_limit_functor(),tol
  );
  
};



};

};


#endif


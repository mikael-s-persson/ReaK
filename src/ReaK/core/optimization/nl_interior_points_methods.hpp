/**
 * \file nl_interior_points_methods.hpp
 *
 * The following library provides methods to perform constrained non-linear optimization using 
 * interior-points methods (or penalty-barrier methods). One option is a trust-region approach 
 * that first solves a quadratic program within a trust-region and then solves an 
 * equality-constrained quadratic program by the projected conjugate gradient method (this 
 * method is essentially an extension of the SQP method to also include inequality constraints).
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

#ifndef REAK_NL_INTERIOR_POINTS_METHODS_HPP
#define REAK_NL_INTERIOR_POINTS_METHODS_HPP

#include "base/defs.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_num_exceptions.hpp"
#include "lin_alg/mat_qr_decomp.hpp"

#include "newton_search_directions.hpp"
#include "limit_functions.hpp"
#include "line_search.hpp"
#include "trust_region_search.hpp"
#include "hessian_approx_update.hpp"

#include "newton_methods.hpp"
#include "quadratic_programs.hpp"
#include "augmented_lagrangian_methods.hpp"

namespace ReaK {
  
  
namespace optim {
  
  
namespace detail {
  
  
  template <typename Function, typename GradFunction, typename HessianFunction, 
            typename Vector, typename EqFunction, typename EqJacFunction, 
	    typename IneqFunction, typename IneqJacFunction, 
	    typename TrustRegionSolver, typename LimitFunction>
  void nl_intpoint_method_tr_impl(Function f, GradFunction df, HessianFunction fill_hessian,  
				  EqFunction g, EqJacFunction fill_g_jac,
				  IneqFunction h, IneqJacFunction fill_h_jac,
				  Vector& x, typename vect_traits<Vector>::value_type max_radius, 
				  typename vect_traits<Vector>::value_type mu, 
				  unsigned int max_iter,
			          TrustRegionSolver solve_step, LimitFunction impose_limits, 
			          typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6), 
				  typename vect_traits<Vector>::value_type kappa = typename vect_traits<Vector>::value_type(1e-4),
				  typename vect_traits<Vector>::value_type tau = typename vect_traits<Vector>::value_type(0.995),
				  typename vect_traits<Vector>::value_type rho = typename vect_traits<Vector>::value_type(1e-4)) {
    typedef typename vect_traits<Vector>::value_type ValueType;
    typedef typename vect_traits<Vector>::size_type SizeType;
    using std::sqrt; using std::fabs; using std::log;
    
    SizeType N = x.size();
    Vector gt_value = g(x);
    SizeType M = gt_value.size();
    Vector ht_value = h(x);
    SizeType K = ht_value.size();
    Vector c; c.resize(M+K);
    c[range(0,M-1)] = gt_value;
    c[range(M,M+K-1)] = ht_value;
    
    if((M == 0) && (K == 0)) { //this means it is an unconstrained problem. TODO change this to dispatch on the type of the fill-hessian functor.
      newton_method_tr_impl(f,df,fill_hessian,x,max_radius,solve_step,impose_limits,tol);
      return;
    };
    
    mat<ValueType, mat_structure::rectangular> Jac_g(M,N);
    fill_g_jac(Jac_g,x,gt_value);
    
    mat<ValueType, mat_structure::rectangular> Jac_h(K,N);
    fill_h_jac(Jac_h,x,ht_value);
    
    mat<ValueType,mat_structure::diagonal> mS(K);
    mat<ValueType,mat_structure::nil> Zero_mk(M,K);
    
    mat_vert_cat< 
      mat_ref_horiz_cat< mat<ValueType, mat_structure::rectangular>, mat<ValueType, mat_structure::nil> >,
      mat_ref_horiz_cat< mat<ValueType, mat_structure::rectangular>, mat<ValueType, mat_structure::diagonal> >
    > Jac_aug(
      mat_ref_horiz_cat< mat<ValueType, mat_structure::rectangular>, mat<ValueType, mat_structure::nil> >(
	Jac_g, Zero_mk
      ),
      mat_ref_horiz_cat< mat<ValueType, mat_structure::rectangular>, mat<ValueType, mat_structure::diagonal> >(
	Jac_h, mS
      )
    );
    
    //compute initial slack vector and roughly adjust x if needed.
    Vector s = ht_value;
    ValueType min_s(0.0);
    for(SizeType i = 0; i < K; ++i)
      if( s[i] < min_s )
	min_s = s[i];
    min_s *= ValueType(-1.5);
    
    if(min_s > ValueType(0.0)) {
      Vector e_s(s);
      for(SizeType i = 0; i < K; ++i)
	e_s[i] = min_s;
      s += e_s;
      x += e_s * Jac_h;
    };
    mS = mat<ValueType,mat_structure::diagonal>(-s);
    c[range(M,M+K-1)] -= s;
    
    ValueType gt_norm = norm_2(c[range(0,M-1)]);
    ValueType ht_norm = norm_2(c[range(M,M+K-1)]);
    
    ValueType radius = ValueType(0.5) * max_radius;
    ValueType x_value = f(x);
    Vector x_grad = df(x);
    mat<ValueType,mat_structure::symmetric> H(mat<ValueType,mat_structure::identity>(x.size()));
    fill_hessian(H,x,x_value,x_grad);
    
    
    Vector xt = x;
    Vector p_x = x;
    Vector st = s;
    Vector p_s = s;
    Vector v; v.resize(N+K);
    Vector p = v;
    Vector p_grad; p_grad.resize(N+K);
    p_grad[range(0,N-1)] = x_grad;
    p_grad[range(N,N+K-1)] = vect_scalar<ValueType>(K,-mu);
    ValueType norm_p = std::numeric_limits<ValueType>::max();
    ValueType norm_v = std::numeric_limits<ValueType>::max();
    Vector r = c;
    //ValueType rho(0.1);
    ValueType tol_mu = 1e-2 + tol;
    
    ValueType log_s(0.0);
    for(SizeType i = 0; i < K; ++i)
      log_s += log(s[i]);
    
    
    ValueType c_norm_star = norm_2(c);
    ValueType abs_c_tol = tol * (c_norm_star + norm_2(s) + norm_2(x));
    ValueType nu(std::numeric_limits<ValueType>::min());
    
    Vector yz; yz.resize(M+K);
    mat_vect_adaptor<Vector> yz_mat(yz);
    vect_ref_view< Vector > y(yz[range(0,M-1)]);
    vect_ref_view< Vector > z(yz[range(M,M+K-1)]);
    linlsq_QR(transpose_view(Jac_aug),yz_mat,mat_vect_adaptor<Vector>(p_grad), abs_c_tol);
    
    ValueType norm_star = norm_2(x_grad) + norm_2(y * Jac_g) + norm_2(z * Jac_h);
    ValueType abs_grad_tol = tol * norm_star;
    Vector l(x_grad - y * Jac_g - z * Jac_h);
    Vector lt = l;
    norm_star = norm_2(l);
    
    mat<ValueType,mat_structure::diagonal> SES(K);
    for(SizeType i = 0; i < K; ++i)
      SES(i,i) = z[i] * s[i];
    
    mat<ValueType,mat_structure::nil> Zero_nk(N,K);
    mat<ValueType,mat_structure::nil> Zero_kn(K,N);
    
    mat_vert_cat< 
      mat_const_ref_horiz_cat< mat<ValueType, mat_structure::symmetric>, mat<ValueType, mat_structure::nil> >,
      mat_const_ref_horiz_cat< mat<ValueType, mat_structure::nil>, mat<ValueType, mat_structure::diagonal> >
    > H_aug(
      mat_ref_horiz_cat< mat<ValueType, mat_structure::symmetric>, mat<ValueType, mat_structure::nil> >(
	H, Zero_nk
      ),
      mat_ref_horiz_cat< mat<ValueType, mat_structure::rectangular>, mat<ValueType, mat_structure::diagonal> >(
	Zero_kn, SES
      )
    );
    
    ValueType Err_value = norm_2( SES * vect_scalar<ValueType>(K,1.0) );
    if(Err_value < norm_star)
      Err_value = norm_star;
    if(Err_value < gt_norm)
      Err_value = gt_norm;
    if(Err_value < ht_norm)
      Err_value = ht_norm;
    
    unsigned int k = 0;
    
    while(Err_value > abs_grad_tol) {
      
      //compute error with mu.
      Err_value = norm_2( SES * vect_scalar<ValueType>(K,1.0) + p_grad[range(N,N+K-1)] );
      if(Err_value < norm_star)
        Err_value = norm_star;
      if(Err_value < gt_norm)
        Err_value = gt_norm;
      if(Err_value < ht_norm)
        Err_value = ht_norm;
      
      ValueType abs_grad_tol_mu = abs_grad_tol * (mu / tol);
      
      while((Err_value > abs_grad_tol_mu) && (++k <= max_iter)) {
        
        solve_step(c,Jac_aug,v,norm_v,ValueType(0.8) * radius, abs_grad_tol_mu);
	for(SizeType i = 0; i < K; ++i)
	  if( v[N + i] < -0.5 * tau )
	    v[N + i] = -0.5 * tau;
	r = c + Jac_aug * v;
	
        projected_CG_method(Jac_aug, r - c, H_aug, p_grad, p, abs_grad_tol_mu);
        norm_p = norm_2(p);
        if(norm_p > radius) {
	  p *= radius / norm_p;
	  norm_p = radius;
        };
	for(SizeType i = 0; i < N; ++i)
	  p_x[i] = s[i] * p[i];
        for(SizeType i = 0; i < K; ++i) {
	  if( p[i+N] < -tau )
	    p[i+N] = -tau;
	  p_s[i] = p[i+N];
	};
	impose_limits(x,p_x);
      
        ValueType pHp = p * (H_aug * p);
	ValueType dm_p = Jac_aug * p;
	ValueType dq_p = p_grad * p + pHp;
        ValueType nu_t;
        if(pHp > ValueType(0.0))
          nu_t = (p_grad * p + ValueType(0.5) * pHp) / ((ValueType(1.0) - rho) * norm_1(c));
        else
	  nu_t = (p_grad * p) / ((ValueType(1.0) - rho) * norm_1(c));
        if(nu < nu_t)
	  nu = nu_t + abs_c_tol;
	ValueType pred = -rho * nu * dm_p;
	
	
        xt = x; xt += p_x;
	st = s; st += p_s;
        ValueType xt_value = f(xt);
        gt_value = g(xt);
        ht_value = h(xt) - st;
	ValueType log_st(0.0);
        for(SizeType i = 0; i < K; ++i)
          log_st += log(st[i]);
	ValueType ct_norm_star = sqrt(norm_2_sqr(gt_value) + norm_2_sqr(ht_value));
	
        ValueType ared = x_value - xt_value - mu * (log_s - log_st) 
	                + nu * (c_norm_star - ct_norm_star);
	
      
        if(ared >= kappa * pred) {
	  //step is accepted. NOTE TODO A lot of shit to update here.
	  x += p_x;
	  x_value = xt_value;
	  x_grad = df(x);
	  p_grad[range(0,N-1)] = x_grad;
	  s += p_s;
	  mS = mat<ValueType,mat_structure::diagonal>(-s);
	  log_s = log_st;
          c[range(0,M-1)] = gt_value;   gt_norm = norm_2(gt_value);
          c[range(M,M+K-1)] = ht_value; ht_norm = norm_2(ht_value);
	  c_norm_star = ct_norm_star;
	  fill_g_jac(Jac_g,x,gt_value);      
	  fill_h_jac(Jac_h,x,ht_value + st);
	  if(norm_p > ValueType(0.8) * radius) {
            radius *= ValueType(2.0);
            if(radius > max_radius)
              radius = max_radius;
          };
          linlsq_QR(transpose_view(Jac_aug),yz_mat,mat_vect_adaptor<Vector>(p_grad), abs_c_tol);
	  lt = x_grad - y * Jac_g - z * Jac_h;
	  norm_star = norm_2(lt);
          fill_hessian(H,x,x_value,x_grad,p_x,lt - l);
	  l = lt;
	  for(SizeType i = 0; i < K; ++i)
            SES(i,i) = z[i] * s[i];
          //compute new error with mu.
	  Err_value = norm_2( SES * vect_scalar<ValueType>(K,1.0) + p_grad[range(N,N+K-1)] );
          if(Err_value < norm_star)
            Err_value = norm_star;
          if(Err_value < gt_norm)
            Err_value = gt_norm;
          if(Err_value < ht_norm)
            Err_value = ht_norm;
        } else {
          radius = ValueType(0.7) * norm_p;
        };
      };
      if(k > max_iter)
	throw maximum_iteration(max_iter);
      
      //decrease mu;
      ValueType sigma(2.0);
      ValueType zeta(1.0);
      ValueType sz_k = s * z / ValueType(K);
      for(SizeType i = 0; i < K; ++i)
	if( SES(i,i) < zeta * sz_k )
	  zeta = SES(i,i) / sz_k;
      if(zeta < 0.5)  //this is just for numerical stability (takes the denominator where it will equalize the quantities involved).
	sigma = (ValueType(1.0) - zeta) / (ValueType(20.0) * zeta);
      else
	sigma = ValueType(0.05) * ((ValueType(1.0) - zeta) / zeta);
      if(sigma > ValueType(2.0))
	sigma = ValueType(2.0);
      sigma *= sigma * sigma * ValueType(0.1);
      mu = sigma * sz_k;
      p_grad[range(N,N+K-1)] = vect_scalar<ValueType>(K,-mu);
      
      //compute error without mu.
      Err_value = norm_2( SES * vect_scalar<ValueType>(K,1.0) );
      if(Err_value < norm_star)
        Err_value = norm_star;
      if(Err_value < gt_norm)
        Err_value = gt_norm;
      if(Err_value < ht_norm)
        Err_value = ht_norm;
      
    };
  };
  
  
  
};





/**
 * This functor is a factory class to construct a Non-Linear Interior-Point optimizer routine 
 * that uses a trust-region approach and incorporates equality and inequality constraints via 
 * the sequential quadratic programming method and a Newton direction. Use make_nlip_newton_tr to
 * construct this without having to specify the template arguments explicitly.
 * This algorithm is roughly as described in Nocedal's Numerical Optimization 
 * book. This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Hess(f)(x), Jac(g)(x) and Jac(h)(x).\n
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
	  typename TrustRegionSolver = tr_solver_right_pinv_dogleg, 
	  typename LimitFunction = no_limit_functor>
struct nlip_newton_tr_factory {
  Function f;
  GradFunction df;
  HessianFunction fill_hessian;
  T max_radius;
  T mu;
  unsigned int max_iter;
  T tol;
  T eta;
  T tau;
  T rho;
  EqFunction g;
  EqJacFunction fill_g_jac;
  IneqFunction h;
  IneqJacFunction fill_h_jac;
  TrustRegionSolver solve_step;
  LimitFunction impose_limits;
  
  typedef nlip_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                                 EqFunction,EqJacFunction,IneqFunction,IneqJacFunction,
				 TrustRegionSolver,LimitFunction> self;
  
  /**
   * Parametrized constructor of the factory object.
   * \param aF The function to minimize.
   * \param aDf The gradient of the function to minimize.
   * \param aFillHessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized. 
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is positive and should start with a rather large value (relative to the scale of the function) and will be progressively decreased by the algorithm as it progresses).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aG The function that computes the equality constraints.
   * \param aFillGJac The function that computes the jacobian matrix of the equality constaints.
   * \param aG The function that computes the inequality constraints.
   * \param aFillGJac The function that computes the jacobian matrix of the inequality constaints.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint (barrier).
   * \param aRho The margin on the sufficient decrease of a step in the trust region.
   * \param aSolveStep The functor that can solve for the step to take within the trust-region.
   * \param aImposeLimits The functor that can impose simple limits on the search domain (e.g. box-constraints or non-negativity).
   */
  nlip_newton_tr_factory(Function aF, GradFunction aDf, HessianFunction aFillHessian, 
			 T aMaxRadius, T aMu, unsigned int aMaxIter,
			 EqFunction aG = EqFunction(), EqJacFunction aFillGJac = EqJacFunction(),
			 IneqFunction aH = EqFunction(), IneqJacFunction aFillHJac = IneqJacFunction(),
			 T aTol = T(1e-6), T aEta = T(1e-4), T aTau = T(0.995), T aRho = T(1e-4),
			 TrustRegionSolver aSolveStep = TrustRegionSolver(),
			 LimitFunction aImposeLimits = LimitFunction()) :
			 f(aF), df(aDf), fill_hessian(aFillHessian),
			 max_radius(aMaxRadius), mu(aMu), max_iter(aMaxIter), tol(aTol), eta(aEta), tau(aTau), rho(aRho),
			 g(aG), fill_g_jac(aFillGJac), 
			 h(aH), fill_h_jac(aFillHJac), 
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
    detail::nl_intpoint_method_tr_impl(
      f, df, hessian_update_dual_exact<HessianFunction>(fill_hessian), 
      g, fill_g_jac, h, fill_h_jac,
      x, max_radius, mu, max_iter, solve_step,
      impose_limits,tol,eta,tau,rho);
  };
  
  /**
   * Sets the maximum trust-region radius to use (i.e. maximum optimization step).
   */
  self& set_max_radius(T aMaxRadius) { max_radius = aMaxRadius; return *this; };
  /**
   * Sets the initial strength of the barrier on the inequalities (initial "barrier parameter"), 
   * this parameter is positive and should start with a rather large value (relative to the scale 
   * of the function and inequalities) and will be progressively decreased by the algorithm as it progresses).
   */
  self& set_initial_barrier_param(T aMu) { mu = aMu; return *this; };
  /**
   * Sets the initial damping of the Hessian matrix.
   */
  self& set_max_iteration(unsigned int aMaxIter) { max_iter = aMaxIter; return *this; };
  /**
   * Sets the relative tolerance on the norm of the step size.
   */
  self& set_tolerance(T aTol) { tol = aTol; return *this; };
  /**
   * Sets the tolerance on the decrease in order to accept a step in the trust region.
   */
  self& set_decrease_tolerance(T aEta) { eta = aEta; return *this; };
  /**
   * Sets the portion (close to 1.0) of a total step to do without coming too close to the inequality 
   * constraint (barrier).
   */
  self& set_barrier_step_margin(T aTau) { tau = aTau; return *this; };
  /**
   * Sets the margin on the sufficient decrease of a step in the trust region.
   */
  self& set_decrease_margin(T aRho) { rho = aRho; return *this; };
    
  /**
   * This function remaps the factory to one which will use a regularized solver within the trust-region.
   * You should regularize the matrix only if there are reasons to expect the Hessian to be near-singular.
   * \param aTau The initial relative damping factor to regularize the Hessian matrix.
   */
  nlip_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                         EqFunction, EqJacFunction, 
                         IneqFunction, IneqJacFunction, 
                         tr_solver_right_pinv_dogleg_reg<T>, LimitFunction>
    regularize(const T& aTau) const {
    return nlip_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                                  EqFunction, EqJacFunction, 
                                  IneqFunction, IneqJacFunction, 
                                  tr_solver_right_pinv_dogleg_reg<T>, LimitFunction>(f,df,fill_hessian,
					                                             max_radius, mu, max_iter,
										     g, fill_g_jac, 
										     h, fill_h_jac, 
										     tol, eta, tau, rho,
										     tr_solver_right_pinv_dogleg_reg<T>(aTau),
										     impose_limits);
  };
    
  /**
   * This function remaps the factory to one which will use the given solver within the trust-region.
   * \tparam NewTrustRegionSolver A new functor type that can solve for a solution step within a trust-region (see trust_region_solver_dogleg for an example).
   * \param new_solver The functor that can solve for the step to take within the trust-region.
   */
  template <typename NewTrustRegionSolver>
  nlip_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                         EqFunction, EqJacFunction, 
                         IneqFunction, IneqJacFunction, 
                         NewTrustRegionSolver, LimitFunction>
    set_tr_solver(NewTrustRegionSolver new_solver) const {
    return nlip_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                                  EqFunction, EqJacFunction, 
                                  IneqFunction, IneqJacFunction, 
                                  NewTrustRegionSolver, LimitFunction>(f,df,fill_hessian,
					                               max_radius, mu, max_iter,
								       g, fill_g_jac, 
								       h, fill_h_jac, 
								       tol, eta, tau, rho,
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
  nlip_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                         EqFunction, EqJacFunction,
                         IneqFunction, IneqJacFunction, 
                         TrustRegionSolver, NewLimitFunction>
    set_limiter(NewLimitFunction new_limits) const {
    return nlip_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                                  EqFunction, EqJacFunction,
                                  IneqFunction, IneqJacFunction, 
                                  TrustRegionSolver, NewLimitFunction>(f,df,fill_hessian,
					                               max_radius, mu, max_iter,
								       g, fill_g_jac, 
								       h, fill_h_jac, 
								       tol, eta, tau, rho,
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
  nlip_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                         NewEqFunction, NewEqJacFunction,
                         IneqFunction, IneqJacFunction, 
                         TrustRegionSolver, LimitFunction>
    set_eq_constraints(NewEqFunction new_g, NewEqJacFunction new_fill_g_jac) const {
    return nlip_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                                  NewEqFunction, NewEqJacFunction,
                                  IneqFunction, IneqJacFunction, 
                                  TrustRegionSolver, LimitFunction>(f,df,fill_hessian,
					                            max_radius, mu, max_iter,
								    new_g, new_fill_g_jac, 
								    h, fill_h_jac, 
								    tol, eta, tau, rho,
								    solve_step,impose_limits);
  };
    
  /**
   * This function remaps the factory to one which will use the given inequality constraints for 
   * the search domain. The inequality constraints must be formulated at h(x) >= 0.
   * \tparam NewIneqFunction The functor type of the inequality constraints function (vector function).
   * \tparam NewIneqJacFunction The functor type of the inequality constraints jacobian function.
   * \param new_h The function that computes the inequality constraints.
   * \param new_fill_h_jac The function that computes the jacobian matrix of the inequality constaints.
   */
  template <typename NewIneqFunction, typename NewIneqJacFunction>
  nlip_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                         EqFunction, EqJacFunction,
                         NewIneqFunction, NewIneqJacFunction, 
                         TrustRegionSolver, LimitFunction>
    set_ineq_constraints(NewIneqFunction new_h, NewIneqJacFunction new_fill_h_jac) const {
    return nlip_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                                  EqFunction, EqJacFunction,
                                  NewIneqFunction, NewIneqJacFunction, 
                                  TrustRegionSolver, LimitFunction>(f,df,fill_hessian,
					                            max_radius, mu, max_iter,
								    g, fill_g_jac, 
								    new_h, new_fill_h_jac, 
								    tol, eta, tau, rho,
								    solve_step,impose_limits);
  };
    
};

/**
 * This function template creates a factory object to construct a Non-Linear 
 * Interior-Point optimizer routine that uses a trust-region approach and 
 * incorporates equality and inequality constraints via the sequential quadratic 
 * programming method and a Newton steps. This algorithm is roughly as described 
 * in Nocedal's Numerical Optimization book. This algorithm solves the following 
 * problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Hess(f)(x), Jac(g)(x) and Jac(h)(x).\n
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized. 
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param mu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is positive and should start with a rather large value (relative to the scale of the function) and will be progressively decreased by the algorithm as it progresses).
 * \param max_iter The maximum number of iterations to perform.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 * \param eta The tolerance on the decrease in order to accept a step in the trust region.
 * \param tau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint (barrier).
 * \param rho The margin on the sufficient decrease of a step in the trust region.
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename T>
nlip_newton_tr_factory<Function,GradFunction,HessianFunction,T> 
  make_nlip_newton_tr(Function f, GradFunction df, HessianFunction fill_hessian, 
		       T max_radius, T mu, unsigned int max_iter,
		       T tol = T(1e-6), T eta = T(1e-4), T tau = T(0.995), T rho = T(1e-4)) {
  return nlip_newton_tr_factory<Function,GradFunction,HessianFunction,T>(f,df,fill_hessian,max_radius,mu,max_iter,no_constraint_functor(),no_constraint_jac_functor(),tol,eta,tau,rho);
};






/**
 * This functor is a factory class to construct a Non-Linear 
 * Interior-Point optimizer routine that uses a trust-region approach and 
 * incorporates equality and inequality constraints via the sequential quadratic 
 * programming method. This version uses a Quasi-Newton search method (by default the Hessian
 * approximation is obtained with a symmetric rank-one update). Use make_bosqp_newton_tr to
 * construct this without having to specify the template arguments explicitly.
 * This algorithm is roughly as described in Nocedal's Numerical Optimization 
 * book. This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Jac(g)(x) and Jac(h)(x).\n
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \tparam EqFunction The functor type of the equality constraints function (vector function).
 * \tparam EqJacFunction The functor type of the equality constraints jacobian function.
 * \tparam IneqFunction The functor type of the inequality constraints function (vector function).
 * \tparam IneqJacFunction The functor type of the inequality constraints jacobian function.
 * \tparam HessianUpdater The functor type to update the Hessian approximation of the function to optimize.
 * \tparam TrustRegionSolver A functor type that can solve for a solution step within a trust-region (see trust_region_solver_dogleg for an example).
 * \tparam LimitFunction A functor type that can impose limits on a proposed solution step (see no_limit_functor or box_limit_function for examples).
 */
template <typename Function, typename GradFunction, typename T,
          typename EqFunction = no_constraint_functor, typename EqJacFunction = no_constraint_jac_functor, 
	  typename IneqFunction = no_constraint_functor, typename IneqJacFunction = no_constraint_jac_functor, 
	  typename HessianUpdater = hessian_update_sr1, typename TrustRegionSolver = tr_solver_right_pinv_dogleg, 
	  typename LimitFunction = no_limit_functor>
struct nlip_quasi_newton_tr_factory {
  Function f;
  GradFunction df;
  T max_radius;
  T mu;
  unsigned int max_iter;
  T tol;
  T eta;
  T tau;
  T rho;
  EqFunction g;
  EqJacFunction fill_g_jac;
  IneqFunction h;
  IneqJacFunction fill_h_jac;
  HessianUpdater update_hessian;
  TrustRegionSolver solve_step;
  LimitFunction impose_limits;
  
  typedef nlip_quasi_newton_tr_factory<Function,GradFunction,T,
                                       EqFunction,EqJacFunction,
				       IneqFunction,IneqJacFunction,
				       HessianUpdater,TrustRegionSolver,
				       LimitFunction> self;
  
  /**
   * Parametrized constructor of the factory object.
   * \param aF The function to minimize.
   * \param aDf The gradient of the function to minimize.
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is positive and should start with a rather large value (relative to the scale of the function) and will be progressively decreased by the algorithm as it progresses).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aG The function that computes the equality constraints.
   * \param aFillGJac The function that computes the jacobian matrix of the equality constaints.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint (barrier).
   * \param aRho The margin on the sufficient decrease of a step in the trust region.
   * \param aUpdateHessian The functor object that can update the approximate Hessian matrix of the function to be optimized. 
   * \param aSolveStep The functor that can solve for the step to take within the trust-region.
   * \param aImposeLimits The functor that can impose simple limits on the search domain (e.g. box-constraints or non-negativity).
   */
  nlip_quasi_newton_tr_factory(Function aF, GradFunction aDf,
			       T aMaxRadius, T aMu, unsigned int aMaxIter,
			       EqFunction aG = EqFunction(), EqJacFunction aFillGJac = EqJacFunction(),
			       IneqFunction aH = IneqFunction(), IneqJacFunction aFillHJac = IneqJacFunction(),
			       T aTol = T(1e-6), T aEta = T(1e-4), T aTau = T(0.995), T aRho = T(1e-4),
			       HessianUpdater aUpdateHessian = HessianUpdater(),
			       TrustRegionSolver aSolveStep = TrustRegionSolver(),
			       LimitFunction aImposeLimits = LimitFunction()) :
			       f(aF), df(aDf), 
			       max_radius(aMaxRadius), mu(aMu), max_iter(aMaxIter), 
			       tol(aTol), eta(aEta), tau(aTau), rho(aRho),
			       g(aG), fill_g_jac(aFillGJac), 
			       h(aH), fill_h_jac(aFillHJac), 
			       update_hessian(aUpdateHessian),
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
    detail::nl_intpoint_method_tr_impl(
      f, df, hessian_update_dual_quasi<HessianUpdater>(update_hessian), g, fill_g_jac, h, fill_h_jac, 
      x, max_radius, mu, max_iter, solve_step, impose_limits,tol,eta,tau,rho);
  };
  
  /**
   * Sets the maximum trust-region radius to use (i.e. maximum optimization step).
   */
  self& set_max_radius(T aMaxRadius) { max_radius = aMaxRadius; return *this; };
  /**
   * Sets the initial strength of the barrier on the inequalities (initial "barrier parameter"), 
   * this parameter is positive and should start with a rather large value (relative to the scale 
   * of the function and inequalities) and will be progressively decreased by the algorithm as it progresses).
   */
  self& set_initial_barrier_param(T aMu) { mu = aMu; return *this; };
  /**
   * Sets the initial damping of the Hessian matrix.
   */
  self& set_max_iteration(unsigned int aMaxIter) { max_iter = aMaxIter; return *this; };
  /**
   * Sets the relative tolerance on the norm of the step size.
   */
  self& set_tolerance(T aTol) { tol = aTol; return *this; };
  /**
   * Sets the tolerance on the decrease in order to accept a step in the trust region.
   */
  self& set_decrease_tolerance(T aEta) { eta = aEta; return *this; };
  /**
   * Sets the portion (close to 1.0) of a total step to do without coming too close to the inequality 
   * constraint (barrier).
   */
  self& set_barrier_step_margin(T aTau) { tau = aTau; return *this; };
  /**
   * Sets the margin on the sufficient decrease of a step in the trust region.
   */
  self& set_decrease_margin(T aRho) { rho = aRho; return *this; };
    
  /**
   * This function remaps the factory to one which will use a regularized solver within the trust-region.
   * You should regularize the matrix only if there are reasons to expect the Hessian to be near-singular.
   * \param tau The initial relative damping factor to regularize the Hessian matrix.
   */
  nlip_quasi_newton_tr_factory<Function,GradFunction,T,
                               EqFunction,EqJacFunction,
			       IneqFunction, IneqJacFunction, HessianUpdater,
                               tr_solver_right_pinv_dogleg_reg<T>,LimitFunction>
    regularize(const T& tau) const {
    return nlip_quasi_newton_tr_factory<Function,GradFunction,T,
                                        EqFunction, EqJacFunction,IneqFunction,IneqJacFunction,HessianUpdater, 
                                        tr_solver_right_pinv_dogleg_reg<T>, LimitFunction>(f,df,
					                                                   max_radius, mu, max_iter,
										           g, fill_g_jac, 
									                   h, fill_h_jac,
										           tol, eta, tau, rho, update_hessian,
										           tr_solver_right_pinv_dogleg_reg<T>(tau),
										           impose_limits);
  };
    
  /**
   * This function remaps the factory to one which will use the given solver within the trust-region.
   * \tparam NewTrustRegionSolver A new functor type that can solve for a solution step within a trust-region (see trust_region_solver_dogleg for an example).
   * \param new_solver The functor that can solve for the step to take within the trust-region.
   */
  template <typename NewTrustRegionSolver>
  nlip_quasi_newton_tr_factory<Function,GradFunction,T,
                               EqFunction, EqJacFunction, 
			       IneqFunction, IneqJacFunction, 
			       HessianUpdater,
                               NewTrustRegionSolver, LimitFunction>
    set_tr_solver(NewTrustRegionSolver new_solver) const {
    return nlip_quasi_newton_tr_factory<Function,GradFunction,T,
                                        EqFunction, EqJacFunction,
			                IneqFunction, IneqJacFunction, 
			                HessianUpdater,
                                        NewTrustRegionSolver, LimitFunction>(f, df, max_radius, mu, max_iter,
								             g, fill_g_jac, 
									     h, fill_h_jac,
									     tol, eta, tau, rho, 
									     update_hessian, new_solver, impose_limits);
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
  nlip_quasi_newton_tr_factory<Function,GradFunction,T,
                               EqFunction,EqJacFunction,
			       IneqFunction, IneqJacFunction, 
			       HessianUpdater,
                               TrustRegionSolver, NewLimitFunction>
    set_limiter(NewLimitFunction new_limits) const {
    return nlip_quasi_newton_tr_factory<Function,GradFunction,T,
                                        EqFunction, EqJacFunction, 
					IneqFunction, IneqJacFunction, 
					HessianUpdater,
                                        TrustRegionSolver, NewLimitFunction>(f, df, max_radius, mu, max_iter,
								             g, fill_g_jac, 
									     h, fill_h_jac,
									     tol, eta, tau, rho, 
									     update_hessian, solve_step, new_limits);
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
  nlip_quasi_newton_tr_factory<Function,GradFunction,T,
                               NewEqFunction,NewEqJacFunction,
			       IneqFunction, IneqJacFunction, 
			       HessianUpdater,
                               TrustRegionSolver,LimitFunction>
    set_eq_constraints(NewEqFunction new_g, NewEqJacFunction new_fill_g_jac) const {
    return nlip_quasi_newton_tr_factory<Function,GradFunction,T,
                                        NewEqFunction, NewEqJacFunction,
					IneqFunction, IneqJacFunction, 
					HessianUpdater,
                                        TrustRegionSolver, LimitFunction>(f, df, max_radius, mu, max_iter,
									  new_g, new_fill_g_jac, 
									  h, fill_h_jac,
									  tol, eta, tau, rho, update_hessian,
									  solve_step,impose_limits);
  };
    
  /**
   * This function remaps the factory to one which will use the given inequality constraints for 
   * the search domain. The inequality constraints must be formulated at h(x) >= 0.
   * \tparam NewIneqFunction The functor type of the inequality constraints function (vector function).
   * \tparam NewIneqJacFunction The functor type of the inequality constraints jacobian function.
   * \param new_h The function that computes the inequality constraints.
   * \param new_fill_h_jac The function that computes the jacobian matrix of the inequality constaints.
   */
  template <typename NewIneqFunction, typename NewIneqJacFunction>
  nlip_quasi_newton_tr_factory<Function,GradFunction,T,
                               EqFunction,EqJacFunction,
			       NewIneqFunction, NewIneqJacFunction, 
			       HessianUpdater,
                               TrustRegionSolver,LimitFunction>
    set_ineq_constraints(NewIneqFunction new_h, NewIneqJacFunction new_fill_h_jac) const {
    return nlip_quasi_newton_tr_factory<Function,GradFunction,T,
                                        EqFunction, EqJacFunction,
					NewIneqFunction, NewIneqJacFunction, 
					HessianUpdater,
                                        TrustRegionSolver, LimitFunction>(f, df, max_radius, mu, max_iter,
									  g, fill_g_jac, 
									  new_h, new_fill_h_jac,
									  tol, eta, tau, rho, update_hessian,
									  solve_step,impose_limits);
  };
    
  /**
   * This function remaps the factory to one which will use the given solver within the trust-region.
   * \tparam NewHessianUpdater A new functor type to update the Hessian approximation of the function to optimize.
   * \param new_update_hessian The new functor object that can update the approximate Hessian matrix of the function to be optimized.
   */
  template <typename NewHessianUpdater>
  nlip_quasi_newton_tr_factory<Function,GradFunction,T,
                               EqFunction, EqJacFunction, 
			       IneqFunction, IneqJacFunction, 
			       NewHessianUpdater,
                               TrustRegionSolver, LimitFunction>
    set_hessian_updater(NewHessianUpdater new_update_hessian) const {
    return nlip_quasi_newton_tr_factory<Function,GradFunction,T,
                                        EqFunction, EqJacFunction, 
					IneqFunction, IneqJacFunction, 
					NewHessianUpdater,
                                        TrustRegionSolver, LimitFunction>(f, df, max_radius, mu, max_iter,
								          g, fill_g_jac, 
									  h, fill_h_jac,
									  tol, eta, tau, rho, 
									  new_update_hessian, 
									  solve_step, impose_limits);
  };
    
};

/**
 * This function template creates a factory object to construct a Non-Linear 
 * Interior-Point optimizer routine that uses a trust-region approach and 
 * incorporates equality and inequality constraints via the sequential quadratic 
 * programming method. This version uses a Quasi-Newton search method (by default the Hessian
 * approximation is obtained with a symmetric rank-one update). 
 * This algorithm is roughly as described in Nocedal's Numerical Optimization 
 * book. This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 *    h(x) >= 0 \n
 * \n
 *  given grad(f)(x), Jac(g)(x) and Jac(h)(x).\n
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param mu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is positive and should start with a rather large value (relative to the scale of the function) and will be progressively decreased by the algorithm as it progresses).
 * \param max_iter The maximum number of iterations to perform.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 * \param eta The tolerance on the decrease in order to accept a step in the trust region.
 * \param tau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint (barrier).
 * \param rho The margin on the sufficient decrease of a step in the trust region.
 */
template <typename Function, typename GradFunction, typename T>
nlip_quasi_newton_tr_factory<Function,GradFunction,T> 
  make_nlip_quasi_newton_tr(Function f, GradFunction df, 
		            T max_radius, T mu, unsigned int max_iter,
		            T tol = T(1e-6), T eta = T(1e-4), T tau = T(0.995), T rho = T(1e-4)) {
  return nlip_quasi_newton_tr_factory<Function,GradFunction,T>(f,df,max_radius,mu,max_iter,no_constraint_functor(),no_constraint_jac_functor(),no_constraint_functor(),no_constraint_jac_functor(),tol,eta,tau,rho);
};



};

};


#endif


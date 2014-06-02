/**
 * \file sequential_qp_methods.hpp
 *
 * The following library provides methods to perform constrained non-linear optimization using 
 * sequential quadratic programming methods (a trust-region approach that first solves a quadratic
 * program within a trust-region and then solves an equality-constrained quadratic program by the 
 * projected conjugate gradient method).
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

#ifndef REAK_SEQUENTIAL_QP_METHODS_HPP
#define REAK_SEQUENTIAL_QP_METHODS_HPP

#include <ReaK/core/base/defs.hpp>

#include <ReaK/core/lin_alg/mat_alg.hpp>
#include <ReaK/core/lin_alg/mat_num_exceptions.hpp>
#include <ReaK/core/lin_alg/mat_qr_decomp.hpp>

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
            typename TrustRegionSolver, typename LimitFunction>
  void byrd_omojokun_sqp_method_tr_impl(Function f, GradFunction df, HessianFunction fill_hessian,  
                                        EqFunction g, EqJacFunction fill_g_jac,
                                        Vector& x, typename vect_traits<Vector>::value_type max_radius,
                                        unsigned int max_iter,
                                        TrustRegionSolver solve_step, LimitFunction impose_limits, 
                                        typename vect_traits<Vector>::value_type abs_tol = typename vect_traits<Vector>::value_type(1e-6), 
                                        typename vect_traits<Vector>::value_type kappa = typename vect_traits<Vector>::value_type(1e-4), 
                                        typename vect_traits<Vector>::value_type rho = typename vect_traits<Vector>::value_type(1e-4)) {
    typedef typename vect_traits<Vector>::value_type ValueType;
    typedef typename vect_traits<Vector>::size_type SizeType;
    using std::sqrt; using std::fabs;
    
    SizeType N = x.size();
    Vector g_value = g(x);
    Vector gt_value = g_value;
    SizeType M = g_value.size();
    
    mat<ValueType, mat_structure::rectangular> Jac_g(M,N);
    fill_g_jac(Jac_g,x,g_value);
    mat_transpose_view< mat<ValueType, mat_structure::rectangular> > Jac_g_t = transpose_view(Jac_g);
    
    if(M == 0) { //this means it is an unconstrained problem.
      newton_method_tr_impl(f,df,fill_hessian,x,max_radius,max_iter,solve_step,impose_limits,abs_tol,abs_tol,kappa);
      return;
    };
    
    ValueType mu = -std::numeric_limits<ValueType>::infinity();
    
    ValueType radius = ValueType(0.5) * max_radius;
    ValueType x_value = f(x);
    Vector x_grad = df(x);
    
    Vector l(M, ValueType(1.0)); //choose initial lambda vector to be 1.0.
    mat_vect_adaptor<Vector> l_mat(l);
    linlsq_QR(Jac_g_t,l_mat,mat_vect_adaptor<Vector>(x_grad),abs_tol);
    
    ValueType c_norm_star = norm_2(g_value);
    Vector lag(x_grad - l * Jac_g);
    Vector lagt = lag;
    ValueType norm_star = norm_2(lag);;
    
    
    mat<ValueType,mat_structure::symmetric> H(mat<ValueType,mat_structure::identity>(x.size()));
    fill_hessian(H,x,x_value,x_grad);
    
    Vector xt = x;
    Vector p = x;
    Vector v = x;
    ValueType norm_p = std::numeric_limits<ValueType>::max();
    ValueType norm_v = std::numeric_limits<ValueType>::max();
    Vector r = g_value;
    
    unsigned int k = 0;
    
    while((norm_star > abs_tol) || (c_norm_star > abs_tol)) {
      
      solve_step(g_value,Jac_g,v,norm_v,ValueType(0.8) * radius,abs_tol);
      r = g_value; r += Jac_g * v;
      
      try {
        null_space_QP_method(Jac_g, r - g_value, H, x_grad, p, abs_tol, radius);
      } catch(singularity_error&) {
        try {
          projected_CG_method(Jac_g, r - g_value, H, x_grad, p, max_iter, abs_tol);
        } catch(maximum_iteration&) { };
      };
      norm_p = norm_2(p);
      if(norm_p > radius) {
        p *= radius / norm_p;
        norm_p = radius;
      };
      impose_limits(x,p);
      
      ValueType pHp = p * (H * p);
      ValueType mu_t;
      if(pHp > ValueType(0.0))
        mu_t = fabs(x_grad * p + ValueType(0.5) * pHp) / ((ValueType(1.0) - rho) * norm_1(g_value));
      else
        mu_t = fabs(x_grad * p) / ((ValueType(1.0) - rho) * norm_1(g_value));
      if(mu < mu_t)
        mu = mu_t + abs_tol;
      
      ValueType pred = mu * (c_norm_star - norm_2(g_value + Jac_g * p))
                       - x_grad * p - pHp;
      xt = x; xt += p;
      ValueType xt_value = f(xt);
      gt_value = g(xt);
      ValueType ared = x_value + mu * c_norm_star - xt_value - mu * norm_2(gt_value);
      
      //RK_NOTICE(1,"  x = " << x << " xt = " << xt << " x_value = " << x_value << " xt_value = " << xt_value);
      //RK_NOTICE(1,"  H = " << H << " Jac_g = " << Jac_g << " p = " << p << " pHp = " << pHp << " grad * p = " << x_grad * p);
      //RK_NOTICE(1,"  pred = " << pred << " ared = " << pred << " mu = " << mu << " c_norm = " << c_norm_star << " ct_norm = " << norm_2(gt_value));
      
      if(ared > kappa * pred) {
        //step is accepted.
        x += p;
        x_value = xt_value;
        x_grad = df(x);
        g_value = gt_value;
        c_norm_star = norm_2(g_value);
        fill_g_jac(Jac_g,x,g_value);
        if((ared > 0.8 * pred) && (norm_p > ValueType(0.8) * radius)) {
          radius *= ValueType(2.0);
          if(radius > max_radius)
            radius = max_radius;
        };
        linlsq_QR(Jac_g_t,l_mat,mat_vect_adaptor<Vector>(x_grad),abs_tol);
        lagt = x_grad - l * Jac_g;
        norm_star = norm_2(lagt);
        fill_hessian(H,x,x_value,x_grad,p,lagt - lag);
        lag = lagt;
      } else {
        radius = ValueType(0.9) * norm_p;
      };
      if(++k > max_iter)
        throw maximum_iteration(max_iter);
      
    };
  };
  
  
  
};




/**
 * This functor is a factory class to construct a Byrd-Omojokun SQP optimizer routine that uses a 
 * trust-region approach and incorporates equality constraints via the sequential quadratic
 * programming method. Use make_bosqp_newton_tr to
 * construct this without having to specify the template arguments explicitly.
 * This algorithm is roughly as described in Nocedal's Numerical Optimization 
 * book. This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 * \n
 *  given grad(f)(x), Hess(f)(x), and Jac(g)(x).\n
 * TEST PASSED, convergence is quite good, close to state-of-the-art methods in commercial packages, however, the nlip method almost always performs better.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \tparam EqFunction The functor type of the equality constraints function (vector function).
 * \tparam EqJacFunction The functor type of the equality constraints jacobian function.
 * \tparam TrustRegionSolver A functor type that can solve for a solution step within a trust-region (see trust_region_solver_dogleg for an example).
 * \tparam LimitFunction A functor type that can impose limits on a proposed solution step (see no_limit_functor or box_limit_function for examples).
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename T,
          typename EqFunction = no_constraint_functor, typename EqJacFunction = no_constraint_jac_functor, 
          typename TrustRegionSolver = tr_solver_right_pinv_dogleg, 
          typename LimitFunction = no_limit_functor>
struct bosqp_newton_tr_factory {
  Function f;
  GradFunction df;
  HessianFunction fill_hessian;
  T max_radius;
  unsigned int max_iter;
  T tol;
  T eta;
  T rho;
  EqFunction g;
  EqJacFunction fill_g_jac;
  TrustRegionSolver solve_step;
  LimitFunction impose_limits;
  
  typedef bosqp_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                                  EqFunction,EqJacFunction,
                                  TrustRegionSolver,LimitFunction> self;
  
  /**
   * Parametrized constructor of the factory object.
   * \param aF The function to minimize.
   * \param aDf The gradient of the function to minimize.
   * \param aFillHessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized. 
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aG The function that computes the equality constraints.
   * \param aFillGJac The function that computes the jacobian matrix of the equality constaints.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aRho The margin on the sufficient decrease of a step in the trust region.
   * \param aSolveStep The functor that can solve for the step to take within the trust-region.
   * \param aImposeLimits The functor that can impose simple limits on the search domain (e.g. box-constraints or non-negativity).
   */
  bosqp_newton_tr_factory(Function aF, GradFunction aDf,
                          HessianFunction aFillHessian, T aMaxRadius, unsigned int aMaxIter,
                          EqFunction aG = EqFunction(), EqJacFunction aFillGJac = EqJacFunction(),
                          T aTol = T(1e-6), T aEta = T(1e-4), T aRho = T(1e-4),
                          TrustRegionSolver aSolveStep = TrustRegionSolver(),
                          LimitFunction aImposeLimits = LimitFunction()) :
                          f(aF), df(aDf), fill_hessian(aFillHessian),
                          max_radius(aMaxRadius), max_iter(aMaxIter), tol(aTol), eta(aEta), rho(aRho),
                          g(aG), fill_g_jac(aFillGJac), 
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
    detail::byrd_omojokun_sqp_method_tr_impl(
      f, df, hessian_update_dual_exact<HessianFunction>(fill_hessian), g, fill_g_jac, 
      x, max_radius, max_iter, solve_step,
      impose_limits,tol,eta,rho);
  };
  
  /**
   * Sets the maximum trust-region radius to use (i.e. maximum optimization step).
   */
  self& set_max_radius(T aMaxRadius) { max_radius = aMaxRadius; return *this; };
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
   * Sets the margin on the sufficient decrease of a step in the trust region.
   */
  self& set_decrease_margin(T aRho) { rho = aRho; return *this; };
    
  /**
   * This function remaps the factory to one which will use a regularized solver within the trust-region.
   * You should regularize the matrix only if there are reasons to expect the Hessian to be near-singular.
   * \param tau The initial relative damping factor to regularize the Hessian matrix.
   */
  bosqp_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                          EqFunction, EqJacFunction, 
                          tr_solver_right_pinv_dogleg_reg<T>, LimitFunction>
    regularize(const T& tau) const {
    return bosqp_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                                   EqFunction, EqJacFunction, 
                                   tr_solver_right_pinv_dogleg_reg<T>, LimitFunction>(f,df,fill_hessian,
                                                                                      max_radius, max_iter,
                                                                                      g, fill_g_jac, 
                                                                                      tol, eta, rho,
                                                                                      tr_solver_right_pinv_dogleg_reg<T>(tau),
                                                                                      impose_limits);
  };
    
  /**
   * This function remaps the factory to one which will use the given solver within the trust-region.
   * \tparam NewTrustRegionSolver A new functor type that can solve for a solution step within a trust-region (see trust_region_solver_dogleg for an example).
   * \param new_solver The functor that can solve for the step to take within the trust-region.
   */
  template <typename NewTrustRegionSolver>
  bosqp_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                               EqFunction, EqJacFunction, 
                               NewTrustRegionSolver, LimitFunction>
    set_tr_solver(NewTrustRegionSolver new_solver) const {
    return bosqp_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                                   EqFunction, EqJacFunction, 
                                   NewTrustRegionSolver, LimitFunction>(f,df,fill_hessian,
                                                                        max_radius, max_iter,
                                                                        g, fill_g_jac, 
                                                                        tol, eta, rho,
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
  bosqp_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                               EqFunction, EqJacFunction,
                               TrustRegionSolver, NewLimitFunction>
    set_limiter(NewLimitFunction new_limits) const {
    return bosqp_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                                   EqFunction, EqJacFunction,
                                   TrustRegionSolver, NewLimitFunction>(f,df,fill_hessian,
                                                                        max_radius, max_iter,
                                                                        g, fill_g_jac, 
                                                                        tol, eta, rho,
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
  bosqp_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                          NewEqFunction, NewEqJacFunction,
                          TrustRegionSolver, LimitFunction>
    set_eq_constraints(NewEqFunction new_g, NewEqJacFunction new_fill_g_jac) const {
    return bosqp_newton_tr_factory<Function,GradFunction,HessianFunction,T,
                                   NewEqFunction, NewEqJacFunction,
                                   TrustRegionSolver, LimitFunction>(f,df,fill_hessian,
                                                                     max_radius, max_iter,
                                                                     new_g, new_fill_g_jac, 
                                                                     tol, eta, rho,
                                                                     solve_step,impose_limits);
  };
    
};

/**
 * This function template creates a factory object to construct a Byrd-Omojokun SQP optimizer 
 * routine that uses a trust-region approach and incorporates equality constraints via the 
 * sequential quadratic programming method. Use make_bosqp_newton_tr to
 * construct this without having to specify the template arguments explicitly.
 * This algorithm is roughly as described in Nocedal's Numerical Optimization 
 * book. This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 * \n
 *  given grad(f)(x), Hess(f)(x), and Jac(g)(x).\n
 * TEST PASSED, convergence is quite good, close to state-of-the-art methods in commercial packages, however, the nlip method almost always performs better.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized. 
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param max_iter The maximum number of iterations to perform.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 * \param eta The tolerance on the decrease in order to accept a step in the trust region.
 * \param rho The margin on the sufficient decrease of a step in the trust region.
 */
template <typename Function, typename GradFunction, typename HessianFunction, typename T>
bosqp_newton_tr_factory<Function,GradFunction,HessianFunction,T> 
  make_bosqp_newton_tr(Function f, GradFunction df, HessianFunction fill_hessian, 
                       T max_radius, unsigned int max_iter,
                       T tol = T(1e-6), T eta = T(1e-4), T rho = T(1e-4)) {
  return bosqp_newton_tr_factory<Function,GradFunction,HessianFunction,T>(f,df,fill_hessian,max_radius,max_iter,no_constraint_functor(),no_constraint_jac_functor(),tol,eta,rho);
};






/**
 * This functor is a factory class to construct a Byrd-Omojokun SQP optimizer routine that uses a 
 * trust-region approach and incorporates equality constraints via the sequential quadratic
 * programming method. This version uses a Quasi-Newton search method (by default the Hessian
 * approximation is obtained with a symmetric rank-one update). Use make_bosqp_newton_tr to
 * construct this without having to specify the template arguments explicitly.
 * This algorithm is roughly as described in Nocedal's Numerical Optimization 
 * book. This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 * \n
 * given grad(f)(x), and Jac(g)(x).\n
 * TEST PASSED, convergence is quite good, close to state-of-the-art methods in commercial packages, however, the nlip method almost always performs better.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianUpdater The functor type to update the Hessian approximation of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \tparam EqFunction The functor type of the equality constraints function (vector function).
 * \tparam EqJacFunction The functor type of the equality constraints jacobian function.
 * \tparam TrustRegionSolver A functor type that can solve for a solution step within a trust-region (see trust_region_solver_dogleg for an example).
 * \tparam LimitFunction A functor type that can impose limits on a proposed solution step (see no_limit_functor or box_limit_function for examples).
 */
template <typename Function, typename GradFunction, typename T,
          typename EqFunction = no_constraint_functor, typename EqJacFunction = no_constraint_jac_functor, 
          typename HessianUpdater = hessian_update_sr1, typename TrustRegionSolver = tr_solver_right_pinv_dogleg, 
          typename LimitFunction = no_limit_functor>
struct bosqp_quasi_newton_tr_factory {
  Function f;
  GradFunction df;
  T max_radius;
  unsigned int max_iter;
  T tol;
  T eta;
  T rho;
  EqFunction g;
  EqJacFunction fill_g_jac;
  HessianUpdater update_hessian;
  TrustRegionSolver solve_step;
  LimitFunction impose_limits;
  
  typedef bosqp_quasi_newton_tr_factory<Function,GradFunction,T,
                                        EqFunction,EqJacFunction,
                                        HessianUpdater,TrustRegionSolver,
                                        LimitFunction> self;
  
  /**
   * Parametrized constructor of the factory object.
   * \param aF The function to minimize.
   * \param aDf The gradient of the function to minimize.
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aG The function that computes the equality constraints.
   * \param aFillGJac The function that computes the jacobian matrix of the equality constaints.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aRho The margin on the sufficient decrease of a step in the trust region.
   * \param aUpdateHessian The functor object that can update the approximate Hessian matrix of the function to be optimized. 
   * \param aSolveStep The functor that can solve for the step to take within the trust-region.
   * \param aImposeLimits The functor that can impose simple limits on the search domain (e.g. box-constraints or non-negativity).
   */
  bosqp_quasi_newton_tr_factory(Function aF, GradFunction aDf,
                                T aMaxRadius, unsigned int aMaxIter,
                                EqFunction aG = EqFunction(), EqJacFunction aFillGJac = EqJacFunction(),
                                T aTol = T(1e-6), T aEta = T(1e-4), T aRho = T(1e-4),
                                HessianUpdater aUpdateHessian = HessianUpdater(),
                                TrustRegionSolver aSolveStep = TrustRegionSolver(),
                                LimitFunction aImposeLimits = LimitFunction()) :
                                f(aF), df(aDf), 
                                max_radius(aMaxRadius), max_iter(aMaxIter), tol(aTol), eta(aEta), rho(aRho),
                                g(aG), fill_g_jac(aFillGJac), 
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
    detail::byrd_omojokun_sqp_method_tr_impl(
      f, df, hessian_update_dual_quasi<HessianUpdater>(update_hessian), g, fill_g_jac, 
      x, max_radius, max_iter, solve_step, impose_limits,tol,eta,rho);
  };
  
  /**
   * Sets the maximum trust-region radius to use (i.e. maximum optimization step).
   */
  self& set_max_radius(T aMaxRadius) { max_radius = aMaxRadius; return *this; };
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
   * Sets the margin on the sufficient decrease of a step in the trust region.
   */
  self& set_decrease_margin(T aRho) { rho = aRho; return *this; };
    
  /**
   * This function remaps the factory to one which will use a regularized solver within the trust-region.
   * You should regularize the matrix only if there are reasons to expect the Hessian to be near-singular.
   * \param tau The initial relative damping factor to regularize the Hessian matrix.
   */
  bosqp_quasi_newton_tr_factory<Function,GradFunction,T,
                                EqFunction,EqJacFunction,HessianUpdater,
                                tr_solver_right_pinv_dogleg_reg<T>,LimitFunction>
    regularize(const T& tau) const {
    return bosqp_quasi_newton_tr_factory<Function,GradFunction,T,
                                         EqFunction, EqJacFunction,HessianUpdater, 
                                         tr_solver_right_pinv_dogleg_reg<T>, LimitFunction>(f,df,
                                                                                            max_radius, max_iter,
                                                                                            g, fill_g_jac, 
                                                                                            tol, eta, rho, update_hessian,
                                                                                            tr_solver_right_pinv_dogleg_reg<T>(tau),
                                                                                            impose_limits);
  };
    
  /**
   * This function remaps the factory to one which will use the given solver within the trust-region.
   * \tparam NewTrustRegionSolver A new functor type that can solve for a solution step within a trust-region (see trust_region_solver_dogleg for an example).
   * \param new_solver The functor that can solve for the step to take within the trust-region.
   */
  template <typename NewTrustRegionSolver>
  bosqp_quasi_newton_tr_factory<Function,GradFunction,T,
                                EqFunction, EqJacFunction, HessianUpdater,
                                NewTrustRegionSolver, LimitFunction>
    set_tr_solver(NewTrustRegionSolver new_solver) const {
    return bosqp_quasi_newton_tr_factory<Function,GradFunction,T,
                                         EqFunction, EqJacFunction, HessianUpdater,
                                         NewTrustRegionSolver, LimitFunction>(f, df, max_radius, max_iter,
                                                                              g, fill_g_jac, tol, eta, rho, 
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
  bosqp_quasi_newton_tr_factory<Function,GradFunction,T,
                                EqFunction,EqJacFunction,HessianUpdater,
                                TrustRegionSolver, NewLimitFunction>
    set_limiter(NewLimitFunction new_limits) const {
    return bosqp_quasi_newton_tr_factory<Function,GradFunction,T,
                                         EqFunction, EqJacFunction, HessianUpdater,
                                         TrustRegionSolver, NewLimitFunction>(f, df, max_radius, max_iter,
                                                                              g, fill_g_jac, tol, eta, rho, 
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
  bosqp_quasi_newton_tr_factory<Function,GradFunction,T,
                                NewEqFunction,NewEqJacFunction,HessianUpdater,
                                TrustRegionSolver,LimitFunction>
    set_eq_constraints(NewEqFunction new_g, NewEqJacFunction new_fill_g_jac) const {
    return bosqp_quasi_newton_tr_factory<Function,GradFunction,T,
                                         NewEqFunction, NewEqJacFunction,HessianUpdater,
                                         TrustRegionSolver, LimitFunction>(f, df, max_radius, max_iter,
                                                                           new_g, new_fill_g_jac, 
                                                                           tol, eta, rho, update_hessian,
                                                                           solve_step,impose_limits);
  };
    
  /**
   * This function remaps the factory to one which will use the given solver within the trust-region.
   * \tparam NewHessianUpdater A new functor type to update the Hessian approximation of the function to optimize.
   * \param new_update_hessian The new functor object that can update the approximate Hessian matrix of the function to be optimized.
   */
  template <typename NewHessianUpdater>
  bosqp_quasi_newton_tr_factory<Function,GradFunction,T,
                                EqFunction, EqJacFunction, NewHessianUpdater,
                                TrustRegionSolver, LimitFunction>
    set_hessian_updater(NewHessianUpdater new_update_hessian) const {
    return bosqp_quasi_newton_tr_factory<Function,GradFunction,T,
                                         EqFunction, EqJacFunction, NewHessianUpdater,
                                         TrustRegionSolver, LimitFunction>(f, df, max_radius, max_iter,
                                                                           g, fill_g_jac, tol, eta, rho, 
                                                                           new_update_hessian, 
                                                                           solve_step, impose_limits);
  };
    
};

/**
 * This function template creates a factory object to construct a Byrd-Omojokun SQP optimizer 
 * routine that uses a trust-region approach and incorporates equality constraints via the 
 * sequential quadratic programming method. This version uses a Quasi-Newton search method 
 * (by default the Hessian approximation is obtained with a symmetric rank-one update). 
 * This algorithm is roughly as described in Nocedal's Numerical Optimization 
 * book. This algorithm solves the following problem:
 * \n
 *   min f(x) \n
 *    g(x) = 0 \n
 * \n
 * given grad(f)(x), and Jac(g)(x).\n
 * TEST PASSED, convergence is quite good, close to state-of-the-art methods in commercial packages, however, the nlip method almost always performs better.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param max_iter The maximum number of iterations to perform.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 * \param eta The tolerance on the decrease in order to accept a step in the trust region.
 * \param rho The margin on the sufficient decrease of a step in the trust region.
 */
template <typename Function, typename GradFunction, typename T>
bosqp_quasi_newton_tr_factory<Function,GradFunction,T> 
  make_bosqp_quasi_newton_tr(Function f, GradFunction df, 
                             T max_radius, unsigned int max_iter,
                             T tol = T(1e-6), T eta = T(1e-4), T rho = T(1e-4)) {
  return bosqp_quasi_newton_tr_factory<Function,GradFunction,T>(f,df,max_radius,max_iter,no_constraint_functor(),no_constraint_jac_functor(),tol,eta,rho);
};




};

};


#endif


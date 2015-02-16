/**
 * \file levenberg_marquardt_method.hpp
 *
 * The following library is a collection of Levenberg-Marquardt methods.
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

#ifndef REAK_LEVENBERG_MARQUARDT_METHOD_HPP
#define REAK_LEVENBERG_MARQUARDT_METHOD_HPP

#include <ReaK/core/base/defs.hpp>

#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/math/lin_alg/mat_norms.hpp>
#include <ReaK/math/lin_alg/mat_cholesky.hpp>
#include <ReaK/math/lin_alg/mat_damped_matrix.hpp>

#include "optim_exceptions.hpp"
#include "limit_functions.hpp"

namespace ReaK {
  
  
namespace optim {


namespace detail { 


template <typename Function, typename JacobianFunction, 
          typename InputVector, typename OutputVector, 
          typename LinearSolver, typename LimitFunction, typename T>
int levenberg_marquardt_nllsq_impl(Function f, JacobianFunction fill_jac, 
                                    InputVector& x, const OutputVector& y, 
                                    LinearSolver lin_solve, LimitFunction impose_limits, unsigned int max_iter, 
                                    T tau, T epsj, T epsx, T epsy)
{
  typedef typename vect_traits<InputVector>::value_type ValueType;
  typedef typename vect_traits<InputVector>::size_type SizeType;
  
  /* Check if the problem is defined properly */
  if (y.size() < x.size())
    throw improper_problem("Levenberg-Marquardt requires M > N!");
  
  mat<ValueType,mat_structure::rectangular> J(y.size(),x.size());
  mat<ValueType,mat_structure::square> JtJ(x.size());
  mat<ValueType,mat_structure::diagonal> diag_JtJ(x.size());
  mat<ValueType,mat_structure::scalar> mu(x.size(),0.0);
  InputVector Jte = x;
  InputVector Dp = x; Dp -= x;
  mat_vect_adaptor<InputVector> Dp_mat(Dp);
  InputVector pDp = x;
  impose_limits(x,Dp); // make sure the initial solution is feasible.
  x += Dp;
  
  if(tau <= 0.0) 
    tau = 1E-03;
  if(epsj <= 0.0) 
    epsj = 1E-17;
  if(epsx <= 0.0) 
    epsx = 1E-17;
  ValueType epsx_sq = epsx * epsx;
  if(epsy <= 0.0) 
    epsy = 1E-17;
  if(max_iter <= 1)
    max_iter = 2;

  /* compute e=x - f(p) and its L2 norm */
  OutputVector y_approx = f(x);
  OutputVector e = y; e -= y_approx;
  OutputVector e_tmp = e;
  ValueType p_eL2 = e * e;
  
  unsigned int nu = 2;
  for(unsigned int k = 0; k < max_iter; ++k) {
    
    if(p_eL2 < epsy)
      return 1;  //residual is too small.

    fill_jac(J,x,y_approx);

    /* J^T J, J^T e */
    for(SizeType i = 0; i < J.get_col_count(); ++i) {
      for(SizeType j = i; j < J.get_col_count(); ++j) {
        ValueType tmp(0.0);
        for(SizeType l = 0; l < J.get_row_count(); ++l)
          tmp += J(l,i) * J(l,j);
        JtJ(i,j) = JtJ(j,i) = tmp;
      };
    };
    Jte = e * J;
    
    ValueType p_L2 = x * x;

    /* check for convergence */
    if( norm_inf(mat_vect_adaptor<InputVector>(Jte)) < epsj) 
      return 2;  //Jacobian is too small.

    /* compute initial damping factor */
    if( k == 0 ) {
      ValueType tmp = std::numeric_limits<ValueType>::min();
      for(SizeType i=0; i < JtJ.get_row_count(); ++i)
        if(JtJ(i,i) > tmp) 
          tmp = JtJ(i,i); /* find max diagonal element */
      mu = mat<ValueType,mat_structure::scalar>(x.size(), tau * tmp);
    };

    /* determine increment using adaptive damping */
    while(true) {

      /* solve augmented equations */
      try {
        lin_solve(make_damped_matrix(JtJ,mu),Dp_mat,mat_vect_adaptor<InputVector>(Jte),epsj);
        
        impose_limits(x,Dp);
        ValueType Dp_L2 = Dp * Dp;
        pDp = x; pDp += Dp;

        if(Dp_L2 < epsx_sq * p_L2) /* relative change in p is small, stop */
          return 3;  //steps are too small.

        if( Dp_L2 >= (p_L2 + epsx) / ( std::numeric_limits<ValueType>::epsilon() * std::numeric_limits<ValueType>::epsilon() ) ) 
          throw 42; //signal to throw a singularity-error (see below).
        
        e_tmp = y; e_tmp -= f(pDp);
        ValueType pDp_eL2 = e_tmp * e_tmp;
        ValueType dL = mu(0,0) * Dp_L2 + Dp * Jte;

        ValueType dF = p_eL2 - pDp_eL2;
        
        if( (dL < 0.0) || (dF < 0.0) ) 
          throw singularity_error("reject inc.");

        // reduction in error, increment is accepted
        ValueType tmp = ( ValueType(2.0) * dF / dL - ValueType(1.0));
        tmp = 1.0 - tmp * tmp * tmp;
        mu *= ( ( tmp >= ValueType(1.0 / 3.0) ) ? tmp : ValueType(1.0 / 3.0) );
        nu = 2;

        x = pDp;
        y_approx = y; y_approx -= e_tmp;
        e = e_tmp;
        p_eL2 = pDp_eL2;
        break; //the step is accepted and the loop is broken.
      } catch(singularity_error&) {
        //the increment must be rejected (either by singularity in damped matrix or no-redux by the step.
        mu *= ValueType(nu);
        nu <<= 1; // 2*nu;
        if( nu == 0 ) /* nu has overflown. */
          throw infeasible_problem("Levenberg-Marquardt method cannot reduce the function further, matrix damping has overflown!");
      } catch(int i) {
        if(i == 42)
          throw singularity_error("Levenberg-Marquardt method has detected a near-singularity in the Jacobian matrix!");
        else
          throw i;  //just in case there might be another integer thrown (very unlikely).
      };

    }; /* inner loop */
  };

  //if this point is reached, it means we have reached the maximum iterations.
  throw maximum_iteration(max_iter);

};


};




/**
 * This functor is a factory class to construct a non-linear least-square optimizer for a vector 
 * function. This algorithm uses the Levenberg-Marquardt method (which is a kind of adaptive 
 * damped-least square method, or trust-region Gauss-Newton). This method performs  
 * very well for most non-linear systems (i.e. a Vanilla solution) and can cope with 
 * some level of rank-deficiency (but not recommended for highly rank-deficient 
 * jacobians), and for very well-conditioned problems, the basic Gauss-Newton method 
 * is often a better choice (see gauss_newton_nllsq).
 * TEST PASSED
 * \tparam Function The functor type of the function to optimize.
 * \tparam JacobianFunction The functor type of the gradient of the function to optimize.
 * \tparam OutputVector The output vector type of the non-linear function provided.
 * \tparam LimitFunction A functor type that can impose limits on a proposed solution step (see no_limit_functor or box_limit_function for examples).
 * \tparam LinearSolver A functor type that can solve a linear system (see Cholesky_linsolver or QR_linlsqsolver for examples).
 * \tparam T The value-type of the field on which the optimization is performed.
 */
template <typename Function, typename JacobianFunction, 
          typename OutputVector, typename T = double,
          typename LimitFunction = no_limit_functor, 
          typename LinearSolver = Cholesky_linsolver>
struct levenberg_marquardt_nllsq_factory {
  Function f;
  JacobianFunction fill_jac;
  OutputVector y;
  unsigned int max_iter;
  T tau;
  T epsj;
  T epsx;
  T epsy;
  LimitFunction impose_limits;
  LinearSolver lin_solve;
  
  typedef levenberg_marquardt_nllsq_factory<
            Function,JacobianFunction,OutputVector,T,
            LimitFunction,LinearSolver> self;
  
  /**
   * Parametrized constructor of the factory object.
   * \param aF The vector-function to fit.
   * \param aFillJac The jacobian of the vector-function to fit.
   * \param aY The desired output vector of the function (to be matched by the algorithm).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aTau The initial, relative damping value to damp the Hessian matrix.
   * \param aEpsj The tolerance on the norm of the gradient.
   * \param aEpsx The tolerance on the norm of the step size.
   * \param aEpsy The tolerance on the norm of the residual.
   * \param aImposeLimits The functor that can impose simple limits on the search domain (i.e. using this boils down to a gradient projection method, for more complex constraints please use a constraint optimization method instead).
   * \param aLinSolver The functor that can solve a linear system.
   */
  levenberg_marquardt_nllsq_factory(Function aF, 
                                    JacobianFunction aFillJac,
                                    OutputVector aY, unsigned int aMaxIter = 100,
                                    T aTau = T(1e-3), T aEpsj = T(1e-6),
                                    T aEpsx = T(1e-6), T aEpsy = T(1e-6),
                                    LimitFunction aImposeLimits = LimitFunction(),
                                    LinearSolver aLinSolver = LinearSolver()) :
                                    f(aF), fill_jac(aFillJac), y(aY), max_iter(aMaxIter),
                                    tau(aTau), epsj(aEpsj), epsx(aEpsx), epsy(aEpsy), 
                                    impose_limits(aImposeLimits),
                                    lin_solve(aLinSolver) { };
  /**
   * This function finds the minimum of a function, given its derivative and Hessian, 
   * using a newton search direction and using a trust-region approach.
   * \tparam Vector The vector type of the independent variable for the function.
   * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
   * \return 1: Stopped by epsy (negl. residual), 2: Stopped by epsj (negl. gradient), 3: Stopped by epsx (negl. steps).
   * \throw improper_problem If the data of the problem is not proper (e.g. M < N).
   * \throw infeasible_problem If it becomes impossible to further damp the Jacobian (cannot make it non-singular), try to re-run with a larger "tau" value.
   * \throw singularity_error If a step is detected to cause a near-infinite jump.
   * \throw maximum_iteration If the maximum number of iterations has been reached.
   * \throw other Exceptions can originate from the functors.
   */
  template <typename Vector>
  int operator()(Vector& x) const {
    return detail::levenberg_marquardt_nllsq_impl(
             f,fill_jac,x,y,lin_solve,impose_limits,
             max_iter,tau,epsj,epsx,epsy);
  };
  
  /**
   * Sets the initial damping of the Hessian matrix.
   */
  self& set_max_iteration(unsigned int aMaxIter) { max_iter = aMaxIter; return *this; };
  /**
   * Sets the initial damping of the Hessian matrix.
   */
  self& set_initial_damping(T aTau) { tau = aTau; return *this; };
  /**
   * Sets the relative tolerance on the norm of the step size.
   */
  self& set_eps_j(T aEpsj) { epsj = aEpsj; return *this; };
  /**
   * Sets the relative tolerance on the norm of the step size.
   */
  self& set_eps_x(T aEpsx) { epsx = aEpsx; return *this; };
  /**
   * Sets the relative tolerance on the norm of the residual.
   */
  self& set_eps_y(T aEpsy) { epsy = aEpsy; return *this; };
  
  /**
   * This function remaps the factory to one which will use the given linear system solver.
   * \tparam NewLinearSolver A functor type that can solve a linear system (see Cholesky_linsolver or QR_linlsqsolver for examples).
   * \param new_lin_solver The functor that can solve a linear system.
   */
  template <typename NewLinearSolver>
  levenberg_marquardt_nllsq_factory<
            Function,JacobianFunction,OutputVector,T, 
            LimitFunction,NewLinearSolver>
    set_lin_solver(NewLinearSolver new_lin_solver) const {
    return levenberg_marquardt_nllsq_factory<
             Function,JacobianFunction,OutputVector,T,
             LimitFunction,NewLinearSolver>(f,fill_jac,y,max_iter,tau,epsj,epsx,epsy,
                                            impose_limits, new_lin_solver);
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
  levenberg_marquardt_nllsq_factory<
            Function,JacobianFunction,OutputVector, 
            NewLimitFunction,LinearSolver,T>
    set_limiter(NewLimitFunction new_limits) const {
    return levenberg_marquardt_nllsq_factory<
             Function,JacobianFunction,OutputVector,T, 
             NewLimitFunction,LinearSolver>(f,fill_jac,y,max_iter,tau,epsj,epsx,epsy,
                                            new_limits, lin_solve);
  };
    
};

/**
 * This function template creates a factory object to construct a Newton-method optimizer routine 
 * that uses a line-search approach. 
 * TEST PASSED
 * \tparam Function The functor type.
 * \tparam JacobianFunction The functor type to fill the Jacobian matrix.
 * \tparam OutputVector The vector type of the dependent variable of the function.
 * \param f The function to match to the fixed output vector.
 * \param fill_jac The Jacobian of the function to minimize.
 * \param y The fixed vector of outputs to be matched in a least-square approximation.
 * \param max_iter The maximum number of iterations to perform.
 * \param tau The initial, relative damping value to damp the Hessian matrix.
 * \param epsj The tolerance on the norm of the jacobian.
 * \param epsx The tolerance on the norm of the estimation steps.
 * \param epsy The tolerance on the norm of the residual vector.
 */
template <typename Function, typename JacobianFunction, typename OutputVector>
levenberg_marquardt_nllsq_factory<Function,JacobianFunction,OutputVector,typename vect_traits<OutputVector>::value_type> 
  make_levenberg_marquardt_nllsq(Function f, JacobianFunction fill_jac, OutputVector y, unsigned int max_iter = 100, 
                                 typename vect_traits<OutputVector>::value_type tau = typename vect_traits<OutputVector>::value_type(1E-03), 
                                 typename vect_traits<OutputVector>::value_type epsj = typename vect_traits<OutputVector>::value_type(1e-6), 
                                 typename vect_traits<OutputVector>::value_type epsx = typename vect_traits<OutputVector>::value_type(1e-6), 
                                 typename vect_traits<OutputVector>::value_type epsy = typename vect_traits<OutputVector>::value_type(1e-6)) {
  return levenberg_marquardt_nllsq_factory<Function,JacobianFunction,OutputVector,typename vect_traits<OutputVector>::value_type>(
           f,fill_jac,y,max_iter,tau,epsj,epsx,epsy);
};



/**
 * This function finds the non-linear least-square solution to a vector function.
 * This function uses the Levenberg-Marquardt method (which is a kind of adaptive 
 * damped-least square method, or trust-region Gauss-Newton). This method performs  
 * very well for most non-linear systems (i.e. a Vanilla solution) and can cope with 
 * some level of rank-deficiency (but not recommended for highly rank-deficient 
 * jacobians), and for very well-conditioned problems, the basic Gauss-Newton method 
 * is often a better choice (see gauss_newton_nllsq).
 * TEST PASSED
 * \tparam Function The functor type.
 * \tparam JacobianFunction The functor type to fill the Jacobian matrix.
 * \tparam InputVector The vector type of the independent variable of the function.
 * \tparam OutputVector The vector type of the dependent variable of the function.
 * \param f The function to match to the fixed output vector.
 * \param fill_jac The Jacobian of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param y The fixed vector of outputs to be matched in a least-square approximation.
 * \param max_iter The maximum number of iterations to perform.
 * \param tau The initial, relative damping value to damp the Hessian matrix.
 * \param epsj The tolerance on the norm of the jacobian.
 * \param epsx The tolerance on the norm of the estimation steps.
 * \param epsy The tolerance on the norm of the residual vector.
 * \return 1: Stopped by epsy (negl. residual), 2: Stopped by epsj (negl. gradient), 3: Stopped by epsx (negl. steps).
 * 
 * \throw improper_problem If the data of the problem is not proper (e.g. M < N).
 * \throw infeasible_problem If it becomes impossible to further damp the Jacobian (cannot make it non-singular), try to re-run with a larger "tau" value.
 * \throw singularity_error If a step is detected to cause a near-infinite jump.
 * \throw maximum_iteration If the maximum number of iterations has been reached.
 */
template <typename Function, typename JacobianFunction, 
          typename InputVector, typename OutputVector>
int levenberg_marquardt_nllsq(Function f, JacobianFunction fill_jac, 
                              InputVector& x, const OutputVector& y, 
                              unsigned int max_iter = 100, 
                              typename vect_traits<OutputVector>::value_type tau = typename vect_traits<OutputVector>::value_type(1E-03), 
                              typename vect_traits<OutputVector>::value_type epsj = typename vect_traits<OutputVector>::value_type(1e-12), 
                              typename vect_traits<OutputVector>::value_type epsx = typename vect_traits<OutputVector>::value_type(1e-12), 
                              typename vect_traits<OutputVector>::value_type epsy = typename vect_traits<OutputVector>::value_type(1e-12)) {
  return detail::levenberg_marquardt_nllsq_impl(f,fill_jac,x,y,Cholesky_linsolver(),no_limit_functor(),max_iter,tau,epsj,epsx,epsy);
};


/**
 * This function finds the non-linear least-square solution to a vector function with
 * limits imposed on the search domain. This function uses the Levenberg-Marquardt method
 * (which is a kind of adaptive damped-least square method, or trust-region Gauss-Newton).
 * This method performs  
 * very well for most non-linear systems (i.e. a Vanilla solution) and can cope with 
 * some level of rank-deficiency (but not recommended for highly rank-deficient 
 * jacobians), and for very well-conditioned problems, the basic Gauss-Newton method 
 * is often a better choice (see gauss_newton_nllsq).
 * TEST PASSED
 * \tparam Function The functor type.
 * \tparam JacobianFunction The functor type to fill the Jacobian matrix.
 * \tparam InputVector The vector type of the independent variable of the function.
 * \tparam OutputVector The vector type of the dependent variable of the function.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to match to the fixed output vector.
 * \param fill_jac The Jacobian of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param y The fixed vector of outputs to be matched in a least-square approximation.
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying independent vector-space.
 * \param max_iter The maximum number of iterations to perform.
 * \param tau The initial, relative damping value to damp the Hessian matrix.
 * \param epsj The tolerance on the norm of the jacobian.
 * \param epsx The tolerance on the norm of the estimation steps.
 * \param epsy The tolerance on the norm of the residual vector.
 * \return 1: Stopped by epsy (negl. residual), 2: Stopped by epsj (negl. gradient), 3: Stopped by epsx (negl. steps).
 * 
 * \throw improper_problem If the data of the problem is not proper (e.g. M < N).
 * \throw infeasible_problem If it becomes impossible to further damp the Jacobian (cannot make it non-singular), try to re-run with a larger "tau" value.
 * \throw singularity_error If a step is detected to cause a near-infinite jump.
 * \throw maximum_iteration If the maximum number of iterations has been reached.
 */
template <typename Function, typename JacobianFunction, 
          typename InputVector, typename OutputVector,
          typename LimitFunction>
int limited_levenberg_marquardt_nllsq(Function f, JacobianFunction fill_jac, 
                                      InputVector& x, const OutputVector& y,
                                      LimitFunction impose_limits,
                                      unsigned int max_iter = 100, 
                                      typename vect_traits<OutputVector>::value_type tau = typename vect_traits<OutputVector>::value_type(1E-03), 
                                      typename vect_traits<OutputVector>::value_type epsj = typename vect_traits<OutputVector>::value_type(1e-12), 
                                      typename vect_traits<OutputVector>::value_type epsx = typename vect_traits<OutputVector>::value_type(1e-12), 
                                      typename vect_traits<OutputVector>::value_type epsy = typename vect_traits<OutputVector>::value_type(1e-12)) {
  return detail::levenberg_marquardt_nllsq_impl(f,fill_jac,x,y,Cholesky_linsolver(),impose_limits,max_iter,tau,epsj,epsx,epsy);
};



};

};

#endif

















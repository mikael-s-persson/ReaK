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

#include "base/defs.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_norms.hpp"
#include "optim_exceptions.hpp"

#include "lin_alg/mat_cholesky.hpp"
#include "lin_alg/mat_damped_matrix.hpp"

#include "limit_functions.hpp"

namespace ReaK {
  
  
namespace optim {


namespace detail { 


template <typename Function, typename JacobianFunction, 
          typename InputVector, typename OutputVector, 
	  typename LimitFunction, typename T>
int levenberg_marquardt_nllsq_impl(Function f, JacobianFunction fill_jac, 
				    InputVector& x, const OutputVector& y, 
				    LimitFunction impose_limits, unsigned int itmax, 
				    T tau, T epsj, T epsx, T epsy)
{
  typedef typename vect_traits<InputVector>::value_type ValueType;
  typedef typename vect_traits<InputVector>::size_type SizeType;
  
  /* Check if the problem is defined properly */
  if ((y.size() < x.size()) ||
      (itmax <= 1))
    throw improper_problem("Levenberg-Marquardt requires M > N!");
  
  mat<ValueType,mat_structure::rectangular> J(y.size(),x.size());
  mat<ValueType,mat_structure::square> JtJ(x.size());
  mat<ValueType,mat_structure::diagonal> diag_JtJ(x.size());
  mat<ValueType,mat_structure::scalar> mu(x.size(),0.0);
  InputVector Jte = x;
  InputVector Dp = x; Dp -= x;
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

  /* compute e=x - f(p) and its L2 norm */
  OutputVector y_approx = f(x);
  OutputVector e = y; e -= y_approx;
  OutputVector e_tmp = e;
  ValueType p_eL2 = e * e;
  
  unsigned int nu = 2;
  for(unsigned int k = 0; k < itmax; ++k) {
    
    if(p_eL2 < epsy){ 
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
        Dp = Jte;
        linsolve_Cholesky(make_damped_matrix(JtJ,mu),Dp,epsj);
	
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
  throw maximum_iteration(itmax);

};


};


/**
 * This function finds the non-linear least-square solution to a vector function.
 * This function uses the Levenberg-Marquardt method (which is a kind of adaptive 
 * damped-least square method, or trust-region Gauss-Newton).
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
  return detail::levenberg_marquardt_nllsq_impl(f,fill_jac,x,y,no_limit_functor(),max_iter,tau,epsj,epsx,epsy);
};


/**
 * This function finds the non-linear least-square solution to a vector function with
 * limits imposed on the search domain. This function uses the Levenberg-Marquardt method
 * (which is a kind of adaptive damped-least square method, or trust-region Gauss-Newton).
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
  return detail::levenberg_marquardt_nllsq_impl(f,fill_jac,x,y,impose_limits,max_iter,tau,epsj,epsx,epsy);
};



};

};

#endif

















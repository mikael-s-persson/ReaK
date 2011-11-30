/**
 * \file gauss_newton_method.hpp
 *
 * The following library provides an implementation of the Gauss-Newton method.
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

#ifndef REAK_GAUSS_NEWTON_METHOD_HPP
#define REAK_GAUSS_NEWTON_METHOD_HPP

#include "base/defs.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_num_exceptions.hpp"

#include "lin_alg/mat_qr_decomp.hpp"

#include "limit_functions.hpp"


namespace ReaK {
  
  
namespace optim {


namespace detail {

  
template <typename Function, typename JacobianFunction, 
          typename InputVector, typename OutputVector, 
	  typename LinearLsqSolver, typename LimitFunction>
void gauss_newton_nllsq_impl(Function f, JacobianFunction fill_jac, InputVector& x, const OutputVector& y, 
			     LinearLsqSolver lin_solve, LimitFunction impose_limits, typename vect_traits<InputVector>::value_type tol = typename vect_traits<InputVector>::value_type(1e-6)) {
  typedef typename vect_traits<InputVector>::value_type ValueType;
  using std::sqrt; using std::fabs;
  
  OutputVector y_approx = f(x);
  OutputVector r = y - y_approx;
  mat<ValueType, mat_structure::rectangular> J(y.size(), x.size());
  fill_jac(J,x,y_approx);
  InputVector e = x;
  mat_vect_adaptor<InputVector> e_mat(e);
  lin_solve(J,e_mat,mat_vect_adaptor<OutputVector>(r),tol);
  impose_limits(x,e);
  while (norm(e) > tol) {
    x += e;
    y_approx = f(x);
    r = y; r -= y_approx;
    fill_jac(J,x,y_approx);
    e = r * J;
    lin_solve(J,e_mat,mat_vect_adaptor<OutputVector>(r),tol);
    impose_limits(x,e);
  };
};

};


/**
 * This function finds the non-linear least-square solution to a vector function.
 * \tparam Function The functor type.
 * \tparam JacobianFunction The functor type to fill the Jacobian matrix.
 * \tparam InputVector The vector type of the independent variable of the function.
 * \tparam OutputVector The vector type of the dependent variable of the function.
 * \param f The function to match to the fixed output vector.
 * \param fill_jac The Jacobian of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param y The fixed vector of outputs to be matched in a least-square approximation.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename JacobianFunction, 
          typename InputVector, typename OutputVector>
void gauss_newton_nllsq(Function f, JacobianFunction fill_jac, InputVector& x, const OutputVector& y, 
			     typename vect_traits<InputVector>::value_type tol = typename vect_traits<InputVector>::value_type(1e-6)) {
  detail::gauss_newton_nllsq_impl(f,fill_jac,x,y,QR_linlsqsolver(),no_limit_functor(),tol);
};

/**
 * This function finds the non-linear least-square solution to a vector function with
 * limits imposed on the search domain.
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
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename JacobianFunction, 
          typename InputVector, typename OutputVector, 
	  typename LimitFunction>
void limited_gauss_newton_nllsq(Function f, JacobianFunction fill_jac, InputVector& x, const OutputVector& y, 
			     LimitFunction impose_limits, typename vect_traits<InputVector>::value_type tol = typename vect_traits<InputVector>::value_type(1e-6)) {
  detail::gauss_newton_nllsq_impl(f,fill_jac,x,y,QR_linlsqsolver(),impose_limits,tol);
};



};

};


#endif


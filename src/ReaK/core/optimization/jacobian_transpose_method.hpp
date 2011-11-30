/**
 * \file jacobian_transpose_method.hpp
 *
 * The following library provides an implementation of the Jacobian-transpose method (or steepest 
 * descent for the non-linear least-square method).
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

#ifndef REAK_JACOBIAN_TRANSPOSE_METHOD_HPP
#define REAK_JACOBIAN_TRANSPOSE_METHOD_HPP

#include "base/defs.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_num_exceptions.hpp"

#include "limit_functions.hpp"


namespace ReaK {
  
  
namespace optim {


namespace detail {

  
template <typename Function, typename GradFunction, 
          typename InputVector, typename OutputVector, 
	  typename LimitFunction>
void jacobian_transpose_nllsq_impl(Function f, GradFunction fill_jac, InputVector& x, const OutputVector& y, 
			           LimitFunction impose_limits, typename vect_traits<InputVector>::value_type tol = typename vect_traits<InputVector>::value_type(1e-6)) {
  typedef typename vect_traits<InputVector>::value_type ValueType;
  using std::sqrt; using std::fabs;
  
  OutputVector y_approx = y;
  OutputVector r = y;
  mat<ValueType, mat_structure::rectangular> J(y.size(), x.size());
  InputVector e = x; e -= x;
  OutputVector Je = y;
  ValueType alpha = ValueType(1.0);
  do {
    x += e;
    y_approx = f(x);
    r = y; r -= y_approx;
    fill_jac(J,x,y_approx);
    e = r * J;
    Je = J * e;
    ValueType alpha = Je * Je;
    if(alpha < tol)
      return;
    alpha = (r * Je) / alpha;
    e *= alpha;
    impose_limits(x,e);
  } while(norm(e) > tol);
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
template <typename Function, typename GradFunction, 
          typename InputVector, typename OutputVector>
void jacobian_transpose_nllsq(Function f, GradFunction fill_jac, InputVector& x, const OutputVector& y, 
			      typename vect_traits<InputVector>::value_type tol = typename vect_traits<InputVector>::value_type(1e-6)) {
  detail::jacobian_transpose_nllsq_impl(f,fill_jac,x,y,no_limit_functor(),tol);
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
template <typename Function, typename GradFunction, 
          typename InputVector, typename OutputVector, 
	  typename LimitFunction>
void limited_jacobian_transpose_nllsq(Function f, GradFunction fill_jac, InputVector& x, const OutputVector& y, 
			     LimitFunction impose_limits, typename vect_traits<InputVector>::value_type tol = typename vect_traits<InputVector>::value_type(1e-6)) {
  detail::jacobian_transpose_nllsq_impl(f,fill_jac,x,y,impose_limits,tol);
};



};

};


#endif


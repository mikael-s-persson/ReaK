/**
 * \file bfgs_methods.hpp
 *
 * The following library provides methods to perform non-linear optimization based on the BFGS method.
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

#ifndef REAK_BFGS_METHODS_HPP
#define REAK_BFGS_METHODS_HPP

#include "base/defs.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_num_exceptions.hpp"

#include "hessian_approx_update.hpp"
#include "line_search.hpp"

#include "quasi_newton_methods.hpp"

namespace ReaK {
  
  
namespace optim {

#if 0
/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search 
 * direction, a BFGS approximation of the Hessian (without restarts), and using a line-search approach 
 * (the line-search is a "sloppy" expand-and-zoom approach that looks to satisfy the strong Wolfe conditions).
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename Vector>
void bfgs_method(Function f, GradFunction df, Vector& x, typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  quasi_newton_line_search(f,df,x,line_search_expand_and_zoom<typename vect_traits<Vector>::value_type>(1e-4,0.9),inv_hessian_update_bfgs(),tol);
  
};

/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search 
 * direction, a DFP approximation of the Hessian (without restarts), and using a line-search approach 
 * (the line-search is a "sloppy" expand-and-zoom approach that looks to satisfy the strong Wolfe conditions).
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename Vector>
void dfp_method(Function f, GradFunction df, Vector& x, typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  quasi_newton_line_search(f,df,x,line_search_expand_and_zoom<typename vect_traits<Vector>::value_type>(1e-4,0.9),inv_hessian_update_dfp(),tol);
  
};

/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search 
 * direction, a Broyden-class approximation of the Hessian (without restarts), and using a line-search approach 
 * (the line-search is a "sloppy" expand-and-zoom approach that looks to satisfy the strong Wolfe conditions).
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param phi The fraction to use in the Broyden-class hessian approximation (0: BFGS only, 1: DFP only).
 * \param tol The tolerance on the norm of the gradient (and thus the step size).
 */
template <typename Function, typename GradFunction, typename Vector>
void broyden_class_method(Function f, GradFunction df, Vector& x, typename vect_traits<Vector>::value_type phi = typename vect_traits<Vector>::value_type(0.5), typename vect_traits<Vector>::value_type tol = typename vect_traits<Vector>::value_type(1e-6)) {
  
  quasi_newton_line_search(f,df,x,line_search_expand_and_zoom<typename vect_traits<Vector>::value_type>(1e-4,0.9),inv_hessian_update_broyden<typename vect_traits<Vector>::value_type>(phi),tol);
  
};
#endif

};

};


#endif


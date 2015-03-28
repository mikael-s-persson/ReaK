/**
 * \file quasi_newton_methods.hpp
 *
 * The following library provides methods to perform non-linear optimization based on Quasi-Newton methods.
 * The methods available in this library include the bgfs_method (recommended as a line-search option),
 * the dfp_method, the broyden_class_method, and the sr1_tr_method (recommended as a trust-region option).
 * All the provided methods are based on two quasi-newton optimization loops (one based on line-search
 * and one based on trust-region), and thus, the user also has the option to use those underlying
 * algorithms directly using various functor types that solve the sub-problems (line-search, trust-region
 * search, approximate hessian updates, limit functions).
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

#ifndef REAK_QUASI_NEWTON_METHODS_HPP
#define REAK_QUASI_NEWTON_METHODS_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/math/lin_alg/mat_num_exceptions.hpp>

#include "trust_region_search.hpp"
#include "line_search.hpp"
#include "hessian_approx_update.hpp"
#include "limit_functions.hpp"


namespace ReaK {


namespace optim {


/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search
 * direction, an approximation of the Hessian (without restarts), and using a line-search approach
 * (satisfying the strong Wolfe conditions). Note that this function is the underlying optimization
 * loop used by most quasi-Newton methods (see bfgs_method, dfp_method, and broyden_class_method).
 * TEST PASSED
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam LineSearcher A functor type that can perform a line-search which satisfy the strong Wolfe conditions.
 * \tparam InvHessianUpdater A functor type that can update the inverse of a Hessian given the changes in the solution
 * and the function gradient.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_iter The maximum number of iterations to perform.
 * \param get_alpha The functor that can perform a line-search which satisfy the strong Wolfe conditions.
 * \param update_inv_hessian The functor that can update the inverse of a Hessian given the changes in the solution and
 * the function gradient.
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying
 * independent vector-space.
 * \param abs_tol The tolerance on the norm of the steps.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template < typename Function, typename GradFunction, typename Vector, typename LineSearcher, typename InvHessianUpdater,
           typename LimitFunction >
void quasi_newton_line_search( Function f, GradFunction df, Vector& x, unsigned int max_iter, LineSearcher get_alpha,
                               InvHessianUpdater update_inv_hessian, LimitFunction impose_limits,
                               typename vect_traits< Vector >::value_type abs_tol
                               = typename vect_traits< Vector >::value_type( 1e-6 ),
                               typename vect_traits< Vector >::value_type abs_grad_tol
                               = typename vect_traits< Vector >::value_type( 1e-6 ) ) {
  typedef typename vect_traits< Vector >::value_type ValueType;
  using std::sqrt;
  using std::fabs;

  ValueType x_value = f( x );
  Vector x_grad = -df( x );

  Vector p = x_grad;
  ValueType alpha = ValueType( 1.0 );
  ValueType pxg = p * x_grad;
  if( ( f( x + p ) > x_value - ValueType( 1e-4 ) * pxg )
      || ( fabs( p * df( x + p ) ) > ValueType( 0.9 ) * fabs( pxg ) ) )
    alpha = get_alpha( f, df, ValueType( 0.0 ), ValueType( 2.0 ), x, p, abs_tol );

  Vector s = alpha * p;
  impose_limits( x, s );
  Vector y = x_grad;
  x += s;
  x_value = f( x );
  x_grad = -df( x );
  y -= x_grad;
  mat< ValueType, mat_structure::symmetric > H(
    mat< ValueType, mat_structure::scalar >( x.size(), ( y * s ) / ( y * y ) ) );

  unsigned int k = 0;

  while( norm_2( x_grad ) > abs_grad_tol ) {
    update_inv_hessian( H, s, y );
    p = H * x_grad;
    // check Wolfe for alpha 1.0
    alpha = ValueType( 1.0 );
    pxg = p * x_grad;
    if( ( f( x + p ) > x_value - ValueType( 1e-4 ) * pxg )
        || ( fabs( p * df( x + p ) ) > ValueType( 0.9 ) * fabs( pxg ) ) )
      alpha = get_alpha( f, df, ValueType( 0.0 ), ValueType( 2.0 ), x, p, abs_tol );
    s = alpha * p;
    impose_limits( x, s );
    if( norm_2( s ) > abs_tol )
      return;
    x += s;

    if( ++k > max_iter )
      throw maximum_iteration( max_iter );

    y = x_grad;
    x_value = f( x );
    x_grad = -df( x );
    y -= x_grad;
  };
};


/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search
 * direction, an approximation of the Hessian (without restarts), and using a trust-region approach.
 * Note that this function is the underlying optimization loop used by most quasi-Newton
 * methods (see bfgs_method, dfp_method, and broyden_class_method).
 * TEST PASSED
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam TrustRegionSolver A functor type that can find the solution to the approximate quadratic over the current
 * trust-region.
 * \tparam HessianUpdater A functor type that can update the Hessian given the changes in the solution and the function
 * gradient.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param max_iter The maximum number of iterations to perform.
 * \param solve_step The functor that can find the solution to the approximate quadratic over the current trust-region.
 * \param update_hessian The functor that can update the Hessian given the changes in the solution and the function
 * gradient.
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying
 * independent vector-space.
 * \param abs_tol The tolerance on the norm of the steps.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 * \param eta The tolerance on the ratio between actual reduction and predicted reduction in order to accept a given
 * step.
 * \param r_tol The tolerance on the alignment of the step and the residual in order to trigger an update of the
 * approximate Hessian matrix.
 */
template < typename Function, typename GradFunction, typename Vector, typename TrustRegionSolver,
           typename HessianUpdater, typename LimitFunction >
void quasi_newton_trust_region(
  Function f, GradFunction df, Vector& x, typename vect_traits< Vector >::value_type max_radius, unsigned int max_iter,
  TrustRegionSolver solve_step, HessianUpdater update_hessian, LimitFunction impose_limits,
  typename vect_traits< Vector >::value_type abs_tol = typename vect_traits< Vector >::value_type( 1e-6 ),
  typename vect_traits< Vector >::value_type abs_grad_tol = typename vect_traits< Vector >::value_type( 1e-6 ),
  typename vect_traits< Vector >::value_type eta = typename vect_traits< Vector >::value_type( 1e-4 ),
  typename vect_traits< Vector >::value_type r_tol = typename vect_traits< Vector >::value_type( 1e-6 ) ) {
  typedef typename vect_traits< Vector >::value_type ValueType;
  using std::sqrt;
  using std::fabs;

  ValueType radius = ValueType( 0.5 ) * max_radius;

  ValueType x_value = f( x );
  Vector x_grad = df( x );

  Vector p = -x_grad;
  ValueType norm_p;
  mat< ValueType, mat_structure::square > B( mat< ValueType, mat_structure::identity >( x.size() ) );
  solve_step( x_grad, B, p, norm_p, radius, abs_tol );
  impose_limits( x, p );
  Vector xt = x;
  xt += p;
  ValueType xt_value = f( xt );
  Vector xt_grad = df( xt );
  ValueType aredux = x_value - xt_value;
  ValueType predux = -( x_grad * p + ValueType( 0.5 ) * ( p * ( B * p ) ) );

  Vector y = xt_grad;
  y -= x_grad;
  B = mat< ValueType, mat_structure::scalar >( x.size(), ( y * p ) / ( p * p ) );
  Vector y_Bp = y;
  y_Bp -= B * p;

  while( true ) {
    if( fabs( p * y_Bp ) >= r_tol * norm_p * norm_2( y_Bp ) )
      update_hessian( B, p, y );

    ValueType ratio = aredux / predux;
    if( ratio > ValueType( 0.75 ) ) {
      if( norm_p > ValueType( 0.8 ) * radius ) {
        radius *= ValueType( 2.0 );
        if( radius > max_radius )
          radius = max_radius;
      };
    } else if( ratio < ValueType( 0.1 ) ) {
      radius *= ValueType( 0.5 );
    };
    if( ratio > eta ) {
      x = xt;
      x_value = xt_value;
      x_grad = xt_grad;
    };
    if( norm_2( x_grad ) < abs_grad_tol )
      return;
    solve_step( x_grad, B, p, norm_p, radius, abs_tol );
    impose_limits( x, p );
    xt = x;
    xt += p;
    if( norm_2( p ) < abs_tol )
      return;
    xt_value = f( xt );
    xt_grad = df( xt );
    aredux = x_value - xt_value;
    predux = -( x_grad * p + ValueType( 0.5 ) * ( p * ( B * p ) ) );
    y = xt_grad;
    y -= x_grad;
    y_Bp = y;
    y_Bp -= B * p;
  };
};


/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search
 * direction, a BFGS approximation of the Hessian (without restarts), and using a line-search approach
 * (the line-search is a "sloppy" expand-and-zoom approach that looks to satisfy the strong Wolfe conditions).
 * TEST PASSED
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance on the norm of the steps.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template < typename Function, typename GradFunction, typename Vector >
void bfgs_method( Function f, GradFunction df, Vector& x, unsigned int max_iter,
                  typename vect_traits< Vector >::value_type abs_tol
                  = typename vect_traits< Vector >::value_type( 1e-6 ),
                  typename vect_traits< Vector >::value_type abs_grad_tol
                  = typename vect_traits< Vector >::value_type( 1e-6 ) ) {

  quasi_newton_line_search( f, df, x, max_iter,
                            line_search_expand_and_zoom< typename vect_traits< Vector >::value_type >( 1e-4, 0.9 ),
                            inv_hessian_update_bfgs(), no_limit_functor(), abs_tol, abs_grad_tol );
};

/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search
 * direction, a BFGS approximation of the Hessian (without restarts), and using a line-search approach
 * (the line-search is a "sloppy" expand-and-zoom approach that looks to satisfy the strong Wolfe conditions).
 * This overload version allows for simple limits to be imposed on the steps.
 * TEST PASSED
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_iter The maximum number of iterations to perform.
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying
 * independent vector-space.
 * \param abs_tol The tolerance on the norm of the steps.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template < typename Function, typename GradFunction, typename Vector, typename LimitFunction >
void limited_bfgs_method( Function f, GradFunction df, Vector& x, unsigned int max_iter, LimitFunction impose_limits,
                          typename vect_traits< Vector >::value_type abs_tol
                          = typename vect_traits< Vector >::value_type( 1e-6 ),
                          typename vect_traits< Vector >::value_type abs_grad_tol
                          = typename vect_traits< Vector >::value_type( 1e-6 ) ) {

  quasi_newton_line_search( f, df, x, max_iter,
                            line_search_expand_and_zoom< typename vect_traits< Vector >::value_type >( 1e-4, 0.9 ),
                            inv_hessian_update_bfgs(), impose_limits, abs_tol, abs_grad_tol );
};

/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search
 * direction, a DFP approximation of the Hessian (without restarts), and using a line-search approach
 * (the line-search is a "sloppy" expand-and-zoom approach that looks to satisfy the strong Wolfe conditions).
 * TEST PASSED
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance on the norm of the steps.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template < typename Function, typename GradFunction, typename Vector >
void dfp_method( Function f, GradFunction df, Vector& x, unsigned int max_iter,
                 typename vect_traits< Vector >::value_type abs_tol
                 = typename vect_traits< Vector >::value_type( 1e-6 ),
                 typename vect_traits< Vector >::value_type abs_grad_tol
                 = typename vect_traits< Vector >::value_type( 1e-6 ) ) {

  quasi_newton_line_search( f, df, x, max_iter,
                            line_search_expand_and_zoom< typename vect_traits< Vector >::value_type >( 1e-4, 0.9 ),
                            inv_hessian_update_dfp(), no_limit_functor(), abs_tol, abs_grad_tol );
};

/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search
 * direction, a DFP approximation of the Hessian (without restarts), and using a line-search approach
 * (the line-search is a "sloppy" expand-and-zoom approach that looks to satisfy the strong Wolfe conditions).
 * This overload version allows for simple limits to be imposed on the steps.
 * TEST PASSED
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_iter The maximum number of iterations to perform.
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying
 * independent vector-space.
 * \param abs_tol The tolerance on the norm of the steps.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template < typename Function, typename GradFunction, typename Vector, typename LimitFunction >
void limited_dfp_method( Function f, GradFunction df, Vector& x, unsigned int max_iter, LimitFunction impose_limits,
                         typename vect_traits< Vector >::value_type abs_tol
                         = typename vect_traits< Vector >::value_type( 1e-6 ),
                         typename vect_traits< Vector >::value_type abs_grad_tol
                         = typename vect_traits< Vector >::value_type( 1e-6 ) ) {

  quasi_newton_line_search( f, df, x, max_iter,
                            line_search_expand_and_zoom< typename vect_traits< Vector >::value_type >( 1e-4, 0.9 ),
                            inv_hessian_update_dfp(), impose_limits, abs_tol, abs_grad_tol );
};

/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search
 * direction, a Broyden-class approximation of the Hessian (without restarts), and using a line-search approach
 * (the line-search is a "sloppy" expand-and-zoom approach that looks to satisfy the strong Wolfe conditions).
 * TEST PASSED
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_iter The maximum number of iterations to perform.
 * \param phi The fraction to use in the Broyden-class hessian approximation (0: BFGS only, 1: DFP only).
 * \param abs_tol The tolerance on the norm of the steps.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template < typename Function, typename GradFunction, typename Vector >
void broyden_class_method(
  Function f, GradFunction df, Vector& x, unsigned int max_iter,
  typename vect_traits< Vector >::value_type phi = typename vect_traits< Vector >::value_type( 0.5 ),
  typename vect_traits< Vector >::value_type abs_tol = typename vect_traits< Vector >::value_type( 1e-6 ),
  typename vect_traits< Vector >::value_type abs_grad_tol = typename vect_traits< Vector >::value_type( 1e-6 ) ) {

  quasi_newton_line_search( f, df, x, max_iter,
                            line_search_expand_and_zoom< typename vect_traits< Vector >::value_type >( 1e-4, 0.9 ),
                            inv_hessian_update_broyden< typename vect_traits< Vector >::value_type >( phi ),
                            no_limit_functor(), abs_tol, abs_grad_tol );
};

/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search
 * direction, a Broyden-class approximation of the Hessian (without restarts), and using a line-search approach
 * (the line-search is a "sloppy" expand-and-zoom approach that looks to satisfy the strong Wolfe conditions).
 * This overload version allows for simple limits to be imposed on the steps.
 * TEST PASSED
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_iter The maximum number of iterations to perform.
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying
 * independent vector-space.
 * \param phi The fraction to use in the Broyden-class hessian approximation (0: BFGS only, 1: DFP only).
 * \param abs_tol The tolerance on the norm of the steps.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template < typename Function, typename GradFunction, typename Vector, typename LimitFunction >
void limited_broyden_class_method(
  Function f, GradFunction df, Vector& x, unsigned int max_iter, LimitFunction impose_limits,
  typename vect_traits< Vector >::value_type phi = typename vect_traits< Vector >::value_type( 0.5 ),
  typename vect_traits< Vector >::value_type abs_tol = typename vect_traits< Vector >::value_type( 1e-6 ),
  typename vect_traits< Vector >::value_type abs_grad_tol = typename vect_traits< Vector >::value_type( 1e-6 ) ) {

  quasi_newton_line_search( f, df, x, max_iter,
                            line_search_expand_and_zoom< typename vect_traits< Vector >::value_type >( 1e-4, 0.9 ),
                            inv_hessian_update_broyden< typename vect_traits< Vector >::value_type >( phi ),
                            impose_limits, abs_tol, abs_grad_tol );
};


/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search
 * direction, a symmetric rank-1 (SR-1) approximation of the Hessian (without restarts), and using a
 * trust-region approach (with Dogleg solver).
 * TEST PASSED
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_radius The maximum radius of the trust-region.
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance on the norm of the steps.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template < typename Function, typename GradFunction, typename Vector >
void sr1_tr_method( Function f, GradFunction df, Vector& x, typename vect_traits< Vector >::value_type max_radius
                                                            = typename vect_traits< Vector >::value_type( 1.0 ),
                    unsigned int max_iter = 100, typename vect_traits< Vector >::value_type abs_tol
                                                 = typename vect_traits< Vector >::value_type( 1e-6 ),
                    typename vect_traits< Vector >::value_type abs_grad_tol
                    = typename vect_traits< Vector >::value_type( 1e-6 ) ) {

  quasi_newton_trust_region( f, df, x, max_radius, max_iter, trust_region_solver_dogleg(), hessian_update_sr1(),
                             no_limit_functor(), abs_tol, abs_grad_tol );
};

/**
 * This function finds the minimum of a function, given its derivative, using a quasi-newton search
 * direction, a symmetric rank-1 (SR-1) approximation of the Hessian (without restarts), and using a
 * trust-region approach (with Dogleg solver).
 * This overload version allows for simple limits to be imposed on the steps.
 * TEST PASSED
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying
 * independent vector-space.
 * \param max_radius The maximum radius of the trust-region.
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance on the norm of the steps.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template < typename Function, typename GradFunction, typename Vector, typename LimitFunction >
void limited_sr1_tr_method( Function f, GradFunction df, Vector& x, LimitFunction impose_limits,
                            typename vect_traits< Vector >::value_type max_radius
                            = typename vect_traits< Vector >::value_type( 1.0 ),
                            unsigned int max_iter = 100, typename vect_traits< Vector >::value_type abs_tol
                                                         = typename vect_traits< Vector >::value_type( 1e-6 ),
                            typename vect_traits< Vector >::value_type abs_grad_tol
                            = typename vect_traits< Vector >::value_type( 1e-6 ) ) {

  quasi_newton_trust_region( f, df, x, max_radius, max_iter, trust_region_solver_dogleg(), hessian_update_sr1(),
                             impose_limits, abs_tol, abs_grad_tol );
};
};
};


#endif

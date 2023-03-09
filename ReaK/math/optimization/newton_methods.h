/**
 * \file newton_methods.h
 *
 * The following library provides methods to perform non-linear optimization based on Newton methods.
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

#ifndef REAK_MATH_OPTIMIZATION_NEWTON_METHODS_H_
#define REAK_MATH_OPTIMIZATION_NEWTON_METHODS_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_num_exceptions.h"

#include "ReaK/math/optimization/limit_functions.h"
#include "ReaK/math/optimization/line_search.h"
#include "ReaK/math/optimization/newton_search_directions.h"
#include "ReaK/math/optimization/trust_region_search.h"

namespace ReaK::optim {

namespace detail {

template <typename Function, typename GradFunction, typename HessianFunction,
          typename Vector, typename LineSearcher, typename LimitFunction,
          typename NewtonDirectioner>
void newton_method_ls_impl(
    Function f, GradFunction df, HessianFunction fill_hessian, Vector& x,
    unsigned int max_iter, LineSearcher get_alpha, LimitFunction impose_limits,
    NewtonDirectioner get_direction,
    vect_value_type_t<Vector> abs_tol = vect_value_type_t<Vector>(1e-6),
    vect_value_type_t<Vector> abs_grad_tol = vect_value_type_t<Vector>(1e-6)) {
  using ValueType = vect_value_type_t<Vector>;
  using std::abs;
  using std::sqrt;

  ValueType x_value = f(x);
  Vector x_grad = df(x);
  ValueType x_grad_norm = norm_2(x_grad);

  Vector p = x_grad;
  mat<ValueType, mat_structure::symmetric> H(
      mat<ValueType, mat_structure::identity>(x.size()));
  fill_hessian(H, x, x_value, x_grad);

  unsigned int k = 0;

  while (x_grad_norm > abs_grad_tol) {
    get_direction(H, x_grad, p, abs_tol);
    // check Wolfe for alpha 1.0
    auto alpha = ValueType(1.0);
    ValueType pxg = p * x_grad;
    if ((f(x + p) > x_value + ValueType(1e-4) * pxg) ||
        (abs(p * df(x + p)) > ValueType(0.8) * abs(pxg))) {
      alpha = get_alpha(f, df, ValueType(0.0), ValueType(2.0), x, p, abs_tol);
    }
    p *= alpha;
    impose_limits(x, p);
    if (norm_2(p) > abs_tol) {
      return;
    }

    x += p;

    if (++k > max_iter) {
      throw maximum_iteration(max_iter);
    }

    x_value = f(x);
    x_grad = df(x);
    x_grad_norm = norm_2(x_grad);
    fill_hessian(H, x, x_value, x_grad);
  }
}

template <typename Function, typename GradFunction, typename HessianFunction,
          typename Vector, typename TrustRegionSolver, typename LimitFunction>
void newton_method_tr_impl(
    Function f, GradFunction df, HessianFunction fill_hessian, Vector& x,
    vect_value_type_t<Vector> max_radius, unsigned int max_iter,
    TrustRegionSolver solve_step, LimitFunction impose_limits,
    vect_value_type_t<Vector> abs_tol = vect_value_type_t<Vector>(1e-6),
    vect_value_type_t<Vector> abs_grad_tol = vect_value_type_t<Vector>(1e-6),
    vect_value_type_t<Vector> eta = vect_value_type_t<Vector>(1e-4)) {
  using ValueType = vect_value_type_t<Vector>;
  using std::abs;
  using std::sqrt;

  ValueType radius = ValueType(0.5) * max_radius;
  ValueType x_value = f(x);
  Vector x_grad = df(x);
  ValueType x_grad_norm = norm_2(x_grad);

  Vector p = x_grad;
  ValueType norm_p = std::numeric_limits<ValueType>::max();
  mat<ValueType, mat_structure::symmetric> H(x.size());
  fill_hessian(H, x, x_value, x_grad);

  Vector xt = x;

  unsigned int k = 0;

  while (x_grad_norm > abs_grad_tol) {
    solve_step(x_grad, H, p, norm_p, radius, abs_tol);
    impose_limits(x, p);
    xt = x;
    xt += p;
    norm_p = norm_2(p);
    ValueType xt_value = f(xt);
    ValueType aredux = x_value - xt_value;
    ValueType predux = -(x_grad * p + ValueType(0.5) * (p * (H * p)));

    ValueType ratio = aredux / predux;
    if (ratio > ValueType(0.75)) {
      if (norm_p > ValueType(0.8) * radius) {
        radius *= ValueType(2.0);
        if (radius > max_radius) {
          radius = max_radius;
        }
      }
    } else if (ratio < ValueType(0.1)) {
      radius *= ValueType(0.5);
    }
    if (ratio > eta) {  // the step is accepted.
      x = xt;
      if (norm_p < abs_tol) {
        return;
      }
      x_value = xt_value;
      x_grad = df(x);
      x_grad_norm = norm_2(x_grad);
      fill_hessian(H, x, x_value, x_grad);
    }
    if (++k > max_iter) {
      throw maximum_iteration(max_iter);
    }
  }
}
}  // namespace detail

/**
 * This functor is a factory class to construct a Newton-method optimizer routine that uses a
 * line-search approach. Use make_newton_method_ls to
 * construct this without having to specify the template arguments explicitly.
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \tparam LimitFunction A functor type that can impose limits on a proposed solution step (see no_limit_functor or
 * box_limit_function for examples).
 * \tparam NewtonDirectioner A functor type that can solve for a search direction (see newton_directioner or
 * regularized_newton_directioner for examples).
 */
template <typename Function, typename GradFunction, typename HessianFunction,
          typename T, typename LimitFunction = no_limit_functor,
          typename NewtonDirectioner = newton_directioner>
struct newton_method_ls_factory {
  Function f;
  GradFunction df;
  HessianFunction fill_hessian;
  unsigned int max_iter;
  T abs_tol;
  T abs_grad_tol;
  LimitFunction impose_limits;
  NewtonDirectioner get_direction;

  /**
   * Parametrized constructor of the factory object.
   * \param aF The function to minimize.
   * \param aDf The gradient of the function to minimize.
   * \param aFillHessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aTol The tolerance on the norm of step size.
   * \param aGradTol The tolerance on the norm of the gradient.
   * \param aImposeLimits The functor that can impose simple limits on the search domain (i.e. using this boils down to
   * a gradient projection method, for more complex constraints please use a constraint optimization method instead).
   * \param aGetDirection The functor that can solve for the search direction.
   */
  newton_method_ls_factory(
      Function aF, GradFunction aDf, HessianFunction aFillHessian,
      unsigned int aMaxIter = 100, T aTol = T(1e-6), T aGradTol = T(1e-6),
      LimitFunction aImposeLimits = LimitFunction(),
      NewtonDirectioner aGetDirection = NewtonDirectioner())
      : f(aF),
        df(aDf),
        fill_hessian(aFillHessian),
        max_iter(aMaxIter),
        abs_tol(aTol),
        abs_grad_tol(aGradTol),
        impose_limits(aImposeLimits),
        get_direction(aGetDirection) {}
  /**
   * This function finds the minimum of a function, given its derivative and Hessian,
   * using a newton search direction and using a trust-region approach.
   * \tparam Vector The vector type of the independent variable for the function.
   * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
   */
  template <typename Vector>
  void operator()(Vector& x) const {
    detail::newton_method_ls_impl(f, df, fill_hessian, x, max_iter,
                                  line_search_expand_and_zoom<T>(1e-4, 0.8),
                                  impose_limits, get_direction, abs_tol,
                                  abs_grad_tol);
  }

  /**
   * This function remaps the factory to one which will use a regularized solver for the search-direction.
   * You should regularize the matrix only if there are reasons to expect the Hessian to be near-singular.
   * \param tau The initial relative damping factor to regularize the Hessian matrix.
   */
  newton_method_ls_factory<Function, GradFunction, HessianFunction, T,
                           LimitFunction, regularized_newton_directioner<T>>
  regularize(const T& tau) const {
    return newton_method_ls_factory<Function, GradFunction, HessianFunction, T,
                                    LimitFunction,
                                    regularized_newton_directioner<T>>(
        f, df, fill_hessian, max_iter, abs_tol, abs_grad_tol, impose_limits,
        regularized_newton_directioner<T>(tau));
  }

  /**
   * This function remaps the factory to one which will use the given solver for the search direction.
   * \tparam NewNewtonDirectioner A functor type that can solve for a search direction (see newton_directioner or
   * regularized_newton_directioner for examples).
   * \param new_directioner The functor that can solve for the search direction.
   */
  template <typename NewNewtonDirectioner>
  newton_method_ls_factory<Function, GradFunction, HessianFunction, T,
                           LimitFunction, NewNewtonDirectioner>
  set_directioner(NewNewtonDirectioner new_directioner) const {
    return newton_method_ls_factory<Function, GradFunction, HessianFunction, T,
                                    LimitFunction, NewNewtonDirectioner>(
        f, df, fill_hessian, max_iter, abs_tol, abs_grad_tol, impose_limits,
        new_directioner);
  }

  /**
   * This function remaps the factory to one which will use the given limit-function for the search domain.
   * Using a limit-function boils down to a gradient projection method, and thus, for more complex
   * constraints (such as non-linear ones or mix of equality-inequality constraints), please use a
   * constraint optimization method instead (see augmented_lagrangian_methods.hpp for example).
   * \tparam NewLimitFunction A new functor type that can impose limits on a proposed solution step (see
   * no_limit_functor or box_limit_function for examples).
   * \param new_limits The functor that can impose simple limits on the search domain (i.e. using this boils down to a
   * gradient projection method, for more complex constraints please use a constraint optimization method instead).
   */
  template <typename NewLimitFunction>
  newton_method_ls_factory<Function, GradFunction, HessianFunction, T,
                           NewLimitFunction, NewtonDirectioner>
  set_limiter(NewLimitFunction new_limits) const {
    return newton_method_ls_factory<Function, GradFunction, HessianFunction, T,
                                    NewLimitFunction, NewtonDirectioner>(
        f, df, fill_hessian, max_iter, abs_tol, abs_grad_tol, new_limits,
        get_direction);
  }
};

/**
 * This function template creates a factory object to construct a Newton-method optimizer routine
 * that uses a line-search approach.
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance on the norm of the step size.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template <typename Function, typename GradFunction, typename HessianFunction>
newton_method_ls_factory<Function, GradFunction, HessianFunction, double>
make_newton_method_ls(Function f, GradFunction df, HessianFunction fill_hessian,
                      unsigned int max_iter = 100, double abs_tol = 1e-6,
                      double abs_grad_tol = 1e-6) {
  return newton_method_ls_factory<Function, GradFunction, HessianFunction,
                                  double>(f, df, fill_hessian, max_iter,
                                          abs_tol, abs_grad_tol);
}

/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a newton search
 * direction and using a line-search approach (satisfying the strong Wolfe conditions).
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance on the norm of the step size.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template <typename Function, typename GradFunction, typename HessianFunction,
          typename Vector>
void newton_method_ls(
    Function f, GradFunction df, HessianFunction fill_hessian, Vector& x,
    unsigned int max_iter,
    vect_value_type_t<Vector> abs_tol = vect_value_type_t<Vector>(1e-6),
    vect_value_type_t<Vector> abs_grad_tol = vect_value_type_t<Vector>(1e-6)) {

  detail::newton_method_ls_impl(
      f, df, fill_hessian, x, max_iter,
      line_search_expand_and_zoom<vect_value_type_t<Vector>>(1e-4, 0.4),
      no_limit_functor(), newton_directioner(), abs_tol, abs_grad_tol);
}

/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a regularized
 * newton search direction and using a line-search approach (satisfying the strong Wolfe conditions).
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_iter The maximum number of iterations to perform.
 * \param tau The initial relative damping factor to regularize the Hessian matrix.
 * \param abs_tol The tolerance on the norm of the step size.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template <typename Function, typename GradFunction, typename HessianFunction,
          typename Vector>
void reg_newton_method_ls(
    Function f, GradFunction df, HessianFunction fill_hessian, Vector& x,
    unsigned int max_iter,
    vect_value_type_t<Vector> tau = vect_value_type_t<Vector>(1e-3),
    vect_value_type_t<Vector> abs_tol = vect_value_type_t<Vector>(1e-6),
    vect_value_type_t<Vector> abs_grad_tol = vect_value_type_t<Vector>(1e-6)) {

  detail::newton_method_ls_impl(
      f, df, fill_hessian, x, max_iter,
      line_search_expand_and_zoom<vect_value_type_t<Vector>>(1e-4, 0.4),
      no_limit_functor(),
      regularized_newton_directioner<vect_value_type_t<Vector>>(tau), abs_tol,
      abs_grad_tol);
}

/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a newton search
 * direction, a function to impose simple limits on the search domain, and a line-search
 * approach (satisfying the strong Wolfe conditions).
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_iter The maximum number of iterations to perform.
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying
 * independent vector-space.
 * \param abs_tol The tolerance on the norm of the step size.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template <typename Function, typename GradFunction, typename HessianFunction,
          typename Vector, typename LimitFunction>
void limited_newton_method_ls(
    Function f, GradFunction df, HessianFunction fill_hessian, Vector& x,
    unsigned int max_iter, LimitFunction impose_limits,
    vect_value_type_t<Vector> abs_tol = vect_value_type_t<Vector>(1e-6),
    vect_value_type_t<Vector> abs_grad_tol = vect_value_type_t<Vector>(1e-6)) {

  detail::newton_method_ls_impl(
      f, df, fill_hessian, x, max_iter,
      line_search_expand_and_zoom<vect_value_type_t<Vector>>(1e-4, 0.4),
      impose_limits, newton_directioner(), abs_tol, abs_grad_tol);
}

/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a regularized
 * newton search direction, a function to impose simple limits on the search domain, and a line-search
 * approach (satisfying the strong Wolfe conditions).
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_iter The maximum number of iterations to perform.
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying
 * independent vector-space.
 * \param tau The initial relative damping factor to regularize the Hessian matrix.
 * \param abs_tol The tolerance on the norm of the step size.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template <typename Function, typename GradFunction, typename HessianFunction,
          typename Vector, typename LimitFunction>
void limited_reg_newton_method_ls(
    Function f, GradFunction df, HessianFunction fill_hessian, Vector& x,
    unsigned int max_iter, LimitFunction impose_limits,
    vect_value_type_t<Vector> tau = vect_value_type_t<Vector>(1e-3),
    vect_value_type_t<Vector> abs_tol = vect_value_type_t<Vector>(1e-6),
    vect_value_type_t<Vector> abs_grad_tol = vect_value_type_t<Vector>(1e-6)) {

  detail::newton_method_ls_impl(
      f, df, fill_hessian, x, max_iter,
      line_search_expand_and_zoom<vect_value_type_t<Vector>>(1e-4, 0.4),
      impose_limits,
      regularized_newton_directioner<vect_value_type_t<Vector>>(tau), abs_tol,
      abs_grad_tol);
}

/**
 * This functor is a factory class to construct a Newton-method optimizer routine that uses a
 * trust-region approach. Use make_newton_method_tr to
 * construct this without having to specify the template arguments explicitly.
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \tparam TrustRegionSolver A functor type that can solve for a solution step within a trust-region (see
 * trust_region_solver_dogleg for an example).
 * \tparam LimitFunction A functor type that can impose limits on a proposed solution step (see no_limit_functor or
 * box_limit_function for examples).
 */
template <typename Function, typename GradFunction, typename HessianFunction,
          typename T, typename TrustRegionSolver = trust_region_solver_dogleg,
          typename LimitFunction = no_limit_functor>
struct newton_method_tr_factory {
  Function f;
  GradFunction df;
  HessianFunction fill_hessian;
  T max_radius;
  unsigned int max_iter;
  T abs_tol;
  T abs_grad_tol;
  T eta;
  TrustRegionSolver solve_step;
  LimitFunction impose_limits;

  /**
   * Parametrized constructor of the factory object.
   * \param aF The function to minimize.
   * \param aDf The gradient of the function to minimize.
   * \param aFillHessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aTol The tolerance on the norm of the step size.
   * \param aGradTol The tolerance on the norm of the gradient.
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aSolveStep The functor that can solve for the step to take within the trust-region.
   * \param aImposeLimits The functor that can impose simple limits on the search domain (i.e. using this boils down to
   * a gradient projection method, for more complex constraints please use a constraint optimization method instead).
   */
  newton_method_tr_factory(Function aF, GradFunction aDf,
                           HessianFunction aFillHessian, T aMaxRadius,
                           unsigned int aMaxIter, T aTol = T(1e-6),
                           T aGradTol = T(1e-6), T aEta = T(1e-4),
                           TrustRegionSolver aSolveStep = TrustRegionSolver(),
                           LimitFunction aImposeLimits = LimitFunction())
      : f(aF),
        df(aDf),
        fill_hessian(aFillHessian),
        max_radius(aMaxRadius),
        max_iter(aMaxIter),
        abs_tol(aTol),
        abs_grad_tol(aGradTol),
        eta(aEta),
        solve_step(aSolveStep),
        impose_limits(aImposeLimits) {}
  /**
   * This function finds the minimum of a function, given its derivative and Hessian,
   * using a newton search direction and using a trust-region approach.
   * \tparam Vector The vector type of the independent variable for the function.
   * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
   */
  template <typename Vector>
  void operator()(Vector& x) const {
    detail::newton_method_tr_impl(f, df, fill_hessian, x, max_radius, max_iter,
                                  solve_step, impose_limits, abs_tol,
                                  abs_grad_tol, eta);
  }

  /**
   * This function remaps the factory to one which will use a regularized solver within the trust-region.
   * You should regularize the matrix only if there are reasons to expect the Hessian to be near-singular.
   * \param tau The initial relative damping factor to regularize the Hessian matrix.
   */
  newton_method_tr_factory<Function, GradFunction, HessianFunction, T,
                           trust_region_solver_dogleg_reg<T>, LimitFunction>
  regularize(const T& tau) const {
    return newton_method_tr_factory<Function, GradFunction, HessianFunction, T,
                                    trust_region_solver_dogleg_reg<T>,
                                    LimitFunction>(
        f, df, fill_hessian, max_radius, max_iter, abs_tol, abs_grad_tol, eta,
        trust_region_solver_dogleg_reg<T>(tau), impose_limits);
  }

  /**
   * This function remaps the factory to one which will use the given solver within the trust-region.
   * \tparam NewTrustRegionSolver A new functor type that can solve for a solution step within a trust-region (see
   * trust_region_solver_dogleg for an example).
   * \param new_solver The functor that can solve for the step to take within the trust-region.
   */
  template <typename NewTrustRegionSolver>
  newton_method_tr_factory<Function, GradFunction, HessianFunction, T,
                           NewTrustRegionSolver, LimitFunction>
  set_tr_solver(NewTrustRegionSolver new_solver) const {
    return newton_method_tr_factory<Function, GradFunction, HessianFunction, T,
                                    NewTrustRegionSolver, LimitFunction>(
        f, df, fill_hessian, max_radius, max_iter, abs_tol, abs_grad_tol, eta,
        new_solver, impose_limits);
  }

  /**
   * This function remaps the factory to one which will use the given limit-function for the search domain.
   * Using a limit-function boils down to a gradient projection method, and thus, for more complex
   * constraints (such as non-linear ones or mix of equality-inequality constraints), please use a
   * constraint optimization method instead (see augmented_lagrangian_methods.hpp for example).
   * \tparam NewLimitFunction A new functor type that can impose limits on a proposed solution step (see
   * no_limit_functor or box_limit_function for examples).
   * \param new_limits The functor that can impose simple limits on the search domain (i.e. using this boils down to a
   * gradient projection method, for more complex constraints please use a constraint optimization method instead).
   */
  template <typename NewLimitFunction>
  newton_method_tr_factory<Function, GradFunction, HessianFunction, T,
                           TrustRegionSolver, NewLimitFunction>
  set_limiter(NewLimitFunction new_limits) const {
    return newton_method_tr_factory<Function, GradFunction, HessianFunction, T,
                                    TrustRegionSolver, NewLimitFunction>(
        f, df, fill_hessian, max_radius, max_iter, abs_tol, abs_grad_tol, eta,
        solve_step, new_limits);
  }
};

/**
 * This function template creates a factory object to construct a Newton-method optimizer routine
 * that uses a trust-region approach.
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
template <typename Function, typename GradFunction, typename HessianFunction,
          typename T>
newton_method_tr_factory<Function, GradFunction, HessianFunction, T>
make_newton_method_tr(Function f, GradFunction df, HessianFunction fill_hessian,
                      T max_radius, unsigned int max_iter = 100,
                      T abs_tol = T(1e-6), T abs_grad_tol = T(1e-6),
                      T eta = T(1e-4)) {
  return newton_method_tr_factory<Function, GradFunction, HessianFunction, T>(
      f, df, fill_hessian, max_radius, max_iter, abs_tol, abs_grad_tol, eta);
}

/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a newton search
 * direction and using a trust-region approach (dogleg solver).
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param max_iter The maximum number of iterations to perform.
 * \param abs_tol The tolerance on the norm of the step size.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template <typename Function, typename GradFunction, typename HessianFunction,
          typename Vector>
void newton_method_tr(
    Function f, GradFunction df, HessianFunction fill_hessian, Vector& x,
    vect_value_type_t<Vector> max_radius, unsigned int max_iter = 100,
    vect_value_type_t<Vector> abs_tol = vect_value_type_t<Vector>(1e-6),
    vect_value_type_t<Vector> abs_grad_tol = vect_value_type_t<Vector>(1e-6)) {

  detail::newton_method_tr_impl(f, df, fill_hessian, x, max_radius, max_iter,
                                trust_region_solver_dogleg(),
                                no_limit_functor(), abs_tol, abs_grad_tol);
}

/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a
 * regularized newton search direction and using a trust-region approach (dogleg solver).
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param max_iter The maximum number of iterations to perform.
 * \param tau The initial relative damping factor to regularize the Hessian matrix.
 * \param abs_tol The tolerance on the norm of the step size.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template <typename Function, typename GradFunction, typename HessianFunction,
          typename Vector>
void reg_newton_method_tr(
    Function f, GradFunction df, HessianFunction fill_hessian, Vector& x,
    vect_value_type_t<Vector> max_radius, unsigned int max_iter = 100,
    vect_value_type_t<Vector> tau = vect_value_type_t<Vector>(1e-3),
    vect_value_type_t<Vector> abs_tol = vect_value_type_t<Vector>(1e-6),
    vect_value_type_t<Vector> abs_grad_tol = vect_value_type_t<Vector>(1e-6)) {

  detail::newton_method_tr_impl(
      f, df, fill_hessian, x, max_radius, max_iter,
      trust_region_solver_dogleg_reg<vect_value_type_t<Vector>>(tau),
      no_limit_functor(), abs_tol, abs_grad_tol);
}

/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a newton search
 * direction, a function to impose simple limits on the search domain and using a trust-region
 * approach (dogleg solver).
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param max_iter The maximum number of iterations to perform.
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying
 * independent vector-space.
 * \param abs_tol The tolerance on the norm of the step size.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template <typename Function, typename GradFunction, typename HessianFunction,
          typename Vector, typename LimitFunction>
void limited_newton_method_tr(
    Function f, GradFunction df, HessianFunction fill_hessian, Vector& x,
    vect_value_type_t<Vector> max_radius, unsigned int max_iter,
    LimitFunction impose_limits,
    vect_value_type_t<Vector> abs_tol = vect_value_type_t<Vector>(1e-6),
    vect_value_type_t<Vector> abs_grad_tol = vect_value_type_t<Vector>(1e-6)) {

  detail::newton_method_tr_impl(f, df, fill_hessian, x, max_radius, max_iter,
                                trust_region_solver_dogleg(), impose_limits,
                                abs_tol, abs_grad_tol);
}

/**
 * This function finds the minimum of a function, given its derivative and Hessian, using a
 * regularized newton search direction, a function to impose simple limits on the search domain
 * and using a trust-region approach (dogleg solver).
 * \test Must create a unit-test for this.
 * \tparam Function The functor type of the function to optimize.
 * \tparam GradFunction The functor type of the gradient of the function to optimize.
 * \tparam HessianFunction The functor type to fill in the Hessian of the function to optimize.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam LimitFunction A functor type that can impose a limit on a proposed step.
 * \param f The function to minimize.
 * \param df The gradient of the function to minimize.
 * \param fill_hessian The functor object that can fill the Hessian symmetric matrix of the function to be optimized.
 * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
 * \param max_radius The maximum trust-region radius to use (i.e. maximum optimization step).
 * \param max_iter The maximum number of iterations to perform.
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying
 * independent vector-space.
 * \param tau The initial relative damping factor to regularize the Hessian matrix.
 * \param abs_tol The tolerance on the norm of the step size.
 * \param abs_grad_tol The tolerance on the norm of the gradient.
 */
template <typename Function, typename GradFunction, typename HessianFunction,
          typename Vector, typename LimitFunction>
void limited_reg_newton_method_tr(
    Function f, GradFunction df, HessianFunction fill_hessian, Vector& x,
    vect_value_type_t<Vector> max_radius, unsigned int max_iter,
    LimitFunction impose_limits,
    vect_value_type_t<Vector> tau = vect_value_type_t<Vector>(1e-3),
    vect_value_type_t<Vector> abs_tol = vect_value_type_t<Vector>(1e-6),
    vect_value_type_t<Vector> abs_grad_tol = vect_value_type_t<Vector>(1e-6)) {

  detail::newton_method_tr_impl(
      f, df, fill_hessian, x, max_radius, max_iter,
      trust_region_solver_dogleg_reg<vect_value_type_t<Vector>>(tau),
      impose_limits, abs_tol, abs_grad_tol);
}

}  // namespace ReaK::optim

#endif  // REAK_MATH_OPTIMIZATION_NEWTON_METHODS_H_

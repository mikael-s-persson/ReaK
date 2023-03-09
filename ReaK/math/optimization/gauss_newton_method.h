/**
 * \file gauss_newton_method.h
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

#ifndef REAK_MATH_OPTIMIZATION_GAUSS_NEWTON_METHOD_H_
#define REAK_MATH_OPTIMIZATION_GAUSS_NEWTON_METHOD_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_num_exceptions.h"
#include "ReaK/math/lin_alg/mat_qr_decomp.h"

#include "ReaK/math/optimization/limit_functions.h"

namespace ReaK::optim {

namespace detail {

template <typename Function, typename JacobianFunction, typename InputVector,
          typename OutputVector, typename LinearLsqSolver,
          typename LimitFunction>
void gauss_newton_nllsq_impl(Function f, JacobianFunction fill_jac,
                             InputVector& x, const OutputVector& y,
                             unsigned int max_iter, LinearLsqSolver lin_solve,
                             LimitFunction impose_limits,
                             vect_value_type_t<InputVector> abs_tol =
                                 vect_value_type_t<InputVector>(1e-6),
                             vect_value_type_t<InputVector> abs_grad_tol =
                                 vect_value_type_t<InputVector>(1e-6)) {
  using ValueType = vect_value_type_t<InputVector>;
  using std::abs;
  using std::sqrt;

  OutputVector y_approx = f(x);
  OutputVector r = y - y_approx;
  mat<ValueType, mat_structure::rectangular> J(y.size(), x.size());
  fill_jac(J, x, y_approx);
  InputVector e = x;
  mat_vect_adaptor<InputVector> e_mat(e);
  lin_solve(J, e_mat, mat_vect_adaptor<OutputVector>(r), abs_grad_tol);
  impose_limits(x, e);
  unsigned int iter = 0;
  while (norm_2(e) > abs_tol) {
    if (++iter > max_iter) {
      throw maximum_iteration(max_iter);
    }
    x += e;
    y_approx = f(x);
    r = y;
    r -= y_approx;
    fill_jac(J, x, y_approx);
    e = r * J;
    lin_solve(J, e_mat, mat_vect_adaptor<OutputVector>(r), abs_grad_tol);
    impose_limits(x, e);
  }
}
}  // namespace detail

/**
 * This functor is a factory class to construct a non-linear least-square optimizer to a vector
 * function, using the Gauss-Newton method. This method performs very well for well-conditioned systems.
 * TEST PASSED
 * \tparam Function The functor type of the function to optimize.
 * \tparam JacobianFunction The functor type of the gradient of the function to optimize.
 * \tparam OutputVector The output vector type of the non-linear function provided.
 * \tparam T The value-type of the field on which the optimization is performed.
 * \tparam LimitFunction A functor type that can impose limits on a proposed solution step (see no_limit_functor or
 * box_limit_function for examples).
 * \tparam LinearSolver A functor type that can solve a linear least-square system (see QR_linlsqsolver for examples).
 */
template <typename Function, typename JacobianFunction, typename OutputVector,
          typename T = double, typename LimitFunction = no_limit_functor,
          typename LinearSolver = QR_linlsqsolver>
struct gauss_newton_nllsq_factory {
  Function f;
  JacobianFunction fill_jac;
  OutputVector y;
  unsigned int max_iter;
  T abs_tol;
  T abs_grad_tol;
  LimitFunction impose_limits;
  LinearSolver lin_solve;

  using self =
      gauss_newton_nllsq_factory<Function, JacobianFunction, OutputVector, T,
                                 LimitFunction, LinearSolver>;

  /**
   * Parametrized constructor of the factory object.
   * \param aF The vector-function to fit.
   * \param aFillJac The jacobian of the vector-function to fit.
   * \param aY The desired output vector of the function (to be matched by the algorithm).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aTol The tolerance on the norm of the step size.
   * \param aGradTol The tolerance on the norm of the gradient.
   * \param aImposeLimits The functor that can impose simple limits on the search domain (i.e. using this boils down to
   * a gradient projection method, for more complex constraints please use a constraint optimization method instead).
   * \param aLinSolver The functor that can solve a linear system.
   */
  gauss_newton_nllsq_factory(Function aF, JacobianFunction aFillJac,
                             OutputVector aY, unsigned int aMaxIter = 100,
                             T aTol = T(1e-6), T aGradTol = T(1e-6),
                             LimitFunction aImposeLimits = LimitFunction(),
                             LinearSolver aLinSolver = LinearSolver())
      : f(aF),
        fill_jac(aFillJac),
        y(aY),
        max_iter(aMaxIter),
        abs_tol(aTol),
        abs_grad_tol(aGradTol),
        impose_limits(aImposeLimits),
        lin_solve(aLinSolver) {}
  /**
   * This function finds the minimum of a function, given its derivative and Hessian,
   * using a newton search direction and using a trust-region approach.
   * \tparam Vector The vector type of the independent variable for the function.
   * \param x The initial guess to the solution, as well as a storage for the result of the algorihm.
   * \throw singularity_error If a step is detected to cause a near-infinite jump.
   * \throw maximum_iteration If the maximum number of iterations has been reached.
   * \throw other Exceptions can originate from the functors.
   */
  template <typename Vector>
  void operator()(Vector& x) const {
    detail::gauss_newton_nllsq_impl(f, fill_jac, x, y, max_iter, lin_solve,
                                    impose_limits, abs_tol, abs_grad_tol);
  }

  /**
   * Sets the initial damping of the Hessian matrix.
   */
  self& set_max_iteration(unsigned int aMaxIter) {
    max_iter = aMaxIter;
    return *this;
  }
  /**
   * Sets the relative tolerance on the norm of the step size.
   */
  self& set_tolerance(T aTol) {
    abs_tol = aTol;
    return *this;
  }
  /**
   * Sets the relative tolerance on the norm of the step size.
   */
  self& set_gradient_tolerance(T aGradTol) {
    abs_grad_tol = aGradTol;
    return *this;
  }

  /**
   * This function remaps the factory to one which will use the given linear system solver.
   * \tparam NewLinearSolver A functor type that can solve a linear system (see Cholesky_linsolver or QR_linlsqsolver
   * for examples).
   * \param new_lin_solver The functor that can solve a linear system.
   */
  template <typename NewLinearSolver>
  gauss_newton_nllsq_factory<Function, JacobianFunction, OutputVector, T,
                             LimitFunction, NewLinearSolver>
  set_lin_solver(NewLinearSolver new_lin_solver) const {
    return gauss_newton_nllsq_factory<Function, JacobianFunction, OutputVector,
                                      T, LimitFunction, NewLinearSolver>(
        f, fill_jac, y, max_iter, abs_tol, abs_grad_tol, impose_limits,
        new_lin_solver);
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
  gauss_newton_nllsq_factory<Function, JacobianFunction, OutputVector,
                             NewLimitFunction, LinearSolver, T>
  set_limiter(NewLimitFunction new_limits) const {
    return gauss_newton_nllsq_factory<Function, JacobianFunction, OutputVector,
                                      T, NewLimitFunction, LinearSolver>(
        f, fill_jac, y, max_iter, abs_tol, abs_grad_tol, new_limits, lin_solve);
  }
};

/**
 * This function template creates a factory object to construct a non-linear least-square
 * optimizer to a vector function. This method performs very well for well-conditioned systems.
 * TEST PASSED
 * \tparam Function The functor type.
 * \tparam JacobianFunction The functor type to fill the Jacobian matrix.
 * \tparam OutputVector The vector type of the dependent variable of the function.
 * \param f The function to match to the fixed output vector.
 * \param fill_jac The Jacobian of the function to minimize.
 * \param y The fixed vector of outputs to be matched in a least-square approximation.
 * \param max_iter The maximum number of iterations to perform.
 * \param tau The initial, relative damping value to damp the Hessian matrix.
 * \param abs_tol The tolerance on the norm of the estimation steps.
 * \param abs_grad_tol The tolerance on the norm of the gradient (or residual vector).
 */
template <typename Function, typename JacobianFunction, typename OutputVector>
gauss_newton_nllsq_factory<Function, JacobianFunction, OutputVector,
                           vect_value_type_t<OutputVector>>
make_gauss_newton_nllsq(Function f, JacobianFunction fill_jac, OutputVector y,
                        unsigned int max_iter = 100,
                        vect_value_type_t<OutputVector> abs_tol =
                            vect_value_type_t<OutputVector>(1E-6),
                        vect_value_type_t<OutputVector> abs_grad_tol =
                            vect_value_type_t<OutputVector>(1E-6)) {
  return gauss_newton_nllsq_factory<Function, JacobianFunction, OutputVector,
                                    vect_value_type_t<OutputVector>>(
      f, fill_jac, y, max_iter, abs_tol, abs_grad_tol);
}

/**
 * This function finds the non-linear least-square solution to a vector function. This method performs very
 * well for well-conditioned systems.
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
 * \param abs_tol The tolerance on the norm of the estimation steps.
 * \param abs_grad_tol The tolerance on the norm of the gradient (or residual vector).
 */
template <typename Function, typename JacobianFunction, typename InputVector,
          typename OutputVector>
void gauss_newton_nllsq(Function f, JacobianFunction fill_jac, InputVector& x,
                        const OutputVector& y, unsigned int max_iter = 100,
                        vect_value_type_t<InputVector> abs_tol =
                            vect_value_type_t<InputVector>(1e-6),
                        vect_value_type_t<InputVector> abs_grad_tol =
                            vect_value_type_t<InputVector>(1e-6)) {
  detail::gauss_newton_nllsq_impl(f, fill_jac, x, y, max_iter,
                                  QR_linlsqsolver(), no_limit_functor(),
                                  abs_tol, abs_grad_tol);
}

/**
 * This function finds the non-linear least-square solution to a vector function with
 * limits imposed on the search domain. This method performs very
 * well for well-conditioned systems.
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
 * \param max_iter The maximum number of iterations to perform.
 * \param impose_limits The functor to use to limit the proposed steps to satisfy some constraints on the underlying
 * independent vector-space.
 * \param abs_tol The tolerance on the norm of the estimation steps.
 * \param abs_grad_tol The tolerance on the norm of the gradient (or residual vector).
 */
template <typename Function, typename JacobianFunction, typename InputVector,
          typename OutputVector, typename LimitFunction>
void limited_gauss_newton_nllsq(Function f, JacobianFunction fill_jac,
                                InputVector& x, const OutputVector& y,
                                unsigned int max_iter,
                                LimitFunction impose_limits,
                                vect_value_type_t<InputVector> abs_tol =
                                    vect_value_type_t<InputVector>(1e-6),
                                vect_value_type_t<InputVector> abs_grad_tol =
                                    vect_value_type_t<InputVector>(1e-6)) {
  detail::gauss_newton_nllsq_impl(f, fill_jac, x, y, max_iter,
                                  QR_linlsqsolver(), impose_limits, abs_tol,
                                  abs_grad_tol);
}

}  // namespace ReaK::optim

#endif  // REAK_MATH_OPTIMIZATION_GAUSS_NEWTON_METHOD_H_

/**
 * \file function_types.h
 *
 * This library declares a number of classes related to creating object-oriented cost function
 * evaluators.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2012
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

#ifndef REAK_MATH_OPTIMIZATION_FUNCTION_TYPES_H_
#define REAK_MATH_OPTIMIZATION_FUNCTION_TYPES_H_

#include <utility>
#include "ReaK/core/base/shared_object.h"
#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/vect_alg.h"

namespace ReaK::optim {

/**
 * This base-class defines the interface for a cost evaluator which lumps the evaluation of
 * the cost value, gradient, and hessian into one interface.
 */
class cost_evaluator : public shared_object {
 public:
  /**
   * This function computes the cost value for a given state vector.
   * \param x The state vector at which to evaluate the cost value.
   * \return The cost value for the given state vector.
   */
  virtual double compute_cost(const vect_n<double>& x) const = 0;

  /**
   * This function computes the cost gradient for a given state vector.
   * \param x The state vector at which to evaluate the cost gradient.
   * \return The cost gradient for the given state vector.
   */
  virtual vect_n<double> compute_cost_grad(const vect_n<double>& x) const = 0;

  /**
   * This function computes the cost hessian for a given state vector.
   * \param H Stores, as output, the hessian of the cost function at the given state vector.
   * \param x The state vector at which to evaluate the cost hessian.
   * \param f The cost value vector at the given state vector.
   * \param x_grad The cost gradient vector at the given state vector.
   */
  virtual void compute_cost_hessian(mat<double, mat_structure::symmetric>& H,
                                    const vect_n<double>& x, double f,
                                    const vect_n<double>& x_grad) const = 0;

  /**
   * Default destructor.
   */
  ~cost_evaluator() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(cost_evaluator, 0xC1500001, 1, "cost_evaluator",
                              shared_object)
};

/**
 * This class creates a quadratic cost evaluator which lumps the evaluation of
 * the cost value, gradient, and hessian for a quadratic function defined by a center
 * vector and a symmetric matrix.
 */
class quadratic_cost_evaluator : public cost_evaluator {
 public:
  /// Holds the center vector of the quadratic cost function.
  vect_n<double> c;
  /// Holds the positive definite matrix of the quadratic cost function.
  mat<double, mat_structure::symmetric> Q;

  /**
   * Parametrized Constructor.
   * \param aC The center vector of the quadratic cost function.
   * \param aQ The positive definite matrix of the quadratic cost function.
   */
  explicit quadratic_cost_evaluator(
      const vect_n<double>& aC = vect_n<double>(),
      const mat<double, mat_structure::symmetric>& aQ =
          (mat<double, mat_structure::symmetric>()))
      : c(aC), Q(aQ) {}

  /**
   * Default destructor.
   */
  ~quadratic_cost_evaluator() override = default;

  double compute_cost(const vect_n<double>& x) const override {
    vect_n<double> tmp = x;
    tmp -= c;
    return 0.5 * (tmp * (Q * tmp));
  }

  vect_n<double> compute_cost_grad(const vect_n<double>& x) const override {
    return Q * (x - c);
  }

  void compute_cost_hessian(mat<double, mat_structure::symmetric>& H,
                            const vect_n<double>& /*x*/, double /*f*/,
                            const vect_n<double>& /*x_grad*/) const override {
    H = Q;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    cost_evaluator::save(A,
                         cost_evaluator::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(c) & RK_SERIAL_SAVE_WITH_NAME(Q);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    cost_evaluator::load(A,
                         cost_evaluator::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(c) & RK_SERIAL_LOAD_WITH_NAME(Q);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(quadratic_cost_evaluator, 0xC1500002, 1,
                              "quadratic_cost_evaluator", cost_evaluator)
};

/**
 * This class creates a cost evaluator which adds the evaluation of two different
 * cost values, gradients, and hessians. Given two cost evaluators, this class
 * adds their results together.
 */
class added_cost_evaluator : public cost_evaluator {
 public:
  /// Points to the first cost evaluator.
  std::shared_ptr<cost_evaluator> first_eval;
  /// Points to the second cost evaluator.
  std::shared_ptr<cost_evaluator> second_eval;

  /**
   * Parametrized Constructor.
   * \param aFirstEval The first cost evaluator.
   * \param aSecondEval The second cost evaluator.
   */
  explicit added_cost_evaluator(std::shared_ptr<cost_evaluator> aFirstEval =
                                    std::shared_ptr<cost_evaluator>(),
                                std::shared_ptr<cost_evaluator> aSecondEval =
                                    std::shared_ptr<cost_evaluator>())
      : first_eval(std::move(aFirstEval)),
        second_eval(std::move(aSecondEval)) {}

  /**
   * Default destructor.
   */
  ~added_cost_evaluator() override = default;

  double compute_cost(const vect_n<double>& x) const override {
    return first_eval->compute_cost(x) + second_eval->compute_cost(x);
  }

  vect_n<double> compute_cost_grad(const vect_n<double>& x) const override {
    return first_eval->compute_cost_grad(x) + second_eval->compute_cost_grad(x);
  }

  void compute_cost_hessian(mat<double, mat_structure::symmetric>& H,
                            const vect_n<double>& x, double f,
                            const vect_n<double>& x_grad) const override {
    first_eval->compute_cost_hessian(H, x, f, x_grad);
    mat<double, mat_structure::symmetric> H_tmp(H.get_row_count());
    second_eval->compute_cost_hessian(H_tmp, x, f, x_grad);
    H += H_tmp;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    cost_evaluator::save(A,
                         cost_evaluator::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(first_eval) &
        RK_SERIAL_SAVE_WITH_NAME(second_eval);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    cost_evaluator::load(A,
                         cost_evaluator::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(first_eval) &
        RK_SERIAL_LOAD_WITH_NAME(second_eval);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(added_cost_evaluator, 0xC1500003, 1,
                              "added_cost_evaluator", cost_evaluator)
};

/**
 * This function wraps an object-oriented cost evaluator pointer and is a function object
 * that can be passed as a cost-function to the optimization algorithms.
 */
struct oop_cost_function {
  /// Points to the OOP cost-evaluator.
  std::shared_ptr<const cost_evaluator> parent;

  /**
   * Parametrized Constructor.
   * \param aParent A pointer to an OOP cost-evaluator object.
   */
  explicit oop_cost_function(std::shared_ptr<const cost_evaluator> aParent)
      : parent(std::move(aParent)) {}

  /**
   * Evaluates the cost function at x.
   * \param x The state-vector at which to evaluate the cost.
   * \return The cost function evaluated at x.
   */
  double operator()(const vect_n<double>& x) const {
    return parent->compute_cost(x);
  }
};

/**
 * This function wraps an object-oriented cost evaluator pointer and is a function object
 * that can be passed as a cost-gradient to the optimization algorithms.
 */
struct oop_cost_grad {
  /// Points to the OOP cost-evaluator.
  std::shared_ptr<const cost_evaluator> parent;

  /**
   * Parametrized Constructor.
   * \param aParent A pointer to an OOP cost-evaluator object.
   */
  explicit oop_cost_grad(std::shared_ptr<const cost_evaluator> aParent)
      : parent(std::move(aParent)) {}

  /**
   * Evaluates the cost gradient at x.
   * \param x The state-vector at which to evaluate the cost gradient.
   * \return The cost gradient evaluated at x.
   */
  vect_n<double> operator()(const vect_n<double>& x) const {
    return parent->compute_cost_grad(x);
  }
};

/**
 * This function wraps an object-oriented cost evaluator pointer and is a function object
 * that can be passed as a cost-hessian to the optimization algorithms.
 */
struct oop_cost_hess {
  /// Points to the OOP cost-evaluator.
  std::shared_ptr<const cost_evaluator> parent;

  /**
   * Parametrized Constructor.
   * \param aParent A pointer to an OOP cost-evaluator object.
   */
  explicit oop_cost_hess(std::shared_ptr<const cost_evaluator> aParent)
      : parent(std::move(aParent)) {}

  /**
   * Evaluates the cost Hessian at x.
   * \param H The Hessian matrix in which to store the output.
   * \param x The state-vector at which to evaluate the cost Hessian.
   * \param f The cost value at x.
   * \param x_grad The cost gradient at x.
   */
  void operator()(mat<double, mat_structure::symmetric>& H,
                  const vect_n<double>& x, double f,
                  const vect_n<double>& x_grad) const {
    parent->compute_cost_hessian(H, x, f, x_grad);
  }
};

}  // namespace ReaK::optim

#endif  // REAK_MATH_OPTIMIZATION_FUNCTION_TYPES_H_

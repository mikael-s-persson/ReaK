/**
 * \file function_types.hpp
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

#ifndef REAK_FUNCTION_TYPES_HPP
#define REAK_FUNCTION_TYPES_HPP

#include "base/defs.hpp"

namespace ReaK {

namespace optim {
  
  
  

/**
 * This base-class defines the interface for a cost evaluator which lumps the evaluation of 
 * the cost value, gradient, and hessian into one interface. 
 */
class cost_evaluator {
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
    virtual void compute_cost_hessian(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double f, const vect_n<double>& x_grad) const = 0;
    
    /**
     * Default destructor.
     */
    virtual ~cost_evaluator() { };
};



/**
 * This class creates a quadratic cost evaluator which lumps the evaluation of 
 * the cost value, gradient, and hessian for a quadratic function defined by a center
 * vector and a symmetric matrix. 
 */
class quadratic_cost_evaluator : public cost_evaluator {
  public:
    vect_n<double> c; ///< Holds the center vector of the quadratic cost function.
    mat<double,mat_structure::symmetric> Q; ///< Holds the positive definite matrix of the quadratic cost function.
    
    /**
     * Parametrized Constructor.
     * \param aC The center vector of the quadratic cost function.
     * \param aQ The positive definite matrix of the quadratic cost function.
     */
    quadratic_cost_evaluator(const vect_n<double>& aC, const mat<double,mat_structure::symmetric>& aQ) : c(aC), Q(aQ) { };
    
    virtual double compute_cost(const vect_n<double>& x) const {
      vect_n<double> tmp = x; tmp -= c;
      return 0.5 * (tmp * (Q * tmp));
    };
    
    virtual vect_n<double> compute_cost_grad(const vect_n<double>& x) const {
      return Q * (x - c);
    };
    
    virtual void compute_cost_hessian(mat<double,mat_structure::symmetric>& H, const vect_n<double>&, double, const vect_n<double>&) const {
      H = Q;
    };
    
};


/**
 * This class creates a cost evaluator which adds the evaluation of two different 
 * cost values, gradients, and hessians. Given two cost evaluators, this class 
 * adds their results together. 
 */
class added_cost_evaluator : public cost_evaluator {
  public:
    shared_ptr<const cost_evaluator> first_eval; ///< Points to the first cost evaluator.
    shared_ptr<const cost_evaluator> second_eval; ///< Points to the second cost evaluator.
    
    /**
     * Parametrized Constructor.
     * \param aFirstEval The first cost evaluator.
     * \param aSecondEval The second cost evaluator.
     */
    added_cost_evaluator(const shared_ptr<const cost_evaluator>& aFirstEval,
			 const shared_ptr<const cost_evaluator>& aSecondEval) :
			 first_eval(aFirstEval),
			 second_eval(aSecondEval) { };
    
    virtual double compute_cost(const vect_n<double>& x) const {
      return first_eval->compute_cost(x) + second_eval->compute_cost(x);
    };
    
    virtual vect_n<double> compute_cost_grad(const vect_n<double>& x) const {
      return first_eval->compute_cost_grad(x) + second_eval->compute_cost_grad(x);
    };
    
    virtual void compute_cost_hessian(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double f, const vect_n<double>& x_grad) const {
      first_eval->compute_cost_hessian(H, x, f, x_grad);
      mat<double,mat_structure::symmetric> H_tmp(H.get_row_count());
      second_eval->compute_cost_hessian(H_tmp, x, f, x_grad);
      H += H_tmp;
    };
};


/**
 * This function wraps an object-oriented cost evaluator pointer and is a function object 
 * that can be passed as a cost-function to the optimization algorithms.
 */
struct oop_cost_function {
  shared_ptr<const cost_evaluator> parent; ///< Points to the OOP cost-evaluator.
      
  /**
   * Parametrized Constructor.
   * \param aParent A pointer to an OOP cost-evaluator object.
   */
  oop_cost_function(const shared_ptr<const cost_evaluator>& aParent) : parent(aParent) { };
  
  /**
   * Evaluates the cost function at x.
   * \param x The state-vector at which to evaluate the cost.
   * \return The cost function evaluated at x.
   */
  double operator()(const vect_n<double>& x) const {
    return parent->compute_cost(x);
  };
};
    
/**
 * This function wraps an object-oriented cost evaluator pointer and is a function object 
 * that can be passed as a cost-gradient to the optimization algorithms.
 */
struct oop_cost_grad {
  shared_ptr<const cost_evaluator> parent; ///< Points to the OOP cost-evaluator.
      
  /**
   * Parametrized Constructor.
   * \param aParent A pointer to an OOP cost-evaluator object.
   */
  oop_cost_grad(const shared_ptr<const cost_evaluator>& aParent) : parent(aParent) { };
  
  /**
   * Evaluates the cost gradient at x.
   * \param x The state-vector at which to evaluate the cost gradient.
   * \return The cost gradient evaluated at x.
   */
  vect_n<double> operator()(const vect_n<double>& x) const {
    return parent->compute_cost_grad(x);
  };
};
    
/**
 * This function wraps an object-oriented cost evaluator pointer and is a function object 
 * that can be passed as a cost-hessian to the optimization algorithms.
 */
struct oop_cost_hess {
  shared_ptr<const cost_evaluator> parent; ///< Points to the OOP cost-evaluator.
      
  /**
   * Parametrized Constructor.
   * \param aParent A pointer to an OOP cost-evaluator object.
   */
  oop_cost_hess(const shared_ptr<const cost_evaluator>& aParent) : parent(aParent) { };
  
  /**
   * Evaluates the cost Hessian at x.
   * \param H The Hessian matrix in which to store the output.
   * \param x The state-vector at which to evaluate the cost Hessian.
   * \param f The cost value at x.
   * \param x_grad The cost gradient at x.
   */
  void operator()(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double f, const vect_n<double>& x_grad) const {
    parent->compute_cost_hessian(H,x,f,x_grad);
  };
};





};


};

#endif








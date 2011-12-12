/**
 * \file newton_search_directions.hpp
 *
 * The following library provides methods to compute Newton search directions.
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

#ifndef REAK_NEWTON_SEARCH_DIRECTIONS_HPP
#define REAK_NEWTON_SEARCH_DIRECTIONS_HPP

#include "base/defs.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_num_exceptions.hpp"

#include "lin_alg/mat_cholesky.hpp"
#include "lin_alg/mat_damped_matrix.hpp"

#include "lin_alg/mat_norms.hpp"

namespace ReaK {
  
  
namespace optim {



/**
 * This function finds the regularized newton search direction from a given Hessian matrix and gradient vector.
 * \test Must create a unit-test for this.
 * \tparam Matrix The Hessian type.
 * \tparam Vector The vector type of the independent variable for the function.
 * \tparam ScalarMatrix The vector type of the independent variable for the function.
 * \param H The Hessian symmetric matrix of the function to be optimized.
 * \param x_grad The gradient of the function being optimized.
 * \param p The resulting search direction.
 * \param mu The damping scalar matrix to damp the Hessian matrix.
 * \param nu The damping multiplication factor.
 * \param abs_tol The tolerance on the singularity-detection.
 */
template <typename Matrix, typename Vector, typename ScalarMatrix>
void regularized_newton_direction(const Matrix& H, const Vector& x_grad, Vector& p, 
				  ScalarMatrix& mu, typename mat_traits<ScalarMatrix>::value_type& nu,
				  const typename mat_traits<Matrix>::value_type& abs_tol) {
  while(true) {
    try {
      p = -x_grad;
      mat_vect_adaptor<Vector> p_mat(p);
      linsolve_Cholesky(make_damped_matrix(H,mu),p_mat,abs_tol);
      mu *= typename mat_traits<ScalarMatrix>::value_type(0.33333);
      nu = 2.0;
      break;
    } catch (singularity_error&) {
      mu *= nu;
      nu *= 2.0;
    };
  };
};


/**
 * This functor can be used to find the regularized newton search direction from a 
 * given Hessian matrix and gradient vector.
 * \test Must create a unit-test for this.
 * \tparam T The value-type.
 */
template <typename T>
struct regularized_newton_directioner {
  mutable mat<T,mat_structure::scalar> mu;
  mutable T nu;
  /**
   * Parametrized Constructor.
   * \param aTau The initial relative damping factor for the damping the Hessian matrix (the actual damping factor is relative to the trace of the Hessian).
   */
  regularized_newton_directioner(T aTau = T(1e-3)) : mu(1,(aTau <= T(0.0) ? T(1e-3) : aTau)), nu(-1.0) { };
  
  /**
   * This function finds the search direction from a given Hessian matrix and gradient vector.
   * \tparam Matrix The Hessian type.
   * \tparam Vector The vector type of the independent variable for the function.
   * \param H The Hessian symmetric matrix of the function to be optimized.
   * \param x_grad The gradient of the function being optimized.
   * \param p The resulting search direction.
   * \param abs_tol The tolerance on the singularity-detection.
   */
  template <typename Matrix, typename Vector>
  void operator()(const Matrix& H, const Vector& x_grad, Vector& p, const T& abs_tol) const {
    if(nu < T(0.0)) {
      mu = mat<T,mat_structure::scalar>(H.get_row_count(),mu(0,0) * frobenius_norm(H));
      nu = 2.0;
    };
    regularized_newton_direction(H,x_grad,p,mu,nu,abs_tol);
  };
};

/**
 * This function finds the newton search direction from a given Hessian matrix and gradient vector.
 * \test Must create a unit-test for this.
 * \tparam Matrix The Hessian type.
 * \tparam Vector The vector type of the independent variable for the function.
 * \param H The Hessian symmetric matrix of the function to be optimized.
 * \param x_grad The gradient of the function being optimized.
 * \param p The resulting search direction.
 * \param abs_tol The tolerance on the singularity-detection.
 */
template <typename Matrix, typename Vector>
void newton_direction(const Matrix& H, const Vector& x_grad, Vector& p, const typename mat_traits<Matrix>::value_type& abs_tol) {
  p = -x_grad;
  mat_vect_adaptor<Vector> p_mat(p);
  linsolve_Cholesky(H,p_mat,abs_tol);
};

/**
 * This functor can be used to find the newton search direction from a 
 * given Hessian matrix and gradient vector.
 * \test Must create a unit-test for this.
 */
struct newton_directioner {
  /**
   * This function finds the search direction from a given Hessian matrix and gradient vector.
   * \tparam Matrix The Hessian type.
   * \tparam Vector The vector type of the independent variable for the function.
   * \param H The Hessian symmetric matrix of the function to be optimized.
   * \param x_grad The gradient of the function being optimized.
   * \param p The resulting search direction.
   * \param abs_tol The tolerance on the singularity-detection.
   */
  template <typename Matrix, typename Vector>
  void operator()(const Matrix& H, const Vector& x_grad, Vector& p, const typename mat_traits<Matrix>::value_type& abs_tol) const {
    newton_direction(H,x_grad,p,abs_tol);
  };
};




};

};


#endif


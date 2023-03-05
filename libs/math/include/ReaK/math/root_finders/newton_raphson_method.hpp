/**
 * \file newton_raphson_method.hpp
 *
 * This library provides a root-finder function that uses the Newton-Raphson method.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
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

#ifndef REAK_NEWTON_RAPHSON_METHOD_HPP
#define REAK_NEWTON_RAPHSON_METHOD_HPP

#include <cmath>
#include <limits>

#include "ReaK/math/lin_alg/mat_num_exceptions.hpp"

namespace ReaK {

/*************************************************************************
                        Newton-Raphson Root Finding Method
*************************************************************************/

/**
 * This function template performs a Newton-Raphson search for the root of a function. This assumes
 * that the function is monotonic and has a unique root between the two given bounds.
 * \tparam T A scalar value type of the independent and dependent value of the function.
 * \tparam RootedFunction A unary functor type.
 * \tparam DerivativeFunction A unary functor type.
 * \param x The starting point of the search for the root, also stores, as output, the root of the function.
 * \param f The functor of which the root is sought.
 * \param df The derivative functor.
 * \param tol The tolerance on the absolute value of the root.
 * \param max_iter The maximum number of iterations to do, this serves as a divergence criteria.
 * \throw maximum_iteration If the maximum number of iterations is reached before convergence.
 * \throw singularity_error If a stationary point is reached.
 */
template <typename T, typename RootedFunction, typename DerivativeFunction>
void newton_raphson_method(T& x, RootedFunction f, DerivativeFunction df,
                           const T& tol = std::numeric_limits<T>::epsilon(),
                           std::size_t max_iter = 50) {
  using std::abs;

  T y_value = f(x);
  T yp_value = df(x);
  std::size_t i = 0;

  while (abs(yp_value) > tol) {

    x = x - y_value / yp_value;

    y_value = f(x);

    if (abs(y_value) < tol) {
      return;
    }

    yp_value = df(x);

    if (++i > max_iter) {
      throw maximum_iteration(max_iter);
    }
  }

  throw singularity_error(
      "Newton-Raphson method failed due to a stationary point!");
}

}  // namespace ReaK

#endif

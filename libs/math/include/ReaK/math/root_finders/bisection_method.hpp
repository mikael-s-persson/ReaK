/**
 * \file bisection_method.hpp
 *
 * This library provides a root-finder function that uses a bisection method.
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

#ifndef REAK_BISECTION_METHOD_HPP
#define REAK_BISECTION_METHOD_HPP

#include <cmath>
#include <limits>

namespace ReaK {

/*************************************************************************
                        Bisection Method Root Finding
*************************************************************************/

/**
 * This function template performs a bisection method search for the root of a function. This assumes
 * that the function is monotonic and has a unique root between the two given bounds. The function
 * ends up narrowing down the bounds to an interval of a span of the given tolerance value.
 * \tparam T A scalar value type of the independent and dependent value of the function.
 * \tparam RootedFunction A unary functor type.
 * \param low The first bound of the search for the root, will also store the independent value at which the root was
 * found.
 * \param hi The second bound of the search for the root, will also store the independent value at which the root was
 * found.
 * \param f The functor of which the root is sought.
 * \param tol The tolerance, i.e., the size of the resulting interval containing the root.
 */
template <typename T, typename RootedFunction>
void bisection_method(T& low, T& hi, RootedFunction f,
                      const T& tol = std::numeric_limits<T>::epsilon()) {
  using std::abs;

  T low_value = f(low);
  T hi_value = f(hi);

  if (low_value * hi_value > 0.0) {
    return;
  }

  while (abs(hi - low) > tol) {

    T mid = 0.5 * (hi + low);
    T mid_value = f(mid);

    if (mid_value * hi_value > 0.0) {
      hi = mid;
      hi_value = mid_value;
    } else {
      low = mid;
      low_value = mid_value;
    }
  }
}

}  // namespace ReaK

#endif

/**
 * \file secant_method.hpp
 *
 * This library provides a root-finder function that uses the secant method.
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

#ifndef REAK_SECANT_METHOD_HPP
#define REAK_SECANT_METHOD_HPP

#include <algorithm>
#include <cmath>
#include <limits>

namespace ReaK {

/*************************************************************************
                        Secant Root Finding Method
*************************************************************************/

/**
 * This function template performs a secant method search for the root of a function.
 * \tparam T A scalar value type of the independent and dependent value of the function.
 * \tparam RootedFunction A unary functor type.
 * \param low_bound The lower bound of the search for the root.
 * \param hi_bound The upper bound of the search for the root.
 * \param f The functor of which the root is sought.
 * \param tol The tolerance on the absolute value of the root.
 * \return The independent variable value at which the root was found, or the bound value that was reached while seeking
 * the root (and diverging).
 */
template <typename T, typename RootedFunction>
T secant_method(const T& low_bound, const T& hi_bound, RootedFunction f,
                const T& tol = std::numeric_limits<T>::epsilon()) {
  using std::abs;

  T x0 = low_bound;
  T x1 = hi_bound;
  T x0_value = f(x0);
  T x1_value = f(x1);

  if (x0_value * x1_value > 0.0) {
    return (abs(x0_value) < abs(x1_value) ? x0 : x1);
  }
  T abs_tol = tol * abs(hi_bound - low_bound);
  T abs_f_tol = tol * (abs(x0_value) + abs(x1_value));
  T cur_interval_size = abs(x0 - x1);

  while (cur_interval_size > abs_tol) {

    //     T x2 = x1 - x1_value * (x1 - x0) / (x1_value - x0_value); // <--- original formula
    // this one is more symmetric and thus, numerically better conditioned:
    T x2 = x0 / (T(1.0) - x0_value / x1_value) +
           x1 / (T(1.0) - x1_value / x0_value);

    // NOTE: added this introspective test to break away from a pathological problem (revert to bisection):
    //  If the calculated new point is too close to an existing bound, then pick the middle-of-interval instead.
    if ((abs(x2 - x0) < 0.05 * cur_interval_size) ||
        (abs(x2 - x1) < 0.05 * cur_interval_size)) {
      x2 = (x0 + x1) * T(0.5);
    }

    T x2_value = f(x2);

    if (abs(x2_value) < abs_f_tol) {
      return x2;
    }

    if (x2_value * x1_value > T(0.0)) {
      // if x2 has the same sign as x1, then we retain x0 and x2:
      x1 = x2;
      x1_value = x2_value;
    } else {
      // else x1 and x2 have opposite signs, then we retain x2 and x1:
      x0 = x2;
      x0_value = x2_value;
    }

    cur_interval_size = abs(x0 - x1);
  }

  return x1;
}

/**
 * This function template performs the Illinois algorithm for the root of a function. This assumes
 * that the function is monotonic and has a unique root between the two given bounds.
 * \tparam T A scalar value type of the independent and dependent value of the function.
 * \tparam RootedFunction A unary functor type.
 * \param low_bound The lower bound of the search for the root, and store as output the lower-bound of the narrowed down
 * interval in which the root exists.
 * \param hi_bound The upper bound of the search for the root, and store as output the upper-bound of the narrowed down
 * interval in which the root exists.
 * \param f The functor of which the root is sought.
 * \param tol The tolerance on the size of the narrowed-down intervale in which the root exists.
 */
template <typename T, typename RootedFunction>
void illinois_method(T& low_bound, T& hi_bound, RootedFunction f,
                     const T& tol = std::numeric_limits<T>::epsilon()) {
  using std::abs;

  T x0_value = f(low_bound);
  T x1_value = f(hi_bound);

  if (x0_value * x1_value > 0.0) {
    return;
  }

  int last_retained_bound = 0;
  T abs_tol = tol * abs(hi_bound - low_bound);
  T abs_f_tol = tol * (abs(x0_value) + abs(x1_value));

  while (abs(low_bound - hi_bound) > abs_tol) {

    //     T x2 = hi_bound - x1_value * (hi_bound - low_bound) / (x1_value - x0_value); // <--- original formula
    // this one is more symmetric and thus, numerically better conditioned:
    T x2 = low_bound / (T(1.0) - x0_value / x1_value) +
           hi_bound / (T(1.0) - x1_value / x0_value);

    T x2_value = f(x2);

    if (x2_value * x1_value > 0.0) {
      hi_bound = x2;
      if (abs(x2_value) < abs_f_tol) {
        low_bound = hi_bound;
        return;
      }
      if (last_retained_bound == -1) {
        x0_value *= 0.5;
      }
      x1_value = x2_value;
      last_retained_bound = -1;
    } else {
      low_bound = x2;
      if (abs(x2_value) < abs_f_tol) {
        hi_bound = low_bound;
        return;
      }
      if (last_retained_bound == 1) {
        x1_value *= 0.5;
      }
      x0_value = x2_value;
      last_retained_bound = 1;
    }
  }
}

/**
 * This function template performs the Ford-3 algorithm for the root of a function. This assumes
 * that the function is monotonic and has a unique root between the two given bounds. The Ford-3
 * method is a Illinois-type method with a refined choice of scaling factor that acheives super-linear
 * convergence.
 * \tparam T A scalar value type of the independent and dependent value of the function.
 * \tparam RootedFunction A unary functor type.
 * \param low_bound The lower bound of the search for the root, and store as output the lower-bound of the narrowed down
 * interval in which the root exists.
 * \param hi_bound The upper bound of the search for the root, and store as output the upper-bound of the narrowed down
 * interval in which the root exists.
 * \param f The functor of which the root is sought.
 * \param tol The tolerance on the size of the narrowed-down intervale in which the root exists.
 */
template <typename T, typename RootedFunction>
void ford3_method(T& low_bound, T& hi_bound, RootedFunction f,
                  const T& tol = std::numeric_limits<T>::epsilon()) {
  using std::abs;

  T x0_value = f(low_bound);
  T x1_value = f(hi_bound);

  if (x0_value * x1_value > 0.0) {
    return;
  }

  int last_retained_bound = 0;
  T abs_tol = tol * abs(hi_bound - low_bound);
  T abs_f_tol = tol * (abs(x0_value) + abs(x1_value));

  while (abs(low_bound - hi_bound) > abs_tol) {

    //     T x2 = hi_bound - x1_value * (hi_bound - low_bound) / (x1_value - x0_value); // <--- original formula
    // this one is more symmetric and thus, numerically better conditioned:
    T x2 = low_bound / (T(1.0) - x0_value / x1_value) +
           hi_bound / (T(1.0) - x1_value / x0_value);

    T x2_value = f(x2);

    if (x2_value * x1_value > 0.0) {
      hi_bound = x2;
      if (abs(x2_value) < abs_f_tol) {
        low_bound = hi_bound;
        return;
      }
      if (last_retained_bound == -1) {
        x0_value *= 1.0 - x2_value / (x1_value * (1.0 - x2_value / x0_value));
      }
      x1_value = x2_value;
      last_retained_bound = -1;
    } else {
      low_bound = x2;
      if (abs(x2_value) < abs_f_tol) {
        hi_bound = low_bound;
        return;
      }
      if (last_retained_bound == 1) {
        x1_value *= 1.0 - x2_value / (x0_value * (1.0 - x2_value / x1_value));
      }
      x0_value = x2_value;
      last_retained_bound = 1;
    }
  }
}

/**
 * This function template performs Brent's method for the root of a function. This assumes
 * that the function is monotonic and has a unique root between the two given bounds.
 * \tparam T A scalar value type of the independent and dependent value of the function.
 * \tparam RootedFunction A unary functor type.
 * \param a The lower bound of the search for the root, and store as output the lower-bound of the narrowed down
 * interval in which the root exists.
 * \param b The upper bound of the search for the root, and store as output the upper-bound of the narrowed down
 * interval in which the root exists.
 * \param f The functor of which the root is sought.
 * \param tol The tolerance on the size of the narrowed-down intervale in which the root exists.
 */
template <typename T, typename RootedFunction>
void brent_method(T& a, T& b, RootedFunction f,
                  const T& tol = std::numeric_limits<T>::epsilon()) {
  using std::abs;
  using std::swap;

  T a_value = f(a);
  T b_value = f(b);

  if (a_value * b_value > 0.0) {
    return;
  }
  if (abs(a_value) < abs(b_value)) {
    swap(a, b);
    swap(a_value, b_value);
  }
  T c = a;
  T c_value = a_value;
  T d = c;

  bool flag = true;
  T abs_tol = tol * abs(a - b);
  T abs_f_tol = tol * abs(a_value - b_value);

  while (abs(b - a) > abs_tol) {
    T s = a;
    if ((abs(a_value - c_value) > tol * abs(a_value)) &&
        (abs(b_value - c_value) > tol * abs(a_value))) {
      s = a * b_value * c_value / ((a_value - b_value) * (a_value - c_value)) +
          b * a_value * c_value / ((b_value - a_value) * (b_value - c_value)) +
          c * a_value * b_value / ((c_value - a_value) * (c_value - b_value));
    } else {
      s = b - b_value * (b - a) / (b_value - a_value);
    }

    T a3_b = 0.25 * (3.0 * a + b);
    if (((abs(s - a3_b) > abs(b - a3_b)) || (abs(s - b) > abs(b - a3_b))) ||
        (flag && (abs(s - b) >= abs(b - c) * 0.5)) ||
        (!flag && (abs(s - b) >= abs(c - d) * 0.5)) ||
        (flag && (abs(b - c) < abs_tol)) || (!flag && (abs(c - d) < abs_tol))) {
      s = 0.5 * (a + b);
      flag = true;
    } else {
      flag = false;
    }

    T s_value = f(s);
    if (abs(s_value) < abs_f_tol) {
      a = b = s;
      return;
    }

    d = c;
    c = b;
    c_value = b_value;

    if (a_value * s_value > 0.0) {
      a = s;
      a_value = s_value;
    } else {
      b = s;
      b_value = s_value;
    }
    if (abs(a_value) < abs(b_value)) {
      swap(a, b);
      swap(a_value, b_value);
    }
  }
}

/**
 * This function template performs Ridder's method for the root of a function. This assumes
 * that the function is monotonic and has a unique root between the two given bounds.
 * \tparam T A scalar value type of the independent and dependent value of the function.
 * \tparam RootedFunction A unary functor type.
 * \param low_bound The lower bound of the search for the root, and store as output the lower-bound of the narrowed down
 * interval in which the root exists.
 * \param hi_bound The upper bound of the search for the root, and store as output the upper-bound of the narrowed down
 * interval in which the root exists.
 * \param f The functor of which the root is sought.
 * \param tol The tolerance on the size of the narrowed-down intervale in which the root exists.
 */
template <typename T, typename RootedFunction>
void ridders_method(T& low_bound, T& hi_bound, RootedFunction f,
                    const T& tol = std::numeric_limits<T>::epsilon()) {
  using std::abs;
  using std::sqrt;

  T x0_value = f(low_bound);
  T x1_value = f(hi_bound);

  if (x0_value * x1_value > 0.0) {
    return;
  }

  T abs_tol = tol * abs(hi_bound - low_bound);
  T abs_f_tol = tol * abs(x0_value - x1_value);

  while (abs(low_bound - hi_bound) > abs_tol) {

    T x2 = 0.5 * (hi_bound + low_bound);
    T x2_value = f(x2);
    if (abs(x2_value) < abs_f_tol) {
      low_bound = hi_bound = x2;
      return;
    }

    T x3 = x2 + (x2 - hi_bound) * (x0_value > 0.0 ? x2_value : -x2_value) /
                    sqrt(x2_value * x2_value - x0_value * x1_value);
    T x3_value = f(x3);
    if (abs(x3_value) < abs_f_tol) {
      low_bound = hi_bound = x3;
      return;
    }

    if (x0_value * x2_value > 0.0) {
      if ((x0_value * x3_value > 0.0) &&
          (abs(x2 - low_bound) < abs(x3 - low_bound))) {
        low_bound = x3;
        x0_value = x3_value;
      } else {
        low_bound = x2;
        x0_value = x2_value;
      }
      if (x1_value * x3_value > 0.0) {
        hi_bound = x3;
        x1_value = x3_value;
      }
    } else {
      if ((x1_value * x3_value > 0.0) &&
          (abs(x2 - hi_bound) < abs(x3 - hi_bound))) {
        hi_bound = x3;
        x1_value = x3_value;
      } else {
        hi_bound = x2;
        x1_value = x2_value;
      }
      if (x0_value * x3_value > 0.0) {
        low_bound = x3;
        x0_value = x3_value;
      }
    }
  }
}

}  // namespace ReaK

#endif

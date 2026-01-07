/**
 * \file line_search.h
 *
 * The following library is a collection of line-search optimization algorithms.
 * For unimodal cost functions (monotonic derivative, i.e. a single global optimum),
 * three classic algorithms are found in this library: Dichotomous search, Golden-Section
 * search and Fibonacci search.
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date September 2010
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

#ifndef REAK_MATH_OPTIMIZATION_LINE_SEARCH_H_
#define REAK_MATH_OPTIMIZATION_LINE_SEARCH_H_


#include <array>
#include <cmath>
#include <vector>

namespace ReaK::optim {

namespace detail {

template <typename T, typename Function>
T dichotomous_search_impl(Function f, T& low_bound, T& up_bound, T delta,
                          T tol) {
  using std::abs;
  while (true) {
    T midpoint = (low_bound + up_bound) * T(0.5);
    if (abs(up_bound - low_bound) < tol) {
      return f(midpoint);
    }
    T f1 = f(midpoint - delta);
    T f2 = f(midpoint + delta);
    if (f1 > f2) {
      low_bound = midpoint - delta;
      delta *= T(0.5);
    } else {
      up_bound = midpoint + delta;
      delta *= T(0.5);
    }
  }
}

const double GoldenRatioPhi = 1.618033988;

template <typename T, typename Function>
T golden_section_search_impl(Function f, T& low_bound, T& up_bound, T tol) {
  using std::abs;

  T mid_value = low_bound + (up_bound - low_bound) / GoldenRatioPhi;
  T mid_cost = f(mid_value);

  while (true) {
    if (abs(low_bound - up_bound) < tol) {
      return f((low_bound + up_bound) * T(0.5));
    }

    T test_value = mid_value + (up_bound - mid_value) / GoldenRatioPhi;
    T test_cost = f(test_value);

    if (test_cost < mid_cost) {
      low_bound = mid_value;
      mid_value = test_value;
      mid_cost = test_cost;
    } else {
      up_bound = low_bound;
      low_bound = test_value;
    }
  }
}

inline std::vector<int>::reverse_iterator get_fibonacci_iter_for_tolerance(
    double low_bound, double up_bound, double tol) {
  using std::abs;
  static std::array<int, 2> first_fib = {0, 1};
  static std::vector<int> fib_seq(first_fib.begin(), first_fib.end());
  while (abs(up_bound - low_bound) > fib_seq.back() * tol) {
    fib_seq.push_back(fib_seq.back() + (*(fib_seq.rbegin() + 1)));
  }
  return fib_seq.rbegin();
}

template <typename T, typename Function>
T fibonacci_search_impl(Function f, T& low_bound, T& up_bound, T tol) {
  auto fib_iter = get_fibonacci_iter_for_tolerance(low_bound, up_bound, tol);

  T Lstar2 = (up_bound - low_bound) * T(*(fib_iter + 2)) / T(*fib_iter);
  T x1 = low_bound + Lstar2;
  T x1_value = f(x1);
  T x2 = up_bound - Lstar2;
  T x2_value = f(x2);
  while (true) {
    ++fib_iter;
    if ((*fib_iter) == 2) {
      return x1_value;
    }
    if (x1_value > x2_value) {
      low_bound = x1;
      x1 = x2;
      x1_value = x2_value;
      x2 =
          up_bound - T(*(fib_iter + 2)) / T(*fib_iter) * (up_bound - low_bound);
      x2_value = f(x2);
    } else {
      up_bound = x2;
      x2 = x1;
      x2_value = x1_value;
      x1 = low_bound +
           T(*(fib_iter + 2)) / T(*fib_iter) * (up_bound - low_bound);
      x1_value = f(x1);
    }
  }
}

template <typename T, typename Function, typename GradFunction>
T backtracking_search_impl(Function f, GradFunction df, T& low_bound,
                           T& up_bound, T tol, T c1, T c2, T tau) {
  using std::abs;
  if (tau >= 1.0) {
    low_bound = up_bound;
    return f(up_bound);
  }
  T x0_value = f(low_bound);
  T x0_grad = df(low_bound);
  T x1_value = f(up_bound);
  T x1_grad = df(up_bound);
  while ((abs(up_bound - low_bound) > tol) &&
         ((x1_value > x0_value + c1 * (up_bound - low_bound) * x0_grad) ||
          (abs(x1_grad) > c2 * abs(x0_grad)))) {
    up_bound = low_bound + tau * (up_bound - low_bound);
    x1_value = f(up_bound);
    x1_grad = df(up_bound);
  }
  low_bound = up_bound;
  return x1_value;
}

template <typename T, typename Function, typename GradFunction>
void expand_and_zoom_zoom_impl(Function f, GradFunction df, T& low_bound,
                               T x0_value, T x0_grad, T& x1, T& x1_value,
                               T& x1_grad, T& x2, T& x2_value, T& x2_grad,
                               T& up_bound, T tol, T c1, T c2) {
  using std::abs;
  using std::sqrt;

  while (abs(x2 - x1) > tol * abs(up_bound - low_bound)) {
    T d1 = x1_grad + x2_grad - T(3.0) * (x1_value - x2_value) / (x1 - x2);
    T d2 = sqrt(d1 * d1 - x1_grad * x2_grad);
    if (x2 < x1) {
      d2 *= T(-1.0);
    }
    T tmp = x2 -
            (x2 - x1) * (x2_grad + d2 - d1) / (x2_grad - x1_grad + T(2.0) * d2);
    T tmp_value = f(tmp);
    T tmp_grad = df(tmp);
    if ((tmp_value > x0_value + c1 * (tmp - low_bound) * x0_grad) ||
        (tmp_value >= x1_value)) {
      x2 = tmp;
      x2_value = tmp_value;
      x2_grad = tmp_grad;
    } else {
      if (abs(tmp_grad) <= -c2 * x0_grad) {
        x1 = tmp;
        x1_value = tmp_value;
        x1_grad = tmp_grad;
        x2 = tmp;
        x2_value = tmp_value;
        x2_grad = tmp_grad;
        return;
      }
      if (tmp_grad * (x2 - x1) >= T(0.0)) {
        x2 = x1;
        x2_value = x1_value;
        x2_grad = x1_grad;
      }
      x1 = tmp;
      x1_value = tmp_value;
      x1_grad = tmp_grad;
    }
  }
}

template <typename T, typename Function, typename GradFunction>
T expand_and_zoom_search_impl(Function f, GradFunction df, T& low_bound,
                              T& up_bound, T tol, T c1, T c2) {
  using std::abs;
  T x0_value = f(low_bound);
  T x0_grad = df(low_bound);
  T x1 = low_bound;
  T x1_value = x0_value;
  T x1_grad = x0_grad;
  T x2 = (low_bound + up_bound) * 0.5;
  while (true) {
    T x2_value = f(x2);
    T x2_grad = df(x2);
    if ((x2_value > x0_value + c1 * (x2 - low_bound)) ||
        ((x1 != low_bound) && (x2_value >= x1_value))) {
      expand_and_zoom_zoom_impl(f, df, low_bound, x0_value, x0_grad, x1,
                                x1_value, x1_grad, x2, x2_value, x2_grad,
                                up_bound, tol, c1, c2);
      low_bound = x2;
      up_bound = x2;
      return x2_value;
    }
    if (abs(x2_grad) <= -c2 * x0_grad) {
      low_bound = x2;
      up_bound = x2;
      return x2_value;
    }
    if (x2_grad >= T(0.0)) {
      expand_and_zoom_zoom_impl(f, df, low_bound, x0_value, x0_grad,  // NOLINT
                                x2, x2_value, x2_grad, x1, x1_value, x1_grad,
                                up_bound, tol, c1, c2);
      low_bound = x1;
      up_bound = x1;
      return x1_value;
    }
    x1 = x2;
    x1_value = x2_value;
    x1_grad = x2_grad;
    x2 = low_bound + (x2 - low_bound) * 1.105;
  }
}

}  // namespace detail

/**
 * This function performs a 1D optimum Dichotomous search on a unimodal cost function to find the value for
 * the lowest cost (both value and cost are returned by this function as a pair). The search guarantees
 * an uncertainty interval of length tol and limits the search to within (low_bound .. up_bound).
 * TEST PASSED
 * \tparam T The value type for the scalar that represents both the independent variable and the cost values.
 * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
 * \param f The cost functor.
 * \param low_bound The lower bound for the search.
 * \param up_bound the upper bound for the search.
 * \param tol the uncertainty on the X value at which the optimization should stop.
 * \return the cost that is the optimum point found, within uncertainty interval (low_bound, up_bound).
 */
template <typename T, typename Function>
T dichotomous_search(Function f, T& low_bound, T& up_bound, T tol) {
  return detail::dichotomous_search_impl(f, low_bound, up_bound,
                                         T((up_bound - low_bound) * 0.1), tol);
}

/**
 * This functor class wraps a call to a Dichotomous search in order to find the minimum of a
 * vector function across a line-search (i.e. with a search direction). This functor only
 * provides bounded searches (with an upper and lower bound provided).
 * TEST PASSED
 */
struct bounded_line_srch_dichotomous {

  /**
   * This overload computes the scalar value (multiplying the search direction vector) which
   * minimizes the given function along the line defined by the vector offset (x0) and the
   * search direction (dx). The solution will lie between a0 and a1 which are assumed to bound
   * the search domain.
   * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
   * \tparam T A scalar type.
   * \tparam U Any vector type which can be linearly composed (i.e. closed under multiplication by a scalar and under
   * addition).
   * \param f The cost function to minimize (should map U -> T).
   * \param a0 The first bound of the search.
   * \param a1 The second bound of the search.
   * \param x0 The starting point in the vector-space (or offset).
   * \param dx The search direction in the vector-space.
   * \param tol The tolerance at which the search is stopped (i.e. when the interval of convergence is narrow enough).
   * \return The optimal scalar value along the line-search.
   */
  template <typename Function, typename T, typename U>
  T operator()(Function f, T a0, T a1, const U& x0, const U& dx,
               const T& tol = T(1e-6)) const {
    dichotomous_search(
        [&f, &x0, &dx](const T& alpha) -> T { return f(x0 + alpha * dx); }, a0,
        a1, tol);
    return (a1 + a0) * T(0.5);
  }

  /**
   * This overload computes the scalar value (multiplying the search direction vector) which
   * minimizes the given function along the line defined by the vector offset (x0) and the
   * search direction (dx). The solution will lie between a0 and a1 which are assumed to bound
   * the search domain.
   * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
   * \tparam GradFunction A functor type for a unary function that computes the gradient of the cost for a given
   * independent variable value.
   * \tparam T A scalar type.
   * \tparam U Any vector type which can be linearly composed (i.e. closed under multiplication by a scalar and under
   * addition).
   * \param f The cost function to minimize (should map U -> T).
   * \param a0 The first bound of the search.
   * \param a1 The second bound of the search.
   * \param x0 The starting point in the vector-space (or offset).
   * \param dx The search direction in the vector-space.
   * \param tol The tolerance at which the search is stopped (i.e. when the interval of convergence is narrow enough).
   * \return The optimal scalar value along the line-search.
   */
  template <typename Function, typename GradFunction, typename T, typename U>
  T operator()(Function f, GradFunction /*unused*/, T a0, T a1, const U& x0,
               const U& dx, const T& tol = T(1e-6)) const {
    dichotomous_search(
        [&f, &x0, &dx](const T& alpha) -> T { return f(x0 + alpha * dx); }, a0,
        a1, tol);
    return (a1 + a0) * T(0.5);
  }
};

/**
 * This function performs a 1D optimum Golden-Section search on a unimodal cost function to find the value for
 * the lowest cost (both value and cost are returned by this function as a pair). The search guarantees
 * an uncertainty interval of length tol and limits the search to within (bound1 .. bound2) (the bounds
 * do not need to be sorted).
 * TEST PASSED
 * \tparam T The value type for the scalar that represents both the independent variable and the cost values.
 * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
 * \param f the cost functor.
 * \param bound1 the first bound for the search.
 * \param bound2 the second bound for the search.
 * \param tol the uncertainty on the X value at which the optimization should stop.
 * \return the cost that is the optimum point found, within uncertainty interval (bound1, bound2).
 */
template <typename T, typename Function>
T golden_section_search(Function f, T& bound1, T& bound2, T tol) {
  return detail::golden_section_search_impl(f, bound1, bound2, tol);
}

/**
 * This functor class wraps a call to a Golden-Section search in order to find the minimum of a
 * vector function across a line-search (i.e. with a search direction). This functor only
 * provides bounded searches (with an upper and lower bound provided).
 * TEST PASSED
 */
struct bounded_line_srch_gold_sect {

  /**
   * This overload computes the scalar value (multiplying the search direction vector) which
   * minimizes the given function along the line defined by the vector offset (x0) and the
   * search direction (dx). The solution will lie between a0 and a1 which are assumed to bound
   * the search domain.
   * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
   * \tparam T A scalar type.
   * \tparam U Any vector type which can be linearly composed (i.e. closed under multiplication by a scalar and under
   * addition).
   * \param f The cost function to minimize (should map U -> T).
   * \param a0 The first bound of the search.
   * \param a1 The second bound of the search.
   * \param x0 The starting point in the vector-space (or offset).
   * \param dx The search direction in the vector-space.
   * \param tol The tolerance at which the search is stopped (i.e. when the interval of convergence is narrow enough).
   * \return The optimal scalar value along the line-search.
   */
  template <typename Function, typename T, typename U>
  T operator()(Function f, T a0, T a1, const U& x0, const U& dx,
               const T& tol = T(1e-6)) const {
    golden_section_search(
        [&f, &x0, &dx](const T& alpha) -> T { return f(x0 + alpha * dx); }, a0,
        a1, tol);
    return (a1 + a0) * T(0.5);
  }

  /**
   * This overload computes the scalar value (multiplying the search direction vector) which
   * minimizes the given function along the line defined by the vector offset (x0) and the
   * search direction (dx). The solution will lie between a0 and a1 which are assumed to bound
   * the search domain.
   * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
   * \tparam GradFunction A functor type for a unary function that computes the gradient of the cost for a given
   * independent variable value.
   * \tparam T A scalar type.
   * \tparam U Any vector type which can be linearly composed (i.e. closed under multiplication by a scalar and under
   * addition).
   * \param f The cost function to minimize (should map U -> T).
   * \param a0 The first bound of the search.
   * \param a1 The second bound of the search.
   * \param x0 The starting point in the vector-space (or offset).
   * \param dx The search direction in the vector-space.
   * \param tol The tolerance at which the search is stopped (i.e. when the interval of convergence is narrow enough).
   * \return The optimal scalar value along the line-search.
   */
  template <typename Function, typename GradFunction, typename T, typename U>
  T operator()(Function f, GradFunction /*unused*/, T a0, T a1, const U& x0,
               const U& dx, const T& tol = T(1e-6)) const {
    golden_section_search(
        [&f, &x0, &dx](const T& alpha) -> T { return f(x0 + alpha * dx); }, a0,
        a1, tol);
    return (a1 + a0) * T(0.5);
  }
};

/**
 * This function performs a 1D optimum Fibonacci search on a unimodal cost function aFunc to find the value for
 * the lowest cost (both value and cost are returned by this function as a pair). The search guarantees
 * an uncertainty interval of length aToleranceX and limits the search to within (aLowBound .. aUpBound).
 * TEST PASSED
 * \tparam T The value type for the scalar that represents both the independent variable and the cost values.
 * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
 * \param f the cost functor.
 * \param low_bound the lower bound for the search.
 * \param up_bound the upper bound for the search.
 * \param tol the uncertainty on the X value at which the optimization should stop.
 * \return the cost that is the optimum point found, within uncertainty interval (low_bound, up_bound).
 */
template <typename T, typename Function>
T fibonacci_search(Function f, T& low_bound, T& up_bound, T tol) {
  using std::swap;
  if (low_bound > up_bound) {
    swap(low_bound, up_bound);
  }

  return detail::fibonacci_search_impl(f, low_bound, up_bound, tol);
}

/**
 * This functor class wraps a call to a Fibonacci search in order to find the minimum of a
 * vector function across a line-search (i.e. with a search direction). This functor only
 * provides bounded searches (with an upper and lower bound provided).
 * TEST PASSED
 */
struct bounded_line_srch_fibonacci {

  /**
   * This overload computes the scalar value (multiplying the search direction vector) which
   * minimizes the given function along the line defined by the vector offset (x0) and the
   * search direction (dx). The solution will lie between a0 and a1 which are assumed to bound
   * the search domain.
   * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
   * \tparam T A scalar type.
   * \tparam U Any vector type which can be linearly composed (i.e. closed under multiplication by a scalar and under
   * addition).
   * \param f The cost function to minimize (should map U -> T).
   * \param a0 The first bound of the search.
   * \param a1 The second bound of the search.
   * \param x0 The starting point in the vector-space (or offset).
   * \param dx The search direction in the vector-space.
   * \param tol The tolerance at which the search is stopped (i.e. when the interval of convergence is narrow enough).
   * \return The optimal scalar value along the line-search.
   */
  template <typename Function, typename T, typename U>
  T operator()(Function f, T a0, T a1, const U& x0, const U& dx,
               const T& tol = T(1e-6)) const {
    fibonacci_search(
        [&f, &x0, &dx](const T& alpha) -> T { return f(x0 + alpha * dx); }, a0,
        a1, tol);
    return (a1 + a0) * T(0.5);
  }

  /**
   * This overload computes the scalar value (multiplying the search direction vector) which
   * minimizes the given function along the line defined by the vector offset (x0) and the
   * search direction (dx). The solution will lie between a0 and a1 which are assumed to bound
   * the search domain.
   * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
   * \tparam GradFunction A functor type for a unary function that computes the gradient of the cost for a given
   * independent variable value.
   * \tparam T A scalar type.
   * \tparam U Any vector type which can be linearly composed (i.e. closed under multiplication by a scalar and under
   * addition).
   * \param f The cost function to minimize (should map U -> T).
   * \param a0 The first bound of the search.
   * \param a1 The second bound of the search.
   * \param x0 The starting point in the vector-space (or offset).
   * \param dx The search direction in the vector-space.
   * \param tol The tolerance at which the search is stopped (i.e. when the interval of convergence is narrow enough).
   * \return The optimal scalar value along the line-search.
   */
  template <typename Function, typename GradFunction, typename T, typename U>
  T operator()(Function f, GradFunction /*unused*/, T a0, T a1, const U& x0,
               const U& dx, const T& tol = T(1e-6)) const {
    fibonacci_search(
        [&f, &x0, &dx](const T& alpha) -> T { return f(x0 + alpha * dx); }, a0,
        a1, tol);
    return (a1 + a0) * T(0.5);
  }
};

/**
 * This function performs a 1D optimum Backtracking search on a unimodal cost function f to find the value for
 * the lowest cost (both value and cost are returned by this function as a pair). The search guarantees
 * that the strong Wolfe conditions are met at the obtained optimal point.
 * TEST PASSED
 * \tparam T The value type for the scalar that represents both the independent variable and the cost values.
 * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
 * \param f the cost functor.
 * \param low_bound the lower bound for the search.
 * \param up_bound the upper bound for the search.
 * \param tol the uncertainty on the X value at which the optimization should stop.
 * \param c1 The first parameter defining the strong Wolfe conditions (Armijo condition) (should be between 0 and 1
 * exclusively).
 * \param c2 The second parameter defining the strong Wolfe conditions (curvature condition) (should be between 0 and 1
 * exclusively, and greater than c1).
 * \param tau The geometric decrease factor for the successive decrease of the search interval (should be between 0 and
 * 1 exclusively).
 * \return the cost that is the optimum point found, within uncertainty interval (low_bound, up_bound).
 */
template <typename T, typename Function, typename GradFunction>
T backtracking_search(Function f, GradFunction df, T& low_bound, T& up_bound,
                      T tol = T(1e-6), T c1 = T(1e-2), T c2 = T(0.9),
                      T tau = T(0.75)) {
  return detail::backtracking_search_impl(f, df, low_bound, up_bound, tol, c1,
                                          c2, tau);
}

/**
 * This functor class wraps a call to a Backtracking search in order to find the minimum of a
 * vector function across a line-search (i.e. with a search direction). This functor only
 * provides bounded searches (with an upper and lower bound provided) and returns the minimum
 * according to the strong Wolfe conditions.
 * TEST PASSED
 * \tparam T A scalar type.
 */
template <typename T>
struct line_search_backtracking {
  T c1, c2, geom_factor;

  /**
   * Parametrized constructor.
   * \param aC1 The first parameter defining the strong Wolfe conditions (Armijo condition) (should be between 0 and 1
   * exclusively).
   * \param aC2 The second parameter defining the strong Wolfe conditions (curvature condition) (should be between 0 and
   * 1 exclusively, and greater than c1).
   * \param aGeomFactor The geometric decrease factor for the successive decrease of the search interval (should be
   * between 0 and 1 exclusively).
   */
  explicit line_search_backtracking(const T& aC1 = T(1e-2),
                                    const T& aC2 = T(0.9),
                                    const T& aGeomFactor = T(0.75))
      : c1(aC1), c2(aC2), geom_factor(aGeomFactor) {}

  /**
   * This overload computes the scalar value (multiplying the search direction vector) which
   * minimizes the given function along the line defined by the vector offset (x0) and the
   * search direction (dx). The solution will lie between a0 and a1 which are assumed to bound
   * the search domain.
   * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
   * \tparam GradFunction A functor type for a unary function that computes the gradient of the cost for a given
   * independent variable value.
   * \tparam U Any vector type which can be linearly composed (i.e. closed under multiplication by a scalar and under
   * addition).
   * \param f The cost function to minimize (should map U -> T).
   * \param df The gradient of the cost function to minimize (should map U -> U).
   * \param a0 The first bound of the search.
   * \param a1 The second bound of the search.
   * \param x0 The starting point in the vector-space (or offset).
   * \param dx The search direction in the vector-space.
   * \param tol The tolerance at which the search is stopped (i.e. when the interval of convergence is narrow enough).
   * \return The optimal scalar value along the line-search.
   */
  template <typename Function, typename GradFunction, typename U>
  T operator()(Function f, GradFunction df, T a0, T a1, const U& x0,
               const U& dx, const T& tol = T(1e-6)) const {
    backtracking_search(
        [&f, &x0, &dx](const T& alpha) -> T { return f(x0 + alpha * dx); },
        [&df, &x0, &dx](const T& alpha) -> T {
          return dx * df(x0 + alpha * dx);
        },
        a0, a1, tol, c1, c2, geom_factor);
    return a1;
  }
};

/**
 * This function performs a 1D optimum Expand-and-Zoom search on a unimodal cost function f to find the value for
 * the lowest cost (both value and cost are returned by this function as a pair). The search guarantees
 * that the strong Wolfe conditions are met at the obtained optimal point.
 * TEST PASSED
 * \tparam T The value type for the scalar that represents both the independent variable and the cost values.
 * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
 * \param f the cost functor.
 * \param low_bound the lower bound for the search.
 * \param up_bound the upper bound for the search.
 * \param tol the uncertainty on the X value at which the optimization should stop.
 * \param c1 The first parameter defining the Wolf conditions (Armijo condition) (should be between 0 and 1
 * exclusively).
 * \param c2 The second parameter defining the Wolf conditions (curvature condition) (should be between 0 and 1
 * exclusively, and greater than c1).
 * \param tau The geometric decrease factor for the successive decrease of the search interval (should be between 0 and
 * 1 exclusively).
 * \return the cost that is the optimum point found, within uncertainty interval (low_bound, up_bound).
 */
template <typename T, typename Function, typename GradFunction>
T expand_and_zoom_search(Function f, GradFunction df, T& low_bound, T& up_bound,
                         T tol = T(1e-6), T c1 = T(1e-2), T c2 = T(0.9)) {
  return detail::expand_and_zoom_search_impl(f, df, low_bound, up_bound, tol,
                                             c1, c2);
}

/**
 * This functor class wraps a call to a Backtracking search in order to find the minimum of a
 * vector function across a line-search (i.e. with a search direction). This functor only
 * provides bounded searches (with an upper and lower bound provided) and returns the minimum
 * according to the strong Wolfe conditions.
 * TEST PASSED
 * \tparam T A scalar type.
 */
template <typename T>
struct line_search_expand_and_zoom {
  T c1, c2;

  /**
   * Parametrized constructor.
   * \param aC1 The first parameter defining the strong Wolfe conditions (Armijo condition) (should be between 0 and 1
   * exclusively).
   * \param aC2 The second parameter defining the strong Wolfe conditions (curvature condition) (should be between 0 and
   * 1 exclusively, and greater than c1).
   */
  explicit line_search_expand_and_zoom(const T& aC1 = T(1e-2),
                                       const T& aC2 = T(0.9))
      : c1(aC1), c2(aC2) {}

  /**
   * This overload computes the scalar value (multiplying the search direction vector) which
   * minimizes the given function along the line defined by the vector offset (x0) and the
   * search direction (dx). The solution will lie between a0 and a1 which are assumed to bound
   * the search domain.
   * \tparam Function A functor type for a unary function that computes the cost for a given independent variable value.
   * \tparam GradFunction A functor type for a unary function that computes the gradient of the cost for a given
   * independent variable value.
   * \tparam U Any vector type which can be linearly composed (i.e. closed under multiplication by a scalar and under
   * addition).
   * \param f The cost function to minimize (should map U -> T).
   * \param df The gradient of the cost function to minimize (should map U -> U).
   * \param a0 The first bound of the search.
   * \param a1 The second bound of the search.
   * \param x0 The starting point in the vector-space (or offset).
   * \param dx The search direction in the vector-space.
   * \param tol The tolerance at which the search is stopped (i.e. when the interval of convergence is narrow enough).
   * \return The optimal scalar value along the line-search.
   */
  template <typename Function, typename GradFunction, typename U>
  T operator()(Function f, GradFunction df, T a0, T a1, const U& x0,
               const U& dx, const T& tol = T(1e-6)) const {
    expand_and_zoom_search(
        [&f, &x0, &dx](const T& alpha) -> T { return f(x0 + alpha * dx); },
        [&df, &x0, &dx](const T& alpha) -> T {
          return dx * df(x0 + alpha * dx);
        },
        a0, a1, tol, c1, c2);
    return a1;
  }
};

}  // namespace ReaK::optim

#endif  // REAK_MATH_OPTIMIZATION_LINE_SEARCH_H_

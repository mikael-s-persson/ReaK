/**
 * \file nelder_mead_method.h
 *
 * The following library provides an implementation of the Nelder-Mead algorithm. The
 * Nelder-Mead algorithm is a gradient-less non-linear optimization method based on
 * maintaining a simplex (simplest polytope) and is thus sometimes referred to as the
 * simplex method (not to be confused with the simplex method in linear programming,
 * see simplex_method.hpp).
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

#ifndef REAK_MATH_OPTIMIZATION_NELDER_MEAD_METHOD_H_
#define REAK_MATH_OPTIMIZATION_NELDER_MEAD_METHOD_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/math/lin_alg/mat_alg.h"

#include <map>
#include <random>
#include <type_traits>

namespace ReaK::optim {

namespace detail {

template <typename T, typename Vector>
void nelder_mead_compute_center_of_gravity(const std::multimap<T, Vector>& pts,
                                           Vector& c) {
  using ValueType = vect_value_type_t<Vector>;
  ValueType factor = ValueType(1.0) / ValueType(pts.size());
  c = c - c;
  for (auto& pt : pts) {
    c += factor * pt.second;
  }
}

template <typename T, typename Vector>
T nelder_mead_compute_std_dev(const std::multimap<T, Vector>& pts,
                              const Vector& c) {
  using ValueType = vect_value_type_t<Vector>;
  using std::sqrt;

  auto result = ValueType(0.0);
  for (auto& pt : pts) {
    result += norm_2_sqr(pt.second - c);
  }
  return sqrt(result / ValueType(pts.size()));
}

template <typename Function, typename Vector, typename T>
void nelder_mead_method_impl(Function f, std::multimap<T, Vector>& pts,
                             Vector& c, T tol, T alpha = T(1.0),
                             T gamma = T(2.0), T rho = T(0.5),
                             T sigma = T(0.5)) {
  nelder_mead_compute_center_of_gravity(pts, c);
  T abs_tol = nelder_mead_compute_std_dev(pts, c) * tol;
  do {
    Vector x_r = c + alpha * (c - pts.rbegin()->second);
    T x_r_value = f(x_r);
    auto it_r = pts.lower_bound(x_r_value);
    if (it_r == pts.begin()) {
      // do expansion.
      Vector x_e = c + gamma * (pts.rbegin()->second - c);
      T x_e_value = f(x_e);
      if (x_e_value < x_r_value) {
        pts.insert(it_r, std::pair<T, Vector>(x_e_value, x_e));
      } else {
        pts.insert(it_r, std::pair<T, Vector>(x_r_value, x_r));
      }
    } else if ((it_r == (++pts.rbegin()).base()) || (it_r == pts.end())) {
      // do contraction (and possibly reduction after).
      Vector x_c = pts.rbegin()->second + rho * (c - pts.rbegin()->second);
      T x_c_value = f(x_c);
      if (x_c_value < pts.rbegin()->first) {
        pts.insert(it_r, std::pair<T, Vector>(x_c_value, x_c));
      } else {
        // do reduction.
        std::multimap<T, Vector> tmp;
        tmp.swap(pts);
        auto it = tmp.begin();
        ++it;
        auto it2 = pts.insert(*tmp.begin());
        for (; it != pts.end(); ++it) {
          Vector x_t =
              tmp.begin()->second + sigma * (it->second - tmp.begin()->second);
          pts.insert(it2, std::pair<T, Vector>(f(x_t), x_t));
          ++it2;
        }
      }
    } else {
      // do reflection.
      pts.insert(it_r, std::pair<T, Vector>(x_r_value, x_r));
    }

    while (pts.size() > c.size() + 1) {
      pts.erase((++pts.rbegin()).base());
    }

    nelder_mead_compute_center_of_gravity(pts, c);
  } while (nelder_mead_compute_std_dev(pts, c) > abs_tol);
}

}  // namespace detail

/**
 * This function performs a Nelder-Mead method to find the minimum of a non-linear multi-dimensional
 * function. The Nelder-Mead algorithm is a gradient-less non-linear optimization method based on
 * maintaining a simplex (simplest polytope) and is thus sometimes referred to as the
 * simplex method (not to be confused with the simplex method in linear programming,
 * see simplex_method.hpp). This overload version takes a simplex that is already provided by the
 * caller through an iterator range.
 * TEST PASSED but bad convergence (expected).
 * \tparam Function The functor type that can map (Vector -> T).
 * \tparam Vector A writable vector type the represents the independent vector.
 * \tparam T A scalar type.
 * \tparam ForwardIter A forward iterator type that can be used to represent a range of vector-points.
 * \param first The start of the point-range which defines the simplex.
 * \param last The end of the point-range which defines the simplex (the range should contain N+1 points).
 * \param f The function to minimize.
 * \param c The initial guess of the minimum point, and stores as output the resulting minimum.
 * \param tol The tolerance on the spread (std-dev) of the final (converged) simplex.
 * \param alpha The reflection coefficient (default 1.0).
 * \param gamma The expansion coefficient (default 2.0).
 * \param rho The contraction coefficient (default 0.5).
 * \param sigma The shrink coefficient (default 0.5).
 */
template <typename Function, typename Vector, typename T, typename ForwardIter>
void nelder_mead_method(ForwardIter first, ForwardIter last, Function f,
                        Vector& c, T tol = T(1e-6), T alpha = T(1.0),
                        T gamma = T(2.0), T rho = T(0.5), T sigma = T(0.5)) {
  static_assert(is_writable_vector_v<Vector>);
  std::multimap<T, Vector> pts;
  for (; first != last; ++first) {
    pts.insert(std::pair<T, Vector>(f(*first), *first));
  }
  detail::nelder_mead_method_impl(f, pts, c, tol, alpha, gamma, rho, sigma);
}

/**
 * This function performs a Nelder-Mead method to find the minimum of a non-linear multi-dimensional
 * function. The Nelder-Mead algorithm is a gradient-less non-linear optimization method based on
 * maintaining a simplex (simplest polytope) and is thus sometimes referred to as the
 * simplex method (not to be confused with the simplex method in linear programming,
 * see simplex_method.hpp). This overload version takes a random-number generator (see Boost.Random)
 * and uses it to generate the initial simplex around the initial guess with a given initial spread (std-dev).
 * TEST PASSED but bad convergence (expected).
 * \tparam Function The functor type that can map (Vector -> T).
 * \tparam Vector A writable vector type the represents the independent vector.
 * \tparam T A scalar type.
 * \tparam RandomNumberGen A uniform random-number generator (see Boost.Random).
 * \param f The function to minimize.
 * \param c The initial guess of the minimum point, and stores as output the resulting minimum.
 * \param init_spread The initial spread (or standard deviation) of the simplex points around the initial guess.
 * \param rng The random-number generator to use to generate the initial simplex.
 * \param tol The tolerance on the spread (std-dev) of the final (converged) simplex.
 * \param alpha The reflection coefficient (default 1.0).
 * \param gamma The expansion coefficient (default 2.0).
 * \param rho The contraction coefficient (default 0.5).
 * \param sigma The shrink coefficient (default 0.5).
 */
template <typename Function, typename Vector, typename T,
          typename RandomNumberGen>
void nelder_mead_method(Function f, Vector& c, T init_spread,
                        RandomNumberGen& rng, T tol = T(1e-6), T alpha = T(1.0),
                        T gamma = T(2.0), T rho = T(0.5), T sigma = T(0.5)) {
  static_assert(is_writable_vector_v<Vector>);
  using ValueType = vect_value_type_t<Vector>;

  std::multimap<T, Vector> pts;
  std::normal_distribution<ValueType> var_rnd;

  for (int i = 0; i <= c.size(); ++i) {
    Vector z = c;
    for (int j = 0; j < z.size(); ++j) {
      z[j] = var_rnd(rng) * init_spread;
    }
    pts.insert(std::pair<T, Vector>(f(z), z));
  }

  detail::nelder_mead_method_impl(f, pts, c, tol, alpha, gamma, rho, sigma);
}

}  // namespace ReaK::optim

#endif  // REAK_MATH_OPTIMIZATION_NELDER_MEAD_METHOD_H_

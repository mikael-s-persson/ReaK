/**
 * \file quat_num.hpp
 *
 * This library implements basic quaternionic numerical methods like interpolations, integration, etc..
 *
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date March 2015
 */

/*
 *    Copyright 2015 Sven Mikael Persson
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

#ifndef REAK_QUAT_NUM_HPP
#define REAK_QUAT_NUM_HPP

#include "ReaK/math/lin_alg/vect_alg.hpp"

#include "ReaK/math/kinetostatics/quat_alg.hpp"

#include <cassert>
#include <cmath>
#include <type_traits>

namespace ReaK::lie_group {

template <typename T>
T ominus(const T& b, const T& a) {
  return b - a;
}

template <typename T>
quat<T> ominus(const quat<T>& b, const quat<T>& a) {
  return T(2) * log(invert(a) * b);
}

template <typename T>
vect<T, 3> ominus(const unit_quat<T>& b, const unit_quat<T>& a) {
  return T(2) * log(invert(a) * b);
}

template <typename T>
T oplus(const T& a, const T& b) {
  return a + b;
}

template <typename T>
quat<T> oplus(const quat<T>& a, const quat<T>& b) {
  return a * exp(b / T(2));
}

template <typename T>
unit_quat<T> oplus(const unit_quat<T>& a, const vect<T, 3>& b) {
  return a * exp(b / T(2));
}

template <typename T>
T otransport(const T& /*unused*/, const T& v, const T& /*unused*/) {
  return v;
}

template <typename T>
quat<T> otransport(const quat<T>& a, const quat<T>& v, const quat<T>& b) {
  quat<T> ab = invert(b) * a;
  return ab * v * invert(ab);
}

template <typename T>
vect<T, 3> otransport(const unit_quat<T>& a, const vect<T, 3>& v,
                      const unit_quat<T>& b) {
  unit_quat<T> ab = invert(b) * a;
  return ab.rotate(v);
}

template <typename Vector>
Vector ocross(const Vector& a, const Vector& b) {
  return a % b;
}

template <typename T>
struct is_commutative : std::true_type {};

template <typename T>
struct is_commutative<quat<T>> : std::false_type {};

template <typename T>
struct is_commutative<unit_quat<T>> : std::false_type {};

template <typename T>
static constexpr bool is_commutative_t = is_commutative<T>::value;

template <typename T, typename Scalar>
T slerp(const T& q0, const Scalar& eta, const T& q1) {
  return oplus(q0, eta * ominus(q1, q0));
}

template <typename FIter, typename Scalar>
auto bezier(FIter qit_first, FIter qit_last, const Scalar& eta) {
  using Value = typename std::iterator_traits<FIter>::value_type;
  std::vector<Value> qb(qit_first, qit_last);
  std::size_t M = qb.size();
  assert(M > 0);
  for (std::size_t i = 1; i < M; ++i) {
    for (std::size_t j = 0; j < M - i; ++j) {
      qb[j] = slerp(qb[j], eta, qb[j + 1]);
    }
  }
  return qb[0];
}

template <typename FIter, typename Scalar>
auto catmull_rom(FIter qit_first, FIter qit_last, const Scalar& eta,
                 const Scalar& alpha = Scalar(1.0)) {
  using std::pow;
  using Value = typename std::iterator_traits<FIter>::value_type;
  std::vector<Value> qb(qit_first, qit_last);
  std::size_t M = qb.size();
  assert(M > 0);
  if (M == 1) {
    return qb[0];
  }
  // first, compute the t-intervals:
  std::vector<Scalar> tv(M, 0.0);
  for (std::size_t i = 1; i < M; ++i) {
    tv[i] = tv[i - 1] + pow(norm_2(ominus(qb[i], qb[i - 1])), alpha);
  }
  // then, compute t from eta [0..1]
  Scalar t = eta * (tv[tv.size() - 1] - tv[0]);
  // and finally, do the catmull-rom folding:
  for (std::size_t i = 1; i < M; ++i) {
    for (std::size_t j = 0; j < M - i; ++j) {
      qb[j] = slerp(qb[j], (t - tv[j]) / (tv[i + j] - tv[j]), qb[j + 1]);
    }
  }
  return qb[0];
}

template <typename FIter>
auto average(FIter qit_first, FIter qit_last) {
  assert(std::distance(qit_first, qit_last) > 0);
  auto result = *qit_first;
  ++qit_first;
  for (std::size_t i = 1; qit_first != qit_last; ++i, ++qit_first) {
    result = slerp(result, 1.0 / (i + 1), *qit_first);
  }
  return result;
}

template <typename FIter, typename WFIter>
auto weighted_avg(FIter qit_first, FIter qit_last, WFIter wit_first,
                  WFIter wit_last) {
  assert(std::distance(qit_first, qit_last) > 0);
  assert(std::distance(qit_first, qit_last) ==
         std::distance(wit_first, wit_last));
  auto result = *qit_first;
  ++qit_first;
  auto w_sum = *wit_first;
  ++wit_first;
  for (; qit_first != qit_last; ++qit_first, ++wit_first) {
    w_sum += *wit_first;
    result = slerp(result, (*wit_first) / w_sum, *qit_first);
  }
  return result;
}

template <typename T, typename Func, typename Scalar>
T forward_euler(T q, Func f, const Scalar& t_start, const Scalar& t_end,
                Scalar dt) {
  assert(t_end > t_start);
  Scalar t = t_start;
  while (t < t_end) {
    if (t + dt > t_end) {
      dt = t_end - t;
    }
    q = oplus(q, dt * f(t, q));
    t += dt;
  }
  return q;
}

template <typename T, typename Func, typename Scalar>
T backward_euler(T q, Func f, const Scalar& t_start, const Scalar& t_end,
                 Scalar dt, const Scalar& tol) {
  assert(t_end > t_start);
  Scalar t = t_start;
  while (t < t_end) {
    if (t + dt > t_end) {
      dt = t_end - t;
    }
    T q0 = q;
    q = oplus(q0, dt * f(t, q0));
    while (true) {
      T q0_back = oplus(q, (-dt) * f(t + dt, q));
      auto q0_err = ominus(q0, q0_back);
      if (norm_2(q0_err) < tol) {
        break;
      }
      q = oplus(q, otransport(q0_back, q0_err, q));
    }
    t += dt;
  }
  return q;
}

template <typename T, typename Func, typename Scalar>
T trapezoidal_rule(T q, Func f, const Scalar& t_start, const Scalar& t_end,
                   Scalar dt, const Scalar& tol) {
  assert(t_end > t_start);
  Scalar t = t_start;
  while (t < t_end) {
    if (t + dt > t_end) {
      dt = t_end - t;
    }
    auto q0 = q;
    auto w0 = f(t, q0);
    auto q0_mid = oplus(q0, (0.5 * dt) * w0);
    q = oplus(q0, dt * w0);
    while (true) {
      auto q1_mid = oplus(q, (-0.5 * dt) * f(t + dt, q));
      auto q_err = ominus(q0_mid, q1_mid);
      if (norm_2(q_err) < tol) {
        break;
      }
      q = oplus(q, otransport(q1_mid, q_err, q));
    }
    t += dt;
  }
  return q;
}

template <typename T, typename Func, typename Scalar>
T runge_kutta_4(T q, Func f, const Scalar& t_start, const Scalar& t_end,
                Scalar dt) {
  assert(t_end > t_start);
  static const std::array<double, 4> w_vals = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0,
                                               1.0 / 6.0};
  Scalar t = t_start;
  while (t < t_end) {
    if (t + dt > t_end) {
      dt = t_end - t;
    }
    auto p0 = q;
    auto k0 = f(t, p0);
    auto p1 = oplus(p0, (0.5 * dt) * k0);
    auto k1 = otransport(p1, f(t + 0.5 * dt, p1), p0);
    auto p2 = oplus(p0, (0.5 * dt) * k1);
    auto k2 = otransport(p2, f(t + 0.5 * dt, p2), p0);
    auto p3 = oplus(p0, dt * k2);
    auto k3 = f(t + 0.5 * dt, p3);
    std::array<T, 4> p_vals = {oplus(p0, dt * k0), oplus(p0, dt * k1), p3,
                               oplus(p0, dt * k3)};
    q = weighted_avg(p_vals.begin(), p_vals.end(), w_vals.begin(),
                     w_vals.end());
    t += dt;
  }
  return q;
}

/* From Munthe-Kaas 1998, with:
 * A =
 *      0
 *      0.5
 *      0 0.5
 *      0 0 1
 * c = 0 0.5 0.5 1
 * d = 0 0 0.25 0.5
 *
 * (m1 m2 m3) ( 0.5 0.25 0   )
 *            ( 0.5 0.25 0.5 )
 *            ( 1   1    1   ) = (1 0 0)
 * m = 2 2 -1
 */
template <typename T, typename Func, typename Scalar>
T corrected_runge_kutta_4(T q, Func f, const Scalar& t_start,
                          const Scalar& t_end, Scalar dt) {
  assert(t_end > t_start);
  Scalar t = t_start;
  while (t < t_end) {
    if (t + dt > t_end) {
      dt = t_end - t;
    }
    auto p0 = q;
    auto k0 = f(t, p0);
    auto p1 = oplus(p0, (0.5 * dt) * k0);
    auto k1 = f(t + 0.5 * dt, p1);
    auto k1_h2 = (0.5 * dt) * k1;
    auto p2 = oplus(p0, k1_h2 + (dt / 12.0) * ocross(k0, k1_h2));
    auto k2 = f(t + 0.5 * dt, p2);
    auto k2_h = dt * k2;
    auto p3 = oplus(p0, k2_h + (dt / 6.0) * ocross(k0, k2_h));
    auto k3 = f(t + 0.5 * dt, p3);
    auto k_2nd = (2.0 / dt) * (k1 + k2 + 0.5 * k3 - 2.5 * k0);
    auto u_avg = (dt / 6.0) * (k0 + 2.0 * k1 + 2.0 * k2 + k3);
    q = oplus(p0, u_avg + (dt / 4.0) * ocross(k0, u_avg) +
                      (dt * dt / 24.0) * ocross(k_2nd, u_avg));
    t += dt;
  }
  return q;
}

}  // namespace ReaK::lie_group

#endif

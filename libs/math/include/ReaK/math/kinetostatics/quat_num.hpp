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

#include <ReaK/math/lin_alg/vect_alg.hpp>

#include "quat_alg.hpp"

#include <cmath>
#include <cassert>

namespace ReaK {
  
namespace lie_group {

template <typename T>
T ominus(const T& b, const T& a) {
  return b - a;
};

template <typename T>
quat<T> ominus(const quat<T>& b, const quat<T>& a) {
  return T(2) * log(invert(a) * b);
};

template <typename T>
unit_quat<T> ominus(const unit_quat<T>& b, const unit_quat<T>& a) {
  return T(2) * log(invert(a) * b);
};


template <typename T>
T oplus(const T& a, const T& b) {
  return a + b;
};

template <typename T>
quat<T> oplus(const quat<T>& a, const quat<T>& b) {
  return a * exp(b / T(2));
};

template <typename T>
unit_quat<T> oplus(const unit_quat<T>& a, const vect<T,3>& b) {
  return a * exp(b / T(2));
};


template <typename T>
struct is_commutative : boost::mpl::true_ {};

template <typename T>
struct is_commutative< quat<T> > : boost::mpl::false_ {};

template <typename T>
struct is_commutative< unit_quat<T> > : boost::mpl::false_ {};


template <typename T, typename Scalar>
T slerp(const T& q0, const Scalar& eta, const T& q1) {
  return oplus(q0, eta * ominus(q1, q0));
};

template <typename T, typename Scalar>
T bezier(std::vector<T> qb, const Scalar& eta) {
  std::size_t M = qb.size();
  assert(M > 0);
  for(std::size_t i = 1; i < M; ++i)
    for(std::size_t j = 0; j < M - i; ++j)
      qb[j] = slerp(qb[j], eta, qb[j+1]);
  return qb[0];
};

template <typename T, typename Scalar>
T catmull_rom(std::vector<T> qb, const Scalar& eta, const Scalar& alpha = Scalar(1.0)) {
  using std::pow;
  std::size_t M = qb.size();
  assert(M > 0);
  if(M == 1)
    return qb[0];
  // first, compute the t-intervals:
  std::vector<Scalar> tv(M, 0.0);
  for(std::size_t i = 1; i < M; ++i)
    tv[i] = tv[i-1] + pow(norm_2(ominus(q[i], q[i-1])), alpha);
  // then, compute t from eta [0..1]
  Scalar t = eta * (tv[tv.size()-1] - tv[0]);
  // and finally, do the catmull-rom folding:
  for(std::size_t i = 1; i < M; ++i)
    for(std::size_t j = 0; j < M - i; ++j)
      qb[j] = slerp(qb[j], (t - tv[j]) / (tv[i + j] - tv[j]), qb[j+1]);
  return qb[0];
};

template <typename T>
T average(const std::vector<T>& qb) {
  assert(qb.size() > 0);
  T result = qb[0];
  for(std::size_t i = 1; i < qb.size(); ++i)
    result = slerp(result, Scalar(1) / (i+1), qb[i]);
  return result;
};

template <typename T, typename Scalar>
T weighted_avg(const std::vector<T>& qb, const std::vector<Scalar>& w) {
  assert(qb.size() > 0);
  assert(qb.size() == w.size());
  T result = qb[0];
  Scalar w_sum = w[0];
  for(std::size_t i = 1; i < qb.size(); ++i) {
    w_sum += w[i];
    result = slerp(result, w[i] / w_sum, qb[i]);
  };
  return result;
};

};

};

#endif



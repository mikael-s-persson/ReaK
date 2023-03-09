/**
 * \file quat_alg.h
 *
 * This library implements quaternionic algebra.
 *
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date June 2011
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

#ifndef REAK_MATH_KINETOSTATICS_QUAT_ALG_H_
#define REAK_MATH_KINETOSTATICS_QUAT_ALG_H_

#include "ReaK/core/serialization/archiver.h"
#include "ReaK/math/lin_alg/vect_alg.h"

#include "ReaK/math/kinetostatics/rotations_3D.h"

#include <cassert>
#include <cmath>

namespace ReaK {

/**
 * This template class defines a quaternion-valued variable (not a unit-quaternion for representing rotations).
 */
template <class T>
class quat {
 public:
  using self = quat<T>;

  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;
  using allocator_type = void;

  using iterator = typename std::array<value_type, 4>::iterator;
  using const_iterator = typename std::array<value_type, 4>::const_iterator;

  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  using scalar_type = T;
  using vector_type = vect<T, 3>;

  static constexpr std::size_t dimensions = 4;

  /// Holds the four values of the quaternion (q[0] is the real part).
  std::array<value_type, 4> q;

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /**
   * Explicit cast from a real number to a quat.
   */
  explicit quat(const scalar_type& aRe = scalar_type()) noexcept {
    q[0] = aRe;
    q[3] = q[2] = q[1] = value_type();
  }

  /**
   * Constructor from cartesian complex values, real and imaginary parts.
   */
  quat(const scalar_type& aRe, const vector_type& aIm) noexcept {
    q[0] = aRe;
    q[1] = aIm[0];
    q[2] = aIm[1];
    q[3] = aIm[2];
  }

  /**
   * Converts a vector into a pure quaternion.
   */
  explicit quat(const vector_type& aIm) noexcept {
    q[0] = scalar_type();
    q[1] = aIm[0];
    q[2] = aIm[1];
    q[3] = aIm[2];
  }

  /**
   * Converts a 4D vector into a quaternion.
   */
  explicit quat(const vect<value_type, 4>& V) noexcept {
    q[0] = V[0];
    q[1] = V[1];
    q[2] = V[2];
    q[3] = V[3];
  }

  /**
   * Constructs a quaternion from 4 components.
   */
  quat(const_reference q0, const_reference q1, const_reference q2,
       const_reference q3) noexcept {
    q[0] = q0;
    q[1] = q1;
    q[2] = q2;
    q[3] = q3;
  }

  // Copy-constructor is default.
  // Assignment operator is default.

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator[](int i) const noexcept {
    assert(i < q.size());
    return q[i];
  }

  /**
   * Array indexing operator, accessor for lvalue.
   * \test PASSED
   */
  reference operator[](int i) noexcept {
    assert(i < q.size());
    return q[i];
  }

  /**
   * Returns the size of the quaternion (viewed as a 4D vector).
   */
  size_type size() const noexcept { return q.size(); }

  /**
   * Returns an iterator to the first element of the quaternion (viewed as a 4D vector).
   */
  iterator begin() noexcept { return q.begin(); }
  /**
   * Returns a const-iterator to the first element of the quaternion (viewed as a 4D vector).
   */
  const_iterator begin() const noexcept { return q.begin(); }
  /**
   * Returns an iterator to the one-past-last element of the quaternion (viewed as a 4D vector).
   */
  iterator end() noexcept { return q.end(); }
  /**
   * Returns a const-iterator to the one-past-last element of the quaternion (viewed as a 4D vector).
   */
  const_iterator end() const noexcept { return q.end(); }

  /*******************************************************************************
                           Assignment Operators
  *******************************************************************************/

  /**
   * Assignment operator.
   */
  self& operator=(const scalar_type& R) noexcept {
    q[0] = R;
    q[3] = q[2] = q[1] = value_type();
    return *this;
  }

  /**
   * Assignment operator.
   */
  self& operator=(const vector_type& I) noexcept {
    q[0] = value_type();
    q[1] = I[0];
    q[2] = I[1];
    q[3] = I[2];
    return *this;
  }

  /** Addition-assignment operator. */
  self& operator+=(const self& C) noexcept {
    q[0] += C.q[0];
    q[1] += C.q[1];
    q[1] += C.q[2];
    q[1] += C.q[3];
    return *this;
  }

  /** Addition-assignment operator with a real value. */
  self& operator+=(const scalar_type& R) noexcept {
    q[0] += R;
    return *this;
  }

  /** Addition-assignment operator with a real value. */
  self& operator+=(const vector_type& V) noexcept {
    q[1] += V[0];
    q[2] += V[1];
    q[3] += V[2];
    return *this;
  }

  /** Substraction-assignment operator. */
  self& operator-=(const self& C) noexcept {
    q[0] -= C.q[0];
    q[1] -= C.q[1];
    q[1] -= C.q[2];
    q[1] -= C.q[3];
    return *this;
  }

  /** Substraction-assignment operator with a real value. */
  self& operator-=(const scalar_type& R) noexcept {
    q[0] -= R;
    return *this;
  }

  /** Substraction-assignment operator with a real value. */
  self& operator-=(const vector_type& V) noexcept {
    q[1] -= V[0];
    q[2] -= V[1];
    q[3] -= V[2];
    return *this;
  }

  /** Multiplication-assignment operator. */
  self& operator*=(const self& C) noexcept { return (*this = ((*this) * C)); }

  /** Multiplication-assignment operator with a real value. */
  self& operator*=(const scalar_type& R) noexcept {
    q[0] *= R;
    q[1] *= R;
    q[2] *= R;
    q[3] *= R;
    return *this;
  }

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /** Addition operator. */
  friend self operator+(const self& C1, const self& C2) noexcept {
    return {C1.q[0] + C2.q[0], C1.q[1] + C2.q[1], C1.q[2] + C2.q[2],
            C1.q[3] + C2.q[3]};
  }

  /** Addition operator. */
  friend self operator+(const self& C1, const scalar_type& C2) noexcept {
    return {C1.q[0] + C2, C1.q[1], C1.q[2], C1.q[3]};
  }

  /** Addition operator. */
  friend self operator+(const scalar_type& C1, const self& C2) noexcept {
    return {C2.q[0] + C1, C2.q[1], C2.q[2], C2.q[3]};
  }

  /** Addition operator. */
  friend self operator+(const self& C1, const vector_type& C2) noexcept {
    return {C1.q[0], C1.q[1] + C2[0], C1.q[2] + C2[1], C1.q[3] + C2[2]};
  }

  /** Addition operator. */
  friend self operator+(const vector_type& C1, const self& C2) noexcept {
    return {C2.q[0], C1[0] + C2.q[1], C1[1] + C2.q[2], C1[2] + C2.q[3]};
  }

  /** Negation operator. */
  friend self operator-(const self& Q) noexcept {
    return {-Q.q[0], -Q.q[1], -Q.q[2], -Q.q[3]};
  }

  /** Substraction operator. */
  friend self operator-(const self& C1, const self& C2) noexcept {
    return {C1.q[0] - C2.q[0], C1.q[1] - C2.q[1], C1.q[2] - C2.q[2],
            C1.q[3] - C2.q[3]};
  }

  /** Substraction operator. */
  friend self operator-(const self& C1, const scalar_type& C2) noexcept {
    return {C1.q[0] - C2, C1.q[1], C1.q[2], C1.q[3]};
  }

  /** Substraction operator. */
  friend self operator-(const scalar_type& C1, const self& C2) noexcept {
    return {C1 - C2.q[0], -C2.q[1], -C2.q[2], -C2.q[3]};
  }

  /** Substraction operator. */
  friend self operator-(const self& C1, const vector_type& C2) noexcept {
    return {C1.q[0], C1.q[1] - C2[0], C1.q[2] - C2[1], C1.q[3] - C2[2]};
  }

  /** Substraction operator. */
  friend self operator-(const vector_type& C1, const self& C2) noexcept {
    return {-C2.q[0], C1[0] - C2.q[1], C1[1] - C2.q[2], C1[2] - C2.q[3]};
  }

  /*******************************************************************************
                             Comparison Operators
  *******************************************************************************/

  /** Equality-comparison operator. */
  friend bool operator==(const self& C1, const self& C2) noexcept {
    return ((C1.q[0] == C2.q[0]) && (C1.q[1] == C2.q[1]) &&
            (C1.q[2] == C2.q[2]) && (C1.q[3] == C2.q[3]));
  }

  /** Equality-comparison operator with a real value. */
  friend bool operator==(const self& C, const scalar_type& R) noexcept {
    return ((C.q[0] == R) && (C.q[1] == value_type()) &&
            (C.q[2] == value_type()) && (C.q[3] == value_type()));
  }

  /** Equality-comparison operator with a real value. */
  friend bool operator==(const scalar_type& R, const self& C) noexcept {
    return ((C.q[0] == R) && (C.q[1] == value_type()) &&
            (C.q[2] == value_type()) && (C.q[3] == value_type()));
  }

  /** Equality-comparison operator with a vector value. */
  friend bool operator==(const self& C, const vector_type& V) noexcept {
    return ((C.q[0] == scalar_type()) && (C.q[1] == V[0]) && (C.q[2] == V[1]) &&
            (C.q[3] == V[2]));
  }

  /** Equality-comparison operator with a vector value. */
  friend bool operator==(const vector_type& V, const self& C) noexcept {
    return ((C.q[0] == scalar_type()) && (C.q[1] == V[0]) && (C.q[2] == V[1]) &&
            (C.q[3] == V[2]));
  }

  /** Inequality-comparison operator. */
  friend bool operator!=(const self& C1, const self& C2) noexcept {
    return ((C1.q[0] != C2.q[0]) || (C1.q[1] != C2.q[1]) ||
            (C1.q[2] != C2.q[2]) || (C1.q[3] != C2.q[3]));
  }

  /** Inequality-comparison operator with a real value. */
  friend bool operator!=(const self& C, const scalar_type& R) noexcept {
    return ((C.q[0] != R) || (C.q[1] != value_type()) ||
            (C.q[2] != value_type()) || (C.q[3] != value_type()));
  }

  /** Inequality-comparison operator with a real value. */
  friend bool operator!=(const scalar_type& R, const self& C) noexcept {
    return ((C.q[0] != R) || (C.q[1] != value_type()) ||
            (C.q[2] != value_type()) || (C.q[3] != value_type()));
  }

  /** Inequality-comparison operator with a vector value. */
  friend bool operator!=(const self& C, const vector_type& V) noexcept {
    return ((C.q[0] == scalar_type()) || (C.q[1] != V[0]) || (C.q[2] != V[1]) ||
            (C.q[3] != V[2]));
  }

  /** Inequality-comparison operator with a vector value. */
  friend bool operator!=(const vector_type& V, const self& C) noexcept {
    return ((C.q[0] != scalar_type()) || (C.q[1] != V[0]) || (C.q[2] != V[1]) ||
            (C.q[3] != V[2]));
  }

  /** Quaternionic conjugate for a quaternion value. */
  friend self conj(const self& x) noexcept {
    return {x.q[0], -x.q[1], -x.q[2], -x.q[3]};
  }

  /**
   * Square magnitude of the quaternion.
   * \test PASSED
   */
  friend value_type norm_2_sqr(const self& v) noexcept {
    return v.q[0] * v.q[0] + v.q[1] * v.q[1] + v.q[2] * v.q[2] +
           v.q[3] * v.q[3];
  }

  /**
   * Magnitude of the quaternion.
   * \test PASSED
   */
  friend value_type norm_2(const self& v) noexcept {
    using std::sqrt;
    return sqrt(norm_2_sqr(v));
  }

  /**
   * Unit quaternion in the same direction.
   * \test PASSED
   */
  friend self unit(const self& v) noexcept { return v * (1.0 / norm_2(v)); }

  /**
   * Checks if two quaternions are colinear.
   * \test PASSED
   */
  friend bool colinear(const self& v1, const self& v2) noexcept {
    using std::abs;
    T tmp_mag1 = norm_2(v1);
    T tmp_mag2 = norm_2(v2);
    T tmp_comb = norm_2(v1 + v2);
    return (
        ((tmp_mag1 + tmp_mag2) * (T(1.0) - std::numeric_limits<T>::epsilon()) <=
         tmp_comb) ||
        (abs(tmp_mag1 - tmp_mag2) *
             (T(1.0) + std::numeric_limits<T>::epsilon()) >=
         tmp_comb));
  }

  // Exponential and logarithmic functions:

  /** Compute exponential function (function), for a quaternion value. */
  friend self exp(const self& x) noexcept {
    using std::cos;
    using std::exp;
    using std::sin;
    using std::sqrt;
    value_type theta =
        sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
    if (theta < std::numeric_limits<value_type>::epsilon()) {
      return self(exp(x.q[0]));
    }
    value_type fact = sin(theta) / theta;
    return exp(x.q[0]) *
           self(cos(theta), fact * x.q[1], fact * x.q[2], fact * x.q[3]);
  }

  /** Compute natural logarithm (function), for a quaternion value. */
  friend self log(const self& x) noexcept {
    using std::abs;
    using std::atan2;
    using std::log;
    using std::sqrt;
    value_type st = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
    if (st < std::numeric_limits<value_type>::epsilon()) {
      return self(log(abs(x.q[0])));
    }
    value_type fact = atan2(st, x.q[0]) / st;
    return self(log(sqrt(x.q[0] * x.q[0] + st * st)), fact * x.q[1],
                fact * x.q[2], fact * x.q[3]);
  }

  // Power functions

  /** Raise to power (function), for a quaternion value.*/
  friend self pow(const self& base, const self& exponent) noexcept {
    return exp(exponent * log(base));
  }

  /** Compute square root (function), for a quaternion value.*/
  friend self sqrt(const self& x) noexcept { return exp(0.5 * log(x)); }

  /**
   * Inverts the quaternion.
   */
  friend self invert(const self& x) noexcept {
    value_type tmp = 1.0 / norm_2_sqr(x);
    return {x.q[0] * tmp, -x.q[1] * tmp, -x.q[2] * tmp, -x.q[3] * tmp};
  }

  // Rounding, absolute value and remainder functions:

  /** Round up value (function), for a quaternion value. */
  friend self ceil(const self& x) noexcept {
    return {ceil(x.q[0]), ceil(x.q[1]), ceil(x.q[2]), ceil(x.q[3])};
  }

  /** Compute absolute value (function), for a quaternion value. */
  friend value_type abs(const self& x) noexcept { return norm_2(x); }

  /** Round down value (function), for a quaternion value.*/
  friend self floor(const self& x) noexcept {
    return {floor(x.q[0]), floor(x.q[1]), floor(x.q[2]), floor(x.q[3])};
  }

  // Trigonometric functions:

  /** Compute cosine (function), for a quaternion value.*/
  friend self cos(const self& x) noexcept {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;
    using std::sqrt;
    value_type theta =
        sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
    if (theta < std::numeric_limits<value_type>::epsilon()) {
      return self(cos(x.q[0]));
    }
    value_type fact = -sin(x.q[0]) * sinh(theta) / theta;
    return self(cos(x.q[0]) * cosh(theta), fact * x.q[1], fact * x.q[2],
                fact * x.q[3]);
  }

  /** Compute sine (function), for a quaternion value.*/
  friend self sin(const self& x) noexcept {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;
    using std::sqrt;
    value_type theta =
        sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
    if (theta < std::numeric_limits<value_type>::epsilon()) {
      return self(sin(x.q[0]));
    }
    value_type fact = cos(x.q[0]) * sinh(theta) / theta;
    return self(sin(x.q[0]) * cosh(theta), fact * x.q[1], fact * x.q[2],
                fact * x.q[3]);
  }

  /** Compute tangent (function), for a quaternion value.*/
  friend self tan(const self& x) noexcept {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;
    using std::sqrt;
    using std::tan;
    value_type theta =
        sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
    if (theta < std::numeric_limits<value_type>::epsilon()) {
      return self(tan(x.q[0]));
    }
    value_type tmp = 1.0 / (cos(2.0 * x.q[0]) + cosh(2.0 * theta));
    value_type fact = sinh(2.0 * theta) * tmp / theta;
    return self(sin(2.0 * x.q[0]) * tmp, fact * x.q[1], fact * x.q[2],
                fact * x.q[3]);
  }

  /** Compute arc cosine (function), for a quaternion value.*/
  friend self acos(const self& x) noexcept {
    using std::acos;
    using std::atan2;
    using std::cos;
    using std::log;
    using std::pow;
    using std::sin;
    using std::sqrt;
    value_type ss_sht =
        sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
    if (ss_sht < std::numeric_limits<value_type>::epsilon()) {
      return self(acos(x.q[0]));
    }

    // complex acos...
    value_type Re1 = x.q[0] * x.q[0] - ss_sht * ss_sht - 1.0;  // x * x - 1.0
    value_type Im1 = 2.0 * x.q[0] * ss_sht;                    // ...........
    value_type angle = atan2(Im1, Re1);      // sqrt( x * x - 1.0 ) + x
    Im1 = pow(Re1 * Re1 + Im1 * Im1, 0.25);  //.......................
    Re1 = cos(0.5 * angle) * Im1 + x.q[0];   //.......................
    Im1 = Im1 * sin(0.5 * angle) + ss_sht;   //.......................

    angle = atan2(Im1, Re1);                  //(-1.0i) * log()
    Im1 = -log(sqrt(Re1 * Re1 + Im1 * Im1));  //...............
    Re1 = angle;                              //...............
    // end complex acos.

    value_type fact = Im1 / ss_sht;
    return {Re1, fact * x.q[1], fact * x.q[2], fact * x.q[3]};
  }

  /** Compute arc sine (function), for a quaternion value.*/
  friend self asin(const self& x) noexcept {
    using std::asin;
    using std::atan2;
    using std::cos;
    using std::log;
    using std::pow;
    using std::sin;
    using std::sqrt;
    value_type cs_sht =
        sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
    if (cs_sht < std::numeric_limits<value_type>::epsilon()) {
      return self(asin(x.q[0]));
    }

    // complex asin...
    value_type Re1 = 1.0 - x.q[0] * x.q[0] + cs_sht * cs_sht;  // 1.0 - x * x
    value_type Im1 = -2.0 * x.q[0] * cs_sht;                   // ...........
    value_type angle = atan2(Im1, Re1);  // sqrt( 1.0 - x * x ) + x.Re i - x.Im
    Im1 = pow(Re1 * Re1 + Im1 * Im1, 0.25);  //.......................
    Re1 = cos(0.5 * angle) * Im1 - cs_sht;   //.......................
    Im1 = Im1 * sin(0.5 * angle) + x.q[0];   //.......................

    angle = atan2(Im1, Re1);                  //(-1.0i) * log()
    Im1 = -log(sqrt(Re1 * Re1 + Im1 * Im1));  //...............
    Re1 = angle;                              //...............
    // end complex asin.

    value_type fact = Im1 / cs_sht;
    return {Re1, fact * x.q[1], fact * x.q[2], fact * x.q[3]};
  }

  /** Compute arc tangent (function), for a quaternion value.*/
  friend self atan(const self& x) noexcept {
    using std::atan;
    using std::atan2;
    using std::cos;
    using std::log;
    using std::pow;
    using std::sin;
    using std::sqrt;
    value_type tmp = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
    if (tmp < std::numeric_limits<value_type>::epsilon()) {
      return self(atan(x.q[0]));
    }

    // complex atan...
    value_type sqr_mag_inv = 1.0 / ((1.0 - tmp) * (1.0 - tmp) +
                                    x.q[0] * x.q[0]);  // complex division
    value_type Re1 =
        sqr_mag_inv * (1.0 - tmp * tmp - x.q[0] * x.q[0]);  //................
    value_type Im1 = -2.0 * sqr_mag_inv * x.q[0];           //................
    value_type angle = atan2(Im1, Re1);                     //(0.5i) * log()
    Im1 = 0.5 * log(sqrt(Re1 * Re1 + Im1 * Im1));           //..............
    Re1 = -0.5 * angle;                                     //..............
    // end complex atan.

    value_type fact = Im1 / tmp;
    return {Re1, fact * x.q[1], fact * x.q[2], fact * x.q[3]};
  }

  /** Compute arc tangent with two parameters (function), for a quaternion value.*/
  friend self atan2(const self& y, const self& x) noexcept {
    return atan(y * invert(x));
  };

  // Hyperbolic functions:

  /** Compute hyperbolic cosine (function), for a quaternion value.*/
  friend self cosh(const self& x) noexcept {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;
    using std::sqrt;
    value_type theta =
        sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
    if (theta < std::numeric_limits<value_type>::epsilon()) {
      return self(cosh(x.q[0]));
    }
    value_type fact = sinh(x.q[0]) * sin(theta) / theta;
    return self(cosh(x.q[0]) * cos(theta), fact * x.q[1], fact * x.q[2],
                fact * x.q[3]);
  }

  /** Compute hyperbolic sine (function), for a quaternion value.*/
  friend self sinh(const self& x) noexcept {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;
    using std::sqrt;
    value_type theta =
        sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
    if (theta < std::numeric_limits<value_type>::epsilon()) {
      return self(sinh(x.q[0]));
    }
    value_type fact = cosh(x.q[0]) * sin(theta) / theta;
    return self(sinh(x.q[0]) * cos(theta), fact * x.q[1], fact * x.q[2],
                fact * x.q[3]);
  }

  /** Compute hyperbolic tangent (function), for a quaternion value.*/
  friend self tanh(const self& x) noexcept {
    using std::cos;
    using std::cosh;
    using std::sin;
    using std::sinh;
    using std::sqrt;
    value_type theta =
        sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
    if (theta < std::numeric_limits<value_type>::epsilon()) {
      return self(tanh(x.q[0]));
    }
    value_type tmp = 1.0 / (cosh(2.0 * x.q[0]) + cos(2.0 * theta));
    value_type fact = sin(2.0 * theta) * tmp / theta;
    return self(sinh(2.0 * x.q[0]) * tmp, fact * x.q[1], fact * x.q[2],
                fact * x.q[3]);
  }
};

/** Multiplication by a quaternion. */
template <typename T>
quat<T> operator*(const quat<T>& Q1, const quat<T>& Q2) noexcept {
  return {Q2[0] * Q1[0] - Q2[1] * Q1[1] - Q2[2] * Q1[2] - Q2[3] * Q1[3],
          Q2[0] * Q1[1] + Q2[3] * Q1[2] - Q2[2] * Q1[3] + Q2[1] * Q1[0],
          Q2[0] * Q1[2] - Q2[3] * Q1[1] + Q2[1] * Q1[3] + Q2[2] * Q1[0],
          Q2[0] * Q1[3] + Q2[2] * Q1[1] - Q2[1] * Q1[2] + Q2[3] * Q1[0]};
}

/** Multiplication by a scalar. */
template <typename T>
quat<T> operator*(const quat<T>& Q1, const T& Q2) noexcept {
  return {Q2 * Q1[0], Q2 * Q1[1], Q2 * Q1[2], Q2 * Q1[3]};
}

/** Multiplication by a scalar. */
template <typename T>
quat<T> operator*(const T& Q1, const quat<T>& Q2) noexcept {
  return {Q1 * Q2[0], Q1 * Q2[1], Q1 * Q2[2], Q1 * Q2[3]};
}

/** Division by a scalar. */
template <typename T>
quat<T> operator/(const quat<T>& Q1, const T& Q2) noexcept {
  return {Q1[0] / Q2, Q1[1] / Q2, Q1[2] / Q2, Q1[3] / Q2};
}

/**
 * Multiplication by a quaternion.
 */
template <typename T>
quat<T> operator*(const quat<T>& Q1, const vect<T, 3>& Q2) noexcept {
  return {-Q2[0] * Q1[1] - Q2[1] * Q1[2] - Q2[2] * Q1[3],
          Q2[2] * Q1[2] - Q2[1] * Q1[3] + Q2[0] * Q1[0],
          -Q2[2] * Q1[1] + Q2[0] * Q1[3] + Q2[1] * Q1[0],
          Q2[1] * Q1[1] - Q2[0] * Q1[2] + Q2[2] * Q1[0]};
}

/**
 * Multiplication by a quaternion.
 */
template <typename T>
quat<T> operator*(const vect<T, 3>& Q1, const quat<T>& Q2) noexcept {
  return {-Q2[1] * Q1[0] - Q2[2] * Q1[1] - Q2[3] * Q1[2],
          Q2[0] * Q1[0] + Q2[3] * Q1[1] - Q2[2] * Q1[2],
          Q2[0] * Q1[1] - Q2[3] * Q1[0] + Q2[1] * Q1[2],
          Q2[0] * Q1[2] + Q2[2] * Q1[0] - Q2[1] * Q1[1]};
}

/**
 * This template class defines a quaternion-valued variable (not a unit-quaternion for representing rotations).
 */
template <class T>
class unit_quat : public quat<T> {
 public:
  using self = unit_quat<T>;

  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;
  using allocator_type = void;

  using iterator = typename quat<T>::iterator;
  using const_iterator = typename quat<T>::const_iterator;

  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  using scalar_type = T;
  using vector_type = vect<T, 3>;

  static constexpr std::size_t dimensions = 4;

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /**
   * Default constructor, always yields (1.0, 0.0, 0.0, 0.0).
   */
  unit_quat() noexcept : quat<T>(scalar_type(1.0)) {}

  /**
   * Constructor from quaternion value.
   */
  explicit unit_quat(const quat<T>& aQ) noexcept : quat<T>(scalar_type(1.0)) {
    using std::sqrt;
    scalar_type factor = sqrt(aQ.q[0] * aQ.q[0] + aQ.q[1] * aQ.q[1] +
                              aQ.q[2] * aQ.q[2] + aQ.q[3] * aQ.q[3]);
    if (factor > std::numeric_limits<scalar_type>::epsilon()) {
      factor = 1.0 / factor;
      this->q[0] = aQ.q[0] * factor;
      this->q[1] = aQ.q[1] * factor;
      this->q[2] = aQ.q[2] * factor;
      this->q[3] = aQ.q[3] * factor;
    }
  }

  /**
   * Converts a 4D vector into a quaternion.
   */
  explicit unit_quat(const vect<value_type, 4>& V) noexcept
      : quat<T>(scalar_type(1.0)) {
    using std::sqrt;
    scalar_type factor = norm_2(V);
    if (factor > std::numeric_limits<scalar_type>::epsilon()) {
      factor = 1.0 / factor;
      this->q[0] = V[0] * factor;
      this->q[1] = V[1] * factor;
      this->q[2] = V[2] * factor;
      this->q[3] = V[3] * factor;
    }
  }

  /**
   * Constructs a quaternion from 4 components.
   */
  unit_quat(const_reference q0, const_reference q1, const_reference q2,
            const_reference q3) noexcept
      : quat<T>(scalar_type(1.0)) {
    using std::sqrt;
    scalar_type factor = sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    if (factor > std::numeric_limits<scalar_type>::epsilon()) {
      factor = 1.0 / factor;
      this->q[0] = q0 * factor;
      this->q[1] = q1 * factor;
      this->q[2] = q2 * factor;
      this->q[3] = q3 * factor;
    }
  }

  explicit unit_quat(const quaternion<value_type>& aQ) noexcept
      : quat<T>(aQ[0], aQ[1], aQ[2], aQ[3]) {}

  // Copy-constructor is default.
  // Assignment operator is default.

  // NOTE hiding the non-const overloads in the base class is intentional here:

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator[](size_type i) const noexcept {
    assert(i < 4);
    return this->q[i];
  }

  /**
   * Returns a const-iterator to the first element of the quaternion (viewed as a 4D vector).
   */
  const_iterator begin() const noexcept { return this->q.begin(); }
  /**
   * Returns a const-iterator to the one-past-last element of the quaternion (viewed as a 4D vector).
   */
  const_iterator end() const noexcept { return this->q.end(); }

  /**
   * Returns the same unit-quaternion, but now as a quaternion class, which is used to represent 3D rotations.
   */
  quaternion<value_type> as_rotation() const noexcept {
    return quaternion<value_type>(this->q[0], this->q[1], this->q[2],
                                  this->q[3]);
  }

  /**
   * Multiplication by a column vector.
   * \test PASSED
   */
  vector_type rotate(const vector_type& V) noexcept {
    std::array<value_type, 9> t;
    t[0] = this->q[0] * this->q[1];
    t[1] = this->q[0] * this->q[2];
    t[2] = this->q[0] * this->q[3];
    t[3] = -this->q[1] * this->q[1];
    t[4] = this->q[1] * this->q[2];
    t[5] = this->q[1] * this->q[3];
    t[6] = -this->q[2] * this->q[2];
    t[7] = this->q[2] * this->q[3];
    t[8] = -this->q[3] * this->q[3];
    return {value_type(2.0) * ((t[6] + t[8]) * V[0] + (t[4] - t[2]) * V[1] +
                               (t[1] + t[5]) * V[2]) +
                V[0],
            value_type(2.0) * ((t[2] + t[4]) * V[0] + (t[3] + t[8]) * V[1] +
                               (t[7] - t[0]) * V[2]) +
                V[1],
            value_type(2.0) * ((t[5] - t[1]) * V[0] + (t[0] + t[7]) * V[1] +
                               (t[3] + t[6]) * V[2]) +
                V[2]};
  }

  /*******************************************************************************
                           Assignment Operators
  *******************************************************************************/

  /**
   * Assignment operator.
   */
  self& operator=(const self& Q) noexcept {
    this->q[0] = Q.q[0];
    this->q[1] = Q.q[1];
    this->q[2] = Q.q[2];
    this->q[3] = Q.q[3];
    return *this;
  }

  /** Multiplication-assignment operator. */
  self& operator*=(const self& C) noexcept { return (*this = ((*this) * C)); }

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /** Negation operator. */
  friend self operator-(self Q) noexcept {
    Q.q[0] = -Q.q[0];
    Q.q[1] = -Q.q[1];
    Q.q[2] = -Q.q[2];
    Q.q[3] = -Q.q[3];
    return Q;
  }

  /** Quaternionic conjugate for a quaternion value. */
  friend self conj(const self& x) noexcept {
    self result(x);
    result.q[1] = -x.q[1];
    result.q[2] = -x.q[2];
    result.q[3] = -x.q[3];
    return result;
  }

  /**
   * Square magnitude of the quaternion.
   * \test PASSED
   */
  friend value_type norm_2_sqr(const self& v) noexcept {
    return value_type(1.0);
  }

  /**
   * Magnitude of the quaternion.
   * \test PASSED
   */
  friend value_type norm_2(const self& v) noexcept { return value_type(1.0); }

  /**
   * Unit quaternion in the same direction.
   * \test PASSED
   */
  friend self unit(const self& v) noexcept { return v; }

  // Exponential and logarithmic functions:

  /** Compute natural logarithm (function), for a quaternion value. */
  friend vector_type log(const self& x) noexcept {
    using std::atan2;
    using std::sqrt;
    value_type st = sqrt(x.q[1] * x.q[1] + x.q[2] * x.q[2] + x.q[3] * x.q[3]);
    if (st < std::numeric_limits<value_type>::epsilon()) {
      return vector_type(value_type(0.0), value_type(0.0), value_type(0.0));
    }
    value_type fact;
    if (x.q[0] < 0.0) {
      fact = -atan2(st, -x.q[0]) / st;
    } else {
      fact = atan2(st, x.q[0]) / st;
    }
    return {fact * x.q[1], fact * x.q[2], fact * x.q[3]};
  }

  // Power functions

  /** Raise to power (function), for a quaternion value.*/
  friend self pow(const self& base, const scalar_type& exponent) noexcept {
    return exp(exponent * log(base));
  }

  /** Compute square root (function), for a quaternion value.*/
  friend self sqrt(const self& x) noexcept { return exp(0.5 * log(x)); }

  /**
   * Inverts the quaternion.
   */
  friend self invert(const self& x) noexcept { return conj(x); }

  // Rounding, absolute value and remainder functions:

  /** Compute absolute value (function), for a quaternion value. */
  friend value_type abs(const self& x) noexcept { return 1.0; }
};

/** Multiplication by a quaternion. */
template <typename T>
unit_quat<T> operator*(const unit_quat<T>& Q1,
                       const unit_quat<T>& Q2) noexcept {
  return {Q2[0] * Q1[0] - Q2[1] * Q1[1] - Q2[2] * Q1[2] - Q2[3] * Q1[3],
          Q2[0] * Q1[1] + Q2[3] * Q1[2] - Q2[2] * Q1[3] + Q2[1] * Q1[0],
          Q2[0] * Q1[2] - Q2[3] * Q1[1] + Q2[1] * Q1[3] + Q2[2] * Q1[0],
          Q2[0] * Q1[3] + Q2[2] * Q1[1] - Q2[1] * Q1[2] + Q2[3] * Q1[0]};
}

/** Multiplication by a quaternion. */
template <typename T>
quat<T> operator*(const unit_quat<T>& Q1, const quat<T>& Q2) noexcept {
  return static_cast<const quat<T>&>(Q1) * Q2;
}

/** Multiplication by a quaternion. */
template <typename T>
quat<T> operator*(const quat<T>& Q1, const unit_quat<T>& Q2) noexcept {
  return Q1 * static_cast<const quat<T>&>(Q2);
}

/** Multiplication by a scalar. */
template <typename T>
quat<T> operator*(const unit_quat<T>& Q1, const T& Q2) noexcept {
  return static_cast<const quat<T>&>(Q1) * Q2;
}

/** Multiplication by a scalar. */
template <typename T>
quat<T> operator*(const T& Q1, const unit_quat<T>& Q2) noexcept {
  return Q1 * static_cast<const quat<T>&>(Q2);
}

/** Division by a scalar. */
template <typename T>
quat<T> operator/(const unit_quat<T>& Q1, const T& Q2) noexcept {
  return {Q1[0] / Q2, Q1[1] / Q2, Q1[2] / Q2, Q1[3] / Q2};
}

/** Multiplication by a quaternion. */
template <typename T>
quat<T> operator*(const unit_quat<T>& Q1, const vect<T, 3>& Q2) noexcept {
  return static_cast<const quat<T>&>(Q1) * Q2;
}

/** Multiplication by a quaternion. */
template <typename T>
quat<T> operator*(const vect<T, 3>& Q1, const unit_quat<T>& Q2) noexcept {
  return Q1 * static_cast<const quat<T>&>(Q2);
}

/** Compute exponential function (function), for a quaternion value. */
template <typename T>
unit_quat<T> exp(const vect<T, 3>& x) noexcept {
  using std::cos;
  using std::exp;
  using std::sin;
  using std::sqrt;
  T theta = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  if (theta < std::numeric_limits<T>::epsilon()) {
    return unit_quat<T>();
  }
  T fact = sin(theta) / theta;
  return {cos(theta), fact * x[0], fact * x[1], fact * x[2]};
}

namespace serialization {

/// Loading a quaternion value.
template <typename T>
iarchive& operator>>(iarchive& in, quat<T>& C) {
  return in >> C.q[0] >> C.q[1] >> C.q[2] >> C.q[3];
}

/// Loading a quaternion value with a name.
template <typename T>
iarchive& operator&(iarchive& in, const std::pair<std::string, quat<T>&>& f) {
  return in & std::pair<std::string, T&>(f.first + "_q0", f.second.q[0]) &
         std::pair<std::string, T&>(f.first + "_q1", f.second.q[1]) &
         std::pair<std::string, T&>(f.first + "_q2", f.second.q[2]) &
         std::pair<std::string, T&>(f.first + "_q3", f.second.q[3]);
}

/// Saving a quaternion value.
template <typename T>
oarchive& operator<<(oarchive& out, const quat<T>& C) {
  return out << C.q[0] << C.q[1] << C.q[2] << C.q[3];
}

/// Saving a quaternion value with a name.
template <typename T>
oarchive& operator&(oarchive& out,
                    const std::pair<std::string, const quat<T>&>& C) {
  return out & std::pair<std::string, T>(C.first + "_q0", C.second.q[0]) &
         std::pair<std::string, T>(C.first + "_q1", C.second.q[1]) &
         std::pair<std::string, T>(C.first + "_q2", C.second.q[2]) &
         std::pair<std::string, T>(C.first + "_q3", C.second.q[3]);
}

}  // namespace serialization

namespace rtti {

template <typename T>
struct get_type_id<quat<T>> {
  static constexpr unsigned int ID = 0x0000002B;
  static constexpr auto type_name = std::string_view{"ReaK::quat"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }

  using save_type = const quat<T>&;
  using load_type = quat<T>&;
};

template <typename T>
struct get_type_id<unit_quat<T>> {
  static constexpr unsigned int ID = 0x0000002D;
  static constexpr auto type_name = std::string_view{"ReaK::unit_quat"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }

  using save_type = const quat<T>&;
  using load_type = quat<T>&;
};
};  // namespace rtti

template <typename T>
struct is_readable_vector<quat<T>> {
  static constexpr bool value = true;
  using type = is_readable_vector<quat<T>>;
};

template <typename T>
struct is_writable_vector<quat<T>> {
  static constexpr bool value = true;
  using type = is_writable_vector<quat<T>>;
};

template <typename T>
struct is_resizable_vector<quat<T>> {
  static constexpr bool value = false;
  using type = is_resizable_vector<quat<T>>;
};

template <typename T>
struct has_allocator_vector<quat<T>> {
  static constexpr bool value = false;
  using type = has_allocator_vector<quat<T>>;
};

template <typename T>
struct is_readable_vector<unit_quat<T>> {
  static constexpr bool value = true;
  using type = is_readable_vector<unit_quat<T>>;
};

template <typename T>
struct is_writable_vector<unit_quat<T>> {
  static constexpr bool value = false;
  using type = is_writable_vector<unit_quat<T>>;
};

template <typename T>
struct is_resizable_vector<unit_quat<T>> {
  static constexpr bool value = false;
  using type = is_resizable_vector<unit_quat<T>>;
};

template <typename T>
struct has_allocator_vector<unit_quat<T>> {
  static constexpr bool value = false;
  using type = has_allocator_vector<unit_quat<T>>;
};

namespace detail {

template <typename T, typename Vector2>
inline void from_vect_impl(unit_quat<T>& lhs, const Vector2& rhs,
                           std::size_t& i) {
  lhs = unit_quat<T>(rhs[i], rhs[i + 1], rhs[i + 2], rhs[i + 3]);
  i += 4;
}

}  // namespace detail

}  // namespace ReaK

#endif  // REAK_MATH_KINETOSTATICS_QUAT_ALG_H_

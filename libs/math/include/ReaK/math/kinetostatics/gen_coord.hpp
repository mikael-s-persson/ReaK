/**
 * \file gen_coord.hpp
 *
 * This library implements the gen_coord class template which can be used to represent a 1D kinetostatic
 * frame of reference, i.e., a generalized coordinate.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2011 (last revision)
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

#ifndef REAK_GEN_COORD_HPP
#define REAK_GEN_COORD_HPP

#include "ReaK/core/base/shared_object.hpp"

namespace ReaK {

/**
 * This class holds the kinematic and dynamic values for a generalized coordinate.
 */
template <typename T>
class gen_coord : public shared_object {
 public:
  using self = gen_coord<T>;
  using value_type = T;

  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;

  value_type q;       ///< Position value of the generalized coordinate.
  value_type q_dot;   ///< Velocity value of the generalized coordinate.
  value_type q_ddot;  ///< Acceleration value of the generalized coordinate.
  value_type f;       ///< Force value of the generalized coordinate.

  /**
   * Default constructor, all is set to zero.
   */
  gen_coord() : shared_object(), q(0.0), q_dot(0.0), q_ddot(0.0), f(0.0) {}

  /**
   * Parametrized constructor, all is set to corresponding parameters.
   */
  gen_coord(const_reference Q, const_reference Q_dot, const_reference Q_ddot,
            const_reference F)
      : shared_object(), q(Q), q_dot(Q_dot), q_ddot(Q_ddot), f(F) {}

  /**
   * Default virtual destructor.
   */
  ~gen_coord() override = default;

  /**
   * Add Q to the position value.
   */
  self& add_Q(const_reference Q) {
    q += Q;
    return *this;
  }

  /**
   * Add Q_dot to the velocity value.
   */
  self& add_Q_dot(const_reference Q_dot) {
    q_dot += Q_dot;
    return *this;
  }

  /**
   * Add Q_ddot to the acceleration value.
   */
  self& add_Q_ddot(const_reference Q_ddot) {
    q_ddot += Q_ddot;
    return *this;
  }

  /**
   * Add F to the force value.
   */
  self& add_F(const_reference F) {
    f += F;
    return *this;
  }

  /**
   * Addition operator.
   */
  friend self operator+(const self& G1, const self& G2) {
    return self(G1.q + G2.q, G1.q_dot + G2.q_dot, G1.q_ddot + G2.q_ddot,
                G1.f + G2.f);
  }

  /**
   * Substraction operator.
   */
  friend self operator-(const self& G1, const self& G2) {
    return self(G1.q - G2.q, G1.q_dot - G2.q_dot, G1.q_ddot - G2.q_ddot,
                G1.f - G2.f);
  }

  /**
   * Negation operator.
   */
  self operator-() const { return self(-q, -q_dot, -q_ddot, -f); }

  /**
   * Addition-assignment operator.
   */
  self& operator+=(const self& G) {
    q += G.q;
    q_dot += G.q_dot;
    q_ddot += G.q_ddot;
    f += G.f;
    return *this;
  }

  /**
   * Substraction-assignment operator.
   */
  self& operator-=(const self& G) {
    q -= G.q;
    q_dot -= G.q_dot;
    q_ddot -= G.q_ddot;
    f -= G.f;
    return *this;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(q) & RK_SERIAL_SAVE_WITH_NAME(q_dot) &
        RK_SERIAL_SAVE_WITH_NAME(q_ddot) & RK_SERIAL_SAVE_WITH_NAME(f);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(q) & RK_SERIAL_LOAD_WITH_NAME(q_dot) &
        RK_SERIAL_LOAD_WITH_NAME(q_ddot) & RK_SERIAL_LOAD_WITH_NAME(f);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x0000000F, 1, "gen_coord", shared_object)
};

template <typename T>
gen_coord<T> gen_coord_pos(const T& value) {
  return gen_coord<T>(value, 0, 0, 0);
}

template <typename T>
gen_coord<T> gen_coord_vel(const T& value) {
  return gen_coord<T>(0, value, 0, 0);
}

template <typename T>
gen_coord<T> gen_coord_acc(const T& value) {
  return gen_coord<T>(0, 0, value, 0);
}

template <typename T>
gen_coord<T> gen_coord_force(const T& value) {
  return gen_coord<T>(0, 0, 0, value);
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const gen_coord<T>& g) {
  out << "(q = " << g.q << "; q_dot = " << g.q_dot << "; q_ddot = " << g.q_ddot
      << "; f = " << g.f << ")";
  return out;
}

}  // namespace ReaK

#endif

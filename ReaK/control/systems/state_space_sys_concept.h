/**
 * \file state_space_sys_concept.h
 *
 * This library defines the traits class for a state-space system and the concept
 * which a state-space system can model. A state-space system is essentially a
 * class which corresponds to the mathematical concept of a state-space system, that is,
 * is can map the current state, time and input to a state-derivative vector as well as
 * compute an output.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2011
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

#ifndef REAK_CONTROL_SYSTEMS_STATE_SPACE_SYS_CONCEPT_H_
#define REAK_CONTROL_SYSTEMS_STATE_SPACE_SYS_CONCEPT_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"

#include <concepts>

namespace ReaK::ctrl {

/**
 * This class template is the traits class that defines the traits that a
 * state-space system should have.
 * \tparam SSSystem The state-space system type for which the traits are sought.
 */
template <typename T>
struct ss_system_traits {
  /** This type is the state-vector type, i.e., a descriptor of the state. */
  using point_type = typename T::point_type;
  /** This type is the state-difference type, i.e., a descriptor of the difference between states. */
  using point_difference_type = typename T::point_difference_type;
  /** This type is the state-derivative type, i.e., a descriptor of the state's derivative. */
  using point_derivative_type = typename T::point_derivative_type;

  /** This type is the time type. */
  using time_type = typename T::time_type;
  /** This type is the time-difference type. */
  using time_difference_type = typename T::time_difference_type;

  /** This type is the input type, i.e., a descriptor of the system's input. */
  using input_type = typename T::input_type;
  /** This type is the output type, i.e., a descriptor of the system's output. */
  using output_type = typename T::output_type;

  /** This constant describes the dimensions of the state-space (0 if not known at compile-time). */
  static constexpr std::size_t dimensions = T::dimensions;
  /** This constant describes the dimensions of the input vector (0 if not known at compile-time). */
  static constexpr std::size_t input_dimensions = T::input_dimensions;
  /** This constant describes the dimensions of the output vector (0 if not known at compile-time). */
  static constexpr std::size_t output_dimensions = T::output_dimensions;
};

template <typename T>
struct ss_system_traits<const T> : ss_system_traits<T> {};

template <typename T>
struct ss_system_traits<T&> : ss_system_traits<T> {};

/**
 * This class template defines the concept for a state-space system as used in the ReaK::ctrl
 * library. In addition to providing the traits defined in ss_system_traits, a state-space
 * system should provide a number of valid expressions.
 *
 * Valid expressions:
 *
 * dp = -dp;  A state-difference can be negated.
 *
 * t = t + dt;  A time-difference can be added to a time value.
 *
 * dp = dp_dt * dt;  A state-derivative times a time-difference yields a state-difference.
 *
 * dp_dt = sys.get_state_derivative(state_space,p,u,t);  The state-space system (sys) can compute the state-derivative
 *given the current state (p), the current input (u) and the current time (t).
 *
 * y = sys.get_output(state_space,p,u,t);  The state-space system (sys) can compute the output (y) given the current
 *state (p), the current input (u) and the current time (t).
 *
 * s = sys.get_state_dimensions();  The state-space system (sys) can deliver the dimensions count (s) for the states of
 *the system.
 *
 * i = sys.get_input_dimensions();  The state-space system (sys) can deliver the dimensions count (i) for the inputs of
 *the system.
 *
 * o = sys.get_output_dimensions();  The state-space system (sys) can deliver the dimensions count (o) for the outputs
 *of the system.
 */
template <typename T, typename StateSpace>
concept SSSystem = pp::Topology<StateSpace>&& requires(
    T sys, StateSpace space, typename ss_system_traits<T>::point_type p,
    typename ss_system_traits<T>::input_type u,
    typename ss_system_traits<T>::time_type t,
    typename ss_system_traits<T>::point_difference_type dp,
    typename ss_system_traits<T>::time_difference_type dt,
    typename ss_system_traits<T>::point_derivative_type dp_dt,
    typename ss_system_traits<T>::output_type y) {
  dp = -dp;
  t = t + dt;
  dp = dp_dt * dt;  // state-space system requirements
  dp_dt = sys.get_state_derivative(space, p, u, t);
  y = sys.get_output(space, p, u, t);
  { sys.get_state_dimensions() } -> std::integral;
  { sys.get_input_dimensions() } -> std::integral;
  { sys.get_output_dimensions() } -> std::integral;
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_SYSTEMS_STATE_SPACE_SYS_CONCEPT_H_

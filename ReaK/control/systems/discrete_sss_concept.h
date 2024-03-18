/**
 * \file discrete_sss_concept.h
 *
 * This library defines the traits class and concept that represent a discrete-time
 * state-space system. This type of systems is the bread-and-butter of state-estimation
 * algorithms, because all useful estimators apply to discrete-time systems.
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

#ifndef REAK_CONTROL_SYSTEMS_DISCRETE_SSS_CONCEPT_H_
#define REAK_CONTROL_SYSTEMS_DISCRETE_SSS_CONCEPT_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"

#include <concepts>

namespace ReaK::ctrl {

/**
 * This traits class defines the characteristics of a discrete-time state-space system.
 */
template <typename T>
struct discrete_sss_traits {
  /** The type which describes the state of the system. */
  using point_type = typename T::point_type;
  /** The type which describes the difference between two states of the system. */
  using point_difference_type = typename T::point_difference_type;

  /** The type which describes the time. */
  using time_type = typename T::time_type;
  /** The type which describes a time difference. */
  using time_difference_type = typename T::time_difference_type;

  /** The type which describes the input vector to the system. */
  using input_type = typename T::input_type;
  /** The type which describes the output of the system. */
  using output_type = typename T::output_type;

  /** This constant describes the dimensions of the state vector (0 if not known at compile-time). */
  static constexpr std::size_t dimensions = T::dimensions;
  /** This constant describes the dimensions of the input vector (0 if not known at compile-time). */
  static constexpr std::size_t input_dimensions = T::input_dimensions;
  /** This constant describes the dimensions of the output vector (0 if not known at compile-time). */
  static constexpr std::size_t output_dimensions = T::output_dimensions;
};

/**
 * This concept class template defines the requirements for a type to be a discrete-time
 * state-space system, as used in ReaK::ctrl.
 *
 * Valid expressions:
 *
 * dt = sys.get_time_step();  The time-step of the discrete-time system can be obtained.
 *
 * p = sys.get_next_state(state_space,p,u,t);  The next state (p) can be obtained from the current state (p), current
 *input (u) and current time (t).
 *
 * y = sys.get_output(state_space,p,u,t);  The system's output (y) can be obtained from the state (p), input (u) and
 *time (t).
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
concept DiscreteSSS = pp::Topology<StateSpace> &&
  requires (const T& sys, const StateSpace& space,
            typename discrete_sss_traits<T>::point_type p,
            typename discrete_sss_traits<T>::input_type u,
            typename discrete_sss_traits<T>::time_type t,
            typename discrete_sss_traits<T>::time_difference_type dt,
            typename discrete_sss_traits<T>::output_type y) {
    dt = sys.get_time_step();
    p = sys.get_next_state(space, p, u, t);
    y = sys.get_output(space, p, u, t);
    { sys.get_state_dimensions() } -> std::integral;
    { sys.get_input_dimensions() } -> std::integral;
    { sys.get_output_dimensions() } -> std::integral;
  };

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_SYSTEMS_DISCRETE_SSS_CONCEPT_H_

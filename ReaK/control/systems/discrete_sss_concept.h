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

#include "boost/concept_check.hpp"

#include "ReaK/topologies/spaces/metric_space_concept.h"

namespace ReaK::ctrl {

/**
 * This traits class defines the characteristics of a discrete-time state-space system.
 * \tparam DiscreteSystem The discrete-time state-space system type for which the traits are sought.
 */
template <typename DiscreteSystem>
struct discrete_sss_traits {
  /** The type which describes the state of the system. */
  using point_type = typename DiscreteSystem::point_type;
  /** The type which describes the difference between two states of the system. */
  using point_difference_type = typename DiscreteSystem::point_difference_type;

  /** The type which describes the time. */
  using time_type = typename DiscreteSystem::time_type;
  /** The type which describes a time difference. */
  using time_difference_type = typename DiscreteSystem::time_difference_type;

  /** The type which describes the input vector to the system. */
  using input_type = typename DiscreteSystem::input_type;
  /** The type which describes the output of the system. */
  using output_type = typename DiscreteSystem::output_type;

  /** This constant describes the dimensions of the state vector (0 if not known at compile-time). */
  static constexpr std::size_t dimensions = DiscreteSystem::dimensions;
  /** This constant describes the dimensions of the input vector (0 if not known at compile-time). */
  static constexpr std::size_t input_dimensions =
      DiscreteSystem::input_dimensions;
  /** This constant describes the dimensions of the output vector (0 if not known at compile-time). */
  static constexpr std::size_t output_dimensions =
      DiscreteSystem::output_dimensions;
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
 *
 * \tparam DiscreteSystem The type to be tested for being a discrete-time state-space system.
 * \tparam StateSpaceType The type of the state-space topology on which the state-space system should be able to act.
 */
template <typename DiscreteSystem, typename StateSpaceType>
struct DiscreteSSSConcept {
  DiscreteSystem sys;
  StateSpaceType state_space;
  typename discrete_sss_traits<DiscreteSystem>::point_type p;
  typename discrete_sss_traits<DiscreteSystem>::time_type t;
  typename discrete_sss_traits<DiscreteSystem>::time_difference_type dt;
  typename discrete_sss_traits<DiscreteSystem>::input_type u;
  typename discrete_sss_traits<DiscreteSystem>::output_type y;

  BOOST_CONCEPT_ASSERT((pp::TopologyConcept<StateSpaceType>));

  BOOST_CONCEPT_USAGE(DiscreteSSSConcept) {
    dt = sys.get_time_step();
    p = sys.get_next_state(state_space, p, u, t);
    y = sys.get_output(state_space, p, u, t);
    std::size_t s = sys.get_state_dimensions();
    RK_UNUSED(s);
    std::size_t i = sys.get_input_dimensions();
    RK_UNUSED(i);
    std::size_t o = sys.get_output_dimensions();
    RK_UNUSED(o);
  }
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_SYSTEMS_DISCRETE_SSS_CONCEPT_H_

/**
 * \file ss_controller_concept.hpp
 *
 * This library defines the concept-check class
 * for state-space controllers. A state-space controller is essentially a
 * class which corresponds to the mathematical concept of a state-space control system, that is,
 * it can map the current state and time of a system into an input vector to apply on it.
 *
 * A state-space controller shares the same traits as a state-space system (see ss_system_traits).
 * A stateful controller is, in essence, also a state-space system whose input vector is the plant
 * state vector, the output is the plant input vector, and the controller state is instrinsic to the
 * controller. For a stateless controller, it is essentially just a function to compute the plant-input
 * given a plant-state.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date August 2012
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

#ifndef REAK_SS_CONTROLLER_CONCEPT_HPP
#define REAK_SS_CONTROLLER_CONCEPT_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/topologies/spaces/metric_space_concept.hpp>

#include <ReaK/control/systems/discrete_linear_sss_concept.hpp>
#include <ReaK/control/systems/discrete_sss_concept.hpp>
#include <ReaK/control/systems/linear_ss_system_concept.hpp>
#include <ReaK/control/systems/state_space_sys_concept.hpp>

#include <boost/concept_check.hpp>

namespace ReaK::ctrl {

/**
 * This class is a constraint specifier for the controller concept classes. This class template
 * specifies that the controller is a stateful system (either discrete-time or continuous-time),
 * whose state space is that given as argument. This specifier essentially makes the control system
 * a complete system in its own right, meaning that additional system concepts may apply, such as
 * linear system concepts (see LinearSSSystemConcept or DiscreteLinearSSSConcept).
 *
 * Valid expressions:
 *
 * The constraints are that the StateSpaceType models the pp::TopologyConcept, and that
 * the system to which this constraint is applied models either the SSSystemConcept or
 * the DiscreteSSSConcept, for the continuous-time or discrete-time cases, respectively.
 *
 * \tparam StateSpaceType The type of the state-space for the controller (or system in general).
 */
template <typename StateSpaceType>
struct Stateful {
  StateSpaceType state_space;

  typename pp::topology_traits<StateSpaceType>::point_type p;

  BOOST_CONCEPT_ASSERT((pp::TopologyConcept<StateSpaceType>));

  template <typename System, typename Output, typename Input, typename Time>
  void ct_constraints(const System& sys, Output& y, const Input& u,
                      const Time& t) {
    BOOST_CONCEPT_ASSERT((SSSystemConcept<System, StateSpaceType>));
    y = sys.get_output(state_space, p, u, t);
  }

  template <typename System, typename Output, typename Input, typename Time>
  void dt_constraints(const System& sys, Output& y, const Input& u,
                      const Time& t) {
    BOOST_CONCEPT_ASSERT((DiscreteSSSConcept<System, StateSpaceType>));
    y = sys.get_output(state_space, p, u, t);
  }
};

/**
 * This class is a constraint specifier for the controller concept classes. This class template
 * specifies that the controller is a stateless system, i.e., just a functional mapping of input
 * to output. The discrete-time and continous-time requirements are the same.
 *
 * Valid expressions:
 *
 * y = sys.get_output(u,t);  The output vector (plant-input) can be obtained from the input vector (plant-state) and
 *time.
 */
struct Stateless {
  template <typename System, typename Output, typename Input, typename Time>
  void ct_constraints(const System& sys, Output& y, const Input& u,
                      const Time& t) {
    y = sys.get_output(u, t);
  }

  template <typename System, typename Output, typename Input, typename Time>
  void dt_constraints(const System& sys, Output& y, const Input& u,
                      const Time& t) {
    y = sys.get_output(u, t);
  }
};

/**
 * This class template defines the concept for a continuous-time state-space controller as used in the ReaK::ctrl
 * library. In addition to providing the traits defined in ss_system_traits, a continuous-time state-space controller
 * should provide a number of valid expressions.
 *
 * Valid expressions:
 *
 * For a 'Stateful' controller:
 *
 * Requires that CtrlSystem models the SSSystemConcept for the given state-space type in the 'Stateful' instance.
 *
 * For a 'Stateless' controller:
 *
 * y = sys.get_output(u,t);  The output vector (plant-input) can be obtained from the input vector (plant-state) and
 *time.
 *
 * \tparam CtrlSystem The state-space control-system type which is tested for modeling the continuous-time state-space
 *controller concept.
 * \tparam PlantSystem The state-space plant system type for which the controller is for.
 * \tparam Statefulness A type specifying the statefulness required of the state-space controller (see Stateless or
 *Stateful).
 */
template <typename CtrlSystem, typename PlantSystem, typename Statefulness>
struct CTSSControllerConcept {
  CtrlSystem ctrl_sys;
  Statefulness statefulness_constraint;
  typename ss_system_traits<CtrlSystem>::time_type t;
  typename ss_system_traits<PlantSystem>::point_type u;
  typename ss_system_traits<PlantSystem>::input_type y;

  BOOST_CONCEPT_USAGE(CTSSControllerConcept) {
    statefulness_constraint.ct_constraints(ctrl_sys, y, u, t);
  }
};

/**
 * This class template defines the concept for a discrete-time state-space controller as used in the ReaK::ctrl
 * library. In addition to providing the traits defined in ss_system_traits, a discrete-time state-space controller
 * should provide a number of valid expressions.
 *
 * Valid expressions:
 *
 * For a 'Stateful' controller:
 *
 * Requires that CtrlSystem models the DiscreteSSSConcept for the given state-space type in the 'Stateful' instance.
 *
 * For a 'Stateless' controller:
 *
 * y = sys.get_output(u,t);  The output vector (plant-input) can be obtained from the input vector (plant-state) and
 *time.
 *
 * \tparam CtrlSystem The state-space control-system type which is tested for modeling the discrete-time state-space
 *controller concept.
 * \tparam PlantSystem The state-space plant system type for which the controller is for.
 * \tparam Statefulness A type specifying the statefulness required of the state-space controller (see Stateless or
 *Stateful).
 */
template <typename CtrlSystem, typename PlantSystem, typename Statefulness>
struct DTSSControllerConcept {
  CtrlSystem ctrl_sys;
  Statefulness statefulness_constraint;
  typename discrete_sss_traits<CtrlSystem>::time_type t;
  typename discrete_sss_traits<PlantSystem>::point_type u;
  typename discrete_sss_traits<PlantSystem>::input_type y;

  BOOST_CONCEPT_USAGE(DTSSControllerConcept) {
    statefulness_constraint.dt_constraints(ctrl_sys, y, u, t);
  }
};

}  // namespace ReaK::ctrl

#endif

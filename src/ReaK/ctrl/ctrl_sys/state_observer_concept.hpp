/**
 * \file state_observer_concept.hpp
 * 
 * This library defines the concept-check class 
 * for state-controller observers. A state-observer is essentially a 
 * class which corresponds to the mathematical concept of a state-observer system, that is,
 * it can map the current output and time of a plant into a state-vector that is an estimate
 * of the state of that plant.
 * 
 * A state-observer shares the same traits as a state-space system (see ss_system_traits).
 * A stateful observer is, in essence, also a state-space system whose input vector is the plant 
 * output vector, the output is the plant state vector, and the observer state is instrinsic to the 
 * observer. For a stateless observer, it is essentially just a function to compute the plant-state
 * given a plant-output.
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

#ifndef REAK_STATE_OBSERVER_CONCEPT_HPP
#define REAK_STATE_OBSERVER_CONCEPT_HPP

#include <ReaK/ctrl/path_planning/metric_space_concept.hpp>

#include "state_space_sys_concept.hpp"
#include "linear_ss_system_concept.hpp"
#include "discrete_sss_concept.hpp"
#include "discrete_linear_sss_concept.hpp"
#include "ss_controller_concept.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Control */
namespace ctrl {


/**
 * This class template defines the concept for a continuous-time state-observer as used in the ReaK::ctrl
 * library. In addition to providing the traits defined in ss_system_traits, a continuous-time state-observer 
 * should provide a number of valid expressions.
 * 
 * Valid expressions:
 * 
 * For a 'Stateful' observer:
 * 
 * Requires that ObserverSystem models the SSSystemConcept for the given state-space type in the 'Stateful' instance.
 * 
 * For a 'Stateless' controller:
 * 
 * y = sys.get_output(u,t);  The output vector (plant-state) can be obtained from the input vector (plant-output) and time.
 * 
 * \tparam ObserverSystem The state-observer system type which is tested for modeling the continuous-time state-observer concept.
 * \tparam PlantSystem The state-space plant system type for which the observer is for.
 * \tparam Statefulness A type specifying the statefulness required of the state-observer (see Stateless or Stateful).
 */
template <typename ObserverSystem, typename PlantSystem, typename Statefulness>
struct CTSSObserverConcept {
  ObserverSystem obs_sys;
  Statefulness statefulness_constraint;
  typename ss_system_traits<ObserverSystem>::time_type t;
  typename ss_system_traits<PlantSystem>::output_type u;
  typename ss_system_traits<PlantSystem>::point_type y;
  
  BOOST_CONCEPT_USAGE(CTSSObserverConcept)
  {
    statefulness_constraint.ct_constraints(obs_sys, y, u, t);
  };
  
};

/**
 * This class template defines the concept for a discrete-time state-observer as used in the ReaK::ctrl
 * library. In addition to providing the traits defined in ss_system_traits, a discrete-time state-observer 
 * should provide a number of valid expressions.
 * 
 * Valid expressions:
 * 
 * For a 'Stateful' observer:
 * 
 * Requires that ObserverSystem models the DiscreteSSSConcept for the given state-space type in the 'Stateful' instance.
 * 
 * For a 'Stateless' observer:
 * 
 * y = sys.get_output(u,t);  The output vector (plant-state) can be obtained from the input vector (plant-output) and time.
 * 
 * \tparam ObserverSystem The state-observer system type which is tested for modeling the discrete-time state-observer concept.
 * \tparam PlantSystem The state-space plant system type for which the observer is for.
 * \tparam Statefulness A type specifying the statefulness required of the state-observer (see Stateless or Stateful).
 */
template <typename ObserverSystem, typename PlantSystem, typename Statefulness>
struct DTSSObserverConcept {
  ObserverSystem obs_sys;
  Statefulness statefulness_constraint;
  typename discrete_sss_traits<ObserverSystem>::time_type t;
  typename discrete_sss_traits<PlantSystem>::output_type u;
  typename discrete_sss_traits<PlantSystem>::point_type y;
  
  BOOST_CONCEPT_USAGE(DTSSObserverConcept)
  {
    statefulness_constraint.dt_constraints(obs_sys, y, u, t);
  };
  
};



};

};

#endif






/**
 * \file discrete_sss_concept.hpp
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

#ifndef REAK_DISCRETE_SSS_CONCEPT_HPP
#define REAK_DISCRETE_SSS_CONCEPT_HPP

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "path_planning/metric_space_concept.hpp"

namespace ReaK {

namespace ctrl {


/**
 * This traits class defines the characteristics of a discrete-time state-space system.
 * \tparam DiscreteSystem The discrete-time state-space system type for which the traits are sought.
 */
template <typename DiscreteSystem>
struct discrete_sss_traits {
  /** The type which describes the state of the system. */
  typedef typename DiscreteSystem::point_type point_type;
  /** The type which describes the difference between two states of the system. */
  typedef typename DiscreteSystem::point_difference_type point_difference_type;
  
  /** The type which describes the time. */
  typedef typename DiscreteSystem::time_type time_type;
  /** The type which describes a time difference. */
  typedef typename DiscreteSystem::time_difference_type time_difference_type;
  
  /** The type which describes the input vector to the system. */
  typedef typename DiscreteSystem::input_type input_type;
  /** The type which describes the output of the system. */
  typedef typename DiscreteSystem::output_type output_type;
  
  /** This constant describes the dimensions of the state vector (0 if not known at compile-time). */
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = DiscreteSystem::dimensions);
  /** This constant describes the dimensions of the input vector (0 if not known at compile-time). */
  BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = DiscreteSystem::input_dimensions);
  /** This constant describes the dimensions of the output vector (0 if not known at compile-time). */
  BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = DiscreteSystem::output_dimensions);
  
};
  
  
/**
 * This concept class template defines the requirements for a type to be a discrete-time
 * state-space system, as used in ReaK::ctrl.
 * 
 * Valid expressions:
 * 
 * dt = sys.get_time_step();  The time-step of the discrete-time system can be obtained.
 * 
 * p = sys.get_next_state(state_space,p,u,t);  The next state (p) can be obtained from the current state (p), current input (u) and current time (t).
 * 
 * y = sys.get_output(state_space,p,u,t);  The system's output (y) can be obtained from the state (p), input (u) and time (t).
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
  
  BOOST_CONCEPT_USAGE(DiscreteSSSConcept)
  { 
    dt = sys.get_time_step();
    p = sys.get_next_state(state_space,p,u,t);
    y = sys.get_output(state_space,p,u,t);
  };
  
};





};

};

#endif











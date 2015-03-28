/**
 * \file steerable_space_concept.hpp
 *
 * This library defines the traits and concepts that pertain to what can be considered
 * a steerable-space, as used in ReaK::pp. Steerable-spaces are based on the Topology concept
 * from the Boost.Graph library, but with additional requirements which are needed
 * in algorithms tailored for path-planning (mostly). Basically, the concept of a
 * steerable-space corresponds to spaces in which non-trivial (steer) motions can
 * be done and their trace be recorded in a steering path or trajectory.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_STEERABLE_SPACE_CONCEPT_HPP
#define REAK_STEERABLE_SPACE_CONCEPT_HPP

#include <ReaK/core/base/defs.hpp>

#include <cmath>
#include <boost/concept_check.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/tuple/tuple.hpp>

#include "metric_space_concept.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {

/**
 * This traits class defines the types and constants associated to a steerable-space.
 * \tparam SteerableSpace The topology type for which the steerable-space traits are sought.
 */
template < typename SteerableSpace >
struct steerable_space_traits {
  /** The type that describes the distance-metric type for the space. */
  typedef typename SteerableSpace::steer_record_type steer_record_type;
};


/**
 * This concept defines the requirements to fulfill in order to model a steerable-space
 * as used in ReaK::pp. A steerable-space is a special kind of topology in which
 * non-trivial (steered) motions can be done and their trace be recorded in a steering path or trajectory.
 *
 * Valid expressions:
 *
 * tie(p2, st_rec) = space.steer_position_toward(p1, d, p3);  A point (p2) and a steer-record (st_rec) can be obtained
 *from attempting to steer a fraction (d) away from one point (p1) to another (p3).
 *
 * \tparam SteerableSpace The topology type to be checked for this concept.
 */
template < typename SteerableSpace >
struct SteerableSpaceConcept {
  typename topology_traits< SteerableSpace >::point_type p1, p2;
  typename steerable_space_traits< SteerableSpace >::steer_record_type st_rec;
  SteerableSpace space;
  double d;

  BOOST_CONCEPT_ASSERT( ( TopologyConcept< SteerableSpace > ) );

  BOOST_CONCEPT_USAGE( SteerableSpaceConcept ) { boost::tie( p1, st_rec ) = space.steer_position_toward( p1, d, p2 ); };
};


template < typename SteerableSpace >
struct is_steerable_space : boost::mpl::false_ {};
};
};


#endif

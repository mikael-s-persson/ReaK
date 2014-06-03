/**
 * \file reachability_space_concept.hpp
 * 
 * This library defines the traits and concept class that represent the idea of a 
 * reachability space (or topology). Essentially, a reachability space is a temporal 
 * space for which the spatial metric and time metric are collapsed into a dual metric
 * which is essentially expressed in a temporal unit that roughly represents the time 
 * it takes to reach a given point. This representation of reachability allows for a 
 * dual sorting based on the reachability with respect to the origin, which can 
 * make nearest reachable neighbor queries very efficient (see reachability_sort.hpp).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2011
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

#ifndef REAK_REACHABILITY_SPACE_CONCEPT_HPP
#define REAK_REACHABILITY_SPACE_CONCEPT_HPP

#include "temporal_space_concept.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace pp {

/**
 * This traits class defines the characteristics of a reachability space (or topology).
 * Since a reachability-space is typically implemented as a thin wrapper for a temporal-space
 * or an extension of a temporal-space, the traits are very similary to those of a temporal-space.
 * \tparam ReachabilityTopology The reachability-space type for which the traits are sought.
 */
template <typename ReachabilityTopology>
struct reachability_topology_traits {
  /** The point type that described the points of the reachability-space. */
  typedef typename ReachabilityTopology::point_type point_type;
  /** The point-difference type that described the difference between points of the reachability-space. */
  typedef typename ReachabilityTopology::point_difference_type point_difference_type;
 
  /** The temporal-topology type in which the points of the reachability-space lie. */
  typedef typename ReachabilityTopology::temporal_space_type temporal_space_type;
  /** The time-topology type in which the time-components of points of the reachability-space lie. */
  typedef typename ReachabilityTopology::time_topology time_topology;
  /** The space-topology type in which the space-components of points of the reachability-space lie. */
  typedef typename ReachabilityTopology::space_topology space_topology;
  /** The distance-metric type for the reachability-space. */
  typedef typename ReachabilityTopology::distance_metric distance_metric;
  
  /** This constant defines the temporal dimensions (0 if unknown at compile-time). */
  BOOST_STATIC_CONSTANT(std::size_t, time_dimensions = time_topology::dimensions);
  /** This constant defines the spatial dimensions (0 if unknown at compile-time). */
  BOOST_STATIC_CONSTANT(std::size_t, space_dimensions = space_topology::point_type::dimensions);
  
};


/**
 * This concept class defines the idea of a 
 * reachability space (or topology). Essentially, a reachability space is a temporal 
 * space for which the spatial metric and time metric are collapsed into a dual metric
 * which is essentially expressed in a temporal unit that roughly represents the time 
 * it takes to reach a given point. This representation of reachability allows for a 
 * dual sorting based on the reachability with respect to the origin, which can 
 * make nearest reachable neighbor queries very efficient (see reachability_sort.hpp).
 * 
 * The dual metric (which is not, strictly a metric as is) has a backward and forward
 * part. Assuming that a point-difference is taken as the difference between two point (A,B),
 * then the backward norm of that difference should, loosely-speaking, represent the 
 * time that can be spent at spatial-point A before point B becomes unreachable. Conversely, 
 * the forward norm should represent the soonest time, from temporal-point A, at which 
 * the spatial-point A can be reached from point B. In other words, if one has a spatial-norm
 * which can describe the minimum travel time required to move a given distance, then the 
 * backward norm can be computed by time-difference minus travel-time and the forward norm
 * can be computed by time-difference plus travel-time, this is what the backward_reachable_norm 
 * and forward_reachable_norm functors compute.
 * 
 * Required concepts:
 * 
 * The topology should model the TemporalSpaceConcept.
 * 
 * Valid expressions:
 * 
 * d = reachable_space.forward_reach(p1);  The forward-norm of a point (p1) with respect to some origin can be obtained.
 * 
 * d = reachable_space.backward_reach(p1);  The backward-norm of a point (p1) with respect to some origin can be obtained.
 * 
 * d = reachable_space.forward_norm(pd);  The forward-norm of a point-difference (pd) can be obtained.
 * 
 * d = reachable_space.backward_norm(pd);  The backward-norm of a point-difference (pd) can be obtained.
 * 
 * \tparam ReachabilityTopology The topology type for which this concept is checked.
 */
template <typename ReachabilityTopology>
struct ReachabilitySpaceConcept {
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept< typename reachability_topology_traits<ReachabilityTopology>::temporal_space_type >));
  
  typename reachability_topology_traits<ReachabilityTopology>::point_type p1;
  typename reachability_topology_traits<ReachabilityTopology>::point_difference_type pd;
  ReachabilityTopology reachable_space;
  double d;
  
  BOOST_CONCEPT_USAGE(ReachabilitySpaceConcept)
  {
    d = reachable_space.forward_reach(p1);
    d = reachable_space.backward_reach(p1);
    d = reachable_space.forward_norm(pd);
    d = reachable_space.backward_norm(pd);
  };
  
};


};

};

#endif




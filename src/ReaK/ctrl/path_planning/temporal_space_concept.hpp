/**
 * \file temporal_space_concept.hpp
 * 
 * This library defines the traits and concepts related to the construction of a 
 * temporal space, i.e. a topology which is constituted of a spatial metric-space 
 * and a time topology (also a metric-space). The temporal space is also a metric-space.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2011
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

#ifndef REAK_TEMPORAL_SPACE_CONCEPT_HPP
#define REAK_TEMPORAL_SPACE_CONCEPT_HPP

#include <boost/config.hpp>
#include <boost/concept_check.hpp>

#include "metric_space_concept.hpp"

namespace ReaK {

namespace pp {


/**
 * This traits class defines the characteristics associated to a temporal-space type.
 * \tparam TemporalTopology The temporal-space type for which the traits are sought.
 */
template <typename TemporalTopology>
struct temporal_topology_traits {
  /** The type that describes a point in the space. */
  typedef typename TemporalTopology::point_type point_type;
  /** The type that describes a difference between points in the space. */
  typedef typename TemporalTopology::point_difference_type point_difference_type;
  
  /** The topology type which describes the space in which the time values reside. */
  typedef typename TemporalTopology::time_topology time_topology;
  /** The topology type which describes the space in which the spatial points reside. */
  typedef typename TemporalTopology::space_topology space_topology;
  
};
  

/**
 * This concept defines the requirements to fulfill in order to model a temporal distance-metric 
 * as used in ReaK::pp. A temporal distance-metric is essentially a callable type that can compute 
 * both the distance between two points and the corresponding norm of a difference between 
 * two points in a temporal-space.
 * 
 * Required concepts:
 * 
 * TemporalTopology should have the traits of a temporal space (temporal_topology_traits).
 * 
 * Valid expressions:
 * 
 * dist = d(p1, p2, t, s);  The distance (dist) can be obtained by calling the distance metric (d) on two points (p1,p2) and providing a const-ref to the time-topology (t) and space-topology (s).
 * 
 * dist = d(pd, t, s);  The distance (dist) can be obtained by calling the distance metric (d) on a point-difference (pd) and providing a const-ref to the time-topology (t) and space-topology (s).
 * 
 * \tparam DistanceMetric The distance metric type to be checked for this concept.
 * \tparam Topology The topology to which the distance metric should apply.
 */
template <typename DistanceMetric, typename TemporalTopology>
struct TemporalDistMetricConcept {
  DistanceMetric d;
  typename temporal_topology_traits<TemporalTopology>::space_topology s;
  typename temporal_topology_traits<TemporalTopology>::time_topology t;
  typename temporal_topology_traits<TemporalTopology>::point_type p1, p2;
  typename temporal_topology_traits<TemporalTopology>::point_difference_type pd;
  double dist;
  BOOST_CONCEPT_USAGE(TemporalDistMetricConcept)
  {
    dist = d(p1, p2, t, s);
    dist = d(pd, t, s);
  };
  
};

/**
 * This concept defines the requirements to fulfill in order to model a temporal-space 
 * as used in ReaK::pp. A temporal space is constituted of a spatial metric-space 
 * and a time topology (also a metric-space). The temporal space is also a metric-space.
 * 
 * Required concepts:
 * 
 * The space-topology should model the MetricSpaceConcept.
 * 
 * The time-topology should model the MetricSpaceConcept.
 * 
 * The temporal-topology should model the MetricSpaceConcept.
 * 
 * Valid expressions:
 * 
 * See MetricSpaceConcept.
 * 
 * \tparam Topology The topology type to be checked for this concept.
 */
template <typename Topology>
struct TemporalSpaceConcept : public MetricSpaceConcept< Topology > {
  BOOST_CONCEPT_ASSERT((MetricSpaceConcept< typename temporal_topology_traits<Topology>::space_topology >));
  BOOST_CONCEPT_ASSERT((MetricSpaceConcept< typename temporal_topology_traits<Topology>::time_topology >));
  
  BOOST_CONCEPT_USAGE(TemporalSpaceConcept)
  {
  };
  
};


};

};

#endif



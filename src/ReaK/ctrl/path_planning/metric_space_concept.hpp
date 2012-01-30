/**
 * \file metric_space_concept.hpp
 * 
 * This library defines the traits and concepts that pertain to what can be considered 
 * a metric-space, as used in ReaK::pp. Metric-spaces are based on the Topology concept 
 * from the Boost.Graph library, but with additional requirements which are needed 
 * in algorithms tailored for a metric-space (see metric_space_search.hpp). Basically,
 * the concept of a metric-space in ReaK::pp corresponds to the mathematical concept of 
 * a metric-space (see wikipedia or any decent math book).
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

#ifndef REAK_METRIC_SPACE_CONCEPT_HPP
#define REAK_METRIC_SPACE_CONCEPT_HPP


#include <boost/config.hpp>
#include <cmath>
#include <boost/concept_check.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {
  
  
  
/**
 * This traits class defines the types and constants associated to a topology.
 * \tparam Topology The topology type for which the topology traits are sought.
 */
template <typename Topology>
struct topology_traits {
  /** The type that describes a point in the space. */
  typedef typename Topology::point_type point_type;
  /** The type that describes a difference between points in the space. */
  typedef typename Topology::point_difference_type point_difference_type;
  
  /** The dimensions of the space (0 if unknown at compile-time). */
  BOOST_STATIC_CONSTANT(std::size_t, dimensions = Topology::dimensions);
  
};


/**
 * This concept defines the requirements to fulfill in order to model a topology 
 * as used in ReaK::pp.
 * 
 * Valid expressions:
 * 
 * dp = space.difference(p1,p2);  The difference (pd) between two points (p1,p2) can be obtained.
 * 
 * p1 = space.origin();  The origin of the space can be obtained.
 * 
 * p2 = space.adjust(p1,dp);  A point-difference can be scaled (d * pd), added / subtracted to another point-difference and added to a point (p1) to obtain an adjusted point.
 * 
 * \tparam Topology The topology type to be checked for this concept.
 */
template <typename Topology>
struct TopologyConcept {
  typename topology_traits<Topology>::point_type p1, p2;
  typename topology_traits<Topology>::point_difference_type dp;
  Topology space;
  
  BOOST_CONCEPT_USAGE(TopologyConcept) 
  {
    dp = space.difference(p1,p2);
    p1 = space.origin();
    p1 = space.adjust(p1,dp);
  };
  
};


/**
 * This concept defines the requirements to fulfill in order to model a Lie Group
 * as used in ReaK::pp. Basically, a Lie Group is a topology on which the point-difference
 * type is an arithmetic type (i.e. vector-space).
 * 
 * Valid expressions:
 * 
 * dp = d * dp + dp - dp;  The differences can be added, subtracted, and multiplied by a scalar.
 * 
 * dp = -dp;  The differences can be reversed.
 * 
 * dp -= dp;  The differences can be subtracted-and-stored.
 * 
 * dp += dp;  The differences can be added-and-stored.
 * 
 * dp *= d;  The differences can be multiplied-and-stored by a scalar.
 * 
 * \tparam LieGroup The Lie Group type to be checked for this concept.
 */
template <typename LieGroup>
struct LieGroupConcept {
  
  BOOST_CONCEPT_ASSERT((TopologyConcept<LieGroup>));
  
  typename topology_traits<LieGroup>::point_difference_type dp;
  double d;
  
  BOOST_CONCEPT_USAGE(LieGroupConcept) 
  {
    dp = d * dp + dp - dp;
    dp = -dp;
    dp -= dp;
    dp += dp;
    dp *= d;
  };
  
};
  
  
/**
 * This traits class defines the types and constants associated to a metric-space.
 * \tparam MetricSpace The topology type for which the metric-space traits are sought.
 */
template <typename MetricSpace>
struct metric_space_traits {
  /** The type that describes the distance-metric type for the space. */
  typedef typename MetricSpace::distance_metric distance_metric;
};

/**
 * This concept defines the requirements to fulfill in order to model a distance-metric 
 * as used in ReaK::pp. A distance-metric is essentially a callable type that can compute 
 * both the distance between two points and the corresponding norm of a difference between 
 * two points.
 * 
 * Required concepts:
 * 
 * Topology should model the TopologyConcept.
 * 
 * Valid expressions:
 * 
 * d = dist(p1, p2, space);  The distance (d) can be obtained by calling the distance metric (dist) on two points (p1,p2) and providing a const-ref to the topology (space).
 * 
 * d = dist(dp, space);  The distance (d) can be obtained by calling the distance metric (dist) on a point-difference (dp) and providing a const-ref to the topology (space).
 * 
 * \tparam DistanceMetric The distance metric type to be checked for this concept.
 * \tparam Topology The topology to which the distance metric should apply.
 */
template <typename DistanceMetric, typename Topology>
struct DistanceMetricConcept {
  DistanceMetric dist;
  Topology space;
  typename topology_traits<Topology>::point_type p1, p2;
  typename topology_traits<Topology>::point_difference_type dp;
  double d;
  
  BOOST_CONCEPT_USAGE(DistanceMetricConcept) 
  {
    d = dist(p1, p2, space);
    d = dist(dp, space);
  };
  
};

/**
 * This tag-type is used to identify (during a "get" call) that the distance-metric object is 
 * to be fetched.
 */
enum distance_metric_t { distance_metric };

/**
 * This concept defines the requirements to fulfill in order to model a metric-space 
 * as used in ReaK::pp. A metric-space is a special kind of topology which has a 
 * distance metric (in theory, satisfying triangular inequality).
 * 
 * Valid expressions:
 * 
 * dist = get(distance_metric, space);  The distance-metric can be obtained by a tagged "get" call on the metric-space.
 * 
 * p1 = space.move_position_toward(p1,d,p2);  A point can be obtained by moving a fraction (d) away from one point (p1) to another (p2).
 * 
 * \tparam MetricSpace The topology type to be checked for this concept.
 */
template <typename MetricSpace>
struct MetricSpaceConcept {
  typename topology_traits<MetricSpace>::point_type p1, p2;
  typename metric_space_traits<MetricSpace>::distance_metric dist;
  MetricSpace space;
  double d;
  
  BOOST_CONCEPT_ASSERT((TopologyConcept<MetricSpace>));
  BOOST_CONCEPT_ASSERT((DistanceMetricConcept<typename metric_space_traits<MetricSpace>::distance_metric>));
  
  BOOST_CONCEPT_USAGE(MetricSpaceConcept) 
  {
    dist = get(distance_metric, space);
    p1 = space.move_position_toward(p1,d,p2);
  };
  
};



};

};


#endif



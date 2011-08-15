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

#ifndef TEMPORAL_SPACES_HPP
#define TEMPORAL_SPACES_HPP

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
  typedef TemporalTopology::point_type point_type;
  /** The type that describes a difference between points in the space. */
  typedef TemporalTopology::point_difference_type point_difference_type;
  
  /** The topology type which describes the space in which the time values reside. */
  typedef TemporalTopology::time_topology time_topology;
  /** The topology type which describes the space in which the spatial points reside. */
  typedef TemporalTopology::space_topology space_topology;
  
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
  void constraints() {
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
 * Valid expressions:
 * 
 * d  = space.distance(p1, p2);  The distance between two points (p1,p2) can be obtained as a double (d).
 * 
 * d  = space.norm(pd);  The norm of the difference (pd) between two points can be obtained as a double (d).
 * 
 * p1 = space.random_point();  A random-point in the metric-space can be obtained.
 * 
 * pd = space.difference(p1,p2);  The difference (pd) between two points (p1,p2) can be obtained.
 * 
 * p1 = space.move_position_toward(p1,d,p2);  A point can be obtained by moving a fraction (d) away from one point (p1) to another (p2).
 * 
 * p1 = space.origin();  The origin of the space can be obtained.
 * 
 * p1 = space.adjust(p1,d * pd);  A point-difference can be scaled (d * pd) and added to a point (p1) to obtain an adjusted point.
 * 
 * pd = -pd;  A point-difference can be negated (reversed).
 * 
 * \tparam Topology The topology type to be checked for this concept.
 */
template <typename Topology>
struct TemporalSpaceConcept {
  typename temporal_topology_traits<Topology>::point_type p1, p2;
  typename temporal_topology_traits<Topology>::point_difference_type pd;
  Topology spacetime;
  double d;
  void constraints() {
    boost::function_requires< MetricSpaceConcept< typename temporal_topology_traits<Topology>::space_topology > >();
    boost::function_requires< MetricSpaceConcept< typename temporal_topology_traits<Topology>::time_topology > >();
    d  = spacetime.distance(p1, p2);
    d  = spacetime.norm(pd);
    p1 = spacetime.random_point();
    pd = spacetime.difference(p1,p2);
    p1 = spacetime.move_position_toward(p1,d,p2);
    p1 = spacetime.origin();
    p1 = spacetime.adjust(p1,d * pd);
    pd = -pd;
  };
  
};


/**
 * This class is a functor type which models the TemporalDistMetricConcept, and computes the 
 * distance based only on the distance in the spatial dimensions (space-topology).
 */
struct spatial_distance_only {
  /**
   * Computes the distance by calling the distance-function of the space-topology (s) on two points (a,b).
   * \tparam Point The point type of points on the temporal-space.
   * \tparam TimeTopology The time-topology type associated to the temporal-space.
   * \tparam SpaceTopology The space-topology type associated to the temporal-space.
   * \param a The first point.
   * \param b The second point.
   * \param s The space-topology.
   * \return the spatial-distance between the two points.
   */
  template <typename Point, typename TimeTopology, typename SpaceTopology>
  double operator()(const Point& a, const Point& b, const TimeTopology&, const SpaceTopology& s) const {
    return s.distance(a.pt, b.pt);
  };
  /**
   * Computes the norm by calling the norm-function of the space-topology (s) on a point-difference (a).
   * \tparam PointDiff The point-difference type of points on the temporal-space.
   * \tparam TimeTopology The time-topology type associated to the temporal-space.
   * \tparam SpaceTopology The space-topology type associated to the temporal-space.
   * \param a The point-difference.
   * \param s The space-topology.
   * \return The spatial-norm of the difference between the two points.
   */
  template <typename PointDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const PointDiff& a, const TimeTopology&, const SpaceTopology& s) const {
    return s.norm(a.pt);
  };
};


/**
 * This class is a functor type which models the TemporalDistMetricConcept, and computes the 
 * distance based only on the distance in the temporal dimensions (time-topology).
 */
struct time_distance_only {
  /**
   * Computes the distance by calling the distance-function of the time-topology (t) on two points (a,b).
   * \tparam Point The point type of points on the temporal-space.
   * \tparam TimeTopology The time-topology type associated to the temporal-space.
   * \tparam SpaceTopology The space-topology type associated to the temporal-space.
   * \param a The first point.
   * \param b The second point.
   * \param t The time-topology.
   * \return the temporal-distance between the two points.
   */
  template <typename Point, typename TimeTopology, typename SpaceTopology>
  double operator()(const Point& a, const Point& b, const TimeTopology& t, const SpaceTopology&) const {
    return t.distance(a.time, b.time);
  };
  /**
   * Computes the norm by calling the norm-function of the time-topology (t) on a point-difference (a).
   * \tparam PointDiff The point-difference type of points on the temporal-space.
   * \tparam TimeTopology The time-topology type associated to the temporal-space.
   * \tparam SpaceTopology The space-topology type associated to the temporal-space.
   * \param a The point-difference.
   * \param t The time-topology.
   * \return The temporal-norm of the difference between the two points.
   */
  template <typename PointDiff, typename TimeTopology, typename SpaceTopology>
  double operator()(const PointDiff& a, const TimeTopology& t, const SpaceTopology&) const {
    return t.norm(a.time);
  };
};


};

};

#endif



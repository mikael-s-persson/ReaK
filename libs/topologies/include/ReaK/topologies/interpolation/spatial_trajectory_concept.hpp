/**
 * \file spatial_trajectory_concept.hpp
 *
 * This library defines the traits and concepts related to a spatial trajectory. A
 * trajectory is simply a continuous curve in a temporal topology (or temporal-space) for which
 * each point is associated to a time (which increases as you travel along the trajectory).
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

#ifndef REAK_SPATIAL_TRAJECTORY_CONCEPT_HPP
#define REAK_SPATIAL_TRAJECTORY_CONCEPT_HPP

#include "ReaK/core/base/defs.hpp"

#include "ReaK/topologies/interpolation/spatial_path_concept.hpp"
#include "ReaK/topologies/spaces/temporal_space_concept.hpp"

#include <cmath>
#include "boost/concept_check.hpp"

namespace ReaK::pp {

/**
 * This traits class defines the traits that characterize a spatial trajectory
 * class.
 * \tparam SpatialTrajectory The trajectory for which the traits are sought.
 */
template <typename SpatialTrajectory>
struct spatial_trajectory_traits {
  /** This type describes a point in the space or topology. */
  using point_type = typename SpatialTrajectory::point_type;
  /** This type describes the difference between two points in the space or topology. */
  using point_difference_type =
      typename SpatialTrajectory::point_difference_type;

  /** This type describes waypoints used by the trajectory to quickly access local parameters of the trajectory. */
  using waypoint_descriptor = typename SpatialTrajectory::waypoint_descriptor;
  /** This type describes const-waypoints used by the trajectory to quickly access local parameters of the trajectory.
   */
  using const_waypoint_descriptor =
      typename SpatialTrajectory::const_waypoint_descriptor;

  /** This type is the temporal-topology type in which the trajectory exists. */
  using topology = typename SpatialTrajectory::topology;
  /** This type is the time-topology type in which the trajectory's time-stamps exist. */
  using time_topology = typename temporal_space_traits<topology>::time_topology;
  /** This type is the space-topology type in which the trajectory's points exist. */
  using space_topology =
      typename temporal_space_traits<topology>::space_topology;
  /** This type is the distance metric type used on the temporal-topology and defining the travel distances along the
   * trajectory. */
  using distance_metric = typename SpatialTrajectory::distance_metric;
};

/**
 * This concept class defines the requirements for a class to model a spatial-trajectory
 * concept as used in ReaK::pp. A spatial trajectory is a path whose points are associated
 * a time-stamp that progresses along the path. Points a given times can be obtained and
 * points can also be obtained at time-differences from a given point.
 *
 * Required concepts:
 *
 * The topology should model the TemporalSpaceConcept.
 *
 * The distance metric should model the TemporalDistMetricConcept over the given temporal-topology.
 *
 * Valid expressions:
 *
 * std::pair< const_waypoint_descriptor, point_type > w_p;  A waypoint-point pair is used to associate a point and a
 *waypoint.
 *
 * pt = p.move_time_diff_from(pt, dt);  A point along the trajectory (p), at a time-difference (dt) from a point (pt),
 *can be obtained.
 *
 * d = p.travel_distance(pt,pt);  The travel distance (as of the distance-metric), along the trajectory (p), between two
 *points (pt,pt), can be obtained.
 *
 * pt = p.get_point_at_time(t);  The point, along the trajectory (p), at a given time (t) can be obtained.
 *
 * w_p = p.move_time_diff_from(w_p, dt);  A waypoint along the trajectory (p), at a time-difference (dt) from a waypoint
 *(w_p), can be obtained.
 *
 * d = p.travel_distance(w_p,w_p);  The travel distance, along the trajectory (p), between two waypoints (w_p,w_p), can
 *be obtained.
 *
 * w_p = p.get_waypoint_at_time(t);  The waypoint, along the trajectory (p), at a given time (t) can be obtained.
 *
 * t = p.get_start_time();  The start time of the trajectory can be obtained.
 *
 * t = p.get_end_time();  The end time of the trajectory can be obtained.
 *
 * const topology& space = p.get_temporal_space();  The temporal space on which the trajectory lies can be obtained.
 *
 * \tparam SpatialTrajectory The trajectory type for which this concept is checked.
 * \tparam Topology The temporal-topology type on which the trajectory should be able to exist.
 */
template <typename SpatialTrajectory, typename Topology>
struct SpatialTrajectoryConcept {
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT(
      (DistanceMetricConcept<typename spatial_trajectory_traits<
                                 SpatialTrajectory>::distance_metric,
                             Topology>));

  SpatialTrajectory* p;
  typename temporal_space_traits<Topology>::point_type pt;
  std::pair<typename spatial_trajectory_traits<
                SpatialTrajectory>::const_waypoint_descriptor,
            typename temporal_space_traits<Topology>::point_type>
      w_p;
  using time_topology = typename temporal_space_traits<Topology>::time_topology;
  typename time_topology::point_difference_type dt;
  typename time_topology::point_type t;
  double d;

  BOOST_CONCEPT_USAGE(SpatialTrajectoryConcept) {
    pt = p->move_time_diff_from(pt, dt);
    d = p->travel_distance(pt, pt);
    pt = p->get_point_at_time(t);
    w_p = p->move_time_diff_from(w_p, dt);
    d = p->travel_distance(w_p, w_p);
    w_p = p->get_waypoint_at_time(t);

    t = p->get_start_time();
    t = p->get_end_time();

    const typename spatial_trajectory_traits<SpatialTrajectory>::topology&
        space = p->get_temporal_space();
    RK_UNUSED(space);
  }
};

}  // namespace ReaK::pp

#endif

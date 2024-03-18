/**
 * \file spatial_trajectory_concept.h
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_SPATIAL_TRAJECTORY_CONCEPT_H_
#define REAK_TOPOLOGIES_INTERPOLATION_SPATIAL_TRAJECTORY_CONCEPT_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/topologies/interpolation/spatial_path_concept.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/temporal_space_concept.h"

#include <cmath>
#include <concepts>

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
  using const_waypoint_descriptor = typename SpatialTrajectory::const_waypoint_descriptor;

  using const_waypoint_pair = std::pair<const_waypoint_descriptor, point_type>;

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
 * The topology should model TemporalSpace.
 *
 * The distance metric should model TemporalDistMetric over the given temporal-topology.
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
 */
template <typename Trajectory, typename Space = typename spatial_trajectory_traits<Trajectory>::topology>
concept SpatialTrajectory = TemporalSpace<Space> && DistanceMetric<typename spatial_trajectory_traits<Trajectory>::distance_metric, Space> &&
  requires (const Trajectory& traj, const topology_point_type_t<Space>& pt,
    topology_point_type_t<typename temporal_space_traits<Space>::time_topology> t,
    topology_point_difference_type_t<typename temporal_space_traits<Space>::time_topology> dt,
    const typename spatial_trajectory_traits<Trajectory>::const_waypoint_pair& wp) {
    { traj.move_time_diff_from(pt, dt) } -> std::convertible_to<topology_point_type_t<Space>>;
    { traj.travel_distance(pt, pt) } -> std::convertible_to<double>;
    { traj.get_point_at_time(t) } -> std::convertible_to<topology_point_type_t<Space>>;
    { traj.move_time_diff_from(wp, dt) } -> std::convertible_to<typename spatial_trajectory_traits<Trajectory>::const_waypoint_pair>;
    { traj.travel_distance(wp, wp) } -> std::convertible_to<double>;
    { traj.get_waypoint_at_time(t) } -> std::convertible_to<typename spatial_trajectory_traits<Trajectory>::const_waypoint_pair>;
    { traj.get_start_time() } -> std::convertible_to<topology_point_type_t<typename temporal_space_traits<Space>::time_topology>>;
    { traj.get_end_time() } -> std::convertible_to<topology_point_type_t<typename temporal_space_traits<Space>::time_topology>>;
    { traj.get_temporal_space() } -> std::convertible_to<const typename spatial_trajectory_traits<Trajectory>::topology&>;
  };

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_INTERPOLATION_SPATIAL_TRAJECTORY_CONCEPT_H_

/**
 * \file spatial_path_concept.h
 *
 * This library defines the traits and concepts related to a spatial path. A
 * path is simply a continuous curve in a topology (or space).
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_SPATIAL_PATH_CONCEPT_H_
#define REAK_TOPOLOGIES_INTERPOLATION_SPATIAL_PATH_CONCEPT_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"

#include <cmath>
#include <concepts>
#include <exception>
#include <string>

namespace ReaK::pp {

/**
 * This traits class defines the traits that characterize a spatial path within a
 * topology.
 * \tparam SpatialPath The spatial path type for which the traits are sought.
 */
template <typename Path>
struct spatial_path_traits {
  /** This type describes a point in the space or topology. */
  using point_type = typename Path::point_type;
  /** This type describes the difference between two points in the space or topology. */
  using point_difference_type = typename Path::point_difference_type;

  /** This type describes waypoints used by the path to quickly access local parameters of the path. */
  using waypoint_descriptor = typename Path::waypoint_descriptor;
  /** This type describes const-waypoints used by the path to quickly access local parameters of the path. */
  using const_waypoint_descriptor = typename Path::const_waypoint_descriptor;

  using const_waypoint_pair = std::pair<const_waypoint_descriptor, point_type>;

  /** This type is the topology type in which the path exists. */
  using topology = typename Path::topology;
  /** This type is the distance metric type used on the topology and defining the travel distances along the path. */
  using distance_metric = typename Path::distance_metric;
};

/**
 * This concept class defines the requirements for a type to model a spatial-path as used in
 * ReaK::pp. A spatial path is a continuous curve within a topology.
 *
 * Required concepts:
 *
 * The topology should model the MetricSpace.
 *
 * The distance-metric should model the DistanceMetric.
 *
 * Valid expressions:
 *
 * std::pair< const_waypoint_descriptor, point_type > w_p;  A waypoint-point pair is used to associate a point and a
 *waypoint.
 *
 * pt = p.move_away_from(pt, d);  A point along the path (p), at a distance (d) from a point (pt), can be obtained.
 *
 * d = p.travel_distance(pt,pt);  The travel distance, along the path (p), between two points (pt,pt), can be obtained.
 *
 * w_p = p.move_away_from(w_p, d);  A waypoint along the path (p), at a distance (d) from a waypoint (w_p), can be
 *obtained.
 *
 * d = p.travel_distance(w_p,w_p);  The travel distance, along the path (p), between two waypoints (w_p,w_p), can be
 *obtained.
 *
 * pt = p.get_start_point();  The starting point of the path can be obtained.
 *
 * w_p = p.get_start_waypoint();  The starting waypoint-pair of the path can be obtained.
 *
 * pt = p.get_end_point();  The end point of the path can be obtained.
 *
 * w_p = p.get_end_waypoint();  The end waypoint-pair of the path can be obtained.
 */
template <typename Path,
          typename Space = typename spatial_path_traits<Path>::topology>
concept SpatialPath = Topology<Space>&&
    DistanceMetric<typename spatial_path_traits<Path>::distance_metric, Space>&&
    requires(const Path& p, const topology_point_type_t<Space>& pt,
             const typename spatial_path_traits<Path>::const_waypoint_pair& wp,
             double d) {
  {
    p.move_away_from(pt, d)
    } -> std::convertible_to<topology_point_type_t<Space>>;
  { p.travel_distance(pt, pt) } -> std::convertible_to<double>;
  {
    p.move_away_from(wp, d)
    } -> std::convertible_to<
        typename spatial_path_traits<Path>::const_waypoint_pair>;
  { p.travel_distance(wp, wp) } -> std::convertible_to<double>;
  { p.get_start_point() } -> std::convertible_to<topology_point_type_t<Space>>;
  { p.get_end_point() } -> std::convertible_to<topology_point_type_t<Space>>;
  {
    p.get_start_waypoint()
    } -> std::convertible_to<
        typename spatial_path_traits<Path>::const_waypoint_pair>;
  {
    p.get_end_waypoint()
    } -> std::convertible_to<
        typename spatial_path_traits<Path>::const_waypoint_pair>;
};

/**
 * This exception class represents an exception occurring when a path is invalid.
 */
struct invalid_path : public std::exception {
  std::string message;
  explicit invalid_path(const std::string& aSender)
      : message(std::string("Invalid path reported! Originating from ") +
                aSender) {}
  ~invalid_path() noexcept override = default;
  const char* what() const noexcept override { return message.c_str(); }
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_INTERPOLATION_SPATIAL_PATH_CONCEPT_H_

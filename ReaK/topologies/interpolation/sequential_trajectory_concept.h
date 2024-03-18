/**
 * \file sequential_trajectory_concept.h
 *
 * This library defines the traits and concepts related to a sequential spatial trajectory. A
 * trajectory is simply a continuous curve in a temporal topology (or time-space) which can be travelled
 * sequentially via either increments in time or in fractions between waypoints.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2013
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_SEQUENTIAL_TRAJECTORY_CONCEPT_H_
#define REAK_TOPOLOGIES_INTERPOLATION_SEQUENTIAL_TRAJECTORY_CONCEPT_H_

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/temporal_space_concept.h"

#include <concepts>

namespace ReaK::pp {

/**
 * This traits class defines the traits that characterize a sequential spatial trajectory within a
 * temporal topology.
 * \tparam SequentialTraj The spatial trajectory type for which the traits are sought.
 */
template <typename SequentialTraj>
struct sequential_trajectory_traits {
  /** This type describes a point in the space or topology. */
  using point_type = typename SequentialTraj::point_type;

  /** This type describes an iterator, corresponding to a point on the trajectory, which can be incremented by time to
   * travel to the next iterator. */
  using point_time_iterator = typename SequentialTraj::point_time_iterator;
  /** This type describes an iterator, corresponding to a point on the trajectory, which can be incremented by a
   * fraction between waypoints to travel to the next iterator. */
  using point_fraction_iterator =
      typename SequentialTraj::point_fraction_iterator;

  /** This type is the topology type in which the path exists. */
  using topology = typename SequentialTraj::topology;
};

/**
 * This concept class defines the requirements for a type to model a sequential trajectory
 * as used in ReaK::pp. A sequential trajectory is a continuous curve within a temporal topology
 * which can be travelled sequentially via either increments in time or in fractions
 * between waypoints.
 *
 * Required concepts:
 *
 * The topology should model the TemporalSpace.
 *
 * Valid expressions:
 *
 * tit = traj.begin_time_travel();  The start of the time-iterator range of the sequential trajectory can be obtained.
 *
 * tit = traj.end_time_travel();  The end of the time-iterator range (one-past-last) of the sequential trajectory can be
 *obtained.
 *
 * pt = *tit;  A point can be obtained from dereferencing a time-iterator.
 *
 * tit = tit + d;
 * tit = d + tit;
 * tit += d;
 * tit = tit - d;
 * tit -= d;  A time-iterator can be incremented by a time (double).
 *
 * b = (tit != tit);
 * b = (tit == tit);  Two time-iterator can be compared for inequality.
 *
 * fit = traj.begin_fraction_travel();  The start of the fraction-iterator range of the sequential trajectory can be
 *obtained.
 *
 * fit = traj.end_fraction_travel();  The end of the fraction-iterator range (one-past-last) of the sequential
 *trajectory can be obtained.
 *
 * pt = *fit;  A point can be obtained from dereferencing a fraction-iterator.
 *
 * fit = fit + f;
 * fit = f + fit;
 * fit += f;
 * fit = fit - f;
 * fit -= f;  A fraction-iterator can be incremented by a fraction (double).
 *
 * b = (fit != fit);
 * b = (fit == fit);  Two fraction-iterator can be compared for equality.
 *
 * d = traj.travel_distance(pt,pt);  The travel distance (as of the distance-metric), along the trajectory (p), between
 *two points (pt,pt), can be obtained.
 */
template <typename Traj, typename Space = typename sequential_trajectory_traits<Traj>::topology>
concept SequentialTrajectory = TemporalSpace<Space> &&
  requires (const Traj& traj, const topology_point_type_t<Space>& pt) {
    { traj.travel_distance(pt, pt) } -> std::convertible_to<double>;
  } &&
  requires (const Traj& traj) {
    { traj.begin_time_travel() } -> std::convertible_to<typename sequential_trajectory_traits<Traj>::point_time_iterator>;
    { traj.end_time_travel() } -> std::convertible_to<typename sequential_trajectory_traits<Traj>::point_time_iterator>;
  } &&
  requires (double d,
            typename sequential_trajectory_traits<Traj>::point_time_iterator& tit) {
    { *tit } -> std::convertible_to<topology_point_type_t<Space>>;
    tit = tit + d;
    tit = d + tit;
    tit += d;
    tit = tit - d;
    tit -= d;
    { tit != tit } -> std::convertible_to<bool>;
    { tit == tit } -> std::convertible_to<bool>;
  } &&
  requires (const Path& path) {
    { path.begin_fraction_travel() } -> std::convertible_to<typename sequential_trajectory_traits<Traj>::point_fraction_iterator>;
    { path.end_fraction_travel() } -> std::convertible_to<typename sequential_trajectory_traits<Traj>::point_fraction_iterator>;
  } &&
  requires (double d,
            typename sequential_trajectory_traits<Traj>::point_fraction_iterator& fit) {
    { *fit } -> std::convertible_to<topology_point_type_t<Space>>;
    fit = fit + d;
    fit = d + fit;
    fit += d;
    fit = fit - d;
    fit -= d;
    { fit != fit } -> std::convertible_to<bool>;
    { fit == fit } -> std::convertible_to<bool>;
  };

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_INTERPOLATION_SEQUENTIAL_TRAJECTORY_CONCEPT_H_

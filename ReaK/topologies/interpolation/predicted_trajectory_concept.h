/**
 * \file predicted_trajectory_concept.h
 *
 * This library defines the concept that represents a predicted trajectory. This concept
 * supplements the SpatialTrajectory with a few additional requirements which
 * characterize a predicted trajectory (see PredictedTrajectory).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2011
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_PREDICTED_TRAJECTORY_CONCEPT_H_
#define REAK_TOPOLOGIES_INTERPOLATION_PREDICTED_TRAJECTORY_CONCEPT_H_


#include "ReaK/topologies/interpolation/spatial_path_concept.h"
#include "ReaK/topologies/interpolation/spatial_trajectory_concept.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/temporal_space_concept.h"

namespace ReaK::pp {

/**
 * This concept represents a predicted trajectory, as used in ReaK::pp. This concept
 * supplements the SpatialTrajectory with a few additional requirements which
 * characterize a predicted trajectory (see PredictedTrajectory).
 *
 * Required concepts:
 *
 * The Trajectory should model the SpatialTrajectory on the given Topology.
 *
 * The Space should model the TemporalSpace.
 *
 * Valid expressions:
 *
 * std::pair< const_waypoint_descriptor, point_type> w_p;  A waypoint is a pair of a const-waypoint descriptor and a
 *point on the topology.
 *
 * p.set_initial_point(pt);  The initial temporal point (pt) can be set to seed the predicted trajectory (p).
 *
 * p.set_initial_point(w_p);  The initial temporal waypoint (w_p) can be set to seed the predicted trajectory (p).
 */
template <typename Trajectory, typename Space>
concept PredictedTrajectory = SpatialTrajectory<Trajectory, Space>&& requires(
    Trajectory& traj, const topology_point_type_t<Space>& pt,
    const typename spatial_trajectory_traits<Trajectory>::const_waypoint_pair&
        wp) {
  traj.set_initial_point(pt);
  traj.set_initial_point(wp);
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_INTERPOLATION_PREDICTED_TRAJECTORY_CONCEPT_H_

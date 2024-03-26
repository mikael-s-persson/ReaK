
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

#ifndef REAK_PLANNING_QT_CRS_PLANNER_DATA_H_
#define REAK_PLANNING_QT_CRS_PLANNER_DATA_H_

#include <atomic>
#include "ReaK/control/systems/satellite_invar_models.h"
#include "ReaK/topologies/interpolation/trajectory_base.h"

class SoTimerSensor;

template <typename Topology>
struct coin_animation_data {
  using temporal_space_type =
      ReaK::pp::temporal_space<Topology, ReaK::pp::time_poisson_topology,
                               ReaK::pp::time_distance_only>;
  using trajectory_type = ReaK::pp::trajectory_base<temporal_space_type>;

  std::shared_ptr<trajectory_type> trajectory;
  SoTimerSensor* animation_timer;
  std::atomic<bool> enabled;
};

using manip_cspace_type = ReaK::pp::hyperbox_topology<ReaK::vect<double, 7>>;
using CRS_sol_anim_data = coin_animation_data<manip_cspace_type>;

using sat_state_space_type = ReaK::pp::se3_1st_order_topology_t<double>;
using CRS_target_anim_data = coin_animation_data<sat_state_space_type>;

#endif  // REAK_PLANNING_QT_CRS_PLANNER_DATA_H_

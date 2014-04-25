
/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_TARGET_PREDICT_DATA_HPP
#define REAK_TARGET_PREDICT_DATA_HPP

#include "CRS_planner_data.hpp"

#include "topologies/temporal_space.hpp"
#include "path_planning/trajectory_base.hpp"

class SoTimerSensor;


struct satellite_predict_data {
  
  typedef ReaK::pp::se3_1st_order_topology<double>::type state_space_type;
  typedef ReaK::pp::temporal_space<state_space_type, ReaK::pp::time_poisson_topology, ReaK::pp::time_distance_only> temp_state_space_type;
  
  typedef ReaK::pp::trajectory_base<temp_state_space_type> state_traj_type;
  
  ReaK::shared_ptr< state_traj_type > trajectory;
  SoTimerSensor* animation_timer;
  volatile bool enabled;
  
};




#endif


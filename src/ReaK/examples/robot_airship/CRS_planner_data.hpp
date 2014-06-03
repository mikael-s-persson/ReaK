
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

#ifndef REAK_CRS_PLANNER_DATA_HPP
#define REAK_CRS_PLANNER_DATA_HPP

#include <ReaK/ctrl/ss_systems/satellite_invar_models.hpp>
#include <ReaK/ctrl/interpolation/trajectory_base.hpp>

class SoTimerSensor;


template <typename Topology>
struct coin_animation_data {
  typedef ReaK::pp::temporal_space<Topology, ReaK::pp::time_poisson_topology, ReaK::pp::time_distance_only> temporal_space_type;
  typedef ReaK::pp::trajectory_base< temporal_space_type > trajectory_type;
  
  ReaK::shared_ptr< trajectory_type > trajectory;
  SoTimerSensor* animation_timer;
  volatile bool enabled;
};


typedef ReaK::pp::hyperbox_topology< ReaK::vect<double,7> > manip_cspace_type;
typedef coin_animation_data<manip_cspace_type> CRS_sol_anim_data;


typedef ReaK::pp::se3_1st_order_topology<double>::type sat_state_space_type;
typedef coin_animation_data<sat_state_space_type> CRS_target_anim_data;





#endif


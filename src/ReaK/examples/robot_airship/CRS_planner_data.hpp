
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

#include "satellite_invar_models.hpp"
#include "interpolation/discrete_point_trajectory.hpp"


class SoTimerSensor;
class SoQtExaminerViewer;
class SoSeparator;
class SoSwitch;
class SoCoordinate3;


namespace ReaK { 

namespace geom {
  class oi_scene_graph;
}; 

};


struct chaser_target_coin_nodes {
  // needed in constructor:
  SoQtExaminerViewer* eviewer;
  SoSeparator* sg_root;
  
  ReaK::geom::oi_scene_graph* sg_chaser_geom;
  SoSwitch* sw_chaser_geom;
  
  ReaK::geom::oi_scene_graph* sg_chaser_kin;
  SoSwitch* sw_chaser_kin;
  
  ReaK::geom::oi_scene_graph* sg_target_geom;
  SoSwitch* sw_target_geom;
  
  ReaK::geom::oi_scene_graph* sg_env_geom;
  SoSwitch* sw_env_geom;
  
  bool trace_motion_graph;
  SoSwitch* sw_motion_graph;
  
  bool trace_solutions;
  SoSwitch* sw_solutions;
};

typedef chaser_target_coin_nodes CRS_coin_nodes;


template <typename Topology>
struct coin_animation_data {
  typedef ReaK::pp::temporal_space<Topology, ReaK::pp::time_poisson_topology, ReaK::pp::time_distance_only> temporal_space_type;
  typedef ReaK::pp::discrete_point_trajectory< temporal_space_type > trajectory_type;
  
  ReaK::shared_ptr< trajectory_type > trajectory;
  SoTimerSensor* animation_timer;
  volatile bool enabled;
};


typedef ReaK::pp::hyperbox_topology< ReaK::vect<double,7> > manip_cspace_type;
typedef coin_animation_data<manip_cspace_type> CRS_sol_anim_data;


typedef ReaK::ctrl::satellite3D_inv_dt_system sat_sys_type;
typedef sat_sys_type::state_space_type sat_state_space_type;
typedef coin_animation_data<sat_state_space_type> CRS_target_anim_data;





#endif


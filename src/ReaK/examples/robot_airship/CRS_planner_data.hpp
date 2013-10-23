
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


#include "CRS_A465_geom_model.hpp"

#include "satellite_invar_models.hpp"
#include "interpolation/discrete_point_trajectory.hpp"
#include "path_planning/path_planner_options.hpp"

#include "serialization/archiver.hpp"

#include "base/chrono_incl.hpp"


class SoTimerSensor;
class SoQtExaminerViewer;
class SoSeparator;
class SoSwitch;
class SoCoordinate3;


namespace ReaK { 

namespace geom {
  class oi_scene_graph;
  class proxy_query_model_3D;
  class proxy_query_pair_3D;
  class oi_scene_graph;
}; 

namespace kte {
  class kte_map_chain;
  class manipulator_kinematics_model;
}; 

};


using namespace ReaK;



typedef ctrl::satellite3D_inv_dt_system sat_sys_type;
typedef sat_sys_type::state_space_type sat_state_space_type;
typedef pp::temporal_space<sat_state_space_type, pp::time_poisson_topology, pp::time_distance_only> sat_temp_space_type;
typedef pp::discrete_point_trajectory< sat_temp_space_type > sat_traj_type;
typedef pp::topology_traits< sat_temp_space_type >::point_type temp_point_type;



struct CRS_coin_nodes {
  // needed in constructor:
  SoQtExaminerViewer * eviewer;
  SoSeparator* sg_root;
  
  geom::oi_scene_graph* sg_robot_geom;
  SoSwitch* sw_robot_geom;
  geom::oi_scene_graph* sg_robot_kin;
  SoSwitch* sw_robot_kin;
  geom::oi_scene_graph* sg_airship_geom;
  SoSwitch* sw_airship_geom;
  geom::oi_scene_graph* sg_lab_geom;
  SoSwitch* sw_lab_geom;
  
  SoSwitch* sw_motion_graph;
  SoSwitch* sw_solutions;
  
  SoSwitch* sw_proxy_geom;
  SoCoordinate3* l_r_proxy_line;
  SoCoordinate3* r_a_proxy_line;
  SoCoordinate3* l_a_proxy_line;
};

struct CRS_model_data {
  robot_airship::CRS_A465_geom_builder builder;
  shared_ptr< kte::kte_map_chain > kin_chain;
  shared_ptr< kte::manipulator_kinematics_model > manip_kin_mdl;
  shared_ptr< pp::joint_limits_collection<double> > manip_jt_limits;
  
  shared_ptr< frame_3D<double> > airship_frame;
  pose_3D<double> target_frame;
  shared_ptr< kte::kte_map_chain > airship_chain;
  
  shared_ptr< geom::proxy_query_model_3D > robot_proxy;
  shared_ptr< geom::proxy_query_model_3D > lab_proxy;
  shared_ptr< geom::proxy_query_model_3D > airship_proxy;
  shared_ptr< geom::proxy_query_pair_3D > robot_lab_proxy;
  shared_ptr< geom::proxy_query_pair_3D > robot_airship_proxy;
  shared_ptr< geom::proxy_query_pair_3D > lab_airship_proxy;
};

struct CRS_sol_anim_data {
  std::vector< vect<double,7> > bestsol_trajectory;
  SoTimerSensor* animation_timer;
  std::chrono::high_resolution_clock::time_point animation_last_render;
  volatile bool enabled;
};

struct CRS_target_anim_data {
  shared_ptr< sat_traj_type > target_trajectory;
  SoTimerSensor* target_anim_timer;
  std::chrono::high_resolution_clock::time_point target_anim_last_render;
  volatile bool enabled;
};



struct CRS_planning_options {
  std::size_t space_order;
  std::size_t interp_id;
  double min_travel;
  double max_travel;
  std::size_t planning_algo;
  std::size_t max_vertices;
  std::size_t prog_interval;
  std::size_t max_results;
  std::size_t planning_options;
  std::size_t store_policy;
  std::size_t knn_method;
  double init_SA_temp;
  double init_relax;
  
  void save(serialization::oarchive& out) const {
    out & RK_SERIAL_SAVE_WITH_NAME(space_order)
        & RK_SERIAL_SAVE_WITH_NAME(interp_id)
        & RK_SERIAL_SAVE_WITH_NAME(min_travel)
        & RK_SERIAL_SAVE_WITH_NAME(max_travel)
        & RK_SERIAL_SAVE_WITH_NAME(planning_algo)
        & RK_SERIAL_SAVE_WITH_NAME(max_vertices)
        & RK_SERIAL_SAVE_WITH_NAME(prog_interval)
        & RK_SERIAL_SAVE_WITH_NAME(max_results)
        & RK_SERIAL_SAVE_WITH_NAME(planning_options)
        & RK_SERIAL_SAVE_WITH_NAME(store_policy)
        & RK_SERIAL_SAVE_WITH_NAME(knn_method)
        & RK_SERIAL_SAVE_WITH_NAME(init_SA_temp)
        & RK_SERIAL_SAVE_WITH_NAME(init_relax);
  };
  void load(serialization::iarchive& in) {
    in & RK_SERIAL_LOAD_WITH_NAME(space_order)
       & RK_SERIAL_LOAD_WITH_NAME(interp_id)
       & RK_SERIAL_LOAD_WITH_NAME(min_travel)
       & RK_SERIAL_LOAD_WITH_NAME(max_travel)
       & RK_SERIAL_LOAD_WITH_NAME(planning_algo)
       & RK_SERIAL_LOAD_WITH_NAME(max_vertices)
       & RK_SERIAL_LOAD_WITH_NAME(prog_interval)
       & RK_SERIAL_LOAD_WITH_NAME(max_results)
       & RK_SERIAL_LOAD_WITH_NAME(planning_options)
       & RK_SERIAL_LOAD_WITH_NAME(store_policy)
       & RK_SERIAL_LOAD_WITH_NAME(knn_method)
       & RK_SERIAL_LOAD_WITH_NAME(init_SA_temp)
       & RK_SERIAL_LOAD_WITH_NAME(init_relax);
  };
  
};



#endif


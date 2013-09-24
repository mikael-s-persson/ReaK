/**
 * \file CRS_planners_utility.hpp
 * 
 * This library defines a class
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2013
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

#ifndef REAK_CRS_PLANNERS_UTILITY_HPP
#define REAK_CRS_PLANNERS_UTILITY_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "shapes/colored_model.hpp"
#include "proximity/proxy_query_model.hpp"
#include "kte_models/manip_P3R3R_arm.hpp"

#include "CRS_workspaces.hpp"
#include "CRS_rrt_planners.hpp"
#include "CRS_rrtstar_planners.hpp"
#include "CRS_prm_planners.hpp"
#include "CRS_fadprm_planners.hpp"
#include "CRS_sbastar_planners.hpp"

namespace ReaK {
  
namespace robot_airship {


class scenario_data : public named_object {
  public:
    shared_ptr< frame_3D<double> >                    chaser_base_frame;
    shared_ptr< kte::manip_P3R3R_kinematics >         chaser_kin_model;
    shared_ptr< pp::joint_limits_collection<double> > chaser_jt_limits;
    shared_ptr< geom::proxy_query_model_3D >          chaser_proxy;
    shared_ptr< geom::colored_model_3D >              chaser_geom_model;
    vect_n<double>                                    chaser_jt_positions;
    
    shared_ptr< kte::kte_map_chain >     target_kin_chain;
    shared_ptr< frame_3D<double> >       target_state;
    shared_ptr< frame_3D<double> >       target_frame;
    shared_ptr< geom::colored_model_3D > target_geom_model;
    
    shared_ptr< geom::proxy_query_pair_3D > chaser_target_proxy;
    
    std::vector< shared_ptr< geom::colored_model_3D > > env_geom_models;
    
    std::vector< shared_ptr< geom::proxy_query_pair_3D > > chaser_env_proxies;
    std::vector< shared_ptr< geom::proxy_query_pair_3D > > target_env_proxies;
    
    scenario_data(const std::string& aChaserFileName);
    
    
    void load_target(const std::string& fileName);
    void load_positions(const std::string& fileName);
    
    void load_environment(const std::string& fileName);
    void clear_environment();
    
  private:
    
    vect_n<double> get_chaser_goal_config();
    
  public:
    
};





/*

typedef pp::manip_pp_traits< kte::manip_P3R3R_kinematics >::rl_o0_jt_space_type CRS3D_jspace_rl_o0_type;
typedef pp::manip_pp_traits< kte::manip_P3R3R_kinematics >::rl_o1_jt_space_type CRS3D_jspace_rl_o1_type;
typedef pp::manip_pp_traits< kte::manip_P3R3R_kinematics >::rl_o2_jt_space_type CRS3D_jspace_rl_o2_type;

typedef pp::manip_pp_traits< kte::manip_P3R3R_kinematics >::o0_jt_space_type CRS3D_jspace_o0_type;
typedef pp::manip_pp_traits< kte::manip_P3R3R_kinematics >::o1_jt_space_type CRS3D_jspace_o1_type;
typedef pp::manip_pp_traits< kte::manip_P3R3R_kinematics >::o2_jt_space_type CRS3D_jspace_o2_type;


typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::linear_interpolation_tag>::rl_o0_workspace_type CRS3D_workspace_o0_i1_type;
typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::linear_interpolation_tag>::rl_o1_workspace_type CRS3D_workspace_o1_i1_type;
typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::linear_interpolation_tag>::rl_o2_workspace_type CRS3D_workspace_o2_i1_type;

typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::cubic_hermite_interpolation_tag>::rl_o1_workspace_type CRS3D_workspace_o1_i3_type;
typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::cubic_hermite_interpolation_tag>::rl_o2_workspace_type CRS3D_workspace_o2_i3_type;

typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::quintic_hermite_interpolation_tag>::rl_o2_workspace_type CRS3D_workspace_o2_i5_type;

typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::svp_Ndof_interpolation_tag>::rl_o1_workspace_type CRS3D_workspace_o1_svp_type;
typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::svp_Ndof_interpolation_tag>::rl_o2_workspace_type CRS3D_workspace_o2_svp_type;

typedef pp::manip_static_workspace< kte::manip_P3R3R_kinematics, pp::sap_Ndof_interpolation_tag>::rl_o2_workspace_type CRS3D_workspace_o2_sap_type;


typedef pp::manip_DK_map< kte::manip_P3R3R_kinematics >::rl_o0_map_type CRS3D_rlDK_o0_type;
typedef pp::manip_DK_map< kte::manip_P3R3R_kinematics >::rl_o1_map_type CRS3D_rlDK_o1_type;
typedef pp::manip_DK_map< kte::manip_P3R3R_kinematics >::rl_o2_map_type CRS3D_rlDK_o2_type;

*/


/*
struct all_robot_info {
  SoQtExaminerViewer * eviewer;
  SoSeparator* sg_root;
  ReaK::geom::oi_scene_graph* sg_robot_geom;
  SoSwitch* sw_robot_geom;
  ReaK::geom::oi_scene_graph* sg_robot_kin;
  SoSwitch* sw_robot_kin;
  ReaK::geom::oi_scene_graph* sg_airship_geom;
  SoSwitch* sw_airship_geom;
  ReaK::geom::oi_scene_graph* sg_lab_geom;
  SoSwitch* sw_lab_geom;
  ReaK::geom::oi_scene_graph* sg_proxy_geom;
  SoSwitch* sw_proxy_geom;
  SoSwitch* sw_motion_graph;
  SoSwitch* sw_solutions;
  ReaK::robot_airship::CRS_A465_geom_builder builder;
  ReaK::shared_ptr< ReaK::kte::kte_map_chain > kin_chain;
  ReaK::shared_ptr< ReaK::kte::manipulator_kinematics_model > manip_kin_mdl;
  ReaK::shared_ptr< ReaK::pp::joint_limits_collection<double> > manip_jt_limits;
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > robot_proxy;
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > lab_proxy;
  ReaK::shared_ptr< ReaK::geom::proxy_query_pair_3D > robot_lab_proxy;
  ReaK::shared_ptr< ReaK::frame_3D<double> > airship_frame;
  ReaK::pose_3D<double> target_frame;
  ReaK::shared_ptr< ReaK::kte::kte_map_chain > airship_chain;
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > airship_proxy;
  ReaK::shared_ptr< ReaK::geom::proxy_query_pair_3D > robot_airship_proxy;
  ReaK::shared_ptr< ReaK::geom::proxy_query_pair_3D > lab_airship_proxy;
  SoCoordinate3* l_r_proxy_line;
  SoCoordinate3* r_a_proxy_line;
  SoCoordinate3* l_a_proxy_line;
  
  std::vector< ReaK::vect<double,7> > bestsol_trajectory;
  SoTimerSensor* animation_timer;
  std::size_t animation_progress;
  std::chrono::high_resolution_clock::time_point animation_last_render;
} r_info;
*/


};

};

#endif


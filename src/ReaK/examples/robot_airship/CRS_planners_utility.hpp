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

#include "CRS_A465_geom_model.hpp"
#include "CRS_A465_models.hpp"
#include "proximity/proxy_query_model.hpp"

#include "kte_models/manip_kinematics_model.hpp"

#include "CRS_workspaces.hpp"
#include "CRS_rrt_planners.hpp"
#include "CRS_rrtstar_planners.hpp"
#include "CRS_prm_planners.hpp"
#include "CRS_fadprm_planners.hpp"
#include "CRS_sbastar_planners.hpp"

#include "serialization/bin_archiver.hpp"
#include "serialization/xml_archiver.hpp"
#include "serialization/protobuf_archiver.hpp"

namespace ReaK {
  
namespace robot_airship {


class scenario_data : public named_object {
  public:
    CRS_A465_geom_builder chaser_builder;
    shared_ptr< kte::manipulator_kinematics_model > chaser_kin_model;
    shared_ptr< pp::joint_limits_collection<double> > chaser_jt_limits;
    
    shared_ptr< kte::kte_map_chain > target_kin_chain;
    shared_ptr< frame_3D<double> > target_state;
    shared_ptr< frame_3D<double> > target_frame;
    
    shared_ptr< geom::proxy_query_pair_3D > chaser_env_proxy;
    shared_ptr< geom::proxy_query_pair_3D > chaser_target_proxy;
    shared_ptr< geom::proxy_query_pair_3D > target_env_proxy;
    
    
    scenario_data() : named_object() {
      setName("robot_airship_scenario_data");
      
      
      chaser_builder.create_geom_from_preset();
      
      shared_ptr< geom::proxy_query_model_3D > chaser_proxy = chaser_builder.get_proximity_model();
      chaser_kin_model = chaser_builder.get_manipulator_kin_model();
      chaser_jt_limits = shared_ptr< pp::joint_limits_collection<double> >(&(chaser_builder.joint_rate_limits), ReaK::null_deleter());
      
      shared_ptr< geom::proxy_query_model_3D > env_proxy;
      shared_ptr< geom::colored_model_3D > env_geom_model;
      {
        serialization::xml_iarchive in("models/MD148_lab_model.xml");
        in >> env_geom_model >> env_proxy;
      };
      
      shared_ptr< geom::proxy_query_model_3D > target_proxy;
      shared_ptr< geom::colored_model_3D > target_geom_model;
      {
        shared_ptr< kte::position_measure_3D > target_position;
        shared_ptr< kte::rotation_measure_3D > target_rotation;
        shared_ptr< kte::driving_actuator_3D > target_actuator;
        shared_ptr< kte::inertia_3D > target_inertia;
        shared_ptr< kte::mass_matrix_calc > target_mass_calc;
        
        serialization::xml_iarchive in("models/airship3D_with_geom.xml");
        in >> target_state
           >> target_position
           >> target_rotation
           >> target_actuator
           >> target_inertia
           >> target_kin_chain
           >> target_mass_calc
           >> target_frame
           >> target_geom_model 
           >> target_proxy;
      };
      
      chaser_env_proxy     = shared_ptr< geom::proxy_query_pair_3D >(new geom::proxy_query_pair_3D("chaser_env_proxy", chaser_proxy, env_proxy));
      chaser_target_proxy  = shared_ptr< geom::proxy_query_pair_3D >(new geom::proxy_query_pair_3D("chaser_target_proxy", chaser_proxy, target_proxy));
      target_env_proxy     = shared_ptr< geom::proxy_query_pair_3D >(new geom::proxy_query_pair_3D("target_env_proxy", env_proxy, target_proxy));
      
      chaser_env_proxy->setModelPair(chaser_proxy, env_proxy);
      chaser_target_proxy->setModelPair(chaser_proxy, target_proxy);
      target_env_proxy->setModelPair(env_proxy, target_proxy);
      
      
    };
    
    
  
  
};

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



};

};

#endif


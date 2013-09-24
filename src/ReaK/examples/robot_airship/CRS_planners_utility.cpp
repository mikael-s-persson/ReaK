
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

#include "CRS_planners_utility.hpp"

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

#include "serialization/archiver_factory.hpp"

namespace ReaK {
  
namespace robot_airship {

  
scenario_data::scenario_data(const std::string& aChaserFileName) : named_object {
  setName("robot_airship_scenario_data");
  
  {
    shared_ptr< serialization::iarchive > p_in = serialization::open_iarchive(aChaserFileName);
    (*p_in) >> chaser_base_frame >> chaser_kin_model >> chaser_jt_limits >> chaser_geom_model >> chaser_proxy;
  };
  
};



void scenario_data::load_target(const std::string& fileName) {
  
  shared_ptr< geom::proxy_query_model_3D > target_proxy;
  {
    shared_ptr< kte::position_measure_3D > target_position;
    shared_ptr< kte::rotation_measure_3D > target_rotation;
    shared_ptr< kte::driving_actuator_3D > target_actuator;
    shared_ptr< kte::inertia_3D > target_inertia;
    shared_ptr< kte::mass_matrix_calc > target_mass_calc;
    
    shared_ptr< serialization::iarchive > p_in = serialization::open_iarchive(fileName);
    (*p_in) >> target_state
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
  
  chaser_target_proxy = shared_ptr< geom::proxy_query_pair_3D >(new geom::proxy_query_pair_3D("chaser_target_proxy", chaser_proxy, target_proxy));
  
};


void scenario_data::load_environment(const std::string& fileName) {
  
  shared_ptr< geom::proxy_query_model_3D > env_proxy;
  shared_ptr< geom::colored_model_3D > env_geom_model;
  {
    shared_ptr< serialization::iarchive > p_in = serialization::open_iarchive(fileName);
    (*p_in) >> env_geom_model >> env_proxy;
  };
  
  env_geom_models.push_back(env_geom_model);
  
  chaser_env_proxies.push_back(
    shared_ptr< geom::proxy_query_pair_3D >(new geom::proxy_query_pair_3D("chaser_env_proxy:" + fileName, chaser_proxy, env_proxy))
  );
  
  target_env_proxies.push_back(
    shared_ptr< geom::proxy_query_pair_3D >(new geom::proxy_query_pair_3D("target_env_proxy:" + fileName, env_proxy, target_proxy))
  );
  
};


void scenario_data::clear_environment() {
  env_geom_models.clear();
  chaser_env_proxies.clear();
  target_env_proxies.clear();
};


void scenario_data::load_positions(const std::string& fileName) {
  
  vect<double,7> chaser_joint_positions;
  vect<double,3> target_position;
  quaternion<double> target_quaternion;
  
  shared_ptr< serialization::iarchive > p_in = serialization::open_iarchive(fileName);
  (*p_in) & RK_SERIAL_LOAD_WITH_NAME(chaser_joint_positions)
          & RK_SERIAL_LOAD_WITH_NAME(target_position)
          & RK_SERIAL_LOAD_WITH_NAME(target_quaternion);
  
  chaser_jt_positions = chaser_joint_positions;
  chaser_kin_model->setJointPositions(chaser_jt_positions);
  chaser_kin_model->doDirectMotion();
  
  target_state->Position = target_position;
  target_state->Quat = target_quaternion;
  target_kin_chain->doMotion();
  
};


namespace {
  
  // without C++11 lambdas, this is the best approximation of a closure for restoring joint positions after an IK attempt.
  struct kin_model_jt_pos_restore {
    kte::direct_kinematics_model& kin_model;
    const vect_n<double>& jt_pos;
    
    kin_model_jt_pos_restore(kte::direct_kinematics_model& aKinModel, const vect_n<double>& aJtPos) :
                             kin_model(aKinModel), jt_pos(aJtPos) { };
    
    ~kin_model_jt_pos_restore() {
      kin_model.setJointPositions(jt_pos);
      kin_model.doDirectMotion();
    };
    
  };
  
};



vect_n<double> scenario_data::get_chaser_goal_config() {
  
  // a simple scope-guard object to restore the chaser configuration upon return (or unwinding):
  kin_model_jt_pos_restore restore_chaser_at_exit( *chaser_kin_model, chaser_jt_positions );
  
  target_kin_chain->doMotion();
  
  *(chaser_kin_model->getDependentFrame3D(0)->mFrame) = *target_frame;
  chaser_kin_model->doInverseMotion();
  
  for( std::vector< shared_ptr< geom::proxy_query_pair_3D > >::const_iterator it = chaser_env_proxies.begin(); it != chaser_env_proxies.end(); ++it) {
    shared_ptr< geom::proximity_finder_3D > tmp = (*it)->findMinimumDistance();
    if((tmp) && (tmp->getLastResult().mDistance < 0.0))
      throw optim::infeasible_problem("The goal configuration for the chaser leads to a collision with the environment!");
  };
  
  {
    shared_ptr< geom::proximity_finder_3D > tmp = chaser_target_proxy->findMinimumDistance();
    if((tmp) && (tmp->getLastResult().mDistance < 0.0))
      throw optim::infeasible_problem("The goal configuration for the chaser leads to a collision with the target geometry!");
  };
  
  return chaser_kin_model->getJointPositions();
};



};

};

#endif


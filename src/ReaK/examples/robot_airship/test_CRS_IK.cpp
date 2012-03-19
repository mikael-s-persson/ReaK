/**
 * \file test_CRS_IK.cpp
 *
 * This application tests the inverse kinematics for the CRS A465 manipulator models.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2012
 */


#include "CRS_A465_models.hpp"
#include "CRS_A465_2D_analog.hpp"

#include "path_planning/manipulator_topo_maps.hpp"

int main() {
  
  ReaK::robot_airship::CRS_A465_model_builder builder;
  
  builder.load_kte_from_file("models/CRS_A465_raw_components.xml");
  builder.load_limits_from_file("models/CRS_A465_limits.xml");
  
  typedef ReaK::robot_airship::CRS_A465_model_builder::joint_space_type JointSpaceType;
  typedef ReaK::robot_airship::CRS_A465_model_builder::end_effector_space_type EESpaceType;
  
  JointSpaceType j_space = builder.get_joint_space();
  EESpaceType ee_space = builder.get_end_effector_space();
  
  ReaK::shared_ptr< ReaK::kte::manipulator_kinematics_model > model = builder.get_manipulator_kin_model();
  
  ReaK::pp::manip_inverse_kin_map ik_map(model,builder.preferred_posture);
  
  typedef ReaK::pp::topology_traits< EESpaceType >::point_type EEPointType;
  typedef ReaK::pp::topology_traits< JointSpaceType >::point_type JointPointType;
  
  EEPointType ee_x;
  set_frame_3D(ee_x, ReaK::frame_3D<double>(ReaK::weak_ptr< ReaK::pose_3D<double> >(),
					    ReaK::vect<double,3>(2.0, 0.3, 0.7),
					    ReaK::quaternion<double>(),
					    ReaK::vect<double,3>(0.0, 0.0, 0.0),
					    ReaK::vect<double,3>(0.0, 0.0, 0.0),
					    ReaK::vect<double,3>(0.0, 0.0, 0.0),
					    ReaK::vect<double,3>(0.0, 0.0, 0.0),
					    ReaK::vect<double,3>(0.0, 0.0, 0.0),
					    ReaK::vect<double,3>(0.0, 0.0, 0.0)));
  
  JointPointType j_x = ik_map.map_to_space(ee_x, ee_space, j_space);
  
  
  ReaK::robot_airship::CRS_A465_2D_model_builder builder2D;
  
  builder2D.load_kte_from_file("models/CRS_A465_2D_raw_components.xml");
  builder2D.load_limits_from_file("models/CRS_A465_2D_limits.xml");
  
  
  return 0;
};








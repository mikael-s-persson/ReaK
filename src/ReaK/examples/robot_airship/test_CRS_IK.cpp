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
  ReaK::pp::manip_direct_kin_map dk_map(model);
  
  typedef ReaK::pp::topology_traits< EESpaceType >::point_type EEPointType;
  typedef ReaK::pp::topology_traits< JointSpaceType >::point_type JointPointType;
  
  EEPointType ee_x;
  set_frame_3D(ee_x, ReaK::frame_3D<double>(ReaK::weak_ptr< ReaK::pose_3D<double> >(),
					    ReaK::vect<double,3>(2.0, 0.5, 0.5),
					    ReaK::axis_angle<double>(M_PI * 0.5,ReaK::vect<double,3>(0.0, 0.0, 1.0)).getQuaternion(),
					    ReaK::vect<double,3>(0.1, 0.1, 0.0),
					    ReaK::vect<double,3>(0.0, 0.0, 0.1),
					    ReaK::vect<double,3>(0.0, 0.0, 0.0),
					    ReaK::vect<double,3>(0.0, 0.0, 0.0),
					    ReaK::vect<double,3>(0.0, 0.0, 0.0),
					    ReaK::vect<double,3>(0.0, 0.0, 0.0)));
  
  JointPointType j_x;
  
  try {
    j_x = ik_map.map_to_space(ee_x, ee_space, j_space);
  } catch(std::exception& e) {
    std::cout << "ERROR: An exception occurred during 3D IK problem: " << e.what() << std::endl;
  };
  
  ReaK::robot_airship::CRS_A465_2D_model_builder builder2D;
  
  builder2D.load_kte_from_file("models/CRS_A465_2D_raw_components.xml");
  builder2D.load_limits_from_file("models/CRS_A465_2D_limits.xml");
  
  
  typedef ReaK::robot_airship::CRS_A465_2D_model_builder::joint_space_type JointSpaceType2D;
  typedef ReaK::robot_airship::CRS_A465_2D_model_builder::end_effector_space_type EESpaceType2D;
  
  JointSpaceType2D j_space2D = builder2D.get_joint_space();
  EESpaceType2D ee_space2D = builder2D.get_end_effector_space();
  
  ReaK::shared_ptr< ReaK::kte::manipulator_kinematics_model > model2D = builder2D.get_manipulator_kin_model();
  
  ReaK::pp::manip_inverse_kin_map ik_map_2D(model2D,builder2D.preferred_posture);
  ReaK::pp::manip_direct_kin_map dk_map_2D(model2D);
  
  typedef ReaK::pp::topology_traits< EESpaceType2D >::point_type EEPointType2D;
  typedef ReaK::pp::topology_traits< JointSpaceType2D >::point_type JointPointType2D;
  
  EEPointType2D ee_x_2D;
  set_frame_2D(ee_x_2D, ReaK::frame_2D<double>(ReaK::weak_ptr< ReaK::pose_2D<double> >(),
					    ReaK::vect<double,2>(2.0, 0.5),
					    ReaK::rot_mat_2D<double>(M_PI * 0.5),
					    ReaK::vect<double,2>(0.1, 0.1),
					    0.1,
					    ReaK::vect<double,2>(0.0, 0.0),
					    0.0,
					    ReaK::vect<double,2>(0.0, 0.0),
					    0.0));
  
  JointPointType2D j_x_2D;
  
  try {
    j_x_2D = ik_map_2D.map_to_space(ee_x_2D, ee_space2D, j_space2D);
  } catch(std::exception& e) {
    std::cout << "ERROR: An exception occurred during 2D IK problem: " << e.what() << std::endl;
  };
  
  
  
  return 0;
};








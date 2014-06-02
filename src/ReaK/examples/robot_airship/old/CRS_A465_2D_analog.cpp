
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



#include "CRS_A465_2D_analog.hpp"

#include "serialization/xml_archiver.hpp"

#include "rtti/typed_primitives.hpp"


namespace ReaK {


namespace robot_airship {


void CRS_A465_2D_model_builder::load_kte_from_file(const std::string& aFileName) {
  
  serialization::xml_iarchive complete_model_input(aFileName);
  serialization::iarchive& input_ref = complete_model_input;
  
  input_ref
   // load the base frame (start of the track.
             & RK_SERIAL_LOAD_WITH_NAME(robot_base)
   // load all the joint coordinates.
             & RK_SERIAL_LOAD_WITH_NAME(track_joint_coord)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_1_coord)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_2_coord)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_3_coord)
   // load the end-effector frame.
             & RK_SERIAL_LOAD_WITH_NAME(arm_EE)
   // load all the intermediate coordinate frames (between links and joints).
             & RK_SERIAL_LOAD_WITH_NAME(track_joint_end)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_1_base)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_1_end)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_2_base)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_2_end)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_3_base)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_3_end)
   // load all the joint jacobian relationships.
             & RK_SERIAL_LOAD_WITH_NAME(track_joint_jacobian)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_1_jacobian)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_2_jacobian)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_3_jacobian)
   // load all the joints.
             & RK_SERIAL_LOAD_WITH_NAME(track_joint)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_1)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_2)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_3)
   // load all the joint-dependency of the link frames.
             & RK_SERIAL_LOAD_WITH_NAME(link_0_dep_frame)
             & RK_SERIAL_LOAD_WITH_NAME(link_1_dep_frame)
             & RK_SERIAL_LOAD_WITH_NAME(link_2_dep_frame)
             & RK_SERIAL_LOAD_WITH_NAME(link_3_dep_frame)
   // load all the links (massless).
             & RK_SERIAL_LOAD_WITH_NAME(link_0)
             & RK_SERIAL_LOAD_WITH_NAME(link_1)
             & RK_SERIAL_LOAD_WITH_NAME(link_2)
             & RK_SERIAL_LOAD_WITH_NAME(link_3)
   // load all the joint inertias (motor inertia).
             & RK_SERIAL_LOAD_WITH_NAME(track_joint_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_1_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_2_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_3_inertia)
   // load all the joint actuators (which apply a driving force).
             & RK_SERIAL_LOAD_WITH_NAME(track_actuator)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_1_actuator)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_2_actuator)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_3_actuator)
   // load all the link inertias (mass information of the links).
             & RK_SERIAL_LOAD_WITH_NAME(link_0_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(link_1_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(link_2_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(link_3_inertia);
             
};


void CRS_A465_2D_model_builder::load_limits_from_file(const std::string& aFileName) {
  serialization::xml_iarchive complete_model_input(aFileName);
  serialization::iarchive& input_ref = complete_model_input;
  
  input_ref & RK_SERIAL_LOAD_WITH_NAME(joint_lower_bounds)
            & RK_SERIAL_LOAD_WITH_NAME(joint_upper_bounds)
            & RK_SERIAL_LOAD_WITH_NAME(joint_rate_limits)
            & RK_SERIAL_LOAD_WITH_NAME(preferred_posture);
};
    
void CRS_A465_2D_model_builder::create_from_preset() {
  
  
  //declare all the intermediate frames.
  robot_base       = shared_ptr< frame_2D<double> >(new frame_2D<double>(), scoped_deleter());
  track_joint_end  = shared_ptr< frame_2D<double> >(new frame_2D<double>(), scoped_deleter());
  arm_joint_1_base = shared_ptr< frame_2D<double> >(new frame_2D<double>(), scoped_deleter());
  arm_joint_1_end  = shared_ptr< frame_2D<double> >(new frame_2D<double>(), scoped_deleter());
  arm_joint_2_base = shared_ptr< frame_2D<double> >(new frame_2D<double>(), scoped_deleter());
  arm_joint_2_end  = shared_ptr< frame_2D<double> >(new frame_2D<double>(), scoped_deleter());
  arm_joint_3_base = shared_ptr< frame_2D<double> >(new frame_2D<double>(), scoped_deleter());
  arm_joint_3_end  = shared_ptr< frame_2D<double> >(new frame_2D<double>(), scoped_deleter());
  arm_EE           = shared_ptr< frame_2D<double> >(new frame_2D<double>(), scoped_deleter());

  //declare all the joint coordinates.
  track_joint_coord = shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter());
  arm_joint_1_coord = shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter());
  arm_joint_2_coord = shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter());
  arm_joint_3_coord = shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter());
  
  //declare all the joint jacobians.
  track_joint_jacobian = shared_ptr< jacobian_gen_2D<double> >(new jacobian_gen_2D<double>(), scoped_deleter());
  arm_joint_1_jacobian = shared_ptr< jacobian_gen_2D<double> >(new jacobian_gen_2D<double>(), scoped_deleter());
  arm_joint_2_jacobian = shared_ptr< jacobian_gen_2D<double> >(new jacobian_gen_2D<double>(), scoped_deleter());
  arm_joint_3_jacobian = shared_ptr< jacobian_gen_2D<double> >(new jacobian_gen_2D<double>(), scoped_deleter());
  
  //set the absolute position of the base and add gravity (z-axis pointing up!) (x-axis pointing forward).
  //normally this frame is set via the feedback from the leg / lower-body motion control of wopa
  //robot_base->Acceleration = vect<double,2>(0.0,0.0,9.81); //put gravity acceleration on base of torso...
  robot_base->Position = vect<double,2>(0.0,0.0); //put the base of the torso somewhere 0.8 m off the ground level.
  
  //create revolute joint
  track_joint = shared_ptr< kte::prismatic_joint_2D >(new kte::prismatic_joint_2D("track_joint",
                                                                              track_joint_coord,
                                                                              vect<double,2>(1.0,0.0),
                                                                              robot_base,
                                                                              track_joint_end,
                                                                              track_joint_jacobian),
                                                       scoped_deleter());
                                              
  //create motor inertia
  shared_ptr< kte::joint_dependent_gen_coord > track_joint_dep_coord(new kte::joint_dependent_gen_coord(track_joint_coord), 
                                                                          scoped_deleter() );
  track_joint_dep_coord->add_joint(track_joint_coord, 
                                   shared_ptr< jacobian_gen_gen<double> >(new jacobian_gen_gen<double>(1.0,0.0), 
                                                                                    scoped_deleter()));
  track_joint_inertia = shared_ptr< kte::inertia_gen >(new kte::inertia_gen("track_joint_inertia",
                                                                        track_joint_dep_coord,
                                                                        1.0), 
                                                        scoped_deleter());
  
  //create force actuator
  track_actuator = shared_ptr< kte::driving_actuator_gen >(new kte::driving_actuator_gen("track_actuator",
                                                                                     track_joint_coord,
                                                                                     track_joint),
                                                            scoped_deleter());
  
  //create link from G to F10.0
  link_0 = shared_ptr< kte::rigid_link_2D >(new kte::rigid_link_2D("link_0",
                                                               track_joint_end,
                                                               arm_joint_1_base,
                                                               pose_2D<double>(weak_ptr<pose_2D<double> >(),
                                                                               vect<double,2>(0.0,0.0),
                                                                               rot_mat_2D<double>())),
                                             scoped_deleter());
  
  //create arm-base inertia of 
  link_0_dep_frame = shared_ptr< kte::joint_dependent_frame_2D >(new kte::joint_dependent_frame_2D(arm_joint_1_base),
                                                                    scoped_deleter());
  link_0_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_0_inertia = shared_ptr< kte::inertia_2D >(new kte::inertia_2D("link_0_inertia",
                                                                 link_0_dep_frame,
                                                                 1.0,
                                                                 1.0),
                                                 scoped_deleter());
  
  //create revolute joint
  arm_joint_1 = shared_ptr< kte::revolute_joint_2D >(new kte::revolute_joint_2D("arm_joint_1",
                                                                            arm_joint_1_coord,
                                                                            arm_joint_1_base,
                                                                            arm_joint_1_end,
                                                                            arm_joint_1_jacobian),
                                                      scoped_deleter());
                                              
  //create motor inertia
  shared_ptr< kte::joint_dependent_gen_coord > arm_joint_1_dep_coord(new kte::joint_dependent_gen_coord(arm_joint_1_coord),
                                                                          scoped_deleter());
  arm_joint_1_dep_coord->add_joint(arm_joint_1_coord,
                                   shared_ptr<jacobian_gen_gen<double> >(new jacobian_gen_gen<double>(1.0,0.0), 
                                                                                   scoped_deleter()));
  arm_joint_1_inertia = shared_ptr< kte::inertia_gen >(new kte::inertia_gen("arm_joint_1_inertia",
                                                                        arm_joint_1_dep_coord,
                                                                        1.0), //~71 kg m^2
                                                        scoped_deleter());
  
  //create force actuator
  arm_joint_1_actuator = shared_ptr< kte::driving_actuator_gen >(new kte::driving_actuator_gen("arm_joint_1_actuator",
                                                                                           arm_joint_1_coord,
                                                                                           arm_joint_1),
                                                                  scoped_deleter());
  
  //create link from F to CM (note that this is very approximate!!!)
  link_1 = shared_ptr< kte::rigid_link_2D >(new kte::rigid_link_2D("link_1",
                                                                arm_joint_1_end,
                                                                arm_joint_2_base,
                                                                pose_2D<double>(weak_ptr<pose_2D<double> >(),
                                                                                vect<double,2>(0.3048,0.0),
                                                                                rot_mat_2D<double>())),
                                             scoped_deleter());
  
  //create link1 inertia 
  link_1_dep_frame = shared_ptr< kte::joint_dependent_frame_2D >(new kte::joint_dependent_frame_2D(arm_joint_2_base),
                                                                    scoped_deleter());
  link_1_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_1_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  link_1_inertia = shared_ptr< kte::inertia_2D >(new kte::inertia_2D("link_1_inertia",
                                                                 link_1_dep_frame,
                                                                 1.0,
                                                                 1.0),
                                                  scoped_deleter());

  //create revolute joint
  arm_joint_2 = shared_ptr< kte::revolute_joint_2D >(new kte::revolute_joint_2D("arm_joint_2",
                                                                             arm_joint_2_coord,
                                                                             arm_joint_2_base,
                                                                             arm_joint_2_end,
                                                                             arm_joint_2_jacobian),
                                                      scoped_deleter());
    
  //create motor inertia
  shared_ptr< kte::joint_dependent_gen_coord > arm_joint_2_dep_coord(new kte::joint_dependent_gen_coord(arm_joint_2_coord),
                                                                          scoped_deleter());
  arm_joint_2_dep_coord->add_joint(arm_joint_2_coord,
                                   shared_ptr<jacobian_gen_gen<double> >(new jacobian_gen_gen<double>(1.0,0.0), 
                                                                                   scoped_deleter()));
  arm_joint_2_inertia = shared_ptr< kte::inertia_gen >(new kte::inertia_gen("arm_joint_2_inertia",
                                                                        arm_joint_2_dep_coord,
                                                                        1.0),
                                                        scoped_deleter());
  
  //create force actuator
  arm_joint_2_actuator = shared_ptr< kte::driving_actuator_gen >(new kte::driving_actuator_gen("arm_joint_2_actuator",
                                                                                           arm_joint_2_coord,
                                                                                           arm_joint_2),
                                                                  scoped_deleter());
  
  //create link 
  link_2 = shared_ptr< kte::rigid_link_2D >(new kte::rigid_link_2D("link_2",
                                                                arm_joint_2_end,
                                                                arm_joint_3_base,
                                                                pose_2D<double>(weak_ptr<pose_2D<double> >(),
                                                                                vect<double,2>(0.3302,0.0),
                                                                                rot_mat_2D<double>())),
                                              scoped_deleter());
  
  //create inertia
  link_2_dep_frame = shared_ptr< kte::joint_dependent_frame_2D >(new kte::joint_dependent_frame_2D(arm_joint_3_base),
                                                                    scoped_deleter());
  link_2_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_2_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  link_2_dep_frame->add_joint(arm_joint_2_coord,arm_joint_2_jacobian);
  link_2_inertia = shared_ptr< kte::inertia_2D >(new kte::inertia_2D("link_2_inertia",
                                                                 link_2_dep_frame,
                                                                 1.0,
                                                                 1.0),
                                                  scoped_deleter());

  //create revolute joint
  arm_joint_3 = shared_ptr< kte::revolute_joint_2D >(new kte::revolute_joint_2D("arm_joint_3",
                                                                            arm_joint_3_coord,
                                                                            arm_joint_3_base,
                                                                            arm_joint_3_end,
                                                                            arm_joint_3_jacobian),
                                                      scoped_deleter());
  
  //create motor inertia
  shared_ptr< kte::joint_dependent_gen_coord > arm_joint_3_dep_coord(new kte::joint_dependent_gen_coord(arm_joint_3_coord),
                                                                          scoped_deleter());
  arm_joint_3_dep_coord->add_joint(arm_joint_3_coord,
                                   shared_ptr<jacobian_gen_gen<double> > (new jacobian_gen_gen<double>(1.0,0.0), 
                                                                                   scoped_deleter()));
  arm_joint_3_inertia = shared_ptr< kte::inertia_gen >(new kte::inertia_gen("arm_joint_3_inertia",
                                                                         arm_joint_3_dep_coord,
                                                                         1.0), 
                                                        scoped_deleter());
  
  //create force actuator
  arm_joint_3_actuator = shared_ptr< kte::driving_actuator_gen >(new kte::driving_actuator_gen("arm_joint_3_actuator",
                                                                                           arm_joint_3_coord,
                                                                                           arm_joint_3),
                                                                  scoped_deleter());
  
  //create link 
  link_3 = shared_ptr< kte::rigid_link_2D >(new kte::rigid_link_2D("link_3",
                                                               arm_joint_3_end,
                                                               arm_EE,
                                                               pose_2D<double>(weak_ptr<pose_2D<double> >(),
                                                                               vect<double,2>(0.0762,0.0),
                                                                               rot_mat_2D<double>())),
                                              scoped_deleter());
  
  //create inertia
  link_3_dep_frame = shared_ptr< kte::joint_dependent_frame_2D >(new kte::joint_dependent_frame_2D(arm_EE),
                                                                    scoped_deleter());
  link_3_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_3_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  link_3_dep_frame->add_joint(arm_joint_2_coord,arm_joint_2_jacobian);
  link_3_dep_frame->add_joint(arm_joint_3_coord,arm_joint_3_jacobian);
  link_3_inertia = shared_ptr< kte::inertia_2D >(new kte::inertia_2D("link_3_inertia",
                                                                 link_3_dep_frame,
                                                                 1.0,
                                                                 1.0),
                                                  scoped_deleter());
  
  
  joint_lower_bounds.resize(4);
  joint_lower_bounds[0] = 0.0;
  joint_lower_bounds[1] = -3.05432619099;
  joint_lower_bounds[2] = -1.91986217719;
  joint_lower_bounds[3] = -1.83259571459;
  
  joint_upper_bounds.resize(4);
  joint_upper_bounds[0] = 3.0;
  joint_upper_bounds[1] = 3.05432619099;
  joint_upper_bounds[2] = 1.91986217719;
  joint_upper_bounds[3] = 1.83259571459;

  joint_rate_limits.gen_speed_limits.resize(4);
  joint_rate_limits.gen_speed_limits[0] = 0.8;
  joint_rate_limits.gen_speed_limits[1] = 3.14159265359;
  joint_rate_limits.gen_speed_limits[2] = 3.14159265359;
  joint_rate_limits.gen_speed_limits[3] = 3.01941960595;
  
  joint_rate_limits.gen_accel_limits.resize(4);
  joint_rate_limits.gen_accel_limits[0] = 3.0;
  joint_rate_limits.gen_accel_limits[1] = 12.5663706144;
  joint_rate_limits.gen_accel_limits[2] = 12.5663706144;
  joint_rate_limits.gen_accel_limits[3] = 24.9582083035;
  
  joint_rate_limits.gen_jerk_limits.resize(4);
  joint_rate_limits.gen_jerk_limits[0] = 12.0;
  joint_rate_limits.gen_jerk_limits[1] = 125.663706144;
  joint_rate_limits.gen_jerk_limits[2] = 125.663706144;
  joint_rate_limits.gen_jerk_limits[3] = 249.582083035;
  
  joint_rate_limits.frame2D_speed_limits.resize(0);
  joint_rate_limits.frame2D_accel_limits.resize(0);
  joint_rate_limits.frame2D_jerk_limits.resize(0);
  joint_rate_limits.frame3D_speed_limits.resize(0);
  joint_rate_limits.frame3D_accel_limits.resize(0);
  joint_rate_limits.frame3D_jerk_limits.resize(0);
  
  preferred_posture.resize(4);
  preferred_posture[0] = 1.5;
  preferred_posture[1] = 0.0;
  preferred_posture[2] = 0.0;
  preferred_posture[3] = 0.0;
  
};
    
void CRS_A465_2D_model_builder::save_kte_to_file(const std::string& aFileName) const {
  serialization::xml_oarchive complete_model_output(aFileName);
  serialization::oarchive& output_ref = complete_model_output;
  
  output_ref
   // save the base frame (start of the track.
             & RK_SERIAL_SAVE_WITH_NAME(robot_base)
   // save all the joint coordinates.
             & RK_SERIAL_SAVE_WITH_NAME(track_joint_coord)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1_coord)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2_coord)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3_coord)
   // save the end-effector frame.
             & RK_SERIAL_SAVE_WITH_NAME(arm_EE)
   // save all the intermediate coordinate frames (between links and joints).
             & RK_SERIAL_SAVE_WITH_NAME(track_joint_end)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1_base)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1_end)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2_base)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2_end)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3_base)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3_end)
   // save all the joint jacobian relationships.
             & RK_SERIAL_SAVE_WITH_NAME(track_joint_jacobian)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1_jacobian)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2_jacobian)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3_jacobian)
   // save all the joints.
             & RK_SERIAL_SAVE_WITH_NAME(track_joint)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3)
   // save all the joint-dependency of the link frames.
             & RK_SERIAL_SAVE_WITH_NAME(link_0_dep_frame)
             & RK_SERIAL_SAVE_WITH_NAME(link_1_dep_frame)
             & RK_SERIAL_SAVE_WITH_NAME(link_2_dep_frame)
             & RK_SERIAL_SAVE_WITH_NAME(link_3_dep_frame)
   // save all the links (massless).
             & RK_SERIAL_SAVE_WITH_NAME(link_0)
             & RK_SERIAL_SAVE_WITH_NAME(link_1)
             & RK_SERIAL_SAVE_WITH_NAME(link_2)
             & RK_SERIAL_SAVE_WITH_NAME(link_3)
   // save all the joint inertias (motor inertia).
             & RK_SERIAL_SAVE_WITH_NAME(track_joint_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3_inertia)
   // save all the joint actuators (which apply a driving force).
             & RK_SERIAL_SAVE_WITH_NAME(track_actuator)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1_actuator)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2_actuator)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3_actuator)
   // save all the link inertias (mass information of the links).
             & RK_SERIAL_SAVE_WITH_NAME(link_0_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(link_1_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(link_2_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(link_3_inertia);
  
};

void CRS_A465_2D_model_builder::save_limits_to_file(const std::string& aFileName) const {
  serialization::xml_oarchive complete_model_output(aFileName);
  serialization::oarchive& output_ref = complete_model_output;
  
  output_ref & RK_SERIAL_SAVE_WITH_NAME(joint_lower_bounds)
             & RK_SERIAL_SAVE_WITH_NAME(joint_upper_bounds)
             & RK_SERIAL_SAVE_WITH_NAME(joint_rate_limits)
             & RK_SERIAL_SAVE_WITH_NAME(preferred_posture);
};





shared_ptr< kte::kte_map_chain > CRS_A465_2D_model_builder::get_kinematics_kte_chain() const {
  shared_ptr< kte::kte_map_chain > CRS_A465_2D_kin_model(new kte::kte_map_chain("CRS_A465_2D_kin_model"), scoped_deleter());
  
  *CRS_A465_2D_kin_model << track_joint
                         << link_0
                         << arm_joint_1
                         << link_1
                         << arm_joint_2
                         << link_2
                         << arm_joint_3
                         << link_3;
  
  return CRS_A465_2D_kin_model;
};
    
    
shared_ptr< kte::kte_map_chain > CRS_A465_2D_model_builder::get_dynamics_kte_chain() const {
  shared_ptr< kte::kte_map_chain > CRS_A465_2D_dyn_model(new kte::kte_map_chain("CRS_A465_2D_dyn_model"), scoped_deleter());
  
  *CRS_A465_2D_dyn_model << track_actuator
                         << track_joint_inertia
                         << track_joint
                         << link_0
                         << link_0_inertia
                         << arm_joint_1_actuator
                         << arm_joint_1_inertia 
                         << arm_joint_1
                         << link_1
                         << link_1_inertia
                         << arm_joint_2_actuator
                         << arm_joint_2_inertia
                         << arm_joint_2
                         << link_2
                         << link_2_inertia
                         << arm_joint_3_actuator
                         << arm_joint_3_inertia
                         << arm_joint_3
                         << link_3
                         << link_3_inertia;
  
  return CRS_A465_2D_dyn_model;
};
    
    
shared_ptr< kte::mass_matrix_calc > CRS_A465_2D_model_builder::get_mass_matrix_calculator( int aInertiaSources ) const {
  shared_ptr< kte::mass_matrix_calc > CRS_A465_2D_M_calc(new kte::mass_matrix_calc("CRS_A465_2D_M_calc"), scoped_deleter());
  
  if( aInertiaSources & link_inertia ) {
    *CRS_A465_2D_M_calc << link_0_inertia
                        << link_1_inertia
                        << link_2_inertia
                        << link_3_inertia;
  };
  if( aInertiaSources & motor_inertia ) {
    *CRS_A465_2D_M_calc << track_joint_inertia
                        << arm_joint_1_inertia
                        << arm_joint_2_inertia
                        << arm_joint_3_inertia;
  };
    
  *CRS_A465_2D_M_calc << track_joint_coord
                      << arm_joint_1_coord
                      << arm_joint_2_coord
                      << arm_joint_3_coord;
  
  return CRS_A465_2D_M_calc; 
};
    
    
shared_ptr< kte::manipulator_kinematics_model > CRS_A465_2D_model_builder::get_manipulator_kin_model( int aDependentFrames ) const {
  
  shared_ptr< kte::manipulator_kinematics_model > CRS_A465_2D_kin_manip(new kte::manipulator_kinematics_model("CRS_A465_2D_kin_manip"), scoped_deleter());
  
  CRS_A465_2D_kin_manip->setModel(get_kinematics_kte_chain());
  
//Register joint coordinates:
  *CRS_A465_2D_kin_manip << track_joint_coord
                         << arm_joint_1_coord
                         << arm_joint_2_coord
                         << arm_joint_3_coord;
//Register joint-dependent link-frames:
  if( aDependentFrames & link_frames ) {
    *CRS_A465_2D_kin_manip << link_0_dep_frame  
                           << link_1_dep_frame
                           << link_2_dep_frame
                           << link_3_dep_frame;
  } else {  // end-effector only.
    *CRS_A465_2D_kin_manip << link_3_dep_frame;
  };
  
  return CRS_A465_2D_kin_manip;
};
    
    
shared_ptr< kte::manipulator_dynamics_model > CRS_A465_2D_model_builder::get_manipulator_dyn_model( int aDependentFrames ) const {
  
  shared_ptr< kte::manipulator_dynamics_model > CRS_A465_2D_dyn_manip(new kte::manipulator_dynamics_model("CRS_A465_2D_dyn_manip"), scoped_deleter());
  
  CRS_A465_2D_dyn_manip->setModel( get_dynamics_kte_chain() );
  
  
//Register joint coordinates:
  *CRS_A465_2D_dyn_manip << track_joint_coord
                         << arm_joint_1_coord
                         << arm_joint_2_coord
                         << arm_joint_3_coord;
//Register joint inertias:
  *CRS_A465_2D_dyn_manip << track_joint_inertia  
                         << arm_joint_1_inertia
                         << arm_joint_2_inertia
                         << arm_joint_3_inertia;
//Register link inertias (and dependent frames):
  *CRS_A465_2D_dyn_manip << link_0_inertia       
                         << link_1_inertia
                         << link_2_inertia
                         << link_3_inertia;
//Register joint actuators:
  *CRS_A465_2D_dyn_manip << track_actuator        
                         << arm_joint_1_actuator
                         << arm_joint_2_actuator
                         << arm_joint_3_actuator;
  
  return CRS_A465_2D_dyn_manip;
};
    




CRS_A465_2D_model_builder::rate_limited_joint_space_type CRS_A465_2D_model_builder::get_rl_joint_space() const {
  return joint_rate_limits.make_rl_joint_space(get_joint_space());
};

CRS_A465_2D_model_builder::joint_space_type CRS_A465_2D_model_builder::get_joint_space() const {
  typedef pp::joint_space_2nd_order<double>::type SingleJointSpace;
  typedef pp::line_segment_topology<double> LinSeg;
  return joint_space_type( arithmetic_tuple< 
                             SingleJointSpace,
                             SingleJointSpace,
                             SingleJointSpace,
                             SingleJointSpace
                           >(
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg,LinSeg,LinSeg
                               >(
                                 LinSeg("track_pos_space",joint_lower_bounds[0],joint_upper_bounds[0]),
                                 LinSeg("track_vel_space",-joint_rate_limits.gen_speed_limits[0],joint_rate_limits.gen_speed_limits[0]),
                                 LinSeg("track_acc_space",-joint_rate_limits.gen_accel_limits[0],joint_rate_limits.gen_accel_limits[0])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg,LinSeg,LinSeg
                               >(
                                 LinSeg("arm_joint_1_pos_space",joint_lower_bounds[1],joint_upper_bounds[1]),
                                 LinSeg("arm_joint_1_vel_space",-joint_rate_limits.gen_speed_limits[1],joint_rate_limits.gen_speed_limits[1]),
                                 LinSeg("arm_joint_1_acc_space",-joint_rate_limits.gen_accel_limits[1],joint_rate_limits.gen_accel_limits[1])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg,LinSeg,LinSeg
                               >(
                                 LinSeg("arm_joint_2_pos_space",joint_lower_bounds[2],joint_upper_bounds[2]),
                                 LinSeg("arm_joint_2_vel_space",-joint_rate_limits.gen_speed_limits[2],joint_rate_limits.gen_speed_limits[2]),
                                 LinSeg("arm_joint_2_acc_space",-joint_rate_limits.gen_accel_limits[2],joint_rate_limits.gen_accel_limits[2])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg,LinSeg,LinSeg
                               >(
                                 LinSeg("arm_joint_3_pos_space",joint_lower_bounds[3],joint_upper_bounds[3]),
                                 LinSeg("arm_joint_3_vel_space",-joint_rate_limits.gen_speed_limits[3],joint_rate_limits.gen_speed_limits[3]),
                                 LinSeg("arm_joint_3_acc_space",-joint_rate_limits.gen_accel_limits[3],joint_rate_limits.gen_accel_limits[3])
                               )
                             )
                           )
                         );
};

CRS_A465_2D_model_builder::end_effector_space_type CRS_A465_2D_model_builder::get_end_effector_space() const {
  
  typedef arithmetic_tuple< pp::hyperbox_topology< vect<double,2> >, 
                            pp::hyperball_topology< vect<double,2> >,
                            pp::hyperball_topology< vect<double,2> >
                          > TranslationTuple;
  typedef arithmetic_tuple< pp::line_segment_topology<double>, 
                            pp::line_segment_topology<double>,
                            pp::line_segment_topology<double>
                          > RotationTuple;
  
  typedef arithmetic_tuple_element<0, end_effector_space_type>::type TranslationDiffSpace;
  typedef arithmetic_tuple_element<1, end_effector_space_type>::type RotationDiffSpace;
  
  typedef arithmetic_tuple<TranslationDiffSpace, RotationDiffSpace> EESpaceTuple;
  
  return end_effector_space_type(
    EESpaceTuple(
      TranslationDiffSpace(
        TranslationTuple(
          pp::hyperbox_topology< vect<double,2> >("EE_pos_space",vect<double,2>(-1.1,-1.1),vect<double,2>(4.5,1.1)),
          pp::hyperball_topology< vect<double,2> >("EE_vel_space",vect<double,2>(0.0,0.0),5.0),
          pp::hyperball_topology< vect<double,2> >("EE_acc_space",vect<double,2>(0.0,0.0),25.0)
        )
      ),
      RotationDiffSpace(
        RotationTuple(
          pp::line_segment_topology<double>("EE_rotation_space",-3.14159265359,3.14159265359),
          pp::line_segment_topology<double>("EE_ang_vel_space",-10.0,10.0),
          pp::line_segment_topology<double>("EE_ang_acc_space",-100.0,100.0)
        )
      )
    )
  );
};







};
  
  
};








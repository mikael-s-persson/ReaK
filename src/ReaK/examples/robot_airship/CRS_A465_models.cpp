
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



#include "CRS_A465_models.hpp"

#include "serialization/xml_archiver.hpp"

#include "rtti/typed_primitives.hpp"

#include "optimization/optim_exceptions.hpp"


namespace ReaK {


namespace robot_airship {
  
  
void CRS_A465_model_builder::load_kte_from_archive(serialization::iarchive& aInput) {
  
  aInput
   // load the base frame (start of the track.
             & RK_SERIAL_LOAD_WITH_NAME(robot_base)
   // load all the joint coordinates.
             & RK_SERIAL_LOAD_WITH_NAME(track_joint_coord)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_1_coord)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_2_coord)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_3_coord)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_4_coord)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_5_coord)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_6_coord)
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
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_4_base)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_4_end)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_5_base)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_5_end)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_6_base)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_6_end)
   // load all the joint jacobian relationships.
             & RK_SERIAL_LOAD_WITH_NAME(track_joint_jacobian)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_1_jacobian)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_2_jacobian)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_3_jacobian)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_4_jacobian)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_5_jacobian)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_6_jacobian)
   // load all the joints.
             & RK_SERIAL_LOAD_WITH_NAME(track_joint)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_1)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_2)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_3)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_4)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_5)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_6)
   // load all the joint-dependency of the link frames.
             & RK_SERIAL_LOAD_WITH_NAME(link_0_dep_frame)
             & RK_SERIAL_LOAD_WITH_NAME(link_1_dep_frame)
             & RK_SERIAL_LOAD_WITH_NAME(link_2_dep_frame)
             & RK_SERIAL_LOAD_WITH_NAME(link_3_dep_frame)
             & RK_SERIAL_LOAD_WITH_NAME(link_4_dep_frame)
             & RK_SERIAL_LOAD_WITH_NAME(link_5_dep_frame)
             & RK_SERIAL_LOAD_WITH_NAME(link_6_dep_frame)
   // load all the links (massless).
             & RK_SERIAL_LOAD_WITH_NAME(link_0)
             & RK_SERIAL_LOAD_WITH_NAME(link_1)
             & RK_SERIAL_LOAD_WITH_NAME(link_2)
             & RK_SERIAL_LOAD_WITH_NAME(link_3)
             & RK_SERIAL_LOAD_WITH_NAME(link_4)
             & RK_SERIAL_LOAD_WITH_NAME(link_5)
             & RK_SERIAL_LOAD_WITH_NAME(link_6)
   // load all the joint inertias (motor inertia).
             & RK_SERIAL_LOAD_WITH_NAME(track_joint_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_1_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_2_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_3_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_4_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_5_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_6_inertia)
   // load all the joint actuators (which apply a driving force).
             & RK_SERIAL_LOAD_WITH_NAME(track_actuator)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_1_actuator)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_2_actuator)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_3_actuator)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_4_actuator)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_5_actuator)
             & RK_SERIAL_LOAD_WITH_NAME(arm_joint_6_actuator)
   // load all the link inertias (mass information of the links).
             & RK_SERIAL_LOAD_WITH_NAME(link_0_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(link_1_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(link_2_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(link_3_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(link_4_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(link_5_inertia)
             & RK_SERIAL_LOAD_WITH_NAME(link_6_inertia);
  
  A465_params.global_to_baseplate = robot_base->Position;
  A465_params.baseplate_to_shoulder_dist = link_1->PoseOffset().Position[2];  // "link_1" offset: vect<double,3>(0.0,0.0,0.3302),
  A465_params.shoulder_to_elbow_dist = link_2->PoseOffset().Position[2];      // "link_2" offset: vect<double,3>(0.0, 0.0, shoulder_to_elbow_dist),
  A465_params.elbow_to_joint_4_dist = link_3->PoseOffset().Position[2];       // "link_3" offset: vect<double,3>(0.0, 0.0, elbow_to_joint_4_dist),
  A465_params.joint_4_to_wrist_dist = link_4->PoseOffset().Position[2];       // "link_4" offset: vect<double,3>(0.0, 0.0, joint_4_to_wrist_dist),
  A465_params.wrist_to_flange_dist = link_5->PoseOffset().Position[2];        // "link_5" offset: vect<double,3>(0.0, 0.0, wrist_to_flange_dist),
};


void CRS_A465_model_builder::save_kte_to_archive(serialization::oarchive& aOutput) const {
  
  aOutput
   // save the base frame (start of the track.
             & RK_SERIAL_SAVE_WITH_NAME(robot_base)
   // save all the joint coordinates.
             & RK_SERIAL_SAVE_WITH_NAME(track_joint_coord)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1_coord)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2_coord)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3_coord)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_4_coord)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_5_coord)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_6_coord)
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
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_4_base)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_4_end)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_5_base)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_5_end)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_6_base)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_6_end)
   // save all the joint jacobian relationships.
             & RK_SERIAL_SAVE_WITH_NAME(track_joint_jacobian)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1_jacobian)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2_jacobian)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3_jacobian)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_4_jacobian)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_5_jacobian)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_6_jacobian)
   // save all the joints.
             & RK_SERIAL_SAVE_WITH_NAME(track_joint)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_4)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_5)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_6)
   // save all the joint-dependency of the link frames.
             & RK_SERIAL_SAVE_WITH_NAME(link_0_dep_frame)
             & RK_SERIAL_SAVE_WITH_NAME(link_1_dep_frame)
             & RK_SERIAL_SAVE_WITH_NAME(link_2_dep_frame)
             & RK_SERIAL_SAVE_WITH_NAME(link_3_dep_frame)
             & RK_SERIAL_SAVE_WITH_NAME(link_4_dep_frame)
             & RK_SERIAL_SAVE_WITH_NAME(link_5_dep_frame)
             & RK_SERIAL_SAVE_WITH_NAME(link_6_dep_frame)
   // save all the links (massless).
             & RK_SERIAL_SAVE_WITH_NAME(link_0)
             & RK_SERIAL_SAVE_WITH_NAME(link_1)
             & RK_SERIAL_SAVE_WITH_NAME(link_2)
             & RK_SERIAL_SAVE_WITH_NAME(link_3)
             & RK_SERIAL_SAVE_WITH_NAME(link_4)
             & RK_SERIAL_SAVE_WITH_NAME(link_5)
             & RK_SERIAL_SAVE_WITH_NAME(link_6)
   // save all the joint inertias (motor inertia).
             & RK_SERIAL_SAVE_WITH_NAME(track_joint_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_4_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_5_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_6_inertia)
   // save all the joint actuators (which apply a driving force).
             & RK_SERIAL_SAVE_WITH_NAME(track_actuator)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1_actuator)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2_actuator)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3_actuator)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_4_actuator)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_5_actuator)
             & RK_SERIAL_SAVE_WITH_NAME(arm_joint_6_actuator)
   // save all the link inertias (mass information of the links).
             & RK_SERIAL_SAVE_WITH_NAME(link_0_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(link_1_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(link_2_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(link_3_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(link_4_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(link_5_inertia)
             & RK_SERIAL_SAVE_WITH_NAME(link_6_inertia);
  
};


void CRS_A465_model_builder::load_kte_from_file(const std::string& aFileName) {
  
  serialization::xml_iarchive complete_model_input(aFileName);
  load_kte_from_archive(complete_model_input);
  
};


void CRS_A465_model_builder::load_limits_from_file(const std::string& aFileName) {
  serialization::xml_iarchive complete_model_input(aFileName);
  serialization::iarchive& input_ref = complete_model_input;
  
  input_ref & RK_SERIAL_LOAD_WITH_NAME(joint_lower_bounds)
            & RK_SERIAL_LOAD_WITH_NAME(joint_upper_bounds)
            & RK_SERIAL_LOAD_WITH_NAME(joint_rate_limits)
            & RK_SERIAL_LOAD_WITH_NAME(preferred_posture);
};
    
void CRS_A465_model_builder::create_from_preset(const CRS_A465_model_parameters& aParams) {
  
  
  //declare all the intermediate frames.
  robot_base       = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());
  track_joint_end  = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());
  arm_joint_1_base = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());
  arm_joint_1_end  = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());
  arm_joint_2_base = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());
  arm_joint_2_end  = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());
  arm_joint_3_base = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());
  arm_joint_3_end  = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());
  arm_joint_4_base = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());
  arm_joint_4_end  = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());
  arm_joint_5_base = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());
  arm_joint_5_end  = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());
  arm_joint_6_base = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());
  arm_joint_6_end  = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());
  arm_EE           = shared_ptr< frame_3D<double> >(new frame_3D<double>(), scoped_deleter());

  //declare all the joint coordinates.
  track_joint_coord = shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter());
  arm_joint_1_coord = shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter());
  arm_joint_2_coord = shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter());
  arm_joint_3_coord = shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter());
  arm_joint_4_coord = shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter());
  arm_joint_5_coord = shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter());
  arm_joint_6_coord = shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter());
  
  //declare all the joint jacobians.
  track_joint_jacobian = shared_ptr< jacobian_gen_3D<double> >(new jacobian_gen_3D<double>(), scoped_deleter());
  arm_joint_1_jacobian = shared_ptr< jacobian_gen_3D<double> >(new jacobian_gen_3D<double>(), scoped_deleter());
  arm_joint_2_jacobian = shared_ptr< jacobian_gen_3D<double> >(new jacobian_gen_3D<double>(), scoped_deleter());
  arm_joint_3_jacobian = shared_ptr< jacobian_gen_3D<double> >(new jacobian_gen_3D<double>(), scoped_deleter());
  arm_joint_4_jacobian = shared_ptr< jacobian_gen_3D<double> >(new jacobian_gen_3D<double>(), scoped_deleter());
  arm_joint_5_jacobian = shared_ptr< jacobian_gen_3D<double> >(new jacobian_gen_3D<double>(), scoped_deleter());
  arm_joint_6_jacobian = shared_ptr< jacobian_gen_3D<double> >(new jacobian_gen_3D<double>(), scoped_deleter());
  
  //set the absolute position of the base and add gravity (z-axis pointing up!) (x-axis pointing forward).
  robot_base->Acceleration = vect<double,3>(0.0,0.0,9.81); //put gravity acceleration on base of the robot
  robot_base->Position = aParams.global_to_baseplate; //put the base of the robot at the near end of the track (global frame is at the far end). Default: robot_base->Position = vect<double,3>(0.0,-3.3,0.3);
  robot_base->Quat = axis_angle<double>(M_PI * 0.5, vect<double,3>(0.0,0.0,1.0)); // align the x-axis along the track.
  
  //create revolute joint
  track_joint = shared_ptr< kte::prismatic_joint_3D >(new kte::prismatic_joint_3D("track_joint",
                                                                              track_joint_coord,
                                                                              vect<double,3>(1.0,0.0,0.0),
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
  link_0 = shared_ptr< kte::rigid_link_3D >(new kte::rigid_link_3D("link_0",
                                                               track_joint_end,
                                                               arm_joint_1_base,
                                                               pose_3D<double>(weak_ptr<pose_3D<double> >(),
                                                                               vect<double,3>(0.0,0.0,0.0),
                                                                               quaternion<double>())),
                                             scoped_deleter());
  
  //create arm-base inertia of 
  link_0_dep_frame = shared_ptr< kte::joint_dependent_frame_3D >(new kte::joint_dependent_frame_3D(arm_joint_1_base),
								    scoped_deleter());
  link_0_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_0_inertia = shared_ptr< kte::inertia_3D >(new kte::inertia_3D("link_0_inertia",
                                                                 link_0_dep_frame,
                                                                 1.0,
                                                                 mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0)),
                                                 scoped_deleter());
  
  //create revolute joint
  arm_joint_1 = shared_ptr< kte::revolute_joint_3D >(new kte::revolute_joint_3D("arm_joint_1",
                                                                            arm_joint_1_coord,
                                                                            vect<double,3>(0.0,0.0,1.0),
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
  link_1 = shared_ptr< kte::rigid_link_3D >(new kte::rigid_link_3D("link_1",
                                                                arm_joint_1_end,
                                                                arm_joint_2_base,
                                                                pose_3D<double>(weak_ptr<pose_3D<double> >(),
                                                                                vect<double,3>(0.0,0.0, aParams.baseplate_to_shoulder_dist /*0.3302*/ ),
                                                                                quaternion<double>())),
                                             scoped_deleter());
  
  //create link1 inertia 
  link_1_dep_frame = shared_ptr< kte::joint_dependent_frame_3D >(new kte::joint_dependent_frame_3D(arm_joint_2_base),
								    scoped_deleter());
  link_1_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_1_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  link_1_inertia = shared_ptr< kte::inertia_3D >(new kte::inertia_3D("link_1_inertia",
                                                                 link_1_dep_frame,
                                                                 1.0,
                                                                 mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0)),
                                                  scoped_deleter());

  //create revolute joint
  arm_joint_2 = shared_ptr< kte::revolute_joint_3D >(new kte::revolute_joint_3D("arm_joint_2",
                                                                             arm_joint_2_coord,
                                                                             vect<double,3>(0.0,-1.0,0.0),
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
  link_2 = shared_ptr< kte::rigid_link_3D >(new kte::rigid_link_3D("link_2",
                                                                arm_joint_2_end,
                                                                arm_joint_3_base,
                                                                pose_3D<double>(weak_ptr<pose_3D<double> >(),
                                                                                vect<double,3>(0.0,0.0, aParams.shoulder_to_elbow_dist /*0.3048*/ ),
                                                                                quaternion<double>())),
                                              scoped_deleter());
  
  //create inertia
  link_2_dep_frame = shared_ptr< kte::joint_dependent_frame_3D >(new kte::joint_dependent_frame_3D(arm_joint_3_base),
								    scoped_deleter());
  link_2_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_2_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  link_2_dep_frame->add_joint(arm_joint_2_coord,arm_joint_2_jacobian);
  link_2_inertia = shared_ptr< kte::inertia_3D >(new kte::inertia_3D("link_2_inertia",
                                                                 link_2_dep_frame,
                                                                 1.0,
                                                                 mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0)),
                                                  scoped_deleter());

  //create revolute joint
  arm_joint_3 = shared_ptr< kte::revolute_joint_3D >(new kte::revolute_joint_3D("arm_joint_3",
                                                                            arm_joint_3_coord,
                                                                            vect<double,3>(0.0,-1.0,0.0),
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
  link_3 = shared_ptr< kte::rigid_link_3D >(new kte::rigid_link_3D("link_3",
                                                               arm_joint_3_end,
                                                               arm_joint_4_base,
                                                               pose_3D<double>(weak_ptr<pose_3D<double> >(),
                                                                               vect<double,3>(0.0,0.0, aParams.elbow_to_joint_4_dist /*0.1500*/ ),
                                                                               quaternion<double>())),
                                              scoped_deleter());
  
  //create inertia
  link_3_dep_frame = shared_ptr< kte::joint_dependent_frame_3D >(new kte::joint_dependent_frame_3D(arm_joint_4_base),
								    scoped_deleter());
  link_3_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_3_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  link_3_dep_frame->add_joint(arm_joint_2_coord,arm_joint_2_jacobian);
  link_3_dep_frame->add_joint(arm_joint_3_coord,arm_joint_3_jacobian);
  link_3_inertia = shared_ptr< kte::inertia_3D >(new kte::inertia_3D("link_3_inertia",
                                                                 link_3_dep_frame,
                                                                 1.0,
                                                                 mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0)),
                                                  scoped_deleter());

  //create revolute joint
  arm_joint_4 = shared_ptr< kte::revolute_joint_3D >(new kte::revolute_joint_3D("arm_joint_4",
                                                                            arm_joint_4_coord,
                                                                            vect<double,3>(0.0,0.0,1.0),
                                                                            arm_joint_4_base,
                                                                            arm_joint_4_end,
                                                                            arm_joint_4_jacobian),
                                                      scoped_deleter());
  
  //create motor inertia
  shared_ptr< kte::joint_dependent_gen_coord > arm_joint_4_dep_coord(new kte::joint_dependent_gen_coord(arm_joint_4_coord),
									  scoped_deleter());
  arm_joint_4_dep_coord->add_joint(arm_joint_4_coord,
                                   shared_ptr<jacobian_gen_gen<double> >(new jacobian_gen_gen<double>(1.0,0.0), 
										   scoped_deleter()));
  arm_joint_4_inertia = shared_ptr< kte::inertia_gen >(new kte::inertia_gen("arm_joint_4_inertia",
                                                                        arm_joint_4_dep_coord,
                                                                        1.0), 
                                                        scoped_deleter());
  
  //create force actuator
  arm_joint_4_actuator = shared_ptr< kte::driving_actuator_gen >(new kte::driving_actuator_gen("arm_joint_4_actuator",
                                                                                           arm_joint_4_coord,
                                                                                           arm_joint_4),
                                                                  scoped_deleter());
  
  //create link 
  link_4 = shared_ptr< kte::rigid_link_3D >(new kte::rigid_link_3D("link_4",
                                                               arm_joint_4_end,
                                                               arm_joint_5_base,
                                                               pose_3D<double>(weak_ptr<pose_3D<double> >(),
                                                                               vect<double,3>(0.0,0.0, aParams.joint_4_to_wrist_dist /*0.1802*/),
                                                                               quaternion<double>())),
                                             scoped_deleter());
  
  //create inertia
  link_4_dep_frame = shared_ptr< kte::joint_dependent_frame_3D >(new kte::joint_dependent_frame_3D(arm_joint_5_base),
								    scoped_deleter());
  link_4_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_4_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  link_4_dep_frame->add_joint(arm_joint_2_coord,arm_joint_2_jacobian);
  link_4_dep_frame->add_joint(arm_joint_3_coord,arm_joint_3_jacobian);
  link_4_dep_frame->add_joint(arm_joint_4_coord,arm_joint_4_jacobian);
  link_4_inertia = shared_ptr< kte::inertia_3D >(new kte::inertia_3D("link_4_inertia",
                                                                 link_4_dep_frame,
                                                                 1.0,
                                                                 mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0)),
                                                  scoped_deleter());

  //create revolute joint
  arm_joint_5 = shared_ptr< kte::revolute_joint_3D >(new kte::revolute_joint_3D("arm_joint_5",
                                                                            arm_joint_5_coord,
                                                                            vect<double,3>(0.0,-1.0,0.0),
                                                                            arm_joint_5_base,
                                                                            arm_joint_5_end,
                                                                            arm_joint_5_jacobian),
                                                      scoped_deleter());
  
  //create motor inertia
  shared_ptr< kte::joint_dependent_gen_coord > arm_joint_5_dep_coord(new kte::joint_dependent_gen_coord(arm_joint_5_coord),
									  scoped_deleter());
  arm_joint_5_dep_coord->add_joint(arm_joint_5_coord,
                                   shared_ptr<jacobian_gen_gen<double> >(new jacobian_gen_gen<double>(1.0,0.0), 
										   scoped_deleter()));
  arm_joint_5_inertia = shared_ptr< kte::inertia_gen >(new kte::inertia_gen("arm_joint_5_inertia",
									arm_joint_5_dep_coord,
									1.0),
                                                        scoped_deleter());
  
  //create force actuator
  arm_joint_5_actuator = shared_ptr< kte::driving_actuator_gen >(new kte::driving_actuator_gen("arm_joint_5_actuator",
											   arm_joint_5_coord,
											   arm_joint_5),
                                                                  scoped_deleter());
  
  //create link 
  link_5 = shared_ptr< kte::rigid_link_3D >(new kte::rigid_link_3D("link_5",
                                                               arm_joint_5_end,
                                                               arm_joint_6_base,
                                                               pose_3D<double>(weak_ptr<pose_3D<double> >(),
                                                                               vect<double,3>(0.0,0.0, aParams.wrist_to_flange_dist /*0.0762*/ ),
                                                                               quaternion<double>())),
                                             scoped_deleter());
  
  //create inertia
  link_5_dep_frame = shared_ptr< kte::joint_dependent_frame_3D >(new kte::joint_dependent_frame_3D(arm_joint_6_base),
								    scoped_deleter());
  link_5_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_5_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  link_5_dep_frame->add_joint(arm_joint_2_coord,arm_joint_2_jacobian);
  link_5_dep_frame->add_joint(arm_joint_3_coord,arm_joint_3_jacobian);
  link_5_dep_frame->add_joint(arm_joint_4_coord,arm_joint_4_jacobian);
  link_5_dep_frame->add_joint(arm_joint_5_coord,arm_joint_5_jacobian);
  link_5_inertia = shared_ptr< kte::inertia_3D >(new kte::inertia_3D("link_5_inertia",
                                                                 link_5_dep_frame,
                                                                 1.0,
                                                                 mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0)),
                                                  scoped_deleter());
  
  //create revolute joint
  arm_joint_6 = shared_ptr< kte::revolute_joint_3D >(new kte::revolute_joint_3D("arm_joint_6",
                                                                            arm_joint_6_coord,
                                                                            vect<double,3>(0.0,0.0,1.0),
                                                                            arm_joint_6_base,
                                                                            arm_joint_6_end,
                                                                            arm_joint_6_jacobian),
                                                      scoped_deleter());
  
  //create motor inertia
  shared_ptr< kte::joint_dependent_gen_coord > arm_joint_6_dep_coord(new kte::joint_dependent_gen_coord(arm_joint_6_coord),
									  scoped_deleter());
  arm_joint_6_dep_coord->add_joint(arm_joint_6_coord,
                                   shared_ptr<jacobian_gen_gen<double> >(new jacobian_gen_gen<double>(1.0,0.0), 
										   scoped_deleter()));
  arm_joint_6_inertia = shared_ptr< kte::inertia_gen >(new kte::inertia_gen("arm_joint_6_inertia",
									arm_joint_6_dep_coord,
                                                                        1.0),
                                                        scoped_deleter());
  
  //create force actuator
  arm_joint_6_actuator = shared_ptr< kte::driving_actuator_gen >(new kte::driving_actuator_gen("arm_joint_6_actuator",
                                                                                           arm_joint_6_coord,
                                                                                           arm_joint_6),
                                                                  scoped_deleter());
  
  //create link 
  link_6 = shared_ptr< kte::rigid_link_3D >(new kte::rigid_link_3D("link_6",
                                                                arm_joint_6_end,
                                                                arm_EE,
                                                                pose_3D<double>(weak_ptr<pose_3D<double> >(),
                                                                                vect<double,3>(0.0,0.0,0.0),
                                                                                quaternion<double>())),
                                             scoped_deleter());
  
  //create inertia
  link_6_dep_frame = shared_ptr< kte::joint_dependent_frame_3D >(new kte::joint_dependent_frame_3D(arm_EE),
								    scoped_deleter());
  link_6_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_6_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  link_6_dep_frame->add_joint(arm_joint_2_coord,arm_joint_2_jacobian);
  link_6_dep_frame->add_joint(arm_joint_3_coord,arm_joint_3_jacobian);
  link_6_dep_frame->add_joint(arm_joint_4_coord,arm_joint_4_jacobian);
  link_6_dep_frame->add_joint(arm_joint_5_coord,arm_joint_5_jacobian);
  link_6_dep_frame->add_joint(arm_joint_6_coord,arm_joint_6_jacobian);
  link_6_inertia = shared_ptr< kte::inertia_3D >(new kte::inertia_3D("link_6_inertia",
                                                                 link_6_dep_frame,
                                                                 1.0,
                                                                 mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0)),
                                                  scoped_deleter());
  
  
  joint_lower_bounds.resize(7);
  joint_lower_bounds[0] = 0.0;
  joint_lower_bounds[1] = -3.05432619099;
  joint_lower_bounds[2] = -1.57079632679;
  joint_lower_bounds[3] = -1.91986217719;
  joint_lower_bounds[4] = -3.14159265359;
  joint_lower_bounds[5] = -1.83259571459;
  joint_lower_bounds[6] = -3.14159265359;
  
  joint_upper_bounds.resize(7);
  joint_upper_bounds[0] = 3.0;
  joint_upper_bounds[1] = 3.05432619099;
  joint_upper_bounds[2] = 1.57079632679;
  joint_upper_bounds[3] = 1.91986217719;
  joint_upper_bounds[4] = 33.14159265359;
  joint_upper_bounds[5] = 1.83259571459;
  joint_upper_bounds[6] = 3.14159265359;

  joint_rate_limits.gen_speed_limits.resize(7);
  joint_rate_limits.gen_speed_limits[0] = 0.8;
  joint_rate_limits.gen_speed_limits[1] = 3.14159265359;
  joint_rate_limits.gen_speed_limits[2] = 3.14159265359;
  joint_rate_limits.gen_speed_limits[3] = 3.14159265359;
  joint_rate_limits.gen_speed_limits[4] = 2.98451302091;
  joint_rate_limits.gen_speed_limits[5] = 3.01941960595;
  joint_rate_limits.gen_speed_limits[6] = 2.98451302091;
  
  joint_rate_limits.gen_accel_limits.resize(7);
  joint_rate_limits.gen_accel_limits[0] = 3.0;
  joint_rate_limits.gen_accel_limits[1] = 12.5663706144;
  joint_rate_limits.gen_accel_limits[2] = 12.5663706144;
  joint_rate_limits.gen_accel_limits[3] = 12.5663706144;
  joint_rate_limits.gen_accel_limits[4] = 24.9582083035;
  joint_rate_limits.gen_accel_limits[5] = 24.9582083035;
  joint_rate_limits.gen_accel_limits[6] = 24.9582083035;
  
  joint_rate_limits.gen_jerk_limits.resize(7);
  joint_rate_limits.gen_jerk_limits[0] = 12.0;
  joint_rate_limits.gen_jerk_limits[1] = 125.663706144;
  joint_rate_limits.gen_jerk_limits[2] = 125.663706144;
  joint_rate_limits.gen_jerk_limits[3] = 125.663706144;
  joint_rate_limits.gen_jerk_limits[4] = 249.582083035;
  joint_rate_limits.gen_jerk_limits[5] = 249.582083035;
  joint_rate_limits.gen_jerk_limits[6] = 249.582083035;
  
  joint_rate_limits.frame2D_speed_limits.resize(0);
  joint_rate_limits.frame2D_accel_limits.resize(0);
  joint_rate_limits.frame2D_jerk_limits.resize(0);
  joint_rate_limits.frame3D_speed_limits.resize(0);
  joint_rate_limits.frame3D_accel_limits.resize(0);
  joint_rate_limits.frame3D_jerk_limits.resize(0);
  
  preferred_posture.resize(7);
  preferred_posture[0] = 1.5;
  preferred_posture[1] = 0.0;
  preferred_posture[2] = 0.0;
  preferred_posture[3] = 1.57079632679;
  preferred_posture[4] = 0.0;
  preferred_posture[5] = 0.0;
  preferred_posture[6] = 0.0;
  
  A465_params = aParams;
  
};
    
void CRS_A465_model_builder::save_kte_to_file(const std::string& aFileName) const {
  
  serialization::xml_oarchive complete_model_output(aFileName);
  save_kte_to_archive(complete_model_output);
  
};

void CRS_A465_model_builder::save_limits_to_file(const std::string& aFileName) const {
  serialization::xml_oarchive complete_model_output(aFileName);
  serialization::oarchive& output_ref = complete_model_output;
  
  output_ref & RK_SERIAL_SAVE_WITH_NAME(joint_lower_bounds)
             & RK_SERIAL_SAVE_WITH_NAME(joint_upper_bounds)
             & RK_SERIAL_SAVE_WITH_NAME(joint_rate_limits)
             & RK_SERIAL_SAVE_WITH_NAME(preferred_posture);
};





shared_ptr< kte::kte_map_chain > CRS_A465_model_builder::get_kinematics_kte_chain() const {
  shared_ptr< kte::kte_map_chain > CRS_A465_kin_model(new kte::kte_map_chain("CRS_A465_kin_model"), scoped_deleter());
  
  *CRS_A465_kin_model << track_joint
                      << link_0
                      << arm_joint_1
                      << link_1
                      << arm_joint_2
                      << link_2
                      << arm_joint_3
                      << link_3
                      << arm_joint_4
                      << link_4
                      << arm_joint_5
                      << link_5
                      << arm_joint_6
                      << link_6;
  
  return CRS_A465_kin_model;
};
    
    
shared_ptr< kte::kte_map_chain > CRS_A465_model_builder::get_dynamics_kte_chain() const {
  shared_ptr< kte::kte_map_chain > CRS_A465_dyn_model(new kte::kte_map_chain("CRS_A465_dyn_model"), scoped_deleter());
  
  *CRS_A465_dyn_model << track_actuator
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
                      << link_3_inertia
                      << arm_joint_4_actuator
                      << arm_joint_4_inertia
                      << arm_joint_4
                      << link_4
                      << link_4_inertia
                      << arm_joint_5_actuator
                      << arm_joint_5_inertia
                      << arm_joint_5
                      << link_5
                      << link_5_inertia
                      << arm_joint_6_actuator
                      << arm_joint_6_inertia
                      << arm_joint_6
                      << link_6
                      << link_6_inertia;
  
  return CRS_A465_dyn_model;
};
    
    
shared_ptr< kte::mass_matrix_calc > CRS_A465_model_builder::get_mass_matrix_calculator( int aInertiaSources ) const {
  shared_ptr< kte::mass_matrix_calc > CRS_A465_M_calc(new kte::mass_matrix_calc("CRS_A465_M_calc"), scoped_deleter());
  
  if( aInertiaSources & link_inertia ) {
    *CRS_A465_M_calc << link_0_inertia
                     << link_1_inertia
                     << link_2_inertia
                     << link_3_inertia
                     << link_4_inertia
                     << link_5_inertia
                     << link_6_inertia;
  };
  if( aInertiaSources & motor_inertia ) {
    *CRS_A465_M_calc << track_joint_inertia
                     << arm_joint_1_inertia
                     << arm_joint_2_inertia
                     << arm_joint_3_inertia
                     << arm_joint_4_inertia
                     << arm_joint_5_inertia
                     << arm_joint_6_inertia;
  };
    
  *CRS_A465_M_calc << track_joint_coord
                   << arm_joint_1_coord
                   << arm_joint_2_coord
                   << arm_joint_3_coord
                   << arm_joint_4_coord
                   << arm_joint_5_coord
                   << arm_joint_6_coord;
  
  return CRS_A465_M_calc; 
};
    
    
shared_ptr< kte::manipulator_kinematics_model > CRS_A465_model_builder::get_manipulator_kin_model( int aDependentFrames ) const {
  
  shared_ptr< kte::manipulator_kinematics_model > CRS_A465_kin_manip(new kte::manipulator_kinematics_model("CRS_A465_kin_manip"), scoped_deleter());
  
  CRS_A465_kin_manip->setModel(get_kinematics_kte_chain());
  
//Register joint coordinates:
  *CRS_A465_kin_manip << track_joint_coord
                      << arm_joint_1_coord
                      << arm_joint_2_coord
                      << arm_joint_3_coord
                      << arm_joint_4_coord
                      << arm_joint_5_coord
                      << arm_joint_6_coord;
//Register joint-dependent link-frames:
  if( aDependentFrames & link_frames ) {
    *CRS_A465_kin_manip << link_0_dep_frame  
                        << link_1_dep_frame
                        << link_2_dep_frame
                        << link_3_dep_frame
                        << link_4_dep_frame
                        << link_5_dep_frame
                        << link_6_dep_frame;
  } else {  // end-effector only.
    *CRS_A465_kin_manip << link_6_dep_frame;
  };
  
  return CRS_A465_kin_manip;
};
    
    
shared_ptr< kte::manipulator_dynamics_model > CRS_A465_model_builder::get_manipulator_dyn_model( int aDependentFrames ) const {
  
  shared_ptr< kte::manipulator_dynamics_model > CRS_A465_dyn_manip(new kte::manipulator_dynamics_model("CRS_A465_dyn_manip"), scoped_deleter());
  
  CRS_A465_dyn_manip->setModel( get_dynamics_kte_chain() );
  
  
//Register joint coordinates:
  *CRS_A465_dyn_manip << track_joint_coord
                      << arm_joint_1_coord
                      << arm_joint_2_coord
                      << arm_joint_3_coord
                      << arm_joint_4_coord
                      << arm_joint_5_coord
                      << arm_joint_6_coord;
//Register joint inertias:
  *CRS_A465_dyn_manip << track_joint_inertia  
                      << arm_joint_1_inertia
                      << arm_joint_2_inertia
                      << arm_joint_3_inertia
                      << arm_joint_4_inertia
                      << arm_joint_5_inertia
                      << arm_joint_6_inertia;
//Register link inertias (and dependent frames):
  *CRS_A465_dyn_manip << link_0_inertia       
                      << link_1_inertia
                      << link_2_inertia
                      << link_3_inertia
                      << link_4_inertia
                      << link_5_inertia
                      << link_6_inertia;
//Register joint actuators:
  *CRS_A465_dyn_manip << track_actuator        
                      << arm_joint_1_actuator
                      << arm_joint_2_actuator
                      << arm_joint_3_actuator
                      << arm_joint_4_actuator
                      << arm_joint_5_actuator
                      << arm_joint_6_actuator;
  
  return CRS_A465_dyn_manip;
};
    




CRS_A465_model_builder::rate_limited_joint_space_type CRS_A465_model_builder::get_rl_joint_space() const {
  return joint_rate_limits.make_rl_joint_space(get_joint_space());
};

CRS_A465_model_builder::joint_space_type CRS_A465_model_builder::get_joint_space() const {
  typedef pp::joint_space_2nd_order<double>::type SingleJointSpace;
  typedef pp::line_segment_topology<double> LinSeg;
  return joint_space_type( arithmetic_tuple< 
                             SingleJointSpace,
                             SingleJointSpace,
                             SingleJointSpace,
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
	                     ),
	                     SingleJointSpace(
	                       arithmetic_tuple<
	                         LinSeg,LinSeg,LinSeg
	                       >(
				 LinSeg("arm_joint_4_pos_space",joint_lower_bounds[4],joint_upper_bounds[4]),
				 LinSeg("arm_joint_4_vel_space",-joint_rate_limits.gen_speed_limits[4],joint_rate_limits.gen_speed_limits[4]),
				 LinSeg("arm_joint_4_acc_space",-joint_rate_limits.gen_accel_limits[4],joint_rate_limits.gen_accel_limits[4])
			       )
	                     ),
	                     SingleJointSpace(
	                       arithmetic_tuple<
	                         LinSeg,LinSeg,LinSeg
	                       >(
				 LinSeg("arm_joint_5_pos_space",joint_lower_bounds[5],joint_upper_bounds[5]),
				 LinSeg("arm_joint_5_vel_space",-joint_rate_limits.gen_speed_limits[5],joint_rate_limits.gen_speed_limits[5]),
				 LinSeg("arm_joint_5_acc_space",-joint_rate_limits.gen_accel_limits[5],joint_rate_limits.gen_accel_limits[5])
			       )
	                     ),
	                     SingleJointSpace(
	                       arithmetic_tuple<
	                         LinSeg,LinSeg,LinSeg
	                       >(
				 LinSeg("arm_joint_6_pos_space",joint_lower_bounds[6],joint_upper_bounds[6]),
				 LinSeg("arm_joint_6_vel_space",-joint_rate_limits.gen_speed_limits[6],joint_rate_limits.gen_speed_limits[6]),
				 LinSeg("arm_joint_6_acc_space",-joint_rate_limits.gen_accel_limits[6],joint_rate_limits.gen_accel_limits[6])
			       )
	                     )
			   )
	                 );
};

CRS_A465_model_builder::end_effector_space_type CRS_A465_model_builder::get_end_effector_space() const {
  
  typedef arithmetic_tuple< pp::hyperbox_topology< vect<double,3> >, 
                            pp::hyperball_topology< vect<double,3> >,
			    pp::hyperball_topology< vect<double,3> >
			  > TranslationTuple;
  typedef arithmetic_tuple< pp::quaternion_topology<double>, 
                            pp::ang_velocity_3D_topology<double>,
			    pp::ang_accel_3D_topology<double>
			  > RotationTuple;
  
  typedef arithmetic_tuple_element<0, end_effector_space_type>::type TranslationDiffSpace;
  typedef arithmetic_tuple_element<1, end_effector_space_type>::type RotationDiffSpace;
  
  typedef arithmetic_tuple<TranslationDiffSpace, RotationDiffSpace> EESpaceTuple;
  
  return end_effector_space_type(
    EESpaceTuple(
      TranslationDiffSpace(
	TranslationTuple(
	  pp::hyperbox_topology< vect<double,3> >("EE_pos_space",vect<double,3>(-1.1,-1.1,0.0),vect<double,3>(4.5,1.1,1.5)),
	  pp::hyperball_topology< vect<double,3> >("EE_vel_space",vect<double,3>(0.0,0.0,0.0),5.0,mat<double,mat_structure::identity>(3)),
	  pp::hyperball_topology< vect<double,3> >("EE_acc_space",vect<double,3>(0.0,0.0,0.0),25.0,mat<double,mat_structure::identity>(3))
	)
      ),
      RotationDiffSpace(
	RotationTuple(
	  pp::quaternion_topology<double>("EE_rotation_space"),
	  pp::ang_velocity_3D_topology<double>("EE_ang_vel_space",10.0),
	  pp::ang_accel_3D_topology<double>("EE_ang_acc_space",100.0)
	)
      )
    )
  );
};


pose_3D<double> CRS_A465_model_builder::compute_direct_kinematics(const vect_n<double>& joint_positions) const {
  
  pose_3D<double> flange;
  
  /* calculate individual rotations */
  //quaternion<double> qm1 = robot_base->Quat;
  quaternion<double>::zrot q1( joint_positions[1]);
  quaternion<double>::yrot q2(-joint_positions[2]);
  quaternion<double>::yrot q3(-joint_positions[3]);
  quaternion<double>::zrot q4( joint_positions[4]);
  quaternion<double>::yrot q5(-joint_positions[5]);
  quaternion<double>::zrot q6( joint_positions[6]);
  
  flange.Position  = robot_base->Position;
  flange.Position += (flange.Quat  = robot_base->Quat) * (A465_params.baseplate_to_shoulder_dist * vect_k);
  flange.Position +=  flange.Quat *                      (joint_positions[0] * vect_i);
  flange.Position += (flange.Quat *= q1 * q2) *          (A465_params.shoulder_to_elbow_dist * vect_k);
  flange.Position += (flange.Quat *= q3) *              ((A465_params.elbow_to_joint_4_dist + A465_params.joint_4_to_wrist_dist) * vect_k);
  flange.Position += (flange.Quat *= q4 * q5 * q6) *     (A465_params.wrist_to_flange_dist * vect_k);
  
  /*
  robot_base->Acceleration = vect<double,3>(0.0,0.0,9.81); //put gravity acceleration on base of the robot
  robot_base->Position = A465_params.global_to_baseplate; //put the base of the robot at the near end of the track (global frame is at the far end).
  robot_base->Quat = axis_angle<double>(M_PI * 0.5, vect<double,3>(0.0,0.0,1.0)); // align the x-axis along the track.
  "track_joint" axis: vect<double,3>(1.0,0.0,0.0),
  "link_0" offset: vect<double,3>(0.0,0.0,0.0),
  "arm_joint_1" axis: vect<double,3>(0.0,0.0,1.0),
  "link_1" offset: vect<double,3>(0.0,0.0,A465_params.baseplate_to_shoulder_dist),
  "arm_joint_2" axis: vect<double,3>(0.0,-1.0,0.0),
  "link_2" offset: vect<double,3>(0.0,0.0,A465_params.shoulder_to_elbow_dist),
  "arm_joint_3" axis: vect<double,3>(0.0,-1.0,0.0),
  "link_3" offset: vect<double,3>(0.0,0.0,A465_params.elbow_to_joint_4_dist),
  "arm_joint_4" axis: vect<double,3>(0.0,0.0,1.0),
  "link_4" offset: vect<double,3>(0.0,0.0,A465_params.joint_4_to_wrist_dist),
  "arm_joint_5" axis: vect<double,3>(0.0,-1.0,0.0),
  "link_5" offset: vect<double,3>(0.0,0.0,A465_params.wrist_to_flange_dist),
  "arm_joint_6" axis: vect<double,3>(0.0,0.0,1.0),
  "link_6" offset: vect<double,3>(0.0,0.0,0.0)
  */
  
  return flange;
};

frame_3D<double> CRS_A465_model_builder::compute_direct_kinematics_with_vel(const vect_n<double>& joint_states) const {
  
  frame_3D<double> flange;
  
  /* calculate individual rotations */
  //quaternion<double> qm1 = robot_base->Quat;
  quaternion<double>::zrot q1( joint_states[1]);
  quaternion<double>::yrot q2(-joint_states[2]);
  quaternion<double>::yrot q3(-joint_states[3]);
  quaternion<double>::zrot q4( joint_states[4]);
  quaternion<double>::yrot q5(-joint_states[5]);
  quaternion<double>::zrot q6( joint_states[6]);
  
  flange.Position  = robot_base->Position;
  flange.Position += (flange.Quat  = robot_base->Quat) * (A465_params.baseplate_to_shoulder_dist * vect_k);
  flange.Position +=  flange.Quat *                      (joint_states[0] * vect_i);
  flange.Position += (flange.Quat *= q1 * q2) *          (A465_params.shoulder_to_elbow_dist * vect_k);
  flange.Position += (flange.Quat *= q3) *              ((A465_params.elbow_to_joint_4_dist + A465_params.joint_4_to_wrist_dist) * vect_k);
  flange.Position += (flange.Quat *= q4 * q5 * q6) *     (A465_params.wrist_to_flange_dist * vect_k);
  
  /*
  robot_base->Acceleration = vect<double,3>(0.0,0.0,9.81); //put gravity acceleration on base of the robot
  robot_base->Position = A465_params.global_to_baseplate; //put the base of the robot at the near end of the track (global frame is at the far end).
  robot_base->Quat = axis_angle<double>(M_PI * 0.5, vect<double,3>(0.0,0.0,1.0)); // align the x-axis along the track.
  "track_joint" axis: vect<double,3>(1.0,0.0,0.0),
  "link_0" offset: vect<double,3>(0.0,0.0,0.0),
  "arm_joint_1" axis: vect<double,3>(0.0,0.0,1.0),
  "link_1" offset: vect<double,3>(0.0,0.0,A465_params.baseplate_to_shoulder_dist),
  "arm_joint_2" axis: vect<double,3>(0.0,-1.0,0.0),
  "link_2" offset: vect<double,3>(0.0,0.0,A465_params.shoulder_to_elbow_dist),
  "arm_joint_3" axis: vect<double,3>(0.0,-1.0,0.0),
  "link_3" offset: vect<double,3>(0.0,0.0,A465_params.elbow_to_joint_4_dist),
  "arm_joint_4" axis: vect<double,3>(0.0,0.0,1.0),
  "link_4" offset: vect<double,3>(0.0,0.0,A465_params.joint_4_to_wrist_dist),
  "arm_joint_5" axis: vect<double,3>(0.0,-1.0,0.0),
  "link_5" offset: vect<double,3>(0.0,0.0,A465_params.wrist_to_flange_dist),
  "arm_joint_6" axis: vect<double,3>(0.0,0.0,1.0),
  "link_6" offset: vect<double,3>(0.0,0.0,0.0)
  */
  
  return flange;
};

static double clamp_to_pi_range(double a) {
  return (a > M_PI ? a - 2.0 * M_PI : (a < -M_PI ? a + 2.0 * M_PI : a) );
};


vect_n<double> CRS_A465_model_builder::compute_inverse_kinematics(const pose_3D<double>& EE_pose) const {
  using std::sin; using std::cos; using std::fabs; 
  using std::atan2; using std::sqrt; using std::pow;
  
  const double pos_epsilon = 1e-3;
  const double extend_epsilon = .06;
  
  quaternion<double>::zrot gl_to_track_rot(-0.5 * M_PI);
  vect<double,3> EE_z_axis = gl_to_track_rot * (EE_pose.Quat * vect_k);
  vect<double,3> EE_y_axis = gl_to_track_rot * (EE_pose.Quat * vect_j);
  vect<double,3> wrist_pos = gl_to_track_rot * (EE_pose.Position - A465_params.global_to_baseplate)
                           - A465_params.wrist_to_flange_dist * EE_z_axis 
                           - A465_params.baseplate_to_shoulder_dist * vect_k;
  double elbow_to_wrist_dist = A465_params.elbow_to_joint_4_dist + A465_params.joint_4_to_wrist_dist;
  double shoulder_to_wrist = A465_params.shoulder_to_elbow_dist + elbow_to_wrist_dist;
  double s2e_dist_sqr = A465_params.shoulder_to_elbow_dist * A465_params.shoulder_to_elbow_dist;
  double e2w_dist_sqr = elbow_to_wrist_dist * elbow_to_wrist_dist;
  
  vect_n<double> solns[8];
  for(std::size_t i = 0; i < 8; ++i)
    solns[i].resize(7,0.0);
  
  /*
   * find the maximum wrist to base distance, x_max and verifies if the required
   * position of the end-effector is within limits
   */
  /*Extended arm*/
  double c2_max = cos(joint_upper_bounds[2]);
  double s2_max = sin(joint_upper_bounds[2]);
  
  double x_max = 0.0;
  if (wrist_pos[2] > (shoulder_to_wrist - extend_epsilon) * c2_max) {
    if (wrist_pos[1] * wrist_pos[1] + wrist_pos[2] * wrist_pos[2] > (shoulder_to_wrist - extend_epsilon) * (shoulder_to_wrist - extend_epsilon))
      throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is outside the cylindrical workspace envelope (fully-extended arm).");
    x_max = sqrt( shoulder_to_wrist * shoulder_to_wrist - wrist_pos[1] * wrist_pos[1] - wrist_pos[2] * wrist_pos[2]);
  } else /*Bent arm*/ {
    // Verifies that the location can be reached, as far as height (z) is concerned
    double max_j23_angle = joint_upper_bounds[2] + joint_upper_bounds[3];
    if(max_j23_angle > M_PI)
      max_j23_angle = M_PI;
    double low_elbow_to_desired_height = wrist_pos[2] - A465_params.shoulder_to_elbow_dist * c2_max;
    double elbow_wrist_eps_dist = elbow_to_wrist_dist - extend_epsilon;
    
    if (low_elbow_to_desired_height < elbow_wrist_eps_dist * cos(max_j23_angle))
      throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is too low for the manipulator to reach.");
    
    elbow_wrist_eps_dist *= elbow_wrist_eps_dist;
    low_elbow_to_desired_height *= low_elbow_to_desired_height;
    
    if ( low_elbow_to_desired_height + pow(fabs(wrist_pos[1]) - A465_params.shoulder_to_elbow_dist * s2_max, 2) > elbow_wrist_eps_dist)
      throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is too far below the manipulator.");
    
    double xy_dist = A465_params.shoulder_to_elbow_dist * s2_max + sqrt(elbow_wrist_eps_dist - low_elbow_to_desired_height);
    x_max = sqrt(xy_dist * xy_dist - wrist_pos[1] * wrist_pos[1]);
  };
  
  /* At this point, the range of joint values for the track is between +- x_max around wrist_pos[0] (to within track-limits). */
  
#if 0 
  /* NOTE: MP: This the original code from the A465 code added by Joel Robert. 
   * However, I am not convinced that this is good or even correct at all.
   * The difference between forward-reach and backward-reach does not imply
   * this arrangement of "in front of" end-effector, it simply means that the 
   * base is flipped around. I believe the best approach is to try and align
   * most of the robot in the direction it needs to point to (EE_z_axis).
   */
  /*Find the best stance possible*/
  if (*direction==FORWARD) {
    *ax = x_max;
    double x_base  = wrist_pos[0] - *ax;
    printf("CASE 1\n");
    if (x_base > joint_upper_bounds[0]) {
      retcode=REACH_TOO_FORWARD;
      printf("xbase retcode=%d x_base%f\n",retcode,x_base);
      return(retcode);
    };
    if (x_base < joint_lower_bounds[0] && wrist_pos[0] > joint_lower_bounds[0] + 0.01) {
      *ax= wrist_pos[0] - (joint_lower_bounds[0] + 0.01);
      printf("CASE 1.1  ax=%f \n",*ax);
    };
  } else { // When the direction is backward, the base will be "in front of" the end-effector
    *ax = -x_max;
    double x_base  = wrist_pos[0] - *ax;
    printf("CASE 2 x=%f  x_base=%f\n", wrist_pos[0], x_base);
    if (x_base < joint_lower_bounds[0]) {
      retcode=REACH_TOO_BACKWARD;
      printf("xbase retcode=%d x_base%f\n",retcode,x_base);
      return(retcode);
    };
    if (x_base > joint_upper_bounds[0] && wrist_pos[0] < joint_upper_bounds[0] - 0.01) {
      *ax=wrist_pos[0] - (joint_upper_bounds[0] - 0.01);
      printf("CASE 2.1  ax=%f \n",*ax);
    };
  };
#endif
  
  
  /* Joint 0 (track) */
  
  // first, check that the limits of the track permit at least some solution:
  if( wrist_pos[0] - x_max > joint_upper_bounds[0] )
    throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is too far ahead of the track (beyond upper track limit).");
  if( wrist_pos[0] + x_max < joint_lower_bounds[0] )
    throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is too far behind the track (beyond lower track limit).");
  
  /* solving for the best track position based on trying to get a right-angle at the elbow */ 
  double R0_rhs = s2e_dist_sqr + e2w_dist_sqr - wrist_pos[1] * wrist_pos[1] - wrist_pos[2] * wrist_pos[2];
  double x_desired;
  if(R0_rhs > 0) {
    // this means that it is possible to obtain an exact right angle at the elbow.
    if(EE_z_axis[0] >= 0.0)  // the EE needs to be pointing in the forward x direction:
      x_desired = wrist_pos[0] - sqrt(R0_rhs);
    else
      x_desired = wrist_pos[0] + sqrt(R0_rhs);
  } else {
    // this means that the best we can do is put the robot-base in x-alignment with the wrist position.
    x_desired = wrist_pos[0];
  };
  // clamp the x-solution to the track's range:
  if( x_desired < joint_lower_bounds[0] )
    x_desired = joint_lower_bounds[0];
  else if( x_desired > joint_upper_bounds[0] )
    x_desired = joint_upper_bounds[0];
  // update the wrist-position vector:
  wrist_pos[0] -= x_desired;
  
  /* reach forward solutions */
  solns[fun_posture][0] = x_desired;
  solns[fuf_posture][0] = x_desired;
  solns[fdn_posture][0] = x_desired;
  solns[fdf_posture][0] = x_desired;
  
  /* reach backward solutions */
  solns[bun_posture][0] = x_desired;
  solns[buf_posture][0] = x_desired;
  solns[bdn_posture][0] = x_desired;
  solns[bdf_posture][0] = x_desired;
  
  
  /* Joint 1 */
  
  if( (fabs(wrist_pos[0]) < pos_epsilon) && (fabs(wrist_pos[1]) < pos_epsilon) ) {
    /* we're in the joint 1 singularity directly above the origin */
    solns[fun_posture][1] = preferred_posture[1];
    solns[bun_posture][1] = clamp_to_pi_range(solns[fun_posture][1] + M_PI);
  } else {
    solns[fun_posture][1] = atan2(wrist_pos[1], wrist_pos[0]);
    solns[bun_posture][1] = clamp_to_pi_range(solns[fun_posture][1] + M_PI);
  };
  
  /* set up some variables for later */
  double s1_F = sin(solns[fun_posture][1]); /* forward */
  double c1_F = cos(solns[fun_posture][1]);
  
#if 0
  /* Joint 3 */ 
  
  double linkvar1 = A465_params.shoulder_to_elbow_dist * A465_params.shoulder_to_elbow_dist
                  + elbow_to_wrist_dist * elbow_to_wrist_dist;
  double linkvar2 = 4.0 * A465_params.shoulder_to_elbow_dist * A465_params.shoulder_to_elbow_dist
                  * elbow_to_wrist_dist * elbow_to_wrist_dist;
  
  double baseplane_dist = c1_F * wrist_pos[0] + s1_F * wrist_pos[1]; /* was f11p  */
  double wrist_dist_sqr = baseplane_dist * baseplane_dist + wrist_pos[2] * wrist_pos[2];
  
  double j3tmp0 = linkvar1 - wrist_dist_sqr;
  double j3tmp1 = linkvar2 - j3tmp0 * j3tmp0;
    
  /* ensure we're within reach */
  if( j3tmp1 < 0.0 )
    throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach!");
  
  /* now determine the Joint 3 angle */
  double j3tmp2 = atan2(-sqrt(j3tmp1), j3tmp0);
  solns[bun_posture][3] = clamp_to_pi_range(0.5 * M_PI + j3tmp2);
  solns[fun_posture][3] = clamp_to_pi_range(0.5 * M_PI - j3tmp2);
  
#else
  /* Joint 3 - modified version (MP: just seems like a more straight forward way to do it) */ 
  
  double baseplane_dist_sqr = wrist_pos[0] * wrist_pos[0] + wrist_pos[1] * wrist_pos[1];
  double wrist_dist_sqr = baseplane_dist_sqr + wrist_pos[2] * wrist_pos[2];
  double j3tmp0 = s2e_dist_sqr + e2w_dist_sqr - wrist_dist_sqr;
  double j3tmp1 = 4.0 * s2e_dist_sqr * e2w_dist_sqr - j3tmp0 * j3tmp0;
    
  /* ensure we're within reach */
  if( j3tmp1 < 0.0 )
    throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Cannot compute an elbow angle that would reach the desired wrist position.");
  
  /* now determine the Joint 3 angle (FUN) */
  solns[fun_posture][3] = -atan2(sqrt(j3tmp1), -j3tmp0);
  solns[bun_posture][3] = -solns[fun_posture][3];
#endif
  
  double inv_wrist_dist_sqr = 1.0 / wrist_dist_sqr;
  double baseplane_dist = sqrt(baseplane_dist_sqr);
  double s2e_projection = A465_params.shoulder_to_elbow_dist * sin(solns[fun_posture][3]);
  
  vect<double,2> to_wrist_along_l1(A465_params.shoulder_to_elbow_dist * cos(solns[fun_posture][3]) + elbow_to_wrist_dist, 
                                   s2e_projection);
  vect<double,2> to_wrist_global(-baseplane_dist, wrist_pos[2]); // reach-forward vector.
  
  
  /* register all solutions for joints 1 and 3 */
  
  /* reach forward solutions */
  solns[fuf_posture][1] = solns[fun_posture][1];
  solns[fdn_posture][1] = solns[fun_posture][1];
  solns[fdf_posture][1] = solns[fun_posture][1];
  
  /* reach backward solutions */
  solns[buf_posture][1] = solns[bun_posture][1];
  solns[bdn_posture][1] = solns[bun_posture][1];
  solns[bdf_posture][1] = solns[bun_posture][1];
  
  /* backward/up and forward/down solutions */
  solns[buf_posture][3] = solns[bun_posture][3];
  solns[fdn_posture][3] = solns[bun_posture][3];
  solns[fdf_posture][3] = solns[bun_posture][3];
  
  /* foward/up and backward/down solutions */
  solns[fuf_posture][3] = solns[fun_posture][3];
  solns[bdn_posture][3] = solns[fun_posture][3];
  solns[bdf_posture][3] = solns[fun_posture][3];
  
  /* We have two solutions for each of J1 and J3 (corresponding
   * to reach forward/back and elbow up/down) so there are now a total of 4
   * stance configurations that we can attain:  FU, BU, FD, BD.  We must find separate
   * solutions for each configuration for all the remaining joints.
   * To simplify this process we will loop through the following calculations 
   * 4 times, once for each of the configurations.
   */
  
  /* Prepare some variables that are constant throughout the iterations. */
  vect<double,2> EE_z_proj_F(c1_F * EE_z_axis[0] + s1_F * EE_z_axis[1], 
                            -s1_F * EE_z_axis[0] + c1_F * EE_z_axis[1]);
  vect<double,2> EE_y_proj_F(c1_F * EE_y_axis[0] + s1_F * EE_y_axis[1], 
                            -s1_F * EE_y_axis[0] + c1_F * EE_y_axis[1]);
  double j4tmp0_F    = -s1_F * EE_z_axis[0] + c1_F * EE_z_axis[1];
  double j4tmp1_F    =  c1_F * EE_z_axis[0] + s1_F * EE_z_axis[1];
  double j6tmp0_F    =  c1_F * EE_y_axis[0] + s1_F * EE_y_axis[1];
  double j6tmp1_F    = -s1_F * EE_y_axis[0] + c1_F * EE_y_axis[1];
  
  /* Because of the way the solution array is organized, we must
   * index through with a step size of two to hit each of the FU, BU
   * FD, and BD solutions.  See A465.h for more info    */
  for(unsigned int index = fun_posture; index <= 6; index += 2) {
    vect<double,3> EE_z_proj(EE_z_proj_F[0], EE_z_proj_F[1], EE_z_axis[2]);
    vect<double,3> EE_y_proj(EE_y_proj_F[0], EE_y_proj_F[1], EE_y_axis[2]);
    if(index & reach_backward) {
      EE_z_proj[0] = -EE_z_proj[0];
      EE_z_proj[1] = -EE_z_proj[1];
      EE_y_proj[0] = -EE_y_proj[0];
      EE_y_proj[1] = -EE_y_proj[1];
      to_wrist_global[0] = baseplane_dist;
      if(index & elbow_down)
        to_wrist_along_l1[1] =  s2e_projection;
      else
        to_wrist_along_l1[1] = -s2e_projection;
    } else {
      to_wrist_global[0] = -baseplane_dist;
      if(index & elbow_down)
        to_wrist_along_l1[1] = -s2e_projection;
      else
        to_wrist_along_l1[1] =  s2e_projection;
    };
    
    /* joint 2 */
    double c23 = (to_wrist_along_l1 % to_wrist_global) * inv_wrist_dist_sqr;
    double s23 = (to_wrist_along_l1 * to_wrist_global) * inv_wrist_dist_sqr;
    solns[index][2]   = atan2(s23, c23) - solns[index][3]; /* noflip solution */
    solns[index+1][2] = solns[index][2];                   /* flip solution */
    
    double c5 = EE_z_proj[0] * s23 + EE_z_proj[2] * c23;
    double s5 = EE_z_proj[0] * c23 - EE_z_proj[2] * s23;
    
    /* Joints 4, 5 & 6 */
    //double j4tmp2 =  c23 * j4tmp1 + s23 * EE_z_axis[2];
    
    /* first check if we're in the J5 singularity (i.e. J5 ~= 0)  */
    if( 1.0 - fabs(c5) < pos_epsilon ) {
      /* we're in the singularity */
      solns[index][5]   = 0.0;  /* noflip solution */
      solns[index+1][5] = 0.0;  /* flip solution */
      
      /* set joint 4 to it's preferred value since position is ambiguous. */
      solns[index][4]   = preferred_posture[4]; /* noflip solution */
      solns[index+1][4] = preferred_posture[4]; /* flip solution */
      double s4 = sin(preferred_posture[4]);
      double c4 = cos(preferred_posture[4]);
      
      //double s6 = -( c4 * ( c23 * j6tmp0 + s23 * EE_y_axis[2]) + s4 * j6tmp1 );
      //double c6 = -s4 * ( c23 * j6tmp0 + s23 * EE_y_axis[2] ) + c4 * j6tmp1;
      //solns[index][6]   = atan2(s6, c6);   /* wrist not flipped */
      //solns[index+1][6] = solns[index][6]; /* in singularity -- flip solution same as noflip solution */
    } else {
      /* we're not singular in jt 5 */
      //solns[index+1][4] = atan2(j4tmp0, j4tmp2);   /* wrist is flipped   */
      solns[index][4]   = clamp_to_pi_range(solns[index+1][4] + M_PI);   /* wrist is not flipped */
      double s4 = sin(solns[index][4]);
      double c4 = cos(solns[index][4]);
      
      solns[index][5]   =  atan2(s5,c5); /* noflip solution */
      solns[index+1][5] = -solns[index][5]; /* flip solution */
      
      //double s6 = -c5 * ( c4 * ( c23 * j6tmp0 + s23 * EE_y_axis[2] ) + s4 * j6tmp1 )
      //           + s5 * ( s23 * j6tmp0 - c23 * EE_y_axis[2] );
      //double c6 = -s4 * ( c23 * j6tmp0 + s23 * EE_y_axis[2] ) + c4 * j6tmp1;
      //solns[index][6]   = atan2(s6, c6);   /* wrist not flipped */
      //solns[index+1][6] = clamp_to_pi_range(solns[index][6] + M_PI);  /* wrist flipped */
    };
  }; /* end of Joints 2-6 for loop*/
  
  return solns[0];
};

vect_n<double> CRS_A465_model_builder::compute_inverse_kinematics(const frame_3D<double>& EE_state) const {
  return vect_n<double>();
};




#if 0

// Below are some code retrieved from the CRS A465 robot software to compute the closed-form inverse kinematics:



/* A465_inverse_kinematics:     Computes the kinematic joint angles 'jointsoln' from 
 *      the basic arm transform 'transform'.
 *
 *      Note:   'old_joints' is an array representing the robots current kinematic
 *                  joint angles.  It is used when a singularity is encountered in joint
 *          1 (wrist directly above origin) or joint 5 (when J5 is exactly 0).
 *          In each case, we appeal to the previous joint angle for a reasonable
 *          approximation of where the next joint angle should be.  Without an 
 *          approach like this the arm may move erratically and dangerously.
 */
static int  A465_inverse_kinematics(double transform[XFORM_SIZE], Joint_Solution jointsoln[NUM_SOLN])
{
        double old_joints[NUM_AXES]={0,0,0,0,0,0};
    double linkvar1, linkvar2;
    double j1soln1, j1soln2;
    double s1, s1_F, s1_B, c1, c1_F, c1_B;   
    double j3tmp1, j3tmp2, j3tmp3, j3tmp4;
    double j3soln1, j3soln2;
    double s3, c3_BUFD, c3_FUBD;
    double threespace_dist_sqr;
    double inv_3space_dist_sqr;
    double baseplane_dist;
    double j2tmp1_BUFD, j2tmp1_FUBD;
    double j2tmp1, j2tmp2;
    double s23, c23;
    double basedist_26;
    double j2soln;
    double j5_singularity;
    double s5, c5, s4, c4, s6, c6;
    double j4tmp1, j4tmp2, j4tmp3;
    double j4soln1, j4soln2;
    double j5soln1, j5soln2;
    double j6tmp1, j6tmp2;
    double j6soln1, j6soln2;
    int index;

                                                                                                
    /* Joint 1 */
#ifndef TESTING 
    if( fabs(transform[Px]) < EPSILON && fabs(transform[Py]) < EPSILON )
    {
        /* we're in the joint 1 singularity directly above the origin */
        j1soln1 = old_joints[J1];
            j1soln2 = j1soln1 + PI;
    }
    else
#endif
    {
            j1soln1 = atan2(transform[Py],transform[Px]);
            j1soln2 = j1soln1 + PI;
    }

    /* reach forward solutions */
    jointsoln[FUN].joint[J1] = j1soln1;
    jointsoln[FUF].joint[J1] = j1soln1;
    jointsoln[FDN].joint[J1] = j1soln1;
    jointsoln[FDF].joint[J1] = j1soln1;

    /* reach backward solutions */
    jointsoln[BUN].joint[J1] = j1soln2;
    jointsoln[BUF].joint[J1] = j1soln2;
    jointsoln[BDN].joint[J1] = j1soln2;
    jointsoln[BDF].joint[J1] = j1soln2;

    /* set up some variables for later */

    s1_F=sin(j1soln1); /* forward */
    c1_F=cos(j1soln1);
    s1_B = -s1_F;       /* backward */
    c1_B = -c1_F;

    
    
    /* Joint 3 */
        
    linkvar1    = A2*A2+D4*D4;
    linkvar2    = 4.0*A2*A2*D4*D4;

    baseplane_dist = c1_F*transform[Px] + s1_F*transform[Py]; /* was f11p  */
    threespace_dist_sqr = baseplane_dist * baseplane_dist + transform[Pz] * transform[Pz];
    j3tmp1 = linkvar1 - threespace_dist_sqr;
    j3tmp2 = linkvar2 - j3tmp1*j3tmp1 ;
    
    /* ensure we're within reach */
        if( j3tmp2 < -linkvar2*EPSILON )/* EPSILON defined in header file */
        {  
                //printf("INVKIN_OUTOFREACH\n");
        return( INVKIN_OUTOFREACH );
        }
    else if( j3tmp2 < 0.0 )
                j3tmp2 = 0.0;
    
    /* now determine the Joint 3 angle */
    j3tmp3 = sqrt( j3tmp2 );
    j3tmp4 = atan2( -j3tmp3,j3tmp1);
    j3soln1 = PID2 + j3tmp4;
    j3soln2 = PID2 - j3tmp4;

    /* backward/up and forward/down solutions */
    jointsoln[BUN].joint[J3] = j3soln1;
    jointsoln[BUF].joint[J3] = j3soln1;
    jointsoln[FDN].joint[J3] = j3soln1;
    jointsoln[FDF].joint[J3] = j3soln1;

    /* foward/up and backward/down solutions */
    jointsoln[FUN].joint[J3] = j3soln2;
    jointsoln[FUF].joint[J3] = j3soln2;
    jointsoln[BDN].joint[J3] = j3soln2;
    jointsoln[BDF].joint[J3] = j3soln2;

    
    /* now prepare some variables for further stages. */
    s3      =sin(j3soln1);  /* was tmp2 */
    c3_BUFD =cos(j3soln1);  /* was tmp1 */
    c3_FUBD = -c3_BUFD;     /* same as cos(j3soln2) */

    /* pre-calculating the inverse here saves cpu cycles later on */
    inv_3space_dist_sqr = 1/threespace_dist_sqr;


    /* Joint 2 --> Joint 6 */

    /* set up joint 2 variables */
    j2tmp1_BUFD = A2*c3_BUFD;
    j2tmp1_FUBD = A2*c3_FUBD;

    j2tmp2 = D4 - A2*s3;

  
    /* We have two solutions for each of J1 and J3 (corresponding
     * to reach forward/back and elbow up/down) so there are now a total of 4
     * stance configurations that we can attain:  FU, BU, FD, BD.  We must find separate
     * solutions for each configuration for all the remaining joints.
     * To simplify this process we will loop through the following calculations 
     * 4 times, once for each of the configurations.
    */


    /* Because of the way the solution array is organized, we must
     * index through with a step size of two to hit each of the FU, BU
     * FD, and BD solutions.  See A465.h for more info    */
    for (index=FUN;index<=6;index+=2)
    {
        switch(index)
        {
            case(FUN):
            {
                c1          = c1_F;
                s1          = s1_F;
                basedist_26 = baseplane_dist;
                j2tmp1      = j2tmp1_FUBD;
                break;
            }
            case(FDN):
            {
                c1          = c1_F;
                s1          = s1_F;
                basedist_26 = baseplane_dist;
                j2tmp1      = j2tmp1_BUFD;  
                break;
            }
            case(BUN):
            {
                c1          = c1_B;
                s1          = s1_B;
                basedist_26 = -baseplane_dist;
                j2tmp1      = j2tmp1_BUFD;
                break;
            }
            case(BDN):
            {
                c1          = c1_B;
                s1          = s1_B;
                basedist_26 = -baseplane_dist;
                j2tmp1      = j2tmp1_FUBD;
                break;
            }
                        default:
                                {
                                //printf("CASE_ERROR\n");
                return CASE_ERROR;
                                }
        }
        
        
        /* joint 2 */
        
        s23 = (j2tmp1*transform[Pz] - j2tmp2*basedist_26) * inv_3space_dist_sqr;
        c23 = ((j2tmp2*transform[Pz]) + (j2tmp1*basedist_26))* inv_3space_dist_sqr; 
        j2soln = atan2(s23, c23) - jointsoln[index].joint[J3];
    
        jointsoln[index].joint[J2]      = j2soln;   /* noflip solution */
        jointsoln[index+1].joint[J2]    = j2soln;   /* flip solution */


        /* Joint 4 */


        j4tmp1 = -s1*transform[r13] + c1*transform[r23];
        j4tmp2 =  c1*transform[r13] + s1*transform[r23];
        j4tmp3 =  c23*j4tmp2 + s23*transform[r33];

    
        /* first check if we're in the J5 singularity (i.e. J5 ~= 0)  */        
#ifndef TESTING
        if( (fabs(j4tmp1) < EPSILON) && (fabs(j4tmp3) < EPSILON) )
            {
            /* we're in the singularity -- set flag and set J5 to exactly 0 */
            j5_singularity = 1;
                jointsoln[index].joint[J5]      = 0.0;  /* noflip solution */
            jointsoln[index+1].joint[J5]    = 0.0;  /* flip solution */
            
            s5 = 0.0;
                c5 = 1.0;
                
            /* set joint 4 to it's old value since we have no way of knowing
                where else it should be (position is ambiguous)  */
            jointsoln[index].joint[J4] = old_joints[J4];   /* noflip solution */
            jointsoln[index+1].joint[J4] = old_joints[J4]; /* flip solution */
            
            s4 = sin(old_joints[J4]);
            c4 = cos(old_joints[J4]);


        }
            else
#endif
        {
            /* we're not singular in jt 5 */
            j5_singularity = 0;
    
                j4soln1 = atan2(j4tmp1,j4tmp3);
            j4soln2 = j4soln1 + PI;

                jointsoln[index+1].joint[J4]    = j4soln1;   /* wrist is flipped   */
                jointsoln[index].joint[J4]      = j4soln2;   /* wrist is not flipped */
            
            s4 = sin(j4soln2);
            c4 = cos(j4soln2);
    
            }


        /* Joint 5 */

        /* now, if we're not in a singularity, we must still solve for J5 */
            if( !j5_singularity )
        {
                s5 = -c4*j4tmp3 - s4*j4tmp1;
                c5 = -s23*j4tmp2 + c23*transform[r33];
            j5soln1 = atan2(s5,c5);
            j5soln2 = -j5soln1;
            jointsoln[index].joint[J5]   = j5soln1; /* noflip solution */
                jointsoln[index+1].joint[J5] = j5soln2; /* flip solution */
        }


        /* Joint 6 */
        
            j6tmp1 = c1*transform[r12] + s1*transform[r22];
        j6tmp2 = -s1*transform[r12] + c1*transform[r22];

            s6   = -c5*(c4*(c23*j6tmp1+s23*transform[r32]) + s4*j6tmp2)
                + s5*(s23*j6tmp1-c23*transform[r32]);
    
        c6   = -s4*(c23*j6tmp1 + s23*transform[r32]) + c4*j6tmp2;

        j6soln1 = atan2(s6, c6);
        j6soln2 = j6soln1 + PI;
    
        jointsoln[index].joint[J6]   = j6soln1;  /* wrist not flipped */

            if( j5_singularity ) 
            /* in singularity -- flip solution same as noflip solution */
            jointsoln[index+1].joint[J6] = j6soln1; 
        else
            jointsoln[index+1].joint[J6] = j6soln2;  /* wrist flipped */        



    } /* end of Joints 2-6 for loop*/

    return(OK);

}   /* end of A465_inverse_kinematics() */



static int A465_find_initial_base_distance(double *ax, double arm_without_terminal[XFORM_SIZE],int *direction)
{
        int retcode=0;
        double x,y,z,l,x_max,xy_dist,x_base,epsilon=.06;
        l=A2+D4;

        /*Obtain the transformation of the wrist*/
        x=arm_without_terminal[Px];
        y=arm_without_terminal[Py];
        z=arm_without_terminal[Pz];

        
        /*
        find the maximum wrist to base distance, x_max and verifies if the required
        position of the end-effector is within limits
        */
        if (z>D1+(l-epsilon)*cos(A465_posjointlim[J2])) /*Extended arm*/
        {
                //printf("Extended Arm\n");
                if (y*y+(z-D1)*(z-D1)>(l-epsilon)*(l-epsilon))
                {
                        retcode=REACH_TOO_FAR;
                        return(retcode);
                }
                x_max=sqrt(l*l-y*y-(z-D1)*(z-D1));
        }
        else /*Bent arm*/
        {
                //printf("Bent Arm\n");
                // Verifies that the location can be reached, as far as height (z) is concerned
                if (z<D1+A2*cos(A465_posjointlim[J2])-(D4-epsilon)*cos(PI-A465_posjointlim[J2]-A465_posjointlim[J3]))
                {
                        retcode=REACH_TOO_LOW;
                        return(retcode);
                }
                else
                {
                        if (pow(D1+A2*cos(A465_posjointlim[J2])-z,2)+pow(fabs(y)-A2*sin(A465_posjointlim[J2]),2)>(D4-epsilon)*(D4-epsilon))
                        {
                                retcode=REACH_TOO_FAR;
                                return(retcode);
                        }
                }
                xy_dist=A2*sin(A465_posjointlim[J2])+sqrt((D4-epsilon)*(D4-epsilon)-pow(D1+A2*cos(A465_posjointlim[J2])-z,2));
                x_max=sqrt(xy_dist*xy_dist-y*y);
        }

        if (x_max < -EPSILON)
        {
                retcode=OUT_OF_REACH;
                return(retcode);
        }
        
        /*Find the best stance possible*/
        if (*direction==FORWARD){
                *ax             = x_max;
                x_base  = x-*ax;
                printf("CASE 1\n");
                if (x_base>A465_posjointlim[J7])
                {
                        retcode=REACH_TOO_FORWARD;
                        printf("xbase retcode=%d x_base%f\n",retcode,x_base);
                        return(retcode);
                }
                if (x_base<A465_negjointlim[J7] && x>A465_negjointlim[J7]+0.01)
                {
                        *ax=x-(A465_negjointlim[J7]+0.01);
                        printf("CASE 1.1  ax=%f \n",*ax);
                }
        }
        else  // When the direction is backward, the base will be "in front of" the end-effector
        {
                *ax             = -x_max;
                x_base  = x-*ax;
                printf("CASE 2 x=%f  x_base=%f\n",x,x_base);
                if (x_base<A465_negjointlim[J7])
                {
                        retcode=REACH_TOO_BACKWARD;
                        printf("xbase retcode=%d x_base%f\n",retcode,x_base);
                        return(retcode);
                }
                if (x_base>A465_posjointlim[J7] && x<A465_posjointlim[J7]-0.01)
                {
                        *ax=x-(A465_posjointlim[J7]-0.01);
                        printf("CASE 2.1  ax=%f \n",*ax);
                }
        }
        return(OK);
}





#endif




};
  
  
};








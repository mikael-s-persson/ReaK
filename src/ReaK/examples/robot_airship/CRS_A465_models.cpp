
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

#include "mbd_kte/inertia.hpp"
#include "mbd_kte/spring.hpp"
#include "mbd_kte/damper.hpp"
#include "mbd_kte/rigid_link.hpp"
#include "mbd_kte/revolute_joint.hpp"
#include "mbd_kte/prismatic_joint.hpp"
#include "mbd_kte/kte_map_chain.hpp"
#include "mbd_kte/force_actuator.hpp"
#include "mbd_kte/joint_friction.hpp"
#include "mbd_kte/driving_actuator.hpp"
#include "mbd_kte/mass_matrix_calculator.hpp"
#include "mbd_kte/jacobian_joint_map.hpp"
#include "kte_models/manip_dynamics_model.hpp"

#include "serialization/xml_archiver.hpp"

#include "rtti/typed_primitives.hpp"

#include "optimization/optim_exceptions.hpp"
#include "lin_alg/mat_qr_decomp.hpp"

#include <utility>
#include <vector>

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
  joint_upper_bounds[4] = 3.14159265359;
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
  return pp::make_Ndof_rl_space<7>(joint_lower_bounds, joint_upper_bounds, 
                                   joint_rate_limits.gen_speed_limits, 
                                   joint_rate_limits.gen_accel_limits, 
                                   joint_rate_limits.gen_jerk_limits);
//   return joint_rate_limits.make_rl_joint_space(get_joint_space());
};

CRS_A465_model_builder::joint_space_type CRS_A465_model_builder::get_joint_space() const {
  return pp::make_Ndof_space<7>(joint_lower_bounds, joint_upper_bounds, 
                                joint_rate_limits.gen_speed_limits, 
                                joint_rate_limits.gen_accel_limits);
  /*
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
	                 );*/
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
	  pp::hyperball_topology< vect<double,3> >("EE_vel_space",vect<double,3>(0.0,0.0,0.0),5.0),
	  pp::hyperball_topology< vect<double,3> >("EE_acc_space",vect<double,3>(0.0,0.0,0.0),25.0)
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



CRS_A465_model_builder::rate_limited_joint_space_1st_type CRS_A465_model_builder::get_rl_joint_space_1st() const {
  return pp::make_Ndof_rl_space<7>(joint_lower_bounds, joint_upper_bounds, 
                                   joint_rate_limits.gen_speed_limits, 
                                   joint_rate_limits.gen_accel_limits);
//   return joint_rate_limits.make_rl_joint_space(get_joint_space_1st());
};

CRS_A465_model_builder::joint_space_1st_type CRS_A465_model_builder::get_joint_space_1st() const {
  return pp::make_Ndof_space<7>(joint_lower_bounds, joint_upper_bounds, 
                                joint_rate_limits.gen_speed_limits);
  /*typedef pp::joint_space_1st_order<double>::type SingleJointSpace;
  typedef pp::line_segment_topology<double> LinSeg;
  return joint_space_1st_type( arithmetic_tuple< 
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
                                 LinSeg,LinSeg
                               >(
                                 LinSeg("track_pos_space",joint_lower_bounds[0],joint_upper_bounds[0]),
                                 LinSeg("track_vel_space",-joint_rate_limits.gen_speed_limits[0],joint_rate_limits.gen_speed_limits[0])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg,LinSeg
                               >(
                                 LinSeg("arm_joint_1_pos_space",joint_lower_bounds[1],joint_upper_bounds[1]),
                                 LinSeg("arm_joint_1_vel_space",-joint_rate_limits.gen_speed_limits[1],joint_rate_limits.gen_speed_limits[1])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg,LinSeg
                               >(
                                 LinSeg("arm_joint_2_pos_space",joint_lower_bounds[2],joint_upper_bounds[2]),
                                 LinSeg("arm_joint_2_vel_space",-joint_rate_limits.gen_speed_limits[2],joint_rate_limits.gen_speed_limits[2])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg,LinSeg
                               >(
                                 LinSeg("arm_joint_3_pos_space",joint_lower_bounds[3],joint_upper_bounds[3]),
                                 LinSeg("arm_joint_3_vel_space",-joint_rate_limits.gen_speed_limits[3],joint_rate_limits.gen_speed_limits[3])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg,LinSeg
                               >(
                                 LinSeg("arm_joint_4_pos_space",joint_lower_bounds[4],joint_upper_bounds[4]),
                                 LinSeg("arm_joint_4_vel_space",-joint_rate_limits.gen_speed_limits[4],joint_rate_limits.gen_speed_limits[4])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg,LinSeg
                               >(
                                 LinSeg("arm_joint_5_pos_space",joint_lower_bounds[5],joint_upper_bounds[5]),
                                 LinSeg("arm_joint_5_vel_space",-joint_rate_limits.gen_speed_limits[5],joint_rate_limits.gen_speed_limits[5])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg,LinSeg
                               >(
                                 LinSeg("arm_joint_6_pos_space",joint_lower_bounds[6],joint_upper_bounds[6]),
                                 LinSeg("arm_joint_6_vel_space",-joint_rate_limits.gen_speed_limits[6],joint_rate_limits.gen_speed_limits[6])
                               )
                             )
                           )
                         );*/
};

CRS_A465_model_builder::end_effector_space_1st_type CRS_A465_model_builder::get_end_effector_space_1st() const {
  
  typedef arithmetic_tuple< pp::hyperbox_topology< vect<double,3> >, 
                            pp::hyperball_topology< vect<double,3> >
                          > TranslationTuple;
  typedef arithmetic_tuple< pp::quaternion_topology<double>, 
                            pp::ang_velocity_3D_topology<double>
                          > RotationTuple;
  
  typedef arithmetic_tuple_element<0, end_effector_space_1st_type>::type TranslationDiffSpace;
  typedef arithmetic_tuple_element<1, end_effector_space_1st_type>::type RotationDiffSpace;
  
  typedef arithmetic_tuple<TranslationDiffSpace, RotationDiffSpace> EESpaceTuple;
  
  return end_effector_space_1st_type(
    EESpaceTuple(
      TranslationDiffSpace(
        TranslationTuple(
          pp::hyperbox_topology< vect<double,3> >("EE_pos_space",vect<double,3>(-1.1,-1.1,0.0),vect<double,3>(4.5,1.1,1.5)),
          pp::hyperball_topology< vect<double,3> >("EE_vel_space",vect<double,3>(0.0,0.0,0.0),5.0)
        )
      ),
      RotationDiffSpace(
        RotationTuple(
          pp::quaternion_topology<double>("EE_rotation_space"),
          pp::ang_velocity_3D_topology<double>("EE_ang_vel_space",10.0)
        )
      )
    )
  );
};




CRS_A465_model_builder::rate_limited_joint_space_0th_type CRS_A465_model_builder::get_rl_joint_space_0th() const {
  return pp::make_Ndof_rl_space<7>(joint_lower_bounds, joint_upper_bounds, 
                                   joint_rate_limits.gen_speed_limits);
//   return joint_rate_limits.make_rl_joint_space(get_joint_space_0th());
};

CRS_A465_model_builder::joint_space_0th_type CRS_A465_model_builder::get_joint_space_0th() const {
  return pp::make_Ndof_space<7>(joint_lower_bounds, joint_upper_bounds);
  /*typedef pp::joint_space_0th_order<double>::type SingleJointSpace;
  typedef pp::line_segment_topology<double> LinSeg;
  return joint_space_0th_type( arithmetic_tuple< 
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
                                 LinSeg
                               >(
                                 LinSeg("track_pos_space",joint_lower_bounds[0],joint_upper_bounds[0])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg
                               >(
                                 LinSeg("arm_joint_1_pos_space",joint_lower_bounds[1],joint_upper_bounds[1])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg
                               >(
                                 LinSeg("arm_joint_2_pos_space",joint_lower_bounds[2],joint_upper_bounds[2])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg
                               >(
                                 LinSeg("arm_joint_3_pos_space",joint_lower_bounds[3],joint_upper_bounds[3])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg
                               >(
                                 LinSeg("arm_joint_4_pos_space",joint_lower_bounds[4],joint_upper_bounds[4])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg
                               >(
                                 LinSeg("arm_joint_5_pos_space",joint_lower_bounds[5],joint_upper_bounds[5])
                               )
                             ),
                             SingleJointSpace(
                               arithmetic_tuple<
                                 LinSeg
                               >(
                                 LinSeg("arm_joint_6_pos_space",joint_lower_bounds[6],joint_upper_bounds[6])
                               )
                             )
                           )
                         );*/
};

CRS_A465_model_builder::end_effector_space_0th_type CRS_A465_model_builder::get_end_effector_space_0th() const {
  
  typedef arithmetic_tuple< pp::hyperbox_topology< vect<double,3> > > TranslationTuple;
  typedef arithmetic_tuple< pp::quaternion_topology<double> > RotationTuple;
  
  typedef arithmetic_tuple_element<0, end_effector_space_0th_type>::type TranslationDiffSpace;
  typedef arithmetic_tuple_element<1, end_effector_space_0th_type>::type RotationDiffSpace;
  
  typedef arithmetic_tuple<TranslationDiffSpace, RotationDiffSpace> EESpaceTuple;
  
  return end_effector_space_0th_type(
    EESpaceTuple(
      TranslationDiffSpace(
        TranslationTuple(
          pp::hyperbox_topology< vect<double,3> >("EE_pos_space",vect<double,3>(-1.1,-1.1,0.0),vect<double,3>(4.5,1.1,1.5))
        )
      ),
      RotationDiffSpace(
        RotationTuple(
          pp::quaternion_topology<double>("EE_rotation_space")
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

frame_3D<double> CRS_A465_model_builder::compute_direct_kinematics_with_vel(const vect_n<double>& joint_states, mat<double,mat_structure::rectangular>* jacobian) const {
  
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
  flange.Quat = robot_base->Quat;
  vect<double,3> e0 = flange.Quat * vect_i;
  flange.Position += joint_states[0] * e0;
  vect<double,3> e1 = flange.Quat * vect_k;
  flange.Position += A465_params.baseplate_to_shoulder_dist * e1;
  vect<double,3> pb = flange.Position;
  flange.Quat *= q1;
  vect<double,3> e2 = -(flange.Quat * vect_j);
  flange.Quat *= q2;
  vect<double,3> a2 = A465_params.shoulder_to_elbow_dist * (flange.Quat * vect_k);
  flange.Position += a2;
  flange.Quat *= q3;
  vect<double,3> a3 = (A465_params.elbow_to_joint_4_dist + A465_params.joint_4_to_wrist_dist) * (flange.Quat * vect_k);
  flange.Position += a3;
  vect<double,3> e4 =  flange.Quat * vect_k;
  flange.Quat *= q4;
  vect<double,3> e5 = -(flange.Quat * vect_j);
  flange.Quat *= q5;
  vect<double,3> e6 =  flange.Quat * vect_k;
  flange.Quat *= q6;
  vect<double,3> a6 = A465_params.wrist_to_flange_dist * e6; 
  flange.Position += a6;
  
  vect<double,3> s2f = flange.Position - pb; // shoulder to flange
  flange.Velocity  = joint_states[7] * e0;
  vect<double,3> w1 = joint_states[8]  * e1;
  flange.AngVelocity  = w1;
  vect<double,3> w2 = joint_states[9]  * e2;
  flange.AngVelocity += w2;
  vect<double,3> w3 = joint_states[10] * e2;
  flange.AngVelocity += w3;
  flange.Velocity += (w1 + w2 + w3) % s2f;
  flange.Velocity += -w3 % a2;
  vect<double,3> w4 = joint_states[11] * e4;
  flange.AngVelocity += w4;
  vect<double,3> w5 = joint_states[12] * e5;
  flange.AngVelocity += w5;
  flange.Velocity += (w4 + w5) % a6;
  flange.AngVelocity += joint_states[13] * e6;
  
  if(jacobian) {
    jacobian->resize(std::make_pair(6,7));
    (*jacobian)(0,0) = e0[0];
    (*jacobian)(1,0) = e0[1];
    (*jacobian)(2,0) = e0[2];
    (*jacobian)(3,0) = 0.0;
    (*jacobian)(4,0) = 0.0;
    (*jacobian)(5,0) = 0.0;
    vect<double,3> v1 = e1 % s2f;
    (*jacobian)(0,1) = v1[0];
    (*jacobian)(1,1) = v1[1];
    (*jacobian)(2,1) = v1[2];
    (*jacobian)(3,1) = e1[0];
    (*jacobian)(4,1) = e1[1];
    (*jacobian)(5,1) = e1[2];
    vect<double,3> v2 = e2 % s2f;
    (*jacobian)(0,2) = v2[0];
    (*jacobian)(1,2) = v2[1];
    (*jacobian)(2,2) = v2[2];
    (*jacobian)(3,2) = e2[0];
    (*jacobian)(4,2) = e2[1];
    (*jacobian)(5,2) = e2[2];
    vect<double,3> v3 = e2 % (s2f - a2);
    (*jacobian)(0,3) = v3[0];
    (*jacobian)(1,3) = v3[1];
    (*jacobian)(2,3) = v3[2];
    (*jacobian)(3,3) = e2[0];
    (*jacobian)(4,3) = e2[1];
    (*jacobian)(5,3) = e2[2];
    vect<double,3> v4 = e4 % a6;
    (*jacobian)(0,4) = v4[0];
    (*jacobian)(1,4) = v4[1];
    (*jacobian)(2,4) = v4[2];
    (*jacobian)(3,4) = e4[0];
    (*jacobian)(4,4) = e4[1];
    (*jacobian)(5,4) = e4[2];
    vect<double,3> v5 = e5 % a6;
    (*jacobian)(0,5) = v5[0];
    (*jacobian)(1,5) = v5[1];
    (*jacobian)(2,5) = v5[2];
    (*jacobian)(3,5) = e5[0];
    (*jacobian)(4,5) = e5[1];
    (*jacobian)(5,5) = e5[2];
    (*jacobian)(0,6) = 0.0;
    (*jacobian)(1,6) = 0.0;
    (*jacobian)(2,6) = 0.0;
    (*jacobian)(3,6) = e6[0];
    (*jacobian)(4,6) = e6[1];
    (*jacobian)(5,6) = e6[2];
  };
  
  return flange;
};


void CRS_A465_model_builder::compute_jacobian_matrix(const vect_n<double>& joint_positions, mat<double,mat_structure::rectangular>& jacobian) const {
  
  /* calculate individual rotations */
  quaternion<double>::zrot q1( joint_positions[1]);
  quaternion<double>::yrot q2(-joint_positions[2]);
  quaternion<double>::yrot q3(-joint_positions[3]);
  quaternion<double>::zrot q4( joint_positions[4]);
  quaternion<double>::yrot q5(-joint_positions[5]);
  quaternion<double>::zrot q6( joint_positions[6]);
  
  quaternion<double> q_accum = robot_base->Quat;
  vect<double,3> e0 = q_accum * vect_i;
  vect<double,3> e1 = q_accum * vect_k;
  vect<double,3> pb = robot_base->Position + joint_positions[0] * e0 + A465_params.baseplate_to_shoulder_dist * e1;
  q_accum *= q1;
  vect<double,3> e2 = -(q_accum * vect_j);
  q_accum *= q2;
  vect<double,3> a2 = A465_params.shoulder_to_elbow_dist * (q_accum * vect_k);
  q_accum *= q3;
  vect<double,3> a3 = (A465_params.elbow_to_joint_4_dist + A465_params.joint_4_to_wrist_dist) * (q_accum * vect_k);
  vect<double,3> e4 =  q_accum * vect_k;
  q_accum *= q4;
  vect<double,3> e5 = -(q_accum * vect_j);
  q_accum *= q5;
  vect<double,3> e6 =  q_accum * vect_k;
  q_accum *= q6;
  vect<double,3> a6 = A465_params.wrist_to_flange_dist * e6; 
  
  vect<double,3> s2f = a2 + a3 + a6; // shoulder to flange
  
  jacobian.resize(std::make_pair(6,7));
  jacobian(0,0) = e0[0];
  jacobian(1,0) = e0[1];
  jacobian(2,0) = e0[2];
  jacobian(3,0) = 0.0;
  jacobian(4,0) = 0.0;
  jacobian(5,0) = 0.0;
  vect<double,3> v1 = e1 % s2f;
  jacobian(0,1) = v1[0];
  jacobian(1,1) = v1[1];
  jacobian(2,1) = v1[2];
  jacobian(3,1) = e1[0];
  jacobian(4,1) = e1[1];
  jacobian(5,1) = e1[2];
  vect<double,3> v2 = e2 % s2f;
  jacobian(0,2) = v2[0];
  jacobian(1,2) = v2[1];
  jacobian(2,2) = v2[2];
  jacobian(3,2) = e2[0];
  jacobian(4,2) = e2[1];
  jacobian(5,2) = e2[2];
  vect<double,3> v3 = e2 % (s2f - a2);
  jacobian(0,3) = v3[0];
  jacobian(1,3) = v3[1];
  jacobian(2,3) = v3[2];
  jacobian(3,3) = e2[0];
  jacobian(4,3) = e2[1];
  jacobian(5,3) = e2[2];
  vect<double,3> v4 = e4 % a6;
  jacobian(0,4) = v4[0];
  jacobian(1,4) = v4[1];
  jacobian(2,4) = v4[2];
  jacobian(3,4) = e4[0];
  jacobian(4,4) = e4[1];
  jacobian(5,4) = e4[2];
  vect<double,3> v5 = e5 % a6;
  jacobian(0,5) = v5[0];
  jacobian(1,5) = v5[1];
  jacobian(2,5) = v5[2];
  jacobian(3,5) = e5[0];
  jacobian(4,5) = e5[1];
  jacobian(5,5) = e5[2];
  jacobian(0,6) = 0.0;
  jacobian(1,6) = 0.0;
  jacobian(2,6) = 0.0;
  jacobian(3,6) = e6[0];
  jacobian(4,6) = e6[1];
  jacobian(5,6) = e6[2];
};


static double clamp_to_pi_range(double a) {
  return (a > M_PI ? a - 2.0 * M_PI : (a < -M_PI ? a + 2.0 * M_PI : a) );
};


vect_n<double> CRS_A465_model_builder::compute_inverse_kinematics(const pose_3D<double>& EE_pose) const {
  using std::sin; using std::cos; using std::fabs; 
  using std::atan2; using std::sqrt; using std::pow;
  
  const double pos_epsilon = 1e-3;
  const double extend_epsilon = .02;
  
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
  
  /* Joint 0 (track) */
  
  // first, check that the limits of the track permit at least some solution:
  if( wrist_pos[0] - x_max > joint_upper_bounds[0] )
    throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is too far ahead of the track (beyond upper track limit).");
  if( wrist_pos[0] + x_max < joint_lower_bounds[0] )
    throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is too far behind the track (beyond lower track limit).");
  
  double x_desired;
#if 0
  /* NOTE: MP: This the original code from the A465 code added by Joel Robert. 
   * However, I am not convinced that this is good or even correct at all.
   * The difference between forward-reach and backward-reach does not imply
   * this arrangement of "in front of" end-effector, it simply means that the 
   * base is flipped around. I believe the best approach is to try and align
   * most of the robot in the direction it needs to point to (EE_z_axis).
   */
  /*Find the best stance possible*/
  if( EE_z_axis[0] > 0.0 ) {
    x_desired = wrist_pos[0] - x_max;
  } else { // When the direction is backward, the base will be "in front of" the end-effector
    x_desired = wrist_pos[0] + x_max;
  };
#else
  /* solving for the best track position based on trying to get a right-angle at the elbow */ 
  double R0_rhs = s2e_dist_sqr + e2w_dist_sqr - wrist_pos[1] * wrist_pos[1] - wrist_pos[2] * wrist_pos[2];
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
#endif
  
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
  
  double inv_wrist_dist_sqr = 1.0 / wrist_dist_sqr;
  double baseplane_dist = sqrt(baseplane_dist_sqr);
  double s2e_projection = A465_params.shoulder_to_elbow_dist * sin(solns[fun_posture][3]);
  
  vect<double,2> to_wrist_along_l1(A465_params.shoulder_to_elbow_dist * cos(solns[fun_posture][3]) + elbow_to_wrist_dist, 
                                   s2e_projection);
  vect<double,2> to_wrist_global(baseplane_dist, wrist_pos[2]); // reach-forward vector.
  
  
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
  
  // solve for both elbow-cases (up or down).
  for(unsigned int i = 0; i < 2; ++i) {
    unsigned int e_o = 2 * i;
    double s2e_projection_i = (i ? -s2e_projection : s2e_projection);
    double c23 = (to_wrist_along_l1[0] * wrist_pos[2] + s2e_projection_i * baseplane_dist) * inv_wrist_dist_sqr;
    double s23 = /*reach_sign * */( s2e_projection_i * wrist_pos[2] - to_wrist_along_l1[0] * baseplane_dist) * inv_wrist_dist_sqr;
    
    solns[fuf_posture + e_o][2] = (solns[fun_posture + e_o][2] = atan2(s23, c23) - solns[fun_posture + e_o][3]);
    solns[buf_posture + e_o][2] = (solns[bun_posture + e_o][2] = atan2(-s23, c23) - solns[bun_posture + e_o][3]);
    
    double jt5_i = s23 * EE_z_axis[2] + c23 * EE_z_proj_F[0];
    
    /* first check if we're in the J5 singularity (i.e. J5 ~= 0)  */
    if( (fabs(jt5_i) < pos_epsilon) && (fabs(EE_z_proj_F[1]) < pos_epsilon) ) {
      double c4 = /*reach_sign * */cos(preferred_posture[4]);
      double s4 = /*reach_sign * */sin(preferred_posture[4]);
      double a4 = preferred_posture[4];  /* F*N or B*F */
      
      solns[buf_posture + e_o][4] = (solns[fun_posture + e_o][4] = a4);
      solns[bun_posture + e_o][4] = (solns[fuf_posture + e_o][4] = clamp_to_pi_range(a4 + M_PI));
      
      solns[bun_posture + e_o][5] = (solns[fun_posture + e_o][5] = 0.0);
      solns[buf_posture + e_o][5] = (solns[fuf_posture + e_o][5] = 0.0);
      
      double c6 = -s4 * (c23 * EE_y_proj_F[0] + s23 * EE_y_axis[2]) 
                 + c4 * EE_y_proj_F[1];
      double s6 = -s4 * EE_y_proj_F[1] 
                 - c4 * (EE_y_axis[2] * s23 + EE_y_proj_F[0] * c23);
      
      double a6 = atan2(s6, c6);  /* wrist not flipped */
      solns[bun_posture + e_o][6] = (solns[fun_posture + e_o][6] = a6);
      solns[buf_posture + e_o][6] = (solns[fuf_posture + e_o][6] = clamp_to_pi_range(a6 + M_PI));
    } else {
      /* we're not singular in jt 5 */
      double a4 = atan2(EE_z_proj_F[1], jt5_i);  /* F*N or B*F */
      double s4 = /* reach_sign *  */ sin(a4);
      double c4 = /* reach_sign *  */ cos(a4);
      
      solns[buf_posture + e_o][4] = (solns[fun_posture + e_o][4] = a4);
      solns[bun_posture + e_o][4] = (solns[fuf_posture + e_o][4] = clamp_to_pi_range(a4 + M_PI));
      
      double c5 =  EE_z_axis[2] * c23 - EE_z_proj_F[0] * s23;
      double s5 = -c4 * jt5_i - s4 * EE_z_proj_F[1]; 
      
      double a5 = atan2(s5, c5);  /* wrist not flipped */
      solns[bun_posture + e_o][5] = (solns[fun_posture + e_o][5] =  a5);
      solns[buf_posture + e_o][5] = (solns[fuf_posture + e_o][5] = -a5);
      
      double c6 = -s4 * (c23 * EE_y_proj_F[0] + s23 * EE_y_axis[2])
                 + c4 * EE_y_proj_F[1];
      double s6 =  s4 * (EE_z_proj_F[1] * (EE_y_axis[2] * c23 - EE_y_proj_F[0] * s23) - c5 * EE_y_proj_F[1])
                 + c4 * (EE_y_axis[2] * EE_z_proj_F[0] - EE_y_proj_F[0] * EE_z_axis[2]);
      
      double a6 = atan2(s6, c6);  /* wrist not flipped */
      solns[bun_posture + e_o][6] = (solns[fun_posture + e_o][6] = a6);
      solns[buf_posture + e_o][6] = (solns[fuf_posture + e_o][6] = clamp_to_pi_range(a6 + M_PI));
    };
  };
  
  
  // Now, we just need to choose a best solution. 
  // We'll choose it based on the most bounds-centric set of joint coordinates:
  double       best_cost = 100.0;
  unsigned int best_soln = 0;
  for(unsigned int i = 0; i < 8; ++i) {
    double cost = 0.0;
    for(unsigned int j = 0; j < 7; ++j) {
      if( ( solns[i][j] > joint_lower_bounds[j] ) && ( solns[i][j] < joint_upper_bounds[j] ) ) {
        cost += fabs( 2.0 * solns[i][j] - joint_lower_bounds[j] - joint_upper_bounds[j] ) / (joint_upper_bounds[j] - joint_lower_bounds[j]);
      } else {
        cost = 100.0; // effectively makes this one of the worst solutions.
        break;
      };
    };
    if(cost < best_cost) {
      best_cost = cost;
      best_soln = i;
    };
  };
  if(best_cost > 99.0)
    throw optim::infeasible_problem("Inverse kinematics problem is infeasible! None of the inverse kinematics solutions respect the joint limits.");
  
  return solns[best_soln];
};

vect_n<double> CRS_A465_model_builder::compute_inverse_kinematics(const frame_3D<double>& EE_state) const {
  
  // first, solve the IK problem for the zeroth-order terms:
  vect_n<double> result = compute_inverse_kinematics(static_cast<const pose_3D<double>&>(EE_state));
  
  // then, use the solution to compute the jacobian matrix:
  mat<double,mat_structure::rectangular> jac;
  compute_jacobian_matrix(result, jac);
  
  // finally, use the jacobian to find a minimum-norm solution for the joint velocities:
  mat<double,mat_structure::rectangular> x(7,1);
  mat<double,mat_structure::rectangular> b(6,1);
  b(0,0) = EE_state.Velocity[0]; b(1,0) = EE_state.Velocity[1]; b(2,0) = EE_state.Velocity[2];
  b(3,0) = EE_state.AngVelocity[0]; b(4,0) = EE_state.AngVelocity[1]; b(5,0) = EE_state.AngVelocity[2];
  minnorm_QR(jac,x,b,1e-4);
  result.resize(13);
  result[7]  = x(0,0); result[8]  = x(1,0); result[9]  = x(2,0);
  result[10] = x(3,0); result[11] = x(4,0); result[12] = x(5,0);
  
  return result;
};







};
  
  
};








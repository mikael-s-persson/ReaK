/**
 * \file build_CRS_A465_model.cpp
 *
 * This application constructs the KTE-based kinematics models for the CRS A465 manipulator.
 * The inertial information is phony and is only there for completeness of the model but do 
 * not reflect the actual inertial information of the CRS A465.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date June 2010
 */

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

#include "kinetostatics/motion_jacobians.hpp"
#include "mbd_kte/jacobian_joint_map.hpp"

#include "mbd_kte/manipulator_model.hpp"


#include "serialization/xml_archiver.hpp"

#include "rtti/typed_primitives.hpp"

using namespace ReaK;
using namespace kte;
using namespace serialization;
using namespace rtti;

int main() {
  
  //declare all the intermediate frames.
  shared_pointer< frame_3D<double> >::type robot_base(new frame_3D<double>(), scoped_deleter());
  shared_pointer< frame_3D<double> >::type track_joint_end(new frame_3D<double>(), scoped_deleter());
  shared_pointer< frame_3D<double> >::type arm_joint_1_base(new frame_3D<double>(), scoped_deleter());
  shared_pointer< frame_3D<double> >::type arm_joint_1_end(new frame_3D<double>(), scoped_deleter());
  shared_pointer< frame_3D<double> >::type arm_joint_2_base(new frame_3D<double>(), scoped_deleter());
  shared_pointer< frame_3D<double> >::type arm_joint_2_end(new frame_3D<double>(), scoped_deleter());
  shared_pointer< frame_3D<double> >::type arm_joint_3_base(new frame_3D<double>(), scoped_deleter());
  shared_pointer< frame_3D<double> >::type arm_joint_3_end(new frame_3D<double>(), scoped_deleter());
  shared_pointer< frame_3D<double> >::type arm_joint_4_base(new frame_3D<double>(), scoped_deleter());
  shared_pointer< frame_3D<double> >::type arm_joint_4_end(new frame_3D<double>(), scoped_deleter());
  shared_pointer< frame_3D<double> >::type arm_joint_5_base(new frame_3D<double>(), scoped_deleter());
  shared_pointer< frame_3D<double> >::type arm_joint_5_end(new frame_3D<double>(), scoped_deleter());
  shared_pointer< frame_3D<double> >::type arm_joint_6_base(new frame_3D<double>(), scoped_deleter());
  shared_pointer< frame_3D<double> >::type arm_joint_6_end(new frame_3D<double>(), scoped_deleter());
  shared_pointer< frame_3D<double> >::type arm_EE(new frame_3D<double>(), scoped_deleter());

  //declare all the joint coordinates.
  shared_pointer<gen_coord<double> >::type track_joint_coord(new gen_coord<double>(), scoped_deleter());
  shared_pointer<gen_coord<double> >::type arm_joint_1_coord(new gen_coord<double>(), scoped_deleter());
  shared_pointer<gen_coord<double> >::type arm_joint_2_coord(new gen_coord<double>(), scoped_deleter());
  shared_pointer<gen_coord<double> >::type arm_joint_3_coord(new gen_coord<double>(), scoped_deleter());
  shared_pointer<gen_coord<double> >::type arm_joint_4_coord(new gen_coord<double>(), scoped_deleter()); 
  shared_pointer<gen_coord<double> >::type arm_joint_5_coord(new gen_coord<double>(), scoped_deleter()); 
  shared_pointer<gen_coord<double> >::type arm_joint_6_coord(new gen_coord<double>(), scoped_deleter()); 
  
  //declare all the joint jacobians.
  shared_pointer<jacobian_gen_3D<double> >::type track_joint_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  shared_pointer<jacobian_gen_3D<double> >::type arm_joint_1_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  shared_pointer<jacobian_gen_3D<double> >::type arm_joint_2_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  shared_pointer<jacobian_gen_3D<double> >::type arm_joint_3_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  shared_pointer<jacobian_gen_3D<double> >::type arm_joint_4_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  shared_pointer<jacobian_gen_3D<double> >::type arm_joint_5_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  shared_pointer<jacobian_gen_3D<double> >::type arm_joint_6_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  
  //set the absolute position of the torso base and add gravity (z-axis pointing up!) (x-axis pointing forward).
  //normally this frame is set via the feedback from the leg / lower-body motion control of wopa
  robot_base->Acceleration = vect<double,3>(0.0,0.0,9.81); //put gravity acceleration on base of torso...
  robot_base->Position = vect<double,3>(0.0,0.0,0.8); //put the base of the torso somewhere 0.8 m off the ground level.
  
  
/******************************************************************************************************************
                        TORSO DYNAMIC MODEL FROM BASE OF TORSO TO THE TWO SHOULDER JOINTS
******************************************************************************************************************/
  
  //create revolute joint
  shared_pointer<prismatic_joint_3D>::type track_joint(new prismatic_joint_3D("track_joint",
                                                                              track_joint_coord,
                                                                              vect<double,3>(1.0,0.0,0.0),
                                                                              robot_base,
                                                                              track_joint_end,
                                                                              track_joint_jacobian),
                                                       scoped_deleter());
                                              
  //create motor inertia
  shared_pointer< joint_dependent_gen_coord >::type track_joint_dep_coord(new joint_dependent_gen_coord(track_joint_coord), 
									  scoped_deleter() );
  track_joint_dep_coord->add_joint(track_joint_coord, 
				   shared_pointer< jacobian_gen_gen<double> >::type(new jacobian_gen_gen<double>(1.0,0.0), 
										    scoped_deleter()));
  shared_pointer<inertia_gen>::type track_joint_inertia(new inertia_gen("track_joint_inertia",
                                                                        track_joint_dep_coord,
                                                                        1.0), 
                                                        scoped_deleter());
  
  //create force actuator
  shared_pointer<driving_actuator_gen>::type track_actuator(new driving_actuator_gen("track_actuator",
                                                                                     track_joint_coord,
                                                                                     track_joint),
                                                            scoped_deleter());
  
  //create link from G to F10.0
  shared_pointer<rigid_link_3D>::type link_0(new rigid_link_3D("link_0",
                                                               track_joint_end,
                                                               arm_joint_1_base,
                                                               pose_3D<double>(weak_pointer<pose_3D<double> >::type(),
                                                                               vect<double,3>(0.0,0.0,0.0),
                                                                               quaternion<double>())),
                                             scoped_deleter());
  
  //create arm-base inertia of 
  shared_pointer< joint_dependent_frame_3D >::type link_0_dep_frame(new joint_dependent_frame_3D(arm_joint_1_base),
								    scoped_deleter());
  link_0_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  shared_pointer<inertia_3D>::type link_0_inertia(new inertia_3D("link_0_inertia",
                                                                 link_0_dep_frame,
                                                                 1.0,
                                                                 mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0)),
                                                 scoped_deleter());
  
  //create revolute joint
  shared_pointer<revolute_joint_3D>::type arm_joint_1(new revolute_joint_3D("arm_joint_1",
                                                                            arm_joint_1_coord,
                                                                            vect<double,3>(0.0,0.0,1.0),
                                                                            arm_joint_1_base,
                                                                            arm_joint_1_end,
                                                                            arm_joint_1_jacobian),
                                                      scoped_deleter());
                                              
  //create motor inertia
  shared_pointer< joint_dependent_gen_coord >::type arm_joint_1_dep_coord(new joint_dependent_gen_coord(arm_joint_1_coord),
									  scoped_deleter());
  arm_joint_1_dep_coord->add_joint(arm_joint_1_coord,
                                   shared_pointer<jacobian_gen_gen<double> >::type(new jacobian_gen_gen<double>(1.0,0.0), 
										   scoped_deleter()));
  shared_pointer<inertia_gen>::type arm_joint_1_inertia(new inertia_gen("arm_joint_1_inertia",
                                                                        arm_joint_1_dep_coord,
                                                                        1.0), //~71 kg m^2
                                                        scoped_deleter());
  
  //create force actuator
  shared_pointer<driving_actuator_gen>::type arm_joint_1_actuator(new driving_actuator_gen("arm_joint_1_actuator",
                                                                                           arm_joint_1_coord,
                                                                                           arm_joint_1),
                                                                  scoped_deleter());
  
  //create link from F to CM (note that this is very approximate!!!)
  shared_pointer<rigid_link_3D>::type link_1(new rigid_link_3D("link_1",
                                                                arm_joint_1_end,
                                                                arm_joint_2_base,
                                                                pose_3D<double>(weak_pointer<pose_3D<double> >::type(),
                                                                                vect<double,3>(0.0,0.0,0.3302),
                                                                                axis_angle<double>(0.5 * M_PI, vect<double,3>(1.0,0.0,0.0)).getQuaternion())),
                                             scoped_deleter());
  
  //create link1 inertia 
  shared_pointer< joint_dependent_frame_3D >::type link_1_dep_frame(new joint_dependent_frame_3D(arm_joint_2_base),
								    scoped_deleter());
  link_1_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_1_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  shared_pointer<inertia_3D>::type link_1_inertia(new inertia_3D("link_1_inertia",
                                                                 link_1_dep_frame,
                                                                 1.0,
                                                                 mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0)),
                                                  scoped_deleter());

  //create revolute joint
  shared_pointer<revolute_joint_3D>::type arm_joint_2(new revolute_joint_3D("arm_joint_2",
                                                                             arm_joint_2_coord,
                                                                             vect<double,3>(0.0,0.0,1.0),
                                                                             arm_joint_2_base,
                                                                             arm_joint_2_end,
                                                                             arm_joint_2_jacobian),
                                                      scoped_deleter());
    
  //create motor inertia
  shared_pointer< joint_dependent_gen_coord >::type arm_joint_2_dep_coord(new joint_dependent_gen_coord(arm_joint_2_coord),
									  scoped_deleter());
  arm_joint_2_dep_coord->add_joint(arm_joint_2_coord,
                                   shared_pointer<jacobian_gen_gen<double> >::type(new jacobian_gen_gen<double>(1.0,0.0), 
										   scoped_deleter()));
  shared_pointer<inertia_gen>::type arm_joint_2_inertia(new inertia_gen("arm_joint_2_inertia",
                                                                        arm_joint_2_dep_coord,
                                                                        1.0),
                                                        scoped_deleter());
  
  //create force actuator
  shared_pointer<driving_actuator_gen>::type arm_joint_2_actuator(new driving_actuator_gen("arm_joint_2_actuator",
                                                                                           arm_joint_2_coord,
                                                                                           arm_joint_2),
                                                                  scoped_deleter());
  
  //create link 
  shared_pointer<rigid_link_3D>::type link_2(new rigid_link_3D("link_2",
                                                                arm_joint_2_end,
                                                                arm_joint_3_base,
                                                                pose_3D<double>(weak_pointer<pose_3D<double> >::type(),
                                                                                vect<double,3>(0.3048,0.0,0.0),
                                                                                quaternion<double>())),
                                              scoped_deleter());
  
  //create inertia
  shared_pointer< joint_dependent_frame_3D >::type link_2_dep_frame(new joint_dependent_frame_3D(arm_joint_3_base),
								    scoped_deleter());
  link_2_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_2_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  link_2_dep_frame->add_joint(arm_joint_2_coord,arm_joint_2_jacobian);
  shared_pointer<inertia_3D>::type link_2_inertia(new inertia_3D("link_2_inertia",
                                                                 link_2_dep_frame,
                                                                 1.0,
                                                                 mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0)),
                                                  scoped_deleter());

  //create revolute joint
  shared_pointer<revolute_joint_3D>::type arm_joint_3(new revolute_joint_3D("arm_joint_3",
                                                                            arm_joint_3_coord,
                                                                            vect<double,3>(0.0,0.0,1.0),
                                                                            arm_joint_3_base,
                                                                            arm_joint_3_end,
                                                                            arm_joint_3_jacobian),
                                                      scoped_deleter());
  
  //create motor inertia
  shared_pointer< joint_dependent_gen_coord >::type arm_joint_3_dep_coord(new joint_dependent_gen_coord(arm_joint_3_coord),
									  scoped_deleter());
  arm_joint_3_dep_coord->add_joint(arm_joint_3_coord,
                                   shared_pointer<jacobian_gen_gen<double> >::type(new jacobian_gen_gen<double>(1.0,0.0), 
										   scoped_deleter()));
  shared_pointer<inertia_gen>::type arm_joint_3_inertia(new inertia_gen("arm_joint_3_inertia",
                                                                         arm_joint_3_dep_coord,
                                                                         1.0), 
                                                        scoped_deleter());
  
  //create force actuator
  shared_pointer<driving_actuator_gen>::type arm_joint_3_actuator(new driving_actuator_gen("arm_joint_3_actuator",
                                                                                           arm_joint_3_coord,
                                                                                           arm_joint_3),
                                                                  scoped_deleter());
  
  //create link 
  shared_pointer<rigid_link_3D>::type link_3(new rigid_link_3D("link_3",
                                                               arm_joint_3_end,
                                                               arm_joint_4_base,
                                                               pose_3D<double>(weak_pointer<pose_3D<double> >::type(),
                                                                               vect<double,3>(0.1500,0.0,0.0),
                                                                               axis_angle<double>(0.5 * M_PI,vect<double,3>(0.0,1.0,0.0)) * axis_angle<double>(0.5 * M_PI,vect<double,3>(0.0,0.0,1.0)))),
                                              scoped_deleter());
  
  //create inertia
  shared_pointer< joint_dependent_frame_3D >::type link_3_dep_frame(new joint_dependent_frame_3D(arm_joint_4_base),
								    scoped_deleter());
  link_3_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_3_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  link_3_dep_frame->add_joint(arm_joint_2_coord,arm_joint_2_jacobian);
  link_3_dep_frame->add_joint(arm_joint_3_coord,arm_joint_3_jacobian);
  shared_pointer<inertia_3D>::type link_3_inertia(new inertia_3D("link_3_inertia",
                                                                 link_3_dep_frame,
                                                                 1.0,
                                                                 mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0)),
                                                  scoped_deleter());

  //create revolute joint
  shared_pointer<revolute_joint_3D>::type arm_joint_4(new revolute_joint_3D("arm_joint_4",
                                                                            arm_joint_4_coord,
                                                                            vect<double,3>(0.0,0.0,1.0),
                                                                            arm_joint_4_base,
                                                                            arm_joint_4_end,
                                                                            arm_joint_4_jacobian),
                                                      scoped_deleter());
  
  //create motor inertia
  shared_pointer< joint_dependent_gen_coord >::type arm_joint_4_dep_coord(new joint_dependent_gen_coord(arm_joint_4_coord),
									  scoped_deleter());
  arm_joint_4_dep_coord->add_joint(arm_joint_4_coord,
                                   shared_pointer<jacobian_gen_gen<double> >::type(new jacobian_gen_gen<double>(1.0,0.0), 
										   scoped_deleter()));
  shared_pointer<inertia_gen>::type arm_joint_4_inertia(new inertia_gen("arm_joint_4_inertia",
                                                                        arm_joint_4_dep_coord,
                                                                        1.0), 
                                                        scoped_deleter());
  
  //create force actuator
  shared_pointer<driving_actuator_gen>::type arm_joint_4_actuator(new driving_actuator_gen("arm_joint_4_actuator",
                                                                                           arm_joint_4_coord,
                                                                                           arm_joint_4),
                                                                  scoped_deleter());
  
  //create link 
  shared_pointer<rigid_link_3D>::type link_4(new rigid_link_3D("link_4",
                                                               arm_joint_4_end,
                                                               arm_joint_5_base,
                                                               pose_3D<double>(weak_pointer<pose_3D<double> >::type(),
                                                                               vect<double,3>(0.0,0.0,0.1802),
                                                                               axis_angle<double>(-0.5*M_PI,vect<double,3>(1.0,0.0,0.0)).getQuaternion())),
                                             scoped_deleter());
  
  //create inertia
  shared_pointer< joint_dependent_frame_3D >::type link_4_dep_frame(new joint_dependent_frame_3D(arm_joint_5_base),
								    scoped_deleter());
  link_4_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_4_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  link_4_dep_frame->add_joint(arm_joint_2_coord,arm_joint_2_jacobian);
  link_4_dep_frame->add_joint(arm_joint_3_coord,arm_joint_3_jacobian);
  link_4_dep_frame->add_joint(arm_joint_4_coord,arm_joint_4_jacobian);
  shared_pointer<inertia_3D>::type link_4_inertia(new inertia_3D("link_4_inertia",
                                                                 link_4_dep_frame,
                                                                 1.0,
                                                                 mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0)),
                                                  scoped_deleter());

  //create revolute joint
  shared_pointer<revolute_joint_3D>::type arm_joint_5(new revolute_joint_3D("arm_joint_5",
                                                                            arm_joint_5_coord,
                                                                            vect<double,3>(0.0,0.0,1.0),
                                                                            arm_joint_5_base,
                                                                            arm_joint_5_end,
                                                                            arm_joint_5_jacobian),
                                                      scoped_deleter());
  
  //create motor inertia
  shared_pointer< joint_dependent_gen_coord >::type arm_joint_5_dep_coord(new joint_dependent_gen_coord(arm_joint_5_coord),
									  scoped_deleter());
  arm_joint_5_dep_coord->add_joint(arm_joint_5_coord,
                                   shared_pointer<jacobian_gen_gen<double> >::type(new jacobian_gen_gen<double>(1.0,0.0), 
										   scoped_deleter()));
  shared_pointer<inertia_gen>::type arm_joint_5_inertia(new inertia_gen("arm_joint_5_inertia",
									arm_joint_5_dep_coord,
									1.0),
                                                        scoped_deleter());
  
  //create force actuator
  shared_pointer<driving_actuator_gen>::type arm_joint_5_actuator(new driving_actuator_gen("arm_joint_5_actuator",
											   arm_joint_5_coord,
											   arm_joint_5),
                                                                  scoped_deleter());
  
  //create link 
  shared_pointer<rigid_link_3D>::type link_5(new rigid_link_3D("link_5",
                                                               arm_joint_5_end,
                                                               arm_joint_6_base,
                                                               pose_3D<double>(weak_pointer<pose_3D<double> >::type(),
                                                                               vect<double,3>(0.0762,0.0,0.0),
                                                                               axis_angle<double>(0.5 * M_PI,vect<double,3>(0.0,1.0,0.0)) * axis_angle<double>(0.5 * M_PI,vect<double,3>(0.0,0.0,1.0)))),
                                             scoped_deleter());
  
  //create inertia
  shared_pointer< joint_dependent_frame_3D >::type link_5_dep_frame(new joint_dependent_frame_3D(arm_joint_6_base),
								    scoped_deleter());
  link_5_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_5_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  link_5_dep_frame->add_joint(arm_joint_2_coord,arm_joint_2_jacobian);
  link_5_dep_frame->add_joint(arm_joint_3_coord,arm_joint_3_jacobian);
  link_5_dep_frame->add_joint(arm_joint_4_coord,arm_joint_4_jacobian);
  link_5_dep_frame->add_joint(arm_joint_5_coord,arm_joint_5_jacobian);
  shared_pointer<inertia_3D>::type link_5_inertia(new inertia_3D("link_5_inertia",
                                                                 link_5_dep_frame,
                                                                 1.0,
                                                                 mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0)),
                                                  scoped_deleter());
  
  //create revolute joint
  shared_pointer<revolute_joint_3D>::type arm_joint_6(new revolute_joint_3D("arm_joint_6",
                                                                            arm_joint_6_coord,
                                                                            vect<double,3>(0.0,0.0,1.0),
                                                                            arm_joint_6_base,
                                                                            arm_joint_6_end,
                                                                            arm_joint_6_jacobian),
                                                      scoped_deleter());
  
  //create motor inertia
  shared_pointer< joint_dependent_gen_coord >::type arm_joint_6_dep_coord(new joint_dependent_gen_coord(arm_joint_6_coord),
									  scoped_deleter());
  arm_joint_6_dep_coord->add_joint(arm_joint_6_coord,
                                   shared_pointer<jacobian_gen_gen<double> >::type(new jacobian_gen_gen<double>(1.0,0.0), 
										   scoped_deleter()));
  shared_pointer<inertia_gen>::type arm_joint_6_inertia(new inertia_gen("arm_joint_6_inertia",
									arm_joint_6_dep_coord,
                                                                        1.0),
                                                        scoped_deleter());
  
  //create force actuator
  shared_pointer<driving_actuator_gen>::type arm_joint_6_actuator(new driving_actuator_gen("arm_joint_6_actuator",
                                                                                           arm_joint_6_coord,
                                                                                           arm_joint_6),
                                                                  scoped_deleter());
  
  //create link 
  shared_pointer<rigid_link_3D>::type link_6(new rigid_link_3D("link_6",
                                                                arm_joint_6_end,
                                                                arm_EE,
                                                                pose_3D<double>(weak_pointer<pose_3D<double> >::type(),
                                                                                vect<double,3>(0.0,0.0,0.0),
                                                                                quaternion<double>())),
                                             scoped_deleter());
  
  //create inertia
  shared_pointer< joint_dependent_frame_3D >::type link_6_dep_frame(new joint_dependent_frame_3D(arm_EE),
								    scoped_deleter());
  link_6_dep_frame->add_joint(track_joint_coord,track_joint_jacobian);
  link_6_dep_frame->add_joint(arm_joint_1_coord,arm_joint_1_jacobian);
  link_6_dep_frame->add_joint(arm_joint_2_coord,arm_joint_2_jacobian);
  link_6_dep_frame->add_joint(arm_joint_3_coord,arm_joint_3_jacobian);
  link_6_dep_frame->add_joint(arm_joint_4_coord,arm_joint_4_jacobian);
  link_6_dep_frame->add_joint(arm_joint_5_coord,arm_joint_5_jacobian);
  link_6_dep_frame->add_joint(arm_joint_6_coord,arm_joint_6_jacobian);
  shared_pointer<inertia_3D>::type link_6_inertia(new inertia_3D("link_6_inertia",
                                                                 link_6_dep_frame,
                                                                 1.0,
                                                                 mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0)),
                                                  scoped_deleter());
  
  shared_pointer<kte_map_chain>::type CRS_A465_dyn_model(new kte_map_chain("CRS_A465_dyn_model"),
							 scoped_deleter());
  
  
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
  
  shared_pointer< manipulator_dynamics_model >::type CRS_A465_dyn_manip(new manipulator_dynamics_model("CRS_A465_dyn_manip"),
									scoped_deleter());
  
  CRS_A465_dyn_manip->setModel(CRS_A465_dyn_model);
  
   
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
                      
  
  shared_pointer< kte_map_chain >::type CRS_A465_kin_model(new kte_map_chain("CRS_A465_kin_model"),
                                                           scoped_deleter());
  
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
		      
  shared_pointer< manipulator_kinematics_model >::type CRS_A465_kin_manip(new manipulator_kinematics_model("CRS_A465_kin_manip"),
									  scoped_deleter());
  
  CRS_A465_kin_manip->setModel(CRS_A465_kin_model);
  
  *CRS_A465_kin_manip 
//Register joint coordinates:
                      << track_joint_coord
                      << arm_joint_1_coord
                      << arm_joint_2_coord
                      << arm_joint_3_coord
                      << arm_joint_4_coord
                      << arm_joint_5_coord
                      << arm_joint_6_coord
//Register joint-dependent link-frames:
                      << link_0_dep_frame  
                      << link_1_dep_frame
                      << link_2_dep_frame
                      << link_3_dep_frame
                      << link_4_dep_frame
                      << link_5_dep_frame
                      << link_6_dep_frame;
  
  shared_pointer< mass_matrix_calc >::type CRS_A465_Mlink_calc(new mass_matrix_calc("CRS_A465_Mlink_calc"),
                                                               scoped_deleter());
  
  *CRS_A465_Mlink_calc << link_0_inertia
                       << link_1_inertia
                       << link_2_inertia
                       << link_3_inertia
                       << link_4_inertia
                       << link_5_inertia
                       << link_6_inertia;
  *CRS_A465_Mlink_calc << track_joint_coord
                       << arm_joint_1_coord
                       << arm_joint_2_coord
                       << arm_joint_3_coord
                       << arm_joint_4_coord
                       << arm_joint_5_coord
                       << arm_joint_6_coord;
		      
  shared_pointer< mass_matrix_calc >::type CRS_A465_Mjoint_calc(new mass_matrix_calc("CRS_A465_Mjoint_calc"),
                                                                scoped_deleter());
  
  *CRS_A465_Mjoint_calc << track_joint_inertia
                        << arm_joint_1_inertia
                        << arm_joint_2_inertia
                        << arm_joint_3_inertia
                        << arm_joint_4_inertia
                        << arm_joint_5_inertia
                        << arm_joint_6_inertia;
  *CRS_A465_Mjoint_calc << track_joint_coord
                        << arm_joint_1_coord
                        << arm_joint_2_coord
                        << arm_joint_3_coord
                        << arm_joint_4_coord
                        << arm_joint_5_coord
                        << arm_joint_6_coord;
                    
  //Now save all to various files.
  
  { //save the complete dynamic model.
    xml_oarchive complete_model_output("models/CRS_A465_dyn_model.xml");
    oarchive& output_ref = complete_model_output;
    
    output_ref & RK_SERIAL_SAVE_WITH_NAME(robot_base)
               & RK_SERIAL_SAVE_WITH_NAME(track_joint_coord)
               & RK_SERIAL_SAVE_WITH_NAME(track_actuator)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1_coord)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1_actuator)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2_coord)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2_actuator)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3_coord)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3_actuator)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_4_coord)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_4_actuator)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_5_coord)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_5_actuator)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_6_coord)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_6_actuator)
               & RK_SERIAL_SAVE_WITH_NAME(arm_EE)
               & RK_SERIAL_SAVE_WITH_NAME(CRS_A465_dyn_model)
               & RK_SERIAL_SAVE_WITH_NAME(CRS_A465_Mlink_calc)
               & RK_SERIAL_SAVE_WITH_NAME(CRS_A465_Mjoint_calc);
  };
  
  { //save the complete kinematic model.
    xml_oarchive complete_model_output("models/CRS_A465_kin_model.xml");
    oarchive& output_ref = complete_model_output;
    
    output_ref & RK_SERIAL_SAVE_WITH_NAME(robot_base)
               & RK_SERIAL_SAVE_WITH_NAME(track_joint_coord)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_1_coord)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_2_coord)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_3_coord)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_4_coord)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_5_coord)
               & RK_SERIAL_SAVE_WITH_NAME(arm_joint_6_coord)
               & RK_SERIAL_SAVE_WITH_NAME(arm_EE)
               & RK_SERIAL_SAVE_WITH_NAME(link_0_dep_frame)
               & RK_SERIAL_SAVE_WITH_NAME(link_1_dep_frame)
               & RK_SERIAL_SAVE_WITH_NAME(link_2_dep_frame)
               & RK_SERIAL_SAVE_WITH_NAME(link_3_dep_frame)
               & RK_SERIAL_SAVE_WITH_NAME(link_4_dep_frame)
               & RK_SERIAL_SAVE_WITH_NAME(link_5_dep_frame)
               & RK_SERIAL_SAVE_WITH_NAME(link_6_dep_frame)
               & RK_SERIAL_SAVE_WITH_NAME(CRS_A465_kin_model);
  };
  
  { //save the complete dynamic model.
    xml_oarchive complete_model_output("models/CRS_A465_dyn_manipulator.xml");
    oarchive& output_ref = complete_model_output;
    
    output_ref & RK_SERIAL_SAVE_WITH_NAME(robot_base)
               & RK_SERIAL_SAVE_WITH_ALIAS("manipulator_model",CRS_A465_dyn_manip);
  };
  
  { //save the complete kinematic model.
    xml_oarchive complete_model_output("models/CRS_A465_kin_manipulator.xml");
    oarchive& output_ref = complete_model_output;
    
    output_ref & RK_SERIAL_SAVE_WITH_NAME(robot_base)
               & RK_SERIAL_SAVE_WITH_ALIAS("manipulator_model",CRS_A465_kin_manip);
  };
  
  
  
};









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



#include "CRS_A465_geom_model.hpp"

#include "serialization/xml_archiver.hpp"

#include "rtti/typed_primitives.hpp"


namespace ReaK {


namespace robot_airship {


void CRS_A465_geom_builder::load_geom_from_archive(serialization::iarchive& aInput) {
  aInput & RK_SERIAL_LOAD_WITH_NAME(geom_model)
         & RK_SERIAL_LOAD_WITH_NAME(proxy_model);
};

void CRS_A465_geom_builder::save_geom_to_archive(serialization::oarchive& aOutput) const {
  aOutput & RK_SERIAL_SAVE_WITH_NAME(geom_model)
          & RK_SERIAL_SAVE_WITH_NAME(proxy_model);
};


void CRS_A465_geom_builder::load_kte_and_geom(const std::string& aFileName) {
  
  serialization::xml_iarchive complete_model_input(aFileName);
  load_kte_from_archive(complete_model_input);
  load_geom_from_archive(complete_model_input);
  
};


void CRS_A465_geom_builder::create_geom_from_preset() {
  
  if(!arm_EE)
    create_from_preset();
  
  global_frame_arrows = shared_ptr< geom::coord_arrows_3D >(new geom::coord_arrows_3D("global_frame_arrows",shared_ptr< pose_3D<double> >(),pose_3D<double>(),0.3));
  robot_base_arrows = shared_ptr< geom::coord_arrows_3D >(new geom::coord_arrows_3D("robot_base_arrows",robot_base,pose_3D<double>(),0.3));
  track_joint_arrows = shared_ptr< geom::coord_arrows_3D >(new geom::coord_arrows_3D("track_joint_arrows",track_joint_end,pose_3D<double>(),0.3));
  arm_joint_1_arrows = shared_ptr< geom::coord_arrows_3D >(new geom::coord_arrows_3D("arm_joint_1_arrows",arm_joint_1_end,pose_3D<double>(),0.3));
  arm_joint_2_arrows = shared_ptr< geom::coord_arrows_3D >(new geom::coord_arrows_3D("arm_joint_2_arrows",arm_joint_2_end,pose_3D<double>(),0.3));
  arm_joint_3_arrows = shared_ptr< geom::coord_arrows_3D >(new geom::coord_arrows_3D("arm_joint_3_arrows",arm_joint_3_end,pose_3D<double>(),0.3));
  arm_joint_4_arrows = shared_ptr< geom::coord_arrows_3D >(new geom::coord_arrows_3D("arm_joint_4_arrows",arm_joint_4_end,pose_3D<double>(),0.3));
  arm_joint_5_arrows = shared_ptr< geom::coord_arrows_3D >(new geom::coord_arrows_3D("arm_joint_5_arrows",arm_joint_5_end,pose_3D<double>(),0.3));
  arm_joint_6_arrows = shared_ptr< geom::coord_arrows_3D >(new geom::coord_arrows_3D("arm_joint_6_arrows",arm_joint_6_end,pose_3D<double>(),0.3));
  
  link1_cyl  = shared_ptr< geom::capped_cylinder >(new geom::capped_cylinder("link1_cyl", track_joint_end, 
    pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0,0.0,0.15),quaternion<double>()), 
    0.3, 0.1));
  joint2_cyl = shared_ptr< geom::capped_cylinder >(new geom::capped_cylinder("joint2_cyl", arm_joint_1_end,
    pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0,0.0,0.3302),axis_angle<double>(M_PI * 0.5, vect<double,3>(1.0,0.0,0.0)).getQuaternion()), 
    0.34, 0.09));
  link2_cyl  = shared_ptr< geom::capped_cylinder >(new geom::capped_cylinder("link2_cyl",arm_joint_2_end,
    pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0,0.0,0.15),quaternion<double>()), 
    0.3, 0.07));
  link3_cyl  = shared_ptr< geom::capped_cylinder >(new geom::capped_cylinder("link3_cyl",arm_joint_3_end,
    pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0,0.0,0.165),quaternion<double>()), 
    0.33, 0.07));
  link5_cyl  = shared_ptr< geom::capped_cylinder >(new geom::capped_cylinder("link5_cyl",arm_joint_5_end,
    pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0,0.0,0.0381),quaternion<double>()), 
    0.0762, 0.05));
  EE_sphere  = shared_ptr< geom::sphere >(new geom::sphere("EE_sphere", arm_joint_6_end, 
    pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(-0.04,0.0,0.05),quaternion<double>()), 
    0.11));
  
  EE_bumblebee = shared_ptr< geom::box >(new geom::box("EE_bumblebee", arm_joint_6_end, 
    pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(-0.105,0.0,0.0483),quaternion<double>()), 
    vect<double,3>(0.04,0.14,0.04)));
  EE_bumblebee_support = shared_ptr< geom::box >(new geom::box("EE_bumblebee_support", arm_joint_6_end, 
    pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(-0.0525,0.0,0.0393),quaternion<double>()), 
    vect<double,3>(0.065,0.06,0.03)));
  EE_gripper_box = shared_ptr< geom::box >(new geom::box("EE_gripper_box", arm_joint_6_end, 
    pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0,0.0,0.0669),quaternion<double>()), 
    vect<double,3>(0.04,0.08,0.1338)));
  EE_gripper_fingers = shared_ptr< geom::box >(new geom::box("EE_gripper_fingers", arm_joint_6_end, 
    pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0,0.0,0.1588),quaternion<double>()), 
    vect<double,3>(0.02,0.01,0.05)));
  
  geom_model = shared_ptr< geom::colored_model_3D >(new geom::colored_model_3D("CRS_A465_model_render"));
  
  (*geom_model)
   .addElement(geom::color(0,0,0),global_frame_arrows)
   .addElement(geom::color(0,0,0),robot_base_arrows)
   .addElement(geom::color(0,0,0),track_joint_arrows)
   .addElement(geom::color(0,0,0),arm_joint_1_arrows)
   .addElement(geom::color(0,0,0),arm_joint_2_arrows)
   .addElement(geom::color(0,0,0),arm_joint_3_arrows)
   .addElement(geom::color(0,0,0),arm_joint_4_arrows)
   .addElement(geom::color(0,0,0),arm_joint_5_arrows)
   .addElement(geom::color(0,0,0),arm_joint_6_arrows)
   .addElement(geom::color(0.2,0.2,0.2),link1_cyl)
   .addElement(geom::color(0.3,0.3,0.3),joint2_cyl)
   .addElement(geom::color(0.7,0.7,0.5),link2_cyl)
   .addElement(geom::color(0.7,0.7,0.5),link3_cyl)
   .addElement(geom::color(0.1,0.1,0.1),link5_cyl)
   //.addElement(geom::color(0.6,0.6,0.0),EE_sphere)
   .addElement(geom::color(0.6,0.6,0.0),EE_bumblebee)
   .addElement(geom::color(0.5,0.5,0.5),EE_bumblebee_support)
   .addElement(geom::color(0.1,0.1,0.1),EE_gripper_box)
   .addElement(geom::color(0.5,0.5,0.5),EE_gripper_fingers);
  
  proxy_model = shared_ptr< geom::proxy_query_model_3D >(new geom::proxy_query_model_3D("CRS_A465_model_proxy"));
  
  (*proxy_model)
   .addShape(joint2_cyl)
   .addShape(link2_cyl)
   .addShape(link3_cyl)
   .addShape(link5_cyl)
   .addShape(EE_sphere);
  
};
    
void CRS_A465_geom_builder::save_kte_and_geom(const std::string& aFileName) const {
  
  serialization::xml_oarchive complete_model_output(aFileName);
  save_kte_to_archive(complete_model_output);
  save_geom_to_archive(complete_model_output);
  
};

shared_ptr< geom::colored_model_3D > CRS_A465_geom_builder::get_geometric_model() const {
  return geom_model;
};

shared_ptr< geom::proxy_query_model_3D > CRS_A465_geom_builder::get_proximity_model() const {
  return proxy_model;
};





};
  
  
};








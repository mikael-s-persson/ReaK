
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


#include "base/defs.hpp"

#include "X8_quadrotor_geom.hpp"

#include "kinetostatics/frame_3D.hpp"
#include "kinetostatics/pose_3D.hpp"
#include "kinetostatics/rotations_3D.hpp"

#include "serialization/xml_archiver.hpp"

#include "rtti/typed_primitives.hpp"

#include "shapes/colored_model.hpp"
#include "shapes/sphere.hpp"
#include "shapes/box.hpp"
#include "shapes/capped_cylinder.hpp"
#include "shapes/cylinder.hpp"
#include "shapes/coord_arrows_3D.hpp"
#include "proximity/proxy_query_model.hpp"


namespace ReaK {


namespace geom {


void X8_quadrotor_geom::create_geom_from_preset(const kte::UAV_kinematics& aModel) {
  
  shared_ptr< frame_3D<double> > out_frame = aModel.getDependentFrame3D(0)->mFrame;
  shared_ptr< coord_arrows_3D > body_frame_arrows(new coord_arrows_3D("body_frame_arrows",out_frame,pose_3D<double>(),0.3));
  
  pose_3D<double> xz_pose(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.0,0.0,0.0),
      quaternion<double>::xrot(M_PI * 0.5).getQuaternion());
  
  pose_3D<double> armLF_pose = xz_pose;
  armLF_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.117,0.0,-0.117),
      quaternion<double>::yrot(M_PI * 0.75).getQuaternion()));
  
  shared_ptr< capped_cylinder > armLF(new capped_cylinder("X8_armLF", out_frame, armLF_pose, 0.331, 0.02));
                                               
  pose_3D<double> armRF_pose = xz_pose;
  armRF_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.117,0.0,0.117),
      quaternion<double>::yrot(M_PI * 0.25).getQuaternion()));
  
  shared_ptr< capped_cylinder > armRF(new capped_cylinder("X8_armRF", out_frame, armRF_pose, 0.331, 0.02));
                                            
  pose_3D<double> armRB_pose = xz_pose;
  armRB_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(-0.117,0.0,0.117),
      quaternion<double>::yrot(-M_PI * 0.25).getQuaternion()));
  
  shared_ptr< capped_cylinder > armRB(new capped_cylinder("X8_armRB", out_frame, armRB_pose, 0.331, 0.02));
                                               
  pose_3D<double> armLB_pose = xz_pose;
  armLB_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(-0.117,0.0,-0.117),
      quaternion<double>::yrot(-M_PI * 0.75).getQuaternion()));
  
  shared_ptr< capped_cylinder > armLB(new capped_cylinder("X8_armLB", out_frame, armLB_pose, 0.331, 0.02));
  
  armLF_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.0,0.0,0.1655),
      quaternion<double>::xrot(M_PI * 0.5).getQuaternion()));
  armRF_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.0,0.0,0.1655),
      quaternion<double>::xrot(M_PI * 0.5).getQuaternion()));
  armRB_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.0,0.0,0.1655),
      quaternion<double>::xrot(M_PI * 0.5).getQuaternion()));
  armLB_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.0,0.0,0.1655),
      quaternion<double>::xrot(M_PI * 0.5).getQuaternion()));
  
  shared_ptr< capped_cylinder > motorLF(new capped_cylinder("X8_motorLF", out_frame, armLF_pose, 0.1, 0.03));
  shared_ptr< capped_cylinder > motorRF(new capped_cylinder("X8_motorRF", out_frame, armRF_pose, 0.1, 0.03));
  shared_ptr< capped_cylinder > motorRB(new capped_cylinder("X8_motorRB", out_frame, armRB_pose, 0.1, 0.03));
  shared_ptr< capped_cylinder > motorLB(new capped_cylinder("X8_motorLB", out_frame, armLB_pose, 0.1, 0.03));
  
  armLF_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.0,0.0,0.05),
      quaternion<double>()));
  armRF_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.0,0.0,0.05),
      quaternion<double>()));
  armRB_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.0,0.0,0.05),
      quaternion<double>()));
  armLB_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.0,0.0,0.05),
      quaternion<double>()));
  
  shared_ptr< cylinder > UrotorLF(new cylinder("X8_UrotorLF", out_frame, armLF_pose, 0.01, 0.202));
  shared_ptr< cylinder > UrotorRF(new cylinder("X8_UrotorRF", out_frame, armRF_pose, 0.01, 0.202));
  shared_ptr< cylinder > UrotorRB(new cylinder("X8_UrotorRB", out_frame, armRB_pose, 0.01, 0.202));
  shared_ptr< cylinder > UrotorLB(new cylinder("X8_UrotorLB", out_frame, armLB_pose, 0.01, 0.202));
  
  armLF_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.0,0.0,-0.1),
      quaternion<double>()));
  armRF_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.0,0.0,-0.1),
      quaternion<double>()));
  armRB_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.0,0.0,-0.1),
      quaternion<double>()));
  armLB_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.0,0.0,-0.1),
      quaternion<double>()));
  
  shared_ptr< cylinder > LrotorLF(new cylinder("X8_LrotorLF", out_frame, armLF_pose, 0.01, 0.19));
  shared_ptr< cylinder > LrotorRF(new cylinder("X8_LrotorRF", out_frame, armRF_pose, 0.01, 0.19));
  shared_ptr< cylinder > LrotorRB(new cylinder("X8_LrotorRB", out_frame, armRB_pose, 0.01, 0.19));
  shared_ptr< cylinder > LrotorLB(new cylinder("X8_LrotorLB", out_frame, armLB_pose, 0.01, 0.19));
  
  
  
  shared_ptr< box > body_case(new box("X8_box", out_frame, pose_3D<double>(), vect<double,3>(0.15,0.15,0.15)));
  
  pose_3D<double> skid_pose = xz_pose;
  skid_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.0,0.25,0.0),
      quaternion<double>::yrot(M_PI * 0.5).getQuaternion()));
  
  skid_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(0.15,0.0,0.0),
      quaternion<double>()));
  
  shared_ptr< capped_cylinder > Lskid(new capped_cylinder("X8_Lskid", out_frame, skid_pose, 0.3, 0.01));
  
  skid_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
      vect<double,3>(-0.3,0.0,0.0),
      quaternion<double>()));
  
  shared_ptr< capped_cylinder > Rskid(new capped_cylinder("X8_Rskid", out_frame, skid_pose, 0.3, 0.01));
  
  
  shared_ptr< sphere > proxy_sphere(new sphere("X8_proxy_sphere", out_frame, pose_3D<double>(), 0.6));
  
  
  geom_model = shared_ptr< colored_model_3D >(new colored_model_3D("X8_model_render"));
  
  (*geom_model)
   //.addElement(color(0,0,0),global_frame_arrows)
   .addElement(color(0.1,0.1,0.1),armLF)
   .addElement(color(0.1,0.1,0.1),armRF)
   .addElement(color(0.1,0.1,0.1),armRB)
   .addElement(color(0.1,0.1,0.1),armLB)
   .addElement(color(0.0,1.0,0.0),motorLF)
   .addElement(color(1.0,0.0,0.0),motorRF)
   .addElement(color(0.1,0.1,0.1),motorRB)
   .addElement(color(0.1,0.1,0.1),motorLB)
   .addElement(color(0,0,0),UrotorLF)
   .addElement(color(0,0,0),UrotorRF)
   .addElement(color(0,0,0),UrotorRB)
   .addElement(color(0,0,0),UrotorLB)
   .addElement(color(0,0,0),LrotorLF)
   .addElement(color(0,0,0),LrotorRF)
   .addElement(color(0,0,0),LrotorRB)
   .addElement(color(0,0,0),LrotorLB)
   .addElement(color(0.1,0.1,0.1),body_case)
   .addElement(color(0.1,0.1,0.1),Lskid)
   .addElement(color(0.1,0.1,0.1),Rskid);
  
  proxy_model = shared_ptr< proxy_query_model_3D >(new proxy_query_model_3D("X8_model_proxy"));
  
  (*proxy_model)
   .addShape(proxy_sphere);
  
};



void RK_CALL X8_quadrotor_geom::save(ReaK::serialization::oarchive& A, unsigned int) const {
  named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(geom_model)
    & RK_SERIAL_SAVE_WITH_NAME(proxy_model);
};

void RK_CALL X8_quadrotor_geom::load(ReaK::serialization::iarchive& A, unsigned int) {
  named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(geom_model)
    & RK_SERIAL_LOAD_WITH_NAME(proxy_model);
};


};
  
  
};








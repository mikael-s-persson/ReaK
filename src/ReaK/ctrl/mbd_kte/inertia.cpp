
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



#include "inertia.hpp"

namespace ReaK {

namespace kte {




void inertia_gen::doMotion(kte_pass_flag aFlag, const shared_pointer<frame_storage>::type& aStorage) {
  if(!mCenterOfMass)
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->gen_coord_mapping[mCenterOfMass]))
      aStorage->gen_coord_mapping[mCenterOfMass] = shared_pointer< ReaK::gen_coord<double> >::type(new ReaK::gen_coord<double>((*mCenterOfMass)),ReaK::scoped_deleter());
    else
      (*(aStorage->gen_coord_mapping[mCenterOfMass])) = (*mCenterOfMass);
  };
};

void inertia_gen::doForce(kte_pass_flag aFlag, const shared_pointer<frame_storage>::type& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if(!mCenterOfMass)
    return;

  mCenterOfMass->f -= mCenterOfMass->q_ddot * mMass; //d'Alembert force

};


void inertia_gen::clearForce() {
  if(mCenterOfMass)
    mCenterOfMass->f = 0.0;
};





void inertia_2D::doMotion(kte_pass_flag aFlag, const shared_pointer<frame_storage>::type& aStorage) {
  if(!mCenterOfMass)
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_2D_mapping[mCenterOfMass]))
      aStorage->frame_2D_mapping[mCenterOfMass] = shared_pointer< frame_2D<double> >::type(new frame_2D<double>((*mCenterOfMass)),scoped_deleter());
    else
      (*(aStorage->frame_2D_mapping[mCenterOfMass])) = (*mCenterOfMass);
  };
};

void inertia_2D::doForce(kte_pass_flag aFlag, const shared_pointer<frame_storage>::type& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if(!mCenterOfMass)
    return;

  frame_2D<double> global_frame = mCenterOfMass->getGlobalFrame(); //get the inertial frame

  mCenterOfMass->Force -= mMass * (global_frame.Acceleration * global_frame.Rotation);
  mCenterOfMass->Torque -= mMomentOfInertia * global_frame.AngAcceleration;

};


void inertia_2D::clearForce() {
  if(mCenterOfMass) {
    mCenterOfMass->Force = vect<double,2>();
    mCenterOfMass->Torque = 0.0;
  };
};




void inertia_3D::doMotion(kte_pass_flag aFlag, const shared_pointer<frame_storage>::type& aStorage) {
  if(!mCenterOfMass)
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_3D_mapping[mCenterOfMass]))
      aStorage->frame_3D_mapping[mCenterOfMass] = shared_pointer< frame_3D<double> >::type(new frame_3D<double>((*mCenterOfMass)),scoped_deleter());
    else
      (*(aStorage->frame_3D_mapping[mCenterOfMass])) = (*mCenterOfMass);
  };
};

void inertia_3D::doForce(kte_pass_flag aFlag, const shared_pointer<frame_storage>::type& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if(!mCenterOfMass)
    return;

  frame_3D<double> global_frame = mCenterOfMass->getGlobalFrame(); //get the inertial frame

  mCenterOfMass->Force -= mMass * (invert(global_frame.Quat) * global_frame.Acceleration);

  mCenterOfMass->Torque -= mInertiaTensor * global_frame.AngAcceleration + global_frame.AngVelocity % (mInertiaTensor * global_frame.AngVelocity);

};


void inertia_3D::clearForce() {
  if(mCenterOfMass) {
    mCenterOfMass->Force = vect<double,3>();
    mCenterOfMass->Torque = vect<double,3>();
  };
};



};


};






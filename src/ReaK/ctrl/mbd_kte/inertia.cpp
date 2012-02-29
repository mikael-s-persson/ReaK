
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




void inertia_gen::doMotion(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) {
  if((!mCenterOfMass) || (!mCenterOfMass->mFrame))
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->gen_coord_mapping[mCenterOfMass->mFrame]))
      aStorage->gen_coord_mapping[mCenterOfMass->mFrame] = shared_ptr< ReaK::gen_coord<double> >(new ReaK::gen_coord<double>((*mCenterOfMass->mFrame)),ReaK::scoped_deleter());
    else
      (*(aStorage->gen_coord_mapping[mCenterOfMass->mFrame])) = (*mCenterOfMass->mFrame);
  };
};

void inertia_gen::doForce(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mCenterOfMass) || (!mCenterOfMass->mFrame))
    return;

  mCenterOfMass->mFrame->f -= mCenterOfMass->mFrame->q_ddot * mMass; //d'Alembert force

};


void inertia_gen::clearForce() {
  if((mCenterOfMass) && (mCenterOfMass->mFrame))
    mCenterOfMass->mFrame->f = 0.0;
};





void inertia_2D::doMotion(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) {
  if((!mCenterOfMass) || (!mCenterOfMass->mFrame))
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_2D_mapping[mCenterOfMass->mFrame]))
      aStorage->frame_2D_mapping[mCenterOfMass->mFrame] = shared_ptr< frame_2D<double> >(new frame_2D<double>((*mCenterOfMass->mFrame)),scoped_deleter());
    else
      (*(aStorage->frame_2D_mapping[mCenterOfMass->mFrame])) = (*mCenterOfMass->mFrame);
  };
};

void inertia_2D::doForce(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mCenterOfMass) || (!mCenterOfMass->mFrame))
    return;

  frame_2D<double> global_frame = mCenterOfMass->mFrame->getGlobalFrame(); //get the inertial frame

  mCenterOfMass->mFrame->Force -= mMass * (global_frame.Acceleration * global_frame.Rotation);
  mCenterOfMass->mFrame->Torque -= mMomentOfInertia * global_frame.AngAcceleration;

};


void inertia_2D::clearForce() {
  if((mCenterOfMass) && (mCenterOfMass->mFrame)) {
    mCenterOfMass->mFrame->Force = vect<double,2>();
    mCenterOfMass->mFrame->Torque = 0.0;
  };
};




void inertia_3D::doMotion(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) {
  if((!mCenterOfMass) || (!mCenterOfMass->mFrame))
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_3D_mapping[mCenterOfMass->mFrame]))
      aStorage->frame_3D_mapping[mCenterOfMass->mFrame] = shared_ptr< frame_3D<double> >(new frame_3D<double>((*mCenterOfMass->mFrame)),scoped_deleter());
    else
      (*(aStorage->frame_3D_mapping[mCenterOfMass->mFrame])) = (*mCenterOfMass->mFrame);
  };
};

void inertia_3D::doForce(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mCenterOfMass) || (!mCenterOfMass->mFrame))
    return;

  frame_3D<double> global_frame = mCenterOfMass->mFrame->getGlobalFrame(); //get the inertial frame

  mCenterOfMass->mFrame->Force -= mMass * (invert(global_frame.Quat) * global_frame.Acceleration);

  mCenterOfMass->mFrame->Torque -= mInertiaTensor * global_frame.AngAcceleration + global_frame.AngVelocity % (mInertiaTensor * global_frame.AngVelocity);

};


void inertia_3D::clearForce() {
  if((mCenterOfMass) && (mCenterOfMass->mFrame)) {
    mCenterOfMass->mFrame->Force = vect<double,3>();
    mCenterOfMass->mFrame->Torque = vect<double,3>();
  };
};



};


};






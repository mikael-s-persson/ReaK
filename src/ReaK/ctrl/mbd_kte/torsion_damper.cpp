
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

#include "torsion_damper.hpp"

namespace ReaK {

namespace kte {




void torsion_damper_2D::doMotion(kte_pass_flag aFlag, const shared_pointer<frame_storage>::type& aStorage) {
  if((!mAnchor1) || (!mAnchor2))
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_2D_mapping[mAnchor1]))
      aStorage->frame_2D_mapping[mAnchor1] = shared_pointer< frame_2D<double> >::type(new frame_2D<double>((*mAnchor1)),scoped_deleter());
    else
      (*(aStorage->frame_2D_mapping[mAnchor1])) = (*mAnchor1);
    if(!(aStorage->frame_2D_mapping[mAnchor2]))
      aStorage->frame_2D_mapping[mAnchor2] = shared_pointer< frame_2D<double> >::type(new frame_2D<double>((*mAnchor2)),scoped_deleter());
    else
      (*(aStorage->frame_2D_mapping[mAnchor2])) = (*mAnchor2);
  };
};

void torsion_damper_2D::doForce(kte_pass_flag aFlag, const shared_pointer<frame_storage>::type& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mAnchor1) || (!mAnchor2))
    return;


  double torque_mag = (mAnchor1->AngVelocity - mAnchor2->AngVelocity) * mDamping;
  mAnchor1->Torque -= torque_mag;
  mAnchor2->Torque += torque_mag;

};


void torsion_damper_2D::clearForce() {
  if(mAnchor1) {
    mAnchor1->Force = vect<double,2>();
    mAnchor1->Torque = 0.0;
  };
  if(mAnchor2) {
    mAnchor2->Force = vect<double,2>();
    mAnchor2->Torque = 0.0;
  };
};






void torsion_damper_3D::doMotion(kte_pass_flag aFlag, const shared_pointer<frame_storage>::type& aStorage) {
  if((!mAnchor1) || (!mAnchor2))
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_3D_mapping[mAnchor1]))
      aStorage->frame_3D_mapping[mAnchor1] = shared_pointer< frame_3D<double> >::type(new frame_3D<double>((*mAnchor1)),scoped_deleter());
    else
      (*(aStorage->frame_3D_mapping[mAnchor1])) = (*mAnchor1);
    if(!(aStorage->frame_3D_mapping[mAnchor2]))
      aStorage->frame_3D_mapping[mAnchor2] = shared_pointer< frame_3D<double> >::type(new frame_3D<double>((*mAnchor2)),scoped_deleter());
    else
      (*(aStorage->frame_3D_mapping[mAnchor2])) = (*mAnchor2);
  };
};

void torsion_damper_3D::doForce(kte_pass_flag aFlag, const shared_pointer<frame_storage>::type& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mAnchor1) || (!mAnchor2))
    return;

  rot_mat_3D<double> R1 = mAnchor1->Quat.getRotMat();
  rot_mat_3D<double> R2 = mAnchor2->Quat.getRotMat();
  vect<double,3> diff = (R1 * mAnchor1->AngVelocity - R2 * mAnchor2->AngVelocity) * mDamping;

  mAnchor1->Torque -= diff * R1;
  mAnchor2->Torque += diff * R2;

};


void torsion_damper_3D::clearForce() {
  if(mAnchor1) {
    mAnchor1->Force = vect<double,3>();
    mAnchor1->Torque = vect<double,3>();
  };
  if(mAnchor2) {
    mAnchor2->Force = vect<double,3>();
    mAnchor2->Torque = vect<double,3>();
  };
};




};

};








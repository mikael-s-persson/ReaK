
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

#include <ReaK/ctrl/mbd_kte/damper.hpp>

namespace ReaK {

namespace kte {



void damper_gen::doMotion(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) {
  if((!mAnchor1) || (!mAnchor2))
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->gen_coord_mapping[mAnchor1]))
      aStorage->gen_coord_mapping[mAnchor1] = shared_ptr< gen_coord<double> >(new gen_coord<double>((*mAnchor1)),scoped_deleter());
    else
      (*(aStorage->gen_coord_mapping[mAnchor1])) = (*mAnchor1);
    if(!(aStorage->gen_coord_mapping[mAnchor2]))
      aStorage->gen_coord_mapping[mAnchor2] = shared_ptr< gen_coord<double> >(new gen_coord<double>((*mAnchor2)),scoped_deleter());
    else
      (*(aStorage->gen_coord_mapping[mAnchor2])) = (*mAnchor2);
  };
};

void damper_gen::doForce(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mAnchor1) || (!mAnchor2))
    return;

  double force_mag = (mAnchor1->q_dot - mAnchor2->q_dot) * mDamping; //Positive if tension, negative if compression.

  mAnchor1->f -= force_mag;
  mAnchor2->f += force_mag;

};


void damper_gen::clearForce() {
  if(mAnchor1) {
    mAnchor1->f = 0.0;
  };
  if(mAnchor2) {
    mAnchor2->f = 0.0;
  };
};




void damper_2D::doMotion(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) {
  if((!mAnchor1) || (!mAnchor2))
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_2D_mapping[mAnchor1]))
      aStorage->frame_2D_mapping[mAnchor1] = shared_ptr< frame_2D<double> >(new frame_2D<double>((*mAnchor1)),scoped_deleter());
    else
      (*(aStorage->frame_2D_mapping[mAnchor1])) = (*mAnchor1);
    if(!(aStorage->frame_2D_mapping[mAnchor2]))
      aStorage->frame_2D_mapping[mAnchor2] = shared_ptr< frame_2D<double> >(new frame_2D<double>((*mAnchor2)),scoped_deleter());
    else
      (*(aStorage->frame_2D_mapping[mAnchor2])) = (*mAnchor2);
  };
};

void damper_2D::doForce(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mAnchor1) || (!mAnchor2))
    return;


  vect<double,2> diff = mAnchor1->Position - mAnchor2->Position;
  double force_sqrmag = norm_2_sqr(diff);
  if(force_sqrmag > 1E-7) {
    diff *= ((mAnchor1->Velocity - mAnchor2->Velocity) * diff) * mDamping / force_sqrmag;

    mAnchor1->Force -= diff * mAnchor1->Rotation;
    mAnchor2->Force += diff * mAnchor2->Rotation;
  };

};


void damper_2D::clearForce() {
  if(mAnchor1) {
    mAnchor1->Force = vect<double,2>();
    mAnchor1->Torque = 0.0;
  };
  if(mAnchor2) {
    mAnchor2->Force = vect<double,2>();
    mAnchor2->Torque = 0.0;
  };
};





void damper_3D::doMotion(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) {
  if((!mAnchor1) || (!mAnchor2))
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_3D_mapping[mAnchor1]))
      aStorage->frame_3D_mapping[mAnchor1] = shared_ptr< frame_3D<double> >(new frame_3D<double>((*mAnchor1)),scoped_deleter());
    else
      (*(aStorage->frame_3D_mapping[mAnchor1])) = (*mAnchor1);
    if(!(aStorage->frame_3D_mapping[mAnchor2]))
      aStorage->frame_3D_mapping[mAnchor2] = shared_ptr< frame_3D<double> >(new frame_3D<double>((*mAnchor2)),scoped_deleter());
    else
      (*(aStorage->frame_3D_mapping[mAnchor2])) = (*mAnchor2);
  };
};

void damper_3D::doForce(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mAnchor1) || (!mAnchor2))
    return;

  vect<double,3> diff = mAnchor1->Position - mAnchor2->Position;
  double force_sqrmag = norm_2_sqr(diff);
  if(force_sqrmag > 1E-7) {
    diff *= ((mAnchor1->Velocity - mAnchor2->Velocity) * diff) * mDamping / force_sqrmag;

    mAnchor1->Force -= invert(mAnchor1->Quat) * diff;
    mAnchor2->Force += invert(mAnchor2->Quat) * diff;
  };

};


void damper_3D::clearForce() {
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






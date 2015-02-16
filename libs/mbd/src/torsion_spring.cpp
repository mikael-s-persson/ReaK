
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

#include <ReaK/mbd/kte/torsion_spring.hpp>

namespace ReaK {

namespace kte {





void torsion_spring_2D::doMotion(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) {
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

void torsion_spring_2D::doForce(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mAnchor1) || (!mAnchor2))
    return;
  
  using std::fabs;


  double angle_diff = (invert(mAnchor1->Rotation) * mAnchor2->Rotation).getAngle() * mStiffness;
  if((mSaturation > 0) && (fabs(angle_diff) > mSaturation)) {
    if(angle_diff > 0) {
      mAnchor1->Torque += mSaturation;
      mAnchor2->Torque -= mSaturation;
    } else {
      mAnchor1->Torque -= mSaturation;
      mAnchor2->Torque += mSaturation;
    };
  } else {
    mAnchor1->Torque += angle_diff;
    mAnchor2->Torque -= angle_diff;
  };

};


void torsion_spring_2D::clearForce() {
  if(mAnchor1) {
    mAnchor1->Force = vect<double,2>();
    mAnchor1->Torque = 0.0;
  };
  if(mAnchor2) {
    mAnchor2->Force = vect<double,2>();
    mAnchor2->Torque = 0.0;
  };
};






void torsion_spring_3D::doMotion(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) {
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

void torsion_spring_3D::doForce(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mAnchor1) || (!mAnchor2))
    return;
  
  using std::fabs;


  axis_angle<double> angle_diff(invert(mAnchor1->Quat) * mAnchor2->Quat);

  double force_mag = mStiffness * angle_diff.angle();
  if((mSaturation > 0) && (fabs(force_mag) > mSaturation)) {
    if(force_mag > 0) {
      mAnchor1->Torque += mSaturation * angle_diff.axis();
      mAnchor2->Torque -= mSaturation * angle_diff.axis();
    } else {
      mAnchor1->Torque -= mSaturation * angle_diff.axis();
      mAnchor2->Torque += mSaturation * angle_diff.axis();
    };
  } else {
    mAnchor1->Torque += mStiffness * angle_diff.angle() * angle_diff.axis();
    mAnchor2->Torque -= mStiffness * angle_diff.angle() * angle_diff.axis();
  };

};


void torsion_spring_3D::clearForce() {
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









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


#include <ReaK/ctrl/mbd_kte/flexible_beam.hpp>

namespace ReaK {

namespace kte {


    
void flexible_beam_2D::doMotion(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) {
  if((!mAnchor1) || (!mAnchor2))
    return;
  
  if(mObjectFrame) {
    mObjectFrame->Position = (mAnchor1->Position + mAnchor2->Position) * 0.5;
    mObjectFrame->Velocity = (mAnchor1->Velocity + mAnchor2->Velocity) * 0.5;
    mObjectFrame->Acceleration = (mAnchor1->Acceleration + mAnchor2->Acceleration) * 0.5;
    
    mObjectFrame->Rotation = rot_mat_2D<double>((mAnchor1->Rotation.getAngle() + mAnchor2->Rotation.getAngle()) * 0.5);
    mObjectFrame->AngVelocity = (mAnchor1->AngVelocity + mAnchor2->AngVelocity) * 0.5;
    mObjectFrame->AngAcceleration = (mAnchor1->AngAcceleration + mAnchor2->AngAcceleration) * 0.5;
  };
  
  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_2D_mapping[mAnchor1]))
      aStorage->frame_2D_mapping[mAnchor1] = shared_ptr< frame_2D<double> >(new frame_2D<double>((*mAnchor1)),scoped_deleter());
    else
      (*(aStorage->frame_2D_mapping[mAnchor1])) = (*mAnchor1);
    if(!(aStorage->frame_2D_mapping[mAnchor2]))
      aStorage->frame_2D_mapping[mAnchor2] = shared_ptr< frame_2D<double> >(new frame_2D<double>((*mAnchor2)),scoped_deleter());
    else
      (*(aStorage->frame_2D_mapping[mAnchor2])) = (*mAnchor2);
    if(mObjectFrame) {
      if(!(aStorage->frame_2D_mapping[mObjectFrame]))
        aStorage->frame_2D_mapping[mObjectFrame] = shared_ptr< frame_2D<double> >(new frame_2D<double>((*mObjectFrame)),scoped_deleter());
      else
        (*(aStorage->frame_2D_mapping[mObjectFrame])) = (*mObjectFrame);
    };
  };
};
    
void flexible_beam_2D::doForce(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) {
  if((!mAnchor1) || (!mAnchor2))
    return;
  
  if(mObjectFrame) {
    pose_2D<double> p1 = mObjectFrame->getPoseRelativeTo(mAnchor1);
    pose_2D<double> p2 = mObjectFrame->getPoseRelativeTo(mAnchor2);
    
    vect<double,2> tmp_force = p1.Rotation * mObjectFrame->Force;
    mAnchor1->Force += tmp_force*0.5;
    mAnchor1->Torque += (mObjectFrame->Torque + p1.Position % tmp_force)*0.5;
        
    tmp_force = p2.Rotation * mObjectFrame->Force;
    mAnchor2->Force += tmp_force*0.5;
    mAnchor2->Torque += (mObjectFrame->Torque + p2.Position % tmp_force)*0.5;
  };
  
  vect<double,2> diff = mAnchor1->Position - mAnchor2->Position;
  vect<double,2> diff_a1 = -diff * mAnchor1->Rotation - vect<double,2>(mRestLength,0.0);
  vect<double,2> diff_a2 = diff * mAnchor2->Rotation + vect<double,2>(mRestLength,0.0);
  double angle_diff = (invert(mAnchor1->Rotation) * mAnchor2->Rotation).getAngle() * mTorsionStiffness;
  
  mAnchor1->Force += mStiffness * diff_a1;
  mAnchor1->Torque += angle_diff;
  mAnchor2->Force += mStiffness * diff_a2;
  mAnchor2->Torque -= angle_diff;
  
  if((aFlag == store_dynamics) && (aStorage)) {
    if(aStorage->frame_2D_mapping[mObjectFrame]) {
      aStorage->frame_2D_mapping[mObjectFrame]->Force = mObjectFrame->Force;
      aStorage->frame_2D_mapping[mObjectFrame]->Torque = mObjectFrame->Torque;
    };
  };  
};


void flexible_beam_2D::clearForce() {
  if(mAnchor1) {
    mAnchor1->Force = vect<double,2>();
    mAnchor1->Torque = 0.0;
  };
  if(mAnchor2) {
    mAnchor2->Force = vect<double,2>();
    mAnchor2->Torque = 0.0;
  };
  if(mObjectFrame) {
    mObjectFrame->Force = vect<double,2>();
    mObjectFrame->Torque = 0.0;
  };
};
    
   

void flexible_beam_3D::doMotion(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) {
  if((!mAnchor1) || (!mAnchor2))
    return;
  
  if(mObjectFrame) {
    mObjectFrame->Position = (mAnchor1->Position + mAnchor2->Position) * 0.5;
    mObjectFrame->Velocity = (mAnchor1->Velocity + mAnchor2->Velocity) * 0.5;
    mObjectFrame->Acceleration = (mAnchor1->Acceleration + mAnchor2->Acceleration) * 0.5;
    
    //mObjectFrame->Rotation = rot_mat_3D<double>((mAnchor1->Rotation.getAngle() + mAnchor2->Rotation.getAngle()) * 0.5);
    axis_angle<double> angle_diff(invert(mAnchor1->Quat) * mAnchor2->Quat);
    angle_diff.angle() *= 0.5;
    mObjectFrame->Quat = mAnchor1->Quat * angle_diff.getQuaternion();
    rot_mat_3D<double> R_o = mObjectFrame->Quat.getRotMat();
    rot_mat_3D<double> R_a1 = mAnchor1->Quat.getRotMat();
    rot_mat_3D<double> R_a2 = mAnchor2->Quat.getRotMat();
    mObjectFrame->AngVelocity = ((R_a1 * mAnchor1->AngVelocity + R_a2 * mAnchor2->AngVelocity) * 0.5) * R_o;
    mObjectFrame->AngAcceleration = ((R_a1 * mAnchor1->AngAcceleration + R_a2 * mAnchor2->AngAcceleration) * 0.5) * R_o;
  };
  
  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_3D_mapping[mAnchor1]))
      aStorage->frame_3D_mapping[mAnchor1] = shared_ptr< frame_3D<double> >(new frame_3D<double>((*mAnchor1)),scoped_deleter());
    else
      (*(aStorage->frame_3D_mapping[mAnchor1])) = (*mAnchor1);
    if(!(aStorage->frame_3D_mapping[mAnchor2]))
      aStorage->frame_3D_mapping[mAnchor2] = shared_ptr< frame_3D<double> >(new frame_3D<double>((*mAnchor2)),scoped_deleter());
    else
      (*(aStorage->frame_3D_mapping[mAnchor2])) = (*mAnchor2);
    if(mObjectFrame) {
      if(!(aStorage->frame_3D_mapping[mObjectFrame]))
        aStorage->frame_3D_mapping[mObjectFrame] = shared_ptr< frame_3D<double> >(new frame_3D<double>((*mObjectFrame)),scoped_deleter());
      else
        (*(aStorage->frame_3D_mapping[mObjectFrame])) = (*mObjectFrame);
    };
  };
};
    
void flexible_beam_3D::doForce(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) {
  if((!mAnchor1) || (!mAnchor2))
    return;
  
  if(mObjectFrame) {
    pose_3D<double> p1 = mObjectFrame->getPoseRelativeTo(mAnchor1);
    pose_3D<double> p2 = mObjectFrame->getPoseRelativeTo(mAnchor2);
    
    rot_mat_3D<double> R(p1.Quat.getRotMat());
    vect<double,3> tmp_force = R * mObjectFrame->Force;
    mAnchor1->Force += tmp_force * 0.5;
    mAnchor1->Torque += (R * mObjectFrame->Torque + p1.Position % tmp_force) * 0.5;
    
    R = p2.Quat.getRotMat();
    tmp_force = R * mObjectFrame->Force;
    mAnchor2->Force += tmp_force * 0.5;
    mAnchor2->Torque += (R * mObjectFrame->Torque + p2.Position % tmp_force) * 0.5;
  };
  
  vect<double,3> diff = mAnchor1->Position - mAnchor2->Position;
  vect<double,3> diff_a1 = invert(mAnchor1->Quat) * (-diff) - vect<double,3>(mRestLength,0.0,0.0);
  vect<double,3> diff_a2 = invert(mAnchor2->Quat) * diff + vect<double,3>(mRestLength,0.0,0.0);
  
  axis_angle<double> angle_diff(invert(mAnchor1->Quat) * mAnchor2->Quat);
  
  mAnchor1->Force += mStiffness * diff_a1;
  mAnchor1->Torque += mTorsionStiffness * angle_diff.angle() * angle_diff.axis();
  mAnchor2->Force += mStiffness * diff_a2;
  mAnchor2->Torque -= mTorsionStiffness * angle_diff.angle() * angle_diff.axis();
  
  if((aFlag == store_dynamics) && (aStorage)) {
    if(aStorage->frame_3D_mapping[mObjectFrame]) {
      aStorage->frame_3D_mapping[mObjectFrame]->Force = mObjectFrame->Force;
      aStorage->frame_3D_mapping[mObjectFrame]->Torque = mObjectFrame->Torque;
    };
  }; 
};


void flexible_beam_3D::clearForce() {
  if(mAnchor1) {
    mAnchor1->Force = vect<double,3>();
    mAnchor1->Torque = vect<double,3>();
  };
  if(mAnchor2) {
    mAnchor2->Force = vect<double,3>();
    mAnchor2->Torque = vect<double,3>();
  };
  if(mObjectFrame) {
    mObjectFrame->Force = vect<double,3>();
    mObjectFrame->Torque = vect<double,3>();
  };
  
};
    
   


};

};






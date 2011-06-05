
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

#include "prismatic_joint.hpp"

namespace ReaK {

namespace kte {

  

    
void prismatic_joint_2D::doMotion(kte_pass_flag aFlag, boost::shared_ptr<frame_storage> aStorage) {
  if((!mEnd) || (!mBase))
    return;
  
  mEnd->Parent = mBase->Parent;
  
  if(!mCoord) {
    mEnd->Position = mBase->Position;
    mEnd->Velocity = mBase->Velocity;
    mEnd->Acceleration = mBase->Acceleration;
  } else {
    vect<double,2> tmp_pos(mCoord->q * mAxis);
    vect<double,2> tmp_vel(mCoord->q_dot * mAxis);
    mEnd->Position = mBase->Position + mBase->Rotation * tmp_pos;
    mEnd->Velocity = mBase->Velocity + mBase->Rotation * ( (mBase->AngVelocity % tmp_pos) + tmp_vel );
    mEnd->Acceleration = mBase->Acceleration + mBase->Rotation * ( (-mBase->AngVelocity * mBase->AngVelocity) * tmp_pos + (2.0 * mBase->AngVelocity) % tmp_vel + mBase->AngAcceleration % tmp_pos + (mCoord->q_ddot * mAxis) );
 
    if(mJacobian) {
      mJacobian->Parent = mEnd;
      
      mJacobian->qd_vel = mAxis;
      mJacobian->qd_avel = 0.0;
      mJacobian->qd_acc = vect<double,2>();
      mJacobian->qd_aacc = 0.0;
    };
  };
  
  mEnd->Rotation = mBase->Rotation;
  mEnd->AngVelocity = mBase->AngVelocity;
  mEnd->AngAcceleration = mBase->AngAcceleration;
  
  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_2D_mapping[mBase]))
      aStorage->frame_2D_mapping[mBase] = boost::shared_ptr< frame_2D<double> >(new frame_2D<double>((*mBase)),scoped_deleter());
    else
      (*(aStorage->frame_2D_mapping[mBase])) = (*mBase);
    if(!(aStorage->frame_2D_mapping[mEnd]))
      aStorage->frame_2D_mapping[mEnd] = boost::shared_ptr< frame_2D<double> >(new frame_2D<double>((*mEnd)),scoped_deleter());
    else
      (*(aStorage->frame_2D_mapping[mEnd])) = (*mEnd);
    if(mCoord) {
      if(!(aStorage->gen_coord_mapping[mCoord]))
        aStorage->gen_coord_mapping[mCoord] = boost::shared_ptr< gen_coord<double> >(new gen_coord<double>((*mCoord)),scoped_deleter());
      else
        (*(aStorage->gen_coord_mapping[mCoord])) = (*mCoord);
    };
  };
};
    
void prismatic_joint_2D::doForce(kte_pass_flag aFlag, boost::shared_ptr<frame_storage> aStorage) {
  if((!mEnd) || (!mBase))
    return;
  
  if(!mCoord) {
    mBase->Force += mEnd->Force;
  } else {
    double tmp_f = mEnd->Force * mAxis;
    mCoord->f += tmp_f;
    mBase->Force += mEnd->Force - tmp_f * mAxis;
  };
  
  mBase->Torque += mEnd->Torque + (mCoord->q * mAxis) % mEnd->Force;
  
  if((aFlag == store_dynamics) && (aStorage)) {
    if(aStorage->frame_2D_mapping[mEnd]) {
      aStorage->frame_2D_mapping[mEnd]->Force = mEnd->Force;
      aStorage->frame_2D_mapping[mEnd]->Torque = mEnd->Torque;
    };
  };
};


void prismatic_joint_2D::clearForce() {
  if(mEnd) {
    mEnd->Force = vect<double,2>();
    mEnd->Torque = 0.0;
  };
  if(mBase) {
    mBase->Force = vect<double,2>();
    mBase->Torque = 0.0;
  };
  if(mCoord) {
    mCoord->f = 0.0;
  };  
};


void prismatic_joint_2D::applyReactionForce(double aForce) {
  if(mCoord)
    mBase->Force -= aForce * mAxis;
};
    




void prismatic_joint_3D::doMotion(kte_pass_flag aFlag, boost::shared_ptr<frame_storage> aStorage) {
  if((!mEnd) || (!mBase))
    return;
  
  mEnd->Parent = mBase->Parent;
  
  if(!mCoord) {
    mEnd->Position = mBase->Position;
    mEnd->Velocity = mBase->Velocity;
    mEnd->Acceleration = mBase->Acceleration;
  } else {
    rot_mat_3D<double> R(mBase->Quat.getRotMat());
    vect<double,3> tmp_pos = mCoord->q * mAxis;
    vect<double,3> tmp_vel = mCoord->q_dot * mAxis;
    mEnd->Position = mBase->Position + R * tmp_pos;
    mEnd->Velocity = mBase->Velocity + R * ( (mBase->AngVelocity % tmp_pos) + tmp_vel );
    mEnd->Acceleration = mBase->Acceleration + R * ( (mBase->AngVelocity % (mBase->AngVelocity % tmp_pos)) + 2.0 * (mBase->AngVelocity % tmp_vel) + (mBase->AngAcceleration % tmp_pos) + (mCoord->q_ddot * mAxis) );
    
    if(mJacobian) {
      mJacobian->Parent = mEnd;
      
      mJacobian->qd_vel = mAxis;
      mJacobian->qd_avel = vect<double,3>();
      mJacobian->qd_acc = vect<double,3>();
      mJacobian->qd_aacc = vect<double,3>();
    };
  };
  
  mEnd->Quat = mBase->Quat;
  mEnd->AngVelocity = mBase->AngVelocity;
  mEnd->AngAcceleration = mBase->AngAcceleration;

  mEnd->UpdateQuatDot();
  
  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_3D_mapping[mBase]))
      aStorage->frame_3D_mapping[mBase] = boost::shared_ptr< frame_3D<double> >(new frame_3D<double>((*mBase)),scoped_deleter());
    else
      (*(aStorage->frame_3D_mapping[mBase])) = (*mBase);
    if(!(aStorage->frame_3D_mapping[mEnd]))
      aStorage->frame_3D_mapping[mEnd] = boost::shared_ptr< frame_3D<double> >(new frame_3D<double>((*mEnd)),scoped_deleter());
    else
      (*(aStorage->frame_3D_mapping[mEnd])) = (*mEnd);
    if(mCoord) {
      if(!(aStorage->gen_coord_mapping[mCoord]))
        aStorage->gen_coord_mapping[mCoord] = boost::shared_ptr< gen_coord<double> >(new gen_coord<double>((*mCoord)),scoped_deleter());
      else
        (*(aStorage->gen_coord_mapping[mCoord])) = (*mCoord);
    };
  };
};
    
void prismatic_joint_3D::doForce(kte_pass_flag aFlag, boost::shared_ptr<frame_storage> aStorage) {
  if((!mEnd) || (!mBase))
    return;
    
  if(!mCoord) {
    mBase->Force += mEnd->Force;
  } else {
    double tmp_f = mEnd->Force * mAxis;
    mCoord->f += tmp_f;
    mBase->Force += mEnd->Force - tmp_f * mAxis;
  };
  
  mBase->Torque += mEnd->Torque + (mCoord->q * mAxis) % mEnd->Force;
  
  if((aFlag == store_dynamics) && (aStorage)) {
    if(aStorage->frame_3D_mapping[mEnd]) {
      aStorage->frame_3D_mapping[mEnd]->Force = mEnd->Force;
      aStorage->frame_3D_mapping[mEnd]->Torque = mEnd->Torque;
    };
  };
};


void prismatic_joint_3D::clearForce() {
  if(mEnd) {
    mEnd->Force = vect<double,3>();
    mEnd->Torque = vect<double,3>();
  };
  if(mBase) {
    mBase->Force = vect<double,3>();
    mBase->Torque = vect<double,3>();
  };
  if(mCoord) {
    mCoord->f = 0.0;
  };  
};


void prismatic_joint_3D::applyReactionForce(double aForce) {
  if(mCoord)
    mBase->Force -= aForce * mAxis;
};
    


};

};






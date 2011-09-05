
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

#include "dry_revolute_joint.hpp"

namespace ReaK {

namespace kte {


    
void dry_revolute_joint_2D::doForce(kte_pass_flag aFlag, const shared_pointer<frame_storage>::type& aStorage) {
  if((!mEnd) || (!mBase))
    return;
  
  using std::fabs;
  
  if(!mAngle) {
    mBase->Force += mEnd->Force;
    mBase->Torque += mEnd->Torque;
  } else {
    vect<double,2> tmp_f = rot_mat_2D<double>(mAngle->q) * mEnd->Force;
    mBase->Force += tmp_f;
    if(fabs(mAngle->q_dot) > mSlipVelocity) 
      mAngle->f += mEnd->Torque - mAngle->q_dot / fabs(mAngle->q_dot) * mSlipCoef * norm(tmp_f);
    else {
      mAngle->f += mEnd->Torque - mAngle->q_dot / mSlipVelocity * mStictionCoef * norm(tmp_f);
    };
  };
  
  if((aFlag == store_dynamics) && (aStorage)) {
    if(aStorage->frame_2D_mapping[mEnd]) {
      aStorage->frame_2D_mapping[mEnd]->Force = mEnd->Force;
      aStorage->frame_2D_mapping[mEnd]->Torque = mEnd->Torque;
    };
  };  
};
  



    
void dry_revolute_joint_3D::doForce(kte_pass_flag aFlag, const shared_pointer<frame_storage>::type& aStorage) {
  if((!mEnd) || (!mBase))
    return;
  
  using std::fabs;
  
  if(!mAngle) {
    mBase->Force += mEnd->Force;
    mBase->Torque += mEnd->Torque;
  } else {
    rot_mat_3D<double> R(axis_angle<double>(mAngle->q,mAxis).getRotMat());
    vect<double,3> tmp_f = R * mEnd->Force;
    mBase->Force += tmp_f;
    double tmp_t = mEnd->Torque * mAxis;
    if(fabs(mAngle->q_dot) > mSlipVelocity) 
      mAngle->f += tmp_t - mAngle->q_dot / fabs(mAngle->q_dot) * mSlipCoef * norm(tmp_f);
    else {
      mAngle->f += tmp_t - mAngle->q_dot / mSlipVelocity * mStictionCoef * norm(tmp_f);
    };
    mBase->Torque += R * ( mEnd->Torque - tmp_t * mAxis );
  };
  
  if((aFlag == store_dynamics) && (aStorage)) {
    if(aStorage->frame_3D_mapping[mEnd]) {
      aStorage->frame_3D_mapping[mEnd]->Force = mEnd->Force;
      aStorage->frame_3D_mapping[mEnd]->Torque = mEnd->Torque;
    };
  };
};



};


};





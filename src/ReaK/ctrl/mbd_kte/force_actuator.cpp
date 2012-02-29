
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

#include "force_actuator.hpp"

namespace ReaK {

namespace kte {


void force_actuator_gen::doMotion(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mFrame) || (!mJoint))
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->gen_coord_mapping[mFrame]))
      aStorage->gen_coord_mapping[mFrame] = shared_ptr< gen_coord<double> >(new gen_coord<double>((*mFrame)),scoped_deleter());
    else
      (*(aStorage->gen_coord_mapping[mFrame])) = (*mFrame);
  };
};

void force_actuator_gen::doForce(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mFrame) || (!mJoint))
    return;

  mJoint->applyReactionForce(-(mFrame->f));

};


void force_actuator_gen::clearForce() {
  if(mFrame) {
    mFrame->f = 0.0;
  };
};









void force_actuator_2D::doMotion(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mFrame) || (!mJoint))
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_2D_mapping[mFrame]))
      aStorage->frame_2D_mapping[mFrame] = shared_ptr< frame_2D<double> >(new frame_2D<double>((*mFrame)),scoped_deleter());
    else
      (*(aStorage->frame_2D_mapping[mFrame])) = (*mFrame);
  };
};

void force_actuator_2D::doForce(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mFrame) || (!mJoint))
    return;

  mJoint->applyReactionForce(-mFrame->Force, -mFrame->Torque);

};


void force_actuator_2D::clearForce() {
  if(mFrame) {
    mFrame->Force = vect<double,2>();
    mFrame->Torque = 0.0;
  };
};










void force_actuator_3D::doMotion(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mFrame) || (!mJoint))
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_3D_mapping[mFrame]))
      aStorage->frame_3D_mapping[mFrame] = shared_ptr< frame_3D<double> >(new frame_3D<double>((*mFrame)),scoped_deleter());
    else
      (*(aStorage->frame_3D_mapping[mFrame])) = (*mFrame);
  };
};

void force_actuator_3D::doForce(kte_pass_flag aFlag, const shared_ptr<frame_storage>& aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mFrame) || (!mJoint))
    return;

  mJoint->applyReactionForce(-mFrame->Force, -mFrame->Torque);

};


void force_actuator_3D::clearForce() {
  if(mFrame) {
    mFrame->Force = vect<double,3>();
    mFrame->Torque = vect<double,3>();
  };
};






};


};






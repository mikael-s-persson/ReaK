
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

#include "driving_actuator.hpp"

namespace ReaK {

namespace kte {


void driving_actuator_gen::doForce(kte_pass_flag aFlag, boost::shared_ptr<frame_storage> aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mFrame) || (!mJoint))
    return;

  mFrame->f += mDriveForce;
  mJoint->applyReactionForce(mDriveForce);

};




void driving_actuator_2D::doForce(kte_pass_flag aFlag, boost::shared_ptr<frame_storage> aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mFrame) || (!mJoint))
    return;

  mFrame->Force += mDriveForce;
  mFrame->Torque += mDriveTorque;
  mJoint->applyReactionForce(mDriveForce, mDriveTorque);

};



void driving_actuator_3D::doForce(kte_pass_flag aFlag, boost::shared_ptr<frame_storage> aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mFrame) || (!mJoint))
    return;

  mFrame->Force += mDriveForce;
  mFrame->Torque += mDriveTorque;
  mJoint->applyReactionForce(mDriveForce, mDriveTorque);

};




};

};







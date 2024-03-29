
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

#include "ReaK/mbd/kte/line_point_mindist.h"

namespace ReaK::kte {

void line_point_mindist_2D::doMotion(
    kte_pass_flag aFlag, const std::shared_ptr<frame_storage>& aStorage) {
  if ((!mBase) || (!mEnd)) {
    return;
  }

  mEnd->Parent = mBase->Parent;
  mEnd->Position = (mBase->Position * mTangent) * mTangent - mOriginMinDist;
  mEnd->Velocity = (mBase->Velocity * mTangent) * mTangent;
  mEnd->Acceleration = (mBase->Acceleration * mTangent) * mTangent;
  mEnd->Rotation = mBase->Rotation;
  mEnd->AngVelocity = mBase->AngVelocity;
  mEnd->AngAcceleration = mBase->AngAcceleration;

  if ((aFlag == store_kinematics) && (aStorage)) {
    if (!(aStorage->frame_2D_mapping[mBase])) {
      aStorage->frame_2D_mapping[mBase] =
          std::make_shared<frame_2D<double>>(*mBase);
    } else {
      (*(aStorage->frame_2D_mapping[mBase])) = (*mBase);
    }
    if (!(aStorage->frame_2D_mapping[mEnd])) {
      aStorage->frame_2D_mapping[mEnd] =
          std::make_shared<frame_2D<double>>(*mEnd);
    } else {
      (*(aStorage->frame_2D_mapping[mEnd])) = (*mEnd);
    }
  }
}

void line_point_mindist_2D::doForce(
    kte_pass_flag aFlag, const std::shared_ptr<frame_storage>& aStorage) {
  if ((!mBase) || (!mEnd)) {
    return;
  }

  if ((aFlag == store_dynamics) && (aStorage)) {
    if (aStorage->frame_2D_mapping[mEnd]) {
      aStorage->frame_2D_mapping[mEnd]->Force = mEnd->Force;
      aStorage->frame_2D_mapping[mEnd]->Torque = mEnd->Torque;
    }
  }
}

void line_point_mindist_2D::clearForce() {
  if (mEnd) {
    mEnd->Force = vect<double, 2>();
    mEnd->Torque = 0.0;
  }
  if (mBase) {
    mBase->Force = vect<double, 2>();
    mBase->Torque = 0.0;
  }
}

void line_point_mindist_3D::doMotion(
    kte_pass_flag aFlag, const std::shared_ptr<frame_storage>& aStorage) {
  if ((!mBase) || (!mEnd)) {
    return;
  }

  mEnd->Parent = mBase->Parent;
  mEnd->Position = (mBase->Position * mTangent) * mTangent - mOriginMinDist;
  mEnd->Velocity = (mBase->Velocity * mTangent) * mTangent;
  mEnd->Acceleration = (mBase->Acceleration * mTangent) * mTangent;
  mEnd->Quat = mBase->Quat;
  mEnd->AngVelocity = mBase->AngVelocity;
  mEnd->AngAcceleration = mBase->AngAcceleration;

  if ((aFlag == store_kinematics) && (aStorage)) {
    if (!(aStorage->frame_3D_mapping[mBase])) {
      aStorage->frame_3D_mapping[mBase] =
          std::make_shared<frame_3D<double>>(*mBase);
    } else {
      (*(aStorage->frame_3D_mapping[mBase])) = (*mBase);
    }
    if (!(aStorage->frame_3D_mapping[mEnd])) {
      aStorage->frame_3D_mapping[mEnd] =
          std::make_shared<frame_3D<double>>(*mEnd);
    } else {
      (*(aStorage->frame_3D_mapping[mEnd])) = (*mEnd);
    }
  }
}

void line_point_mindist_3D::doForce(
    kte_pass_flag aFlag, const std::shared_ptr<frame_storage>& aStorage) {
  if ((!mBase) || (!mEnd)) {
    return;
  }

  if ((aFlag == store_dynamics) && (aStorage)) {
    if (aStorage->frame_3D_mapping[mEnd]) {
      aStorage->frame_3D_mapping[mEnd]->Force = mEnd->Force;
      aStorage->frame_3D_mapping[mEnd]->Torque = mEnd->Torque;
    }
  }
}

void line_point_mindist_3D::clearForce() {
  if (mEnd) {
    mEnd->Force = vect<double, 3>();
    mEnd->Torque = vect<double, 3>();
  }
  if (mBase) {
    mBase->Force = vect<double, 3>();
    mBase->Torque = vect<double, 3>();
  }
}
}  // namespace ReaK::kte


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

#include <ReaK/mbd/kte/inertial_beam.hpp>

namespace ReaK::kte {

void inertial_beam_2D::doMotion(
    kte_pass_flag aFlag, const std::shared_ptr<frame_storage>& aStorage) {
  if ((!mAnchor1) || (!mAnchor2)) {
    return;
  }

  if ((aFlag == store_kinematics) && (aStorage)) {
    if (!(aStorage->frame_2D_mapping[mAnchor1])) {
      aStorage->frame_2D_mapping[mAnchor1] =
          std::make_shared<frame_2D<double>>(*mAnchor1);
    } else {
      (*(aStorage->frame_2D_mapping[mAnchor1])) = (*mAnchor1);
    }
    if (!(aStorage->frame_2D_mapping[mAnchor2])) {
      aStorage->frame_2D_mapping[mAnchor2] =
          std::make_shared<frame_2D<double>>(*mAnchor2);
    } else {
      (*(aStorage->frame_2D_mapping[mAnchor2])) = (*mAnchor2);
    }
  }
}

void inertial_beam_2D::doForce(kte_pass_flag aFlag,
                               const std::shared_ptr<frame_storage>& aStorage) {
  RK_UNUSED(aFlag);
  RK_UNUSED(aStorage);
  if ((!mAnchor1) || (!mAnchor2)) {
    return;
  }

  vect<double, 2> diff = mAnchor1->Position - mAnchor2->Position;

  mAnchor1->Force -= 0.5 * mMass * mAnchor1->Acceleration * mAnchor1->Rotation;
  mAnchor1->Torque +=
      0.5 * mMass * mAnchor1->AngAcceleration * norm_2_sqr(diff);
  mAnchor2->Force -= 0.5 * mMass * mAnchor2->Acceleration * mAnchor2->Rotation;
  mAnchor2->Torque +=
      0.5 * mMass * mAnchor2->AngAcceleration * norm_2_sqr(diff);
}

void inertial_beam_2D::clearForce() {
  if (mAnchor1) {
    mAnchor1->Force = vect<double, 2>();
    mAnchor1->Torque = 0.0;
  }
  if (mAnchor2) {
    mAnchor2->Force = vect<double, 2>();
    mAnchor2->Torque = 0.0;
  }
}

void inertial_beam_3D::doMotion(
    kte_pass_flag aFlag, const std::shared_ptr<frame_storage>& aStorage) {
  if ((!mAnchor1) || (!mAnchor2)) {
    return;
  }

  if ((aFlag == store_kinematics) && (aStorage)) {
    if (!(aStorage->frame_3D_mapping[mAnchor1])) {
      aStorage->frame_3D_mapping[mAnchor1] =
          std::make_shared<frame_3D<double>>(*mAnchor1);
    } else {
      (*(aStorage->frame_3D_mapping[mAnchor1])) = (*mAnchor1);
    }
    if (!(aStorage->frame_3D_mapping[mAnchor2])) {
      aStorage->frame_3D_mapping[mAnchor2] =
          std::make_shared<frame_3D<double>>(*mAnchor2);
    } else {
      (*(aStorage->frame_3D_mapping[mAnchor2])) = (*mAnchor2);
    }
  }
}

void inertial_beam_3D::doForce(kte_pass_flag aFlag,
                               const std::shared_ptr<frame_storage>& aStorage) {
  RK_UNUSED(aFlag);
  RK_UNUSED(aStorage);
  if ((!mAnchor1) || (!mAnchor2)) {
    return;
  }

  vect<double, 3> diff = mAnchor1->Position - mAnchor2->Position;
  vect<double, 3> diff_a1 = invert(mAnchor1->Quat) * -diff;
  vect<double, 3> diff_a2 = invert(mAnchor2->Quat) * diff;

  mAnchor1->Force -=
      0.5 * mMass * (invert(mAnchor1->Quat) * mAnchor1->Acceleration);
  mAnchor1->Torque -=
      0.5 * mMass * diff_a1 % (mAnchor1->AngAcceleration % diff_a1);
  mAnchor2->Force -=
      0.5 * mMass * (invert(mAnchor2->Quat) * mAnchor2->Acceleration);
  mAnchor2->Torque -=
      0.5 * mMass * diff_a2 % (mAnchor2->AngAcceleration % diff_a2);
}

void inertial_beam_3D::clearForce() {
  if (mAnchor1) {
    mAnchor1->Force = vect<double, 3>();
    mAnchor1->Torque = vect<double, 3>();
  }
  if (mAnchor2) {
    mAnchor2->Force = vect<double, 3>();
    mAnchor2->Torque = vect<double, 3>();
  }
}
}  // namespace ReaK::kte

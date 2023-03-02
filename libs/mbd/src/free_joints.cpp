
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

#include <ReaK/mbd/kte/free_joints.hpp>

namespace ReaK::kte {

void free_joint_2D::doMotion(kte_pass_flag aFlag,
                             const std::shared_ptr<frame_storage>& aStorage) {
  if ((!mEnd) || (!mBase)) {
    return;
  }

  if (!mCoord) {
    *mEnd = *mBase;
  } else {
    *mEnd = (*mBase) * (*mCoord);
    mEnd->Parent = mBase->Parent;

    if (mJacobian) {
      mJacobian->Parent = mEnd;

      mJacobian->vel_vel = vect<vect<double, 2>, 2>(vect<double, 2>(1.0, 0.0),
                                                    vect<double, 2>(0.0, 1.0));
      mJacobian->vel_avel = vect<double, 2>();
      mJacobian->avel_vel = vect<double, 2>();
      mJacobian->avel_avel = 1.0;
      mJacobian->vel_acc = vect<vect<double, 2>, 2>();
      mJacobian->vel_aacc = vect<double, 2>();
      mJacobian->avel_acc = vect<double, 2>();
      mJacobian->avel_aacc = 0.0;
    }
  }

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
    if (mCoord) {
      if (!(aStorage->frame_2D_mapping[mCoord])) {
        aStorage->frame_2D_mapping[mCoord] =
            std::make_shared<frame_2D<double>>(*mCoord);
      } else {
        (*(aStorage->frame_2D_mapping[mCoord])) = (*mCoord);
      }
    }
  }
}

void free_joint_2D::doForce(kte_pass_flag aFlag,
                            const std::shared_ptr<frame_storage>& aStorage) {
  if ((!mEnd) || (!mBase)) {
    return;
  }

  if (!mCoord) {
    mBase->Force += mEnd->Force;
    mBase->Torque += mEnd->Torque;
  } else {
    mCoord->Force += mEnd->Force;
    mCoord->Torque += mEnd->Torque;
  }

  if ((aFlag == store_dynamics) && (aStorage)) {
    if (aStorage->frame_2D_mapping[mEnd]) {
      aStorage->frame_2D_mapping[mEnd]->Force = mEnd->Force;
      aStorage->frame_2D_mapping[mEnd]->Torque = mEnd->Torque;
    }
  }
}

void free_joint_2D::clearForce() {
  if (mEnd) {
    mEnd->Force = vect<double, 2>();
    mEnd->Torque = 0.0;
  }
  if (mBase) {
    mBase->Force = vect<double, 2>();
    mBase->Torque = 0.0;
  }
  if (mCoord) {
    mCoord->Force = vect<double, 2>();
    mCoord->Torque = 0.0;
  }
}

void free_joint_2D::applyReactionForce(vect<double, 2> aForce, double aTorque) {
  if (mBase) {
    mBase->Force -= aForce;
    mBase->Torque -= aTorque;
  }
}

void free_joint_3D::doMotion(kte_pass_flag aFlag,
                             const std::shared_ptr<frame_storage>& aStorage) {
  if ((!mEnd) || (!mBase)) {
    return;
  }

  if (!mCoord) {
    *mEnd = *mBase;
  } else {
    *mEnd = (*mBase) * (*mCoord);
    mEnd->Parent = mBase->Parent;
    mEnd->UpdateQuatDot();

    if (mJacobian) {
      mJacobian->Parent = mEnd;

      mJacobian->vel_vel = vect<vect<double, 3>, 3>(
          vect<double, 3>(1.0, 0.0, 0.0), vect<double, 3>(0.0, 1.0, 0.0),
          vect<double, 3>(0.0, 0.0, 1.0));
      mJacobian->vel_avel = vect<vect<double, 3>, 3>();
      mJacobian->avel_vel = vect<vect<double, 3>, 3>();
      mJacobian->avel_avel = vect<vect<double, 3>, 3>(
          vect<double, 3>(1.0, 0.0, 0.0), vect<double, 3>(0.0, 1.0, 0.0),
          vect<double, 3>(0.0, 0.0, 1.0));
      mJacobian->vel_acc = vect<vect<double, 3>, 3>();
      mJacobian->vel_aacc = vect<vect<double, 3>, 3>();
      mJacobian->avel_acc = vect<vect<double, 3>, 3>();
      mJacobian->avel_aacc = vect<vect<double, 3>, 3>();
    }
  }

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
    if (mCoord) {
      if (!(aStorage->frame_3D_mapping[mCoord])) {
        aStorage->frame_3D_mapping[mCoord] =
            std::make_shared<frame_3D<double>>(*mCoord);
      } else {
        (*(aStorage->frame_3D_mapping[mCoord])) = (*mCoord);
      }
    }
  }
}

void free_joint_3D::doForce(kte_pass_flag aFlag,
                            const std::shared_ptr<frame_storage>& aStorage) {
  if ((!mEnd) || (!mBase)) {
    return;
  }

  if (!mCoord) {
    mBase->Force += mEnd->Force;
    mBase->Torque += mEnd->Torque;
  } else {
    mCoord->Force += mEnd->Force;
    mCoord->Torque += mEnd->Torque;
  }

  if ((aFlag == store_dynamics) && (aStorage)) {
    if (aStorage->frame_3D_mapping[mEnd]) {
      aStorage->frame_3D_mapping[mEnd]->Force = mEnd->Force;
      aStorage->frame_3D_mapping[mEnd]->Torque = mEnd->Torque;
    }
  }
}

void free_joint_3D::clearForce() {
  if (mEnd) {
    mEnd->Force = vect<double, 3>();
    mEnd->Torque = vect<double, 3>();
  }
  if (mBase) {
    mBase->Force = vect<double, 3>();
    mBase->Torque = vect<double, 3>();
  }
  if (mCoord) {
    mCoord->Force = vect<double, 3>();
    mCoord->Torque = vect<double, 3>();
  }
}

void free_joint_3D::applyReactionForce(vect<double, 3> aForce,
                                       vect<double, 3> aTorque) {
  if (mBase) {
    mBase->Force -= aForce;
    mBase->Torque -= aTorque;
  }
}
}  // namespace ReaK::kte

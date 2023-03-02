
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

#include <ReaK/mbd/kte/joint_friction.hpp>

namespace ReaK::kte {

void joint_dry_microslip_gen::doMotion(
    kte_pass_flag aFlag, const std::shared_ptr<frame_storage>& aStorage) {
  if (!mAnchor) {
    return;
  }

  if ((aFlag == store_kinematics) && (aStorage)) {
    if (!(aStorage->gen_coord_mapping[mAnchor])) {
      aStorage->gen_coord_mapping[mAnchor] =
          std::make_shared<gen_coord<double>>(*mAnchor);
    } else {
      (*(aStorage->gen_coord_mapping[mAnchor])) = (*mAnchor);
    }
  }
}

void joint_dry_microslip_gen::doForce(
    kte_pass_flag aFlag, const std::shared_ptr<frame_storage>& aStorage) {
  RK_UNUSED(aFlag);
  RK_UNUSED(aStorage);
  if (!mAnchor) {
    return;
  }

  using std::abs;

  double tmp_speed = abs(mAnchor->q_dot);

  if (tmp_speed <= mStictionVelocity) {
    mAnchor->f -= mAnchor->q_dot * mStictionCoef / mStictionVelocity;
  } else if (tmp_speed < mSlipVelocity) {
    mAnchor->f -= (mAnchor->q_dot > 0.0 ? 1.0 : -1.0) *
                  (mStictionCoef + (mSlipCoef - mStictionCoef) *
                                       (tmp_speed - mStictionVelocity) /
                                       (mSlipVelocity - mStictionVelocity));
  } else {
    mAnchor->f -= (mAnchor->q_dot > 0.0 ? 1.0 : -1.0) * mSlipCoef;
  }
}

void joint_dry_microslip_gen::clearForce() {
  if (mAnchor) {
    mAnchor->f = 0.0;
  }
}

void joint_viscosity_gen::doMotion(
    kte_pass_flag aFlag, const std::shared_ptr<frame_storage>& aStorage) {
  if (!mAnchor) {
    return;
  }

  if ((aFlag == store_kinematics) && (aStorage)) {
    if (!(aStorage->gen_coord_mapping[mAnchor])) {
      aStorage->gen_coord_mapping[mAnchor] =
          std::make_shared<gen_coord<double>>(*mAnchor);
    } else {
      (*(aStorage->gen_coord_mapping[mAnchor])) = (*mAnchor);
    }
  }
}

void joint_viscosity_gen::doForce(
    kte_pass_flag aFlag, const std::shared_ptr<frame_storage>& aStorage) {
  RK_UNUSED(aFlag);
  RK_UNUSED(aStorage);
  if (!mAnchor) {
    return;
  }

  mAnchor->f -= mAnchor->q_dot * mViscosity;
}

void joint_viscosity_gen::clearForce() {
  if (mAnchor) {
    mAnchor->f = 0.0;
  }
}
}  // namespace ReaK::kte

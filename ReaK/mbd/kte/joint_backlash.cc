
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

#include "ReaK/mbd/kte/joint_backlash.h"

namespace ReaK::kte {

void joint_backlash_gen::doMotion(
    kte_pass_flag aFlag, const std::shared_ptr<frame_storage>& aStorage) {
  if ((!mBase) || (!mEnd)) {
    return;
  }

  if (mBase->q_dot >= mEnd->q_dot) {
    if (mEnd->q <= mBase->q - 0.5 * mGapSize) {
      // Connection is rigid from input to output
      mEnd->q = mBase->q - 0.5 * mGapSize;
      mEnd->q_dot = mBase->q_dot;
      if (mBase->q_ddot > mEnd->q_ddot) {
        mEnd->q_ddot = mBase->q_ddot;
      }
    }  // Otherwise, connection is loose and mEnd is left unchanged (free coordinate, can be integrated or not).
  } else {
    if (mEnd->q >= mBase->q + 0.5 * mGapSize) {
      // Connection is rigid from input to output
      mEnd->q = mBase->q + 0.5 * mGapSize;
      mEnd->q_dot = mBase->q_dot;
      if (mBase->q_ddot < mEnd->q_ddot) {
        mEnd->q_ddot = mBase->q_ddot;
      }
    }  // Otherwise, connection is loose and mEnd is left unchanged (free coordinate, can be integrated or not).
  }

  if ((aFlag == store_kinematics) && (aStorage)) {
    if (!(aStorage->gen_coord_mapping[mBase])) {
      aStorage->gen_coord_mapping[mBase] =
          std::make_shared<gen_coord<double>>(*mBase);
    } else {
      (*(aStorage->gen_coord_mapping[mBase])) = (*mBase);
    }
    if (!(aStorage->gen_coord_mapping[mEnd])) {
      aStorage->gen_coord_mapping[mEnd] =
          std::make_shared<gen_coord<double>>(*mEnd);
    } else {
      (*(aStorage->gen_coord_mapping[mEnd])) = (*mEnd);
    }
  }
}

void joint_backlash_gen::doForce(
    kte_pass_flag aFlag, const std::shared_ptr<frame_storage>& aStorage) {
  if ((!mBase) || (!mEnd)) {
    return;
  }

  if (mBase->q_dot >= mEnd->q_dot) {
    if ((mEnd->q <= mBase->q - 0.5 * mGapSize) && (mEnd->f < 0.0)) {
      // Connection is rigid from input to output
      mBase->f += mEnd->f;
    }  // Otherwise, connection is loose and mBase is left unchanged (free coordinate, can be integrated or not).
  } else {
    if ((mEnd->q >= mBase->q - 0.5 * mGapSize) && (mEnd->f > 0.0)) {
      // Connection is rigid from input to output
      mBase->f += mEnd->f;
    }  // Otherwise, connection is loose and mEnd is left unchanged (free coordinate, can be integrated or not).
  }

  if ((aFlag == store_dynamics) && (aStorage)) {
    if (aStorage->gen_coord_mapping[mEnd]) {
      aStorage->gen_coord_mapping[mEnd]->f = mEnd->f;
    }
  }
}

void joint_backlash_gen::clearForce() {
  if (mBase) {
    mBase->f = 0.0;
  }
  if (mEnd) {
    mEnd->f = 0.0;
  }
}
}  // namespace ReaK::kte

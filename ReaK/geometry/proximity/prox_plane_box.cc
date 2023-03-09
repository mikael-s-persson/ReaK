
/*
 *    Copyright 2012 Sven Mikael Persson
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

#include "ReaK/geometry/proximity/prox_plane_box.h"

#include <cmath>

namespace ReaK::geom {

proximity_record_3D compute_proximity(const plane& aPlane,
                                      const shape_3D_precompute_pack& aPack1,
                                      const box& aBox,
                                      const shape_3D_precompute_pack& aPack2) {
  using ReaK::unit;
  using std::abs;
  using std::sqrt;
  proximity_record_3D result;

  const pose_3D<double>& bx_pose = aPack1.global_pose;
  const pose_3D<double>& pl_pose = aPack2.global_pose;

  vect<double, 3> bx_c = bx_pose.Position;
  vect<double, 3> bx_x = pl_pose.rotateFromGlobal(
      bx_pose.rotateToGlobal(vect<double, 3>(1.0, 0.0, 0.0)));
  vect<double, 3> bx_y = pl_pose.rotateFromGlobal(
      bx_pose.rotateToGlobal(vect<double, 3>(1.0, 0.0, 0.0)));
  vect<double, 3> bx_z = pl_pose.rotateFromGlobal(
      bx_pose.rotateToGlobal(vect<double, 3>(1.0, 0.0, 0.0)));

  const vect<double, 3> bx_dim = aBox.getDimensions();

  if (bx_x[2] > 0.0) {
    bx_x = -bx_x;
  }
  if (bx_y[2] > 0.0) {
    bx_y = -bx_y;
  }
  if (bx_z[2] > 0.0) {
    bx_z = -bx_z;
  }

  vect<double, 3> bx_c_rel = pl_pose.transformFromGlobal(bx_c);
  vect<double, 3> bx_pt_rel =
      bx_c_rel + 0.5 * (bx_dim[0] * bx_x + bx_dim[1] * bx_y + bx_dim[2] * bx_z);

  result.mPoint1 = pl_pose.transformToGlobal(
      vect<double, 3>(bx_pt_rel[0], bx_pt_rel[1], 0.0));
  result.mPoint2 = pl_pose.transformToGlobal(bx_pt_rel);
  result.mDistance = bx_pt_rel[2];
  return result;
}

proximity_record_3D compute_proximity(const box& aBox,
                                      const shape_3D_precompute_pack& aPack1,
                                      const plane& aPlane,
                                      const shape_3D_precompute_pack& aPack2) {
  using std::swap;
  proximity_record_3D result =
      compute_proximity(aPlane, aPack2, aBox, aPack1);  // NOLINT
  swap(result.mPoint1, result.mPoint2);
  return result;
}

proximity_record_3D prox_plane_box::computeProximity(
    const shape_3D_precompute_pack& aPack1,
    const shape_3D_precompute_pack& aPack2) {
  if ((mBox == nullptr) || (mPlane == nullptr)) {
    return {};
  }

  if (aPack1.parent == mPlane) {
    return compute_proximity(*mPlane, aPack1, *mBox, aPack2);
  }
  return compute_proximity(*mBox, aPack1, *mPlane, aPack2);
}

prox_plane_box::prox_plane_box(const plane* aPlane, const box* aBox)
    : mPlane(aPlane), mBox(aBox) {}

}  // namespace ReaK::geom


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

#include "ReaK/geometry/proximity/prox_ccylinder_box.hpp"

#include "ReaK/geometry/proximity/prox_fundamentals_3D.hpp"

namespace ReaK::geom {

proximity_record_3D compute_proximity(const capped_cylinder& aCCylinder,
                                      const shape_3D_precompute_pack& aPack1,
                                      const box& aBox,
                                      const shape_3D_precompute_pack& aPack2) {
  using std::abs;
  using std::sqrt;

  const pose_3D<double>& cy_pose = aPack1.global_pose;
  const pose_3D<double>& bx_pose = aPack2.global_pose;

  const vect<double, 3> cy_c = cy_pose.Position;
  const vect<double, 3> cy_t =
      cy_pose.rotateToGlobal(vect<double, 3>(0.0, 0.0, 1.0));
  const double cy_len = aCCylinder.getLength();
  const double cy_rad = aCCylinder.getRadius();

  proximity_record_3D bxln_result =
      findProximityBoxToLine(aBox, bx_pose, cy_c, cy_t, 0.5 * cy_len);

  // add a sphere-sweep around the point-box solution.
  vect<double, 3> diff_v = bxln_result.mPoint1 - bxln_result.mPoint2;
  double diff_d = norm_2(diff_v);
  proximity_record_3D result;
  if (bxln_result.mDistance < 0.0) {
    result.mPoint1 = bxln_result.mPoint2 - (cy_rad / diff_d) * diff_v;
  } else {
    result.mPoint1 = bxln_result.mPoint2 + (cy_rad / diff_d) * diff_v;
  }
  result.mPoint2 = bxln_result.mPoint1;
  result.mDistance = bxln_result.mDistance - cy_rad;
  return result;
}

proximity_record_3D compute_proximity(const box& aBox,
                                      const shape_3D_precompute_pack& aPack1,
                                      const capped_cylinder& aCCylinder,
                                      const shape_3D_precompute_pack& aPack2) {
  using std::swap;
  proximity_record_3D result =
      compute_proximity(aCCylinder, aPack2, aBox, aPack1);  // NOLINT
  swap(result.mPoint1, result.mPoint2);
  return result;
}

proximity_record_3D prox_ccylinder_box::computeProximity(
    const shape_3D_precompute_pack& aPack1,
    const shape_3D_precompute_pack& aPack2) {
  if ((mCCylinder == nullptr) || (mBox == nullptr)) {
    return {};
  }

  if (aPack1.parent == mCCylinder) {
    return compute_proximity(*mCCylinder, aPack1, *mBox, aPack2);
  }
  return compute_proximity(*mBox, aPack1, *mCCylinder, aPack2);
}

prox_ccylinder_box::prox_ccylinder_box(const capped_cylinder* aCCylinder,
                                       const box* aBox)
    : mCCylinder(aCCylinder), mBox(aBox) {}

}  // namespace ReaK::geom

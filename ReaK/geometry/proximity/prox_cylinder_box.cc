
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

#include "ReaK/geometry/proximity/prox_cylinder_box.h"

#include "ReaK/geometry/proximity/prox_fundamentals_3D.h"

namespace ReaK::geom {

proximity_record_3D compute_proximity(const cylinder& aCylinder,
                                      const shape_3D_precompute_pack& aPack1,
                                      const box& aBox,
                                      const shape_3D_precompute_pack& aPack2) {

  return findProximityByGJKEPA(
      cylinder_support_func(aCylinder, aPack1.global_pose),
      box_support_func(aBox, aPack2.global_pose));
}

proximity_record_3D compute_proximity(const box& aBox,
                                      const shape_3D_precompute_pack& aPack1,
                                      const cylinder& aCylinder,
                                      const shape_3D_precompute_pack& aPack2) {
  using std::swap;
  proximity_record_3D result =
      compute_proximity(aCylinder, aPack2, aBox, aPack1);  // NOLINT
  swap(result.mPoint1, result.mPoint2);
  return result;
}

proximity_record_3D prox_cylinder_box::computeProximity(
    const shape_3D_precompute_pack& aPack1,
    const shape_3D_precompute_pack& aPack2) {
  if ((mCylinder == nullptr) || (mBox == nullptr)) {
    return {};
  }

  if (aPack1.parent == mCylinder) {
    return compute_proximity(*mCylinder, aPack1, *mBox, aPack2);
  }
  return compute_proximity(*mBox, aPack1, *mCylinder, aPack2);
}

prox_cylinder_box::prox_cylinder_box(const cylinder* aCylinder, const box* aBox)
    : mCylinder(aCylinder), mBox(aBox) {}

}  // namespace ReaK::geom

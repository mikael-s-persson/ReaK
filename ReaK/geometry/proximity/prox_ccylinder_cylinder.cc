
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

#include "ReaK/geometry/proximity/prox_ccylinder_cylinder.h"

#include "ReaK/geometry/proximity/prox_fundamentals_3D.h"

#include <cmath>

namespace ReaK::geom {

proximity_record_3D compute_proximity(const capped_cylinder& aCCylinder,
                                      const shape_3D_precompute_pack& aPack1,
                                      const cylinder& aCylinder,
                                      const shape_3D_precompute_pack& aPack2) {

  return findProximityByGJKEPA(
      ccylinder_support_func(aCCylinder, aPack1.global_pose),
      cylinder_support_func(aCylinder, aPack2.global_pose));
}

proximity_record_3D compute_proximity(const cylinder& aCylinder,
                                      const shape_3D_precompute_pack& aPack1,
                                      const capped_cylinder& aCCylinder,
                                      const shape_3D_precompute_pack& aPack2) {
  using std::swap;
  proximity_record_3D result =
      compute_proximity(aCCylinder, aPack2, aCylinder, aPack1);  // NOLINT
  swap(result.mPoint1, result.mPoint2);
  return result;
}

proximity_record_3D prox_ccylinder_cylinder::computeProximity(
    const shape_3D_precompute_pack& aPack1,
    const shape_3D_precompute_pack& aPack2) {
  if ((mCCylinder == nullptr) || (mCylinder == nullptr)) {
    return {};
  }

  if (aPack1.parent == mCCylinder) {
    return compute_proximity(*mCCylinder, aPack1, *mCylinder, aPack2);
  }
  return compute_proximity(*mCylinder, aPack1, *mCCylinder, aPack2);
}

prox_ccylinder_cylinder::prox_ccylinder_cylinder(
    const capped_cylinder* aCCylinder, const cylinder* aCylinder)
    : mCCylinder(aCCylinder), mCylinder(aCylinder) {}

}  // namespace ReaK::geom

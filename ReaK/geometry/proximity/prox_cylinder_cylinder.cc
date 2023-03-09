
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

#include "ReaK/geometry/proximity/prox_cylinder_cylinder.h"

#include "ReaK/geometry/proximity/prox_fundamentals_3D.h"

#include <cmath>

namespace ReaK::geom {

proximity_record_3D compute_proximity(const cylinder& aCylinder1,
                                      const shape_3D_precompute_pack& aPack1,
                                      const cylinder& aCylinder2,
                                      const shape_3D_precompute_pack& aPack2) {

  return findProximityByGJKEPA(
      cylinder_support_func(aCylinder1, aPack1.global_pose),
      cylinder_support_func(aCylinder2, aPack2.global_pose));
}

proximity_record_3D prox_cylinder_cylinder::computeProximity(
    const shape_3D_precompute_pack& aPack1,
    const shape_3D_precompute_pack& aPack2) {
  if ((mCylinder1 == nullptr) || (mCylinder2 == nullptr)) {
    return {};
  }

  if (aPack1.parent == mCylinder1) {
    return compute_proximity(*mCylinder1, aPack1, *mCylinder2, aPack2);
  }
  return compute_proximity(*mCylinder2, aPack1, *mCylinder1, aPack2);
}

prox_cylinder_cylinder::prox_cylinder_cylinder(const cylinder* aCylinder1,
                                               const cylinder* aCylinder2)
    : mCylinder1(aCylinder1), mCylinder2(aCylinder2) {}

}  // namespace ReaK::geom


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

#include "ReaK/geometry/proximity/prox_sphere_sphere.hpp"

namespace ReaK::geom {

proximity_record_3D compute_proximity(const sphere& aSphere1,
                                      const shape_3D_precompute_pack& aPack1,
                                      const sphere& aSphere2,
                                      const shape_3D_precompute_pack& aPack2) {
  proximity_record_3D result;

  vect<double, 3> c1 = aPack1.global_pose.Position;
  vect<double, 3> c2 = aPack2.global_pose.Position;

  const double s1_rad = aSphere1.getRadius();
  const double s2_rad = aSphere2.getRadius();

  vect<double, 3> diff_cc = c2 - c1;
  double dist_cc = norm_2(diff_cc);

  result.mDistance = dist_cc - s1_rad - s2_rad;
  result.mPoint1 = c1 + (s1_rad / dist_cc) * diff_cc;
  result.mPoint2 = c2 - (s2_rad / dist_cc) * diff_cc;
  return result;
}

proximity_record_3D prox_sphere_sphere::computeProximity(
    const shape_3D_precompute_pack& aPack1,
    const shape_3D_precompute_pack& aPack2) {
  if ((mSphere1 == nullptr) || (mSphere2 == nullptr)) {
    return {};
  }

  if (aPack1.parent == mSphere1) {
    return compute_proximity(*mSphere1, aPack1, *mSphere2, aPack2);
  }
  return compute_proximity(*mSphere2, aPack1, *mSphere1, aPack2);
}

prox_sphere_sphere::prox_sphere_sphere(const sphere* aSphere1,
                                       const sphere* aSphere2)
    : mSphere1(aSphere1), mSphere2(aSphere2) {}

}  // namespace ReaK::geom


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

#include <ReaK/geometry/proximity/prox_circle_circle.hpp>

namespace ReaK::geom {

proximity_record_2D compute_proximity(const circle& aCircle1,
                                      const shape_2D_precompute_pack& aPack1,
                                      const circle& aCircle2,
                                      const shape_2D_precompute_pack& aPack2) {
  proximity_record_2D result;

  const vect<double, 2> c1 = aPack1.global_pose.Position;
  const vect<double, 2> c2 = aPack2.global_pose.Position;

  const double c1_rad = aCircle1.getRadius();
  const double c2_rad = aCircle2.getRadius();

  const vect<double, 2> diff_cc = c2 - c1;
  const double dist_cc = norm_2(diff_cc);

  result.mDistance = dist_cc - c1_rad - c2_rad;
  result.mPoint1 = c1 + (c1_rad / dist_cc) * diff_cc;
  result.mPoint2 = c2 - (c2_rad / dist_cc) * diff_cc;

  return result;
}

proximity_record_2D prox_circle_circle::computeProximity(
    const shape_2D_precompute_pack& aPack1,
    const shape_2D_precompute_pack& aPack2) {
  if ((mCircle1 == nullptr) || (mCircle2 == nullptr)) {
    return {};
  }

  if (aPack1.parent == mCircle1) {
    return compute_proximity(*mCircle1, aPack1, *mCircle2, aPack2);
  }
  return compute_proximity(*mCircle2, aPack1, *mCircle1, aPack2);
}

prox_circle_circle::prox_circle_circle(const circle* aCircle1,
                                       const circle* aCircle2)
    : mCircle1(aCircle1), mCircle2(aCircle2) {}

}  // namespace ReaK::geom

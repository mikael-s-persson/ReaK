
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

#include <ReaK/geometry/proximity/prox_circle_rectangle.hpp>

namespace ReaK::geom {

proximity_record_2D compute_proximity(const circle& aCircle,
                                      const shape_2D_precompute_pack& aPack1,
                                      const rectangle& aRectangle,
                                      const shape_2D_precompute_pack& aPack2) {
  using std::abs;
  proximity_record_2D result;

  const pose_2D<double>& ci_pose = aPack1.global_pose;
  const pose_2D<double>& re_pose = aPack2.global_pose;

  const vect<double, 2> ci_c = ci_pose.Position;
  const vect<double, 2> ci_c_rel = re_pose.transformFromGlobal(ci_c);

  const double ci_rad = aCircle.getRadius();
  const vect<double, 2> re_dim = aRectangle.getDimensions();

  bool in_x_range =
      ((ci_c_rel[0] > -0.5 * re_dim[0]) && (ci_c_rel[0] < 0.5 * re_dim[0]));
  bool in_y_range =
      ((ci_c_rel[1] > -0.5 * re_dim[1]) && (ci_c_rel[1] < 0.5 * re_dim[1]));

  if (in_x_range && in_y_range) {
    // The circle is inside the rectangle.
    const vect<double, 2> bound_dists(0.5 * re_dim[0] - abs(ci_c_rel[0]),
                                      0.5 * re_dim[1] - abs(ci_c_rel[1]));
    if (bound_dists[0] <= bound_dists[1]) {
      in_x_range = false;
    } else {
      in_y_range = false;
    }
  }

  vect<double, 2> corner_pt = 0.5 * re_dim;
  if (in_x_range) {
    corner_pt[0] = ci_c_rel[0];
  } else if (ci_c_rel[0] < 0.0) {
    corner_pt[0] = -corner_pt[0];
  }
  if (in_y_range) {
    corner_pt[1] = ci_c_rel[1];
  } else if (ci_c_rel[1] < 0.0) {
    corner_pt[1] = -corner_pt[1];
  }
  result.mPoint2 = re_pose.transformToGlobal(corner_pt);
  const vect<double, 2> diff_v = result.mPoint2 - ci_c;
  const double diff_d = norm_2(diff_v);
  result.mPoint1 = ci_c + (ci_rad / diff_d) * diff_v;
  result.mDistance = diff_d - ci_rad;
  return result;
}

proximity_record_2D compute_proximity(const rectangle& aRectangle,
                                      const shape_2D_precompute_pack& aPack1,
                                      const circle& aCircle,
                                      const shape_2D_precompute_pack& aPack2) {
  using std::swap;
  proximity_record_2D result =
      compute_proximity(aCircle, aPack2, aRectangle, aPack1);  // NOLINT
  swap(result.mPoint1, result.mPoint2);
  return result;
}

proximity_record_2D prox_circle_rectangle::computeProximity(
    const shape_2D_precompute_pack& aPack1,
    const shape_2D_precompute_pack& aPack2) {
  if ((mCircle == nullptr) || (mRectangle == nullptr)) {
    return {};
  }

  if (aPack1.parent == mCircle) {
    return compute_proximity(*mCircle, aPack1, *mRectangle, aPack2);
  }
  return compute_proximity(*mRectangle, aPack1, *mCircle, aPack2);
}

prox_circle_rectangle::prox_circle_rectangle(const circle* aCircle,
                                             const rectangle* aRectangle)
    : mCircle(aCircle), mRectangle(aRectangle) {}

}  // namespace ReaK::geom

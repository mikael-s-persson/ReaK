
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

#include <cmath>

#include <ReaK/geometry/proximity/prox_rectangle_rectangle.hpp>

namespace ReaK::geom {

void compute_proximity_of_point(const rectangle& aRectangle,
                                const pose_2D<double>& aRecGblPose,
                                const vect<double, 2>& aPoint,
                                vect<double, 2>& aPointRec, double& aDistance) {
  using std::abs;

  vect<double, 2> pt_rel = aRecGblPose.transformFromGlobal(aPoint);

  const vect<double, 2> re_dim = aRectangle.getDimensions();

  bool in_x_range =
      ((pt_rel[0] > -0.5 * re_dim[0]) && (pt_rel[0] < 0.5 * re_dim[0]));
  bool in_y_range =
      ((pt_rel[1] > -0.5 * re_dim[1]) && (pt_rel[1] < 0.5 * re_dim[1]));

  if (in_x_range && in_y_range) {
    // The circle is inside the rectangle.
    vect<double, 2> bound_dists = vect<double, 2>(
        0.5 * re_dim[0] - abs(pt_rel[0]), 0.5 * re_dim[1] - abs(pt_rel[1]));
    if (bound_dists[0] <= bound_dists[1]) {
      in_x_range = false;
    } else {
      in_y_range = false;
    }
  }

  vect<double, 2> corner_pt = 0.5 * re_dim;
  if (in_x_range) {
    corner_pt[0] = pt_rel[0];
  } else if (pt_rel[0] < 0.0) {
    corner_pt[0] = -corner_pt[0];
  }
  if (in_y_range) {
    corner_pt[1] = pt_rel[1];
  } else if (pt_rel[1] < 0.0) {
    corner_pt[1] = -corner_pt[1];
  }
  aPointRec = aRecGblPose.transformToGlobal(corner_pt);
  aDistance = norm_2(aPointRec - aPoint);
}

proximity_record_2D compute_proximity(const rectangle& aRectangle1,
                                      const shape_2D_precompute_pack& aPack1,
                                      const rectangle& aRectangle2,
                                      const shape_2D_precompute_pack& aPack2) {
  proximity_record_2D result;
  result.mDistance = std::numeric_limits<double>::infinity();

  const pose_2D<double>& r1_pose = aPack1.global_pose;
  const pose_2D<double>& r2_pose = aPack2.global_pose;

  vect<double, 2> temp_pt;
  double temp_dist = NAN;

  vect<double, 2> corner = 0.5 * aRectangle2.getDimensions();
  vect<double, 2> corner_gbl = r2_pose.transformToGlobal(corner);
  compute_proximity_of_point(aRectangle1, r1_pose, corner_gbl, temp_pt,
                             temp_dist);
  if (temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint1 = temp_pt;
    result.mPoint2 = corner_gbl;
  }

  corner[1] = -corner[1];
  corner_gbl = r2_pose.transformToGlobal(corner);
  compute_proximity_of_point(aRectangle1, r1_pose, corner_gbl, temp_pt,
                             temp_dist);
  if (temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint1 = temp_pt;
    result.mPoint2 = corner_gbl;
  }

  corner[0] = -corner[0];
  corner_gbl = r2_pose.transformToGlobal(corner);
  compute_proximity_of_point(aRectangle1, r1_pose, corner_gbl, temp_pt,
                             temp_dist);
  if (temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint1 = temp_pt;
    result.mPoint2 = corner_gbl;
  }

  corner[1] = -corner[1];
  corner_gbl = r2_pose.transformToGlobal(corner);
  compute_proximity_of_point(aRectangle1, r1_pose, corner_gbl, temp_pt,
                             temp_dist);
  if (temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint1 = temp_pt;
    result.mPoint2 = corner_gbl;
  }

  corner = 0.5 * aRectangle1.getDimensions();
  corner_gbl = r1_pose.transformToGlobal(corner);
  compute_proximity_of_point(aRectangle2, r2_pose, corner_gbl, temp_pt,
                             temp_dist);
  if (temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint2 = temp_pt;
    result.mPoint1 = corner_gbl;
  }

  corner[1] = -corner[1];
  corner_gbl = r1_pose.transformToGlobal(corner);
  compute_proximity_of_point(aRectangle2, r2_pose, corner_gbl, temp_pt,
                             temp_dist);
  if (temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint2 = temp_pt;
    result.mPoint1 = corner_gbl;
  }

  corner[0] = -corner[0];
  corner_gbl = r1_pose.transformToGlobal(corner);
  compute_proximity_of_point(aRectangle2, r2_pose, corner_gbl, temp_pt,
                             temp_dist);
  if (temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint2 = temp_pt;
    result.mPoint1 = corner_gbl;
  }

  corner[1] = -corner[1];
  corner_gbl = r1_pose.transformToGlobal(corner);
  compute_proximity_of_point(aRectangle2, r2_pose, corner_gbl, temp_pt,
                             temp_dist);
  if (temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint2 = temp_pt;
    result.mPoint1 = corner_gbl;
  }

  return result;
}

proximity_record_2D prox_rectangle_rectangle::computeProximity(
    const shape_2D_precompute_pack& aPack1,
    const shape_2D_precompute_pack& aPack2) {
  if ((mRectangle1 == nullptr) || (mRectangle2 == nullptr)) {
    return {};
  }

  if (aPack1.parent == mRectangle1) {
    return compute_proximity(*mRectangle1, aPack1, *mRectangle2, aPack2);
  }
  return compute_proximity(*mRectangle2, aPack1, *mRectangle1, aPack2);
}

prox_rectangle_rectangle::prox_rectangle_rectangle(const rectangle* aRectangle1,
                                                   const rectangle* aRectangle2)
    : mRectangle1(aRectangle1), mRectangle2(aRectangle2) {}

}  // namespace ReaK::geom

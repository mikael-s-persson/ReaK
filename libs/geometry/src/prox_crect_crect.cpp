
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

#include "ReaK/geometry/proximity/prox_crect_crect.hpp"

namespace ReaK::geom {

proximity_record_2D compute_proximity(const capped_rectangle& aCRect1,
                                      const shape_2D_precompute_pack& aPack1,
                                      const capped_rectangle& aCRect2,
                                      const shape_2D_precompute_pack& aPack2) {
  using std::abs;
  proximity_record_2D result;

  const pose_2D<double>& cr1_pose = aPack1.global_pose;
  const pose_2D<double>& cr2_pose = aPack2.global_pose;

  const vect<double, 2> cr2_c = cr2_pose.Position;
  const vect<double, 2> cr2_t =
      cr2_pose.rotateToGlobal(vect<double, 2>(1.0, 0.0));

  const vect<double, 2> cr2_c_rel = cr1_pose.transformFromGlobal(cr2_c);
  const vect<double, 2> cr2_t_rel = cr1_pose.rotateFromGlobal(cr2_t);

  const vect<double, 2> cr1_dim = aCRect1.getDimensions();
  const vect<double, 2> cr2_dim = aCRect2.getDimensions();

  if (abs(cr2_t_rel[1]) < 1e-5) {
    // The capped-rectangles are parallel.
    if ((cr2_c_rel[0] + 0.5 * cr2_dim[0] > -0.5 * cr1_dim[0]) ||
        (cr2_c_rel[0] - 0.5 * cr2_dim[0] < 0.5 * cr1_dim[0])) {
      // there is an overlap between the capped-rectangle sides.
      const double max_x_rel =
          ((cr2_c_rel[0] + 0.5 * cr2_dim[0] < 0.5 * cr1_dim[0])
               ? (cr2_c_rel[0] + 0.5 * cr2_dim[0])
               : (0.5 * cr1_dim[0]));
      const double min_x_rel =
          ((cr2_c_rel[0] - 0.5 * cr2_dim[0] > -0.5 * cr1_dim[0])
               ? (cr2_c_rel[0] - 0.5 * cr2_dim[0])
               : (-0.5 * cr1_dim[0]));
      const double avg_x_rel = (max_x_rel + min_x_rel) * 0.5;
      vect<double, 2> cr2_r_rel(0.0, 1.0);
      if (cr2_c_rel[1] < 0.0) {
        cr2_r_rel[1] = -1.0;
      }
      result.mPoint1 = cr1_pose.transformToGlobal(
          vect<double, 2>(avg_x_rel, 0.5 * cr1_dim[1] * cr2_r_rel[1]));
      result.mPoint2 = cr1_pose.transformToGlobal(vect<double, 2>(
          avg_x_rel, cr2_c_rel[1] - 0.5 * cr2_dim[1] * cr2_r_rel[1]));
      result.mDistance =
          abs(cr2_c_rel[1]) - 0.5 * cr1_dim[1] - 0.5 * cr2_dim[1];
      return result;
    }
    // there is no overlap, and thus, this boils down to a circle-circle problem.
    vect<double, 2> cr1_cic_rel(0.0, 0.0);
    vect<double, 2> cr2_cic_rel = cr2_c_rel;
    if (cr2_c_rel[0] < 0.0) {
      cr1_cic_rel[0] -= 0.5 * cr1_dim[0];
      cr2_cic_rel[0] += 0.5 * cr2_dim[0];
    } else {
      cr1_cic_rel[0] += 0.5 * cr1_dim[0];
      cr2_cic_rel[0] -= 0.5 * cr2_dim[0];
    }
    const vect<double, 2> diff_v_rel = cr2_cic_rel - cr1_cic_rel;
    const double dist_v_rel = norm_2(diff_v_rel);
    result.mPoint1 = cr1_pose.transformToGlobal(
        cr1_cic_rel + (0.5 * cr1_dim[1] / dist_v_rel) * diff_v_rel);
    result.mPoint2 = cr1_pose.transformToGlobal(
        cr2_cic_rel - (0.5 * cr2_dim[1] / dist_v_rel) * diff_v_rel);
    result.mDistance = dist_v_rel - 0.5 * cr1_dim[1] - 0.5 * cr2_dim[1];
    return result;
  }

  // Line-Line solution:
  const double d = cr2_t_rel * cr2_c_rel;
  const double denom = 1.0 - cr2_t_rel[0] * cr2_t_rel[0];
  double s_c = (cr2_t_rel[0] * cr2_c_rel[0] - d) / denom;
  double t_c = (cr2_c_rel[0] - cr2_t_rel[0] * d) / denom;

  // Segment-Segment solution:
  if (s_c < -0.5 * cr2_dim[0]) {
    s_c = -0.5 * cr2_dim[0];
    t_c = cr2_c_rel[0] - 0.5 * cr2_dim[0] * cr2_t_rel[0];
  } else if (s_c > 0.5 * cr2_dim[0]) {
    s_c = 0.5 * cr2_dim[0];
    t_c = cr2_c_rel[0] + 0.5 * cr2_dim[0] * cr2_t_rel[0];
  }

  if (t_c < -0.5 * cr1_dim[0]) {
    t_c = -0.5 * cr1_dim[0];
    s_c = -0.5 * cr1_dim[0] * cr2_t_rel[0] - d;
  } else if (t_c > 0.5 * cr1_dim[0]) {
    t_c = 0.5 * cr1_dim[0];
    s_c = 0.5 * cr1_dim[0] * cr2_t_rel[0] - d;
  }

  if (s_c < -0.5 * cr2_dim[0]) {
    s_c = -0.5 * cr2_dim[0];
  } else if (s_c > 0.5 * cr2_dim[0]) {
    s_c = 0.5 * cr2_dim[0];
  }

  // we have parameters s and t for the min-dist points on the center segments.

  // just apply a circle-sweep on the line-segments.
  const vect<double, 2> cr1_ptc(t_c, 0.0);
  const vect<double, 2> cr2_ptc = cr2_c_rel + s_c * cr2_t_rel;

  const vect<double, 2> diff_v_rel = cr2_ptc - cr1_ptc;
  const double dist_v_rel = norm_2(diff_v_rel);
  result.mPoint1 = cr1_pose.transformToGlobal(
      cr1_ptc + (0.5 * cr1_dim[1] / dist_v_rel) * diff_v_rel);
  result.mPoint2 = cr1_pose.transformToGlobal(
      cr2_ptc - (0.5 * cr2_dim[1] / dist_v_rel) * diff_v_rel);
  result.mDistance = dist_v_rel - 0.5 * cr1_dim[1] - 0.5 * cr2_dim[1];
  return result;
}

proximity_record_2D prox_crect_crect::computeProximity(
    const shape_2D_precompute_pack& aPack1,
    const shape_2D_precompute_pack& aPack2) {
  if ((mCRect1 == nullptr) || (mCRect2 == nullptr)) {
    return {};
  }

  if (aPack1.parent == mCRect1) {
    return compute_proximity(*mCRect1, aPack1, *mCRect2, aPack2);
  }
  return compute_proximity(*mCRect2, aPack1, *mCRect1, aPack2);
}

prox_crect_crect::prox_crect_crect(const capped_rectangle* aCRect1,
                                   const capped_rectangle* aCRect2)
    : mCRect1(aCRect1), mCRect2(aCRect2) {}

}  // namespace ReaK::geom


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

#include <ReaK/geometry/proximity/prox_crect_rectangle.hpp>

namespace ReaK::geom {

proximity_record_2D compute_proximity_of_line(const rectangle& aRectangle,
                                              const pose_2D<double>& aGblPose,
                                              const vect<double, 2>& ln_c,
                                              const vect<double, 2>& ln_t,
                                              double half_length) {
  using std::abs;
  proximity_record_2D result;

  const vect<double, 2> re_dim = aRectangle.getDimensions();

  vect<double, 2> ln_c_rel = aGblPose.transformFromGlobal(ln_c);
  vect<double, 2> ln_t_rel = aGblPose.rotateFromGlobal(ln_t);

  if (abs(ln_t_rel[0]) < 1e-5) {
    // this means the line is vertical.
    if ((ln_c_rel[1] + half_length > -0.5 * re_dim[1]) ||
        (ln_c_rel[1] - half_length < 0.5 * re_dim[1])) {
      // there is an overlap between the rectangle side and the line.
      double max_y_rel = ((ln_c_rel[1] + half_length < 0.5 * re_dim[1])
                              ? (ln_c_rel[1] + half_length)
                              : (0.5 * re_dim[1]));
      double min_y_rel = ((ln_c_rel[1] - half_length > -0.5 * re_dim[1])
                              ? (ln_c_rel[1] - half_length)
                              : (-0.5 * re_dim[1]));
      double avg_y_rel = (max_y_rel + min_y_rel) * 0.5;
      vect<double, 2> ln_r_rel(1.0, 0.0);
      if (ln_c_rel[0] < 0.0) {
        ln_r_rel[0] = -1.0;
      }
      result.mPoint1 =
          aGblPose.transformToGlobal(vect<double, 2>(ln_c_rel[0], avg_y_rel));
      result.mPoint2 = aGblPose.transformToGlobal(
          vect<double, 2>(0.5 * re_dim[0] * ln_r_rel[0], avg_y_rel));
      result.mDistance = abs(ln_c_rel[0]) - 0.5 * re_dim[0];
      return result;
    }

    // there is no overlap, and thus, this boils down to a point-point problem.
    vect<double, 2> re_pt_rel(0.0, 0.0);
    vect<double, 2> ln_pt_rel = ln_c_rel;
    if (ln_c_rel[0] < 0.0) {
      re_pt_rel[0] -= 0.5 * re_dim[0];
    } else {
      re_pt_rel[0] += 0.5 * re_dim[0];
    }
    if (ln_c_rel[1] < 0.0) {
      re_pt_rel[1] -= 0.5 * re_dim[1];
      ln_pt_rel[1] += half_length;
    } else {
      re_pt_rel[1] += 0.5 * re_dim[1];
      ln_pt_rel[1] -= half_length;
    }

    vect<double, 2> diff_v_rel = ln_pt_rel - re_pt_rel;
    double dist_v_rel = norm_2(diff_v_rel);
    result.mPoint1 = aGblPose.transformToGlobal(ln_pt_rel);
    result.mPoint2 = aGblPose.transformToGlobal(re_pt_rel);
    result.mDistance = dist_v_rel;
    return result;
  }

  if (abs(ln_t_rel[1]) < 1e-5) {
    // this means the line is horizontal.
    if ((ln_c_rel[0] + half_length > -0.5 * re_dim[0]) ||
        (ln_c_rel[0] - half_length < 0.5 * re_dim[0])) {
      // there is an overlap between the rectangle side and the line.
      double max_x_rel = ((ln_c_rel[0] + half_length < 0.5 * re_dim[0])
                              ? (ln_c_rel[0] + half_length)
                              : (0.5 * re_dim[0]));
      double min_x_rel = ((ln_c_rel[0] - half_length > -0.5 * re_dim[0])
                              ? (ln_c_rel[0] - half_length)
                              : (-0.5 * re_dim[0]));
      double avg_x_rel = (max_x_rel + min_x_rel) * 0.5;
      vect<double, 2> ln_r_rel(0.0, 1.0);
      if (ln_c_rel[1] < 0.0) {
        ln_r_rel[1] = -1.0;
      }
      result.mPoint1 =
          aGblPose.transformToGlobal(vect<double, 2>(avg_x_rel, ln_c_rel[1]));
      result.mPoint2 = aGblPose.transformToGlobal(
          vect<double, 2>(avg_x_rel, 0.5 * re_dim[1] * ln_r_rel[1]));
      result.mDistance = abs(ln_c_rel[1]) - 0.5 * re_dim[1];
      return result;
    }

    // there is no overlap, and thus, this boils down to a point-point problem.
    vect<double, 2> re_pt_rel(0.0, 0.0);
    vect<double, 2> ln_pt_rel = ln_c_rel;
    if (ln_c_rel[1] < 0.0) {
      re_pt_rel[1] -= 0.5 * re_dim[1];
    } else {
      re_pt_rel[1] += 0.5 * re_dim[1];
    }
    if (ln_c_rel[0] < 0.0) {
      re_pt_rel[0] -= 0.5 * re_dim[0];
      ln_pt_rel[0] += half_length;
    } else {
      re_pt_rel[0] += 0.5 * re_dim[0];
      ln_pt_rel[0] -= half_length;
    }

    vect<double, 2> diff_v_rel = ln_pt_rel - re_pt_rel;
    double dist_v_rel = norm_2(diff_v_rel);
    result.mPoint1 = aGblPose.transformToGlobal(ln_pt_rel);
    result.mPoint2 = aGblPose.transformToGlobal(re_pt_rel);
    result.mDistance = dist_v_rel;
    return result;
  }

  // in any other case, we have to resort to the segment-point test.

  vect<double, 2> ln_n_rel = 1.0 % ln_t_rel;
  if (ln_n_rel * ln_c_rel < 0.0) {
    ln_n_rel = -ln_n_rel;
  }
  vect<double, 2> corner_pt(-0.5 * re_dim[0], -0.5 * re_dim[1]);
  if (ln_n_rel[0] > 0.0) {
    corner_pt[0] = 0.5 * re_dim[0];
  }
  if (ln_n_rel[1] > 0.0) {
    corner_pt[1] = 0.5 * re_dim[1];
  }

  vect<double, 2> corner_pt_diff = (ln_c_rel - corner_pt);
  double dist_tmp = corner_pt_diff * ln_n_rel;
  double t_tmp = -(corner_pt_diff * ln_t_rel);
  if (abs(t_tmp) > half_length) {
    if (t_tmp < 0.0) {
      t_tmp = -half_length;
    } else {
      t_tmp = half_length;
    }
    vect<double, 2> ln_pt_rel = ln_c_rel + t_tmp * ln_t_rel;

    double in_x_range = abs(ln_pt_rel[0]) - 0.5 * re_dim[0];
    double in_y_range = abs(ln_pt_rel[1]) - 0.5 * re_dim[1];

    if ((in_x_range < 0.0) && (in_y_range > in_x_range)) {
      corner_pt[0] = ln_pt_rel[0];
      dist_tmp = abs(ln_pt_rel[1]) - 0.5 * re_dim[1];
    } else if ((in_y_range < 0.0) && (in_x_range > in_y_range)) {
      corner_pt[1] = ln_pt_rel[1];
      dist_tmp = abs(ln_pt_rel[0]) - 0.5 * re_dim[0];
    } else {
      if (ln_pt_rel[0] < 0.0) {
        corner_pt[0] = -0.5 * re_dim[0];
      } else {
        corner_pt[0] = 0.5 * re_dim[0];
      }
      if (ln_pt_rel[1] < 0.0) {
        corner_pt[1] = -0.5 * re_dim[1];
      } else {
        corner_pt[1] = 0.5 * re_dim[1];
      }
      dist_tmp = norm_2(ln_pt_rel - corner_pt);
    }

    result.mPoint1 = aGblPose.transformToGlobal(ln_pt_rel);
    result.mPoint2 = aGblPose.transformToGlobal(corner_pt);
    result.mDistance = dist_tmp;
  } else {
    result.mPoint1 =
        aGblPose.transformToGlobal(corner_pt + dist_tmp * ln_n_rel);
    result.mPoint2 = aGblPose.transformToGlobal(corner_pt);
    result.mDistance = dist_tmp;
  }

  return result;
}

proximity_record_2D compute_proximity(const capped_rectangle& aCRect,
                                      const shape_2D_precompute_pack& aPack1,
                                      const rectangle& aRectangle,
                                      const shape_2D_precompute_pack& aPack2) {
  proximity_record_2D result;

  const pose_2D<double>& re_pose = aPack1.global_pose;
  const pose_2D<double>& cr_pose = aPack2.global_pose;

  vect<double, 2> cr_c = cr_pose.Position;
  vect<double, 2> cr_t = cr_pose.rotateToGlobal(vect<double, 2>(1.0, 0.0));

  const vect<double, 2> cr_dim = aCRect.getDimensions();

  result = compute_proximity_of_line(aRectangle, re_pose, cr_c, cr_t,
                                     0.5 * cr_dim[0]);

  // add a circle-sweep around the line-rectangle solution.
  vect<double, 2> diff_v = result.mPoint2 - result.mPoint1;
  double diff_d = norm_2(diff_v);
  if (result.mDistance < 0.0) {
    result.mPoint1 -= (0.5 * cr_dim[1] / diff_d) * diff_v;
  } else {
    result.mPoint1 += (0.5 * cr_dim[1] / diff_d) * diff_v;
  }
  result.mDistance -= 0.5 * cr_dim[1];
  return result;
}

proximity_record_2D compute_proximity(const rectangle& aRectangle,
                                      const shape_2D_precompute_pack& aPack1,
                                      const capped_rectangle& aCRect,
                                      const shape_2D_precompute_pack& aPack2) {
  using std::swap;
  proximity_record_2D result =
      compute_proximity(aCRect, aPack2, aRectangle, aPack1);  // NOLINT
  swap(result.mPoint1, result.mPoint2);
  return result;
}

proximity_record_2D prox_crect_rectangle::computeProximity(
    const shape_2D_precompute_pack& aPack1,
    const shape_2D_precompute_pack& aPack2) {
  if ((mCRect == nullptr) || (mRectangle == nullptr)) {
    return {};
  }

  if (aPack1.parent == mCRect) {
    return compute_proximity(*mCRect, aPack1, *mRectangle, aPack2);
  }
  return compute_proximity(*mRectangle, aPack1, *mCRect, aPack2);
}

prox_crect_rectangle::prox_crect_rectangle(const capped_rectangle* aCRect,
                                           const rectangle* aRectangle)
    : mCRect(aCRect), mRectangle(aRectangle) {}

}  // namespace ReaK::geom

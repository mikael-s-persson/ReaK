
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

#include <ReaK/geometry/proximity/prox_circle_crect.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


proximity_record_2D compute_proximity(const circle& aCircle, 
                                      const shape_2D_precompute_pack& aPack1,
                                      const capped_rectangle& aCRect, 
                                      const shape_2D_precompute_pack& aPack2) {
  using std::fabs;
  proximity_record_2D result;
  
  const pose_2D<double>& ci_pose = aPack1.global_pose;
  const pose_2D<double>& re_pose = aPack2.global_pose;
  
  const vect<double,2> ci_c = ci_pose.Position;
  const vect<double,2> ci_c_rel = re_pose.transformFromGlobal(ci_c);
  
  const double ci_rad = aCircle.getRadius();
  const vect<double,2> cr_dim = aCRect.getDimensions();
  
  const bool in_x_range = ((ci_c_rel[0] > -0.5 * cr_dim[0]) &&
                           (ci_c_rel[0] <  0.5 * cr_dim[0]));
  
  if(in_x_range) {
    if(ci_c_rel[1] > 0.0) {
      result.mPoint1 = re_pose.transformToGlobal(vect<double,2>(ci_c_rel[0], ci_c_rel[1] - ci_rad));
      result.mPoint2 = re_pose.transformToGlobal(vect<double,2>(ci_c_rel[0], 0.5 * cr_dim[1]));
      result.mDistance = ci_c_rel[1] - ci_rad - 0.5 * cr_dim[1];
    } else {
      result.mPoint1 = re_pose.transformToGlobal(vect<double,2>(ci_c_rel[0], ci_c_rel[1] + ci_rad));
      result.mPoint2 = re_pose.transformToGlobal(vect<double,2>(ci_c_rel[0], -0.5 * cr_dim[1]));
      result.mDistance = -0.5 * cr_dim[1] - ci_c_rel[1] - ci_rad;
    };
    return result;
  };
  
  // this boils down to a circle-circle test.
  vect<double,2> re_endc(0.0,0.0);
  if(ci_c_rel[0] > 0.0)
    re_endc[0] += 0.5 * cr_dim[0];
  else
    re_endc[0] -= 0.5 * cr_dim[0];
  const vect<double,2> diff_v_rel = ci_c_rel - re_endc;
  const double diff_d_rel = norm_2(diff_v_rel);
  result.mPoint1 = re_pose.transformToGlobal(ci_c_rel - (ci_rad / diff_d_rel) * diff_v_rel);
  result.mPoint2 = re_pose.transformToGlobal(re_endc + (0.5 * cr_dim[1] / diff_d_rel) * diff_v_rel);
  result.mDistance = diff_d_rel - 0.5 * cr_dim[1] - ci_rad;
  return result;
};

proximity_record_2D compute_proximity(const capped_rectangle& aCRect, 
                                      const shape_2D_precompute_pack& aPack1,
                                      const circle& aCircle, 
                                      const shape_2D_precompute_pack& aPack2) {
  using std::swap;
  proximity_record_2D result = compute_proximity(aCircle, aPack2, aCRect, aPack1);
  swap(result.mPoint1,result.mPoint2);
  return result;
};

void prox_circle_crect::computeProximity(const shape_2D_precompute_pack& aPack1, 
                                         const shape_2D_precompute_pack& aPack2) {
  if((!mCircle) || (!mCRect)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,2>(0.0,0.0);
    mLastResult.mPoint2 = vect<double,2>(0.0,0.0);
    return;
  };
  
  if( aPack1.parent == mCircle )
    mLastResult = compute_proximity(*mCircle,aPack1,*mCRect,aPack2);
  else
    mLastResult = compute_proximity(*mCRect,aPack1,*mCircle,aPack2);
  
};


prox_circle_crect::prox_circle_crect(const circle* aCircle,
                                     const capped_rectangle* aCRect) :
                                     proximity_finder_2D(),
                                     mCircle(aCircle),
                                     mCRect(aCRect) { };


};

};


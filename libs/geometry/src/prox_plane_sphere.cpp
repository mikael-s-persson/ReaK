
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

#include <ReaK/geometry/proximity/prox_plane_sphere.hpp>

#include <cmath>

namespace ReaK::geom {

proximity_record_3D compute_proximity(const plane& aPlane,
                                      const shape_3D_precompute_pack& aPack1,
                                      const sphere& aSphere,
                                      const shape_3D_precompute_pack& aPack2) {
  proximity_record_3D result;

  const pose_3D<double>& sp_pose = aPack1.global_pose;
  const pose_3D<double>& pl_pose = aPack2.global_pose;

  vect<double, 3> sp_c =
      sp_pose.transformToGlobal(vect<double, 3>(0.0, 0.0, 0.0));
  vect<double, 3> sp_c_rel = pl_pose.transformFromGlobal(sp_c);

  const double sp_rad = aSphere.getRadius();

  result.mPoint1 =
      pl_pose.transformToGlobal(vect<double, 3>(sp_c_rel[0], sp_c_rel[1], 0.0));
  result.mPoint2 = pl_pose.transformToGlobal(
      vect<double, 3>(sp_c_rel[0], sp_c_rel[1], sp_c_rel[2] - sp_rad));
  result.mDistance = sp_c_rel[2] - sp_rad;
  return result;
}

proximity_record_3D compute_proximity(const sphere& aSphere,
                                      const shape_3D_precompute_pack& aPack1,
                                      const plane& aPlane,
                                      const shape_3D_precompute_pack& aPack2) {
  using std::swap;
  proximity_record_3D result =
      compute_proximity(aPlane, aPack2, aSphere, aPack1);  // NOLINT
  swap(result.mPoint1, result.mPoint2);
  return result;
}

/*
// this version assumes a finite plane for proximity purposes.
void prox_plane_sphere::computeProximity(const shape_3D_precompute_pack& aPack1,
                                         const shape_3D_precompute_pack& aPack2) {
  proximity_record_3D result;
  if((!mSphere) || (!mPlane))
    return result;

  using std::abs; using std::sqrt;

  const pose_3D<double>& sp_pose = (aPack1.parent == mSphere ?
                                    aPack1.global_pose : aPack2.global_pose);
  const pose_3D<double>& pl_pose = (aPack1.parent == mSphere ?
                                    aPack2.global_pose : aPack1.global_pose);

  vect<double,3> sp_c = sp_pose.transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  vect<double,3> pl_c = pl_pose.transformToGlobal(vect<double,3>(0.0,0.0,0.0));

  vect<double,3> sp_c_rel = pl_pose.transformFromGlobal(sp_c);

  if((sp_c_rel[0] > -0.5 * mPlane->getDimensions()[0]) &&
     (sp_c_rel[0] <  0.5 * mPlane->getDimensions()[0]) &&
     (sp_c_rel[1] > -0.5 * mPlane->getDimensions()[1]) &&
     (sp_c_rel[1] <  0.5 * mPlane->getDimensions()[1])) {
    // The sphere is within the dimensions of the plane.

    double fact = 1.0;
    if(sp_c_rel[2] < 0.0)
      fact = -1.0;

    result.mPoint1 = pl_pose.transformToGlobal(vect<double,3>(sp_c_rel[0],sp_c_rel[1],0.0));
    result.mPoint2 = pl_pose.transformToGlobal(vect<double,3>(sp_c_rel[0],sp_c_rel[1],sp_c_rel[2] - fact *
mSphere->getRadius()));
    result.mDistance = fact * sp_c_rel[2] - mSphere->getRadius();
  } else {
    vect<double,3> rim_pt;
    if((sp_c_rel[0] > -0.5 * mPlane->getDimensions()[0]) &&
       (sp_c_rel[0] <  0.5 * mPlane->getDimensions()[0])) {
      // The sphere is on the above or below the plane (y-axis).
      double fact = 1.0;
      if(sp_c_rel[1] < 0.0)
        fact = -1.0;
      vect<double,3> rim_pt = vect<double,3>(sp_c_rel[0],fact * 0.5 * mPlane->getDimensions()[1],0.0);
      result.mPoint1 = pl_pose.transformToGlobal(rim_pt);
    } else if((sp_c_rel[1] > -0.5 * mPlane->getDimensions()[1]) &&
              (sp_c_rel[1] <  0.5 * mPlane->getDimensions()[1])) {
      // The sphere is on the right or left of the plane (x-axis).
      double fact = 1.0;
      if(sp_c_rel[0] < 0.0)
        fact = -1.0;
      vect<double,3> rim_pt = vect<double,3>(fact * 0.5 * mPlane->getDimensions()[0],sp_c_rel[1],0.0);
      result.mPoint1 = pl_pose.transformToGlobal(rim_pt);
    } else {
      // The sphere is outside one of the corners of the plane.
      vect<double,3> rim_pt = vect<double,3>(0.5 * mPlane->getDimensions()[0],0.5 * mPlane->getDimensions()[1],0.0);
      if(sp_c_rel[0] < 0.0)
        rim_pt[0] = -rim_pt[0];
      if(sp_c_rel[1] < 0.0)
        rim_pt[1] = -rim_pt[1];
      result.mPoint1 = pl_pose.transformToGlobal(rim_pt);
    };
    vect<double,3> diff = result.mPoint1 - sp_c;
    double diff_d = norm_2(diff);
    result.mPoint2 = sp_c + (mSphere->getRadius() / diff_d) * diff;
    result.mDistance = diff_d - mSphere->getRadius();
  };
  return result;
};*/

proximity_record_3D prox_plane_sphere::computeProximity(
    const shape_3D_precompute_pack& aPack1,
    const shape_3D_precompute_pack& aPack2) {
  if ((mSphere == nullptr) || (mPlane == nullptr)) {
    return {};
  }

  if (aPack1.parent == mPlane) {
    return compute_proximity(*mPlane, aPack1, *mSphere, aPack2);
  }
  return compute_proximity(*mSphere, aPack1, *mPlane, aPack2);
}

prox_plane_sphere::prox_plane_sphere(const plane* aPlane, const sphere* aSphere)
    : mPlane(aPlane), mSphere(aSphere) {}

}  // namespace ReaK::geom

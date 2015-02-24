
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

#include <ReaK/geometry/proximity/prox_plane_ccylinder.hpp>

#include <cmath>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


proximity_record_3D compute_proximity(const plane& aPlane, 
                                      const shape_3D_precompute_pack& aPack1,
                                      const capped_cylinder& aCCylinder, 
                                      const shape_3D_precompute_pack& aPack2) {
  using std::fabs; using std::sqrt; using ReaK::unit;
  proximity_record_3D result;
  
  const pose_3D<double>& cy_pose = aPack1.global_pose;
  const pose_3D<double>& pl_pose = aPack2.global_pose;
  
  vect<double,3> cy_c = cy_pose.Position;
  vect<double,3> cy_t = cy_pose.rotateToGlobal(vect<double,3>(0.0,0.0,1.0));
  
  vect<double,3> cy_c_rel = pl_pose.transformFromGlobal(cy_c);
  vect<double,3> cy_t_rel = pl_pose.rotateFromGlobal(cy_t);
  
  const double cy_len = aCCylinder.getLength();
  const double cy_rad = aCCylinder.getRadius();
  
  if(fabs(cy_t_rel[2]) < 1e-6) {
    // The capped-cylinder is sitting flat (on its side) on the plane.
    result.mPoint1 = pl_pose.transformToGlobal(vect<double,3>(cy_c_rel[0],cy_c_rel[1],0.0));
    result.mPoint2 = pl_pose.transformToGlobal(vect<double,3>(cy_c_rel[0],cy_c_rel[1],cy_c_rel[2] - cy_rad));
    result.mDistance = cy_c_rel[2] - cy_rad;
  } else {
    // The capped-cylinder is at an angle to the plane.
    if(cy_t_rel[2] > 0.0)
      cy_t_rel = -cy_t_rel;
    vect<double,3> cypt_rel = cy_c_rel + (0.5 * cy_len) * cy_t_rel + vect<double,3>(0.0,0.0,-cy_rad);
    result.mPoint1 = pl_pose.transformToGlobal(vect<double,3>(cypt_rel[0],cypt_rel[1],0.0));
    result.mPoint2 = pl_pose.transformToGlobal(cypt_rel);
    result.mDistance = cypt_rel[2];
  };
  return result;
};

proximity_record_3D compute_proximity(const capped_cylinder& aCCylinder, 
                                      const shape_3D_precompute_pack& aPack1,
                                      const plane& aPlane, 
                                      const shape_3D_precompute_pack& aPack2) {
  using std::swap;
  proximity_record_3D result = compute_proximity(aPlane, aPack2, aCCylinder, aPack1);
  swap(result.mPoint1,result.mPoint2);
  return result;
};

proximity_record_3D prox_plane_ccylinder::computeProximity(const shape_3D_precompute_pack& aPack1, 
                                                           const shape_3D_precompute_pack& aPack2) {
  if((!mCCylinder) || (!mPlane)) 
    return proximity_record_3D();
  
  if( aPack1.parent == mPlane ) 
    return compute_proximity(*mPlane,aPack1,*mCCylinder,aPack2);
  else
    return compute_proximity(*mCCylinder,aPack1,*mPlane,aPack2);
};


prox_plane_ccylinder::prox_plane_ccylinder(const plane* aPlane,
                                           const capped_cylinder* aCCylinder) :
                                           proximity_finder_3D(),
                                           mPlane(aPlane),
                                           mCCylinder(aCCylinder) { };


};

};


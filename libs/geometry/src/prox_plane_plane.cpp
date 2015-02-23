
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

#include <ReaK/geometry/proximity/prox_plane_plane.hpp>

#include <cmath>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


void compute_proximity_of_point(const plane& aPlane,
                                const pose_3D<double>& aPlGblPose,
                                const vect<double,3>& aPoint, 
                                vect<double,3>& aPointRec, 
                                double& aDistance) {
  using std::fabs; using std::sqrt;
  
  vect<double,3> pt_rel = aPlGblPose.transformFromGlobal(aPoint);
  
  if((pt_rel[0] > -0.5 * aPlane.getDimensions()[0]) &&
     (pt_rel[0] <  0.5 * aPlane.getDimensions()[0]) &&
     (pt_rel[1] > -0.5 * aPlane.getDimensions()[1]) &&
     (pt_rel[1] <  0.5 * aPlane.getDimensions()[1])) {
    // The sphere is within the dimensions of the plane.
    
    double fact = 1.0;
    if(pt_rel[2] < 0.0)
      fact = -1.0;
    
    aPointRec = aPlGblPose.transformToGlobal(vect<double,3>(pt_rel[0],pt_rel[1],0.0));
    aDistance = fact * pt_rel[2];
  } else {
    vect<double,3> rim_pt;
    if((pt_rel[0] > -0.5 * aPlane.getDimensions()[0]) &&
       (pt_rel[0] <  0.5 * aPlane.getDimensions()[0])) {
      // The sphere is on the above or below the plane (y-axis).
      double fact = 1.0;
      if(pt_rel[1] < 0.0)
        fact = -1.0;
      vect<double,3> rim_pt = vect<double,3>(pt_rel[0],fact * 0.5 * aPlane.getDimensions()[1],0.0);
      aPointRec = aPlGblPose.transformToGlobal(rim_pt);
    } else if((pt_rel[1] > -0.5 * aPlane.getDimensions()[1]) &&
              (pt_rel[1] <  0.5 * aPlane.getDimensions()[1])) {
      // The sphere is on the right or left of the plane (x-axis).
      double fact = 1.0;
      if(pt_rel[0] < 0.0)
        fact = -1.0;
      vect<double,3> rim_pt = vect<double,3>(fact * 0.5 * aPlane.getDimensions()[0],pt_rel[1],0.0);
      aPointRec = aPlGblPose.transformToGlobal(rim_pt);
    } else {
      // The sphere is outside one of the corners of the plane.
      vect<double,3> rim_pt = vect<double,3>(0.5 * aPlane.getDimensions()[0],
                                             0.5 * aPlane.getDimensions()[1],0.0);
      if(pt_rel[0] < 0.0)
        rim_pt[0] = -rim_pt[0];
      if(pt_rel[1] < 0.0)
        rim_pt[1] = -rim_pt[1];
      aPointRec = aPlGblPose.transformToGlobal(rim_pt);
    };
    aDistance = norm_2(aPointRec - aPoint);
  };
  
};

proximity_record_3D compute_proximity(const plane& aPlane1, 
                                      const shape_3D_precompute_pack& aPack1,
                                      const plane& aPlane2, 
                                      const shape_3D_precompute_pack& aPack2) {
  proximity_record_3D result;
  result.mDistance = std::numeric_limits<double>::infinity();
  
  const pose_3D<double>& p1_pose = aPack1.global_pose;
  const pose_3D<double>& p2_pose = aPack2.global_pose;
  
  const vect<double,2> p1_dim = aPlane1.getDimensions();
  const vect<double,2> p2_dim = aPlane2.getDimensions();
  
  vect<double,3> temp_pt;
  double temp_dist;
  
  vect<double,3> corner = vect<double,3>(0.5 * p2_dim[0],
                                         0.5 * p2_dim[1], 0.0);
  vect<double,3> corner_gbl = p2_pose.transformToGlobal(corner);
  compute_proximity_of_point(aPlane1, p1_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint1 = temp_pt;
    result.mPoint2 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = p2_pose.transformToGlobal(corner);
  compute_proximity_of_point(aPlane1, p1_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint1 = temp_pt;
    result.mPoint2 = corner_gbl;
  };
  
  corner[0] = -corner[0];
  corner_gbl = p2_pose.transformToGlobal(corner);
  compute_proximity_of_point(aPlane1, p1_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint1 = temp_pt;
    result.mPoint2 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = p2_pose.transformToGlobal(corner);
  compute_proximity_of_point(aPlane1, p1_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint1 = temp_pt;
    result.mPoint2 = corner_gbl;
  };
  
  
  corner = vect<double,3>(0.5 * p1_dim[0],
                          0.5 * p1_dim[1], 0.0);
  corner_gbl = p1_pose.transformToGlobal(corner);
  compute_proximity_of_point(aPlane2, p2_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint2 = temp_pt;
    result.mPoint1 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = p1_pose.transformToGlobal(corner);
  compute_proximity_of_point(aPlane2, p2_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint2 = temp_pt;
    result.mPoint1 = corner_gbl;
  };
  
  corner[0] = -corner[0];
  corner_gbl = p1_pose.transformToGlobal(corner);
  compute_proximity_of_point(aPlane2, p2_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint2 = temp_pt;
    result.mPoint1 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = p1_pose.transformToGlobal(corner);
  compute_proximity_of_point(aPlane2, p2_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < result.mDistance) {
    result.mDistance = temp_dist;
    result.mPoint2 = temp_pt;
    result.mPoint1 = corner_gbl;
  };
  
  return result;
};

    
void prox_plane_plane::computeProximity(const shape_3D_precompute_pack& aPack1, 
                                        const shape_3D_precompute_pack& aPack2) {
  mLastResult.mDistance = std::numeric_limits<double>::infinity();
  mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
  mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
  
  if((!mPlane1) || (!mPlane2))
    return;
  
  if( aPack1.parent == mPlane1 ) 
    mLastResult = compute_proximity(*mPlane1,aPack1,*mPlane2,aPack2);
  else
    mLastResult = compute_proximity(*mPlane2,aPack1,*mPlane1,aPack2);
  
};


prox_plane_plane::prox_plane_plane(const plane* aPlane1,
                                   const plane* aPlane2) :
                                   proximity_finder_3D(),
                                   mPlane1(aPlane1),
                                   mPlane2(aPlane2) { };


};

};


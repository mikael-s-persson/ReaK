
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

#include <ReaK/geometry/proximity/prox_rectangle_rectangle.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


void prox_rectangle_rectangle::computeProximityOfPoint(const rectangle& aRectangle, 
                                                       const pose_2D<double>& aRecGblPose,
                                                       const vect<double,2>& aPoint, 
                                                       vect<double,2>& aPointRec, 
                                                       double& aDistance) {
  using std::fabs;
  
  vect<double,2> pt_rel = aRecGblPose.transformFromGlobal(aPoint);
  
  bool in_x_range = ((pt_rel[0] > -0.5 * aRectangle.getDimensions()[0]) &&
                     (pt_rel[0] <  0.5 * aRectangle.getDimensions()[0]));
  bool in_y_range = ((pt_rel[1] > -0.5 * aRectangle.getDimensions()[1]) &&
                     (pt_rel[1] <  0.5 * aRectangle.getDimensions()[1]));
  
  if(in_x_range && in_y_range) {
    // The circle is inside the rectangle.
    vect<double,2> bound_dists = vect<double,2>(0.5 * aRectangle.getDimensions()[0] - fabs(pt_rel[0]),
                                                0.5 * aRectangle.getDimensions()[1] - fabs(pt_rel[1]));
    if(bound_dists[0] <= bound_dists[1]) {
      in_x_range = false;
    } else {
      in_y_range = false;
    };
  };
  
  vect<double,2> corner_pt = 0.5 * aRectangle.getDimensions();
  if(in_x_range)
    corner_pt[0] = pt_rel[0];
  else if(pt_rel[0] < 0.0)
    corner_pt[0] = -corner_pt[0];
  if(in_y_range)
    corner_pt[1] = pt_rel[1];
  else if(pt_rel[1] < 0.0)
    corner_pt[1] = -corner_pt[1];
  aPointRec = aRecGblPose.transformToGlobal(corner_pt);
  aDistance = norm_2(aPointRec - aPoint);
};
    
void prox_rectangle_rectangle::computeProximity(const shape_2D_precompute_pack& aPack1, 
                                                const shape_2D_precompute_pack& aPack2) {
  mLastResult.mDistance = std::numeric_limits<double>::infinity();
  mLastResult.mPoint1 = vect<double,2>(0.0,0.0);
  mLastResult.mPoint2 = vect<double,2>(0.0,0.0);
  if((!mRectangle1) || (!mRectangle2))
    return;
  
  const pose_2D<double>& r1_pose = (aPack1.parent == mRectangle1 ? 
                                    aPack1.global_pose : aPack2.global_pose);
  const pose_2D<double>& r2_pose = (aPack1.parent == mRectangle1 ? 
                                    aPack2.global_pose : aPack1.global_pose);
  
  vect<double,2> temp_pt;
  double temp_dist;
  
  vect<double,2> corner = 0.5 * mRectangle2->getDimensions();
  vect<double,2> corner_gbl = r2_pose.transformToGlobal(corner);
  computeProximityOfPoint(*mRectangle1, r1_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint1 = temp_pt;
    mLastResult.mPoint2 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = r2_pose.transformToGlobal(corner);
  computeProximityOfPoint(*mRectangle1, r1_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint1 = temp_pt;
    mLastResult.mPoint2 = corner_gbl;
  };
  
  corner[0] = -corner[0];
  corner_gbl = r2_pose.transformToGlobal(corner);
  computeProximityOfPoint(*mRectangle1, r1_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint1 = temp_pt;
    mLastResult.mPoint2 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = r2_pose.transformToGlobal(corner);
  computeProximityOfPoint(*mRectangle1, r1_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint1 = temp_pt;
    mLastResult.mPoint2 = corner_gbl;
  };
  
  
  corner = 0.5 * mRectangle1->getDimensions();
  corner_gbl = r1_pose.transformToGlobal(corner);
  computeProximityOfPoint(*mRectangle2, r2_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint2 = temp_pt;
    mLastResult.mPoint1 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = r1_pose.transformToGlobal(corner);
  computeProximityOfPoint(*mRectangle2, r2_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint2 = temp_pt;
    mLastResult.mPoint1 = corner_gbl;
  };
  
  corner[0] = -corner[0];
  corner_gbl = r1_pose.transformToGlobal(corner);
  computeProximityOfPoint(*mRectangle2, r2_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint2 = temp_pt;
    mLastResult.mPoint1 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = r1_pose.transformToGlobal(corner);
  computeProximityOfPoint(*mRectangle2, r2_pose, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint2 = temp_pt;
    mLastResult.mPoint1 = corner_gbl;
  };
  
};


prox_rectangle_rectangle::prox_rectangle_rectangle(const rectangle* aRectangle1,
                                                   const rectangle* aRectangle2) :
                                                   proximity_finder_2D(),
                                                   mRectangle1(aRectangle1),
                                                   mRectangle2(aRectangle2) { };


};

};


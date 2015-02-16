
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


shared_ptr< shape_3D > prox_plane_plane::getShape1() const {
  return mPlane1;
};

shared_ptr< shape_3D > prox_plane_plane::getShape2() const {
  return mPlane2;
};


void prox_plane_plane::computeProximityOfPoint(const shared_ptr< plane >& aPlane, 
                                               const vect<double,3>& aPoint, 
                                               vect<double,3>& aPointRec, 
                                               double& aDistance) {
  using std::fabs; using std::sqrt;
  
  vect<double,3> pt_rel = aPlane->getPose().transformFromGlobal(aPoint);
  
  if((pt_rel[0] > -0.5 * aPlane->getDimensions()[0]) &&
     (pt_rel[0] <  0.5 * aPlane->getDimensions()[0]) &&
     (pt_rel[1] > -0.5 * aPlane->getDimensions()[1]) &&
     (pt_rel[1] <  0.5 * aPlane->getDimensions()[1])) {
    // The sphere is within the dimensions of the plane.
    
    double fact = 1.0;
    if(pt_rel[2] < 0.0)
      fact = -1.0;
    
    aPointRec = aPlane->getPose().transformToGlobal(vect<double,3>(pt_rel[0],pt_rel[1],0.0));
    aDistance = fact * pt_rel[2];
  } else {
    vect<double,3> rim_pt;
    if((pt_rel[0] > -0.5 * aPlane->getDimensions()[0]) &&
       (pt_rel[0] <  0.5 * aPlane->getDimensions()[0])) {
      // The sphere is on the above or below the plane (y-axis).
      double fact = 1.0;
      if(pt_rel[1] < 0.0)
        fact = -1.0;
      vect<double,3> rim_pt = vect<double,3>(pt_rel[0],fact * 0.5 * aPlane->getDimensions()[1],0.0);
      aPointRec = aPlane->getPose().transformToGlobal(rim_pt);
    } else if((pt_rel[1] > -0.5 * aPlane->getDimensions()[1]) &&
              (pt_rel[1] <  0.5 * aPlane->getDimensions()[1])) {
      // The sphere is on the right or left of the plane (x-axis).
      double fact = 1.0;
      if(pt_rel[0] < 0.0)
        fact = -1.0;
      vect<double,3> rim_pt = vect<double,3>(fact * 0.5 * aPlane->getDimensions()[0],pt_rel[1],0.0);
      aPointRec = aPlane->getPose().transformToGlobal(rim_pt);
    } else {
      // The sphere is outside one of the corners of the plane.
      vect<double,3> rim_pt = vect<double,3>(0.5 * aPlane->getDimensions()[0],
                                             0.5 * aPlane->getDimensions()[1],0.0);
      if(pt_rel[0] < 0.0)
        rim_pt[0] = -rim_pt[0];
      if(pt_rel[1] < 0.0)
        rim_pt[1] = -rim_pt[1];
      aPointRec = aPlane->getPose().transformToGlobal(rim_pt);
    };
    aDistance = norm_2(aPointRec - aPoint);
  };
  
};

    
void prox_plane_plane::computeProximity() {
  mLastResult.mDistance = std::numeric_limits<double>::infinity();
  mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
  mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
  
  if((!mPlane1) || (!mPlane2))
    return;
  
  vect<double,3> temp_pt;
  double temp_dist;
  
  vect<double,3> corner = vect<double,3>(0.5 * mPlane2->getDimensions()[0],
                                         0.5 * mPlane2->getDimensions()[1], 0.0);
  vect<double,3> corner_gbl = mPlane2->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mPlane1, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint1 = temp_pt;
    mLastResult.mPoint2 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = mPlane2->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mPlane1, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint1 = temp_pt;
    mLastResult.mPoint2 = corner_gbl;
  };
  
  corner[0] = -corner[0];
  corner_gbl = mPlane2->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mPlane1, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint1 = temp_pt;
    mLastResult.mPoint2 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = mPlane2->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mPlane1, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint1 = temp_pt;
    mLastResult.mPoint2 = corner_gbl;
  };
  
  
  corner = vect<double,3>(0.5 * mPlane1->getDimensions()[0],
                          0.5 * mPlane1->getDimensions()[1], 0.0);
  corner_gbl = mPlane1->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mPlane2, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint2 = temp_pt;
    mLastResult.mPoint1 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = mPlane1->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mPlane2, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint2 = temp_pt;
    mLastResult.mPoint1 = corner_gbl;
  };
  
  corner[0] = -corner[0];
  corner_gbl = mPlane1->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mPlane2, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint2 = temp_pt;
    mLastResult.mPoint1 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = mPlane1->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mPlane2, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint2 = temp_pt;
    mLastResult.mPoint1 = corner_gbl;
  };
  
};


prox_plane_plane::prox_plane_plane(const shared_ptr< plane >& aPlane1,
                                   const shared_ptr< plane >& aPlane2) :
                                   proximity_finder_3D(),
                                   mPlane1(aPlane1),
                                   mPlane2(aPlane2) { };
    
    
void RK_CALL prox_plane_plane::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_3D::save(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mPlane1)
    & RK_SERIAL_SAVE_WITH_NAME(mPlane2);
};
    
void RK_CALL prox_plane_plane::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_3D::load(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mPlane1)
    & RK_SERIAL_LOAD_WITH_NAME(mPlane2);
};


};

};











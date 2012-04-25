
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

#include "prox_rectangle_rectangle.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_2D > prox_rectangle_rectangle::getShape1() const {
  return mRectangle1;
};

shared_ptr< shape_2D > prox_rectangle_rectangle::getShape2() const {
  return mRectangle2;
};


void prox_rectangle_rectangle::computeProximityOfPoint(const shared_ptr< rectangle >& aRectangle, 
						       const vect<double,2>& aPoint, 
						       vect<double,2>& aPointRec, 
						       double& aDistance) {
  using std::fabs;
  
  vect<double,2> pt_rel = aRectangle->getPose().transformFromGlobal(aPoint);
  
  bool in_x_range = ((pt_rel[0] > -0.5 * aRectangle->getDimensions()[0]) &&
                     (pt_rel[0] <  0.5 * aRectangle->getDimensions()[0]));
  bool in_y_range = ((pt_rel[1] > -0.5 * aRectangle->getDimensions()[1]) &&
                     (pt_rel[1] <  0.5 * aRectangle->getDimensions()[1]));
  
  if(in_x_range && in_y_range) {
    // The circle is inside the rectangle.
    vect<double,2> bound_dists = vect<double,2>(0.5 * aRectangle->getDimensions()[0] - fabs(pt_rel[0]),
						0.5 * aRectangle->getDimensions()[1] - fabs(pt_rel[1]));
    if(bound_dists[0] <= bound_dists[1]) {
      in_x_range = false;
    } else {
      in_y_range = false;
    };
  };
  
  vect<double,2> corner_pt = 0.5 * aRectangle->getDimensions();
  if(in_x_range)
    corner_pt[0] = pt_rel[0];
  else if(pt_rel[0] < 0.0)
    corner_pt[0] = -corner_pt[0];
  if(in_y_range)
    corner_pt[1] = pt_rel[1];
  else if(pt_rel[1] < 0.0)
    corner_pt[1] = -corner_pt[1];
  aPointRec = aRectangle->getPose().transformToGlobal(corner_pt);
  aDistance = norm_2(aPointRec - aPoint);
};
    
void prox_rectangle_rectangle::computeProximity() {
  mLastResult.mDistance = std::numeric_limits<double>::infinity();
  mLastResult.mPoint1 = vect<double,2>(0.0,0.0);
  mLastResult.mPoint2 = vect<double,2>(0.0,0.0);
  if((!mRectangle1) || (!mRectangle2))
    return;
  
  vect<double,2> temp_pt;
  double temp_dist;
  
  vect<double,2> corner = 0.5 * mRectangle2->getDimensions();
  vect<double,2> corner_gbl = mRectangle2->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mRectangle1, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint1 = temp_pt;
    mLastResult.mPoint2 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = mRectangle2->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mRectangle1, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint1 = temp_pt;
    mLastResult.mPoint2 = corner_gbl;
  };
  
  corner[0] = -corner[0];
  corner_gbl = mRectangle2->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mRectangle1, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint1 = temp_pt;
    mLastResult.mPoint2 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = mRectangle2->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mRectangle1, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint1 = temp_pt;
    mLastResult.mPoint2 = corner_gbl;
  };
  
  
  corner = 0.5 * mRectangle1->getDimensions();
  corner_gbl = mRectangle1->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mRectangle2, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint2 = temp_pt;
    mLastResult.mPoint1 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = mRectangle1->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mRectangle2, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint2 = temp_pt;
    mLastResult.mPoint1 = corner_gbl;
  };
  
  corner[0] = -corner[0];
  corner_gbl = mRectangle1->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mRectangle2, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint2 = temp_pt;
    mLastResult.mPoint1 = corner_gbl;
  };
  
  corner[1] = -corner[1];
  corner_gbl = mRectangle1->getPose().transformToGlobal(corner);
  computeProximityOfPoint(mRectangle2, corner_gbl, temp_pt, temp_dist);
  if(temp_dist < mLastResult.mDistance) {
    mLastResult.mDistance = temp_dist;
    mLastResult.mPoint2 = temp_pt;
    mLastResult.mPoint1 = corner_gbl;
  };
  
};


prox_rectangle_rectangle::prox_rectangle_rectangle(const shared_ptr< rectangle >& aRectangle1,
						   const shared_ptr< rectangle >& aRectangle2) :
						   proximity_finder_2D(),
						   mRectangle1(aRectangle1),
						   mRectangle2(aRectangle2) { };
    
    
void RK_CALL prox_rectangle_rectangle::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_2D::save(A,proximity_finder_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mRectangle1)
    & RK_SERIAL_SAVE_WITH_NAME(mRectangle2);
};
    
void RK_CALL prox_rectangle_rectangle::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_2D::load(A,proximity_finder_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mRectangle1)
    & RK_SERIAL_LOAD_WITH_NAME(mRectangle2);
};


};

};











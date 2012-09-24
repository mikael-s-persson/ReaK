
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

#include "prox_plane_sphere.hpp"

#include <cmath>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_3D > prox_plane_sphere::getShape1() const {
  return mPlane;
};

shared_ptr< shape_3D > prox_plane_sphere::getShape2() const {
  return mSphere;
};
    
/*
// this version assumes a finite plane for proximity purposes.
void prox_plane_sphere::computeProximity() {
  if((!mSphere) || (!mPlane)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  using std::fabs; using std::sqrt;
  
  vect<double,3> sp_c = mSphere->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  vect<double,3> pl_c = mPlane->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  
  vect<double,3> sp_c_rel = mPlane->getPose().transformFromGlobal(sp_c);
  
  if((sp_c_rel[0] > -0.5 * mPlane->getDimensions()[0]) &&
     (sp_c_rel[0] <  0.5 * mPlane->getDimensions()[0]) &&
     (sp_c_rel[1] > -0.5 * mPlane->getDimensions()[1]) &&
     (sp_c_rel[1] <  0.5 * mPlane->getDimensions()[1])) {
    // The sphere is within the dimensions of the plane.
    
    double fact = 1.0;
    if(sp_c_rel[2] < 0.0)
      fact = -1.0;
    
    mLastResult.mPoint1 = mPlane->getPose().transformToGlobal(vect<double,3>(sp_c_rel[0],sp_c_rel[1],0.0));
    mLastResult.mPoint2 = mPlane->getPose().transformToGlobal(vect<double,3>(sp_c_rel[0],sp_c_rel[1],sp_c_rel[2] - fact * mSphere->getRadius()));
    mLastResult.mDistance = fact * sp_c_rel[2] - mSphere->getRadius();
  } else {
    vect<double,3> rim_pt;
    if((sp_c_rel[0] > -0.5 * mPlane->getDimensions()[0]) &&
       (sp_c_rel[0] <  0.5 * mPlane->getDimensions()[0])) {
      // The sphere is on the above or below the plane (y-axis).
      double fact = 1.0;
      if(sp_c_rel[1] < 0.0)
	fact = -1.0;
      vect<double,3> rim_pt = vect<double,3>(sp_c_rel[0],fact * 0.5 * mPlane->getDimensions()[1],0.0);
      mLastResult.mPoint1 = mPlane->getPose().transformToGlobal(rim_pt);
    } else if((sp_c_rel[1] > -0.5 * mPlane->getDimensions()[1]) &&
              (sp_c_rel[1] <  0.5 * mPlane->getDimensions()[1])) {
      // The sphere is on the right or left of the plane (x-axis).
      double fact = 1.0;
      if(sp_c_rel[0] < 0.0)
	fact = -1.0;
      vect<double,3> rim_pt = vect<double,3>(fact * 0.5 * mPlane->getDimensions()[0],sp_c_rel[1],0.0);
      mLastResult.mPoint1 = mPlane->getPose().transformToGlobal(rim_pt);
    } else {
      // The sphere is outside one of the corners of the plane.
      vect<double,3> rim_pt = vect<double,3>(0.5 * mPlane->getDimensions()[0],0.5 * mPlane->getDimensions()[1],0.0);
      if(sp_c_rel[0] < 0.0)
	rim_pt[0] = -rim_pt[0];
      if(sp_c_rel[1] < 0.0)
	rim_pt[1] = -rim_pt[1];
      mLastResult.mPoint1 = mPlane->getPose().transformToGlobal(rim_pt);
    };
    vect<double,3> diff = mLastResult.mPoint1 - sp_c;
    double diff_d = norm_2(diff);
    mLastResult.mPoint2 = sp_c + (mSphere->getRadius() / diff_d) * diff;
    mLastResult.mDistance = diff_d - mSphere->getRadius();
  };
};*/

void prox_plane_sphere::computeProximity() {
  if((!mSphere) || (!mPlane)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  
  vect<double,3> sp_c = mSphere->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  vect<double,3> pl_c = mPlane->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  
  vect<double,3> sp_c_rel = mPlane->getPose().transformFromGlobal(sp_c);
  
  mLastResult.mPoint1 = mPlane->getPose().transformToGlobal(vect<double,3>(sp_c_rel[0],sp_c_rel[1],0.0));
  mLastResult.mPoint2 = mPlane->getPose().transformToGlobal(vect<double,3>(sp_c_rel[0],sp_c_rel[1],sp_c_rel[2] - mSphere->getRadius()));
  mLastResult.mDistance = sp_c_rel[2] - mSphere->getRadius();
};


prox_plane_sphere::prox_plane_sphere(const shared_ptr< plane >& aPlane,
				     const shared_ptr< sphere >& aSphere) :
				     proximity_finder_3D(),
				     mPlane(aPlane),
				     mSphere(aSphere) { };
    
    
void RK_CALL prox_plane_sphere::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_3D::save(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mPlane)
    & RK_SERIAL_SAVE_WITH_NAME(mSphere);
};
    
void RK_CALL prox_plane_sphere::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_3D::load(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mPlane)
    & RK_SERIAL_LOAD_WITH_NAME(mSphere);
};


};

};











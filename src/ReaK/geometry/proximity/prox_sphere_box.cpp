
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

#include <ReaK/geometry/proximity/prox_sphere_box.hpp>

#include <ReaK/geometry/proximity/prox_fundamentals_3D.hpp>

#include <cmath>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_3D > prox_sphere_box::getShape1() const {
  return mSphere;
};

shared_ptr< shape_3D > prox_sphere_box::getShape2() const {
  return mBox;
};
    
void prox_sphere_box::computeProximity() {
  if((!mSphere) || (!mBox)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  using std::fabs; using std::sqrt;
  
  vect<double,3> sp_c = mSphere->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  vect<double,3> bx_c = mBox->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  
  proximity_record_3D bxpt_result = findProximityBoxToPoint(mBox, sp_c);
  
  // add a sphere-sweep around the point-box solution.
  vect<double,3> diff_v = bxpt_result.mPoint1 - bxpt_result.mPoint2;
  double diff_d = norm_2(diff_v);
  if(bxpt_result.mDistance < 0.0)
    mLastResult.mPoint1 = bxpt_result.mPoint2 - (mSphere->getRadius() / diff_d) * diff_v;
  else
    mLastResult.mPoint1 = bxpt_result.mPoint2 + (mSphere->getRadius() / diff_d) * diff_v;
  mLastResult.mPoint2 = bxpt_result.mPoint1;
  mLastResult.mDistance = bxpt_result.mDistance - mSphere->getRadius();
  return;
};


prox_sphere_box::prox_sphere_box(const shared_ptr< sphere >& aSphere,
                                 const shared_ptr< box >& aBox) :
                                 proximity_finder_3D(),
                                 mSphere(aSphere),
                                 mBox(aBox) { };
    
    
void RK_CALL prox_sphere_box::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_3D::save(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mSphere)
    & RK_SERIAL_SAVE_WITH_NAME(mBox);
};
    
void RK_CALL prox_sphere_box::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_3D::load(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mSphere)
    & RK_SERIAL_LOAD_WITH_NAME(mBox);
};


};

};











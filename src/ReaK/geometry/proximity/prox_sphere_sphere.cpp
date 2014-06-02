
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

#include <ReaK/geometry/proximity/prox_sphere_sphere.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_3D > prox_sphere_sphere::getShape1() const {
  return mSphere1;
};

shared_ptr< shape_3D > prox_sphere_sphere::getShape2() const {
  return mSphere2;
};
    
void prox_sphere_sphere::computeProximity() {
  if((!mSphere1) || (!mSphere2)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  vect<double,3> c1 = mSphere1->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  vect<double,3> c2 = mSphere2->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  
  vect<double,3> diff_cc = c2 - c1;
  double dist_cc = norm_2(diff_cc);
  
  mLastResult.mDistance = dist_cc - mSphere1->getRadius() - mSphere2->getRadius();
  mLastResult.mPoint1 = c1 + (mSphere1->getRadius() / dist_cc) * diff_cc;
  mLastResult.mPoint2 = c2 - (mSphere2->getRadius() / dist_cc) * diff_cc;
  
};


prox_sphere_sphere::prox_sphere_sphere(const shared_ptr< sphere >& aSphere1,
                                       const shared_ptr< sphere >& aSphere2) :
                                       proximity_finder_3D(),
                                       mSphere1(aSphere1),
                                       mSphere2(aSphere2) { };
    
    
void RK_CALL prox_sphere_sphere::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_3D::save(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mSphere1)
    & RK_SERIAL_SAVE_WITH_NAME(mSphere2);
};
    
void RK_CALL prox_sphere_sphere::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_3D::load(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mSphere1)
    & RK_SERIAL_LOAD_WITH_NAME(mSphere2);
};


};

};











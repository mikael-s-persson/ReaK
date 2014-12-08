
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

#include <ReaK/geometry/proximity/prox_ccylinder_cylinder.hpp>

#include <ReaK/geometry/proximity/prox_fundamentals_3D.hpp>

#include <cmath>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_3D > prox_ccylinder_cylinder::getShape1() const {
  return mCCylinder;
};

shared_ptr< shape_3D > prox_ccylinder_cylinder::getShape2() const {
  return mCylinder;
};

void prox_ccylinder_cylinder::computeProximity() {
  if((!mCCylinder) || (!mCylinder)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  
  mLastResult = findProximityByGJKEPA(
    ccylinder_support_func(mCCylinder), 
    cylinder_support_func(mCylinder));
  
};


prox_ccylinder_cylinder::prox_ccylinder_cylinder(const shared_ptr< capped_cylinder >& aCCylinder,
                                                 const shared_ptr< cylinder >& aCylinder) :
                                                 proximity_finder_3D(),
                                                 mCCylinder(aCCylinder),
                                                 mCylinder(aCylinder) { };


void RK_CALL prox_ccylinder_cylinder::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_3D::save(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mCCylinder)
    & RK_SERIAL_SAVE_WITH_NAME(mCylinder);
};

void RK_CALL prox_ccylinder_cylinder::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_3D::load(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mCCylinder)
    & RK_SERIAL_LOAD_WITH_NAME(mCylinder);
};


};

};











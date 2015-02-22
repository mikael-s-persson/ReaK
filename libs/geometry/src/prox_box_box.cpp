
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

#include <ReaK/geometry/proximity/prox_box_box.hpp>

#include <ReaK/geometry/proximity/prox_fundamentals_3D.hpp>


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_3D > prox_box_box::getShape1() const {
  return mBox1;
};

shared_ptr< shape_3D > prox_box_box::getShape2() const {
  return mBox2;
};



void prox_box_box::computeProximity(const shape_3D_precompute_pack& aPack1, 
                                    const shape_3D_precompute_pack& aPack2) {
  if((!mBox1) || (!mBox2)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  
  const pose_3D<double>& b1_pose = (aPack1.parent == mBox1.get() ? 
                                    aPack1.global_pose : aPack2.global_pose);
  const pose_3D<double>& b2_pose = (aPack1.parent == mBox1.get() ? 
                                    aPack2.global_pose : aPack1.global_pose);
  
  mLastResult = findProximityByGJKEPA(
    box_support_func(*mBox1, b1_pose), 
    box_support_func(*mBox2, b2_pose));
  
};


prox_box_box::prox_box_box(const shared_ptr< box >& aBox1,
                           const shared_ptr< box >& aBox2) :
                           proximity_finder_3D(),
                           mBox1(aBox1),
                           mBox2(aBox2) { };
    
    
void RK_CALL prox_box_box::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_3D::save(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mBox1)
    & RK_SERIAL_SAVE_WITH_NAME(mBox2);
};
    
void RK_CALL prox_box_box::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_3D::load(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mBox1)
    & RK_SERIAL_LOAD_WITH_NAME(mBox2);
};



};

};











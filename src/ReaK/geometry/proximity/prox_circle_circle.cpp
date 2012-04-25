
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

#include "prox_circle_circle.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_2D > prox_circle_circle::getShape1() const {
  return mCircle1;
};

shared_ptr< shape_2D > prox_circle_circle::getShape2() const {
  return mCircle2;
};
    
void prox_circle_circle::computeProximity() {
  if((!mCircle1) || (!mCircle2)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,2>(0.0,0.0);
    mLastResult.mPoint2 = vect<double,2>(0.0,0.0);
    return;
  };
  vect<double,2> c1 = mCircle1->getPose().transformToGlobal(vect<double,2>(0.0,0.0));
  vect<double,2> c2 = mCircle2->getPose().transformToGlobal(vect<double,2>(0.0,0.0));
  
  vect<double,2> diff_cc = c2 - c1;
  double dist_cc = norm_2(diff_cc);
  
  mLastResult.mDistance = dist_cc - mCircle1->getRadius() - mCircle2->getRadius();
  mLastResult.mPoint1 = c1 + (mCircle1->getRadius() / dist_cc) * diff_cc;
  mLastResult.mPoint2 = c2 - (mCircle2->getRadius() / dist_cc) * diff_cc;
  
};


prox_circle_circle::prox_circle_circle(const shared_ptr< circle >& aCircle1,
				       const shared_ptr< circle >& aCircle2) :
				       proximity_finder_2D(),
				       mCircle1(aCircle1),
				       mCircle2(aCircle2) { };
    
    
void RK_CALL prox_circle_circle::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_2D::save(A,proximity_finder_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mCircle1)
    & RK_SERIAL_SAVE_WITH_NAME(mCircle2);
};
    
void RK_CALL prox_circle_circle::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_2D::load(A,proximity_finder_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mCircle1)
    & RK_SERIAL_LOAD_WITH_NAME(mCircle2);
};


};

};











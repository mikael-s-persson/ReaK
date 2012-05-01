
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

#include "prox_circle_crect.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_2D > prox_circle_crect::getShape1() const {
  return mCircle;
};

shared_ptr< shape_2D > prox_circle_crect::getShape2() const {
  return mCRect;
};
    
void prox_circle_crect::computeProximity() {
  if((!mCircle) || (!mCRect)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,2>(0.0,0.0);
    mLastResult.mPoint2 = vect<double,2>(0.0,0.0);
    return;
  };
  
  using std::fabs;
  
  vect<double,2> ci_c = mCircle->getPose().transformToGlobal(vect<double,2>(0.0,0.0));
  vect<double,2> re_c = mCRect->getPose().transformToGlobal(vect<double,2>(0.0,0.0));
  
  vect<double,2> ci_c_rel = mCRect->getPose().transformFromGlobal(ci_c);
  
  bool in_x_range = ((ci_c_rel[0] > -0.5 * mCRect->getDimensions()[0]) &&
                     (ci_c_rel[0] <  0.5 * mCRect->getDimensions()[0]));
  bool in_y_range = ((ci_c_rel[1] > -0.5 * mCRect->getDimensions()[1]) &&
                     (ci_c_rel[1] <  0.5 * mCRect->getDimensions()[1]));
  
  if(in_x_range && in_y_range) {
    // The circle is inside the rectangle.
    vect<double,2> bound_dists = vect<double,2>(0.5 * mCRect->getDimensions()[0] - fabs(ci_c_rel[0]),
						0.5 * mCRect->getDimensions()[1] - fabs(ci_c_rel[1]));
    if(bound_dists[0] <= bound_dists[1]) {
      in_x_range = false;
    } else {
      in_y_range = false;
    };
  };
  
  vect<double,2> corner_pt = 0.5 * mCRect->getDimensions();
  if(in_x_range)
    corner_pt[0] = ci_c_rel[0];
  else if(ci_c_rel[0] < 0.0)
    corner_pt[0] = -corner_pt[0];
  if(in_y_range)
    corner_pt[1] = ci_c_rel[1];
  else if(ci_c_rel[1] < 0.0)
    corner_pt[1] = -corner_pt[1];
  mLastResult.mPoint2 = mCRect->getPose().transformToGlobal(corner_pt);
  vect<double,2> diff_v = mLastResult.mPoint2 - ci_c;
  double diff_d = norm_2(diff_v);
  mLastResult.mPoint1 = ci_c + (mCircle->getRadius() / diff_d) * diff_v;
  mLastResult.mDistance = diff_d - mCircle->getRadius();
};


prox_circle_crect::prox_circle_crect(const shared_ptr< circle >& aCircle,
				     const shared_ptr< capped_rectangle >& aCRect) :
				     proximity_finder_2D(),
				     mCircle(aCircle),
				     mCRect(aCRect) { };
    
    
void RK_CALL prox_circle_crect::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_2D::save(A,proximity_finder_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mCircle)
    & RK_SERIAL_SAVE_WITH_NAME(mCRect);
};
    
void RK_CALL prox_circle_crect::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_2D::load(A,proximity_finder_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mCircle)
    & RK_SERIAL_LOAD_WITH_NAME(mCRect);
};


};

};












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

#include <ReaK/geometry/proximity/prox_circle_rectangle.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_2D > prox_circle_rectangle::getShape1() const {
  return mCircle;
};

shared_ptr< shape_2D > prox_circle_rectangle::getShape2() const {
  return mRectangle;
};
    
void prox_circle_rectangle::computeProximity(const shape_2D_precompute_pack& aPack1, 
                                             const shape_2D_precompute_pack& aPack2) {
  if((!mCircle) || (!mRectangle)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,2>(0.0,0.0);
    mLastResult.mPoint2 = vect<double,2>(0.0,0.0);
    return;
  };
  
  using std::fabs;
  
  const pose_2D<double>& ci_pose = (aPack1.parent == mCircle.get() ? 
                                    aPack1.global_pose : aPack2.global_pose);
  const pose_2D<double>& re_pose = (aPack1.parent == mCircle.get() ? 
                                    aPack2.global_pose : aPack1.global_pose);
  
  vect<double,2> ci_c = ci_pose.Position;
  
  vect<double,2> ci_c_rel = re_pose.transformFromGlobal(ci_c);
  
  bool in_x_range = ((ci_c_rel[0] > -0.5 * mRectangle->getDimensions()[0]) &&
                     (ci_c_rel[0] <  0.5 * mRectangle->getDimensions()[0]));
  bool in_y_range = ((ci_c_rel[1] > -0.5 * mRectangle->getDimensions()[1]) &&
                     (ci_c_rel[1] <  0.5 * mRectangle->getDimensions()[1]));
  
  if(in_x_range && in_y_range) {
    // The circle is inside the rectangle.
    vect<double,2> bound_dists = vect<double,2>(0.5 * mRectangle->getDimensions()[0] - fabs(ci_c_rel[0]),
                                                0.5 * mRectangle->getDimensions()[1] - fabs(ci_c_rel[1]));
    if(bound_dists[0] <= bound_dists[1]) {
      in_x_range = false;
    } else {
      in_y_range = false;
    };
  };
  
  vect<double,2> corner_pt = 0.5 * mRectangle->getDimensions();
  if(in_x_range)
    corner_pt[0] = ci_c_rel[0];
  else if(ci_c_rel[0] < 0.0)
    corner_pt[0] = -corner_pt[0];
  if(in_y_range)
    corner_pt[1] = ci_c_rel[1];
  else if(ci_c_rel[1] < 0.0)
    corner_pt[1] = -corner_pt[1];
  mLastResult.mPoint2 = re_pose.transformToGlobal(corner_pt);
  vect<double,2> diff_v = mLastResult.mPoint2 - ci_c;
  double diff_d = norm_2(diff_v);
  mLastResult.mPoint1 = ci_c + (mCircle->getRadius() / diff_d) * diff_v;
  mLastResult.mDistance = diff_d - mCircle->getRadius();
};


prox_circle_rectangle::prox_circle_rectangle(const shared_ptr< circle >& aCircle,
                                             const shared_ptr< rectangle >& aRectangle) :
                                             proximity_finder_2D(),
                                             mCircle(aCircle),
                                             mRectangle(aRectangle) { };
    
    
void RK_CALL prox_circle_rectangle::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_2D::save(A,proximity_finder_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mCircle)
    & RK_SERIAL_SAVE_WITH_NAME(mRectangle);
};
    
void RK_CALL prox_circle_rectangle::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_2D::load(A,proximity_finder_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mCircle)
    & RK_SERIAL_LOAD_WITH_NAME(mRectangle);
};


};

};











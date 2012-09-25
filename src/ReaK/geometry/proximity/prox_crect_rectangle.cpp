
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

#include "prox_crect_rectangle.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_2D > prox_crect_rectangle::getShape1() const {
  return mCRect;
};

shared_ptr< shape_2D > prox_crect_rectangle::getShape2() const {
  return mRectangle;
};


void prox_crect_rectangle::computeProximityOfLine(const shared_ptr< rectangle >& aRectangle, const vect<double,2>& ln_c, const vect<double,2>& ln_t, double half_length, proximity_record_2D& result) {
  
  using std::fabs;
  
  vect<double,2> ln_c_rel = aRectangle->getPose().transformFromGlobal(ln_c);
  vect<double,2> ln_t_rel = aRectangle->getPose().rotateFromGlobal(ln_t);
  
  if(fabs(ln_t_rel[0]) < 1e-5) {
    // this means the line is vertical.
    if((ln_c_rel[1] + half_length > -0.5 * aRectangle->getDimensions()[1]) || 
       (ln_c_rel[1] - half_length <  0.5 * aRectangle->getDimensions()[1])) {
      // there is an overlap between the rectangle side and the line.
      double max_y_rel = ((ln_c_rel[1] + half_length <  0.5 * aRectangle->getDimensions()[1]) ? (ln_c_rel[1] + half_length) : ( 0.5 * aRectangle->getDimensions()[1]));
      double min_y_rel = ((ln_c_rel[1] - half_length > -0.5 * aRectangle->getDimensions()[1]) ? (ln_c_rel[1] - half_length) : (-0.5 * aRectangle->getDimensions()[1]));
      double avg_y_rel = (max_y_rel + min_y_rel) * 0.5;
      vect<double,2> ln_r_rel(1.0,0.0);
      if(ln_c_rel[0] < 0.0)
        ln_r_rel[0] = -1.0;
      result.mPoint1 = aRectangle->getPose().transformToGlobal(vect<double,2>(ln_c_rel[0], avg_y_rel));
      result.mPoint2 = aRectangle->getPose().transformToGlobal(vect<double,2>(0.5 * aRectangle->getDimensions()[0] * ln_r_rel[0], avg_y_rel));
      result.mDistance = fabs(ln_c_rel[0]) - 0.5 * aRectangle->getDimensions()[0];
      return;
    };
    
    // there is no overlap, and thus, this boils down to a point-point problem.
    vect<double,2> re_pt_rel(0.0,0.0);
    vect<double,2> ln_pt_rel = ln_c_rel;
    if(ln_c_rel[0] < 0.0)
      re_pt_rel[0] -= 0.5 * aRectangle->getDimensions()[0];
    else
      re_pt_rel[0] += 0.5 * aRectangle->getDimensions()[0];
    if(ln_c_rel[1] < 0.0) {
      re_pt_rel[1] -= 0.5 * aRectangle->getDimensions()[1];
      ln_pt_rel[1] += half_length;
    } else {
      re_pt_rel[1] += 0.5 * aRectangle->getDimensions()[1];
      ln_pt_rel[1] -= half_length;
    };
    
    vect<double,2> diff_v_rel = ln_pt_rel - re_pt_rel;
    double dist_v_rel = norm_2(diff_v_rel);
    result.mPoint1 = aRectangle->getPose().transformToGlobal(ln_pt_rel);
    result.mPoint2 = aRectangle->getPose().transformToGlobal(re_pt_rel);
    result.mDistance = dist_v_rel;
    return;
  };
  
  if(fabs(ln_t_rel[1]) < 1e-5) {
    // this means the line is horizontal.
    if((ln_c_rel[0] + half_length > -0.5 * aRectangle->getDimensions()[0]) || 
       (ln_c_rel[0] - half_length <  0.5 * aRectangle->getDimensions()[0])) {
      // there is an overlap between the rectangle side and the line.
      double max_x_rel = ((ln_c_rel[0] + half_length <  0.5 * aRectangle->getDimensions()[0]) ? (ln_c_rel[0] + half_length) : ( 0.5 * aRectangle->getDimensions()[0]));
      double min_x_rel = ((ln_c_rel[0] - half_length > -0.5 * aRectangle->getDimensions()[0]) ? (ln_c_rel[0] - half_length) : (-0.5 * aRectangle->getDimensions()[0]));
      double avg_x_rel = (max_x_rel + min_x_rel) * 0.5;
      vect<double,2> ln_r_rel(0.0,1.0);
      if(ln_c_rel[1] < 0.0)
        ln_r_rel[1] = -1.0;
      result.mPoint1 = aRectangle->getPose().transformToGlobal(vect<double,2>(avg_x_rel, ln_c_rel[1]));
      result.mPoint2 = aRectangle->getPose().transformToGlobal(vect<double,2>(avg_x_rel, 0.5 * aRectangle->getDimensions()[1] * ln_r_rel[1]));
      result.mDistance = fabs(ln_c_rel[1]) - 0.5 * aRectangle->getDimensions()[1];
      return;
    };
    
    // there is no overlap, and thus, this boils down to a point-point problem.
    vect<double,2> re_pt_rel(0.0,0.0);
    vect<double,2> ln_pt_rel = ln_c_rel;
    if(ln_c_rel[1] < 0.0)
      re_pt_rel[1] -= 0.5 * aRectangle->getDimensions()[1];
    else
      re_pt_rel[1] += 0.5 * aRectangle->getDimensions()[1];
    if(ln_c_rel[0] < 0.0) {
      re_pt_rel[0] -= 0.5 * aRectangle->getDimensions()[0];
      ln_pt_rel[0] += half_length;
    } else {
      re_pt_rel[0] += 0.5 * aRectangle->getDimensions()[0];
      ln_pt_rel[0] -= half_length;
    };
    
    vect<double,2> diff_v_rel = ln_pt_rel - re_pt_rel;
    double dist_v_rel = norm_2(diff_v_rel);
    result.mPoint1 = aRectangle->getPose().transformToGlobal(ln_pt_rel);
    result.mPoint2 = aRectangle->getPose().transformToGlobal(re_pt_rel);
    result.mDistance = dist_v_rel;
    return;
  };
  
  // in any other case, we have to resort to the segment-point test.
  
  vect<double,2> ln_n_rel = 1.0 % ln_t_rel;
  if(ln_n_rel * ln_c_rel < 0.0)
    ln_n_rel = -ln_n_rel;
  vect<double,2> corner_pt(-0.5 * aRectangle->getDimensions()[0], -0.5 * aRectangle->getDimensions()[1]);
  if(ln_n_rel[0] > 0.0)
    corner_pt[0] = 0.5 * aRectangle->getDimensions()[0];
  if(ln_n_rel[1] > 0.0)
    corner_pt[1] = 0.5 * aRectangle->getDimensions()[1];
  
  vect<double,2> corner_pt_diff = (ln_c_rel - corner_pt);
  double dist_tmp = corner_pt_diff * ln_n_rel;
  double t_tmp = -(corner_pt_diff * ln_t_rel);
  if(fabs(t_tmp) > half_length) {
    if(t_tmp < 0.0)
      t_tmp = -half_length;
    else
      t_tmp = half_length;
    vect<double,2> ln_pt_rel = ln_c_rel + t_tmp * ln_t_rel;
    
    double in_x_range = fabs(ln_pt_rel[0]) - 0.5 * aRectangle->getDimensions()[0];
    double in_y_range = fabs(ln_pt_rel[1]) - 0.5 * aRectangle->getDimensions()[1];
    
    if((in_x_range < 0.0) && (in_y_range > in_x_range)) {
      corner_pt[0] = ln_pt_rel[0];
      dist_tmp = fabs(ln_pt_rel[1]) - 0.5 * aRectangle->getDimensions()[1];
    } else if((in_y_range < 0.0) && (in_x_range > in_y_range)) {
      corner_pt[1] = ln_pt_rel[1];
      dist_tmp = fabs(ln_pt_rel[0]) - 0.5 * aRectangle->getDimensions()[0];
    } else {
      if(ln_pt_rel[0] < 0.0)
        corner_pt[0] = -0.5 * aRectangle->getDimensions()[0];
      else
        corner_pt[0] =  0.5 * aRectangle->getDimensions()[0];
      if(ln_pt_rel[1] < 0.0)
        corner_pt[1] = -0.5 * aRectangle->getDimensions()[1];
      else
        corner_pt[1] =  0.5 * aRectangle->getDimensions()[1];
      dist_tmp = norm_2(ln_pt_rel - corner_pt);
    };
    
    result.mPoint1 = aRectangle->getPose().transformToGlobal(ln_pt_rel);
    result.mPoint2 = aRectangle->getPose().transformToGlobal(corner_pt);
    result.mDistance = dist_tmp;
  } else {
    result.mPoint1 = aRectangle->getPose().transformToGlobal(corner_pt + dist_tmp * ln_n_rel);
    result.mPoint2 = aRectangle->getPose().transformToGlobal(corner_pt);
    result.mDistance = dist_tmp;
  };
  
  return;
};


void prox_crect_rectangle::computeProximity() {
  mLastResult.mDistance = std::numeric_limits<double>::infinity();
  mLastResult.mPoint1 = vect<double,2>(0.0,0.0);
  mLastResult.mPoint2 = vect<double,2>(0.0,0.0);
  if((!mCRect) || (!mRectangle))
    return;
  
  
  vect<double,2> re_c = mRectangle->getPose().transformToGlobal(vect<double,2>(0.0,0.0));
  vect<double,2> cr_c = mCRect->getPose().transformToGlobal(vect<double,2>(0.0,0.0));
  vect<double,2> cr_t = mCRect->getPose().rotateToGlobal(vect<double,2>(1.0,0.0));
  
  computeProximityOfLine(mRectangle, cr_c, cr_t, 0.5 * mCRect->getDimensions()[0], mLastResult); 
  
  // add a circle-sweep around the line-rectangle solution.
  vect<double,2> diff_v = mLastResult.mPoint2 - mLastResult.mPoint1;
  double diff_d = norm_2(diff_v);
  if(mLastResult.mDistance < 0.0)
    mLastResult.mPoint1 -= (0.5 * mCRect->getDimensions()[1] / diff_d) * diff_v;
  else
    mLastResult.mPoint1 += (0.5 * mCRect->getDimensions()[1] / diff_d) * diff_v;
  mLastResult.mDistance -= 0.5 * mCRect->getDimensions()[1];
  return;
};


prox_crect_rectangle::prox_crect_rectangle(const shared_ptr< capped_rectangle >& aCRect,
                                           const shared_ptr< rectangle >& aRectangle) :
                                           proximity_finder_2D(),
                                           mCRect(aCRect),
                                           mRectangle(aRectangle) { };
    
    
void RK_CALL prox_crect_rectangle::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_2D::save(A,proximity_finder_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mCRect)
    & RK_SERIAL_SAVE_WITH_NAME(mRectangle);
};
    
void RK_CALL prox_crect_rectangle::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_2D::load(A,proximity_finder_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mCRect)
    & RK_SERIAL_LOAD_WITH_NAME(mRectangle);
};


};

};











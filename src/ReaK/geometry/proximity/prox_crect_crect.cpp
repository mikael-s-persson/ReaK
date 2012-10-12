
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

#include "prox_crect_crect.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_2D > prox_crect_crect::getShape1() const {
  return mCRect1;
};

shared_ptr< shape_2D > prox_crect_crect::getShape2() const {
  return mCRect2;
};
    
void prox_crect_crect::computeProximity() {
  mLastResult.mDistance = std::numeric_limits<double>::infinity();
  mLastResult.mPoint1 = vect<double,2>(0.0,0.0);
  mLastResult.mPoint2 = vect<double,2>(0.0,0.0);
  if((!mCRect1) || (!mCRect2))
    return;
  
  
  vect<double,2> cr1_c = mCRect1->getPose().transformToGlobal(vect<double,2>(0.0,0.0));
  vect<double,2> cr2_c = mCRect2->getPose().transformToGlobal(vect<double,2>(0.0,0.0));
  vect<double,2> cr2_t = mCRect2->getPose().rotateToGlobal(vect<double,2>(1.0,0.0));
  
  vect<double,2> cr2_c_rel = mCRect1->getPose().transformFromGlobal(cr2_c);
  vect<double,2> cr2_t_rel = mCRect1->getPose().rotateFromGlobal(cr2_t);
  
  
  if(fabs(cr2_t_rel[1]) < 1e-5) {
    // The capped-rectangles are parallel.
    if((cr2_c_rel[0] + 0.5 * mCRect2->getDimensions()[0] > -0.5 * mCRect1->getDimensions()[0]) || 
       (cr2_c_rel[0] - 0.5 * mCRect2->getDimensions()[0] <  0.5 * mCRect1->getDimensions()[0])) {
      // there is an overlap between the capped-rectangle sides.
      double max_x_rel = ((cr2_c_rel[0] + 0.5 * mCRect2->getDimensions()[0] <  0.5 * mCRect1->getDimensions()[0]) ? (cr2_c_rel[0] + 0.5 * mCRect2->getDimensions()[0]) : ( 0.5 * mCRect1->getDimensions()[0]));
      double min_x_rel = ((cr2_c_rel[0] - 0.5 * mCRect2->getDimensions()[0] > -0.5 * mCRect1->getDimensions()[0]) ? (cr2_c_rel[0] - 0.5 * mCRect2->getDimensions()[0]) : (-0.5 * mCRect1->getDimensions()[0]));
      double avg_x_rel = (max_x_rel + min_x_rel) * 0.5;
      vect<double,2> cr2_r_rel(0.0,1.0);
      if(cr2_c_rel[1] < 0.0)
        cr2_r_rel[1] = -1.0;
      mLastResult.mPoint1 = mCRect1->getPose().transformToGlobal(vect<double,2>(avg_x_rel, 0.5 * mCRect1->getDimensions()[1] * cr2_r_rel[1]));
      mLastResult.mPoint2 = mCRect1->getPose().transformToGlobal(vect<double,2>(avg_x_rel, cr2_c_rel[1] - 0.5 * mCRect2->getDimensions()[1] * cr2_r_rel[1]));
      mLastResult.mDistance = fabs(cr2_c_rel[1]) - 0.5 * mCRect1->getDimensions()[1] - 0.5 * mCRect2->getDimensions()[1];
      return;
    };
    // there is no overlap, and thus, this boils down to a circle-circle problem.
    vect<double,2> cr1_cic_rel(0.0,0.0);
    vect<double,2> cr2_cic_rel = cr2_c_rel;
    if(cr2_c_rel[0] < 0.0) {
      cr1_cic_rel[0] -= 0.5 * mCRect1->getDimensions()[0];
      cr2_cic_rel[0] += 0.5 * mCRect2->getDimensions()[0];
    } else {
      cr1_cic_rel[0] += 0.5 * mCRect1->getDimensions()[0];
      cr2_cic_rel[0] -= 0.5 * mCRect2->getDimensions()[0];
    };
    vect<double,2> diff_v_rel = cr2_cic_rel - cr1_cic_rel;
    double dist_v_rel = norm_2(diff_v_rel);
    mLastResult.mPoint1 = mCRect1->getPose().transformToGlobal(cr1_cic_rel + (0.5 * mCRect1->getDimensions()[1] / dist_v_rel) * diff_v_rel);
    mLastResult.mPoint2 = mCRect1->getPose().transformToGlobal(cr2_cic_rel - (0.5 * mCRect2->getDimensions()[1] / dist_v_rel) * diff_v_rel);
    mLastResult.mDistance = dist_v_rel - 0.5 * mCRect1->getDimensions()[1] - 0.5 * mCRect2->getDimensions()[1];
    return;
  };
  
  
  // Line-Line solution:
  double d = cr2_t_rel * cr2_c_rel;
  double denom = 1.0 - cr2_t_rel[0] * cr2_t_rel[0];
  double s_c = (cr2_t_rel[0] * cr2_c_rel[0] - d) / denom;
  double t_c = (cr2_c_rel[0] - cr2_t_rel[0] * d) / denom;
  
  // Segment-Segment solution:
  if(s_c < -0.5 * mCRect2->getDimensions()[0]) {
    s_c = -0.5 * mCRect2->getDimensions()[0];
    t_c = cr2_c_rel[0] - 0.5 * mCRect2->getDimensions()[0] * cr2_t_rel[0];
  } else if(s_c > 0.5 * mCRect2->getDimensions()[0]) {
    s_c = 0.5 * mCRect2->getDimensions()[0];
    t_c = cr2_c_rel[0] + 0.5 * mCRect2->getDimensions()[0] * cr2_t_rel[0];
  };
  
  if(t_c < -0.5 * mCRect1->getDimensions()[0]) {
    t_c = -0.5 * mCRect1->getDimensions()[0];
    s_c = -0.5 * mCRect1->getDimensions()[0] * cr2_t_rel[0] - d;
  } else if(t_c > 0.5 * mCRect1->getDimensions()[0]) {
    t_c = 0.5 * mCRect1->getDimensions()[0];
    s_c = 0.5 * mCRect1->getDimensions()[0] * cr2_t_rel[0] - d;
  };
  
  if(s_c < -0.5 * mCRect2->getDimensions()[0])
    s_c = -0.5 * mCRect2->getDimensions()[0];
  else if(s_c > 0.5 * mCRect2->getDimensions()[0])
    s_c = 0.5 * mCRect2->getDimensions()[0];
  
  // we have parameters s and t for the min-dist points on the center segments.
  
  // just apply a circle-sweep on the line-segments.
  vect<double,2> cr1_ptc(t_c,0.0);
  vect<double,2> cr2_ptc = cr2_c_rel + s_c * cr2_t_rel;
  
  vect<double,2> diff_v_rel = cr2_ptc - cr1_ptc;
  double dist_v_rel = norm_2(diff_v_rel);
  mLastResult.mPoint1 = mCRect1->getPose().transformToGlobal(cr1_ptc + (0.5 * mCRect1->getDimensions()[1] / dist_v_rel) * diff_v_rel);
  mLastResult.mPoint2 = mCRect1->getPose().transformToGlobal(cr2_ptc - (0.5 * mCRect2->getDimensions()[1] / dist_v_rel) * diff_v_rel);
  mLastResult.mDistance = dist_v_rel - 0.5 * mCRect1->getDimensions()[1] - 0.5 * mCRect2->getDimensions()[1];
  
};


prox_crect_crect::prox_crect_crect(const shared_ptr< capped_rectangle >& aCRect1,
                                   const shared_ptr< capped_rectangle >& aCRect2) :
                                   proximity_finder_2D(),
                                   mCRect1(aCRect1),
                                   mCRect2(aCRect2) { };
    
    
void RK_CALL prox_crect_crect::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_2D::save(A,proximity_finder_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mCRect1)
    & RK_SERIAL_SAVE_WITH_NAME(mCRect2);
};
    
void RK_CALL prox_crect_crect::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_2D::load(A,proximity_finder_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mCRect1)
    & RK_SERIAL_LOAD_WITH_NAME(mCRect2);
};


};

};











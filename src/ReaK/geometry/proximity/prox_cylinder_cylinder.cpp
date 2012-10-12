
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

#include "prox_cylinder_cylinder.hpp"

#include <cmath>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_3D > prox_cylinder_cylinder::getShape1() const {
  return mCylinder1;
};

shared_ptr< shape_3D > prox_cylinder_cylinder::getShape2() const {
  return mCylinder2;
};

void prox_cylinder_cylinder::computeProximity() {
  if((!mCylinder1) || (!mCylinder2)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  using std::fabs; using std::sqrt; using ReaK::unit; using ReaK::norm_2;
  
  vect<double,3> cy1_c = mCylinder1->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  vect<double,3> cy2_c = mCylinder2->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  vect<double,3> cy2_t = mCylinder2->getPose().rotateToGlobal(vect<double,3>(0.0,0.0,1.0));
  
  vect<double,3> cy2_c_rel = mCylinder1->getPose().transformFromGlobal(cy2_c);
  vect<double,3> cy2_t_rel = mCylinder1->getPose().rotateFromGlobal(cy2_t);
  
  if(sqrt(cy2_t_rel[0] * cy2_t_rel[0] + cy2_t_rel[1] * cy2_t_rel[1]) < 1e-5) {
    // The capped-cylinders are parallel.
    if((cy2_c_rel[2] + 0.5 * mCylinder2->getLength() > -0.5 * mCylinder1->getLength()) || 
       (cy2_c_rel[2] - 0.5 * mCylinder2->getLength() <  0.5 * mCylinder1->getLength())) {
      // there is an overlap between the capped-cylinder sides.
      double max_z_rel = ((cy2_c_rel[2] + 0.5 * mCylinder2->getLength() <  0.5 * mCylinder1->getLength()) ? (cy2_c_rel[2] + 0.5 * mCylinder2->getLength()) : ( 0.5 * mCylinder1->getLength()));
      double min_z_rel = ((cy2_c_rel[2] - 0.5 * mCylinder2->getLength() > -0.5 * mCylinder1->getLength()) ? (cy2_c_rel[2] - 0.5 * mCylinder2->getLength()) : (-0.5 * mCylinder1->getLength()));
      double avg_z_rel = (max_z_rel + min_z_rel) * 0.5;
      vect<double,3> cy2_r_rel(cy2_c_rel[0],cy2_c_rel[1],0.0);
      double cy2_r_mag = norm_2(cy2_r_rel);
      mLastResult.mPoint1 = mCylinder1->getPose().transformToGlobal(vect<double,3>(mCylinder1->getRadius() * cy2_r_rel[0] / cy2_r_mag, mCylinder1->getRadius() * cy2_r_rel[1] / cy2_r_mag, avg_z_rel));
      mLastResult.mPoint2 = mCylinder1->getPose().transformToGlobal(vect<double,3>(cy2_c_rel[0] - mCylinder2->getRadius() * cy2_r_rel[0] / cy2_r_mag, cy2_c_rel[1] - mCylinder2->getRadius() * cy2_r_rel[1] / cy2_r_mag, avg_z_rel));
      mLastResult.mDistance = cy2_r_mag - mCylinder1->getRadius() - mCylinder2->getRadius();
      if((mLastResult.mDistance < 0.0) && (mLastResult.mDistance > min_z_rel - max_z_rel)) {
        // this means that the collision is mostly on the top/bottom sides
        mLastResult.mDistance = min_z_rel - max_z_rel;
        if(cy2_c_rel[2] < 0.0) {
          mLastResult.mPoint1 = mCylinder1->getPose().transformToGlobal(vect<double,3>(0.5 * cy2_r_rel[0], 0.5 * cy2_r_rel[1], -0.5 * mCylinder1->getLength()));
          mLastResult.mPoint2 = mCylinder1->getPose().transformToGlobal(vect<double,3>(0.5 * cy2_r_rel[0], 0.5 * cy2_r_rel[1], cy2_c_rel[2] + 0.5 * mCylinder2->getLength()));
        } else {
          mLastResult.mPoint1 = mCylinder1->getPose().transformToGlobal(vect<double,3>(0.5 * cy2_r_rel[0], 0.5 * cy2_r_rel[1],  0.5 * mCylinder1->getLength()));
          mLastResult.mPoint2 = mCylinder1->getPose().transformToGlobal(vect<double,3>(0.5 * cy2_r_rel[0], 0.5 * cy2_r_rel[1], cy2_c_rel[2] - 0.5 * mCylinder2->getLength()));
        };
      };
      return;
    };
    // there is no overlap, and thus, this boils down to a disk-disk problem.
    vect<double,3> cy1_spc_rel(0.0,0.0,0.0);
    vect<double,3> cy2_spc_rel = cy2_c_rel;
    if(cy2_c_rel[2] < 0.0) {
      cy1_spc_rel[2] -= 0.5 * mCylinder1->getLength();
      cy2_spc_rel[2] += 0.5 * mCylinder2->getLength();
    } else {
      cy1_spc_rel[2] += 0.5 * mCylinder1->getLength();
      cy2_spc_rel[2] -= 0.5 * mCylinder2->getLength();
    };
    vect<double,3> diff_v_rel = cy2_spc_rel - cy1_spc_rel;
    mLastResult.mDistance = fabs(diff_v_rel[2]);
    diff_v_rel[2] = 0.0;
    double dist_v_rel = norm_2(diff_v_rel);
    if(dist_v_rel > mCylinder1->getRadius() + mCylinder2->getRadius()) {
      mLastResult.mPoint1 = mCylinder1->getPose().transformToGlobal(cy1_spc_rel + (mCylinder1->getRadius() / dist_v_rel) * diff_v_rel);
      mLastResult.mPoint2 = mCylinder1->getPose().transformToGlobal(cy2_spc_rel - (mCylinder2->getRadius() / dist_v_rel) * diff_v_rel);
    } else {
      double d_offset = 0.5 * (dist_v_rel - mCylinder1->getRadius() - mCylinder2->getRadius());
      mLastResult.mPoint1 = mCylinder1->getPose().transformToGlobal(cy1_spc_rel + ((mCylinder1->getRadius() + d_offset) / dist_v_rel) * diff_v_rel);
      mLastResult.mPoint2 = mCylinder1->getPose().transformToGlobal(cy2_spc_rel - ((mCylinder2->getRadius() + d_offset) / dist_v_rel) * diff_v_rel);
    };
    return;
  };
  
  // Line-Line solution:
  double d = cy2_t_rel * cy2_c_rel;
  double denom = 1.0 - cy2_t_rel[2] * cy2_t_rel[2];
  double s_c = (cy2_t_rel[2] * cy2_c_rel[2] - d) / denom;
  double t_c = (cy2_c_rel[2] - cy2_t_rel[2] * d) / denom;
  double s_m = 0.0;
  double t_m = 0.0;
  
  // Segment-Segment solution:
  if(s_c < -0.5 * mCylinder2->getLength()) {
    s_c = -0.5 * mCylinder2->getLength();
    s_m = -1.0; // lower-bound.
    t_c = cy2_c_rel[2] - 0.5 * mCylinder2->getLength() * cy2_t_rel[2];
  } else if(s_c > 0.5 * mCylinder2->getLength()) {
    s_c = 0.5 * mCylinder2->getLength();
    s_m = 1.0; // upper-bound.
    t_c = cy2_c_rel[2] + 0.5 * mCylinder2->getLength() * cy2_t_rel[2];
  };
  
  if(t_c < -0.5 * mCylinder1->getLength()) {
    t_c = -0.5 * mCylinder1->getLength();
    t_m = -1.0; // lower-bound.
    s_m = 0.0; // reset.
    s_c = -0.5 * mCylinder1->getLength() * cy2_t_rel[2] - d;
  } else if(t_c > 0.5 * mCylinder1->getLength()) {
    t_c = 0.5 * mCylinder1->getLength();
    t_m = 1.0; // upper-bound.
    s_m = 0.0; // reset.
    s_c = 0.5 * mCylinder1->getLength() * cy2_t_rel[2] - d;
  };
  
  if(s_c < -0.5 * mCylinder2->getLength()) {
    s_c = -0.5 * mCylinder2->getLength();
    s_m = -1.0; // lower-bound.
  } else if(s_c > 0.5 * mCylinder2->getLength()) {
    s_c = 0.5 * mCylinder2->getLength();
    s_m = 1.0; // upper-bound.
  };
  
  // we have parameters s and t for the min-dist points on the center segments.
  
  // just apply a sphere-sweep on the line-segments.
  vect<double,3> cy1_ptc(0.0,0.0,t_c);
  vect<double,3> cy2_ptc = cy2_c_rel + s_c * cy2_t_rel;
  vect<double,3> diff_v_rel = cy2_ptc - cy1_ptc;
  
  if((fabs(s_m) < 0.5) && (fabs(t_m) < 0.5)) {
    // this means we have a side-to-side proximity.
    double dist_v_rel = norm_2(diff_v_rel);
    mLastResult.mPoint1 = mCylinder1->getPose().transformToGlobal(cy1_ptc + (mCylinder1->getRadius() / dist_v_rel) * diff_v_rel);
    mLastResult.mPoint2 = mCylinder1->getPose().transformToGlobal(cy2_ptc - (mCylinder2->getRadius() / dist_v_rel) * diff_v_rel);
    mLastResult.mDistance = dist_v_rel - mCylinder1->getRadius() - mCylinder2->getRadius();
    return;
  };
  
  // at this point, we have to deal with more complicated cases.
  
  if((fabs(s_m) > 0.5) && (fabs(t_m) > 0.5)) {
    // nominally, we can expect the cylinders to be end near end.
    
  };
  
};


prox_cylinder_cylinder::prox_cylinder_cylinder(const shared_ptr< cylinder >& aCylinder1,
                                               const shared_ptr< cylinder >& aCylinder2) :
                                               proximity_finder_3D(),
                                               mCylinder1(aCylinder1),
                                               mCylinder2(aCylinder2) { };


void RK_CALL prox_cylinder_cylinder::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_3D::save(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mCylinder1)
    & RK_SERIAL_SAVE_WITH_NAME(mCylinder2);
};

void RK_CALL prox_cylinder_cylinder::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_3D::load(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mCylinder1)
    & RK_SERIAL_LOAD_WITH_NAME(mCylinder2);
};


};

};












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

#include <ReaK/geometry/proximity/prox_ccylinder_ccylinder.hpp>

#include <cmath>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_3D > prox_ccylinder_ccylinder::getShape1() const {
  return mCCylinder1;
};

shared_ptr< shape_3D > prox_ccylinder_ccylinder::getShape2() const {
  return mCCylinder2;
};

void prox_ccylinder_ccylinder::computeProximity() {
  if((!mCCylinder1) || (!mCCylinder2)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  using std::fabs; using std::sqrt; using ReaK::unit; using ReaK::norm_2;
  
  vect<double,3> cy1_c = mCCylinder1->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  vect<double,3> cy2_c = mCCylinder2->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  vect<double,3> cy2_t = mCCylinder2->getPose().rotateToGlobal(vect<double,3>(0.0,0.0,1.0));
  
  vect<double,3> cy2_c_rel = mCCylinder1->getPose().transformFromGlobal(cy2_c);
  vect<double,3> cy2_t_rel = mCCylinder1->getPose().rotateFromGlobal(cy2_t);
  
  if(sqrt(cy2_t_rel[0] * cy2_t_rel[0] + cy2_t_rel[1] * cy2_t_rel[1]) < 1e-5) {
    // The capped-cylinders are parallel.
    if((cy2_c_rel[2] + 0.5 * mCCylinder2->getLength() > -0.5 * mCCylinder1->getLength()) || 
       (cy2_c_rel[2] - 0.5 * mCCylinder2->getLength() <  0.5 * mCCylinder1->getLength())) {
      // there is an overlap between the capped-cylinder sides.
      double max_z_rel = ((cy2_c_rel[2] + 0.5 * mCCylinder2->getLength() <  0.5 * mCCylinder1->getLength()) ? (cy2_c_rel[2] + 0.5 * mCCylinder2->getLength()) : ( 0.5 * mCCylinder1->getLength()));
      double min_z_rel = ((cy2_c_rel[2] - 0.5 * mCCylinder2->getLength() > -0.5 * mCCylinder1->getLength()) ? (cy2_c_rel[2] - 0.5 * mCCylinder2->getLength()) : (-0.5 * mCCylinder1->getLength()));
      double avg_z_rel = (max_z_rel + min_z_rel) * 0.5;
      vect<double,3> cy2_r_rel = unit(vect<double,3>(cy2_c_rel[0],cy2_c_rel[1],0.0));
      mLastResult.mPoint1 = mCCylinder1->getPose().transformToGlobal(vect<double,3>(mCCylinder1->getRadius() * cy2_r_rel[0], mCCylinder1->getRadius() * cy2_r_rel[1], avg_z_rel));
      mLastResult.mPoint2 = mCCylinder1->getPose().transformToGlobal(vect<double,3>(cy2_c_rel[0] - mCCylinder2->getRadius() * cy2_r_rel[0], cy2_c_rel[1] - mCCylinder2->getRadius() * cy2_r_rel[1], avg_z_rel));
      mLastResult.mDistance = sqrt(cy2_c_rel[0] * cy2_c_rel[0] + cy2_c_rel[1] * cy2_c_rel[1]) - mCCylinder1->getRadius() - mCCylinder2->getRadius();
      return;
    };
    // there is no overlap, and thus, this boils down to a sphere-sphere problem.
    vect<double,3> cy1_spc_rel(0.0,0.0,0.0);
    vect<double,3> cy2_spc_rel = cy2_c_rel;
    if(cy2_c_rel[2] < 0.0) {
      cy1_spc_rel[2] -= 0.5 * mCCylinder1->getLength();
      cy2_spc_rel[2] += 0.5 * mCCylinder2->getLength();
    } else {
      cy1_spc_rel[2] += 0.5 * mCCylinder1->getLength();
      cy2_spc_rel[2] -= 0.5 * mCCylinder2->getLength();
    };
    vect<double,3> diff_v_rel = cy2_spc_rel - cy1_spc_rel;
    double dist_v_rel = norm_2(diff_v_rel);
    mLastResult.mPoint1 = mCCylinder1->getPose().transformToGlobal(cy1_spc_rel + (mCCylinder1->getRadius() / dist_v_rel) * diff_v_rel);
    mLastResult.mPoint2 = mCCylinder1->getPose().transformToGlobal(cy2_spc_rel - (mCCylinder2->getRadius() / dist_v_rel) * diff_v_rel);
    mLastResult.mDistance = dist_v_rel - mCCylinder1->getRadius() - mCCylinder2->getRadius();
    return;
  };
  
  // Line-Line solution:
  double d = cy2_t_rel * cy2_c_rel;
  double denom = 1.0 - cy2_t_rel[2] * cy2_t_rel[2];
  double s_c = (cy2_t_rel[2] * cy2_c_rel[2] - d) / denom;
  double t_c = (cy2_c_rel[2] - cy2_t_rel[2] * d) / denom;
  
  // Segment-Segment solution:
  if(s_c < -0.5 * mCCylinder2->getLength()) {
    s_c = -0.5 * mCCylinder2->getLength();
    t_c = cy2_c_rel[2] - 0.5 * mCCylinder2->getLength() * cy2_t_rel[2];
  } else if(s_c > 0.5 * mCCylinder2->getLength()) {
    s_c = 0.5 * mCCylinder2->getLength();
    t_c = cy2_c_rel[2] + 0.5 * mCCylinder2->getLength() * cy2_t_rel[2];
  };
  
  if(t_c < -0.5 * mCCylinder1->getLength()) {
    t_c = -0.5 * mCCylinder1->getLength();
    s_c = -0.5 * mCCylinder1->getLength() * cy2_t_rel[2] - d;
  } else if(t_c > 0.5 * mCCylinder1->getLength()) {
    t_c = 0.5 * mCCylinder1->getLength();
    s_c = 0.5 * mCCylinder1->getLength() * cy2_t_rel[2] - d;
  };
  
  if(s_c < -0.5 * mCCylinder2->getLength())
    s_c = -0.5 * mCCylinder2->getLength();
  else if(s_c > 0.5 * mCCylinder2->getLength())
    s_c = 0.5 * mCCylinder2->getLength();
  
  // we have parameters s and t for the min-dist points on the center segments.
  
  // just apply a sphere-sweep on the line-segments.
  vect<double,3> cy1_ptc(0.0,0.0,t_c);
  vect<double,3> cy2_ptc = cy2_c_rel + s_c * cy2_t_rel;
  
  vect<double,3> diff_v_rel = cy2_ptc - cy1_ptc;
  double dist_v_rel = norm_2(diff_v_rel);
  mLastResult.mPoint1 = mCCylinder1->getPose().transformToGlobal(cy1_ptc + (mCCylinder1->getRadius() / dist_v_rel) * diff_v_rel);
  mLastResult.mPoint2 = mCCylinder1->getPose().transformToGlobal(cy2_ptc - (mCCylinder2->getRadius() / dist_v_rel) * diff_v_rel);
  mLastResult.mDistance = dist_v_rel - mCCylinder1->getRadius() - mCCylinder2->getRadius();
};


prox_ccylinder_ccylinder::prox_ccylinder_ccylinder(const shared_ptr< capped_cylinder >& aCCylinder1,
                                                   const shared_ptr< capped_cylinder >& aCCylinder2) :
                                                   proximity_finder_3D(),
                                                   mCCylinder1(aCCylinder1),
                                                   mCCylinder2(aCCylinder2) { };


void RK_CALL prox_ccylinder_ccylinder::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_3D::save(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mCCylinder1)
    & RK_SERIAL_SAVE_WITH_NAME(mCCylinder2);
};

void RK_CALL prox_ccylinder_ccylinder::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_3D::load(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mCCylinder1)
    & RK_SERIAL_LOAD_WITH_NAME(mCCylinder2);
};


};

};











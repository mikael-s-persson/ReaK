
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

#include "prox_ccylinder_cylinder.hpp"

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
  using std::fabs; using std::sqrt; using ReaK::unit; using ReaK::norm_2;
  
  vect<double,3> cy1_c = mCylinder->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  vect<double,3> cy2_c = mCCylinder->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  vect<double,3> cy2_t = mCCylinder->getPose().rotateToGlobal(vect<double,3>(0.0,0.0,1.0));
  
  vect<double,3> cy2_c_rel = mCylinder->getPose().transformFromGlobal(cy2_c);
  vect<double,3> cy2_t_rel = mCylinder->getPose().rotateFromGlobal(cy2_t);
  
  if(sqrt(cy2_t_rel[0] * cy2_t_rel[0] + cy2_t_rel[1] * cy2_t_rel[1]) < 1e-5) {
    // The capped-cylinders are parallel.
    if((cy2_c_rel[2] + 0.5 * mCCylinder->getLength() > -0.5 * mCylinder->getLength()) || 
       (cy2_c_rel[2] - 0.5 * mCCylinder->getLength() <  0.5 * mCylinder->getLength())) {
      // there is an overlap between the capped-cylinder sides.
      double max_z_rel = ((cy2_c_rel[2] + 0.5 * mCCylinder->getLength() <  0.5 * mCylinder->getLength()) ? (cy2_c_rel[2] + 0.5 * mCCylinder->getLength()) : ( 0.5 * mCCylinder1->getLength()));
      double min_z_rel = ((cy2_c_rel[2] - 0.5 * mCCylinder->getLength() > -0.5 * mCylinder->getLength()) ? (cy2_c_rel[2] - 0.5 * mCCylinder->getLength()) : (-0.5 * mCCylinder1->getLength()));
      double avg_z_rel = (max_z_rel + min_z_rel) * 0.5;
      vect<double,3> cy2_r_rel = unit(vect<double,3>(cy2_c_rel[0],cy2_c_rel[1],0.0));
      mLastResult.mPoint1 = mCylinder->getPose().transformToGlobal(vect<double,3>(mCylinder->getRadius() * cy2_r_rel[0], mCylinder->getRadius() * cy2_r_rel[1], avg_z_rel));
      mLastResult.mPoint2 = mCylinder->getPose().transformToGlobal(vect<double,3>(cy2_c_rel[0] - mCCylinder->getRadius() * cy2_r_rel[0], cy2_c_rel[1] - mCCylinder->getRadius() * cy2_r_rel[1], avg_z_rel));
      mLastResult.mDistance = sqrt(cy2_c_rel[0] * cy2_c_rel[0] + cy2_c_rel[1] * cy2_c_rel[1]) - mCylinder->getRadius() - mCCylinder->getRadius();
      return;
    };
    // there is no overlap, and thus, this boils down to a sphere-disk problem.
    vect<double,3> cy1_spc_rel(0.0,0.0,0.0);
    vect<double,3> cy2_spc_rel = cy2_c_rel;
    if(cy2_c_rel[2] < 0.0) {
      cy1_spc_rel[2] -= 0.5 * mCylinder->getLength();
      cy2_spc_rel[2] += 0.5 * mCCylinder->getLength();
    } else {
      cy1_spc_rel[2] += 0.5 * mCylinder->getLength();
      cy2_spc_rel[2] -= 0.5 * mCCylinder->getLength();
    };
    double dist_v_rel = sqrt(cy2_spc_rel[0] * cy2_spc_rel[0] + cy2_spc_rel[1] * cy2_spc_rel[1]);
    if(dist_v_rel < mCylinder->getRadius()) {
      cy1_spc_rel[0] = cy2_spc_rel[0];
      cy1_spc_rel[1] = cy2_spc_rel[1];
    } else {
      cy1_spc_rel[0] = (mCylinder->getRadius() / dist_v_rel) * cy2_spc_rel[0];
      cy1_spc_rel[1] = (mCylinder->getRadius() / dist_v_rel) * cy2_spc_rel[1];
    };
    vect<double,3> diff_v_rel = cy2_spc_rel - cy1_spc_rel;
    dist_v_rel = norm_2(diff_v_rel);
    mLastResult.mPoint1 = mCylinder->getPose().transformToGlobal(cy2_spc_rel - (mCCylinder->getRadius() / dist_v_rel) * diff_v_rel);
    mLastResult.mPoint2 = mCylinder->getPose().transformToGlobal(cy1_spc_rel);
    mLastResult.mDistance = dist_v_rel - mCCylinder->getRadius();
    return;
  };
  
  
  // NOTE: must resort to a non-linear solver.
  
  
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











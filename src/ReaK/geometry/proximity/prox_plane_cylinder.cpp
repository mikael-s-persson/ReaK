
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

#include <ReaK/geometry/proximity/prox_plane_cylinder.hpp>

#include <cmath>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_3D > prox_plane_cylinder::getShape1() const {
  return mPlane;
};

shared_ptr< shape_3D > prox_plane_cylinder::getShape2() const {
  return mCylinder;
};

void prox_plane_cylinder::computeProximity() {
  if((!mCylinder) || (!mPlane)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  using std::fabs; using std::sqrt; using ReaK::unit;
  
  vect<double,3> cy_c = mCylinder->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  vect<double,3> cy_t = mCylinder->getPose().rotateToGlobal(vect<double,3>(0.0,0.0,1.0));
  vect<double,3> pl_c = mPlane->getPose().transformToGlobal(vect<double,3>(0.0,0.0,0.0));
  
  vect<double,3> cy_c_rel = mPlane->getPose().transformFromGlobal(cy_c);
  vect<double,3> cy_t_rel = mPlane->getPose().rotateFromGlobal(cy_t);
  
  if(fabs(cy_t_rel[2]) < 1e-6) {
    // The cylinder is sitting flat (on round side) on the plane.
    mLastResult.mPoint1 = mPlane->getPose().transformToGlobal(vect<double,3>(cy_c_rel[0],cy_c_rel[1],0.0));
    mLastResult.mPoint2 = mPlane->getPose().transformToGlobal(vect<double,3>(cy_c_rel[0],cy_c_rel[1],cy_c_rel[2] - mCylinder->getRadius()));
    mLastResult.mDistance = cy_c_rel[2] - mCylinder->getRadius();
  } else if(sqrt(cy_t_rel[0] * cy_t_rel[0] + cy_t_rel[1] * cy_t_rel[1]) < 1e-6) {
    // The cylinder is sitting flat (on flat ends) on the plane.
    mLastResult.mPoint1 = mPlane->getPose().transformToGlobal(vect<double,3>(cy_c_rel[0],cy_c_rel[1],0.0));
    mLastResult.mPoint2 = mPlane->getPose().transformToGlobal(vect<double,3>(cy_c_rel[0],cy_c_rel[1],cy_c_rel[2] - 0.5 * mCylinder->getLength()));
    mLastResult.mDistance = cy_c_rel[2] - 0.5 * mCylinder->getLength();
  } else {
    // The cylinder is at an angle to the plane.
    if(cy_t_rel[2] > 0.0)
      cy_t_rel = -cy_t_rel;
    vect<double,3> cy_r_rel = unit(vect<double,3>(0.0,0.0,-1.0) + cy_t_rel[2] * cy_t_rel);
    vect<double,3> cypt_rel = cy_c_rel + (0.5 * mCylinder->getLength()) * cy_t_rel + mCylinder->getRadius() * cy_r_rel;
    mLastResult.mPoint1 = mPlane->getPose().transformToGlobal(vect<double,3>(cypt_rel[0],cypt_rel[1],0.0));
    mLastResult.mPoint2 = mPlane->getPose().transformToGlobal(cypt_rel);
    mLastResult.mDistance = cypt_rel[2];
  };
};


prox_plane_cylinder::prox_plane_cylinder(const shared_ptr< plane >& aPlane,
                                         const shared_ptr< cylinder >& aCylinder) :
                                         proximity_finder_3D(),
                                         mPlane(aPlane),
                                         mCylinder(aCylinder) { };


void RK_CALL prox_plane_cylinder::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_3D::save(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mPlane)
    & RK_SERIAL_SAVE_WITH_NAME(mCylinder);
};

void RK_CALL prox_plane_cylinder::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_3D::load(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mPlane)
    & RK_SERIAL_LOAD_WITH_NAME(mCylinder);
};


};

};











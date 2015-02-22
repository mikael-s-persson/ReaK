
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

#include <ReaK/geometry/proximity/prox_plane_box.hpp>

#include <cmath>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_3D > prox_plane_box::getShape1() const {
  return mPlane;
};

shared_ptr< shape_3D > prox_plane_box::getShape2() const {
  return mBox;
};

void prox_plane_box::computeProximity(const shape_3D_precompute_pack& aPack1, 
                                      const shape_3D_precompute_pack& aPack2) {
  if((!mBox) || (!mPlane)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  using std::fabs; using std::sqrt; using ReaK::unit;
  
  const pose_3D<double>& bx_pose = (aPack1.parent == mBox.get() ? 
                                    aPack1.global_pose : aPack2.global_pose);
  const pose_3D<double>& pl_pose = (aPack1.parent == mBox.get() ? 
                                    aPack2.global_pose : aPack1.global_pose);
  
  vect<double,3> bx_c = bx_pose.Position;
  vect<double,3> bx_x = pl_pose.rotateFromGlobal(bx_pose.rotateToGlobal(vect<double,3>(1.0,0.0,0.0)));
  vect<double,3> bx_y = pl_pose.rotateFromGlobal(bx_pose.rotateToGlobal(vect<double,3>(1.0,0.0,0.0)));
  vect<double,3> bx_z = pl_pose.rotateFromGlobal(bx_pose.rotateToGlobal(vect<double,3>(1.0,0.0,0.0)));
  
  if(bx_x[2] > 0.0)
    bx_x = -bx_x;
  if(bx_y[2] > 0.0)
    bx_y = -bx_y;
  if(bx_z[2] > 0.0)
    bx_z = -bx_z;
  
  vect<double,3> bx_c_rel = pl_pose.transformFromGlobal(bx_c);
  vect<double,3> bx_pt_rel = bx_c_rel + 0.5 * (mBox->getDimensions()[0] * bx_x + mBox->getDimensions()[1] * bx_y + mBox->getDimensions()[2] * bx_z);
  
  mLastResult.mPoint1 = pl_pose.transformToGlobal(vect<double,3>(bx_pt_rel[0],bx_pt_rel[1],0.0));
  mLastResult.mPoint2 = pl_pose.transformToGlobal(bx_pt_rel);
  mLastResult.mDistance = bx_pt_rel[2];
};


prox_plane_box::prox_plane_box(const shared_ptr< plane >& aPlane,
                               const shared_ptr< box >& aBox) :
                               proximity_finder_3D(),
                               mPlane(aPlane),
                               mBox(aBox) { };


void RK_CALL prox_plane_box::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_3D::save(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mPlane)
    & RK_SERIAL_SAVE_WITH_NAME(mBox);
};

void RK_CALL prox_plane_box::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_3D::load(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mPlane)
    & RK_SERIAL_LOAD_WITH_NAME(mBox);
};


};

};











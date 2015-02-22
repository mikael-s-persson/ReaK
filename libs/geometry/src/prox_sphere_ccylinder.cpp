
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

#include <ReaK/geometry/proximity/prox_sphere_ccylinder.hpp>

#include <cmath>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


shared_ptr< shape_3D > prox_sphere_ccylinder::getShape1() const {
  return mSphere;
};

shared_ptr< shape_3D > prox_sphere_ccylinder::getShape2() const {
  return mCCylinder;
};
    
void prox_sphere_ccylinder::computeProximity(const shape_3D_precompute_pack& aPack1, 
                                             const shape_3D_precompute_pack& aPack2) {
  if((!mSphere) || (!mCCylinder)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  using std::fabs; using std::sqrt;
  
  const pose_3D<double>& sp_pose = (aPack1.parent == mSphere.get() ? 
                                    aPack1.global_pose : aPack2.global_pose);
  const pose_3D<double>& cy_pose = (aPack1.parent == mSphere.get() ? 
                                    aPack2.global_pose : aPack1.global_pose);
  
  vect<double,3> sp_c = sp_pose.Position;
  
  vect<double,3> sp_c_rel = cy_pose.transformFromGlobal(sp_c);
  //double sp_c_rel_rad = sqrt(sp_c_rel[0] * sp_c_rel[0] + sp_c_rel[1] * sp_c_rel[1]);
  
  if(fabs(sp_c_rel[2]) <= 0.5 * mCCylinder->getLength()) {
    // The sphere is around the round side of the cylinder.
    //  this means the min-dist point is on the round shell of the cylinder in the direction of sphere center.
    vect<double,3> sp_c_proj = vect<double,3>(sp_c_rel[0],sp_c_rel[1],0.0);
    double sp_c_proj_d = norm_2(sp_c_proj);
    mLastResult.mPoint2 = cy_pose.transformToGlobal(vect<double,3>(0.0,0.0,sp_c_rel[2]) + sp_c_proj * (mCCylinder->getRadius() / sp_c_proj_d));
    mLastResult.mPoint1 = cy_pose.transformToGlobal(sp_c_rel - sp_c_proj * (mSphere->getRadius() / sp_c_proj_d));
    mLastResult.mDistance = sp_c_proj_d - mSphere->getRadius() - mCCylinder->getRadius();
  } else {
    // The sphere is above or below the capped cylinder.
    //  this boils down to a simple sphere-sphere proximity.
    double fact = 1.0;
    if(sp_c_rel[2] < 0.0)
      fact = -1.0;
    
    vect<double,3> cy_c2 = cy_pose.transformToGlobal(vect<double,3>(0.0,0.0,fact * 0.5 * mCCylinder->getLength()));
    vect<double,3> diff_cc = cy_c2 - sp_c;
    double dist_cc = norm_2(diff_cc);
    
    mLastResult.mDistance = dist_cc - mSphere->getRadius() - mCCylinder->getRadius();
    mLastResult.mPoint1 = sp_c + (mSphere->getRadius() / dist_cc) * diff_cc;
    mLastResult.mPoint2 = cy_c2 - (mCCylinder->getRadius() / dist_cc) * diff_cc;
  };
};


prox_sphere_ccylinder::prox_sphere_ccylinder(const shared_ptr< sphere >& aSphere,
                                             const shared_ptr< capped_cylinder >& aCCylinder) :
                                             proximity_finder_3D(),
                                             mSphere(aSphere),
                                             mCCylinder(aCCylinder) { };


void RK_CALL prox_sphere_ccylinder::save(ReaK::serialization::oarchive& A, unsigned int) const {
  proximity_finder_3D::save(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mSphere)
    & RK_SERIAL_SAVE_WITH_NAME(mCCylinder);
};

void RK_CALL prox_sphere_ccylinder::load(ReaK::serialization::iarchive& A, unsigned int) {
  proximity_finder_3D::load(A,proximity_finder_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mSphere)
    & RK_SERIAL_LOAD_WITH_NAME(mCCylinder);
};


};

};











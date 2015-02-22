
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

#include <ReaK/geometry/proximity/prox_sphere_cylinder.hpp>

#include <cmath>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


void prox_sphere_cylinder::computeProximity(const shape_3D_precompute_pack& aPack1, 
                                            const shape_3D_precompute_pack& aPack2) {
  if((!mSphere) || (!mCylinder)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  using std::fabs; using std::sqrt;
  
  const pose_3D<double>& sp_pose = (aPack1.parent == mSphere ? 
                                    aPack1.global_pose : aPack2.global_pose);
  const pose_3D<double>& cy_pose = (aPack1.parent == mSphere ? 
                                    aPack2.global_pose : aPack1.global_pose);
  
  vect<double,3> sp_c = sp_pose.Position;
  
  vect<double,3> sp_c_rel = cy_pose.transformFromGlobal(sp_c);
  double sp_c_rel_rad = sqrt(sp_c_rel[0] * sp_c_rel[0] + sp_c_rel[1] * sp_c_rel[1]);
  
  const double sp_rad = mSphere->getRadius();
  const double cy_len = mCylinder->getLength();
  const double cy_rad = mCylinder->getRadius();
  
  if(fabs(sp_c_rel[2]) <= 0.5 * cy_len) {
    // The sphere is around the round side of the cylinder.
    //  this means the min-dist point is on the round shell of the cylinder in the direction of sphere center.
    vect<double,3> sp_c_proj = vect<double,3>(sp_c_rel[0],sp_c_rel[1],0.0);
    double sp_c_proj_d = norm_2(sp_c_proj);
    mLastResult.mPoint2 = cy_pose.transformToGlobal(vect<double,3>(0.0,0.0,sp_c_rel[2]) + sp_c_proj * (cy_rad / sp_c_proj_d));
    mLastResult.mPoint1 = cy_pose.transformToGlobal(sp_c_rel - sp_c_proj * (sp_rad / sp_c_proj_d));
    mLastResult.mDistance = sp_c_proj_d - sp_rad - cy_rad;
  } else if(sp_c_rel_rad < cy_rad) {
    // The sphere is above or below the cylinder.
    //  this boils down to a simple plane-sphere proximity.
    double fact = 1.0;
    if(sp_c_rel[2] < 0.0)
      fact = -1.0;
    mLastResult.mPoint2 = cy_pose.transformToGlobal(vect<double,3>(sp_c_rel[0],sp_c_rel[1],fact * 0.5 * cy_len));
    mLastResult.mPoint1 = cy_pose.transformToGlobal(vect<double,3>(sp_c_rel[0],sp_c_rel[1],sp_c_rel[2] - fact * sp_rad));
    mLastResult.mDistance = fact * sp_c_rel[2] - 0.5 * cy_len - sp_rad;
  } else {
    // The sphere is outside the rims of the cylinder.
    //  this means the min-dist point is on the rim of the cylinder, in the direction of the sphere center.
    vect<double,3> sp_c_proj = vect<double,3>(sp_c_rel[0],sp_c_rel[1],0.0);
    double sp_c_proj_d = norm_2(sp_c_proj);
    double fact = 1.0;
    if(sp_c_rel[2] < 0.0)
      fact = -1.0;
    vect<double,3> rim_pt = (cy_rad / sp_c_proj_d) * sp_c_proj + vect<double,3>(0.0,0.0,fact * 0.5 * cy_len);
    mLastResult.mPoint2 = cy_pose.transformToGlobal(rim_pt);
    sp_c_proj = mLastResult.mPoint2 - sp_c;
    sp_c_proj_d = norm_2(sp_c_proj);
    mLastResult.mPoint1 = sp_c + (sp_rad / sp_c_proj_d) * sp_c_proj;
    mLastResult.mDistance = sp_c_proj_d - sp_rad;
  };
};


prox_sphere_cylinder::prox_sphere_cylinder(const sphere* aSphere,
                                           const cylinder* aCylinder) :
                                           proximity_finder_3D(),
                                           mSphere(aSphere),
                                           mCylinder(aCylinder) { };


};

};


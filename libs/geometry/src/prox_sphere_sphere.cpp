
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

#include <ReaK/geometry/proximity/prox_sphere_sphere.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


void prox_sphere_sphere::computeProximity(const shape_3D_precompute_pack& aPack1, 
                                          const shape_3D_precompute_pack& aPack2) {
  if((!mSphere1) || (!mSphere2)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  vect<double,3> c1 = (aPack1.parent == mSphere1 ? 
                       aPack1.global_pose.Position : aPack2.global_pose.Position);
  vect<double,3> c2 = (aPack1.parent == mSphere1 ? 
                       aPack2.global_pose.Position : aPack1.global_pose.Position);
  
  const double s1_rad = mSphere1->getRadius();
  const double s2_rad = mSphere2->getRadius();
  
  vect<double,3> diff_cc = c2 - c1;
  double dist_cc = norm_2(diff_cc);
  
  mLastResult.mDistance = dist_cc - s1_rad - s2_rad;
  mLastResult.mPoint1 = c1 + (s1_rad / dist_cc) * diff_cc;
  mLastResult.mPoint2 = c2 - (s2_rad / dist_cc) * diff_cc;
  
};


prox_sphere_sphere::prox_sphere_sphere(const sphere* aSphere1,
                                       const sphere* aSphere2) :
                                       proximity_finder_3D(),
                                       mSphere1(aSphere1),
                                       mSphere2(aSphere2) { };


};

};


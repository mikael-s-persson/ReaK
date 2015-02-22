
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

#include <ReaK/geometry/proximity/prox_sphere_box.hpp>

#include <ReaK/geometry/proximity/prox_fundamentals_3D.hpp>

#include <cmath>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


void prox_sphere_box::computeProximity(const shape_3D_precompute_pack& aPack1, 
                                       const shape_3D_precompute_pack& aPack2) {
  if((!mSphere) || (!mBox)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  using std::fabs; using std::sqrt;
  
  const pose_3D<double>& sp_pose = (aPack1.parent == mSphere ? 
                                    aPack1.global_pose : aPack2.global_pose);
  const pose_3D<double>& bx_pose = (aPack1.parent == mSphere ? 
                                    aPack2.global_pose : aPack1.global_pose);
  
  vect<double,3> sp_c = sp_pose.Position;
  
  proximity_record_3D bxpt_result = findProximityBoxToPoint(*mBox, bx_pose, sp_c);
  
  // add a sphere-sweep around the point-box solution.
  vect<double,3> diff_v = bxpt_result.mPoint1 - bxpt_result.mPoint2;
  double diff_d = norm_2(diff_v);
  const double sp_rad = mSphere->getRadius();
  if(bxpt_result.mDistance < 0.0)
    mLastResult.mPoint1 = bxpt_result.mPoint2 - (sp_rad / diff_d) * diff_v;
  else
    mLastResult.mPoint1 = bxpt_result.mPoint2 + (sp_rad / diff_d) * diff_v;
  mLastResult.mPoint2 = bxpt_result.mPoint1;
  mLastResult.mDistance = bxpt_result.mDistance - sp_rad;
  return;
};


prox_sphere_box::prox_sphere_box(const sphere* aSphere,
                                 const box* aBox) :
                                 proximity_finder_3D(),
                                 mSphere(aSphere),
                                 mBox(aBox) { };


};

};


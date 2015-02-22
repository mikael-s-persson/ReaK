
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

#include <ReaK/geometry/proximity/prox_ccylinder_box.hpp>

#include <ReaK/geometry/proximity/prox_fundamentals_3D.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


void prox_ccylinder_box::computeProximity(const shape_3D_precompute_pack& aPack1, 
                                          const shape_3D_precompute_pack& aPack2) {
  if((!mCCylinder) || (!mBox)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,3>(0.0,0.0,0.0);
    mLastResult.mPoint2 = vect<double,3>(0.0,0.0,0.0);
    return;
  };
  using std::fabs; using std::sqrt;
  
  const pose_3D<double>& cy_pose = (aPack1.parent == mCCylinder ? 
                                    aPack1.global_pose : aPack2.global_pose);
  const pose_3D<double>& bx_pose = (aPack1.parent == mCCylinder ? 
                                    aPack2.global_pose : aPack1.global_pose);
  
  
  vect<double,3> cy_c = cy_pose.Position;
  vect<double,3> cy_t = cy_pose.rotateToGlobal(vect<double,3>(0.0,0.0,1.0));
  
  proximity_record_3D bxln_result = findProximityBoxToLine(*mBox, bx_pose, cy_c, cy_t, 0.5 * mCCylinder->getLength());
  
  // add a sphere-sweep around the point-box solution.
  vect<double,3> diff_v = bxln_result.mPoint1 - bxln_result.mPoint2;
  double diff_d = norm_2(diff_v);
  if(bxln_result.mDistance < 0.0)
    mLastResult.mPoint1 = bxln_result.mPoint2 - (mCCylinder->getRadius() / diff_d) * diff_v;
  else
    mLastResult.mPoint1 = bxln_result.mPoint2 + (mCCylinder->getRadius() / diff_d) * diff_v;
  mLastResult.mPoint2 = bxln_result.mPoint1;
  mLastResult.mDistance = bxln_result.mDistance - mCCylinder->getRadius();
  return;
};


prox_ccylinder_box::prox_ccylinder_box(const capped_cylinder* aCCylinder,
                                       const box* aBox) :
                                       proximity_finder_3D(),
                                       mCCylinder(aCCylinder),
                                       mBox(aBox) { };

};

};



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

#include <ReaK/geometry/proximity/prox_circle_circle.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


void prox_circle_circle::computeProximity(const shape_2D_precompute_pack& aPack1, 
                                          const shape_2D_precompute_pack& aPack2) {
  if((!mCircle1) || (!mCircle2)) {
    mLastResult.mDistance = std::numeric_limits<double>::infinity();
    mLastResult.mPoint1 = vect<double,2>(0.0,0.0);
    mLastResult.mPoint2 = vect<double,2>(0.0,0.0);
    return;
  };
  
  const pose_2D<double>& c1_pose = (aPack1.parent == mCircle1 ? 
                                    aPack1.global_pose : aPack2.global_pose);
  const pose_2D<double>& c2_pose = (aPack1.parent == mCircle1 ? 
                                    aPack2.global_pose : aPack1.global_pose);
  
  vect<double,2> c1 = c1_pose.Position;
  vect<double,2> c2 = c2_pose.Position;
  
  vect<double,2> diff_cc = c2 - c1;
  double dist_cc = norm_2(diff_cc);
  
  mLastResult.mDistance = dist_cc - mCircle1->getRadius() - mCircle2->getRadius();
  mLastResult.mPoint1 = c1 + (mCircle1->getRadius() / dist_cc) * diff_cc;
  mLastResult.mPoint2 = c2 - (mCircle2->getRadius() / dist_cc) * diff_cc;
  
};


prox_circle_circle::prox_circle_circle(const circle* aCircle1,
                                       const circle* aCircle2) :
                                       proximity_finder_2D(),
                                       mCircle1(aCircle1),
                                       mCircle2(aCircle2) { };

};

};


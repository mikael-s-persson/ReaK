
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

#include <ReaK/geometry/proximity/prox_fundamentals_3D.hpp>

#include <ReaK/core/optimization/line_search.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


proximity_record_3D findProximityBoxToPoint(const shared_ptr< box >& aBox, const vect<double,3>& aPoint) {
  vect<double,3> pt_rel = aBox->getPose().transformFromGlobal(aPoint);
  
  bool in_x_range = ((pt_rel[0] > -0.5 * aBox->getDimensions()[0]) &&
                     (pt_rel[0] <  0.5 * aBox->getDimensions()[0]));
  bool in_y_range = ((pt_rel[1] > -0.5 * aBox->getDimensions()[1]) &&
                     (pt_rel[1] <  0.5 * aBox->getDimensions()[1]));
  bool in_z_range = ((pt_rel[2] > -0.5 * aBox->getDimensions()[2]) &&
                     (pt_rel[2] <  0.5 * aBox->getDimensions()[2]));
  bool is_inside = (in_x_range && in_y_range && in_z_range);
  if(is_inside) {
    // The point is inside the box.
    vect<double,3> bound_dists = vect<double,3>(0.5 * aBox->getDimensions()[0] - fabs(pt_rel[0]),
                                                0.5 * aBox->getDimensions()[1] - fabs(pt_rel[1]),
                                                0.5 * aBox->getDimensions()[2] - fabs(pt_rel[2]));
    if((bound_dists[0] <= bound_dists[1]) &&
       (bound_dists[0] <= bound_dists[2])) {
      in_x_range = false;
    } else if((bound_dists[1] <= bound_dists[0]) &&
              (bound_dists[1] <= bound_dists[2])) {
      in_y_range = false;
    } else {
      in_z_range = false;
    };
  }
  
  vect<double,3> corner_pt = 0.5 * aBox->getDimensions();
  if(in_x_range)
    corner_pt[0] = pt_rel[0];
  else if(pt_rel[0] < 0.0)
    corner_pt[0] = -corner_pt[0];
  if(in_y_range)
    corner_pt[1] = pt_rel[1];
  else if(pt_rel[1] < 0.0)
    corner_pt[1] = -corner_pt[1];
  if(in_z_range)
    corner_pt[2] = pt_rel[2];
  else if(pt_rel[2] < 0.0)
    corner_pt[2] = -corner_pt[2];
  
  proximity_record_3D result;
  result.mPoint1 = aBox->getPose().transformToGlobal(corner_pt);
  double diff_d = norm_2(corner_pt - pt_rel);
  result.mPoint2 = aPoint;
  result.mDistance = (is_inside ? -diff_d : diff_d);
  return result;
};


namespace detail {
  
  struct ProxBoxToLineFunctor {
    shared_ptr< box > mBox;
    vect<double,3> mCenter;
    vect<double,3> mTangent;
    proximity_record_3D* mResult;
    
    ProxBoxToLineFunctor(const shared_ptr< box >& aBox, 
                         const vect<double,3>& aCenter, 
                         const vect<double,3>& aTangent, 
                         proximity_record_3D& aResult) :
                         mBox(aBox), mCenter(aCenter), mTangent(aTangent), mResult(&aResult) { };
    
    double operator()(double t) const {
      (*mResult) = findProximityBoxToPoint(mBox, mCenter + mTangent * t);
      return mResult->mDistance;
    };
    
  };
  
};


proximity_record_3D findProximityBoxToLine(const shared_ptr< box >& aBox, const vect<double,3>& aCenter, const vect<double,3>& aTangent, double aHalfLength) {
  proximity_record_3D result;
  detail::ProxBoxToLineFunctor fct(aBox, aCenter, aTangent, result);
  double lb = -aHalfLength;
  double ub = aHalfLength;
  optim::golden_section_search(fct, lb, ub, 1e-3 * aHalfLength);
  return result;  // the result of the search should be found in the result object (filled in by 'fct').
};

  


};

};


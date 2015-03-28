
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


proximity_record_3D compute_proximity( const sphere& aSphere, const shape_3D_precompute_pack& aPack1,
                                       const capped_cylinder& aCCylinder, const shape_3D_precompute_pack& aPack2 ) {
  using std::fabs;
  using std::sqrt;
  proximity_record_3D result;

  const pose_3D< double >& sp_pose = aPack1.global_pose;
  const pose_3D< double >& cy_pose = aPack2.global_pose;

  vect< double, 3 > sp_c = sp_pose.Position;

  vect< double, 3 > sp_c_rel = cy_pose.transformFromGlobal( sp_c );
  // double sp_c_rel_rad = sqrt(sp_c_rel[0] * sp_c_rel[0] + sp_c_rel[1] * sp_c_rel[1]);

  const double sp_rad = aSphere.getRadius();
  const double cy_len = aCCylinder.getLength();
  const double cy_rad = aCCylinder.getRadius();

  if( fabs( sp_c_rel[2] ) <= 0.5 * cy_len ) {
    // The sphere is around the round side of the cylinder.
    //  this means the min-dist point is on the round shell of the cylinder in the direction of sphere center.
    vect< double, 3 > sp_c_proj = vect< double, 3 >( sp_c_rel[0], sp_c_rel[1], 0.0 );
    double sp_c_proj_d = norm_2( sp_c_proj );
    result.mPoint2
      = cy_pose.transformToGlobal( vect< double, 3 >( 0.0, 0.0, sp_c_rel[2] ) + sp_c_proj * ( cy_rad / sp_c_proj_d ) );
    result.mPoint1 = cy_pose.transformToGlobal( sp_c_rel - sp_c_proj * ( sp_rad / sp_c_proj_d ) );
    result.mDistance = sp_c_proj_d - sp_rad - cy_rad;
  } else {
    // The sphere is above or below the capped cylinder.
    //  this boils down to a simple sphere-sphere proximity.
    double fact = 1.0;
    if( sp_c_rel[2] < 0.0 )
      fact = -1.0;

    vect< double, 3 > cy_c2 = cy_pose.transformToGlobal( vect< double, 3 >( 0.0, 0.0, fact * 0.5 * cy_len ) );
    vect< double, 3 > diff_cc = cy_c2 - sp_c;
    double dist_cc = norm_2( diff_cc );

    result.mDistance = dist_cc - sp_rad - cy_rad;
    result.mPoint1 = sp_c + ( sp_rad / dist_cc ) * diff_cc;
    result.mPoint2 = cy_c2 - ( cy_rad / dist_cc ) * diff_cc;
  };
  return result;
};

proximity_record_3D compute_proximity( const capped_cylinder& aCCylinder, const shape_3D_precompute_pack& aPack1,
                                       const sphere& aSphere, const shape_3D_precompute_pack& aPack2 ) {
  using std::swap;
  proximity_record_3D result = compute_proximity( aSphere, aPack2, aCCylinder, aPack1 );
  swap( result.mPoint1, result.mPoint2 );
  return result;
};

proximity_record_3D prox_sphere_ccylinder::computeProximity( const shape_3D_precompute_pack& aPack1,
                                                             const shape_3D_precompute_pack& aPack2 ) {
  if( ( !mSphere ) || ( !mCCylinder ) )
    return proximity_record_3D();

  if( aPack1.parent == mSphere )
    return compute_proximity( *mSphere, aPack1, *mCCylinder, aPack2 );
  else
    return compute_proximity( *mCCylinder, aPack1, *mSphere, aPack2 );
};


prox_sphere_ccylinder::prox_sphere_ccylinder( const sphere* aSphere, const capped_cylinder* aCCylinder )
    : proximity_finder_3D(), mSphere( aSphere ), mCCylinder( aCCylinder ){};
};
};

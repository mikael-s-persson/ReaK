
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

#include <ReaK/math/optimization/line_search.hpp>

#include <vector>
#include <queue>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {

proximity_record_3D findProximityBoxToPoint( const box& aBox, const pose_3D< double >& aBoxGblPose,
                                             const vect< double, 3 >& aPoint ) {
  const vect< double, 3 > pt_rel = aBoxGblPose.transformFromGlobal( aPoint );
  const vect< double, 3 > bx_dim = aBox.getDimensions();

  bool in_x_range = ( ( pt_rel[0] > -0.5 * bx_dim[0] ) && ( pt_rel[0] < 0.5 * bx_dim[0] ) );
  bool in_y_range = ( ( pt_rel[1] > -0.5 * bx_dim[1] ) && ( pt_rel[1] < 0.5 * bx_dim[1] ) );
  bool in_z_range = ( ( pt_rel[2] > -0.5 * bx_dim[2] ) && ( pt_rel[2] < 0.5 * bx_dim[2] ) );
  bool is_inside = ( in_x_range && in_y_range && in_z_range );
  if( is_inside ) {
    // The point is inside the box.
    vect< double, 3 > bound_dists = vect< double, 3 >(
      0.5 * bx_dim[0] - fabs( pt_rel[0] ), 0.5 * bx_dim[1] - fabs( pt_rel[1] ), 0.5 * bx_dim[2] - fabs( pt_rel[2] ) );
    if( ( bound_dists[0] <= bound_dists[1] ) && ( bound_dists[0] <= bound_dists[2] ) ) {
      in_x_range = false;
    } else if( ( bound_dists[1] <= bound_dists[0] ) && ( bound_dists[1] <= bound_dists[2] ) ) {
      in_y_range = false;
    } else {
      in_z_range = false;
    };
  }

  vect< double, 3 > corner_pt = 0.5 * bx_dim;
  if( in_x_range )
    corner_pt[0] = pt_rel[0];
  else if( pt_rel[0] < 0.0 )
    corner_pt[0] = -corner_pt[0];
  if( in_y_range )
    corner_pt[1] = pt_rel[1];
  else if( pt_rel[1] < 0.0 )
    corner_pt[1] = -corner_pt[1];
  if( in_z_range )
    corner_pt[2] = pt_rel[2];
  else if( pt_rel[2] < 0.0 )
    corner_pt[2] = -corner_pt[2];

  proximity_record_3D result;
  result.mPoint1 = aBoxGblPose.transformToGlobal( corner_pt );
  double diff_d = norm_2( corner_pt - pt_rel );
  result.mPoint2 = aPoint;
  result.mDistance = ( is_inside ? -diff_d : diff_d );
  return result;
};


namespace detail {

struct ProxBoxToLineFunctor {
  const box* mBox;
  const pose_3D< double >* mBoxGblPose;
  vect< double, 3 > mCenter;
  vect< double, 3 > mTangent;
  proximity_record_3D* mResult;

  ProxBoxToLineFunctor( const box& aBox, const pose_3D< double >& aBoxGblPose, const vect< double, 3 >& aCenter,
                        const vect< double, 3 >& aTangent, proximity_record_3D& aResult )
      : mBox( &aBox ), mBoxGblPose( &aBoxGblPose ), mCenter( aCenter ), mTangent( aTangent ), mResult( &aResult ){};

  double operator()( double t ) const {
    ( *mResult ) = findProximityBoxToPoint( *mBox, *mBoxGblPose, mCenter + mTangent * t );
    return mResult->mDistance;
  };
};
};


proximity_record_3D findProximityBoxToLine( const box& aBox, const pose_3D< double >& aBoxGblPose,
                                            const vect< double, 3 >& aCenter, const vect< double, 3 >& aTangent,
                                            double aHalfLength ) {
  proximity_record_3D result;
  detail::ProxBoxToLineFunctor fct( aBox, aBoxGblPose, aCenter, aTangent, result );
  double lb = -aHalfLength;
  double ub = aHalfLength;
  optim::golden_section_search( fct, lb, ub, 1e-3 * aHalfLength );
  return result; // the result of the search should be found in the result object (filled in by 'fct').
};


vect< double, 3 > cylinder_boundary_func::operator()( vect< double, 3 > v ) {
  using std::sqrt;
  using std::fabs;
  v = mGblPose->rotateFromGlobal( v );
  const double v_radius = sqrt( v[0] * v[0] + v[1] * v[1] );
  const double cy_len = mCylinder->getLength();
  const double cy_rad = mCylinder->getRadius();
  if( v_radius * 0.5 * cy_len > cy_rad * fabs( v[2] ) ) {
    v *= cy_rad / v_radius;
  } else {
    v *= 0.5 * cy_len / fabs( v[2] );
  };
  return mGblPose->transformToGlobal( v );
};

mat< double, mat_structure::square > cylinder_boundary_jac::operator()( vect< double, 3 > v ) {
  using std::sqrt;
  using std::fabs;
  mat< double, mat_structure::square > R = mGblPose->Quat.getMat();
  const double cy_len = mCylinder->getLength();
  const double cy_rad = mCylinder->getRadius();
  // Jac_gbl = R * Jac_local * R^T

  // Jac_local:
  v = mGblPose->rotateFromGlobal( v );
  double v_radius = sqrt( v[0] * v[0] + v[1] * v[1] );
  mat< double, mat_structure::square > Jac_local( 3, 0.0 );
  if( v_radius * 0.5 * cy_len > cy_rad * fabs( v[2] ) ) {
    // v *= cy_rad / v_radius;
    double v_rad3 = v_radius * v_radius * v_radius;
    Jac_local( 0, 0 ) = cy_rad * v[1] * v[1] / v_rad3;
    Jac_local( 0, 1 ) = -cy_rad * v[0] * v[1] / v_rad3;
    Jac_local( 1, 0 ) = -cy_rad * v[0] * v[1] / v_rad3;
    Jac_local( 1, 1 ) = cy_rad * v[0] * v[0] / v_rad3;
    Jac_local( 2, 0 ) = -cy_rad * v[0] * v[2] / v_rad3;
    Jac_local( 2, 1 ) = -cy_rad * v[1] * v[2] / v_rad3;
    Jac_local( 2, 2 ) = cy_rad / v_radius;
  } else {
    // v *= 0.5 * cy_len / fabs(v[2]);
    double v_z3 = v[2] * v[2] * fabs( v[2] );
    Jac_local( 0, 0 ) = 0.5 * cy_len / fabs( v[2] );
    Jac_local( 0, 2 ) = -0.5 * cy_len * v[0] * v[2] / v_z3;
    Jac_local( 1, 1 ) = 0.5 * cy_len / fabs( v[2] );
    Jac_local( 1, 2 ) = -0.5 * cy_len * v[1] * v[2] / v_z3;
    //       Jac_local(2,2) =  0.0;
  };

  // Jac_gbl:
  return mat< double, mat_structure::square >( R * Jac_local * transpose_view( R ) );
};

vect< double, 3 > ccylinder_boundary_func::operator()( vect< double, 3 > v ) {
  using std::sqrt;
  using std::fabs;
  const double cy_len = mCCylinder->getLength();
  const double cy_rad = mCCylinder->getRadius();
  v = mGblPose->rotateFromGlobal( v );
  double v_radius = sqrt( v[0] * v[0] + v[1] * v[1] );
  if( v_radius * 0.5 * cy_len > cy_rad * fabs( v[2] ) ) {
    v *= cy_rad / v_radius;
  } else {
    double v_l2 = v_radius * v_radius + v[2] * v[2];
    double A
      = sqrt( ( cy_rad * cy_rad - 0.25 * cy_len * cy_len ) * v_radius * v_radius + cy_rad * cy_rad * v[2] * v[2] );
    v *= ( 0.5 * cy_len * fabs( v[2] ) + A ) / v_l2;
  };
  return mGblPose->transformToGlobal( v );
};

mat< double, mat_structure::square > ccylinder_boundary_jac::operator()( vect< double, 3 > v ) {
  using std::sqrt;
  using std::fabs;
  const double cy_len = mCCylinder->getLength();
  const double cy_rad = mCCylinder->getRadius();
  mat< double, mat_structure::square > R = mGblPose->Quat.getMat();
  // Jac_gbl = R * Jac_local * R^T

  // Jac_local:
  v = mGblPose->rotateFromGlobal( v );
  double v_radius = sqrt( v[0] * v[0] + v[1] * v[1] );
  mat< double, mat_structure::square > Jac_local( 3, 0.0 );
  if( v_radius * 0.5 * cy_len > cy_rad * fabs( v[2] ) ) {
    // v *= cy_rad / v_radius;
    double v_rad3 = v_radius * v_radius * v_radius;
    Jac_local( 0, 0 ) = cy_rad * v[1] * v[1] / v_rad3;
    Jac_local( 0, 1 ) = -cy_rad * v[0] * v[1] / v_rad3;
    Jac_local( 1, 0 ) = -cy_rad * v[0] * v[1] / v_rad3;
    Jac_local( 1, 1 ) = cy_rad * v[0] * v[0] / v_rad3;
    Jac_local( 2, 0 ) = -cy_rad * v[0] * v[2] / v_rad3;
    Jac_local( 2, 1 ) = -cy_rad * v[1] * v[2] / v_rad3;
    Jac_local( 2, 2 ) = cy_rad / v_radius;
  } else {
    double v_l2 = v_radius * v_radius + v[2] * v[2];
    double R2 = cy_rad * cy_rad;
    double R_L = ( R2 - 0.25 * cy_len * cy_len );
    double A = sqrt( R_L * v_radius * v_radius + R2 * v[2] * v[2] );
    double Lv_A = 0.5 * cy_len * fabs( v[2] ) + A;
    double F1 = R_L / A - 2.0 * Lv_A / v_l2;
    double F2 = R2 / A + 0.5 * cy_len / fabs( v[2] ) - 2.0 * Lv_A / v_l2;
    //       v *= Lv_A / v_l2;

    Jac_local( 0, 0 ) = v[0] * v[0] / v_l2 * F1 + Lv_A / v_l2;
    Jac_local( 0, 1 ) = v[0] * v[1] / v_l2 * F1;
    Jac_local( 0, 2 ) = v[0] * v[2] / v_l2 * F2;

    Jac_local( 1, 0 ) = v[1] * v[0] / v_l2 * F1;
    Jac_local( 1, 1 ) = v[1] * v[1] / v_l2 * F1 + Lv_A / v_l2;
    Jac_local( 1, 2 ) = v[1] * v[2] / v_l2 * F2;

    Jac_local( 2, 0 ) = v[2] * v[0] / v_l2 * F1;
    Jac_local( 2, 1 ) = v[2] * v[1] / v_l2 * F1;
    Jac_local( 2, 2 ) = v[2] * v[2] / v_l2 * F2 + Lv_A / v_l2;
  };

  // Jac_gbl:
  return mat< double, mat_structure::square >( R * Jac_local * transpose_view( R ) );
};

vect< double, 3 > box_boundary_func::operator()( vect< double, 3 > v ) {
  using std::fabs;
  const vect< double, 3 > bx_dim = mBox->getDimensions();
  v = mGblPose->rotateFromGlobal( v );
  // there are three possible quadrants (facet where the support lays):
  if( ( fabs( v[1] ) * bx_dim[0] < fabs( v[0] ) * bx_dim[1] )
      && ( fabs( v[2] ) * bx_dim[0] < fabs( v[0] ) * bx_dim[2] ) ) {
    // on the x-normal plane:
    v *= 0.5 * bx_dim[0] / fabs( v[0] );
  } else if( ( fabs( v[0] ) * bx_dim[1] < fabs( v[1] ) * bx_dim[0] )
             && ( fabs( v[2] ) * bx_dim[1] < fabs( v[1] ) * bx_dim[2] ) ) {
    // on the y-normal plane:
    v *= 0.5 * bx_dim[1] / fabs( v[1] );
  } else {
    // on the z-normal plane:
    v *= 0.5 * bx_dim[2] / fabs( v[2] );
  };
  return mGblPose->transformToGlobal( v );
};

mat< double, mat_structure::square > box_boundary_jac::operator()( vect< double, 3 > v ) {
  using std::sqrt;
  using std::fabs;
  const vect< double, 3 > bx_dim = mBox->getDimensions();
  pose_3D< double > gbl_pose = mGblPose->getGlobalPose();
  mat< double, mat_structure::square > R = gbl_pose.Quat.getMat();
  // Jac_gbl = R * Jac_local * R^T

  // Jac_local:
  v = gbl_pose.rotateFromGlobal( v );
  mat< double, mat_structure::square > Jac_local( 3, 0.0 );
  if( ( fabs( v[1] ) * bx_dim[0] < fabs( v[0] ) * bx_dim[1] )
      && ( fabs( v[2] ) * bx_dim[0] < fabs( v[0] ) * bx_dim[2] ) ) {
    // on the x-normal plane:
    //       v *= 0.5 * bx_dim[0] / fabs(v[0]);
    double v_x3 = v[0] * v[0] * fabs( v[0] );
    //       Jac_local(0,0) =  0.0;
    Jac_local( 1, 0 ) = -0.5 * bx_dim[0] * v[1] * v[0] / v_x3;
    Jac_local( 1, 1 ) = 0.5 * bx_dim[0] / fabs( v[0] );
    Jac_local( 2, 0 ) = -0.5 * bx_dim[0] * v[2] * v[0] / v_x3;
    Jac_local( 2, 2 ) = 0.5 * bx_dim[0] / fabs( v[0] );
  } else if( ( fabs( v[0] ) * bx_dim[1] < fabs( v[1] ) * bx_dim[0] )
             && ( fabs( v[2] ) * bx_dim[1] < fabs( v[1] ) * bx_dim[2] ) ) {
    // on the y-normal plane:
    //       v *= 0.5 * bx_dim[1] / fabs(v[1]);
    double v_y3 = v[1] * v[1] * fabs( v[1] );
    Jac_local( 0, 0 ) = 0.5 * bx_dim[1] / fabs( v[1] );
    Jac_local( 0, 1 ) = -0.5 * bx_dim[1] * v[0] * v[1] / v_y3;
    //       Jac_local(1,1) =  0.0;
    Jac_local( 2, 1 ) = -0.5 * bx_dim[1] * v[2] * v[1] / v_y3;
    Jac_local( 2, 2 ) = 0.5 * bx_dim[1] / fabs( v[1] );
  } else {
    // on the z-normal plane:
    //       v *= 0.5 * bx_dim[2] / fabs(v[2]);
    double v_z3 = v[2] * v[2] * fabs( v[2] );
    Jac_local( 0, 0 ) = 0.5 * bx_dim[2] / fabs( v[2] );
    Jac_local( 0, 2 ) = -0.5 * bx_dim[2] * v[0] * v[2] / v_z3;
    Jac_local( 1, 1 ) = 0.5 * bx_dim[2] / fabs( v[2] );
    Jac_local( 1, 2 ) = -0.5 * bx_dim[2] * v[1] * v[2] / v_z3;
    //       Jac_local(2,2) =  0.0;
  };

  // Jac_gbl:
  return mat< double, mat_structure::square >( R * Jac_local * transpose_view( R ) );
};


#define REAK_PROX_SUPPORTS_USE_DISCRETE_APPROX

vect< double, 3 > cylinder_support_func::operator()( vect< double, 3 > v ) const {
  using std::sqrt;
  using std::fabs;
  const double cy_len = mCylinder->getLength();
  const double cy_rad = mCylinder->getRadius();
  v = mGblPose->rotateFromGlobal( v );
  double v_radius = sqrt( v[0] * v[0] + v[1] * v[1] );
  if( v_radius > 1e-5 ) {
#ifndef REAK_PROX_SUPPORTS_USE_DISCRETE_APPROX
    v[0] *= cy_rad / v_radius;
    v[1] *= cy_rad / v_radius;
#else
    using std::atan2;
    using std::sin;
    using std::cos;
    using std::round;
    // octogonal facets: (octogon vs. circle is pretty OK approximation, with 8-pts vs. infinite).
    double theta = round( atan2( v[1], v[0] ) * 8.0 / M_PI ) * M_PI / 8.0;
    v[0] = cos( theta ) * cy_rad;
    v[1] = sin( theta ) * cy_rad;
#endif
    if( fabs( v[2] ) > 1e-5 ) {
      v[2] *= 0.5 * cy_len / fabs( v[2] );
    } else {
      v[2] = 0.0; // at the center of cylinder.
    };
  } else {
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] *= 0.5 * cy_len / fabs( v[2] );
  };
  return mGblPose->transformToGlobal( v );
};

vect< double, 3 > ccylinder_support_func::operator()( vect< double, 3 > v ) const {
  using std::sqrt;
  using std::fabs;
  const double cy_len = mCCylinder->getLength();
  const double cy_rad = mCCylinder->getRadius();
  v = mGblPose->rotateFromGlobal( v );
  double v_radius = sqrt( v[0] * v[0] + v[1] * v[1] );
  if( v_radius > 1e-5 ) {
#ifndef REAK_PROX_SUPPORTS_USE_DISCRETE_APPROX
    v[0] *= cy_rad / v_radius;
    v[1] *= cy_rad / v_radius;
#else
    using std::atan2;
    using std::sin;
    using std::cos;
    using std::round;
    // octogonal facets: (octogon vs. circle is pretty OK approximation, with 8-pts vs. infinite).
    double theta = round( atan2( v[1], v[0] ) * 8.0 / M_PI ) * M_PI / 8.0;
    v[0] = cos( theta ) * cy_rad;
    v[1] = sin( theta ) * cy_rad;
#endif
    if( fabs( v[2] ) > 1e-5 ) {
      double phi = atan2( v[2], v_radius );
#ifdef REAK_PROX_SUPPORTS_USE_DISCRETE_APPROX
      phi = round( phi * 8.0 / M_PI ) * M_PI / 8.0;
#endif
      v[0] *= cos( phi );
      v[1] *= cos( phi );
      v[2] *= 0.5 * cy_len / fabs( v[2] );
      v[2] += sin( phi ) * cy_rad;
    } else {
      v[2] = 0.0; // at the center of cylinder.
    };
  } else {
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] *= 0.5 * cy_len / fabs( v[2] );
  };
  return mGblPose->transformToGlobal( v );
};

vect< double, 3 > box_support_func::operator()( vect< double, 3 > v ) const {
  using std::fabs;
  const vect< double, 3 > bx_dim = mBox->getDimensions();
  v = mGblPose->rotateFromGlobal( v );

  if( fabs( v[0] ) > 1e-5 )
    v[0] *= 0.5 * bx_dim[0] / fabs( v[0] );
  else
    v[0] = 0.0;
  if( fabs( v[1] ) > 1e-5 )
    v[1] *= 0.5 * bx_dim[1] / fabs( v[1] );
  else
    v[1] = 0.0;
  if( fabs( v[2] ) > 1e-5 )
    v[2] *= 0.5 * bx_dim[2] / fabs( v[2] );
  else
    v[2] = 0.0;

  return mGblPose->transformToGlobal( v );
};


namespace {

struct GJK_simplex_point {
  vect< double, 3 > dir;
  vect< double, 3 > s1;
  vect< double, 3 > s2;
  vect< double, 3 > cso;
  double dist;

  void init( const vect< double, 3 >& aDir, const support_func_base& aSFunc1, const support_func_base& aSFunc2 ) {
    dir = aDir;
    s1 = aSFunc1( aDir );
    s2 = aSFunc2( -aDir );
    cso = s1 - s2;
    dist = norm_2( cso );
  };
};

enum EPA_color { EPA_black, EPA_grey, EPA_white };

struct EPA_edge_info {
  int v[2];
  int f[2];
  vect< double, 3 > min_pt;
  double min_dist;
  int min_vert;
  EPA_color color;

  EPA_edge_info( const std::vector< GJK_simplex_point >& aVec, int aV1, int aV2 ) {
    v[0] = aV1;
    v[1] = aV2;
    f[0] = -1;
    f[1] = -1;
    color = EPA_white;
    vect< double, 3 > line_v = aVec[aV1].cso - aVec[aV2].cso;
    double line_length = norm_2( line_v );
    if( line_length < 1e-5 ) {
      min_pt = aVec[aV1].cso;
      min_dist = aVec[aV1].dist;
      min_vert = 0;
      return;
    };
    line_v *= 1.0 / line_length;
    double i_proj = aVec[aV1].cso * line_v;
    if( i_proj > line_length ) {
      min_pt = aVec[aV2].cso;
      min_dist = aVec[aV2].dist;
      min_vert = 1;
    } else if( i_proj < 0.0 ) {
      min_pt = aVec[aV1].cso;
      min_dist = aVec[aV1].dist;
      min_vert = 0;
    } else {
      min_pt = aVec[aV1].cso - i_proj * line_v;
      min_dist = norm_2( min_pt );
      min_vert = -1; // on the edge
    };
  };
};

vect< double, 3 > proj_barycentric_coord( const vect< double, 3 >& v1, vect< double, 3 > v2, vect< double, 3 > v3 ) {
  vect< double, 3 > b;
  v2 -= v1;
  v3 -= v1;
  vect< double, 3 > n = v2 % v3;
  double fact = 1.0 / ( n * n );
  b[2] = ( ( v1 % v2 ) * n ) * fact;
  b[1] = ( ( v3 % v1 ) * n ) * fact;
  b[0] = 1.0 - b[1] - b[2];
  return b;
};

struct EPA_face_info {
  int v[3];
  int e[3];
  vect< double, 3 > min_pt;
  double min_dist;
  int min_edge;
  EPA_color color;

  EPA_face_info( int aId, const std::vector< GJK_simplex_point >& aVec, std::vector< EPA_edge_info >& aEd, int aE1,
                 int aE2, int aE3 ) {
    e[0] = aE1;
    e[1] = aE2;
    e[2] = aE3;
    color = EPA_white;

    for( int i = 0; i < 3; ++i ) {
      if( aEd[e[i]].f[0] < 0 )
        aEd[e[i]].f[0] = aId;
      else
        aEd[e[i]].f[1] = aId;
    };

    int aV3 = -1;
    int aVCommon = -1;
    if( ( aEd[aE2].v[0] != aEd[aE1].v[0] ) && ( aEd[aE2].v[0] != aEd[aE1].v[1] ) ) {
      aV3 = aEd[aE2].v[0];
      aVCommon = aEd[aE2].v[1];
    } else {
      aV3 = aEd[aE2].v[1];
      aVCommon = aEd[aE2].v[0];
    };
    v[0] = aEd[aE1].v[0];
    v[1] = aEd[aE1].v[1];
    v[2] = aV3;

    //       std::cout << " - Adding face(" << v[0] << " ," << v[1] << " ," << v[2] << ")..." << std::endl;
    //       std::cout << "   with edge min-pt = " << aEd[aE1].min_pt
    //                 << " and sommet = " << aVec[aV3].cso << std::endl;

    vect< double, 3 > b = proj_barycentric_coord( aVec[v[0]].cso, aVec[v[1]].cso, aVec[v[2]].cso );

    if( b[0] < 1e-5 ) {
      // on edge opposite to v[0]:
      if( aVCommon == v[0] ) {
        min_pt = aEd[aE3].min_pt;
        min_dist = aEd[aE3].min_dist;
        min_edge = 2;
      } else {
        min_pt = aEd[aE2].min_pt;
        min_dist = aEd[aE2].min_dist;
        min_edge = 1;
      };
    } else if( b[1] < 1e-5 ) {
      // on edge opposite to v[1]:
      if( aVCommon == v[1] ) {
        min_pt = aEd[aE3].min_pt;
        min_dist = aEd[aE3].min_dist;
        min_edge = 2;
      } else {
        min_pt = aEd[aE2].min_pt;
        min_dist = aEd[aE2].min_dist;
        min_edge = 1;
      };
    } else if( b[2] < 1e-5 ) {
      // on edge opposite to v[2]:
      min_pt = aEd[aE1].min_pt;
      min_dist = aEd[aE1].min_dist;
      min_edge = 0;
    } else {
      // inside of the triangle:
      min_pt = b[0] * aVec[v[0]].cso + b[1] * aVec[v[1]].cso + b[2] * aVec[v[2]].cso;
      min_dist = norm_2( min_pt );
      min_edge = -1;
    };
  };
};

int EPA_add_face( std::vector< EPA_face_info >& epa_faces, std::vector< int >& graveyard,
                  const std::vector< GJK_simplex_point >& epa_pts, std::vector< EPA_edge_info >& epa_edges, int e1,
                  int e2, int e3 ) {
  int result;
  if( graveyard.empty() ) {
    // create a new face:
    result = static_cast<int>(epa_faces.size());
    epa_faces.push_back( EPA_face_info( result, epa_pts, epa_edges, e1, e2, e3 ) );
  } else {
    // resurrect a face from the graveyard:
    result = graveyard.back();
    epa_faces[result] = EPA_face_info( result, epa_pts, epa_edges, e1, e2, e3 );
    graveyard.pop_back();
  };
  return result;
};

struct EPA_face_entry {
  double dist;
  int f;
  EPA_face_entry( double aDist, int aF ) : dist( aDist ), f( aF ){};
  bool operator<( const EPA_face_entry& rhs ) const {
    return dist > rhs.dist; // less-than is greater-than, for min-heap.
  };
};
};


proximity_record_3D findProximityByGJKEPA( const support_func_base& aSFunc1, const support_func_base& aSFunc2 ) {
  using std::fabs;
  GJK_simplex_point simplex[5];
  int cur_size = 0;

  // First, do GJK:
  simplex[cur_size].init( vect< double, 3 >( 1.0, 0.0, 0.0 ), aSFunc1, aSFunc2 );
  ++cur_size;
  vect< double, 3 > closest;
  double closest_dist = std::numeric_limits< double >::infinity();
  int closest_case = -1;

  double last_closest_dist = std::numeric_limits< double >::infinity();
  int last_cur_size = 1;

  while( true ) {
    // Find the closest point (test each point, line, triangle of the simplex):
    closest_dist = std::numeric_limits< double >::infinity();
    closest_case = -1;
    vect< double, 3 > facet_proj_pts[4];
    int facet_count = 0;
    // cases are: 0, 1, 2, 3 for closest point index;
    //            16, 32, 36, 48, 52, 56 for closest point on edge;
    //            2304, 3328, 3584, 3648 for closest point in triangle.
    for( int i = 0; i < cur_size; ++i ) {
      // Test point:
      if( simplex[i].dist < closest_dist ) {
        closest = simplex[i].cso;
        closest_dist = simplex[i].dist;
        closest_case = i;
      };
      for( int j = i + 1; j < cur_size; ++j ) {
        // Test line:
        vect< double, 3 > line_v = simplex[i].cso - simplex[j].cso;
        double line_length = norm_2( line_v );
        if( line_length < 1e-5 )
          continue;
        line_v *= 1.0 / line_length;
        double i_proj = simplex[i].cso * line_v;
        if( ( i_proj > line_length ) || ( i_proj < 0.0 ) )
          continue; // closest point not on this line.
        vect< double, 3 > closest_on_line = simplex[i].cso - i_proj * line_v;
        double dist_from_line = norm_2( closest_on_line );
        if( dist_from_line < closest_dist ) {
          closest = closest_on_line;
          closest_dist = dist_from_line;
          closest_case = ( i + ( j << 2 ) ) << 2;
        };
        for( int k = j + 1; k < cur_size; ++k ) {
          // Test facet:
          vect< double, 3 > b = proj_barycentric_coord( simplex[i].cso, simplex[j].cso, simplex[k].cso );

          if( ( b[0] < 1e-5 ) || ( b[1] < 1e-5 ) || ( b[2] < 1e-5 ) )
            continue; // closest point not on this facet.

          vect< double, 3 > closest_on_facet = b[0] * simplex[i].cso + b[1] * simplex[j].cso + b[2] * simplex[k].cso;
          facet_proj_pts[facet_count] = closest_on_facet;
          ++facet_count;
          double dist_from_facet = norm_2( closest_on_facet );
          if( dist_from_facet < closest_dist ) {
            closest = closest_on_facet;
            closest_dist = dist_from_facet;
            closest_case = ( i + ( j << 2 ) + ( k << 4 ) ) << 6;
          };
        };
      };
    };

    // Check if the origin is inside the current simplex:
    if( facet_count >= 3 ) { // it's impossible for the origin to be inside the simplex if < 2 facets.
      // Test tetrahedron:
      for( int i = 0; i < facet_count; ++i ) {
        for( int j = i + 1; j < facet_count; ++j ) {
          double d = facet_proj_pts[i] * facet_proj_pts[j];
          if( d < 0.0 ) {
            // this means that the origin is inside the tetrahedron.
            closest = vect< double, 3 >( 0.0, 0.0, 0.0 );
            closest_dist = 0.0;
            closest_case = 1 << 12;
          };
        };
      };
    };

    // Check if the origin is the closest point:
    if( closest_dist < 1e-5 ) {
      // Colliding, exit GJK
      break;
    };

    bool something_changed = ( fabs( last_closest_dist - closest_dist ) > 1e-5 );
    // Eliminate any points not part of the closest point / line / facet.
    if( closest_case >= 64 )
      closest_case >>= 6;
    else if( closest_case >= 4 )
      closest_case >>= 2;
    for( int i = 0; i < cur_size; ++i ) {
      int j = closest_case & 3;
      if( i != j ) {
        simplex[i] = simplex[j];
        something_changed = true;
      };
      closest_case >>= 2;
      if( closest_case == 0 ) {
        cur_size = i + 1;
        break;
      };
    };

    simplex[cur_size].init( -closest, aSFunc1, aSFunc2 );

    something_changed = ( something_changed || ( last_cur_size != cur_size + 1 ) );

    //     std::cout << "- Current closest: " <<  closest << std::endl;
    //     std::cout << "    at distance: " <<  closest_dist << std::endl;
    //     std::cout << "- Current support simplex: " << std::endl;
    //     for(int i = 0; i <= cur_size; ++i)
    //       std::cout << "-    " << simplex[i].cso << std::endl;

    if( !something_changed || ( closest - simplex[cur_size].cso ) * closest < 1e-5 ) {
      // Non-colliding, exit GJK
      break;
    };

    ++cur_size;

    last_closest_dist = closest_dist;
    last_cur_size = cur_size;
  };

  proximity_record_3D result;

  // If we have a non-zero distance, then it means that we have no collision
  // and therefore, the simplex should be the last primitive (point, line, facet)
  // touching the closest distance, and should be on the edge of the minkowski sum/diff.
  // If there is only a single point coincident with the origin, then it's the point of collision.
  // The closest point represents the min-dist vector:
  if( cur_size == 1 ) {
    result.mDistance = closest_dist;
    result.mPoint1 = simplex[0].s1;
    result.mPoint2 = simplex[0].s2;
    return result;
  };
  if( ( closest_dist > 1e-5 ) && ( cur_size == 2 ) ) {
    result.mDistance = closest_dist;
    result.mPoint1 = simplex[0].s1;
    result.mPoint2 = simplex[0].s2;
    vect< double, 3 > line_v = simplex[1].cso - simplex[0].cso;
    double i_proj = -( simplex[0].cso * line_v ) / ( line_v * line_v );
    result.mPoint1 += i_proj * ( simplex[1].s1 - result.mPoint1 );
    result.mPoint2 += i_proj * ( simplex[1].s2 - result.mPoint2 );
    return result;
  };
  if( ( closest_dist > 1e-5 ) && ( cur_size == 3 ) ) {
    result.mDistance = closest_dist;
    vect< double, 3 > b = proj_barycentric_coord( simplex[0].cso, simplex[1].cso, simplex[2].cso );
    result.mPoint1 = b[0] * simplex[0].s1 + b[1] * simplex[1].s1 + b[2] * simplex[2].s1;
    result.mPoint2 = b[0] * simplex[0].s2 + b[1] * simplex[1].s2 + b[2] * simplex[2].s2;
    return result;
  };
  //   if( (closest_dist > 1e-5) && (cur_size == 4) ) {
  //     std::cout << "ERROR THIS POINT SHOULD NOT HAPPEN!!!" << std::endl;
  //   };

  // Otherwise, do EPA to find the penetration depth, normal and deepest points:
  std::vector< GJK_simplex_point > epa_pts( simplex, simplex + cur_size );
  std::vector< EPA_edge_info > epa_edges;
  std::vector< EPA_face_info > epa_faces;


  // First, initialize the polytope with the simplex:
  if( cur_size == 4 ) {
    // check degeneracy of the simplex:
    vect< double, 3 > d1 = simplex[1].cso - simplex[0].cso;
    vect< double, 3 > d2 = simplex[2].cso - simplex[0].cso;
    vect< double, 3 > d3 = simplex[3].cso - simplex[0].cso;
    double volume = ( d1 % d2 ) * d3; // triple-product
    if( fabs( volume ) < 1e-5 ) {
      // find a point to eliminate:
      // TODO

      epa_pts.pop_back();
      cur_size = 3; // go to case 3
    } else {
      // Form the tetrahedron edges and faces:
      epa_edges.push_back( EPA_edge_info( epa_pts, 0, 1 ) );
      epa_edges.push_back( EPA_edge_info( epa_pts, 0, 2 ) );
      epa_edges.push_back( EPA_edge_info( epa_pts, 0, 3 ) );
      epa_edges.push_back( EPA_edge_info( epa_pts, 1, 2 ) );
      epa_edges.push_back( EPA_edge_info( epa_pts, 1, 3 ) );
      epa_edges.push_back( EPA_edge_info( epa_pts, 2, 3 ) );

      epa_faces.push_back( EPA_face_info( 0, epa_pts, epa_edges, 0, 1, 2 ) );
      epa_faces.push_back( EPA_face_info( 1, epa_pts, epa_edges, 0, 1, 3 ) );
      epa_faces.push_back( EPA_face_info( 2, epa_pts, epa_edges, 0, 2, 3 ) );
      epa_faces.push_back( EPA_face_info( 3, epa_pts, epa_edges, 1, 2, 3 ) );
    };
  };
  if( cur_size == 3 ) {
    // find the normal to the plane (in or out is not important):
    vect< double, 3 > line_n = ( simplex[1].cso - simplex[0].cso ) % ( simplex[2].cso - simplex[0].cso );
    double n_norm = norm_2( line_n );
    if( n_norm < 1e-5 ) {
      // degenerate facet.
      // find the redundant point:
      if( fabs( simplex[1].dist - simplex[0].dist ) < 0.1 * fabs( simplex[2].dist - simplex[0].dist ) ) {
        // v 1 is very close to v 0
        epa_pts[1] = epa_pts[2];
      };
      epa_pts.pop_back();
      cur_size = 2; // go to case 2
    } else {
      line_n /= n_norm;

      // add three points:
      epa_pts.resize( cur_size + 2 );
      epa_pts[3].init( line_n, aSFunc1, aSFunc2 );
      epa_pts[4].init( -line_n, aSFunc1, aSFunc2 );

      // Form the 'diamond' volume around the origin with edges and faces:
      epa_edges.push_back( EPA_edge_info( epa_pts, 3, 0 ) );
      epa_edges.push_back( EPA_edge_info( epa_pts, 3, 1 ) );
      epa_edges.push_back( EPA_edge_info( epa_pts, 3, 2 ) );
      epa_edges.push_back( EPA_edge_info( epa_pts, 4, 0 ) );
      epa_edges.push_back( EPA_edge_info( epa_pts, 4, 1 ) );
      epa_edges.push_back( EPA_edge_info( epa_pts, 4, 2 ) );
      epa_edges.push_back( EPA_edge_info( epa_pts, 0, 1 ) );
      epa_edges.push_back( EPA_edge_info( epa_pts, 0, 2 ) );
      epa_edges.push_back( EPA_edge_info( epa_pts, 1, 2 ) );

      epa_faces.push_back( EPA_face_info( 0, epa_pts, epa_edges, 0, 6, 1 ) );
      epa_faces.push_back( EPA_face_info( 1, epa_pts, epa_edges, 1, 8, 2 ) );
      epa_faces.push_back( EPA_face_info( 2, epa_pts, epa_edges, 0, 7, 2 ) );
      epa_faces.push_back( EPA_face_info( 3, epa_pts, epa_edges, 3, 6, 4 ) );
      epa_faces.push_back( EPA_face_info( 4, epa_pts, epa_edges, 4, 8, 5 ) );
      epa_faces.push_back( EPA_face_info( 5, epa_pts, epa_edges, 5, 7, 3 ) );
    };
  };
  if( cur_size == 2 ) {
    // find orthogonal basis normal to line segment:
    vect< double, 3 > line_v = simplex[0].cso - simplex[1].cso;
    vect< double, 3 > b1;
    if( fabs( line_v[0] ) > 1e-5 )
      b1 = line_v % vect< double, 3 >( 0.0, 1.0, 1.0 );
    else if( fabs( line_v[1] ) > 1e-5 )
      b1 = line_v % vect< double, 3 >( 1.0, 0.0, 1.0 );
    else
      b1 = line_v % vect< double, 3 >( 1.0, 1.0, 0.0 );
    vect< double, 3 > b2 = b1 % line_v;
    b1 /= norm_2( b1 );
    b2 /= norm_2( b2 );
    // add three points:
    epa_pts.resize( cur_size + 3 );
    epa_pts[2].init( b1, aSFunc1, aSFunc2 );
    epa_pts[3].init( -0.5 * b1 + 0.866 * b2, aSFunc1, aSFunc2 );
    epa_pts[4].init( -0.5 * b1 - 0.866 * b2, aSFunc1, aSFunc2 );

    // Form the 'diamond' volume around the origin with edges and faces:
    epa_edges.push_back( EPA_edge_info( epa_pts, 0, 2 ) );
    epa_edges.push_back( EPA_edge_info( epa_pts, 0, 3 ) );
    epa_edges.push_back( EPA_edge_info( epa_pts, 0, 4 ) );
    epa_edges.push_back( EPA_edge_info( epa_pts, 1, 2 ) );
    epa_edges.push_back( EPA_edge_info( epa_pts, 1, 3 ) );
    epa_edges.push_back( EPA_edge_info( epa_pts, 1, 4 ) );
    epa_edges.push_back( EPA_edge_info( epa_pts, 2, 3 ) );
    epa_edges.push_back( EPA_edge_info( epa_pts, 2, 4 ) );
    epa_edges.push_back( EPA_edge_info( epa_pts, 3, 4 ) );

    epa_faces.push_back( EPA_face_info( 0, epa_pts, epa_edges, 0, 6, 1 ) );
    epa_faces.push_back( EPA_face_info( 1, epa_pts, epa_edges, 1, 8, 2 ) );
    epa_faces.push_back( EPA_face_info( 2, epa_pts, epa_edges, 0, 7, 2 ) );
    epa_faces.push_back( EPA_face_info( 3, epa_pts, epa_edges, 3, 6, 4 ) );
    epa_faces.push_back( EPA_face_info( 4, epa_pts, epa_edges, 4, 8, 5 ) );
    epa_faces.push_back( EPA_face_info( 5, epa_pts, epa_edges, 5, 7, 3 ) );
  };

  // Then, initialize the priority-queue of faces to expand:
  std::priority_queue< EPA_face_entry > Q;
  std::vector< int > graveyard;
  std::vector< int > I;
  for( int i = 0; i < int( epa_faces.size() ); ++i )
    Q.push( EPA_face_entry( epa_faces[i].min_dist, i ) );
  double last_min_dist = std::numeric_limits< double >::infinity();
  int last_f_sol = -1;

  // Run the EPA iterations until you have no change in distance:
  while( Q.size()
         && !( ( last_f_sol == Q.top().f ) && ( fabs( last_min_dist - epa_faces[Q.top().f].min_dist ) < 1e-5 ) ) ) {

    int f_sol = Q.top().f;
    if( ( epa_faces[f_sol].color == EPA_black ) || ( epa_faces[f_sol].min_dist != Q.top().dist ) ) {
      Q.pop(); // skip any old inconsistent faces in the queue.
      continue;
    };
    last_min_dist = epa_faces[f_sol].min_dist;
    Q.pop();
    epa_faces[f_sol].color = EPA_grey; // being removed

#if 0
    for(std::size_t i = 0; i < epa_pts.size(); ++i) {
      std::cout << " - epa_pts[" << i << "].s1 = " << epa_pts[i].s1 
                << " - epa_pts[" << i << "].s2 = " << epa_pts[i].s2 
                << " - epa_pts[" << i << "].cso = " << epa_pts[i].cso << std::endl;
    };
    for(std::size_t i = 0; i < epa_edges.size(); ++i) {
      std::cout << " - edge(" << epa_edges[i].v[0] << " ," << epa_edges[i].v[1] 
                << ")->min_dist = " << epa_edges[i].min_dist << std::endl;
    };
    for(std::size_t i = 0; i < epa_faces.size(); ++i) {
      if( epa_faces[i].color == EPA_black )
        continue;
      std::cout << " - face(" << epa_faces[i].v[0] << " ," << epa_faces[i].v[1] << " ," << epa_faces[i].v[2] 
                << ")->min_dist = " << epa_faces[i].min_dist << std::endl;
    };
    std::cout << " - Top queue has face(" << epa_faces[f_sol].v[0] 
              << " ," << epa_faces[f_sol].v[1] 
              << " ," << epa_faces[f_sol].v[2] 
              << ")->min_dist = " << epa_faces[f_sol].min_dist << std::endl;
    std::cout << "   with normal = " << epa_faces[f_sol].min_pt << std::endl;
    vect<double,3> real_normal = (epa_pts[epa_faces[f_sol].v[1]].cso - epa_pts[epa_faces[f_sol].v[0]].cso)
                               % (epa_pts[epa_faces[f_sol].v[2]].cso - epa_pts[epa_faces[f_sol].v[0]].cso);
    std::cout << "   with real normal = " << real_normal << " or " << (-real_normal) << std::endl;
#endif

    I.push_back( epa_faces[f_sol].e[0] );              //
    epa_edges[epa_faces[f_sol].e[0]].color = EPA_grey; //
    I.push_back( epa_faces[f_sol].e[1] );              // edges must be re-evaluated
    epa_edges[epa_faces[f_sol].e[1]].color = EPA_grey; //
    I.push_back( epa_faces[f_sol].e[2] );              //
    epa_edges[epa_faces[f_sol].e[2]].color = EPA_grey; //

    // Add a support point in the normal direction of the min-face:
    epa_pts.resize( epa_pts.size() + 1 );
    epa_pts.back().init( epa_faces[f_sol].min_pt, aSFunc1, aSFunc2 );
    //     std::cout << "   giving new supports s1 = "
    //               << epa_pts.back().s1 << " s2 = "
    //               << epa_pts.back().s2 << " cso = " << epa_pts.back().cso << std::endl;

    // Add edges to the new support point:
    int e1 = static_cast<int>(epa_edges.size());
    epa_edges.push_back( EPA_edge_info( epa_pts, epa_faces[f_sol].v[0], static_cast<int>(epa_pts.size() - 1) ) );
    int e2 = e1 + 1;
    epa_edges.push_back(EPA_edge_info(epa_pts, epa_faces[f_sol].v[1], static_cast<int>(epa_pts.size() - 1)));
    int e3 = e2 + 1;
    epa_edges.push_back(EPA_edge_info(epa_pts, epa_faces[f_sol].v[2], static_cast<int>(epa_pts.size() - 1)));

    double best_new_dist = std::numeric_limits< double >::infinity();

    // Add faces:
    for( int i = 0; i < 3; ++i ) {
      if( epa_edges[epa_faces[f_sol].e[i]].f[0] == f_sol )
        epa_edges[epa_faces[f_sol].e[i]].f[0] = -1;
      else
        epa_edges[epa_faces[f_sol].e[i]].f[1] = -1;
    };
    int f1, f2, f3;
    f1 = EPA_add_face( epa_faces, graveyard, epa_pts, epa_edges, e1, e2, epa_faces[f_sol].e[0] );
    if( ( epa_edges[epa_faces[f_sol].e[1]].v[0] == epa_faces[f_sol].v[0] )
        || ( epa_edges[epa_faces[f_sol].e[1]].v[1] == epa_faces[f_sol].v[0] ) ) {
      // e-1 is connected to v-0:
      f2 = EPA_add_face( epa_faces, graveyard, epa_pts, epa_edges, e1, e3, epa_faces[f_sol].e[1] );
      f3 = EPA_add_face( epa_faces, graveyard, epa_pts, epa_edges, e2, e3, epa_faces[f_sol].e[2] );
    } else {
      f2 = EPA_add_face( epa_faces, graveyard, epa_pts, epa_edges, e1, e3, epa_faces[f_sol].e[2] );
      f3 = EPA_add_face( epa_faces, graveyard, epa_pts, epa_edges, e2, e3, epa_faces[f_sol].e[1] );
    };
    Q.push( EPA_face_entry( epa_faces[f1].min_dist, f1 ) );
    Q.push( EPA_face_entry( epa_faces[f2].min_dist, f2 ) );
    Q.push( EPA_face_entry( epa_faces[f3].min_dist, f3 ) );

    if( epa_faces[f1].min_dist < best_new_dist )
      best_new_dist = epa_faces[f1].min_dist;
    if( epa_faces[f2].min_dist < best_new_dist )
      best_new_dist = epa_faces[f2].min_dist;
    if( epa_faces[f3].min_dist < best_new_dist )
      best_new_dist = epa_faces[f3].min_dist;

    // Visit all possibly inconsistent edges to repair convexity of the polytope / hull:
    for( std::size_t i = 0; i < I.size(); ++i ) {
      if( epa_edges[I[i]].color == EPA_black )
        continue;
      int e_m = I[i];

      // take the two faces around the edge and check convexity:
      int f_1 = epa_edges[e_m].f[0];
      int f_2 = epa_edges[e_m].f[1];
      int v_A, v_B, v_m1, v_m2;
      v_m1 = epa_edges[e_m].v[0];
      v_m2 = epa_edges[e_m].v[1];
      if( ( epa_faces[f_1].v[0] != v_m1 ) && ( epa_faces[f_1].v[0] != v_m2 ) )
        v_A = epa_faces[f_1].v[0];
      else if( ( epa_faces[f_1].v[1] != v_m1 ) && ( epa_faces[f_1].v[1] != v_m2 ) )
        v_A = epa_faces[f_1].v[1];
      else
        v_A = epa_faces[f_1].v[2];
      if( ( epa_faces[f_2].v[0] != v_m1 ) && ( epa_faces[f_2].v[0] != v_m2 ) )
        v_B = epa_faces[f_2].v[0];
      else if( ( epa_faces[f_2].v[1] != v_m1 ) && ( epa_faces[f_2].v[1] != v_m2 ) )
        v_B = epa_faces[f_2].v[1];
      else
        v_B = epa_faces[f_2].v[2];

      if( true ) {
        //       vect<double,3> n = (epa_pts[v_m1].cso - epa_pts[v_A].cso) % (epa_pts[v_m2].cso - epa_pts[v_A].cso);
        //       double dA = epa_pts[v_A].cso * n;
        //       double dB = epa_pts[v_B].cso * n;
        //       if( fabs(dB) < 1.00001 * fabs(dA) ) {
        // If B is closer to the origin than A, then the edge is convex.
        // Close the edge and move on:
        epa_edges[e_m].color = EPA_black;
        continue;
      };

      // if concave, flip that middle edge to the other diagonal between the 4 points.
      //  and add the four exterior edges to the incons list, unless black already.
      int e_1A, e_2A, e_1B, e_2B;
      for( int j = 0; j < 3; ++j ) {
        if( epa_faces[f_1].e[j] != e_m ) {
          if( ( epa_edges[epa_faces[f_1].e[j]].v[0] == v_m1 ) || ( epa_edges[epa_faces[f_1].e[j]].v[1] == v_m1 ) )
            e_1A = epa_faces[f_1].e[j];
          else
            e_2A = epa_faces[f_1].e[j];
        };
        if( epa_faces[f_2].e[j] != e_m ) {
          if( ( epa_edges[epa_faces[f_2].e[j]].v[0] == v_m1 ) || ( epa_edges[epa_faces[f_2].e[j]].v[1] == v_m1 ) )
            e_1B = epa_faces[f_2].e[j];
          else
            e_2B = epa_faces[f_2].e[j];
        };
      };
      if( epa_edges[e_1A].f[0] == f_1 )
        epa_edges[e_1A].f[0] = -1;
      else
        epa_edges[e_1A].f[1] = -1;
      if( epa_edges[e_2A].f[0] == f_1 )
        epa_edges[e_2A].f[0] = -1;
      else
        epa_edges[e_2A].f[1] = -1;
      if( epa_edges[e_1B].f[0] == f_2 )
        epa_edges[e_1B].f[0] = -1;
      else
        epa_edges[e_1B].f[1] = -1;
      if( epa_edges[e_2B].f[0] == f_2 )
        epa_edges[e_2B].f[0] = -1;
      else
        epa_edges[e_2B].f[1] = -1;

      epa_edges[e_m] = EPA_edge_info( epa_pts, v_A, v_B );

      epa_faces[f_1] = EPA_face_info( f_1, epa_pts, epa_edges, e_1A, e_1B, e_m );
      epa_faces[f_2] = EPA_face_info( f_2, epa_pts, epa_edges, e_2A, e_2B, e_m );

      Q.push( EPA_face_entry( epa_faces[f_1].min_dist, f_1 ) );
      Q.push( EPA_face_entry( epa_faces[f_2].min_dist, f_2 ) );

      if( epa_faces[f_1].min_dist < best_new_dist )
        best_new_dist = epa_faces[f_1].min_dist;
      if( epa_faces[f_2].min_dist < best_new_dist )
        best_new_dist = epa_faces[f_2].min_dist;

      epa_edges[e_m].color = EPA_black;
      epa_edges[e_1A].color = EPA_grey;
      I.push_back( e_1A );
      epa_edges[e_2A].color = EPA_grey;
      I.push_back( e_2A );
      epa_edges[e_1B].color = EPA_grey;
      I.push_back( e_1B );
      epa_edges[e_2B].color = EPA_grey;
      I.push_back( e_2B );
    };

    // Re-color all the edges:
    for( std::size_t i = 0; i < I.size(); ++i )
      epa_edges[I[i]].color = EPA_white;
    I.clear();

    // Finally remove the min-face:
    epa_faces[f_sol].color = EPA_black; // removed
    graveyard.push_back( f_sol );

    if( best_new_dist < ( 1.0 + 1e-5 ) * last_min_dist ) {
      // if the best new face is not further way from the current face,
      // then that's the optimal (should be on top of queue).
      break;
    };

#if 0
    for(std::size_t i = 0; i < epa_faces.size(); ++i) {
      if( epa_faces[i].color == EPA_black )
        continue;
      std::cout << " - face(" << epa_faces[i].v[0] << " ," << epa_faces[i].v[1] << " ," << epa_faces[i].v[2] 
                << ")->min_dist = " << epa_faces[i].min_dist << std::endl;
    };
#endif
  };

  // Solution is on the top face:
  int f_sol = Q.top().f;

  // Fill in the colliding contact information and exit.
  result.mDistance = -epa_faces[f_sol].min_dist;
  result.mPoint1 = epa_pts[epa_faces[f_sol].v[0]].s1;
  result.mPoint2 = epa_pts[epa_faces[f_sol].v[0]].s2;
  vect< double, 3 > cur_cso = epa_pts[epa_faces[f_sol].v[0]].cso;
  for( int i = 1; i < 3; ++i ) {
    vect< double, 3 > line_v = epa_pts[epa_faces[f_sol].v[i]].cso - cur_cso;
    double i_proj = -( cur_cso * line_v ) / ( line_v * line_v );
    result.mPoint1 += i_proj * ( epa_pts[epa_faces[f_sol].v[i]].s1 - result.mPoint1 );
    result.mPoint2 += i_proj * ( epa_pts[epa_faces[f_sol].v[i]].s2 - result.mPoint2 );
    cur_cso += i_proj * line_v;
  };
  return result;
};
};
};

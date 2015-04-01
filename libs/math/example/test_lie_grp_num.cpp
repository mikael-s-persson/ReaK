
/*
 *    Copyright 2015 Sven Mikael Persson
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

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/lin_alg/complex_math.hpp>
#include <ReaK/math/kinetostatics/quat_alg.hpp>
#include <ReaK/math/kinetostatics/quat_num.hpp>

#include <iostream>
#include <fstream>
#include <cstdio>

using namespace ReaK;

struct const_ang_vel_sys {
  vect< double, 3 > ang_velocity;

  const_ang_vel_sys( const vect< double, 3 >& aAngVel ) : ang_velocity( aAngVel ){};

  vect< double, 3 > operator()( double, const unit_quat< double >& ) const { return ang_velocity; };

  vect< double, 4 > operator()( double, const vect< double, 4 >& q ) const {
    auto qw = quat< double >( q ) * ang_velocity;
    return vect< double, 4 >( 0.5 * qw[0], 0.5 * qw[1], 0.5 * qw[2], 0.5 * qw[3] );
  };
};


int main() {
  using namespace ReaK;

#if 0
  /* --------- interpolation tests ---------- */
  
  unit_quat< double > qx_05(0.877582562, 0.4794255386, 0.0, 0.0);
  unit_quat< double > qy_05(0.877582562, 0.0, 0.4794255386, 0.0);
  unit_quat< double > qz_05(0.877582562, 0.0, 0.0, 0.4794255386);

  std::vector< unit_quat< double > > ctrl_pts = {
    unit_quat< double >()
    , qx_05
//     , qx_05 * qy_05
//     , qz_05 * qy_05
    , qz_05 * qx_05 * qy_05
//     , qx_05 * qx_05 * qz_05 * qy_05
    , qz_05 * qx_05 * qx_05 * qz_05 * qy_05
  };
  std::vector< vect< double, 4 > > ctrl_pts_v4;
  for( auto q : ctrl_pts ) 
    ctrl_pts_v4.emplace_back( q[0], q[1], q[2], q[3] );
  
  // run interpolators on both kinds of ctrl-pts vectors:
  for(double eta = 0.0; eta < 1.005; eta += 0.01) {
//     unit_quat< double > q_mid  = lie_group::bezier( ctrl_pts.begin(), ctrl_pts.end(), eta );
//     vect< double, 4 > qv_mid   = lie_group::bezier( ctrl_pts_v4.begin(), ctrl_pts_v4.end(), eta );
    unit_quat< double > q_mid  = lie_group::catmull_rom( ctrl_pts.begin(), ctrl_pts.end(), eta, 0.0 );
    vect< double, 4 > qv_mid   = lie_group::catmull_rom( ctrl_pts_v4.begin(), ctrl_pts_v4.end(), eta, 0.0 );
    vect< double, 4 > qv_mid_u = unit(qv_mid);
    std::cout << " " << eta 
              << " " << q_mid[0] << " " << q_mid[1] << " " << q_mid[2] << " " << q_mid[3] 
              << " " << qv_mid[0] << " " << qv_mid[1] << " " << qv_mid[2] << " " << qv_mid[3] 
              << " " << qv_mid_u[0] << " " << qv_mid_u[1] << " " << qv_mid_u[2] << " " << qv_mid_u[3] << std::endl;
  };
#endif

  /* ----- numerical integration tests ------ */
  unit_quat< double > q0;
  vect< double, 4 > qv0( q0[0], q0[1], q0[2], q0[3] );
  const double t_start = 0.0;
  const double t_end = 100.0;
  const double dt_samp = 0.1;
  const double dt_int = 0.01;
  const_ang_vel_sys f_sys( vect< double, 3 >( 2, 4, 1 ) );

  for( double t = t_start; t < t_end + dt_samp * 0.5; t += dt_samp ) {
    //     q0  = lie_group::forward_euler( q0, f_sys, t, t + dt_samp, dt_int );
    //     qv0 = lie_group::forward_euler( qv0, f_sys, t, t + dt_samp, dt_int );
    q0 = lie_group::backward_euler( q0, f_sys, t, t + dt_samp, dt_int, 0.01 );
    qv0 = lie_group::backward_euler( qv0, f_sys, t, t + dt_samp, dt_int, 0.01 );
    vect< double, 4 > qv0_u = unit( qv0 );
    std::cout << " " << t << " " << q0[0] << " " << q0[1] << " " << q0[2] << " " << q0[3] << " " << qv0[0] << " "
              << qv0[1] << " " << qv0[2] << " " << qv0[3] << " " << qv0_u[0] << " " << qv0_u[1] << " " << qv0_u[2]
              << " " << qv0_u[3] << std::endl;
    qv0 = unit( qv0 ); // with normalization at sample steps.
  };

  return 0;
};


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

#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/math/lin_alg/mat_cholesky.hpp>

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


struct so3_state {
  unit_quat<double> q;
  vect<double, 3> w;
  explicit so3_state(const unit_quat<double>& aQ = (unit_quat<double>()), 
                     const vect<double,3>& aW = (vect<double,3>(0.0,0.0,0.0))) : q(aQ), w(aW) {};
};

struct so3_state_tangent {
  vect<double,3> dq;
  vect<double,3> dw;
  explicit so3_state_tangent(const vect<double,3>& aDQ = (vect<double,3>(0.0,0.0,0.0)), 
                             const vect<double,3>& aDW = (vect<double,3>(0.0,0.0,0.0))) : dq(aDQ), dw(aDW) {};
  
};

so3_state_tangent operator*( double s, const so3_state_tangent& x ) {
  return so3_state_tangent( s * x.dq, s * x.dw );
};

so3_state_tangent operator*( const so3_state_tangent& x, double s ) {
  return so3_state_tangent( s * x.dq, s * x.dw );
};

so3_state_tangent operator+( const so3_state_tangent& lhs, const so3_state_tangent& rhs ) {
  return so3_state_tangent( lhs.dq + rhs.dq, lhs.dw + rhs.dw );
};
so3_state_tangent operator-( const so3_state_tangent& lhs, const so3_state_tangent& rhs ) {
  return so3_state_tangent( lhs.dq - rhs.dq, lhs.dw - rhs.dw );
};
so3_state_tangent operator-( const so3_state_tangent& lhs ) {
  return so3_state_tangent( -lhs.dq, -lhs.dw );
};


so3_state oplus(const so3_state& lhs, const so3_state_tangent& rhs) {
  using lie_group::oplus;
  return so3_state( oplus(lhs.q, rhs.dq), lhs.w + rhs.dw );
};

so3_state_tangent ominus(const so3_state& lhs, const so3_state& rhs) {
  using lie_group::ominus;
  return so3_state_tangent( ominus(lhs.q, rhs.q), lhs.w - rhs.w );
};

so3_state_tangent otransport( const so3_state& a, const so3_state_tangent& v, const so3_state& b ) {
  unit_quat< double > ab = invert( b.q ) * a.q;
  return so3_state_tangent( ab.rotate( v.dq ), ab.rotate( v.dw ) );
};

so3_state_tangent ocross(const so3_state_tangent& lhs, const so3_state_tangent& rhs) {
  return so3_state_tangent(lhs.dq % rhs.dq);
};


struct const_ang_mom_sys {
  mat< double, mat_structure::symmetric > I;
  mat< double, mat_structure::symmetric > Iinv;
  
  const_ang_mom_sys( const mat< double, mat_structure::symmetric >& aInertiaMoment )
      : I( aInertiaMoment ) {
    invert_Cholesky( I, Iinv );
  };
  
  so3_state_tangent operator()( double, const so3_state& x ) const {
    vect< double, 3 > w = x.w;
    vect< double, 3 > aacc = Iinv * ( ( I * w ) % w );
    
    return so3_state_tangent( w, aacc );
//     vect< double, 3 > v = get_velocity( x );
//     return point_derivative_type(
//       make_arithmetic_tuple( v, vect< double, 3 >( 0.0, 0.0, 9.81 ) ),
//       make_arithmetic_tuple( w, aacc ) );
  };
  
  vect< double, 7 > operator()( double, const vect< double, 7 >& x ) const {
    quat< double > q( vect< double, 4 >( x[0], x[1], x[2], x[3] ) );
    vect< double, 3 > w( x[4], x[5], x[6] );
    auto qw   = q * w;
    auto aacc = Iinv * ( ( I * w ) % w );
    
    return vect< double, 7 >( 0.5 * qw[0], 0.5 * qw[1], 0.5 * qw[2], 0.5 * qw[3], aacc[0], aacc[1], aacc[2] );
  };
};



struct compare_quat_lex_less {
  template < typename Vector >
  bool operator()( const Vector& u, const Vector& v, int i = 0 ) const {
    if( u[i] < v[i] )
      return true;
    else if( u[i] > v[i] )
      return false;
    else {
      if( i == 3 )
        return false;
      return ( *this )( u, v, i + 1 ); // tail-call
    };
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

#if 0
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
//     q0 = lie_group::backward_euler( q0, f_sys, t, t + dt_samp, dt_int, 0.01 );
//     qv0 = lie_group::backward_euler( qv0, f_sys, t, t + dt_samp, dt_int, 0.01 );
//     q0 = lie_group::trapezoidal_rule( q0, f_sys, t, t + dt_samp, dt_int, 0.01 );
//     qv0 = lie_group::trapezoidal_rule( qv0, f_sys, t, t + dt_samp, dt_int, 0.01 );
//     q0 = lie_group::runge_kutta_4( q0, f_sys, t, t + dt_samp, dt_int );
//     qv0 = lie_group::runge_kutta_4( qv0, f_sys, t, t + dt_samp, dt_int );
    q0 = lie_group::corrected_runge_kutta_4( q0, f_sys, t, t + dt_samp, dt_int );
//     qv0 = lie_group::corrected_runge_kutta_4( qv0, f_sys, t, t + dt_samp, dt_int );
    vect< double, 4 > qv0_u = unit( qv0 );
    std::cout << " " << t << " " << q0[0] << " " << q0[1] << " " << q0[2] << " " << q0[3] << " " << qv0[0] << " "
              << qv0[1] << " " << qv0[2] << " " << qv0[3] << " " << qv0_u[0] << " " << qv0_u[1] << " " << qv0_u[2]
              << " " << qv0_u[3] << std::endl;
    qv0 = unit( qv0 ); // with normalization at sample steps.
  };
#endif

#if 0
  /* --------- interpolation tests ---------- */

  unit_quat< double > qx_05( 0.877582562, 0.4794255386, 0.0, 0.0 );
  unit_quat< double > qy_05( 0.877582562, 0.0, 0.4794255386, 0.0 );
  unit_quat< double > qz_05( 0.877582562, 0.0, 0.0, 0.4794255386 );

  std::vector< unit_quat< double > > avg_pts
    = {unit_quat< double >(), qx_05,                         qx_05 * qy_05,                        qz_05 * qy_05,
       qz_05 * qx_05 * qy_05, qx_05 * qx_05 * qz_05 * qy_05, qz_05 * qx_05 * qx_05 * qz_05 * qy_05};
  std::vector< vect< double, 4 > > avg_pts_v4;
  for( auto q : avg_pts )
    avg_pts_v4.emplace_back( q[0], q[1], q[2], q[3] );

  std::sort( avg_pts.begin(), avg_pts.end(), compare_quat_lex_less() );
  std::sort( avg_pts_v4.begin(), avg_pts_v4.end(), compare_quat_lex_less() );

  // run averaging on both kinds of ctrl-pts vectors:
  do {
    unit_quat< double > q_avg = lie_group::average( avg_pts.begin(), avg_pts.end() );
    vect< double, 4 > qv_avg = lie_group::average( avg_pts_v4.begin(), avg_pts_v4.end() );
    vect< double, 4 > qv_avg_u = unit( qv_avg );
    std::cout << " " << q_avg[0] << " " << q_avg[1] << " " << q_avg[2] << " " << q_avg[3] << " " << qv_avg[0] << " "
              << qv_avg[1] << " " << qv_avg[2] << " " << qv_avg[3] << " " << qv_avg_u[0] << " " << qv_avg_u[1] << " "
              << qv_avg_u[2] << " " << qv_avg_u[3] << std::endl;
    std::next_permutation( avg_pts.begin(), avg_pts.end(), compare_quat_lex_less() );
  } while( std::next_permutation( avg_pts_v4.begin(), avg_pts_v4.end(), compare_quat_lex_less() ) );
#endif

#if 1
  /* ----- numerical integration tests on spinning rigid body ------ */
  unit_quat< double > q0;
  vect< double, 3 > w0(1, 2, 0.5);
  so3_state s0( q0, w0 );
  vect< double, 7 > sv0( q0[0], q0[1], q0[2], q0[3], w0[0], w0[1], w0[2] );
  const double t_start = 0.0;
  const double t_end = 100.0;
  const double dt_samp = 0.1;
  const double dt_int = 0.01;
  const_ang_mom_sys f_sys( mat< double, mat_structure::symmetric >( 2.0, 0.5, -0.25, 3.0, 1.25, 4.5 ) );
  
  for( double t = t_start; t < t_end + dt_samp * 0.5; t += dt_samp ) {
//     s0  = lie_group::forward_euler( s0, f_sys, t, t + dt_samp, dt_int );
//     sv0 = lie_group::forward_euler( sv0, f_sys, t, t + dt_samp, dt_int );
//     s0 = lie_group::backward_euler( s0, f_sys, t, t + dt_samp, dt_int, 0.01 );
//     sv0 = lie_group::backward_euler( sv0, f_sys, t, t + dt_samp, dt_int, 0.01 );
//     s0 = lie_group::trapezoidal_rule( s0, f_sys, t, t + dt_samp, dt_int, 0.01 );
//     sv0 = lie_group::trapezoidal_rule( sv0, f_sys, t, t + dt_samp, dt_int, 0.01 );
//     s0 = lie_group::runge_kutta_4( s0, f_sys, t, t + dt_samp, dt_int );
//     sv0 = lie_group::runge_kutta_4( sv0, f_sys, t, t + dt_samp, dt_int );
    s0 = lie_group::corrected_runge_kutta_4( s0, f_sys, t, t + dt_samp, dt_int );
//     sv0 = lie_group::corrected_runge_kutta_4( sv0, f_sys, t, t + dt_samp, dt_int );
    vect< double, 4 > qv0_u = unit( vect< double, 4 >( sv0[0], sv0[1], sv0[2], sv0[3] ) );
    vect< double, 7 > sv0_u( qv0_u[0], qv0_u[1], qv0_u[2], qv0_u[3], sv0[0], sv0[1], sv0[2] );
    std::cout << " " << t << " " 
              << s0.q[0] << " " << s0.q[1] << " " << s0.q[2] << " " << s0.q[3] << " "
              << s0.w[0] << " " << s0.w[1] << " " << s0.w[2] << " " 
              << sv0[0] << " " << sv0[1] << " " << sv0[2] << " " << sv0[3] << " " 
              << sv0[4] << " " << sv0[5] << " " << sv0[6] << " " 
              << qv0_u[0] << " " << qv0_u[1] << " " << qv0_u[2] << " " << qv0_u[3] << std::endl;
    sv0 = sv0_u; // with normalization at sample steps.
  };
#endif

  return 0;
};

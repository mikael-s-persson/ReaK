
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
#include <ReaK/math/kinetostatics/quat_alg.hpp>
#include <ReaK/math/kinetostatics/quat_num.hpp>
#include <ReaK/math/lin_alg/complex_math.hpp>

#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/math/lin_alg/mat_cholesky.hpp>

#include <cstdio>
#include <fstream>
#include <iostream>

using namespace ReaK;

struct const_ang_vel_sys {
  vect<double, 3> ang_velocity;

  explicit const_ang_vel_sys(const vect<double, 3>& aAngVel)
      : ang_velocity(aAngVel){}

  vect<double, 3> operator()(double /*unused*/,
                             const unit_quat<double>& /*unused*/) const {
    return ang_velocity;
  }

  vect<double, 4> operator()(double /*unused*/,
                             const vect<double, 4>& q) const {
    auto qw = quat<double>(q) * ang_velocity;
    return {0.5 * qw[0], 0.5 * qw[1], 0.5 * qw[2], 0.5 * qw[3]};
  }
};

struct so3_state {
  unit_quat<double> q;
  vect<double, 3> w;
  explicit so3_state(const unit_quat<double>& aQ = (unit_quat<double>()),
                     const vect<double, 3>& aW = (vect<double, 3>(0.0, 0.0,
                                                                  0.0)))
      : q(aQ), w(aW){}
};

struct so3_state_tangent {
  vect<double, 3> dq;
  vect<double, 3> dw;
  explicit so3_state_tangent(
      const vect<double, 3>& aDQ = (vect<double, 3>(0.0, 0.0, 0.0)),
      const vect<double, 3>& aDW = (vect<double, 3>(0.0, 0.0, 0.0)))
      : dq(aDQ), dw(aDW){}
};

so3_state_tangent operator*(double s, const so3_state_tangent& x) {
  return so3_state_tangent(s * x.dq, s * x.dw);
}

so3_state_tangent operator*(const so3_state_tangent& x, double s) {
  return so3_state_tangent(s * x.dq, s * x.dw);
}

so3_state_tangent operator+(const so3_state_tangent& lhs,
                            const so3_state_tangent& rhs) {
  return so3_state_tangent(lhs.dq + rhs.dq, lhs.dw + rhs.dw);
}
so3_state_tangent operator-(const so3_state_tangent& lhs,
                            const so3_state_tangent& rhs) {
  return so3_state_tangent(lhs.dq - rhs.dq, lhs.dw - rhs.dw);
}
so3_state_tangent operator-(const so3_state_tangent& lhs) {
  return so3_state_tangent(-lhs.dq, -lhs.dw);
}

so3_state oplus(const so3_state& lhs, const so3_state_tangent& rhs) {
  using lie_group::oplus;
  return so3_state(oplus(lhs.q, rhs.dq), lhs.w + rhs.dw);
}

so3_state_tangent ominus(const so3_state& lhs, const so3_state& rhs) {
  using lie_group::ominus;
  return so3_state_tangent(ominus(lhs.q, rhs.q), lhs.w - rhs.w);
}

so3_state_tangent otransport(const so3_state& a, const so3_state_tangent& v,
                             const so3_state& b) {
  unit_quat<double> ab = invert(b.q) * a.q;
  return so3_state_tangent(ab.rotate(v.dq), ab.rotate(v.dw));
}

so3_state_tangent ocross(const so3_state_tangent& lhs,
                         const so3_state_tangent& rhs) {
  return so3_state_tangent(lhs.dq % rhs.dq);
}

double norm_2(const so3_state_tangent& x) {
  return norm_2(x.dq) + norm_2(x.dw);
}

struct const_ang_mom_sys {
  mat<double, mat_structure::symmetric> I;
  mat<double, mat_structure::symmetric> Iinv;

  explicit const_ang_mom_sys(
      const mat<double, mat_structure::symmetric>& aInertiaMoment)
      : I(aInertiaMoment) {
    invert_Cholesky(I, Iinv);
  }

  so3_state_tangent operator()(double /*unused*/, const so3_state& x) const {
    vect<double, 3> w = x.w;
    vect<double, 3> aacc = Iinv * ((I * w) % w);

    return so3_state_tangent(w, aacc);
  }

  vect<double, 7> operator()(double /*unused*/,
                             const vect<double, 7>& x) const {
    quat<double> q(vect<double, 4>(x[0], x[1], x[2], x[3]));
    vect<double, 3> w(x[4], x[5], x[6]);
    auto qw = q * w;
    auto aacc = Iinv * ((I * w) % w);

    return {0.5 * qw[0], 0.5 * qw[1], 0.5 * qw[2], 0.5 * qw[3],
                           aacc[0], aacc[1], aacc[2]};
  }
};

so3_state trapm_rule(so3_state q, const_ang_mom_sys f, double t_start,
                     double t_end, double dt, double tol) {
  using lie_group::ominus;
  using lie_group::oplus;
  using lie_group::otransport;
  assert(t_end > t_start);
  double t = t_start;
  while (t < t_end) {
    if (t + dt > t_end) {
      dt = t_end - t;
    }
    unit_quat<double> q0 = q.q;
    vect<double, 3> w0 = q.w;
    unit_quat<double> q0_mid = oplus(q0, (0.5 * dt) * w0);
    q.q = oplus(q0, dt * w0);
    while (true) {
      q.w = f.Iinv * otransport(q0, f.I * w0, q.q);
      unit_quat<double> q1_mid = oplus(q.q, (-0.5 * dt) * q.w);
      auto q_err = ominus(q0_mid, q1_mid);
      if (norm_2(q_err) < tol) {
        break;
      }
      q.q = oplus(q0_mid, otransport(q.q, (0.5 * dt) * q.w, q0_mid));
    }
    t += dt;
  }
  return q;
}

struct compare_quat_lex_less {
  template <typename Vector>
  bool operator()(const Vector& u, const Vector& v, int i = 0) const {
    if (u[i] < v[i]) {
      return true;
    }
    if (u[i] > v[i]) {
      return false;
    }
    if (i == 3) {
      return false;
    }
    return (*this)(u, v, i + 1);  // tail-call
  }
};

inline void save_state_vectors(std::ostream& out, double t, const so3_state& s,
                               const vect<double, 7>& sv) {
  out << " " << t << " " << s.q[0] << " " << s.q[1] << " " << s.q[2] << " "
      << s.q[3] << " " << s.w[0] << " " << s.w[1] << " " << s.w[2] << " "
      << sv[0] << " " << sv[1] << " " << sv[2] << " " << sv[3] << " " << sv[4]
      << " " << sv[5] << " " << sv[6] << std::endl;
}

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
  }
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
  }
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
  }
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
  unit_quat<double> q0;
  vect<double, 3> w0(1, 2, 0.5);
  so3_state s0(q0, w0);
  vect<double, 7> sv0(q0[0], q0[1], q0[2], q0[3], w0[0], w0[1], w0[2]);
  so3_state s0_fe = s0;
  so3_state s0_be = s0;
  so3_state s0_tr = s0;
  so3_state s0_rk = s0;
  so3_state s0_crk = s0;
  so3_state s0_trm = s0;
  vect<double, 7> sv0_fe = sv0;
  vect<double, 7> sv0_be = sv0;
  vect<double, 7> sv0_tr = sv0;
  vect<double, 7> sv0_rk = sv0;
  const double t_start = 0.0;
  const double t_end = 100.0;
  const double dt_samp = 1;
  const double dt_int = 0.1;
  const_ang_mom_sys f_sys(
      mat<double, mat_structure::symmetric>(2.0, 0.5, -0.25, 3.0, 1.25, 4.5));

  std::ofstream f_out_fe("test_lie_grp_num_result/cm_fe.ssv");
  std::ofstream f_out_be("test_lie_grp_num_result/cm_be.ssv");
  std::ofstream f_out_tr("test_lie_grp_num_result/cm_tr.ssv");
  std::ofstream f_out_rk("test_lie_grp_num_result/cm_rk.ssv");
  std::ofstream f_out_crk("test_lie_grp_num_result/cm_crk.ssv");
  std::ofstream f_out_trm("test_lie_grp_num_result/cm_trm.ssv");

  for (double t = t_start; t < t_end + dt_samp * 0.5; t += dt_samp) {
    s0_fe = lie_group::forward_euler(s0_fe, f_sys, t, t + dt_samp, dt_int);
    sv0_fe = lie_group::forward_euler(sv0_fe, f_sys, t, t + dt_samp, dt_int);
    s0_be =
        lie_group::backward_euler(s0_be, f_sys, t, t + dt_samp, dt_int, 0.0001);
    sv0_be = lie_group::backward_euler(sv0_be, f_sys, t, t + dt_samp, dt_int,
                                       0.0001);
    s0_tr = lie_group::trapezoidal_rule(s0_tr, f_sys, t, t + dt_samp, dt_int,
                                        0.0001);
    sv0_tr = lie_group::trapezoidal_rule(sv0_tr, f_sys, t, t + dt_samp, dt_int,
                                         0.0001);
    s0_rk = lie_group::runge_kutta_4(s0_rk, f_sys, t, t + dt_samp, dt_int);
    sv0_rk = lie_group::runge_kutta_4(sv0_rk, f_sys, t, t + dt_samp, dt_int);
    s0_crk = lie_group::corrected_runge_kutta_4(s0_crk, f_sys, t, t + dt_samp,
                                                dt_int);
    s0_trm = trapm_rule(s0_trm, f_sys, t, t + dt_samp, dt_int, 0.0001);
    vect<double, 4> qv0_fe_u =
        unit(vect<double, 4>(sv0_fe[0], sv0_fe[1], sv0_fe[2], sv0_fe[3]));
    vect<double, 7> sv0_fe_u(qv0_fe_u[0], qv0_fe_u[1], qv0_fe_u[2], qv0_fe_u[3],
                             sv0_fe[4], sv0_fe[5], sv0_fe[6]);
    vect<double, 4> qv0_be_u =
        unit(vect<double, 4>(sv0_be[0], sv0_be[1], sv0_be[2], sv0_be[3]));
    vect<double, 7> sv0_be_u(qv0_be_u[0], qv0_be_u[1], qv0_be_u[2], qv0_be_u[3],
                             sv0_be[4], sv0_be[5], sv0_be[6]);
    vect<double, 4> qv0_tr_u =
        unit(vect<double, 4>(sv0_tr[0], sv0_tr[1], sv0_tr[2], sv0_tr[3]));
    vect<double, 7> sv0_tr_u(qv0_tr_u[0], qv0_tr_u[1], qv0_tr_u[2], qv0_tr_u[3],
                             sv0_tr[4], sv0_tr[5], sv0_tr[6]);
    vect<double, 4> qv0_rk_u =
        unit(vect<double, 4>(sv0_rk[0], sv0_rk[1], sv0_rk[2], sv0_rk[3]));
    vect<double, 7> sv0_rk_u(qv0_rk_u[0], qv0_rk_u[1], qv0_rk_u[2], qv0_rk_u[3],
                             sv0_rk[4], sv0_rk[5], sv0_rk[6]);
    save_state_vectors(f_out_fe, t, s0_fe, sv0_fe_u);
    save_state_vectors(f_out_be, t, s0_be, sv0_be_u);
    save_state_vectors(f_out_tr, t, s0_tr, sv0_tr_u);
    save_state_vectors(f_out_rk, t, s0_rk, sv0_rk_u);
    save_state_vectors(f_out_crk, t, s0_crk, sv0_rk_u);
    save_state_vectors(f_out_trm, t, s0_trm, sv0_tr_u);
    sv0_fe = sv0_fe_u;  // with normalization at sample steps.
    sv0_be = sv0_be_u;  // with normalization at sample steps.
    sv0_tr = sv0_tr_u;  // with normalization at sample steps.
    sv0_rk = sv0_rk_u;  // with normalization at sample steps.
  }
#endif

  return 0;
}

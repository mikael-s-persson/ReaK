
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

#include <cmath>

#include <ReaK/topologies/interpolation/sustained_acceleration_pulse_Ndof_detail.hpp>

#include <ReaK/math/root_finders/bisection_method.hpp>
#include <ReaK/math/root_finders/secant_method.hpp>

#include <ReaK/core/sorting/insertion_sort.hpp>

#include <ReaK/math/optimization/optim_exceptions.hpp>

namespace ReaK::pp::detail {

static void sap_Ndof_compute_interpolated_values_balanced(
    double start_position, double end_position, double start_velocity,
    double end_velocity, double peak_velocity, double max_velocity,
    double max_acceleration, double dt, double dt_total, double& result_pos,
    double& result_vel, double& result_acc, double& result_jerk) {

  using std::sqrt;

  double dt_vp1_1st = peak_velocity - start_velocity;
  double sgn_vp1 = 1.0;
  if (peak_velocity < start_velocity) {
    sgn_vp1 = -1.0;
    dt_vp1_1st = -dt_vp1_1st;
  }
  // we know that dt_vp_2nd = dt_vp_1st + dt_amax
  double dt_vp1 = dt_vp1_1st - max_acceleration;
  double dt_ap1 = max_acceleration;
  if (dt_vp1 < 0.0) {
    // means that we don't have time to reach the maximum acceleration:
    dt_vp1 = 0.0;
    dt_ap1 = sqrt(max_acceleration * dt_vp1_1st);
  }
  double dp_ramp1 = 0.5 * (dt_vp1 + 2.0 * dt_ap1) *
                    (start_velocity + peak_velocity) / max_velocity;

  double dt_vp2_1st = end_velocity - peak_velocity;
  double sgn_vp2 = 1.0;
  if (end_velocity < peak_velocity) {
    sgn_vp2 = -1.0;
    dt_vp2_1st = -dt_vp2_1st;
  }
  // we know that dt_vp_2nd = dt_vp_1st + dt_amax
  double dt_vp2 = dt_vp2_1st - max_acceleration;
  double dt_ap2 = max_acceleration;
  if (dt_vp2 < 0.0) {
    // means that we don't have time to reach the maximum acceleration:
    dt_vp2 = 0.0;
    dt_ap2 = sqrt(max_acceleration * dt_vp2_1st);
  }
  double dp_ramp2 = 0.5 * (dt_vp2 + 2.0 * dt_ap2) *
                    (end_velocity + peak_velocity) / max_velocity;

  dt_total -= dt_vp2 + 2.0 * dt_ap2 + dt_vp1 + 2.0 * dt_ap1;

  if (dt <= 0.0) {

    result_pos = start_position;
    result_vel = start_velocity;
    result_acc = 0.0;
    result_jerk = 0.0;

  } else if (dt < (2.0 * dt_ap1 + dt_vp1)) {

    if (dt < dt_ap1) {

      // Segment 1: in the jerk-up phase of velocity ramp-up.
      result_pos = start_position + (start_velocity + sgn_vp1 / 6.0 * dt * dt /
                                                          max_acceleration) *
                                        dt / max_velocity;
      result_vel = start_velocity + 0.5 * sgn_vp1 * dt * dt / max_acceleration;
      result_acc = dt * sgn_vp1;
      result_jerk = sgn_vp1;

    } else if (dt < dt_ap1 + dt_vp1) {
      // Segment 2: in the constant accel phase of velocity ramp-up.
      dt -= dt_ap1;

      result_pos =
          start_position +
          ((start_velocity + 0.5 * dt * sgn_vp1) * (max_acceleration + dt) +
           max_acceleration * max_acceleration * sgn_vp1 / 6.0) /
              max_velocity;
      result_vel = start_velocity + (dt + 0.5 * max_acceleration) * sgn_vp1;
      result_acc = max_acceleration * sgn_vp1;
      result_jerk = 0.0;

    } else {
      // Segment 3: in the jerk-down phase of velocity ramp-up.
      double mdt = dt - 2.0 * dt_ap1 - dt_vp1;

      result_pos =
          start_position + dp_ramp1 +
          (peak_velocity - sgn_vp1 / 6.0 * mdt * mdt / max_acceleration) * mdt /
              max_velocity;
      result_vel = peak_velocity - 0.5 * sgn_vp1 * mdt * mdt / max_acceleration;
      result_acc = mdt * sgn_vp1;
      result_jerk = -sgn_vp1;
    }

  } else if (dt < (2.0 * dt_ap1 + dt_vp1 + dt_total)) {

    // Segment 4: in the cruise phase.
    double idt = dt - 2.0 * dt_ap1 - dt_vp1;
    double pis = start_position + dp_ramp1;
    double pie = end_position - dp_ramp2;
    result_pos = idt / dt_total * pie + (1.0 - idt / dt_total) * pis;
    result_vel = peak_velocity;
    result_acc = 0.0;
    result_jerk = 0.0;

  } else if (dt < (2.0 * dt_ap1 + dt_vp1 + dt_total + 2.0 * dt_ap2 + dt_vp2)) {

    if (dt < (2.0 * dt_ap1 + dt_vp1 + dt_total + dt_ap2)) {

      // Segment 5: in the jerk-up phase of velocity ramp-down.
      dt -= 2.0 * dt_ap1 + dt_vp1 + dt_total;

      result_pos =
          end_position - dp_ramp2 +
          (peak_velocity + sgn_vp2 / 6.0 * dt * dt / max_acceleration) * dt /
              max_velocity;
      result_vel = peak_velocity + 0.5 * sgn_vp2 * dt * dt / max_acceleration;
      result_acc = dt * sgn_vp2;
      result_jerk = sgn_vp2;

    } else if (dt < (2.0 * dt_ap1 + dt_vp1 + dt_total + dt_ap2 + dt_vp2)) {

      // Segment 6: in the constant accel phase of velocity ramp-down.
      dt -= 2.0 * dt_ap1 + dt_vp1 + dt_total + dt_ap2;

      result_pos = end_position - dp_ramp2 +
                   (peak_velocity * (max_acceleration + dt) +
                    (max_acceleration * max_acceleration / 6.0 +
                     0.5 * max_acceleration * dt + 0.5 * dt * dt) *
                        sgn_vp2) /
                       max_velocity;
      result_vel = peak_velocity + (dt + 0.5 * max_acceleration) * sgn_vp2;
      result_acc = max_acceleration * sgn_vp2;
      result_jerk = 0.0;

    } else {

      // Segment 7: in the jerk-down phase of velocity ramp-down.
      double mdt =
          dt - 2.0 * dt_ap1 - dt_vp1 - dt_total - 2.0 * dt_ap2 - dt_vp2;

      result_pos = end_position + (end_velocity - sgn_vp2 / 6.0 * mdt * mdt /
                                                      max_acceleration) *
                                      mdt / max_velocity;
      result_vel = end_velocity - 0.5 * sgn_vp2 * mdt * mdt / max_acceleration;
      result_acc = mdt * sgn_vp2;
      result_jerk = -sgn_vp2;
    }

  } else {

    result_pos = end_position;
    result_vel = end_velocity;
    result_acc = 0.0;
    result_jerk = 0.0;
  }
}

void sap_Ndof_compute_interpolated_values(
    double start_position, double end_position, double start_velocity,
    double end_velocity, double peak_velocity, double max_velocity,
    double max_acceleration, double dt, double dt_total, double& result_pos,
    double& result_vel, double& result_acc, double& result_jerk) {

  sap_Ndof_compute_interpolated_values_balanced(
      start_position, end_position, start_velocity, end_velocity, peak_velocity,
      max_velocity, max_acceleration, dt, dt_total, result_pos, result_vel,
      result_acc, result_jerk);
}

inline void sap_Ndof_compute_ramp_dist_and_time(double v1, double v2,
                                                double vmax, double amax,
                                                double& d_pos, double& dt) {

  using std::abs;
  using std::sqrt;

  if (abs(v2 - v1) >= amax) {
    dt = abs(v2 - v1) + amax;
    d_pos = 0.5 * dt * (v1 + v2) / vmax;
  } else {
    dt = 2.0 * sqrt(amax * abs(v2 - v1));
    d_pos = 0.5 * dt * (v1 + v2) / vmax;
  }
}

struct sap_Ndof_no_cruise_calculator {
  double dp, v1, v2, vmax, amax;

  sap_Ndof_no_cruise_calculator(double a_dp, double a_v1, double a_v2,
                                double a_vmax, double a_amax)
      : dp(a_dp), v1(a_v1), v2(a_v2), vmax(a_vmax), amax(a_amax) {}

  double operator()(double vp) const {
    double dp1 = NAN;
    double dt1 = NAN;
    double dp2 = NAN;
    double dt2 = NAN;
    sap_Ndof_compute_ramp_dist_and_time(v1, vp, vmax, amax, dp1, dt1);
    sap_Ndof_compute_ramp_dist_and_time(vp, v2, vmax, amax, dp2, dt2);

    if (dp < 0.0) {
      return dp1 + dp2 - dp;
    }
    return dp - dp1 - dp2;
  }
};

static double sap_Ndof_compute_min_delta_time_numsolve(
    double start_position, double end_position, double start_velocity,
    double end_velocity, double& peak_velocity, double max_velocity,
    double max_acceleration) {
  using std::abs;
  using std::sqrt;

  if ((abs(end_position - start_position) < 1e-6 * max_velocity) &&
      (abs(end_velocity - start_velocity) < 1e-6 * max_acceleration)) {
    peak_velocity = start_velocity;
    return 0.0;
  }

  if ((abs(start_velocity) > max_velocity) ||
      (abs(end_velocity) > max_velocity)) {
    peak_velocity = 0.0;
    return std::numeric_limits<double>::infinity();
  }

  double sign_p1_p0 = 1.0;
  if (start_position > end_position) {
    sign_p1_p0 = -1.0;
  }

  double dt_ramp1 = 0.0;
  double dp_ramp1 = 0.0;
  double dt_ramp2 = 0.0;
  double dp_ramp2 = 0.0;

  sap_Ndof_no_cruise_calculator nc_calc(end_position - start_position,
                                        start_velocity, end_velocity,
                                        max_velocity, max_acceleration);
  double peak_vel_low = -sign_p1_p0 * max_velocity;
  double peak_vel_hi = sign_p1_p0 * max_velocity;

  double delta_first_order = nc_calc(peak_vel_hi);
  if (delta_first_order > 0.0) {
    peak_velocity = peak_vel_hi;
    sap_Ndof_compute_ramp_dist_and_time(start_velocity, peak_velocity,
                                        max_velocity, max_acceleration,
                                        dp_ramp1, dt_ramp1);
    sap_Ndof_compute_ramp_dist_and_time(peak_velocity, end_velocity,
                                        max_velocity, max_acceleration,
                                        dp_ramp2, dt_ramp2);
    return delta_first_order + dt_ramp1 + dt_ramp2;
  }

  delta_first_order = nc_calc(peak_vel_low);
  if (delta_first_order < 0.0) {
    peak_velocity = peak_vel_low;
    sap_Ndof_compute_ramp_dist_and_time(start_velocity, peak_velocity,
                                        max_velocity, max_acceleration,
                                        dp_ramp1, dt_ramp1);
    sap_Ndof_compute_ramp_dist_and_time(peak_velocity, end_velocity,
                                        max_velocity, max_acceleration,
                                        dp_ramp2, dt_ramp2);
    return -delta_first_order + dt_ramp1 + dt_ramp2;
  }

  brent_method(peak_vel_low, peak_vel_hi, nc_calc, 1e-8);

  peak_velocity = peak_vel_hi;
  sap_Ndof_compute_ramp_dist_and_time(start_velocity, peak_velocity,
                                      max_velocity, max_acceleration, dp_ramp1,
                                      dt_ramp1);
  sap_Ndof_compute_ramp_dist_and_time(peak_velocity, end_velocity, max_velocity,
                                      max_acceleration, dp_ramp2, dt_ramp2);
  return abs(end_position - start_position - dp_ramp1 - dp_ramp2) + dt_ramp1 +
         dt_ramp2;
}

double sap_Ndof_compute_min_delta_time(
    double start_position, double end_position, double start_velocity,
    double end_velocity, double& peak_velocity, double max_velocity,
    double max_acceleration) {

  double min_delta_time = sap_Ndof_compute_min_delta_time_numsolve(
      start_position, end_position, start_velocity, end_velocity, peak_velocity,
      max_velocity, max_acceleration);

  return min_delta_time;
}

struct sap_Ndof_pos_diff_calculator {
  double dp;
  double v1;
  double v2;
  double vmax;
  double amax;
  double dt;

  sap_Ndof_pos_diff_calculator(double a_dp, double a_v1, double a_v2,
                               double a_vmax, double a_amax, double a_dt)
      : dp(a_dp), v1(a_v1), v2(a_v2), vmax(a_vmax), amax(a_amax), dt(a_dt) {}

  double operator()(double vp) const {
    double dp1 = 0.0;
    double dt1 = 0.0;
    double dp2 = 0.0;
    double dt2 = 0.0;
    sap_Ndof_compute_ramp_dist_and_time(v1, vp, vmax, amax, dp1, dt1);
    sap_Ndof_compute_ramp_dist_and_time(vp, v2, vmax, amax, dp2, dt2);
    return dp - dp1 - dp2 - vp / vmax * (dt - dt1 - dt2);
  }

  [[nodiscard]] double get_delta_time_diff(double vp) const {
    double dp1 = NAN;
    double dt1 = NAN;
    double dp2 = NAN;
    double dt2 = NAN;
    sap_Ndof_compute_ramp_dist_and_time(v1, vp, vmax, amax, dp1, dt1);
    sap_Ndof_compute_ramp_dist_and_time(vp, v2, vmax, amax, dp2, dt2);
    return dt - dt1 - dt2;
  }
};

inline void sap_Ndof_comp_vp_num2_find_sign_change(
    double& vp_lo, double& vp_hi, double tol,
    const sap_Ndof_pos_diff_calculator& pd_calc) {
  double pd_lo = pd_calc(vp_lo);
  double pd_hi = pd_calc(vp_hi);
  double dt_lo = pd_calc.get_delta_time_diff(vp_lo);
  double dt_hi = pd_calc.get_delta_time_diff(vp_hi);

  while ((std::abs(vp_lo - vp_hi) > tol) && (pd_lo * pd_hi > 0.0) &&
         (dt_lo * dt_hi < 0.0)) {
    double vp_mid = (vp_lo + vp_hi) * 0.5;
    double pd_mid = pd_calc(vp_mid);
    double dt_mid = pd_calc.get_delta_time_diff(vp_mid);

    if (dt_mid * dt_hi < 0.0) {
      vp_lo = vp_mid;
      pd_lo = pd_mid;
      dt_lo = dt_mid;
    } else {
      vp_hi = vp_mid;
      pd_hi = pd_mid;
      dt_hi = dt_mid;
    }
  }
}

static void sap_Ndof_compute_peak_velocity_numsolve(
    double start_position, double end_position, double start_velocity,
    double end_velocity, double& peak_velocity, double max_velocity,
    double max_acceleration, double delta_time) {

  using std::abs;

  if ((abs(end_position - start_position) < 1e-6 * max_velocity) &&
      (abs(end_velocity - start_velocity) < 1e-6 * max_acceleration)) {
    peak_velocity = start_velocity;
    return;
  };

  if ((abs(start_velocity) > max_velocity) ||
      (abs(end_velocity) > max_velocity)) {
    peak_velocity = 0.0;
    throw optim::infeasible_problem(
        "Violation of the velocity bounds on invocation of the SAP "
        "peak-velocity solver!");
  }

  double sign_p1_p0 = 1.0;
  if (start_position > end_position) {
    sign_p1_p0 = -1.0;
  }

  sap_Ndof_pos_diff_calculator pd_calc(
      end_position - start_position, start_velocity, end_velocity, max_velocity,
      max_acceleration, delta_time);

  // Here are all the intervals that could exist:
  // vp == vmax;                      // upper boundary.
  // vmax > vp >= (v1 + amax);        // between upper boundary and lowest all-trapezoidal-ramps peak-velocity.
  // (v1 + amax) > vp >= v1;          // when the first ramp is a triangular ramp UP.
  // v1 >= vp > (v1 - amax);          // when the first ramp is a triangular ramp DOWN.
  // (v1 - amax) > vp >= (v2 + amax); // in the regime when both ramps are trapezoidal, and going DOWN and UP.
  // (v2 + amax) > vp >= v2;          // when the second ramp is a triangular ramp DOWN.
  // v2 > vp > (v2 - amax);           // when the second ramp is a triangular ramp UP.
  // (v2 - amax) >= vp > -vmax;       // when the second ramp is a trapezoidal ramp UP.
  // vp == -vmax;                     // lower boundary.

  // This means that the points of interest are:
  //   vmax, v1 + amax, v1, v1 - amax, v2 + amax, v2, v2 - amax, -vmax

  // Lets try to just record all the points of interest and solve between them:
  std::array<double, 8> interest_pts = {
      max_velocity,
      sign_p1_p0 * start_velocity + max_acceleration,
      sign_p1_p0 * start_velocity,
      sign_p1_p0 * start_velocity - max_acceleration,
      sign_p1_p0 * end_velocity + max_acceleration,
      sign_p1_p0 * end_velocity,
      sign_p1_p0 * end_velocity - max_acceleration,
      -max_velocity};
  // using insertion sort because it is the fastest method for such a small array:
  sorting::insertion_sort(interest_pts.begin(), interest_pts.end(),
                          std::greater<>());

  for (int i = 0; i < 7; ++i) {
    if (interest_pts[i] > 1.001 * max_velocity) {
      continue;
    }
    if (interest_pts[i + 1] < -1.001 * max_velocity) {
      break;
    }
    if (interest_pts[i] - interest_pts[i + 1] < 1e-8 * max_velocity) {
      double vp_sol = sign_p1_p0 * interest_pts[i];
      double pd_sol = pd_calc(vp_sol);
      if ((pd_calc.get_delta_time_diff(vp_sol) >= -1e-3 * max_velocity) &&
          (abs(pd_sol) < 1e-3 * max_velocity)) {
        peak_velocity = vp_sol;
        return;
      }
      continue;
    }

    double vp_lo = sign_p1_p0 * interest_pts[i];
    double vp_hi = sign_p1_p0 * interest_pts[i + 1];

    double pd_lo = pd_calc(vp_lo);
    double pd_hi = pd_calc(vp_hi);

    double dt_lo = pd_calc.get_delta_time_diff(vp_lo);
    double dt_hi = pd_calc.get_delta_time_diff(vp_hi);

    if ((dt_lo >= -1e-3 * max_velocity) && (abs(pd_lo) < 1e-3 * max_velocity)) {
      peak_velocity = vp_lo;
      return;
    }
    if ((dt_hi >= -1e-3 * max_velocity) && (abs(pd_hi) < 1e-3 * max_velocity)) {
      peak_velocity = vp_hi;
      return;
    }
    if ((dt_lo < -1e-3 * max_velocity) && (dt_hi < -1e-3 * max_velocity)) {
      continue;
    }

    if (pd_lo * pd_hi < 0.0) {
      //         ford3_method(vp_lo, vp_hi, pd_calc, 1e-9);
      brent_method(vp_lo, vp_hi, pd_calc, 1e-8);
      //         bisection_method(vp_lo, vp_hi, pd_calc, 1e-9 * max_velocity);
      double vp_sol = (vp_lo + vp_hi) * 0.5;
      double pd_sol = pd_calc(vp_sol);
      if ((pd_calc.get_delta_time_diff(vp_sol) >= -1e-3 * max_velocity) &&
          (abs(pd_sol) < 1e-3 * max_velocity)) {
        peak_velocity = vp_sol;
        return;
      }
    } else {
      double vp_mid = 0.5 * (vp_hi + vp_lo);
      double pd_mid = pd_calc(vp_mid);

      if (abs(pd_lo + pd_hi - 2.0 * pd_mid) > 1e-3 * max_velocity) {
        if ((dt_lo * dt_hi > 0.0)) {
          continue;
        }

        sap_Ndof_comp_vp_num2_find_sign_change(vp_lo, vp_hi,
                                               1e-6 * max_velocity, pd_calc);
        //           ford3_method(vp_lo, vp_hi, pd_calc, 1e-9);
        brent_method(vp_lo, vp_hi, pd_calc, 1e-8);
        //           bisection_method(vp_lo, vp_hi, pd_calc, 1e-9 * max_velocity);
        double vp_sol = (vp_lo + vp_hi) * 0.5;
        double pd_sol = pd_calc(vp_sol);
        if ((pd_calc.get_delta_time_diff(vp_sol) >= -1e-3 * max_velocity) &&
            (abs(pd_sol) < 1e-3 * max_velocity)) {
          peak_velocity = vp_sol;
          return;
        }
      }
    }
  }

  peak_velocity = 0.0;
  throw optim::infeasible_problem(
      "The SAP peak-velocity solver could not find a solution for the given "
      "boundary conditions!");
}

void sap_Ndof_compute_peak_velocity(double start_position, double end_position,
                                    double start_velocity, double end_velocity,
                                    double& peak_velocity, double max_velocity,
                                    double max_acceleration,
                                    double delta_time) {
  sap_Ndof_compute_peak_velocity_numsolve(
      start_position, end_position, start_velocity, end_velocity, peak_velocity,
      max_velocity, max_acceleration, delta_time);
}

}  // namespace ReaK::pp::detail

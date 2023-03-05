/**
 * \file sustained_velocity_pulse_detail.hpp
 *
 * This library contains the implementation details of the rate-limited sustained velocity
 * pulse (SVP) interpolation.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
 */

/*
 *    Copyright 2011 Sven Mikael Persson
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

#ifndef REAK_SUSTAINED_VELOCITY_PULSE_DETAIL_HPP
#define REAK_SUSTAINED_VELOCITY_PULSE_DETAIL_HPP

#include "ReaK/math/lin_alg/arithmetic_tuple.hpp"

#include "ReaK/topologies/spaces/bounded_space_concept.hpp"
#include "ReaK/topologies/spaces/tangent_bundle_concept.hpp"

#include <functional>
#include <type_traits>

namespace ReaK::pp::detail {

double svp_solve_for_min_dt_beta(double beta, double norm_delta,
                                 const std::array<double, 5>& coefs,
                                 double num_tol);

bool svp_min_dt_predicate(double beta, double norm_delta,
                          const std::array<double, 5>& coefs, double num_tol,
                          double& result);

double svp_solve_for_no_slack_beta(double beta, double norm_delta,
                                   const std::array<double, 5>& coefs,
                                   double num_tol, double delta_time);

bool svp_no_slack_predicate(double beta, double norm_delta,
                            const std::array<double, 5>& coefs, double num_tol,
                            double& slack, double delta_time);

template <int Order, typename PointType, typename PointDiff1,
          typename DiffSpace, typename TimeSpace>
inline void svp_constant_accel_motion_impl(PointType& result,
                                           const PointDiff1& descended_accel,
                                           const DiffSpace& space,
                                           const TimeSpace& t_space,
                                           double dt) {
  static_assert(Order >= 1);

  if constexpr (Order > 2) {
    svp_constant_accel_motion_impl<Order - 1>(result, descended_accel, space,
                                              t_space, dt);
    get<Order>(result) = get_space<Order>(space, t_space).origin();
    return;
  }

  get<0>(result) =
      get_space<0>(space, t_space)
          .adjust(get<0>(result),
                  descend_to_space<0>(
                      get_space<1>(space, t_space)
                          .adjust(get<1>(result), 0.5 * dt * descended_accel),
                      dt, space, t_space));

  get<1>(result) =
      get_space<1>(space, t_space).adjust(get<1>(result), dt * descended_accel);

  if constexpr (Order > 1) {
    get<2>(result) = lift_to_space<2>(descended_accel, 1.0, space, t_space);
  }
}

template <int Order, typename PointType, typename PointDiff0,
          typename DiffSpace, typename TimeSpace>
inline void svp_constant_vel_motion_impl(PointType& result,
                                         const PointDiff0& descended_vel,
                                         const DiffSpace& space,
                                         const TimeSpace& t_space, double dt) {
  static_assert(Order >= 0);

  if constexpr (Order > 1) {
    svp_constant_vel_motion_impl<Order - 1>(result, descended_vel, space,
                                            t_space, dt);
    get<Order>(result) = get_space<Order>(space, t_space).origin();
    return;
  }

  get<0>(result) =
      get_space<0>(space, t_space).adjust(get<0>(result), dt * descended_vel);

  if constexpr (Order > 0) {
    get<1>(result) = lift_to_space<1>(descended_vel, 1.0, space, t_space);
  }
}

template <int Order, typename PointType, typename PointDiff0,
          typename PointType1, typename DiffSpace, typename TimeSpace>
inline void svp_interpolate_impl(
    PointType& result, const PointType& start_point, const PointType& end_point,
    const PointDiff0& delta_first_order, const PointType1& peak_velocity,
    const DiffSpace& space, const TimeSpace& t_space, double dt,
    double dt_total) {
  static_assert(Order >= 1);

  double dt1 = get(distance_metric, get_space<1>(space, t_space))(
      get<1>(start_point), peak_velocity, get_space<1>(space, t_space));
  double dt2 = get(distance_metric, get_space<1>(space, t_space))(
      peak_velocity, get<1>(end_point), get_space<1>(space, t_space));
  dt_total -= dt1 + dt2;

  result = start_point;

  // Phase 1: constant acceleration to the peak-velocity:

  if (dt1 > std::numeric_limits<double>::epsilon()) {
    svp_constant_accel_motion_impl<Order>(
        result,
        (1.0 / dt1) * get_space<1>(space, t_space)
                          .difference(peak_velocity, get<1>(start_point)),
        space, t_space, (dt > dt1 ? dt1 : dt));
  }
  dt -= dt1;
  if (dt < 0.0) {
    return;
  }

  // Phase 2: constant velocity (or cruise phase):

  if (dt_total > std::numeric_limits<double>::epsilon()) {
    svp_constant_vel_motion_impl<Order>(
        result, descend_to_space<0>(peak_velocity, 1.0, space, t_space), space,
        t_space, (dt > dt_total ? dt_total : dt));
  }
  dt -= dt_total;
  if (dt < 0.0) {
    return;
  }

  // Phase 3: constant acceleration to end-velocity:

  if (dt2 > std::numeric_limits<double>::epsilon()) {
    svp_constant_accel_motion_impl<Order>(
        result,
        (1.0 / dt2) * get_space<1>(space, t_space)
                          .difference(get<1>(end_point), peak_velocity),
        space, t_space, (dt > dt2 ? dt2 : dt));
  }
}

struct svp_update_delta_first_order {
  template <typename PointType, typename PointDiff0, typename PointType1,
            typename DiffSpace, typename TimeSpace>
  void operator()(const PointType& start_point, const PointType& end_point,
                  PointDiff0& delta_first_order, PointType1& peak_velocity,
                  double& norm_delta, const double& beta,
                  const DiffSpace& space, const TimeSpace& t_space,
                  double num_tol = 1E-6) {

    PointDiff0 descended_peak_velocity =
        descend_to_space<0>(peak_velocity, 1.0, space, t_space);
    delta_first_order =
        get_space<0>(space, t_space)
            .difference(
                get_space<0>(space, t_space)
                    .adjust(get<0>(end_point),
                            (-0.5 *
                             get(distance_metric, get_space<1>(space, t_space))(
                                 peak_velocity, get<1>(end_point),
                                 get_space<1>(space, t_space))) *
                                (descended_peak_velocity +
                                 descend_to_space<0>(get<1>(end_point), 1.0,
                                                     space, t_space))),
                get_space<0>(space, t_space)
                    .adjust(get<0>(start_point),
                            (0.5 *
                             get(distance_metric, get_space<1>(space, t_space))(
                                 peak_velocity, get<1>(start_point),
                                 get_space<1>(space, t_space))) *
                                (descended_peak_velocity +
                                 descend_to_space<0>(get<1>(start_point), 1.0,
                                                     space, t_space))));
    norm_delta = get(distance_metric, get_space<0>(space, t_space))(
        delta_first_order, get_space<0>(space, t_space));
    if (norm_delta > num_tol) {
      peak_velocity = lift_to_space<1>(delta_first_order * beta, norm_delta,
                                       space, t_space);
    }
    if (get(distance_metric, get_space<0>(space, t_space))(
            descended_peak_velocity - delta_first_order,
            get_space<0>(space, t_space)) >
        get(distance_metric, get_space<0>(space, t_space))(
            descended_peak_velocity, get_space<0>(space, t_space))) {
      norm_delta = -norm_delta;
    }
  }
};

template <typename PointType, typename PointDiff0, typename PointType1,
          typename DiffSpace, typename TimeSpace, typename UpdateDeltaPFunctor,
          typename SolveForBetaFunctor, typename AddedPredicate>
void svp_peak_velocity_iteration(
    const PointType& start_point, const PointType& end_point,
    PointDiff0& delta_first_order, PointType1& peak_velocity,
    double& norm_delta, double& beta, const DiffSpace& space,
    const TimeSpace& t_space, UpdateDeltaPFunctor update_delta_p,
    SolveForBetaFunctor get_next_beta, AddedPredicate pred,
    double num_tol = 1E-6, unsigned int max_iter = 20) {

  using Space1 = derived_N_order_space_t<DiffSpace, TimeSpace, 1>;

  const Space1& space_1 = get_space<1>(space, t_space);

  std::array<double, 5> coefs = {
      get(distance_metric, space_1)(get<1>(start_point), space_1.origin(),
                                    space_1),
      0.0,
      get(distance_metric, space_1)(get<1>(end_point), space_1.origin(),
                                    space_1),
      0.0, space_1.get_radius()};

  update_delta_p(start_point, end_point, delta_first_order, peak_velocity,
                 norm_delta, beta, space, t_space, num_tol);

  for (unsigned int i = 0; i < max_iter; ++i) {
    PointType1 prev_peak_velocity = peak_velocity;

    peak_velocity = space_1.adjust(peak_velocity,
                                   space_1.get_diff_to_boundary(peak_velocity));

    double temp = get(distance_metric, space_1)(get<1>(start_point),
                                                peak_velocity, space_1);
    coefs[1] = (coefs[0] * coefs[0] + coefs[4] * coefs[4] - temp * temp) /
               (2.0 * coefs[4]);
    temp = get(distance_metric, space_1)(get<1>(end_point), peak_velocity,
                                         space_1);
    coefs[3] = (coefs[2] * coefs[2] + coefs[4] * coefs[4] - temp * temp) /
               (2.0 * coefs[4]);

    beta = get_next_beta(beta, norm_delta, coefs, num_tol);

    peak_velocity = space_1.adjust(
        space_1.origin(),
        beta * space_1.difference(peak_velocity, space_1.origin()));

    update_delta_p(start_point, end_point, delta_first_order, peak_velocity,
                   norm_delta, beta, space, t_space, num_tol);

    if (pred(beta, norm_delta, coefs, num_tol) && (beta > num_tol) &&
        (beta <= 1.0) &&
        (get(distance_metric, space_1)(prev_peak_velocity, peak_velocity,
                                       space_1) <
         1e-6 * space_1.get_radius())) {
      break;
    }
  }
}

template <typename PointType, typename PointDiff0, typename PointType1,
          typename DiffSpace, typename TimeSpace>
double svp_compute_min_delta_time(const PointType& start_point,
                                  const PointType& end_point,
                                  PointDiff0& delta_first_order,
                                  PointType1& peak_velocity, double& norm_delta,
                                  double& beta, const DiffSpace& space,
                                  const TimeSpace& t_space,
                                  double num_tol = 1E-6,
                                  unsigned int max_iter = 20) {

  double result = 0.0;

  svp_peak_velocity_iteration(
      start_point, end_point, delta_first_order, peak_velocity, norm_delta,
      beta, space, t_space, svp_update_delta_first_order(),
      svp_solve_for_min_dt_beta,
      [&](double beta, double norm_delta, const std::array<double, 5>& coefs,
          double num_tol) {
        return svp_min_dt_predicate(beta, norm_delta, coefs, num_tol, result);
      },
      num_tol, max_iter);

  return result;
}

template <typename PointType, typename PointDiff0, typename PointType1,
          typename DiffSpace, typename TimeSpace>
double svp_compute_peak_velocity(
    const PointType& start_point, const PointType& end_point,
    PointDiff0& delta_first_order, PointType1& peak_velocity,
    double& norm_delta, double& beta, double delta_time, const DiffSpace& space,
    const TimeSpace& t_space, double num_tol = 1E-6,
    unsigned int max_iter = 20) {

  double slack = 0.0;

  svp_peak_velocity_iteration(
      start_point, end_point, delta_first_order, peak_velocity, norm_delta,
      beta, space, t_space, svp_update_delta_first_order(),
      [&](double beta, double norm_delta, const std::array<double, 5>& coefs,
          double num_tol) {
        return svp_solve_for_no_slack_beta(beta, norm_delta, coefs, num_tol,
                                           delta_time);
      },
      [&](double beta, double norm_delta, const std::array<double, 5>& coefs,
          double num_tol) {
        return svp_no_slack_predicate(beta, norm_delta, coefs, num_tol, slack,
                                      delta_time);
      },
      num_tol, max_iter);

  return slack;
}

template <typename PointType, typename DiffSpace, typename TimeSpace>
double svp_compute_interpolation_data_impl(
    const PointType& start_point, const PointType& end_point,
    topology_point_difference_type_t<
        derived_N_order_space_t<DiffSpace, TimeSpace, 0>>& delta_first_order,
    topology_point_type_t<derived_N_order_space_t<DiffSpace, TimeSpace, 1>>&
        peak_velocity,
    const DiffSpace& space, const TimeSpace& t_space, double delta_time = 0.0,
    topology_point_type_t<derived_N_order_space_t<DiffSpace, TimeSpace, 1>>*
        best_peak_velocity = nullptr,
    double num_tol = 1e-6, unsigned int max_iter = 20) {
  delta_first_order = get_space<0>(space, t_space)
                          .difference(get<0>(end_point), get<0>(start_point));
  double norm_delta = get(distance_metric, get_space<0>(space, t_space))(
      delta_first_order, get_space<0>(space, t_space));
  double beta = 0.0;
  peak_velocity = get_space<1>(space, t_space).origin();

  double min_delta_time = svp_compute_min_delta_time(
      start_point, end_point, delta_first_order, peak_velocity, norm_delta,
      beta, space, t_space, 1e-6, 20);
  if (best_peak_velocity) {
    *best_peak_velocity = peak_velocity;
  }

  if (min_delta_time > delta_time) {
    return min_delta_time;
  }

  beta = beta * min_delta_time / delta_time;
  peak_velocity =
      get_space<1>(space, t_space)
          .adjust(get_space<1>(space, t_space).origin(),
                  (min_delta_time / delta_time) *
                      get_space<1>(space, t_space)
                          .difference(peak_velocity,
                                      get_space<1>(space, t_space).origin()));

  svp_compute_peak_velocity(start_point, end_point, delta_first_order,
                            peak_velocity, norm_delta, beta, delta_time, space,
                            t_space, 1e-6, 100);

  return min_delta_time;
}

}  // namespace ReaK::pp::detail

#endif

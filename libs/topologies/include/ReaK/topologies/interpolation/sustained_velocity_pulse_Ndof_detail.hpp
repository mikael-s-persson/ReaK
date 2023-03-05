/**
 * \file sustained_velocity_pulse_Ndof_detail.hpp
 *
 * This library contains the implementation details of the rate-limited sustained velocity
 * pulse (SVP) interpolation.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2013
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

#ifndef REAK_SUSTAINED_VELOCITY_PULSE_NDOF_DETAIL_HPP
#define REAK_SUSTAINED_VELOCITY_PULSE_NDOF_DETAIL_HPP

#include "ReaK/math/lin_alg/arithmetic_tuple.hpp"
#include "ReaK/math/optimization/optim_exceptions.hpp"

#include "ReaK/topologies/spaces/tangent_bundle_concept.hpp"

#include <cmath>
#include <type_traits>

namespace ReaK::pp::detail {

double svp_Ndof_compute_min_delta_time(
    double start_position, double end_position, double start_velocity,
    double end_velocity, double& peak_velocity, double max_velocity);

void svp_Ndof_compute_peak_velocity(double start_position, double end_position,
                                    double start_velocity, double end_velocity,
                                    double& peak_velocity, double max_velocity,
                                    double delta_time);

void svp_Ndof_compute_interpolated_values(
    double start_position, double end_position, double start_velocity,
    double end_velocity, double peak_velocity, double max_velocity, double dt,
    double dt_total, double& result_pos, double& result_vel,
    double& result_acc);

template <int Order, int MaxOrder, typename PointType>
inline void svp_Ndof_zero_out_HOT_impl(PointType& result, std::size_t i) {
  if constexpr (Order <= MaxOrder) {
    get<Order>(result)[i] = 0.0;
    svp_Ndof_zero_out_HOT_impl<Order + 1, MaxOrder>(result, i);
  }
}

template <int Order, typename PointType, typename PointType1,
          typename DiffSpace, typename TimeSpace>
inline void svp_Ndof_interpolate_impl(
    PointType& result, const PointType& start_point, const PointType& end_point,
    const PointType1& peak_velocity, const DiffSpace& space,
    const TimeSpace& t_space, double dt, double dt_total) {
  static_assert(Order >= 1);

  PointType1 max_velocity = get_space<1>(space, t_space).get_upper_corner();

  for (std::size_t i = 0; i < peak_velocity.size(); ++i) {
    double result_pos = 0.0;
    double result_vel = 0.0;
    double result_acc = 0.0;

    svp_Ndof_compute_interpolated_values(
        get<0>(start_point)[i], get<0>(end_point)[i], get<1>(start_point)[i],
        get<1>(end_point)[i], peak_velocity[i], max_velocity[i], dt, dt_total,
        result_pos, result_vel, result_acc);

    get<0>(result)[i] = result_pos;
    get<1>(result)[i] = result_vel;
    if constexpr (Order > 1) {
      const auto& max_acceleration =
          get_space<2>(space, t_space).get_upper_corner();
      get<2>(result)[i] = result_acc * max_acceleration[i];
    }
    svp_Ndof_zero_out_HOT_impl<3, Order>(result, i);
  }
}

template <typename PointType, typename DiffSpace, typename TimeSpace>
double svp_compute_Ndof_interpolation_data_impl(
    const PointType& start_point, const PointType& end_point,
    topology_point_type_t<derived_N_order_space_t<DiffSpace, TimeSpace, 1>>&
        peak_velocity,
    const DiffSpace& space, const TimeSpace& t_space, double delta_time = 0.0,
    topology_point_type_t<derived_N_order_space_t<DiffSpace, TimeSpace, 1>>*
        best_peak_velocity = nullptr) {
  auto max_velocity = get_space<1>(space, t_space).get_upper_corner();
  peak_velocity = max_velocity;
  double min_dt_final = 0.0;

  for (std::size_t i = 0; i < peak_velocity.size(); ++i) {
    double vp = 0.0;
    double min_delta_time = svp_Ndof_compute_min_delta_time(
        get<0>(start_point)[i], get<0>(end_point)[i], get<1>(start_point)[i],
        get<1>(end_point)[i], vp, max_velocity[i]);
    peak_velocity[i] = vp;

    if (min_dt_final < min_delta_time) {
      min_dt_final = min_delta_time;
    }
  }

  if (best_peak_velocity) {
    *best_peak_velocity = peak_velocity;
  }

  if (min_dt_final > delta_time) {
    delta_time = min_dt_final;
  }

  for (std::size_t i = 0; i < peak_velocity.size(); ++i) {
    double vp = 0.0;
    svp_Ndof_compute_peak_velocity(get<0>(start_point)[i], get<0>(end_point)[i],
                                   get<1>(start_point)[i], get<1>(end_point)[i],
                                   vp, max_velocity[i], delta_time);
    peak_velocity[i] = vp;
  }

  return min_dt_final;
}

}  // namespace ReaK::pp::detail

#endif

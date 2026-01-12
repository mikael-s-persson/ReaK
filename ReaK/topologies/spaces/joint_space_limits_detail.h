/**
 * \file joint_space_limits_detail.h
 *
 * This library contains the implementation details for classes in the joint_space_limits.hpp file.
 * \note This is not a stand-alone header file, it is simply to be included in the joint_space_limits.hpp header. This
 *"detail" header just serves to tuck away the details (TMP magic).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2012
 */

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

#ifndef REAK_TOPOLOGIES_SPACES_JOINT_SPACE_LIMITS_DETAIL_H_
#define REAK_TOPOLOGIES_SPACES_JOINT_SPACE_LIMITS_DETAIL_H_

#include "ReaK/math/lin_alg/vect_concepts.h"
#include "ReaK/topologies/spaces/differentiable_space.h"
#include "ReaK/topologies/spaces/joint_space_topologies.h"
#include "ReaK/topologies/spaces/ndof_spaces.h"
#include "ReaK/topologies/spaces/se2_topologies.h"
#include "ReaK/topologies/spaces/se3_topologies.h"

#include <type_traits>

namespace ReaK::pp::detail {

/*****************************************************************************************************
                                 FUNCTIONS TO CREATE RATE-LIMITED JOINT-SPACES
******************************************************************************************************/

template <typename OutSpace, typename InSpace, typename RateLimitMap>
void create_rl_joint_space_impl(OutSpace& space_out, const InSpace& space_in,
                                const RateLimitMap& j_limits,
                                std::size_t& gen_i, std::size_t& f2d_i,
                                std::size_t& f3d_i) {
  if constexpr (is_rate_limited_joint_space_v<OutSpace>) {
    constexpr int InOrder = max_derivation_order_v<InSpace, time_topology>;
    using ValueType = typename RateLimitMap::value_type;

    line_segment_topology<ValueType> topo_0(
        get<0>(space_in).get_name() + "_rl",
        (get<0>(space_in).origin() - get<0>(space_in).get_radius()) /
            j_limits.gen_speed_limits[gen_i],
        (get<0>(space_in).origin() + get<0>(space_in).get_radius()) /
            j_limits.gen_speed_limits[gen_i]);
    if constexpr (InOrder > 0) {
      line_segment_topology<ValueType> topo_1(
          get<1>(space_in).get_name() + "_rl",
          (get<1>(space_in).origin() - get<1>(space_in).get_radius()) /
              j_limits.gen_accel_limits[gen_i],
          (get<1>(space_in).origin() + get<1>(space_in).get_radius()) /
              j_limits.gen_accel_limits[gen_i]);
      reach_time_differentiation diff_0(j_limits.gen_speed_limits[gen_i] /
                                        j_limits.gen_accel_limits[gen_i]);
      if constexpr (InOrder > 1) {
        line_segment_topology<ValueType> topo_2(
            get<2>(space_in).get_name() + "_rl",
            (get<2>(space_in).origin() - get<2>(space_in).get_radius()) /
                j_limits.gen_jerk_limits[gen_i],
            (get<2>(space_in).origin() + get<2>(space_in).get_radius()) /
                j_limits.gen_jerk_limits[gen_i]);
        reach_time_differentiation diff_1(j_limits.gen_accel_limits[gen_i] /
                                          j_limits.gen_jerk_limits[gen_i]);

        space_out =
            OutSpace(arithmetic_tuple(std::move(topo_0), std::move(topo_1),
                                      std::move(topo_2)),
                     euclidean_tuple_distance(),
                     arithmetic_tuple(std::move(diff_0), std::move(diff_1)));
      } else {
        space_out = OutSpace(
            arithmetic_tuple(std::move(topo_0), std::move(topo_1)),
            euclidean_tuple_distance(), arithmetic_tuple(std::move(diff_0)));
      }
    } else {
      space_out = OutSpace(arithmetic_tuple(std::move(topo_0)));
    }

    ++gen_i;
  } else if constexpr (is_Ndof_space_v<InSpace>) {
    constexpr int InOrder = max_derivation_order_v<InSpace, time_topology>;
    using BoxTopo = derived_N_order_space_t<InSpace, time_topology, 0>;
    using VectorType = typename BoxTopo::point_type;
    VectorType lower_bnd = get<0>(space_in).get_lower_corner();
    VectorType upper_bnd = get<0>(space_in).get_upper_corner();
    VectorType speed_lim = get<0>(space_in).get_upper_corner();
    VectorType accel_lim = get<0>(space_in).get_upper_corner();
    if constexpr (InOrder > 0) {
      speed_lim = get<1>(space_in).get_upper_corner();
    }
    if constexpr (InOrder > 1) {
      accel_lim = get<2>(space_in).get_upper_corner();
    }
    for (std::size_t i = 0; i < lower_bnd.size(); ++i) {
      lower_bnd[i] /= j_limits.gen_speed_limits[gen_i];
      upper_bnd[i] /= j_limits.gen_speed_limits[gen_i];
      if constexpr (InOrder > 0) {
        speed_lim[i] /= j_limits.gen_accel_limits[gen_i];
      }
      if constexpr (InOrder > 1) {
        accel_lim[i] /= j_limits.gen_jerk_limits[gen_i];
      }
      ++gen_i;
    }
    if constexpr (InOrder == 0) {
      space_out = OutSpace(arithmetic_tuple(
          hyperbox_topology<VectorType, inf_norm_distance_metric>(
              get<0>(space_in).get_name() + "_rl", lower_bnd, upper_bnd)));
    } else if constexpr (InOrder == 1) {
      space_out = OutSpace(
          arithmetic_tuple(
              hyperbox_topology<VectorType, inf_norm_distance_metric>(
                  get<0>(space_in).get_name() + "_rl", lower_bnd, upper_bnd),
              hyperbox_topology<VectorType, inf_norm_distance_metric>(
                  get<1>(space_in).get_name() + "_rl", -speed_lim, speed_lim)),
          manhattan_tuple_distance(),
          arithmetic_tuple(
              Ndof_reach_time_differentiation<VectorType>(speed_lim)));
    } else {
      space_out = OutSpace(
          arithmetic_tuple(
              hyperbox_topology<VectorType, inf_norm_distance_metric>(
                  get<0>(space_in).get_name() + "_rl", lower_bnd, upper_bnd),
              hyperbox_topology<VectorType, inf_norm_distance_metric>(
                  get<1>(space_in).get_name() + "_rl", -speed_lim, speed_lim),
              hyperbox_topology<VectorType, inf_norm_distance_metric>(
                  get<2>(space_in).get_name() + "_rl", -accel_lim, accel_lim)),
          manhattan_tuple_distance(),
          arithmetic_tuple(
              Ndof_reach_time_differentiation<VectorType>(speed_lim),
              Ndof_reach_time_differentiation<VectorType>(accel_lim)));
    }
  } else if constexpr (is_rate_limited_se2_space_v<OutSpace>) {
    constexpr int InOrder =
        max_derivation_order_v<arithmetic_tuple_element_t<0, InSpace>,
                               time_topology>;

    using ValueType = typename RateLimitMap::value_type;
    using LineSegTopo = line_segment_topology<ValueType>;
    using BoxTopo = hyperbox_topology<vect<ValueType, 2>>;
    using BallTopo = hyperball_topology<vect<ValueType, 2>>;
    using PosTopoOutType = arithmetic_tuple_element_t<0, OutSpace>;
    using RotTopoOutType = arithmetic_tuple_element_t<1, OutSpace>;

    BoxTopo pos_0(get<0>(get<0>(space_in)).get_name() + "_rl",
                  get<0>(get<0>(space_in)).get_lower_corner() *
                      (ValueType(1.0) / j_limits.frame2D_speed_limits[f2d_i]),
                  get<0>(get<0>(space_in)).get_upper_corner() *
                      (ValueType(1.0) / j_limits.frame2D_speed_limits[f2d_i]));
    LineSegTopo rot_0(get<0>(get<1>(space_in)).get_name() + "_rl",
                      (get<0>(get<1>(space_in)).origin() -
                       get<0>(get<1>(space_in)).get_radius()) /
                          j_limits.frame2D_speed_limits[f2d_i + 1],
                      (get<0>(get<1>(space_in)).origin() +
                       get<0>(get<1>(space_in)).get_radius()) /
                          j_limits.frame2D_speed_limits[f2d_i + 1]);
    if constexpr (InOrder > 0) {
      BallTopo pos_1(
          get<1>(get<0>(space_in)).get_name() + "_rl",
          get<1>(get<0>(space_in)).origin() *
              (ValueType(1.0) / j_limits.frame2D_accel_limits[f2d_i]),
          get<1>(get<0>(space_in)).get_radius() /
              j_limits.frame2D_accel_limits[f2d_i]);
      LineSegTopo rot_1(get<1>(get<1>(space_in)).get_name() + "_rl",
                        (get<1>(get<1>(space_in)).origin() -
                         get<1>(get<1>(space_in)).get_radius()) /
                            j_limits.frame2D_accel_limits[f2d_i + 1],
                        (get<1>(get<1>(space_in)).origin() +
                         get<1>(get<1>(space_in)).get_radius()) /
                            j_limits.frame2D_accel_limits[f2d_i + 1]);
      reach_time_differentiation pos_diff_0(
          j_limits.frame2D_speed_limits[f2d_i] /
          j_limits.frame2D_accel_limits[f2d_i]);
      reach_time_differentiation rot_diff_0(
          j_limits.frame2D_speed_limits[f2d_i + 1] /
          j_limits.frame2D_accel_limits[f2d_i + 1]);

      if constexpr (InOrder > 1) {
        BallTopo pos_2(
            get<2>(get<0>(space_in)).get_name() + "_rl",
            get<2>(get<0>(space_in)).origin() *
                (ValueType(1.0) / j_limits.frame2D_jerk_limits[f2d_i]),
            get<2>(get<0>(space_in)).get_radius() /
                j_limits.frame2D_jerk_limits[f2d_i]);
        LineSegTopo rot_2(get<2>(get<1>(space_in)).get_name() + "_rl",
                          (get<2>(get<1>(space_in)).origin() -
                           get<2>(get<1>(space_in)).get_radius()) /
                              j_limits.frame2D_jerk_limits[f2d_i + 1],
                          (get<2>(get<1>(space_in)).origin() +
                           get<2>(get<1>(space_in)).get_radius()) /
                              j_limits.frame2D_jerk_limits[f2d_i + 1]);
        reach_time_differentiation pos_diff_1(
            j_limits.frame2D_accel_limits[f2d_i] /
            j_limits.frame2D_jerk_limits[f2d_i]);
        reach_time_differentiation rot_diff_1(
            j_limits.frame2D_accel_limits[f2d_i + 1] /
            j_limits.frame2D_jerk_limits[f2d_i + 1]);

        space_out = OutSpace(arithmetic_tuple(
            PosTopoOutType(
                arithmetic_tuple(std::move(pos_0), std::move(pos_1),
                                 std::move(pos_2)),
                euclidean_tuple_distance(),
                arithmetic_tuple(std::move(pos_diff_0), std::move(pos_diff_1))),
            RotTopoOutType(arithmetic_tuple(std::move(rot_0), std::move(rot_1),
                                            std::move(rot_2)),
                           euclidean_tuple_distance(),
                           arithmetic_tuple(std::move(rot_diff_0),
                                            std::move(rot_diff_1)))));
      } else {
        space_out = OutSpace(arithmetic_tuple(
            PosTopoOutType(arithmetic_tuple(std::move(pos_0), std::move(pos_1)),
                           euclidean_tuple_distance(),
                           arithmetic_tuple(std::move(pos_diff_0))),
            RotTopoOutType(arithmetic_tuple(std::move(rot_0), std::move(rot_1)),
                           euclidean_tuple_distance(),
                           arithmetic_tuple(std::move(rot_diff_0)))));
      }
    } else {
      space_out = OutSpace(
          arithmetic_tuple(PosTopoOutType(arithmetic_tuple(std::move(pos_0))),
                           RotTopoOutType(arithmetic_tuple(std::move(rot_0)))));
    }

    f2d_i += 2;
  } else if constexpr (is_rate_limited_se3_space_v<OutSpace>) {
    constexpr int InOrder =
        max_derivation_order_v<arithmetic_tuple_element_t<0, InSpace>,
                               time_topology>;

    using ValueType = typename RateLimitMap::value_type;
    using QuatTopo = rate_limited_quat_space<ValueType>;
    using AngVelTopo = ang_velocity_3D_topology<ValueType>;
    using AngAccTopo = ang_accel_3D_topology<ValueType>;
    using BoxTopo = hyperbox_topology<vect<ValueType, 3>>;
    using BallTopo = hyperball_topology<vect<ValueType, 3>>;
    using PosTopoOutType = arithmetic_tuple_element_t<0, OutSpace>;
    using RotTopoOutType = arithmetic_tuple_element_t<1, OutSpace>;

    BoxTopo pos_0(get<0>(get<0>(space_in)).get_name() + "_rl",
                  get<0>(get<0>(space_in)).get_lower_corner() *
                      (ValueType(1.0) / j_limits.frame3D_speed_limits[f3d_i]),
                  get<0>(get<0>(space_in)).get_upper_corner() *
                      (ValueType(1.0) / j_limits.frame3D_speed_limits[f3d_i]));
    QuatTopo rot_0(get<0>(get<1>(space_in)).get_name() + "_rl",
                   j_limits.frame3D_speed_limits[f3d_i + 1]);

    if constexpr (InOrder > 0) {
      BallTopo pos_1(
          get<1>(get<0>(space_in)).get_name() + "_rl",
          get<1>(get<0>(space_in)).origin() *
              (ValueType(1.0) / j_limits.frame3D_accel_limits[f3d_i]),
          get<1>(get<0>(space_in)).get_radius() /
              j_limits.frame3D_accel_limits[f3d_i]);
      AngVelTopo rot_1(get<1>(get<1>(space_in)).get_name() + "_rl",
                       get<1>(get<1>(space_in)).get_radius() /
                           j_limits.frame3D_accel_limits[f3d_i + 1]);
      reach_time_differentiation pos_diff_0(
          j_limits.frame3D_speed_limits[f3d_i] /
          j_limits.frame3D_accel_limits[f3d_i]);
      reach_time_differentiation rot_diff_0(
          j_limits.frame3D_speed_limits[f3d_i + 1] /
          j_limits.frame3D_accel_limits[f3d_i + 1]);

      if constexpr (InOrder > 1) {
        BallTopo pos_2(
            get<2>(get<0>(space_in)).get_name() + "_rl",
            get<2>(get<0>(space_in)).origin() *
                (ValueType(1.0) / j_limits.frame3D_jerk_limits[f3d_i]),
            get<2>(get<0>(space_in)).get_radius() /
                j_limits.frame3D_jerk_limits[f3d_i]);
        AngAccTopo rot_2(get<2>(get<1>(space_in)).get_name() + "_rl",
                         get<2>(get<1>(space_in)).get_radius() /
                             j_limits.frame3D_jerk_limits[f3d_i + 1]);
        reach_time_differentiation pos_diff_1(
            j_limits.frame3D_accel_limits[f3d_i] /
            j_limits.frame3D_jerk_limits[f3d_i]);
        reach_time_differentiation rot_diff_1(
            j_limits.frame3D_accel_limits[f3d_i + 1] /
            j_limits.frame3D_jerk_limits[f3d_i + 1]);

        space_out = OutSpace(arithmetic_tuple(
            PosTopoOutType(
                arithmetic_tuple(std::move(pos_0), std::move(pos_1),
                                 std::move(pos_2)),
                euclidean_tuple_distance(),
                arithmetic_tuple(std::move(pos_diff_0), std::move(pos_diff_1))),
            RotTopoOutType(arithmetic_tuple(std::move(rot_0), std::move(rot_1),
                                            std::move(rot_2)),
                           euclidean_tuple_distance(),
                           arithmetic_tuple(std::move(rot_diff_0),
                                            std::move(rot_diff_1)))));
      } else {
        space_out = OutSpace(arithmetic_tuple(
            PosTopoOutType(arithmetic_tuple(std::move(pos_0), std::move(pos_1)),
                           euclidean_tuple_distance(),
                           arithmetic_tuple(std::move(pos_diff_0))),
            RotTopoOutType(arithmetic_tuple(std::move(rot_0), std::move(rot_1)),
                           euclidean_tuple_distance(),
                           arithmetic_tuple(std::move(rot_diff_0)))));
      }
    } else {
      space_out = OutSpace(
          arithmetic_tuple(PosTopoOutType(arithmetic_tuple(std::move(pos_0))),
                           RotTopoOutType(arithmetic_tuple(std::move(rot_0)))));
    }

    f3d_i += 2;
  } else {
    arithmetic_tuple_details::tuple_for_each(
        space_out, space_in,
        [&](auto& space_elem_out, const auto& space_elem_in) {
          create_rl_joint_space_impl(space_elem_out, space_elem_in, j_limits,
                                     gen_i, f2d_i, f3d_i);
          return int{};
        });
  }
}

template <typename OutSpace, typename InSpace, typename RateLimitMap>
void create_rl_joint_spaces_impl(OutSpace& space_out, const InSpace& space_in,
                                 const RateLimitMap& j_limits) {
  std::size_t gen_i = 0;
  std::size_t f2d_i = 0;
  std::size_t f3d_i = 0;
  create_rl_joint_space_impl(space_out, space_in, j_limits, gen_i, f2d_i,
                             f3d_i);
}

/***********************************************************************************************************
                                 FUNCTIONS TO CREATE NORMAL JOINT-SPACES
************************************************************************************************************/

template <typename OutSpace, typename InSpace, typename RateLimitMap>
void create_normal_joint_space_impl(OutSpace& space_out,
                                    const InSpace& space_in,
                                    const RateLimitMap& j_limits,
                                    std::size_t& gen_i, std::size_t& f2d_i,
                                    std::size_t& f3d_i) {
  if constexpr (is_rate_limited_joint_space_v<InSpace>) {
    constexpr int InOrder = max_derivation_order_v<InSpace, time_topology>;
    using ValueType = typename RateLimitMap::value_type;

    line_segment_topology<ValueType> topo_0(
        get<0>(space_in).get_name() + "_non_rl",
        (get<0>(space_in).origin() - get<0>(space_in).get_radius()) *
            j_limits.gen_speed_limits[gen_i],
        (get<0>(space_in).origin() + get<0>(space_in).get_radius()) *
            j_limits.gen_speed_limits[gen_i]);

    if constexpr (InOrder > 0) {
      line_segment_topology<ValueType> topo_1(
          get<1>(space_in).get_name() + "_non_rl",
          (get<1>(space_in).origin() - get<1>(space_in).get_radius()) *
              j_limits.gen_accel_limits[gen_i],
          (get<1>(space_in).origin() + get<1>(space_in).get_radius()) *
              j_limits.gen_accel_limits[gen_i]);
      if constexpr (InOrder > 1) {
        line_segment_topology<ValueType> topo_2(
            get<2>(space_in).get_name() + "_non_rl",
            (get<2>(space_in).origin() - get<2>(space_in).get_radius()) *
                j_limits.gen_jerk_limits[gen_i],
            (get<2>(space_in).origin() + get<2>(space_in).get_radius()) *
                j_limits.gen_jerk_limits[gen_i]);
        space_out = OutSpace(arithmetic_tuple(
            std::move(topo_0), std::move(topo_1), std::move(topo_2)));
      } else {
        space_out =
            OutSpace(arithmetic_tuple(std::move(topo_0), std::move(topo_1)));
      }
    } else {
      space_out = OutSpace(arithmetic_tuple(std::move(topo_0)));
    }

    ++gen_i;
  } else if constexpr (is_Ndof_rl_space_v<InSpace>) {
    constexpr int InOrder = max_derivation_order_v<InSpace, time_topology>;

    using BoxTopo = derived_N_order_space_t<InSpace, time_topology, 0>;
    using VectorType = typename BoxTopo::point_type;
    VectorType lower_bnd = get<0>(space_in).get_lower_corner();
    VectorType upper_bnd = get<0>(space_in).get_upper_corner();
    VectorType speed_lim = get<0>(space_in).get_upper_corner();
    VectorType accel_lim = get<0>(space_in).get_upper_corner();
    if constexpr (InOrder > 0) {
      speed_lim = get<1>(space_in).get_upper_corner();
    }
    if constexpr (InOrder > 1) {
      accel_lim = get<2>(space_in).get_upper_corner();
    }
    for (std::size_t i = 0; i < lower_bnd.size(); ++i) {
      lower_bnd[i] *= j_limits.gen_speed_limits[gen_i];
      upper_bnd[i] *= j_limits.gen_speed_limits[gen_i];
      if constexpr (InOrder > 0) {
        speed_lim[i] *= j_limits.gen_accel_limits[gen_i];
      }
      if constexpr (InOrder > 1) {
        accel_lim[i] *= j_limits.gen_jerk_limits[gen_i];
      }
      ++gen_i;
    }
    hyperbox_topology<VectorType, manhattan_distance_metric> topo_0(
        get<0>(space_in).get_name() + "_non_rl", lower_bnd, upper_bnd);
    if constexpr (InOrder > 0) {
      hyperbox_topology<VectorType, manhattan_distance_metric> topo_1(
          get<1>(space_in).get_name() + "_non_rl", -speed_lim, speed_lim);
      if constexpr (InOrder > 1) {
        hyperbox_topology<VectorType, manhattan_distance_metric> topo_2(
            get<2>(space_in).get_name() + "_non_rl", -accel_lim, accel_lim);
        space_out = OutSpace(arithmetic_tuple(
            std::move(topo_0), std::move(topo_1), std::move(topo_2)));
      } else {
        space_out =
            OutSpace(arithmetic_tuple(std::move(topo_0), std::move(topo_1)));
      }
    } else {
      space_out = OutSpace(arithmetic_tuple(std::move(topo_0)));
    }
  } else if constexpr (is_rate_limited_se2_space_v<InSpace>) {
    constexpr int InOrder =
        max_derivation_order_v<arithmetic_tuple_element_t<0, InSpace>,
                               time_topology>;

    using ValueType = typename RateLimitMap::value_type;
    using LineSegTopo = line_segment_topology<ValueType>;
    using BoxTopo = hyperbox_topology<vect<ValueType, 2>>;
    using BallTopo = hyperball_topology<vect<ValueType, 2>>;
    using PosTopoOutType = arithmetic_tuple_element_t<0, OutSpace>;
    using RotTopoOutType = arithmetic_tuple_element_t<1, OutSpace>;

    BoxTopo pos_0(get<0>(get<0>(space_in)).get_name() + "_non_rl",
                  get<0>(get<0>(space_in)).get_lower_corner() *
                      j_limits.frame2D_speed_limits[f2d_i],
                  get<0>(get<0>(space_in)).get_upper_corner() *
                      j_limits.frame2D_speed_limits[f2d_i]);
    LineSegTopo rot_0(get<0>(get<1>(space_in)).get_name() + "_non_rl",
                      (get<0>(get<1>(space_in)).origin() -
                       get<0>(get<1>(space_in)).get_radius()) *
                          j_limits.frame2D_speed_limits[f2d_i + 1],
                      (get<0>(get<1>(space_in)).origin() +
                       get<0>(get<1>(space_in)).get_radius()) *
                          j_limits.frame2D_speed_limits[f2d_i + 1]);

    if constexpr (InOrder > 0) {
      BallTopo pos_1(get<1>(get<0>(space_in)).get_name() + "_non_rl",
                     get<1>(get<0>(space_in)).origin() *
                         j_limits.frame2D_accel_limits[f2d_i],
                     get<1>(get<0>(space_in)).get_radius() *
                         j_limits.frame2D_accel_limits[f2d_i]);
      LineSegTopo rot_1(get<1>(get<1>(space_in)).get_name() + "_non_rl",
                        (get<1>(get<1>(space_in)).origin() -
                         get<1>(get<1>(space_in)).get_radius()) *
                            j_limits.frame2D_accel_limits[f2d_i + 1],
                        (get<1>(get<1>(space_in)).origin() +
                         get<1>(get<1>(space_in)).get_radius()) *
                            j_limits.frame2D_accel_limits[f2d_i + 1]);
      if constexpr (InOrder > 1) {
        BallTopo pos_2(get<2>(get<0>(space_in)).get_name() + "_non_rl",
                       get<2>(get<0>(space_in)).origin() *
                           j_limits.frame2D_jerk_limits[f2d_i],
                       get<2>(get<0>(space_in)).get_radius() *
                           j_limits.frame2D_jerk_limits[f2d_i]);
        LineSegTopo rot_2(get<2>(get<1>(space_in)).get_name() + "_non_rl",
                          (get<2>(get<1>(space_in)).origin() -
                           get<2>(get<1>(space_in)).get_radius()) *
                              j_limits.frame2D_jerk_limits[f2d_i + 1],
                          (get<2>(get<1>(space_in)).origin() +
                           get<2>(get<1>(space_in)).get_radius()) *
                              j_limits.frame2D_jerk_limits[f2d_i + 1]);
        space_out = OutSpace(arithmetic_tuple(
            PosTopoOutType(arithmetic_tuple(std::move(pos_0), std::move(pos_1),
                                            std::move(pos_2))),
            RotTopoOutType(arithmetic_tuple(std::move(rot_0), std::move(rot_1),
                                            std::move(rot_2)))));
      } else {
        space_out = OutSpace(arithmetic_tuple(
            PosTopoOutType(
                arithmetic_tuple(std::move(pos_0), std::move(pos_1))),
            RotTopoOutType(
                arithmetic_tuple(std::move(rot_0), std::move(rot_1)))));
      }
    } else {
      space_out = OutSpace(
          arithmetic_tuple(PosTopoOutType(arithmetic_tuple(std::move(pos_0))),
                           RotTopoOutType(arithmetic_tuple(std::move(rot_0)))));
    }

    f2d_i += 2;
  } else if constexpr (is_rate_limited_se3_space_v<InSpace>) {
    constexpr int InOrder =
        max_derivation_order_v<arithmetic_tuple_element_t<0, InSpace>,
                               time_topology>;

    using ValueType = typename RateLimitMap::value_type;
    using QuatTopo = quaternion_topology<ValueType>;
    using AngVelTopo = ang_velocity_3D_topology<ValueType>;
    using AngAccTopo = ang_accel_3D_topology<ValueType>;
    using BoxTopo = hyperbox_topology<vect<ValueType, 3>>;
    using BallTopo = hyperball_topology<vect<ValueType, 3>>;
    using PosTopoOutType = arithmetic_tuple_element_t<0, OutSpace>;
    using RotTopoOutType = arithmetic_tuple_element_t<1, OutSpace>;

    BoxTopo pos_0(get<0>(get<0>(space_in)).get_name() + "_non_rl",
                  get<0>(get<0>(space_in)).get_lower_corner() *
                      j_limits.frame3D_speed_limits[f3d_i],
                  get<0>(get<0>(space_in)).get_upper_corner() *
                      j_limits.frame3D_speed_limits[f3d_i]);
    QuatTopo rot_0(get<0>(get<1>(space_in)).get_name() + "_non_rl");

    if constexpr (InOrder > 0) {
      BallTopo pos_1(get<1>(get<0>(space_in)).get_name() + "_non_rl",
                     get<1>(get<0>(space_in)).origin() *
                         j_limits.frame3D_accel_limits[f3d_i],
                     get<1>(get<0>(space_in)).get_radius() *
                         j_limits.frame3D_accel_limits[f3d_i]);
      AngVelTopo rot_1(get<1>(get<1>(space_in)).get_name() + "_non_rl",
                       get<1>(get<1>(space_in)).get_radius() *
                           j_limits.frame3D_accel_limits[f3d_i + 1]);
      if constexpr (InOrder > 1) {
        BallTopo pos_2(get<2>(get<0>(space_in)).get_name() + "_non_rl",
                       get<2>(get<0>(space_in)).origin() *
                           j_limits.frame3D_jerk_limits[f3d_i],
                       get<2>(get<0>(space_in)).get_radius() *
                           j_limits.frame3D_jerk_limits[f3d_i]);
        AngAccTopo rot_2(get<2>(get<1>(space_in)).get_name() + "_non_rl",
                         get<2>(get<1>(space_in)).get_radius() *
                             j_limits.frame3D_jerk_limits[f3d_i + 1]);
        space_out = OutSpace(arithmetic_tuple(
            PosTopoOutType(arithmetic_tuple(std::move(pos_0), std::move(pos_1),
                                            std::move(pos_2))),
            RotTopoOutType(arithmetic_tuple(std::move(rot_0), std::move(rot_1),
                                            std::move(rot_2)))));
      } else {
        space_out = OutSpace(arithmetic_tuple(
            PosTopoOutType(
                arithmetic_tuple(std::move(pos_0), std::move(pos_1))),
            RotTopoOutType(
                arithmetic_tuple(std::move(rot_0), std::move(rot_1)))));
      }
    } else {
      space_out = OutSpace(
          arithmetic_tuple(PosTopoOutType(arithmetic_tuple(std::move(pos_0))),
                           RotTopoOutType(arithmetic_tuple(std::move(rot_0)))));
    }

    f3d_i += 2;
  } else {
    arithmetic_tuple_details::tuple_for_each(
        space_out, space_in,
        [&](auto& space_elem_out, const auto& space_elem_in) {
          create_normal_joint_space_impl(space_elem_out, space_elem_in,
                                         j_limits, gen_i, f2d_i, f3d_i);
          return int{};
        });
  }
}

template <typename OutSpace, typename InSpace, typename RateLimitMap>
void create_normal_joint_spaces_impl(OutSpace& space_out,
                                     const InSpace& space_in,
                                     const RateLimitMap& j_limits) {
  std::size_t gen_i = 0;
  std::size_t f2d_i = 0;
  std::size_t f3d_i = 0;
  create_normal_joint_space_impl(space_out, space_in, j_limits, gen_i, f2d_i,
                                 f3d_i);
}

/**************************************************************************************************
                                 FUNCTIONS TO CREATE RATE-LIMITED JOINT-SPACE VECTORS
***************************************************************************************************/

template <typename RateLimitMap, typename Arg0, typename... Args>
void create_rl_joint_vector_impl(arithmetic_tuple<Arg0, Args...>& result,
                                 const arithmetic_tuple<Arg0, Args...>& pt,
                                 const RateLimitMap& j_limits,
                                 std::size_t& gen_i, std::size_t& f2d_i,
                                 std::size_t& f3d_i) {
  using ValueType = typename RateLimitMap::value_type;
  if constexpr (std::is_same_v<Arg0, ValueType>) {
    constexpr int Order = sizeof...(Args);
    get<0>(result) = get<0>(pt) / j_limits.gen_speed_limits[gen_i];
    if constexpr (Order > 0) {
      get<1>(result) = get<1>(pt) / j_limits.gen_accel_limits[gen_i];
    }
    if constexpr (Order > 1) {
      get<2>(result) = get<2>(pt) / j_limits.gen_jerk_limits[gen_i];
    }
    ++gen_i;
  } else if constexpr (WritableVector<Arg0>) {
    constexpr int Order = sizeof...(Args);
    for (std::size_t i = 0; i < get<0>(pt).size(); ++i) {
      get<0>(result)[i] = get<0>(pt)[i] / j_limits.gen_speed_limits[gen_i];
      if constexpr (Order > 0) {
        get<1>(result)[i] = get<1>(pt)[i] / j_limits.gen_accel_limits[gen_i];
      }
      if constexpr (Order > 1) {
        get<2>(result)[i] = get<2>(pt)[i] / j_limits.gen_jerk_limits[gen_i];
      }
      ++gen_i;
    }
  } else if constexpr (is_arithmetic_tuple_v<Arg0>) {
    constexpr int Order = arithmetic_tuple_size_v<Arg0> - 1;
    using SubArg0 = std::decay_t<decltype(get<0>(get<0>(pt)))>;
    if constexpr (std::is_same_v<SubArg0, vect<ValueType, 2>>) {
      // SE2
      get<0>(get<0>(result)) =
          get<0>(get<0>(pt)) * (1.0 / j_limits.frame2D_speed_limits[f2d_i]);
      get<0>(get<1>(result)) =
          get<0>(get<1>(pt)) / j_limits.frame2D_speed_limits[f2d_i + 1];
      if constexpr (Order > 0) {
        get<1>(get<0>(result)) =
            get<1>(get<0>(pt)) * (1.0 / j_limits.frame2D_accel_limits[f2d_i]);
        get<1>(get<1>(result)) =
            get<1>(get<1>(pt)) / j_limits.frame2D_accel_limits[f2d_i + 1];
      }
      if constexpr (Order > 1) {
        get<2>(get<0>(result)) =
            get<2>(get<0>(pt)) * (1.0 / j_limits.frame2D_jerk_limits[f2d_i]);
        get<2>(get<1>(result)) =
            get<2>(get<1>(pt)) / j_limits.frame2D_jerk_limits[f2d_i + 1];
      }
      f2d_i += 2;
    } else if constexpr (std::is_same_v<SubArg0, vect<ValueType, 3>>) {
      // SE3
      get<0>(get<0>(result)) =
          get<0>(get<0>(pt)) * (1.0 / j_limits.frame3D_speed_limits[f3d_i]);
      get<0>(get<1>(result)) = get<0>(get<1>(pt));
      if constexpr (Order > 0) {
        get<1>(get<0>(result)) =
            get<1>(get<0>(pt)) * (1.0 / j_limits.frame3D_accel_limits[f3d_i]);
        get<1>(get<1>(result)) =
            get<1>(get<1>(pt)) *
            (1.0 / j_limits.frame3D_accel_limits[f3d_i + 1]);
      }
      if constexpr (Order > 1) {
        get<2>(get<0>(result)) =
            get<2>(get<0>(pt)) * (1.0 / j_limits.frame3D_jerk_limits[f3d_i]);
        get<2>(get<1>(result)) =
            get<2>(get<1>(pt)) *
            (1.0 / j_limits.frame3D_jerk_limits[f3d_i + 1]);
      }
      f3d_i += 2;
    } else {
      // Unknown, just recurse.
      arithmetic_tuple_details::tuple_for_each(
          result, pt, [&](auto& elem_out, const auto& elem_in) {
            create_rl_joint_vector_impl(elem_out, elem_in, j_limits, gen_i,
                                        f2d_i, f3d_i);
            return int{};
          });
    }
  } else {
    // Unknown, just recurse (will probably fail).
    arithmetic_tuple_details::tuple_for_each(
        result, pt, [&](auto& elem_out, const auto& elem_in) {
          create_rl_joint_vector_impl(elem_out, elem_in, j_limits, gen_i, f2d_i,
                                      f3d_i);
          return int{};
        });
  }
}

template <typename OutPoint, typename InPoint, typename RateLimitMap>
void create_rl_joint_vectors_impl(OutPoint& result, const InPoint& pt,
                                  const RateLimitMap& j_limits) {
  std::size_t gen_i = 0;
  std::size_t f2d_i = 0;
  std::size_t f3d_i = 0;
  create_rl_joint_vector_impl(result, pt, j_limits, gen_i, f2d_i, f3d_i);
}

/*******************************************************************************************************************
                                 FUNCTIONS TO CREATE NORMAL JOINT-SPACE VECTORS
*******************************************************************************************************************/

template <typename RateLimitMap, typename Arg0, typename... Args>
void create_normal_joint_vector_impl(arithmetic_tuple<Arg0, Args...>& result,
                                     const arithmetic_tuple<Arg0, Args...>& pt,
                                     const RateLimitMap& j_limits,
                                     std::size_t& gen_i, std::size_t& f2d_i,
                                     std::size_t& f3d_i) {
  using ValueType = typename RateLimitMap::value_type;
  if constexpr (std::is_same_v<Arg0, ValueType>) {
    constexpr int Order = sizeof...(Args);
    get<0>(result) = get<0>(pt) * j_limits.gen_speed_limits[gen_i];
    if constexpr (Order > 0) {
      get<1>(result) = get<1>(pt) * j_limits.gen_accel_limits[gen_i];
    }
    if constexpr (Order > 1) {
      get<2>(result) = get<2>(pt) * j_limits.gen_jerk_limits[gen_i];
    }
    ++gen_i;
  } else if constexpr (WritableVector<Arg0>) {
    constexpr int Order = sizeof...(Args);
    for (std::size_t i = 0; i < get<0>(pt).size(); ++i) {
      get<0>(result)[i] = get<0>(pt)[i] * j_limits.gen_speed_limits[gen_i];
      if constexpr (Order > 0) {
        get<1>(result)[i] = get<1>(pt)[i] * j_limits.gen_accel_limits[gen_i];
      }
      if constexpr (Order > 1) {
        get<2>(result)[i] = get<2>(pt)[i] * j_limits.gen_jerk_limits[gen_i];
      }
      ++gen_i;
    }
  } else if constexpr (is_arithmetic_tuple_v<Arg0>) {
    constexpr int Order = arithmetic_tuple_size_v<Arg0> - 1;
    using SubArg0 = std::decay_t<decltype(get<0>(get<0>(pt)))>;
    if constexpr (std::is_same_v<SubArg0, vect<ValueType, 2>>) {
      // SE2
      get<0>(get<0>(result)) =
          get<0>(get<0>(pt)) * j_limits.frame2D_speed_limits[f2d_i];
      get<0>(get<1>(result)) =
          get<0>(get<1>(pt)) * j_limits.frame2D_speed_limits[f2d_i + 1];
      if constexpr (Order > 0) {
        get<1>(get<0>(result)) =
            get<1>(get<0>(pt)) * j_limits.frame2D_accel_limits[f2d_i];
        get<1>(get<1>(result)) =
            get<1>(get<1>(pt)) * j_limits.frame2D_accel_limits[f2d_i + 1];
      }
      if constexpr (Order > 1) {
        get<2>(get<0>(result)) =
            get<2>(get<0>(pt)) * j_limits.frame2D_jerk_limits[f2d_i];
        get<2>(get<1>(result)) =
            get<2>(get<1>(pt)) * j_limits.frame2D_jerk_limits[f2d_i + 1];
      }
      f2d_i += 2;
    } else if constexpr (std::is_same_v<SubArg0, vect<ValueType, 3>>) {
      // SE3
      get<0>(get<0>(result)) =
          get<0>(get<0>(pt)) * j_limits.frame3D_speed_limits[f3d_i];
      get<0>(get<1>(result)) = get<0>(get<1>(pt));
      if constexpr (Order > 0) {
        get<1>(get<0>(result)) =
            get<1>(get<0>(pt)) * j_limits.frame3D_accel_limits[f3d_i];
        get<1>(get<1>(result)) =
            get<1>(get<1>(pt)) * j_limits.frame3D_accel_limits[f3d_i + 1];
      }
      if constexpr (Order > 1) {
        get<2>(get<0>(result)) =
            get<2>(get<0>(pt)) * j_limits.frame3D_jerk_limits[f3d_i];
        get<2>(get<1>(result)) =
            get<2>(get<1>(pt)) * j_limits.frame3D_jerk_limits[f3d_i + 1];
      }
      f3d_i += 2;
    } else {
      // Unknown, just recurse.
      arithmetic_tuple_details::tuple_for_each(
          result, pt, [&](auto& elem_out, const auto& elem_in) {
            create_normal_joint_vector_impl(elem_out, elem_in, j_limits, gen_i,
                                            f2d_i, f3d_i);
            return int{};
          });
    }
  } else {
    // Unknown, just recurse (will probably fail).
    arithmetic_tuple_details::tuple_for_each(
        result, pt, [&](auto& elem_out, const auto& elem_in) {
          create_normal_joint_vector_impl(elem_out, elem_in, j_limits, gen_i,
                                          f2d_i, f3d_i);
          return int{};
        });
  }
}

template <typename OutPoint, typename InPoint, typename RateLimitMap>
void create_normal_joint_vectors_impl(OutPoint& result, const InPoint& pt,
                                      const RateLimitMap& j_limits) {
  std::size_t gen_i = 0;
  std::size_t f2d_i = 0;
  std::size_t f3d_i = 0;
  create_normal_joint_vector_impl(result, pt, j_limits, gen_i, f2d_i, f3d_i);
}

}  // namespace ReaK::pp::detail

#endif  // REAK_TOPOLOGIES_SPACES_JOINT_SPACE_LIMITS_DETAIL_H_

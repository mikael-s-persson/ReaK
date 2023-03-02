/**
 * \file se2_topologies.hpp
 *
 * This library provides classes that define topologies on SE(2) (2D rigid-body motion).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2012
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

#ifndef REAK_SE2_TOPOLOGIES_HPP
#define REAK_SE2_TOPOLOGIES_HPP

#include <ReaK/core/base/defs.hpp>

#include "differentiable_space.hpp"
#include "hyperball_topology.hpp"
#include "hyperbox_topology.hpp"
#include "line_topology.hpp"
#include "metric_space_tuple.hpp"
#include "rate_limited_spaces.hpp"

#include <ReaK/math/kinetostatics/frame_2D.hpp>
#include <ReaK/math/lin_alg/arithmetic_tuple.hpp>
#include <ReaK/math/lin_alg/vect_alg.hpp>

namespace ReaK {
namespace pp {

/**
 * This meta-function defines the type for a 0th order SE(2) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se2_0th_order_topology {
  using type = metric_space_tuple<
      arithmetic_tuple<
          differentiable_space<time_topology,
                               arithmetic_tuple<hyperbox_topology<vect<T, 2>>>,
                               DistanceMetric>,
          differentiable_space<time_topology,
                               arithmetic_tuple<line_segment_topology<T>>,
                               DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using se2_0th_order_topology_t =
    typename se2_0th_order_topology<T, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 1st order SE(2) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se2_1st_order_topology {
  using type = metric_space_tuple<
      arithmetic_tuple<
          differentiable_space<time_topology,
                               arithmetic_tuple<hyperbox_topology<vect<T, 2>>,
                                                hyperball_topology<vect<T, 2>>>,
                               DistanceMetric>,
          differentiable_space<time_topology,
                               arithmetic_tuple<line_segment_topology<T>,
                                                line_segment_topology<T>>,
                               DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using se2_1st_order_topology_t =
    typename se2_1st_order_topology<T, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 2nd order SE(2) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se2_2nd_order_topology {
  using type = metric_space_tuple<
      arithmetic_tuple<
          differentiable_space<time_topology,
                               arithmetic_tuple<hyperbox_topology<vect<T, 2>>,
                                                hyperball_topology<vect<T, 2>>,
                                                hyperball_topology<vect<T, 2>>>,
                               DistanceMetric>,
          differentiable_space<time_topology,
                               arithmetic_tuple<line_segment_topology<T>,
                                                line_segment_topology<T>,
                                                line_segment_topology<T>>,
                               DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using se2_2nd_order_topology_t =
    typename se2_2nd_order_topology<T, DistanceMetric>::type;

template <typename T, int Order,
          typename DistanceMetric = euclidean_tuple_distance>
struct se2_topology {
  using type = std::conditional_t<
      (Order == 0), se2_0th_order_topology_t<T, DistanceMetric>,
      std::conditional_t<(Order == 1),
                         se2_1st_order_topology_t<T, DistanceMetric>,
                         se2_2nd_order_topology_t<T, DistanceMetric>>>;
};

template <typename T, int Order,
          typename DistanceMetric = euclidean_tuple_distance>
using se2_topology_t = typename se2_topology<T, Order, DistanceMetric>::type;

template <typename SE2Space>
struct is_se2_space : std::false_type {};

template <typename SE2Space>
static constexpr bool is_se2_space_v = is_se2_space<SE2Space>::value;

template <typename T, typename DistanceMetric>
struct is_se2_space<metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<time_topology,
                             arithmetic_tuple<hyperbox_topology<vect<T, 2>>>,
                             DistanceMetric>,
        differentiable_space<time_topology,
                             arithmetic_tuple<line_segment_topology<T>>,
                             DistanceMetric>>,
    DistanceMetric>> : std::true_type {};

template <typename T, typename DistanceMetric>
struct is_se2_space<metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<time_topology,
                             arithmetic_tuple<hyperbox_topology<vect<T, 2>>,
                                              hyperball_topology<vect<T, 2>>>,
                             DistanceMetric>,
        differentiable_space<time_topology,
                             arithmetic_tuple<line_segment_topology<T>,
                                              line_segment_topology<T>>,
                             DistanceMetric>>,
    DistanceMetric>> : std::true_type {};

template <typename T, typename DistanceMetric>
struct is_se2_space<metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<time_topology,
                             arithmetic_tuple<hyperbox_topology<vect<T, 2>>,
                                              hyperball_topology<vect<T, 2>>,
                                              hyperball_topology<vect<T, 2>>>,
                             DistanceMetric>,
        differentiable_space<
            time_topology,
            arithmetic_tuple<line_segment_topology<T>, line_segment_topology<T>,
                             line_segment_topology<T>>,
            DistanceMetric>>,
    DistanceMetric>> : std::true_type {};

/**
 * This meta-function defines the type for a 0th order SE(2) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se2_0th_order_rl_topology {
  using type = metric_space_tuple<
      arithmetic_tuple<
          reach_time_diff_space<time_topology,
                                arithmetic_tuple<hyperbox_topology<vect<T, 2>>>,
                                DistanceMetric>,
          reach_time_diff_space<time_topology,
                                arithmetic_tuple<line_segment_topology<T>>,
                                DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using se2_0th_order_rl_topology_t =
    typename se2_0th_order_rl_topology<T, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 1st order SE(2) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se2_1st_order_rl_topology {
  using type = metric_space_tuple<
      arithmetic_tuple<
          reach_time_diff_space<
              time_topology,
              arithmetic_tuple<hyperbox_topology<vect<T, 2>>,
                               hyperball_topology<vect<T, 2>>>,
              DistanceMetric>,
          reach_time_diff_space<time_topology,
                                arithmetic_tuple<line_segment_topology<T>,
                                                 line_segment_topology<T>>,
                                DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using se2_1st_order_rl_topology_t =
    typename se2_1st_order_rl_topology<T, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 2nd order SE(2) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se2_2nd_order_rl_topology {
  using type = metric_space_tuple<
      arithmetic_tuple<
          reach_time_diff_space<
              time_topology,
              arithmetic_tuple<hyperbox_topology<vect<T, 2>>,
                               hyperball_topology<vect<T, 2>>,
                               hyperball_topology<vect<T, 2>>>,
              DistanceMetric>,
          reach_time_diff_space<time_topology,
                                arithmetic_tuple<line_segment_topology<T>,
                                                 line_segment_topology<T>,
                                                 line_segment_topology<T>>,
                                DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using se2_2nd_order_rl_topology_t =
    typename se2_2nd_order_rl_topology<T, DistanceMetric>::type;

template <typename T, int Order,
          typename DistanceMetric = euclidean_tuple_distance>
struct se2_rl_topology {
  using type = std::conditional_t<
      (Order == 0), se2_0th_order_rl_topology_t<T, DistanceMetric>,
      std::conditional_t<(Order == 1),
                         se2_1st_order_rl_topology_t<T, DistanceMetric>,
                         se2_2nd_order_rl_topology_t<T, DistanceMetric>>>;
};

template <typename T, int Order,
          typename DistanceMetric = euclidean_tuple_distance>
using se2_rl_topology_t =
    typename se2_rl_topology<T, Order, DistanceMetric>::type;

template <typename SE2Space>
struct is_rate_limited_se2_space : std::false_type {};

template <typename SE2Space>
static constexpr bool is_rate_limited_se2_space_v =
    is_rate_limited_se2_space<SE2Space>::value;

template <typename T, typename DistanceMetric>
struct is_rate_limited_se2_space<metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<hyperbox_topology<vect<T, 2>>>,
                              DistanceMetric>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<line_segment_topology<T>>,
                              DistanceMetric>>,
    DistanceMetric>> : std::true_type {};

template <typename T, typename DistanceMetric>
struct is_rate_limited_se2_space<metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<hyperbox_topology<vect<T, 2>>,
                                               hyperball_topology<vect<T, 2>>>,
                              DistanceMetric>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<line_segment_topology<T>,
                                               line_segment_topology<T>>,
                              DistanceMetric>>,
    DistanceMetric>> : std::true_type {};

template <typename T, typename DistanceMetric>
struct is_rate_limited_se2_space<metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<hyperbox_topology<vect<T, 2>>,
                                               hyperball_topology<vect<T, 2>>,
                                               hyperball_topology<vect<T, 2>>>,
                              DistanceMetric>,
        reach_time_diff_space<
            time_topology,
            arithmetic_tuple<line_segment_topology<T>, line_segment_topology<T>,
                             line_segment_topology<T>>,
            DistanceMetric>>,
    DistanceMetric>> : std::true_type {};

}  // namespace pp

// Because of ADL rules, the get functions for the arithmetic-tuple types that represent SE(3) states should be in the
// ReaK namespace.

template <typename T>
frame_2D<T> get_frame_2D(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>, vect<T, 2>>,
                           arithmetic_tuple<T, T, T>>& pt) {
  return frame_2D<T>(std::weak_ptr<pose_2D<T>>(), get<0>(get<0>(pt)),
                     rot_mat_2D<T>(get<0>(get<1>(pt))), get<1>(get<0>(pt)),
                     get<1>(get<1>(pt)), get<2>(get<0>(pt)), get<2>(get<1>(pt)),
                     vect<T, 2>(0.0, 0.0), 0.0);
}

template <typename T>
frame_2D<T> get_frame_2D(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>>,
                           arithmetic_tuple<T, T>>& pt) {
  return frame_2D<T>(std::weak_ptr<pose_2D<T>>(), get<0>(get<0>(pt)),
                     rot_mat_2D<T>(get<0>(get<1>(pt))), get<1>(get<0>(pt)),
                     get<1>(get<1>(pt)), vect<T, 2>(0.0, 0.0), 0.0,
                     vect<T, 2>(0.0, 0.0), 0.0);
}

template <typename T>
frame_2D<T> get_frame_2D(const arithmetic_tuple<arithmetic_tuple<vect<T, 2>>,
                                                arithmetic_tuple<T>>& pt) {
  return frame_2D<T>(std::weak_ptr<pose_2D<T>>(), get<0>(get<0>(pt)),
                     rot_mat_2D<T>(get<0>(get<1>(pt))), vect<T, 2>(0.0, 0.0),
                     0.0, vect<T, 2>(0.0, 0.0), 0.0, vect<T, 2>(0.0, 0.0), 0.0);
}

template <typename T>
void set_frame_2D(
    arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>, vect<T, 2>>,
                     arithmetic_tuple<T, T, T>>& pt,
    const frame_2D<T>& p) {
  get<0>(get<0>(pt)) = p.Position;
  get<0>(get<1>(pt)) = p.Rotation.getAngle();
  get<1>(get<0>(pt)) = p.Velocity;
  get<1>(get<1>(pt)) = p.AngVelocity;
  get<2>(get<0>(pt)) = p.Acceleration;
  get<2>(get<1>(pt)) = p.AngAcceleration;
}

template <typename T>
void set_frame_2D(arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>>,
                                   arithmetic_tuple<T, T>>& pt,
                  const frame_2D<T>& p) {
  get<0>(get<0>(pt)) = p.Position;
  get<0>(get<1>(pt)) = p.Rotation.getAngle();
  get<1>(get<0>(pt)) = p.Velocity;
  get<1>(get<1>(pt)) = p.AngVelocity;
}

template <typename T>
void set_frame_2D(
    arithmetic_tuple<arithmetic_tuple<vect<T, 2>>, arithmetic_tuple<T>>& pt,
    const frame_2D<T>& p) {
  get<0>(get<0>(pt)) = p.Position;
  get<0>(get<1>(pt)) = p.Rotation.getAngle();
}

template <typename T>
pose_2D<T> get_pose_2D(const arithmetic_tuple<arithmetic_tuple<vect<T, 2>>,
                                              arithmetic_tuple<T>>& pt) {
  return pose_2D<T>(std::weak_ptr<pose_2D<T>>(), get<0>(get<0>(pt)),
                    rot_mat_2D<T>(get<0>(get<1>(pt))));
}

template <typename T>
void set_pose_2D(
    arithmetic_tuple<arithmetic_tuple<vect<T, 2>>, arithmetic_tuple<T>>& pt,
    const pose_2D<T>& p) {
  get<0>(get<0>(pt)) = p.Position;
  get<0>(get<1>(pt)) = p.Rotation.getAngle();
}

template <typename T>
const T& get_rotation(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>, vect<T, 2>>,
                           arithmetic_tuple<T, T, T>>& pt) {
  return get<0>(get<1>(pt));
}

template <typename T>
const T& get_rotation(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>>,
                           arithmetic_tuple<T, T>>& pt) {
  return get<0>(get<1>(pt));
}

template <typename T>
const T& get_rotation(const arithmetic_tuple<arithmetic_tuple<vect<T, 2>>,
                                             arithmetic_tuple<T>>& pt) {
  return get<0>(get<1>(pt));
}

template <typename T>
void set_rotation(
    arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>, vect<T, 2>>,
                     arithmetic_tuple<T, T, T>>& pt,
    const T& q) {
  get<0>(get<1>(pt)) = q;
}

template <typename T>
void set_rotation(arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>>,
                                   arithmetic_tuple<T, T>>& pt,
                  const T& q) {
  get<0>(get<1>(pt)) = q;
}

template <typename T>
void set_rotation(
    arithmetic_tuple<arithmetic_tuple<vect<T, 2>>, arithmetic_tuple<T>>& pt,
    const T& q) {
  get<0>(get<1>(pt)) = q;
}

template <typename T>
const vect<T, 2>& get_position(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>, vect<T, 2>>,
                           arithmetic_tuple<T, T, T>>& pt) {
  return get<0>(get<0>(pt));
}

template <typename T>
const vect<T, 2>& get_position(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>>,
                           arithmetic_tuple<T, T>>& pt) {
  return get<0>(get<0>(pt));
}

template <typename T>
const vect<T, 2>& get_position(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 2>>, arithmetic_tuple<T>>&
        pt) {
  return get<0>(get<0>(pt));
}

template <typename T>
void set_position(
    arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>, vect<T, 2>>,
                     arithmetic_tuple<T, T, T>>& pt,
    const vect<T, 2>& p) {
  get<0>(get<0>(pt)) = p;
}

template <typename T>
void set_position(arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>>,
                                   arithmetic_tuple<T, T>>& pt,
                  const vect<T, 2>& p) {
  get<0>(get<0>(pt)) = p;
}

template <typename T>
void set_position(
    arithmetic_tuple<arithmetic_tuple<vect<T, 2>>, arithmetic_tuple<T>>& pt,
    const vect<T, 2>& p) {
  get<0>(get<0>(pt)) = p;
}

template <typename T>
const T& get_ang_velocity(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>, vect<T, 2>>,
                           arithmetic_tuple<T, T, T>>& pt) {
  return get<1>(get<1>(pt));
}

template <typename T>
const T& get_ang_velocity(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>>,
                           arithmetic_tuple<T, T>>& pt) {
  return get<1>(get<1>(pt));
}

template <typename T>
void set_ang_velocity(
    arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>, vect<T, 2>>,
                     arithmetic_tuple<T, T, T>>& pt,
    const T& p) {
  get<1>(get<1>(pt)) = p;
}

template <typename T>
void set_ang_velocity(arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>>,
                                       arithmetic_tuple<T, T>>& pt,
                      const T& p) {
  get<1>(get<1>(pt)) = p;
}

template <typename T>
const vect<T, 2>& get_velocity(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>, vect<T, 2>>,
                           arithmetic_tuple<T, T, T>>& pt) {
  return get<1>(get<0>(pt));
}

template <typename T>
const vect<T, 2>& get_velocity(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>>,
                           arithmetic_tuple<T, T>>& pt) {
  return get<1>(get<0>(pt));
}

template <typename T>
void set_velocity(
    arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>, vect<T, 2>>,
                     arithmetic_tuple<T, T, T>>& pt,
    const vect<T, 2>& p) {
  get<1>(get<0>(pt)) = p;
}

template <typename T>
void set_velocity(arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>>,
                                   arithmetic_tuple<T, T>>& pt,
                  const vect<T, 2>& p) {
  get<1>(get<0>(pt)) = p;
}

template <typename T>
const T& get_ang_acceleration(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>, vect<T, 2>>,
                           arithmetic_tuple<T, T, T>>& pt) {
  return get<2>(get<1>(pt));
}

template <typename T>
void set_ang_acceleration(
    arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>, vect<T, 2>>,
                     arithmetic_tuple<T, T, T>>& pt,
    const T& p) {
  get<2>(get<1>(pt)) = p;
}

template <typename T>
const vect<T, 2>& get_acceleration(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>, vect<T, 2>>,
                           arithmetic_tuple<T, T, T>>& pt) {
  return get<2>(get<0>(pt));
}

template <typename T>
void set_acceleration(
    arithmetic_tuple<arithmetic_tuple<vect<T, 2>, vect<T, 2>, vect<T, 2>>,
                     arithmetic_tuple<T, T, T>>& pt,
    const vect<T, 2>& p) {
  get<2>(get<0>(pt)) = p;
}

}  // namespace ReaK

#include "reachability_space.hpp"
#include "temporal_space.hpp"
#include "time_poisson_topology.hpp"

#include "joint_space_limits.hpp"

namespace ReaK::pp {

// se2_0th_order_topology
extern template class metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<
            time_topology, arithmetic_tuple<hyperbox_topology<vect<double, 2>>>,
            euclidean_tuple_distance>,
        differentiable_space<time_topology,
                             arithmetic_tuple<line_segment_topology<double>>,
                             euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se2_1st_order_topology
extern template class metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<
            time_topology,
            arithmetic_tuple<hyperbox_topology<vect<double, 2>>,
                             hyperball_topology<vect<double, 2>>>,
            euclidean_tuple_distance>,
        differentiable_space<time_topology,
                             arithmetic_tuple<line_segment_topology<double>,
                                              line_segment_topology<double>>,
                             euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se2_2nd_order_topology
extern template class metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<
            time_topology,
            arithmetic_tuple<hyperbox_topology<vect<double, 2>>,
                             hyperball_topology<vect<double, 2>>,
                             hyperball_topology<vect<double, 2>>>,
            euclidean_tuple_distance>,
        differentiable_space<time_topology,
                             arithmetic_tuple<line_segment_topology<double>,
                                              line_segment_topology<double>,
                                              line_segment_topology<double>>,
                             euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se2_0th_order_rl_topology
extern template class metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<
            time_topology, arithmetic_tuple<hyperbox_topology<vect<double, 2>>>,
            euclidean_tuple_distance>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<line_segment_topology<double>>,
                              euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se2_1st_order_rl_topology
extern template class metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<
            time_topology,
            arithmetic_tuple<hyperbox_topology<vect<double, 2>>,
                             hyperball_topology<vect<double, 2>>>,
            euclidean_tuple_distance>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<line_segment_topology<double>,
                                               line_segment_topology<double>>,
                              euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se2_2nd_order_rl_topology
extern template class metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<
            time_topology,
            arithmetic_tuple<hyperbox_topology<vect<double, 2>>,
                             hyperball_topology<vect<double, 2>>,
                             hyperball_topology<vect<double, 2>>>,
            euclidean_tuple_distance>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<line_segment_topology<double>,
                                               line_segment_topology<double>,
                                               line_segment_topology<double>>,
                              euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se2_0th_order_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            differentiable_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 2>>>,
                euclidean_tuple_distance>,
            differentiable_space<
                time_topology, arithmetic_tuple<line_segment_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;

// se2_1st_order_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            differentiable_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 2>>,
                                 hyperball_topology<vect<double, 2>>>,
                euclidean_tuple_distance>,
            differentiable_space<
                time_topology,
                arithmetic_tuple<line_segment_topology<double>,
                                 line_segment_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;

// se2_2nd_order_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            differentiable_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 2>>,
                                 hyperball_topology<vect<double, 2>>,
                                 hyperball_topology<vect<double, 2>>>,
                euclidean_tuple_distance>,
            differentiable_space<
                time_topology,
                arithmetic_tuple<line_segment_topology<double>,
                                 line_segment_topology<double>,
                                 line_segment_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;

// se2_0th_order_rl_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 2>>>,
                euclidean_tuple_distance>,
            reach_time_diff_space<
                time_topology, arithmetic_tuple<line_segment_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;

// se2_1st_order_rl_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 2>>,
                                 hyperball_topology<vect<double, 2>>>,
                euclidean_tuple_distance>,
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<line_segment_topology<double>,
                                 line_segment_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;

// se2_2nd_order_rl_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 2>>,
                                 hyperball_topology<vect<double, 2>>,
                                 hyperball_topology<vect<double, 2>>>,
                euclidean_tuple_distance>,
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<line_segment_topology<double>,
                                 line_segment_topology<double>,
                                 line_segment_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;

// se2_0th_order_rl_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 2>>>,
                euclidean_tuple_distance>,
            reach_time_diff_space<
                time_topology, arithmetic_tuple<line_segment_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, reach_plus_time_metric>;

// se2_1st_order_rl_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 2>>,
                                 hyperball_topology<vect<double, 2>>>,
                euclidean_tuple_distance>,
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<line_segment_topology<double>,
                                 line_segment_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, reach_plus_time_metric>;

// se2_2nd_order_rl_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 2>>,
                                 hyperball_topology<vect<double, 2>>,
                                 hyperball_topology<vect<double, 2>>>,
                euclidean_tuple_distance>,
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<line_segment_topology<double>,
                                 line_segment_topology<double>,
                                 line_segment_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, reach_plus_time_metric>;

extern template metric_space_array_t<se2_0th_order_rl_topology_t<double>, 1>
joint_limits_mapping<double>::make_rl_joint_space(
    const metric_space_array_t<se2_0th_order_topology_t<double>, 1>&) const;
extern template metric_space_array_t<se2_1st_order_rl_topology_t<double>, 1>
joint_limits_mapping<double>::make_rl_joint_space(
    const metric_space_array_t<se2_1st_order_topology_t<double>, 1>&) const;
extern template metric_space_array_t<se2_2nd_order_rl_topology_t<double>, 1>
joint_limits_mapping<double>::make_rl_joint_space(
    const metric_space_array_t<se2_2nd_order_topology_t<double>, 1>&) const;

extern template metric_space_array_t<se2_0th_order_topology_t<double>, 1>
joint_limits_mapping<double>::make_normal_joint_space(
    const metric_space_array_t<se2_0th_order_rl_topology_t<double>, 1>&) const;
extern template metric_space_array_t<se2_1st_order_topology_t<double>, 1>
joint_limits_mapping<double>::make_normal_joint_space(
    const metric_space_array_t<se2_1st_order_rl_topology_t<double>, 1>&) const;
extern template metric_space_array_t<se2_2nd_order_topology_t<double>, 1>
joint_limits_mapping<double>::make_normal_joint_space(
    const metric_space_array_t<se2_2nd_order_rl_topology_t<double>, 1>&) const;

extern template topology_point_type_t<
    metric_space_array_t<se2_0th_order_rl_topology_t<double>, 1>>
joint_limits_mapping<double>::map_to_space(
    const topology_point_type_t<
        metric_space_array_t<se2_0th_order_topology_t<double>, 1>>& pt,
    const metric_space_array_t<se2_0th_order_topology_t<double>, 1>&,
    const metric_space_array_t<se2_0th_order_rl_topology_t<double>, 1>&) const;
extern template topology_point_type_t<
    metric_space_array_t<se2_1st_order_rl_topology_t<double>, 1>>
joint_limits_mapping<double>::map_to_space(
    const topology_point_type_t<
        metric_space_array_t<se2_1st_order_topology_t<double>, 1>>& pt,
    const metric_space_array_t<se2_1st_order_topology_t<double>, 1>&,
    const metric_space_array_t<se2_1st_order_rl_topology_t<double>, 1>&) const;
extern template topology_point_type_t<
    metric_space_array_t<se2_2nd_order_rl_topology_t<double>, 1>>
joint_limits_mapping<double>::map_to_space(
    const topology_point_type_t<
        metric_space_array_t<se2_2nd_order_topology_t<double>, 1>>& pt,
    const metric_space_array_t<se2_2nd_order_topology_t<double>, 1>&,
    const metric_space_array_t<se2_2nd_order_rl_topology_t<double>, 1>&) const;

extern template topology_point_type_t<
    metric_space_array_t<se2_0th_order_topology_t<double>, 1>>
joint_limits_mapping<double>::map_to_space(
    const topology_point_type_t<
        metric_space_array_t<se2_0th_order_rl_topology_t<double>, 1>>& pt,
    const metric_space_array_t<se2_0th_order_rl_topology_t<double>, 1>&,
    const metric_space_array_t<se2_0th_order_topology_t<double>, 1>&) const;
extern template topology_point_type_t<
    metric_space_array_t<se2_1st_order_topology_t<double>, 1>>
joint_limits_mapping<double>::map_to_space(
    const topology_point_type_t<
        metric_space_array_t<se2_1st_order_rl_topology_t<double>, 1>>& pt,
    const metric_space_array_t<se2_1st_order_rl_topology_t<double>, 1>&,
    const metric_space_array_t<se2_1st_order_topology_t<double>, 1>&) const;
extern template topology_point_type_t<
    metric_space_array_t<se2_2nd_order_topology_t<double>, 1>>
joint_limits_mapping<double>::map_to_space(
    const topology_point_type_t<
        metric_space_array_t<se2_2nd_order_rl_topology_t<double>, 1>>& pt,
    const metric_space_array_t<se2_2nd_order_rl_topology_t<double>, 1>&,
    const metric_space_array_t<se2_2nd_order_topology_t<double>, 1>&) const;

}  // namespace ReaK::pp

#endif

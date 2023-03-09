/**
 * \file se3_topologies.h
 *
 * This library provides classes that define topologies on SE(3) (3D rigid-body motion).
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

#ifndef REAK_TOPOLOGIES_SPACES_SE3_TOPOLOGIES_H_
#define REAK_TOPOLOGIES_SPACES_SE3_TOPOLOGIES_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/topologies/spaces/so3_topologies.h"

#include "ReaK/topologies/spaces/differentiable_space.h"
#include "ReaK/topologies/spaces/metric_space_tuple.h"
#include "ReaK/topologies/spaces/rate_limited_spaces.h"

#include "ReaK/topologies/spaces/hyperball_topology.h"
#include "ReaK/topologies/spaces/hyperbox_topology.h"
#include "ReaK/topologies/spaces/line_topology.h"

#include "ReaK/math/kinetostatics/frame_3D.h"
#include "ReaK/math/lin_alg/arithmetic_tuple.h"
#include "ReaK/math/lin_alg/vect_alg.h"

namespace ReaK {

namespace pp {

/**
 * This meta-function defines the type for a 0th order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se3_0th_order_topology {
  using type = metric_space_tuple<
      arithmetic_tuple<
          differentiable_space<time_topology,
                               arithmetic_tuple<hyperbox_topology<vect<T, 3>>>,
                               DistanceMetric>,
          differentiable_space<time_topology,
                               arithmetic_tuple<quaternion_topology<T>>,
                               DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using se3_0th_order_topology_t =
    typename se3_0th_order_topology<T, DistanceMetric>::type;

template <typename T>
auto make_se3_space(const std::string& aName, const vect<T, 3>& aMinCorner,
                    const vect<T, 3>& aMaxCorner) {
  return se3_0th_order_topology_t<T>{arithmetic_tuple(
      differentiable_space<time_topology,
                           arithmetic_tuple<hyperbox_topology<vect<T, 3>>>,
                           euclidean_tuple_distance>(
          arithmetic_tuple(hyperbox_topology<vect<T, 3>>(
              aName + "_pos", aMinCorner, aMaxCorner))),
      differentiable_space<time_topology,
                           arithmetic_tuple<quaternion_topology<T>>,
                           euclidean_tuple_distance>(
          arithmetic_tuple(quaternion_topology<T>(aName + "_quat"))))};
}

template <typename TupleDistanceMetric, typename T>
auto make_se3_space(const std::string& aName, const vect<T, 3>& aMinCorner,
                    const vect<T, 3>& aMaxCorner) {

  return se3_0th_order_topology_t<T, TupleDistanceMetric>{arithmetic_tuple(
      differentiable_space<time_topology,
                           arithmetic_tuple<hyperbox_topology<vect<T, 3>>>,
                           TupleDistanceMetric>(
          arithmetic_tuple(hyperbox_topology<vect<T, 3>>(
              aName + "_pos", aMinCorner, aMaxCorner))),
      differentiable_space<time_topology,
                           arithmetic_tuple<quaternion_topology<T>>,
                           TupleDistanceMetric>(
          arithmetic_tuple(quaternion_topology<T>(aName + "_quat"))))};
}

/**
 * This meta-function defines the type for a 1st order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se3_1st_order_topology {
  using type = metric_space_tuple<
      arithmetic_tuple<
          differentiable_space<time_topology,
                               arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                                hyperball_topology<vect<T, 3>>>,
                               DistanceMetric>,
          differentiable_space<time_topology,
                               arithmetic_tuple<quaternion_topology<T>,
                                                ang_velocity_3D_topology<T>>,
                               DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using se3_1st_order_topology_t =
    typename se3_1st_order_topology<T, DistanceMetric>::type;

template <typename T>
auto make_se3_space(const std::string& aName, const vect<T, 3>& aMinCorner,
                    const vect<T, 3>& aMaxCorner, const T& aMaxSpeed,
                    const T& aMaxAngularSpeed) {

  return se3_1st_order_topology_t<T>{arithmetic_tuple(
      differentiable_space<time_topology,
                           arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                            hyperball_topology<vect<T, 3>>>,
                           euclidean_tuple_distance>(arithmetic_tuple(
          hyperbox_topology<vect<T, 3>>(aName + "_pos", aMinCorner, aMaxCorner),
          hyperball_topology<vect<T, 3>>(
              aName + "_vel", vect<T, 3>(0.0, 0.0, 0.0), aMaxSpeed))),
      differentiable_space<
          time_topology,
          arithmetic_tuple<quaternion_topology<T>, ang_velocity_3D_topology<T>>,
          euclidean_tuple_distance>(arithmetic_tuple(
          quaternion_topology<T>(aName + "_quat"),
          ang_velocity_3D_topology<T>(aName + "_ang_vel", aMaxAngularSpeed))))};
}

template <typename TupleDistanceMetric, typename T>
auto make_se3_space(const std::string& aName, const vect<T, 3>& aMinCorner,
                    const vect<T, 3>& aMaxCorner, const T& aMaxSpeed,
                    const T& aMaxAngularSpeed) {

  return se3_1st_order_topology_t<T, TupleDistanceMetric>{arithmetic_tuple(
      differentiable_space<time_topology,
                           arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                            hyperball_topology<vect<T, 3>>>,
                           TupleDistanceMetric>(arithmetic_tuple(
          hyperbox_topology<vect<T, 3>>(aName + "_pos", aMinCorner, aMaxCorner),
          hyperball_topology<vect<T, 3>>(
              aName + "_vel", vect<T, 3>(0.0, 0.0, 0.0), aMaxSpeed))),
      differentiable_space<
          time_topology,
          arithmetic_tuple<quaternion_topology<T>, ang_velocity_3D_topology<T>>,
          TupleDistanceMetric>(arithmetic_tuple(
          quaternion_topology<T>(aName + "_quat"),
          ang_velocity_3D_topology<T>(aName + "_ang_vel", aMaxAngularSpeed))))};
}

/**
 * This meta-function defines the type for a 2nd order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se3_2nd_order_topology {
  using type = metric_space_tuple<
      arithmetic_tuple<
          differentiable_space<time_topology,
                               arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                                hyperball_topology<vect<T, 3>>,
                                                hyperball_topology<vect<T, 3>>>,
                               DistanceMetric>,
          differentiable_space<time_topology,
                               arithmetic_tuple<quaternion_topology<T>,
                                                ang_velocity_3D_topology<T>,
                                                ang_accel_3D_topology<T>>,
                               DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using se3_2nd_order_topology_t =
    typename se3_2nd_order_topology<T, DistanceMetric>::type;

template <typename T>
auto make_se3_space(const std::string& aName, const vect<T, 3>& aMinCorner,
                    const vect<T, 3>& aMaxCorner, const T& aMaxSpeed,
                    const T& aMaxAngularSpeed, const T& aMaxAcceleration,
                    const T& aMaxAngularAccel) {

  return se3_2nd_order_topology_t<T>{arithmetic_tuple(
      differentiable_space<time_topology,
                           arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                            hyperball_topology<vect<T, 3>>,
                                            hyperball_topology<vect<T, 3>>>,
                           euclidean_tuple_distance>(arithmetic_tuple(
          hyperbox_topology<vect<T, 3>>(aName + "_pos", aMinCorner, aMaxCorner),
          hyperball_topology<vect<T, 3>>(aName + "_vel",
                                         vect<T, 3>(0.0, 0.0, 0.0), aMaxSpeed),
          hyperball_topology<vect<T, 3>>(
              aName + "_acc", vect<T, 3>(0.0, 0.0, 0.0), aMaxAcceleration))),
      differentiable_space<
          time_topology,
          arithmetic_tuple<quaternion_topology<T>, ang_velocity_3D_topology<T>,
                           ang_accel_3D_topology<T>>,
          euclidean_tuple_distance>(arithmetic_tuple(
          quaternion_topology<T>(aName + "_quat"),
          ang_velocity_3D_topology<T>(aName + "_ang_vel", aMaxAngularSpeed),
          ang_accel_3D_topology<T>(aName + "_ang_acc", aMaxAngularAccel))))};
}

template <typename TupleDistanceMetric, typename T>
auto make_se3_space(const std::string& aName, const vect<T, 3>& aMinCorner,
                    const vect<T, 3>& aMaxCorner, const T& aMaxSpeed,
                    const T& aMaxAngularSpeed, const T& aMaxAcceleration,
                    const T& aMaxAngularAccel) {

  return se3_2nd_order_topology_t<T, TupleDistanceMetric>{arithmetic_tuple(
      differentiable_space<time_topology,
                           arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                            hyperball_topology<vect<T, 3>>,
                                            hyperball_topology<vect<T, 3>>>,
                           TupleDistanceMetric>(arithmetic_tuple(
          hyperbox_topology<vect<T, 3>>(aName + "_pos", aMinCorner, aMaxCorner),
          hyperball_topology<vect<T, 3>>(aName + "_vel",
                                         vect<T, 3>(0.0, 0.0, 0.0), aMaxSpeed),
          hyperball_topology<vect<T, 3>>(
              aName + "_acc", vect<T, 3>(0.0, 0.0, 0.0), aMaxAcceleration))),
      differentiable_space<
          time_topology,
          arithmetic_tuple<quaternion_topology<T>, ang_velocity_3D_topology<T>,
                           ang_accel_3D_topology<T>>,
          TupleDistanceMetric>(arithmetic_tuple(
          quaternion_topology<T>(aName + "_quat"),
          ang_velocity_3D_topology<T>(aName + "_ang_vel", aMaxAngularSpeed),
          ang_accel_3D_topology<T>(aName + "_ang_acc", aMaxAngularAccel))))};
}

template <typename T, int Order,
          typename DistanceMetric = euclidean_tuple_distance>
struct se3_topology {
  using type = std::conditional_t<
      (Order == 0), se3_0th_order_topology_t<T, DistanceMetric>,
      std::conditional_t<(Order == 1),
                         se3_1st_order_topology_t<T, DistanceMetric>,
                         se3_2nd_order_topology_t<T, DistanceMetric>>>;
};

template <typename T, int Order,
          typename DistanceMetric = euclidean_tuple_distance>
using se3_topology_t = typename se3_topology<T, Order, DistanceMetric>::type;

template <typename SE3Space>
struct is_se3_space : std::false_type {};

template <typename SE3Space>
static constexpr bool is_se3_space_v = is_se3_space<SE3Space>::value;

template <typename T, typename DistanceMetric>
struct is_se3_space<metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<time_topology,
                             arithmetic_tuple<hyperbox_topology<vect<T, 3>>>,
                             DistanceMetric>,
        differentiable_space<time_topology,
                             arithmetic_tuple<quaternion_topology<T>>,
                             DistanceMetric>>,
    DistanceMetric>> : std::true_type {};

template <typename T, typename DistanceMetric>
struct is_se3_space<metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<time_topology,
                             arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                              hyperball_topology<vect<T, 3>>>,
                             DistanceMetric>,
        differentiable_space<time_topology,
                             arithmetic_tuple<quaternion_topology<T>,
                                              ang_velocity_3D_topology<T>>,
                             DistanceMetric>>,
    DistanceMetric>> : std::true_type {};

template <typename T, typename DistanceMetric>
struct is_se3_space<metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<time_topology,
                             arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                              hyperball_topology<vect<T, 3>>,
                                              hyperball_topology<vect<T, 3>>>,
                             DistanceMetric>,
        differentiable_space<time_topology,
                             arithmetic_tuple<quaternion_topology<T>,
                                              ang_velocity_3D_topology<T>,
                                              ang_accel_3D_topology<T>>,
                             DistanceMetric>>,
    DistanceMetric>> : std::true_type {};

/**
 * This meta-function defines the type for a 0th order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se3_0th_order_rl_topology {
  using type = metric_space_tuple<
      arithmetic_tuple<
          reach_time_diff_space<time_topology,
                                arithmetic_tuple<hyperbox_topology<vect<T, 3>>>,
                                DistanceMetric>,
          reach_time_diff_space<time_topology,
                                arithmetic_tuple<rate_limited_quat_space<T>>,
                                DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using se3_0th_order_rl_topology_t =
    typename se3_0th_order_rl_topology<T, DistanceMetric>::type;

template <typename T>
auto make_rl_se3_space(const std::string& aName, const vect<T, 3>& aMinCorner,
                       const vect<T, 3>& aMaxCorner, const T& aMaxSpeed,
                       const T& aMaxAngularSpeed) {

  return se3_0th_order_rl_topology_t<T>{arithmetic_tuple(
      reach_time_diff_space<time_topology,
                            arithmetic_tuple<hyperbox_topology<vect<T, 3>>>,
                            euclidean_tuple_distance>(
          arithmetic_tuple(hyperbox_topology<vect<T, 3>>(
              aName + "_pos", aMinCorner * (1.0 / aMaxSpeed),
              aMaxCorner * (1.0 / aMaxSpeed)))),
      reach_time_diff_space<time_topology,
                            arithmetic_tuple<rate_limited_quat_space<T>>,
                            euclidean_tuple_distance>(arithmetic_tuple(
          rate_limited_quat_space<T>(aName + "_quat", aMaxAngularSpeed))))};
}

template <typename TupleDistanceMetric, typename T>
auto make_rl_se3_space(const std::string& aName, const vect<T, 3>& aMinCorner,
                       const vect<T, 3>& aMaxCorner, const T& aMaxSpeed,
                       const T& aMaxAngularSpeed) {

  return se3_0th_order_rl_topology_t<T, TupleDistanceMetric>{arithmetic_tuple(
      reach_time_diff_space<time_topology,
                            arithmetic_tuple<hyperbox_topology<vect<T, 3>>>,
                            TupleDistanceMetric>(
          arithmetic_tuple(hyperbox_topology<vect<T, 3>>(
              aName + "_pos", aMinCorner * (1.0 / aMaxSpeed),
              aMaxCorner * (1.0 / aMaxSpeed)))),
      reach_time_diff_space<time_topology,
                            arithmetic_tuple<rate_limited_quat_space<T>>,
                            TupleDistanceMetric>(arithmetic_tuple(
          rate_limited_quat_space<T>(aName + "_quat", aMaxAngularSpeed))))};
}

/**
 * This meta-function defines the type for a 1st order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se3_1st_order_rl_topology {
  using type = metric_space_tuple<
      arithmetic_tuple<
          reach_time_diff_space<
              time_topology,
              arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                               hyperball_topology<vect<T, 3>>>,
              DistanceMetric>,
          reach_time_diff_space<time_topology,
                                arithmetic_tuple<rate_limited_quat_space<T>,
                                                 ang_velocity_3D_topology<T>>,
                                DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using se3_1st_order_rl_topology_t =
    typename se3_1st_order_rl_topology<T, DistanceMetric>::type;

template <typename T>
auto make_rl_se3_space(const std::string& aName, const vect<T, 3>& aMinCorner,
                       const vect<T, 3>& aMaxCorner, const T& aMaxSpeed,
                       const T& aMaxAngularSpeed, const T& aMaxAcceleration,
                       const T& aMaxAngularAccel) {

  return se3_1st_order_rl_topology_t<T>{arithmetic_tuple(
      reach_time_diff_space<time_topology,
                            arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                             hyperball_topology<vect<T, 3>>>,
                            euclidean_tuple_distance>(
          arithmetic_tuple(hyperbox_topology<vect<T, 3>>(
                               aName + "_pos", aMinCorner * (1.0 / aMaxSpeed),
                               aMaxCorner * (1.0 / aMaxSpeed)),
                           hyperball_topology<vect<T, 3>>(
                               aName + "_vel", vect<T, 3>(0.0, 0.0, 0.0),
                               aMaxSpeed / aMaxAcceleration)),
          euclidean_tuple_distance(),
          arithmetic_tuple(
              reach_time_differentiation(aMaxSpeed / aMaxAcceleration))),
      reach_time_diff_space<time_topology,
                            arithmetic_tuple<rate_limited_quat_space<T>,
                                             ang_velocity_3D_topology<T>>,
                            euclidean_tuple_distance>(
          arithmetic_tuple(
              rate_limited_quat_space<T>(aName + "_quat", aMaxAngularSpeed),
              ang_velocity_3D_topology<T>(aName + "_ang_vel",
                                          aMaxAngularSpeed / aMaxAngularAccel)),
          euclidean_tuple_distance(),
          arithmetic_tuple(reach_time_differentiation(aMaxAngularSpeed /
                                                      aMaxAngularAccel))))};
}

template <typename TupleDistanceMetric, typename T>
auto make_rl_se3_space(const std::string& aName, const vect<T, 3>& aMinCorner,
                       const vect<T, 3>& aMaxCorner, const T& aMaxSpeed,
                       const T& aMaxAngularSpeed, const T& aMaxAcceleration,
                       const T& aMaxAngularAccel) {

  return se3_1st_order_rl_topology_t<T, TupleDistanceMetric>{arithmetic_tuple(
      reach_time_diff_space<time_topology,
                            arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                             hyperball_topology<vect<T, 3>>>,
                            TupleDistanceMetric>(
          arithmetic_tuple(hyperbox_topology<vect<T, 3>>(
                               aName + "_pos", aMinCorner * (1.0 / aMaxSpeed),
                               aMaxCorner * (1.0 / aMaxSpeed)),
                           hyperball_topology<vect<T, 3>>(
                               aName + "_vel", vect<T, 3>(0.0, 0.0, 0.0),
                               aMaxSpeed / aMaxAcceleration)),
          TupleDistanceMetric(),
          arithmetic_tuple(
              reach_time_differentiation(aMaxSpeed / aMaxAcceleration))),
      reach_time_diff_space<time_topology,
                            arithmetic_tuple<rate_limited_quat_space<T>,
                                             ang_velocity_3D_topology<T>>,
                            TupleDistanceMetric>(
          arithmetic_tuple(
              rate_limited_quat_space<T>(aName + "_quat", aMaxAngularSpeed),
              ang_velocity_3D_topology<T>(aName + "_ang_vel",
                                          aMaxAngularSpeed / aMaxAngularAccel)),
          TupleDistanceMetric(),
          arithmetic_tuple(reach_time_differentiation(aMaxAngularSpeed /
                                                      aMaxAngularAccel))))};
}

/**
 * This meta-function defines the type for a 2nd order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se3_2nd_order_rl_topology {
  using type = metric_space_tuple<
      arithmetic_tuple<
          reach_time_diff_space<
              time_topology,
              arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                               hyperball_topology<vect<T, 3>>,
                               hyperball_topology<vect<T, 3>>>,
              DistanceMetric>,
          reach_time_diff_space<time_topology,
                                arithmetic_tuple<rate_limited_quat_space<T>,
                                                 ang_velocity_3D_topology<T>,
                                                 ang_accel_3D_topology<T>>,
                                DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using se3_2nd_order_rl_topology_t =
    typename se3_2nd_order_rl_topology<T, DistanceMetric>::type;

template <typename T>
auto make_rl_se3_space(const std::string& aName, const vect<T, 3>& aMinCorner,
                       const vect<T, 3>& aMaxCorner, const T& aMaxSpeed,
                       const T& aMaxAngularSpeed, const T& aMaxAcceleration,
                       const T& aMaxAngularAccel, const T& aMaxJerk,
                       const T& aMaxAngularJerk) {

  return se3_2nd_order_rl_topology_t<T>{arithmetic_tuple(
      reach_time_diff_space<time_topology,
                            arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                             hyperball_topology<vect<T, 3>>,
                                             hyperball_topology<vect<T, 3>>>,
                            euclidean_tuple_distance>(
          arithmetic_tuple(hyperbox_topology<vect<T, 3>>(
                               aName + "_pos", aMinCorner * (1.0 / aMaxSpeed),
                               aMaxCorner * (1.0 / aMaxSpeed)),
                           hyperball_topology<vect<T, 3>>(
                               aName + "_vel", vect<T, 3>(0.0, 0.0, 0.0),
                               aMaxSpeed / aMaxAcceleration),
                           hyperball_topology<vect<T, 3>>(
                               aName + "_acc", vect<T, 3>(0.0, 0.0, 0.0),
                               aMaxAcceleration / aMaxJerk)),
          euclidean_tuple_distance(),
          arithmetic_tuple(
              reach_time_differentiation(aMaxSpeed / aMaxAcceleration),
              reach_time_differentiation(aMaxAcceleration / aMaxJerk))),
      reach_time_diff_space<time_topology,
                            arithmetic_tuple<rate_limited_quat_space<T>,
                                             ang_velocity_3D_topology<T>,
                                             ang_accel_3D_topology<T>>,
                            euclidean_tuple_distance>(
          arithmetic_tuple(
              rate_limited_quat_space<T>(aName + "_quat", aMaxAngularSpeed),
              ang_velocity_3D_topology<T>(aName + "_ang_vel",
                                          aMaxAngularSpeed / aMaxAngularAccel),
              ang_accel_3D_topology<T>(aName + "_ang_acc",
                                       aMaxAngularAccel / aMaxAngularJerk)),
          euclidean_tuple_distance(),
          arithmetic_tuple(
              reach_time_differentiation(aMaxAngularSpeed / aMaxAngularAccel),
              reach_time_differentiation(aMaxAngularAccel /
                                         aMaxAngularJerk))))};
}

template <typename TupleDistanceMetric, typename T>
auto make_rl_se3_space(const std::string& aName, const vect<T, 3>& aMinCorner,
                       const vect<T, 3>& aMaxCorner, const T& aMaxSpeed,
                       const T& aMaxAngularSpeed, const T& aMaxAcceleration,
                       const T& aMaxAngularAccel, const T& aMaxJerk,
                       const T& aMaxAngularJerk) {

  return se3_2nd_order_rl_topology_t<T, TupleDistanceMetric>{arithmetic_tuple(
      reach_time_diff_space<time_topology,
                            arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                             hyperball_topology<vect<T, 3>>,
                                             hyperball_topology<vect<T, 3>>>,
                            TupleDistanceMetric>(
          arithmetic_tuple(hyperbox_topology<vect<T, 3>>(
                               aName + "_pos", aMinCorner * (1.0 / aMaxSpeed),
                               aMaxCorner * (1.0 / aMaxSpeed)),
                           hyperball_topology<vect<T, 3>>(
                               aName + "_vel", vect<T, 3>(0.0, 0.0, 0.0),
                               aMaxSpeed / aMaxAcceleration),
                           hyperball_topology<vect<T, 3>>(
                               aName + "_acc", vect<T, 3>(0.0, 0.0, 0.0),
                               aMaxAcceleration / aMaxJerk)),
          TupleDistanceMetric(),
          arithmetic_tuple(
              reach_time_differentiation(aMaxSpeed / aMaxAcceleration),
              reach_time_differentiation(aMaxAcceleration / aMaxJerk))),
      reach_time_diff_space<time_topology,
                            arithmetic_tuple<rate_limited_quat_space<T>,
                                             ang_velocity_3D_topology<T>,
                                             ang_accel_3D_topology<T>>,
                            TupleDistanceMetric>(
          arithmetic_tuple(
              rate_limited_quat_space<T>(aName + "_quat", aMaxAngularSpeed),
              ang_velocity_3D_topology<T>(aName + "_ang_vel",
                                          aMaxAngularSpeed / aMaxAngularAccel),
              ang_accel_3D_topology<T>(aName + "_ang_acc",
                                       aMaxAngularAccel / aMaxAngularJerk)),
          TupleDistanceMetric(),
          arithmetic_tuple(
              reach_time_differentiation(aMaxAngularSpeed / aMaxAngularAccel),
              reach_time_differentiation(aMaxAngularAccel /
                                         aMaxAngularJerk))))};
}

template <typename T, int Order,
          typename DistanceMetric = euclidean_tuple_distance>
struct se3_rl_topology {
  using type = std::conditional_t<
      (Order == 0), se3_0th_order_rl_topology_t<T, DistanceMetric>,
      std::conditional_t<(Order == 1),
                         se3_1st_order_rl_topology_t<T, DistanceMetric>,
                         se3_2nd_order_rl_topology_t<T, DistanceMetric>>>;
};

template <typename T, int Order,
          typename DistanceMetric = euclidean_tuple_distance>
using se3_rl_topology_t =
    typename se3_rl_topology<T, Order, DistanceMetric>::type;

template <typename SE3Space>
struct is_rate_limited_se3_space : std::false_type {};

template <typename SE3Space>
static constexpr bool is_rate_limited_se3_space_v =
    is_rate_limited_se3_space<SE3Space>::value;

template <typename T, typename DistanceMetric>
struct is_rate_limited_se3_space<metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<hyperbox_topology<vect<T, 3>>>,
                              DistanceMetric>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<rate_limited_quat_space<T>>,
                              DistanceMetric>>,
    DistanceMetric>> : std::true_type {};

template <typename T, typename DistanceMetric>
struct is_rate_limited_se3_space<metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                               hyperball_topology<vect<T, 3>>>,
                              DistanceMetric>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<rate_limited_quat_space<T>,
                                               ang_velocity_3D_topology<T>>,
                              DistanceMetric>>,
    DistanceMetric>> : std::true_type {};

template <typename T, typename DistanceMetric>
struct is_rate_limited_se3_space<metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                               hyperball_topology<vect<T, 3>>,
                                               hyperball_topology<vect<T, 3>>>,
                              DistanceMetric>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<rate_limited_quat_space<T>,
                                               ang_velocity_3D_topology<T>,
                                               ang_accel_3D_topology<T>>,
                              DistanceMetric>>,
    DistanceMetric>> : std::true_type {};

template <typename T, typename DistanceMetric>
struct get_rate_limited_space<metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<time_topology,
                             arithmetic_tuple<hyperbox_topology<vect<T, 3>>>,
                             DistanceMetric>,
        differentiable_space<time_topology,
                             arithmetic_tuple<quaternion_topology<T>>,
                             DistanceMetric>>,
    DistanceMetric>> {
  using type = metric_space_tuple<
      arithmetic_tuple<
          reach_time_diff_space<time_topology,
                                arithmetic_tuple<hyperbox_topology<vect<T, 3>>>,
                                DistanceMetric>,
          reach_time_diff_space<time_topology,
                                arithmetic_tuple<rate_limited_quat_space<T>>,
                                DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric>
struct get_rate_limited_space<metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<time_topology,
                             arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                              hyperball_topology<vect<T, 3>>>,
                             DistanceMetric>,
        differentiable_space<time_topology,
                             arithmetic_tuple<quaternion_topology<T>,
                                              ang_velocity_3D_topology<T>>,
                             DistanceMetric>>,
    DistanceMetric>> {
  using type = metric_space_tuple<
      arithmetic_tuple<
          reach_time_diff_space<
              time_topology,
              arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                               hyperball_topology<vect<T, 3>>>,
              DistanceMetric>,
          reach_time_diff_space<time_topology,
                                arithmetic_tuple<rate_limited_quat_space<T>,
                                                 ang_velocity_3D_topology<T>>,
                                DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric>
struct get_rate_limited_space<metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<time_topology,
                             arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                              hyperball_topology<vect<T, 3>>,
                                              hyperball_topology<vect<T, 3>>>,
                             DistanceMetric>,
        differentiable_space<time_topology,
                             arithmetic_tuple<quaternion_topology<T>,
                                              ang_velocity_3D_topology<T>,
                                              ang_accel_3D_topology<T>>,
                             DistanceMetric>>,
    DistanceMetric>> {
  using type = metric_space_tuple<
      arithmetic_tuple<
          reach_time_diff_space<
              time_topology,
              arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                               hyperball_topology<vect<T, 3>>,
                               hyperball_topology<vect<T, 3>>>,
              DistanceMetric>,
          reach_time_diff_space<time_topology,
                                arithmetic_tuple<rate_limited_quat_space<T>,
                                                 ang_velocity_3D_topology<T>,
                                                 ang_accel_3D_topology<T>>,
                                DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric>
struct get_rate_illimited_space<metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<hyperbox_topology<vect<T, 3>>>,
                              DistanceMetric>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<rate_limited_quat_space<T>>,
                              DistanceMetric>>,
    DistanceMetric>> {
  using type = metric_space_tuple<
      arithmetic_tuple<
          differentiable_space<time_topology,
                               arithmetic_tuple<hyperbox_topology<vect<T, 3>>>,
                               DistanceMetric>,
          differentiable_space<time_topology,
                               arithmetic_tuple<quaternion_topology<T>>,
                               DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric>
struct get_rate_illimited_space<metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                               hyperball_topology<vect<T, 3>>>,
                              DistanceMetric>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<rate_limited_quat_space<T>,
                                               ang_velocity_3D_topology<T>>,
                              DistanceMetric>>,
    DistanceMetric>> {
  using type = metric_space_tuple<
      arithmetic_tuple<
          differentiable_space<time_topology,
                               arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                                hyperball_topology<vect<T, 3>>>,
                               DistanceMetric>,
          differentiable_space<time_topology,
                               arithmetic_tuple<quaternion_topology<T>,
                                                ang_velocity_3D_topology<T>>,
                               DistanceMetric>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric>
struct get_rate_illimited_space<metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                               hyperball_topology<vect<T, 3>>,
                                               hyperball_topology<vect<T, 3>>>,
                              DistanceMetric>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<rate_limited_quat_space<T>,
                                               ang_velocity_3D_topology<T>,
                                               ang_accel_3D_topology<T>>,
                              DistanceMetric>>,
    DistanceMetric>> {
  using type = metric_space_tuple<
      arithmetic_tuple<
          differentiable_space<time_topology,
                               arithmetic_tuple<hyperbox_topology<vect<T, 3>>,
                                                hyperball_topology<vect<T, 3>>,
                                                hyperball_topology<vect<T, 3>>>,
                               DistanceMetric>,
          differentiable_space<time_topology,
                               arithmetic_tuple<quaternion_topology<T>,
                                                ang_velocity_3D_topology<T>,
                                                ang_accel_3D_topology<T>>,
                               DistanceMetric>>,
      DistanceMetric>;
};

}  // namespace pp

// Because of ADL rules, the get functions for the arithmetic-tuple types that represent SE(3) states should be in the
// ReaK namespace.

template <typename T>
frame_3D<T> get_frame_3D(
    const arithmetic_tuple<
        arithmetic_tuple<vect<T, 3>, vect<T, 3>, vect<T, 3>>,
        arithmetic_tuple<unit_quat<T>, vect<T, 3>, vect<T, 3>>>& pt) {
  return frame_3D<T>(std::weak_ptr<pose_3D<T>>(), get<0>(get<0>(pt)),
                     quaternion<T>(get<0>(get<1>(pt))), get<1>(get<0>(pt)),
                     get<1>(get<1>(pt)), get<2>(get<0>(pt)), get<2>(get<1>(pt)),
                     vect<T, 3>(0.0, 0.0, 0.0), vect<T, 3>(0.0, 0.0, 0.0));
};

template <typename T>
frame_3D<T> get_frame_3D(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>>,
                           arithmetic_tuple<unit_quat<T>, vect<T, 3>>>& pt) {
  return frame_3D<T>(std::weak_ptr<pose_3D<T>>(), get<0>(get<0>(pt)),
                     quaternion<T>(get<0>(get<1>(pt))), get<1>(get<0>(pt)),
                     get<1>(get<1>(pt)), vect<T, 3>(0.0, 0.0, 0.0),
                     vect<T, 3>(0.0, 0.0, 0.0), vect<T, 3>(0.0, 0.0, 0.0),
                     vect<T, 3>(0.0, 0.0, 0.0));
};

template <typename T>
frame_3D<T> get_frame_3D(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 3>>,
                           arithmetic_tuple<unit_quat<T>>>& pt) {
  return frame_3D<T>(std::weak_ptr<pose_3D<T>>(), get<0>(get<0>(pt)),
                     quaternion<T>(get<0>(get<1>(pt))),
                     vect<T, 3>(0.0, 0.0, 0.0), vect<T, 3>(0.0, 0.0, 0.0),
                     vect<T, 3>(0.0, 0.0, 0.0), vect<T, 3>(0.0, 0.0, 0.0),
                     vect<T, 3>(0.0, 0.0, 0.0), vect<T, 3>(0.0, 0.0, 0.0));
};

template <typename T>
void set_frame_3D(
    arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>, vect<T, 3>>,
                     arithmetic_tuple<unit_quat<T>, vect<T, 3>, vect<T, 3>>>&
        pt,
    const frame_3D<T>& p) {
  get<0>(get<0>(pt)) = p.Position;
  get<0>(get<1>(pt)) = unit_quat<T>(p.Quat);
  get<1>(get<0>(pt)) = p.Velocity;
  get<1>(get<1>(pt)) = p.AngVelocity;
  get<2>(get<0>(pt)) = p.Acceleration;
  get<2>(get<1>(pt)) = p.AngAcceleration;
};

template <typename T>
void set_frame_3D(
    arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>>,
                     arithmetic_tuple<unit_quat<T>, vect<T, 3>>>& pt,
    const frame_3D<T>& p) {
  get<0>(get<0>(pt)) = p.Position;
  get<0>(get<1>(pt)) = unit_quat<T>(p.Quat);
  get<1>(get<0>(pt)) = p.Velocity;
  get<1>(get<1>(pt)) = p.AngVelocity;
};

template <typename T>
void set_frame_3D(arithmetic_tuple<arithmetic_tuple<vect<T, 3>>,
                                   arithmetic_tuple<unit_quat<T>>>& pt,
                  const frame_3D<T>& p) {
  get<0>(get<0>(pt)) = p.Position;
  get<0>(get<1>(pt)) = unit_quat<T>(p.Quat);
};

template <typename T>
pose_3D<T> get_pose_3D(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 3>>,
                           arithmetic_tuple<unit_quat<T>>>& pt) {
  return pose_3D<T>(std::weak_ptr<pose_3D<T>>(), get<0>(get<0>(pt)),
                    quaternion<T>(get<0>(get<1>(pt))));
};

template <typename T>
void set_pose_3D(arithmetic_tuple<arithmetic_tuple<vect<T, 3>>,
                                  arithmetic_tuple<unit_quat<T>>>& pt,
                 const pose_3D<T>& p) {
  get<0>(get<0>(pt)) = p.Position;
  get<0>(get<1>(pt)) = unit_quat<T>(p.Quat);
};

template <typename T>
const unit_quat<T>& get_quaternion(
    const arithmetic_tuple<
        arithmetic_tuple<vect<T, 3>, vect<T, 3>, vect<T, 3>>,
        arithmetic_tuple<unit_quat<T>, vect<T, 3>, vect<T, 3>>>& pt) {
  return get<0>(get<1>(pt));
};

template <typename T>
const unit_quat<T>& get_quaternion(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>>,
                           arithmetic_tuple<unit_quat<T>, vect<T, 3>>>& pt) {
  return get<0>(get<1>(pt));
};

template <typename T>
const unit_quat<T>& get_quaternion(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 3>>,
                           arithmetic_tuple<unit_quat<T>>>& pt) {
  return get<0>(get<1>(pt));
};

template <typename T>
void set_quaternion(
    arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>, vect<T, 3>>,
                     arithmetic_tuple<unit_quat<T>, vect<T, 3>, vect<T, 3>>>&
        pt,
    const unit_quat<T>& q) {
  get<0>(get<1>(pt)) = q;
};

template <typename T>
void set_quaternion(
    arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>>,
                     arithmetic_tuple<unit_quat<T>, vect<T, 3>>>& pt,
    const unit_quat<T>& q) {
  get<0>(get<1>(pt)) = q;
};

template <typename T>
void set_quaternion(arithmetic_tuple<arithmetic_tuple<vect<T, 3>>,
                                     arithmetic_tuple<unit_quat<T>>>& pt,
                    const unit_quat<T>& q) {
  get<0>(get<1>(pt)) = q;
};

template <typename T>
const vect<T, 3>& get_position(
    const arithmetic_tuple<
        arithmetic_tuple<vect<T, 3>, vect<T, 3>, vect<T, 3>>,
        arithmetic_tuple<unit_quat<T>, vect<T, 3>, vect<T, 3>>>& pt) {
  return get<0>(get<0>(pt));
};

template <typename T>
const vect<T, 3>& get_position(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>>,
                           arithmetic_tuple<unit_quat<T>, vect<T, 3>>>& pt) {
  return get<0>(get<0>(pt));
};

template <typename T>
const vect<T, 3>& get_position(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 3>>,
                           arithmetic_tuple<unit_quat<T>>>& pt) {
  return get<0>(get<0>(pt));
};

template <typename T>
void set_position(
    arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>, vect<T, 3>>,
                     arithmetic_tuple<unit_quat<T>, vect<T, 3>, vect<T, 3>>>&
        pt,
    const vect<T, 3>& p) {
  get<0>(get<0>(pt)) = p;
};

template <typename T>
void set_position(
    arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>>,
                     arithmetic_tuple<unit_quat<T>, vect<T, 3>>>& pt,
    const vect<T, 3>& p) {
  get<0>(get<0>(pt)) = p;
};

template <typename T>
void set_position(arithmetic_tuple<arithmetic_tuple<vect<T, 3>>,
                                   arithmetic_tuple<unit_quat<T>>>& pt,
                  const vect<T, 3>& p) {
  get<0>(get<0>(pt)) = p;
};

template <typename T>
const vect<T, 3>& get_ang_velocity(
    const arithmetic_tuple<
        arithmetic_tuple<vect<T, 3>, vect<T, 3>, vect<T, 3>>,
        arithmetic_tuple<unit_quat<T>, vect<T, 3>, vect<T, 3>>>& pt) {
  return get<1>(get<1>(pt));
};

template <typename T>
const vect<T, 3>& get_ang_velocity(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>>,
                           arithmetic_tuple<unit_quat<T>, vect<T, 3>>>& pt) {
  return get<1>(get<1>(pt));
};

template <typename T>
void set_ang_velocity(
    arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>, vect<T, 3>>,
                     arithmetic_tuple<unit_quat<T>, vect<T, 3>, vect<T, 3>>>&
        pt,
    const vect<T, 3>& p) {
  get<1>(get<1>(pt)) = p;
};

template <typename T>
void set_ang_velocity(
    arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>>,
                     arithmetic_tuple<unit_quat<T>, vect<T, 3>>>& pt,
    const vect<T, 3>& p) {
  get<1>(get<1>(pt)) = p;
};

template <typename T>
const vect<T, 3>& get_velocity(
    const arithmetic_tuple<
        arithmetic_tuple<vect<T, 3>, vect<T, 3>, vect<T, 3>>,
        arithmetic_tuple<unit_quat<T>, vect<T, 3>, vect<T, 3>>>& pt) {
  return get<1>(get<0>(pt));
};

template <typename T>
const vect<T, 3>& get_velocity(
    const arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>>,
                           arithmetic_tuple<unit_quat<T>, vect<T, 3>>>& pt) {
  return get<1>(get<0>(pt));
};

template <typename T>
void set_velocity(
    arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>, vect<T, 3>>,
                     arithmetic_tuple<unit_quat<T>, vect<T, 3>, vect<T, 3>>>&
        pt,
    const vect<T, 3>& p) {
  get<1>(get<0>(pt)) = p;
};

template <typename T>
void set_velocity(
    arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>>,
                     arithmetic_tuple<unit_quat<T>, vect<T, 3>>>& pt,
    const vect<T, 3>& p) {
  get<1>(get<0>(pt)) = p;
};

template <typename T>
const vect<T, 3>& get_ang_acceleration(
    const arithmetic_tuple<
        arithmetic_tuple<vect<T, 3>, vect<T, 3>, vect<T, 3>>,
        arithmetic_tuple<unit_quat<T>, vect<T, 3>, vect<T, 3>>>& pt) {
  return get<2>(get<1>(pt));
};

template <typename T>
void set_ang_acceleration(
    arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>, vect<T, 3>>,
                     arithmetic_tuple<unit_quat<T>, vect<T, 3>, vect<T, 3>>>&
        pt,
    const vect<T, 3>& p) {
  get<2>(get<1>(pt)) = p;
};

template <typename T>
const vect<T, 3>& get_acceleration(
    const arithmetic_tuple<
        arithmetic_tuple<vect<T, 3>, vect<T, 3>, vect<T, 3>>,
        arithmetic_tuple<unit_quat<T>, vect<T, 3>, vect<T, 3>>>& pt) {
  return get<2>(get<0>(pt));
};

template <typename T>
void set_acceleration(
    arithmetic_tuple<arithmetic_tuple<vect<T, 3>, vect<T, 3>, vect<T, 3>>,
                     arithmetic_tuple<unit_quat<T>, vect<T, 3>, vect<T, 3>>>&
        pt,
    const vect<T, 3>& p) {
  get<2>(get<0>(pt)) = p;
};

}  // namespace ReaK

#include "ReaK/topologies/spaces/reachability_space.h"
#include "ReaK/topologies/spaces/temporal_space.h"
#include "ReaK/topologies/spaces/time_poisson_topology.h"

#include "ReaK/topologies/spaces/joint_space_limits.h"

namespace ReaK::pp {

// se3_0th_order_topology
extern template class metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<
            time_topology, arithmetic_tuple<hyperbox_topology<vect<double, 3>>>,
            euclidean_tuple_distance>,
        differentiable_space<time_topology,
                             arithmetic_tuple<quaternion_topology<double>>,
                             euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se3_1st_order_topology
extern template class metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<
            time_topology,
            arithmetic_tuple<hyperbox_topology<vect<double, 3>>,
                             hyperball_topology<vect<double, 3>>>,
            euclidean_tuple_distance>,
        differentiable_space<time_topology,
                             arithmetic_tuple<quaternion_topology<double>,
                                              ang_velocity_3D_topology<double>>,
                             euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se3_2nd_order_topology
extern template class metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<
            time_topology,
            arithmetic_tuple<hyperbox_topology<vect<double, 3>>,
                             hyperball_topology<vect<double, 3>>,
                             hyperball_topology<vect<double, 3>>>,
            euclidean_tuple_distance>,
        differentiable_space<time_topology,
                             arithmetic_tuple<quaternion_topology<double>,
                                              ang_velocity_3D_topology<double>,
                                              ang_accel_3D_topology<double>>,
                             euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se3_0th_order_rl_topology
extern template class metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<
            time_topology, arithmetic_tuple<hyperbox_topology<vect<double, 3>>>,
            euclidean_tuple_distance>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<rate_limited_quat_space<double>>,
                              euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se3_1st_order_rl_topology
extern template class metric_space_tuple<
    arithmetic_tuple<reach_time_diff_space<
                         time_topology,
                         arithmetic_tuple<hyperbox_topology<vect<double, 3>>,
                                          hyperball_topology<vect<double, 3>>>,
                         euclidean_tuple_distance>,
                     reach_time_diff_space<
                         time_topology,
                         arithmetic_tuple<rate_limited_quat_space<double>,
                                          ang_velocity_3D_topology<double>>,
                         euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se3_2nd_order_rl_topology
extern template class metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<
            time_topology,
            arithmetic_tuple<hyperbox_topology<vect<double, 3>>,
                             hyperball_topology<vect<double, 3>>,
                             hyperball_topology<vect<double, 3>>>,
            euclidean_tuple_distance>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<rate_limited_quat_space<double>,
                                               ang_velocity_3D_topology<double>,
                                               ang_accel_3D_topology<double>>,
                              euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se3_0th_order_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            differentiable_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 3>>>,
                euclidean_tuple_distance>,
            differentiable_space<time_topology,
                                 arithmetic_tuple<quaternion_topology<double>>,
                                 euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;

// se3_1st_order_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            differentiable_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 3>>,
                                 hyperball_topology<vect<double, 3>>>,
                euclidean_tuple_distance>,
            differentiable_space<
                time_topology,
                arithmetic_tuple<quaternion_topology<double>,
                                 ang_velocity_3D_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;

// se3_2nd_order_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            differentiable_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 3>>,
                                 hyperball_topology<vect<double, 3>>,
                                 hyperball_topology<vect<double, 3>>>,
                euclidean_tuple_distance>,
            differentiable_space<
                time_topology,
                arithmetic_tuple<quaternion_topology<double>,
                                 ang_velocity_3D_topology<double>,
                                 ang_accel_3D_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;

// se3_0th_order_rl_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 3>>>,
                euclidean_tuple_distance>,
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<rate_limited_quat_space<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;

// se3_1st_order_rl_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 3>>,
                                 hyperball_topology<vect<double, 3>>>,
                euclidean_tuple_distance>,
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<rate_limited_quat_space<double>,
                                 ang_velocity_3D_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;

// se3_2nd_order_rl_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 3>>,
                                 hyperball_topology<vect<double, 3>>,
                                 hyperball_topology<vect<double, 3>>>,
                euclidean_tuple_distance>,
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<rate_limited_quat_space<double>,
                                 ang_velocity_3D_topology<double>,
                                 ang_accel_3D_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;

// se3_0th_order_rl_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 3>>>,
                euclidean_tuple_distance>,
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<rate_limited_quat_space<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, reach_plus_time_metric>;

// se3_1st_order_rl_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 3>>,
                                 hyperball_topology<vect<double, 3>>>,
                euclidean_tuple_distance>,
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<rate_limited_quat_space<double>,
                                 ang_velocity_3D_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, reach_plus_time_metric>;

// se3_2nd_order_rl_topology
extern template class temporal_space<
    metric_space_tuple<
        arithmetic_tuple<
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<hyperbox_topology<vect<double, 3>>,
                                 hyperball_topology<vect<double, 3>>,
                                 hyperball_topology<vect<double, 3>>>,
                euclidean_tuple_distance>,
            reach_time_diff_space<
                time_topology,
                arithmetic_tuple<rate_limited_quat_space<double>,
                                 ang_velocity_3D_topology<double>,
                                 ang_accel_3D_topology<double>>,
                euclidean_tuple_distance>>,
        euclidean_tuple_distance>,
    time_poisson_topology, reach_plus_time_metric>;

extern template topology_point_type_t<
    metric_space_array_t<se3_0th_order_rl_topology_t<double>, 1>>
joint_limits_mapping<double>::map_to_space(
    const topology_point_type_t<
        metric_space_array_t<se3_0th_order_topology_t<double>, 1>>& pt,
    const metric_space_array_t<se3_0th_order_topology_t<double>, 1>&,
    const metric_space_array_t<se3_0th_order_rl_topology_t<double>, 1>&) const;
extern template topology_point_type_t<
    metric_space_array_t<se3_1st_order_rl_topology_t<double>, 1>>
joint_limits_mapping<double>::map_to_space(
    const topology_point_type_t<
        metric_space_array_t<se3_1st_order_topology_t<double>, 1>>& pt,
    const metric_space_array_t<se3_1st_order_topology_t<double>, 1>&,
    const metric_space_array_t<se3_1st_order_rl_topology_t<double>, 1>&) const;
extern template topology_point_type_t<
    metric_space_array_t<se3_2nd_order_rl_topology_t<double>, 1>>
joint_limits_mapping<double>::map_to_space(
    const topology_point_type_t<
        metric_space_array_t<se3_2nd_order_topology_t<double>, 1>>& pt,
    const metric_space_array_t<se3_2nd_order_topology_t<double>, 1>&,
    const metric_space_array_t<se3_2nd_order_rl_topology_t<double>, 1>&) const;

extern template topology_point_type_t<
    metric_space_array_t<se3_0th_order_topology_t<double>, 1>>
joint_limits_mapping<double>::map_to_space(
    const topology_point_type_t<
        metric_space_array_t<se3_0th_order_rl_topology_t<double>, 1>>& pt,
    const metric_space_array_t<se3_0th_order_rl_topology_t<double>, 1>&,
    const metric_space_array_t<se3_0th_order_topology_t<double>, 1>&) const;
extern template topology_point_type_t<
    metric_space_array_t<se3_1st_order_topology_t<double>, 1>>
joint_limits_mapping<double>::map_to_space(
    const topology_point_type_t<
        metric_space_array_t<se3_1st_order_rl_topology_t<double>, 1>>& pt,
    const metric_space_array_t<se3_1st_order_rl_topology_t<double>, 1>&,
    const metric_space_array_t<se3_1st_order_topology_t<double>, 1>&) const;
extern template topology_point_type_t<
    metric_space_array_t<se3_2nd_order_topology_t<double>, 1>>
joint_limits_mapping<double>::map_to_space(
    const topology_point_type_t<
        metric_space_array_t<se3_2nd_order_rl_topology_t<double>, 1>>& pt,
    const metric_space_array_t<se3_2nd_order_rl_topology_t<double>, 1>&,
    const metric_space_array_t<se3_2nd_order_topology_t<double>, 1>&) const;

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_SE3_TOPOLOGIES_H_

/**
 * \file joint_space_topologies.hpp
 *
 * This library provides classes that define topologies in joint-space (generalized coordinates).
 * All the topologies included here are for a single joint. To allow for more joints, use the
 * meta-function metric_space_array to generate a metric_space_tuple (or form the metric_space_tuple
 * manually).
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

#ifndef REAK_TOPOLOGIES_SPACES_JOINT_SPACE_TOPOLOGIES_H_
#define REAK_TOPOLOGIES_SPACES_JOINT_SPACE_TOPOLOGIES_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/topologies/spaces/differentiable_space.h"
#include "ReaK/topologies/spaces/metric_space_tuple.h"
#include "ReaK/topologies/spaces/rate_limited_spaces.h"

#include "ReaK/topologies/spaces/line_topology.h"

#include "ReaK/math/kinetostatics/gen_coord.h"
#include "ReaK/math/lin_alg/arithmetic_tuple.h"

namespace ReaK {

namespace pp {

/**
 * This meta-function defines the type for a 0th order single-joint space (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct joint_space_0th_order {
  using type = differentiable_space<time_topology,
                                    arithmetic_tuple<line_segment_topology<T>>,
                                    DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using joint_space_0th_order_t =
    typename joint_space_0th_order<T, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 0th order multi-joint space (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam N The number of degrees of freedom of the joint space.
 * \tparam DistanceMetric The distance metric to apply to the tuple (top-level tuple).
 */
template <typename T, std::size_t N,
          typename DistanceMetric = euclidean_tuple_distance>
struct Ndof_0th_order_space {
  using type =
      metric_space_array_t<joint_space_0th_order_t<T>, N, DistanceMetric>;
};

template <typename T, std::size_t N,
          typename DistanceMetric = euclidean_tuple_distance>
using Ndof_0th_order_space_t =
    typename Ndof_0th_order_space<T, N, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 1st order single-joint space.
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct joint_space_1st_order {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<line_segment_topology<T>, line_segment_topology<T>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using joint_space_1st_order_t =
    typename joint_space_1st_order<T, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 1st order multi-joint space (a once-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam N The number of degrees of freedom of the joint space.
 * \tparam DistanceMetric The distance metric to apply to the tuple (top-level tuple).
 */
template <typename T, std::size_t N,
          typename DistanceMetric = euclidean_tuple_distance>
struct Ndof_1st_order_space {
  using type =
      metric_space_array_t<joint_space_1st_order_t<T>, N, DistanceMetric>;
};

template <typename T, std::size_t N,
          typename DistanceMetric = euclidean_tuple_distance>
using Ndof_1st_order_space_t =
    typename Ndof_1st_order_space<T, N, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 2nd order single-joint space.
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct joint_space_2nd_order {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<line_segment_topology<T>, line_segment_topology<T>,
                       line_segment_topology<T>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using joint_space_2nd_order_t =
    typename joint_space_2nd_order<T, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 2nd order multi-joint space (a once-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam N The number of degrees of freedom of the joint space.
 * \tparam DistanceMetric The distance metric to apply to the tuple (top-level tuple).
 */
template <typename T, std::size_t N,
          typename DistanceMetric = euclidean_tuple_distance>
struct Ndof_2nd_order_space {
  using type =
      metric_space_array_t<joint_space_2nd_order_t<T>, N, DistanceMetric>;
};

template <typename T, std::size_t N,
          typename DistanceMetric = euclidean_tuple_distance>
using Ndof_2nd_order_space_t =
    typename Ndof_2nd_order_space<T, N, DistanceMetric>::type;

template <typename JointSpace>
struct is_normal_joint_space : std::false_type {};

template <typename JointSpace>
static constexpr bool is_normal_joint_space_v =
    is_normal_joint_space<JointSpace>::value;

template <typename T, typename DistanceMetric>
struct is_normal_joint_space<differentiable_space<
    time_topology, arithmetic_tuple<line_segment_topology<T>>, DistanceMetric>>
    : std::true_type {};

template <typename T, typename DistanceMetric>
struct is_normal_joint_space<differentiable_space<
    time_topology,
    arithmetic_tuple<line_segment_topology<T>, line_segment_topology<T>>,
    DistanceMetric>> : std::true_type {};

template <typename T, typename DistanceMetric>
struct is_normal_joint_space<differentiable_space<
    time_topology,
    arithmetic_tuple<line_segment_topology<T>, line_segment_topology<T>,
                     line_segment_topology<T>>,
    DistanceMetric>> : std::true_type {};

/**
 * This meta-function defines the type for a rate-limited 0th order single-joint space (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct rl_joint_space_0th_order {
  using type = reach_time_diff_space<time_topology,
                                     arithmetic_tuple<line_segment_topology<T>>,
                                     DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using rl_joint_space_0th_order_t =
    typename rl_joint_space_0th_order<T, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 0th order multi-joint rate-limited space (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam N The number of degrees of freedom of the joint space.
 * \tparam DistanceMetric The distance metric to apply to the tuple (top-level tuple).
 */
template <typename T, std::size_t N,
          typename DistanceMetric = euclidean_tuple_distance>
struct Ndof_0th_order_rl_space {
  using type =
      metric_space_array_t<rl_joint_space_0th_order_t<T>, N, DistanceMetric>;
};

template <typename T, std::size_t N,
          typename DistanceMetric = euclidean_tuple_distance>
using Ndof_0th_order_rl_space_t =
    typename Ndof_0th_order_rl_space<T, N, DistanceMetric>::type;

/**
 * This meta-function defines the type for a rate-limited 1st order single-joint space.
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct rl_joint_space_1st_order {
  using type = reach_time_diff_space<
      time_topology,
      arithmetic_tuple<line_segment_topology<T>, line_segment_topology<T>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using rl_joint_space_1st_order_t =
    typename rl_joint_space_1st_order<T, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 1st order multi-joint rate-limited space (a once-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam N The number of degrees of freedom of the joint space.
 * \tparam DistanceMetric The distance metric to apply to the tuple (top-level tuple).
 */
template <typename T, std::size_t N,
          typename DistanceMetric = euclidean_tuple_distance>
struct Ndof_1st_order_rl_space {
  using type =
      metric_space_array_t<rl_joint_space_1st_order_t<T>, N, DistanceMetric>;
};

template <typename T, std::size_t N,
          typename DistanceMetric = euclidean_tuple_distance>
using Ndof_1st_order_rl_space_t =
    typename Ndof_1st_order_rl_space<T, N, DistanceMetric>::type;

/**
 * This meta-function defines the type for a rate-limited 2nd order single-joint space.
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct rl_joint_space_2nd_order {
  using type = reach_time_diff_space<
      time_topology,
      arithmetic_tuple<line_segment_topology<T>, line_segment_topology<T>,
                       line_segment_topology<T>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using rl_joint_space_2nd_order_t =
    typename rl_joint_space_2nd_order<T, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 2nd order multi-joint rate-limited space (a twice-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam N The number of degrees of freedom of the joint space.
 * \tparam DistanceMetric The distance metric to apply to the tuple (top-level tuple).
 */
template <typename T, std::size_t N,
          typename DistanceMetric = euclidean_tuple_distance>
struct Ndof_2nd_order_rl_space {
  using type =
      metric_space_array_t<rl_joint_space_2nd_order_t<T>, N, DistanceMetric>;
};

template <typename T, std::size_t N,
          typename DistanceMetric = euclidean_tuple_distance>
using Ndof_2nd_order_rl_space_t =
    typename Ndof_2nd_order_rl_space<T, N, DistanceMetric>::type;

template <typename JointSpace>
struct is_rate_limited_joint_space : std::false_type {};

template <typename JointSpace>
static constexpr bool is_rate_limited_joint_space_v =
    is_rate_limited_joint_space<JointSpace>::value;

template <typename T, typename DistanceMetric>
struct is_rate_limited_joint_space<reach_time_diff_space<
    time_topology, arithmetic_tuple<line_segment_topology<T>>, DistanceMetric>>
    : std::true_type {};

template <typename T, typename DistanceMetric>
struct is_rate_limited_joint_space<reach_time_diff_space<
    time_topology,
    arithmetic_tuple<line_segment_topology<T>, line_segment_topology<T>>,
    DistanceMetric>> : std::true_type {};

template <typename T, typename DistanceMetric>
struct is_rate_limited_joint_space<reach_time_diff_space<
    time_topology,
    arithmetic_tuple<line_segment_topology<T>, line_segment_topology<T>,
                     line_segment_topology<T>>,
    DistanceMetric>> : std::true_type {};

}  // namespace pp

template <typename T>
gen_coord<T> get_gen_coord(const arithmetic_tuple<T, T, T>& pt) {
  return gen_coord<T>(get<0>(pt), get<1>(pt), get<2>(pt), 0.0);
}

template <typename T>
gen_coord<T> get_gen_coord(const arithmetic_tuple<T, T>& pt) {
  return gen_coord<T>(get<0>(pt), get<1>(pt), 0.0, 0.0);
}

template <typename T>
gen_coord<T> get_gen_coord(const arithmetic_tuple<T>& pt) {
  return gen_coord<T>(get<0>(pt), 0.0, 0.0, 0.0);
}

template <typename T>
void set_gen_coord(arithmetic_tuple<T, T, T>& pt, const gen_coord<T>& p) {
  get<0>(pt) = p.q;
  get<1>(pt) = p.q_dot;
  get<2>(pt) = p.q_ddot;
}

template <typename T>
void set_gen_coord(arithmetic_tuple<T, T>& pt, const gen_coord<T>& p) {
  get<0>(pt) = p.q;
  get<1>(pt) = p.q_dot;
}

template <typename T>
void set_gen_coord(arithmetic_tuple<T>& pt, const gen_coord<T>& p) {
  get<0>(pt) = p.q;
}

template <typename T>
const T& get_position(const arithmetic_tuple<T, T, T>& pt) {
  return get<0>(pt);
}

template <typename T>
const T& get_position(const arithmetic_tuple<T, T>& pt) {
  return get<0>(pt);
}

template <typename T>
const T& get_position(const arithmetic_tuple<T>& pt) {
  return get<0>(pt);
}

template <typename T>
void set_position(arithmetic_tuple<T, T, T>& pt, const T& p) {
  get<0>(pt) = p;
}

template <typename T>
void set_position(arithmetic_tuple<T, T>& pt, const T& p) {
  get<0>(pt) = p;
}

template <typename T>
void set_position(arithmetic_tuple<T>& pt, const T& p) {
  get<0>(pt) = p;
}

template <typename T>
const T& get_velocity(const arithmetic_tuple<T, T, T>& pt) {
  return get<1>(pt);
}

template <typename T>
const T& get_velocity(const arithmetic_tuple<T, T>& pt) {
  return get<1>(pt);
}

template <typename T>
void set_velocity(arithmetic_tuple<T, T, T>& pt, const T& p) {
  get<1>(pt) = p;
}

template <typename T>
void set_velocity(arithmetic_tuple<T, T>& pt, const T& p) {
  get<1>(pt) = p;
}

template <typename T>
const T& get_acceleration(const arithmetic_tuple<T, T, T>& pt) {
  return get<2>(pt);
}

template <typename T>
void set_acceleration(arithmetic_tuple<T, T, T>& pt, const T& p) {
  get<2>(pt) = p;
}

}  // namespace ReaK

#include "ReaK/topologies/spaces/joint_space_topologies_ext.h"

#endif  // REAK_TOPOLOGIES_SPACES_JOINT_SPACE_TOPOLOGIES_H_

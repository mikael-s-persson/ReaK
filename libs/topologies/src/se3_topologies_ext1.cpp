
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

#include <ReaK/core/base/defs.hpp>

#include <ReaK/topologies/spaces/joint_space_limits.tpp>
#include <ReaK/topologies/spaces/se3_topologies.hpp>

namespace ReaK::pp {

// se3_0th_order_topology
template class metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<
            time_topology, arithmetic_tuple<hyperbox_topology<vect<double, 3>>>,
            euclidean_tuple_distance>,
        differentiable_space<time_topology,
                             arithmetic_tuple<quaternion_topology<double>>,
                             euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se3_1st_order_topology
template class metric_space_tuple<
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
template class metric_space_tuple<
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
template class metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<
            time_topology, arithmetic_tuple<hyperbox_topology<vect<double, 3>>>,
            euclidean_tuple_distance>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<rate_limited_quat_space<double>>,
                              euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se3_1st_order_rl_topology
template class metric_space_tuple<
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
template class metric_space_tuple<
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
template class temporal_space<
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
template class temporal_space<
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
template class temporal_space<
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
template class temporal_space<
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
template class temporal_space<
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
template class temporal_space<
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
template class temporal_space<
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
template class temporal_space<
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
template class temporal_space<
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

template topology_traits<
    metric_space_array_t<se3_0th_order_rl_topology_t<double>, 1>>::point_type
joint_limits_mapping<double>::map_to_space(
    const topology_traits<metric_space_array_t<se3_0th_order_topology_t<double>,
                                               1>>::point_type& pt,
    const metric_space_array_t<se3_0th_order_topology_t<double>, 1>&,
    const metric_space_array_t<se3_0th_order_rl_topology_t<double>, 1>&) const;
template topology_traits<
    metric_space_array_t<se3_1st_order_rl_topology_t<double>, 1>>::point_type
joint_limits_mapping<double>::map_to_space(
    const topology_traits<metric_space_array_t<se3_1st_order_topology_t<double>,
                                               1>>::point_type& pt,
    const metric_space_array_t<se3_1st_order_topology_t<double>, 1>&,
    const metric_space_array_t<se3_1st_order_rl_topology_t<double>, 1>&) const;
template topology_traits<
    metric_space_array_t<se3_2nd_order_rl_topology_t<double>, 1>>::point_type
joint_limits_mapping<double>::map_to_space(
    const topology_traits<metric_space_array_t<se3_2nd_order_topology_t<double>,
                                               1>>::point_type& pt,
    const metric_space_array_t<se3_2nd_order_topology_t<double>, 1>&,
    const metric_space_array_t<se3_2nd_order_rl_topology_t<double>, 1>&) const;

template topology_traits<
    metric_space_array_t<se3_0th_order_topology_t<double>, 1>>::point_type
joint_limits_mapping<double>::map_to_space(
    const topology_traits<metric_space_array_t<
        se3_0th_order_rl_topology_t<double>, 1>>::point_type& pt,
    const metric_space_array_t<se3_0th_order_rl_topology_t<double>, 1>&,
    const metric_space_array_t<se3_0th_order_topology_t<double>, 1>&) const;
template topology_traits<
    metric_space_array_t<se3_1st_order_topology_t<double>, 1>>::point_type
joint_limits_mapping<double>::map_to_space(
    const topology_traits<metric_space_array_t<
        se3_1st_order_rl_topology_t<double>, 1>>::point_type& pt,
    const metric_space_array_t<se3_1st_order_rl_topology_t<double>, 1>&,
    const metric_space_array_t<se3_1st_order_topology_t<double>, 1>&) const;
template topology_traits<
    metric_space_array_t<se3_2nd_order_topology_t<double>, 1>>::point_type
joint_limits_mapping<double>::map_to_space(
    const topology_traits<metric_space_array_t<
        se3_2nd_order_rl_topology_t<double>, 1>>::point_type& pt,
    const metric_space_array_t<se3_2nd_order_rl_topology_t<double>, 1>&,
    const metric_space_array_t<se3_2nd_order_topology_t<double>, 1>&) const;

}  // namespace ReaK::pp

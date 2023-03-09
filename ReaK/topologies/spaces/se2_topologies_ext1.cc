
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

#include "ReaK/core/base/defs.h"

#include "ReaK/topologies/spaces/joint_space_limits.inc"
#include "ReaK/topologies/spaces/se2_topologies.h"

namespace ReaK::pp {

// se2_0th_order_topology
template class metric_space_tuple<
    arithmetic_tuple<
        differentiable_space<
            time_topology, arithmetic_tuple<hyperbox_topology<vect<double, 2>>>,
            euclidean_tuple_distance>,
        differentiable_space<time_topology,
                             arithmetic_tuple<line_segment_topology<double>>,
                             euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se2_1st_order_topology
template class metric_space_tuple<
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
template class metric_space_tuple<
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
template class metric_space_tuple<
    arithmetic_tuple<
        reach_time_diff_space<
            time_topology, arithmetic_tuple<hyperbox_topology<vect<double, 2>>>,
            euclidean_tuple_distance>,
        reach_time_diff_space<time_topology,
                              arithmetic_tuple<line_segment_topology<double>>,
                              euclidean_tuple_distance>>,
    euclidean_tuple_distance>;

// se2_1st_order_rl_topology
template class metric_space_tuple<
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
template class metric_space_tuple<
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
template class temporal_space<
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
template class temporal_space<
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
template class temporal_space<
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
template class temporal_space<
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
template class temporal_space<
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
template class temporal_space<
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
template class temporal_space<
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
template class temporal_space<
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
template class temporal_space<
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

template metric_space_array_t<se2_0th_order_rl_topology_t<double>, 1>
joint_limits_mapping<double>::make_rl_joint_space(
    const metric_space_array_t<se2_0th_order_topology_t<double>, 1>&) const;
template metric_space_array_t<se2_1st_order_rl_topology_t<double>, 1>
joint_limits_mapping<double>::make_rl_joint_space(
    const metric_space_array_t<se2_1st_order_topology_t<double>, 1>&) const;
template metric_space_array_t<se2_2nd_order_rl_topology_t<double>, 1>
joint_limits_mapping<double>::make_rl_joint_space(
    const metric_space_array_t<se2_2nd_order_topology_t<double>, 1>&) const;

template metric_space_array_t<se2_0th_order_topology_t<double>, 1>
joint_limits_mapping<double>::make_normal_joint_space(
    const metric_space_array_t<se2_0th_order_rl_topology_t<double>, 1>&) const;
template metric_space_array_t<se2_1st_order_topology_t<double>, 1>
joint_limits_mapping<double>::make_normal_joint_space(
    const metric_space_array_t<se2_1st_order_rl_topology_t<double>, 1>&) const;
template metric_space_array_t<se2_2nd_order_topology_t<double>, 1>
joint_limits_mapping<double>::make_normal_joint_space(
    const metric_space_array_t<se2_2nd_order_rl_topology_t<double>, 1>&) const;

template topology_traits<
    metric_space_array_t<se2_0th_order_rl_topology_t<double>, 1>>::point_type
joint_limits_mapping<double>::map_to_space(
    const topology_traits<metric_space_array_t<se2_0th_order_topology_t<double>,
                                               1>>::point_type& pt,
    const metric_space_array_t<se2_0th_order_topology_t<double>, 1>&,
    const metric_space_array_t<se2_0th_order_rl_topology_t<double>, 1>&) const;
template topology_traits<
    metric_space_array_t<se2_1st_order_rl_topology_t<double>, 1>>::point_type
joint_limits_mapping<double>::map_to_space(
    const topology_traits<metric_space_array_t<se2_1st_order_topology_t<double>,
                                               1>>::point_type& pt,
    const metric_space_array_t<se2_1st_order_topology_t<double>, 1>&,
    const metric_space_array_t<se2_1st_order_rl_topology_t<double>, 1>&) const;
template topology_traits<
    metric_space_array_t<se2_2nd_order_rl_topology_t<double>, 1>>::point_type
joint_limits_mapping<double>::map_to_space(
    const topology_traits<metric_space_array_t<se2_2nd_order_topology_t<double>,
                                               1>>::point_type& pt,
    const metric_space_array_t<se2_2nd_order_topology_t<double>, 1>&,
    const metric_space_array_t<se2_2nd_order_rl_topology_t<double>, 1>&) const;

template topology_traits<
    metric_space_array_t<se2_0th_order_topology_t<double>, 1>>::point_type
joint_limits_mapping<double>::map_to_space(
    const topology_traits<metric_space_array_t<
        se2_0th_order_rl_topology_t<double>, 1>>::point_type& pt,
    const metric_space_array_t<se2_0th_order_rl_topology_t<double>, 1>&,
    const metric_space_array_t<se2_0th_order_topology_t<double>, 1>&) const;
template topology_traits<
    metric_space_array_t<se2_1st_order_topology_t<double>, 1>>::point_type
joint_limits_mapping<double>::map_to_space(
    const topology_traits<metric_space_array_t<
        se2_1st_order_rl_topology_t<double>, 1>>::point_type& pt,
    const metric_space_array_t<se2_1st_order_rl_topology_t<double>, 1>&,
    const metric_space_array_t<se2_1st_order_topology_t<double>, 1>&) const;
template topology_traits<
    metric_space_array_t<se2_2nd_order_topology_t<double>, 1>>::point_type
joint_limits_mapping<double>::map_to_space(
    const topology_traits<metric_space_array_t<
        se2_2nd_order_rl_topology_t<double>, 1>>::point_type& pt,
    const metric_space_array_t<se2_2nd_order_rl_topology_t<double>, 1>&,
    const metric_space_array_t<se2_2nd_order_topology_t<double>, 1>&) const;

}  // namespace ReaK::pp

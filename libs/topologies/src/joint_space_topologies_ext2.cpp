
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

#include "ReaK/core/base/defs.hpp"

#include "ReaK/topologies/spaces/joint_space_topologies.hpp"

namespace ReaK::pp {

// Multiple degrees of freedom (1-10) for 0th-order space:

template class metric_space_tuple<
    arithmetic_tuple<joint_space_0th_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<joint_space_0th_order_t<double>,
                     joint_space_0th_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<joint_space_0th_order_t<double>,
                     joint_space_0th_order_t<double>,
                     joint_space_0th_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>,
        joint_space_0th_order_t<double>, joint_space_0th_order_t<double>>,
    euclidean_tuple_distance>;

// Multiple degrees of freedom (1-10) for 1st-order space:

template class metric_space_tuple<
    arithmetic_tuple<joint_space_1st_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<joint_space_1st_order_t<double>,
                     joint_space_1st_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<joint_space_1st_order_t<double>,
                     joint_space_1st_order_t<double>,
                     joint_space_1st_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>,
        joint_space_1st_order_t<double>, joint_space_1st_order_t<double>>,
    euclidean_tuple_distance>;

// Multiple degrees of freedom (1-10) for 1st-order space:

template class metric_space_tuple<
    arithmetic_tuple<joint_space_2nd_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<joint_space_2nd_order_t<double>,
                     joint_space_2nd_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<joint_space_2nd_order_t<double>,
                     joint_space_2nd_order_t<double>,
                     joint_space_2nd_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>>,
    euclidean_tuple_distance>;

template class metric_space_tuple<
    arithmetic_tuple<
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>,
        joint_space_2nd_order_t<double>, joint_space_2nd_order_t<double>>,
    euclidean_tuple_distance>;

// joint_space_0th_order
template class temporal_space<Ndof_0th_order_space_t<double, 1>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_0th_order_space_t<double, 2>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_0th_order_space_t<double, 3>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_0th_order_space_t<double, 4>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_0th_order_space_t<double, 5>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_0th_order_space_t<double, 6>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_0th_order_space_t<double, 7>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_0th_order_space_t<double, 8>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_0th_order_space_t<double, 9>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_0th_order_space_t<double, 10>,
                              time_poisson_topology, spatial_distance_only>;

// joint_space_1st_order
template class temporal_space<Ndof_1st_order_space_t<double, 1>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_1st_order_space_t<double, 2>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_1st_order_space_t<double, 3>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_1st_order_space_t<double, 4>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_1st_order_space_t<double, 5>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_1st_order_space_t<double, 6>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_1st_order_space_t<double, 7>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_1st_order_space_t<double, 8>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_1st_order_space_t<double, 9>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_1st_order_space_t<double, 10>,
                              time_poisson_topology, spatial_distance_only>;

// joint_space_2nd_order
template class temporal_space<Ndof_2nd_order_space_t<double, 1>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_2nd_order_space_t<double, 2>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_2nd_order_space_t<double, 3>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_2nd_order_space_t<double, 4>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_2nd_order_space_t<double, 5>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_2nd_order_space_t<double, 6>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_2nd_order_space_t<double, 7>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_2nd_order_space_t<double, 8>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_2nd_order_space_t<double, 9>,
                              time_poisson_topology, spatial_distance_only>;
template class temporal_space<Ndof_2nd_order_space_t<double, 10>,
                              time_poisson_topology, spatial_distance_only>;

}  // namespace ReaK::pp

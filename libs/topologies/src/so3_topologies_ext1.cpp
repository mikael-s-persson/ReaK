
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

#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

#include <ReaK/topologies/spaces/so3_topologies.hpp>

namespace ReaK {

namespace pp {

template class quaternion_topology< double >;
template class rate_limited_quat_space< double >;
template class ang_velocity_3D_topology< double >;
template class ang_accel_3D_topology< double >;

template class temporal_space< quaternion_topology< double >, time_poisson_topology, spatial_distance_only >;
template class temporal_space< rate_limited_quat_space< double >, time_poisson_topology, spatial_distance_only >;
template class temporal_space< rate_limited_quat_space< double >, time_poisson_topology, reach_plus_time_metric >;
template class temporal_space< ang_velocity_3D_topology< double >, time_poisson_topology, spatial_distance_only >;
template class temporal_space< ang_accel_3D_topology< double >, time_poisson_topology, spatial_distance_only >;

// so3_0th_order_topology
template class differentiable_space< time_topology, arithmetic_tuple< quaternion_topology< double > >,
                                     euclidean_tuple_distance >;
// so3_1st_order_topology
template class differentiable_space< time_topology, arithmetic_tuple< quaternion_topology< double >,
                                                                      ang_velocity_3D_topology< double > >,
                                     euclidean_tuple_distance >;
// so3_2nd_order_topology
template class differentiable_space< time_topology, arithmetic_tuple< quaternion_topology< double >,
                                                                      ang_velocity_3D_topology< double >,
                                                                      ang_accel_3D_topology< double > >,
                                     euclidean_tuple_distance >;

// so3_0th_order_topology
template class temporal_space< differentiable_space< time_topology, arithmetic_tuple< quaternion_topology< double > >,
                                                     euclidean_tuple_distance >,
                               time_poisson_topology, spatial_distance_only >;
// so3_1st_order_topology
template class temporal_space< differentiable_space< time_topology,
                                                     arithmetic_tuple< quaternion_topology< double >,
                                                                       ang_velocity_3D_topology< double > >,
                                                     euclidean_tuple_distance >,
                               time_poisson_topology, spatial_distance_only >;
// so3_2nd_order_topology
template class temporal_space< differentiable_space< time_topology,
                                                     arithmetic_tuple< quaternion_topology< double >,
                                                                       ang_velocity_3D_topology< double >,
                                                                       ang_accel_3D_topology< double > >,
                                                     euclidean_tuple_distance >,
                               time_poisson_topology, spatial_distance_only >;


// so3_0th_order_rl_topology
template class reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space< double > >,
                                      euclidean_tuple_distance >;
// so3_1st_order_rl_topology
template class reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space< double >,
                                                                       ang_velocity_3D_topology< double > >,
                                      euclidean_tuple_distance >;
// so3_2nd_order_rl_topology
template class reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space< double >,
                                                                       ang_velocity_3D_topology< double >,
                                                                       ang_accel_3D_topology< double > >,
                                      euclidean_tuple_distance >;

// so3_0th_order_rl_topology
template class temporal_space< reach_time_diff_space< time_topology,
                                                      arithmetic_tuple< rate_limited_quat_space< double > >,
                                                      euclidean_tuple_distance >,
                               time_poisson_topology, spatial_distance_only >;
// so3_1st_order_rl_topology
template class temporal_space< reach_time_diff_space< time_topology,
                                                      arithmetic_tuple< rate_limited_quat_space< double >,
                                                                        ang_velocity_3D_topology< double > >,
                                                      euclidean_tuple_distance >,
                               time_poisson_topology, spatial_distance_only >;
// so3_2nd_order_rl_topology
template class temporal_space< reach_time_diff_space< time_topology,
                                                      arithmetic_tuple< rate_limited_quat_space< double >,
                                                                        ang_velocity_3D_topology< double >,
                                                                        ang_accel_3D_topology< double > >,
                                                      euclidean_tuple_distance >,
                               time_poisson_topology, spatial_distance_only >;

// so3_0th_order_rl_topology
template class temporal_space< reach_time_diff_space< time_topology,
                                                      arithmetic_tuple< rate_limited_quat_space< double > >,
                                                      euclidean_tuple_distance >,
                               time_poisson_topology, reach_plus_time_metric >;
// so3_1st_order_rl_topology
template class temporal_space< reach_time_diff_space< time_topology,
                                                      arithmetic_tuple< rate_limited_quat_space< double >,
                                                                        ang_velocity_3D_topology< double > >,
                                                      euclidean_tuple_distance >,
                               time_poisson_topology, reach_plus_time_metric >;
// so3_2nd_order_rl_topology
template class temporal_space< reach_time_diff_space< time_topology,
                                                      arithmetic_tuple< rate_limited_quat_space< double >,
                                                                        ang_velocity_3D_topology< double >,
                                                                        ang_accel_3D_topology< double > >,
                                                      euclidean_tuple_distance >,
                               time_poisson_topology, reach_plus_time_metric >;
};
};

#else

namespace ReaK {

namespace pp {

void dummy_so3_topologies_externs_1_symbol(){};
};
};

#endif

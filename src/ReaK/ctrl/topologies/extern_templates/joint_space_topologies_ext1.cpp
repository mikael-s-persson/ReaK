
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

#include "base/defs.hpp"

#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

#include "topologies/joint_space_topologies.hpp"

namespace ReaK {

namespace pp {

// joint_space_0th_order
template class differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance >;
// joint_space_1st_order
template class differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >;
// joint_space_2nd_order
template class differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >;

// joint_space_0th_order
template class temporal_space< differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;
// joint_space_1st_order
template class temporal_space< differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;
// joint_space_2nd_order
template class temporal_space< differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;


// rl_joint_space_0th_order
template class reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance >;
// rl_joint_space_1st_order
template class reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >;
// rl_joint_space_2nd_order
template class reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >;

// rl_joint_space_0th_order
template class temporal_space< reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;
// rl_joint_space_1st_order
template class temporal_space< reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;
// rl_joint_space_2nd_order
template class temporal_space< reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;

// rl_joint_space_0th_order
template class temporal_space< reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, reach_plus_time_metric>;
// rl_joint_space_1st_order
template class temporal_space< reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, reach_plus_time_metric>;
// rl_joint_space_2nd_order
template class temporal_space< reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, reach_plus_time_metric>;


};

};

#else

namespace ReaK {

namespace pp {

void dummy_joint_space_topologies_externs_1_symbol() { };

};

};

#endif















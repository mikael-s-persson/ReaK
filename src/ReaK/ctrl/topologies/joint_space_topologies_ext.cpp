
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

#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

#include "joint_space_topologies.hpp"

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




// Multiple degrees of freedom (1-10) for 0th-order space:

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;


  
// Multiple degrees of freedom (1-10) for 1st-order space:

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

  
// Multiple degrees of freedom (1-10) for 1st-order space:

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

  


// Multiple degrees of freedom (1-10) for 0th-order space:

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;


  
// Multiple degrees of freedom (1-10) for 1st-order space:

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

  
// Multiple degrees of freedom (1-10) for 1st-order space:

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;





// Multiple degrees of freedom (1-10) for 0th-order space:

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;


  
// Multiple degrees of freedom (1-10) for 1st-order space:

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

  
// Multiple degrees of freedom (1-10) for 1st-order space:

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;



  
  
// joint_space_0th_order
template class temporal_space< Ndof_0th_order_space<double, 1>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_space<double, 2>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_space<double, 3>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_space<double, 4>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_space<double, 5>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_space<double, 6>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_space<double, 7>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_space<double, 8>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_space<double, 9>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_space<double,10>::type, time_poisson_topology, spatial_distance_only>;

// joint_space_1st_order
template class temporal_space< Ndof_1st_order_space<double, 1>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_space<double, 2>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_space<double, 3>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_space<double, 4>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_space<double, 5>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_space<double, 6>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_space<double, 7>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_space<double, 8>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_space<double, 9>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_space<double,10>::type, time_poisson_topology, spatial_distance_only>;

// joint_space_2nd_order
template class temporal_space< Ndof_2nd_order_space<double, 1>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_space<double, 2>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_space<double, 3>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_space<double, 4>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_space<double, 5>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_space<double, 6>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_space<double, 7>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_space<double, 8>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_space<double, 9>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_space<double,10>::type, time_poisson_topology, spatial_distance_only>;

  
  
// rl_joint_space_0th_order
template class temporal_space< Ndof_0th_order_rl_space<double, 1>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 2>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 3>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 4>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 5>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 6>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 7>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 8>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 9>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double,10>::type, time_poisson_topology, spatial_distance_only>;

// rl_joint_space_1st_order
template class temporal_space< Ndof_1st_order_rl_space<double, 1>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 2>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 3>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 4>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 5>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 6>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 7>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 8>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 9>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double,10>::type, time_poisson_topology, spatial_distance_only>;

// rl_joint_space_2nd_order
template class temporal_space< Ndof_2nd_order_rl_space<double, 1>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 2>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 3>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 4>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 5>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 6>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 7>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 8>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 9>::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double,10>::type, time_poisson_topology, spatial_distance_only>;


  
// rl_joint_space_0th_order
template class temporal_space< Ndof_0th_order_rl_space<double, 1>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 2>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 3>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 4>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 5>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 6>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 7>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 8>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 9>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double,10>::type, time_poisson_topology, reach_plus_time_metric>;

// rl_joint_space_1st_order
template class temporal_space< Ndof_1st_order_rl_space<double, 1>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 2>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 3>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 4>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 5>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 6>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 7>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 8>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 9>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double,10>::type, time_poisson_topology, reach_plus_time_metric>;

// rl_joint_space_2nd_order
template class temporal_space< Ndof_2nd_order_rl_space<double, 1>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 2>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 3>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 4>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 5>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 6>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 7>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 8>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 9>::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double,10>::type, time_poisson_topology, reach_plus_time_metric>;



// rl_joint_space_0th_order
template class temporal_space< Ndof_0th_order_rl_space<double, 1, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 2, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 3, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 4, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 5, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 6, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 7, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 8, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double, 9, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_0th_order_rl_space<double,10, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;

// rl_joint_space_1st_order
template class temporal_space< Ndof_1st_order_rl_space<double, 1, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 2, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 3, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 4, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 5, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 6, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 7, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 8, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double, 9, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_1st_order_rl_space<double,10, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;

// rl_joint_space_2nd_order
template class temporal_space< Ndof_2nd_order_rl_space<double, 1, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 2, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 3, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 4, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 5, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 6, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 7, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 8, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 9, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;
template class temporal_space< Ndof_2nd_order_rl_space<double,10, inf_norm_tuple_distance >::type, time_poisson_topology, spatial_distance_only>;


  
// rl_joint_space_0th_order
template class temporal_space< Ndof_0th_order_rl_space<double, 1, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 2, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 3, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 4, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 5, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 6, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 7, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 8, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double, 9, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_0th_order_rl_space<double,10, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;

// rl_joint_space_1st_order
template class temporal_space< Ndof_1st_order_rl_space<double, 1, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 2, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 3, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 4, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 5, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 6, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 7, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 8, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double, 9, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_1st_order_rl_space<double,10, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;

// rl_joint_space_2nd_order
template class temporal_space< Ndof_2nd_order_rl_space<double, 1, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 2, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 3, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 4, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 5, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 6, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 7, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 8, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double, 9, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;
template class temporal_space< Ndof_2nd_order_rl_space<double,10, inf_norm_tuple_distance >::type, time_poisson_topology, reach_plus_time_metric>;





};

};

#else

namespace ReaK {

namespace pp {

void dummy_joint_space_topologies_externs_symbol() { };

};

};

#endif















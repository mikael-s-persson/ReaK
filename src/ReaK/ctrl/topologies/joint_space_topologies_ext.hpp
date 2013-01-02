/**
 * \file joint_space_topologies_ext.hpp
 * 
 * This library provides extern template declarations for classes defined in the joint_space_topologies.hpp header file.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2012
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

#ifndef REAK_JOINT_SPACE_TOPOLOGIES_EXT_HPP
#define REAK_JOINT_SPACE_TOPOLOGIES_EXT_HPP

#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

#include "base/defs.hpp"
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "differentiable_space.hpp"
#include "metric_space_tuple.hpp"
#include "rate_limited_spaces.hpp"

#include "line_topology.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "kinetostatics/gen_coord.hpp"

#include "time_poisson_topology.hpp"
#include "temporal_space.hpp"
#include "reachability_space.hpp"

namespace ReaK {

namespace pp {

// joint_space_0th_order
extern template class differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance >;
// joint_space_1st_order
extern template class differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >;
// joint_space_2nd_order
extern template class differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >;

// joint_space_0th_order
extern template class temporal_space< differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;
// joint_space_1st_order
extern template class temporal_space< differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;
// joint_space_2nd_order
extern template class temporal_space< differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;


// rl_joint_space_0th_order
extern template class reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance >;
// rl_joint_space_1st_order
extern template class reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >;
// rl_joint_space_2nd_order
extern template class reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >;

// rl_joint_space_0th_order
extern template class temporal_space< reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;
// rl_joint_space_1st_order
extern template class temporal_space< reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;
// rl_joint_space_2nd_order
extern template class temporal_space< reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;

// rl_joint_space_0th_order
extern template class temporal_space< reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, reach_plus_time_metric>;
// rl_joint_space_1st_order
extern template class temporal_space< reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, reach_plus_time_metric>;
// rl_joint_space_2nd_order
extern template class temporal_space< reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance >, time_poisson_topology, reach_plus_time_metric>;



// Multiple degrees of freedom (1-10) for 0th-order space:

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type,
  joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
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

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type,
  joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
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

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type,
  joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
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

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
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

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
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

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, euclidean_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
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

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type,
  rl_joint_space_0th_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
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

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type,
  rl_joint_space_1st_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
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

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type,
  rl_joint_space_2nd_order<double>::type>, inf_norm_tuple_distance >;

extern template class metric_space_tuple< arithmetic_tuple< 
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
extern template class temporal_space< Ndof_0th_order_space<double, 1>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_space<double, 2>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_space<double, 3>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_space<double, 4>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_space<double, 5>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_space<double, 6>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_space<double, 7>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_space<double, 8>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_space<double, 9>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_space<double,10>::type, time_poisson_topology, spatial_distance_only>;

// joint_space_1st_order
extern template class temporal_space< Ndof_1st_order_space<double, 1>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_space<double, 2>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_space<double, 3>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_space<double, 4>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_space<double, 5>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_space<double, 6>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_space<double, 7>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_space<double, 8>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_space<double, 9>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_space<double,10>::type, time_poisson_topology, spatial_distance_only>;

// joint_space_2nd_order
extern template class temporal_space< Ndof_2nd_order_space<double, 1>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_space<double, 2>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_space<double, 3>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_space<double, 4>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_space<double, 5>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_space<double, 6>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_space<double, 7>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_space<double, 8>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_space<double, 9>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_space<double,10>::type, time_poisson_topology, spatial_distance_only>;

  
  
// rl_joint_space_0th_order
extern template class temporal_space< Ndof_0th_order_rl_space<double, 1>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 2>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 3>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 4>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 5>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 6>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 7>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 8>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 9>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_0th_order_rl_space<double,10>::type, time_poisson_topology, spatial_distance_only>;

// rl_joint_space_1st_order
extern template class temporal_space< Ndof_1st_order_rl_space<double, 1>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 2>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 3>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 4>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 5>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 6>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 7>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 8>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 9>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_1st_order_rl_space<double,10>::type, time_poisson_topology, spatial_distance_only>;

// rl_joint_space_2nd_order
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 1>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 2>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 3>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 4>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 5>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 6>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 7>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 8>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 9>::type, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double,10>::type, time_poisson_topology, spatial_distance_only>;

  
// rl_joint_space_0th_order
extern template class temporal_space< Ndof_0th_order_rl_space<double, 1>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 2>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 3>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 4>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 5>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 6>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 7>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 8>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_0th_order_rl_space<double, 9>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_0th_order_rl_space<double,10>::type, time_poisson_topology, reach_plus_time_metric>;

// rl_joint_space_1st_order
extern template class temporal_space< Ndof_1st_order_rl_space<double, 1>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 2>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 3>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 4>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 5>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 6>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 7>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 8>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_1st_order_rl_space<double, 9>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_1st_order_rl_space<double,10>::type, time_poisson_topology, reach_plus_time_metric>;

// rl_joint_space_2nd_order
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 1>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 2>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 3>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 4>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 5>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 6>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 7>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 8>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double, 9>::type, time_poisson_topology, reach_plus_time_metric>;
extern template class temporal_space< Ndof_2nd_order_rl_space<double,10>::type, time_poisson_topology, reach_plus_time_metric>;

  
  

};

};

#endif


#endif









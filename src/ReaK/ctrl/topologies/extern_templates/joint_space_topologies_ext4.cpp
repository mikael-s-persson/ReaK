
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

#include <ReaK/ctrl/topologies/joint_space_topologies.hpp>

namespace ReaK {

namespace pp {


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

void dummy_joint_space_topologies_externs_4_symbol() { };

};

};

#endif















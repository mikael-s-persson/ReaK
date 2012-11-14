
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

#include "se2_topologies.hpp"

namespace ReaK {

namespace pp {
  
  
// se2_0th_order_topology
template class metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;

// se2_1st_order_topology
template class metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> >, hyperball_topology< vect<double,2> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;

// se2_2nd_order_topology
template class metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> >, hyperball_topology< vect<double,2> >, hyperball_topology< vect<double,2> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;



// se2_0th_order_rl_topology
template class metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;

// se2_1st_order_rl_topology
template class metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> >, hyperball_topology< vect<double,2> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;

// se2_2nd_order_rl_topology
template class metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> >, hyperball_topology< vect<double,2> >, hyperball_topology< vect<double,2> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;


};


template frame_2D<double> get_frame_2D(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>,    vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double, double > >& pt);

template frame_2D<double> get_frame_2D(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double > >& pt);

template frame_2D<double> get_frame_2D(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                          arithmetic_tuple< double > >& pt);


template void set_frame_2D(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>,    vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double, double > >& pt,
  const frame_2D<double>& p);

template void set_frame_2D(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double > >& pt,
  const frame_2D<double>& p);

template void set_frame_2D(
  arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                    arithmetic_tuple< double > >& pt,
  const frame_2D<double>& p);

template pose_2D<double> get_pose_2D(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                          arithmetic_tuple< double > >& pt);

template void set_pose_2D(
  arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                    arithmetic_tuple< double > >& pt,
  const pose_2D<double>& p);

template const double& get_rotation(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double, double > >& pt);

template const double& get_rotation(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double > >& pt);

template const double& get_rotation(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                          arithmetic_tuple< double > >& pt);

template void set_rotation(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double, double > >& pt,
  const double& q);

template void set_rotation(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double > >& pt,
  const double& q);

template void set_rotation(
  arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                    arithmetic_tuple< double > >& pt,
  const double& q);

template const vect<double,2>& get_position(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double, double > >& pt);

template const vect<double,2>& get_position(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double > >& pt);

template const vect<double,2>& get_position(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                          arithmetic_tuple< double > >& pt);

template void set_position(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double, double > >& pt,
  const vect<double,2>& p);

template void set_position(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double > >& pt,
  const vect<double,2>& p);

template void set_position(
  arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                    arithmetic_tuple< double > >& pt,
  const vect<double,2>& p);

template const double& get_ang_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double, double > >& pt);

template const double& get_ang_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double > >& pt);

template void set_ang_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double, double > >& pt,
  const double& p);

template void set_ang_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double > >& pt,
  const double& p);

template const vect<double,2>& get_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double, double > >& pt);

template const vect<double,2>& get_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double > >& pt);

template void set_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double, double > >& pt,
  const vect<double,2>& p);

template void set_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double > >& pt,
  const vect<double,2>& p);

template const double& get_ang_acceleration(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double, double > >& pt);

template void set_ang_acceleration(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double, double > >& pt,
  const double& p);

template const vect<double,2>& get_acceleration(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double, double > >& pt);

template void set_acceleration(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double, double > >& pt,
  const vect<double,2>& p);


};

#else

namespace ReaK {

namespace pp {

void dummy_se2_topologies_externs_symbol() { };

};

};

#endif















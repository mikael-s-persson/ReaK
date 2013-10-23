/**
 * \file se2_topologies.hpp
 * 
 * This library provides classes that define topologies on SE(2) (2D rigid-body motion). 
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

#ifndef REAK_SE2_TOPOLOGIES_HPP
#define REAK_SE2_TOPOLOGIES_HPP


#include "base/defs.hpp"
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "differentiable_space.hpp"
#include "rate_limited_spaces.hpp"
#include "metric_space_tuple.hpp"
#include "hyperbox_topology.hpp"
#include "hyperball_topology.hpp"
#include "line_topology.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "lin_alg/vect_alg.hpp"
#include "kinetostatics/frame_2D.hpp"

namespace ReaK {

namespace pp {


  

/**
 * This meta-function defines the type for a 0th order SE(2) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se2_0th_order_topology {
  typedef 
    metric_space_tuple< arithmetic_tuple<
      differentiable_space< 
        time_topology, 
        arithmetic_tuple< hyperbox_topology< vect<T,2> > >, 
        DistanceMetric 
      >,
      differentiable_space< 
        time_topology, 
        arithmetic_tuple< line_segment_topology<T> >, 
        DistanceMetric 
      > >,
      DistanceMetric 
    > type;
};

/**
 * This meta-function defines the type for a 1st order SE(2) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se2_1st_order_topology {
  typedef 
    metric_space_tuple< arithmetic_tuple<
      differentiable_space< 
        time_topology, 
        arithmetic_tuple< 
          hyperbox_topology< vect<T,2> >,
          hyperball_topology< vect<T,2> >
        >, 
        DistanceMetric 
      >,
      differentiable_space< 
        time_topology, 
        arithmetic_tuple< 
          line_segment_topology<T>,
          line_segment_topology<T>
        >, 
        DistanceMetric 
      > >,
      DistanceMetric 
    > type;
};

/**
 * This meta-function defines the type for a 2nd order SE(2) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se2_2nd_order_topology {
  typedef 
    metric_space_tuple< arithmetic_tuple<
      differentiable_space< 
        time_topology, 
        arithmetic_tuple< 
          hyperbox_topology< vect<T,2> >,
          hyperball_topology< vect<T,2> >,
          hyperball_topology< vect<T,2> >
        >, 
        DistanceMetric 
      >,
      differentiable_space< 
        time_topology, 
        arithmetic_tuple< 
          line_segment_topology<T>,
          line_segment_topology<T>,
          line_segment_topology<T>
        >, 
        DistanceMetric 
      > >,
      DistanceMetric 
    > type;
};



template <typename T, int Order, typename DistanceMetric = euclidean_tuple_distance>
struct se2_topology {
  typedef typename boost::mpl::if_<
    boost::mpl::equal_to< boost::mpl::int_<0>, boost::mpl::int_<Order> >,
    typename se2_0th_order_topology<T,DistanceMetric>::type,
    typename boost::mpl::if_< 
      boost::mpl::equal_to< boost::mpl::int_<1>, boost::mpl::int_<Order> >,
      typename se2_1st_order_topology<T,DistanceMetric>::type,
      typename se2_2nd_order_topology<T,DistanceMetric>::type
    >::type
  >::type type;
};



template <typename SE2Space>
struct is_se2_space : boost::mpl::false_ { };


template <typename T, typename DistanceMetric>
struct is_se2_space< 
    metric_space_tuple< arithmetic_tuple<
      differentiable_space< 
        time_topology, 
        arithmetic_tuple< hyperbox_topology< vect<T,2> > >, 
        DistanceMetric 
      >,
      differentiable_space< 
        time_topology, 
        arithmetic_tuple< line_segment_topology<T> >, 
        DistanceMetric 
      > >,
      DistanceMetric 
    > > : boost::mpl::true_ { };

template <typename T, typename DistanceMetric>
struct is_se2_space< 
    metric_space_tuple< arithmetic_tuple<
      differentiable_space< 
        time_topology, 
        arithmetic_tuple< 
          hyperbox_topology< vect<T,2> >,
          hyperball_topology< vect<T,2> >
        >, 
        DistanceMetric 
      >,
      differentiable_space< 
        time_topology, 
        arithmetic_tuple< 
          line_segment_topology<T>,
          line_segment_topology<T>
        >, 
        DistanceMetric 
      > >,
      DistanceMetric 
    > > : boost::mpl::true_ { };

template <typename T, typename DistanceMetric>
struct is_se2_space< 
    metric_space_tuple< arithmetic_tuple<
      differentiable_space< 
        time_topology, 
        arithmetic_tuple< 
          hyperbox_topology< vect<T,2> >,
          hyperball_topology< vect<T,2> >,
          hyperball_topology< vect<T,2> >
        >, 
        DistanceMetric 
      >,
      differentiable_space< 
        time_topology, 
        arithmetic_tuple< 
          line_segment_topology<T>,
          line_segment_topology<T>,
          line_segment_topology<T>
        >, 
        DistanceMetric 
      > >,
      DistanceMetric 
    > > : boost::mpl::true_ { };



/**
 * This meta-function defines the type for a 0th order SE(2) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se2_0th_order_rl_topology {
  typedef 
    metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< 
        time_topology, 
        arithmetic_tuple< hyperbox_topology< vect<T,2> > >, 
        DistanceMetric 
      >,
      reach_time_diff_space< 
        time_topology, 
        arithmetic_tuple< line_segment_topology<T> >, 
        DistanceMetric 
      > >,
      DistanceMetric 
    > type;
};

/**
 * This meta-function defines the type for a 1st order SE(2) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se2_1st_order_rl_topology {
  typedef 
    metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< 
        time_topology, 
        arithmetic_tuple< 
          hyperbox_topology< vect<T,2> >,
          hyperball_topology< vect<T,2> >
        >, 
        DistanceMetric 
      >,
      reach_time_diff_space< 
        time_topology, 
        arithmetic_tuple< 
          line_segment_topology<T>,
          line_segment_topology<T>
        >, 
        DistanceMetric 
      > >,
      DistanceMetric 
    > type;
};

/**
 * This meta-function defines the type for a 2nd order SE(2) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se2_2nd_order_rl_topology {
  typedef 
    metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< 
        time_topology, 
        arithmetic_tuple< 
          hyperbox_topology< vect<T,2> >,
          hyperball_topology< vect<T,2> >,
          hyperball_topology< vect<T,2> >
        >, 
        DistanceMetric 
      >,
      reach_time_diff_space< 
        time_topology, 
        arithmetic_tuple< 
          line_segment_topology<T>,
          line_segment_topology<T>,
          line_segment_topology<T>
        >, 
        DistanceMetric 
      > >,
      DistanceMetric 
    > type;
};



template <typename T, int Order, typename DistanceMetric = euclidean_tuple_distance>
struct se2_rl_topology {
  typedef typename boost::mpl::if_<
    boost::mpl::equal_to< boost::mpl::int_<0>, boost::mpl::int_<Order> >,
    typename se2_0th_order_rl_topology<T,DistanceMetric>::type,
    typename boost::mpl::if_< 
      boost::mpl::equal_to< boost::mpl::int_<1>, boost::mpl::int_<Order> >,
      typename se2_1st_order_rl_topology<T,DistanceMetric>::type,
      typename se2_2nd_order_rl_topology<T,DistanceMetric>::type
    >::type
  >::type type;
};




template <typename SE2Space>
struct is_rate_limited_se2_space : boost::mpl::false_ { };


template <typename T, typename DistanceMetric>
struct is_rate_limited_se2_space< 
    metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< 
        time_topology, 
        arithmetic_tuple< hyperbox_topology< vect<T,2> > >, 
        DistanceMetric 
      >,
      reach_time_diff_space< 
        time_topology, 
        arithmetic_tuple< line_segment_topology<T> >, 
        DistanceMetric 
      > >,
      DistanceMetric 
    > > : boost::mpl::true_ { };

template <typename T, typename DistanceMetric>
struct is_rate_limited_se2_space< 
    metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< 
        time_topology, 
        arithmetic_tuple< 
          hyperbox_topology< vect<T,2> >,
          hyperball_topology< vect<T,2> >
        >, 
        DistanceMetric 
      >,
      reach_time_diff_space< 
        time_topology, 
        arithmetic_tuple< 
          line_segment_topology<T>,
          line_segment_topology<T>
        >, 
        DistanceMetric 
      > >,
      DistanceMetric 
    > > : boost::mpl::true_ { };

template <typename T, typename DistanceMetric>
struct is_rate_limited_se2_space< 
    metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< 
        time_topology, 
        arithmetic_tuple<  
          hyperbox_topology< vect<T,2> >,
          hyperball_topology< vect<T,2> >,
          hyperball_topology< vect<T,2> >
        >, 
        DistanceMetric 
      >,
      reach_time_diff_space< 
        time_topology, 
        arithmetic_tuple< 
          line_segment_topology<T>,
          line_segment_topology<T>,
          line_segment_topology<T>
        >, 
        DistanceMetric 
      > >,
      DistanceMetric 
    > > : boost::mpl::true_ { };






};





// Because of ADL rules, the get functions for the arithmetic-tuple types that represent SE(3) states should be in the ReaK namespace.


template <typename T>
frame_2D<T> get_frame_2D(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2>,    vect<T,2>, vect<T,2> >,
                          arithmetic_tuple< T, T, T > >& pt) {
  return frame_2D<T>(weak_ptr< pose_2D<T> >(),
                     get<0>(get<0>(pt)), 
                     rot_mat_2D<T>(get<0>(get<1>(pt))), 
                     get<1>(get<0>(pt)), 
                     get<1>(get<1>(pt)), 
                     get<2>(get<0>(pt)), 
                     get<2>(get<1>(pt)),
                     vect<T,2>(0.0,0.0), 
                     0.0);
};

template <typename T>
frame_2D<T> get_frame_2D(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2> >,
                          arithmetic_tuple< T, T > >& pt) {
  return frame_2D<T>(weak_ptr< pose_2D<T> >(),
                     get<0>(get<0>(pt)), 
                     rot_mat_2D<T>(get<0>(get<1>(pt))), 
                     get<1>(get<0>(pt)), 
                     get<1>(get<1>(pt)), 
                     vect<T,2>(0.0,0.0), 
                     0.0,
                     vect<T,2>(0.0,0.0), 
                     0.0);
};

template <typename T>
frame_2D<T> get_frame_2D(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2> >,
                          arithmetic_tuple< T > >& pt) {
  return frame_2D<T>(weak_ptr< pose_2D<T> >(),
                     get<0>(get<0>(pt)), 
                     rot_mat_2D<T>(get<0>(get<1>(pt))), 
                     vect<T,2>(0.0,0.0), 
                     0.0,
                     vect<T,2>(0.0,0.0), 
                     0.0, 
                     vect<T,2>(0.0,0.0), 
                     0.0);
};


template <typename T>
void set_frame_2D(
  arithmetic_tuple< arithmetic_tuple< vect<T,2>,    vect<T,2>, vect<T,2> >,
                    arithmetic_tuple< T, T, T > >& pt,
  const frame_2D<T>& p) {
  get<0>(get<0>(pt)) = p.Position;
  get<0>(get<1>(pt)) = p.Rotation.getAngle();
  get<1>(get<0>(pt)) = p.Velocity;
  get<1>(get<1>(pt)) = p.AngVelocity;
  get<2>(get<0>(pt)) = p.Acceleration;
  get<2>(get<1>(pt)) = p.AngAcceleration;
};

template <typename T>
void set_frame_2D(
  arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2> >,
                    arithmetic_tuple< T, T > >& pt,
  const frame_2D<T>& p) {
  get<0>(get<0>(pt)) = p.Position;
  get<0>(get<1>(pt)) = p.Rotation.getAngle();
  get<1>(get<0>(pt)) = p.Velocity;
  get<1>(get<1>(pt)) = p.AngVelocity;
};

template <typename T>
void set_frame_2D(
  arithmetic_tuple< arithmetic_tuple< vect<T,2> >,
                    arithmetic_tuple< T > >& pt,
  const frame_2D<T>& p) {
  get<0>(get<0>(pt)) = p.Position;
  get<0>(get<1>(pt)) = p.Rotation.getAngle();
};



template <typename T>
pose_2D<T> get_pose_2D(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2> >,
                          arithmetic_tuple< T > >& pt) {
  return pose_2D<T>(weak_ptr< pose_2D<T> >(),
                    get<0>(get<0>(pt)), 
                    rot_mat_2D<T>(get<0>(get<1>(pt))));
};

template <typename T>
void set_pose_2D(
  arithmetic_tuple< arithmetic_tuple< vect<T,2> >,
                    arithmetic_tuple< T > >& pt,
  const pose_2D<T>& p) {
  get<0>(get<0>(pt)) = p.Position;
  get<0>(get<1>(pt)) = p.Rotation.getAngle();
};





template <typename T>
const T& get_rotation(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2>, vect<T,2> >,
                          arithmetic_tuple< T, T, T > >& pt) {
  return get<0>(get<1>(pt));
};

template <typename T>
const T& get_rotation(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2> >,
                          arithmetic_tuple< T, T > >& pt) {
  return get<0>(get<1>(pt));
};

template <typename T>
const T& get_rotation(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2> >,
                          arithmetic_tuple< T > >& pt) {
  return get<0>(get<1>(pt));
};

template <typename T>
void set_rotation(
  arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2>, vect<T,2> >,
                    arithmetic_tuple< T, T, T > >& pt,
  const T& q) {
  get<0>(get<1>(pt)) = q;
};

template <typename T>
void set_rotation(
  arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2> >,
                    arithmetic_tuple< T, T > >& pt,
  const T& q) {
  get<0>(get<1>(pt)) = q;
};

template <typename T>
void set_rotation(
  arithmetic_tuple< arithmetic_tuple< vect<T,2> >,
                    arithmetic_tuple< T > >& pt,
  const T& q) {
  get<0>(get<1>(pt)) = q;
};



template <typename T>
const vect<T,2>& get_position(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2>, vect<T,2> >,
                          arithmetic_tuple< T, T, T > >& pt) {
  return get<0>(get<0>(pt));
};

template <typename T>
const vect<T,2>& get_position(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2> >,
                          arithmetic_tuple< T, T > >& pt) {
  return get<0>(get<0>(pt));
};

template <typename T>
const vect<T,2>& get_position(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2> >,
                          arithmetic_tuple< T > >& pt) {
  return get<0>(get<0>(pt));
};

template <typename T>
void set_position(
  arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2>, vect<T,2> >,
                    arithmetic_tuple< T, T, T > >& pt,
  const vect<T,2>& p) {
  get<0>(get<0>(pt)) = p;
};

template <typename T>
void set_position(
  arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2> >,
                    arithmetic_tuple< T, T > >& pt,
  const vect<T,2>& p) {
  get<0>(get<0>(pt)) = p;
};

template <typename T>
void set_position(
  arithmetic_tuple< arithmetic_tuple< vect<T,2> >,
                    arithmetic_tuple< T > >& pt,
  const vect<T,2>& p) {
  get<0>(get<0>(pt)) = p;
};



template <typename T>
const T& get_ang_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2>, vect<T,2> >,
                          arithmetic_tuple< T, T, T > >& pt) {
  return get<1>(get<1>(pt));
};

template <typename T>
const T& get_ang_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2> >,
                          arithmetic_tuple< T, T > >& pt) {
  return get<1>(get<1>(pt));
};

template <typename T>
void set_ang_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2>, vect<T,2> >,
                    arithmetic_tuple< T, T, T > >& pt,
  const T& p) {
  get<1>(get<1>(pt)) = p;
};

template <typename T>
void set_ang_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2> >,
                    arithmetic_tuple< T, T > >& pt,
  const T& p) {
  get<1>(get<1>(pt)) = p;
};



template <typename T>
const vect<T,2>& get_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2>, vect<T,2> >,
                          arithmetic_tuple< T, T, T > >& pt) {
  return get<1>(get<0>(pt));
};

template <typename T>
const vect<T,2>& get_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2> >,
                          arithmetic_tuple< T, T > >& pt) {
  return get<1>(get<0>(pt));
};

template <typename T>
void set_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2>, vect<T,2> >,
                    arithmetic_tuple< T, T, T > >& pt,
  const vect<T,2>& p) {
  get<1>(get<0>(pt)) = p;
};

template <typename T>
void set_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2> >,
                    arithmetic_tuple< T, T > >& pt,
  const vect<T,2>& p) {
  get<1>(get<0>(pt)) = p;
};




template <typename T>
const T& get_ang_acceleration(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2>, vect<T,2> >,
                          arithmetic_tuple< T, T, T > >& pt) {
  return get<2>(get<1>(pt));
};

template <typename T>
void set_ang_acceleration(
  arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2>, vect<T,2> >,
                    arithmetic_tuple< T, T, T > >& pt,
  const T& p) {
  get<2>(get<1>(pt)) = p;
};



template <typename T>
const vect<T,2>& get_acceleration(
  const arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2>, vect<T,2> >,
                          arithmetic_tuple< T, T, T > >& pt) {
  return get<2>(get<0>(pt));
};

template <typename T>
void set_acceleration(
  arithmetic_tuple< arithmetic_tuple< vect<T,2>, vect<T,2>, vect<T,2> >,
                    arithmetic_tuple< T, T, T > >& pt,
  const vect<T,2>& p) {
  get<2>(get<0>(pt)) = p;
};









};




#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

#include "time_poisson_topology.hpp"
#include "temporal_space.hpp"
#include "reachability_space.hpp"

#include "joint_space_limits.hpp"


namespace ReaK {

namespace pp {
  
  
// se2_0th_order_topology
extern template class metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;

// se2_1st_order_topology
extern template class metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> >, hyperball_topology< vect<double,2> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;

// se2_2nd_order_topology
extern template class metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> >, hyperball_topology< vect<double,2> >, hyperball_topology< vect<double,2> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;



// se2_0th_order_rl_topology
extern template class metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;

// se2_1st_order_rl_topology
extern template class metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> >, hyperball_topology< vect<double,2> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;

// se2_2nd_order_rl_topology
extern template class metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> >, hyperball_topology< vect<double,2> >, hyperball_topology< vect<double,2> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >;

      
      
      
// se2_0th_order_topology
extern template class temporal_space< metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;

// se2_1st_order_topology
extern template class temporal_space< metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> >, hyperball_topology< vect<double,2> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;

// se2_2nd_order_topology
extern template class temporal_space< metric_space_tuple< arithmetic_tuple<
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> >, hyperball_topology< vect<double,2> >, hyperball_topology< vect<double,2> > >, euclidean_tuple_distance >,
      differentiable_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;


// se2_0th_order_rl_topology
extern template class temporal_space< metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;

// se2_1st_order_rl_topology
extern template class temporal_space< metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> >, hyperball_topology< vect<double,2> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;

// se2_2nd_order_rl_topology
extern template class temporal_space< metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> >, hyperball_topology< vect<double,2> >, hyperball_topology< vect<double,2> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, spatial_distance_only>;

      
// se2_0th_order_rl_topology
extern template class temporal_space< metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, reach_plus_time_metric>;

// se2_1st_order_rl_topology
extern template class temporal_space< metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> >, hyperball_topology< vect<double,2> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, reach_plus_time_metric>;

// se2_2nd_order_rl_topology
extern template class temporal_space< metric_space_tuple< arithmetic_tuple<
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect<double,2> >, hyperball_topology< vect<double,2> >, hyperball_topology< vect<double,2> > >, euclidean_tuple_distance >,
      reach_time_diff_space< time_topology, arithmetic_tuple< line_segment_topology<double>, line_segment_topology<double>, line_segment_topology<double> >, euclidean_tuple_distance > >, euclidean_tuple_distance >, time_poisson_topology, reach_plus_time_metric>;





extern template metric_space_array< se2_0th_order_rl_topology<double>::type, 1>::type joint_limits_collection<double>::make_rl_joint_space(const metric_space_array< se2_0th_order_topology<double>::type, 1>::type&) const;
extern template metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type joint_limits_collection<double>::make_rl_joint_space(const metric_space_array< se2_1st_order_topology<double>::type, 1>::type&) const;
extern template metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type joint_limits_collection<double>::make_rl_joint_space(const metric_space_array< se2_2nd_order_topology<double>::type, 1>::type&) const;

extern template metric_space_array< se2_0th_order_topology<double>::type, 1>::type joint_limits_collection<double>::make_normal_joint_space(const metric_space_array< se2_0th_order_rl_topology<double>::type, 1>::type&) const;
extern template metric_space_array< se2_1st_order_topology<double>::type, 1>::type joint_limits_collection<double>::make_normal_joint_space(const metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type&) const;
extern template metric_space_array< se2_2nd_order_topology<double>::type, 1>::type joint_limits_collection<double>::make_normal_joint_space(const metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type&) const;

extern template 
topology_traits< metric_space_array< se2_0th_order_rl_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se2_0th_order_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se2_0th_order_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se2_0th_order_rl_topology<double>::type, 1>::type& ) const;
extern template 
topology_traits< metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se2_1st_order_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se2_1st_order_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type& ) const;
extern template 
topology_traits< metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se2_2nd_order_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se2_2nd_order_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type& ) const;

extern template 
topology_traits< metric_space_array< se2_0th_order_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se2_0th_order_rl_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se2_0th_order_rl_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se2_0th_order_topology<double>::type, 1>::type& ) const;
extern template 
topology_traits< metric_space_array< se2_1st_order_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se2_1st_order_topology<double>::type, 1>::type& ) const;
extern template 
topology_traits< metric_space_array< se2_2nd_order_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se2_2nd_order_topology<double>::type, 1>::type& ) const;




};


extern template frame_2D<double> get_frame_2D(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>,    vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double, double > >& pt);

extern template frame_2D<double> get_frame_2D(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double > >& pt);

extern template frame_2D<double> get_frame_2D(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                          arithmetic_tuple< double > >& pt);


extern template void set_frame_2D(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>,    vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double, double > >& pt,
  const frame_2D<double>& p);

extern template void set_frame_2D(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double > >& pt,
  const frame_2D<double>& p);

extern template void set_frame_2D(
  arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                    arithmetic_tuple< double > >& pt,
  const frame_2D<double>& p);

extern template pose_2D<double> get_pose_2D(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                          arithmetic_tuple< double > >& pt);

extern template void set_pose_2D(
  arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                    arithmetic_tuple< double > >& pt,
  const pose_2D<double>& p);

extern template const double& get_rotation(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double, double > >& pt);

extern template const double& get_rotation(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double > >& pt);

extern template const double& get_rotation(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                          arithmetic_tuple< double > >& pt);

extern template void set_rotation(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double, double > >& pt,
  const double& q);

extern template void set_rotation(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double > >& pt,
  const double& q);

extern template void set_rotation(
  arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                    arithmetic_tuple< double > >& pt,
  const double& q);

extern template const vect<double,2>& get_position(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double, double > >& pt);

extern template const vect<double,2>& get_position(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double > >& pt);

extern template const vect<double,2>& get_position(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                          arithmetic_tuple< double > >& pt);

extern template void set_position(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double, double > >& pt,
  const vect<double,2>& p);

extern template void set_position(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double > >& pt,
  const vect<double,2>& p);

extern template void set_position(
  arithmetic_tuple< arithmetic_tuple< vect<double,2> >,
                    arithmetic_tuple< double > >& pt,
  const vect<double,2>& p);

extern template const double& get_ang_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double, double > >& pt);

extern template const double& get_ang_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double > >& pt);

extern template void set_ang_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double, double > >& pt,
  const double& p);

extern template void set_ang_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double > >& pt,
  const double& p);

extern template const vect<double,2>& get_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double, double > >& pt);

extern template const vect<double,2>& get_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double > >& pt);

extern template void set_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double, double > >& pt,
  const vect<double,2>& p);

extern template void set_velocity(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double > >& pt,
  const vect<double,2>& p);

extern template const double& get_ang_acceleration(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double, double > >& pt);

extern template void set_ang_acceleration(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double, double > >& pt,
  const double& p);

extern template const vect<double,2>& get_acceleration(
  const arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                          arithmetic_tuple< double, double, double > >& pt);

extern template void set_acceleration(
  arithmetic_tuple< arithmetic_tuple< vect<double,2>, vect<double,2>, vect<double,2> >,
                    arithmetic_tuple< double, double, double > >& pt,
  const vect<double,2>& p);

};

#endif



#endif









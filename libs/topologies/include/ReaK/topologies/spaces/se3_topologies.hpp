/**
 * \file se3_topologies.hpp
 *
 * This library provides classes that define topologies on SE(3) (3D rigid-body motion).
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

#ifndef REAK_SE3_TOPOLOGIES_HPP
#define REAK_SE3_TOPOLOGIES_HPP


#include <ReaK/core/base/defs.hpp>

#include "so3_topologies.hpp"

#include "differentiable_space.hpp"
#include "rate_limited_spaces.hpp"
#include "metric_space_tuple.hpp"

#include "hyperbox_topology.hpp"
#include "hyperball_topology.hpp"
#include "line_topology.hpp"

#include <ReaK/math/lin_alg/arithmetic_tuple.hpp>
#include <ReaK/math/lin_alg/vect_alg.hpp>
#include <ReaK/math/kinetostatics/frame_3D.hpp>

namespace ReaK {

namespace pp {


/**
 * This meta-function defines the type for a 0th order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template < typename T, typename DistanceMetric = euclidean_tuple_distance >
struct se3_0th_order_topology {
  typedef metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                 3 > > >,
                                                                      DistanceMetric >,
                                                differentiable_space< time_topology,
                                                                      arithmetic_tuple< quaternion_topology< T > >,
                                                                      DistanceMetric > >,
                              DistanceMetric > type;
};


template < typename T >
typename se3_0th_order_topology< T >::type make_se3_space( const std::string& aName, const vect< T, 3 >& aMinCorner,
                                                           const vect< T, 3 >& aMaxCorner ) {

  return metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                     arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                3 > > >,
                                                                     euclidean_tuple_distance >,
                                               differentiable_space< time_topology,
                                                                     arithmetic_tuple< quaternion_topology< T > >,
                                                                     euclidean_tuple_distance > >,
                             euclidean_tuple_distance >(
    arithmetic_tuple< differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > > >,
                                            euclidean_tuple_distance >,
                      differentiable_space< time_topology, arithmetic_tuple< quaternion_topology< T > >,
                                            euclidean_tuple_distance > >(
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > > >,
                            euclidean_tuple_distance >( arithmetic_tuple< hyperbox_topology< vect< T, 3 > > >(
        hyperbox_topology< vect< T, 3 > >( aName + "_pos", aMinCorner, aMaxCorner ) ) ),
      differentiable_space< time_topology, arithmetic_tuple< quaternion_topology< T > >, euclidean_tuple_distance >(
        arithmetic_tuple< quaternion_topology< T > >( quaternion_topology< T >( aName + "_quat" ) ) ) ) );
};


template < typename TupleDistanceMetric, typename T >
typename se3_0th_order_topology< T, TupleDistanceMetric >::type
  make_se3_space( const std::string& aName, const vect< T, 3 >& aMinCorner, const vect< T, 3 >& aMaxCorner ) {

  return metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                     arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                3 > > >,
                                                                     TupleDistanceMetric >,
                                               differentiable_space< time_topology,
                                                                     arithmetic_tuple< quaternion_topology< T > >,
                                                                     TupleDistanceMetric > >,
                             TupleDistanceMetric >(
    arithmetic_tuple< differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > > >,
                                            TupleDistanceMetric >,
                      differentiable_space< time_topology, arithmetic_tuple< quaternion_topology< T > >,
                                            TupleDistanceMetric > >(
      differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > > >, TupleDistanceMetric >(
        arithmetic_tuple< hyperbox_topology< vect< T, 3 > > >(
          hyperbox_topology< vect< T, 3 > >( aName + "_pos", aMinCorner, aMaxCorner ) ) ),
      differentiable_space< time_topology, arithmetic_tuple< quaternion_topology< T > >, TupleDistanceMetric >(
        arithmetic_tuple< quaternion_topology< T > >( quaternion_topology< T >( aName + "_quat" ) ) ) ) );
};


/**
 * This meta-function defines the type for a 1st order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template < typename T, typename DistanceMetric = euclidean_tuple_distance >
struct se3_1st_order_topology {
  typedef metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                 3 > >,
                                                                                        hyperball_topology< vect< T,
                                                                                                                  3 > > >,
                                                                      DistanceMetric >,
                                                differentiable_space< time_topology,
                                                                      arithmetic_tuple< quaternion_topology< T >,
                                                                                        ang_velocity_3D_topology< T > >,
                                                                      DistanceMetric > >,
                              DistanceMetric > type;
};


template < typename T >
typename se3_1st_order_topology< T >::type make_se3_space( const std::string& aName, const vect< T, 3 >& aMinCorner,
                                                           const vect< T, 3 >& aMaxCorner, const T& aMaxSpeed,
                                                           const T& aMaxAngularSpeed ) {

  return metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                     arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                3 > >,
                                                                                       hyperball_topology< vect< T,
                                                                                                                 3 > > >,
                                                                     euclidean_tuple_distance >,
                                               differentiable_space< time_topology,
                                                                     arithmetic_tuple< quaternion_topology< T >,
                                                                                       ang_velocity_3D_topology< T > >,
                                                                     euclidean_tuple_distance > >,
                             euclidean_tuple_distance >(
    arithmetic_tuple< differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > >,
                                                                             hyperball_topology< vect< T, 3 > > >,
                                            euclidean_tuple_distance >,
                      differentiable_space< time_topology,
                                            arithmetic_tuple< quaternion_topology< T >, ang_velocity_3D_topology< T > >,
                                            euclidean_tuple_distance > >(
      differentiable_space< time_topology,
                            arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > > >,
                            euclidean_tuple_distance >(
        arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > > >(
          hyperbox_topology< vect< T, 3 > >( aName + "_pos", aMinCorner, aMaxCorner ),
          hyperball_topology< vect< T, 3 > >( aName + "_vel", vect< T, 3 >( 0.0, 0.0, 0.0 ), aMaxSpeed ) ) ),
      differentiable_space< time_topology, arithmetic_tuple< quaternion_topology< T >, ang_velocity_3D_topology< T > >,
                            euclidean_tuple_distance >(
        arithmetic_tuple< quaternion_topology< T >, ang_velocity_3D_topology< T > >(
          quaternion_topology< T >( aName + "_quat" ),
          ang_velocity_3D_topology< T >( aName + "_ang_vel", aMaxAngularSpeed ) ) ) ) );
};


template < typename TupleDistanceMetric, typename T >
typename se3_1st_order_topology< T, TupleDistanceMetric >::type
  make_se3_space( const std::string& aName, const vect< T, 3 >& aMinCorner, const vect< T, 3 >& aMaxCorner,
                  const T& aMaxSpeed, const T& aMaxAngularSpeed ) {

  return metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                     arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                3 > >,
                                                                                       hyperball_topology< vect< T,
                                                                                                                 3 > > >,
                                                                     TupleDistanceMetric >,
                                               differentiable_space< time_topology,
                                                                     arithmetic_tuple< quaternion_topology< T >,
                                                                                       ang_velocity_3D_topology< T > >,
                                                                     TupleDistanceMetric > >,
                             TupleDistanceMetric >(
    arithmetic_tuple< differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > >,
                                                                             hyperball_topology< vect< T, 3 > > >,
                                            TupleDistanceMetric >,
                      differentiable_space< time_topology,
                                            arithmetic_tuple< quaternion_topology< T >, ang_velocity_3D_topology< T > >,
                                            TupleDistanceMetric > >(
      differentiable_space< time_topology,
                            arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > > >,
                            TupleDistanceMetric >(
        arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > > >(
          hyperbox_topology< vect< T, 3 > >( aName + "_pos", aMinCorner, aMaxCorner ),
          hyperball_topology< vect< T, 3 > >( aName + "_vel", vect< T, 3 >( 0.0, 0.0, 0.0 ), aMaxSpeed ) ) ),
      differentiable_space< time_topology, arithmetic_tuple< quaternion_topology< T >, ang_velocity_3D_topology< T > >,
                            TupleDistanceMetric >(
        arithmetic_tuple< quaternion_topology< T >, ang_velocity_3D_topology< T > >(
          quaternion_topology< T >( aName + "_quat" ),
          ang_velocity_3D_topology< T >( aName + "_ang_vel", aMaxAngularSpeed ) ) ) ) );
};


/**
 * This meta-function defines the type for a 2nd order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template < typename T, typename DistanceMetric = euclidean_tuple_distance >
struct se3_2nd_order_topology {
  typedef metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                 3 > >,
                                                                                        hyperball_topology< vect< T,
                                                                                                                  3 > >,
                                                                                        hyperball_topology< vect< T,
                                                                                                                  3 > > >,
                                                                      DistanceMetric >,
                                                differentiable_space< time_topology,
                                                                      arithmetic_tuple< quaternion_topology< T >,
                                                                                        ang_velocity_3D_topology< T >,
                                                                                        ang_accel_3D_topology< T > >,
                                                                      DistanceMetric > >,
                              DistanceMetric > type;
};


template < typename T >
typename se3_2nd_order_topology< T >::type make_se3_space( const std::string& aName, const vect< T, 3 >& aMinCorner,
                                                           const vect< T, 3 >& aMaxCorner, const T& aMaxSpeed,
                                                           const T& aMaxAngularSpeed, const T& aMaxAcceleration,
                                                           const T& aMaxAngularAccel ) {

  return metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                     arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                3 > >,
                                                                                       hyperball_topology< vect< T,
                                                                                                                 3 > >,
                                                                                       hyperball_topology< vect< T,
                                                                                                                 3 > > >,
                                                                     euclidean_tuple_distance >,
                                               differentiable_space< time_topology,
                                                                     arithmetic_tuple< quaternion_topology< T >,
                                                                                       ang_velocity_3D_topology< T >,
                                                                                       ang_accel_3D_topology< T > >,
                                                                     euclidean_tuple_distance > >,
                             euclidean_tuple_distance >(
    arithmetic_tuple< differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > >,
                                                                             hyperball_topology< vect< T, 3 > >,
                                                                             hyperball_topology< vect< T, 3 > > >,
                                            euclidean_tuple_distance >,
                      differentiable_space< time_topology,
                                            arithmetic_tuple< quaternion_topology< T >, ang_velocity_3D_topology< T >,
                                                              ang_accel_3D_topology< T > >,
                                            euclidean_tuple_distance > >(
      differentiable_space< time_topology,
                            arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > >,
                                              hyperball_topology< vect< T, 3 > > >,
                            euclidean_tuple_distance >(
        arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > >,
                          hyperball_topology< vect< T, 3 > > >(
          hyperbox_topology< vect< T, 3 > >( aName + "_pos", aMinCorner, aMaxCorner ),
          hyperball_topology< vect< T, 3 > >( aName + "_vel", vect< T, 3 >( 0.0, 0.0, 0.0 ), aMaxSpeed ),
          hyperball_topology< vect< T, 3 > >( aName + "_acc", vect< T, 3 >( 0.0, 0.0, 0.0 ), aMaxAcceleration ) ) ),
      differentiable_space< time_topology, arithmetic_tuple< quaternion_topology< T >, ang_velocity_3D_topology< T >,
                                                             ang_accel_3D_topology< T > >,
                            euclidean_tuple_distance >(
        arithmetic_tuple< quaternion_topology< T >, ang_velocity_3D_topology< T >, ang_accel_3D_topology< T > >(
          quaternion_topology< T >( aName + "_quat" ),
          ang_velocity_3D_topology< T >( aName + "_ang_vel", aMaxAngularSpeed ),
          ang_accel_3D_topology< T >( aName + "_ang_acc", aMaxAngularAccel ) ) ) ) );
};


template < typename TupleDistanceMetric, typename T >
typename se3_2nd_order_topology< T, TupleDistanceMetric >::type
  make_se3_space( const std::string& aName, const vect< T, 3 >& aMinCorner, const vect< T, 3 >& aMaxCorner,
                  const T& aMaxSpeed, const T& aMaxAngularSpeed, const T& aMaxAcceleration,
                  const T& aMaxAngularAccel ) {

  return metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                     arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                3 > >,
                                                                                       hyperball_topology< vect< T,
                                                                                                                 3 > >,
                                                                                       hyperball_topology< vect< T,
                                                                                                                 3 > > >,
                                                                     TupleDistanceMetric >,
                                               differentiable_space< time_topology,
                                                                     arithmetic_tuple< quaternion_topology< T >,
                                                                                       ang_velocity_3D_topology< T >,
                                                                                       ang_accel_3D_topology< T > >,
                                                                     TupleDistanceMetric > >,
                             TupleDistanceMetric >(
    arithmetic_tuple< differentiable_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > >,
                                                                             hyperball_topology< vect< T, 3 > >,
                                                                             hyperball_topology< vect< T, 3 > > >,
                                            TupleDistanceMetric >,
                      differentiable_space< time_topology,
                                            arithmetic_tuple< quaternion_topology< T >, ang_velocity_3D_topology< T >,
                                                              ang_accel_3D_topology< T > >,
                                            TupleDistanceMetric > >(
      differentiable_space< time_topology,
                            arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > >,
                                              hyperball_topology< vect< T, 3 > > >,
                            TupleDistanceMetric >(
        arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > >,
                          hyperball_topology< vect< T, 3 > > >(
          hyperbox_topology< vect< T, 3 > >( aName + "_pos", aMinCorner, aMaxCorner ),
          hyperball_topology< vect< T, 3 > >( aName + "_vel", vect< T, 3 >( 0.0, 0.0, 0.0 ), aMaxSpeed ),
          hyperball_topology< vect< T, 3 > >( aName + "_acc", vect< T, 3 >( 0.0, 0.0, 0.0 ), aMaxAcceleration ) ) ),
      differentiable_space< time_topology, arithmetic_tuple< quaternion_topology< T >, ang_velocity_3D_topology< T >,
                                                             ang_accel_3D_topology< T > >,
                            TupleDistanceMetric >(
        arithmetic_tuple< quaternion_topology< T >, ang_velocity_3D_topology< T >, ang_accel_3D_topology< T > >(
          quaternion_topology< T >( aName + "_quat" ),
          ang_velocity_3D_topology< T >( aName + "_ang_vel", aMaxAngularSpeed ),
          ang_accel_3D_topology< T >( aName + "_ang_acc", aMaxAngularAccel ) ) ) ) );
};


template < typename T, int Order, typename DistanceMetric = euclidean_tuple_distance >
struct se3_topology {
  typedef typename boost::mpl::
    if_< boost::mpl::equal_to< boost::mpl::int_< 0 >, boost::mpl::int_< Order > >,
         typename se3_0th_order_topology< T, DistanceMetric >::type,
         typename boost::mpl::if_< boost::mpl::equal_to< boost::mpl::int_< 1 >, boost::mpl::int_< Order > >,
                                   typename se3_1st_order_topology< T, DistanceMetric >::type,
                                   typename se3_2nd_order_topology< T, DistanceMetric >::type >::type >::type type;
};


template < typename SE3Space >
struct is_se3_space : boost::mpl::false_ {};


template < typename T, typename DistanceMetric >
struct
  is_se3_space< metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                            arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                       3 > > >,
                                                                            DistanceMetric >,
                                                      differentiable_space< time_topology,
                                                                            arithmetic_tuple< quaternion_topology< T > >,
                                                                            DistanceMetric > >,
                                    DistanceMetric > > : boost::mpl::true_ {};

template < typename T, typename DistanceMetric >
struct
  is_se3_space< metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                            arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                       3 > >,
                                                                                              hyperball_topology< vect< T,
                                                                                                                        3 > > >,
                                                                            DistanceMetric >,
                                                      differentiable_space< time_topology,
                                                                            arithmetic_tuple< quaternion_topology< T >,
                                                                                              ang_velocity_3D_topology< T > >,
                                                                            DistanceMetric > >,
                                    DistanceMetric > > : boost::mpl::true_ {};

template < typename T, typename DistanceMetric >
struct
  is_se3_space< metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                            arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                       3 > >,
                                                                                              hyperball_topology< vect< T,
                                                                                                                        3 > >,
                                                                                              hyperball_topology< vect< T,
                                                                                                                        3 > > >,
                                                                            DistanceMetric >,
                                                      differentiable_space< time_topology,
                                                                            arithmetic_tuple< quaternion_topology< T >,
                                                                                              ang_velocity_3D_topology< T >,
                                                                                              ang_accel_3D_topology< T > >,
                                                                            DistanceMetric > >,
                                    DistanceMetric > > : boost::mpl::true_ {};


/**
 * This meta-function defines the type for a 0th order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template < typename T, typename DistanceMetric = euclidean_tuple_distance >
struct se3_0th_order_rl_topology {
  typedef metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                       arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                  3 > > >,
                                                                       DistanceMetric >,
                                                reach_time_diff_space< time_topology,
                                                                       arithmetic_tuple< rate_limited_quat_space< T > >,
                                                                       DistanceMetric > >,
                              DistanceMetric > type;
};

template < typename T >
typename se3_0th_order_rl_topology< T >::type
  make_rl_se3_space( const std::string& aName, const vect< T, 3 >& aMinCorner, const vect< T, 3 >& aMaxCorner,
                     const T& aMaxSpeed, const T& aMaxAngularSpeed ) {

  return metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                 3 > > >,
                                                                      euclidean_tuple_distance >,
                                               reach_time_diff_space< time_topology,
                                                                      arithmetic_tuple< rate_limited_quat_space< T > >,
                                                                      euclidean_tuple_distance > >,
                             euclidean_tuple_distance >(
    arithmetic_tuple< reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > > >,
                                             euclidean_tuple_distance >,
                      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space< T > >,
                                             euclidean_tuple_distance > >(
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > > >,
                             euclidean_tuple_distance >(
        arithmetic_tuple< hyperbox_topology< vect< T, 3 > > >( hyperbox_topology< vect< T, 3 > >(
          aName + "_pos", aMinCorner * ( 1.0 / aMaxSpeed ), aMaxCorner * ( 1.0 / aMaxSpeed ) ) ) ),
      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space< T > >,
                             euclidean_tuple_distance >( arithmetic_tuple< rate_limited_quat_space< T > >(
        rate_limited_quat_space< T >( aName + "_quat", aMaxAngularSpeed ) ) ) ) );
};


template < typename TupleDistanceMetric, typename T >
typename se3_0th_order_rl_topology< T, TupleDistanceMetric >::type
  make_rl_se3_space( const std::string& aName, const vect< T, 3 >& aMinCorner, const vect< T, 3 >& aMaxCorner,
                     const T& aMaxSpeed, const T& aMaxAngularSpeed ) {

  return metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                 3 > > >,
                                                                      TupleDistanceMetric >,
                                               reach_time_diff_space< time_topology,
                                                                      arithmetic_tuple< rate_limited_quat_space< T > >,
                                                                      TupleDistanceMetric > >,
                             TupleDistanceMetric >(
    arithmetic_tuple< reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > > >,
                                             TupleDistanceMetric >,
                      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space< T > >,
                                             TupleDistanceMetric > >(
      reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > > >,
                             TupleDistanceMetric >(
        arithmetic_tuple< hyperbox_topology< vect< T, 3 > > >( hyperbox_topology< vect< T, 3 > >(
          aName + "_pos", aMinCorner * ( 1.0 / aMaxSpeed ), aMaxCorner * ( 1.0 / aMaxSpeed ) ) ) ),
      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space< T > >, TupleDistanceMetric >(
        arithmetic_tuple< rate_limited_quat_space< T > >(
          rate_limited_quat_space< T >( aName + "_quat", aMaxAngularSpeed ) ) ) ) );
};


/**
 * This meta-function defines the type for a 1st order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template < typename T, typename DistanceMetric = euclidean_tuple_distance >
struct se3_1st_order_rl_topology {
  typedef metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                       arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                  3 > >,
                                                                                         hyperball_topology< vect< T,
                                                                                                                   3 > > >,
                                                                       DistanceMetric >,
                                                reach_time_diff_space< time_topology,
                                                                       arithmetic_tuple< rate_limited_quat_space< T >,
                                                                                         ang_velocity_3D_topology< T > >,
                                                                       DistanceMetric > >,
                              DistanceMetric > type;
};


template < typename T >
typename se3_1st_order_rl_topology< T >::type
  make_rl_se3_space( const std::string& aName, const vect< T, 3 >& aMinCorner, const vect< T, 3 >& aMaxCorner,
                     const T& aMaxSpeed, const T& aMaxAngularSpeed, const T& aMaxAcceleration,
                     const T& aMaxAngularAccel ) {

  return metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                 3 > >,
                                                                                        hyperball_topology< vect< T,
                                                                                                                  3 > > >,
                                                                      euclidean_tuple_distance >,
                                               reach_time_diff_space< time_topology,
                                                                      arithmetic_tuple< rate_limited_quat_space< T >,
                                                                                        ang_velocity_3D_topology< T > >,
                                                                      euclidean_tuple_distance > >,
                             euclidean_tuple_distance >(
    arithmetic_tuple< reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > >,
                                                                              hyperball_topology< vect< T, 3 > > >,
                                             euclidean_tuple_distance >,
                      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space< T >,
                                                                              ang_velocity_3D_topology< T > >,
                                             euclidean_tuple_distance > >(
      reach_time_diff_space< time_topology,
                             arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > > >,
                             euclidean_tuple_distance >(
        arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > > >(
          hyperbox_topology< vect< T, 3 > >( aName + "_pos", aMinCorner * ( 1.0 / aMaxSpeed ),
                                             aMaxCorner * ( 1.0 / aMaxSpeed ) ),
          hyperball_topology< vect< T, 3 > >( aName + "_vel", vect< T, 3 >( 0.0, 0.0, 0.0 ),
                                              aMaxSpeed / aMaxAcceleration ) ),
        euclidean_tuple_distance(),
        arithmetic_tuple< reach_time_differentiation >( reach_time_differentiation( aMaxSpeed / aMaxAcceleration ) ) ),
      reach_time_diff_space< time_topology,
                             arithmetic_tuple< rate_limited_quat_space< T >, ang_velocity_3D_topology< T > >,
                             euclidean_tuple_distance >(
        arithmetic_tuple< rate_limited_quat_space< T >, ang_velocity_3D_topology< T > >(
          rate_limited_quat_space< T >( aName + "_quat", aMaxAngularSpeed ),
          ang_velocity_3D_topology< T >( aName + "_ang_vel", aMaxAngularSpeed / aMaxAngularAccel ) ),
        euclidean_tuple_distance(), arithmetic_tuple< reach_time_differentiation >(
                                      reach_time_differentiation( aMaxAngularSpeed / aMaxAngularAccel ) ) ) ) );
};


template < typename TupleDistanceMetric, typename T >
typename se3_1st_order_rl_topology< T, TupleDistanceMetric >::type
  make_rl_se3_space( const std::string& aName, const vect< T, 3 >& aMinCorner, const vect< T, 3 >& aMaxCorner,
                     const T& aMaxSpeed, const T& aMaxAngularSpeed, const T& aMaxAcceleration,
                     const T& aMaxAngularAccel ) {

  return metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                 3 > >,
                                                                                        hyperball_topology< vect< T,
                                                                                                                  3 > > >,
                                                                      TupleDistanceMetric >,
                                               reach_time_diff_space< time_topology,
                                                                      arithmetic_tuple< rate_limited_quat_space< T >,
                                                                                        ang_velocity_3D_topology< T > >,
                                                                      TupleDistanceMetric > >,
                             TupleDistanceMetric >(
    arithmetic_tuple< reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > >,
                                                                              hyperball_topology< vect< T, 3 > > >,
                                             TupleDistanceMetric >,
                      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space< T >,
                                                                              ang_velocity_3D_topology< T > >,
                                             TupleDistanceMetric > >(
      reach_time_diff_space< time_topology,
                             arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > > >,
                             TupleDistanceMetric >(
        arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > > >(
          hyperbox_topology< vect< T, 3 > >( aName + "_pos", aMinCorner * ( 1.0 / aMaxSpeed ),
                                             aMaxCorner * ( 1.0 / aMaxSpeed ) ),
          hyperball_topology< vect< T, 3 > >( aName + "_vel", vect< T, 3 >( 0.0, 0.0, 0.0 ),
                                              aMaxSpeed / aMaxAcceleration ) ),
        TupleDistanceMetric(),
        arithmetic_tuple< reach_time_differentiation >( reach_time_differentiation( aMaxSpeed / aMaxAcceleration ) ) ),
      reach_time_diff_space< time_topology,
                             arithmetic_tuple< rate_limited_quat_space< T >, ang_velocity_3D_topology< T > >,
                             TupleDistanceMetric >(
        arithmetic_tuple< rate_limited_quat_space< T >, ang_velocity_3D_topology< T > >(
          rate_limited_quat_space< T >( aName + "_quat", aMaxAngularSpeed ),
          ang_velocity_3D_topology< T >( aName + "_ang_vel", aMaxAngularSpeed / aMaxAngularAccel ) ),
        TupleDistanceMetric(), arithmetic_tuple< reach_time_differentiation >(
                                 reach_time_differentiation( aMaxAngularSpeed / aMaxAngularAccel ) ) ) ) );
};


/**
 * This meta-function defines the type for a 2nd order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template < typename T, typename DistanceMetric = euclidean_tuple_distance >
struct se3_2nd_order_rl_topology {
  typedef metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                       arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                  3 > >,
                                                                                         hyperball_topology< vect< T,
                                                                                                                   3 > >,
                                                                                         hyperball_topology< vect< T,
                                                                                                                   3 > > >,
                                                                       DistanceMetric >,
                                                reach_time_diff_space< time_topology,
                                                                       arithmetic_tuple< rate_limited_quat_space< T >,
                                                                                         ang_velocity_3D_topology< T >,
                                                                                         ang_accel_3D_topology< T > >,
                                                                       DistanceMetric > >,
                              DistanceMetric > type;
};


template < typename T >
typename se3_2nd_order_rl_topology< T >::type
  make_rl_se3_space( const std::string& aName, const vect< T, 3 >& aMinCorner, const vect< T, 3 >& aMaxCorner,
                     const T& aMaxSpeed, const T& aMaxAngularSpeed, const T& aMaxAcceleration,
                     const T& aMaxAngularAccel, const T& aMaxJerk, const T& aMaxAngularJerk ) {

  return metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                 3 > >,
                                                                                        hyperball_topology< vect< T,
                                                                                                                  3 > >,
                                                                                        hyperball_topology< vect< T,
                                                                                                                  3 > > >,
                                                                      euclidean_tuple_distance >,
                                               reach_time_diff_space< time_topology,
                                                                      arithmetic_tuple< rate_limited_quat_space< T >,
                                                                                        ang_velocity_3D_topology< T >,
                                                                                        ang_accel_3D_topology< T > >,
                                                                      euclidean_tuple_distance > >,
                             euclidean_tuple_distance >(
    arithmetic_tuple< reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > >,
                                                                              hyperball_topology< vect< T, 3 > >,
                                                                              hyperball_topology< vect< T, 3 > > >,
                                             euclidean_tuple_distance >,
                      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space< T >,
                                                                              ang_velocity_3D_topology< T >,
                                                                              ang_accel_3D_topology< T > >,
                                             euclidean_tuple_distance > >(
      reach_time_diff_space< time_topology,
                             arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > >,
                                               hyperball_topology< vect< T, 3 > > >,
                             euclidean_tuple_distance >(
        arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > >,
                          hyperball_topology< vect< T, 3 > > >(
          hyperbox_topology< vect< T, 3 > >( aName + "_pos", aMinCorner * ( 1.0 / aMaxSpeed ),
                                             aMaxCorner * ( 1.0 / aMaxSpeed ) ),
          hyperball_topology< vect< T, 3 > >( aName + "_vel", vect< T, 3 >( 0.0, 0.0, 0.0 ),
                                              aMaxSpeed / aMaxAcceleration ),
          hyperball_topology< vect< T, 3 > >( aName + "_acc", vect< T, 3 >( 0.0, 0.0, 0.0 ),
                                              aMaxAcceleration / aMaxJerk ) ),
        euclidean_tuple_distance(), arithmetic_tuple< reach_time_differentiation, reach_time_differentiation >(
                                      reach_time_differentiation( aMaxSpeed / aMaxAcceleration ),
                                      reach_time_differentiation( aMaxAcceleration / aMaxJerk ) ) ),
      reach_time_diff_space< time_topology,
                             arithmetic_tuple< rate_limited_quat_space< T >, ang_velocity_3D_topology< T >,
                                               ang_accel_3D_topology< T > >,
                             euclidean_tuple_distance >(
        arithmetic_tuple< rate_limited_quat_space< T >, ang_velocity_3D_topology< T >, ang_accel_3D_topology< T > >(
          rate_limited_quat_space< T >( aName + "_quat", aMaxAngularSpeed ),
          ang_velocity_3D_topology< T >( aName + "_ang_vel", aMaxAngularSpeed / aMaxAngularAccel ),
          ang_accel_3D_topology< T >( aName + "_ang_acc", aMaxAngularAccel / aMaxAngularJerk ) ),
        euclidean_tuple_distance(), arithmetic_tuple< reach_time_differentiation, reach_time_differentiation >(
                                      reach_time_differentiation( aMaxAngularSpeed / aMaxAngularAccel ),
                                      reach_time_differentiation( aMaxAngularAccel / aMaxAngularJerk ) ) ) ) );
};


template < typename TupleDistanceMetric, typename T >
typename se3_2nd_order_rl_topology< T, TupleDistanceMetric >::type
  make_rl_se3_space( const std::string& aName, const vect< T, 3 >& aMinCorner, const vect< T, 3 >& aMaxCorner,
                     const T& aMaxSpeed, const T& aMaxAngularSpeed, const T& aMaxAcceleration,
                     const T& aMaxAngularAccel, const T& aMaxJerk, const T& aMaxAngularJerk ) {

  return metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                 3 > >,
                                                                                        hyperball_topology< vect< T,
                                                                                                                  3 > >,
                                                                                        hyperball_topology< vect< T,
                                                                                                                  3 > > >,
                                                                      TupleDistanceMetric >,
                                               reach_time_diff_space< time_topology,
                                                                      arithmetic_tuple< rate_limited_quat_space< T >,
                                                                                        ang_velocity_3D_topology< T >,
                                                                                        ang_accel_3D_topology< T > >,
                                                                      TupleDistanceMetric > >,
                             TupleDistanceMetric >(
    arithmetic_tuple< reach_time_diff_space< time_topology, arithmetic_tuple< hyperbox_topology< vect< T, 3 > >,
                                                                              hyperball_topology< vect< T, 3 > >,
                                                                              hyperball_topology< vect< T, 3 > > >,
                                             TupleDistanceMetric >,
                      reach_time_diff_space< time_topology, arithmetic_tuple< rate_limited_quat_space< T >,
                                                                              ang_velocity_3D_topology< T >,
                                                                              ang_accel_3D_topology< T > >,
                                             TupleDistanceMetric > >(
      reach_time_diff_space< time_topology,
                             arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > >,
                                               hyperball_topology< vect< T, 3 > > >,
                             TupleDistanceMetric >(
        arithmetic_tuple< hyperbox_topology< vect< T, 3 > >, hyperball_topology< vect< T, 3 > >,
                          hyperball_topology< vect< T, 3 > > >(
          hyperbox_topology< vect< T, 3 > >( aName + "_pos", aMinCorner * ( 1.0 / aMaxSpeed ),
                                             aMaxCorner * ( 1.0 / aMaxSpeed ) ),
          hyperball_topology< vect< T, 3 > >( aName + "_vel", vect< T, 3 >( 0.0, 0.0, 0.0 ),
                                              aMaxSpeed / aMaxAcceleration ),
          hyperball_topology< vect< T, 3 > >( aName + "_acc", vect< T, 3 >( 0.0, 0.0, 0.0 ),
                                              aMaxAcceleration / aMaxJerk ) ),
        TupleDistanceMetric(), arithmetic_tuple< reach_time_differentiation, reach_time_differentiation >(
                                 reach_time_differentiation( aMaxSpeed / aMaxAcceleration ),
                                 reach_time_differentiation( aMaxAcceleration / aMaxJerk ) ) ),
      reach_time_diff_space< time_topology,
                             arithmetic_tuple< rate_limited_quat_space< T >, ang_velocity_3D_topology< T >,
                                               ang_accel_3D_topology< T > >,
                             TupleDistanceMetric >(
        arithmetic_tuple< rate_limited_quat_space< T >, ang_velocity_3D_topology< T >, ang_accel_3D_topology< T > >(
          rate_limited_quat_space< T >( aName + "_quat", aMaxAngularSpeed ),
          ang_velocity_3D_topology< T >( aName + "_ang_vel", aMaxAngularSpeed / aMaxAngularAccel ),
          ang_accel_3D_topology< T >( aName + "_ang_acc", aMaxAngularAccel / aMaxAngularJerk ) ),
        TupleDistanceMetric(), arithmetic_tuple< reach_time_differentiation, reach_time_differentiation >(
                                 reach_time_differentiation( aMaxAngularSpeed / aMaxAngularAccel ),
                                 reach_time_differentiation( aMaxAngularAccel / aMaxAngularJerk ) ) ) ) );
};


template < typename T, int Order, typename DistanceMetric = euclidean_tuple_distance >
struct se3_rl_topology {
  typedef typename boost::mpl::
    if_< boost::mpl::equal_to< boost::mpl::int_< 0 >, boost::mpl::int_< Order > >,
         typename se3_0th_order_rl_topology< T, DistanceMetric >::type,
         typename boost::mpl::if_< boost::mpl::equal_to< boost::mpl::int_< 1 >, boost::mpl::int_< Order > >,
                                   typename se3_1st_order_rl_topology< T, DistanceMetric >::type,
                                   typename se3_2nd_order_rl_topology< T, DistanceMetric >::type >::type >::type type;
};


template < typename SE3Space >
struct is_rate_limited_se3_space : boost::mpl::false_ {};


template < typename T, typename DistanceMetric >
struct
  is_rate_limited_se3_space< metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                                          arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                                     3 > > >,
                                                                                          DistanceMetric >,
                                                                   reach_time_diff_space< time_topology,
                                                                                          arithmetic_tuple< rate_limited_quat_space< T > >,
                                                                                          DistanceMetric > >,
                                                 DistanceMetric > > : boost::mpl::true_ {};

template < typename T, typename DistanceMetric >
struct
  is_rate_limited_se3_space< metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                                          arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                                     3 > >,
                                                                                                            hyperball_topology< vect< T,
                                                                                                                                      3 > > >,
                                                                                          DistanceMetric >,
                                                                   reach_time_diff_space< time_topology,
                                                                                          arithmetic_tuple< rate_limited_quat_space< T >,
                                                                                                            ang_velocity_3D_topology< T > >,
                                                                                          DistanceMetric > >,
                                                 DistanceMetric > > : boost::mpl::true_ {};

template < typename T, typename DistanceMetric >
struct
  is_rate_limited_se3_space< metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                                          arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                                     3 > >,
                                                                                                            hyperball_topology< vect< T,
                                                                                                                                      3 > >,
                                                                                                            hyperball_topology< vect< T,
                                                                                                                                      3 > > >,
                                                                                          DistanceMetric >,
                                                                   reach_time_diff_space< time_topology,
                                                                                          arithmetic_tuple< rate_limited_quat_space< T >,
                                                                                                            ang_velocity_3D_topology< T >,
                                                                                                            ang_accel_3D_topology< T > >,
                                                                                          DistanceMetric > >,
                                                 DistanceMetric > > : boost::mpl::true_ {};


template < typename T, typename DistanceMetric >
struct
  get_rate_limited_space< metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                                 3 > > >,
                                                                                      DistanceMetric >,
                                                                differentiable_space< time_topology,
                                                                                      arithmetic_tuple< quaternion_topology< T > >,
                                                                                      DistanceMetric > >,
                                              DistanceMetric > > {
  typedef metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                       arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                  3 > > >,
                                                                       DistanceMetric >,
                                                reach_time_diff_space< time_topology,
                                                                       arithmetic_tuple< rate_limited_quat_space< T > >,
                                                                       DistanceMetric > >,
                              DistanceMetric > type;
};

template < typename T, typename DistanceMetric >
struct
  get_rate_limited_space< metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                                 3 > >,
                                                                                                        hyperball_topology< vect< T,
                                                                                                                                  3 > > >,
                                                                                      DistanceMetric >,
                                                                differentiable_space< time_topology,
                                                                                      arithmetic_tuple< quaternion_topology< T >,
                                                                                                        ang_velocity_3D_topology< T > >,
                                                                                      DistanceMetric > >,
                                              DistanceMetric > > {
  typedef metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                       arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                  3 > >,
                                                                                         hyperball_topology< vect< T,
                                                                                                                   3 > > >,
                                                                       DistanceMetric >,
                                                reach_time_diff_space< time_topology,
                                                                       arithmetic_tuple< rate_limited_quat_space< T >,
                                                                                         ang_velocity_3D_topology< T > >,
                                                                       DistanceMetric > >,
                              DistanceMetric > type;
};

template < typename T, typename DistanceMetric >
struct
  get_rate_limited_space< metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                                 3 > >,
                                                                                                        hyperball_topology< vect< T,
                                                                                                                                  3 > >,
                                                                                                        hyperball_topology< vect< T,
                                                                                                                                  3 > > >,
                                                                                      DistanceMetric >,
                                                                differentiable_space< time_topology,
                                                                                      arithmetic_tuple< quaternion_topology< T >,
                                                                                                        ang_velocity_3D_topology< T >,
                                                                                                        ang_accel_3D_topology< T > >,
                                                                                      DistanceMetric > >,
                                              DistanceMetric > > {
  typedef metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                       arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                  3 > >,
                                                                                         hyperball_topology< vect< T,
                                                                                                                   3 > >,
                                                                                         hyperball_topology< vect< T,
                                                                                                                   3 > > >,
                                                                       DistanceMetric >,
                                                reach_time_diff_space< time_topology,
                                                                       arithmetic_tuple< rate_limited_quat_space< T >,
                                                                                         ang_velocity_3D_topology< T >,
                                                                                         ang_accel_3D_topology< T > >,
                                                                       DistanceMetric > >,
                              DistanceMetric > type;
};


template < typename T, typename DistanceMetric >
struct
  get_rate_illimited_space< metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                                         arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                                    3 > > >,
                                                                                         DistanceMetric >,
                                                                  reach_time_diff_space< time_topology,
                                                                                         arithmetic_tuple< rate_limited_quat_space< T > >,
                                                                                         DistanceMetric > >,
                                                DistanceMetric > > {
  typedef metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                 3 > > >,
                                                                      DistanceMetric >,
                                                differentiable_space< time_topology,
                                                                      arithmetic_tuple< quaternion_topology< T > >,
                                                                      DistanceMetric > >,
                              DistanceMetric > type;
};

template < typename T, typename DistanceMetric >
struct
  get_rate_illimited_space< metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                                         arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                                    3 > >,
                                                                                                           hyperball_topology< vect< T,
                                                                                                                                     3 > > >,
                                                                                         DistanceMetric >,
                                                                  reach_time_diff_space< time_topology,
                                                                                         arithmetic_tuple< rate_limited_quat_space< T >,
                                                                                                           ang_velocity_3D_topology< T > >,
                                                                                         DistanceMetric > >,
                                                DistanceMetric > > {
  typedef metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                 3 > >,
                                                                                        hyperball_topology< vect< T,
                                                                                                                  3 > > >,
                                                                      DistanceMetric >,
                                                differentiable_space< time_topology,
                                                                      arithmetic_tuple< quaternion_topology< T >,
                                                                                        ang_velocity_3D_topology< T > >,
                                                                      DistanceMetric > >,
                              DistanceMetric > type;
};

template < typename T, typename DistanceMetric >
struct
  get_rate_illimited_space< metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                                         arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                                    3 > >,
                                                                                                           hyperball_topology< vect< T,
                                                                                                                                     3 > >,
                                                                                                           hyperball_topology< vect< T,
                                                                                                                                     3 > > >,
                                                                                         DistanceMetric >,
                                                                  reach_time_diff_space< time_topology,
                                                                                         arithmetic_tuple< rate_limited_quat_space< T >,
                                                                                                           ang_velocity_3D_topology< T >,
                                                                                                           ang_accel_3D_topology< T > >,
                                                                                         DistanceMetric > >,
                                                DistanceMetric > > {
  typedef metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                      arithmetic_tuple< hyperbox_topology< vect< T,
                                                                                                                 3 > >,
                                                                                        hyperball_topology< vect< T,
                                                                                                                  3 > >,
                                                                                        hyperball_topology< vect< T,
                                                                                                                  3 > > >,
                                                                      DistanceMetric >,
                                                differentiable_space< time_topology,
                                                                      arithmetic_tuple< quaternion_topology< T >,
                                                                                        ang_velocity_3D_topology< T >,
                                                                                        ang_accel_3D_topology< T > >,
                                                                      DistanceMetric > >,
                              DistanceMetric > type;
};
};


// Because of ADL rules, the get functions for the arithmetic-tuple types that represent SE(3) states should be in the
// ReaK namespace.


template < typename T >
frame_3D< T >
  get_frame_3D( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 >, vect< T, 3 > >,
                                        arithmetic_tuple< unit_quat< T >, vect< T, 3 >, vect< T, 3 > > >& pt ) {
  return frame_3D< T >( weak_ptr< pose_3D< T > >(), get< 0 >( get< 0 >( pt ) ),
                        quaternion< T >( get< 0 >( get< 1 >( pt ) ) ), get< 1 >( get< 0 >( pt ) ),
                        get< 1 >( get< 1 >( pt ) ), get< 2 >( get< 0 >( pt ) ), get< 2 >( get< 1 >( pt ) ),
                        vect< T, 3 >( 0.0, 0.0, 0.0 ), vect< T, 3 >( 0.0, 0.0, 0.0 ) );
};

template < typename T >
frame_3D< T > get_frame_3D( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 > >,
                                                    arithmetic_tuple< unit_quat< T >, vect< T, 3 > > >& pt ) {
  return frame_3D< T >( weak_ptr< pose_3D< T > >(), get< 0 >( get< 0 >( pt ) ),
                        quaternion< T >( get< 0 >( get< 1 >( pt ) ) ), get< 1 >( get< 0 >( pt ) ),
                        get< 1 >( get< 1 >( pt ) ), vect< T, 3 >( 0.0, 0.0, 0.0 ), vect< T, 3 >( 0.0, 0.0, 0.0 ),
                        vect< T, 3 >( 0.0, 0.0, 0.0 ), vect< T, 3 >( 0.0, 0.0, 0.0 ) );
};

template < typename T >
frame_3D< T >
  get_frame_3D( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 > >, arithmetic_tuple< unit_quat< T > > >& pt ) {
  return frame_3D< T >( weak_ptr< pose_3D< T > >(), get< 0 >( get< 0 >( pt ) ),
                        quaternion< T >( get< 0 >( get< 1 >( pt ) ) ), vect< T, 3 >( 0.0, 0.0, 0.0 ),
                        vect< T, 3 >( 0.0, 0.0, 0.0 ), vect< T, 3 >( 0.0, 0.0, 0.0 ), vect< T, 3 >( 0.0, 0.0, 0.0 ),
                        vect< T, 3 >( 0.0, 0.0, 0.0 ), vect< T, 3 >( 0.0, 0.0, 0.0 ) );
};


template < typename T >
void set_frame_3D( arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 >, vect< T, 3 > >,
                                     arithmetic_tuple< unit_quat< T >, vect< T, 3 >, vect< T, 3 > > >& pt,
                   const frame_3D< T >& p ) {
  get< 0 >( get< 0 >( pt ) ) = p.Position;
  get< 0 >( get< 1 >( pt ) ) = unit_quat< T >( p.Quat );
  get< 1 >( get< 0 >( pt ) ) = p.Velocity;
  get< 1 >( get< 1 >( pt ) ) = p.AngVelocity;
  get< 2 >( get< 0 >( pt ) ) = p.Acceleration;
  get< 2 >( get< 1 >( pt ) ) = p.AngAcceleration;
};

template < typename T >
void set_frame_3D( arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 > >,
                                     arithmetic_tuple< unit_quat< T >, vect< T, 3 > > >& pt,
                   const frame_3D< T >& p ) {
  get< 0 >( get< 0 >( pt ) ) = p.Position;
  get< 0 >( get< 1 >( pt ) ) = unit_quat< T >( p.Quat );
  get< 1 >( get< 0 >( pt ) ) = p.Velocity;
  get< 1 >( get< 1 >( pt ) ) = p.AngVelocity;
};

template < typename T >
void set_frame_3D( arithmetic_tuple< arithmetic_tuple< vect< T, 3 > >, arithmetic_tuple< unit_quat< T > > >& pt,
                   const frame_3D< T >& p ) {
  get< 0 >( get< 0 >( pt ) ) = p.Position;
  get< 0 >( get< 1 >( pt ) ) = unit_quat< T >( p.Quat );
};


template < typename T >
pose_3D< T >
  get_pose_3D( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 > >, arithmetic_tuple< unit_quat< T > > >& pt ) {
  return pose_3D< T >( weak_ptr< pose_3D< T > >(), get< 0 >( get< 0 >( pt ) ),
                       quaternion< T >( get< 0 >( get< 1 >( pt ) ) ) );
};

template < typename T >
void set_pose_3D( arithmetic_tuple< arithmetic_tuple< vect< T, 3 > >, arithmetic_tuple< unit_quat< T > > >& pt,
                  const pose_3D< T >& p ) {
  get< 0 >( get< 0 >( pt ) ) = p.Position;
  get< 0 >( get< 1 >( pt ) ) = unit_quat< T >( p.Quat );
};


template < typename T >
const unit_quat< T >&
  get_quaternion( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 >, vect< T, 3 > >,
                                          arithmetic_tuple< unit_quat< T >, vect< T, 3 >, vect< T, 3 > > >& pt ) {
  return get< 0 >( get< 1 >( pt ) );
};

template < typename T >
const unit_quat< T >& get_quaternion( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 > >,
                                                              arithmetic_tuple< unit_quat< T >, vect< T, 3 > > >& pt ) {
  return get< 0 >( get< 1 >( pt ) );
};

template < typename T >
const unit_quat< T >&
  get_quaternion( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 > >, arithmetic_tuple< unit_quat< T > > >& pt ) {
  return get< 0 >( get< 1 >( pt ) );
};

template < typename T >
void set_quaternion( arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 >, vect< T, 3 > >,
                                       arithmetic_tuple< unit_quat< T >, vect< T, 3 >, vect< T, 3 > > >& pt,
                     const unit_quat< T >& q ) {
  get< 0 >( get< 1 >( pt ) ) = q;
};

template < typename T >
void set_quaternion( arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 > >,
                                       arithmetic_tuple< unit_quat< T >, vect< T, 3 > > >& pt,
                     const unit_quat< T >& q ) {
  get< 0 >( get< 1 >( pt ) ) = q;
};

template < typename T >
void set_quaternion( arithmetic_tuple< arithmetic_tuple< vect< T, 3 > >, arithmetic_tuple< unit_quat< T > > >& pt,
                     const unit_quat< T >& q ) {
  get< 0 >( get< 1 >( pt ) ) = q;
};


template < typename T >
const vect< T, 3 >&
  get_position( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 >, vect< T, 3 > >,
                                        arithmetic_tuple< unit_quat< T >, vect< T, 3 >, vect< T, 3 > > >& pt ) {
  return get< 0 >( get< 0 >( pt ) );
};

template < typename T >
const vect< T, 3 >& get_position( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 > >,
                                                          arithmetic_tuple< unit_quat< T >, vect< T, 3 > > >& pt ) {
  return get< 0 >( get< 0 >( pt ) );
};

template < typename T >
const vect< T, 3 >&
  get_position( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 > >, arithmetic_tuple< unit_quat< T > > >& pt ) {
  return get< 0 >( get< 0 >( pt ) );
};

template < typename T >
void set_position( arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 >, vect< T, 3 > >,
                                     arithmetic_tuple< unit_quat< T >, vect< T, 3 >, vect< T, 3 > > >& pt,
                   const vect< T, 3 >& p ) {
  get< 0 >( get< 0 >( pt ) ) = p;
};

template < typename T >
void set_position( arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 > >,
                                     arithmetic_tuple< unit_quat< T >, vect< T, 3 > > >& pt,
                   const vect< T, 3 >& p ) {
  get< 0 >( get< 0 >( pt ) ) = p;
};

template < typename T >
void set_position( arithmetic_tuple< arithmetic_tuple< vect< T, 3 > >, arithmetic_tuple< unit_quat< T > > >& pt,
                   const vect< T, 3 >& p ) {
  get< 0 >( get< 0 >( pt ) ) = p;
};


template < typename T >
const vect< T, 3 >&
  get_ang_velocity( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 >, vect< T, 3 > >,
                                            arithmetic_tuple< unit_quat< T >, vect< T, 3 >, vect< T, 3 > > >& pt ) {
  return get< 1 >( get< 1 >( pt ) );
};

template < typename T >
const vect< T, 3 >& get_ang_velocity( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 > >,
                                                              arithmetic_tuple< unit_quat< T >, vect< T, 3 > > >& pt ) {
  return get< 1 >( get< 1 >( pt ) );
};

template < typename T >
void set_ang_velocity( arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 >, vect< T, 3 > >,
                                         arithmetic_tuple< unit_quat< T >, vect< T, 3 >, vect< T, 3 > > >& pt,
                       const vect< T, 3 >& p ) {
  get< 1 >( get< 1 >( pt ) ) = p;
};

template < typename T >
void set_ang_velocity( arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 > >,
                                         arithmetic_tuple< unit_quat< T >, vect< T, 3 > > >& pt,
                       const vect< T, 3 >& p ) {
  get< 1 >( get< 1 >( pt ) ) = p;
};


template < typename T >
const vect< T, 3 >&
  get_velocity( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 >, vect< T, 3 > >,
                                        arithmetic_tuple< unit_quat< T >, vect< T, 3 >, vect< T, 3 > > >& pt ) {
  return get< 1 >( get< 0 >( pt ) );
};

template < typename T >
const vect< T, 3 >& get_velocity( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 > >,
                                                          arithmetic_tuple< unit_quat< T >, vect< T, 3 > > >& pt ) {
  return get< 1 >( get< 0 >( pt ) );
};

template < typename T >
void set_velocity( arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 >, vect< T, 3 > >,
                                     arithmetic_tuple< unit_quat< T >, vect< T, 3 >, vect< T, 3 > > >& pt,
                   const vect< T, 3 >& p ) {
  get< 1 >( get< 0 >( pt ) ) = p;
};

template < typename T >
void set_velocity( arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 > >,
                                     arithmetic_tuple< unit_quat< T >, vect< T, 3 > > >& pt,
                   const vect< T, 3 >& p ) {
  get< 1 >( get< 0 >( pt ) ) = p;
};


template < typename T >
const vect< T, 3 >&
  get_ang_acceleration( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 >, vect< T, 3 > >,
                                                arithmetic_tuple< unit_quat< T >, vect< T, 3 >, vect< T, 3 > > >& pt ) {
  return get< 2 >( get< 1 >( pt ) );
};

template < typename T >
void set_ang_acceleration( arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 >, vect< T, 3 > >,
                                             arithmetic_tuple< unit_quat< T >, vect< T, 3 >, vect< T, 3 > > >& pt,
                           const vect< T, 3 >& p ) {
  get< 2 >( get< 1 >( pt ) ) = p;
};


template < typename T >
const vect< T, 3 >&
  get_acceleration( const arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 >, vect< T, 3 > >,
                                            arithmetic_tuple< unit_quat< T >, vect< T, 3 >, vect< T, 3 > > >& pt ) {
  return get< 2 >( get< 0 >( pt ) );
};

template < typename T >
void set_acceleration( arithmetic_tuple< arithmetic_tuple< vect< T, 3 >, vect< T, 3 >, vect< T, 3 > >,
                                         arithmetic_tuple< unit_quat< T >, vect< T, 3 >, vect< T, 3 > > >& pt,
                       const vect< T, 3 >& p ) {
  get< 2 >( get< 0 >( pt ) ) = p;
};
};


#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

#include "time_poisson_topology.hpp"
#include "temporal_space.hpp"
#include "reachability_space.hpp"

#include "joint_space_limits.hpp"

namespace ReaK {

namespace pp {


// se3_0th_order_topology
extern template class
  metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                              arithmetic_tuple< hyperbox_topology< vect< double,
                                                                                                         3 > > >,
                                                              euclidean_tuple_distance >,
                                        differentiable_space< time_topology,
                                                              arithmetic_tuple< quaternion_topology< double > >,
                                                              euclidean_tuple_distance > >,
                      euclidean_tuple_distance >;

extern template se3_0th_order_topology< double >::type
  make_se3_space( const std::string& aName, const vect< double, 3 >& aMinCorner, const vect< double, 3 >& aMaxCorner );

// se3_1st_order_topology
extern template class
  metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                              arithmetic_tuple< hyperbox_topology< vect< double, 3 > >,
                                                                                hyperball_topology< vect< double,
                                                                                                          3 > > >,
                                                              euclidean_tuple_distance >,
                                        differentiable_space< time_topology,
                                                              arithmetic_tuple< quaternion_topology< double >,
                                                                                ang_velocity_3D_topology< double > >,
                                                              euclidean_tuple_distance > >,
                      euclidean_tuple_distance >;

extern template se3_1st_order_topology< double >::type
  make_se3_space( const std::string& aName, const vect< double, 3 >& aMinCorner, const vect< double, 3 >& aMaxCorner,
                  const double& aMaxSpeed, const double& aMaxAngularSpeed );

// se3_2nd_order_topology
extern template class
  metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                              arithmetic_tuple< hyperbox_topology< vect< double, 3 > >,
                                                                                hyperball_topology< vect< double, 3 > >,
                                                                                hyperball_topology< vect< double,
                                                                                                          3 > > >,
                                                              euclidean_tuple_distance >,
                                        differentiable_space< time_topology,
                                                              arithmetic_tuple< quaternion_topology< double >,
                                                                                ang_velocity_3D_topology< double >,
                                                                                ang_accel_3D_topology< double > >,
                                                              euclidean_tuple_distance > >,
                      euclidean_tuple_distance >;

extern template se3_2nd_order_topology< double >::type
  make_se3_space( const std::string& aName, const vect< double, 3 >& aMinCorner, const vect< double, 3 >& aMaxCorner,
                  const double& aMaxSpeed, const double& aMaxAngularSpeed, const double& aMaxAcceleration,
                  const double& aMaxAngularAccel );

// se3_0th_order_rl_topology
extern template class
  metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                               arithmetic_tuple< hyperbox_topology< vect< double,
                                                                                                          3 > > >,
                                                               euclidean_tuple_distance >,
                                        reach_time_diff_space< time_topology,
                                                               arithmetic_tuple< rate_limited_quat_space< double > >,
                                                               euclidean_tuple_distance > >,
                      euclidean_tuple_distance >;

extern template se3_0th_order_rl_topology< double >::type
  make_rl_se3_space( const std::string& aName, const vect< double, 3 >& aMinCorner, const vect< double, 3 >& aMaxCorner,
                     const double& aMaxSpeed, const double& aMaxAngularSpeed );

// se3_1st_order_rl_topology
extern template class
  metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                               arithmetic_tuple< hyperbox_topology< vect< double, 3 > >,
                                                                                 hyperball_topology< vect< double,
                                                                                                           3 > > >,
                                                               euclidean_tuple_distance >,
                                        reach_time_diff_space< time_topology,
                                                               arithmetic_tuple< rate_limited_quat_space< double >,
                                                                                 ang_velocity_3D_topology< double > >,
                                                               euclidean_tuple_distance > >,
                      euclidean_tuple_distance >;

extern template se3_1st_order_rl_topology< double >::type
  make_rl_se3_space( const std::string& aName, const vect< double, 3 >& aMinCorner, const vect< double, 3 >& aMaxCorner,
                     const double& aMaxSpeed, const double& aMaxAngularSpeed, const double& aMaxAcceleration,
                     const double& aMaxAngularAccel );

// se3_2nd_order_rl_topology
extern template class
  metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                               arithmetic_tuple< hyperbox_topology< vect< double, 3 > >,
                                                                                 hyperball_topology< vect< double,
                                                                                                           3 > >,
                                                                                 hyperball_topology< vect< double,
                                                                                                           3 > > >,
                                                               euclidean_tuple_distance >,
                                        reach_time_diff_space< time_topology,
                                                               arithmetic_tuple< rate_limited_quat_space< double >,
                                                                                 ang_velocity_3D_topology< double >,
                                                                                 ang_accel_3D_topology< double > >,
                                                               euclidean_tuple_distance > >,
                      euclidean_tuple_distance >;

extern template se3_2nd_order_rl_topology< double >::type
  make_rl_se3_space( const std::string& aName, const vect< double, 3 >& aMinCorner, const vect< double, 3 >& aMaxCorner,
                     const double& aMaxSpeed, const double& aMaxAngularSpeed, const double& aMaxAcceleration,
                     const double& aMaxAngularAccel, const double& aMaxJerk, const double& aMaxAngularJerk );


// se3_0th_order_topology
extern template class
  temporal_space< metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                              arithmetic_tuple< hyperbox_topology< vect< double,
                                                                                                                         3 > > >,
                                                                              euclidean_tuple_distance >,
                                                        differentiable_space< time_topology,
                                                                              arithmetic_tuple< quaternion_topology< double > >,
                                                                              euclidean_tuple_distance > >,
                                      euclidean_tuple_distance >,
                  time_poisson_topology, spatial_distance_only >;

// se3_1st_order_topology
extern template class
  temporal_space< metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                              arithmetic_tuple< hyperbox_topology< vect< double,
                                                                                                                         3 > >,
                                                                                                hyperball_topology< vect< double,
                                                                                                                          3 > > >,
                                                                              euclidean_tuple_distance >,
                                                        differentiable_space< time_topology,
                                                                              arithmetic_tuple< quaternion_topology< double >,
                                                                                                ang_velocity_3D_topology< double > >,
                                                                              euclidean_tuple_distance > >,
                                      euclidean_tuple_distance >,
                  time_poisson_topology, spatial_distance_only >;

// se3_2nd_order_topology
extern template class
  temporal_space< metric_space_tuple< arithmetic_tuple< differentiable_space< time_topology,
                                                                              arithmetic_tuple< hyperbox_topology< vect< double,
                                                                                                                         3 > >,
                                                                                                hyperball_topology< vect< double,
                                                                                                                          3 > >,
                                                                                                hyperball_topology< vect< double,
                                                                                                                          3 > > >,
                                                                              euclidean_tuple_distance >,
                                                        differentiable_space< time_topology,
                                                                              arithmetic_tuple< quaternion_topology< double >,
                                                                                                ang_velocity_3D_topology< double >,
                                                                                                ang_accel_3D_topology< double > >,
                                                                              euclidean_tuple_distance > >,
                                      euclidean_tuple_distance >,
                  time_poisson_topology, spatial_distance_only >;


// se3_0th_order_rl_topology
extern template class
  temporal_space< metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                               arithmetic_tuple< hyperbox_topology< vect< double,
                                                                                                                          3 > > >,
                                                                               euclidean_tuple_distance >,
                                                        reach_time_diff_space< time_topology,
                                                                               arithmetic_tuple< rate_limited_quat_space< double > >,
                                                                               euclidean_tuple_distance > >,
                                      euclidean_tuple_distance >,
                  time_poisson_topology, spatial_distance_only >;

// se3_1st_order_rl_topology
extern template class
  temporal_space< metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                               arithmetic_tuple< hyperbox_topology< vect< double,
                                                                                                                          3 > >,
                                                                                                 hyperball_topology< vect< double,
                                                                                                                           3 > > >,
                                                                               euclidean_tuple_distance >,
                                                        reach_time_diff_space< time_topology,
                                                                               arithmetic_tuple< rate_limited_quat_space< double >,
                                                                                                 ang_velocity_3D_topology< double > >,
                                                                               euclidean_tuple_distance > >,
                                      euclidean_tuple_distance >,
                  time_poisson_topology, spatial_distance_only >;

// se3_2nd_order_rl_topology
extern template class
  temporal_space< metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                               arithmetic_tuple< hyperbox_topology< vect< double,
                                                                                                                          3 > >,
                                                                                                 hyperball_topology< vect< double,
                                                                                                                           3 > >,
                                                                                                 hyperball_topology< vect< double,
                                                                                                                           3 > > >,
                                                                               euclidean_tuple_distance >,
                                                        reach_time_diff_space< time_topology,
                                                                               arithmetic_tuple< rate_limited_quat_space< double >,
                                                                                                 ang_velocity_3D_topology< double >,
                                                                                                 ang_accel_3D_topology< double > >,
                                                                               euclidean_tuple_distance > >,
                                      euclidean_tuple_distance >,
                  time_poisson_topology, spatial_distance_only >;


// se3_0th_order_rl_topology
extern template class
  temporal_space< metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                               arithmetic_tuple< hyperbox_topology< vect< double,
                                                                                                                          3 > > >,
                                                                               euclidean_tuple_distance >,
                                                        reach_time_diff_space< time_topology,
                                                                               arithmetic_tuple< rate_limited_quat_space< double > >,
                                                                               euclidean_tuple_distance > >,
                                      euclidean_tuple_distance >,
                  time_poisson_topology, reach_plus_time_metric >;

// se3_1st_order_rl_topology
extern template class
  temporal_space< metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                               arithmetic_tuple< hyperbox_topology< vect< double,
                                                                                                                          3 > >,
                                                                                                 hyperball_topology< vect< double,
                                                                                                                           3 > > >,
                                                                               euclidean_tuple_distance >,
                                                        reach_time_diff_space< time_topology,
                                                                               arithmetic_tuple< rate_limited_quat_space< double >,
                                                                                                 ang_velocity_3D_topology< double > >,
                                                                               euclidean_tuple_distance > >,
                                      euclidean_tuple_distance >,
                  time_poisson_topology, reach_plus_time_metric >;

// se3_2nd_order_rl_topology
extern template class
  temporal_space< metric_space_tuple< arithmetic_tuple< reach_time_diff_space< time_topology,
                                                                               arithmetic_tuple< hyperbox_topology< vect< double,
                                                                                                                          3 > >,
                                                                                                 hyperball_topology< vect< double,
                                                                                                                           3 > >,
                                                                                                 hyperball_topology< vect< double,
                                                                                                                           3 > > >,
                                                                               euclidean_tuple_distance >,
                                                        reach_time_diff_space< time_topology,
                                                                               arithmetic_tuple< rate_limited_quat_space< double >,
                                                                                                 ang_velocity_3D_topology< double >,
                                                                                                 ang_accel_3D_topology< double > >,
                                                                               euclidean_tuple_distance > >,
                                      euclidean_tuple_distance >,
                  time_poisson_topology, reach_plus_time_metric >;


extern template metric_space_array< se3_0th_order_rl_topology< double >::type, 1 >::type
  joint_limits_mapping< double >::make_rl_joint_space(
    const metric_space_array< se3_0th_order_topology< double >::type, 1 >::type& ) const;
extern template metric_space_array< se3_1st_order_rl_topology< double >::type, 1 >::type
  joint_limits_mapping< double >::make_rl_joint_space(
    const metric_space_array< se3_1st_order_topology< double >::type, 1 >::type& ) const;
extern template metric_space_array< se3_2nd_order_rl_topology< double >::type, 1 >::type
  joint_limits_mapping< double >::make_rl_joint_space(
    const metric_space_array< se3_2nd_order_topology< double >::type, 1 >::type& ) const;

extern template metric_space_array< se3_0th_order_topology< double >::type, 1 >::type
  joint_limits_mapping< double >::make_normal_joint_space(
    const metric_space_array< se3_0th_order_rl_topology< double >::type, 1 >::type& ) const;
extern template metric_space_array< se3_1st_order_topology< double >::type, 1 >::type
  joint_limits_mapping< double >::make_normal_joint_space(
    const metric_space_array< se3_1st_order_rl_topology< double >::type, 1 >::type& ) const;
extern template metric_space_array< se3_2nd_order_topology< double >::type, 1 >::type
  joint_limits_mapping< double >::make_normal_joint_space(
    const metric_space_array< se3_2nd_order_rl_topology< double >::type, 1 >::type& ) const;

extern template topology_traits< metric_space_array< se3_0th_order_rl_topology< double >::type, 1 >::type >::point_type
  joint_limits_mapping< double >::map_to_space(
    const topology_traits< metric_space_array< se3_0th_order_topology< double >::type, 1 >::type >::point_type& pt,
    const metric_space_array< se3_0th_order_topology< double >::type, 1 >::type&,
    const metric_space_array< se3_0th_order_rl_topology< double >::type, 1 >::type& ) const;
extern template topology_traits< metric_space_array< se3_1st_order_rl_topology< double >::type, 1 >::type >::point_type
  joint_limits_mapping< double >::map_to_space(
    const topology_traits< metric_space_array< se3_1st_order_topology< double >::type, 1 >::type >::point_type& pt,
    const metric_space_array< se3_1st_order_topology< double >::type, 1 >::type&,
    const metric_space_array< se3_1st_order_rl_topology< double >::type, 1 >::type& ) const;
extern template topology_traits< metric_space_array< se3_2nd_order_rl_topology< double >::type, 1 >::type >::point_type
  joint_limits_mapping< double >::map_to_space(
    const topology_traits< metric_space_array< se3_2nd_order_topology< double >::type, 1 >::type >::point_type& pt,
    const metric_space_array< se3_2nd_order_topology< double >::type, 1 >::type&,
    const metric_space_array< se3_2nd_order_rl_topology< double >::type, 1 >::type& ) const;

extern template topology_traits< metric_space_array< se3_0th_order_topology< double >::type, 1 >::type >::point_type
  joint_limits_mapping< double >::map_to_space(
    const topology_traits< metric_space_array< se3_0th_order_rl_topology< double >::type, 1 >::type >::point_type& pt,
    const metric_space_array< se3_0th_order_rl_topology< double >::type, 1 >::type&,
    const metric_space_array< se3_0th_order_topology< double >::type, 1 >::type& ) const;
extern template topology_traits< metric_space_array< se3_1st_order_topology< double >::type, 1 >::type >::point_type
  joint_limits_mapping< double >::map_to_space(
    const topology_traits< metric_space_array< se3_1st_order_rl_topology< double >::type, 1 >::type >::point_type& pt,
    const metric_space_array< se3_1st_order_rl_topology< double >::type, 1 >::type&,
    const metric_space_array< se3_1st_order_topology< double >::type, 1 >::type& ) const;
extern template topology_traits< metric_space_array< se3_2nd_order_topology< double >::type, 1 >::type >::point_type
  joint_limits_mapping< double >::map_to_space(
    const topology_traits< metric_space_array< se3_2nd_order_rl_topology< double >::type, 1 >::type >::point_type& pt,
    const metric_space_array< se3_2nd_order_rl_topology< double >::type, 1 >::type&,
    const metric_space_array< se3_2nd_order_topology< double >::type, 1 >::type& ) const;
};


extern template frame_3D< double > get_frame_3D(
  const arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 >, vect< double, 3 > >,
                          arithmetic_tuple< unit_quat< double >, vect< double, 3 >, vect< double, 3 > > >& pt );

extern template frame_3D< double >
  get_frame_3D( const arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 > >,
                                        arithmetic_tuple< unit_quat< double >, vect< double, 3 > > >& pt );

extern template frame_3D< double > get_frame_3D(
  const arithmetic_tuple< arithmetic_tuple< vect< double, 3 > >, arithmetic_tuple< unit_quat< double > > >& pt );

extern template void
  set_frame_3D( arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 >, vect< double, 3 > >,
                                  arithmetic_tuple< unit_quat< double >, vect< double, 3 >, vect< double, 3 > > >& pt,
                const frame_3D< double >& p );

extern template void set_frame_3D( arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 > >,
                                                     arithmetic_tuple< unit_quat< double >, vect< double, 3 > > >& pt,
                                   const frame_3D< double >& p );

extern template void
  set_frame_3D( arithmetic_tuple< arithmetic_tuple< vect< double, 3 > >, arithmetic_tuple< unit_quat< double > > >& pt,
                const frame_3D< double >& p );

extern template pose_3D< double > get_pose_3D(
  const arithmetic_tuple< arithmetic_tuple< vect< double, 3 > >, arithmetic_tuple< unit_quat< double > > >& pt );

extern template void
  set_pose_3D( arithmetic_tuple< arithmetic_tuple< vect< double, 3 > >, arithmetic_tuple< unit_quat< double > > >& pt,
               const pose_3D< double >& p );

extern template const unit_quat< double >& get_quaternion(
  const arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 >, vect< double, 3 > >,
                          arithmetic_tuple< unit_quat< double >, vect< double, 3 >, vect< double, 3 > > >& pt );

extern template const unit_quat< double >&
  get_quaternion( const arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 > >,
                                          arithmetic_tuple< unit_quat< double >, vect< double, 3 > > >& pt );

extern template const unit_quat< double >& get_quaternion(
  const arithmetic_tuple< arithmetic_tuple< vect< double, 3 > >, arithmetic_tuple< unit_quat< double > > >& pt );

extern template void
  set_quaternion( arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 >, vect< double, 3 > >,
                                    arithmetic_tuple< unit_quat< double >, vect< double, 3 >, vect< double, 3 > > >& pt,
                  const unit_quat< double >& q );

extern template void set_quaternion( arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 > >,
                                                       arithmetic_tuple< unit_quat< double >, vect< double, 3 > > >& pt,
                                     const unit_quat< double >& q );

extern template void set_quaternion(
  arithmetic_tuple< arithmetic_tuple< vect< double, 3 > >, arithmetic_tuple< unit_quat< double > > >& pt,
  const unit_quat< double >& q );

extern template const vect< double, 3 >& get_position(
  const arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 >, vect< double, 3 > >,
                          arithmetic_tuple< unit_quat< double >, vect< double, 3 >, vect< double, 3 > > >& pt );

extern template const vect< double, 3 >&
  get_position( const arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 > >,
                                        arithmetic_tuple< unit_quat< double >, vect< double, 3 > > >& pt );

extern template const vect< double, 3 >& get_position(
  const arithmetic_tuple< arithmetic_tuple< vect< double, 3 > >, arithmetic_tuple< unit_quat< double > > >& pt );

extern template void
  set_position( arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 >, vect< double, 3 > >,
                                  arithmetic_tuple< unit_quat< double >, vect< double, 3 >, vect< double, 3 > > >& pt,
                const vect< double, 3 >& p );

extern template void set_position( arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 > >,
                                                     arithmetic_tuple< unit_quat< double >, vect< double, 3 > > >& pt,
                                   const vect< double, 3 >& p );

extern template void
  set_position( arithmetic_tuple< arithmetic_tuple< vect< double, 3 > >, arithmetic_tuple< unit_quat< double > > >& pt,
                const vect< double, 3 >& p );

extern template const vect< double, 3 >& get_ang_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 >, vect< double, 3 > >,
                          arithmetic_tuple< unit_quat< double >, vect< double, 3 >, vect< double, 3 > > >& pt );

extern template const vect< double, 3 >&
  get_ang_velocity( const arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 > >,
                                            arithmetic_tuple< unit_quat< double >, vect< double, 3 > > >& pt );

extern template void set_ang_velocity(
  arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 >, vect< double, 3 > >,
                    arithmetic_tuple< unit_quat< double >, vect< double, 3 >, vect< double, 3 > > >& pt,
  const vect< double, 3 >& p );

extern template void
  set_ang_velocity( arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 > >,
                                      arithmetic_tuple< unit_quat< double >, vect< double, 3 > > >& pt,
                    const vect< double, 3 >& p );

extern template const vect< double, 3 >& get_velocity(
  const arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 >, vect< double, 3 > >,
                          arithmetic_tuple< unit_quat< double >, vect< double, 3 >, vect< double, 3 > > >& pt );

extern template const vect< double, 3 >&
  get_velocity( const arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 > >,
                                        arithmetic_tuple< unit_quat< double >, vect< double, 3 > > >& pt );

extern template void
  set_velocity( arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 >, vect< double, 3 > >,
                                  arithmetic_tuple< unit_quat< double >, vect< double, 3 >, vect< double, 3 > > >& pt,
                const vect< double, 3 >& p );

extern template void set_velocity( arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 > >,
                                                     arithmetic_tuple< unit_quat< double >, vect< double, 3 > > >& pt,
                                   const vect< double, 3 >& p );


extern template const vect< double, 3 >& get_ang_acceleration(
  const arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 >, vect< double, 3 > >,
                          arithmetic_tuple< unit_quat< double >, vect< double, 3 >, vect< double, 3 > > >& pt );

extern template void set_ang_acceleration(
  arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 >, vect< double, 3 > >,
                    arithmetic_tuple< unit_quat< double >, vect< double, 3 >, vect< double, 3 > > >& pt,
  const vect< double, 3 >& p );


extern template const vect< double, 3 >& get_acceleration(
  const arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 >, vect< double, 3 > >,
                          arithmetic_tuple< unit_quat< double >, vect< double, 3 >, vect< double, 3 > > >& pt );

extern template void set_acceleration(
  arithmetic_tuple< arithmetic_tuple< vect< double, 3 >, vect< double, 3 >, vect< double, 3 > >,
                    arithmetic_tuple< unit_quat< double >, vect< double, 3 >, vect< double, 3 > > >& pt,
  const vect< double, 3 >& p );
};


#endif


#endif

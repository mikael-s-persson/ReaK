/**
 * \file joint_space_topologies.hpp
 * 
 * This library provides classes that define topologies in joint-space (generalized coordinates). 
 * All the topologies included here are for a single joint. To allow for more joints, use the 
 * meta-function metric_space_array to generate a metric_space_tuple (or form the metric_space_tuple 
 * manually).
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

#ifndef REAK_JOINT_SPACE_TOPOLOGIES_HPP
#define REAK_JOINT_SPACE_TOPOLOGIES_HPP


#include "base/defs.hpp"

#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "so3_topologies.hpp"

#include "differentiable_space.hpp"
#include "metric_space_tuple.hpp"

#include "hyperbox_topology.hpp"
#include "hyperball_topology.hpp"
#include "line_topology.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "lin_alg/vect_alg.hpp"

#include "kinetostatics/gen_coord.hpp"
#include "rate_limited_spaces.hpp"

namespace ReaK {

namespace pp {




/**
 * This meta-function defines the type for a 0th order single-joint space (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct joint_space_0th_order {
  typedef 
    differentiable_space< 
      time_topology, 
      arithmetic_tuple< line_segment_topology<T> >, 
      DistanceMetric 
    > type;
};

/**
 * This meta-function defines the type for a 1st order single-joint space.
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct joint_space_1st_order {
  typedef 
    differentiable_space< 
      time_topology, 
      arithmetic_tuple< 
        line_segment_topology<T>,
        line_segment_topology<T>
      >, 
      DistanceMetric 
    > type;
};

/**
 * This meta-function defines the type for a 2nd order single-joint space.
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct joint_space_2nd_order {
  typedef 
    differentiable_space< 
      time_topology, 
      arithmetic_tuple< 
        line_segment_topology<T>,
        line_segment_topology<T>,
        line_segment_topology<T>
      >, 
      DistanceMetric 
    > type;
};




/**
 * This meta-function defines the type for a rate-limited 0th order single-joint space (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct rl_joint_space_0th_order {
  typedef 
    reach_time_diff_space< 
      time_topology, 
      arithmetic_tuple< 
        line_segment_topology<T> 
      >, 
      DistanceMetric 
    > type;
};

/**
 * This meta-function defines the type for a rate-limited 1st order single-joint space.
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct rl_joint_space_1st_order {
  typedef 
    reach_time_diff_space< 
      time_topology, 
      arithmetic_tuple< 
        line_segment_topology<T>,
        line_segment_topology<T>
      >, 
      DistanceMetric 
    > type;
};

/**
 * This meta-function defines the type for a rate-limited 2nd order single-joint space.
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct rl_joint_space_2nd_order {
  typedef 
    reach_time_diff_space< 
      time_topology, 
      arithmetic_tuple< 
        line_segment_topology<T>,
        line_segment_topology<T>,
        line_segment_topology<T>
      >, 
      DistanceMetric 
    > type;
};




};



template <typename T>
gen_coord<T> get_gen_coord(
  const arithmetic_tuple< T, T, T >& pt) {
  return gen_coord<T>(get<0>(pt), 
		      get<1>(pt), 
		      get<2>(pt), 
		      0.0);
};

template <typename T>
gen_coord<T> get_gen_coord(
  const arithmetic_tuple< T, T >& pt) {
  return gen_coord<T>(get<0>(pt), 
		      get<1>(pt), 
		      0.0, 
		      0.0);
};

template <typename T>
gen_coord<T> get_gen_coord(
  const arithmetic_tuple< T >& pt) {
  return gen_coord<T>(get<0>(pt), 
		      0.0, 
		      0.0, 
		      0.0);
};

template <typename T>
void set_gen_coord(
  arithmetic_tuple< T, T, T >& pt,
  const gen_coord<T>& p) {
  get<0>(pt) = p.q;
  get<1>(pt) = p.q_dot;
  get<2>(pt) = p.q_ddot;
};

template <typename T>
void set_gen_coord(
  arithmetic_tuple< T, T >& pt,
  const gen_coord<T>& p) {
  get<0>(pt) = p.q;
  get<1>(pt) = p.q_dot;
};

template <typename T>
void set_gen_coord(
  arithmetic_tuple< T >& pt,
  const gen_coord<T>& p) {
  get<0>(pt) = p.q;
};



template <typename T>
const T& get_position(
  const arithmetic_tuple< T, T, T >& pt) {
  return get<0>(pt);
};

template <typename T>
const T& get_position(
  const arithmetic_tuple< T, T >& pt) {
  return get<0>(pt);
};

template <typename T>
const T& get_position(
  const arithmetic_tuple< T >& pt) {
  return get<0>(pt);
};

template <typename T>
void set_position(
  arithmetic_tuple< T, T, T >& pt,
  const T& p) {
  get<0>(pt) = p;
};

template <typename T>
void set_position(
  arithmetic_tuple< T, T >& pt,
  const T& p) {
  get<0>(pt) = p;
};

template <typename T>
void set_position(
  arithmetic_tuple< T >& pt,
  const T& p) {
  get<0>(pt) = p;
};



template <typename T>
const T& get_velocity(
  const arithmetic_tuple< T, T, T >& pt) {
  return get<1>(pt);
};

template <typename T>
const T& get_velocity(
  const arithmetic_tuple< T, T >& pt) {
  return get<1>(pt);
};

template <typename T>
void set_velocity(
  arithmetic_tuple< T, T, T >& pt,
  const T& p) {
  get<1>(pt) = p;
};

template <typename T>
void set_velocity(
  arithmetic_tuple< T, T >& pt,
  const T& p) {
  get<1>(pt) = p;
};


template <typename T>
const T& get_acceleration(
  const arithmetic_tuple< T, T, T >& pt) {
  return get<2>(pt);
};

template <typename T>
void set_acceleration(
  arithmetic_tuple< T, T, T >& pt,
  const T& p) {
  get<2>(pt) = p;
};







};

#endif









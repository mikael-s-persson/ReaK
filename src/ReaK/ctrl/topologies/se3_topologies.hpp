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


#include "base/defs.hpp"

#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "so3_topologies.hpp"

#include "differentiable_space.hpp"
#include "metric_space_tuple.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "lin_alg/vect_alg.hpp"

namespace ReaK {

namespace pp {




/**
 * This meta-function defines the type for a 0th order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se3_0th_order_topology {
  typedef 
    metric_space_tuple< arithmetic_tuple<
      differentiable_space< 
        time_topology, 
	arithmetic_tuple< vector_topology< vect<T,3> > >, 
	DistanceMetric 
      >,
      differentiable_space< 
        time_topology, 
	arithmetic_tuple< quaternion_topology<T> >, 
	DistanceMetric 
      > >,
      DistanceMetric 
    > type;
};

/**
 * This meta-function defines the type for a 1st order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se3_1st_order_topology {
  typedef 
    metric_space_tuple< arithmetic_tuple<
      differentiable_space< 
        time_topology, 
	arithmetic_tuple< 
	  vector_topology< vect<T,3> >,
	  vector_topology< vect<T,3> >
	>, 
	DistanceMetric 
      >,
      differentiable_space< 
        time_topology, 
	arithmetic_tuple< 
	  quaternion_topology<T>,
	  ang_velocity_3D_topology<T>
	>, 
	DistanceMetric 
      > >,
      DistanceMetric 
    > type;
};

/**
 * This meta-function defines the type for a 2nd order SE(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct se3_2nd_order_topology {
  typedef 
    metric_space_tuple< arithmetic_tuple<
      differentiable_space< 
        time_topology, 
	arithmetic_tuple< 
	  vector_topology< vect<T,3> >,
	  vector_topology< vect<T,3> >,
	  vector_topology< vect<T,3> >
	>, 
	DistanceMetric 
      >,
      differentiable_space< 
        time_topology, 
	arithmetic_tuple< 
	  quaternion_topology<T>,
	  ang_velocity_3D_topology<T>,
	  ang_accel_3D_topology<T>
	>, 
	DistanceMetric 
      > >,
      DistanceMetric 
    > type;
};






};

};

#endif









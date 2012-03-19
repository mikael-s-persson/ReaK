/**
 * \file joint_space_limits.hpp
 * 
 * This library provides classes to help create and manipulate joint-space topologies in over 
 * a joint-space with limits (speed, acceleration, and jerk limits).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2012
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

#ifndef REAK_JOINT_SPACE_LIMITS_HPP
#define REAK_JOINT_SPACE_LIMITS_HPP


#include "base/defs.hpp"

#include "base/named_object.hpp"

#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "joint_space_topologies.hpp"
#include "se2_topologies.hpp"
#include "se3_topologies.hpp"


#include <boost/mpl/and.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/prior.hpp>
#include <boost/mpl/less.hpp>


namespace ReaK {

namespace pp {

  
  
namespace detail {
  
  
  
  
/*******************************************************************************************************************
                                 FUNCTIONS TO CREATE RATE-LIMITED JOINT-SPACES
*******************************************************************************************************************/ 
  
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_joint_space< OutSpace >,
      boost::mpl::equal_to<
        max_derivation_order< InSpace, time_topology >,
	boost::mpl::size_t<0>
      >
    >,
  void >::type create_rl_joint_space_impl(OutSpace& space_out,
					  const InSpace& space_in,
					  const RateLimitMap& j_limits,
					  std::size_t& gen_i, std::size_t&, std::size_t&) {
    
    space_out = OutSpace(arithmetic_tuple<
                           line_segment_topology< typename RateLimitMap::value_type >
                         >(line_segment_topology< typename RateLimitMap::value_type >(
			     get<0>(space_in).getName() + "_rl",
			     ( get<0>(space_in).origin() - get<0>(space_in).get_radius() ) / j_limits.gen_speed_limits[gen_i],
		             ( get<0>(space_in).origin() + get<0>(space_in).get_radius() ) / j_limits.gen_speed_limits[gen_i]
			   )
			 )
		);
    ++gen_i;
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_joint_space< OutSpace >,
      boost::mpl::equal_to<
        max_derivation_order< InSpace, time_topology >,
	boost::mpl::size_t<1>
      >
    >,
  void >::type create_rl_joint_space_impl(OutSpace& space_out,
					  const InSpace& space_in,
					  const RateLimitMap& j_limits,
					  std::size_t& gen_i, std::size_t&, std::size_t&) {
    
    space_out = OutSpace(arithmetic_tuple<
                           line_segment_topology< typename RateLimitMap::value_type >,
                           line_segment_topology< typename RateLimitMap::value_type >
                         >(line_segment_topology< typename RateLimitMap::value_type >(
			     get<0>(space_in).getName() + "_rl",
			     ( get<0>(space_in).origin() - get<0>(space_in).get_radius() ) / j_limits.gen_speed_limits[gen_i],
		             ( get<0>(space_in).origin() + get<0>(space_in).get_radius() ) / j_limits.gen_speed_limits[gen_i]
			   ),
			   line_segment_topology< typename RateLimitMap::value_type >(
			     get<1>(space_in).getName() + "_rl",
			     ( get<1>(space_in).origin() - get<1>(space_in).get_radius() ) / j_limits.gen_accel_limits[gen_i],
		             ( get<1>(space_in).origin() + get<1>(space_in).get_radius() ) / j_limits.gen_accel_limits[gen_i]
			   )
			 ),
	                 euclidean_tuple_distance(),
	                 reach_time_differentiation_tuple< 1 >::type(
	                   reach_time_differentiation( j_limits.gen_speed_limits[gen_i] / j_limits.gen_accel_limits[gen_i] )
	                 )
		);
    ++gen_i;
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_joint_space< OutSpace >,
      boost::mpl::equal_to<
        max_derivation_order< InSpace, time_topology >,
	boost::mpl::size_t<2>
      >
    >,
  void >::type create_rl_joint_space_impl(OutSpace& space_out,
					  const InSpace& space_in,
					  const RateLimitMap& j_limits,
					  std::size_t& gen_i, std::size_t&, std::size_t&) {
    
    space_out = OutSpace(arithmetic_tuple<
                           line_segment_topology< typename RateLimitMap::value_type >,
                           line_segment_topology< typename RateLimitMap::value_type >,
                           line_segment_topology< typename RateLimitMap::value_type >
                         >(line_segment_topology< typename RateLimitMap::value_type >(
			     get<0>(space_in).getName() + "_rl",
			     ( get<0>(space_in).origin() - get<0>(space_in).get_radius() ) / j_limits.gen_speed_limits[gen_i],
		             ( get<0>(space_in).origin() + get<0>(space_in).get_radius() ) / j_limits.gen_speed_limits[gen_i]
			   ),
			   line_segment_topology< typename RateLimitMap::value_type >(
			     get<1>(space_in).getName() + "_rl",
			     ( get<1>(space_in).origin() - get<1>(space_in).get_radius() ) / j_limits.gen_accel_limits[gen_i],
		             ( get<1>(space_in).origin() + get<1>(space_in).get_radius() ) / j_limits.gen_accel_limits[gen_i]
			   ),
			   line_segment_topology< typename RateLimitMap::value_type >(
			     get<2>(space_in).getName() + "_rl",
			     ( get<2>(space_in).origin() - get<2>(space_in).get_radius() ) / j_limits.gen_jerk_limits[gen_i],
		             ( get<2>(space_in).origin() + get<2>(space_in).get_radius() ) / j_limits.gen_jerk_limits[gen_i]
			   )
			 ),
	                 euclidean_tuple_distance(),
	                 reach_time_differentiation_tuple< 2 >::type(
	                   reach_time_differentiation( j_limits.gen_speed_limits[gen_i] / j_limits.gen_accel_limits[gen_i] ),
	                   reach_time_differentiation( j_limits.gen_accel_limits[gen_i] / j_limits.gen_jerk_limits[gen_i] )
	                 )
		);
    ++gen_i;
  };
  
  
  
  
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_se2_space< OutSpace >,
      boost::mpl::equal_to<
        max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type, time_topology >,
	boost::mpl::size_t<0>
      >
    >,
  void >::type create_rl_joint_space_impl(OutSpace& space_out,
					  const InSpace& space_in,
					  const RateLimitMap& j_limits,
					  std::size_t&, std::size_t& f2d_i, std::size_t&) {
    
    typedef typename RateLimitMap::value_type ValueType;
    typedef line_segment_topology< ValueType > LineSegTopo;
    typedef hyperbox_topology< vect<ValueType, 2> > BoxTopo;
    typedef arithmetic_tuple< LineSegTopo > RotTopoTuple;
    typedef arithmetic_tuple< BoxTopo > PosTopoTuple;
    typedef typename arithmetic_tuple_element<0, OutSpace>::type PosTopoOutType;
    typedef typename arithmetic_tuple_element<1, OutSpace>::type RotTopoOutType;
    
    space_out = OutSpace(arithmetic_tuple< PosTopoTuple, RotTopoTuple >(
                           PosTopoOutType(
			     PosTopoTuple(
			       BoxTopo(
			         get<0>(get<0>(space_in)).getName() + "_rl",
			         get<0>(get<0>(space_in)).get_lower_corner() * ( ValueType(1.0) / j_limits.frame2D_speed_limits[f2d_i]),
			         get<0>(get<0>(space_in)).get_upper_corner() * ( ValueType(1.0) / j_limits.frame2D_speed_limits[f2d_i])
			       )
			     )
			   ),
			   RotTopoOutType(
			     RotTopoTuple(
			       LineSegTopo(
			         get<0>(get<1>(space_in)).getName() + "_rl",
			         ( get<0>(get<1>(space_in)).origin() - get<0>(get<1>(space_in)).get_radius() ) / j_limits.frame2D_speed_limits[f2d_i + 1],
		                 ( get<0>(get<1>(space_in)).origin() + get<0>(get<1>(space_in)).get_radius() ) / j_limits.frame2D_speed_limits[f2d_i + 1]
			       )
			     )
			   )
			 )
		);
    f2d_i += 2;
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_se2_space< OutSpace >,
      boost::mpl::equal_to<
        max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type, time_topology >,
	boost::mpl::size_t<1>
      >
    >,
  void >::type create_rl_joint_space_impl(OutSpace& space_out,
					  const InSpace& space_in,
					  const RateLimitMap& j_limits,
					  std::size_t&, std::size_t& f2d_i, std::size_t&) {
    
    typedef typename RateLimitMap::value_type ValueType;
    typedef line_segment_topology< ValueType > LineSegTopo;
    typedef hyperbox_topology< vect<ValueType, 2> > BoxTopo;
    typedef hyperball_topology< vect<ValueType, 2> > BallTopo;
    typedef arithmetic_tuple< LineSegTopo, LineSegTopo > RotTopoTuple;
    typedef arithmetic_tuple< BoxTopo, BallTopo > PosTopoTuple;
    typedef typename arithmetic_tuple_element<0, OutSpace>::type PosTopoOutType;
    typedef typename arithmetic_tuple_element<1, OutSpace>::type RotTopoOutType;
    
    space_out = OutSpace(arithmetic_tuple< PosTopoTuple, RotTopoTuple >(
                           PosTopoOutType(
			     PosTopoTuple(
			       BoxTopo(
			         get<0>(get<0>(space_in)).getName() + "_rl",
			         get<0>(get<0>(space_in)).get_lower_corner() * ( ValueType(1.0) / j_limits.frame2D_speed_limits[f2d_i]),
			         get<0>(get<0>(space_in)).get_upper_corner() * ( ValueType(1.0) / j_limits.frame2D_speed_limits[f2d_i])
			       ),
			       BallTopo(
			         get<1>(get<0>(space_in)).getName() + "_rl",
			         get<1>(get<0>(space_in)).origin() * ( ValueType(1.0) / j_limits.frame2D_accel_limits[f2d_i]),
			         get<1>(get<0>(space_in)).get_radius() / j_limits.frame2D_accel_limits[f2d_i]
			       )
			     ),
		             euclidean_tuple_distance(),
			     reach_time_differentiation_tuple< 1 >::type(
	                       reach_time_differentiation( j_limits.frame2D_speed_limits[f2d_i] / j_limits.frame2D_accel_limits[f2d_i] )
	                     )
			   ),
			   RotTopoOutType(
			     RotTopoTuple(
			       LineSegTopo(
			         get<0>(get<1>(space_in)).getName() + "_rl",
			         ( get<0>(get<1>(space_in)).origin() - get<0>(get<1>(space_in)).get_radius() ) / j_limits.frame2D_speed_limits[f2d_i + 1],
		                 ( get<0>(get<1>(space_in)).origin() + get<0>(get<1>(space_in)).get_radius() ) / j_limits.frame2D_speed_limits[f2d_i + 1]
			       ),
			       LineSegTopo(
			         get<1>(get<1>(space_in)).getName() + "_rl",
			         ( get<1>(get<1>(space_in)).origin() - get<1>(get<1>(space_in)).get_radius() ) / j_limits.frame2D_accel_limits[f2d_i + 1],
		                 ( get<1>(get<1>(space_in)).origin() + get<1>(get<1>(space_in)).get_radius() ) / j_limits.frame2D_accel_limits[f2d_i + 1]
			       )
			     ),
		             euclidean_tuple_distance(),
			     reach_time_differentiation_tuple< 1 >::type(
	                       reach_time_differentiation( j_limits.frame2D_speed_limits[f2d_i + 1] / j_limits.frame2D_accel_limits[f2d_i + 1] )
	                     )
			   )
			 )
		);
    f2d_i += 2;
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_se2_space< OutSpace >,
      boost::mpl::equal_to<
        max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type, time_topology >,
	boost::mpl::size_t<2>
      >
    >,
  void >::type create_rl_joint_space_impl(OutSpace& space_out,
					  const InSpace& space_in,
					  const RateLimitMap& j_limits,
					  std::size_t&, std::size_t& f2d_i, std::size_t&) {
    
    typedef typename RateLimitMap::value_type ValueType;
    typedef line_segment_topology< ValueType > LineSegTopo;
    typedef hyperbox_topology< vect<ValueType, 2> > BoxTopo;
    typedef hyperball_topology< vect<ValueType, 2> > BallTopo;
    typedef arithmetic_tuple< LineSegTopo, LineSegTopo, LineSegTopo > RotTopoTuple;
    typedef arithmetic_tuple< BoxTopo, BallTopo, BallTopo > PosTopoTuple;
    typedef typename arithmetic_tuple_element<0, OutSpace>::type PosTopoOutType;
    typedef typename arithmetic_tuple_element<1, OutSpace>::type RotTopoOutType;
    
    space_out = OutSpace(arithmetic_tuple< PosTopoTuple, RotTopoTuple >(
                           PosTopoOutType(
			     PosTopoTuple(
			       BoxTopo(
			         get<0>(get<0>(space_in)).getName() + "_rl",
			         get<0>(get<0>(space_in)).get_lower_corner() * ( ValueType(1.0) / j_limits.frame2D_speed_limits[f2d_i]),
			         get<0>(get<0>(space_in)).get_upper_corner() * ( ValueType(1.0) / j_limits.frame2D_speed_limits[f2d_i])
			       ),
			       BallTopo(
			         get<1>(get<0>(space_in)).getName() + "_rl",
			         get<1>(get<0>(space_in)).origin() * ( ValueType(1.0) / j_limits.frame2D_accel_limits[f2d_i]),
			         get<1>(get<0>(space_in)).get_radius() / j_limits.frame2D_accel_limits[f2d_i]
			       ),
			       BallTopo(
			         get<2>(get<0>(space_in)).getName() + "_rl",
			         get<2>(get<0>(space_in)).origin() * ( ValueType(1.0) / j_limits.frame2D_jerk_limits[f2d_i]),
			         get<2>(get<0>(space_in)).get_radius() / j_limits.frame2D_jerk_limits[f2d_i]
			       )
			     ),
		             euclidean_tuple_distance(),
			     reach_time_differentiation_tuple< 2 >::type(
	                       reach_time_differentiation( j_limits.frame2D_speed_limits[f2d_i] / j_limits.frame2D_accel_limits[f2d_i] ),
	                       reach_time_differentiation( j_limits.frame2D_accel_limits[f2d_i] / j_limits.frame2D_jerk_limits[f2d_i] )
	                     )
			   ),
			   RotTopoOutType(
			     RotTopoTuple(
			       LineSegTopo(
			         get<0>(get<1>(space_in)).getName() + "_rl",
			         ( get<0>(get<1>(space_in)).origin() - get<0>(get<1>(space_in)).get_radius() ) / j_limits.frame2D_speed_limits[f2d_i + 1],
		                 ( get<0>(get<1>(space_in)).origin() + get<0>(get<1>(space_in)).get_radius() ) / j_limits.frame2D_speed_limits[f2d_i + 1]
			       ),
			       LineSegTopo(
			         get<1>(get<1>(space_in)).getName() + "_rl",
			         ( get<1>(get<1>(space_in)).origin() - get<1>(get<1>(space_in)).get_radius() ) / j_limits.frame2D_accel_limits[f2d_i + 1],
		                 ( get<1>(get<1>(space_in)).origin() + get<1>(get<1>(space_in)).get_radius() ) / j_limits.frame2D_accel_limits[f2d_i + 1]
			       ),
			       LineSegTopo(
			         get<2>(get<1>(space_in)).getName() + "_rl",
			         ( get<2>(get<1>(space_in)).origin() - get<2>(get<1>(space_in)).get_radius() ) / j_limits.frame2D_jerk_limits[f2d_i + 1],
		                 ( get<2>(get<1>(space_in)).origin() + get<2>(get<1>(space_in)).get_radius() ) / j_limits.frame2D_jerk_limits[f2d_i + 1]
			       )
			     ),
		             euclidean_tuple_distance(),
			     reach_time_differentiation_tuple< 2 >::type(
	                       reach_time_differentiation( j_limits.frame2D_speed_limits[f2d_i + 1] / j_limits.frame2D_accel_limits[f2d_i + 1] ),
	                       reach_time_differentiation( j_limits.frame2D_accel_limits[f2d_i + 1] / j_limits.frame2D_jerk_limits[f2d_i + 1] )
	                     )
			   )
			 )
		);
    f2d_i += 2;
  };
  
  
  
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_se3_space< OutSpace >,
      boost::mpl::equal_to<
        max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type, time_topology >,
	boost::mpl::size_t<0>
      >
    >,
  void >::type create_rl_joint_space_impl(OutSpace& space_out,
					  const InSpace& space_in,
					  const RateLimitMap& j_limits,
					  std::size_t&, std::size_t&, std::size_t& f3d_i) {
    
    typedef typename RateLimitMap::value_type ValueType;
    typedef rate_limited_quat_space< ValueType > QuatTopo;
    typedef hyperbox_topology< vect<ValueType, 3> > BoxTopo;
    typedef arithmetic_tuple< QuatTopo > RotTopoTuple;
    typedef arithmetic_tuple< BoxTopo > PosTopoTuple;
    typedef typename arithmetic_tuple_element<0, OutSpace>::type PosTopoOutType;
    typedef typename arithmetic_tuple_element<1, OutSpace>::type RotTopoOutType;
    
    space_out = OutSpace(arithmetic_tuple< PosTopoTuple, RotTopoTuple >(
                           PosTopoOutType(
			     PosTopoTuple(
			       BoxTopo(
			         get<0>(get<0>(space_in)).getName() + "_rl",
			         get<0>(get<0>(space_in)).get_lower_corner() * ( ValueType(1.0) / j_limits.frame3D_speed_limits[f3d_i]),
			         get<0>(get<0>(space_in)).get_upper_corner() * ( ValueType(1.0) / j_limits.frame3D_speed_limits[f3d_i])
			       )
			     )
			   ),
			   RotTopoOutType(
			     RotTopoTuple(
			       QuatTopo(
			         get<0>(get<1>(space_in)).getName() + "_rl",
			         j_limits.frame3D_speed_limits[f3d_i + 1]
			       )
			     )
			   )
			 )
		);
    f3d_i += 2;
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_se3_space< OutSpace >,
      boost::mpl::equal_to<
        max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type, time_topology >,
	boost::mpl::size_t<1>
      >
    >,
  void >::type create_rl_joint_space_impl(OutSpace& space_out,
					  const InSpace& space_in,
					  const RateLimitMap& j_limits,
					  std::size_t&, std::size_t&, std::size_t& f3d_i) {
    
    typedef typename RateLimitMap::value_type ValueType;
    typedef rate_limited_quat_space< ValueType > QuatTopo;
    typedef ang_velocity_3D_topology< ValueType > AngVelTopo;
    typedef hyperbox_topology< vect<ValueType, 3> > BoxTopo;
    typedef hyperball_topology< vect<ValueType, 3> > BallTopo;
    typedef arithmetic_tuple< QuatTopo, AngVelTopo > RotTopoTuple;
    typedef arithmetic_tuple< BoxTopo, BallTopo > PosTopoTuple;
    typedef typename arithmetic_tuple_element<0, OutSpace>::type PosTopoOutType;
    typedef typename arithmetic_tuple_element<1, OutSpace>::type RotTopoOutType;
    
    space_out = OutSpace(arithmetic_tuple< PosTopoTuple, RotTopoTuple >(
                           PosTopoOutType(
			     PosTopoTuple(
			       BoxTopo(
			         get<0>(get<0>(space_in)).getName() + "_rl",
			         get<0>(get<0>(space_in)).get_lower_corner() * ( ValueType(1.0) / j_limits.frame3D_speed_limits[f3d_i]),
			         get<0>(get<0>(space_in)).get_upper_corner() * ( ValueType(1.0) / j_limits.frame3D_speed_limits[f3d_i])
			       ),
			       BallTopo(
			         get<1>(get<0>(space_in)).getName() + "_rl",
			         get<1>(get<0>(space_in)).origin() * ( ValueType(1.0) / j_limits.frame3D_accel_limits[f3d_i]),
			         get<1>(get<0>(space_in)).get_radius() / j_limits.frame3D_accel_limits[f3d_i]
			       )
			     ),
		             euclidean_tuple_distance(),
			     reach_time_differentiation_tuple< 1 >::type(
	                       reach_time_differentiation( j_limits.frame3D_speed_limits[f3d_i] / j_limits.frame3D_accel_limits[f3d_i] )
	                     )
			   ),
			   RotTopoOutType(
			     RotTopoTuple(
			       QuatTopo(
			         get<0>(get<1>(space_in)).getName() + "_rl",
			         j_limits.frame3D_speed_limits[f3d_i + 1]
			       ),
			       AngVelTopo(
			         get<1>(get<1>(space_in)).getName() + "_rl",
			         get<1>(get<1>(space_in)).get_radius() / j_limits.frame3D_accel_limits[f3d_i + 1]
		               )
			     ),
		             euclidean_tuple_distance(),
			     reach_time_differentiation_tuple< 1 >::type(
	                       reach_time_differentiation( j_limits.frame3D_speed_limits[f3d_i + 1] / j_limits.frame3D_accel_limits[f3d_i + 1] )
	                     )
			   )
			 )
		);
    f3d_i += 2;
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_se3_space< OutSpace >,
      boost::mpl::equal_to<
        max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type, time_topology >,
	boost::mpl::size_t<2>
      >
    >,
  void >::type create_rl_joint_space_impl(OutSpace& space_out,
					  const InSpace& space_in,
					  const RateLimitMap& j_limits,
					  std::size_t&, std::size_t&, std::size_t& f3d_i) {
    
    typedef typename RateLimitMap::value_type ValueType;
    typedef rate_limited_quat_space< ValueType > QuatTopo;
    typedef ang_velocity_3D_topology< ValueType > AngVelTopo;
    typedef ang_accel_3D_topology< ValueType > AngAccTopo;
    typedef hyperbox_topology< vect<ValueType, 3> > BoxTopo;
    typedef hyperball_topology< vect<ValueType, 3> > BallTopo;
    typedef arithmetic_tuple< QuatTopo, AngVelTopo, AngAccTopo > RotTopoTuple;
    typedef arithmetic_tuple< BoxTopo, BallTopo, BallTopo > PosTopoTuple;
    typedef typename arithmetic_tuple_element<0, OutSpace>::type PosTopoOutType;
    typedef typename arithmetic_tuple_element<1, OutSpace>::type RotTopoOutType;
    
    space_out = OutSpace(arithmetic_tuple< PosTopoTuple, RotTopoTuple >(
                           PosTopoOutType(
			     PosTopoTuple(
			       BoxTopo(
			         get<0>(get<0>(space_in)).getName() + "_rl",
			         get<0>(get<0>(space_in)).get_lower_corner() * ( ValueType(1.0) / j_limits.frame3D_speed_limits[f3d_i]),
			         get<0>(get<0>(space_in)).get_upper_corner() * ( ValueType(1.0) / j_limits.frame3D_speed_limits[f3d_i])
			       ),
			       BallTopo(
			         get<1>(get<0>(space_in)).getName() + "_rl",
			         get<1>(get<0>(space_in)).origin() * ( ValueType(1.0) / j_limits.frame3D_accel_limits[f3d_i]),
			         get<1>(get<0>(space_in)).get_radius() / j_limits.frame3D_accel_limits[f3d_i]
			       ),
			       BallTopo(
			         get<2>(get<0>(space_in)).getName() + "_rl",
			         get<2>(get<0>(space_in)).origin() * ( ValueType(1.0) / j_limits.frame3D_jerk_limits[f3d_i]),
			         get<2>(get<0>(space_in)).get_radius() / j_limits.frame3D_jerk_limits[f3d_i]
			       )
			     ),
		             euclidean_tuple_distance(),
			     reach_time_differentiation_tuple< 2 >::type(
	                       reach_time_differentiation( j_limits.frame3D_speed_limits[f3d_i] / j_limits.frame3D_accel_limits[f3d_i] ),
	                       reach_time_differentiation( j_limits.frame3D_accel_limits[f3d_i] / j_limits.frame3D_jerk_limits[f3d_i] )
	                     )
			   ),
			   RotTopoOutType(
			     RotTopoTuple(
			       QuatTopo(
			         get<0>(get<1>(space_in)).getName() + "_rl",
			         j_limits.frame3D_speed_limits[f3d_i + 1]
			       ),
			       AngVelTopo(
			         get<1>(get<1>(space_in)).getName() + "_rl",
			         get<1>(get<1>(space_in)).get_radius() / j_limits.frame3D_accel_limits[f3d_i + 1]
		               ),
			       AngAccTopo(
			         get<2>(get<1>(space_in)).getName() + "_rl",
			         get<2>(get<1>(space_in)).get_radius() / j_limits.frame3D_jerk_limits[f3d_i + 1]
		               )
			     ),
		             euclidean_tuple_distance(),
			     reach_time_differentiation_tuple< 2 >::type(
	                       reach_time_differentiation( j_limits.frame3D_speed_limits[f3d_i + 1] / j_limits.frame3D_accel_limits[f3d_i + 1] ),
	                       reach_time_differentiation( j_limits.frame3D_accel_limits[f3d_i + 1] / j_limits.frame3D_jerk_limits[f3d_i + 1] )
	                     )
			   )
			 )
		);
    f3d_i += 2;
  };
  
  // declarations only:
  template <typename Idx, typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_rl_joint_spaces_impl(OutSpace& space_out,
					   const InSpace& space_in,
					   const RateLimitMap& j_limits,
					   std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i);
  
  // declarations only:
  template <typename Idx, typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_rl_joint_spaces_impl(OutSpace& space_out,
					   const InSpace& space_in,
					   const RateLimitMap& j_limits,
					   std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i);
  
  
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_rate_limited_joint_space< OutSpace >,
      is_rate_limited_se2_space< OutSpace >,
      is_rate_limited_se3_space< OutSpace >
    >,
  void >::type create_rl_joint_space_impl(OutSpace& space_out,
					  const InSpace& space_in,
					  const RateLimitMap& j_limits,
					  std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i) {
    create_rl_joint_spaces_impl< typename boost::mpl::prior< arithmetic_tuple_size< OutSpace > >::type >(space_out, space_in, j_limits, gen_i, f2d_i, f3d_i);
  };
  
  
  
  template <typename Idx, typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_rl_joint_spaces_impl(OutSpace& space_out,
					   const InSpace& space_in,
					   const RateLimitMap& j_limits,
					   std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i) {
    
    create_rl_joint_spaces_impl< typename boost::mpl::prior<Idx>::type >(space_out,space_in,j_limits,gen_i,f2d_i,f3d_i);
    
    create_rl_joint_space_impl(get< Idx::type::value >(space_out), 
                               get< Idx::type::value >(space_in),
			       j_limits, gen_i, f2d_i, f3d_i);
    
  };
  
  template <typename Idx, typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_rl_joint_spaces_impl(OutSpace& space_out,
					   const InSpace& space_in,
					   const RateLimitMap& j_limits,
					   std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i) {
    
    create_rl_joint_space_impl(get< 0 >(space_out), 
                               get< 0 >(space_in),
			       j_limits, gen_i, f2d_i, f3d_i);
    
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_rate_limited_joint_space< OutSpace >,
      is_rate_limited_se2_space< OutSpace >,
      is_rate_limited_se3_space< OutSpace >
    >,
  void >::type create_rl_joint_spaces_impl(OutSpace& space_out,
				   const InSpace& space_in,
				   const RateLimitMap& j_limits) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    create_rl_joint_spaces_impl< typename boost::mpl::prior< arithmetic_tuple_size< OutSpace > >::type >(space_out, 
											  space_in,
											  j_limits, gen_i, f2d_i, f3d_i);
    
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::or_<
      is_rate_limited_joint_space< OutSpace >,
      is_rate_limited_se2_space< OutSpace >,
      is_rate_limited_se3_space< OutSpace >
    >,
  void >::type create_rl_joint_spaces_impl(OutSpace& space_out,
				   const InSpace& space_in,
				   const RateLimitMap& j_limits) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    create_rl_joint_space_impl(space_out, space_in,
		               j_limits, gen_i, f2d_i, f3d_i);
    
  };
  
  
  
  
  
  
  
/*******************************************************************************************************************
                                 FUNCTIONS TO CREATE NORMAL JOINT-SPACES
*******************************************************************************************************************/ 
  
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_joint_space< InSpace >,
      boost::mpl::equal_to<
        max_derivation_order< InSpace, time_topology >,
	boost::mpl::size_t<0>
      >
    >,
  void >::type create_normal_joint_space_impl(OutSpace& space_out,
					      const InSpace& space_in,
					      const RateLimitMap& j_limits,
					      std::size_t& gen_i, std::size_t&, std::size_t&) {
    
    space_out = OutSpace(arithmetic_tuple<
                           line_segment_topology< typename RateLimitMap::value_type >
                         >(line_segment_topology< typename RateLimitMap::value_type >(
			     get<0>(space_in).getName() + "_non_rl",
			     ( get<0>(space_in).origin() - get<0>(space_in).get_radius() ) * j_limits.gen_speed_limits[gen_i],
		             ( get<0>(space_in).origin() + get<0>(space_in).get_radius() ) * j_limits.gen_speed_limits[gen_i]
			   )
			 )
		);
    ++gen_i;
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_joint_space< InSpace >,
      boost::mpl::equal_to<
        max_derivation_order< InSpace, time_topology >,
	boost::mpl::size_t<1>
      >
    >,
  void >::type create_normal_joint_space_impl(OutSpace& space_out,
					      const InSpace& space_in,
					      const RateLimitMap& j_limits,
					      std::size_t& gen_i, std::size_t&, std::size_t&) {
    
    space_out = OutSpace(arithmetic_tuple<
                           line_segment_topology< typename RateLimitMap::value_type >,
                           line_segment_topology< typename RateLimitMap::value_type >
                         >(line_segment_topology< typename RateLimitMap::value_type >(
			     get<0>(space_in).getName() + "_non_rl",
			     ( get<0>(space_in).origin() - get<0>(space_in).get_radius() ) * j_limits.gen_speed_limits[gen_i],
		             ( get<0>(space_in).origin() + get<0>(space_in).get_radius() ) * j_limits.gen_speed_limits[gen_i]
			   ),
			   line_segment_topology< typename RateLimitMap::value_type >(
			     get<1>(space_in).getName() + "_non_rl",
			     ( get<1>(space_in).origin() - get<1>(space_in).get_radius() ) * j_limits.gen_accel_limits[gen_i],
		             ( get<1>(space_in).origin() + get<1>(space_in).get_radius() ) * j_limits.gen_accel_limits[gen_i]
			   )
			 )
		);
    ++gen_i;
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_joint_space< InSpace >,
      boost::mpl::equal_to<
        max_derivation_order< InSpace, time_topology >,
	boost::mpl::size_t<2>
      >
    >,
  void >::type create_normal_joint_space_impl(OutSpace& space_out,
					      const InSpace& space_in,
					      const RateLimitMap& j_limits,
					      std::size_t& gen_i, std::size_t&, std::size_t&) {
    
    space_out = OutSpace(arithmetic_tuple<
                           line_segment_topology< typename RateLimitMap::value_type >,
                           line_segment_topology< typename RateLimitMap::value_type >,
                           line_segment_topology< typename RateLimitMap::value_type >
                         >(line_segment_topology< typename RateLimitMap::value_type >(
			     get<0>(space_in).getName() + "_non_rl",
			     ( get<0>(space_in).origin() - get<0>(space_in).get_radius() ) * j_limits.gen_speed_limits[gen_i],
		             ( get<0>(space_in).origin() + get<0>(space_in).get_radius() ) * j_limits.gen_speed_limits[gen_i]
			   ),
			   line_segment_topology< typename RateLimitMap::value_type >(
			     get<1>(space_in).getName() + "_non_rl",
			     ( get<1>(space_in).origin() - get<1>(space_in).get_radius() ) * j_limits.gen_accel_limits[gen_i],
		             ( get<1>(space_in).origin() + get<1>(space_in).get_radius() ) * j_limits.gen_accel_limits[gen_i]
			   ),
			   line_segment_topology< typename RateLimitMap::value_type >(
			     get<2>(space_in).getName() + "_non_rl",
			     ( get<2>(space_in).origin() - get<2>(space_in).get_radius() ) * j_limits.gen_jerk_limits[gen_i],
		             ( get<2>(space_in).origin() + get<2>(space_in).get_radius() ) * j_limits.gen_jerk_limits[gen_i]
			   )
			 )
		);
    ++gen_i;
  };
  
  
  
  
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_se2_space< InSpace >,
      boost::mpl::equal_to<
        max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type, time_topology >,
	boost::mpl::size_t<0>
      >
    >,
  void >::type create_normal_joint_space_impl(OutSpace& space_out,
					      const InSpace& space_in,
					      const RateLimitMap& j_limits,
					      std::size_t&, std::size_t& f2d_i, std::size_t&) {
    
    typedef typename RateLimitMap::value_type ValueType;
    typedef line_segment_topology< ValueType > LineSegTopo;
    typedef hyperbox_topology< vect<ValueType, 2> > BoxTopo;
    typedef arithmetic_tuple< LineSegTopo > RotTopoTuple;
    typedef arithmetic_tuple< BoxTopo > PosTopoTuple;
    typedef typename arithmetic_tuple_element<0, OutSpace>::type PosTopoOutType;
    typedef typename arithmetic_tuple_element<1, OutSpace>::type RotTopoOutType;
    
    space_out = OutSpace(arithmetic_tuple< PosTopoTuple, RotTopoTuple >(
                           PosTopoOutType(
			     PosTopoTuple(
			       BoxTopo(
			         get<0>(get<0>(space_in)).getName() + "_non_rl",
			         get<0>(get<0>(space_in)).get_lower_corner() * j_limits.frame2D_speed_limits[f2d_i],
			         get<0>(get<0>(space_in)).get_upper_corner() * j_limits.frame2D_speed_limits[f2d_i]
			       )
			     )
			   ),
			   RotTopoOutType(
			     RotTopoTuple(
			       LineSegTopo(
			         get<0>(get<1>(space_in)).getName() + "_non_rl",
			         ( get<0>(get<1>(space_in)).origin() - get<0>(get<1>(space_in)).get_radius() ) * j_limits.frame2D_speed_limits[f2d_i + 1],
		                 ( get<0>(get<1>(space_in)).origin() + get<0>(get<1>(space_in)).get_radius() ) * j_limits.frame2D_speed_limits[f2d_i + 1]
			       )
			     )
			   )
			 )
		);
    f2d_i += 2;
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_se2_space< InSpace >,
      boost::mpl::equal_to<
        max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type, time_topology >,
	boost::mpl::size_t<1>
      >
    >,
  void >::type create_normal_joint_space_impl(OutSpace& space_out,
					      const InSpace& space_in,
					      const RateLimitMap& j_limits,
					      std::size_t&, std::size_t& f2d_i, std::size_t&) {
    
    typedef typename RateLimitMap::value_type ValueType;
    typedef line_segment_topology< ValueType > LineSegTopo;
    typedef hyperbox_topology< vect<ValueType, 2> > BoxTopo;
    typedef hyperball_topology< vect<ValueType, 2> > BallTopo;
    typedef arithmetic_tuple< LineSegTopo, LineSegTopo > RotTopoTuple;
    typedef arithmetic_tuple< BoxTopo, BallTopo > PosTopoTuple;
    typedef typename arithmetic_tuple_element<0, OutSpace>::type PosTopoOutType;
    typedef typename arithmetic_tuple_element<1, OutSpace>::type RotTopoOutType;
    
    space_out = OutSpace(arithmetic_tuple< PosTopoTuple, RotTopoTuple >(
                           PosTopoOutType(
			     PosTopoTuple(
			       BoxTopo(
			         get<0>(get<0>(space_in)).getName() + "_non_rl",
			         get<0>(get<0>(space_in)).get_lower_corner() * j_limits.frame2D_speed_limits[f2d_i],
			         get<0>(get<0>(space_in)).get_upper_corner() * j_limits.frame2D_speed_limits[f2d_i]
			       ),
			       BallTopo(
			         get<1>(get<0>(space_in)).getName() + "_non_rl",
			         get<1>(get<0>(space_in)).origin() * j_limits.frame2D_accel_limits[f2d_i],
			         get<1>(get<0>(space_in)).get_radius() * j_limits.frame2D_accel_limits[f2d_i]
			       )
			     )
			   ),
			   RotTopoOutType(
			     RotTopoTuple(
			       LineSegTopo(
			         get<0>(get<1>(space_in)).getName() + "_non_rl",
			         ( get<0>(get<1>(space_in)).origin() - get<0>(get<1>(space_in)).get_radius() ) * j_limits.frame2D_speed_limits[f2d_i + 1],
		                 ( get<0>(get<1>(space_in)).origin() + get<0>(get<1>(space_in)).get_radius() ) * j_limits.frame2D_speed_limits[f2d_i + 1]
			       ),
			       LineSegTopo(
			         get<1>(get<1>(space_in)).getName() + "_non_rl",
			         ( get<1>(get<1>(space_in)).origin() - get<1>(get<1>(space_in)).get_radius() ) * j_limits.frame2D_accel_limits[f2d_i + 1],
		                 ( get<1>(get<1>(space_in)).origin() + get<1>(get<1>(space_in)).get_radius() ) * j_limits.frame2D_accel_limits[f2d_i + 1]
			       )
			     )
			   )
			 )
		);
    f2d_i += 2;
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_se2_space< InSpace >,
      boost::mpl::equal_to<
        max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type, time_topology >,
	boost::mpl::size_t<2>
      >
    >,
  void >::type create_normal_joint_space_impl(OutSpace& space_out,
					      const InSpace& space_in,
					      const RateLimitMap& j_limits,
					      std::size_t&, std::size_t& f2d_i, std::size_t&) {
    
    typedef typename RateLimitMap::value_type ValueType;
    typedef line_segment_topology< ValueType > LineSegTopo;
    typedef hyperbox_topology< vect<ValueType, 2> > BoxTopo;
    typedef hyperball_topology< vect<ValueType, 2> > BallTopo;
    typedef arithmetic_tuple< LineSegTopo, LineSegTopo, LineSegTopo > RotTopoTuple;
    typedef arithmetic_tuple< BoxTopo, BallTopo, BallTopo > PosTopoTuple;
    typedef typename arithmetic_tuple_element<0, OutSpace>::type PosTopoOutType;
    typedef typename arithmetic_tuple_element<1, OutSpace>::type RotTopoOutType;
    
    space_out = OutSpace(arithmetic_tuple< PosTopoTuple, RotTopoTuple >(
                           PosTopoOutType(
			     PosTopoTuple(
			       BoxTopo(
			         get<0>(get<0>(space_in)).getName() + "_non_rl",
			         get<0>(get<0>(space_in)).get_lower_corner() * j_limits.frame2D_speed_limits[f2d_i],
			         get<0>(get<0>(space_in)).get_upper_corner() * j_limits.frame2D_speed_limits[f2d_i]
			       ),
			       BallTopo(
			         get<1>(get<0>(space_in)).getName() + "_non_rl",
			         get<1>(get<0>(space_in)).origin() * j_limits.frame2D_accel_limits[f2d_i],
			         get<1>(get<0>(space_in)).get_radius() * j_limits.frame2D_accel_limits[f2d_i]
			       ),
			       BallTopo(
			         get<2>(get<0>(space_in)).getName() + "_non_rl",
			         get<2>(get<0>(space_in)).origin() * j_limits.frame2D_jerk_limits[f2d_i],
			         get<2>(get<0>(space_in)).get_radius() * j_limits.frame2D_jerk_limits[f2d_i]
			       )
			     )
			   ),
			   RotTopoOutType(
			     RotTopoTuple(
			       LineSegTopo(
			         get<0>(get<1>(space_in)).getName() + "_non_rl",
			         ( get<0>(get<1>(space_in)).origin() - get<0>(get<1>(space_in)).get_radius() ) * j_limits.frame2D_speed_limits[f2d_i + 1],
		                 ( get<0>(get<1>(space_in)).origin() + get<0>(get<1>(space_in)).get_radius() ) * j_limits.frame2D_speed_limits[f2d_i + 1]
			       ),
			       LineSegTopo(
			         get<1>(get<1>(space_in)).getName() + "_non_rl",
			         ( get<1>(get<1>(space_in)).origin() - get<1>(get<1>(space_in)).get_radius() ) * j_limits.frame2D_accel_limits[f2d_i + 1],
		                 ( get<1>(get<1>(space_in)).origin() + get<1>(get<1>(space_in)).get_radius() ) * j_limits.frame2D_accel_limits[f2d_i + 1]
			       ),
			       LineSegTopo(
			         get<2>(get<1>(space_in)).getName() + "_non_rl",
			         ( get<2>(get<1>(space_in)).origin() - get<2>(get<1>(space_in)).get_radius() ) * j_limits.frame2D_jerk_limits[f2d_i + 1],
		                 ( get<2>(get<1>(space_in)).origin() + get<2>(get<1>(space_in)).get_radius() ) * j_limits.frame2D_jerk_limits[f2d_i + 1]
			       )
			     )
			   )
			 )
		);
    f2d_i += 2;
  };
  
  
  
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_se3_space< InSpace >,
      boost::mpl::equal_to<
        max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type, time_topology >,
	boost::mpl::size_t<0>
      >
    >,
  void >::type create_normal_joint_space_impl(OutSpace& space_out,
					      const InSpace& space_in,
					      const RateLimitMap& j_limits,
					      std::size_t&, std::size_t&, std::size_t& f3d_i) {
    
    typedef typename RateLimitMap::value_type ValueType;
    typedef quaternion_topology< ValueType > QuatTopo;
    typedef hyperbox_topology< vect<ValueType, 3> > BoxTopo;
    typedef arithmetic_tuple< QuatTopo > RotTopoTuple;
    typedef arithmetic_tuple< BoxTopo > PosTopoTuple;
    typedef typename arithmetic_tuple_element<0, OutSpace>::type PosTopoOutType;
    typedef typename arithmetic_tuple_element<1, OutSpace>::type RotTopoOutType;
    
    space_out = OutSpace(arithmetic_tuple< PosTopoTuple, RotTopoTuple >(
                           PosTopoOutType(
			     PosTopoTuple(
			       BoxTopo(
			         get<0>(get<0>(space_in)).getName() + "_non_rl",
			         get<0>(get<0>(space_in)).get_lower_corner() * j_limits.frame3D_speed_limits[f3d_i],
			         get<0>(get<0>(space_in)).get_upper_corner() * j_limits.frame3D_speed_limits[f3d_i]
			       )
			     )
			   ),
			   RotTopoOutType(
			     RotTopoTuple(
			       QuatTopo(
			         get<0>(get<1>(space_in)).getName() + "_non_rl"
			       )
			     )
			   )
			 )
		);
    f3d_i += 2;
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_se3_space< InSpace >,
      boost::mpl::equal_to<
        max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type, time_topology >,
	boost::mpl::size_t<1>
      >
    >,
  void >::type create_normal_joint_space_impl(OutSpace& space_out,
					      const InSpace& space_in,
					      const RateLimitMap& j_limits,
					      std::size_t&, std::size_t&, std::size_t& f3d_i) {
    
    typedef typename RateLimitMap::value_type ValueType;
    typedef quaternion_topology< ValueType > QuatTopo;
    typedef ang_velocity_3D_topology< ValueType > AngVelTopo;
    typedef hyperbox_topology< vect<ValueType, 3> > BoxTopo;
    typedef hyperball_topology< vect<ValueType, 3> > BallTopo;
    typedef arithmetic_tuple< QuatTopo, AngVelTopo > RotTopoTuple;
    typedef arithmetic_tuple< BoxTopo, BallTopo > PosTopoTuple;
    typedef typename arithmetic_tuple_element<0, OutSpace>::type PosTopoOutType;
    typedef typename arithmetic_tuple_element<1, OutSpace>::type RotTopoOutType;
    
    space_out = OutSpace(arithmetic_tuple< PosTopoTuple, RotTopoTuple >(
                           PosTopoOutType(
			     PosTopoTuple(
			       BoxTopo(
			         get<0>(get<0>(space_in)).getName() + "_non_rl",
			         get<0>(get<0>(space_in)).get_lower_corner() * j_limits.frame3D_speed_limits[f3d_i],
			         get<0>(get<0>(space_in)).get_upper_corner() * j_limits.frame3D_speed_limits[f3d_i]
			       ),
			       BallTopo(
			         get<1>(get<0>(space_in)).getName() + "_non_rl",
			         get<1>(get<0>(space_in)).origin() * j_limits.frame3D_accel_limits[f3d_i],
			         get<1>(get<0>(space_in)).get_radius() * j_limits.frame3D_accel_limits[f3d_i]
			       )
			     )
			   ),
			   RotTopoOutType(
			     RotTopoTuple(
			       QuatTopo(
			         get<0>(get<1>(space_in)).getName() + "_non_rl"
			       ),
			       AngVelTopo(
			         get<1>(get<1>(space_in)).getName() + "_non_rl",
			         get<1>(get<1>(space_in)).get_radius() * j_limits.frame3D_accel_limits[f3d_i + 1]
		               )
			     )
			   )
			 )
		);
    f3d_i += 2;
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if<
    boost::mpl::and_<
      is_rate_limited_se3_space< InSpace >,
      boost::mpl::equal_to<
        max_derivation_order< typename arithmetic_tuple_element< 0, InSpace >::type, time_topology >,
	boost::mpl::size_t<2>
      >
    >,
  void >::type create_normal_joint_space_impl(OutSpace& space_out,
					      const InSpace& space_in,
					      const RateLimitMap& j_limits,
					      std::size_t&, std::size_t&, std::size_t& f3d_i) {
    
    typedef typename RateLimitMap::value_type ValueType;
    typedef quaternion_topology< ValueType > QuatTopo;
    typedef ang_velocity_3D_topology< ValueType > AngVelTopo;
    typedef ang_accel_3D_topology< ValueType > AngAccTopo;
    typedef hyperbox_topology< vect<ValueType, 3> > BoxTopo;
    typedef hyperball_topology< vect<ValueType, 3> > BallTopo;
    typedef arithmetic_tuple< QuatTopo, AngVelTopo, AngAccTopo > RotTopoTuple;
    typedef arithmetic_tuple< BoxTopo, BallTopo, BallTopo > PosTopoTuple;
    typedef typename arithmetic_tuple_element<0, OutSpace>::type PosTopoOutType;
    typedef typename arithmetic_tuple_element<1, OutSpace>::type RotTopoOutType;
    
    space_out = OutSpace(arithmetic_tuple< PosTopoTuple, RotTopoTuple >(
                           PosTopoOutType(
			     PosTopoTuple(
			       BoxTopo(
			         get<0>(get<0>(space_in)).getName() + "_non_rl",
			         get<0>(get<0>(space_in)).get_lower_corner() * j_limits.frame3D_speed_limits[f3d_i],
			         get<0>(get<0>(space_in)).get_upper_corner() * j_limits.frame3D_speed_limits[f3d_i]
			       ),
			       BallTopo(
			         get<1>(get<0>(space_in)).getName() + "_non_rl",
			         get<1>(get<0>(space_in)).origin() * j_limits.frame3D_accel_limits[f3d_i],
			         get<1>(get<0>(space_in)).get_radius() * j_limits.frame3D_accel_limits[f3d_i]
			       ),
			       BallTopo(
			         get<2>(get<0>(space_in)).getName() + "_non_rl",
			         get<2>(get<0>(space_in)).origin() * j_limits.frame3D_jerk_limits[f3d_i],
			         get<2>(get<0>(space_in)).get_radius() * j_limits.frame3D_jerk_limits[f3d_i]
			       )
			     )
			   ),
			   RotTopoOutType(
			     RotTopoTuple(
			       QuatTopo(
			         get<0>(get<1>(space_in)).getName() + "_non_rl"
			       ),
			       AngVelTopo(
			         get<1>(get<1>(space_in)).getName() + "_non_rl",
			         get<1>(get<1>(space_in)).get_radius() * j_limits.frame3D_accel_limits[f3d_i + 1]
		               ),
			       AngAccTopo(
			         get<2>(get<1>(space_in)).getName() + "_non_rl",
			         get<2>(get<1>(space_in)).get_radius() * j_limits.frame3D_jerk_limits[f3d_i + 1]
		               )
			     )
			   )
			 )
		);
    f3d_i += 2;
  };
  
  
  // declaration only:
  template <typename Idx, typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_normal_joint_spaces_impl(OutSpace& space_out,
					       const InSpace& space_in,
					       const RateLimitMap& j_limits,
					       std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i);
  
  // declaration only:
  template <typename Idx, typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_normal_joint_spaces_impl(OutSpace& space_out,
					       const InSpace& space_in,
					       const RateLimitMap& j_limits,
					       std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i);
  
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_rate_limited_joint_space< InSpace >,
      is_rate_limited_se2_space< InSpace >,
      is_rate_limited_se3_space< InSpace >
    >,
  void >::type create_normal_joint_space_impl(OutSpace& space_out,
					      const InSpace& space_in,
					      const RateLimitMap& j_limits,
					      std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i) {
    create_normal_joint_spaces_impl< typename boost::mpl::prior< arithmetic_tuple_size< OutSpace > >::type >(space_out, space_in, j_limits, gen_i, f2d_i, f3d_i);
  };
  
  
  
  template <typename Idx, typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_normal_joint_spaces_impl(OutSpace& space_out,
					       const InSpace& space_in,
					       const RateLimitMap& j_limits,
					       std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i) {
    
    create_normal_joint_spaces_impl< typename boost::mpl::prior<Idx>::type >(space_out,space_in,j_limits,gen_i,f2d_i,f3d_i);
    
    create_normal_joint_space_impl(get< Idx::type::value >(space_out), 
                                   get< Idx::type::value >(space_in),
			           j_limits, gen_i, f2d_i, f3d_i);
    
  };
  
  template <typename Idx, typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_normal_joint_spaces_impl(OutSpace& space_out,
					       const InSpace& space_in,
					       const RateLimitMap& j_limits,
					       std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i) {
    
    create_normal_joint_space_impl(get< 0 >(space_out), 
                                   get< 0 >(space_in),
			           j_limits, gen_i, f2d_i, f3d_i);
    
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_rate_limited_joint_space< InSpace >,
      is_rate_limited_se2_space< InSpace >,
      is_rate_limited_se3_space< InSpace >
    >,
  void >::type create_normal_joint_spaces_impl(OutSpace& space_out,
				       const InSpace& space_in,
				       const RateLimitMap& j_limits) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    create_normal_joint_spaces_impl< typename boost::mpl::prior< arithmetic_tuple_size< OutSpace > >::type >(space_out, 
											      space_in,
											      j_limits, gen_i, f2d_i, f3d_i);
    
  };
  
  template <typename OutSpace, typename InSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::or_<
      is_rate_limited_joint_space< InSpace >,
      is_rate_limited_se2_space< InSpace >,
      is_rate_limited_se3_space< InSpace >
    >,
  void >::type create_normal_joint_spaces_impl(OutSpace& space_out,
				       const InSpace& space_in,
				       const RateLimitMap& j_limits) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    create_normal_joint_space_impl(space_out, space_in, j_limits, gen_i, f2d_i, f3d_i);
  };
  
  
  
  
/*******************************************************************************************************************
                                 FUNCTIONS TO CREATE RATE-LIMITED JOINT-SPACE VECTORS
*******************************************************************************************************************/ 
  

  template <typename RateLimitMap>
  void create_rl_joint_vector_impl(arithmetic_tuple< typename RateLimitMap::value_type >& result,
				   const arithmetic_tuple< typename RateLimitMap::value_type >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t& gen_i, std::size_t&, std::size_t&) {
    get< 0 >(result) = get< 0 >(pt) / j_limits.gen_speed_limits[gen_i];
    ++gen_i;
  };

  template <typename RateLimitMap>
  void create_rl_joint_vector_impl(arithmetic_tuple< 
                                     typename RateLimitMap::value_type,
                                     typename RateLimitMap::value_type
                                   >& result,
				   const arithmetic_tuple< 
				     typename RateLimitMap::value_type,
				     typename RateLimitMap::value_type
				   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t& gen_i, std::size_t&, std::size_t&) {
    get< 0 >(result) = get< 0 >(pt) / j_limits.gen_speed_limits[gen_i];
    get< 1 >(result) = get< 1 >(pt) / j_limits.gen_accel_limits[gen_i];
    ++gen_i;
  };

  template <typename RateLimitMap>
  void create_rl_joint_vector_impl(arithmetic_tuple< 
                                     typename RateLimitMap::value_type,
                                     typename RateLimitMap::value_type,
                                     typename RateLimitMap::value_type
                                   >& result,
				   const arithmetic_tuple< 
				     typename RateLimitMap::value_type,
				     typename RateLimitMap::value_type,
				     typename RateLimitMap::value_type
				   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t& gen_i, std::size_t&, std::size_t&) {
    get< 0 >(result) = get< 0 >(pt) / j_limits.gen_speed_limits[gen_i];
    get< 1 >(result) = get< 1 >(pt) / j_limits.gen_accel_limits[gen_i];
    get< 2 >(result) = get< 2 >(pt) / j_limits.gen_jerk_limits[gen_i];
    ++gen_i;
  };
  
  
  template <typename RateLimitMap>
  void create_rl_joint_vector_impl(arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,2> 
                                     >,
                                     arithmetic_tuple< 
                                       typename RateLimitMap::value_type 
                                     > 
                                   >& result,
				   const arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,2> 
                                     >,
                                     arithmetic_tuple< 
                                       typename RateLimitMap::value_type 
                                     > 
                                   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t&, std::size_t& f2d_i, std::size_t&) {
    get< 0 >(get< 0 >(result)) = get< 0 >(get< 0 >(pt)) * (1.0 / j_limits.frame2D_speed_limits[f2d_i]);
    get< 0 >(get< 1 >(result)) = get< 0 >(get< 1 >(pt)) / j_limits.frame2D_speed_limits[f2d_i + 1];
    f2d_i += 2;
  };

  template <typename RateLimitMap>
  void create_rl_joint_vector_impl(arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,2>,
                                       vect<typename RateLimitMap::value_type,2> 
                                     >,
                                     arithmetic_tuple< 
                                       typename RateLimitMap::value_type,
                                       typename RateLimitMap::value_type 
                                     > 
                                   >& result,
				   const arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,2>,
                                       vect<typename RateLimitMap::value_type,2> 
                                     >,
                                     arithmetic_tuple< 
                                       typename RateLimitMap::value_type,
                                       typename RateLimitMap::value_type 
                                     > 
                                   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t&, std::size_t& f2d_i, std::size_t&) {
    get< 0 >(get< 0 >(result)) = get< 0 >(get< 0 >(pt)) * (1.0 / j_limits.frame2D_speed_limits[f2d_i]);
    get< 0 >(get< 1 >(result)) = get< 0 >(get< 1 >(pt)) / j_limits.frame2D_speed_limits[f2d_i + 1];
    get< 1 >(get< 0 >(result)) = get< 1 >(get< 0 >(pt)) * (1.0 / j_limits.frame2D_accel_limits[f2d_i]);
    get< 1 >(get< 1 >(result)) = get< 1 >(get< 1 >(pt)) / j_limits.frame2D_accel_limits[f2d_i + 1];
    f2d_i += 2;
  };

  template <typename RateLimitMap>
  void create_rl_joint_vector_impl(arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,2>,
                                       vect<typename RateLimitMap::value_type,2>,
                                       vect<typename RateLimitMap::value_type,2> 
                                     >,
                                     arithmetic_tuple< 
                                       typename RateLimitMap::value_type,
                                       typename RateLimitMap::value_type,
                                       typename RateLimitMap::value_type
                                     > 
                                   >& result,
				   const arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,2>,
                                       vect<typename RateLimitMap::value_type,2>,
                                       vect<typename RateLimitMap::value_type,2> 
                                     >,
                                     arithmetic_tuple< 
                                       typename RateLimitMap::value_type,
                                       typename RateLimitMap::value_type,
                                       typename RateLimitMap::value_type
                                     > 
                                   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t&, std::size_t& f2d_i, std::size_t&) {
    get< 0 >(get< 0 >(result)) = get< 0 >(get< 0 >(pt)) * (1.0 / j_limits.frame2D_speed_limits[f2d_i]);
    get< 0 >(get< 1 >(result)) = get< 0 >(get< 1 >(pt)) / j_limits.frame2D_speed_limits[f2d_i + 1];
    get< 1 >(get< 0 >(result)) = get< 1 >(get< 0 >(pt)) * (1.0 / j_limits.frame2D_accel_limits[f2d_i]);
    get< 1 >(get< 1 >(result)) = get< 1 >(get< 1 >(pt)) / j_limits.frame2D_accel_limits[f2d_i + 1];
    get< 2 >(get< 0 >(result)) = get< 2 >(get< 0 >(pt)) * (1.0 / j_limits.frame2D_jerk_limits[f2d_i]);
    get< 2 >(get< 1 >(result)) = get< 2 >(get< 1 >(pt)) / j_limits.frame2D_jerk_limits[f2d_i + 1];
    f2d_i += 2;
  };
  
  
  template <typename RateLimitMap>
  void create_rl_joint_vector_impl(arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,3> 
                                     >,
                                     arithmetic_tuple< 
                                       unit_quat< typename RateLimitMap::value_type >
                                     > 
                                   >& result,
				   const arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,3> 
                                     >,
                                     arithmetic_tuple< 
                                       unit_quat< typename RateLimitMap::value_type >
                                     > 
                                   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t&, std::size_t&, std::size_t& f3d_i) {
    get< 0 >(get< 0 >(result)) = get< 0 >(get< 0 >(pt)) * (1.0 / j_limits.frame3D_speed_limits[f3d_i]);
    get< 0 >(get< 1 >(result)) = get< 0 >(get< 1 >(pt));
    f3d_i += 2;
  };

  template <typename RateLimitMap>
  void create_rl_joint_vector_impl(arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3> 
                                     >,
                                     arithmetic_tuple< 
                                       unit_quat< typename RateLimitMap::value_type >,
                                       vect<typename RateLimitMap::value_type,3>
                                     > 
                                   >& result,
				   const arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3> 
                                     >,
                                     arithmetic_tuple< 
                                       unit_quat< typename RateLimitMap::value_type >,
                                       vect<typename RateLimitMap::value_type,3> 
                                     > 
                                   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t&, std::size_t&, std::size_t& f3d_i) {
    get< 0 >(get< 0 >(result)) = get< 0 >(get< 0 >(pt)) * (1.0 / j_limits.frame3D_speed_limits[f3d_i]);
    get< 0 >(get< 1 >(result)) = get< 0 >(get< 1 >(pt));
    get< 1 >(get< 0 >(result)) = get< 1 >(get< 0 >(pt)) * (1.0 / j_limits.frame3D_accel_limits[f3d_i]);
    get< 1 >(get< 1 >(result)) = get< 1 >(get< 1 >(pt)) * (1.0 / j_limits.frame3D_accel_limits[f3d_i + 1]);
    f3d_i += 2;
  };

  template <typename RateLimitMap>
  void create_rl_joint_vector_impl(arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3> 
                                     >,
                                     arithmetic_tuple< 
                                       unit_quat< typename RateLimitMap::value_type >,
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3>
                                     > 
                                   >& result,
				   const arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3> 
                                     >,
                                     arithmetic_tuple< 
                                       unit_quat< typename RateLimitMap::value_type >,
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3>
                                     > 
                                   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t&, std::size_t&, std::size_t& f3d_i) {
    get< 0 >(get< 0 >(result)) = get< 0 >(get< 0 >(pt)) * (1.0 / j_limits.frame3D_speed_limits[f3d_i]);
    get< 0 >(get< 1 >(result)) = get< 0 >(get< 1 >(pt));
    get< 1 >(get< 0 >(result)) = get< 1 >(get< 0 >(pt)) * (1.0 / j_limits.frame3D_accel_limits[f3d_i]);
    get< 1 >(get< 1 >(result)) = get< 1 >(get< 1 >(pt)) * (1.0 / j_limits.frame3D_accel_limits[f3d_i + 1]);
    get< 2 >(get< 0 >(result)) = get< 2 >(get< 0 >(pt)) * (1.0 / j_limits.frame3D_jerk_limits[f3d_i]);
    get< 2 >(get< 1 >(result)) = get< 2 >(get< 1 >(pt)) * (1.0 / j_limits.frame3D_jerk_limits[f3d_i + 1]);
    f3d_i += 2;
  };
  
  
  // declaration only:
  template <typename Idx, typename OutPoint, typename InPoint, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_rl_joint_vectors_impl(OutPoint& result,
				            const InPoint& pt,
					    const RateLimitMap& j_limits,
					    std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i);
  
  // declaration only:
  template <typename Idx, typename OutPoint, typename InPoint, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_rl_joint_vectors_impl(OutPoint& result,
				            const InPoint& pt,
					    const RateLimitMap& j_limits,
					    std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i);
  

  template <typename OutPoint, typename InPoint, typename RateLimitMap>
  void create_rl_joint_vector_impl(OutPoint& result,
				   const InPoint& pt,
				   const RateLimitMap& j_limits,
				   std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i) {
    create_rl_joint_vectors_impl< typename boost::mpl::prior< arithmetic_tuple_size< OutPoint > >::type >(result, pt, j_limits, gen_i, f2d_i, f3d_i);
  };
  
  
  
  template <typename Idx, typename OutPoint, typename InPoint, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_rl_joint_vectors_impl(OutPoint& result,
				            const InPoint& pt,
					    const RateLimitMap& j_limits,
					    std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i) {
    
    create_rl_joint_vectors_impl< typename boost::mpl::prior<Idx>::type >(result,pt,j_limits,gen_i,f2d_i,f3d_i);
    
    create_rl_joint_vector_impl(get< Idx::type::value >(result), 
                                get< Idx::type::value >(pt),
			        j_limits, gen_i, f2d_i, f3d_i);
    
  };
  
  template <typename Idx, typename OutPoint, typename InPoint, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_rl_joint_vectors_impl(OutPoint& result,
				            const InPoint& pt,
					    const RateLimitMap& j_limits,
					    std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i) {
    
    create_rl_joint_vector_impl(get< 0 >(result), 
                                get< 0 >(pt),
			        j_limits, gen_i, f2d_i, f3d_i);
    
  };
  
  template <typename OutPoint, typename InPoint, typename RateLimitMap>
  void create_rl_joint_vectors_impl(OutPoint& result,
				    const InPoint& pt,
				    const RateLimitMap& j_limits) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    create_rl_joint_vectors_impl< typename boost::mpl::prior< arithmetic_tuple_size< OutPoint > >::type >(result, pt, j_limits, gen_i, f2d_i, f3d_i);
  };
  
  
  
/*******************************************************************************************************************
                                 FUNCTIONS TO CREATE NORMAL JOINT-SPACE VECTORS
*******************************************************************************************************************/ 
  
  

  template <typename RateLimitMap>
  void create_normal_joint_vector_impl(arithmetic_tuple< typename RateLimitMap::value_type >& result,
				   const arithmetic_tuple< typename RateLimitMap::value_type >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t& gen_i, std::size_t&, std::size_t&) {
    get< 0 >(result) = get< 0 >(pt) * j_limits.gen_speed_limits[gen_i];
    ++gen_i;
  };

  template <typename RateLimitMap>
  void create_normal_joint_vector_impl(arithmetic_tuple< 
                                     typename RateLimitMap::value_type,
                                     typename RateLimitMap::value_type
                                   >& result,
				   const arithmetic_tuple< 
				     typename RateLimitMap::value_type,
				     typename RateLimitMap::value_type
				   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t& gen_i, std::size_t&, std::size_t&) {
    get< 0 >(result) = get< 0 >(pt) * j_limits.gen_speed_limits[gen_i];
    get< 1 >(result) = get< 1 >(pt) * j_limits.gen_accel_limits[gen_i];
    ++gen_i;
  };

  template <typename RateLimitMap>
  void create_normal_joint_vector_impl(arithmetic_tuple< 
                                     typename RateLimitMap::value_type,
                                     typename RateLimitMap::value_type,
                                     typename RateLimitMap::value_type
                                   >& result,
				   const arithmetic_tuple< 
				     typename RateLimitMap::value_type,
				     typename RateLimitMap::value_type,
				     typename RateLimitMap::value_type
				   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t& gen_i, std::size_t&, std::size_t&) {
    get< 0 >(result) = get< 0 >(pt) * j_limits.gen_speed_limits[gen_i];
    get< 1 >(result) = get< 1 >(pt) * j_limits.gen_accel_limits[gen_i];
    get< 2 >(result) = get< 2 >(pt) * j_limits.gen_jerk_limits[gen_i];
    ++gen_i;
  };
  
  
  template <typename RateLimitMap>
  void create_normal_joint_vector_impl(arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,2> 
                                     >,
                                     arithmetic_tuple< 
                                       typename RateLimitMap::value_type 
                                     > 
                                   >& result,
				   const arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,2> 
                                     >,
                                     arithmetic_tuple< 
                                       typename RateLimitMap::value_type 
                                     > 
                                   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t&, std::size_t& f2d_i, std::size_t&) {
    get< 0 >(get< 0 >(result)) = get< 0 >(get< 0 >(pt)) * j_limits.frame2D_speed_limits[f2d_i];
    get< 0 >(get< 1 >(result)) = get< 0 >(get< 1 >(pt)) * j_limits.frame2D_speed_limits[f2d_i + 1];
    f2d_i += 2;
  };

  template <typename RateLimitMap>
  void create_normal_joint_vector_impl(arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,2>,
                                       vect<typename RateLimitMap::value_type,2> 
                                     >,
                                     arithmetic_tuple< 
                                       typename RateLimitMap::value_type,
                                       typename RateLimitMap::value_type 
                                     > 
                                   >& result,
				   const arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,2>,
                                       vect<typename RateLimitMap::value_type,2> 
                                     >,
                                     arithmetic_tuple< 
                                       typename RateLimitMap::value_type,
                                       typename RateLimitMap::value_type 
                                     > 
                                   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t&, std::size_t& f2d_i, std::size_t&) {
    get< 0 >(get< 0 >(result)) = get< 0 >(get< 0 >(pt)) * j_limits.frame2D_speed_limits[f2d_i];
    get< 0 >(get< 1 >(result)) = get< 0 >(get< 1 >(pt)) * j_limits.frame2D_speed_limits[f2d_i + 1];
    get< 1 >(get< 0 >(result)) = get< 1 >(get< 0 >(pt)) * j_limits.frame2D_accel_limits[f2d_i];
    get< 1 >(get< 1 >(result)) = get< 1 >(get< 1 >(pt)) * j_limits.frame2D_accel_limits[f2d_i + 1];
    f2d_i += 2;
  };

  template <typename RateLimitMap>
  void create_normal_joint_vector_impl(arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,2>,
                                       vect<typename RateLimitMap::value_type,2>,
                                       vect<typename RateLimitMap::value_type,2> 
                                     >,
                                     arithmetic_tuple< 
                                       typename RateLimitMap::value_type,
                                       typename RateLimitMap::value_type,
                                       typename RateLimitMap::value_type
                                     > 
                                   >& result,
				   const arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,2>,
                                       vect<typename RateLimitMap::value_type,2>,
                                       vect<typename RateLimitMap::value_type,2> 
                                     >,
                                     arithmetic_tuple< 
                                       typename RateLimitMap::value_type,
                                       typename RateLimitMap::value_type,
                                       typename RateLimitMap::value_type
                                     > 
                                   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t&, std::size_t& f2d_i, std::size_t&) {
    get< 0 >(get< 0 >(result)) = get< 0 >(get< 0 >(pt)) * j_limits.frame2D_speed_limits[f2d_i];
    get< 0 >(get< 1 >(result)) = get< 0 >(get< 1 >(pt)) * j_limits.frame2D_speed_limits[f2d_i + 1];
    get< 1 >(get< 0 >(result)) = get< 1 >(get< 0 >(pt)) * j_limits.frame2D_accel_limits[f2d_i];
    get< 1 >(get< 1 >(result)) = get< 1 >(get< 1 >(pt)) * j_limits.frame2D_accel_limits[f2d_i + 1];
    get< 2 >(get< 0 >(result)) = get< 2 >(get< 0 >(pt)) * j_limits.frame2D_jerk_limits[f2d_i];
    get< 2 >(get< 1 >(result)) = get< 2 >(get< 1 >(pt)) * j_limits.frame2D_jerk_limits[f2d_i + 1];
    f2d_i += 2;
  };
  
  
  template <typename RateLimitMap>
  void create_normal_joint_vector_impl(arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,3> 
                                     >,
                                     arithmetic_tuple< 
                                       unit_quat< typename RateLimitMap::value_type >
                                     > 
                                   >& result,
				   const arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,3> 
                                     >,
                                     arithmetic_tuple< 
                                       unit_quat< typename RateLimitMap::value_type >
                                     > 
                                   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t&, std::size_t&, std::size_t& f3d_i) {
    get< 0 >(get< 0 >(result)) = get< 0 >(get< 0 >(pt)) * j_limits.frame3D_speed_limits[f3d_i];
    get< 0 >(get< 1 >(result)) = get< 0 >(get< 1 >(pt));
    f3d_i += 2;
  };

  template <typename RateLimitMap>
  void create_normal_joint_vector_impl(arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3> 
                                     >,
                                     arithmetic_tuple< 
                                       unit_quat< typename RateLimitMap::value_type >,
                                       vect<typename RateLimitMap::value_type,3>
                                     > 
                                   >& result,
				   const arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3> 
                                     >,
                                     arithmetic_tuple< 
                                       unit_quat< typename RateLimitMap::value_type >,
                                       vect<typename RateLimitMap::value_type,3> 
                                     > 
                                   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t&, std::size_t&, std::size_t& f3d_i) {
    get< 0 >(get< 0 >(result)) = get< 0 >(get< 0 >(pt)) * j_limits.frame3D_speed_limits[f3d_i];
    get< 0 >(get< 1 >(result)) = get< 0 >(get< 1 >(pt));
    get< 1 >(get< 0 >(result)) = get< 1 >(get< 0 >(pt)) * j_limits.frame3D_accel_limits[f3d_i];
    get< 1 >(get< 1 >(result)) = get< 1 >(get< 1 >(pt)) * j_limits.frame3D_accel_limits[f3d_i + 1];
    f3d_i += 2;
  };

  template <typename RateLimitMap>
  void create_normal_joint_vector_impl(arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3> 
                                     >,
                                     arithmetic_tuple< 
                                       unit_quat< typename RateLimitMap::value_type >,
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3>
                                     > 
                                   >& result,
				   const arithmetic_tuple< 
                                     arithmetic_tuple< 
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3> 
                                     >,
                                     arithmetic_tuple< 
                                       unit_quat< typename RateLimitMap::value_type >,
                                       vect<typename RateLimitMap::value_type,3>,
                                       vect<typename RateLimitMap::value_type,3>
                                     > 
                                   >& pt,
				   const RateLimitMap& j_limits,
				   std::size_t&, std::size_t&, std::size_t& f3d_i) {
    get< 0 >(get< 0 >(result)) = get< 0 >(get< 0 >(pt)) * j_limits.frame3D_speed_limits[f3d_i];
    get< 0 >(get< 1 >(result)) = get< 0 >(get< 1 >(pt));
    get< 1 >(get< 0 >(result)) = get< 1 >(get< 0 >(pt)) * j_limits.frame3D_accel_limits[f3d_i];
    get< 1 >(get< 1 >(result)) = get< 1 >(get< 1 >(pt)) * j_limits.frame3D_accel_limits[f3d_i + 1];
    get< 2 >(get< 0 >(result)) = get< 2 >(get< 0 >(pt)) * j_limits.frame3D_jerk_limits[f3d_i];
    get< 2 >(get< 1 >(result)) = get< 2 >(get< 1 >(pt)) * j_limits.frame3D_jerk_limits[f3d_i + 1];
    f3d_i += 2;
  };
  
  
  
  // declaration only:
  template <typename Idx, typename OutPoint, typename InPoint, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_normal_joint_vectors_impl(OutPoint& result,
				            const InPoint& pt,
					    const RateLimitMap& j_limits,
					    std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i);
  
  // declaration only:
  template <typename Idx, typename OutPoint, typename InPoint, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_normal_joint_vectors_impl(OutPoint& result,
				            const InPoint& pt,
					    const RateLimitMap& j_limits,
					    std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i);
  
  

  template <typename OutPoint, typename InPoint, typename RateLimitMap>
  void create_normal_joint_vector_impl(OutPoint& result,
				   const InPoint& pt,
				   const RateLimitMap& j_limits,
				   std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i) {
    create_normal_joint_vectors_impl< typename boost::mpl::prior< arithmetic_tuple_size< OutPoint > >::type >(result, pt, j_limits, gen_i, f2d_i, f3d_i);
  };
  
  
  
  
  template <typename Idx, typename OutPoint, typename InPoint, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_normal_joint_vectors_impl(OutPoint& result,
				            const InPoint& pt,
					    const RateLimitMap& j_limits,
					    std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i) {
    
    create_normal_joint_vectors_impl< typename boost::mpl::prior<Idx>::type >(result,pt,j_limits,gen_i,f2d_i,f3d_i);
    
    create_normal_joint_vector_impl(get< Idx::type::value >(result), 
                                get< Idx::type::value >(pt),
			        j_limits, gen_i, f2d_i, f3d_i);
    
  };
  
  template <typename Idx, typename OutPoint, typename InPoint, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_normal_joint_vectors_impl(OutPoint& result,
				            const InPoint& pt,
					    const RateLimitMap& j_limits,
					    std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i) {
    
    create_normal_joint_vectors_impl(get< 0 >(result), 
                                get< 0 >(pt),
			        j_limits, gen_i, f2d_i, f3d_i);
    
  };
  
  template <typename OutPoint, typename InPoint, typename RateLimitMap>
  void create_normal_joint_vectors_impl(OutPoint& result,
				    const InPoint& pt,
				    const RateLimitMap& j_limits) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    create_rl_joint_vectors_impl< typename boost::mpl::prior< arithmetic_tuple_size< OutPoint > >::type >(result, pt, j_limits, gen_i, f2d_i, f3d_i);
  };

  
  
  
};





/**
 * This class template stores a set of vectors to represent the rate-limits on the joints 
 * of a manipulator. Basically, this class is just a POD class, but it also provides functions
 * to construct a rate-limited joint-space from a normal joint-space, or vice-versa. Also, 
 * it can act as a mapping between rate-limited joint coordinates and normal joint coordinates.
 * \tparam T The value type of the underlying joint-space.
 */
template <typename T>
struct joint_limits_collection : public named_object {
  /** Holds the speed limit for all generalized coordinates. */
  vect_n<T> gen_speed_limits;
  /** Holds the acceleration limit for all generalized coordinates. */
  vect_n<T> gen_accel_limits;
  /** Holds the jerk limit for all generalized coordinates. */
  vect_n<T> gen_jerk_limits;
  /** Holds the speed limit for all 2D frames (alternating velocity limit and angular velocity limit). */
  vect_n<T> frame2D_speed_limits;
  /** Holds the acceleration limit for all 2D frames (alternating acceleration limit and angular acceleration limit). */
  vect_n<T> frame2D_accel_limits;
  /** Holds the jerk limit for all 2D frames (alternating jerk limit and angular jerk limit). */
  vect_n<T> frame2D_jerk_limits;
  /** Holds the speed limit for all 3D frames (alternating velocity limit and angular velocity limit). */
  vect_n<T> frame3D_speed_limits;
  /** Holds the acceleration limit for all 3D frames (alternating acceleration limit and angular acceleration limit). */
  vect_n<T> frame3D_accel_limits;
  /** Holds the jerk limit for all 3D frames (alternating jerk limit and angular jerk limit). */
  vect_n<T> frame3D_jerk_limits;
  
  typedef T value_type;
  typedef joint_limits_collection<T> self;
  
  /**
   * Default constructor.
   */
  joint_limits_collection(const std::string& aName = "") : named_object() {
    this->setName(aName);
  };
  
  /**
   * This function constructs a rate-limited joint-space out of the given normal joint-space.
   * \tparam NormalSpaceType The topology type of the joint-space.
   * \param j_space The normal joint-space.
   * \return A rate-limited joint-space corresponding to given joint-space and the stored limit values.
   */
  template <typename NormalSpaceType>
  typename get_rate_limited_space< NormalSpaceType >::type make_rl_joint_space(const NormalSpaceType& j_space) const {
    typename get_rate_limited_space< NormalSpaceType >::type result;
    detail::create_rl_joint_spaces_impl(result, j_space, *this);
    return result;
  };
  
  /**
   * This function constructs a normal joint-space out of the given rate-limited joint-space.
   * \tparam RateLimitedSpaceType The topology type of the rate-limited joint-space.
   * \param j_space The rate-limited joint-space.
   * \return A normal joint-space corresponding to given rate-limited joint-space and the stored limit values.
   */
  template <typename RateLimitedSpaceType>
  typename get_rate_illimited_space< RateLimitedSpaceType >::type make_normal_joint_space(const RateLimitedSpaceType& j_space) const {
    typename get_rate_illimited_space< RateLimitedSpaceType >::type result;
    detail::create_normal_joint_spaces_impl(result, j_space, *this);
    return result;
  };
  
  /**
   * This function maps a set of normal joint coordinates into a set of rate-limited joint coordinates.
   * \tparam NormalSpaceType The topology type of the joint-space.
   * \param pt A point in the normal joint-space.
   * \param j_space The normal joint-space.
   * \param rl_j_space The rate-limited joint-space (in which the output lies).
   * \return A set of rate-limited joint coordinates corresponding to given normal joint coordinates and the stored limit values.
   */
  template <typename NormalSpaceType>
  typename topology_traits< typename get_rate_limited_space< NormalSpaceType >::type >::point_type map_to_space(
      const typename topology_traits< NormalSpaceType >::point_type& pt,
      const NormalSpaceType& , const typename get_rate_limited_space< NormalSpaceType >::type& ) const {
    typename topology_traits< typename get_rate_limited_space< NormalSpaceType >::type >::point_type result;
    detail::create_rl_joint_vectors_impl(result, pt, *this);
    return result;
  };
  
  
  /**
   * This function maps a set of rate-limited joint coordinates into a set of normal joint coordinates.
   * \tparam RateLimitedSpaceType The topology type of the rate-limited joint-space.
   * \param pt A point in the rate-limited joint-space.
   * \param j_space The rate-limited joint-space.
   * \param rl_j_space The normal joint-space (in which the output lies).
   * \return A set of normal joint coordinates corresponding to given rate-limited joint coordinates and the stored limit values.
   */
  template <typename RateLimitedSpaceType>
  typename topology_traits< typename get_rate_illimited_space< RateLimitedSpaceType >::type >::point_type map_to_space(
      const typename topology_traits< RateLimitedSpaceType >::point_type& pt,
      const RateLimitedSpaceType& , const typename get_rate_illimited_space< RateLimitedSpaceType >::type& ) const {
    typename topology_traits< typename get_rate_illimited_space< RateLimitedSpaceType >::type >::point_type result;
    detail::create_normal_joint_vectors_impl(result, pt, *this);
    return result;
  };
  
  
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(gen_speed_limits)
        & RK_SERIAL_SAVE_WITH_NAME(gen_accel_limits)
        & RK_SERIAL_SAVE_WITH_NAME(gen_jerk_limits)
        & RK_SERIAL_SAVE_WITH_NAME(frame2D_speed_limits)
        & RK_SERIAL_SAVE_WITH_NAME(frame2D_accel_limits)
        & RK_SERIAL_SAVE_WITH_NAME(frame2D_jerk_limits)
        & RK_SERIAL_SAVE_WITH_NAME(frame3D_speed_limits)
        & RK_SERIAL_SAVE_WITH_NAME(frame3D_accel_limits)
        & RK_SERIAL_SAVE_WITH_NAME(frame3D_jerk_limits);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(gen_speed_limits)
        & RK_SERIAL_LOAD_WITH_NAME(gen_accel_limits)
        & RK_SERIAL_LOAD_WITH_NAME(gen_jerk_limits)
        & RK_SERIAL_LOAD_WITH_NAME(frame2D_speed_limits)
        & RK_SERIAL_LOAD_WITH_NAME(frame2D_accel_limits)
        & RK_SERIAL_LOAD_WITH_NAME(frame2D_jerk_limits)
        & RK_SERIAL_LOAD_WITH_NAME(frame3D_speed_limits)
        & RK_SERIAL_LOAD_WITH_NAME(frame3D_accel_limits)
        & RK_SERIAL_LOAD_WITH_NAME(frame3D_jerk_limits);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400011,1,"joint_limits_collection",named_object)
    
  
};




};



};

#endif









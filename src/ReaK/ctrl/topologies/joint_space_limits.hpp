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

#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "joint_space_topologies.hpp"

namespace ReaK {

namespace pp {

  
  
namespace detail {
  
  
  
  template <typename Idx, typename T, std::size_t N, typename DistanceMetric>
  typename boost::disable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_0th_rl_joint_space_impl(
    typename metric_space_array<
      rl_joint_space_0th_order<T>::type,
      N,
      DistanceMetric >::type& result,
    const typename metric_space_array<
      joint_space_0th_order<T>::type,
      N,
      DistanceMetric >::type& j_space,
    const vect<T,N>& speed_limits) {
    
    create_Oth_rl_joint_space_impl< boost::mpl::prior<Idx>, T, N, DistanceMetric >(result,j_space,speed_limits);
    
    get< Idx::type::value >(result) 
      = reach_time_diff_space< 
          time_topology, 
          arithmetic_tuple< 
            line_segment_topology<T>
          >, 
          euclidean_tuple_distance 
        >(arithmetic_tuple< 
            line_segment_topology<T>
          >(line_segment_topology<T>(
	      get<0>(get< Idx::type::value >(j_space)).getName() + "_rl",
	      ( get<0>(get< Idx::type::value >(j_space)).origin() - get<0>(get< Idx::type::value >(j_space)).get_radius() ) / speed_limits[Idx::type::value],
	      ( get<0>(get< Idx::type::value >(j_space)).origin() + get<0>(get< Idx::type::value >(j_space)).get_radius() ) / speed_limits[Idx::type::value]
	    )
	  )
 	);
    
  };
  
  template <typename Idx, typename T, std::size_t N, typename DistanceMetric>
  typename boost::enable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_0th_rl_joint_space_impl(
    typename metric_space_array<
      rl_joint_space_0th_order<T>::type,
      N,
      DistanceMetric >::type& result,
    const typename metric_space_array<
      joint_space_0th_order<T>::type,
      N,
      DistanceMetric >::type& j_space,
    const vect<T,N>& speed_limits) {
    
    get< 0 >(result) 
      = reach_time_diff_space< 
          time_topology, 
          arithmetic_tuple< 
            line_segment_topology<T>
          >, 
          euclidean_tuple_distance 
        >(arithmetic_tuple< 
            line_segment_topology<T>
          >(line_segment_topology<T>(
	      get<0>(get< 0 >(j_space)).getName() + "_rl",
	      ( get<0>(get< 0 >(j_space)).origin() - get<0>(get< 0 >(j_space)).get_radius() ) / speed_limits[0],
	      ( get<0>(get< 0 >(j_space)).origin() + get<0>(get< 0 >(j_space)).get_radius() ) / speed_limits[0]
	    )
	  )
 	);
    
  };
  
  
  
  
  
  template <typename Idx, typename T, std::size_t N, typename DistanceMetric>
  typename boost::disable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_1st_rl_joint_space_impl(
    typename metric_space_array<
      rl_joint_space_1st_order<T>::type,
      N,
      DistanceMetric >::type& result,
    const typename metric_space_array<
      joint_space_1st_order<T>::type,
      N,
      DistanceMetric >::type& j_space,
    const vect<T,N>& speed_limits,
    const vect<T,N>& accel_limits) {
    
    create_1st_rl_joint_space_impl< boost::mpl::prior<Idx>, T, N, DistanceMetric >(result,j_space,speed_limits,accel_limits);
    
    get< Idx::type::value >(result) 
      = reach_time_diff_space< 
          time_topology, 
          arithmetic_tuple< 
            line_segment_topology<T>,
            line_segment_topology<T>
          >, 
          euclidean_tuple_distance 
        >(arithmetic_tuple< 
            line_segment_topology<T>,
            line_segment_topology<T>
          >(line_segment_topology<T>(
	      get<0>(get< Idx::type::value >(j_space)).getName() + "_rl",
	      ( get<0>(get< Idx::type::value >(j_space)).origin() - get<0>(get< Idx::type::value >(j_space)).get_radius() ) / speed_limits[Idx::type::value],
	      ( get<0>(get< Idx::type::value >(j_space)).origin() + get<0>(get< Idx::type::value >(j_space)).get_radius() ) / speed_limits[Idx::type::value]
	    ),
	    line_segment_topology<T>(
	      get<1>(get< Idx::type::value >(j_space)).getName() + "_rl",
	      ( get<1>(get< Idx::type::value >(j_space)).origin() - get<1>(get< Idx::type::value >(j_space)).get_radius() ) / accel_limits[Idx::type::value],
	      ( get<1>(get< Idx::type::value >(j_space)).origin() + get<1>(get< Idx::type::value >(j_space)).get_radius() ) / accel_limits[Idx::type::value]
	    )
	  ),
	  euclidean_tuple_distance(),
	  reach_time_differentiation_tuple< 1 >::type(
	    reach_time_differentiation( speed_limits[ Idx::type::value ] / accel_limits[ Idx::type::value ] )
	  )
 	);
    
  };
  
  template <typename Idx, typename T, std::size_t N, typename DistanceMetric>
  typename boost::enable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_1st_rl_joint_space_impl(
    typename metric_space_array<
      rl_joint_space_1st_order<T>::type,
      N,
      DistanceMetric >::type& result,
    const typename metric_space_array<
      joint_space_1st_order<T>::type,
      N,
      DistanceMetric >::type& j_space,
    const vect<T,N>& speed_limits,
    const vect<T,N>& accel_limits) {
    
    get< 0 >(result) 
      = reach_time_diff_space< 
          time_topology, 
          arithmetic_tuple< 
            line_segment_topology<T>,
            line_segment_topology<T>
          >, 
          euclidean_tuple_distance 
        >(arithmetic_tuple< 
            line_segment_topology<T>,
            line_segment_topology<T>
          >(line_segment_topology<T>(
	      get<0>(get< 0 >(j_space)).getName() + "_rl",
	      ( get<0>(get< 0 >(j_space)).origin() - get<0>(get< 0 >(j_space)).get_radius() ) / speed_limits[0],
	      ( get<0>(get< 0 >(j_space)).origin() + get<0>(get< 0 >(j_space)).get_radius() ) / speed_limits[0]
	    ),
	    line_segment_topology<T>(
	      get<1>(get< 0 >(j_space)).getName() + "_rl",
	      ( get<1>(get< 0 >(j_space)).origin() - get<1>(get< 0 >(j_space)).get_radius() ) / accel_limits[0],
	      ( get<1>(get< 0 >(j_space)).origin() + get<1>(get< 0 >(j_space)).get_radius() ) / accel_limits[0]
	    )
	  ),
	  euclidean_tuple_distance(),
	  reach_time_differentiation_tuple< 1 >::type(
	    reach_time_differentiation( speed_limits[ 0 ] / accel_limits[ 0 ] )
	  )
 	);
    
  };
  
  
  
  
  
  
  
  template <typename Idx, typename T, std::size_t N, typename DistanceMetric>
  typename boost::disable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_2nd_rl_joint_space_impl(
    typename metric_space_array<
      rl_joint_space_2nd_order<T>::type,
      N,
      DistanceMetric >::type& result,
    const typename metric_space_array<
      joint_space_2nd_order<T>::type,
      N,
      DistanceMetric >::type& j_space,
    const vect<T,N>& speed_limits,
    const vect<T,N>& accel_limits,
    const vect<T,N>& jerk_limits) {
    
    create_2nd_rl_joint_space_impl< boost::mpl::prior<Idx>, T, N, DistanceMetric >(result,j_space,speed_limits,accel_limits,jerk_limits);
    
    get< Idx::type::value >(result) 
      = reach_time_diff_space< 
          time_topology, 
          arithmetic_tuple< 
            line_segment_topology<T>,
            line_segment_topology<T>,
            line_segment_topology<T>
          >, 
          euclidean_tuple_distance 
        >(arithmetic_tuple< 
            line_segment_topology<T>,
            line_segment_topology<T>,
            line_segment_topology<T>
          >(line_segment_topology<T>(
	      get<0>(get< Idx::type::value >(j_space)).getName() + "_rl",
	      ( get<0>(get< Idx::type::value >(j_space)).origin() - get<0>(get< Idx::type::value >(j_space)).get_radius() ) / speed_limits[Idx::type::value],
	      ( get<0>(get< Idx::type::value >(j_space)).origin() + get<0>(get< Idx::type::value >(j_space)).get_radius() ) / speed_limits[Idx::type::value]
	    ),
	    line_segment_topology<T>(
	      get<1>(get< Idx::type::value >(j_space)).getName() + "_rl",
	      ( get<1>(get< Idx::type::value >(j_space)).origin() - get<1>(get< Idx::type::value >(j_space)).get_radius() ) / accel_limits[Idx::type::value],
	      ( get<1>(get< Idx::type::value >(j_space)).origin() + get<1>(get< Idx::type::value >(j_space)).get_radius() ) / accel_limits[Idx::type::value]
	    ),
	    line_segment_topology<T>(
	      get<2>(get< Idx::type::value >(j_space)).getName() + "_rl",
	      ( get<2>(get< Idx::type::value >(j_space)).origin() - get<2>(get< Idx::type::value >(j_space)).get_radius() ) / jerk_limits[Idx::type::value],
	      ( get<2>(get< Idx::type::value >(j_space)).origin() + get<2>(get< Idx::type::value >(j_space)).get_radius() ) / jerk_limits[Idx::type::value]
	    )
	  ),
	  euclidean_tuple_distance(),
	  reach_time_differentiation_tuple< 2 >::type(
	    reach_time_differentiation( speed_limits[ Idx::type::value ] / accel_limits[ Idx::type::value ] ),
	    reach_time_differentiation( accel_limits[ Idx::type::value ] / jerk_limits[ Idx::type::value ] )
	  )
 	);
    
  };
  
  template <typename Idx, typename T, std::size_t N, typename DistanceMetric>
  typename boost::enable_if< 
    boost::mpl::less<
      Idx,
      boost::mpl::size_t<1>
    >,
  void >::type create_2nd_rl_joint_space_impl(
    typename metric_space_array<
      rl_joint_space_2nd_order<T>::type,
      N,
      DistanceMetric >::type& result,
    const typename metric_space_array<
      joint_space_2nd_order<T>::type,
      N,
      DistanceMetric >::type& j_space,
    const vect<T,N>& speed_limits,
    const vect<T,N>& accel_limits,
    const vect<T,N>& jerk_limits) {
    
    get< 0 >(result) 
      = reach_time_diff_space< 
          time_topology, 
          arithmetic_tuple< 
            line_segment_topology<T>,
            line_segment_topology<T>,
            line_segment_topology<T>
          >, 
          euclidean_tuple_distance 
        >(arithmetic_tuple< 
            line_segment_topology<T>,
            line_segment_topology<T>,
            line_segment_topology<T>
          >(line_segment_topology<T>(
	      get<0>(get< 0 >(j_space)).getName() + "_rl",
	      ( get<0>(get< 0 >(j_space)).origin() - get<0>(get< 0 >(j_space)).get_radius() ) / speed_limits[0],
	      ( get<0>(get< 0 >(j_space)).origin() + get<0>(get< 0 >(j_space)).get_radius() ) / speed_limits[0]
	    ),
	    line_segment_topology<T>(
	      get<1>(get< 0 >(j_space)).getName() + "_rl",
	      ( get<1>(get< 0 >(j_space)).origin() - get<1>(get< 0 >(j_space)).get_radius() ) / accel_limits[0],
	      ( get<1>(get< 0 >(j_space)).origin() + get<1>(get< 0 >(j_space)).get_radius() ) / accel_limits[0]
	    ),
	    line_segment_topology<T>(
	      get<2>(get< 0 >(j_space)).getName() + "_rl",
	      ( get<2>(get< 0 >(j_space)).origin() - get<2>(get< 0 >(j_space)).get_radius() ) / jerk_limits[0],
	      ( get<2>(get< 0 >(j_space)).origin() + get<2>(get< 0 >(j_space)).get_radius() ) / jerk_limits[0]
	    )
	  ),
	  euclidean_tuple_distance(),
	  reach_time_differentiation_tuple< 2 >::type(
	    reach_time_differentiation( speed_limits[ 0 ] / accel_limits[ 0 ] ),
	    reach_time_differentiation( accel_limits[ 0 ] / jerk_limits[ 0 ] )
	  )
 	);
    
  };
  
  
  
  
};




template <typename T, std::size_t N, typename DistanceMetric = inf_norm_tuple_distance>
struct joint_limits_1st_order {
  vect<T, N> speed_limits;
  
  joint_limits_1st_order(const vect<T,N>& aSpeedLimits) : speed_limits(aSpeedLimits) { };
  
  typename metric_space_array<
    rl_joint_space_0th_order<T>::type,
    N,
    DistanceMetric >::type make_rl_joint_space(
      const typename metric_space_array<
        joint_space_0th_order<T>::type,
	N,
	DistanceMetric >::type& j_space
    ) {
    typename metric_space_array<
      rl_joint_space_0th_order<T>::type,
      N,
      DistanceMetric >::type result;
    detail::create_0th_rl_joint_space_impl<boost::mpl::size_t<N-1>, T, N, DistanceMetric >(result, j_space, speed_limits);
    return result;
  };
  
};




template <typename T, std::size_t N, typename DistanceMetric = inf_norm_tuple_distance>
struct joint_limits_2nd_order {
  vect<T, N> speed_limits;
  vect<T, N> accel_limits;
  
  joint_limits_2nd_order(const vect<T,N>& aSpeedLimits, 
                         const vect<T,N>& aAccelLimits) : 
                         speed_limits(aSpeedLimits), 
                         accel_limits(aAccelLimits) { };
  
  typename metric_space_array<
    rl_joint_space_1st_order<T>::type,
    N,
    DistanceMetric >::type make_rl_joint_space(
      const typename metric_space_array<
        joint_space_1st_order<T>::type,
	N,
	DistanceMetric >::type& j_space
    ) {
    typename metric_space_array<
      rl_joint_space_1st_order<T>::type,
      N,
      DistanceMetric >::type result;
    detail::create_1st_rl_joint_space_impl<boost::mpl::size_t<N-1>, T, N, DistanceMetric >(result, j_space, speed_limits);
    return result;
  };
  
};







template <typename T, std::size_t N, typename DistanceMetric = inf_norm_tuple_distance>
struct joint_limits_3rd_order {
  vect<T, N> speed_limits;
  vect<T, N> accel_limits;
  vect<T, N> jerk_limits;
  
  joint_limits_3rd_order(const vect<T,N>& aSpeedLimits, 
                         const vect<T,N>& aAccelLimits, 
                         const vect<T,N>& aJerkLimits) : 
                         speed_limits(aSpeedLimits), 
                         accel_limits(aAccelLimits), 
                         jerk_limits(aJerkLimits) { };
  
  typename metric_space_array<
    rl_joint_space_2nd_order<T>::type,
    N,
    DistanceMetric >::type make_rl_joint_space(
      const typename metric_space_array<
        joint_space_2nd_order<T>::type,
	N,
	DistanceMetric >::type& j_space
    ) {
    typename metric_space_array<
      rl_joint_space_2nd_order<T>::type,
      N,
      DistanceMetric >::type result;
    detail::create_2nd_rl_joint_space_impl<boost::mpl::size_t<N-1>, T, N, DistanceMetric >(result, j_space, speed_limits);
    return result;
  };
  
};



};



};

#endif









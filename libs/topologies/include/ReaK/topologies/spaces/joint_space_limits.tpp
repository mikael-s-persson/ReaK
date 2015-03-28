/**
 * \file joint_space_limits.tpp
 *
 * This library provides classes to help create and manipulate joint-space topologies in over
 * a joint-space with limits (speed, acceleration, and jerk limits).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2012
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

#ifndef REAK_JOINT_SPACE_LIMITS_TPP
#define REAK_JOINT_SPACE_LIMITS_TPP

#include "joint_space_limits.hpp"
#include "joint_space_limits_detail.hpp"


namespace ReaK {

namespace pp {


template < typename T >
template < typename NormalSpaceType >
typename get_rate_limited_space< NormalSpaceType >::type
  joint_limits_mapping< T >::make_rl_joint_space( const NormalSpaceType& j_space ) const {
  typename get_rate_limited_space< NormalSpaceType >::type result;
  detail::create_rl_joint_spaces_impl( result, j_space, *this->limits );
  return result;
};

template < typename T >
template < typename RateLimitedSpaceType >
typename get_rate_illimited_space< RateLimitedSpaceType >::type
  joint_limits_mapping< T >::make_normal_joint_space( const RateLimitedSpaceType& j_space ) const {
  typename get_rate_illimited_space< RateLimitedSpaceType >::type result;
  detail::create_normal_joint_spaces_impl( result, j_space, *this->limits );
  return result;
};

template < typename T >
template < typename NormalSpaceType >
typename topology_traits< typename get_rate_limited_space< NormalSpaceType >::type >::point_type
  joint_limits_mapping< T >::map_to_space( const typename topology_traits< NormalSpaceType >::point_type& pt,
                                           const NormalSpaceType&,
                                           const typename get_rate_limited_space< NormalSpaceType >::type& ) const {
  typename topology_traits< typename get_rate_limited_space< NormalSpaceType >::type >::point_type result;
  detail::create_rl_joint_vectors_impl( result, pt, *this->limits );
  return result;
};


template < typename T >
template < typename RateLimitedSpaceType >
typename topology_traits< typename get_rate_illimited_space< RateLimitedSpaceType >::type >::point_type
  joint_limits_mapping< T >::map_to_space(
    const typename topology_traits< RateLimitedSpaceType >::point_type& pt, const RateLimitedSpaceType&,
    const typename get_rate_illimited_space< RateLimitedSpaceType >::type& ) const {
  typename topology_traits< typename get_rate_illimited_space< RateLimitedSpaceType >::type >::point_type result;
  detail::create_normal_joint_vectors_impl( result, pt, *this->limits );
  return result;
};
};
};


#endif

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

#include "ReaK/topologies/spaces/joint_space_limits.hpp"
#include "ReaK/topologies/spaces/joint_space_limits_detail.hpp"

namespace ReaK::pp {

template <typename T>
template <typename NormalSpaceType>
get_rate_limited_space_t<NormalSpaceType>
joint_limits_mapping<T>::make_rl_joint_space(
    const NormalSpaceType& j_space) const {
  get_rate_limited_space_t<NormalSpaceType> result;
  detail::create_rl_joint_spaces_impl(result, j_space, *this->limits);
  return result;
}

template <typename T>
template <typename RateLimitedSpaceType>
get_rate_illimited_space_t<RateLimitedSpaceType>
joint_limits_mapping<T>::make_normal_joint_space(
    const RateLimitedSpaceType& j_space) const {
  get_rate_illimited_space_t<RateLimitedSpaceType> result;
  detail::create_normal_joint_spaces_impl(result, j_space, *this->limits);
  return result;
}

template <typename T>
template <typename NormalSpaceType>
topology_point_type_t<get_rate_limited_space_t<NormalSpaceType>>
joint_limits_mapping<T>::map_to_space(
    const topology_point_type_t<NormalSpaceType>& pt, const NormalSpaceType& /*unused*/,
    const get_rate_limited_space_t<NormalSpaceType>& /*unused*/) const {
  topology_point_type_t<get_rate_limited_space_t<NormalSpaceType>> result;
  detail::create_rl_joint_vectors_impl(result, pt, *this->limits);
  return result;
}

template <typename T>
template <typename RateLimitedSpaceType>
topology_point_type_t<get_rate_illimited_space_t<RateLimitedSpaceType>>
joint_limits_mapping<T>::map_to_space(
    const topology_point_type_t<RateLimitedSpaceType>& pt,
    const RateLimitedSpaceType& /*unused*/,
    const get_rate_illimited_space_t<RateLimitedSpaceType>& /*unused*/) const {
  topology_point_type_t<get_rate_illimited_space_t<RateLimitedSpaceType>>
      result;
  detail::create_normal_joint_vectors_impl(result, pt, *this->limits);
  return result;
}

}  // namespace ReaK::pp

#endif

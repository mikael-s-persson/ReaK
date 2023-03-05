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

#ifndef REAK_JOINT_SPACE_LIMITS_MAPPING_HPP
#define REAK_JOINT_SPACE_LIMITS_MAPPING_HPP

#include "ReaK/core/base/defs.hpp"
#include "ReaK/core/base/named_object.hpp"

#include "ReaK/mbd/models/joint_space_limits.hpp"

#include "ReaK/topologies/spaces/metric_space_concept.hpp"
#include "ReaK/topologies/spaces/rate_limited_space_metamaps.hpp"

namespace ReaK::pp {

/**
 * This class template stores a set of vectors to represent the rate-limits on the joints
 * of a manipulator. Basically, this class is just a POD class, but it also provides functions
 * to construct a rate-limited joint-space from a normal joint-space, or vice-versa. Also,
 * it can act as a mapping between rate-limited joint coordinates and normal joint coordinates.
 * \tparam T The value type of the underlying joint-space.
 */
template <typename T>
struct joint_limits_mapping : public named_object {
  /** Holds the all limits. */
  std::shared_ptr<kte::joint_limits_collection<T>> limits;

  using value_type = T;
  using self = joint_limits_mapping<T>;

  /**
   * Default constructor.
   */
  explicit joint_limits_mapping(
      const std::shared_ptr<kte::joint_limits_collection<T>>& aLimits = {})
      : named_object(), limits(aLimits) {
    if (limits) {
      this->setName(limits->getName() + "_mapping");
    } else {
      this->setName("joint_limits_mapping");
    }
  }

  /**
   * This function constructs a rate-limited joint-space out of the given normal joint-space.
   * \tparam NormalSpaceType The topology type of the joint-space.
   * \param j_space The normal joint-space.
   * \return A rate-limited joint-space corresponding to given joint-space and the stored limit values.
   */
  template <typename NormalSpaceType>
  get_rate_limited_space_t<NormalSpaceType> make_rl_joint_space(
      const NormalSpaceType& j_space) const;

  /**
   * This function constructs a normal joint-space out of the given rate-limited joint-space.
   * \tparam RateLimitedSpaceType The topology type of the rate-limited joint-space.
   * \param j_space The rate-limited joint-space.
   * \return A normal joint-space corresponding to given rate-limited joint-space and the stored limit values.
   */
  template <typename RateLimitedSpaceType>
  get_rate_illimited_space_t<RateLimitedSpaceType> make_normal_joint_space(
      const RateLimitedSpaceType& j_space) const;

  /**
   * This function maps a set of normal joint coordinates into a set of rate-limited joint coordinates.
   * \tparam NormalSpaceType The topology type of the joint-space.
   * \param pt A point in the normal joint-space.
   * \param j_space The normal joint-space.
   * \param rl_j_space The rate-limited joint-space (in which the output lies).
   * \return A set of rate-limited joint coordinates corresponding to given normal joint coordinates and the stored
   * limit values.
   */
  template <typename NormalSpaceType>
  topology_point_type_t<get_rate_limited_space_t<NormalSpaceType>> map_to_space(
      const topology_point_type_t<NormalSpaceType>& pt,
      const NormalSpaceType& /*unused*/,
      const get_rate_limited_space_t<NormalSpaceType>& /*unused*/) const;

  /**
   * This function maps a set of rate-limited joint coordinates into a set of normal joint coordinates.
   * \tparam RateLimitedSpaceType The topology type of the rate-limited joint-space.
   * \param pt A point in the rate-limited joint-space.
   * \param j_space The rate-limited joint-space.
   * \param rl_j_space The normal joint-space (in which the output lies).
   * \return A set of normal joint coordinates corresponding to given rate-limited joint coordinates and the stored
   * limit values.
   */
  template <typename RateLimitedSpaceType>
  topology_point_type_t<get_rate_illimited_space_t<RateLimitedSpaceType>>
  map_to_space(
      const topology_point_type_t<RateLimitedSpaceType>& pt,
      const RateLimitedSpaceType& /*unused*/,
      const get_rate_illimited_space_t<RateLimitedSpaceType>& /*unused*/) const;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(limits);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(limits);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240002F, 1, "joint_limits_mapping",
                              named_object)
};

}  // namespace ReaK::pp

#endif

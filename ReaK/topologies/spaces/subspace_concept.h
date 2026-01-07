/**
 * \file subspace_concept.h
 *
 * This library defines the traits and concepts that pertain to what can be considered
 * a sub-space, as used in ReaK::pp. Sub-spaces are topologies which are embedded in a
 * larger topology (or super-space). A typical example is the free configuration space (e.g. C-free)
 * embedded in the overall configuration space. It is useful to be able to relate a sub-space
 * to its super-space, for example, if a motion planner plans a trajectory through C-free, it
 * is useful to cast that trajectory onto the configuration space such that collision detections
 * are avoided when executing the trajectory.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
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

#ifndef REAK_TOPOLOGIES_SPACES_SUBSPACE_CONCEPT_H_
#define REAK_TOPOLOGIES_SPACES_SUBSPACE_CONCEPT_H_

#include "ReaK/core/base/shared_object.h"

#include <concepts>

/** Main namespace for ReaK */
namespace ReaK::pp {

/**
 * This traits class defines the types and constants associated to a sub-space.
 */
template <Topology Space>
struct subspace_traits {
  /** The type of the super-space in which this sub-space is embedded. */
  using super_space_type = typename Space::super_space_type;
};

template <Topology Space>
using super_space_t = typename subspace_traits<Space>::super_space_type;

/**
 * This concept defines the requirements to fulfill in order to model a sub-space
 * as used in ReaK::pp.
 *
 * Valid expressions:
 *
 * super_space = space.get_super_space();  The super-space can be obtained (at least, by const-reference).
 *
 * \tparam Topology The topology type to be checked for this concept.
 */
template <typename Space>
concept SubSpace = Topology<Space>&& requires(const Space& space) {
  {
    space.get_super_space()
    } -> std::convertible_to<const super_space_t<Space>&>;
};

struct subspace_map : public shared_object {
  using self = subspace_map;

  subspace_map() = default;

  template <typename PointType, Topology SpaceIn>
  PointType map_to_space(const PointType& p_in, const SpaceIn& /*unused*/,
                         const super_space_t<SpaceIn>& /*unused*/) const {
    return p_in;
  }

  template <typename PointType, Topology SpaceOut>
  PointType map_to_space(const PointType& p_in,
                         const super_space_t<SpaceOut>& /*unused*/,
                         const SpaceOut& /*unused*/) const {
    return p_in;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240002C, 1, "subspace_map",
                              shared_object)
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_SUBSPACE_CONCEPT_H_

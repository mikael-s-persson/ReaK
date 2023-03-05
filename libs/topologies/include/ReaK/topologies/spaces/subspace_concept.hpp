/**
 * \file subspace_concept.hpp
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

#ifndef REAK_SUBSPACE_CONCEPT_HPP
#define REAK_SUBSPACE_CONCEPT_HPP

#include "ReaK/core/base/defs.hpp"
#include "ReaK/core/base/shared_object.hpp"

#include "boost/concept_check.hpp"

/** Main namespace for ReaK */
namespace ReaK::pp {

/**
 * This traits class defines the types and constants associated to a sub-space.
 * \tparam Topology The topology type for which the topology traits are sought.
 */
template <typename Topology>
struct subspace_traits {
  /** The type of the super-space in which this sub-space is embedded. */
  using super_space_type = typename Topology::super_space_type;
};

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
template <typename Topology>
struct SubSpaceConcept {
  Topology space;

  BOOST_CONCEPT_USAGE(SubSpaceConcept) {
    const typename subspace_traits<Topology>::super_space_type& super_space =
        space.get_super_space();
    RK_UNUSED(super_space);
  }
};

struct subspace_map : public shared_object {
  using self = subspace_map;

  subspace_map() = default;
  ;

  template <typename PointType, typename SpaceIn>
  PointType map_to_space(
      const PointType& p_in, const SpaceIn& /*unused*/,
      const typename subspace_traits<SpaceIn>::super_space_type& /*unused*/)
      const {
    return p_in;
  }

  template <typename PointType, typename SpaceOut>
  PointType map_to_space(
      const PointType& p_in,
      const typename subspace_traits<SpaceOut>::super_space_type& /*unused*/,
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

#endif

/**
 * \file augmented_to_state_mapping.hpp
 *
 * This library provides a class template which can map augmented-state points into state-points of a
 * compatibly state-space topology.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_AUGMENTED_TO_STATE_MAPPING_HPP
#define REAK_AUGMENTED_TO_STATE_MAPPING_HPP

#include <ReaK/core/base/named_object.hpp>

#include <ReaK/topologies/spaces/metric_space_concept.hpp>
#include <ReaK/topologies/spaces/temporal_space_concept.hpp>

#include <ReaK/math/lin_alg/arithmetic_tuple.hpp>

#include <type_traits>

namespace ReaK::ctrl {

/**
 * This class is a bijection mapping from a belief-space to a state-space (topology)
 * by making the assumption of maximum likelihood. In other words, this mapping reduces
 * belief-states into their most likely value only.
 *
 * Models: BijectionConcept between a belief-state topology and a compatible state topology.
 */
struct augmented_to_state_map : public named_object {

  using self = augmented_to_state_map;

  augmented_to_state_map() { setName("augmented_to_state_map"); }

  /**
   * This function extracts the state from a augmented-state.
   * \tparam AugStateSpace A augmented-state topoloogy whose state-space is compatible with the given state-space
   * topology.
   * \tparam StateSpaceOut A state-space topology whose point-types are compatible with the state type that the
   * augmented-state type would produce.
   * \param b The augmented-state from which the state is sought.
   * \return The state of the augmented-state.
   */
  template <typename AugStateSpace, typename StateSpaceOut>
  pp::topology_point_type_t<StateSpaceOut> map_to_space(
      const pp::topology_point_type_t<AugStateSpace>& b,
      const AugStateSpace& /*unused*/, const StateSpaceOut& /*unused*/) const {
    if constexpr (pp::is_temporal_space_v<AugStateSpace>) {
      return {b.time, get<0>(b.pt)};
    } else {
      return {get<0>(b)};
    }
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2300017, 1, "augmented_to_state_map",
                              named_object)
};

}  // namespace ReaK::ctrl

#endif

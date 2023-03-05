/**
 * \file maximum_likelihood_mapping.hpp
 *
 * This library provides a class template which can map belief-state points into points of a
 * compatibly state-space topology.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2011
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

#ifndef REAK_MAXIMUM_LIKELIHOOD_MAPPING_HPP
#define REAK_MAXIMUM_LIKELIHOOD_MAPPING_HPP

#include "ReaK/core/base/named_object.hpp"

#include "ReaK/topologies/spaces/metric_space_concept.hpp"
#include "ReaK/topologies/spaces/temporal_space_concept.hpp"

#include "ReaK/control/estimators/belief_state_concept.hpp"

#include <type_traits>

namespace ReaK::ctrl {

/**
 * This class is a bijection mapping from a belief-space to a state-space (topology)
 * by making the assumption of maximum likelihood. In other words, this mapping reduces
 * belief-states into their most likely value only.
 *
 * Models: BijectionConcept between a belief-state topology and a compatible state topology.
 */
struct maximum_likelihood_map : public named_object {

  using self = maximum_likelihood_map;

  maximum_likelihood_map() { setName("maximum_likelihood_map"); }

  /**
   * This function extracts the most-likely-value from a belief-state.
   * \tparam BeliefSpace A belief-state topoloogy whose state-space is compatible with the given state-space topology.
   * \tparam StateSpaceOut A state-space topology whose point-types are compatible with the state type that the
   * belief-state type would produce.
   * \param b The belief-state from which the maximum likelihood value is sought.
   * \return The maximum likelihood value of the belief-state.
   */
  template <typename BeliefSpace, typename StateSpaceOut>
  pp::topology_point_type_t<StateSpaceOut> map_to_space(
      const pp::topology_point_type_t<BeliefSpace>& b,
      const BeliefSpace& /*unused*/, const StateSpaceOut& /*unused*/) const {
    if constexpr (pp::is_temporal_space_v<BeliefSpace>) {
      return {b.time, b.pt.get_most_likely_state()};
    } else {
      return {b.get_most_likely_state()};
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

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2300014, 1, "maximum_likelihood_map",
                              named_object)
};

}  // namespace ReaK::ctrl

#endif

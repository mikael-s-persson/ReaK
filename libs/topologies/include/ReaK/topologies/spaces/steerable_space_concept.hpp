/**
 * \file steerable_space_concept.hpp
 *
 * This library defines the traits and concepts that pertain to what can be considered
 * a steerable-space, as used in ReaK::pp. Steerable-spaces are based on the Topology concept
 * from the Boost.Graph library, but with additional requirements which are needed
 * in algorithms tailored for path-planning (mostly). Basically, the concept of a
 * steerable-space corresponds to spaces in which non-trivial (steer) motions can
 * be done and their trace be recorded in a steering path or trajectory.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_STEERABLE_SPACE_CONCEPT_HPP
#define REAK_STEERABLE_SPACE_CONCEPT_HPP

#include "ReaK/core/base/defs.hpp"

#include <cmath>
#include <tuple>
#include "boost/concept_check.hpp"

#include "ReaK/topologies/spaces/metric_space_concept.hpp"

/** Main namespace for ReaK */
namespace ReaK::pp {

template <typename SteerableSpace, typename = void>
struct steerable_space_steer_record {
  using type = void;
};

// Get steer_record_type.
template <typename SteerableSpace>
struct steerable_space_steer_record<
    SteerableSpace,
    std::void_t<decltype(std::get<1>(
        std::declval<SteerableSpace>().steer_position_toward(
            std::declval<topology_point_type_t<SteerableSpace>>(), 0.0,
            std::declval<topology_point_type_t<SteerableSpace>>())))>> {
  using type = std::decay_t<decltype(std::get<1>(
      std::declval<std::add_const_t<SteerableSpace>>().steer_position_toward(
          std::declval<topology_point_type_t<SteerableSpace>>(), 0.0,
          std::declval<topology_point_type_t<SteerableSpace>>())))>;
};

template <typename SteerableSpace>
using steerable_space_steer_record_t =
    typename steerable_space_steer_record<SteerableSpace>::type;

/**
 * This concept defines the requirements to fulfill in order to model a steerable-space
 * as used in ReaK::pp. A steerable-space is a special kind of topology in which
 * non-trivial (steered) motions can be done and their trace be recorded in a steering path or trajectory.
 *
 * Valid expressions:
 *
 * tie(p2, st_rec) = space.steer_position_toward(p1, d, p3);  A point (p2) and a steer-record (st_rec) can be obtained
 *from attempting to steer a fraction (d) away from one point (p1) to another (p3).
 *
 * \tparam SteerableSpace The topology type to be checked for this concept.
 */
template <typename SteerableSpace>
struct SteerableSpaceConcept {
  topology_point_type_t<SteerableSpace> p1, p2;
  steerable_space_steer_record_t<SteerableSpace> st_rec;
  SteerableSpace space;
  double d;

  BOOST_CONCEPT_ASSERT((TopologyConcept<SteerableSpace>));

  BOOST_CONCEPT_USAGE(SteerableSpaceConcept) {
    std::tie(p1, st_rec) = space.steer_position_toward(p1, d, p2);
  }
};

template <typename SteerableSpace>
struct is_steerable_space : std::false_type {};

template <typename SteerableSpace>
static constexpr bool is_steerable_space_v =
    is_steerable_space<SteerableSpace>::value;

}  // namespace ReaK::pp

#endif

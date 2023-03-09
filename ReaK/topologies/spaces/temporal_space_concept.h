/**
 * \file temporal_space_concept.h
 *
 * This library defines the traits and concepts related to the construction of a
 * temporal space, i.e. a topology which is constituted of a spatial metric-space
 * and a time topology (also a metric-space). The temporal space is also a metric-space.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2011
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

#ifndef REAK_TOPOLOGIES_SPACES_TEMPORAL_SPACE_CONCEPT_H_
#define REAK_TOPOLOGIES_SPACES_TEMPORAL_SPACE_CONCEPT_H_

#include "ReaK/core/base/defs.h"

#include "boost/concept_check.hpp"

#include "ReaK/topologies/spaces/metric_space_concept.h"

namespace ReaK::pp {

/**
 * This traits class defines the characteristics associated to a temporal-space type.
 * \tparam TemporalSpace The temporal-space type for which the traits are sought.
 */
template <typename TemporalSpace>
struct temporal_space_traits {
  /** The type that describes a point in the space. */
  using point_type = typename TemporalSpace::point_type;
  /** The type that describes a difference between points in the space. */
  using point_difference_type = typename TemporalSpace::point_difference_type;

  /** The topology type which describes the space in which the time values reside. */
  using time_topology = typename TemporalSpace::time_topology;
  /** The topology type which describes the space in which the spatial points reside. */
  using space_topology = typename TemporalSpace::space_topology;
};

/**
 * This concept defines the requirements to fulfill in order to model a temporal-space
 * as used in ReaK::pp. A temporal space is constituted of a spatial topology
 * and a time topology. The temporal space is also a topology.
 *
 * Required concepts:
 *
 * The space-topology should model the TopologyConcept.
 *
 * The time-topology should model the TopologyConcept.
 *
 * The temporal-space should model the TopologyConcept.
 *
 * Valid expressions:
 *
 * s_space = space.get_space_topology();  The spatial-topology (s_space) can be obtained from the temporal-space
 *(space).
 *
 * t_space = space.get_time_topology();  The time-topology (t_space) can be obtained from the temporal-space (space).
 *
 * \tparam TemporalSpace The topology type to be checked for this concept.
 */
template <typename TemporalSpace>
struct TemporalSpaceConcept : public TopologyConcept<TemporalSpace> {
  BOOST_CONCEPT_ASSERT(
      (TopologyConcept<
          typename temporal_space_traits<TemporalSpace>::space_topology>));
  BOOST_CONCEPT_ASSERT(
      (TopologyConcept<
          typename temporal_space_traits<TemporalSpace>::time_topology>));

  BOOST_CONCEPT_USAGE(TemporalSpaceConcept) {
    const typename temporal_space_traits<TemporalSpace>::space_topology&
        cs_space = this->space.get_space_topology();
    RK_UNUSED(cs_space);
    const typename temporal_space_traits<TemporalSpace>::time_topology&
        ct_space = this->space.get_time_topology();
    RK_UNUSED(ct_space);
  }
};

template <typename TemporalSpace>
struct is_temporal_space : std::false_type {};

template <typename TemporalSpace>
static constexpr bool is_temporal_space_v =
    is_temporal_space<TemporalSpace>::value;

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_TEMPORAL_SPACE_CONCEPT_H_

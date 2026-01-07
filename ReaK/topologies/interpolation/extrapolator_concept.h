/**
 * \file extrapolator_concept.h
 *
 * This library defines the traits and concepts related to an extrapolator and extrapolator
 * factories. An extrapolator is simply an object that can generate a trajectory
 * segment from a starting point. The basic scheme used in ReaK is to have one extrapolator
 * factory object that serves to generate fresh extrapolators for a given trajectory segment.
 * In other words, from a user's perspective, the factory can be used to store certain data
 * that all extrapolators will need (e.g. starting point or general parameters of the extrapolation
 * scheme). The intended use from within the library is to store one extrapolator factory
 * object (set by the user) and use that factory object to create extrapolators for each
 * segment of a trajectory.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date August 2012
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_EXTRAPOLATOR_CONCEPT_H_
#define REAK_TOPOLOGIES_INTERPOLATION_EXTRAPOLATOR_CONCEPT_H_


#include "ReaK/topologies/spaces/temporal_space_concept.h"

#include <concepts>

namespace ReaK::pp {

/**
 * This traits class defines the traits that characterize an extrapolator factory
 * class.
 * \tparam ExtrapolatorFactory The trajectory for which the traits are sought.
 */
template <typename ExtrapolatorFactory>
struct extrapolator_factory_traits {
  /** This type describes a point in the temporal space or topology. */
  using point_type = typename ExtrapolatorFactory::point_type;

  /** This type is the extrapolator type that is generated from the factory. */
  using extrapolator_type = typename ExtrapolatorFactory::extrapolator_type;
};

/**
 * This concept class defines the requirements for a class to model an extrapolator
 * concept as used in ReaK::pp. An extrapolator is simply an object that stores some
 * representation of an extrapolated segment from a starting point.
 *
 * Required concepts:
 *
 * The topology should model the TemporalSpace.
 *
 * Valid expressions:
 *
 * extrap.set_start_point(&pt);  The start point of the extrapolated segment can be set (as a const pointers to a
 *point).
 *
 * const point_type* ppt = extrap.get_start_point();  The starting point of the extrapolated segment can be obtained as
 *a const pointer.
 *
 * pt = extrap.get_point_at_time(t);  The point, along the extrapolated segment (extrap), at a given time (t) can be
 *obtained.
 */
template <typename Extrap, typename Space>
concept Extrapolator = TemporalSpace<Space>&& requires(
    Extrap& extrap, const topology_point_type_t<Space>& pt,
    const topology_point_type_t<
        typename temporal_space_traits<Space>::time_topology>& t) {
  extrap.set_start_point(&pt);
  {
    extrap.get_start_point()
    } -> std::convertible_to<const topology_point_type_t<Space>*>;
  {
    extrap.get_point_at_time(t)
    } -> std::convertible_to<topology_point_type_t<Space>>;
};

/**
 * This concept class defines the requirements for a class to model an extrapolator factory
 * concept as used in ReaK::pp. An extrapolator is simply an object that stores some
 * representation of an extrapolated segment from a start point. The basic scheme used
 * in ReaK is to have one extrapolator factory object that serves to generate fresh extrapolators
 * for a given trajectory segment. In other words, from a user's perspective, the factory can
 * be used to store certain data that all extrapolators will need (e.g. information about motion
 * limits or general parameters of the extrapolation scheme). The intended use from within the
 * library is to store one extrapolator factory object (set by the user) and use that factory
 * object to create extrapolators for each segment of a trajectory.
 *
 * Required concepts:
 *
 * The topology should model the TemporalSpace.
 *
 * The extrapolator type should model the Extrapolator on the given topology.
 *
 * Valid expressions:
 *
 * extrap = extrap_fact.create_extrapolator(&pt);  An extrapolator object (extrap) can be created from the factory
 *object (extrap_fact) given a start point as pointers (&pt).
 *
 * extrap_fact.set_temporal_space(pspace);  The temporal space object used by the extrapolators can be set as a const
 *shared-pointer to a topology.
 */
template <typename Factory, typename Space>
concept ExtrapolatorFactoryConcept = TemporalSpace<Space>&& requires(
    Factory& fact, const topology_point_type_t<Space>& pt,
    std::shared_ptr<Space> pspace) {
  fact.set_temporal_space(pspace);
  { fact.create_extrapolator(&pt) } -> Extrapolator<Space>;
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_INTERPOLATION_EXTRAPOLATOR_CONCEPT_H_

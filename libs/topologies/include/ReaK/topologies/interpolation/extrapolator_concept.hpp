/**
 * \file extrapolator_concept.hpp
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

#ifndef REAK_EXTRAPOLATOR_CONCEPT_HPP
#define REAK_EXTRAPOLATOR_CONCEPT_HPP

#include <ReaK/core/base/defs.hpp>

#include <ReaK/topologies/spaces/temporal_space_concept.hpp>

#include <boost/concept_check.hpp>

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
 * The topology should model the TemporalSpaceConcept.
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
 *
 * \tparam Extrapolator The trajectory type for which this concept is checked.
 * \tparam Topology The temporal-topology type on which the trajectory should be able to exist.
 */
template <typename Extrapolator, typename Topology>
struct ExtrapolatorConcept {
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));

  Extrapolator extrap;
  topology_point_type_t<Topology> pt;
  const topology_point_type_t<Topology>* ppt;
  using time_topology = typename temporal_space_traits<Topology>::time_topology;
  topology_point_type_t<time_topology> t;

  BOOST_CONCEPT_USAGE(ExtrapolatorConcept) {
    extrap.set_start_point(&pt);
    ppt = extrap.get_start_point();
    pt = extrap.get_point_at_time(t);
  }
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
 * The topology should model the TemporalSpaceConcept.
 *
 * The extrapolator type should model the ExtrapolatorConcept on the given topology.
 *
 * Valid expressions:
 *
 * extrap = extrap_fact.create_extrapolator(&pt);  An extrapolator object (extrap) can be created from the factory
 *object (extrap_fact) given a start point as pointers (&pt).
 *
 * extrap_fact.set_temporal_space(pspace);  The temporal space object used by the extrapolators can be set as a const
 *shared-pointer to a topology.
 *
 * \tparam ExtrapolatorFactory The extrapolator factory type for which this concept is checked.
 * \tparam Topology The temporal-topology type on which the extrapolated segments should exist.
 */
template <typename ExtrapolatorFactory, typename Topology>
struct ExtrapolatorFactoryConcept {
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT(
      (ExtrapolatorConcept<typename extrapolator_factory_traits<
                               ExtrapolatorFactory>::extrapolator_type,
                           Topology>));

  ExtrapolatorFactory extrap_fact;
  typename extrapolator_factory_traits<ExtrapolatorFactory>::extrapolator_type
      extrap;
  typename extrapolator_factory_traits<ExtrapolatorFactory>::point_type pt;
  std::shared_ptr<Topology> pspace;

  BOOST_CONCEPT_USAGE(ExtrapolatorFactoryConcept) {
    extrap = extrap_fact.create_extrapolator(&pt);
    extrap_fact.set_temporal_space(pspace);
  }
};

}  // namespace ReaK::pp

#endif

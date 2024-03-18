/**
 * \file random_sampler_concept.h
 *
 * This library defines the traits and concepts that pertain to a random-sampler over
 * a topology, as used in ReaK::pp. Random-samplers are used in ReaK to generate
 * points within a topology, often with the hopes of connecting them into a planned path
 * or to propagate or evaluate them in estimation or optimization methods.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2012
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

#ifndef REAK_TOPOLOGIES_SPACES_RANDOM_SAMPLER_CONCEPT_H_
#define REAK_TOPOLOGIES_SPACES_RANDOM_SAMPLER_CONCEPT_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"

#include <concepts>

/** Main namespace for ReaK */
namespace ReaK::pp {

/**
 * This concept defines the requirements to fulfill in order to model a random-sampler
 * as used in ReaK::pp. A random-sampler is essentially a callable type that can produce
 * a point within a given topology.
 *
 * Valid expressions:
 *
 * p = sampler(space);  A random point (p) can be obtained by calling the random-sampler (sampler) by
 *providing a const-ref to the topology (space).
 */
template <typename Sampler, typename Space>
concept RandomSampler = Topology<Space> && requires (const Space& space, Sampler& sampler) {
  { sampler(space) } -> std::convertible_to<topology_point_type_t<Space>>;
};

/**
 * This tag-type is used to identify (during a "get" call) that the random-sampler object is
 * to be fetched, if one is associated to a given topology (usually the default one).
 */
enum random_sampler_t { random_sampler };

/**
 * This traits class defines the types and constants associated to a point distribution.
 * \tparam PointDistribution The topology type for which the point distribution traits are sought.
 */
template <typename PointDistribution>
struct point_distribution_traits {
  /** The type that describes the random-sampler type for the distribution. */
  using random_sampler_type = typename PointDistribution::random_sampler_type;
};

template <typename PointDistribution>
struct point_distribution_random_sampler {
  using type = typename point_distribution_traits<
      PointDistribution>::random_sampler_type;
};
template <typename PointDistribution>
using point_distribution_random_sampler_t =
    typename point_distribution_random_sampler<PointDistribution>::type;

/**
 * This concept defines the requirements to fulfill in order to model a point distribution
 * as used in ReaK::pp. A point distribution is a special kind of topology which has a
 * random-sampler over it.
 *
 * Required concepts:
 *
 * The Space should model Topology.
 *
 * The random-sampler should model RandomSampler.
 *
 * Valid expressions:
 *
 * rand_sampler = get(random_sampler,space);  The random-sampler can be obtained by a tagged "get" call on the
 *point-distribution.
 */
template <typename Space>
concept PointDistribution = Topology<Space> && RandomSampler<point_distribution_random_sampler_t<Space>, Space> &&
  requires (const Space& space) {
    { get(random_sampler, space) } -> std::convertible_to<point_distribution_random_sampler_t<Space>>;
  };

template <typename T>
static constexpr bool is_point_distribution_v = PointDistribution<T>;

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_RANDOM_SAMPLER_CONCEPT_H_

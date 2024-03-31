/**
 * \file generic_sampler_factory.h
 *
 * This library provides an implementation of a generic sampler factory.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2012
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

#ifndef REAK_TOPOLOGIES_SPACES_GENERIC_SAMPLER_FACTORY_H_
#define REAK_TOPOLOGIES_SPACES_GENERIC_SAMPLER_FACTORY_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/math/lin_alg/arithmetic_tuple.h"
#include "ReaK/topologies/spaces/metric_space_tuple.h"
#include "ReaK/topologies/spaces/tangent_bundle_concept.h"
#include "ReaK/topologies/spaces/temporal_space_concept.h"

#include <cmath>
#include <utility>

namespace ReaK::pp {

namespace detail {
namespace {

template <typename Sampler, typename PointType, typename SpaceType,
          typename Factory>
void generate_sample_impl(const Sampler& sampler, PointType& result,
                          const SpaceType& space, const Factory& factory) {
  result = sampler(space, factory);
}

template <typename Sampler, typename PointType, typename SpaceType>
void generate_sample_impl(const Sampler& sampler, PointType& result,
                          const SpaceType& space) {
  result = sampler(space);
}

template <typename Sampler, typename PointType, typename SpaceTuple,
          typename TupleDistMetric, typename Factory>
void generate_sample_impl(
    const Sampler& sampler, PointType& result,
    const metric_space_tuple<SpaceTuple, TupleDistMetric>& space,
    const Factory& factory);

template <typename Sampler, typename PointType, typename SpaceTuple,
          typename TupleDistMetric>
void generate_sample_impl(
    const Sampler& sampler, PointType& result,
    const metric_space_tuple<SpaceTuple, TupleDistMetric>& space);

template <typename Sampler, typename Factory = void>
struct tuple_sample_generator {
  Sampler* p_sampler;
  const Factory* p_factory;
  tuple_sample_generator(Sampler& sampler, const Factory& factory)
      : p_sampler(&sampler), p_factory(&factory) {}
  template <typename PointType, typename SpaceType>
  int operator()(PointType& result, const SpaceType& space) const {
    generate_sample_impl(*p_sampler, result, space, *p_factory);
    return int{};
  }
};

template <typename Sampler>
struct tuple_sample_generator<Sampler, void> {
  Sampler* p_sampler;
  explicit tuple_sample_generator(Sampler& sampler) : p_sampler(&sampler) {}
  template <typename PointType, typename SpaceType>
  int operator()(PointType& result, const SpaceType& space) const {
    generate_sample_impl(*p_sampler, result, space);
    return int{};
  }
};

template <typename Sampler, typename PointType, typename SpaceTuple,
          typename TupleDistMetric, typename Factory>
void generate_sample_impl(
    const Sampler& sampler, PointType& result,
    const metric_space_tuple<SpaceTuple, TupleDistMetric>& space,
    const Factory& factory) {
  arithmetic_tuple_details::tuple_for_each(
      result, space,
      tuple_sample_generator<Sampler, Factory>(sampler, factory));
}

template <typename Sampler, typename PointType, typename SpaceTuple,
          typename TupleDistMetric>
void generate_sample_impl(
    const Sampler& sampler, PointType& result,
    const metric_space_tuple<SpaceTuple, TupleDistMetric>& space) {
  arithmetic_tuple_details::tuple_for_each(
      result, space, tuple_sample_generator<Sampler>(sampler));
}

}  // namespace
}  // namespace detail

/**
 * This functor class implements a generic sampler in a topology.
 * \tparam Factory The sampler factory type which created this generic sampler (this is a fly-weight class to store
 * whatever global parameters are needed by the sampler implementation).
 * \tparam Sampler The sampler implementation which is used to perform the specifics of the actual sampling.
 */
template <typename Sampler, typename SpaceType, typename Factory = void>
class generic_sampler {
 public:
  using self = generic_sampler<Sampler, SpaceType, Factory>;
  using point_type = typename topology_traits<SpaceType>::point_type;
  using topology = SpaceType;

 private:
  const Factory* parent;
  Sampler sampler;

 public:
  /**
   * Default constructor.
   */
  explicit generic_sampler(Sampler aSampler = Sampler(),
                           const Factory* aParent = nullptr)
      : parent(aParent), sampler(aSampler) {}

  template <typename OtherFactory>
  point_type operator()(const topology& space,
                        const OtherFactory& aNewParent) const {
    point_type result;
    detail::generate_sample_impl(sampler, result, space, aNewParent);
    return result;
  }

  point_type operator()(const topology& space) const {
    point_type result;
    detail::generate_sample_impl(sampler, result, space, *parent);
    return result;
  }
};

template <typename Sampler, typename SpaceType>
class generic_sampler<Sampler, SpaceType, void> {
 public:
  using self = generic_sampler<Sampler, SpaceType, void>;
  using point_type = typename topology_traits<SpaceType>::point_type;
  using topology = SpaceType;

 private:
  Sampler sampler;

 public:
  /**
   * Default constructor.
   */
  explicit generic_sampler(Sampler aSampler = Sampler())
      : sampler(std::move(aSampler)) {}

  template <typename OtherFactory>
  point_type operator()(const topology& space,
                        const OtherFactory& aParent) const {
    point_type result;
    detail::generate_sample_impl(sampler, result, space, aParent);
    return result;
  }

  point_type operator()(const topology& space) const {
    point_type result;
    detail::generate_sample_impl(sampler, result, space);
    return result;
  }
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_GENERIC_SAMPLER_FACTORY_H_

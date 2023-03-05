/**
 * \file tangent_bundle_concept.hpp
 *
 * This library defines the meta-functions, traits and concepts related to the construction of a
 * tangent bundle, i.e. a topology which bundles a sequence of derivative topologies can be
 * obtained given a topology against which the derivative is taken (e.g. a time-topology).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date September 2011
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

#ifndef REAK_TANGENT_BUNDLE_CONCEPT_HPP
#define REAK_TANGENT_BUNDLE_CONCEPT_HPP

#include "ReaK/core/base/defs.hpp"

#include "boost/concept_check.hpp"

#include "ReaK/topologies/spaces/metric_space_concept.hpp"

#include <type_traits>

namespace ReaK::pp {

/* Just a prototype. */
template <int Idx, typename SpaceType, typename IndependentSpace>
const SpaceType& get_space(const SpaceType& s,
                           const IndependentSpace& /*unused*/) {
  return s;
}

/* Just a prototype. */
template <int Idx, typename SpaceType, typename IndependentSpace>
SpaceType& get_space(SpaceType& s, const IndependentSpace& /*unused*/) {
  return s;
}

/* Just a prototype. */
template <int Idx, typename PointDiffType, typename TimeDiffType,
          typename SpaceType, typename IndependentSpace>
PointDiffType lift_to_space(const PointDiffType& dp, const TimeDiffType& dt,
                            const SpaceType& /*unused*/,
                            const IndependentSpace& /*unused*/) {
  return dp / dt;
}

/* Just a prototype. */
template <int Idx, typename PointType, typename TimeDiffType,
          typename SpaceType, typename IndependentSpace>
PointType descend_to_space(const PointType& v, const TimeDiffType& dt,
                           const SpaceType& /*unused*/,
                           const IndependentSpace& /*unused*/) {
  return v * dt;
}

/**
 * This meta-function provides an integral-constant type with the maximum differential
 * order that can be provided by the tangent bundle against a given independent space.
 * \tparam TangentBundle The tangent bundle whose maximum differential order is sought.
 * \tparam IndependentSpace The independent space against which the derivation is taken.
 */
template <typename TangentBundle, typename IndependentSpace>
struct max_derivation_order : std::integral_constant<int, 0> {};

template <typename TangentBundle, typename IndependentSpace>
static constexpr int max_derivation_order_v =
    max_derivation_order<TangentBundle, IndependentSpace>::value;

/**
 * This meta-function provides the type of N-order differential space with a tangent bundle.
 * \tparam TangentBundle The tangent bundle for which the N-order differential space is sought.
 * \tparam IndependentSpace The independent space against which the differentiation is applied.
 * \tparam Order The order of differentiation of the differential space type that is sought.
 */
template <typename TangentBundle, typename IndependentSpace, std::size_t Order>
struct derived_N_order_space {
  using type = TangentBundle;
};

template <typename TangentBundle, typename IndependentSpace, std::size_t Order>
using derived_N_order_space_t =
    typename derived_N_order_space<TangentBundle, IndependentSpace,
                                   Order>::type;

/**
 * This concept defines the requirements to fulfill in order to model a differential relation
 * as used in ReaK::pp. A differentiable relation serves to map a spatial topology
 * its derivative topology with respect to a given independent space (e.g. time). The mapping
 * also serves to map elements between to two spaces (tangent lifts and descents).
 *
 * Required concepts:
 *
 * All topologies involved should model the TopologyConcept.
 *
 * Valid expressions:
 *
 * space = get_space<0..N>(diff_space,t_space);  The metric space (space) corresponding to the 0 to Nth order derivative
 *space can be obtained given an independent space (e.g. time topology, t_space).
 *
 * v = lift_to_space<1..N>(dp,dt,diff_space,t_space);  A derivative-point (v) can be obtained from lifting a
 *point-difference (dp) from the space via a difference-point on the independent space (dt). This expression is
 *analogous to v = dp / dt.
 *
 * dp = descend_to_space<0..N-1>(v,dt,diff_space,t_space);  A point-difference (dp) can be obtained from descending a
 *derivative-point (v) to the space via a difference-point on the independent space (dt). This expression is analogous
 *to dp = v * dt.
 *
 * \tparam TangentBundle The topology type to be checked for this concept.
 * \tparam Order The maximum order of differentiation of the tangent bundle type.
 * \tparam IndependentSpace The topology type to be checked for this concept.
 */
template <typename TangentBundle, std::size_t Order, typename IndependentSpace>
struct TangentBundleConcept
    : TangentBundleConcept<TangentBundle, Order - 1, IndependentSpace> {

  static_assert(max_derivation_order_v<TangentBundle, IndependentSpace> >=
                Order);

  using base_space_type =
      derived_N_order_space_t<TangentBundle, IndependentSpace, Order - 1>;
  using derived_space_type =
      derived_N_order_space_t<TangentBundle, IndependentSpace, Order>;

  BOOST_CONCEPT_ASSERT((TopologyConcept<derived_space_type>));

  topology_point_difference_type_t<base_space_type> dp;

  topology_point_type_t<derived_space_type> v;

  topology_point_difference_type_t<IndependentSpace> dt;

  BOOST_CONCEPT_USAGE(TangentBundleConcept) {
    const derived_space_type& space =
        get_space<Order>(this->diff_space, this->t_space);
    RK_UNUSED(space);
    v = lift_to_space<Order>(dp, dt, this->diff_space, this->t_space);
    dp = descend_to_space<Order - 1>(v, dt, this->diff_space, this->t_space);
  }
};

template <typename TangentBundle, typename IndependentSpace>
struct TangentBundleConcept<TangentBundle, 0, IndependentSpace> {

  using base_space_type =
      derived_N_order_space_t<TangentBundle, IndependentSpace, 0>;

  BOOST_CONCEPT_ASSERT((TopologyConcept<base_space_type>));
  BOOST_CONCEPT_ASSERT((TopologyConcept<IndependentSpace>));

  TangentBundle diff_space;
  IndependentSpace t_space;

  BOOST_CONCEPT_USAGE(TangentBundleConcept) {
    const base_space_type& space =
        get_space<0>(this->diff_space, this->t_space);
    RK_UNUSED(space);
  }
};

}  // namespace ReaK::pp

#endif

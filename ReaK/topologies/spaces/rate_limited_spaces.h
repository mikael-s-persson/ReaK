/**
 * \file rate_limited_spaces.h
 *
 * This library provides classes that define a number of reach-time-space class templates. A reach-time-space is
 * a transformation on a tangent bundle such that points in the spaces and their differences
 * represent the time-to-reach (or reach-time) based on the N-order rate limitations on the tangent spaces.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2011
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

#ifndef REAK_TOPOLOGIES_SPACES_RATE_LIMITED_SPACES_H_
#define REAK_TOPOLOGIES_SPACES_RATE_LIMITED_SPACES_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/core/serialization/serializable.h"

#include "ReaK/math/lin_alg/arithmetic_tuple.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/prob_distribution_concept.h"

#include "ReaK/topologies/spaces/default_random_sampler.h"
#include "ReaK/topologies/spaces/differentiable_space.h"
#include "ReaK/topologies/spaces/tuple_distance_metrics.h"

#include "ReaK/topologies/spaces/rate_limited_space_metamaps.h"

#include <type_traits>

namespace ReaK::pp {

/**
 * This class defines the differentiation rule to apply either to lift a
 * point-difference (e.g. finite-difference) to the tangent space, or to descend
 * a tangent vector to a point-difference, for topologies whose point-difference
 * vectors are expressed as reach-time values.
 */
struct reach_time_differentiation : public serializable {

  double max_rate_reach_time;

  explicit reach_time_differentiation(double aMaxRateReachTime = 1.0)
      : max_rate_reach_time(aMaxRateReachTime) {}

  /**
   * This function will lift a point-difference vector into its corresponding tangent vector.
   * This function performs a simple division, dp * (max_rate_reach_time / dt).
   * \tparam T The destination type, a point in the tangent space.
   * \tparam U The source type, a point-difference in the base space.
   * \tparam V A type representing the independent variable's difference (e.g. time-difference).
   * \tparam TSpace The type of the independent space (e.g. time-space).
   * \param v The resulting point in the tangent space.
   * \param dp The point-difference that is being lifted.
   * \param dt The time-difference value (i.e. the difference in the independent variable).
   */
  template <typename T, typename U, typename V, typename TSpace>
  void lift(T& v, const U& dp, const V& dt, const TSpace& /*unused*/) const {
    v = dp * (max_rate_reach_time / dt);
  }
  /**
   * This function will descend a tangent vector into its corresponding point-difference vector.
   * This function performs a simple multiplication, v * (dt / max_rate_reach_time).
   * \tparam T The destination type, a point-difference in the base space.
   * \tparam U The source type, a point in the tangent space.
   * \tparam V A type representing the independent variable's difference (e.g. time-difference).
   * \tparam TSpace The type of the independent space (e.g. time-space).
   * \param dp The resulting point-difference in the base space.
   * \param v The point in the tangent space that is being descended.
   * \param dt The time-difference value (i.e. the difference in the independent variable).
   */
  template <typename T, typename U, typename V, typename TSpace>
  void descend(T& dp, const U& v, const V& dt, const TSpace& /*unused*/) const {
    dp = v * (dt / max_rate_reach_time);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(max_rate_reach_time);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(max_rate_reach_time);
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(reach_time_differentiation, 0xC2420001, 1,
                              "reach_time_differentiation", serializable)
};

namespace detail {

template <int Idx, typename ForwardIter, typename ReachTimeDiffTuple>
void assign_reach_time_diff_rules(ForwardIter it,
                                  ReachTimeDiffTuple& diff_rule) {
  if constexpr (Idx < arithmetic_tuple_size_v<ReachTimeDiffTuple>) {
    get<Idx>(diff_rule) = reach_time_differentiation{*it};
    assign_reach_time_diff_rules<Idx + 1>(++it, diff_rule);
  }
}

template <typename ReachTimeDiffTuple, typename ForwardIter>
ReachTimeDiffTuple construct_reach_time_diff_rules(ForwardIter first,
                                                   ForwardIter last) {
  ReachTimeDiffTuple result;
  assign_reach_time_diff_rules<0>(first, result);
  return result;
}

}  // namespace detail

/**
 * This class template can be used to glue together a number of spaces by a reach-time differentiation
 * relationship, where each differentiation / integration operation (or more formally speaking, each
 * lift and descent through the tangent bundle) is governed by its own reach-time differentiation rule.
 * This class template models Topology (if all underlying spaces do as well), and models
 * TangentBundle for as high an order as there are differentiation rules and spaces to support.
 *
 * \tparam IndependentSpace The type of the independent-space against which the differentiation is
 *                          taken (e.g. time_topology). There are no formal requirements on this type,
 *                          it is merely used as a placeholder by this class (although the differentiation
 *                          rules might require more of this type).
 * \tparam SpaceTuple A tuple type (e.g. arithmetic_tuple) which provides a set of spaces that are arranged
 *                    in sequence of differentiation levels (e.g. space 0 -- diff --> space 1 -- diff --> space 2 ...).
 * \tparam TupleDistanceMetric A distance metric type which models the DistanceMetric and operates on a
 *                             space-tuple (e.g. arithmetic_tuple).
 */
template <Topology IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric = manhattan_tuple_distance>
class reach_time_diff_space
    : public differentiable_space<IndependentSpace, SpaceTuple,
                                  TupleDistanceMetric,
                                  reach_time_differentiation> {
 public:
  using self =
      reach_time_diff_space<IndependentSpace, SpaceTuple, TupleDistanceMetric>;
  using base_type =
      differentiable_space<IndependentSpace, SpaceTuple, TupleDistanceMetric,
                           reach_time_differentiation>;
  using diff_rule_tuple = typename base_type::diff_rule_tuple;

  static constexpr std::size_t dimensions = base_type::dimensions;

  /**
   * Parametrized and default constructor.
   * \param aSpaces The space tuple to initialize the spaces with.
   * \param aDist The distance metric functor on the space-tuple.
   * \param aDiffRules The differentiation rule tuple to initialize the diff-rule functors with.
   */
  explicit reach_time_diff_space(
      const SpaceTuple& aSpaces = SpaceTuple{},
      const TupleDistanceMetric& aDist = TupleDistanceMetric{},
      const diff_rule_tuple& aDiffRules = diff_rule_tuple{})
      : base_type(aSpaces, aDist, aDiffRules) {}

  /**
   * Parametrized constructor which creates the reach-time differentiation rules from a list
   * of rate-limits on spaces of the tangent bundle (element 0 should be the velocity limit,
   * then the acceleration limit, and so on, so forth).
   * \param first An iterator to the first rate-limit used for differentiation rules.
   * \param last An iterator to the one-past-last rate-limit used for differentiation rules.
   * \param aSpaces The space tuple to initialize the spaces with.
   * \param aDist The distance metric functor on the space-tuple.
   */
  template <typename ForwardIter>
  reach_time_diff_space(
      ForwardIter first, ForwardIter last,
      const SpaceTuple& aSpaces = SpaceTuple{},
      const TupleDistanceMetric& aDist = TupleDistanceMetric{})
      : base_type(aSpaces, aDist,
                  detail::construct_reach_time_diff_rules<diff_rule_tuple>(
                      first, last)) {}

  /**
   * Parametrized constructor which creates the reach-time differentiation rules from a list
   * of rate-limits on spaces of the tangent bundle (element 0 should be the velocity limit,
   * then the acceleration limit, and so on, so forth).
   * \param aList An initializer-list to the rate-limits used for differentiation rules.
   * \param aSpaces The space tuple to initialize the spaces with.
   * \param aDist The distance metric functor on the space-tuple.
   */
  reach_time_diff_space(
      std::initializer_list<double> aList,
      const SpaceTuple& aSpaces = SpaceTuple{},
      const TupleDistanceMetric& aDist = TupleDistanceMetric{})
      : base_type(aSpaces, aDist,
                  detail::construct_reach_time_diff_rules<diff_rule_tuple>(
                      aList.begin(), aList.end())) {}

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_type::save(A, base_type::getStaticObjectType()->TypeVersion());
  }
  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base_type::load(A, base_type::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400010, 1, "reach_time_diff_space",
                              base_type)
};

template <typename IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric>
struct is_metric_symmetric<
    reach_time_diff_space<IndependentSpace, SpaceTuple, TupleDistanceMetric>>
    : is_metric_symmetric<typename reach_time_diff_space<
          IndependentSpace, SpaceTuple, TupleDistanceMetric>::base_type> {};

template <typename SpaceTuple, typename TupleDistanceMetric,
          typename IndependentSpace, typename IndependentSpace2>
struct max_derivation_order<
    reach_time_diff_space<IndependentSpace, SpaceTuple, TupleDistanceMetric>,
    IndependentSpace2>
    : max_derivation_order<
          typename reach_time_diff_space<IndependentSpace, SpaceTuple,
                                         TupleDistanceMetric>::base_type,
          IndependentSpace2> {};

template <typename SpaceTuple, typename TupleDistanceMetric,
          typename IndependentSpace, typename IndependentSpace2,
          std::size_t Order>
struct derived_N_order_space<
    reach_time_diff_space<IndependentSpace, SpaceTuple, TupleDistanceMetric>,
    IndependentSpace2, Order> {
  using type = arithmetic_tuple_element_t<Order, SpaceTuple>;
};

/**
 * This function returns the space at a given differential order against a given independent-space.
 * \tparam Idx The differential order (e.g. 0: position, 1: velocity, 2: acceleration).
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric, typename IndependentSpace2>
const auto& get_space(const reach_time_diff_space<IndependentSpace, SpaceTuple,
                                                  TupleDistanceMetric>& s,
                      const IndependentSpace2& t) {
  return get_space<Idx>(
      static_cast<const differentiable_space<IndependentSpace, SpaceTuple,
                                             TupleDistanceMetric,
                                             reach_time_differentiation>&>(s),
      t);
}

/**
 * This function returns the space at a given differential order against a given independent-space.
 * \tparam Idx The differential order (e.g. 0: position, 1: velocity, 2: acceleration).
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric, typename IndependentSpace2>
auto& get_space(
    reach_time_diff_space<IndependentSpace, SpaceTuple, TupleDistanceMetric>& s,
    const IndependentSpace2& t) {
  return get_space<Idx>(
      static_cast<differentiable_space<IndependentSpace, SpaceTuple,
                                       TupleDistanceMetric,
                                       reach_time_differentiation>&>(s),
      t);
}

/**
 * This function returns the differentiation functor at a given order against a given independent-space.
 * \tparam Idx The differential order (e.g. 0: position/velocity, 1: velocity/acceleration).
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric, typename IndependentSpace2>
const reach_time_differentiation& get_diff_rule(
    const reach_time_diff_space<IndependentSpace, SpaceTuple,
                                TupleDistanceMetric>& s,
    const IndependentSpace2& /*unused*/) {
  return s.template get_diff_rule_impl<Idx>();
}

/**
 * This function returns the differentiation functor at a given order against a given independent-space.
 * \tparam Idx The differential order (e.g. 0: position/velocity, 1: velocity/acceleration).
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric, typename IndependentSpace2>
reach_time_differentiation& get_diff_rule(
    reach_time_diff_space<IndependentSpace, SpaceTuple, TupleDistanceMetric>& s,
    const IndependentSpace2& /*unused*/) {
  return s.template get_diff_rule_impl<Idx>();
}

/**
 * This function lifts a point-difference in space Idx-1 into a point in space Idx.
 * \tparam Idx The differential order of the destination space.
 * \param dp The point-difference in the space Idx-1.
 * \param dt The point-difference in the independent-space (e.g. time).
 * \param space The differentiable-space.
 * \param t_space The independent-space.
 * \return The point in space Idx which is the tangential lift of the point-difference in space Idx-1.
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric, typename IndependentSpace2>
arithmetic_tuple_element_t<
    Idx, typename reach_time_diff_space<IndependentSpace, SpaceTuple,
                                        TupleDistanceMetric>::point_type>
lift_to_space(const arithmetic_tuple_element_t<
                  Idx - 1, typename reach_time_diff_space<
                               IndependentSpace, SpaceTuple,
                               TupleDistanceMetric>::point_difference_type>& dp,
              const topology_point_difference_type_t<IndependentSpace2>& dt,
              const reach_time_diff_space<IndependentSpace, SpaceTuple,
                                          TupleDistanceMetric>& space,
              const IndependentSpace2& t_space) {
  return space.template lift_to_space<Idx>(dp, dt, t_space);
}

/**
 * This function descends a point in space Idx+1 into a point-difference in space Idx.
 * \tparam Idx The differential order of the destination space.
 * \param v The point in the space Idx+1.
 * \param dt The point-difference in the independent-space (e.g. time).
 * \param space The differentiable-space.
 * \param t_space The independent-space.
 * \return The point-difference in space Idx which is the tangential descent of the point in space Idx+1.
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric, typename IndependentSpace2>
arithmetic_tuple_element_t<Idx, typename reach_time_diff_space<
                                    IndependentSpace, SpaceTuple,
                                    TupleDistanceMetric>::point_difference_type>
descend_to_space(
    const arithmetic_tuple_element_t<
        Idx + 1,
        typename reach_time_diff_space<IndependentSpace, SpaceTuple,
                                       TupleDistanceMetric>::point_type>& v,
    const topology_point_difference_type_t<IndependentSpace2>& dt,
    const reach_time_diff_space<IndependentSpace, SpaceTuple,
                                TupleDistanceMetric>& space,
    const IndependentSpace2& t_space) {
  return space.template descend_to_space<Idx>(v, dt, t_space);
}

namespace detail {

template <std::size_t Size, typename SpaceTuple>
struct get_rate_illimited_space_tuple_impl {
  // static_assert(false);
};

template <std::size_t Size, typename... Spaces>
struct get_rate_illimited_space_tuple_impl<Size, std::tuple<Spaces...>> {
  using type = arithmetic_tuple<get_rate_illimited_space_t<Spaces>...>;
};

template <std::size_t Size, typename... Spaces>
struct get_rate_illimited_space_tuple_impl<Size, arithmetic_tuple<Spaces...>> {
  using type = arithmetic_tuple<get_rate_illimited_space_t<Spaces>...>;
};

template <typename SpaceTuple>
struct get_rate_illimited_space_tuple
    : get_rate_illimited_space_tuple_impl<arithmetic_tuple_size_v<SpaceTuple>,
                                          SpaceTuple> {};

template <typename SpaceTuple>
using get_rate_illimited_space_tuple_t =
    typename get_rate_illimited_space_tuple<SpaceTuple>::type;

template <std::size_t Size, typename SpaceTuple>
struct get_rate_limited_space_tuple_impl {
  // static_assert(false);
};

template <std::size_t Size, typename... Spaces>
struct get_rate_limited_space_tuple_impl<Size, std::tuple<Spaces...>> {
  using type = arithmetic_tuple<get_rate_limited_space_t<Spaces>...>;
};

template <std::size_t Size, typename... Spaces>
struct get_rate_limited_space_tuple_impl<Size, arithmetic_tuple<Spaces...>> {
  using type = arithmetic_tuple<get_rate_limited_space_t<Spaces>...>;
};

template <typename SpaceTuple>
struct get_rate_limited_space_tuple
    : get_rate_limited_space_tuple_impl<arithmetic_tuple_size_v<SpaceTuple>,
                                        SpaceTuple> {};

template <typename SpaceTuple>
using get_rate_limited_space_tuple_t =
    typename get_rate_limited_space_tuple<SpaceTuple>::type;

}  // namespace detail

template <typename IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric>
struct get_rate_illimited_space<
    reach_time_diff_space<IndependentSpace, SpaceTuple, TupleDistanceMetric>> {
  using type =
      differentiable_space<IndependentSpace,
                           detail::get_rate_illimited_space_tuple_t<SpaceTuple>,
                           TupleDistanceMetric>;
};

template <typename IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric>
struct get_rate_limited_space<
    differentiable_space<IndependentSpace, SpaceTuple, TupleDistanceMetric>> {
  using type =
      reach_time_diff_space<IndependentSpace,
                            detail::get_rate_limited_space_tuple_t<SpaceTuple>,
                            TupleDistanceMetric>;
};

template <typename SpaceTuple, typename TupleDistanceMetric>
struct get_rate_illimited_space<
    metric_space_tuple<SpaceTuple, TupleDistanceMetric>> {
  using type =
      metric_space_tuple<detail::get_rate_illimited_space_tuple_t<SpaceTuple>,
                         TupleDistanceMetric>;
};

template <typename SpaceTuple, typename TupleDistanceMetric>
struct get_rate_limited_space<
    metric_space_tuple<SpaceTuple, TupleDistanceMetric>> {
  using type =
      metric_space_tuple<detail::get_rate_limited_space_tuple_t<SpaceTuple>,
                         TupleDistanceMetric>;
};

}  // namespace ReaK::pp

namespace ReaK {

/* Specialization, see general template docs. */
template <typename IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric>
struct arithmetic_tuple_size<pp::reach_time_diff_space<
    IndependentSpace, SpaceTuple, TupleDistanceMetric>>
    : arithmetic_tuple_size<SpaceTuple> {};

/* Specialization, see general template docs. */
template <int Idx, typename IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric>
struct arithmetic_tuple_element<
    Idx, pp::reach_time_diff_space<IndependentSpace, SpaceTuple,
                                   TupleDistanceMetric>> {
  using type = arithmetic_tuple_element_t<Idx, SpaceTuple>;
};

/* Specialization, see general template docs. */
template <int Idx, typename IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric>
struct arithmetic_tuple_element<
    Idx, const pp::reach_time_diff_space<IndependentSpace, SpaceTuple,
                                         TupleDistanceMetric>> {
  using type = arithmetic_tuple_element_t<Idx, const SpaceTuple>;
};

/* Specialization, see general template docs. */
template <int Idx, typename IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric>
struct arithmetic_tuple_element<
    Idx, volatile pp::reach_time_diff_space<IndependentSpace, SpaceTuple,
                                            TupleDistanceMetric>> {
  using type = arithmetic_tuple_element_t<Idx, volatile SpaceTuple>;
};

/* Specialization, see general template docs. */
template <int Idx, typename IndependentSpace, typename SpaceTuple,
          typename TupleDistanceMetric>
struct arithmetic_tuple_element<
    Idx, const volatile pp::reach_time_diff_space<IndependentSpace, SpaceTuple,
                                                  TupleDistanceMetric>> {
  using type = arithmetic_tuple_element_t<Idx, const volatile SpaceTuple>;
};

}  // namespace ReaK

#endif  // REAK_TOPOLOGIES_SPACES_RATE_LIMITED_SPACES_H_

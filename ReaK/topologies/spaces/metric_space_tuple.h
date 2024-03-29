/**
 * \file metric_space_tuple.h
 *
 * This library provides classes that define a metric-space tuple class template. A metric-space tuple is
 * a simple association of several topologies (metric-spaces) which, in turn, also models a metric-space
 * (conditional upon each underlying spaces being a metric-space as well).
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

#ifndef REAK_TOPOLOGIES_SPACES_METRIC_SPACE_TUPLE_H_
#define REAK_TOPOLOGIES_SPACES_METRIC_SPACE_TUPLE_H_

#include <utility>
#include "ReaK/core/base/defs.h"
#include "ReaK/core/base/shared_object.h"
#include "ReaK/core/serialization/serializable.h"
#include "ReaK/math/lin_alg/arithmetic_tuple.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/metric_space_tuple_fwd.h"

#include "ReaK/topologies/spaces/default_random_sampler.h"
#include "ReaK/topologies/spaces/tuple_distance_metrics.h"

namespace ReaK::pp {

namespace metric_space_tuple_details {

template <typename SpaceTuple>
struct topology_traits_tuple_impl {
  // static_assert(false);
};

template <typename... Spaces>
struct topology_traits_tuple_impl<arithmetic_tuple<Spaces...>> {
  using point_type = arithmetic_tuple<topology_point_type_t<Spaces>...>;
  using point_difference_type =
      arithmetic_tuple<topology_point_difference_type_t<Spaces>...>;

  static constexpr std::size_t dimensions = 0;  // TODO
};

template <typename SpaceTuple>
struct topology_traits_tuple : topology_traits_tuple_impl<SpaceTuple> {};

struct mst_random_point_computer {
  template <typename Space>
  auto operator()(const Space& s) const {
    return get(random_sampler, s)(s);
  }
};

struct mst_difference_computer {
  template <typename Space, typename Point>
  auto operator()(const Space& s, const Point& p1, const Point& p2) const {
    return s.difference(p1, p2);
  }
};

struct mst_move_position_toward_computer {
  double d;
  explicit mst_move_position_toward_computer(double aD) : d(aD) {}
  template <typename Space, typename Point>
  auto operator()(const Space& s, const Point& p1, const Point& p2) const {
    return s.move_position_toward(p1, d, p2);
  }
};

struct mst_move_position_back_to_computer {
  double d;
  explicit mst_move_position_back_to_computer(double aD) : d(aD) {}
  template <typename Space, typename Point>
  auto operator()(const Space& s, const Point& p1, const Point& p2) const {
    return s.move_position_back_to(p1, d, p2);
  }
};

struct mst_origin_computer {
  template <typename Space>
  auto operator()(const Space& s) const {
    return s.origin();
  }
};

struct mst_adjust_computer {
  template <typename Space, typename Point, typename PointDiff>
  auto operator()(const Space& s, const Point& p, const PointDiff& dp) const {
    return s.adjust(p, dp);
  }
};

struct mst_bring_point_in_bounds_computer {
  template <typename Space, typename Point>
  int operator()(const Space& s, Point& p) const {
    s.bring_point_in_bounds(p);
    return int{};
  }
};

struct mst_get_diff_to_boundary_computer {
  template <typename Space, typename Point>
  auto operator()(const Space& s, const Point& p) const {
    return s.get_diff_to_boundary(p);
  }
};

struct mst_is_in_bounds_computer {
  bool* p_result;
  explicit mst_is_in_bounds_computer(bool& result) : p_result(&result) {}
  template <typename Space, typename Point>
  int operator()(const Space& s, const Point& p) const {
    (*p_result) = ((*p_result) && s.is_in_bounds(p));
    return int{};
  }
};

}  // namespace metric_space_tuple_details

/**
 * This class template can be used to glue together a number of spaces into a tuple. Depending on the models
 * supported by the underlying spaces included in the tuple, this class template models
 * the TopologyConcept, the MetricSpaceConcept, the PointDistributionConcept, the LieGroupConcept, and
 * the BoundedSpaceConcept (i.e. the metric-space tuple class will model all the concepts which are also
 * modeled by all the spaces it includes). This class is also a tuple class (meaning
 * that the std::tuple meta-functions, as well as the ReaK.Arithmetic-tuple meta-functions
 * will work on this class, with the usual semantics).
 *
 * \tparam SpaceTuple A tuple type (e.g. arithmetic_tuple) which provides a set of spaces to glue together.
 * \tparam TupleDistanceMetric A distance metric type which models the DistanceMetricConcept and operates on a
 *space-tuple (e.g. arithmetic_tuple).
 */
template <typename SpaceTuple, typename TupleDistanceMetric>
class metric_space_tuple : public shared_object {
 protected:
  SpaceTuple m_spaces;
  TupleDistanceMetric m_dist;

 public:
  using self = metric_space_tuple<SpaceTuple, TupleDistanceMetric>;
  using self_traits =
      metric_space_tuple_details::topology_traits_tuple_impl<SpaceTuple>;

  using point_type = typename self_traits::point_type;
  using point_difference_type = typename self_traits::point_difference_type;

  using distance_metric_type = TupleDistanceMetric;
  using random_sampler_type = default_random_sampler;

  static constexpr std::size_t dimensions = self_traits::dimensions;

  /**
   * Parametrized and default constructor.
   * \param aSpaces The space tuple to initialize the spaces with.
   * \param aDist The distance metric functor on the space-tuple.
   */
  explicit metric_space_tuple(SpaceTuple aSpaces = SpaceTuple{},
                              TupleDistanceMetric aDist = TupleDistanceMetric{})
      : m_spaces(std::move(aSpaces)), m_dist(std::move(aDist)) {}

  /*************************************************************************
  *                             MetricSpaceConcept
  * **********************************************************************/

  /**
   * Returns the distance between two points.
   */
  double distance(const point_type& p1, const point_type& p2) const {
    return m_dist(p1, p2, m_spaces);
  }

  /**
   * Returns the norm of the difference between two points.
   */
  double norm(const point_difference_type& dp) const {
    return m_dist(dp, m_spaces);
  }

  friend TupleDistanceMetric& get(distance_metric_t /*unused*/, self& space) {
    return space.m_dist;
  }

  friend const TupleDistanceMetric& get(distance_metric_t /*unused*/,
                                        const self& space) {
    return space.m_dist;
  }

  /*************************************************************************
   *                             PointDistributionConcept
   * **********************************************************************/

  /**
   * Generates a random point in the space, uniformly distributed.
   */
  point_type random_point() const {
    return arithmetic_tuple_details::tuple_for_each(
        m_spaces, metric_space_tuple_details::mst_random_point_computer());
  }

  /*************************************************************************
   *                             TopologyConcept
   * **********************************************************************/

  /**
   * Returns the difference between two points (a - b).
   */
  point_difference_type difference(const point_type& p1,
                                   const point_type& p2) const {
    return arithmetic_tuple_details::tuple_for_each(
        m_spaces, p1, p2,
        metric_space_tuple_details::mst_difference_computer());
  }

  /**
   * Returns the origin of the space (the lower-limit).
   */
  point_type origin() const {
    return arithmetic_tuple_details::tuple_for_each(
        m_spaces, metric_space_tuple_details::mst_origin_computer());
  }

  /**
   * Returns the addition of a point-difference to a point.
   */
  point_type adjust(const point_type& p1,
                    const point_difference_type& dp) const {
    return arithmetic_tuple_details::tuple_for_each(
        m_spaces, p1, dp, metric_space_tuple_details::mst_adjust_computer());
  }

  /*************************************************************************
  *                             LieGroupConcept
  * **********************************************************************/

  /**
   * Returns a point which is at a fraction between two points a to b.
   */
  point_type move_position_toward(const point_type& p1, double d,
                                  const point_type& p2) const {
    return arithmetic_tuple_details::tuple_for_each(
        m_spaces, p1, p2,
        metric_space_tuple_details::mst_move_position_toward_computer(d));
  }

  /**
   * Returns a point which is at a fraction between two points a to b.
   */
  point_type move_position_back_to(const point_type& p1, double d,
                                   const point_type& p2) const {
    return arithmetic_tuple_details::tuple_for_each(
        m_spaces, p1, p2,
        metric_space_tuple_details::mst_move_position_back_to_computer(d));
  }

  /*************************************************************************
  *                             BoundedSpaceConcept
  * **********************************************************************/

  /**
   * Brings a given point back with the bounds of the space.
   */
  void bring_point_in_bounds(point_type& p1) const {
    arithmetic_tuple_details::tuple_for_each(
        m_spaces, p1,
        metric_space_tuple_details::mst_bring_point_in_bounds_computer());
  }

  /**
   * Returns the difference between two points (a - b).
   */
  point_difference_type get_diff_to_boundary(const point_type& p1) const {
    return arithmetic_tuple_details::tuple_for_each(
        m_spaces, p1,
        metric_space_tuple_details::mst_get_diff_to_boundary_computer());
  }

  /**
   * Returns the addition of a point-difference to a point.
   */
  bool is_in_bounds(const point_type& p1) const {
    bool result = true;
    arithmetic_tuple_details::tuple_for_each(
        m_spaces, p1,
        metric_space_tuple_details::mst_is_in_bounds_computer(result));
    return result;
  }

  /**
   * This function returns the space at a given index.
   * \tparam Idx The index of the space.
   */
  template <int Idx>
  const arithmetic_tuple_element_t<Idx, SpaceTuple>& get_space_impl() const {
    return get<Idx>(m_spaces);
  }

  /**
   * This function returns the space at a given index.
   * \tparam Idx The index of the space.
   */
  template <int Idx>
  arithmetic_tuple_element_t<Idx, SpaceTuple>& get_space_impl() {
    return get<Idx>(m_spaces);
  }

  /**
   * This function returns the space at a given index.
   * \tparam Idx The index of the space.
   */
  template <int Idx>
  const arithmetic_tuple_element_t<Idx, SpaceTuple>& get_impl() const {
    return get<Idx>(m_spaces);
  }

  /**
   * This function returns the space at a given index.
   * \tparam Idx The index of the space.
   */
  template <int Idx>
  arithmetic_tuple_element_t<Idx, SpaceTuple>& get_impl() {
    return get<Idx>(m_spaces);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(m_spaces) & RK_SERIAL_SAVE_WITH_NAME(m_dist);
  }
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(m_spaces) & RK_SERIAL_LOAD_WITH_NAME(m_dist);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240000A, 1, "metric_space_tuple",
                              shared_object)
};

/**
 * This function returns the space at a given index.
 * \tparam Idx The index of the space.
 */
template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
const auto& get_space(
    const metric_space_tuple<SpaceTuple, TupleDistanceMetric>& s) {
  return s.template get_space_impl<Idx>();
}

/**
 * This function returns the space at a given index.
 * \tparam Idx The index of the space.
 */
template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
auto& get_space(metric_space_tuple<SpaceTuple, TupleDistanceMetric>& s) {
  return s.template get_space_impl<Idx>();
}

/**
 * This function returns the space at a given index.
 * \tparam Idx The index of the space.
 */
template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
const auto& get(const metric_space_tuple<SpaceTuple, TupleDistanceMetric>& s) {
  return s.template get_impl<Idx>();
}

/**
 * This function returns the space at a given index.
 * \tparam Idx The index of the space.
 */
template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
auto& get(metric_space_tuple<SpaceTuple, TupleDistanceMetric>& s) {
  return s.template get_impl<Idx>();
}

/**
 * This meta-function can be used to glue together a number of spaces of the same type into a tuple.
 * This class will generate a metric_space_tuple class, which has N spaces of type SpaceType.
 *
 * \tparam SpaceType The type of the spaces to glue together as a metric-space tuple.
 * \tparam N The number of spaces to glue together as a metric-space tuple.
 * \tparam TupleDistanceMetric A distance metric type which models the DistanceMetricConcept and operates on a
 *space-tuple (e.g. arithmetic_tuple).
 */
template <typename SpaceType, std::size_t N, typename TupleDistanceMetric>
struct metric_space_array {
  char cannot_instantiation_the_general_template[0];  // NOLINT
};

template <typename SpaceType, std::size_t N,
          typename TupleDistanceMetric = manhattan_tuple_distance>
using metric_space_array_t =
    typename metric_space_array<SpaceType, N, TupleDistanceMetric>::type;

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType, 1, TupleDistanceMetric> {
  using type =
      metric_space_tuple<arithmetic_tuple<SpaceType>, TupleDistanceMetric>;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType, 2, TupleDistanceMetric> {
  using type = metric_space_tuple<arithmetic_tuple<SpaceType, SpaceType>,
                                  TupleDistanceMetric>;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType, 3, TupleDistanceMetric> {
  using type =
      metric_space_tuple<arithmetic_tuple<SpaceType, SpaceType, SpaceType>,
                         TupleDistanceMetric>;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType, 4, TupleDistanceMetric> {
  using type = metric_space_tuple<
      arithmetic_tuple<SpaceType, SpaceType, SpaceType, SpaceType>,
      TupleDistanceMetric>;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType, 5, TupleDistanceMetric> {
  using type = metric_space_tuple<
      arithmetic_tuple<SpaceType, SpaceType, SpaceType, SpaceType, SpaceType>,
      TupleDistanceMetric>;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType, 6, TupleDistanceMetric> {
  using type =
      metric_space_tuple<arithmetic_tuple<SpaceType, SpaceType, SpaceType,
                                          SpaceType, SpaceType, SpaceType>,
                         TupleDistanceMetric>;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType, 7, TupleDistanceMetric> {
  using type = metric_space_tuple<
      arithmetic_tuple<SpaceType, SpaceType, SpaceType, SpaceType, SpaceType,
                       SpaceType, SpaceType>,
      TupleDistanceMetric>;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType, 8, TupleDistanceMetric> {
  using type = metric_space_tuple<
      arithmetic_tuple<SpaceType, SpaceType, SpaceType, SpaceType, SpaceType,
                       SpaceType, SpaceType, SpaceType>,
      TupleDistanceMetric>;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType, 9, TupleDistanceMetric> {
  using type = metric_space_tuple<
      arithmetic_tuple<SpaceType, SpaceType, SpaceType, SpaceType, SpaceType,
                       SpaceType, SpaceType, SpaceType, SpaceType>,
      TupleDistanceMetric>;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType, 10, TupleDistanceMetric> {
  using type = metric_space_tuple<
      arithmetic_tuple<SpaceType, SpaceType, SpaceType, SpaceType, SpaceType,
                       SpaceType, SpaceType, SpaceType, SpaceType, SpaceType>,
      TupleDistanceMetric>;
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_METRIC_SPACE_TUPLE_H_

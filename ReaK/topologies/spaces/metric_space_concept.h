/**
 * \file metric_space_concept.h
 *
 * This library defines the traits and concepts that pertain to what can be considered
 * a metric-space, as used in ReaK::pp. Metric-spaces are based on the Topology concept
 * from the Boost.Graph library, but with additional requirements which are needed
 * in algorithms tailored for a metric-space (see metric_space_search.hpp). Basically,
 * the concept of a metric-space in ReaK::pp corresponds to the mathematical concept of
 * a metric-space (see wikipedia or any decent math book).
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

#ifndef REAK_TOPOLOGIES_SPACES_METRIC_SPACE_CONCEPT_H_
#define REAK_TOPOLOGIES_SPACES_METRIC_SPACE_CONCEPT_H_

#include "ReaK/core/serialization/serializable.h"

#include <cmath>
#include <type_traits>

/** Main namespace for ReaK */
namespace ReaK::pp {

/**
 * This traits class defines the types and constants associated to a topology.
 */
template <typename Space>
struct topology_traits {
  using space_type = std::decay_t<Space>;
  /** The type that describes a point in the space. */
  using point_type = typename space_type::point_type;
  /** The type that describes a difference between points in the space. */
  using point_difference_type = typename space_type::point_difference_type;

  /** The dimensions of the space (0 if unknown at compile-time). */
  static constexpr std::size_t dimensions = space_type::dimensions;
};

template <typename Space>
struct topology_point_type {
  using type = typename topology_traits<Space>::point_type;
};
template <typename Space>
using topology_point_type_t = typename topology_point_type<Space>::type;

template <typename Space>
struct topology_point_difference_type {
  using type = typename topology_traits<Space>::point_difference_type;
};
template <typename Space>
using topology_point_difference_type_t =
    typename topology_point_difference_type<Space>::type;

/**
 * This concept defines the requirements to fulfill in order to model a topology
 * as used in ReaK::pp.
 *
 * Valid expressions:
 *
 * dp = space.difference(p1,p2);  The difference (pd) between two points (p1,p2) can be obtained.
 *
 * p1 = space.origin();  The origin of the space can be obtained.
 *
 * p2 = space.adjust(p1,dp);  A point-difference can be scaled (d * pd), added / subtracted to another point-difference
 *and added to a point (p1) to obtain an adjusted point.
 */
template <typename Space>
concept Topology = requires(const Space& space) {
  space.origin();
  space.difference(space.origin(), space.origin());
  {
    space.adjust(space.origin(),
                 space.difference(space.origin(), space.origin()))
    } -> std::convertible_to<std::decay_t<decltype(space.origin())>>;
};

/**
 * This concept defines the requirements to fulfill in order to model a Lie Group
 * as used in ReaK::pp. Basically, a Lie Group is a topology on which the point-difference
 * type is an arithmetic type (i.e. vector-space).
 *
 * Valid expressions:
 *
 * dp = d * dp + dp - dp;  The differences can be added, subtracted, and multiplied by a scalar.
 *
 * dp = -dp;  The differences can be reversed.
 *
 * dp -= dp;  The differences can be subtracted-and-stored.
 *
 * dp += dp;  The differences can be added-and-stored.
 *
 * dp *= d;  The differences can be multiplied-and-stored by a scalar.
 */
template <typename Space>
concept LieGroup = Topology<Space>&& requires(
    topology_point_difference_type_t<Space>& dp, double d) {
  {
    d* dp + dp - dp
    } -> std::convertible_to<topology_point_difference_type_t<Space>>;
  { -dp } -> std::convertible_to<topology_point_difference_type_t<Space>>;
  dp -= dp;
  dp += dp;
  dp *= d;
};

/**
 * This tag-type is used to identify (during a "get" call) that the distance-metric object is
 * to be fetched.
 */
enum distance_metric_t { distance_metric };

/**
 * This concept defines the requirements to fulfill in order to model a distance-metric
 * as used in ReaK::pp. A distance-metric is essentially a callable type that can compute
 * both the distance between two points and the corresponding norm of a difference between
 * two points.
 *
 * Required concepts:
 *
 * Space should model the Topology.
 *
 * Valid expressions:
 *
 * d = dist(p1, p2, space);  The distance (d) can be obtained by calling the distance metric (dist) on two points
 *(p1,p2) and providing a const-ref to the topology (space).
 *
 * d = dist(dp, space);  The distance (d) can be obtained by calling the distance metric (dist) on a point-difference
 *(dp) and providing a const-ref to the topology (space).
 */
template <typename Metric, typename Space>
concept DistanceMetric = Topology<Space>&& requires(
    const Metric& dist, const Space& space,
    const topology_point_type_t<Space>& p,
    const topology_point_difference_type_t<Space>& d) {
  { dist(p, p, space) } -> std::convertible_to<double>;
  { dist(d, space) } -> std::convertible_to<double>;
};

/**
 * This traits class defines the types and constants associated to a metric-space.
 * \tparam MetricSpace The topology type for which the metric-space traits are sought.
 */
template <typename Space>
struct metric_space_traits {
  using space_type = std::decay_t<Space>;
  /** The type that describes the distance-metric type for the space. */
  using distance_metric_type = typename space_type::distance_metric_type;
};

template <typename Space>
struct metric_space_distance_metric {
  using type = typename metric_space_traits<Space>::distance_metric_type;
};

template <typename Space>
using metric_space_distance_metric_t =
    typename metric_space_distance_metric<Space>::type;

/**
 * This concept defines the requirements to fulfill in order to model a metric-space
 * as used in ReaK::pp. A metric-space is a special kind of topology which has a
 * distance metric (in theory, satisfying triangular inequality).
 *
 * Valid expressions:
 *
 * dist = get(distance_metric,space);  The distance-metric can be obtained by a tagged "get" call on the metric-space.
 *
 * p1 = space.move_position_toward(p1,d,p2);  A point can be obtained by moving a fraction (d) away from one point (p1)
 *to another (p2).
 *
 * \tparam MetricSpace The topology type to be checked for this concept.
 */
template <typename Space>
concept MetricSpace = Topology<Space>&& requires(
    const Space& space, const topology_point_type_t<Space>& p, double d) {
  { get(distance_metric, space) } -> DistanceMetric<Space>;
  {
    space.move_position_toward(p, d, p)
    } -> std::convertible_to<topology_point_type_t<Space>>;
};

template <typename Space>
static constexpr bool is_metric_space_v = MetricSpace<Space>;

template <typename Space>
struct is_metric_symmetric : std::true_type {};

template <typename Space>
static constexpr bool is_metric_symmetric_v = is_metric_symmetric<Space>::value;

/**
 * This class is the default distance metric functor which models the DistanceMetric.
 * This class will simply rely on the distance and norm functions included in the
 * given topology (assuming it models the MetricSpaceConcept).
 * \note Do not use this distance metric to define a topology, because it will be cyclic (infinite recursion).
 */
struct default_distance_metric : public serializable {

  default_distance_metric() = default;

  template <Topology Space>
  explicit default_distance_metric(const Space& /*unused*/) {}

  /**
   * This function returns the distance between two points on a topology.
   * \tparam Point The point-type.
   * \tparam Space The topology.
   * \param a The first point.
   * \param b The second point.
   * \param s The topology or space on which the points lie.
   * \return The distance between two points on a topology.
   */
  template <typename Point, Topology Space>
  double operator()(const Point& a, const Point& b, const Space& s) const {
    return s.distance(a, b);
  }

  /**
   * This function returns the norm of a difference between two points on a topology.
   * \tparam PointDiff The point-difference-type.
   * \tparam Space The topology.
   * \param a The point-difference.
   * \param s The topology or space on which the points lie.
   * \return The norm of the difference between two points on a topology.
   */
  template <typename PointDiff, Topology Space>
  double operator()(const PointDiff& a, const Space& s) const {
    return s.norm(a);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {}

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {}

  RK_RTTI_MAKE_ABSTRACT_1BASE(default_distance_metric, 0xC2410000, 1,
                              "default_distance_metric", serializable)
};

template <typename Space>
auto get(distance_metric_t /*unused*/, const Space& s) {
  return metric_space_distance_metric_t<Space>{s};
}

/**
 * This class is the default distance metric functor which models the DistanceMetric.
 * This class will simply rely on the distance and norm functions included in the
 * given topology (assuming it models the MetricSpaceConcept).
 * \note Do not use this distance metric to define a topology, because it will be cyclic (infinite recursion).
 */
template <typename Metric>
struct symmetrized_metric : public serializable {
  using self = symmetrized_metric<Metric>;

  Metric unsym_distance;

  explicit symmetrized_metric(Metric aUnsymDistance = Metric{})
      : unsym_distance(aUnsymDistance) {}

  /**
   * This function returns the distance between two points on a topology.
   * \tparam Point The point-type.
   * \tparam Topology The topology.
   * \param a The first point.
   * \param b The second point.
   * \param s The topology or space on which the points lie.
   * \return The distance between two points on a topology.
   */
  template <typename Point, Topology Space>
  requires DistanceMetric<Metric, Space> double operator()(
      const Point& a, const Point& b, const Space& s) const {
    double d_left = unsym_distance(a, b, s);
    double d_right = unsym_distance(b, a, s);
    return (d_left < d_right ? d_left : d_right);
  }
  /**
   * This function returns the norm of a difference between two points on a topology.
   * \tparam PointDiff The point-difference-type.
   * \tparam Topology The topology.
   * \param a The point-difference.
   * \param s The topology or space on which the points lie.
   * \return The norm of the difference between two points on a topology.
   */
  template <typename PointDiff, Topology Space>
  requires DistanceMetric<Metric, Space> double operator()(
      const PointDiff& a, const Space& s) const {
    return unsym_distance(a, s);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(unsym_distance);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(unsym_distance);
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC241000B, 1, "symmetrized_metric",
                              serializable)
};

enum unsymmetrized_metric_t { unsymmetrized_metric };

template <typename Metric>
struct unsymmetrize {
  using type = Metric;
};

template <typename Metric>
using unsymmetrize_t = typename unsymmetrize<Metric>::type;

template <typename Metric>
struct unsymmetrize<symmetrized_metric<Metric>> {
  using type = Metric;
};

template <typename Metric>
const Metric& get(unsymmetrized_metric_t /*unused*/, const Metric& d) {
  return d;
}

template <typename Metric>
Metric get(unsymmetrized_metric_t /*unused*/,
           const symmetrized_metric<Metric>& d) {
  return d.unsym_distance;
}

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_METRIC_SPACE_CONCEPT_H_

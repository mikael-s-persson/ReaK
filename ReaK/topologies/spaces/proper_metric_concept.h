/**
 * \file proper_metric_concept.h
 *
 * This library defines the traits and concepts that pertain to what can be considered
 * a proper metric-space, as used in ReaK::pp. Proper metric-spaces are based on the
 * metric-space concept, but with additional requirements tailored for algorithms that
 * need to have the guarantee that the distance values are always finite (never NaN or Inf)
 * and that the triangular inequality is guaranteed to hold true. Basically,
 * the concept of a proper metric-space in ReaK::pp corresponds to the mathematical concept of
 * a metric-space with strict obedience to the four rules on the metric: non-negative, identity
 * of indiscernibles, symmetry and triangular inequality. The more relaxed concept of metric-space
 * used in ReaK::pp does not require symmetry or triangular inequality to hold true (although, it
 * is usually the case and always preferred), this is because many algorithms can work just as well
 * with these relaxations, but some cannot, and thus, need to rely on the concepts and traits defined
 * in this header file to appeal to a stricter compliance to the aforementioned rules.
 *
 * Note that, by default, metric-spaces in ReaK are considered to be proper metric-spaces
 * and their distance metrics are used as the proper metrics. However, the functions from
 * this header provides a path for the implementers to introduce a new distance metric in
 * cases where the plain distance metric is not proper.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2013
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

#ifndef REAK_TOPOLOGIES_SPACES_PROPER_METRIC_CONCEPT_H_
#define REAK_TOPOLOGIES_SPACES_PROPER_METRIC_CONCEPT_H_

#include "ReaK/core/serialization/serializable.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"

#include <concepts>
#include <type_traits>

/** Main namespace for ReaK */
namespace ReaK::pp {

/**
 * This tag-type is used to identify (during a "get" call) that the proper distance-metric object is
 * to be fetched.
 */
enum proper_metric_t { proper_metric };

/**
 * This meta-function gets the type of the proper metric associated to a metric-space and its distance metric.
 */
template <typename Space>
struct get_proper_metric {
  using type = metric_space_distance_metric_t<Space>;
};

template <typename Space>
using get_proper_metric_t = typename get_proper_metric<Space>::type;

/**
 * This meta-function turns a distance-metric type to a proper distance metric type.
 */
template <typename Metric>
struct get_proper_metric_from_metric {
  using type = Metric;
};

template <typename Metric>
using get_proper_metric_from_metric_t =
    typename get_proper_metric_from_metric<Metric>::type;

/**
 * This concept defines the requirements to fulfill in order to model a proper metric-space
 * as used in ReaK::pp. A proper metric-space is a special kind of metric-space which has a
 * proper distance metric (must satisfy triangular inequality and always yield a valid distance
 * value (never NaN or Inf)).
 *
 * Valid expressions:
 *
 * pdist = get(proper_metric,space);  The proper-metric can be obtained by a tagged "get" call on the metric-space.
 *
 * \tparam MetricSpace The topology type to be checked for this concept.
 */
template <typename Space>
concept ProperMetricSpace =
    MetricSpace<Space>&& DistanceMetric<get_proper_metric_t<Space>, Space>&&
    requires(const Space& space) {
  {
    get(proper_metric, space)
    } -> std::convertible_to<get_proper_metric_t<Space>>;
};

template <MetricSpace Space>
auto get(proper_metric_t /*unused*/, const Space& space) {
  return get_proper_metric_t<Space>{get(distance_metric, space)};
}

/**
 * This class is the default proper distance metric functor which models the DistanceMetricConcept.
 * This class will simply rely on the proper_distance and proper_norm functions included in the
 * given topology (assuming it models the ProperMetricSpaceConcept).
 * \note Do not use this distance metric to define a topology, because it will be cyclic (infinite recursion).
 */
struct default_proper_metric : public serializable {

  default_proper_metric() = default;

  template <typename SpaceOrMetric>
  explicit default_proper_metric(const SpaceOrMetric& /*unused*/) {}

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
  double operator()(const Point& a, const Point& b, const Space& s) const {
    return s.proper_distance(a, b);
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
  double operator()(const PointDiff& a, const Space& s) const {
    return s.proper_norm(a);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {}

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {}

  RK_RTTI_MAKE_ABSTRACT_1BASE(default_proper_metric, 0xC241000E, 1,
                              "default_proper_metric", serializable)
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_PROPER_METRIC_CONCEPT_H_

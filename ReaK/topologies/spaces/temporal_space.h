/**
 * \file temporal_space.h
 *
 * This library provides an implementation of a temporal-space which augments a
 * topology with a temporal dimension (time-stamp).
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

#ifndef REAK_TOPOLOGIES_SPACES_TEMPORAL_SPACE_H_
#define REAK_TOPOLOGIES_SPACES_TEMPORAL_SPACE_H_

#include <utility>
#include "ReaK/core/base/defs.h"
#include "ReaK/core/base/named_object.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/proper_metric_concept.h"
#include "ReaK/topologies/spaces/reversible_space_concept.h"
#include "ReaK/topologies/spaces/temporal_space_concept.h"

#include "ReaK/topologies/spaces/default_random_sampler.h"
#include "ReaK/topologies/spaces/temporal_distance_metrics.h"

namespace ReaK::pp {

/**
 * This type represents the points of a temporal-space.
 */
template <typename SpacePoint, typename TimePoint = double>
struct temporal_point : public serializable {
  using self = temporal_point<SpacePoint, TimePoint>;

  TimePoint time;  ///< Holds the time.
  SpacePoint pt;   ///< Holds the spatial-point.

  /**
   * Default constructor.
   */
  temporal_point() : time(), pt() {}
  /**
   * Constructor from a time and spatial point.
   * \param aTime The time associated to the space-time point.
   * \param aPt The spatial-point associated to the space-time point.
   */
  temporal_point(const TimePoint& aTime, SpacePoint aPt)
      : time(aTime), pt(std::move(aPt)) {}

  temporal_point(const self& rhs) : time(rhs.time), pt(rhs.pt) {}

  self& operator=(const self& rhs) {
    time = rhs.time;
    pt = rhs.pt;
    return *this;
  }

  temporal_point(self&& rhs) noexcept
      : time(std::move(rhs.time)), pt(std::move(rhs.pt)) {}

  self& operator=(self&& rhs) noexcept {
    time = std::move(rhs.time);
    pt = std::move(rhs.pt);
    return *this;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(time) & RK_SERIAL_SAVE_WITH_NAME(pt);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(time) & RK_SERIAL_LOAD_WITH_NAME(pt);
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0x0000002E, 1, "temporal_point",
                              serializable)
};

struct temporal_point_time_ordering {
  template <typename TempPoint>
  bool operator()(const TempPoint& lhs, const TempPoint& rhs) const {
    return lhs.time < rhs.time;
  }
};

template <typename Vector, typename SpacePoint, typename TimePoint>
void to_vect_impl(Vector& lhs,
                  const temporal_point<SpacePoint, TimePoint>& rhs) {
  using ReaK::detail::to_vect_impl;  // for ADL
  to_vect_impl(lhs, rhs.time);
  to_vect_impl(lhs, rhs.pt);
}

template <typename SpacePoint, typename TimePoint, typename Vector2>
void from_vect_impl(temporal_point<SpacePoint, TimePoint>& lhs,
                    const Vector2& rhs, std::size_t& i) {
  using ReaK::detail::from_vect_impl;  // for ADL
  from_vect_impl(lhs.time, rhs, i);
  from_vect_impl(lhs.pt, rhs, i);
}

/**
 * This nested type represents the difference between two points of the temporal-space.
 */
template <typename SpaceDiff, typename TimeDiff = double>
struct temporal_point_difference : public serializable {
  using self = temporal_point_difference<SpaceDiff, TimeDiff>;

  TimeDiff time;  ///< Holds the time-difference.
  SpaceDiff pt;   ///< Holds the spatial-difference.

  /**
   * Default constructor.
   */
  temporal_point_difference() : time(), pt() {}
  /**
   * Constructor from a time and space difference.
   * \param aTime The time difference.
   * \param aPt The spatial-difference.
   */
  temporal_point_difference(const TimeDiff& aTime, SpaceDiff aPt)
      : time(aTime), pt(std::move(aPt)) {}

  temporal_point_difference(const self& rhs) : time(rhs.time), pt(rhs.pt) {}

  self& operator=(const self& rhs) {
    time = rhs.time;
    pt = rhs.pt;
    return *this;
  }

  temporal_point_difference(self&& rhs) noexcept
      : time(std::move(rhs.time)), pt(std::move(rhs.pt)) {}

  self& operator=(self&& rhs) noexcept {
    time = std::move(rhs.time);
    pt = std::move(rhs.pt);
    return *this;
  }

  self operator-() const { return self(-time, -pt); }

  friend self operator*(const self& a, double b) {
    return self(a.time * b, a.pt * b);
  }

  friend self operator*(double a, const self& b) {
    return self(a * b.time, a * b.pt);
  }

  self& operator+=(const self& b) {
    time += b.time;
    pt += b.pt;
    return *this;
  }

  self& operator-=(const self& b) {
    time -= b.time;
    pt -= b.pt;
    return *this;
  }

  friend self operator+(const self& a, const self& b) {
    self result;
    result.time = a.time + b.time;
    result.pt = a.pt + b.pt;
    return result;
  }

  friend self operator-(const self& a, const self& b) {
    self result;
    result.time = a.time - b.time;
    result.pt = a.pt - b.pt;
    return result;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(time) & RK_SERIAL_SAVE_WITH_NAME(pt);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(time) & RK_SERIAL_LOAD_WITH_NAME(pt);
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0x0000002F, 1, "temporal_point_difference",
                              serializable)
};

template <typename Vector, typename SpaceDiff, typename TimeDiff>
void to_vect_impl(Vector& lhs,
                  const temporal_point_difference<SpaceDiff, TimeDiff>& rhs) {
  using ReaK::detail::to_vect_impl;  // for ADL
  to_vect_impl(lhs, rhs.time);
  to_vect_impl(lhs, rhs.pt);
}

template <typename SpaceDiff, typename TimeDiff, typename Vector2>
void from_vect_impl(temporal_point_difference<SpaceDiff, TimeDiff>& lhs,
                    const Vector2& rhs, std::size_t& i) {
  using ReaK::detail::from_vect_impl;  // for ADL
  from_vect_impl(lhs.time, rhs, i);
  from_vect_impl(lhs.pt, rhs, i);
}

/**
 * This class template can be used to transform a spatial topological map into a
 * topological map over a temporal space.
 * \tparam SpatialMap A spatial topological map type. (held by value)
 */
template <typename SpatialMap>
class temporal_topo_map : public serializable {
 private:
  SpatialMap spatial_map;

 public:
  using self = temporal_topo_map<SpatialMap>;

  /**
   * Parametric and default constructor.
   * \param aSMap The underlying spatial map.
   */
  explicit temporal_topo_map(const SpatialMap& aSMap = SpatialMap())
      : spatial_map(aSMap) {}

  /**
   * This function template maps a given temporal point in an input temporal space
   * to an output temporal point in the output temporal space, given that the
   * spatial map type of this class template is able to perform the mapping in
   * the space topology level.
   * \tparam PointType The point-type of the input temporal space.
   * \tparam InSpace The type of the input temporal space.
   * \tparam OutSpace The type of the output temporal space.
   * \param pt The temporal point in the input temporal space.
   * \param space_in The input temporal space.
   * \param space_out The output temporal space.
   * \return A temporal point in the output temporal space.
   */
  template <typename PointType, typename InSpace, typename OutSpace>
  topology_point_type_t<OutSpace> map_to_space(
      const PointType& pt, const InSpace& space_in,
      const OutSpace& space_out) const {
    return {pt.time,
            spatial_map.map_to_space(pt.pt, space_in.get_space_topology(),
                                     space_out.get_space_topology())};
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(spatial_map);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(spatial_map);
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2400038, 1, "temporal_topo_map",
                              serializable)
};

/**
 * This class can be used to transform a time-space points into space-only points, i.e.,
 * it extracts the spatial component only of the time-space point.
 */
class extract_spatial_component : public serializable {
 public:
  using self = extract_spatial_component;

  /**
   * Parametric and default constructor.
   */
  extract_spatial_component() = default;

  /**
   * This function template maps a given temporal point in an input temporal space
   * to an output spatial point in the output space.
   * \tparam PointType The point-type of the input temporal space.
   * \tparam InSpace The type of the input temporal space.
   * \tparam OutSpace The type of the output temporal space.
   * \param pt The temporal point in the input temporal space.
   * \return A spatial point in the output space.
   */
  template <typename PointType, typename InSpace, typename OutSpace>
  topology_point_type_t<OutSpace> map_to_space(
      const PointType& pt, const InSpace& /*unused*/,
      const OutSpace& /*unused*/) const {
    return pt.pt;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/
  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {}
  void load(serialization::iarchive& A, unsigned int /*Version*/) override {}
  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC240003E, 1, "extract_spatial_component",
                              serializable)
};

/**
 * This class implementats a temporal-space which augments a
 * topology with a temporal dimension (time-stamp). The time-dimension resides on a line-segment
 * topology (see line_segment_topology), while the spatial topology and distance-metric
 * is provided by the user. Models TemporalSpace.
 * \tparam Space The topology type which represents the spatial dimensions.
 * \tparam TimeTopology The topology type which represents the time dimension.
 * \tparam TemporalDistanceMetric The distance metric type for the temporal-space, should model the
 * TemporalDistMetric.
 */
template <Topology Space, Topology TimeTopology,
          typename TemporalDistanceMetric = spatial_distance_only>
class temporal_space : public named_object {
 public:
  using time_topology = TimeTopology;
  using space_topology = Space;

  using distance_metric_type = TemporalDistanceMetric;
  using random_sampler_type = default_random_sampler;

  using self = temporal_space<Space, TimeTopology, TemporalDistanceMetric>;

 protected:
  space_topology space;
  time_topology time;
  distance_metric_type dist;

 public:
  /**
   * Parametrized constructor.
   * \param aSpace The space topology to be used.
   * \param aTime The time topology to be used.
   * \param aDist The temporal distance metric functor to use.
   */
  explicit temporal_space(
      const std::string& aName = "",
      const space_topology& aSpace = space_topology(),             // NOLINT
      const time_topology& aTime = time_topology(),                // NOLINT
      const distance_metric_type& aDist = distance_metric_type())  // NOLINT
      : named_object(), space(aSpace), time(aTime), dist(aDist) {
    this->setName(aName);
  }

  using point_type = temporal_point<topology_point_type_t<space_topology>,
                                    topology_point_type_t<time_topology>>;
  using point_difference_type = temporal_point_difference<
      topology_point_difference_type_t<space_topology>,
      topology_point_difference_type_t<time_topology>>;

  static constexpr std::size_t dimensions =
      topology_traits<space_topology>::dimensions +
      topology_traits<time_topology>::dimensions;

  /** Returns the underlying space topology. */
  const space_topology& get_space_topology() const { return space; }
  /** Returns the underlying time topology. */
  const time_topology& get_time_topology() const { return time; }
  /** Returns the temporal distance metric functor used. */
  const distance_metric_type& get_distance_metric() const { return dist; }

  friend const TemporalDistanceMetric& get(distance_metric_t /*unused*/,
                                           const self& t_space) {
    return t_space.dist;
  }

  /** Returns the underlying space topology. */
  space_topology& get_space_topology() { return space; }
  /** Returns the underlying time topology. */
  time_topology& get_time_topology() { return time; }
  /** Returns the temporal distance metric functor used. */
  distance_metric_type& get_distance_metric() { return dist; }

  friend TemporalDistanceMetric& get(distance_metric_t /*unused*/,
                                     self& t_space) {
    return t_space.dist;
  }

  /*************************************************************************
   *                             Topology
   * **********************************************************************/

  /**
   * Returns the difference between two points (a - b).
   * \param a The first point.
   * \param b The second point.
   * \return The difference between the two points.
   */
  point_difference_type difference(const point_type& a,
                                   const point_type& b) const {
    return point_difference_type(time.difference(a.time, b.time),
                                 space.difference(a.pt, b.pt));
  }

  /**
   * Returns the addition of a point-difference to a point.
   * \param a The starting point.
   * \param delta The point-difference.
   * \return The addition of a point-difference to a point.
   */
  point_type adjust(const point_type& a,
                    const point_difference_type& delta) const {
    return point_type(time.adjust(a.time, delta.time),
                      space.adjust(a.pt, delta.pt));
  }

  /**
   * Returns the origin of the temporal-space.
   * \return The origin of the temporal-space.
   */
  point_type origin() const {
    return point_type(time.origin(), space.origin());
  }

  /**
   * Tests if a given point is within the boundary of this space.
   */
  bool is_in_bounds(const point_type& a) const {
    return time.is_in_bounds(a.time) && space.is_in_bounds(a.pt);
  }

  /*************************************************************************
  *                             MetricSpace
  * **********************************************************************/

  /**
   * Computes the distance between two points in the temporal-space.
   * \param a The first point.
   * \param b The second point.
   * \return The distance between a and b.
   */
  double distance(const point_type& a, const point_type& b) const {
    return dist(a, b, *this);
  }

  /**
   * Computes the norm of a difference between two points.
   * \param a The difference between two points.
   * \return The norm of a difference between two points.
   */
  double norm(const point_difference_type& a) const { return dist(a, *this); }

  /*************************************************************************
  *                             LieGroup
  * **********************************************************************/

  /**
   * Returns a point which is at a fraction between two points a to b.
   * \param a The first point.
   * \param fraction The fraction between the two points (0 to 1).
   * \param b The second point.
   * \return The point which is at a fraction between two points.
   */
  point_type move_position_toward(const point_type& a, double fraction,
                                  const point_type& b) const {
    return point_type(time.move_position_toward(a.time, fraction, b.time),
                      space.move_position_toward(a.pt, fraction, b.pt));
  }

  /**
   * Returns a point which is at a backward fraction between two points a to b.
   * \param a The first point.
   * \param fraction The backward fraction between the two points (0 to 1).
   * \param b The second point.
   * \return The point which is at a backward fraction between two points.
   */
  point_type move_position_back_to(const point_type& a, double fraction,
                                   const point_type& b) const {
    return point_type(time.move_position_back_to(a.time, fraction, b.time),
                      space.move_position_back_to(a.pt, fraction, b.pt));
  }

  /*************************************************************************
  *                             PointDistribution
  * **********************************************************************/

  /**
   * Returns a random point within the temporal-space.
   * \return A random point within the temporal-space.
   */
  point_type random_point() const {
    return point_type(get(random_sampler, time)(time),
                      get(random_sampler, space)(space));
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    ReaK::named_object::save(
        A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(space) & RK_SERIAL_SAVE_WITH_NAME(time) &
        RK_SERIAL_SAVE_WITH_NAME(dist);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    ReaK::named_object::load(
        A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(space) & RK_SERIAL_LOAD_WITH_NAME(time) &
        RK_SERIAL_LOAD_WITH_NAME(dist);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400004, 1, "temporal_space",
                              named_object)
};

template <typename Topology, typename TimeTopology,
          typename TemporalDistanceMetric>
struct is_metric_symmetric<
    temporal_space<Topology, TimeTopology, TemporalDistanceMetric>>
    : std::integral_constant<
          bool, is_metric_symmetric_v<Topology> &&
                    is_metric_symmetric_v<TimeTopology> &&
                    is_metric_symmetric_v<TemporalDistanceMetric>> {};

template <typename Topology, typename TimeTopology,
          typename TemporalDistanceMetric>
struct get_proper_metric<
    temporal_space<Topology, TimeTopology, TemporalDistanceMetric>>
    : get_proper_metric_from_metric<TemporalDistanceMetric> {};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_TEMPORAL_SPACE_H_

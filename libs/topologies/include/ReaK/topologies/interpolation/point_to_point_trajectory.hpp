/**
 * \file point_to_point_trajectory.hpp
 *
 * This library provides an implementation of a point-to-point trajectory within a temporal topology.
 * The trajectory is represented by a set of waypoints and all intermediate points
 * are computed with the topology's functions.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
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

#ifndef REAK_POINT_TO_POINT_TRAJECTORY_HPP
#define REAK_POINT_TO_POINT_TRAJECTORY_HPP

#include <ReaK/core/base/defs.hpp>

#include "spatial_trajectory_concept.hpp"

#include "waypoint_container.hpp"

#include <boost/concept_check.hpp>

#include <cmath>
#include <limits>
#include <utility>

namespace ReaK::pp {

/**
 * This class implements a point-to-point trajectory within a temporal topology (interpolated by whatever method is used
 * in
 * the topology's move function). The path is represented by a set of waypoints and all intermediate points
 * are computed with the topology's move_position_toward function based on the distance or fraction of travel.
 * This class models the SpatialTrajectoryConcept and SequentialTrajectoryConcept.
 * \tparam Topology The topology type on which the points and the path can reside, should model the
 * TemporalSpaceConcept.
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the
 * DistanceMetricConcept.
 */
template <typename Topology,
          typename DistanceMetric =
              typename metric_space_traits<Topology>::distance_metric_type>
class point_to_point_trajectory
    : public waypoint_container<Topology, DistanceMetric> {
 public:
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));

  using self = point_to_point_trajectory<Topology, DistanceMetric>;
  using base_class_type = waypoint_container<Topology, DistanceMetric>;

  using const_waypoint_descriptor =
      typename base_class_type::const_waypoint_descriptor;
  using const_waypoint_bounds = typename base_class_type::const_waypoint_bounds;
  using point_type = typename base_class_type::point_type;
  using topology = typename base_class_type::topology;
  using distance_metric = typename base_class_type::distance_metric;

  using waypoint_pair = typename base_class_type::waypoint_pair;

  using point_time_iterator = typename base_class_type::point_time_iterator;
  using point_fraction_iterator =
      typename base_class_type::point_fraction_iterator;

 private:
  using container_type = typename base_class_type::container_type;

  const container_type& get_waypoints() const { return this->waypoints; }
  const distance_metric& get_dist() const { return this->dist; }
  const topology& get_space() const { return *(this->space); }

  double travel_distance_impl(const point_type& a,
                              const point_type& b) const override {
    if (a.time > b.time) {
      return travel_distance_impl(b, a);
    }
    const_waypoint_bounds wpb_a = this->get_waypoint_bounds(a.time);
    const_waypoint_bounds wpb_b = this->get_waypoint_bounds(b.time);

    double sum = 0;
    if ((wpb_a.first == wpb_b.first) && (wpb_a.first == wpb_b.first)) {
      // this means that a and b are in the same segment.
      return this->dist(a, b, *(this->space));
    }

    sum += this->dist(a, wpb_a.second->second, *(this->space));

    auto it = wpb_a.second;
    auto it_prev = it;
    while (it++ != wpb_b.first) {
      sum += this->dist((it_prev++)->second, it->second, *(this->space));
    }

    sum += this->dist(wpb_b.first->second, b, *(this->space));

    return sum;
  }

  waypoint_pair get_point_at_time_impl(
      double t, const const_waypoint_bounds& wpb_a) const {
    using std::abs;

    if (wpb_a.first == wpb_a.second) {
      // one way or another, the point is at the boundary:
      waypoint_pair result(wpb_a.first, wpb_a.first->second);
      result.second.time = t;
      return result;
    }

    double t_total = wpb_a.second->first - wpb_a.first->first;
    double dt = t - wpb_a.first->first;
    if (t_total < 1e-6 * abs(wpb_a.second->first)) {
      waypoint_pair result(wpb_a.first, wpb_a.first->second);
      result.second.time = t;
      return result;
    }

    return waypoint_pair(wpb_a.first, this->space->move_position_toward(
                                          wpb_a.first->second, dt / t_total,
                                          wpb_a.second->second));
  }

  waypoint_pair move_time_diff_from_impl(const point_type& a,
                                         const const_waypoint_bounds& wpb_a,
                                         double dt) const override {
    if ((a.time + dt >= wpb_a.first->first) &&
        (a.time + dt <= wpb_a.second->first)) {
      return get_point_at_time_impl(a.time + dt, wpb_a);
    }
    return get_point_at_time_impl(a.time + dt,
                                  this->get_waypoint_bounds(a.time + dt));
  }

 public:
  /**
   * Constructs the trajectory from a space, assumes the start and end are at the origin
   * of the space.
   * \param aSpace The space on which the trajectory is.
   * \param aDist The distance metric functor that the trajectory should use.
   */
  explicit point_to_point_trajectory(
      const std::shared_ptr<topology>& aSpace = std::make_shared<topology>(),
      const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aDist) {}

  /**
   * Constructs the trajectory from a space, the start and end points.
   * \param aSpace The space on which the trajectory is.
   * \param aStart The start point of the trajectory.
   * \param aEnd The end-point of the trajectory.
   * \param aDist The distance metric functor that the trajectory should use.
   */
  point_to_point_trajectory(const std::shared_ptr<topology>& aSpace,
                            const point_type& aStart, const point_type& aEnd,
                            const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aStart, aEnd, aDist) {}

  /**
   * Constructs the trajectory from a range of points and their space.
   * \tparam ForwardIter A forward-iterator type for getting points to initialize the trajectory with.
   * \param aBegin An iterator to the first point of the trajectory.
   * \param aEnd An iterator to the on-past-last point of the trajectory.
   * \param aSpace The space on which the trajectory is.
   * \param aDist The distance metric functor that the trajectory should use.
   */
  template <typename ForwardIter>
  point_to_point_trajectory(ForwardIter aBegin, ForwardIter aEnd,
                            const std::shared_ptr<topology>& aSpace,
                            const distance_metric& aDist = distance_metric())
      : base_class_type(aBegin, aEnd, aSpace, aDist) {}

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(static_cast<base_class_type&>(lhs),
         static_cast<base_class_type&>(rhs));
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_class_type::save(
        A, base_class_type::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base_class_type::load(
        A, base_class_type::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2440017, 1, "point_to_point_trajectory",
                              base_class_type)
};

}  // namespace ReaK::pp

#endif

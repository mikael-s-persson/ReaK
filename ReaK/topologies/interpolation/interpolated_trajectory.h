/**
 * \file interpolated_trajectory.h
 *
 * This library provides an implementation of an interpolated trajectory within a temporal topology.
 * The trajectory is represented by a set of waypoints and all intermediate points
 * are computed with an interpolation functor provided as a template argument.
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_INTERPOLATED_TRAJECTORY_H_
#define REAK_TOPOLOGIES_INTERPOLATION_INTERPOLATED_TRAJECTORY_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/math/lin_alg/mat_num_exceptions.h"

#include "ReaK/topologies/interpolation/interpolator_concept.h"
#include "ReaK/topologies/interpolation/spatial_trajectory_concept.h"

#include "ReaK/topologies/interpolation/waypoint_container.h"

#include <cmath>
#include "boost/concept_check.hpp"

#include <limits>
#include <list>
#include <map>

namespace ReaK::pp {

/**
 * This class implements a trajectory in a temporal and once-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points
 * are computed with a linear interpolation. This class models the SpatialTrajectoryConcept.
 * \tparam Space The topology type on which the points and the path can reside.
 * \tparam InterpFactory The interpolation factory type which can create interpolators for on the given topology.
 * \tparam Metric The distance metric used to assess the distance between points in the path.
 */
template <TemporalSpace Space, typename InterpFactory,
          DistanceMetric<Space> Metric =
              typename metric_space_traits<Space>::distance_metric_type>
requires InterpolatorFactory<InterpFactory, Space,
                             Metric> class interpolated_trajectory
    : public waypoint_container<Space, Metric> {
 public:
  using self = interpolated_trajectory<Space, InterpFactory, Metric>;
  using base_class_type = waypoint_container<Space, Metric>;

  using interpolator_factory_type = InterpFactory;
  using interpolator_type = typename interpolator_factory_traits<
      interpolator_factory_type>::interpolator_type;

  using const_waypoint_descriptor =
      typename base_class_type::const_waypoint_descriptor;
  using const_waypoint_bounds = typename base_class_type::const_waypoint_bounds;
  using point_type = typename base_class_type::point_type;
  using topology = typename base_class_type::topology;
  using distance_metric = typename base_class_type::distance_metric;

  using waypoint_pair = typename base_class_type::waypoint_pair;

  using interpolator_map_type = std::map<const point_type*, interpolator_type>;

 protected:
  interpolator_factory_type interp_fact;
  mutable interpolator_map_type interp_segments;

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
      interpolator_type seg = interp_fact.create_interpolator(&a, &b);
      return seg.travel_distance_to(b, this->dist);
    }
    if (wpb_a.first == wpb_a.second) {
      // this means that a is before the first point in the trajectory.
      interpolator_type seg =
          interp_fact.create_interpolator(&a, &(wpb_a.second->second));
      sum += seg.travel_distance_from(a, this->dist);
    } else {
      auto it_int = interp_segments.find(&(wpb_a.first->second));
      if (it_int == interp_segments.end()) {
        interp_segments[&(wpb_a.first->second)] =
            interp_fact.create_interpolator(&(wpb_a.first->second),
                                            &(wpb_a.second->second));
        it_int = interp_segments.find(&(wpb_a.first->second));
      } else if (it_int->second.get_end_point() != &(wpb_a.second->second)) {
        it_int->second.set_segment(&(wpb_a.first->second),
                                   &(wpb_a.second->second));
      }
      sum += it_int->second.travel_distance_from(a, this->dist);
    }
    auto it = wpb_a.second;
    auto it_prev = it;
    while (it++ != wpb_b.first) {
      auto it_int = interp_segments.find(&(it_prev->second));
      if (it_int == interp_segments.end()) {
        interp_segments[&(it_prev->second)] =
            interp_fact.create_interpolator(&(it_prev->second), &(it->second));
        it_int = interp_segments.find(&(it_prev->second));
      } else if (it_int->second.get_end_point() != &(it->second)) {
        it_int->second.set_segment(&(it_prev->second), &(it->second));
      }
      sum += it_int->second.travel_distance_from(it_prev->second, this->dist);
      ++it_prev;
    }
    {
      auto it_int = interp_segments.find(&(wpb_b.first->second));
      if (it_int == interp_segments.end()) {
        if (wpb_b.first == wpb_b.second) {
          interpolator_type seg =
              interp_fact.create_interpolator(&(wpb_b.first->second), &b);
          sum += seg.travel_distance_to(b, this->dist);
        } else {
          interp_segments[&(wpb_b.first->second)] =
              interp_fact.create_interpolator(&(wpb_b.first->second),
                                              &(wpb_b.second->second));
          sum += interp_segments[&(wpb_b.first->second)].travel_distance_to(
              b, this->dist);
        }
      } else if (it_int->second.get_end_point() != &(wpb_b.second->second)) {
        it_int->second.set_segment(&(wpb_b.first->second),
                                   &(wpb_b.second->second));
        sum += it_int->second.travel_distance_to(b, this->dist);
      } else {
        sum += it_int->second.travel_distance_to(b, this->dist);
      }
    }
    return sum;
  }

  waypoint_pair get_point_at_time_impl(
      double t, const const_waypoint_bounds& wpb_a) const {

    if (wpb_a.first == wpb_a.second) {
      // one way or another, the point is at the boundary:
      waypoint_pair result(wpb_a.first, wpb_a.first->second);
      result.second.time = t;
      return result;
    }

    auto it_int = interp_segments.find(&(wpb_a.first->second));
    if (it_int == interp_segments.end()) {
      return waypoint_pair(
          wpb_a.first,
          (interp_segments[&(wpb_a.first->second)] =
               interp_fact.create_interpolator(&(wpb_a.first->second),
                                               &(wpb_a.second->second)))
              .get_point_at_time(t));
    }
    if (it_int->second.get_end_point() != &(wpb_a.second->second)) {
      it_int->second.set_segment(&(wpb_a.first->second),
                                 &(wpb_a.second->second));
    }
    return waypoint_pair(wpb_a.first, it_int->second.get_point_at_time(t));
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
   * \param aInterp The interpolator functor that the trajectory should use.
   */
  explicit interpolated_trajectory(
      const std::shared_ptr<topology>& aSpace = std::make_shared<topology>(),
      const distance_metric& aDist = distance_metric(),
      const interpolator_factory_type& aInterpFactory =
          interpolator_factory_type())
      : base_class_type(aSpace, aDist),
        interp_fact(aInterpFactory),
        interp_segments() {
    interp_fact.set_temporal_space(this->space);
  }

  /**
   * Constructs the trajectory from a space, the start and end points.
   * \param aSpace The space on which the trajectory is.
   * \param aStart The start point of the trajectory.
   * \param aEnd The end-point of the trajectory.
   * \param aDist The distance metric functor that the trajectory should use.
   * \param aInterp The interpolator functor that the trajectory should use.
   */
  interpolated_trajectory(const std::shared_ptr<topology>& aSpace,
                          const point_type& aStart, const point_type& aEnd,
                          const distance_metric& aDist = distance_metric(),
                          const interpolator_factory_type& aInterpFactory =
                              interpolator_factory_type())
      : base_class_type(aSpace, aStart, aEnd, aDist),
        interp_fact(aInterpFactory),
        interp_segments() {
    interp_fact.set_temporal_space(this->space);
  }

  /**
   * Constructs the trajectory from a range of points and their space.
   * \tparam ForwardIter A forward-iterator type for getting points to initialize the trajectory with.
   * \param aBegin An iterator to the first point of the trajectory.
   * \param aEnd An iterator to the on-past-last point of the trajectory.
   * \param aSpace The space on which the trajectory is.
   * \param aDist The distance metric functor that the trajectory should use.
   * \param aInterp The interpolator functor that the trajectory should use.
   */
  template <typename ForwardIter>
  interpolated_trajectory(ForwardIter aBegin, ForwardIter aEnd,
                          const std::shared_ptr<topology>& aSpace,
                          const distance_metric& aDist = distance_metric(),
                          const interpolator_factory_type& aInterpFactory =
                              interpolator_factory_type())
      : base_class_type(aBegin, aEnd, aSpace, aDist),
        interp_fact(aInterpFactory),
        interp_segments() {
    interp_fact.set_temporal_space(this->space);
  }

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) noexcept {
    using std::swap;
    swap(static_cast<base_class_type&>(lhs),
         static_cast<base_class_type&>(rhs));
    swap(lhs.interp_fact, rhs.interp_fact);
    swap(lhs.interp_segments, rhs.interp_segments);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_class_type::save(
        A, base_class_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(interp_fact);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base_class_type::load(
        A, base_class_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(interp_fact);
    interp_segments.clear();
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2440002, 1, "interpolated_trajectory",
                              base_class_type)
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_INTERPOLATION_INTERPOLATED_TRAJECTORY_H_

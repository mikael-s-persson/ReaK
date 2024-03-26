/**
 * \file sustained_acceleration_pulse.h
 *
 * This library provides an implementation of a trajectory within a temporal topology.
 * The path is represented by a set of waypoints and all intermediate points
 * are computed with a rate-limited sustained acceleration pulse (SAP) interpolation.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_SUSTAINED_ACCELERATION_PULSE_H_
#define REAK_TOPOLOGIES_INTERPOLATION_SUSTAINED_ACCELERATION_PULSE_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/topologies/spaces/bounded_space_concept.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/rate_limited_spaces.h"
#include "ReaK/topologies/spaces/tangent_bundle_concept.h"
#include "ReaK/topologies/spaces/temporal_space_concept.h"

#include "ReaK/topologies/interpolation/interpolated_trajectory.h"
#include "ReaK/topologies/spaces/generic_interpolator_factory.h"

#include "ReaK/topologies/interpolation/sustained_acceleration_pulse_detail.h"

#include <limits>
#include <type_traits>

namespace ReaK {
namespace pp {

/**
 * Use this tag type for some class templates that are parametrized in terms of the interpolation method used overall.
 */
struct sap_interpolation_tag {};

template <>
struct is_metric_symmetric<sap_interpolation_tag> : std::false_type {};

/**
 * This function template computes a Sustained Acceleration Pulse (SAP) interpolation between two points in a
 * temporal and twice-differentiable topology.
 * \tparam PointType The point type on the temporal and twice-differentiable topology.
 * \tparam Topology The temporal and twice-differentiable topology type.
 * \param a The starting point of the interpolation.
 * \param b The ending point of the interpolation.
 * \param t The time value at which the interpolated point is sought.
 * \param space The space on which the points reside.
 * \return The interpolated point at time t, between a and b.
 */
template <typename PointType, TemporalSpace Space>
PointType sap_interpolate(const PointType& a, const PointType& b, double t,
                          const Space& space) {
  using SpaceType = typename temporal_space_traits<Space>::space_topology;
  using TimeSpaceType = typename temporal_space_traits<Space>::time_topology;
  static_assert(MetricSpace<SpaceType>);
  static_assert(TangentBundle<SpaceType, TimeSpaceType, 0, 1, 2>);
  static_assert(
      SphereBoundedSpace<derived_N_order_space_t<SpaceType, TimeSpaceType, 1>>);
  static_assert(
      SphereBoundedSpace<derived_N_order_space_t<SpaceType, TimeSpaceType, 2>>);

  using Space0 = derived_N_order_space_t<SpaceType, TimeSpaceType, 0>;
  using PointDiff0 = topology_point_difference_type_t<Space0>;

  using Space1 = derived_N_order_space_t<SpaceType, TimeSpaceType, 1>;
  using PointType1 = topology_point_type_t<Space1>;

  using Space2 = derived_N_order_space_t<SpaceType, TimeSpaceType, 2>;

  static_assert(MetricSpace<Space0>);
  static_assert(MetricSpace<Space1>);
  static_assert(MetricSpace<Space2>);
  static_assert(LieGroup<Space0>);
  static_assert(LieGroup<Space1>);
  static_assert(LieGroup<Space2>);

  if (t <= a.time) {
    return a;
  }
  if (t >= b.time) {
    return b;
  }

  PointDiff0 delta_first_order;
  PointType1 peak_velocity;
  double delta_time = b.time - a.time;

  double min_delta_time = detail::sap_compute_interpolation_data_impl(
      a.pt, b.pt, delta_first_order, peak_velocity, space.get_space_topology(),
      space.get_time_topology(), delta_time, nullptr, 1e-6, 60);

  if (min_delta_time > delta_time) {
    delta_time = min_delta_time;
  }
  double dt = t - a.time;

  PointType result;
  result.time = t;

  detail::sap_interpolate_impl<
      max_derivation_order_v<SpaceType, TimeSpaceType>>(
      result.pt, a.pt, b.pt, delta_first_order, peak_velocity,
      space.get_space_topology(), space.get_time_topology(), dt, delta_time);

  return result;
}

template <BoundedSpace Space, typename TimeSpace>
bool sap_is_in_bounds(const topology_point_type_t<Space>& pt,
                      const Space& space, const TimeSpace& t_space) {
  using Space0 = derived_N_order_space_t<Space, TimeSpace, 0>;
  using Space1 = derived_N_order_space_t<Space, TimeSpace, 1>;
  using Space2 = derived_N_order_space_t<Space, TimeSpace, 2>;
  using PointType = topology_point_type_t<Space>;
  using Point2 = topology_point_type_t<Space2>;
  using PointDiff1 = topology_point_difference_type_t<Space1>;
  using PointDiff2 = topology_point_difference_type_t<Space2>;

  static_assert(LieGroup<Space0>);
  static_assert(LieGroup<Space1>);
  static_assert(LieGroup<Space2>);
  static_assert(SphereBoundedSpace<Space1>);
  static_assert(SphereBoundedSpace<Space2>);

  const Space1& s1 = get_space<1>(space, t_space);
  const Space2& s2 = get_space<2>(space, t_space);
  const auto& get_dist1 = get(distance_metric, s1);
  const auto& get_dist2 = get(distance_metric, s2);

  // Get the descended acceleration to the point of zero velocity (origin).
  PointDiff1 dp1 = s1.difference(s1.origin(), get<1>(pt));
  double dt1 = get_dist1(s1.origin(), get<1>(pt), s1);
  // Get the corresponding acceleration.
  Point2 p2 = lift_to_space<2>(dp1, dt1, space, t_space);
  // Get the descended jerk to reach that acceleration.
  PointDiff2 dp2 = s2.difference(p2, s2.origin());
  double dt2 = get_dist2(p2, s2.origin(), s2);

  // Check if we can safely ramp-up to that acceleration:
  PointType result_a = pt;
  detail::sap_constant_jerk_motion_impl<max_derivation_order<Space, TimeSpace>>(
      result_a, dp2, space, t_space, dt2);
  if (!space.is_in_bounds(result_a)) {
    return false;  // reject the sample.
  }

  if (dt1 > get_dist1(get<1>(result_a), get<1>(pt), s1)) {
    // This means, we didn't cross the zero-velocity during the jerk-down.

    // Get the new descended acceleration to the point of zero velocity (origin).
    dp1 = s1.difference(s1.origin(), get<1>(result_a));
    double dt1a = get_dist1(s1.origin(), get<1>(result_a), s1);

    // Check if we can safely stop before the boundary:
    detail::svp_constant_accel_motion_impl<
        max_derivation_order<Space, TimeSpace>>(result_a, dp1, space, t_space,
                                                dt1a);
    if (!space.is_in_bounds(result_a)) {
      return false;  // reject the sample.
    }
  }

  // Check if we could have ramped-down from that inverse acceleration:
  PointType result_b = pt;
  detail::sap_constant_jerk_motion_impl<max_derivation_order<Space, TimeSpace>>(
      result_b, -dp2, space, t_space, -dt2);
  if (!space.is_in_bounds(result_b)) {
    return false;  // reject the sample.
  }

  if (dt1 > get_dist1(get<1>(pt), get<1>(result_b), s1)) {
    // This means, the zero-velocity point is not in the wake of the last jerk-down.

    // Get the new descended acceleration to the point of zero velocity (origin).
    dp1 = s1.difference(s1.origin(), get<1>(result_b));
    double dt1b = get_dist1(s1.origin(), get<1>(result_b), s1);

    // Check if we could have ramped up from within the boundary:
    detail::svp_constant_accel_motion_impl<
        max_derivation_order<Space, TimeSpace>>(result_b, -dp1, space, t_space,
                                                -dt1b);
    if (!space.is_in_bounds(result_b)) {
      return false;  // reject the sample.
    }
  }

  // if this point is reached it means that the sample is acceptable:
  return true;
}

/**
 * This functor class implements a sustained acceleration pulse (SAP) interpolation in a temporal
 * and twice-differentiable topology.
 * \tparam SpaceType The topology on which the interpolation is done, should model MetricSpace and
 * DifferentiableSpace once against time.
 */
template <MetricSpace SpaceType, Topology TimeSpaceType>
class sap_interpolator {
 public:
  using self = sap_interpolator<SpaceType, TimeSpaceType>;
  using point_type = topology_point_type_t<SpaceType>;

  using Space0 = derived_N_order_space_t<SpaceType, TimeSpaceType, 0>;
  using PointType0 = topology_point_type_t<Space0>;
  using PointDiff0 = topology_point_difference_type_t<Space0>;
  using Space1 = derived_N_order_space_t<SpaceType, TimeSpaceType, 1>;
  using PointType1 = topology_point_type_t<Space1>;
  using PointDiff1 = topology_point_difference_type_t<Space1>;
  using Space2 = derived_N_order_space_t<SpaceType, TimeSpaceType, 2>;
  using PointType2 = topology_point_type_t<Space1>;
  using PointDiff2 = topology_point_difference_type_t<Space1>;

  static_assert(TangentBundle<SpaceType, TimeSpaceType, 0, 1, 2>);
  static_assert(
      SphereBoundedSpace<derived_N_order_space_t<SpaceType, TimeSpaceType, 1>>);
  static_assert(
      SphereBoundedSpace<derived_N_order_space_t<SpaceType, TimeSpaceType, 2>>);

  static_assert(MetricSpace<Space0>);
  static_assert(MetricSpace<Space1>);
  static_assert(MetricSpace<Space2>);
  static_assert(LieGroup<Space0>);
  static_assert(LieGroup<Space1>);
  static_assert(LieGroup<Space2>);

 private:
  PointDiff0 delta_first_order;
  PointType1 peak_velocity;
  double min_delta_time;
  PointType1 best_peak_velocity;

 public:
  /**
   * Default constructor.
   */
  sap_interpolator()
      : min_delta_time(std::numeric_limits<double>::infinity()) {}

  /**
   * Constructs the interpolator with its start and end points.
   * \tparam Factory The factory type that can be used to store fly-weight parameters used by the interpolator.
   * \param start_point The start point of the interpolation.
   * \param end_point The end point of the interpolation.
   * \param space The metric space on which the interpolation resides.
   * \param t_space The time-space against which the interpolation is done.
   * \param factory The factory object that stores relevant fly-weight parameters for the interpolator.
   */
  template <typename Factory>
  sap_interpolator(const point_type& start_point, const point_type& end_point,
                   double dt, const SpaceType& space,
                   const TimeSpaceType& t_space, const Factory& factory) {
    initialize(start_point, end_point, dt, space, t_space, factory);
  }

  /**
   * Initializes the interpolator with its start and end points.
   * \tparam Factory The factory type that can be used to store fly-weight parameters used by the interpolator.
   * \param start_point The start point of the interpolation.
   * \param end_point The end point of the interpolation.
   * \param dt The time difference between the start point to the end point of the interpolation.
   * \param space The metric space on which the interpolation resides.
   * \param t_space The time-space against which the interpolation is done.
   * \param factory The factory object that stores relevant fly-weight parameters for the interpolator.
   */
  template <typename Factory>
  void initialize(const point_type& start_point, const point_type& end_point,
                  double dt, const SpaceType& space,
                  const TimeSpaceType& t_space, const Factory& factory) {

    min_delta_time = detail::sap_compute_interpolation_data_impl(
        start_point, end_point, delta_first_order, peak_velocity, space,
        t_space, dt, &best_peak_velocity, factory.tolerance,
        factory.maximum_iterations);
  }

  /**
   * Computes the point at a given delta-time from the start-point.
   * \tparam Factory The factory type that can be used to store fly-weight parameters used by the interpolator.
   * \param result The result point of the interpolation.
   * \param start_point The start point of the interpolation.
   * \param end_point The end point of the interpolation.
   * \param space The metric space on which the interpolation resides.
   * \param t_space The time-space against which the interpolation is done.
   * \param dt The time difference from the start-point to the resulting interpolated point.
   * \param dt_total The time difference from the start-point to the end point.
   * \param factory The factory object that stores relevant fly-weight parameters for the interpolator.
   */
  template <typename Factory>
  void compute_point(point_type& result, const point_type& start_point,
                     const point_type& end_point, const SpaceType& space,
                     const TimeSpaceType& t_space, double dt, double dt_total,
                     const Factory& factory) const {
    if (dt <= 0.0) {
      result = start_point;
      return;
    }
    if (dt >= dt_total) {
      result = end_point;
      return;
    }

    detail::sap_interpolate_impl<
        max_derivation_order_v<SpaceType, TimeSpaceType>>(
        result, start_point, end_point, delta_first_order, peak_velocity, space,
        t_space, dt, dt_total);
  }

  /**
   * Returns the minimum travel time between the initialized start and end points.
   * \return The minimum travel time between the initialized start and end points.
   */
  double get_minimum_travel_time() const { return min_delta_time; }
};

/**
 * This class is a factory class for sustained acceleration pulse (SAP) interpolators on a temporal
 * differentiable space.
 * \tparam TemporalTopology The temporal topology on which the interpolation is done, should model TemporalSpace,
 *                          with a spatial topology that is twice-differentiable (see DifferentiableSpace) and
 *                          whose 1-order and 2-order derivative space has a spherical bound (see
 * SphereBoundedSpace).
 */
template <TemporalSpace Space>
class sap_interpolator_factory : public serializable {
 public:
  using self = sap_interpolator_factory<Space>;
  using topology = Space;
  using point_type = topology_point_type_t<Space>;

  using interpolator_type = generic_interpolator<self, sap_interpolator>;

 private:
  std::shared_ptr<Space> space;

 public:
  double tolerance;
  unsigned int maximum_iterations;

  explicit sap_interpolator_factory(const std::shared_ptr<Space>& aSpace = {},
                                    double aTolerance = 1e-6,
                                    unsigned int aMaxIter = 60)
      : space(aSpace), tolerance(aTolerance), maximum_iterations(aMaxIter) {}

  void set_temporal_space(const std::shared_ptr<Space>& aSpace) {
    space = aSpace;
  }
  const std::shared_ptr<Space>& get_temporal_space() const { return space; }

  interpolator_type create_interpolator(const point_type* pp1,
                                        const point_type* pp2) const {
    return interpolator_type(this, pp1, pp2);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(space) & RK_SERIAL_SAVE_WITH_NAME(tolerance) &
        RK_SERIAL_SAVE_WITH_NAME(maximum_iterations);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(space) & RK_SERIAL_LOAD_WITH_NAME(tolerance) &
        RK_SERIAL_LOAD_WITH_NAME(maximum_iterations);
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2430005, 1, "sap_interpolator_factory",
                              serializable)
};

/**
 * This class implements a trajectory in a temporal and twice-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points
 * are computed with a sustained acceleration pulse interpolation (limited). This class models
 * the SpatialTrajectory.
 * \tparam Topology The temporal topology on which the interpolation is done, should model TemporalSpace,
 *                  with a spatial topology that is twice-differentiable (see DifferentiableSpace) and
 *                  whose 1-order and 2-order derivative spaces have a spherical bound (see SphereBoundedSpace).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the
 * DistanceMetric.
 */
template <TemporalSpace Space,
          DistanceMetric<Space> Metric =
              typename metric_space_traits<Space>::distance_metric_type>
class sap_interp_traj
    : public interpolated_trajectory<Space, sap_interpolator_factory<Space>,
                                     Metric> {
 public:
  using self = sap_interp_traj<Space, Metric>;
  using base_class_type =
      interpolated_trajectory<Space, sap_interpolator_factory<Space>, Metric>;

  using point_type = typename base_class_type::point_type;
  using topology = typename base_class_type::topology;
  using distance_metric = typename base_class_type::distance_metric;

  /**
   * Constructs the path from a space, assumes the start and end are at the origin
   * of the space.
   * \param aSpace The space on which the path is.
   * \param aDist The distance metric functor that the path should use.
   */
  explicit sap_interp_traj(
      const std::shared_ptr<topology>& aSpace = std::make_shared<topology>(),
      const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aDist,
                        sap_interpolator_factory<Space>(aSpace)) {}

  /**
   * Constructs the path from a space, the start and end points.
   * \param aSpace The space on which the path is.
   * \param aStart The start point of the path.
   * \param aEnd The end-point of the path.
   * \param aDist The distance metric functor that the path should use.
   */
  sap_interp_traj(const std::shared_ptr<topology>& aSpace,
                  const point_type& aStart, const point_type& aEnd,
                  const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aStart, aEnd, aDist,
                        sap_interpolator_factory<Space>(aSpace)) {}

  /**
   * Constructs the path from a range of points and their space.
   * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
   * \param aBegin An iterator to the first point of the path.
   * \param aEnd An iterator to the second point of the path.
   * \param aSpace The space on which the path is.
   * \param aDist The distance metric functor that the path should use.
   */
  template <typename ForwardIter>
  sap_interp_traj(ForwardIter aBegin, ForwardIter aEnd,
                  const std::shared_ptr<topology>& aSpace,
                  const distance_metric& aDist = distance_metric())
      : base_class_type(aBegin, aEnd, aSpace, aDist,
                        sap_interpolator_factory<Space>(aSpace)) {}

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

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2440007, 1, "sap_interp_traj",
                              base_class_type)
};

}  // namespace pp

namespace rtti {

template <>
struct get_type_id<pp::sap_interpolation_tag> {
  static constexpr unsigned int ID = 4;
  static constexpr auto type_name = std::string_view{"sap_interpolation_tag"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

}  // namespace rtti
}  // namespace ReaK

#endif  // REAK_TOPOLOGIES_INTERPOLATION_SUSTAINED_ACCELERATION_PULSE_H_

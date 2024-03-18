/**
 * \file sustained_velocity_pulse.h
 *
 * This library provides an implementation of a trajectory within a temporal topology.
 * The path is represented by a set of waypoints and all intermediate points
 * are computed with a rate-limited sustained velocity pulse (SVP) interpolation.
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_SUSTAINED_VELOCITY_PULSE_H_
#define REAK_TOPOLOGIES_INTERPOLATION_SUSTAINED_VELOCITY_PULSE_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/topologies/spaces/bounded_space_concept.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/rate_limited_spaces.h"
#include "ReaK/topologies/spaces/tangent_bundle_concept.h"
#include "ReaK/topologies/spaces/temporal_space_concept.h"

#include "ReaK/topologies/spaces/generic_interpolator_factory.h"
#include "ReaK/topologies/interpolation/interpolated_trajectory.h"

#include "ReaK/topologies/interpolation/sustained_velocity_pulse_detail.h"

#include <limits>
#include <type_traits>

namespace ReaK {

namespace pp {

/**
 * Use this tag type for some class templates that are parametrized in terms of the interpolation method used overall.
 */
struct svp_interpolation_tag {};

template <>
struct is_metric_symmetric<svp_interpolation_tag> : std::false_type {};

/**
 * This function template computes a Sustained Velocity Pulse (SVP) interpolation between two points in a
 * temporal and once-differentiable topology.
 * \tparam PointType The point type on the temporal and once-differentiable topology.
 * \tparam Topology The temporal and once-differentiable topology type.
 * \param a The starting point of the interpolation.
 * \param b The ending point of the interpolation.
 * \param t The time value at which the interpolated point is sought.
 * \param space The space on which the points reside.
 * \return The interpolated point at time t, between a and b.
 */
template <typename PointType, TemporalSpace Space>
PointType svp_interpolate(const PointType& a, const PointType& b, double t,
                          const Space& space, double tolerance = 1e-6,
                          unsigned int maximum_iterations = 60) {
  using SpaceType = typename temporal_space_traits<Space>::space_topology;
  using TimeSpaceType = typename temporal_space_traits<Space>::time_topology;
  static_assert(TangentBundle<SpaceType, TimeSpaceType, 0, 1>);
  static_assert(SphereBoundedSpace<
                        derived_N_order_space_t<SpaceType, TimeSpaceType, 1>>);

  using Space0 = derived_N_order_space_t<SpaceType, TimeSpaceType, 0>;
  using PointDiff0 = topology_point_difference_type_t<Space0>;

  using Space1 = derived_N_order_space_t<SpaceType, TimeSpaceType, 1>;
  using PointType1 = topology_point_type_t<Space1>;

  static_assert(MetricSpace<Space0>);
  static_assert(MetricSpace<Space1>);
  static_assert(LieGroup<Space0>);
  static_assert(LieGroup<Space1>);

  if (t <= a.time) {
    return a;
  }
  //   if(t >= b.time)
  //     return b;

  PointDiff0 delta_first_order;
  PointType1 peak_velocity;

  double delta_time = b.time - a.time;

  double slack = detail::svp_compute_interpolation_data_impl(
      a.pt, b.pt, delta_first_order, peak_velocity, space.get_space_topology(),
      space.get_time_topology(), delta_time, nullptr, tolerance,
      maximum_iterations);

  if (slack < 0.0) {
    delta_time -= slack;
  }
  double dt = t - a.time;

  PointType result;
  result.time = t;

  detail::svp_interpolate_impl<
      max_derivation_order_v<SpaceType, TimeSpaceType>>(
      result.pt, a.pt, b.pt, delta_first_order, peak_velocity,
      space.get_space_topology(), space.get_time_topology(), dt, delta_time);

  return result;
}

template <BoundedSpace Space, Topology TimeSpace>
bool svp_is_in_bounds(const topology_point_type_t<Space>& pt,
                      const Space& space, const TimeSpace& t_space) {
  using Space0 = derived_N_order_space_t<Space, TimeSpace, 0>;
  using Space1 = derived_N_order_space_t<Space, TimeSpace, 1>;
  using PointType = topology_point_type_t<Space>;
  using PointDiff1 = topology_point_difference_type_t<Space1>;

  static_assert(LieGroup<Space0>);
  static_assert(LieGroup<Space1>);
  static_assert(SphereBoundedSpace<Space1>);

  const Space1& s1 = get_space<1>(space, t_space);
  auto& get_vel_dist = get(distance_metric, s1);

  PointDiff1 dp1 = s1.difference(s1.origin(), get<1>(pt));
  double dt = get_vel_dist(s1.origin(), get<1>(pt), s1);

  // Check if we can stop the motion before the boundary.
  PointType result = pt;
  detail::svp_constant_accel_motion_impl<
      max_derivation_order<Space, TimeSpace>>(result, dp1, space, t_space,
                                                    dt);
  if (!space.is_in_bounds(result)) {
    return false;  // reject the sample.
  }

  // Check if we could have initiated the motion from within the boundary.
  result = pt;
  detail::svp_constant_accel_motion_impl<
      max_derivation_order<Space, TimeSpace>>(result, -dp1, space,
                                                    t_space, -dt);
  return space.is_in_bounds(result);
}

/**
 * This functor class implements a sustained velocity pulse (SVP) interpolation in a temporal and once-differentiable
 * topology.
 * \tparam Space The topology on which the interpolation is done.
 * \tparam TimeSpace The time topology.
 */
template <Topology Space, Topology TimeSpace>
class svp_interpolator {
 public:
  using self = svp_interpolator<Space, TimeSpace>;
  using point_type = topology_point_type_t<Space>;

  using Space0 = derived_N_order_space_t<Space, TimeSpace, 0>;
  using PointType0 = topology_point_type_t<Space0>;
  using PointDiff0 = topology_point_difference_type_t<Space0>;
  using Space1 = derived_N_order_space_t<Space, TimeSpace, 1>;
  using PointType1 = topology_point_type_t<Space1>;
  using PointDiff1 = topology_point_difference_type_t<Space1>;

  static_assert(MetricSpace<Space0>);
  static_assert(MetricSpace<Space1>);
  static_assert(LieGroup<Space0>);
  static_assert(LieGroup<Space1>);
  static_assert(TangentBundle<Space, TimeSpace, 0, 1>);
  static_assert(SphereBoundedSpace<Space1>);

 private:
  PointDiff0 delta_first_order;
  PointType1 peak_velocity;
  double min_delta_time;
  PointType1 best_peak_velocity;

 public:
  /**
   * Default constructor.
   */
  svp_interpolator()
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
  svp_interpolator(const point_type& start_point, const point_type& end_point,
                   double dt, const Space& space,
                   const TimeSpace& t_space, const Factory& factory) {
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
                  double dt, const Space& space,
                  const TimeSpace& t_space, const Factory& factory) {

    min_delta_time = detail::svp_compute_interpolation_data_impl(
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
                     const point_type& end_point, const Space& space,
                     const TimeSpace& t_space, double dt, double dt_total,
                     const Factory& factory) const {
    if (dt <= 0.0) {
      result = start_point;
      return;
    }

    detail::svp_interpolate_impl<
        max_derivation_order_v<Space, TimeSpace>>(
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
 * This class is a factory class for sustained velocity pulse (SVP) interpolators on a temporal differentiable space.
 * \tparam TemporalTopology The temporal topology on which the interpolation is done, should model TemporalSpace,
 *                          with a spatial topology that is once-differentiable (see DifferentiableSpace) and
 *                          whose 1-order derivative space has a spherical bound (see SphereBoundedSpace).
 */
template <TemporalSpace Space>
class svp_interpolator_factory : public serializable {
 public:
  using self = svp_interpolator_factory<Space>;
  using topology = Space;
  using point_type = topology_point_type_t<Space>;

  using interpolator_type = generic_interpolator<self, svp_interpolator>;

 private:
  std::shared_ptr<topology> space;

 public:
  double tolerance;
  unsigned int maximum_iterations;

  explicit svp_interpolator_factory(
      const std::shared_ptr<topology>& aSpace = std::shared_ptr<topology>(),
      double aTolerance = 1e-6, unsigned int aMaxIter = 60)
      : space(aSpace), tolerance(aTolerance), maximum_iterations(aMaxIter){};

  void set_temporal_space(const std::shared_ptr<topology>& aSpace) {
    space = aSpace;
  }
  const std::shared_ptr<topology>& get_temporal_space() const { return space; }

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

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2430004, 1, "svp_interpolator_factory",
                              serializable)
};

/**
 * This class implements a trajectory in a temporal and once-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points
 * are computed with a sustained velocity pulse interpolation (limited). This class models
 * the SpatialTrajectoryConcept.
 * \tparam Topology The temporal topology on which the interpolation is done, should model TemporalSpaceConcept,
 *                  with a spatial topology that is once-differentiable (see DifferentiableSpaceConcept) and
 *                  whose 1-order derivative space has a spherical bound (see SphereBoundedSpaceConcept).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the
 * DistanceMetricConcept.
 */
template <TemporalSpace Space,
          DistanceMetric<Space> Metric =
              typename metric_space_traits<Space>::distance_metric_type>
class svp_interp_traj
    : public interpolated_trajectory<
          Space, svp_interpolator_factory<Space>, Metric> {
 public:
  using self = svp_interp_traj<Space, Metric>;
  using base_class_type =
      interpolated_trajectory<Space, svp_interpolator_factory<Space>, Metric>;

  using point_type = typename base_class_type::point_type;
  using topology = typename base_class_type::topology;
  using distance_metric = typename base_class_type::distance_metric;

  /**
   * Constructs the path from a space, assumes the start and end are at the origin
   * of the space.
   * \param aSpace The space on which the path is.
   * \param aDist The distance metric functor that the path should use.
   */
  explicit svp_interp_traj(
      const std::shared_ptr<topology>& aSpace = std::make_shared<topology>(),
      const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aDist,
                        svp_interpolator_factory<Space>(aSpace)) {}

  /**
   * Constructs the path from a space, the start and end points.
   * \param aSpace The space on which the path is.
   * \param aStart The start point of the path.
   * \param aEnd The end-point of the path.
   * \param aDist The distance metric functor that the path should use.
   */
  svp_interp_traj(const std::shared_ptr<topology>& aSpace,
                  const point_type& aStart, const point_type& aEnd,
                  const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aStart, aEnd, aDist,
                        svp_interpolator_factory<Space>(aSpace)) {}

  /**
   * Constructs the path from a range of points and their space.
   * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
   * \param aBegin An iterator to the first point of the path.
   * \param aEnd An iterator to the second point of the path.
   * \param aSpace The space on which the path is.
   * \param aDist The distance metric functor that the path should use.
   */
  template <typename ForwardIter>
  svp_interp_traj(ForwardIter aBegin, ForwardIter aEnd,
                  const std::shared_ptr<topology>& aSpace,
                  const distance_metric& aDist = distance_metric())
      : base_class_type(aBegin, aEnd, aSpace, aDist,
                        svp_interpolator_factory<Space>(aSpace)) {}

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

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2440006, 1, "svp_interp_traj",
                              base_class_type)
};

}  // namespace pp

namespace rtti {

template <>
struct get_type_id<pp::svp_interpolation_tag> {
  static constexpr unsigned int ID = 5;
  static constexpr auto type_name = std::string_view{"svp_interpolation_tag"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

}  // namespace rtti
}  // namespace ReaK

#endif  // REAK_TOPOLOGIES_INTERPOLATION_SUSTAINED_VELOCITY_PULSE_H_

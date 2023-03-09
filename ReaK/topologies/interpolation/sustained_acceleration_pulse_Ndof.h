/**
 * \file sustained_acceleration_pulse_Ndof.h
 *
 * This library provides an implementation of a trajectory within a temporal N-dof topology.
 * The path is represented by a set of waypoints and all intermediate points
 * are computed with a rate-limited sustained acceleration pulse (SAP) interpolation.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2013
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_SUSTAINED_ACCELERATION_PULSE_NDOF_H_
#define REAK_TOPOLOGIES_INTERPOLATION_SUSTAINED_ACCELERATION_PULSE_NDOF_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/topologies/spaces/bounded_space_concept.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/rate_limited_spaces.h"
#include "ReaK/topologies/spaces/tangent_bundle_concept.h"
#include "ReaK/topologies/spaces/temporal_space_concept.h"

#include "ReaK/topologies/interpolation/generic_interpolator_factory.h"
#include "ReaK/topologies/interpolation/interpolated_trajectory.h"

#include "ReaK/topologies/interpolation/sustained_acceleration_pulse_Ndof_detail.h"

#include "boost/concept_check.hpp"

#include <limits>
#include <type_traits>

namespace ReaK {

namespace pp {

/**
 * Use this tag type for some class templates that are parametrized in terms of the interpolation method used overall.
 */
struct sap_Ndof_interpolation_tag {};

template <>
struct is_metric_symmetric<sap_Ndof_interpolation_tag> : std::false_type {};

/**
 * This function template computes a Sustained Acceleration Pulse (SAP) interpolation between two points in a
 * temporal and twice-differentiable N-dof topology.
 * \tparam PointType The point type on the temporal and twice-differentiable topology.
 * \tparam Topology The temporal and twice-differentiable topology type.
 * \param a The starting point of the interpolation.
 * \param b The ending point of the interpolation.
 * \param t The time value at which the interpolated point is sought.
 * \param space The space on which the points reside.
 * \return The interpolated point at time t, between a and b.
 */
template <typename PointType, typename Topology>
PointType sap_Ndof_interpolate(const PointType& a, const PointType& b, double t,
                               const Topology& space) {
  using SpaceType = typename temporal_space_traits<Topology>::space_topology;
  using TimeSpaceType = typename temporal_space_traits<Topology>::time_topology;
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((MetricSpaceConcept<SpaceType>));
  BOOST_CONCEPT_ASSERT((TangentBundleConcept<SpaceType, 2, TimeSpaceType>));

  using Space0 = derived_N_order_space_t<SpaceType, TimeSpaceType, 0>;
  using PointType0 = topology_point_type_t<Space0>;

  using Space1 = derived_N_order_space_t<SpaceType, TimeSpaceType, 1>;
  using PointType1 = topology_point_type_t<Space1>;

  using Space2 = derived_N_order_space_t<SpaceType, TimeSpaceType, 2>;
  using PointType2 = topology_point_type_t<Space2>;

  BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Space0>));
  BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Space1>));
  BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Space2>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space1>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space2>));
  BOOST_CONCEPT_ASSERT((BoxBoundedSpaceConcept<Space1>));
  BOOST_CONCEPT_ASSERT((BoxBoundedSpaceConcept<Space2>));
  BOOST_CONCEPT_ASSERT((WritableVectorConcept<PointType0>));
  BOOST_CONCEPT_ASSERT((WritableVectorConcept<PointType1>));
  BOOST_CONCEPT_ASSERT((WritableVectorConcept<PointType2>));

  if (t <= a.time) {
    return a;
  }
  //   if(t >= b.time)
  //     return b;

  PointType1 peak_velocity;
  double delta_time = b.time - a.time;

  double min_delta_time = detail::sap_compute_Ndof_interpolation_data_impl(
      a.pt, b.pt, peak_velocity, space.get_space_topology(),
      space.get_time_topology(), delta_time, nullptr);

  if (min_delta_time > delta_time) {
    delta_time = min_delta_time;
  }
  double dt = t - a.time;

  PointType result;
  result.time = t;

  detail::sap_Ndof_interpolate_impl<
      max_derivation_order_v<SpaceType, TimeSpaceType>>(
      result.pt, a.pt, b.pt, peak_velocity, space.get_space_topology(),
      space.get_time_topology(), dt, delta_time);

  return result;
}

template <typename Topology, typename TimeTopology>
bool sap_Ndof_is_in_bounds(const topology_point_type_t<Topology>& pt,
                           const Topology& space, const TimeTopology& t_space) {
  BOOST_CONCEPT_ASSERT((TopologyConcept<Topology>));
  BOOST_CONCEPT_ASSERT((BoundedSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((TangentBundleConcept<Topology, 2, TimeTopology>));

  using Space0 = derived_N_order_space_t<Topology, TimeTopology, 0>;
  using Space1 = derived_N_order_space_t<Topology, TimeTopology, 1>;
  using Space2 = derived_N_order_space_t<Topology, TimeTopology, 2>;
  using Point0 = topology_point_type_t<Space0>;
  using Point1 = topology_point_type_t<Space1>;
  using Point2 = topology_point_type_t<Space2>;

  using std::abs;

  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space1>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space2>));
  BOOST_CONCEPT_ASSERT((BoundedSpaceConcept<Space0>));
  BOOST_CONCEPT_ASSERT((BoxBoundedSpaceConcept<Space1>));
  BOOST_CONCEPT_ASSERT((BoxBoundedSpaceConcept<Space2>));
  BOOST_CONCEPT_ASSERT((WritableVectorConcept<Point0>));
  BOOST_CONCEPT_ASSERT((WritableVectorConcept<Point1>));
  BOOST_CONCEPT_ASSERT((WritableVectorConcept<Point2>));

  const Space0& s0 = get_space<0>(space, t_space);
  const Space1& s1 = get_space<1>(space, t_space);
  const Space2& s2 = get_space<2>(space, t_space);

  if (!s0.is_in_bounds(get<0>(pt)) || !s1.is_in_bounds(get<1>(pt)) ||
      !s2.is_in_bounds(get<2>(pt))) {
    return false;
  }

  Point1 max_velocity = s1.get_upper_corner();
  Point2 max_acceleration = s2.get_upper_corner();

  Point0 stopping_point = get<0>(pt);
  Point0 starting_point = get<0>(pt);

  for (std::size_t i = 0; i < max_velocity.size(); ++i) {
    double dt_vp1_1st = abs(get<1>(pt)[i]);
    if (dt_vp1_1st <= std::numeric_limits<double>::epsilon()) {
      continue;
    }
    // we know that dt_vp_2nd = dt_vp_1st + dt_amax
    double dt_vp1 = dt_vp1_1st - max_acceleration[i];
    if (dt_vp1 < 0.0) {
      // means that we don't have time to reach the maximum acceleration:
      double integ = get<1>(pt)[i] * dt_vp1_1st / max_velocity[i];
      stopping_point[i] += integ;
      starting_point[i] -= integ;
    } else {
      double integ = 0.5 * get<1>(pt)[i] * (dt_vp1_1st + max_acceleration[i]) /
                     max_velocity[i];
      stopping_point[i] += integ;
      starting_point[i] -= integ;
    }
  }

  // Check if we could have initiated the motion from within the boundary or if we can stop the motion before the
  // boundary.
  return (s0.is_in_bounds(stopping_point) && s0.is_in_bounds(starting_point));
}

/**
 * This functor class implements a sustained acceleration pulse (SAP) interpolation in a temporal
 * and twice-differentiable topology.
 * \tparam SpaceType The topology on which the interpolation is done, should model MetricSpaceConcept and
 * DifferentiableSpaceConcept once against time.
 * \tparam TimeSpaceType The time topology.
 */
template <typename SpaceType, typename TimeSpaceType>
class sap_Ndof_interpolator {
 public:
  using self = sap_Ndof_interpolator<SpaceType, TimeSpaceType>;
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

  BOOST_CONCEPT_ASSERT((MetricSpaceConcept<SpaceType>));
  BOOST_CONCEPT_ASSERT((TangentBundleConcept<SpaceType, 2, TimeSpaceType>));
  BOOST_CONCEPT_ASSERT((BoxBoundedSpaceConcept<Space1>));
  BOOST_CONCEPT_ASSERT((BoxBoundedSpaceConcept<Space2>));

  BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Space0>));
  BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Space1>));
  BOOST_CONCEPT_ASSERT((MetricSpaceConcept<Space2>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space1>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space2>));
  BOOST_CONCEPT_ASSERT((WritableVectorConcept<PointType0>));
  BOOST_CONCEPT_ASSERT((WritableVectorConcept<PointType1>));
  BOOST_CONCEPT_ASSERT((WritableVectorConcept<PointType2>));

 private:
  PointType1 peak_velocity;
  double min_delta_time;
  PointType1 best_peak_velocity;

 public:
  /**
   * Default constructor.
   */
  sap_Ndof_interpolator()
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
  sap_Ndof_interpolator(const point_type& start_point,
                        const point_type& end_point, double dt,
                        const SpaceType& space, const TimeSpaceType& t_space,
                        const Factory& factory) {
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
    RK_UNUSED(factory);
    try {
      min_delta_time = detail::sap_compute_Ndof_interpolation_data_impl(
          start_point, end_point, peak_velocity, space, t_space, dt,
          &best_peak_velocity);
    } catch (optim::infeasible_problem& e) {
      RK_UNUSED(e);
      min_delta_time = std::numeric_limits<double>::infinity();
    }
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
    RK_UNUSED(factory);
    if (min_delta_time == std::numeric_limits<double>::infinity()) {
      return;
    }
    if (dt <= 0.0) {
      result = start_point;
      return;
    }

    detail::sap_Ndof_interpolate_impl<
        max_derivation_order_v<SpaceType, TimeSpaceType>>(
        result, start_point, end_point, peak_velocity, space, t_space, dt,
        dt_total);
  }

  /**
   * Returns the minimum travel time between the initialized start and end points.
   * \return The minimum travel time between the initialized start and end points.
   */
  double get_minimum_travel_time() const { return min_delta_time; }
};

/**
 * This class is a factory class for sustained acceleration pulse (SAP) interpolators on a temporal
 * differentiable N-dof space.
 * \tparam TemporalTopology The temporal topology on which the interpolation is done, should model TemporalSpaceConcept,
 *                          with a spatial topology that is twice-differentiable (see DifferentiableSpaceConcept) and
 *                          whose 1-order and 2-order derivative space has a spherical bound (see
 * SphereBoundedSpaceConcept).
 */
template <typename TemporalTopology>
class sap_Ndof_interp_factory : public serializable {
 public:
  using self = sap_Ndof_interp_factory<TemporalTopology>;
  using topology = TemporalTopology;
  using point_type = topology_point_type_t<TemporalTopology>;

  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<topology>));

  using interpolator_type = generic_interpolator<self, sap_Ndof_interpolator>;

 private:
  std::shared_ptr<topology> space;

 public:
  explicit sap_Ndof_interp_factory(const std::shared_ptr<topology>& aSpace)
      : space(aSpace) {}

  sap_Ndof_interp_factory()
      : sap_Ndof_interp_factory(std::shared_ptr<topology>()) {}

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
    A& RK_SERIAL_SAVE_WITH_NAME(space);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(space);
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2430007, 1, "sap_Ndof_interp_factory",
                              serializable)
};

/**
 * This class implements a trajectory in a temporal and twice-differentiable N-dof topology.
 * The trajectory is represented by a set of waypoints and all intermediate points
 * are computed with a sustained acceleration pulse interpolation (limited). This class models
 * the SpatialTrajectoryConcept.
 * \tparam Topology The temporal topology on which the interpolation is done, should model TemporalSpaceConcept,
 *                  with a spatial topology that is twice-differentiable (see DifferentiableSpaceConcept) and
 *                  whose 1-order and 2-order derivative spaces have a spherical bound (see SphereBoundedSpaceConcept).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the
 * DistanceMetricConcept.
 */
template <typename Topology,
          typename DistanceMetric =
              typename metric_space_traits<Topology>::distance_metric_type>
class sap_Ndof_interp_traj
    : public interpolated_trajectory<
          Topology, sap_Ndof_interp_factory<Topology>, DistanceMetric> {
 public:
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));

  using self = sap_Ndof_interp_traj<Topology, DistanceMetric>;
  using base_class_type =
      interpolated_trajectory<Topology, sap_Ndof_interp_factory<Topology>,
                              DistanceMetric>;

  using point_type = typename base_class_type::point_type;
  using topology = typename base_class_type::topology;
  using distance_metric = typename base_class_type::distance_metric;

  /**
   * Constructs the path from a space, assumes the start and end are at the origin
   * of the space.
   * \param aSpace The space on which the path is.
   * \param aDist The distance metric functor that the path should use.
   */
  explicit sap_Ndof_interp_traj(
      const std::shared_ptr<topology>& aSpace = std::make_shared<topology>(),
      const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aDist,
                        sap_Ndof_interp_factory<Topology>(aSpace)) {}

  /**
   * Constructs the path from a space, the start and end points.
   * \param aSpace The space on which the path is.
   * \param aStart The start point of the path.
   * \param aEnd The end-point of the path.
   * \param aDist The distance metric functor that the path should use.
   */
  sap_Ndof_interp_traj(const std::shared_ptr<topology>& aSpace,
                       const point_type& aStart, const point_type& aEnd,
                       const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aStart, aEnd, aDist,
                        sap_Ndof_interp_factory<Topology>(aSpace)) {}

  /**
   * Constructs the path from a range of points and their space.
   * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
   * \param aBegin An iterator to the first point of the path.
   * \param aEnd An iterator to the second point of the path.
   * \param aSpace The space on which the path is.
   * \param aDist The distance metric functor that the path should use.
   */
  template <typename ForwardIter>
  sap_Ndof_interp_traj(ForwardIter aBegin, ForwardIter aEnd,
                       const std::shared_ptr<topology>& aSpace,
                       const distance_metric& aDist = distance_metric())
      : base_class_type(aBegin, aEnd, aSpace, aDist,
                        sap_Ndof_interp_factory<Topology>(aSpace)) {}

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

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC244000F, 1, "sap_Ndof_interp_traj",
                              base_class_type)
};

}  // namespace pp

namespace rtti {

template <>
struct get_type_id<pp::sap_Ndof_interpolation_tag> {
  static constexpr unsigned int ID = 4;
  static constexpr auto type_name =
      std::string_view{"sap_Ndof_interpolation_tag"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

}  // namespace rtti
}  // namespace ReaK

#endif  // REAK_TOPOLOGIES_INTERPOLATION_SUSTAINED_ACCELERATION_PULSE_NDOF_H_

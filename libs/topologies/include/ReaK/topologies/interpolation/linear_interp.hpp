/**
 * \file linear_interp.hpp
 *
 * This library provides an implementation of a trajectory within a temporal topology.
 * The path is represented by a set of waypoints and all intermediate points
 * are computed with a linear interpolation.
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

#ifndef REAK_LINEAR_INTERP_HPP
#define REAK_LINEAR_INTERP_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/lin_alg/arithmetic_tuple.hpp>
#include <ReaK/math/lin_alg/mat_num_exceptions.hpp>

#include <ReaK/topologies/spaces/tangent_bundle_concept.hpp>
#include "spatial_trajectory_concept.hpp"

#include "generic_interpolator_factory.hpp"
#include "interpolated_trajectory.hpp"

#include <boost/concept_check.hpp>
#include <cmath>

#include <limits>
#include <list>
#include <map>
#include <type_traits>

namespace ReaK {

namespace pp {

/**
 * Use this tag type for some class templates that are parametrized in terms of the interpolation method used overall.
 */
struct linear_interpolation_tag {};

namespace detail {

template <int Order, typename PointType, typename PointDiff0,
          typename DiffSpace, typename TimeSpace>
void linear_interpolate_impl(PointType& result, const PointType& a,
                             const PointDiff0& dp1p0, const DiffSpace& space,
                             const TimeSpace& t_space, double t_factor,
                             double t_normal) {
  if constexpr (Order > 1) {
    linear_interpolate_impl<(Order - 1)>(result, a, dp1p0, space, t_space,
                                         t_factor, t_normal);
    get<Order>(result) = get_space<Order>(space, t_space).origin();
    return;
  }

  get<0>(result) =
      get_space<0>(space, t_space).adjust(get<0>(a), t_normal * dp1p0);

  if constexpr (Order > 0) {
    get<1>(result) = lift_to_space<1>(dp1p0, t_factor, space, t_space);
  }
}

}  // namespace detail

/**
 * This function template computes a linear interpolation between two points in a
 * temporal and zero-differentiable topology.
 * \tparam PointType The point type on the temporal and zero-differentiable topology.
 * \tparam Topology The temporal and zero-differentiable topology type.
 * \param a The starting point of the interpolation.
 * \param b The ending point of the interpolation.
 * \param t The time value at which the interpolated point is sought.
 * \param space The space on which the points reside.
 * \return The interpolated point at time t, between a and b.
 */
template <typename PointType, typename Topology>
PointType linear_interpolate(const PointType& a, const PointType& b, double t,
                             const Topology& space) {
  using SpaceType = typename temporal_space_traits<Topology>::space_topology;
  using TimeSpaceType = typename temporal_space_traits<Topology>::time_topology;

  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((TangentBundleConcept<SpaceType, 0, TimeSpaceType>));

  double t_factor = b.time - a.time;
  if (std::abs(t_factor) < std::numeric_limits<double>::epsilon()) {
    throw singularity_error(
        "Normalizing factor in cubic Hermite spline is zero!");
  }
  double t_normal = (t - a.time) / (b.time - a.time);

  PointType result;
  result.time = t;

  using Space0 =
      typename derived_N_order_space<SpaceType, TimeSpaceType, 0>::type;

  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));

  using PointDiff0 = topology_point_difference_type_t<Space0>;

  PointDiff0 dp1p0 =
      get_space<0>(space.get_space_topology(), space.get_time_topology())
          .difference(get<0>(b.pt), get<0>(a.pt));

  detail::linear_interpolate_impl<
      max_derivation_order_v<SpaceType, TimeSpaceType>>(
      result.pt, a.pt, dp1p0, space.get_space_topology(),
      space.get_time_topology(), t_factor, t_normal);

  return result;
}

/**
 * This functor class implements a linear interpolation in a temporal and zero-differentiable
 * topology.
 * \tparam SpaceType The topology on which the interpolation is done, should model MetricSpaceConcept and
 * DifferentiableSpaceConcept zero times against time.
 * \tparam TimeSpaceType The time topology.
 */

template <typename SpaceType, typename TimeSpaceType>
class linear_interpolator {
 public:
  using self = linear_interpolator<SpaceType, TimeSpaceType>;
  using point_type = topology_point_type_t<SpaceType>;

  using Space0 =
      typename derived_N_order_space<SpaceType, TimeSpaceType, 0>::type;
  using PointType0 = topology_point_type_t<Space0>;
  using PointDiff0 = topology_point_difference_type_t<Space0>;

  BOOST_CONCEPT_ASSERT((TopologyConcept<SpaceType>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));
  BOOST_CONCEPT_ASSERT((TangentBundleConcept<SpaceType, 0, TimeSpaceType>));

 private:
  PointDiff0 delta_first_order;

 public:
  /**
   * Default constructor.
   */
  linear_interpolator() = default;

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
  linear_interpolator(const point_type& start_point,
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
   * \param space The metric space on which the interpolation resides.
   * \param t_space The time-space against which the interpolation is done.
   * \param factory The factory object that stores relevant fly-weight parameters for the interpolator.
   */
  template <typename Factory>
  void initialize(const point_type& start_point, const point_type& end_point,
                  double /*unused*/, const SpaceType& space,
                  const TimeSpaceType& t_space, const Factory& factory) {
    delta_first_order = get_space<0>(space, t_space)
                            .difference(get<0>(end_point), get<0>(start_point));
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
    if (std::abs(dt_total) < std::numeric_limits<double>::epsilon()) {
      throw singularity_error(
          "Normalizing factor in linear interpolation is zero!");
    }
    double t_normal = dt / dt_total;

    detail::linear_interpolate_impl<
        max_derivation_order_v<SpaceType, TimeSpaceType>>(
        result, start_point, delta_first_order, space, t_space, dt_total,
        t_normal);
  }

  /**
   * Returns the minimum travel time between the initialized start and end points.
   * \return The minimum travel time between the initialized start and end points.
   */
  double get_minimum_travel_time() const { return 0.0; }
};

/**
 * This class is a factory class for linear interpolators on a temporal differentiable space.
 * \tparam TemporalTopology The temporal topology on which the interpolation is done, should model TemporalSpaceConcept.
 */
template <typename TemporalTopology>
class linear_interpolator_factory : public serializable {
 public:
  using self = linear_interpolator_factory<TemporalTopology>;
  using topology = TemporalTopology;
  using point_type = topology_point_type_t<TemporalTopology>;
  using interpolator_type = generic_interpolator<self, linear_interpolator>;

  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<topology>));

 private:
  std::shared_ptr<topology> space;

 public:
  explicit linear_interpolator_factory(const std::shared_ptr<topology>& aSpace)
      : space(aSpace) {}

  linear_interpolator_factory()
      : linear_interpolator_factory(std::shared_ptr<topology>()) {}

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

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2430001, 1,
                              "linear_interpolator_factory", serializable)
};

template <typename SpaceType, typename TimeTopology>
struct get_tagged_spatial_interpolator<linear_interpolation_tag, SpaceType,
                                       TimeTopology> {
  using type = detail::generic_interpolator_impl<linear_interpolator, SpaceType,
                                                 TimeTopology>;
  using pseudo_factory_type = void*;
};

template <typename TemporalSpaceType>
struct get_tagged_temporal_interpolator<linear_interpolation_tag,
                                        TemporalSpaceType> {
  using type =
      generic_interpolator<linear_interpolator_factory<TemporalSpaceType>,
                           linear_interpolator>;
};

/**
 * This class implements a trajectory in a temporal and zero-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points
 * are computed with a linear interpolation. This class models the SpatialTrajectoryConcept.
 * \tparam Topology The topology type on which the points and the path can reside, should model the TemporalSpaceConcept
 * and the DifferentiableSpaceConcept (order 1 with space against time).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the
 * DistanceMetricConcept.
 */
template <typename Topology,
          typename DistanceMetric =
              typename metric_space_traits<Topology>::distance_metric_type>
class linear_interp_traj
    : public interpolated_trajectory<
          Topology, linear_interpolator_factory<Topology>, DistanceMetric> {
 public:
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT(
      (TangentBundleConcept<
          typename temporal_space_traits<Topology>::space_topology, 1,
          typename temporal_space_traits<Topology>::time_topology>));

  using self = linear_interp_traj<Topology, DistanceMetric>;
  using base_class_type =
      interpolated_trajectory<Topology, linear_interpolator_factory<Topology>,
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
  explicit linear_interp_traj(const std::shared_ptr<topology>& aSpace,
                              const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aDist,
                        linear_interpolator_factory<Topology>(aSpace)) {}

  linear_interp_traj() : linear_interp_traj(std::make_shared<topology>()) {}

  /**
   * Constructs the path from a space, the start and end points.
   * \param aSpace The space on which the path is.
   * \param aStart The start point of the path.
   * \param aEnd The end-point of the path.
   * \param aDist The distance metric functor that the path should use.
   */
  linear_interp_traj(const std::shared_ptr<topology>& aSpace,
                     const point_type& aStart, const point_type& aEnd,
                     const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aStart, aEnd, aDist,
                        linear_interpolator_factory<Topology>(aSpace)) {}

  /**
   * Constructs the path from a range of points and their space.
   * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
   * \param aBegin An iterator to the first point of the path.
   * \param aEnd An iterator to the second point of the path.
   * \param aSpace The space on which the path is.
   * \param aDist The distance metric functor that the path should use.
   */
  template <typename ForwardIter>
  linear_interp_traj(ForwardIter aBegin, ForwardIter aEnd,
                     const std::shared_ptr<topology>& aSpace,
                     const distance_metric& aDist = distance_metric())
      : base_class_type(aBegin, aEnd, aSpace, aDist,
                        linear_interpolator_factory<Topology>(aSpace)) {}

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

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2440003, 1, "linear_interp_traj",
                              base_class_type)
};

}  // namespace pp

namespace rtti {

template <>
struct get_type_id<pp::linear_interpolation_tag> {
  static constexpr unsigned int ID = 1;
  static constexpr auto type_name =
      std::string_view{"linear_interpolation_tag"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

}  // namespace rtti
}  // namespace ReaK

#endif

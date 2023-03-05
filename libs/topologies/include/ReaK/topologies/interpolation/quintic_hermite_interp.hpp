/**
 * \file quintic_hermite_interp.hpp
 *
 * This library provides an implementation of a trajectory within a temporal and twice-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points
 * are computed with a quintic Hermite interpolation (quintic hermite spline, or qspline).
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

#ifndef REAK_QUINTIC_HERMITE_INTERP_HPP
#define REAK_QUINTIC_HERMITE_INTERP_HPP

#include "ReaK/core/base/defs.hpp"
#include "ReaK/math/lin_alg/arithmetic_tuple.hpp"
#include "ReaK/math/lin_alg/mat_num_exceptions.hpp"

#include "ReaK/topologies/interpolation/spatial_trajectory_concept.hpp"
#include "ReaK/topologies/spaces/tangent_bundle_concept.hpp"

#include "ReaK/topologies/interpolation/generic_interpolator_factory.hpp"
#include "ReaK/topologies/interpolation/interpolated_trajectory.hpp"

#include <cmath>
#include "boost/concept_check.hpp"

#include <limits>
#include <list>
#include <map>
#include <type_traits>

namespace ReaK {

namespace pp {

/**
 * Use this tag type for some class templates that are parametrized in terms of the interpolation method used overall.
 */
struct quintic_hermite_interpolation_tag {};

namespace detail {

template <int Order, typename PointType, typename PointDiff0,
          typename PointDiff1, typename PointDiff2, typename DiffSpace,
          typename TimeSpace>
void quintic_hermite_interpolate_impl(
    PointType& result, const PointType& a, const PointType& b,
    const PointDiff0& dp1p0, const PointDiff1& dv1v0,
    const PointDiff1& d_ldp1p0_v0, const PointDiff1& d_v1_ldp1p0,
    const PointDiff1& i_a0, const PointDiff1& i_a1, const PointDiff2& da1a0,
    const PointDiff2& da_term1, const PointDiff2& da_term2,
    const DiffSpace& space, const TimeSpace& t_space, double t_factor,
    double t_normal) {
  static_assert(Order >= 2);

  if constexpr (Order > 5) {
    quintic_hermite_interpolate_impl<Order - 1>(
        result, a, b, dp1p0, dv1v0, d_ldp1p0_v0, d_v1_ldp1p0, i_a0, i_a1, da1a0,
        da_term1, da_term2, space, t_space, t_factor, t_normal);
    get<Order>(result) = get_space<Order>(space, t_space).origin();
    return;
  }

  double t2 = t_normal * t_normal;
  double t3 = t_normal * t2;
  double t4 = t2 * t2;
  double t5 = t4 * t_normal;

  get<0>(result) =
      get_space<0>(space, t_space)
          .adjust(
              get<0>(a),
              (10.0 * t3 - 15.0 * t4 + 6.0 * t5) *
                      (dp1p0 - 0.5 * (descend_to_space<0>(get<1>(a), t_factor,
                                                          space, t_space) +
                                      descend_to_space<0>(get<1>(b), t_factor,
                                                          space, t_space))) +
                  t_normal *
                      descend_to_space<0>(
                          get_space<1>(space, t_space)
                              .adjust(get<1>(a),
                                      (t2 - 0.5 * t3) * dv1v0 +
                                          (0.5 * (t_normal - 3.0 * t2 +
                                                  3.0 * t3 - t4)) *
                                              i_a0 +
                                          (0.5 * (t2 - 2.0 * t3 + t4)) * i_a1),
                          t_factor, space, t_space));

  get<1>(result) =
      get_space<1>(space, t_space)
          .adjust(get<1>(a),
                  (15.0 * (t2 - 2.0 * t3 + t4)) * (d_ldp1p0_v0 - d_v1_ldp1p0) +
                      (3.0 * t2 - 2.0 * t3) * dv1v0 +
                      (t_normal - 4.5 * t2 + 6.0 * t3 - 2.5 * t4) * i_a0 +
                      (1.5 * t2 - 4.0 * t3 + 2.5 * t4) * i_a1);

  get<2>(result) =
      get_space<2>(space, t_space)
          .adjust(get<2>(a),
                  (5.0 * (t_normal - 3.0 * t2 + 2.0 * t3)) * da_term1 +
                      (3.0 * (t_normal - t2)) * da_term2 + t_normal * da1a0);

  if constexpr (Order > 2) {
    get<3>(result) = lift_to_space<3>(
        (5.0 - (30.0 - 30.0 * t_normal) * t_normal) * da_term1 +
            (3.0 - 6.0 * t_normal) * da_term2 + da1a0,
        t_factor, space, t_space);
  }

  if constexpr (Order > 3) {
    get<4>(result) = lift_to_space<4>(
        get_space<3>(space, t_space)
            .difference(
                lift_to_space<3>((60.0 * t_normal - 30.0) * da_term1, t_factor,
                                 space, t_space),
                lift_to_space<3>(6.0 * da_term2, t_factor, space, t_space)),
        t_factor, space, t_space);
  }

  if constexpr (Order > 4) {
    get<5>(result) = lift_to_space<5>(
        get_space<4>(space, t_space)
            .difference(
                lift_to_space<4>(
                    60.0 *
                        get_space<3>(space, t_space)
                            .difference(lift_to_space<3>(da_term1, t_factor,
                                                         space, t_space),
                                        get_space<3>(space, t_space).origin()),
                    t_factor, space, t_space),
                get_space<4>(space, t_space).origin()),
        t_factor, space, t_space);
  }
}

}  // namespace detail

/**
 * This function template computes a quintic Hermite interpolation between two points in a
 * temporal and twice-differentiable topology.
 * \tparam PointType The point type on the temporal and twice-differentiable topology.
 * \tparam Topology The temporal and twice-differentiable topology type.
 * \param a The starting point of the interpolation.
 * \param b The ending point of the interpolation.
 * \param t The time value at which the interpolated point is sought.
 * \param space The space on which the points reside.
 * \return The interpolated point at time t, between a and b.
 */
template <typename PointType, typename Topology>
PointType quintic_hermite_interpolate(const PointType& a, const PointType& b,
                                      double t, const Topology& space) {
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  using SpaceType = typename temporal_space_traits<Topology>::space_topology;
  using TimeSpaceType = typename temporal_space_traits<Topology>::time_topology;
  BOOST_CONCEPT_ASSERT((TangentBundleConcept<SpaceType, 2, TimeSpaceType>));

  double t_factor = b.time - a.time;
  if (std::abs(t_factor) < std::numeric_limits<double>::epsilon()) {
    throw singularity_error(
        "Normalizing factor in cubic Hermite spline is zero!");
  }
  double t_normal = (t - a.time) / (b.time - a.time);

  PointType result;
  result.time = t;

  using Space0 = derived_N_order_space_t<SpaceType, TimeSpaceType, 0>;
  using Space1 = derived_N_order_space_t<SpaceType, TimeSpaceType, 1>;
  using Space2 = derived_N_order_space_t<SpaceType, TimeSpaceType, 2>;

  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space1>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space2>));

  using PointType1 = topology_point_type_t<Space1>;
  using PointType2 = topology_point_type_t<Space2>;

  using PointDiff0 = topology_point_difference_type_t<Space0>;
  using PointDiff1 = topology_point_difference_type_t<Space1>;
  using PointDiff2 = topology_point_difference_type_t<Space2>;

  PointDiff0 dp1p0 =
      get_space<0>(space.get_space_topology(), space.get_time_topology())
          .difference(get<0>(b.pt), get<0>(a.pt));
  PointDiff1 dv1v0 =
      get_space<1>(space.get_space_topology(), space.get_time_topology())
          .difference(get<1>(b.pt), get<1>(a.pt));
  PointType1 ldp1p0 = lift_to_space<1>(
      dp1p0, t_factor, space.get_space_topology(), space.get_time_topology());
  PointDiff1 d_ldp1p0_v0 =
      get_space<1>(space.get_space_topology(), space.get_time_topology())
          .difference(ldp1p0, get<1>(a.pt));
  PointDiff1 d_v1_ldp1p0 =
      get_space<1>(space.get_space_topology(), space.get_time_topology())
          .difference(get<1>(b.pt), ldp1p0);
  PointDiff1 i_a0 =
      descend_to_space<1>(get<2>(a.pt), t_factor, space.get_space_topology(),
                          space.get_time_topology());
  PointDiff1 i_a1 =
      descend_to_space<1>(get<2>(b.pt), t_factor, space.get_space_topology(),
                          space.get_time_topology());

  PointType2 l6d_ldp1p0_v0 =
      lift_to_space<2>(6.0 * d_ldp1p0_v0, t_factor, space.get_space_topology(),
                       space.get_time_topology());
  PointType2 l6d_v1_ldp1p0 =
      lift_to_space<2>(6.0 * d_v1_ldp1p0, t_factor, space.get_space_topology(),
                       space.get_time_topology());
  PointDiff2 d_l6d_ldp1p0_v0_a0 =
      get_space<2>(space.get_space_topology(), space.get_time_topology())
          .difference(l6d_ldp1p0_v0, get<2>(a.pt));
  PointDiff2 d_a1_l6d_v1_ldp1p0 =
      get_space<2>(space.get_space_topology(), space.get_time_topology())
          .difference(get<2>(b.pt), l6d_v1_ldp1p0);
  PointType2 ldv1v0 = lift_to_space<2>(
      dv1v0, t_factor, space.get_space_topology(), space.get_time_topology());
  PointDiff2 d_ldv1v0_a0 =
      get_space<2>(space.get_space_topology(), space.get_time_topology())
          .difference(ldv1v0, get<2>(a.pt));
  PointDiff2 d_a1_ldv1v0 =
      get_space<2>(space.get_space_topology(), space.get_time_topology())
          .difference(get<2>(b.pt), ldv1v0);
  PointDiff2 da1a0 =
      get_space<2>(space.get_space_topology(), space.get_time_topology())
          .difference(get<2>(b.pt), get<2>(a.pt));

  PointDiff2 da_term1 = d_l6d_ldp1p0_v0_a0 + d_a1_l6d_v1_ldp1p0;
  PointDiff2 da_term2 = d_ldv1v0_a0 - d_a1_ldv1v0;

  detail::quintic_hermite_interpolate_impl<
      max_derivation_order_v<SpaceType, TimeSpaceType>>(
      result.pt, a.pt, b.pt, dp1p0, dv1v0, d_ldp1p0_v0, d_v1_ldp1p0, i_a0, i_a1,
      da1a0, da_term1, da_term2, space.get_space_topology(),
      space.get_time_topology(), t_factor, t_normal);

  return result;
}

/**
 * This functor class implements a quintic Hermite interpolation in a temporal and twice-differentiable
 * topology.
 * \tparam SpaceType The topology on which the interpolation is done, should model MetricSpaceConcept and
 * DifferentiableSpaceConcept once against time.
 * \tparam TimeSpaceType The time topology.
 */
template <typename SpaceType, typename TimeSpaceType>
class quintic_hermite_interpolator {
 public:
  using self = quintic_hermite_interpolator<SpaceType, TimeSpaceType>;
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

  BOOST_CONCEPT_ASSERT((TopologyConcept<SpaceType>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space0>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space1>));
  BOOST_CONCEPT_ASSERT((LieGroupConcept<Space2>));
  BOOST_CONCEPT_ASSERT((TangentBundleConcept<SpaceType, 2, TimeSpaceType>));

 private:
  PointDiff0 delta_first_order;
  PointDiff1 delta_second_order;
  PointDiff1 delta_lifted_first_and_second;
  PointDiff1 delta_second_and_lifted_first;
  PointDiff1 delta_integral_third_0;
  PointDiff1 delta_integral_third_1;
  PointDiff2 delta_third_order;
  PointDiff2 third_order_term1;
  PointDiff2 third_order_term2;

 public:
  /**
   * Default constructor.
   */
  quintic_hermite_interpolator() = default;

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
  quintic_hermite_interpolator(const point_type& start_point,
                               const point_type& end_point, double dt,
                               const SpaceType& space,
                               const TimeSpaceType& t_space,
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
    delta_first_order = get_space<0>(space, t_space)
                            .difference(get<0>(end_point), get<0>(start_point));
    delta_second_order =
        get_space<1>(space, t_space)
            .difference(get<1>(end_point), get<1>(start_point));

    PointType1 ldp1p0 = lift_to_space<1>(delta_first_order, dt, space, t_space);
    delta_lifted_first_and_second =
        get_space<1>(space, t_space).difference(ldp1p0, get<1>(start_point));
    delta_second_and_lifted_first =
        get_space<1>(space, t_space).difference(get<1>(end_point), ldp1p0);

    delta_integral_third_0 =
        descend_to_space<1>(get<2>(start_point), dt, space, t_space);
    delta_integral_third_1 =
        descend_to_space<1>(get<2>(end_point), dt, space, t_space);

    PointType2 l6d_ldp1p0_v0 = lift_to_space<2>(
        6.0 * delta_lifted_first_and_second, dt, space, t_space);
    PointType2 l6d_v1_ldp1p0 = lift_to_space<2>(
        6.0 * delta_second_and_lifted_first, dt, space, t_space);
    PointType2 ldv1v0 =
        lift_to_space<2>(delta_second_order, dt, space, t_space);

    delta_third_order = get_space<2>(space, t_space)
                            .difference(get<2>(end_point), get<2>(start_point));
    third_order_term1 = get_space<2>(space, t_space)
                            .difference(l6d_ldp1p0_v0, get<2>(start_point)) +
                        get_space<2>(space, t_space)
                            .difference(get<2>(end_point), l6d_v1_ldp1p0);
    third_order_term2 =
        get_space<2>(space, t_space).difference(ldv1v0, get<2>(start_point)) -
        get_space<2>(space, t_space).difference(get<2>(end_point), ldv1v0);
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
          "Normalizing factor in quintic Hermite spline is zero!");
    }
    double t_normal = dt / dt_total;

    detail::quintic_hermite_interpolate_impl<
        max_derivation_order_v<SpaceType, TimeSpaceType>>(
        result, start_point, end_point, delta_first_order, delta_second_order,
        delta_lifted_first_and_second, delta_second_and_lifted_first,
        delta_integral_third_0, delta_integral_third_1, delta_third_order,
        third_order_term1, third_order_term2, space, t_space, dt_total,
        t_normal);
  }

  /**
   * Returns the minimum travel time between the initialized start and end points.
   * \return The minimum travel time between the initialized start and end points.
   */
  double get_minimum_travel_time() const { return 0.0; }
};

/**
 * This class is a factory class for quintic Hermite interpolators on a temporal differentiable space.
 * \tparam TemporalTopology The temporal topology on which the interpolation is done, should model TemporalSpaceConcept.
 */
template <typename TemporalTopology>
class quintic_hermite_interp_factory : public serializable {
 public:
  using self = quintic_hermite_interp_factory<TemporalTopology>;
  using topology = TemporalTopology;
  using point_type = topology_point_type_t<TemporalTopology>;
  using interpolator_type =
      generic_interpolator<self, quintic_hermite_interpolator>;

  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<TemporalTopology>));

 private:
  std::shared_ptr<topology> space;

 public:
  explicit quintic_hermite_interp_factory(
      const std::shared_ptr<topology>& aSpace)
      : space(aSpace) {}

  quintic_hermite_interp_factory()
      : quintic_hermite_interp_factory(std::shared_ptr<topology>()) {}

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

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2430003, 1,
                              "quintic_hermite_interp_factory", serializable)
};

template <typename SpaceType, typename TimeTopology>
struct get_tagged_spatial_interpolator<quintic_hermite_interpolation_tag,
                                       SpaceType, TimeTopology> {
  using type = detail::generic_interpolator_impl<quintic_hermite_interpolator,
                                                 SpaceType, TimeTopology>;
  using pseudo_factory_type = void*;
};

template <typename TemporalSpaceType>
struct get_tagged_temporal_interpolator<quintic_hermite_interpolation_tag,
                                        TemporalSpaceType> {
  using type =
      generic_interpolator<quintic_hermite_interp_factory<TemporalSpaceType>,
                           quintic_hermite_interpolator>;
};

/**
 * This class implements a trajectory in a temporal and twice-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points
 * are computed with a quintic Hermite interpolation. This class models the SpatialTrajectoryConcept.
 * \tparam Topology The topology type on which the points and the path can reside, should model the TemporalSpaceConcept
 * and the DifferentiableSpaceConcept (order 1 with space against time).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the
 * DistanceMetricConcept.
 */
template <typename Topology,
          typename DistanceMetric =
              typename metric_space_traits<Topology>::distance_metric_type>
class quintic_hermite_interp_traj
    : public interpolated_trajectory<
          Topology, quintic_hermite_interp_factory<Topology>, DistanceMetric> {
 public:
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT(
      (TangentBundleConcept<
          typename temporal_space_traits<Topology>::space_topology, 2,
          typename temporal_space_traits<Topology>::time_topology>));

  using self = quintic_hermite_interp_traj<Topology, DistanceMetric>;
  using base_class_type = interpolated_trajectory<
      Topology, quintic_hermite_interp_factory<Topology>, DistanceMetric>;

  using topology = typename base_class_type::topology;
  using distance_metric = typename base_class_type::distance_metric;
  using point_type = typename base_class_type::point_type;

  /**
   * Constructs the path from a space, assumes the start and end are at the origin
   * of the space.
   * \param aSpace The space on which the path is.
   * \param aDist The distance metric functor that the path should use.
   */
  explicit quintic_hermite_interp_traj(
      const std::shared_ptr<topology>& aSpace = std::make_shared<topology>(),
      const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aDist,
                        quintic_hermite_interp_factory<Topology>(aSpace)) {}

  /**
   * Constructs the path from a space, the start and end points.
   * \param aSpace The space on which the path is.
   * \param aStart The start point of the path.
   * \param aEnd The end-point of the path.
   * \param aDist The distance metric functor that the path should use.
   */
  quintic_hermite_interp_traj(const std::shared_ptr<topology>& aSpace,
                              const point_type& aStart, const point_type& aEnd,
                              const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aStart, aEnd, aDist,
                        quintic_hermite_interp_factory<Topology>(aSpace)) {}

  /**
   * Constructs the path from a range of points and their space.
   * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
   * \param aBegin An iterator to the first point of the path.
   * \param aEnd An iterator to the second point of the path.
   * \param aSpace The space on which the path is.
   * \param aDist The distance metric functor that the path should use.
   */
  template <typename ForwardIter>
  quintic_hermite_interp_traj(ForwardIter aBegin, ForwardIter aEnd,
                              const std::shared_ptr<topology>& aSpace,
                              const distance_metric& aDist = distance_metric())
      : base_class_type(aBegin, aEnd, aSpace, aDist,
                        quintic_hermite_interp_factory<Topology>(aSpace)) {}

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

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2440005, 1,
                              "quintic_hermite_interp_traj", base_class_type)
};

}  // namespace pp

namespace rtti {

template <>
struct get_type_id<pp::quintic_hermite_interpolation_tag> {
  static constexpr unsigned int ID = 3;
  static constexpr auto type_name =
      std::string_view{"quintic_hermite_interpolation_tag"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }
};

}  // namespace rtti
}  // namespace ReaK

#endif

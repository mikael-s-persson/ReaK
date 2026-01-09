/**
 * \file cubic_hermite_interp.h
 *
 * This library provides an implementation of a trajectory within a temporal and once-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points
 * are computed with a cubic Hermite interpolation (cubic hermite spline, or cspline).
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_CUBIC_HERMITE_INTERP_H_
#define REAK_TOPOLOGIES_INTERPOLATION_CUBIC_HERMITE_INTERP_H_

#include "ReaK/math/lin_alg/arithmetic_tuple.h"
#include "ReaK/math/lin_alg/mat_num_exceptions.h"

#include "ReaK/topologies/interpolation/spatial_trajectory_concept.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/tangent_bundle_concept.h"

#include "ReaK/topologies/interpolation/interpolated_trajectory.h"
#include "ReaK/topologies/spaces/generic_interpolator_factory.h"

#include "ReaK/topologies/spaces/temporal_space_concept.h"

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
struct cubic_hermite_interpolation_tag {};

namespace detail {

template <int Order, typename PointType, typename PointDiff0,
          typename PointDiff1, typename DiffSpace, typename TimeSpace>
inline void cubic_hermite_interpolate_impl(
    PointType& result, const PointType& a, const PointType& b,
    const PointDiff0& dp1p0, const PointDiff1& dv1v0,
    const PointDiff1 d_ldp1p0_v0, const DiffSpace& space,
    const TimeSpace& t_space, double t_factor, double t_normal) {
  static_assert(Order >= 1);

  if constexpr (Order > 3) {
    cubic_hermite_interpolate_impl<Order - 1>(result, a, b, dp1p0, dv1v0,
                                              d_ldp1p0_v0, space, t_space,
                                              t_factor, t_normal);
    get<Order>(result) = get_space<Order>(space, t_space).origin();
    return;
  }

  double t2 = t_normal * t_normal;
  double t3 = t_normal * t2;

  get<0>(result) =
      get_space<0>(space, t_space)
          .adjust(
              get<0>(a),
              (3.0 * t2 - 2.0 * t3) * dp1p0 +
                  (t_normal - t2 * 2.0 + t3) *
                      descend_to_space<0>(get<1>(a), t_factor, space, t_space) +
                  (t3 - t2) *
                      descend_to_space<0>(get<1>(b), t_factor, space, t_space));

  get<1>(result) =
      get_space<1>(space, t_space)
          .adjust(get<1>(a), ((t_normal - t2) * 6.0) * d_ldp1p0_v0 -
                                 (2.0 * t_normal - 3.0 * t2) * dv1v0);

  if constexpr (Order > 1) {
    get<2>(result) =
        get_space<2>(space, t_space)
            .adjust(lift_to_space<2>(dv1v0, t_factor, space, t_space),
                    (6.0 - 12.0 * t_normal) *
                        get_space<2>(space, t_space)
                            .difference(lift_to_space<2>(d_ldp1p0_v0, t_factor,
                                                         space, t_space),
                                        lift_to_space<2>(0.5 * dv1v0, t_factor,
                                                         space, t_space)));
  }

  if constexpr (Order > 2) {
    get<3>(result) = lift_to_space<3>(
        -12.0 *
            get_space<2>(space, t_space)
                .difference(
                    lift_to_space<2>(d_ldp1p0_v0, t_factor, space, t_space),
                    lift_to_space<2>(0.5 * dv1v0, t_factor, space, t_space)),
        t_factor, space, t_space);
  }
}

}  // namespace detail

/**
 * This function template computes a cubic Hermite interpolation between two points in a
 * temporal and once-differentiable topology.
 * \param a The starting point of the interpolation.
 * \param b The ending point of the interpolation.
 * \param t The time value at which the interpolated point is sought.
 * \param space The space on which the points reside.
 * \return The interpolated point at time t, between a and b.
 */
template <typename PointType, TemporalSpace Space>
PointType cubic_hermite_interpolate(const PointType& a, const PointType& b,
                                    double t, const Space& space) {
  using SpaceType = typename temporal_space_traits<Space>::space_topology;
  using TimeSpaceType = typename temporal_space_traits<Space>::time_topology;
  static_assert(TangentBundle<SpaceType, TimeSpaceType, 0, 1>);

  double t_factor = b.time - a.time;
  if (std::abs(t_factor) < std::numeric_limits<double>::epsilon()) {
    throw singularity_error(
        "Normalizing factor in cubic Hermite spline is zero!");
  }
  double t_normal = (t - a.time) / (b.time - a.time);

  PointType result;
  result.time = t;

  using Space0 = derived_N_order_space_t<SpaceType, TimeSpaceType, 0>;
  using Space1 = derived_N_order_space<SpaceType, TimeSpaceType, 1>;

  static_assert(LieGroup<Space0>);
  static_assert(LieGroup<Space1>);

  using PointDiff0 = topology_point_difference_type_t<Space0>;
  using PointDiff1 = topology_point_difference_type_t<Space1>;

  PointDiff0 dp1p0 =
      get_space<0>(space.get_space_topology(), space.get_time_topology())
          .difference(get<0>(b.pt), get<0>(a.pt));
  PointDiff1 dv1v0 =
      get_space<1>(space.get_space_topology(), space.get_time_topology())
          .difference(get<1>(b.pt), get<1>(a.pt));
  PointDiff1 d_ldp1p0_v0 =
      get_space<1>(space.get_space_topology(), space.get_time_topology())
          .difference(
              lift_to_space<1>(dp1p0, t_factor, space.get_space_topology(),
                               space.get_time_topology()),
              get<1>(a.pt));

  detail::cubic_hermite_interpolate_impl<
      max_derivation_order_v<SpaceType, TimeSpaceType>>(
      result.pt, a.pt, b.pt, dp1p0, dv1v0, d_ldp1p0_v0,
      space.get_space_topology(), space.get_time_topology(), t_factor,
      t_normal);

  return result;
}

/**
 * This functor class implements a cubic Hermite interpolation in a temporal and once-differentiable
 * topology.
 * \tparam SpaceType The topology on which the interpolation is done.
 * \tparam TimeSpaceType The time topology.
 */
template <MetricSpace SpaceType, Topology TimeSpaceType>
class cubic_hermite_interpolator {
 public:
  using self = cubic_hermite_interpolator<SpaceType, TimeSpaceType>;
  using point_type = topology_point_type_t<SpaceType>;

  using Space0 = derived_N_order_space_t<SpaceType, TimeSpaceType, 0>;
  using PointType0 = topology_point_type_t<Space0>;
  using PointDiff0 = topology_point_difference_type_t<Space0>;
  using Space1 = derived_N_order_space_t<SpaceType, TimeSpaceType, 1>;
  using PointType1 = topology_point_type_t<Space1>;
  using PointDiff1 = topology_point_difference_type_t<Space1>;

  static_assert(LieGroup<Space0>);
  static_assert(LieGroup<Space1>);
  static_assert(TangentBundle<SpaceType, TimeSpaceType, 0, 1>);

 private:
  PointDiff0 delta_first_order;
  PointDiff1 delta_second_order;
  PointDiff1 delta_lifted_first_and_second;

 public:
  /**
   * Default constructor.
   */
  cubic_hermite_interpolator() = default;

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
  cubic_hermite_interpolator(const point_type& start_point,
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
    delta_lifted_first_and_second =
        get_space<1>(space, t_space)
            .difference(lift_to_space<1>(delta_first_order, dt, space, t_space),
                        get<1>(start_point));
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
          "Normalizing factor in cubic Hermite spline is zero!");
    }
    double t_normal = dt / dt_total;

    detail::cubic_hermite_interpolate_impl<
        max_derivation_order_v<SpaceType, TimeSpaceType>>(
        result, start_point, end_point, delta_first_order, delta_second_order,
        delta_lifted_first_and_second, space, t_space, dt_total, t_normal);
  }

  /**
   * Returns the minimum travel time between the initialized start and end points.
   * \return The minimum travel time between the initialized start and end points.
   */
  double get_minimum_travel_time() const { return 0.0; }
};

/**
 * This class is a factory class for cubic Hermite interpolators on a temporal differentiable space.
 */
template <TemporalSpace Space>
class cubic_hermite_interp_factory : public serializable {
 public:
  using self = cubic_hermite_interp_factory<Space>;
  using topology = Space;
  using point_type = topology_point_type_t<Space>;
  using interpolator_type =
      generic_interpolator<self, cubic_hermite_interpolator>;

 private:
  std::shared_ptr<topology> space;

 public:
  explicit cubic_hermite_interp_factory(
      const std::shared_ptr<topology>& aSpace = {})
      : space(aSpace) {}

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

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2430002, 1,
                              "cubic_hermite_interp_factory", serializable)
};

template <typename SpaceType, typename TimeTopology>
struct get_tagged_spatial_interpolator<cubic_hermite_interpolation_tag,
                                       SpaceType, TimeTopology> {
  using type = detail::generic_interpolator_impl<cubic_hermite_interpolator,
                                                 SpaceType, TimeTopology>;
  using pseudo_factory_type = void*;
};

template <typename TemporalSpaceType>
struct get_tagged_temporal_interpolator<cubic_hermite_interpolation_tag,
                                        TemporalSpaceType> {
  using type =
      generic_interpolator<cubic_hermite_interp_factory<TemporalSpaceType>,
                           cubic_hermite_interpolator>;
};

/**
 * This class implements a trajectory in a temporal and once-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points
 * are computed with a cubic Hermite interpolation. This class models the SpatialTrajectoryConcept.
 * \tparam Space The topology type on which the points and the path can reside.
 * \tparam Metric The distance metric used to assess the distance between points in the path.
 */
template <TemporalSpace Space,
          DistanceMetric<Space> Metric =
              typename metric_space_traits<Space>::distance_metric_type>
class cubic_hermite_interp_traj
    : public interpolated_trajectory<Space, cubic_hermite_interp_factory<Space>,
                                     Metric> {
 public:
  static_assert(TangentBundle<
                typename temporal_space_traits<Space>::space_topology,
                typename temporal_space_traits<Space>::time_topology, 0, 1>);

  using self = cubic_hermite_interp_traj<Space, Metric>;
  using base_class_type =
      interpolated_trajectory<Space, cubic_hermite_interp_factory<Space>,
                              Metric>;

  using topology = typename base_class_type::topology;
  using distance_metric = typename base_class_type::distance_metric;
  using point_type = typename base_class_type::point_type;

  /**
   * Constructs the path from a space, assumes the start and end are at the origin
   * of the space.
   * \param aSpace The space on which the path is.
   * \param aDist The distance metric functor that the path should use.
   */
  explicit cubic_hermite_interp_traj(
      const std::shared_ptr<topology>& aSpace = std::make_shared<topology>(),
      const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aDist,
                        cubic_hermite_interp_factory<Space>(aSpace)) {}

  /**
   * Constructs the path from a space, the start and end points.
   * \param aSpace The space on which the path is.
   * \param aStart The start point of the path.
   * \param aEnd The end-point of the path.
   * \param aDist The distance metric functor that the path should use.
   */
  cubic_hermite_interp_traj(const std::shared_ptr<topology>& aSpace,
                            const point_type& aStart, const point_type& aEnd,
                            const distance_metric& aDist = distance_metric())
      : base_class_type(aSpace, aStart, aEnd, aDist,
                        cubic_hermite_interp_factory<Space>(aSpace)) {}

  /**
   * Constructs the path from a range of points and their space.
   * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
   * \param aBegin An iterator to the first point of the path.
   * \param aEnd An iterator to the second point of the path.
   * \param aSpace The space on which the path is.
   * \param aDist The distance metric functor that the path should use.
   */
  template <typename ForwardIter>
  cubic_hermite_interp_traj(ForwardIter aBegin, ForwardIter aEnd,
                            const std::shared_ptr<topology>& aSpace,
                            const distance_metric& aDist = distance_metric())
      : base_class_type(aBegin, aEnd, aSpace, aDist,
                        cubic_hermite_interp_factory<Space>(aSpace)) {}

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_class_type::save(A,
                          base_class_type::get_static_object_type()->version());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base_class_type::load(A,
                          base_class_type::get_static_object_type()->version());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2440004, 1, "cubic_hermite_interp_traj",
                              base_class_type)
};

}  // namespace pp

namespace rtti {

template <>
struct get_type_id<pp::cubic_hermite_interpolation_tag> {
  static constexpr unsigned int id = 2;
  static constexpr auto type_name =
      std::string_view{"cubic_hermite_interpolation_tag"};
  static construct_ptr create_ptr() noexcept { return nullptr; }
};

}  // namespace rtti

}  // namespace ReaK

#endif  // REAK_TOPOLOGIES_INTERPOLATION_CUBIC_HERMITE_INTERP_H_

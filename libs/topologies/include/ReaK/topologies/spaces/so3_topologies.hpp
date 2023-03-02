/**
 * \file so3_topologies.hpp
 *
 * This library provides classes that define topologies on SO(3) (3D rotation). A quaternion-topology is
 * a simple metric-space where the points are unit quaternion values. Higher-order differential spaces
 * in SO(3) are just normal vector-spaces.
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

#ifndef REAK_SO3_TOPOLOGIES_HPP
#define REAK_SO3_TOPOLOGIES_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/global_rng.hpp>
#include <ReaK/core/base/named_object.hpp>

#include "tuple_distance_metrics.hpp"

#include "differentiable_space.hpp"
#include "hyperball_topology.hpp"
#include "rate_limited_spaces.hpp"

#include <ReaK/math/kinetostatics/quat_alg.hpp>
#include <ReaK/math/lin_alg/mat_num_exceptions.hpp>
#include <ReaK/math/lin_alg/vect_alg.hpp>

#include <random>
#include <type_traits>

namespace ReaK::pp {

/**
 * This class implements a quaternion-topology. Because quaternions are constrained on the unit
 * hyper-sphere, this topology is indeed bounded (yet infinite at the same time). This class
 * models the MetricSpaceConcept, the LieGroupConcept, and the PointDistributionConcept.
 * \tparam T The value-type for the topology.
 */
template <typename T>
class quaternion_topology : public named_object {
 public:
  using self = quaternion_topology<T>;

  using point_type = unit_quat<T>;
  using point_difference_type = vect<T, 3>;

  using distance_metric_type = default_distance_metric;
  using random_sampler_type = default_random_sampler;

  static constexpr std::size_t dimensions = 3;

  /**
   * Default constructor.
   */
  explicit quaternion_topology(const std::string& aName = "quaternion_topology")
      : named_object() {
    setName(aName);
  }

  /*************************************************************************
  *                             MetricSpaceConcept
  * **********************************************************************/

  /**
   * Returns the distance between two points.
   */
  virtual double distance(const point_type& a, const point_type& b) const {
    return ReaK::norm_2(this->difference(b, a));
  }

  /**
   * Returns the norm of the difference between two points.
   */
  virtual double norm(const point_difference_type& delta) const {
    return ReaK::norm_2(delta);
  }

  /*************************************************************************
   *                         for PointDistributionConcept
   * **********************************************************************/

  /**
   * Generates a random point in the space, uniformly distributed.
   */
  point_type random_point() const {
    global_rng_type& rng = get_global_rng();
    std::normal_distribution<vect_value_type_t<point_difference_type>> var_rnd;
    return point_type(var_rnd(rng), var_rnd(rng), var_rnd(rng), var_rnd(rng));
    // According to most sources, normalizing a vector of Normal-distributed components will yield a uniform
    // distribution of the components on the unit hyper-sphere. N.B., the normalization happens in the
    // constructor of the unit-quaternion (when constructed from four values).
  }

  /*************************************************************************
   *                             TopologyConcept
   * **********************************************************************/

  /**
   * Returns the difference between two points (analogous to a - b, but implemented in SO(3) Lie algebra).
   */
  virtual point_difference_type difference(const point_type& a,
                                           const point_type& b) const {
    return T(2.0) * log(conj(b) * a);
  }

  /**
   * Returns the addition of a point-difference to a point.
   */
  virtual point_type adjust(const point_type& a,
                            const point_difference_type& delta) const {
    return a * exp(0.5 * delta);
  }

  /**
   * Returns the origin of the space (the lower-limit).
   */
  point_type origin() const {
    return point_type();  // this just creates the "no-rotation" quaternion.
  }

  /**
   * Tests if a given point is within the boundary of this space.
   */
  bool is_in_bounds(const point_type& a) const { return true; }

  point_difference_type get_diff_to_boundary(
      const point_type& /*unused*/) const {
    return point_difference_type();
  }

  void bring_point_in_bounds(point_type& /*unused*/) const {}

  /*************************************************************************
  *                             LieGroupConcept
  * **********************************************************************/

  /**
   * Returns a point which is at a fraction between two points a to b. This function uses SLERP.
   */
  point_type move_position_toward(const point_type& a, double fraction,
                                  const point_type& b) const {
    return unit(a * pow(conj(a) * b, fraction));
  }

  /**
   * Returns a point which is at a backward fraction between two points a to b. This function uses SLERP.
   */
  point_type move_position_back_to(const point_type& a, double fraction,
                                   const point_type& b) const {
    return unit(a * pow(conj(a) * b, (1.0 - fraction)));
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    ReaK::named_object::save(
        A, named_object::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    ReaK::named_object::load(
        A, named_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240000C, 1, "quaternion_topology",
                              named_object)
};

template <typename T>
struct is_metric_space<quaternion_topology<T>> : std::true_type {};

template <typename T>
struct is_reversible_space<quaternion_topology<T>> : std::true_type {};

template <typename T>
struct is_point_distribution<quaternion_topology<T>> : std::true_type {};

/**
 * This class implements a rate-limited quaternion-space. This topology will produce distance and norm
 * values which are expressed in a time-unit and represent the time needed to travel along a given
 * quaternion difference vector (for the norm function, and the distance function gives the time
 * needed to travel between two given quaternions) assuming that the entire travel is done at the
 * specified maximum angular speed. This class models the MetricSpaceConcept, the LieGroupConcept,
 * and the PointDistributionConcept.
 * \tparam T The value-type for the topology.
 */
template <typename T>
class rate_limited_quat_space : public quaternion_topology<T> {
 public:
  using self = rate_limited_quat_space<T>;
  using base = quaternion_topology<T>;

  using point_type = typename base::point_type;
  using point_difference_type = typename base::point_difference_type;

  using distance_metric_type = default_distance_metric;
  using random_sampler_type = typename base::random_sampler_type;

  static constexpr std::size_t dimensions = 3;

  double max_angular_speed;

  /**
   * Default constructor.
   * \param aName The name of the topology.
   * \param aMaxAngSpeed The maximum angular speed that limits the rate of motion between SO(3) configurations.
   */
  explicit rate_limited_quat_space(
      const std::string& aName = "rate_limited_quat_space",
      double aMaxAngSpeed = 1.0)
      : base(aName), max_angular_speed(aMaxAngSpeed) {
    if (max_angular_speed < std::numeric_limits<double>::epsilon()) {
      throw singularity_error(
          "Maximum angular speed cannot be less-than or equal to 0!");
    }
  }

  /*************************************************************************
  *                             MetricSpaceConcept
  * **********************************************************************/

  /**
   * Returns the distance between two points.
   */
  double distance(const point_type& a, const point_type& b) const override {
    return ReaK::norm_2(this->difference(b, a)) / max_angular_speed;
  }

  /**
   * Returns the norm of the difference between two points.
   */
  double norm(const point_difference_type& delta) const override {
    return ReaK::norm_2(delta) / max_angular_speed;
  }

  /*************************************************************************
   *                             TopologyConcept
   * **********************************************************************/

  /**
   * Returns the difference between two points (analogous to a - b, but implemented in SO(3) Lie algebra).
   */
  point_difference_type difference(const point_type& a,
                                   const point_type& b) const override {
    return T(2.0 / max_angular_speed) * log(conj(b) * a);
  }

  /**
   * Returns the addition of a point-difference to a point.
   */
  point_type adjust(const point_type& a,
                    const point_difference_type& delta) const override {
    return a * exp((0.5 * max_angular_speed) * delta);
  }

  /**
   * Tests if a given point is within the boundary of this space.
   */
  bool is_in_bounds(const point_type& a) const { return true; }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base::save(A, base::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(max_angular_speed);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base::load(A, base::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(max_angular_speed);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240000F, 1, "rate_limited_quat_space",
                              base)
};

template <typename T>
struct is_metric_space<rate_limited_quat_space<T>> : std::true_type {};

template <typename T>
struct is_reversible_space<rate_limited_quat_space<T>> : std::true_type {};

template <typename T>
struct is_point_distribution<rate_limited_quat_space<T>> : std::true_type {};

/**
 * This class implements an angular velocity topology (for SO(3)). The angular velocities are constrained
 * to within a hyper-ball of a given maximum radius (max angular speed), this topology models
 * the MetricSpaceConcept, and is bounded spherically (models BoundedSpaceConcept and SphereBoundedSpaceConcept).
 * \tparam T The value-type for the topology.
 */
template <typename T>
class ang_velocity_3D_topology : public hyperball_topology<vect<T, 3>> {
 public:
  using self = ang_velocity_3D_topology<T>;
  using base = hyperball_topology<vect<T, 3>>;

  using point_type = typename base::point_type;
  using point_difference_type = typename base::point_type;

  using distance_metric_type = typename base::distance_metric_type;
  using random_sampler_type = typename base::random_sampler_type;

  static constexpr std::size_t dimensions = 3;

  /**
   * Default constructor.
   * \param aMaxAngSpeed The maximum (scalar) angular velocity that bounds this hyper-ball topology.
   */
  explicit ang_velocity_3D_topology(
      const std::string& aName = "ang_velocity_3D_topology",
      double aMaxAngSpeed = 1.0)
      : base(aName, vect<T, 3>(), aMaxAngSpeed) {}

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base::save(A, base::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base::load(A, base::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240000D, 1, "ang_velocity_3D_topology",
                              base)
};

template <typename T>
struct is_metric_space<ang_velocity_3D_topology<T>> : std::true_type {};

template <typename T>
struct is_reversible_space<ang_velocity_3D_topology<T>> : std::true_type {};

template <typename T>
struct is_point_distribution<ang_velocity_3D_topology<T>> : std::true_type {};

/**
 * This class implements an angular acceleration topology (for SO(3)). The angular accelerations are constrained
 * to within a hyper-ball of a given maximum radius (max angular acceleration), this topology models
 * the MetricSpaceConcept, and is bounded spherically (models BoundedSpaceConcept and SphereBoundedSpaceConcept).
 * \tparam T The value-type for the topology.
 */
template <typename T>
class ang_accel_3D_topology : public hyperball_topology<vect<T, 3>> {
 public:
  using self = ang_accel_3D_topology<T>;
  using base = hyperball_topology<vect<T, 3>>;

  using point_type = typename base::point_type;
  using point_difference_type = typename base::point_type;

  using distance_metric_type = typename base::distance_metric_type;
  using random_sampler_type = typename base::random_sampler_type;

  static constexpr std::size_t dimensions = 3;

  /**
   * Default constructor.
   * \param aMaxAngSpeed The maximum (scalar) angular acceleration that bounds this hyper-ball topology.
   */
  explicit ang_accel_3D_topology(
      const std::string& aName = "ang_accel_3D_topology",
      double aMaxAngAcc = 1.0)
      : base(aName, vect<T, 3>(), aMaxAngAcc) {}

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base::save(A, base::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base::load(A, base::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240000E, 1, "ang_accel_3D_topology",
                              base)
};

template <typename T>
struct is_metric_space<ang_accel_3D_topology<T>> : std::true_type {};

template <typename T>
struct is_reversible_space<ang_accel_3D_topology<T>> : std::true_type {};

template <typename T>
struct is_point_distribution<ang_accel_3D_topology<T>> : std::true_type {};

/**
 * This meta-function defines the type for a 0th order SO(3) topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct so3_0th_order_topology {
  using type = differentiable_space<
      time_topology, arithmetic_tuple<quaternion_topology<T>>, DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using so3_0th_order_topology_t =
    typename so3_0th_order_topology<T, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 1st order SO(3) topology (a once-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct so3_1st_order_topology {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<quaternion_topology<T>, ang_velocity_3D_topology<T>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using so3_1st_order_topology_t =
    typename so3_1st_order_topology<T, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 2nd order SO(3) topology (a twice-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct so3_2nd_order_topology {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<quaternion_topology<T>, ang_velocity_3D_topology<T>,
                       ang_accel_3D_topology<T>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using so3_2nd_order_topology_t =
    typename so3_2nd_order_topology<T, DistanceMetric>::type;

template <typename T, int Order,
          typename DistanceMetric = euclidean_tuple_distance>
struct so3_topology {
  using type = std::conditional_t<
      (Order == 0), so3_0th_order_topology_t<T, DistanceMetric>,
      std::conditional_t<(Order == 1),
                         so3_1st_order_topology_t<T, DistanceMetric>,
                         so3_2nd_order_topology_t<T, DistanceMetric>>>;
};

template <typename T, int Order,
          typename DistanceMetric = euclidean_tuple_distance>
using so3_topology_t = typename so3_topology<T, Order, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 0th order SO(3) rate-limited topology (a zero-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct so3_0th_order_rl_topology {
  using type =
      reach_time_diff_space<time_topology,
                            arithmetic_tuple<rate_limited_quat_space<T>>,
                            DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using so3_0th_order_rl_topology_t =
    typename so3_0th_order_rl_topology<T, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 1st order SO(3) rate-limited topology (a once-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct so3_1st_order_rl_topology {
  using type = reach_time_diff_space<
      time_topology,
      arithmetic_tuple<rate_limited_quat_space<T>, ang_velocity_3D_topology<T>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using so3_1st_order_rl_topology_t =
    typename so3_1st_order_rl_topology<T, DistanceMetric>::type;

/**
 * This meta-function defines the type for a 2nd order SO(3) rate-limited topology (a twice-differentiable space).
 * \tparam T The value type for the topology.
 * \tparam DistanceMetric The distance metric to apply to the tuple.
 */
template <typename T, typename DistanceMetric = euclidean_tuple_distance>
struct so3_2nd_order_rl_topology {
  using type = reach_time_diff_space<
      time_topology,
      arithmetic_tuple<rate_limited_quat_space<T>, ang_velocity_3D_topology<T>,
                       ang_accel_3D_topology<T>>,
      DistanceMetric>;
};

template <typename T, typename DistanceMetric = euclidean_tuple_distance>
using so3_2nd_order_rl_topology_t =
    typename so3_2nd_order_rl_topology<T, DistanceMetric>::type;

template <typename T, int Order,
          typename DistanceMetric = euclidean_tuple_distance>
struct so3_rl_topology {
  using type = std::conditional_t<
      (Order == 0), so3_0th_order_rl_topology_t<T, DistanceMetric>,
      std::conditional_t<(Order == 1),
                         so3_1st_order_rl_topology_t<T, DistanceMetric>,
                         so3_2nd_order_rl_topology_t<T, DistanceMetric>>>;
};

template <typename T, int Order,
          typename DistanceMetric = euclidean_tuple_distance>
using so3_rl_topology_t =
    typename so3_rl_topology<T, Order, DistanceMetric>::type;

}  // namespace ReaK::pp

#include "reachability_space.hpp"
#include "temporal_space.hpp"
#include "time_poisson_topology.hpp"

namespace ReaK::pp {

extern template class quaternion_topology<double>;
extern template class rate_limited_quat_space<double>;
extern template class ang_velocity_3D_topology<double>;
extern template class ang_accel_3D_topology<double>;

extern template class temporal_space<
    quaternion_topology<double>, time_poisson_topology, spatial_distance_only>;
extern template class temporal_space<rate_limited_quat_space<double>,
                                     time_poisson_topology,
                                     spatial_distance_only>;
extern template class temporal_space<rate_limited_quat_space<double>,
                                     time_poisson_topology,
                                     reach_plus_time_metric>;
extern template class temporal_space<ang_velocity_3D_topology<double>,
                                     time_poisson_topology,
                                     spatial_distance_only>;
extern template class temporal_space<ang_accel_3D_topology<double>,
                                     time_poisson_topology,
                                     spatial_distance_only>;

// so3_0th_order_topology
extern template class differentiable_space<
    time_topology, arithmetic_tuple<quaternion_topology<double>>,
    euclidean_tuple_distance>;
// so3_1st_order_topology
extern template class differentiable_space<
    time_topology,
    arithmetic_tuple<quaternion_topology<double>,
                     ang_velocity_3D_topology<double>>,
    euclidean_tuple_distance>;
// so3_2nd_order_topology
extern template class differentiable_space<
    time_topology,
    arithmetic_tuple<quaternion_topology<double>,
                     ang_velocity_3D_topology<double>,
                     ang_accel_3D_topology<double>>,
    euclidean_tuple_distance>;

// so3_0th_order_topology
extern template class temporal_space<
    differentiable_space<time_topology,
                         arithmetic_tuple<quaternion_topology<double>>,
                         euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;
// so3_1st_order_topology
extern template class temporal_space<
    differentiable_space<time_topology,
                         arithmetic_tuple<quaternion_topology<double>,
                                          ang_velocity_3D_topology<double>>,
                         euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;
// so3_2nd_order_topology
extern template class temporal_space<
    differentiable_space<time_topology,
                         arithmetic_tuple<quaternion_topology<double>,
                                          ang_velocity_3D_topology<double>,
                                          ang_accel_3D_topology<double>>,
                         euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;

// so3_0th_order_rl_topology
extern template class reach_time_diff_space<
    time_topology, arithmetic_tuple<rate_limited_quat_space<double>>,
    euclidean_tuple_distance>;
// so3_1st_order_rl_topology
extern template class reach_time_diff_space<
    time_topology,
    arithmetic_tuple<rate_limited_quat_space<double>,
                     ang_velocity_3D_topology<double>>,
    euclidean_tuple_distance>;
// so3_2nd_order_rl_topology
extern template class reach_time_diff_space<
    time_topology,
    arithmetic_tuple<rate_limited_quat_space<double>,
                     ang_velocity_3D_topology<double>,
                     ang_accel_3D_topology<double>>,
    euclidean_tuple_distance>;

// so3_0th_order_rl_topology
extern template class temporal_space<
    reach_time_diff_space<time_topology,
                          arithmetic_tuple<rate_limited_quat_space<double>>,
                          euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;
// so3_1st_order_rl_topology
extern template class temporal_space<
    reach_time_diff_space<time_topology,
                          arithmetic_tuple<rate_limited_quat_space<double>,
                                           ang_velocity_3D_topology<double>>,
                          euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;
// so3_2nd_order_rl_topology
extern template class temporal_space<
    reach_time_diff_space<time_topology,
                          arithmetic_tuple<rate_limited_quat_space<double>,
                                           ang_velocity_3D_topology<double>,
                                           ang_accel_3D_topology<double>>,
                          euclidean_tuple_distance>,
    time_poisson_topology, spatial_distance_only>;

// so3_0th_order_rl_topology
extern template class temporal_space<
    reach_time_diff_space<time_topology,
                          arithmetic_tuple<rate_limited_quat_space<double>>,
                          euclidean_tuple_distance>,
    time_poisson_topology, reach_plus_time_metric>;
// so3_1st_order_rl_topology
extern template class temporal_space<
    reach_time_diff_space<time_topology,
                          arithmetic_tuple<rate_limited_quat_space<double>,
                                           ang_velocity_3D_topology<double>>,
                          euclidean_tuple_distance>,
    time_poisson_topology, reach_plus_time_metric>;
// so3_2nd_order_rl_topology
extern template class temporal_space<
    reach_time_diff_space<time_topology,
                          arithmetic_tuple<rate_limited_quat_space<double>,
                                           ang_velocity_3D_topology<double>,
                                           ang_accel_3D_topology<double>>,
                          euclidean_tuple_distance>,
    time_poisson_topology, reach_plus_time_metric>;

}  // namespace ReaK::pp

#endif

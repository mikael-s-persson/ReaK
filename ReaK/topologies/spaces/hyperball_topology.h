/**
 * \file hyperball_topology.h
 *
 * This library provides classes that define a hyper-ball vector-topology. A hyper-ball vector-topology is
 * a vector-topology where the points are vector values and the boundary is a hyper-ellipsoid.
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

#ifndef REAK_TOPOLOGIES_SPACES_HYPERBALL_TOPOLOGY_H_
#define REAK_TOPOLOGIES_SPACES_HYPERBALL_TOPOLOGY_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/core/base/global_rng.h"
#include "ReaK/core/base/named_object.h"
#include "ReaK/math/lin_alg/vect_alg.h"
#include "ReaK/math/lin_alg/vect_concepts.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"

#include "ReaK/topologies/spaces/default_random_sampler.h"
#include "ReaK/topologies/spaces/vector_topology.h"

#include <cmath>
#include <random>
#include <type_traits>

namespace ReaK::pp {

/**
 * This library provides classes that define a hyper-ball vector-topology. A hyper-ball vector-topology is
 * a vector-topology where the points are vector values and the boundary is a hyper-ellipsoid.
 * This class models the MetricSpaceConcept, the LieGroupConcept, the BoundedSpaceConcept,
 * the SphereBoundedSpaceConcept, and the PointDistributionConcept.
 * \tparam Vector The vector-type for the topology, should model an Arithmetic concept and WritableVectorConcept.
 */
template <typename Vector>
class hyperball_topology : public vector_topology<Vector> {
 public:
  using self = hyperball_topology<Vector>;

  using point_type = Vector;
  using point_difference_type = Vector;

  using distance_metric_type = default_distance_metric;
  using random_sampler_type = default_random_sampler;

  static constexpr std::size_t dimensions = vect_traits<Vector>::dimensions;

 protected:
  point_type center_point;
  double radius_value;

 public:
  explicit hyperball_topology(const std::string& aName = "hyperball_topology",
                              const point_type& aOrigin = point_type(),
                              double aRadius = 1.0)
      : vector_topology<Vector>(aName),
        center_point(aOrigin),
        radius_value(aRadius) {}

  /*************************************************************************
   *                             MetricSpaceConcept
   * **********************************************************************/

  /**
   * Returns the distance between two points.
   */
  double distance(const point_type& a, const point_type& b) const {
    return this->norm(this->difference(b, a));
  }

  /**
   * Returns the norm of the difference between two points.
   */
  double norm(const point_difference_type& delta) const {
    using std::sqrt;
    double result = sqrt(delta * delta);
    return result;
  }

  /*************************************************************************
   *                         for PointDistributionConcept
   * **********************************************************************/

  /**
   * Generates a random point in the space, uniformly distributed.
   */
  point_type random_point() const {
    using std::pow;
    using std::sqrt;

    auto dp = this->difference(center_point, center_point);
    if (dp.size() == 0) {
      return center_point;
    }

    auto radial_dim_correction = double(dp.size());

    global_rng_type& rng = get_global_rng();
    std::normal_distribution<std::decay_t<decltype(dp[0])>> var_rnd;
    for (int i = 0; i < dp.size(); ++i) {
      dp[i] = var_rnd(rng);
    }

    std::uniform_real_distribution<double> uniform_rng;
    double factor = pow(uniform_rng(rng), 1.0 / radial_dim_correction) *
                    radius_value / sqrt(dp * dp);

    return this->adjust(center_point, factor * dp);
  }

  /*************************************************************************
   *                             BoundedSpaceConcept
   * **********************************************************************/

  /**
   * Takes a point and clips it to within this hyperball space.
   */
  void bring_point_in_bounds(point_type& a) const {
    a = this->adjust(a, this->get_diff_to_boundary(a));
  }

  /**
   * Returns the distance to the boundary of the space.
   */
  double distance_from_boundary(const point_type& a) const {
    using std::abs;
    auto c2a = this->difference(a, center_point);
    return abs(radius_value - this->norm(c2a));
  }

  /**
   * Returns the difference to the closest boundary.
   */
  point_difference_type get_diff_to_boundary(const point_type& a) const {
    auto c2a = this->difference(a, center_point);
    return (radius_value - this->norm(c2a)) * c2a;
  }

  /**
   * Tests if a given point is within the boundary of this space.
   */
  bool is_in_bounds(const point_type& a) const override {
    auto c2a = this->difference(a, center_point);
    return radius_value >= this->norm(c2a);
  }

  /**
   * Returns the origin of the space (the lower-limit).
   */
  point_type origin() const override { return center_point; }

  /*************************************************************************
   *                             SphereBoundedSpaceConcept
   * **********************************************************************/

  /**
   * Returns the radius of the space.
   */
  double get_radius() const { return radius_value; }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    ReaK::named_object::save(
        A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(center_point) &
        RK_SERIAL_SAVE_WITH_NAME(radius_value);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    ReaK::named_object::load(
        A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(center_point) &
        RK_SERIAL_LOAD_WITH_NAME(radius_value);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400008, 1, "hyperball_topology",
                              vector_topology<Vector>)
};

template <typename Vector>
struct is_metric_space<hyperball_topology<Vector>> : std::true_type {};

template <typename Vector>
struct is_reversible_space<hyperball_topology<Vector>> : std::true_type {};

template <typename Vector>
struct is_point_distribution<hyperball_topology<Vector>> : std::true_type {};

extern template class hyperball_topology<vect<double, 2>>;
extern template class hyperball_topology<vect<double, 3>>;
extern template class hyperball_topology<vect<double, 4>>;
extern template class hyperball_topology<vect<double, 6>>;
extern template class hyperball_topology<vect_n<double>>;

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_HYPERBALL_TOPOLOGY_H_

/**
 * \file hyperbox_topology.h
 *
 * This library provides classes that define a hyper-box vector-topology. A hyper-box vector-topology is
 * a vector-topology where the points are vector values and the boundary is a hyper-box.
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

#ifndef REAK_TOPOLOGIES_SPACES_HYPERBOX_TOPOLOGY_H_
#define REAK_TOPOLOGIES_SPACES_HYPERBOX_TOPOLOGY_H_

#include "ReaK/core/base/global_rng.h"
#include "ReaK/core/base/named_object.h"
#include "ReaK/math/lin_alg/vect_alg.h"
#include "ReaK/math/lin_alg/vect_concepts.h"

#include "ReaK/topologies/spaces/default_random_sampler.h"
#include "ReaK/topologies/spaces/vect_distance_metrics.h"
#include "ReaK/topologies/spaces/vector_topology.h"

#include <cmath>
#include <random>
#include <type_traits>

namespace ReaK::pp {

/**
 * This class defines a hyper-box vector-topology. A hyper-box vector-topology is
 * a vector-topology where the points are vector values and the boundary is a hyper-box.
 * This class models MetricSpace, LieGroup, BoundedSpace, and PointDistribution.
 * \tparam Vector The vector-type for the topology, should model an Arithmetic concept.
 */
template <WritableVector Vector, typename Metric = euclidean_distance_metric>
class hyperbox_topology : public vector_topology<Vector> {
 public:
  using self = hyperbox_topology<Vector, Metric>;

  using point_type = Vector;
  using point_difference_type = Vector;

  using distance_metric_type = Metric;
  using random_sampler_type = default_random_sampler;

  static constexpr std::size_t dimensions = vect_traits<Vector>::dimensions;

 protected:
  point_type lower_corner;
  point_type upper_corner;

 public:
  const point_type& get_lower_corner() const { return lower_corner; }
  const point_type& get_upper_corner() const { return upper_corner; }

  explicit hyperbox_topology(const std::string& aName = "hyperbox_topology",
                             const point_type& aLowerCorner = point_type(),
                             const point_type& aUpperCorner = point_type())
      : vector_topology<Vector>(aName),
        lower_corner(aLowerCorner),
        upper_corner(aUpperCorner) {}

  /*************************************************************************
  *                         for PointDistributionConcept
  * **********************************************************************/

  /**
   * Generates a random point in the space, uniformly distributed.
   */
  point_type random_point() const {
    global_rng_type& rng = get_global_rng();
    std::uniform_real_distribution<double> var_rnd;

    point_type p = lower_corner;
    for (int i = 0; i < p.size(); ++i) {
      p[i] += var_rnd(rng) * (upper_corner[i] - lower_corner[i]);
    }

    return p;
  }

  /*************************************************************************
  *                             BoundedSpace
  * **********************************************************************/

  /**
   * Takes a point and clips it to within this line-segment space.
   */
  void bring_point_in_bounds(point_type& a) const {
    for (int i = 0; i < a.size(); ++i) {
      if (lower_corner[i] < upper_corner[i]) {
        if (a[i] < lower_corner[i]) {
          a[i] = lower_corner[i];
        } else if (a[i] > upper_corner[i]) {
          a[i] = upper_corner[i];
        }
      } else {
        if (a[i] > lower_corner[i]) {
          a[i] = lower_corner[i];
        } else if (a[i] < upper_corner[i]) {
          a[i] = upper_corner[i];
        }
      }
    }
  }

  /**
   * Returns the distance to the boundary of the space.
   */
  double distance_from_boundary(const point_type& a) const {
    double dist = std::numeric_limits<std::decay_t<decltype(a[0])>>::max();
    using std::abs;
    for (int i = 0; i < a.size(); ++i) {
      if (dist > abs(a[i] - lower_corner[i])) {
        dist = abs(a[i] - lower_corner[i]);
      }
      if (dist > abs(a[i] - upper_corner[i])) {
        dist = abs(a[i] - upper_corner[i]);
      }
    }
    return dist;
  }

  /**
   * Returns the difference to the closest boundary.
   */
  point_difference_type get_diff_to_boundary(const point_type& a) const {
    double dist = std::numeric_limits<std::decay_t<decltype(a[0])>>::max();
    int j = a.size();
    bool at_upper = false;
    using std::abs;
    for (int i = 0; i < a.size(); ++i) {
      if (dist > abs(a[i] - lower_corner[i])) {
        j = i;
        at_upper = false;
        dist = abs(a[i] - lower_corner[i]);
      }
      if (dist > abs(a[i] - upper_corner[i])) {
        j = i;
        at_upper = true;
        dist = abs(a[i] - upper_corner[i]);
      }
    }
    auto dp = this->difference(a, a);
    if (j == a.size()) {
      return dp;
    }
    if (at_upper) {
      dp[j] = upper_corner[j] - a[j];
    } else {
      dp[j] = lower_corner[j] - a[j];
    }
    return dp;
  }

  /**
   * Tests if a given point is within the boundary of this space.
   */
  bool is_in_bounds(const point_type& a) const override {
    for (int i = 0; i < a.size(); ++i) {
      if (lower_corner[i] < upper_corner[i]) {
        if ((a[i] < lower_corner[i]) || (a[i] > upper_corner[i])) {
          return false;
        }
      } else {
        if ((a[i] > lower_corner[i]) || (a[i] < upper_corner[i])) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Returns the origin of the space (the center).
   */
  point_type origin() const override {
    return this->adjust(lower_corner,
                        0.5 * this->difference(upper_corner, lower_corner));
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    ReaK::named_object::save(A,
                             named_object::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(lower_corner) &
        RK_SERIAL_SAVE_WITH_NAME(upper_corner);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    ReaK::named_object::load(A,
                             named_object::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(lower_corner) &
        RK_SERIAL_LOAD_WITH_NAME(upper_corner);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400009, 1, "hyperbox_topology",
                              vector_topology<Vector>)
};

template <typename Vector, typename Metric>
struct is_metric_symmetric<hyperbox_topology<Vector, Metric>>
    : is_metric_symmetric<Metric> {};

extern template class hyperbox_topology<vect<double, 2>>;
extern template class hyperbox_topology<vect<double, 3>>;
extern template class hyperbox_topology<vect<double, 4>>;
extern template class hyperbox_topology<vect<double, 6>>;
extern template class hyperbox_topology<vect_n<double>>;

extern template class hyperbox_topology<vect<double, 2>,
                                        manhattan_distance_metric>;
extern template class hyperbox_topology<vect<double, 3>,
                                        manhattan_distance_metric>;
extern template class hyperbox_topology<vect<double, 4>,
                                        manhattan_distance_metric>;
extern template class hyperbox_topology<vect<double, 6>,
                                        manhattan_distance_metric>;
extern template class hyperbox_topology<vect_n<double>,
                                        manhattan_distance_metric>;

extern template class hyperbox_topology<vect<double, 2>,
                                        inf_norm_distance_metric>;
extern template class hyperbox_topology<vect<double, 3>,
                                        inf_norm_distance_metric>;
extern template class hyperbox_topology<vect<double, 4>,
                                        inf_norm_distance_metric>;
extern template class hyperbox_topology<vect<double, 6>,
                                        inf_norm_distance_metric>;
extern template class hyperbox_topology<vect_n<double>,
                                        inf_norm_distance_metric>;

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_HYPERBOX_TOPOLOGY_H_

/**
 * \file line_topology.h
 *
 * This library provides classes that define a line-topology. A line-topology is
 * a simple metric-space where the points are real values (doubles) along a 1D
 * space (line-segment).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2011
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

#ifndef REAK_TOPOLOGIES_SPACES_LINE_TOPOLOGY_H_
#define REAK_TOPOLOGIES_SPACES_LINE_TOPOLOGY_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/core/base/global_rng.h"
#include "ReaK/core/base/named_object.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/reversible_space_concept.h"

#include "ReaK/topologies/spaces/default_random_sampler.h"

#include <cmath>

namespace ReaK::pp {

/**
 * This class implements an infinite line topology. This class models Topology,
 * LieGroup, and MetricSpace.
 * \tparam T The value-type for the topology (should be an arithmetic type that is implicitly convertable to double).
 */
template <typename T = double>
class line_topology : public named_object {
 public:
  using self = line_topology<T>;

  using point_type = T;
  using point_difference_type = T;

  using distance_metric_type = default_distance_metric;

  static constexpr std::size_t dimensions = 1;

  explicit line_topology(const std::string& aName = "line_topology")
      : named_object() {
    setName(aName);
  }

  /*************************************************************************
   *                             Topology
   * **********************************************************************/

  /**
   * Returns the difference between two points (a - b).
   */
  point_difference_type difference(const point_type& a,
                                   const point_type& b) const {
    return a - b;
  }

  /**
   * Returns the addition of a point-difference to a point.
   */
  point_type adjust(const point_type& a,
                    const point_difference_type& delta) const {
    return a + delta;
  }

  /**
   * Returns the origin of the space (the lower-limit).
   */
  virtual point_type origin() const { return point_type(0); }

  /**
   * Tests if a given point is within the boundary of this space.
   */
  virtual bool is_in_bounds(const point_type& a) const { return true; }

  /*************************************************************************
  *                             MetricSpace
  * **********************************************************************/

  /**
   * Returns the distance between two points.
   */
  double distance(const point_type& a, const point_type& b) const {
    using std::abs;
    return abs(b - a);
  }

  /**
   * Returns the norm of the difference between two points.
   */
  double norm(const point_difference_type& delta) const {
    using std::abs;
    return abs(delta);
  }

  /**
   * Returns the volume of the difference between two points.
   */
  double volume(const point_difference_type& delta) const {
    using std::abs;
    return abs(delta);
  }

  /*************************************************************************
  *                             LieGroup
  * **********************************************************************/

  /**
   * Returns a point which is at a fraction between two points a to b.
   */
  point_type move_position_toward(const point_type& a, double fraction,
                                  const point_type& b) const {
    return a + (b - a) * fraction;
  }

  /**
   * Returns a point which is at a backward fraction between two points a to b.
   */
  point_type move_position_back_to(const point_type& a, double fraction,
                                   const point_type& b) const {
    return b + (a - b) * fraction;
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

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400001, 1, "line_topology",
                              named_object)
};

/**
 * This class implements a line-segment topology. The space extends from the minimum value up to some
 * maximum value. Models Topology, LieGroup, MetricSpace,
 * BoundedSpace, and SphereBoundedSpace, and also provides a distance-metric and a
 * random-sampler (default, uniform sampler).
 *
 * \tparam T A value-type (scalar value).
 */
template <typename T = double>
class line_segment_topology : public line_topology<T> {
  using rand_t = std::uniform_real_distribution<T>;

 public:
  using self = line_segment_topology<T>;

  using point_type = typename line_topology<T>::point_type;
  using point_difference_type =
      typename line_topology<T>::point_difference_type;

  using distance_metric_type = typename line_topology<T>::distance_metric_type;
  using random_sampler_type = default_random_sampler;

  static constexpr std::size_t dimensions = line_topology<T>::dimensions;

  /**
   * Default constructor.
   * \param aStart The minimum bound of the line-segment.
   * \param aEnd The maximum bound of the line-segment.
   */
  explicit line_segment_topology(
      const std::string& aName = "line_segment_topology",
      point_type aStart = point_type(0.0), point_type aEnd = point_type(1.0))
      : line_topology<T>(aName), start_pt(aStart), end_pt(aEnd) {}

  /*************************************************************************
  *                         for PointDistribution
  * **********************************************************************/

  /**
   * Generates a random point in the space, uniformly distributed.
   */
  point_type random_point() const {
    return rand_t()(get_global_rng()) * (end_pt - start_pt) + start_pt;
  }

  /*************************************************************************
  *                             BoundedSpace
  * **********************************************************************/

  /**
   * Takes a point and clips it to within this line-segment space.
   */
  void bring_point_in_bounds(point_type& a) const {
    if (end_pt > start_pt) {
      if (a > end_pt) {
        a = end_pt;
      } else if (a < start_pt) {
        a = start_pt;
      }
    } else {
      if (a < end_pt) {
        a = end_pt;
      } else if (a > start_pt) {
        a = start_pt;
      }
    }
  }

  /**
   * Returns the distance to the boundary of the space.
   */
  double distance_from_boundary(point_type a) const {
    using std::abs;
    double dist = abs(end_pt - a);
    if (abs(a - start_pt) < dist) {
      return abs(a - start_pt);
    }
    return dist;
  }

  /**
   * Returns the difference to the closest boundary.
   */
  point_difference_type get_diff_to_boundary(point_type a) const {
    using std::abs;
    double dist = abs(end_pt - a);
    if (abs(a - start_pt) < dist) {
      return start_pt - a;
    }
    return end_pt - a;
  }

  /**
   * Tests if a given point is within the boundary of this space.
   */
  bool is_in_bounds(const point_type& a) const override {
    if (end_pt > start_pt) {
      if ((a > end_pt) || (a < start_pt)) {
        return false;
      }
    } else {
      if ((a < end_pt) || (a > start_pt)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Returns the origin of the space (the lower-limit).
   */
  point_type origin() const override {
    return (end_pt - start_pt) * 0.5 + start_pt;
  }

  /*************************************************************************
  *                             SphereBoundedSpace
  * **********************************************************************/

  /**
   * Returns the radius of the space.
   */
  double get_radius() const {
    using std::abs;
    return abs(end_pt - start_pt) * 0.5;
  }

 private:
  point_type start_pt;
  point_type end_pt;

 public:
  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    line_topology<T>::save(
        A, line_topology<T>::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(start_pt) & RK_SERIAL_SAVE_WITH_NAME(end_pt);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    line_topology<T>::load(
        A, line_topology<T>::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(start_pt) & RK_SERIAL_LOAD_WITH_NAME(end_pt);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400006, 1, "line_segment_topology",
                              line_topology<T>)
};

extern template class line_topology<double>;
extern template class line_segment_topology<double>;

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_LINE_TOPOLOGY_H_

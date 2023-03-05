/**
 * \file time_topology.hpp
 *
 * This library provides classes that define a time-topology. A time-topology is
 * a simple metric-space where the points are real values (doubles) along a 1D
 * space. However, because time is unlimited, this topology, although modeling the
 * MetricSpaceConcept, is not strictly a metric-space since random-points cannot be
 * generated.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date September 2011
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

#ifndef REAK_TIME_TOPOLOGY_HPP
#define REAK_TIME_TOPOLOGY_HPP

#include "ReaK/core/base/defs.hpp"
#include "ReaK/core/base/named_object.hpp"

#include "ReaK/topologies/spaces/metric_space_concept.hpp"
#include "ReaK/topologies/spaces/reversible_space_concept.hpp"

#include <cmath>

namespace ReaK::pp {

/**
 * This class implements an infinite time-topology. Since the space is
 * infinite, there is no way to generate random points from it, and thus,
 * this class does not strictly model the topology concepts, but defines all
 * the functions required to provide the full model of a MetricSpaceConcept.
 */
class time_topology : public named_object {
 public:
  using point_type = double;
  using point_difference_type = double;

  using distance_metric_type = default_distance_metric;

  static constexpr std::size_t dimensions = 1;

  explicit time_topology(const std::string& aName = "time_topology") {
    setName(aName);
  }

  /*************************************************************************
   *                             TopologyConcept
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
   * Returns the origin of the space.
   */
  point_type origin() const { return 0.0; }

  /*************************************************************************
  *                             MetricSpaceConcept
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
  *                             LieGroupConcept
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

  /**
   * Returns the norm of the difference between two points.
   */
  point_type pointwise_min(const point_type& a, const point_type& b) const {
    return std::min(a, b);
  }

  /**
   * Returns the norm of the difference between two points.
   */
  point_type pointwise_max(const point_type& a, const point_type& b) const {
    return std::max(a, b);
  }

  /**
   * Tests if a given point is within the boundary of this space.
   */
  virtual bool is_in_bounds(const point_type& a) const { return true; }

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

  RK_RTTI_MAKE_CONCRETE_1BASE(time_topology, 0xC2400002, 1, "time_topology",
                              named_object)
};

template <>
struct is_metric_space<time_topology> : std::true_type {};

template <>
struct is_reversible_space<time_topology> : std::true_type {};

}  // namespace ReaK::pp

#endif

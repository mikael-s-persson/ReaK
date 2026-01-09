/**
 * \file vector_topology.h
 *
 * This library provides classes that define a vector-topology. A vector-topology is
 * a simple metric-space where the points are vector values.
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

#ifndef REAK_TOPOLOGIES_SPACES_VECTOR_TOPOLOGY_H_
#define REAK_TOPOLOGIES_SPACES_VECTOR_TOPOLOGY_H_

#include "ReaK/core/base/named_object.h"

#include "ReaK/math/lin_alg/vect_concepts.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/reversible_space_concept.h"

namespace ReaK::pp {

/**
 * This class implements an infinite vector topology. This class
 * models Topology and LieGroup.
 * \tparam Vector The vector-type for the topology, should model an Arithmetic concept and possess a norm() function.
 */
template <typename Vector>
class vector_topology : public named_object {
 public:
  using self = vector_topology<Vector>;

  using point_type = Vector;
  using point_difference_type = Vector;

  static constexpr std::size_t dimensions = 0;

  explicit vector_topology(const std::string& aName = "vector_topology")
      : named_object() {
    set_name(aName);
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
   * Returns the addition of a point-difference to a point.
   */
  virtual point_type origin() const { return point_type(); }

  /**
   * Tests if a given point is within the boundary of this space.
   */
  virtual bool is_in_bounds(const point_type& a) const { return true; }

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
    ReaK::named_object::save(A,
                             named_object::get_static_object_type()->version());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    ReaK::named_object::load(A,
                             named_object::get_static_object_type()->version());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400007, 1, "vector_topology",
                              named_object)
};

}  // namespace ReaK::pp

#include "ReaK/math/lin_alg/vect_alg.h"

namespace ReaK::pp {

extern template class vector_topology<vect<double, 2>>;
extern template class vector_topology<vect<double, 3>>;
extern template class vector_topology<vect<double, 4>>;
extern template class vector_topology<vect<double, 6>>;
extern template class vector_topology<vect_n<double>>;

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_VECTOR_TOPOLOGY_H_

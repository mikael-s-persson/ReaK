/**
 * \file bounded_space_concept.h
 *
 * This library defines the concepts that pertain to what can be considered
 * a bounded metric-space, as used in ReaK::pp. A bounded metric space can be
 * bounded
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

#ifndef REAK_TOPOLOGIES_SPACES_BOUNDED_SPACE_CONCEPT_H_
#define REAK_TOPOLOGIES_SPACES_BOUNDED_SPACE_CONCEPT_H_

#include "ReaK/topologies/spaces/metric_space_concept.h"

#include "boost/concept_check.hpp"

namespace ReaK::pp {

/**
 * This concept defines the requirements to fulfill in order to model a bounded metric-space
 * as used in ReaK::pp.
 *
 * Required Models:
 *
 * The topology should model the TopologyConcept.
 *
 * Valid expressions:
 *
 * pd = space.get_diff_to_boundary(p);  The difference (pd) between a point (p) and the closest boundary-point of the
 *metric-space can be obtained.
 *
 * bool b = space.is_in_bounds(p);  A point can be checked for being with the bounds of the space.
 *
 * space.bring_point_in_bounds(p);  A point (non-const ref) can be brought within bounds.
 *
 * \tparam BoundedSpace The topology type to be checked for this concept.
 */
template <typename BoundedSpace>
struct BoundedSpaceConcept {
  typename topology_traits<BoundedSpace>::point_type p;
  typename topology_traits<BoundedSpace>::point_difference_type pd;
  BoundedSpace space;

  BOOST_CONCEPT_ASSERT((TopologyConcept<BoundedSpace>));

  BOOST_CONCEPT_USAGE(BoundedSpaceConcept) {
    pd = space.get_diff_to_boundary(p);
    bool b = space.is_in_bounds(p);
    RK_UNUSED(b);
    space.bring_point_in_bounds(p);
  }
};

/**
 * This concept defines the requirements to fulfill in order to model a box-bounded metric-space
 * as used in ReaK::pp.
 *
 * Required Models:
 *
 * The topology should model the BoundedSpaceConcept.
 *
 * Valid expressions:
 *
 * p = space.get_upper_corner(p);  The point (p) at the upper boundary of the metric-space can be obtained.
 *
 * p = space.get_lower_corner(p);  The point (p) at the lower boundary of the metric-space can be obtained.
 *
 * \tparam BoundedSpace The topology type to be checked for this concept.
 */
template <typename BoundedSpace>
struct BoxBoundedSpaceConcept : BoundedSpaceConcept<BoundedSpace> {
  BOOST_CONCEPT_USAGE(BoxBoundedSpaceConcept) {
    this->p = this->space.get_upper_corner();
    this->p = this->space.get_lower_corner();
  }
};

/**
 * This concept defines the requirements to fulfill in order to model a sphere-bounded metric-space
 * as used in ReaK::pp.
 *
 * Required Models:
 *
 * The topology should model the BoundedSpaceConcept.
 *
 * Valid expressions:
 *
 * d = space.get_radius();  The radius (d) can be obtained (maximum distance between a point and the origin of the
 *metric-space).
 *
 * \tparam BoundedSpace The topology type to be checked for this concept.
 */
template <typename BoundedSpace>
struct SphereBoundedSpaceConcept : BoundedSpaceConcept<BoundedSpace> {
  double d;

  BOOST_CONCEPT_USAGE(SphereBoundedSpaceConcept) {
    d = this->space.get_radius();
    RK_UNUSED(d);
  }
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_BOUNDED_SPACE_CONCEPT_H_

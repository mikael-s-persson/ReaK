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

#include <concepts>

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
 */
template <typename Space>
concept BoundedSpace = Topology<Space>&& requires(
    const Space& space, const topology_point_type_t<Space>& p_in,
    topology_point_type_t<Space>& p_out) {
  {
    space.get_diff_to_boundary(p_in)
    } -> std::convertible_to<topology_point_difference_type_t<Space>>;
  { space.is_in_bounds(p_in) } -> std::convertible_to<bool>;
  space.bring_point_in_bounds(p_out);
};

/**
 * This concept defines the requirements to fulfill in order to model a box-bounded metric-space
 * as used in ReaK::pp.
 *
 * Required Models:
 *
 * The topology should model the BoundedSpace.
 *
 * Valid expressions:
 *
 * p = space.get_upper_corner(p);  The point (p) at the upper boundary of the metric-space can be obtained.
 *
 * p = space.get_lower_corner(p);  The point (p) at the lower boundary of the metric-space can be obtained.
 */
template <typename Space>
concept BoxBoundedSpace = BoundedSpace<Space>&& requires(const Space& space) {
  {
    space.get_lower_corner()
    } -> std::convertible_to<topology_point_type_t<Space>>;
  {
    space.get_upper_corner()
    } -> std::convertible_to<topology_point_type_t<Space>>;
};

/**
 * This concept defines the requirements to fulfill in order to model a sphere-bounded metric-space
 * as used in ReaK::pp.
 *
 * Required Models:
 *
 * The topology should model the BoundedSpace.
 *
 * Valid expressions:
 *
 * d = space.get_radius();  The radius (d) can be obtained (maximum distance between a point and the origin of the
 *metric-space).
 */
template <typename Space>
concept SphereBoundedSpace =
    BoundedSpace<Space>&& requires(const Space& space) {
  { space.get_radius() } -> std::convertible_to<double>;
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_BOUNDED_SPACE_CONCEPT_H_

/**
 * \file reversible_space_concept.h
 *
 * This library defines the traits and concepts that pertain to what can be considered
 * a reversible-space, as used in ReaK::pp. Reversible-spaces are based on the Topology concept
 * from the Boost.Graph library, but with additional requirements which are needed
 * in algorithms tailored for path-planning (mostly). Basically, the concept of a
 * reversible-space corresponds to spaces in which motions can
 * be done in reverse, whether steerable or not.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_TOPOLOGIES_SPACES_REVERSIBLE_SPACE_CONCEPT_H_
#define REAK_TOPOLOGIES_SPACES_REVERSIBLE_SPACE_CONCEPT_H_


#include <cmath>
#include <concepts>
#include <tuple>

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/steerable_space_concept.h"

/** Main namespace for ReaK */
namespace ReaK::pp {

/**
 * This concept defines the requirements to fulfill in order to model a reversible-space
 * as used in ReaK::pp. A reversible-space is a special kind of topology in which
 * motions can be done in reverse, whether steerable or not.
 *
 * Valid expressions:
 *
 * if(is_metric_space<ReversibleSpace>)
 *   p2 = space.move_position_back_to(p1, d, p3);  A point (p2) can be obtained by moving a fraction (d) back from a
 *final point (p3) to previous point (p1).
 *
 * if(is_steerable_space<ReversibleSpace>)
 *   tie(p2, st_rec) = space.steer_position_back_to(p1, d, p3);  A point (p2) and a steer-record (st_rec) can be
 *obtained from attempting to steer a fraction (d) back from a final point (p3) to previous point (p1).
 *
 * \tparam ReversibleSpace The topology type to be checked for this concept.
 */
template <typename Space>
concept ReversibleSpace =
    (MetricSpace<Space> || SteerableSpace<Space>)&&(
        !MetricSpace<Space> ||
        requires(const Space& s, const topology_point_type_t<Space>& p,
                 double d) {
          {
            s.move_position_back_to(p, d, p)
            } -> std::convertible_to<topology_point_type_t<Space>>;
        }) &&
    (!SteerableSpace<Space> ||
     requires(topology_point_type_t<Space> & p_out,
              steerable_space_steer_record_t<Space>& st_rec, const Space& s,
              const topology_point_type_t<Space>& p, double d) {
       std::tie(p_out, st_rec) = s.steer_position_back_to(p, d, p);
     });

template <typename Space>
static constexpr bool is_reversible_space_v = ReversibleSpace<Space>;

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_REVERSIBLE_SPACE_CONCEPT_H_

/**
 * \file motion_graph_structures.h
 *
 * This library defines a class to solve path planning problems using the
 * Rapidly-exploring Random Tree Star (RRT*) algorithm (or one of its variants).
 * Given a C_free (configuration space restricted to non-colliding points) and a
 * result reporting policy, this class will probabilistically construct a motion-graph
 * that will connect a starting point and a goal point with a path through C-free
 * that is as close as possible to the optimal path in terms of distance.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date August 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_PLANNING_PATH_PLANNING_MOTION_GRAPH_STRUCTURES_H_
#define REAK_PLANNING_PATH_PLANNING_MOTION_GRAPH_STRUCTURES_H_

// #define RK_PLANNERS_ENABLE_VEBL_TREE
// #define RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT


#include "bagl/adjacency_list.h"
#include "bagl/bfl_d_ary_tree.h"

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE
#include "bagl/vebl_d_ary_tree.h"
#endif

#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT
#include "ReaK/planning/path_planning/dvp_layout_adjacency_list.h"
#endif

#include "ReaK/topologies/spaces/metric_space_concept.h"

#include <type_traits>

namespace ReaK::pp {

template <typename FreeSpaceType>
struct motion_segment_directionality {
  using type = std::conditional_t<is_metric_symmetric_v<FreeSpaceType>,
                                  bagl::undirected_s, bagl::bidirectional_s>;
};

template <typename FreeSpaceType>
using motion_segment_directionality_t =
    typename motion_segment_directionality<FreeSpaceType>::type;

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_MOTION_GRAPH_STRUCTURES_H_

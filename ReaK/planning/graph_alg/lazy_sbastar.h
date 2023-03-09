/**
 * \file lazy_sbastar.h
 *
 * This library provides function templates and concepts that implement a Lazy Sampling-based A* search
 * algorithm. A Lazy-SBA* uses the A* search algorithm to drive the expansion of a roadmap into the free-space
 * in order to connect a start and goal location. This algorithm has many customization points because there
 * are many choices to be made in the method, such as how to find nearest neighbors for attempting to
 * connect them through free-space, how to expand vertices, when to stop the algorithm, etc.
 * All these customization points are left to the user to implement, some are defined by the
 * SBAStarVisitorConcept (random-walk, edge-added, etc.).
 *
 * The Lazy-SBA* algorithm is a generalization of the A* algorithm where the neighborhood of a given node of
 * the motion graph is not defined as a fixed set of neighbors (as in a classic A* over a fixed graph),
 * but rather as a region from which samples can be drawn (biased or not). In an ordinary A* algorithm,
 * vertices are closed when their entire neighborhood has been explored. In an SBA* algorithm, the same
 * criteria cannot apply since samples could be drawn ad infinitum, so, instead, this concept of the
 * neighborhood being fully explored is derived from the expected information gained (or conversely, the
 * "surprisal") from drawing a new sample in the neighborhood. In this lazy version, the computation of the
 * edge weights as a cost-to-go through the free-space is tentatively replaced by the cost-to-go in the
 * configuration space (without obstacles), and collision along the path is only performed once the edge
 * has relaxed (identified as a segment of the local optimal path).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_PLANNING_GRAPH_ALG_LAZY_SBASTAR_H_
#define REAK_PLANNING_GRAPH_ALG_LAZY_SBASTAR_H_

#include "ReaK/planning/graph_alg/sbastar_search.h"

namespace ReaK::graph {}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_LAZY_SBASTAR_H_

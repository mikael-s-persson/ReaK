/**
 * \file tree_organizer_concept.h
 *
 * This library provides a concept class that defines a tree-organizing visitor. Such a visitor is meant
 * to be used to lump insertion / deletion callbacks that maintain a tree structure.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2012
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

#ifndef REAK_PLANNING_GRAPH_ALG_TREE_ORGANIZER_CONCEPT_H_
#define REAK_PLANNING_GRAPH_ALG_TREE_ORGANIZER_CONCEPT_H_

#include "ReaK/planning/graph_alg/simple_graph_traits.h"

namespace ReaK::graph {

/**
 * This concept defines the requirements to fulfill in order to model a tree-organizing visitor concept
 * as used in ReaK::graph.
 *
 * Valid expressions (vertex v; vertex_iterator cv_it, cv_it_end; TreeType tree):
 *
 * vis.remove_vertex(v, tree);  The visitor can perform the removal of the vertex (v) from the tree.
 *
 * v = vis.add_vertex(tree_vp, tree);  The visitor can perform the addition a vertex (v) with property (vp) to the tree.
 *
 * v = vis.add_vertex(std::move(tree_vp), tree);  The visitor can perform the addition a vertex (v) by moving in the
 *property (vp) to the tree.
 *
 * \tparam TreeOrganizerVisitor The visitor type to be checked for this concept.
 * \tparam TreeType The tree type on which the visitor must operate.
 */
template <typename Visitor, typename TreeType>
concept TreeOrganizerVisitor = requires(Visitor vis, TreeType tree,
                                        graph_vertex_t<TreeType> v,
                                        graph_vertex_property_t<TreeType> vp) {
  vis.remove_vertex(v, tree);
  vis.add_vertex(vp, tree);
  vis.add_vertex(std::move(vp), tree);
};

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_TREE_ORGANIZER_CONCEPT_H_

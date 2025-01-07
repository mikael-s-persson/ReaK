/**
 * \file rr_graph.h
 *
 * This library contains the Rapidly-Exploring Random Graph generation algorithm.
 * This is a method to create a random graph that will span over a non-convex space
 * as rapidly as possible. The method relies on a simple randomized insertion algorithm.
 * At each step, a random point is picked from the underlying topology (i.e. configuration
 * space in path-planning terms). Then, the points in the current graph that are nearest
 * to the random point are picked for expansion. Finally, edges (of a maximum length) are
 * added to the vertex of the graph towards the random point while it is still possible to
 * add such an edge without leaving the free space (the part of the configuration space which
 * is not occupied by an obstacle). The algorithm will stop when either the number of vertices
 * in the tree has reached a maximum or when the user callback signals the stop.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2012
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

#ifndef REAK_PLANNING_GRAPH_ALG_RR_GRAPH_H_
#define REAK_PLANNING_GRAPH_ALG_RR_GRAPH_H_

#include <iterator>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include "bagl/graph_concepts.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/random_sampler_concept.h"

#include "ReaK/planning/graph_alg/node_generators.h"
#include "ReaK/planning/graph_alg/sbmp_visitor_concepts.h"

namespace ReaK::graph {

namespace rrg_detail {

template <typename Graph, typename Topology, typename RRGVisitor,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          typename NodeGenerator, typename NcSelector>
void generate_rrg_loop(Graph& g, const Topology& space, RRGVisitor vis,
                       PositionMap position, NodeGenerator node_generator_func,
                       NcSelector select_neighborhood,
                       std::size_t max_vertex_count) {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;

  while ((num_vertices(g) < max_vertex_count) && (vis.keep_going())) {

    auto [x_near, p_new, eprop] = node_generator_func(
        g, vis,
        bagl::composite_property_map(position, get(bagl::vertex_all, g)));

    std::vector<Vertex> pred;
    std::vector<Vertex> succ;
    if constexpr (bagl::is_undirected_graph_v<Graph>) {
      select_neighborhood(
          p_new, std::back_inserter(pred), g, space,
          bagl::composite_property_map(position, get(bagl::vertex_all, g)));
    } else {
      select_neighborhood(
          p_new, std::back_inserter(pred), std::back_inserter(succ), g, space,
          bagl::composite_property_map(position, get(bagl::vertex_all, g)));
    }

    bagl::vertex_property_type<Graph> xp_new;
    put(position, xp_new, p_new);
    Vertex x_new = add_vertex(g, std::move(xp_new));
    vis.vertex_added(x_new, g);

    {
      auto [e_new, e_new_exists] = add_edge(x_near, x_new, g, std::move(eprop));
      if (e_new_exists) {
        vis.edge_added(e_new, g);
      }
    }

    for (Vertex u : pred) {
      if (u == x_near) {
        continue;
      }

      auto [can_connect, eprop_new] = vis.can_be_connected(u, x_new, g);
      if (!can_connect) {
        continue;
      }

      auto [e_new, e_new_exists] = add_edge(u, x_new, g, std::move(eprop_new));
      if (e_new_exists) {
        vis.edge_added(e_new, g);
      }
    }

    if constexpr (!bagl::is_undirected_graph_v<Graph>) {
      for (Vertex u : succ) {
        auto [can_connect, eprop_new] = vis.can_be_connected(x_new, u, g);
        if (!can_connect) {
          continue;
        }

        auto [e_new, e_new_exists] =
            add_edge(x_new, u, g, std::move(eprop_new));
        if (e_new_exists) {
          vis.edge_added(e_new, g);
        }
      }
    }
  }
}

}  // namespace rrg_detail

/**
  * This function template is the unidirectional version of the RRG algorithm (refer to rr_graph.hpp dox).
  * \tparam Graph A mutable graph type that will represent the generated tree, should model
  *bagl::concepts::VertexListGraph and bagl::concepts::MutableGraph
  * \tparam Topology A topology type that will represent the space in which the configurations (or positions) exist,
  *should model BGL's Topology concept
  * \tparam RRGVisitor An RRG visitor type that implements the customizations to this RRG algorithm, should model the
  *RRGVisitorConcept.
  * \tparam PositionMap A property-map type that can store the configurations (or positions) of the vertices.
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \tparam NcSelector A functor type which can perform a neighborhood search of a point to a graph in the topology (see
  *topological_search.hpp).
  * \param g A mutable graph that should initially store the starting and goal
  *        vertex and will store the generated graph once the algorithm has finished.
  * \param space A topology (as defined by the Born Again Graph Library). Note
  *        that it is not required to generate only random points in
  *        the free-space.
  * \param vis A RRG visitor implementing the RRGVisitorConcept. This is the
  *        main point of customization and recording of results that the
  *        user can implement.
  * \param position A mapping that implements the MutablePropertyMap Concept. Also,
  *        the value_type of this map should be the same type as the topology's
  *        value_type.
  * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
  * \param select_neighborhood A callable object (functor) which can perform a
  *        nearest neighbor search of a point to a graph in the topology. (see star_neighborhood)
  * \param max_vertex_count The maximum number of vertices beyond which the algorithm
  *        should stop regardless of whether the resulting tree is satisfactory or not.
  *
  */
template <bagl::concepts::MutableGraph Graph, typename Topology,
          RRGVisitorConcept<Graph, Topology> RRGVisitor,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          pp::RandomSamplerConcept<Topology> RandomSampler, typename NcSelector>
void generate_rrg(Graph& g, const Topology& space, RRGVisitor vis,
                  PositionMap position, RandomSampler get_sample,
                  NcSelector select_neighborhood,
                  std::size_t max_vertex_count) {
  if (num_vertices(g) == 0) {
    auto p = get_sample(space);
    while (!vis.is_position_free(p)) {
      p = get_sample(space);
    }

    bagl::vertex_property_type<Graph> up;
    put(position, up, p);
    auto u = add_vertex(g, std::move(up));
    vis.vertex_added(u, g);
  }

  rrg_detail::generate_rrg_loop(
      g, space, vis, position,
      rrg_node_generator(&space, get_sample, select_neighborhood),
      select_neighborhood, max_vertex_count);
}

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_RR_GRAPH_H_

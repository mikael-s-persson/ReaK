/**
 * \file probabilistic_roadmap.h
 *
 * This library provides function templates and concepts that implement the probabilistic roadmap (PRM)
 * algorithm (as of "Geraerts and Overmars, 2002"). A PRM is generated in two phases. First, a number of
 * vertices are generated at random from the free-space (configuration-space which is not occupied by
 * obstacles) and bi-directional (or undirected) connections between those vertices are established, as
 * much as possible (if they can be directly connected through free-space). Second, difficult areas
 * of the roadmap are determined using some density function (heuristic) and a priority-queue of vertices
 * with high density (whatever that means for the heuristic chosen). Then, high-density vertices are
 * selected for expansion, in which case a number of neighboring vertices are generated coming out of
 * the selected vertex into free-space. The first phase is called the construction phase, and the
 * second is called the expansion phase. This algorithm has many customization points because there
 * are many choices to be made in the method, such as how to find nearest neighbors for attempting to
 * connect them through free-space, how to expand vertices, how to measure density, how to compare
 * densities, when to stop the algorithm, what proportion of constructed vs. expanded vertices should
 * be generated, etc. All these customization points are left to the user to implement, some are
 * defined by the PRMVisitorConcept (expand-vertex and compute-density) while others are provided
 * as functors to the function template that generates the PRM (generate_prm).
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

#ifndef REAK_PLANNING_GRAPH_ALG_PROBABILISTIC_ROADMAP_H_
#define REAK_PLANNING_GRAPH_ALG_PROBABILISTIC_ROADMAP_H_

#include <limits>
#include <set>
#include <stack>
#include <tuple>
#include <utility>

#include "bagl/d_ary_heap.h"
#include "bagl/graph_concepts.h"
#include "bagl/more_property_maps.h"
#include "bagl/property_map.h"

#include "ReaK/core/base/global_rng.h"
#include "ReaK/planning/graph_alg/prm_connector.h"
#include "ReaK/planning/graph_alg/sbmp_visitor_concepts.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/random_sampler_concept.h"

namespace ReaK::graph {

namespace prm_detail {

template <typename PRMVisitor, typename UpdatableQueue, typename IndexInHeapMap,
          typename PositionMap, typename CCRootMap>
struct prm_conn_visitor {

  using PositionValue = bagl::property_traits_value_t<PositionMap>;
  using CCRootValue = bagl::property_traits_value_t<CCRootMap>;

  prm_conn_visitor(PRMVisitor vis, UpdatableQueue& q,
                   IndexInHeapMap index_in_heap, PositionMap pos,
                   CCRootMap cc_root, std::set<CCRootValue>& cc_set)
      : vis_(vis),
        q_(q),
        index_in_heap_(index_in_heap),
        position_(pos),
        cc_root_(cc_root),
        cc_set_(cc_set) {}

  template <class Graph>
  auto create_vertex(const PositionValue& p, Graph& g) const {
    using Vertex = bagl::graph_vertex_descriptor_t<Graph>;
    using VertexProp = bagl::vertex_property_type<Graph>;

    VertexProp up;
    put(position_, up, p);
    Vertex u = add_vertex(g, std::move(up));
    put(cc_root_, u, u);
    cc_set_.insert(u);
    vis_.vertex_added(u, g);
    put(index_in_heap_, u, std::numeric_limits<std::size_t>::max());
    return u;
  }

  template <typename Vertex, typename Graph>
  void shortcut_cc_root(Vertex u, Graph& g) const {
    std::stack<Vertex> u_trace;
    u_trace.push(u);
    while (get(cc_root_, u) != u) {
      u = get(cc_root_, u);
      u_trace.push(u);
    }
    while (!u_trace.empty()) {
      put(cc_root_, u_trace.top(), u);
      u_trace.pop();
    }
  }

  template <typename Edge, typename Graph>
  void edge_added(Edge e, Graph& g) const {
    using Vertex = bagl::graph_vertex_descriptor_t<Graph>;

    vis_.edge_added(e, g);

    Vertex u = source(e, g);
    Vertex v = target(e, g);
    shortcut_cc_root(u, g);
    shortcut_cc_root(v, g);
    if (get(cc_root_, v) != get(cc_root_, u)) {
      Vertex r1 = get(cc_root_, u);
      Vertex r2 = get(cc_root_, v);
      put(cc_root_, r2, r1);
      put(cc_root_, v, r1);
      cc_set_.erase(r2);

      if (cc_set_.size() < 2) {
        vis_.publish_path(g);
      }
    }
  }
  template <typename Vertex, typename Graph>
  void vertex_added(Vertex u, Graph& g) const {
    vis_.vertex_added(u, g);
  }

  template <typename Vertex, typename Graph>
  void travel_explored(Vertex u, Vertex v, Graph& g) const {
    vis_.travel_explored(u, v, g);
  }

  template <typename Vertex, typename Graph>
  void travel_succeeded(Vertex u, Vertex v, Graph& g) const {
    vis_.travel_succeeded(u, v, g);
  }

  template <typename Vertex, typename Graph>
  void travel_failed(Vertex u, Vertex v, Graph& g) const {
    vis_.travel_failed(u, v, g);
  }

  template <typename Vertex, typename Graph>
  void requeue_vertex(Vertex u, Graph& g) const {
    vis_.affected_vertex(u, g);
    q_.push_or_update(u);
  }
  template <typename Vertex, typename Graph>
  void affected_vertex(Vertex u, Graph& g) const {
    // same function, different name.
    requeue_vertex(u, g);
  }

  bool keep_going() const { return vis_.keep_going(); }

  template <typename Vertex, typename Graph>
  auto can_be_connected(Vertex u, Vertex v, Graph& g) const {
    return vis_.can_be_connected(u, v, g);
  }

  template <typename Vertex, typename Graph>
  auto random_walk(Vertex u, Graph& g) const {
    return vis_.random_walk(u, g);
  }

  bool is_position_free(const PositionValue& p) const {
    return vis_.is_position_free(p);
  }

  PRMVisitor vis_;
  UpdatableQueue& q_;
  IndexInHeapMap index_in_heap_;

  PositionMap position_;
  CCRootMap cc_root_;

  std::set<CCRootValue>& cc_set_;
};

template <typename Graph, pp::MetricSpace Space, typename PRMConnVisitor,
          bagl::concepts::ReadableVPropMemberMap<Graph> PositionMap,
          pp::RandomSampler<Space> Sampler, typename MutableQueue,
          typename NodeConnector, typename NcSelector>
void generate_prm_impl(Graph& g, const Space& super_space, PRMConnVisitor& vis,
                       PositionMap position, Sampler get_sample,
                       MutableQueue& q, NodeConnector connect_vertex,
                       const NcSelector& select_neighborhood,
                       double expand_probability) {
  using EdgeProp = bagl::edge_property_type<Graph>;

  std::uniform_real_distribution<double> unit_dist{};

  while (vis.keep_going()) {
    // Generate random-number between 0 and 1.
    double rand_value = unit_dist(ReaK::get_global_rng());

    if (rand_value > expand_probability) {
      // Construction node:
      auto p_rnd = get_sample(super_space);
      while (!vis.is_position_free(p_rnd)) {
        p_rnd = get_sample(super_space);
      }
      EdgeProp ep;
      connect_vertex(p_rnd, bagl::graph_traits<Graph>::null_vertex(), ep, g,
                     super_space, vis, position, select_neighborhood);
    } else {
      // Expansion node:
      // use the priority queue to get the vertices that need expansion.
      auto v = q.top();
      auto [p_rnd, expanding_worked, ep] = vis.random_walk(v, g);
      if (expanding_worked) {
        connect_vertex(p_rnd, v, ep, g, super_space, vis, position,
                       select_neighborhood);
      } else {
        // if one cannot expand from this vertex, then leave it out of future attempts.
        q.pop();
      }
    }
  }
}

}  // namespace prm_detail

/**
 * This function is the basic PRM algorithm (as of "Geraerts and Overmars, 2002"). Note that
 * all the design decisions (how to generate vertices, how to select the neighborhood, how to
 * check the collision) have been externalized via the provided free_space topology, visitor
 * object, and neighborhood selection functor.
 *
 * \tparam Graph A mutable graph type that can store the roadmap (either bidirectional or
 *         undirected graph, the algorithm will deal with either cases as appropriate).
 * \tparam Space A topology type on which the vertex positions lie.
 * \tparam Visitor A PRM visitor type.
 * \tparam PositionMap A property-map type that can store the position of each vertex.
 * \tparam DensityMap A property-map type that can store the density-measures for each vertex.
 * \tparam NcSelector A functor type that can select a list of vertices of the graph that are
 *         the nearest-neighbors of a given vertex (or some other heuristic to select the neighbors).
 *         See classes in the topological_search.hpp header-file.
 * \tparam CompareFunction A functor type that can be used to compare density values (strict weak-ordering).
 * \param g A mutable graph that should initially store the starting
 *        vertex (if not it will be randomly generated) and will store
 *        the generated graph once the algorithm has finished.
 * \param super_space A topology (as defined by ReaK) that represents the configuration-space
 *        used for the planning (a space with no obstacles).
 * \param vis A PRM visitor implementing the PRMVisitor. This is the
 *        main point of customization and recording of results that the
 *        user can implement.
 * \param position A mapping that implements the MutablePropertyMap Concept. Also,
 *        the value_type of this map should be the same type as the topology's
 *        value_type.
 * \param get_sample A random sampler of positions in the super-space (obstacle-free topology).
 * \param density A property-map that provides the density values assiciated to each vertex.
 * \param select_neighborhood A callable object (functor) that can select a list of
 *        vertices of the graph that ought to be connected to a new
 *        vertex. The list should be sorted in order of increasing "distance".
 * \param expand_probability The probability (between 0 and 1) that an expansion will be
 *        performed as opposed to a general sampling from the free-space.
 */
template <bagl::concepts::VertexListGraph Graph, pp::MetricSpace Space,
          PRMVisitor<Graph, Space> Visitor,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          pp::RandomSampler<Space> Sampler,
          bagl::concepts::ReadableVPropMemberMap<Graph> DensityMap,
          typename NcSelector>
inline void generate_prm(Graph& g, const Space& super_space, Visitor vis,
                         PositionMap position, Sampler get_sample,
                         DensityMap density,
                         const NcSelector& select_neighborhood,
                         double expand_probability) {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;
  using VertexProp = bagl::vertex_property_type<Graph>;

  auto cc_root =
      bagl::vector_property_map(num_vertices(g), get(bagl::vertex_index, g),
                                bagl::graph_traits<Graph>::null_vertex());
  std::set<Vertex> cc_set;

  auto index_in_heap =
      bagl::vector_property_map(num_vertices(g), get(bagl::vertex_index, g),
                                std::numeric_limits<std::size_t>::max());

  auto g_density =
      bagl::composite_property_map(density, get(bagl::vertex_all, g));

  // priority queue holding the "expandable" vertices.
  auto q = bagl::make_d_ary_heap_indirect<Vertex, 4>(
      g_density, index_in_heap.ref(), std::less<>());

  if (num_vertices(g) == 0) {
    auto p = get_sample(super_space);
    while (!vis.is_position_free(p)) {
      p = get_sample(super_space);
    }
    VertexProp up;
    put(position, up, p);
    Vertex u = add_vertex(g, std::move(up));
    put(cc_root, u, u);
    vis.vertex_added(u, g);
    cc_set.insert(u);
  } else {
    for (auto u : vertices(g)) {
      vis.affected_vertex(u, g);
      q.push(u);
      put(cc_root, u, u);
    }

    for (auto u : vertices(g)) {
      if (get(cc_root, u) != u) {
        continue;
      }
      cc_set.insert(u);
      std::stack<Vertex> u_trace;
      u_trace.push(u);
      while (!u_trace.empty()) {
        Vertex v = u_trace.top();
        u_trace.pop();
        if (get(cc_root, v) == u) {
          continue;
        }
        put(cc_root, v, u);
        for (auto e : out_edges(v, g)) {
          if (target(e, g) == v) {
            u_trace.push(source(e, g));
          } else {
            u_trace.push(target(e, g));
          }
        }
      }
    }
  }

  prm_detail::prm_conn_visitor prm_conn_vis(vis, q, index_in_heap.ref(),
                                            position, cc_root.ref(), cc_set);

  prm_detail::generate_prm_impl(g, super_space, prm_conn_vis, position,
                                get_sample, q, prm_node_connector(),
                                select_neighborhood, expand_probability);
}

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_PROBABILISTIC_ROADMAP_H_

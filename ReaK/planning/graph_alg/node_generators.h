/**
 * \file node_generators.h
 *
 * This library contains a node generator that can be used in RRT-like algorithms.
 * Essentially, the node generator is a callable object that will perform a "Voronoi-pull"
 * operation, characteristic of RRT-style algorithms.
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

#ifndef REAK_PLANNING_GRAPH_ALG_NODE_GENERATORS_H_
#define REAK_PLANNING_GRAPH_ALG_NODE_GENERATORS_H_

#include <iterator>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/random_sampler_concept.h"

namespace ReaK::graph {

namespace node_gen_detail {

template <typename Graph>
struct rrg_node_puller {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;
  using EdgeProp = bagl::edge_property_type<Graph>;

  template <typename PositionValue, typename RRGVisitor,
            typename PredecessorMap = std::nullptr_t>
  static auto expand_to_nearest(PositionValue& p_new,
                                const std::vector<Vertex>& pred, Graph& g,
                                RRGVisitor vis,
                                PredecessorMap predecessor = nullptr) {
    const Vertex null_v = bagl::graph_traits<Graph>::null_vertex();
    for (Vertex u : pred) {
      if constexpr (!std::is_same_v<PredecessorMap, std::nullptr_t>) {
        if (get(predecessor, get_property(g, u)) == null_v) {
          continue;
        }
      }
      auto [p_tmp, expand_succeeded, ep] =
          vis.steer_towards_position(p_new, u, g);
      if (expand_succeeded) {
        p_new = p_tmp;
        return std::tuple(u, true, ep);
      }
    }
    return std::tuple(null_v, false, EdgeProp());
  }

  template <typename PositionValue, typename RRGVisitor,
            typename SuccessorMap = std::nullptr_t>
  static auto retract_from_nearest(PositionValue& p_new,
                                   const std::vector<Vertex>& succ, Graph& g,
                                   RRGVisitor vis,
                                   SuccessorMap successor = nullptr) {
    const Vertex null_v = bagl::graph_traits<Graph>::null_vertex();
    for (Vertex u : succ) {
      if constexpr (!std::is_same_v<SuccessorMap, std::nullptr_t>) {
        if (get(successor, get_property(g, u)) == null_v) {
          continue;
        }
      }
      auto [p_tmp, retract_succeeded, ep] =
          vis.steer_back_to_position(p_new, u, g);
      if (retract_succeeded) {
        p_new = p_tmp;
        return std::tuple(u, true, ep);
      }
    }
    return std::tuple(null_v, false, EdgeProp());
  }
};

}  // namespace node_gen_detail

/**
 * This node generator that can be used in RRT-like algorithms.
 * Essentially, the node generator is a callable object that will perform a "Voronoi-pull"
 * operation, characteristic of RRT-style algorithms.
 * \tparam Topology The topology type on which the planning is performed (i.e., the configuration space type).
 * \tparam RandomSampler The type of the random-sampler that can generate random points in the configuration space,
 * should model ReaK::pp::RandomSamplerConcept.
 * \tparam NcSelector The type of a functor that can be used to perform nearest-neighbor queries.
 */
template <typename Topology, typename RandomSampler, typename NcSelector>
struct rrg_node_generator {

  const Topology* space;
  RandomSampler get_sample;
  NcSelector select_neighborhood;

  rrg_node_generator(const Topology* s, RandomSampler gs, NcSelector sn)
      : space(s), get_sample(gs), select_neighborhood(sn) {}

  template <typename Graph, typename RRGVisitor,
            bagl::concepts::ReadableVertexPropertyMap<Graph> PositionMap>
  auto operator()(Graph& g, RRGVisitor vis, PositionMap g_position) const {
    using Vertex = bagl::graph_vertex_descriptor_t<Graph>;
    using NodePuller = node_gen_detail::rrg_node_puller<Graph>;

    std::size_t i = 0;
    while (true) {
      auto p_new = get_sample(*space);

      std::vector<Vertex> pred;
      std::vector<Vertex> succ;
      if constexpr (bagl::is_undirected_graph_v<Graph>) {
        select_neighborhood(p_new, std::back_inserter(pred), g, *space,
                            g_position);
      } else {
        select_neighborhood(p_new, std::back_inserter(pred),
                            std::back_inserter(succ), g, *space, g_position);
      }

      auto [x_near, was_expanded, ep] =
          NodePuller::expand_to_nearest(p_new, pred, g, vis);
      if (was_expanded) {
        return std::tuple(x_near, p_new, ep);
      }
      if (i >= 10) {
        return std::tuple(bagl::graph_traits<Graph>::null_vertex(), p_new, ep);
      }
      ++i;
    }
  }
};

/**
 * This bidirectional node generator that can be used in RRT-like algorithms.
 * Essentially, the node generator is a callable object that will perform a "Voronoi-pull"
 * operation, characteristic of RRT-style algorithms.
 * \tparam Topology The topology type on which the planning is performed (i.e., the configuration space type).
 * \tparam RandomSampler The type of the random-sampler that can generate random points in the configuration space,
 * should model ReaK::pp::RandomSamplerConcept.
 * \tparam NcSelector The type of a functor that can be used to perform nearest-neighbor queries.
 */
template <typename Topology, typename RandomSampler, typename NcSelector,
          typename PredecessorMap, typename SuccessorMap>
struct rrg_bidir_generator {

  const Topology* space;
  RandomSampler get_sample;
  NcSelector select_neighborhood;
  PredecessorMap predecessor;
  SuccessorMap successor;

  rrg_bidir_generator(const Topology* s, RandomSampler gs, NcSelector sn,
                      PredecessorMap pred, SuccessorMap succ)
      : space(s),
        get_sample(gs),
        select_neighborhood(sn),
        predecessor(pred),
        successor(succ) {}

  template <typename Graph, typename RRGVisitor,
            bagl::concepts::ReadableVertexPropertyMap<Graph> PositionMap>
  auto operator()(Graph& g, RRGVisitor vis, PositionMap g_position) const {
    using Vertex = bagl::graph_vertex_descriptor_t<Graph>;
    using NodePuller = node_gen_detail::rrg_node_puller<Graph>;

    const Vertex null_v = bagl::graph_traits<Graph>::null_vertex();

    int i = 0;
    while (true) {
      auto p_pred = get_sample(*space);
      auto p_succ = p_pred;

      std::vector<Vertex> pred;
      std::vector<Vertex> succ;
      if constexpr (bagl::is_undirected_graph_v<Graph>) {
        select_neighborhood(p_pred, std::back_inserter(pred), g, *space,
                            g_position);
      } else {
        select_neighborhood(p_pred, std::back_inserter(pred),
                            std::back_inserter(succ), g, *space, g_position);
      }

      auto [x_pred, was_expanded, ep_pred] =
          NodePuller::expand_to_nearest(p_pred, pred, g, vis, predecessor);
      if constexpr (bagl::is_undirected_graph_v<Graph>) {
        pred.swap(succ);
      }
      auto [x_succ, was_retracted, ep_succ] =
          NodePuller::retract_from_nearest(p_succ, succ, g, vis, successor);
      if (was_expanded || was_retracted) {
        return std::tuple(x_pred, p_pred, ep_pred, x_succ, p_succ, ep_succ);
      }
      if (i >= 10) {
        return std::tuple(null_v, p_pred, ep_pred, null_v, p_succ, ep_succ);
      }
      ++i;
    }
  }
};

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_NODE_GENERATORS_H_

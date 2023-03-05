/**
 * \file node_generators.hpp
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

#ifndef REAK_NODE_GENERATORS_HPP
#define REAK_NODE_GENERATORS_HPP

#include <iterator>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include "boost/graph/graph_concepts.hpp"
#include "boost/property_map/property_map.hpp"

#include "ReaK/topologies/spaces/metric_space_concept.hpp"
#include "ReaK/topologies/spaces/random_sampler_concept.hpp"

#include "ReaK/planning/graph_alg/simple_graph_traits.hpp"

namespace ReaK::graph {

namespace detail {
namespace {

template <typename Graph>
struct rrg_node_puller {
  using Vertex = graph_vertex_t<Graph>;
  using EdgeProp = graph_edge_bundle_t<Graph>;

  template <typename PositionValue, typename RRGVisitor,
            typename PredecessorMap = std::nullptr_t>
  static auto expand_to_nearest(PositionValue& p_new,
                                const std::vector<Vertex>& Pred, Graph& g,
                                RRGVisitor vis,
                                PredecessorMap predecessor = nullptr) {
    const Vertex null_v = boost::graph_traits<Graph>::null_vertex();
    for (Vertex u : Pred) {
      if constexpr (!std::is_same_v<PredecessorMap, std::nullptr_t>) {
        if (get(predecessor, g[u]) == null_v) {
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
                                   const std::vector<Vertex>& Succ, Graph& g,
                                   RRGVisitor vis,
                                   SuccessorMap successor = nullptr) {
    const Vertex null_v = boost::graph_traits<Graph>::null_vertex();
    for (Vertex u : Succ) {
      if constexpr (!std::is_same_v<SuccessorMap, std::nullptr_t>) {
        if (get(successor, g[u]) == null_v) {
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
}  // namespace
}  // namespace detail

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

  rrg_node_generator(const Topology* aSpace, RandomSampler aGetSample,
                     NcSelector aSelectNeighborhood)
      : space(aSpace),
        get_sample(aGetSample),
        select_neighborhood(aSelectNeighborhood) {}

  template <typename Graph, typename RRGVisitor, typename PositionMap>
  auto operator()(Graph& g, RRGVisitor vis, PositionMap g_position) const {
    using Vertex = graph_vertex_t<Graph>;
    using NodePuller = detail::rrg_node_puller<Graph>;
    using std::back_inserter;

    std::size_t i = 0;
    while (true) {
      auto p_new = get_sample(*space);

      std::vector<Vertex> Pred;
      std::vector<Vertex> Succ;
      if constexpr (boost::is_undirected_graph<Graph>::value) {
        select_neighborhood(p_new, back_inserter(Pred), g, *space, g_position);
      } else {
        select_neighborhood(p_new, back_inserter(Pred), back_inserter(Succ), g,
                            *space, g_position);
      }

      auto [x_near, was_expanded, ep] =
          NodePuller::expand_to_nearest(p_new, Pred, g, vis);
      if (was_expanded) {
        return std::tuple(x_near, p_new, ep);
      }
      if (i >= 10) {
        return std::tuple(boost::graph_traits<Graph>::null_vertex(), p_new, ep);
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

  rrg_bidir_generator(const Topology* aSpace, RandomSampler aGetSample,
                      NcSelector aSelectNeighborhood, PredecessorMap aPred,
                      SuccessorMap aSucc)
      : space(aSpace),
        get_sample(aGetSample),
        select_neighborhood(aSelectNeighborhood),
        predecessor(aPred),
        successor(aSucc) {}

  template <typename Graph, typename RRGVisitor, typename PositionMap>
  auto operator()(Graph& g, RRGVisitor vis, PositionMap g_position) const {
    using Vertex = graph_vertex_t<Graph>;
    using NodePuller = detail::rrg_node_puller<Graph>;
    using std::back_inserter;

    const Vertex null_v = boost::graph_traits<Graph>::null_vertex();

    int i = 0;
    while (true) {
      auto p_pred = get_sample(*space);
      auto p_succ = p_pred;

      std::vector<Vertex> Pred;
      std::vector<Vertex> Succ;
      if constexpr (boost::is_undirected_graph<Graph>::value) {
        select_neighborhood(p_pred, back_inserter(Pred), g, *space, g_position);
      } else {
        select_neighborhood(p_pred, back_inserter(Pred), back_inserter(Succ), g,
                            *space, g_position);
      }

      auto [x_pred, was_expanded, ep_pred] =
          NodePuller::expand_to_nearest(p_pred, Pred, g, vis, predecessor);
      if constexpr (boost::is_undirected_graph<Graph>::value) {
        Pred.swap(Succ);
      }
      auto [x_succ, was_retracted, ep_succ] =
          NodePuller::retract_from_nearest(p_succ, Succ, g, vis, successor);
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

#endif

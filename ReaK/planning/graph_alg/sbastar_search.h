/**
 * \file sbastar_search.h
 *
 * This library provides function templates and concepts that implement a Sampling-based A* search
 * algorithm. A SBA* uses the A* search algorithm to drive the expansion of a roadmap into the free-space
 * in order to connect a start and goal location. This algorithm has many customization points because there
 * are many choices to be made in the method, such as how to find nearest neighbors for attempting to
 * connect them through free-space, how to expand vertices, when to stop the algorithm, etc.
 * All these customization points are left to the user to implement, some are defined by the
 * SBAStarVisitorConcept (random-walk, edge-added, etc.).
 *
 * The SBA* algorithm is a generalization of the A* algorithm where the neighborhood of a given node of
 * the motion graph is not defined as a fixed set of neighbors (as in a classic A* over a fixed graph),
 * but rather as a region from which samples can be drawn (biased or not). In an ordinary A* algorithm,
 * vertices are closed when their entire neighborhood has been explored. In an SBA* algorithm, the same
 * criteria cannot apply since samples could be drawn ad infinitum, so, instead, this concept of the
 * neighborhood being fully explored is derived from the expected information gained (or conversely, the
 * "surprisal") from drawing a new sample in the neighborhood.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2013
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

#ifndef REAK_PLANNING_GRAPH_ALG_SBASTAR_SEARCH_H_
#define REAK_PLANNING_GRAPH_ALG_SBASTAR_SEARCH_H_

#include <limits>
#include <tuple>
#include <utility>

#include "ReaK/topologies/spaces/metric_space_concept.h"

#include "adj_list_tree_overlay.h"
#include "bagl/d_ary_heap.h"
#include "bagl/graph_concepts.h"
#include "bagl/more_property_maps.h"
#include "bagl/properties.h"

#include "ReaK/planning/graph_alg/branch_and_bound_connector.h"
#include "ReaK/planning/graph_alg/lazy_connector.h"
#include "ReaK/planning/graph_alg/pruned_connector.h"
#include "ReaK/planning/graph_alg/sbmp_visitor_concepts.h"

namespace ReaK::graph {

/**
  * This concept class defines the valid expressions required of a class to be used as a visitor
  * class for the SBA* algorithm. A visitor class is essentially a class that regroups a number of
  * callback functions that can be used to inject customization into the SBA* algorithm. In other
  * words, the visitor pattern in generic programming is an implementation of IoC
  * (Inversion of Control), since the SBA* algorithm is in control of execution, but custom behavior can
  * be injected in several places, even blocking the algorithm if needed.
  *
  * Required concepts:
  *
  * The visitor class should model SBMPVisitor, NodePushingVisitor, NodeReConnectVisitor,
  * NeighborhoodTrackingVisitor, and NodeExploringVisitor.
  *
  * Valid expressions:
  *
  * vis.publish_path(g);  A function to notify the visitor that at least one A* round has completed and its resulting
  *path (partial or complete) can be published (the path is encoded in the predecessor property-map).
  */
template <typename Visitor, typename Graph, typename Space>
concept SBAStarVisitor =
    SBMPVisitor<Visitor, Graph>&& NodePushingVisitor<Visitor, Graph, Space>&&
        NodeReConnectVisitor<Visitor, Graph>&& NeighborhoodTrackingVisitor<
            Visitor, Graph>&& NodeExploringVisitor<Visitor, Graph>&&
        requires(Visitor vis, Graph g) {
  vis.publish_path(g);
};

/**
 * This class is simply an archetype visitor for the SBA* algorithm.
 */
template <typename Space>
struct sbastar_visitor_archetype : sbmp_visitor_archetype,
                                   node_pushing_visitor_archetype<Space>,
                                   node_reconnect_visitor_archetype,
                                   neighborhood_tracking_visitor_archetype,
                                   node_exploring_visitor_archetype {
  template <typename Graph>
  void publish_path(const Graph& /*unused*/) const {}
};

/**
  * This concept class defines the valid expressions required of a class to be used as a visitor
  * class for the bi-directional SBA* algorithm. A visitor class is essentially a class that regroups a number of
  * callback functions that can be used to inject customization into the bi-directional SBA* algorithm.
  * In other words, the visitor pattern in generic programming is an implementation of IoC
  * (Inversion of Control), since the bi-directional SBA* algorithm is in control of execution,
  * but custom behavior can be injected in several places, even blocking the algorithm if needed.
  *
  * Required concepts:
  *
  * The visitor class should model SBAStarVisitor and NodeBackPushingVisitor.
  */
template <typename Visitor, typename Graph, typename Space>
concept SBAStarBidirVisitor = SBAStarVisitor<Visitor, Graph, Space>&&
    NodeBackPushingVisitor<Visitor, Graph, Space>;

/**
 * This class is simply an archetype visitor for the bi-directional SBA* algorithm.
 */
template <typename Space>
struct sbastar_bidir_visitor_archetype
    : sbastar_visitor_archetype<Space>,
      node_back_pushing_visitor_archetype<Space> {};

namespace sbastar_detail {

template <typename Graph, typename UniformCostVisitor, typename UpdatableQueue,
          bagl::concepts::ReadWriteVertexPropertyMap<Graph> IndexInHeapMap,
          bagl::concepts::ReadWriteVertexPropertyMap<Graph> KeyMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
          bagl::concepts::ReadableVPropMemberMap<Graph> DensityMap,
          bagl::concepts::ReadableVPropMemberMap<Graph> ConstrictionMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> FwdDistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> SuccessorMap =
              null_vertex_prop_map<Graph>>
struct sbastar_bfs_visitor {

  sbastar_bfs_visitor(UniformCostVisitor vis, UpdatableQueue& q,
                      IndexInHeapMap index_in_heap, KeyMap key, PositionMap pos,
                      WeightMap weight, DensityMap density,
                      ConstrictionMap constriction, DistanceMap dist,
                      PredecessorMap pred, FwdDistanceMap fwd_dist,
                      SuccessorMap succ = SuccessorMap())
      : vis_(vis),
        q_(q),
        index_in_heap_(index_in_heap),
        key_(key),
        position_(pos),
        weight_(weight),
        density_(density),
        constriction_(constriction),
        distance_(dist),
        predecessor_(pred),
        fwd_distance_(fwd_dist),
        successor_(succ) {}

  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;
  using Edge = bagl::graph_edge_descriptor_t<Graph>;

  template <typename PositionValue>
  Vertex create_vertex(const PositionValue& p, Graph& g) const {
    bagl::vertex_bundle_type<Graph> up;
    put(position_, up, p);
    put(distance_, up, std::numeric_limits<double>::infinity());
    put(predecessor_, up, bagl::graph_traits<Graph>::null_vertex());
    put(fwd_distance_, up, std::numeric_limits<double>::infinity());
    put(successor_, up, bagl::graph_traits<Graph>::null_vertex());
    Vertex u = add_vertex(g, std::move(up));
    vis_.vertex_added(u, g);
    put(index_in_heap_, u, std::numeric_limits<std::size_t>::max());
    put(key_, u, 0.0);
    return u;
  }

  template <typename PositionValue>
  auto steer_towards_position(const PositionValue& p, Vertex u,
                              const Graph& g) const {
    return vis_.steer_towards_position(p, u, g);
  }

  template <typename PositionValue>
  auto steer_back_to_position(const PositionValue& p, Vertex u,
                              const Graph& g) const {
    return vis_.steer_back_to_position(p, u, g);
  }

  auto can_be_connected(Vertex u, Vertex v, Graph& g) const {
    return vis_.can_be_connected(u, v, g);
  }

  auto random_walk(Vertex u, Graph& g) const { return vis_.random_walk(u, g); }

  auto random_back_walk(Vertex u, Graph& g) const {
    return vis_.random_back_walk(u, g);
  }

  void vertex_to_be_removed(Vertex u, Graph& g) const {
    put(key_, u, -std::numeric_limits<double>::infinity());
    q_.push_or_update(u);
    q_.pop();
    vis_.vertex_to_be_removed(u, g);
  }

  void vertex_added(Vertex v, Graph& g) const { vis_.vertex_added(v, g); }

  void edge_added(Edge e, Graph& g) const {
    vis_.edge_added(e, g);
    vis_.examine_edge(e, g);
  }

  void travel_explored(Vertex u, Vertex v, Graph& g) const {
    vis_.travel_explored(u, v, g);
  }

  void travel_succeeded(Vertex u, Vertex v, Graph& g) const {
    vis_.travel_succeeded(u, v, g);
  }

  void travel_failed(Vertex u, Vertex v, Graph& g) const {
    vis_.travel_failed(u, v, g);
  }

  void requeue_vertex(Vertex u, Graph& g) const {
    update_key(u, g);
    if (!vis_.should_close(u, g)) {
      q_.push_or_update(u);
      vis_.discover_vertex(u, g);
    }
  }
  void affected_vertex(Vertex u, Graph& g) const {
    // same function, different name.
    requeue_vertex(u, g);
  }

  void examine_vertex(Vertex u, Graph& g) const { vis_.examine_vertex(u, g); }

  void examine_edge(Edge e, Graph& g) const { vis_.examine_edge(e, g); }

  void publish_path(const Graph& g) const { vis_.publish_path(g); }
  bool keep_going() const { return vis_.keep_going(); }
  bool has_search_potential(Vertex u, const Graph& g) const {
    return vis_.has_search_potential(u, g);
  }
  bool should_close(Vertex u, const Graph& g) const {
    return vis_.should_close(u, g);
  }

  void update_key(Vertex u, Graph& g) const {
    vis_.affected_vertex(u, g);
    double g_u = get(distance_, get_property(g, u));
    double h_u = get(fwd_distance_, get_property(g, u));
    // Key-value for the min-heap (priority-queue):
    put(key_, u,
        (g_u + h_u) / (1.0 - get(constriction_, get_property(g, u))) /
            (1.0 - get(density_, get_property(g, u))));
  }

  UniformCostVisitor vis_;
  UpdatableQueue& q_;
  IndexInHeapMap index_in_heap_;

  KeyMap key_;
  PositionMap position_;
  WeightMap weight_;
  DensityMap density_;
  ConstrictionMap constriction_;
  DistanceMap distance_;
  PredecessorMap predecessor_;
  FwdDistanceMap fwd_distance_;
  SuccessorMap successor_;
};

template <typename Graph, typename V, typename Q, typename... Maps>
sbastar_bfs_visitor<Graph, V, Q, Maps...> make_sbastar_bfs_visitor(V v, Q& q,
                                                                   Maps... m) {
  return {v, q, m...};
}

struct sba_node_generator {

  template <typename Graph, typename SBAVisitor, typename PositionMap>
  auto operator()(bagl::graph_vertex_descriptor_t<Graph> u, Graph& g,
                  const SBAVisitor& sba_vis, PositionMap /*unused*/) const {
    auto [p_new, was_expanded, ep_new] = sba_vis.random_walk(u, g);
    auto u_exp = (was_expanded ? u : bagl::graph_traits<Graph>::null_vertex());
    return std::tuple(u_exp, p_new, ep_new);
  }
};

struct sba_bidir_node_generator {

  template <typename Graph, typename SBAVisitor, typename PositionMap>
  auto operator()(bagl::graph_vertex_descriptor_t<Graph> u, Graph& g,
                  const SBAVisitor& sba_vis, PositionMap /*unused*/) const {
    using std::get;
    using PositionValue =
        std::decay_t<decltype(get<0>(sba_vis.random_walk(u, g)))>;
    using EdgeProp = bagl::edge_property_type<Graph>;

    PositionValue p_exp;
    PositionValue p_ret;
    bool was_expanded = false;
    bool was_retracted = false;
    EdgeProp ep_exp;
    EdgeProp ep_ret;
    if (get(sba_vis.predecessor_, get_property(g, u)) !=
        bagl::graph_traits<Graph>::null_vertex()) {
      std::tie(p_exp, was_expanded, ep_exp) = sba_vis.random_walk(u, g);
    }
    if (get(sba_vis.successor_, get_property(g, u)) !=
        bagl::graph_traits<Graph>::null_vertex()) {
      std::tie(p_ret, was_retracted, ep_ret) = sba_vis.random_back_walk(u, g);
    }

    auto u_exp = (was_expanded ? u : bagl::graph_traits<Graph>::null_vertex());
    auto u_ret = (was_retracted ? u : bagl::graph_traits<Graph>::null_vertex());
    return std::tuple(u_exp, p_exp, ep_exp, u_ret, p_ret, ep_ret);
  }
};

template <typename Graph, pp::MetricSpace Space, typename Visitor,
          typename MotionGraphConnector, typename SBANodeGenerator,
          typename MutableQueue, typename NcSelector>
void sbastar_search_loop(Graph& g, const Space& super_space, Visitor& sba_vis,
                         MotionGraphConnector connect_vertex,
                         SBANodeGenerator sba_generate_node, MutableQueue& q,
                         NcSelector select_neighborhood) {
  while (!q.empty() && sba_vis.keep_going()) {
    auto u = q.top();
    q.pop();

    // stop if the best nodes do not meet the potential threshold.
    while (!sba_vis.has_search_potential(u, g)) {
      if (q.empty()) {
        u = bagl::graph_traits<Graph>::null_vertex();
        break;
      }
      u = q.top();
      q.pop();
    }
    if (u == bagl::graph_traits<Graph>::null_vertex()) {
      break;  // no more nodes with search potential.
    }

    sba_vis.examine_vertex(u, g);

    auto [x_near, p_new, eprop] =
        sba_generate_node(u, g, sba_vis, sba_vis.position_);

    // then push it back on the OPEN queue.
    if ((x_near != bagl::graph_traits<Graph>::null_vertex()) || (q.empty())) {
      sba_vis.requeue_vertex(u, g);
    }

    if (x_near != bagl::graph_traits<Graph>::null_vertex()) {
      connect_vertex(p_new, x_near, eprop, g, super_space, sba_vis,
                     sba_vis.position_, sba_vis.distance_, sba_vis.predecessor_,
                     sba_vis.weight_, select_neighborhood);
    }

  }  // while
}

template <typename Graph, pp::MetricSpace Space, typename Visitor,
          typename MotionGraphConnector, typename SBANodeGenerator,
          typename MutableQueue, typename NcSelector>
void sbastar_bidir_loop(Graph& g, const Space& super_space, Visitor& sba_vis,
                        MotionGraphConnector connect_vertex,
                        SBANodeGenerator sba_generate_node, MutableQueue& q,
                        NcSelector select_neighborhood) {
  while (!q.empty() && sba_vis.keep_going()) {
    auto u = q.top();
    q.pop();

    // stop if the best node does not meet the potential threshold.
    while (!sba_vis.has_search_potential(u, g)) {
      if (q.empty()) {
        u = bagl::graph_traits<Graph>::null_vertex();
        break;
      }
      u = q.top();
      q.pop();
    }
    if (u == bagl::graph_traits<Graph>::null_vertex()) {
      break;  // no more nodes with search potential.
    }

    sba_vis.examine_vertex(u, g);

    auto [x_near_pred, p_new_pred, ep_pred, x_near_succ, p_new_succ, ep_succ] =
        sba_generate_node(u, g, sba_vis, sba_vis.position_);

    // then push it back on the OPEN queue.
    if ((x_near_pred != bagl::graph_traits<Graph>::null_vertex()) ||
        (x_near_succ != bagl::graph_traits<Graph>::null_vertex()) ||
        (q.empty())) {
      sba_vis.requeue_vertex(u, g);
    }

    if (x_near_pred != bagl::graph_traits<Graph>::null_vertex()) {
      auto x_near_other = bagl::graph_traits<Graph>::null_vertex();
      bagl::edge_property_type<Graph> ep_other;
      connect_vertex(p_new_pred, x_near_pred, ep_pred, x_near_other, ep_other,
                     g, super_space, sba_vis, sba_vis.position_,
                     sba_vis.distance_, sba_vis.predecessor_,
                     sba_vis.fwd_distance_, sba_vis.successor_, sba_vis.weight_,
                     select_neighborhood);
    }
    if (x_near_succ != bagl::graph_traits<Graph>::null_vertex()) {
      auto x_near_other = bagl::graph_traits<Graph>::null_vertex();
      bagl::edge_property_type<Graph> ep_other;
      connect_vertex(p_new_succ, x_near_other, ep_other, x_near_succ, ep_succ,
                     g, super_space, sba_vis, sba_vis.position_,
                     sba_vis.distance_, sba_vis.predecessor_,
                     sba_vis.fwd_distance_, sba_vis.successor_, sba_vis.weight_,
                     select_neighborhood);
    }

  }  // while
}

template <bagl::concepts::VertexListGraph Graph, pp::MetricSpace Space,
          SBAStarVisitor<Graph, Space> Visitor, typename NodeConnector,
          typename KeyMap, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, typename DistanceMap,
          typename PredecessorMap, typename FwdDistanceMap, typename NcSelector>
void generate_sbastar_no_init_impl(
    Graph& g, bagl::graph_vertex_descriptor_t<Graph> start_vertex,
    const Space& super_space, Visitor vis, NodeConnector connect_vertex,
    KeyMap key, PositionMap position, WeightMap weight, DensityMap density,
    ConstrictionMap constriction, DistanceMap distance,
    PredecessorMap predecessor, FwdDistanceMap fwd_distance,
    NcSelector select_neighborhood) {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;

  auto index_in_heap =
      bagl::vector_property_map(num_vertices(g), get(bagl::vertex_index, g),
                                std::numeric_limits<std::size_t>::max());

  // priority queue holding the OPEN set.
  using KeyCompareType = std::less<>;  // <---- this is a min-heap.
  auto q = bagl::make_d_ary_heap_indirect<Vertex, 4>(key, index_in_heap.ref(),
                                                     KeyCompareType());

  auto sba_bfs_vis = make_sbastar_bfs_visitor<Graph>(
      vis, q, index_in_heap.ref(), key, position, weight, density, constriction,
      distance, predecessor, fwd_distance);

  put(distance, get_property(g, start_vertex), 0.0);
  put(predecessor, get_property(g, start_vertex), start_vertex);
  sba_bfs_vis.requeue_vertex(start_vertex, g);

  sbastar_search_loop(g, super_space, sba_bfs_vis, connect_vertex,
                      sba_node_generator(), q, select_neighborhood);
}

template <bagl::concepts::VertexListGraph Graph, pp::MetricSpace Space,
          SBAStarBidirVisitor<Graph, Space> Visitor, typename NodeConnector,
          typename KeyMap, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, typename DistanceMap,
          typename PredecessorMap, typename FwdDistanceMap,
          typename SuccessorMap, typename NcSelector>
void generate_sbastar_bidir_no_init_impl(
    Graph& g, bagl::graph_vertex_descriptor_t<Graph> start_vertex,
    bagl::graph_vertex_descriptor_t<Graph> goal_vertex,
    const Space& super_space, Visitor vis, NodeConnector connect_vertex,
    KeyMap key, PositionMap position, WeightMap weight, DensityMap density,
    ConstrictionMap constriction, DistanceMap distance,
    PredecessorMap predecessor, FwdDistanceMap fwd_distance,
    SuccessorMap successor, NcSelector select_neighborhood) {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;

  auto index_in_heap =
      bagl::vector_property_map(num_vertices(g), get(bagl::vertex_index, g),
                                std::numeric_limits<std::size_t>::max());

  // priority queue holding the OPEN set.
  using KeyCompareType = std::less<>;  // <---- this is a min-heap.
  auto q = bagl::make_d_ary_heap_indirect<Vertex, 4>(key, index_in_heap.ref(),
                                                     KeyCompareType());

  auto sba_bfs_vis = make_sbastar_bfs_visitor<Graph>(
      vis, q, index_in_heap.ref(), key, position, weight, density, constriction,
      distance, predecessor, fwd_distance, successor);

  put(distance, get_property(g, start_vertex), 0.0);
  put(predecessor, get_property(g, start_vertex), start_vertex);
  sba_bfs_vis.requeue_vertex(start_vertex, g);
  if (goal_vertex != bagl::graph_traits<Graph>::null_vertex()) {
    put(fwd_distance, get_property(g, goal_vertex), 0.0);
    put(successor, get_property(g, goal_vertex), goal_vertex);
    sba_bfs_vis.requeue_vertex(goal_vertex, g);
  }

  sbastar_bidir_loop(g, super_space, sba_bfs_vis, connect_vertex,
                     sba_bidir_node_generator(), q, select_neighborhood);
}

template <bagl::concepts::VertexListGraph Graph, typename Visitor,
          typename KeyMap, typename DistanceMap, typename PredecessorMap>
void initialize_sbastar_nodes(Graph& g, Visitor vis, KeyMap key,
                              DistanceMap distance,
                              PredecessorMap predecessor) {
  for (auto u : vertices(g)) {
    put(key, u, 0.0);
    put(distance, get_property(g, u), std::numeric_limits<double>::infinity());
    put(predecessor, get_property(g, u),
        bagl::graph_traits<Graph>::null_vertex());
    vis.initialize_vertex(u, g);
  }
}

template <bagl::concepts::VertexListGraph Graph, typename Visitor,
          typename KeyMap, typename DistanceMap, typename PredecessorMap,
          typename FwdDistanceMap, typename SuccessorMap>
void initialize_sbastar_nodes(Graph& g, Visitor vis, KeyMap key,
                              DistanceMap distance, PredecessorMap predecessor,
                              FwdDistanceMap fwd_distance,
                              SuccessorMap successor) {
  for (auto u : vertices(g)) {
    put(key, u, 0.0);
    put(distance, get_property(g, u), std::numeric_limits<double>::infinity());
    put(predecessor, get_property(g, u),
        bagl::graph_traits<Graph>::null_vertex());
    put(fwd_distance, get_property(g, u),
        std::numeric_limits<double>::infinity());
    put(successor, get_property(g, u),
        bagl::graph_traits<Graph>::null_vertex());
    vis.initialize_vertex(u, g);
  }
}

}  // namespace sbastar_detail

template <typename Graph, pp::MetricSpace Space,
          SBAStarVisitor<Graph, Space> Visitor, typename NcSelector,
          typename KeyMap, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, typename DistanceMap,
          typename PredecessorMap, typename FwdDistanceMap,
          typename SuccessorMap = null_vertex_prop_map<Graph>>
struct sbastar_bundle {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;

  Graph* g_;
  Vertex start_vertex_;
  Vertex goal_vertex_;
  const Space* super_space_;
  Visitor vis_;
  NcSelector select_neighborhood_;
  KeyMap key_;
  PositionMap position_;
  WeightMap weight_;
  DensityMap density_;
  ConstrictionMap constriction_;
  DistanceMap distance_;
  PredecessorMap predecessor_;
  FwdDistanceMap fwd_distance_;
  SuccessorMap successor_;

  sbastar_bundle(Graph& g, Vertex start_vertex, const Space& super_space,
                 Visitor vis, NcSelector select_neighborhood, KeyMap key,
                 PositionMap position, WeightMap weight, DensityMap density,
                 ConstrictionMap constriction, DistanceMap distance,
                 PredecessorMap predecessor, FwdDistanceMap fwd_distance,
                 SuccessorMap successor = SuccessorMap())
      : g_(&g),
        start_vertex_(start_vertex),
        goal_vertex_(bagl::graph_traits<Graph>::null_vertex()),
        super_space_(&super_space),
        vis_(vis),
        select_neighborhood_(select_neighborhood),
        key_(key),
        position_(position),
        weight_(weight),
        density_(density),
        constriction_(constriction),
        distance_(distance),
        predecessor_(predecessor),
        fwd_distance_(fwd_distance),
        successor_(successor) {}

  sbastar_bundle(Graph& g, Vertex start_vertex, Vertex goal_vertex,
                 const Space& super_space, Visitor vis,
                 NcSelector select_neighborhood, KeyMap key,
                 PositionMap position, WeightMap weight, DensityMap density,
                 ConstrictionMap constriction, DistanceMap distance,
                 PredecessorMap predecessor, FwdDistanceMap fwd_distance,
                 SuccessorMap successor = SuccessorMap())
      : g_(&g),
        start_vertex_(start_vertex),
        goal_vertex_(goal_vertex),
        super_space_(&super_space),
        vis_(vis),
        select_neighborhood_(select_neighborhood),
        key_(key),
        position_(position),
        weight_(weight),
        density_(density),
        constriction_(constriction),
        distance_(distance),
        predecessor_(predecessor),
        fwd_distance_(fwd_distance),
        successor_(successor) {}
};

/**
  * This function template creates a bundle of parameters to be fed to any of the
  * SBA* algorithms. This is mainly to simply the interface and the code of all these
  * different variants of the SBA* algorithm.
  * \tparam Graph The graph type that can store the generated roadmap, should model
  *         BidirectionalGraph and MutableGraph.
  * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
  * \tparam Space The topology type that represents the free-space, should model BGL's Topology concept.
  * \tparam Visitor The type of the SBA* visitor to be used.
  * \tparam AStarHeuristicMap This property-map type is used to obtain the heuristic-function values
  *         for each vertex in the graph.
  * \tparam PositionMap A property-map type that can store the position of each vertex-property object.
  * \tparam WeightMap This property-map type is used to store the weights of the edge-properties of the
  *         graph (cost of travel along an edge).
  * \tparam DensityMap A property-map type that can store the probability-measure of the expected common information
  *         between a new sample and the current neighborhood for each vertex-property object.
  * \tparam ConstrictionMap A property-map type that can store the probability-measure of sampling a colliding point
  *         for each vertex-property object.
  * \tparam DistanceMap This property-map type is used to store the estimated distance of each vertex-property object
  *         to the goal.
  * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting
  *         vertex-property object together with its optimal predecessor.
  * \tparam KeyMap This property-map type is used to store the priority-keys of the vertices of the
  *         graph (cost of travel along an edge).
  * \tparam NcSelector A functor type that can select a list of vertices of the graph that are
  *         the nearest-neighbors of a given vertex (or some other heuristic to select the neighbors).
  *         See classes in the topological_search.hpp header-file.
  *
  * \param g A mutable graph that should initially store the starting
  *        vertex (if not it will be randomly generated) and will store
  *        the generated graph once the algorithm has finished.
  * \param start_vertex The starting point of the algorithm, on the graph.
  * \param super_space A topology (as defined by the Boost Graph Library). This topology
  *        should not include collision checking in its distance metric.
  * \param vis A SBA* visitor implementing the FADPRMVisitorConcept. This is the
  *        main point of customization and recording of results that the
  *        user can implement.
  * \param select_neighborhood A callable object (functor) that can select a list of
  *        vertices of the graph that ought to be connected to a new
  *        vertex. The list should be sorted in order of increasing "distance".
  * \param key The property-map which stores the AD* key-values associated to each vertex.
  * \param position A mapping that implements the MutablePropertyMap Concept. Also,
  *        the value_type of this map should be the same type as the topology's
  *        value_type.
  * \param weight The property-map which stores the weight of each edge-property object (the cost of travel
  *        along the edge).
  * \param density A property-map that provides the expected common information associated with a sample drawn near
  *        to a vertex w.r.t. the current neighborhood of that vertex.
  * \param constriction A property-map that provides the probability of a collision when a sample is drawn near to a
  *        vertex (i.e., that a sample near this vertex will not be in the free-space).
  * \param distance The property-map which stores the estimated distance of each vertex to the goal.
  * \param predecessor The property-map which will store the resulting path by connecting
  *        vertices together with their optimal predecessor (follow in reverse to discover the
  *        complete path).
  * \param fwd_distance The property-map which stores the estimated distance of each vertex to the goal.
  */
template <bagl::concepts::VertexListGraph Graph, pp::MetricSpace Space,
          SBAStarVisitor<Graph, Space> Visitor, typename NcSelector,
          typename KeyMap, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, typename DistanceMap,
          typename PredecessorMap, typename FwdDistanceMap>
requires bagl::concepts::MutableGraph<Graph> sbastar_bundle<
    Graph, Space, Visitor, NcSelector, KeyMap, PositionMap, WeightMap,
    DensityMap, ConstrictionMap, DistanceMap, PredecessorMap, FwdDistanceMap>
make_sbastar_bundle(Graph& g,
                    bagl::graph_vertex_descriptor_t<Graph> start_vertex,
                    const Space& super_space, Visitor vis,
                    NcSelector select_neighborhood, KeyMap key,
                    PositionMap position, WeightMap weight, DensityMap density,
                    ConstrictionMap constriction, DistanceMap distance,
                    PredecessorMap predecessor, FwdDistanceMap fwd_distance) {
  return {g,        start_vertex, super_space, vis,     select_neighborhood,
          key,      position,     weight,      density, constriction,
          distance, predecessor,  fwd_distance};
}

template <bagl::concepts::VertexListGraph Graph, pp::MetricSpace Space,
          SBAStarVisitor<Graph, Space> Visitor, typename NcSelector,
          typename KeyMap, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, typename DistanceMap,
          typename PredecessorMap, typename FwdDistanceMap>
requires bagl::concepts::MutableGraph<Graph> sbastar_bundle<
    Graph, Space, Visitor, NcSelector, KeyMap, PositionMap, WeightMap,
    DensityMap, ConstrictionMap, DistanceMap, PredecessorMap, FwdDistanceMap>
make_sbastar_bundle(Graph& g,
                    bagl::graph_vertex_descriptor_t<Graph> start_vertex,
                    bagl::graph_vertex_descriptor_t<Graph> goal_vertex,
                    const Space& super_space, Visitor vis,
                    NcSelector select_neighborhood, KeyMap key,
                    PositionMap position, WeightMap weight, DensityMap density,
                    ConstrictionMap constriction, DistanceMap distance,
                    PredecessorMap predecessor, FwdDistanceMap fwd_distance) {
  return {g,           start_vertex, goal_vertex,
          super_space, vis,          select_neighborhood,
          key,         position,     weight,
          density,     constriction, distance,
          predecessor, fwd_distance};
}

/**
  * This function template creates a bundle of parameters to be fed to any of the
  * SBA* algorithms. This is mainly to simply the interface and the code of all these
  * different variants of the SBA* algorithm.
  * \tparam Graph The graph type that can store the generated roadmap, should model
  *         BidirectionalGraph and MutableGraph.
  * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
  * \tparam Space The topology type that represents the free-space, should model BGL's Topology concept.
  * \tparam Visitor The type of the SBA* visitor to be used.
  * \tparam AStarHeuristicMap This property-map type is used to obtain the heuristic-function values
  *         for each vertex in the graph.
  * \tparam PositionMap A property-map type that can store the position of each vertex-property object.
  * \tparam WeightMap This property-map type is used to store the weights of the edge-properties of the
  *         graph (cost of travel along an edge).
  * \tparam DensityMap A property-map type that can store the probability-measure of the expected common information
  *         between a new sample and the current neighborhood for each vertex-property object.
  * \tparam ConstrictionMap A property-map type that can store the probability-measure of sampling a colliding point
  *         for each vertex-property object.
  * \tparam DistanceMap This property-map type is used to store the estimated distance of each vertex-property object
  *         to the goal.
  * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting
  *         vertex-property object together with its optimal predecessor.
  * \tparam KeyMap This property-map type is used to store the priority-keys of the vertices of the
  *         graph (cost of travel along an edge).
  * \tparam NcSelector A functor type that can select a list of vertices of the graph that are
  *         the nearest-neighbors of a given vertex (or some other heuristic to select the neighbors).
  *         See classes in the topological_search.hpp header-file.
  *
  * \param g A mutable graph that should initially store the starting
  *        vertex (if not it will be randomly generated) and will store
  *        the generated graph once the algorithm has finished.
  * \param start_vertex The starting point of the algorithm, on the graph.
  * \param goal_vertex The goal point of the algorithm, on the graph.
  * \param super_space A topology (as defined by the Boost Graph Library). This topology
  *        should not include collision checking in its distance metric.
  * \param vis A SBA* visitor implementing the FADPRMVisitorConcept. This is the
  *        main point of customization and recording of results that the
  *        user can implement.
  * \param select_neighborhood A callable object (functor) that can select a list of
  *        vertices of the graph that ought to be connected to a new
  *        vertex. The list should be sorted in order of increasing "distance".
  * \param key The property-map which stores the AD* key-values associated to each vertex.
  * \param position A mapping that implements the MutablePropertyMap Concept. Also,
  *        the value_type of this map should be the same type as the topology's
  *        value_type.
  * \param weight The property-map which stores the weight of each edge-property object (the cost of travel
  *        along the edge).
  * \param density A property-map that provides the expected common information associated with a sample drawn near
  *        to a vertex w.r.t. the current neighborhood of that vertex.
  * \param constriction A property-map that provides the probability of a collision when a sample is drawn near to a
  *        vertex (i.e., that a sample near this vertex will not be in the free-space).
  * \param distance The property-map which stores the estimated distance of each vertex from the start.
  * \param predecessor The property-map which will store the resulting path by connecting
  *        vertices together with their optimal predecessor (follow in reverse to discover the
  *        complete path).
  * \param fwd_distance The property-map which stores the estimated distance of each vertex to the goal.
  * \param successor The property-map which will store the resulting path by connecting
  *        vertices together with their optimal successor (follow in order to discover the
  *        remaining path to the goal).
  */

template <bagl::concepts::VertexListGraph Graph, pp::MetricSpace Space,
          SBAStarBidirVisitor<Graph, Space> Visitor, typename NcSelector,
          typename KeyMap, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, typename DistanceMap,
          typename PredecessorMap, typename FwdDistanceMap,
          typename SuccessorMap>
requires bagl::concepts::MutableGraph<Graph>
    sbastar_bundle<Graph, Space, Visitor, NcSelector, KeyMap, PositionMap,
                   WeightMap, DensityMap, ConstrictionMap, DistanceMap,
                   PredecessorMap, FwdDistanceMap, SuccessorMap>
    make_sbastar_bundle(Graph& g,
                        bagl::graph_vertex_descriptor_t<Graph> start_vertex,
                        bagl::graph_vertex_descriptor_t<Graph> goal_vertex,
                        const Space& super_space, Visitor vis,
                        NcSelector select_neighborhood, KeyMap key,
                        PositionMap position, WeightMap weight,
                        DensityMap density, ConstrictionMap constriction,
                        DistanceMap distance, PredecessorMap predecessor,
                        FwdDistanceMap fwd_distance, SuccessorMap successor) {
  return {g,           start_vertex, goal_vertex,
          super_space, vis,          select_neighborhood,
          key,         position,     weight,
          density,     constriction, distance,
          predecessor, fwd_distance, successor};
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the SBA* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  */
template <typename SBAStarBundle>
void generate_sbastar_no_init(const SBAStarBundle& bdl) {
  sbastar_detail::generate_sbastar_no_init_impl(
      *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
      pruned_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
      bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
      bdl.fwd_distance_, bdl.select_neighborhood_);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the SBA* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  */
template <typename SBAStarBundle>
void generate_sbastar(const SBAStarBundle& bdl) {
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_);
  generate_sbastar_no_init(bdl);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional SBA* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  */
template <typename SBAStarBundle>
void generate_sbastar_bidir_no_init(const SBAStarBundle& bdl) {
  sbastar_detail::generate_sbastar_bidir_no_init_impl(
      *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_, bdl.vis_,
      pruned_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
      bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
      bdl.fwd_distance_, bdl.successor_, bdl.select_neighborhood_);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional SBA* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  */
template <typename SBAStarBundle>
void generate_sbastar_bidir(const SBAStarBundle& bdl) {
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_,
                                           bdl.fwd_distance_, bdl.successor_);

  generate_sbastar_bidir_no_init(bdl);
}

/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Lazy-SBA* algorithm, without initialization of the existing graph.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_sbastar_no_init(const SBAStarBundle& bdl) {
  sbastar_detail::generate_sbastar_no_init_impl(
      *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
      lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_, bdl.density_,
      bdl.constriction_, bdl.distance_, bdl.predecessor_, bdl.fwd_distance_,
      bdl.select_neighborhood_);
}

/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the search.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_sbastar(const SBAStarBundle& bdl) {
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_);
  generate_lazy_sbastar_no_init(bdl);
}

/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Bi-directional Lazy-SBA* algorithm, without initialization of the existing graph.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_sbastar_bidir_no_init(const SBAStarBundle& bdl) {
  sbastar_detail::generate_sbastar_bidir_no_init_impl(
      *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_, bdl.vis_,
      lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_, bdl.density_,
      bdl.constriction_, bdl.distance_, bdl.predecessor_, bdl.fwd_distance_,
      bdl.successor_, bdl.select_neighborhood_);
}

/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Bi-directional Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the search.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_sbastar_bidir(const SBAStarBundle& bdl) {
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_,
                                           bdl.fwd_distance_, bdl.successor_);

  generate_lazy_sbastar_bidir_no_init(bdl);
}

/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Lazy-SBA* algorithm, without initialization of the existing graph.
 * This function uses a branch-and-bound heuristic to limit the number of nodes.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_bnb_sbastar_no_init(const SBAStarBundle& bdl) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  if (bdl.goal_vertex_ == bagl::graph_traits<Graph>::null_vertex()) {
    sbastar_detail::generate_sbastar_no_init_impl(
        *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
        lazy_node_connector(), bdl.key_, bdl.position_, bdl.weight_,
        bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
        bdl.fwd_distance_, bdl.select_neighborhood_);
  } else {
    bnb_ordering_data<Graph> bnb_data(*bdl.g_, bdl.start_vertex_,
                                      bdl.goal_vertex_);
    sbastar_detail::generate_sbastar_no_init_impl(
        *bdl.g_, bdl.start_vertex_, *bdl.super_space_, bdl.vis_,
        bnb_connector<Graph>(bnb_data), bdl.key_, bdl.position_, bdl.weight_,
        bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
        bdl.fwd_distance_, bdl.select_neighborhood_);
  }
}

/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the search.
 * This function uses a branch-and-bound heuristic to limit the number of nodes.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_bnb_sbastar(const SBAStarBundle& bdl) {
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_);
  generate_lazy_bnb_sbastar_no_init(bdl);
}

/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Bi-directional Lazy-SBA* algorithm, without initialization of the existing graph.
 * This function uses a branch-and-bound heuristic to limit the number of nodes.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_bnb_sbastar_bidir_no_init(const SBAStarBundle& bdl) {
  using Graph = std::decay_t<decltype(*bdl.g_)>;
  bnb_ordering_data<Graph> bnb_data(*bdl.g_, bdl.start_vertex_,
                                    bdl.goal_vertex_);
  sbastar_detail::generate_sbastar_bidir_no_init_impl(
      *bdl.g_, bdl.start_vertex_, bdl.goal_vertex_, *bdl.super_space_, bdl.vis_,
      bnb_connector<Graph>(bnb_data), bdl.key_, bdl.position_, bdl.weight_,
      bdl.density_, bdl.constriction_, bdl.distance_, bdl.predecessor_,
      bdl.fwd_distance_, bdl.successor_, bdl.select_neighborhood_);
}

/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Bi-directional Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the search.
 * This function uses a branch-and-bound heuristic to limit the number of nodes.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_bnb_sbastar_bidir(const SBAStarBundle& bdl) {
  sbastar_detail::initialize_sbastar_nodes(*bdl.g_, bdl.vis_, bdl.key_,
                                           bdl.distance_, bdl.predecessor_,
                                           bdl.fwd_distance_, bdl.successor_);
  generate_lazy_bnb_sbastar_bidir_no_init(bdl);
}

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_SBASTAR_SEARCH_H_

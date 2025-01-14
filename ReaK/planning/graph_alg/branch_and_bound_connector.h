/**
 * \file branch_and_bound_connector.h
 *
 * This library provides a class template and concept that implement a Lazy Branch-and-bound Motion-graph Connector.
 * A Lazy Branch-and-bound Connector uses the accumulated distance to assess the local optimality of the wirings
 * on a motion-graph and prunes away any node that cannot yield a better path than the current best path.
 * This algorithm has many customization points because it can be used in many different sampling-based
 * motion-planners.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2013
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

#ifndef REAK_PLANNING_GRAPH_ALG_BRANCH_AND_BOUND_CONNECTOR_H_
#define REAK_PLANNING_GRAPH_ALG_BRANCH_AND_BOUND_CONNECTOR_H_

#include <functional>
#include <limits>
#include <stack>
#include <tuple>
#include <type_traits>

#include "ReaK/topologies/spaces/metric_space_concept.h"

#include "ReaK/planning/graph_alg/lazy_connector.h"

#include "bagl/d_ary_heap.h"
#include "bagl/graph_concepts.h"
#include "bagl/more_property_maps.h"
#include "bagl/properties.h"

namespace ReaK::graph {

/**
  * This concept class defines the valid expressions required of a class to be used as a visitor
  * class for a branch-and-bound connection strategies. A visitor class is essentially a class that regroups a number of
  * callback functions that can be used to inject customization into a branch-and-bound connection strategy. In other
  * words, the visitor pattern in generic programming is an implementation of IoC
  * (Inversion of Control), since the connection strategy is in control of execution, but custom behavior can
  * be injected in several places, even blocking the algorithm if needed.
  *
  * Required concepts:
  *
  * The visitor class should model the MotionGraphConnectorVisitorConcept.
  *
  * Valid expressions:
  *
  * conn_vis.vertex_to_be_removed(u, g);  This function is called just before a vertex is removed from the motion-graph
  *(i.e., pruned by the branch-and-bound heuristic).
  */
template <typename Visitor, typename Graph, typename Space>
concept BNBConnectorVisitor =
    MotionGraphConnectorVisitor<Visitor, Graph, Space>&& requires(
        Visitor vis, Graph g, bagl::graph_vertex_descriptor_t<Graph> u) {
  vis.vertex_to_be_removed(u, g);
};

template <typename Graph>
struct bnb_ordering_data {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;
  using IndexInHeapMap = bagl::vector_property_map<std::size_t>;
  using KeyMap = bagl::vector_property_map<double>;
  using KeyCompareType = std::greater<>;  // <---- this is a max-heap.
  using MutableQueue =
      bagl::d_ary_heap_indirect<Vertex, 4,
                                bagl::property_map_ref<IndexInHeapMap>,
                                bagl::property_map_ref<KeyMap>, KeyCompareType>;

  Vertex start_vertex;
  Vertex goal_vertex;
  IndexInHeapMap index_in_heap;
  KeyMap key;
  MutableQueue q;  // priority queue holding the OPEN set.

  bnb_ordering_data(Graph& g, Vertex aStart, Vertex aGoal)
      : start_vertex(aStart),
        goal_vertex(aGoal),
        q(key.ref(), index_in_heap.ref(), KeyCompareType()) {
    for (auto u : vertices(g)) {
      put(index_in_heap, u, std::numeric_limits<std::size_t>::max());
    }
  }

  template <SBMPPruningVisitor<Graph> Visitor>
  void prune_worst_nodes(double threshold, Graph& g, const Visitor& conn_vis) {
    // prune all the worst nodes:
    double top_value = get(key, q.top());
    while (top_value > threshold) {
      conn_vis.vertex_to_be_removed(q.top(), g);
      clear_vertex(q.top(), g);
      remove_vertex(q.top(), g);
      q.pop();
      top_value = get(key, q.top());
    }
  }
};

/**
 * This callable class template implements a Lazy Branch-and-bound Motion-graph Connector.
 * A Lazy Branch-and-bound Connector uses the accumulated distance to assess the local optimality of the wirings
 * on a motion-graph and prunes away any node that cannot yield a better path than the current best path.
 * The call operator accepts a visitor object to provide customized behavior because it can be used in many
 * different sampling-based motion-planners. The visitor must model the BNBConnectorVisitorConcept concept.
 */
template <typename Graph>
struct bnb_connector {

  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;
  using EdgeProp = bagl::edge_property_type<Graph>;

  bnb_ordering_data<Graph>* data;

  explicit bnb_connector(bnb_ordering_data<Graph>& d) : data(&d){};

  template <pp::MetricSpace Space, BNBConnectorVisitor<Graph, Space> Visitor,
            typename PositionMap, typename DistanceMap, typename PredecessorMap,
            typename WeightMap>
  void update_successors(Vertex v, Graph& g, const Space& super_space,
                         const Visitor& conn_vis, PositionMap position,
                         DistanceMap distance, PredecessorMap predecessor,
                         WeightMap weight) const {

    // need to update all the children of the v node:
    std::stack<Vertex> incons;
    incons.push(v);
    while (!incons.empty()) {
      Vertex s = incons.top();
      incons.pop();
      for (auto e : out_edges(s, g)) {
        Vertex t = target(e, g);
        if (t == s) {
          t = source(e, g);
        }
        if (s != get(predecessor, get_property(g, t))) {
          continue;
        }
        put(distance, get_property(g, t),
            get(distance, get_property(g, s)) +
                get(weight, get_property(g, e)));

        conn_vis.affected_vertex(t, g);  // affected by changed distance value.

        put(data->key, t,
            get(distance, get_property(g, t)) +
                get(ReaK::pp::distance_metric, super_space)(
                    get(position, get_property(g, t)),
                    get(position, get_property(g, data->goal_vertex)),
                    super_space));
        data->q.push_or_update(t);

        incons.push(t);
      }
    }

    if (get(predecessor, get_property(g, data->goal_vertex)) !=
        bagl::graph_traits<Graph>::null_vertex()) {
      data->prune_worst_nodes(get(distance, get_property(g, data->goal_vertex)),
                              g, conn_vis);
    }
  }

  template <pp::MetricSpace Space, BNBConnectorVisitor<Graph, Space> Visitor,
            typename PositionMap, typename FwdDistanceMap,
            typename SuccessorMap, typename WeightMap>
  void update_predecessors(Vertex v, Graph& g, const Space& super_space,
                           const Visitor& conn_vis, PositionMap position,
                           FwdDistanceMap fwd_distance, SuccessorMap successor,
                           WeightMap weight) const {

    // need to update all the children of the v node:
    std::stack<Vertex> incons;
    incons.push(v);
    while (!incons.empty()) {
      Vertex t = incons.top();
      incons.pop();
      for (auto e : in_edges(t, g)) {
        Vertex s = source(e, g);
        if (t == s) {
          s = target(e, g);
        }
        if (t != get(successor, get_property(g, s))) {
          continue;
        }
        put(fwd_distance, get_property(g, s),
            get(fwd_distance, get_property(g, t)) +
                get(weight, get_property(g, e)));

        conn_vis.affected_vertex(s, g);  // affected by changed distance value.

        put(data->key, s,
            get(fwd_distance, get_property(g, s)) +
                get(ReaK::pp::distance_metric, super_space)(
                    get(position, get_property(g, data->start_vertex)),
                    get(position, get_property(g, s)), super_space));
        data->q.push_or_update(s);

        incons.push(s);
      }
    }

    if (get(successor, get_property(g, data->start_vertex)) !=
        bagl::graph_traits<Graph>::null_vertex()) {
      data->prune_worst_nodes(
          get(fwd_distance, get_property(g, data->start_vertex)), g, conn_vis);
    }
  }

  /**
   * This call operator takes a position value, the predecessor from which the new position was obtained,
   * the travel-record (as an edge property) that can do the travel from the predecessor to the new position,
   * and the other objects needed for motion planning, and it creates a new vertex for the new position and
   * connects that new vertex to the motion-graph using a lazy and pruned strategy.
   * \note This version applies to a undirected graph (and undirected / symmetric distance metric).
   *
   * \tparam Graph2 The graph type that can store the generated roadmap, should model
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Space The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam Visitor The type of the node-connector visitor to be used.
   * \tparam PositionMap A property-map type that can store the position of each vertex.
   * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting
   *         vertex together with its optimal predecessor.
   * \tparam WeightMap This property-map type is used to store the weights of the edge-properties of the
   *         graph (cost of travel along an edge).
   * \tparam NcSelector A functor type that can select a list of vertices of the graph that are
   *         the nearest-neighbors of a given vertex (or some other heuristic to select the neighbors).
   *         See classes in the topological_search.hpp header-file.
   *
   * \param p The position of the new vertex to be added and connected to the motion-graph.
   * \param x_near The predecessor from which the new vertex was generated (e.g., expanded from, random-walk, etc.).
   * \param eprop The edge-property corresponding to the travel from x_near to p.
   * \param g A mutable graph that should initially store the starting
   *        vertex (if not it will be randomly generated) and will store
   *        the generated graph once the algorithm has finished.
   * \param super_space A topology (as defined by the Born Again Graph Library). This topology
   *        should not include collision checking in its distance metric.
   * \param conn_vis A node-connector visitor implementing the BNBConnectorVisitorConcept. This is the
   *        main point of customization and recording of results that the user can implement.
   * \param position A mapping that implements the MutablePropertyMap Concept. Also,
   *        the value_type of this map should be the same type as the topology's point_type.
   * \param distance The property-map which stores the accumulated distance of each vertex to the root.
   * \param predecessor The property-map which will store the resulting path by connecting
   *        vertices together with their optimal predecessor (follow in reverse to discover the
   *        complete path).
   * \param weight The property-map which stores the weight of each edge-property object (the cost of travel
   *        along the edge).
   * \param select_neighborhood A callable object (functor) that can select a list of
   *        vertices of the graph that ought to be connected to a new
   *        vertex. The list should be sorted in order of increasing "distance".
   */
  template <typename Graph2, pp::MetricSpace Space,
            BNBConnectorVisitor<Graph2, Space> Visitor,
            bagl::concepts::ReadableVPropMemberMap<Graph2> PositionMap,
            bagl::concepts::ReadWriteVPropMemberMap<Graph2> DistanceMap,
            bagl::concepts::ReadWriteVPropMemberMap<Graph2> PredecessorMap,
            bagl::concepts::ReadableEPropMemberMap<Graph2> WeightMap,
            typename NcSelector>
  void operator()(const bagl::property_traits_value_t<PositionMap>& p,
                  Vertex& x_near, EdgeProp& eprop, Graph2& g,
                  const Space& super_space, const Visitor& conn_vis,
                  PositionMap position, DistanceMap distance,
                  PredecessorMap predecessor, WeightMap weight,
                  NcSelector select_neighborhood) const {
    double dist_from_start = get(ReaK::pp::distance_metric, super_space)(
        get(position, get_property(g, data->start_vertex)), p, super_space);
    double dist_to_goal = get(ReaK::pp::distance_metric, super_space)(
        p, get(position, get_property(g, data->goal_vertex)), super_space);

    if ((get(predecessor, get_property(g, data->goal_vertex)) !=
         bagl::graph_traits<Graph>::null_vertex()) &&
        (dist_from_start + dist_to_goal >
         get(distance, get_property(g, data->goal_vertex)))) {
      return;
    }

    using std::back_inserter;
    std::vector<Vertex> Pred;
    std::vector<Vertex> Succ;
    if constexpr (bagl::is_undirected_graph_v<Graph2>) {
      select_neighborhood(
          p, back_inserter(Pred), g, super_space,
          bagl::composite_property_map(position, get(bagl::vertex_all, g)));
    } else {
      select_neighborhood(
          p, back_inserter(Pred), back_inserter(Succ), g, super_space,
          bagl::composite_property_map(position, get(bagl::vertex_all, g)));
    }

    Vertex v = conn_vis.create_vertex(p, g);
    put(data->index_in_heap, v, std::numeric_limits<std::size_t>::max());

    if (x_near != bagl::graph_traits<Graph>::null_vertex()) {
      conn_vis.travel_explored(x_near, v, g);
      conn_vis.travel_succeeded(x_near, v, g);
      conn_vis.affected_vertex(x_near, g);
    }

    lazy_node_connector::connect_best_predecessor(
        v, x_near, eprop, g, super_space, conn_vis, position, distance,
        predecessor, weight, Pred);
    pruned_node_connector::create_pred_edge(v, x_near, std::move(eprop), g,
                                            conn_vis, distance, predecessor,
                                            weight);

    if ((get(predecessor, get_property(g, data->goal_vertex)) !=
         bagl::graph_traits<Graph>::null_vertex()) &&
        (get(distance, get_property(g, v)) + dist_to_goal >
         get(distance, get_property(g, data->goal_vertex)))) {
      conn_vis.vertex_to_be_removed(v, g);
      clear_vertex(v, g);
      remove_vertex(v, g);
      return;
    }

    put(data->key, v, get(distance, get_property(g, v)) + dist_to_goal);
    data->q.push(v);

    if constexpr (bagl::is_undirected_graph_v<Graph2>) {
      Pred.swap(Succ);
    }
    lazy_node_connector::connect_successors(v, x_near, g, super_space, conn_vis,
                                            position, distance, predecessor,
                                            weight, Succ);
    update_successors(v, g, super_space, conn_vis, position, distance,
                      predecessor, weight);
  }

  template <typename Graph2, pp::MetricSpace Space,
            BNBConnectorVisitor<Graph2, Space> Visitor,
            bagl::concepts::ReadableVPropMemberMap<Graph2> PositionMap,
            bagl::concepts::ReadWriteVPropMemberMap<Graph2> DistanceMap,
            bagl::concepts::ReadWriteVPropMemberMap<Graph2> PredecessorMap,
            bagl::concepts::ReadWriteVPropMemberMap<Graph2> FwdDistanceMap,
            bagl::concepts::ReadWriteVPropMemberMap<Graph2> SuccessorMap,
            bagl::concepts::ReadableEPropMemberMap<Graph2> WeightMap,
            typename NcSelector>
  void operator()(const bagl::property_traits_value_t<PositionMap>& p,
                  Vertex& x_pred, EdgeProp& eprop_pred, Vertex& x_succ,
                  EdgeProp& eprop_succ, Graph2& g, const Space& super_space,
                  const Visitor& conn_vis, PositionMap position,
                  DistanceMap distance, PredecessorMap predecessor,
                  FwdDistanceMap fwd_distance, SuccessorMap successor,
                  WeightMap weight, NcSelector select_neighborhood) const {
    using std::back_inserter;

    double dist_from_start = get(ReaK::pp::distance_metric, super_space)(
        get(position, get_property(g, data->start_vertex)), p, super_space);
    double dist_to_goal = get(ReaK::pp::distance_metric, super_space)(
        p, get(position, get_property(g, data->goal_vertex)), super_space);

    if ((get(predecessor, get_property(g, data->goal_vertex)) !=
         bagl::graph_traits<Graph>::null_vertex()) &&
        (dist_from_start + dist_to_goal >
         get(distance, get_property(g, data->goal_vertex)))) {
      return;
    }

    using std::back_inserter;
    std::vector<Vertex> Pred;
    std::vector<Vertex> Succ;
    if constexpr (bagl::is_undirected_graph_v<Graph2>) {
      select_neighborhood(
          p, back_inserter(Pred), g, super_space,
          bagl::composite_property_map(position, get(bagl::vertex_all, g)));
    } else {
      select_neighborhood(
          p, back_inserter(Pred), back_inserter(Succ), g, super_space,
          bagl::composite_property_map(position, get(bagl::vertex_all, g)));
    }

    Vertex v = conn_vis.create_vertex(p, g);
    put(data->index_in_heap, v, std::numeric_limits<std::size_t>::max());

    if (x_pred != bagl::graph_traits<Graph>::null_vertex()) {
      conn_vis.travel_explored(x_pred, v, g);
      conn_vis.travel_succeeded(x_pred, v, g);
      conn_vis.affected_vertex(x_pred, g);
    }
    if (x_succ != bagl::graph_traits<Graph>::null_vertex()) {
      conn_vis.travel_explored(v, x_succ, g);
      conn_vis.travel_succeeded(v, x_succ, g);
      conn_vis.affected_vertex(x_succ, g);
    }

    lazy_node_connector::connect_best_predecessor(
        v, x_pred, eprop_pred, g, super_space, conn_vis, position, distance,
        predecessor, weight, Pred);
    if constexpr (bagl::is_undirected_graph_v<Graph2>) {
      Pred.swap(Succ);
    }
    lazy_node_connector::connect_best_successor(
        v, x_succ, eprop_succ, g, super_space, conn_vis, position, fwd_distance,
        successor, weight, Succ);

    if ((x_pred == bagl::graph_traits<Graph>::null_vertex()) &&
        (x_succ == bagl::graph_traits<Graph>::null_vertex())) {
      conn_vis.vertex_to_be_removed(v, g);
      clear_vertex(v, g);
      remove_vertex(v, g);
      return;
    }

    if (x_pred != bagl::graph_traits<Graph>::null_vertex()) {
      pruned_node_connector::create_pred_edge(v, x_pred, std::move(eprop_pred),
                                              g, conn_vis, distance,
                                              predecessor, weight);
    }
    if (x_succ != bagl::graph_traits<Graph>::null_vertex()) {
      pruned_node_connector::create_succ_edge(v, x_succ, std::move(eprop_succ),
                                              g, conn_vis, fwd_distance,
                                              successor, weight);
    }

    if ((get(predecessor, get_property(g, data->goal_vertex)) !=
         bagl::graph_traits<Graph>::null_vertex()) &&
        (get(distance, get_property(g, v)) + dist_to_goal >
         get(distance, get_property(g, data->goal_vertex)))) {
      conn_vis.vertex_to_be_removed(v, g);
      clear_vertex(v, g);
      remove_vertex(v, g);
      return;
    }
    put(data->key, v, get(distance, get_property(g, v)) + dist_to_goal);
    data->q.push(v);

    lazy_node_connector::connect_successors(v, x_pred, g, super_space, conn_vis,
                                            position, distance, predecessor,
                                            weight, Succ, successor);
    update_successors(v, g, super_space, conn_vis, position, distance,
                      predecessor, weight);
    if constexpr (bagl::is_undirected_graph_v<Graph2>) {
      Pred.swap(Succ);
    }
    lazy_node_connector::connect_predecessors(
        v, x_succ, g, super_space, conn_vis, position, fwd_distance, successor,
        weight, Pred, predecessor);
    update_predecessors(v, g, super_space, conn_vis, position, fwd_distance,
                        successor, weight);
  }
};

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_BRANCH_AND_BOUND_CONNECTOR_H_

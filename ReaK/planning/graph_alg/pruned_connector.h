/**
 * \file pruned_connector.h
 *
 * This library provides a class template and concept that implement a Pruned Motion-graph Connector.
 * A Pruned-Connector uses the accumulated distance to assess the local optimality of the wirings on a motion-graph.
 * This algorithm has many customization points because it can be used in many different sampling-based
 * motion-planners.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2014
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

#ifndef REAK_PLANNING_GRAPH_ALG_PRUNED_CONNECTOR_H_
#define REAK_PLANNING_GRAPH_ALG_PRUNED_CONNECTOR_H_

#include <iterator>
#include <tuple>
#include <type_traits>
#include <utility>

#include "ReaK/topologies/spaces/metric_space_concept.h"

#include "ReaK/planning/graph_alg/sbmp_visitor_concepts.h"

#include "bagl/graph_concepts.h"
#include "bagl/more_property_maps.h"
#include "bagl/properties.h"
#include "bagl/single_property_map.h"

#include <stack>
#include <vector>

namespace ReaK::graph {

template <typename Graph>
struct null_vertex_prop_map : bagl::single_writable_property_map<
                                  bagl::graph_vertex_descriptor_t<Graph>> {
  null_vertex_prop_map()
      : bagl::single_writable_property_map<
            bagl::graph_vertex_descriptor_t<Graph>>(
            bagl::graph_traits<Graph>::null_vertex()) {}
};

struct infinite_double_value_prop_map
    : bagl::single_writable_property_map<double> {
  infinite_double_value_prop_map()
      : bagl::single_writable_property_map<double>(
            std::numeric_limits<double>::infinity()) {}
};

/**
 * This callable class template implements a Pruned Motion-graph Connector.
 * A Pruned-Connector uses the accumulated distance to assess the local optimality of the wirings on a motion-graph.
 * The call operator accepts a visitor object to provide customized behavior because it can be used in many
 * different sampling-based motion-planners. The visitor must model the MotionGraphConnectorVisitorConcept concept.
 */
struct pruned_node_connector {

  template <typename Graph, pp::MetricSpace Space, typename Visitor,
            typename PositionMap, typename DistanceMap, typename PredecessorMap,
            typename WeightMap>
  requires NodeReConnectVisitor<Visitor, Graph>&&
      NeighborhoodTrackingVisitor<Visitor, Graph> static void
      connect_best_predecessor(
          bagl::graph_vertex_descriptor_t<Graph> v,
          bagl::graph_vertex_descriptor_t<Graph>& x_near,
          bagl::edge_property_type<Graph>& eprop, Graph& g,
          const Space& super_space, const Visitor& conn_vis,
          PositionMap position, DistanceMap distance,
          PredecessorMap predecessor, WeightMap weight,
          std::vector<bagl::graph_vertex_descriptor_t<Graph>>& Pred) {
    const auto null_v = bagl::graph_traits<Graph>::null_vertex();

    auto x_near_original = x_near;
    double d_near = std::numeric_limits<double>::infinity();
    if (x_near != null_v) {
      d_near = get(distance, get_property(g, x_near)) + get(weight, eprop);
    }

    for (auto u : Pred) {
      if ((u == x_near_original) ||
          (get(predecessor, get_property(g, u)) == null_v)) {
        continue;
      }

      auto [can_connect, eprop_new] = conn_vis.can_be_connected(u, v, g);
      conn_vis.travel_explored(u, v, g);
      if (can_connect) {
        conn_vis.travel_succeeded(u, v, g);
        double d_out =
            get(weight, eprop_new) + get(distance, get_property(g, u));
        if (d_out < d_near) {
          // edge will be useful as an in-edge to v.
          x_near = u;
          d_near = d_out;
          eprop = std::move(eprop_new);
        }
      } else {
        conn_vis.travel_failed(u, v, g);
      }
      conn_vis.affected_vertex(u, g);  // affected by travel attempts.
    }
    // affected by travel attempts and new in-going edge.
    conn_vis.affected_vertex(v, g);
  }

  template <typename Graph, pp::MetricSpace Space, typename Visitor,
            typename PositionMap, typename FwdDistanceMap,
            typename SuccessorMap, typename WeightMap>
  requires NodeReConnectVisitor<Visitor, Graph>&&
      NeighborhoodTrackingVisitor<Visitor, Graph> static void
      connect_best_successor(
          bagl::graph_vertex_descriptor_t<Graph> v,
          bagl::graph_vertex_descriptor_t<Graph>& x_near,
          bagl::edge_property_type<Graph>& eprop, Graph& g,
          const Space& super_space, const Visitor& conn_vis,
          PositionMap position, FwdDistanceMap fwd_distance,
          SuccessorMap successor, WeightMap weight,
          std::vector<bagl::graph_vertex_descriptor_t<Graph>>& succ) {
    const auto null_v = bagl::graph_traits<Graph>::null_vertex();

    auto x_near_original = x_near;
    double d_near = std::numeric_limits<double>::infinity();
    if (x_near != null_v) {
      d_near = get(fwd_distance, get_property(g, x_near)) + get(weight, eprop);
    }

    for (auto u : succ) {
      if ((u == x_near_original) ||
          (get(successor, get_property(g, u)) == null_v)) {
        continue;
      }

      auto [can_connect, eprop_new] = conn_vis.can_be_connected(v, u, g);
      conn_vis.travel_explored(v, u, g);
      double d_in =
          get(weight, eprop_new) + get(fwd_distance, get_property(g, u));
      if (can_connect) {
        conn_vis.travel_succeeded(v, u, g);
        if (d_in < d_near) {
          // edge will be useful as an in-edge to v.
          x_near = u;
          d_near = d_in;
          eprop = std::move(eprop_new);
        }
      } else {
        conn_vis.travel_failed(v, u, g);
      }
      conn_vis.affected_vertex(u, g);  // affected by travel attempts.
    }
    // affected by travel attempts and new in-going edge.
    conn_vis.affected_vertex(v, g);
  }

  template <typename Graph, pp::MetricSpace Space, typename Visitor,
            typename PositionMap, typename FwdDistanceMap,
            typename SuccessorMap, typename WeightMap,
            typename PredecessorMap = null_vertex_prop_map<Graph>>
  requires NodeReConnectVisitor<Visitor, Graph>&& NeighborhoodTrackingVisitor<
      Visitor, Graph>&& SBMPVisitor<Visitor, Graph> static void
  connect_predecessors(
      bagl::graph_vertex_descriptor_t<Graph> v,
      bagl::graph_vertex_descriptor_t<Graph> x_near, Graph& g,
      const Space& super_space, const Visitor& conn_vis, PositionMap position,
      FwdDistanceMap fwd_distance, SuccessorMap successor, WeightMap weight,
      std::vector<bagl::graph_vertex_descriptor_t<Graph>>& pred,
      PredecessorMap predecessor = PredecessorMap()) {
    const auto null_v = bagl::graph_traits<Graph>::null_vertex();

    for (auto u : pred) {
      if ((u == x_near) || (get(predecessor, get_property(g, u)) != null_v)) {
        continue;
      }

      auto [can_connect, eprop_new] = conn_vis.can_be_connected(u, v, g);
      conn_vis.travel_explored(u, v, g);
      if (can_connect) {
        conn_vis.travel_succeeded(u, v, g);
        double d_in =
            get(weight, eprop_new) + get(fwd_distance, get_property(g, v));
        if (d_in < get(fwd_distance, get_property(g, u))) {
          // edge is useful as an in-edge to (u).
          auto [e_new, e_new_exists] = add_edge(u, v, g, std::move(eprop_new));
          if (e_new_exists) {
            put(fwd_distance, get_property(g, u), d_in);
            auto old_succ = get(successor, get_property(g, u));
            put(successor, get_property(g, u), v);
            conn_vis.edge_added(e_new, g);
            if ((old_succ != u) && (old_succ != null_v)) {
              remove_edge(u, old_succ, g);
            }
          }
        }
      } else {
        conn_vis.travel_failed(v, u, g);
      }
      conn_vis.affected_vertex(u, g);  // affected by travel attempts.
    }
    // affected by travel attempts and new out-going edges.
    conn_vis.affected_vertex(v, g);
  }

  template <typename Graph, pp::MetricSpace Space, typename Visitor,
            typename PositionMap, typename DistanceMap, typename PredecessorMap,
            typename WeightMap,
            typename SuccessorMap = null_vertex_prop_map<Graph>>
  requires NodeReConnectVisitor<Visitor, Graph>&& NeighborhoodTrackingVisitor<
      Visitor, Graph>&& SBMPVisitor<Visitor, Graph> static void
  connect_successors(bagl::graph_vertex_descriptor_t<Graph> v,
                     bagl::graph_vertex_descriptor_t<Graph> x_near, Graph& g,
                     const Space& super_space, const Visitor& conn_vis,
                     PositionMap position, DistanceMap distance,
                     PredecessorMap predecessor, WeightMap weight,
                     std::vector<bagl::graph_vertex_descriptor_t<Graph>>& succ,
                     SuccessorMap successor = SuccessorMap()) {
    const auto null_v = bagl::graph_traits<Graph>::null_vertex();

    for (auto u : succ) {
      if ((u == x_near) || (get(successor, get_property(g, u)) != null_v)) {
        continue;
      }

      auto [can_connect, eprop_new] = conn_vis.can_be_connected(v, u, g);
      conn_vis.travel_explored(v, u, g);
      if (can_connect) {
        conn_vis.travel_succeeded(v, u, g);
        double d_in =
            get(weight, eprop_new) + get(distance, get_property(g, v));
        if (d_in < get(distance, get_property(g, u))) {
          // edge is useful as an in-edge to (u).
          auto [e_new, e_new_exists] = add_edge(v, u, g, std::move(eprop_new));
          if (e_new_exists) {
            put(distance, get_property(g, u), d_in);
            auto old_pred = get(predecessor, get_property(g, u));
            put(predecessor, get_property(g, u), v);
            conn_vis.edge_added(e_new, g);
            if ((old_pred != u) && (old_pred != null_v)) {
              remove_edge(old_pred, u, g);
            }
          }
        }
      } else {
        conn_vis.travel_failed(u, v, g);
      }
      conn_vis.affected_vertex(u, g);  // affected by travel attempts.
    }
    // affected by travel attempts and new out-going edges.
    conn_vis.affected_vertex(v, g);
  }

  template <typename Graph, NeighborhoodTrackingVisitor<Graph> Visitor,
            typename DistanceMap, typename PredecessorMap, typename WeightMap>
  static void update_successors(bagl::graph_vertex_descriptor_t<Graph> v,
                                Graph& g, const Visitor& conn_vis,
                                DistanceMap distance,
                                PredecessorMap predecessor, WeightMap weight) {
    // need to update all the children of the v node:
    std::stack<bagl::graph_vertex_descriptor_t<Graph>> incons;
    incons.push(v);
    while (!incons.empty()) {
      auto s = incons.top();
      incons.pop();
      for (auto e : out_edges(s, g)) {
        auto t = target(e, g);
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
        incons.push(t);
      }
    }
  }

  template <typename Graph, NeighborhoodTrackingVisitor<Graph> Visitor,
            typename FwdDistanceMap, typename SuccessorMap, typename WeightMap>
  static void update_predecessors(bagl::graph_vertex_descriptor_t<Graph> v,
                                  Graph& g, const Visitor& conn_vis,
                                  FwdDistanceMap fwd_distance,
                                  SuccessorMap successor, WeightMap weight) {
    // need to update all the children of the v node:
    std::stack<bagl::graph_vertex_descriptor_t<Graph>> incons;
    incons.push(v);
    while (!incons.empty()) {
      auto t = incons.top();
      incons.pop();
      for (auto e : in_edges(t, g)) {
        auto s = source(e, g);
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
        incons.push(s);
      }
    }
  }

  template <typename Graph, SBMPVisitor<Graph> Visitor, typename DistanceMap,
            typename PredecessorMap, typename WeightMap>
  static void create_pred_edge(bagl::graph_vertex_descriptor_t<Graph> v,
                               bagl::graph_vertex_descriptor_t<Graph>& x_near,
                               bagl::edge_property_type<Graph>&& eprop,
                               Graph& g, const Visitor& conn_vis,
                               DistanceMap distance, PredecessorMap predecessor,
                               WeightMap weight) {
    double d_near = get(weight, eprop) + get(distance, get_property(g, x_near));
    auto [e_new, e_new_exists] = add_edge(x_near, v, g, std::move(eprop));
    if (e_new_exists) {
      put(distance, get_property(g, v), d_near);
      put(predecessor, get_property(g, v), x_near);
      conn_vis.edge_added(e_new, g);
    }
  }

  template <typename Graph, SBMPVisitor<Graph> Visitor, typename FwdDistanceMap,
            typename SuccessorMap, typename WeightMap>
  static void create_succ_edge(bagl::graph_vertex_descriptor_t<Graph> v,
                               bagl::graph_vertex_descriptor_t<Graph>& x_near,
                               bagl::edge_property_type<Graph>&& eprop,
                               Graph& g, const Visitor& conn_vis,
                               FwdDistanceMap fwd_distance,
                               SuccessorMap successor, WeightMap weight) {
    double d_near =
        get(weight, eprop) + get(fwd_distance, get_property(g, x_near));
    auto [e_new, e_new_exists] = add_edge(v, x_near, g, std::move(eprop));
    if (e_new_exists) {
      put(fwd_distance, get_property(g, v), d_near);
      put(successor, get_property(g, v), x_near);
      conn_vis.edge_added(e_new, g);
    }
  }

  /**
   * This call operator takes a position value, the predecessor from which the new position was obtained,
   * the travel-record (as an edge property) that can do the travel from the predecessor to the new position,
   * and the other objects needed for motion planning, and it creates a new vertex for the new position and
   * connects that new vertex to the motion-graph using a lazy and pruned strategy.
   * \note This version applies to a undirected graph (and undirected / symmetric distance metric).
   *
   * \tparam Graph The graph type that can store the generated roadmap, should model
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Space The topology type that represents the free-space, should model MetricSpace.
   * \tparam Visitor The type of the node-connector visitor to be used, should model the
   *MotionGraphConnectorVisitorConcept.
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
   * \param conn_vis A node-connector visitor implementing the MotionGraphConnectorVisitorConcept. This is the
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
  template <typename Graph, pp::MetricSpace Space,
            MotionGraphConnectorVisitor<Graph, Space> Visitor,
            bagl::concepts::ReadableVPropMemberMap<Graph> PositionMap,
            bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
            bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap,
            bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
            typename NcSelector>
  void operator()(const bagl::property_traits_value_t<PositionMap>& p,
                  bagl::graph_vertex_descriptor_t<Graph>& x_near,
                  bagl::edge_property_type<Graph>& eprop, Graph& g,
                  const Space& super_space, const Visitor& conn_vis,
                  PositionMap position, DistanceMap distance,
                  PredecessorMap predecessor, WeightMap weight,
                  NcSelector select_neighborhood) const {
    using Vertex = bagl::graph_vertex_descriptor_t<Graph>;

    std::vector<Vertex> pred;
    std::vector<Vertex> succ;
    if constexpr (bagl::is_undirected_graph_v<Graph>) {
      select_neighborhood(
          p, std::back_inserter(pred), g, super_space,
          bagl::composite_property_map(position, get(bagl::vertex_all, g)));
    } else {
      select_neighborhood(
          p, std::back_inserter(pred), std::back_inserter(succ), g, super_space,
          bagl::composite_property_map(position, get(bagl::vertex_all, g)));
    }

    Vertex v = conn_vis.create_vertex(p, g);

    if (x_near != bagl::graph_traits<Graph>::null_vertex()) {
      conn_vis.travel_explored(x_near, v, g);
      conn_vis.travel_succeeded(x_near, v, g);
      conn_vis.affected_vertex(x_near, g);
    }

    connect_best_predecessor(v, x_near, eprop, g, super_space, conn_vis,
                             position, distance, predecessor, weight, pred);

    if (x_near == bagl::graph_traits<Graph>::null_vertex()) {
      conn_vis.vertex_to_be_removed(v, g);
      clear_vertex(v, g);
      remove_vertex(v, g);
      return;
    }

    create_pred_edge(v, x_near, std::move(eprop), g, conn_vis, distance,
                     predecessor, weight);
    if constexpr (bagl::is_undirected_graph_v<Graph>) {
      pred.swap(succ);
    }
    connect_successors(v, x_near, g, super_space, conn_vis, position, distance,
                       predecessor, weight, succ);
    update_successors(v, g, conn_vis, distance, predecessor, weight);
  }

  template <typename Graph, pp::MetricSpace Space,
            MotionGraphConnectorVisitor<Graph, Space> Visitor,
            bagl::concepts::ReadableVPropMemberMap<Graph> PositionMap,
            bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
            bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap,
            bagl::concepts::ReadWriteVPropMemberMap<Graph> FwdDistanceMap,
            bagl::concepts::ReadWriteVPropMemberMap<Graph> SuccessorMap,
            bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
            typename NcSelector>
  void operator()(const bagl::property_traits_value_t<PositionMap>& p,
                  bagl::graph_vertex_descriptor_t<Graph>& x_pred,
                  bagl::edge_property_type<Graph>& eprop_pred,
                  bagl::graph_vertex_descriptor_t<Graph>& x_succ,
                  bagl::edge_property_type<Graph>& eprop_succ, Graph& g,
                  const Space& super_space, const Visitor& conn_vis,
                  PositionMap position, DistanceMap distance,
                  PredecessorMap predecessor, FwdDistanceMap fwd_distance,
                  SuccessorMap successor, WeightMap weight,
                  NcSelector select_neighborhood) const {
    using Vertex = bagl::graph_vertex_descriptor_t<Graph>;
    const Vertex null_v = bagl::graph_traits<Graph>::null_vertex();

    std::vector<Vertex> pred;
    std::vector<Vertex> succ;
    if constexpr (bagl::is_undirected_graph_v<Graph>) {
      select_neighborhood(
          p, std::back_inserter(pred), g, super_space,
          bagl::composite_property_map(position, get(bagl::vertex_all, g)));
    } else {
      select_neighborhood(
          p, std::back_inserter(pred), std::back_inserter(succ), g, super_space,
          bagl::composite_property_map(position, get(bagl::vertex_all, g)));
    }

    Vertex v = conn_vis.create_vertex(p, g);

    if (x_pred != null_v) {
      conn_vis.travel_explored(x_pred, v, g);
      conn_vis.travel_succeeded(x_pred, v, g);
      conn_vis.affected_vertex(x_pred, g);
    }
    if (x_succ != null_v) {
      conn_vis.travel_explored(v, x_succ, g);
      conn_vis.travel_succeeded(v, x_succ, g);
      conn_vis.affected_vertex(x_succ, g);
    }

    connect_best_predecessor(v, x_pred, eprop_pred, g, super_space, conn_vis,
                             position, distance, predecessor, weight, pred);
    if constexpr (bagl::is_undirected_graph_v<Graph>) {
      pred.swap(succ);
    }
    connect_best_successor(v, x_succ, eprop_succ, g, super_space, conn_vis,
                           position, fwd_distance, successor, weight, succ);

    if ((x_pred == null_v) && (x_succ == null_v)) {
      conn_vis.vertex_to_be_removed(v, g);
      clear_vertex(v, g);
      remove_vertex(v, g);
      return;
    }

    if (x_pred != null_v) {
      create_pred_edge(v, x_pred, std::move(eprop_pred), g, conn_vis, distance,
                       predecessor, weight);
    }
    if (x_succ != null_v) {
      create_succ_edge(v, x_succ, std::move(eprop_succ), g, conn_vis,
                       fwd_distance, successor, weight);
    }

    connect_successors(v, x_pred, g, super_space, conn_vis, position, distance,
                       predecessor, weight, succ, successor);
    update_successors(v, g, conn_vis, distance, predecessor, weight);
    if constexpr (bagl::is_undirected_graph_v<Graph>) {
      pred.swap(succ);
    }
    connect_predecessors(v, x_succ, g, super_space, conn_vis, position,
                         fwd_distance, successor, weight, pred, predecessor);
    update_predecessors(v, g, conn_vis, fwd_distance, successor, weight);
  }
};

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_PRUNED_CONNECTOR_H_

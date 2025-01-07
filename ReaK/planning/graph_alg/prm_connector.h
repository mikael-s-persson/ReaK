/**
 * \file prm_connector.h
 *
 * This library provides a class template and concept that implement a PRM Motion-graph Connector.
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

#ifndef REAK_PLANNING_GRAPH_ALG_PRM_CONNECTOR_H_
#define REAK_PLANNING_GRAPH_ALG_PRM_CONNECTOR_H_

#include <functional>
#include <tuple>
#include <type_traits>

#include "ReaK/topologies/spaces/metric_space_concept.h"

#include "bagl/graph_concepts.h"
#include "bagl/more_property_maps.h"
#include "bagl/properties.h"

#include "ReaK/planning/graph_alg/motion_graph_connector.h"

#include <stack>

namespace ReaK::graph {

/**
 * This callable class template implements a Motion-graph Connector.
 * A connector uses the accumulated distance to assess the local optimality of the wirings on a motion-graph.
 * The call operator accepts a visitor object to provide customized behavior because it can be used in many
 * different sampling-based motion-planners. The visitor must model the MotionGraphConnectorVisitor concept.
 */
struct prm_node_connector {

  template <typename Graph, NeighborhoodTrackingVisitor<Graph> Visitor>
  requires SBMPVisitor<Visitor, Graph> void connect_to_first_pred(
      bagl::graph_vertex_descriptor_t<Graph> v,
      bagl::graph_vertex_descriptor_t<Graph> x_near,
      bagl::edge_property_type<Graph>& eprop, Graph& g, Visitor& conn_vis) {
    conn_vis.travel_explored(x_near, v, g);
    conn_vis.travel_succeeded(x_near, v, g);
    auto [e_new, e_new_exists] = add_edge(x_near, v, g, std::move(eprop));
    if (e_new_exists) {
      conn_vis.edge_added(e_new, g);
      conn_vis.affected_vertex(x_near, g);
    }
  }

  /**
   * This call operator takes a position value, the predecessor from which the new position was obtained,
   * the travel-record (as an edge property) that can do the travel from the predecessor to the new position,
   * and the other objects needed for motion planning, and it creates a new vertex for the new position and
   * connects that new vertex to the motion-graph while keeping track of the connected-component.
   * \note This version applies to a directed graph (and directed / asymmetric distance metric).
   *
   * \tparam Graph The graph type that can store the generated roadmap, should model
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Space The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam Visitor The type of the node-connector visitor to be used.
   * \tparam PositionMap A property-map type that can store the position of each vertex.
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
   * \param conn_vis A node-connector visitor implementing the MotionGraphConnectorVisitor. This is the
   *        main point of customization and recording of results that the user can implement.
   * \param position A mapping that implements the MutablePropertyMap Concept. Also,
   *        the value_type of this map should be the same type as the topology's point_type.
   * \param select_neighborhood A callable object (functor) that can select a list of
   *        vertices of the graph that ought to be connected to a new
   *        vertex. The list should be sorted in order of increasing "distance".
   */
  template <typename Graph, pp::Topology Space,
            MotionGraphConnectorVisitor<Graph, Space> Visitor,
            bagl::concepts::ReadableVPropMemberMap<Graph> PositionMap,
            typename NcSelector>
  void operator()(const bagl::property_traits_value_t<PositionMap>& p,
                  bagl::graph_vertex_descriptor_t<Graph> x_near,
                  bagl::edge_property_type<Graph>& eprop, Graph& g,
                  const Space& super_space, Visitor& conn_vis,
                  PositionMap position, const NcSelector& select_neighborhood) {
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
      connect_to_first_pred(v, x_near, eprop, g, conn_vis);
    }
    conn_vis.affected_vertex(v, g);

    for (Vertex u : pred) {
      if (u == x_near) {
        continue;
      }

      auto [can_connect, eprop_new] = conn_vis.can_be_connected(u, v, g);
      conn_vis.travel_explored(u, v, g);
      if (can_connect) {
        conn_vis.travel_succeeded(u, v, g);
        auto [e_new, e_new_exists] = add_edge(u, v, g, std::move(eprop_new));
        if (e_new_exists) {
          conn_vis.edge_added(e_new, g);
        }
      } else {
        conn_vis.travel_failed(u, v, g);
      }
      conn_vis.affected_vertex(u, g);  // affected by travel attempts.
    }

    if constexpr (!bagl::is_undirected_graph_v<Graph>) {
      for (Vertex u : succ) {
        auto [can_connect, eprop_new] = conn_vis.can_be_connected(v, u, g);
        conn_vis.travel_explored(v, u, g);
        if (can_connect) {
          conn_vis.travel_succeeded(v, u, g);
          auto [e_new, e_new_exists] = add_edge(v, u, g, std::move(eprop_new));
          if (e_new_exists) {
            conn_vis.edge_added(e_new, g);
          }
        } else {
          conn_vis.travel_failed(v, u, g);
        }
        conn_vis.affected_vertex(u, g);  // affected by travel attempts.
      }
    }

    conn_vis.affected_vertex(v, g);
  }
};

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_PRM_CONNECTOR_H_

/**
 * \file motion_graph_connector.hpp
 *
 * This library provides a class template and concept that implement a Motion-graph Connector.
 * A connector uses the accumulated distance to assess the local optimality of the wirings on a motion-graph.
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

#ifndef REAK_MOTION_GRAPH_CONNECTOR_HPP
#define REAK_MOTION_GRAPH_CONNECTOR_HPP

#include <iterator>
#include <tuple>
#include <type_traits>
#include <utility>

#include <ReaK/topologies/spaces/metric_space_concept.hpp>

#include "sbmp_visitor_concepts.hpp"
#include "simple_graph_traits.hpp"

#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>

// BGL-Extra includes:
#include <boost/graph/more_property_maps.hpp>
#include <boost/graph/more_property_tags.hpp>

#include <stack>
#include <vector>

namespace ReaK::graph {

/**
 * This callable class template implements a Motion-graph Connector.
 * A connector uses the accumulated distance to assess the local optimality of the wirings on a motion-graph.
 * The call operator accepts a visitor object to provide customized behavior because it can be used in many
 * different sampling-based motion-planners. The visitor must model the MotionGraphConnectorVisitorConcept concept.
 */
struct motion_graph_connector {

  /**
   * This call operator takes a position value, the predecessor from which the new position was obtained,
   * the travel-record (as an edge property) that can do the travel from the predecessor to the new position,
   * and the other objects needed for motion planning, and it creates a new vertex for the new position and
   * connects that new vertex to the motion-graph using a lazy and pruned strategy.
   * \note This version applies to a undirected graph (and undirected / symmetric distance metric).
   *
   * \tparam Graph The graph type that can store the generated roadmap, should model
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam SBAStarVisitor The type of the node-connector visitor to be used, should model the
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
   * \param super_space A topology (as defined by the Boost Graph Library). This topology
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
  template <typename Graph, typename Topology, typename ConnectorVisitor,
            typename PositionMap, typename DistanceMap, typename PredecessorMap,
            typename WeightMap, typename NcSelector>
  void operator()(const property_value_t<PositionMap>& p,
                  graph_vertex_t<Graph>& x_near,
                  graph_edge_bundle_t<Graph>& eprop, Graph& g,
                  const Topology& super_space, const ConnectorVisitor& conn_vis,
                  PositionMap position, DistanceMap distance,
                  PredecessorMap predecessor, WeightMap weight,
                  NcSelector select_neighborhood) const {
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((
        MotionGraphConnectorVisitorConcept<ConnectorVisitor, Graph, Topology>));

    using Vertex = graph_vertex_t<Graph>;
    using std::back_inserter;

    std::vector<Vertex> Pred;
    std::vector<Vertex> Succ;
    if constexpr (boost::is_undirected_graph<Graph>::value) {
      select_neighborhood(p, back_inserter(Pred), g, super_space,
                          boost::bundle_prop_to_vertex_prop(position, g));
    } else {
      select_neighborhood(p, back_inserter(Pred), back_inserter(Succ), g,
                          super_space,
                          boost::bundle_prop_to_vertex_prop(position, g));
    }

    Vertex v = conn_vis.create_vertex(p, g);

    Vertex x_near_original = x_near;
    conn_vis.travel_explored(x_near, v, g);
    conn_vis.travel_succeeded(x_near, v, g);
    double d_near = get(distance, g[x_near]) + get(weight, eprop);
    {
      auto [e_new, e_new_exists] = add_edge(x_near, v, std::move(eprop), g);
      if (e_new_exists) {
        put(distance, g[v], d_near);
        put(predecessor, g[v], x_near);
        conn_vis.edge_added(e_new, g);
        conn_vis.affected_vertex(x_near, g);
      }
    }
    // affected by travel attempts and new in-going edge.
    conn_vis.affected_vertex(v, g);

    for (Vertex u : Pred) {
      if (u == x_near_original) {
        continue;
      }

      auto [can_connect, eprop2] = conn_vis.can_be_connected(v, u, g);
      conn_vis.travel_explored(v, u, g);
      if (can_connect) {
        conn_vis.travel_succeeded(v, u, g);
        double d_in = get(distance, g[u]) + get(weight, eprop2);
        auto [e_new, e_new_exists] = add_edge(v, u, std::move(eprop2), g);
        if (e_new_exists) {
          if (d_in < d_near) {
            d_near = d_in;
            x_near = u;
            put(distance, g[v], d_in);
            put(predecessor, g[v], u);
          }
          conn_vis.edge_added(e_new, g);
        }
      } else {
        conn_vis.travel_failed(v, u, g);
      }
      conn_vis.affected_vertex(u, g);  // affected by travel attempts.
    }
    // affected by travel attempts and new edges.
    conn_vis.affected_vertex(v, g);

    if constexpr (boost::is_directed_graph<Graph>::value) {
      for (Vertex u : Succ) {
        auto [can_connect, eprop2] = conn_vis.can_be_connected(v, u, g);
        conn_vis.travel_explored(v, u, g);
        if (can_connect) {
          conn_vis.travel_succeeded(v, u, g);
          auto [e_new, e_new_exists] = add_edge(v, u, std::move(eprop2), g);
          if (e_new_exists) {
            conn_vis.edge_added(e_new, g);
          }
        } else {
          conn_vis.travel_failed(v, u, g);
        }
        conn_vis.affected_vertex(u, g);  // affected by travel attempts.
      }
    }

    // need to update all the children of the v node:
    std::stack<Vertex> incons;
    incons.push(v);
    while (!incons.empty()) {
      Vertex s = incons.top();
      incons.pop();
      for (auto [eo, eo_end] = out_edges(s, g); eo != eo_end; ++eo) {
        Vertex t = target(*eo, g);
        if (t == s) {
          t = source(*eo, g);
        }
        if (s != get(predecessor, g[t])) {
          continue;
        }
        put(distance, g[t], get(distance, g[s]) + get(weight, g[*eo]));
        conn_vis.affected_vertex(t, g);  // affected by changed distance value.
        incons.push(t);
      }
    }
  }
};

}  // namespace ReaK::graph

#endif

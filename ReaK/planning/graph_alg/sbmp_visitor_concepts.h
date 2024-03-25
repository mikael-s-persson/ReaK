/**
 * \file sbmp_visitor_concepts.h
 *
 * This library contains the
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

#ifndef REAK_PLANNING_GRAPH_ALG_SBMP_VISITOR_CONCEPTS_H_
#define REAK_PLANNING_GRAPH_ALG_SBMP_VISITOR_CONCEPTS_H_

#include "ReaK/planning/graph_alg/simple_graph_traits.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"

#include <concepts>

namespace ReaK::graph {

/**
  * This concept class defines what is required of a class to serve as a visitor to
  * a Sampling-based Motion-planning algorithm.
  *
  * Required concepts:
  *
  * Valid expressions:
  *
  * vis.vertex_added(u, g);  This function is called whenever a new vertex (u) has been added to the graph (g), but not
  *yet connected.
  *
  * vis.edge_added(e, g);  This function is called whenever a new edge (e) has been created between the last created
  *vertex and its neighbor in the graph (g).
  *
  * bool b = vis.keep_going();  This function is called at each attempt to expand the graph to verify that the user
  *still wants more vertices to be generated in the motion-graph.
  */
template <typename Visitor, typename Graph>
concept SBMPVisitor = requires(Visitor vis, Graph g, graph_vertex_t<Graph> u,
                               graph_edge_t<Graph> e) {
  vis.vertex_added(u, g);
  vis.edge_added(e, g);
  { vis.keep_going() } -> std::convertible_to<bool>;
};

/**
  * This class is simply an archetype visitor for SBMPVisitor.
  */
struct sbmp_visitor_archetype {
  template <typename Vertex, typename Graph>
  void vertex_added(Vertex /*unused*/, Graph& /*unused*/) const {}
  template <typename Edge, typename Graph>
  void edge_added(Edge /*unused*/, Graph& /*unused*/) const {}
  bool keep_going() const { return true; }
};

/**
  * This concept class defines what is required of a class to serve as a visitor to
  * a Sampling-based Motion-planning algorithm that involves pruning of the motion-graph.
  *
  * Required concepts:
  *
  * the visitor should model SBMPVisitor.
  *
  * Valid expressions:
  *
  * vis.vertex_to_be_removed(u, g);  This function is called just before a vertex is removed from the motion-graph
  *(i.e., pruned).
  */
template <typename Visitor, typename Graph>
concept SBMPPruningVisitor = SBMPVisitor<Visitor, Graph>&& requires(
    Visitor vis, Graph g, graph_vertex_t<Graph> u) {
  vis.vertex_to_be_removed(u, g);
};

/**
  * This class is simply an archetype visitor for SBMPPruningVisitor.
  */
struct sbmp_pruning_visitor_archetype : sbmp_visitor_archetype {
  template <typename Vertex, typename Graph>
  void vertex_to_be_removed(Vertex /*unused*/, Graph& /*unused*/) const {}
};

/**
  * This concept class defines what is required of a class to serve as a visitor to
  * a Sampling-based Motion-planning algorithm that involves finding joining-points of
  * multiple motion-graphs (e.g., bi-directional algorithms).
  *
  * Required concepts:
  *
  * the visitor should model SBMPVisitor.
  *
  * Valid expressions:
  *
  * vis.joining_vertex_found(u1, u2, g1, g2);  This function is called by the algorithm when two vertices (u1,u2) of two
  *graphs (g1,g2) is found that meet each other.
  */
template <typename Visitor, typename Graph>
concept SBMPJoiningVisitor = SBMPVisitor<Visitor, Graph>&& requires(
    Visitor vis, Graph g, graph_vertex_t<Graph> u) {
  vis.joining_vertex_found(u, u, g, g);
};

/**
  * This class is simply an archetype visitor for SBMPJoiningVisitor.
  */
struct sbmp_joining_visitor_archetype : sbmp_visitor_archetype {
  template <typename Vertex, typename Graph>
  void joining_vertex_found(Vertex /*unused*/, Vertex /*unused*/,
                            Graph& /*unused*/, Graph& /*unused*/) const {}
};

/**
  * This concept class defines what is required of a class to serve as a node-creator for a motion-planning
  * algorithm. This is typical used internally (in intermediate visitors) for most motion-planners.
  *
  * Required concepts:
  *
  * Valid expressions:
  *
  * u = conn_vis.create_vertex(p, g);  This function is called to request the visitor to create a new vertex (u) in the
  *graph (g), with position (p).
  */
template <typename Visitor, typename Graph, typename Space>
concept NodeCreatorVisitor = pp::Topology<Space>&& requires(
    Visitor vis, Graph g, pp::topology_point_type_t<Space> p) {
  { vis.create_vertex(p, g) } -> std::convertible_to<graph_vertex_t<Graph>>;
};

/**
  * This class is simply an archetype visitor for NodeCreatorVisitor.
  */
struct node_creator_visitor_archetype {
  template <typename Position, typename Graph>
  auto create_vertex(const Position& /*unused*/, Graph& /*unused*/) const {
    return boost::graph_traits<Graph>::null_vertex();
  }
};

/**
  * This concept class defines what is required of a class to serve as a visitor to a motion-planning
  * algorithm that pulls nodes towards some sample-point (i.e., attempt to steer, at least, partially
  * towards a point). This is typical of RRT-style motion-planners.
  *
  * Required concepts:
  *
  * Valid expressions:
  *
  * tie(p,b,ep) = vis.steer_towards_position(p,u,g);  This function is called to attempt to steer from vertex u to
  *position p, it returns a std::pair with the position that could be reached and a boolean value to indicate whether
  *any significant motion occurred (collision-free).
  */
template <typename Visitor, typename Graph, typename Space>
concept NodePullingVisitor = pp::Topology<Space>&& requires(
    Visitor vis, Graph g, graph_vertex_t<Graph> u,
    pp::topology_point_type_t<Space> p, bool b, graph_edge_bundle_t<Graph> ep) {
  std::tie(p, b, ep) = vis.steer_towards_position(p, u, g);
};

/**
  * This class is simply an archetype visitor for NodePullingVisitor.
  */
struct node_pulling_visitor_archetype {
  template <typename Position, typename Vertex, typename Graph>
  auto steer_towards_position(const Position& p, Vertex u, Graph& g) const {
    return std::tuple(p, false, graph_edge_bundle_t<Graph>());
  }
};

/**
  * This concept class defines what is required of a class to serve as a visitor to a motion-planning
  * algorithm that pulls nodes backwards towards some sample-point (i.e., attempt to steer (reverse-time),
  * at least, partially towards a point). This is typical of RRT-style motion-planners.
  *
  * Required concepts:
  *
  * Valid expressions:
  *
  * tie(p,b,ep) = vis.steer_back_to_position(p,u,g);  This function is called to attempt to steer backwards in time from
  *vertex u to position p, it returns a std::pair with the position that could be reached and a boolean value to
  *indicate whether any significant motion occurred (collision-free).
  */
template <typename Visitor, typename Graph, typename Space>
concept NodeBackPullingVisitor = pp::Topology<Space>&& requires(
    Visitor vis, Graph g, graph_vertex_t<Graph> u,
    pp::topology_point_type_t<Space> p, bool b, graph_edge_bundle_t<Graph> ep) {
  std::tie(p, b, ep) = vis.steer_back_to_position(p, u, g);
};

/**
  * This class is simply an archetype visitor for NodeBackPullingVisitorConcept.
  */
struct node_back_pulling_visitor_archetype {
  template <typename Position, typename Vertex, typename Graph>
  auto steer_back_to_position(const Position& p, Vertex u, Graph& g) const {
    return std::tuple(p, false, graph_edge_bundle_t<Graph>());
  }
};

/**
  * This concept class defines what is required of a class to serve as a visitor to a motion-planning
  * algorithm that expands existing nodes towards some unspecified direction (i.e., attempt to steer,
  * at least, partially towards some random or guided direction). This is typical of PRM-style
  * motion-planners and of what one could call "pushed" expansion algorithms (SBA*, FADPRM, etc.).
  *
  * Required concepts:
  *
  * Valid expressions:
  *
  * tie(p,b,ep) = vis.random_walk(u, g);  This function is called to perform the expansion of the roadmap from a given
  *vertex (u) in the graph (g). This function returns a newly generated position value that is a candidate to be added
  *to the graph, and with an edge-property object associated with a new edge.
  */
template <typename Visitor, typename Graph, typename Space>
concept NodePushingVisitor = pp::Topology<Space>&& requires(
    Visitor vis, Graph g, graph_vertex_t<Graph> u,
    pp::topology_point_type_t<Space> p, bool b, graph_edge_bundle_t<Graph> ep) {
  std::tie(p, b, ep) = vis.random_walk(u, g);
};

/**
  * This class is simply an archetype visitor for NodePushingVisitor.
  */
template <pp::Topology Space>
struct node_pushing_visitor_archetype {
  using PointType = pp::topology_point_type_t<Space>;
  template <typename Vertex, typename Graph>
  auto random_walk(Vertex u, Graph& g) const {
    return std::tuple(PointType(), false, graph_edge_bundle_t<Graph>());
  }
};

/**
  * This concept class defines what is required of a class to serve as a visitor to a motion-planning
  * algorithm that retracts existing nodes backwards towards some unspecified direction (i.e., attempt to steer,
  * at least, partially, backwards towards some random or guided direction). This is typical of PRM-style
  * motion-planners and of what one could call "pushed" expansion algorithms (SBA*, FADPRM, etc.).
  *
  * Required concepts:
  *
  * Valid expressions:
  *
  * tie(p,b,ep) = vis.random_back_walk(u, g);  This function is called to perform the retraction of the roadmap from a
  *given vertex (u) in the graph (g). This function returns a newly generated position value that is a candidate to be
  *added to the graph, and with an edge-property object associated with a new edge.
  */
template <typename Visitor, typename Graph, typename Space>
concept NodeBackPushingVisitor = pp::Topology<Space>&& requires(
    Visitor vis, Graph g, graph_vertex_t<Graph> u,
    pp::topology_point_type_t<Space> p, bool b, graph_edge_bundle_t<Graph> ep) {
  std::tie(p, b, ep) = vis.random_back_walk(u, g);
};

/**
  * This class is simply an archetype visitor for NodeBackPushingVisitor.
  */
template <pp::Topology Space>
struct node_back_pushing_visitor_archetype {
  using PointType = pp::topology_point_type_t<Space>;
  template <typename Vertex, typename Graph>
  auto random_back_walk(Vertex u, Graph& g) const {
    return std::tuple(PointType(), false, graph_edge_bundle_t<Graph>());
  }
};

/**
  * This concept class defines what is required of a class to serve as a visitor to a motion-planning
  * algorithm that attempts re-connections of existing nodes (i.e., attempt to steer,
  * fully, from one node to another). This is typical of PRM-style motion-planners and of
  * probabilistically optimal connection strategies which must involve re-wirings (e.g., RRT*, SBA*, etc.).
  *
  * Required concepts:
  *
  * Valid expressions:
  *
  * tie(b, ep) = vis.can_be_connected(u,v,g);  This function is called to attempt to steer from vertex u to vertex v, it
  *returns true if a local path exists and is collision-free, and it also returns the edge-property of the edge that
  *could connect those two vertices.
  */
template <typename Visitor, typename Graph>
concept NodeReConnectVisitor = requires(Visitor vis, Graph g,
                                        graph_vertex_t<Graph> u, bool b,
                                        graph_edge_bundle_t<Graph> ep) {
  std::tie(b, ep) = vis.can_be_connected(u, u, g);
};

/**
  * This class is simply an archetype visitor for NodeReConnectVisitor.
  */
struct node_reconnect_visitor_archetype {
  template <typename Vertex, typename Graph>
  auto can_be_connected(Vertex u, Vertex /*unused*/, Graph& g) const {
    return std::pair(false, graph_edge_bundle_t<Graph>());
  }
};

/**
  * This concept class defines what is required of a class to serve as a visitor to a motion-planning
  * algorithm that requires the tracking of the neighborhood of the nodes (e.g., keeping tabs of the
  * density of the neighborhood). This is typical of PRM-style motion-planners and of
  * other "pushed" expansion algorithms (e.g., SBA*, FADPRM, etc.).
  *
  * Required concepts:
  *
  * Valid expressions:
  *
  * vis.travel_explored(u, v, g);  This function is called whenever a source-destination pair of vertices have been
  *matched for a potential edge between them.
  *
  * vis.travel_succeeded(u, v, g);  This function is called whenever a source-destination pair of vertices led to a
  *successful travel.
  *
  * vis.travel_failed(u, v, g);  This function is called whenever a source-destination pair of vertices led to a failure
  *to travel (e.g., a collision occurred).
  *
  * vis.affected_vertex(u,g);  This function is called to notify the visitor that something about the vertex's
  *neighborhood might have changed (possibly nothing). For example, this would be the place where to re-compute a
  *density metric.
  */
template <typename Visitor, typename Graph>
concept NeighborhoodTrackingVisitor = requires(Visitor vis, Graph g,
                                               graph_vertex_t<Graph> u) {
  vis.travel_explored(u, u, g);
  vis.travel_succeeded(u, u, g);
  vis.travel_failed(u, u, g);
  vis.affected_vertex(u, g);
};

/**
  * This class is simply an archetype visitor for NeighborhoodTrackingVisitor.
  */
struct neighborhood_tracking_visitor_archetype {
  template <typename Vertex, typename Graph>
  void travel_explored(Vertex /*unused*/, Vertex /*unused*/,
                       Graph& /*unused*/) const {}
  template <typename Vertex, typename Graph>
  void travel_succeeded(Vertex /*unused*/, Vertex /*unused*/,
                        Graph& /*unused*/) const {}
  template <typename Vertex, typename Graph>
  void travel_failed(Vertex /*unused*/, Vertex /*unused*/,
                     Graph& /*unused*/) const {}
  template <typename Vertex, typename Graph>
  void affected_vertex(Vertex /*unused*/, Graph& /*unused*/) const {}
};

/**
  * This concept class defines what is required of a class to serve as a visitor to check collision of
  * a point during a motion-planning algorithm.
  *
  * Valid expressions:
  *
  * bool b = vis.is_position_free(p);  This function is called to query whether a particular configuration (position, p)
  *is free.
  *
  * \tparam Visitor The visitor class to be checked for modeling this concept.
  * \tparam Topology The topology that provides the positions with which the visitor class is required to work.
  */
template <typename Visitor, typename Space>
concept CollisionCheckingVisitor = pp::Topology<Space>&& requires(
    Visitor vis, pp::topology_point_type_t<Space> p) {
  { vis.is_position_free(p) } -> std::convertible_to<bool>;
};

/**
  * This class is simply an archetype visitor for CollisionCheckingVisitor.
  */
struct collision_checking_visitor_archetype {
  template <typename Position>
  bool is_position_free(const Position& /*unused*/) const {
    return true;
  }
};

/**
  * This concept class defines the valid expressions required of a class to be used as a
  * node-exploring visitor. This is used in algorithms that look through a set of vertices
  * of the motion-graph (either to expand from them or simply doing a shortest-path algorithm).
  *
  * Required concepts:
  *
  * Valid expressions:
  *
  * vis.initialize_vertex(u, g);  A function that gets called whenever a vertex (u) is first initialized before the
  *search.
  *
  * vis.discover_vertex(u, g);  A function that gets called whenever a vertex (u) is added to the OPEN set (or updated
  *in the OPEN set).
  *
  * vis.examine_vertex(u, g);  A function that gets called whenever a vertex (u) is taken out of the OPEN set to be
  *examined, this is called before it gets expanded.
  *
  * vis.examine_edge(e, g);  A function that is called whenever an edge (e) is examined (usually after its source vertex
  *has been examined).
  *
  * b = vis.has_search_potential(u, g);  This function is called to check if the vertex is deemed minimally useful to
  *examine by some measure.
  *
  * b = vis.should_close(u, g);  This function is called to check if a vertex should be closed according to some
  *measure, i.e., it should no longer be elligible to be examined in the future.
  */
template <typename Visitor, typename Graph>
concept NodeExploringVisitor = requires(Visitor vis, Graph g,
                                        graph_vertex_t<Graph> u,
                                        graph_edge_t<Graph> e) {
  vis.initialize_vertex(u, g);
  vis.discover_vertex(u, g);
  vis.examine_vertex(u, g);
  vis.examine_edge(e, g);
  { vis.has_search_potential(u, g) } -> std::convertible_to<bool>;
  { vis.should_close(u, g) } -> std::convertible_to<bool>;
};

/**
  * This class is simply an archetype visitor for NodeExploringVisitor.
  */
struct node_exploring_visitor_archetype {

  template <typename Vertex, typename Graph>
  void initialize_vertex(Vertex /*unused*/, const Graph& /*unused*/) const {}

  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex /*unused*/, const Graph& /*unused*/) const {}

  template <typename Vertex, typename Graph>
  void examine_vertex(Vertex /*unused*/, const Graph& /*unused*/) const {}

  template <typename Edge, typename Graph>
  void examine_edge(Edge /*unused*/, const Graph& /*unused*/) const {}

  template <typename Vertex, typename Graph>
  bool has_search_potential(Vertex /*unused*/, const Graph& /*unused*/) const {
    return true;
  }

  template <typename Vertex, typename Graph>
  bool should_close(Vertex /*unused*/, const Graph& /*unused*/) const {
    return false;
  }
};

/**
  * This concept class defines what is required of a class to serve as a visitor to the RRT algorithm.
  *
  * Required concepts:
  *
  * the visitor should model SBMPVisitor, NodePullingVisitor, and CollisionCheckingVisitor.
  */
template <typename Visitor, typename Graph, typename Space>
concept RRTVisitor =
    SBMPVisitor<Visitor, Graph>&& pp::Topology<Space>&& NodePullingVisitor<
        Visitor, Graph, Space>&& CollisionCheckingVisitor<Visitor, Space>;

/**
  * This class is simply an archetype visitor for the RRT algorithm.
  */
struct rrt_visitor_archetype : sbmp_visitor_archetype,
                               node_pulling_visitor_archetype,
                               collision_checking_visitor_archetype {};

/**
  * This concept class defines what is required of a class to serve as a visitor to the Bi-RRT algorithm.
  *
  * Required concepts:
  *
  * the visitor should model RRTVisitor and SBMPJoiningVisitor.
  */
template <typename Visitor, typename Graph, typename Space>
concept BiRRTVisitor =
    RRTVisitor<Visitor, Graph, Space>&& SBMPJoiningVisitor<Visitor, Graph>;

/**
  * This class is simply an archetype visitor for the Bi-RRT algorithm.
  */
struct birrt_visitor_archetype : rrt_visitor_archetype,
                                 sbmp_joining_visitor_archetype {};

/**
  * This concept class defines what is required of a class to serve as a visitor to the PRM algorithm.
  *
  * Required concepts:
  *
  * the visitor should model SBMPVisitor, NodePushingVisitor, NodeReConnectVisitor,
  * NeighborhoodTrackingVisitor, and CollisionCheckingVisitor.
  */
template <typename Visitor, typename Graph, typename Space>
concept PRMVisitor =
    SBMPVisitor<Visitor, Graph>&& NodePushingVisitor<Visitor, Graph, Space>&&
        NodeReConnectVisitor<Visitor, Graph>&& NeighborhoodTrackingVisitor<
            Visitor, Graph>&& CollisionCheckingVisitor<Visitor, Space>;

/**
  * This class is simply an archetype visitor for the PRM algorithm.
  */
template <typename Space>
struct prm_visitor_archetype : sbmp_visitor_archetype,
                               node_pushing_visitor_archetype<Space>,
                               node_reconnect_visitor_archetype,
                               neighborhood_tracking_visitor_archetype,
                               collision_checking_visitor_archetype {};

/**
  * This concept class defines what is required of a class to serve as a visitor to the RRG algorithm.
  *
  * Required concepts:
  *
  * the visitor should model RRTVisitor, and NodeReConnectVisitor.
  */
template <typename Visitor, typename Graph, typename Space>
concept RRGVisitor =
    RRTVisitor<Visitor, Graph, Space>&& NodeReConnectVisitor<Visitor, Graph>;

/**
  * This class is simply an archetype visitor for RRGVisitor.
  */
struct rrg_visitor_archetype : rrt_visitor_archetype,
                               node_reconnect_visitor_archetype {};

/**
  * This concept class defines what is required of a class to serve as a visitor to the bi-directional RRG algorithm.
  *
  * Required concepts:
  *
  * the visitor should model RRGVisitorConcept, and NodeBackPullingVisitorConcept.
  */
template <typename Visitor, typename Graph, typename Space>
concept RRGBidirVisitor = RRGVisitor<Visitor, Graph, Space>&&
    NodeBackPullingVisitor<Visitor, Graph, Space>;

/**
  * This class is simply an archetype visitor for RRGBidirVisitor.
  */
struct rrg_bidir_visitor_archetype : rrg_visitor_archetype,
                                     node_back_pulling_visitor_archetype {};

/**
  * This concept class defines the valid expressions required of a class to be used as a visitor
  * class for the connection strategies.
  *
  * Required concepts:
  *
  * The visitor class should model SBMPVisitor, NodeCreatorVisitor, NodeReConnectVisitor,
  * and NeighborhoodTrackingVisitor.
  */
template <typename ConnectorVisitor, typename Graph, typename Space>
concept MotionGraphConnectorVisitor = SBMPVisitor<ConnectorVisitor, Graph>&&
    NodeCreatorVisitor<ConnectorVisitor, Graph, Space>&&
        NodeReConnectVisitor<ConnectorVisitor, Graph>&&
            NeighborhoodTrackingVisitor<ConnectorVisitor, Graph>;

/**
  * This class is simply an archetype visitor for MotionGraphConnectorVisitor.
  */
struct mg_connector_visitor_archetype
    : sbmp_visitor_archetype,
      node_creator_visitor_archetype,
      node_reconnect_visitor_archetype,
      neighborhood_tracking_visitor_archetype {};

/**
  * This concept class defines the valid expressions required of a class to be used as an anytime-heuristic
  * visitor.
  *
  * Valid expressions:
  *
  * d = vis.adjust_relaxation(d, g);  This function should return a new value for the relaxation factor used in the
  *anytime-heuristic algorithm.
  */
template <typename Visitor, typename Graph>
concept AnytimeHeuristicVisitor = requires(Visitor vis, Graph g, double d) {
  { vis.adjust_relaxation(d, g) } -> std::convertible_to<double>;
};

/**
  * This class is simply an archetype visitor for AnytimeHeuristicVisitor.
  */
struct anytime_heuristic_visitor_archetype {
  template <typename Graph>
  double adjust_relaxation(double d, Graph& /*unused*/) const {
    return d;
  }
};

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_SBMP_VISITOR_CONCEPTS_H_

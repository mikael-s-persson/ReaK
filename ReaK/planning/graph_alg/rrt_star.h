/**
 * \file rrt_star.h
 *
 * This library contains the Rapidly-Exploring Random Tree generation algorithm with an
 * A* rewiring policy.
 * This is a method to create a random tree that will span over a non-convex space
 * as rapidly as possible. The method relies on a simple randomized insertion algorithm.
 * At each step, a random point is picked from the underlying topology (i.e. configuration
 * space in path-planning terms). Then, the points in the current graph that are nearest
 * to the random point are picked for expansion. Then, edges (of a maximum length) are
 * added to the nearest vertex towards the random point while it is still possible to
 * add such an edge without leaving the free space (the part of the configuration space which
 * is not occupied by an obstacle). Finally, the last point of the expansion is used to attempt
 * the same expansion with the other nearest neighbors and the accumulated cost-to-go of all
 * these alternate paths are compared and the edge that leads to the shortest accumulated cost-to-go
 * is selected as a new edge in the tree. The algorithm will stop when either the number of vertices
 * in the tree has reached a maximum or when the user callback signals the stop.
 *
 * This library also provides the bidirectional version of the RRT* algorithm. In this version,
 * two trees are generated. Typically, one tree is initialized with the starting vertex and
 * the other is initialized with the goal vertex. The algorithm works to try and join the
 * two graphs as quickly as possible and with the most direct path. The algorithm alternates
 * between the two graphs. It first uses the normal procedure (as in the unidirectional variant)
 * to try and add a vertex to one graph. If it does not succeed (i.e. there was no free-space in
 * the expanded direction), it will try to expand the other graph. If it does succeed,
 * than the last vertex that was added to the tree will be the point towards which the other tree
 * will be expanded (if free-space permits). In other words, any successful vertex addition to one
 * tree causes the other to attempt an expansion in that direction, and any failure at adding a vertex
 * to one tree causes the other to attempt to expand towards a random point. This version of the
 * algorithm will also notify the user (via a visitor's callback) whenever one tree was successfully
 * expanded towards the other tree to the point that they meet at two vertices (one on each graph).
 * The user can thus record successful connections as paths and decide whether it's worth continuing
 * with the generation the Bi-RRT in the hopes of finding a better path with richer trees.
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

#ifndef REAK_PLANNING_GRAPH_ALG_RRT_STAR_H_
#define REAK_PLANNING_GRAPH_ALG_RRT_STAR_H_

#include <limits>
#include <tuple>
#include <utility>
#include "bagl/graph_concepts.h"
#include "bagl/property_map.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/random_sampler_concept.h"

#include "ReaK/planning/graph_alg/neighborhood_functors.h"
#include "ReaK/planning/graph_alg/node_generators.h"
#include "ReaK/planning/graph_alg/sbmp_visitor_concepts.h"

#include "ReaK/planning/graph_alg/branch_and_bound_connector.h"
#include "ReaK/planning/graph_alg/lazy_connector.h"
#include "ReaK/planning/graph_alg/pruned_connector.h"

#include "ReaK/math/optimization/optim_exceptions.h"

namespace ReaK::graph {

namespace rrt_detail {

template <typename Graph, typename RRGVisitor, typename PositionMap,
          typename WeightMap, typename DistanceMap, typename PredecessorMap,
          typename FwdDistanceMap = infinite_double_value_prop_map,
          typename SuccessorMap = null_vertex_prop_map<Graph>>
struct rrt_conn_visitor {

  rrt_conn_visitor(RRGVisitor vis, PositionMap pos, WeightMap weight,
                   DistanceMap dist, PredecessorMap pred,
                   FwdDistanceMap fwd_dist = FwdDistanceMap(),
                   SuccessorMap succ = SuccessorMap())
      : vis_(vis),
        position_(pos),
        weight_(weight),
        distance_(dist),
        predecessor_(pred),
        fwd_distance_(fwd_dist),
        successor_(succ){};

  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;
  using Edge = bagl::graph_edge_descriptor_t<Graph>;

  template <typename PositionValue>
  Vertex create_vertex(const PositionValue& p, Graph& g) const {
    bagl::vertex_property_type<Graph> up;
    put(position_, up, p);
    put(distance_, up, std::numeric_limits<double>::infinity());
    put(predecessor_, up, bagl::graph_traits<Graph>::null_vertex());
    put(fwd_distance_, up, std::numeric_limits<double>::infinity());
    put(successor_, up, bagl::graph_traits<Graph>::null_vertex());
    Vertex u = add_vertex(g, std::move(up));
    vis_.vertex_added(u, g);
    return u;
  }

  void vertex_to_be_removed(Vertex u, Graph& g) const {
    vis_.vertex_to_be_removed(u, g);
  }

  void vertex_added(Vertex v, Graph& g) const { vis_.vertex_added(v, g); }
  void edge_added(Edge e, Graph& g) const { vis_.edge_added(e, g); }

  void travel_explored(Vertex /*unused*/, Vertex /*unused*/,
                       Graph& /*unused*/) const {}
  void travel_succeeded(Vertex /*unused*/, Vertex /*unused*/,
                        Graph& /*unused*/) const {}
  void travel_failed(Vertex /*unused*/, Vertex /*unused*/,
                     Graph& /*unused*/) const {}
  void affected_vertex(Vertex /*unused*/, Graph& /*unused*/) const {}

  bool keep_going() const { return vis_.keep_going(); }

  auto can_be_connected(Vertex u, Vertex v, Graph& g) const {
    return vis_.can_be_connected(u, v, g);
  }

  template <typename PositionValue>
  auto steer_towards_position(const PositionValue& p, Vertex u,
                              Graph& g) const {
    return vis_.steer_towards_position(p, u, g);
  }

  template <typename PositionValue>
  auto steer_back_to_position(const PositionValue& p, Vertex u,
                              Graph& g) const {
    return vis_.steer_back_to_position(p, u, g);
  }

  RRGVisitor vis_;
  PositionMap position_;
  // needed by generate_rrt_star_loop (given to connector call)
  WeightMap weight_;
  DistanceMap distance_;
  PredecessorMap predecessor_;
  FwdDistanceMap fwd_distance_;
  SuccessorMap successor_;
};

template <typename Graph, pp::MetricSpace Space, typename RRTStarConnVisitor,
          typename MotionGraphConnector, typename PositionMap,
          typename NodeGenerator, typename NcSelector>
void generate_rrt_star_loop(Graph& g, const Space& super_space,
                            RRTStarConnVisitor conn_vis,
                            MotionGraphConnector connect_vertex,
                            PositionMap position,
                            NodeGenerator node_generator_func,
                            NcSelector select_neighborhood) {
  while (conn_vis.keep_going()) {
    auto [x_near, p_new, eprop] = node_generator_func(
        g, conn_vis,
        bagl::composite_property_map(position, get(bagl::vertex_all, g)));

    if ((x_near != bagl::graph_traits<Graph>::null_vertex()) &&
        (get(conn_vis.distance_, get_property(g, x_near)) !=
         std::numeric_limits<double>::infinity())) {
      connect_vertex(p_new, x_near, eprop, g, super_space, conn_vis,
                     conn_vis.position_, conn_vis.distance_,
                     conn_vis.predecessor_, conn_vis.weight_,
                     select_neighborhood);
    }
  }
}

template <typename Graph, pp::MetricSpace Space, typename RRTStarConnVisitor,
          typename MotionGraphConnector, typename PositionMap,
          typename NodeGenerator, typename NcSelector>
void generate_rrt_star_bidir_loop(Graph& g, const Space& super_space,
                                  RRTStarConnVisitor conn_vis,
                                  MotionGraphConnector connect_vertex,
                                  PositionMap position,
                                  NodeGenerator node_generator_func,
                                  NcSelector select_neighborhood) {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;
  using EdgeProp = bagl::edge_bundle_type<Graph>;

  while (conn_vis.keep_going()) {
    auto [x_near_pred, p_new_pred, ep_pred, x_near_succ, p_new_succ, ep_succ] =
        node_generator_func(
            g, conn_vis,
            bagl::composite_property_map(position, get(bagl::vertex_all, g)));

    if (x_near_pred != bagl::graph_traits<Graph>::null_vertex()) {
      Vertex x_near_other = bagl::graph_traits<Graph>::null_vertex();
      EdgeProp ep_other;
      connect_vertex(p_new_pred, x_near_pred, ep_pred, x_near_other, ep_other,
                     g, super_space, conn_vis, conn_vis.position_,
                     conn_vis.distance_, conn_vis.predecessor_,
                     conn_vis.fwd_distance_, conn_vis.successor_,
                     conn_vis.weight_, select_neighborhood);
    }
    if (x_near_succ != bagl::graph_traits<Graph>::null_vertex()) {
      Vertex x_near_other = bagl::graph_traits<Graph>::null_vertex();
      EdgeProp ep_other;
      connect_vertex(p_new_succ, x_near_other, ep_other, x_near_succ, ep_succ,
                     g, super_space, conn_vis, conn_vis.position_,
                     conn_vis.distance_, conn_vis.predecessor_,
                     conn_vis.fwd_distance_, conn_vis.successor_,
                     conn_vis.weight_, select_neighborhood);
    }
  }
}

template <typename Graph, pp::MetricSpace Space, typename RRGVisitor,
          typename PositionMap, typename DistanceMap, typename PredecessorMap,
          typename WeightMap, typename NodeGenerator, typename NcSelector>
void generate_rrt_star_loop(Graph& g, const Space& super_space, RRGVisitor vis,
                            PositionMap position, DistanceMap distance,
                            PredecessorMap pred, WeightMap weight,
                            NodeGenerator node_generator_func,
                            NcSelector select_neighborhood) {
  rrt_conn_visitor<Graph, RRGVisitor, PositionMap, WeightMap, DistanceMap,
                   PredecessorMap>
      conn_vis(vis, position, weight, distance, pred);

  generate_rrt_star_loop(g, super_space, conn_vis, lazy_node_connector(),
                         position, node_generator_func, select_neighborhood);
}

template <typename Graph, pp::MetricSpace Space, typename RRGVisitor,
          typename PositionMap, typename DistanceMap, typename PredecessorMap,
          typename FwdDistanceMap, typename SuccessorMap, typename WeightMap,
          typename NodeGenerator, typename NcSelector>
void generate_rrt_star_bidir_loop(Graph& g, const Space& super_space,
                                  RRGVisitor vis, PositionMap position,
                                  DistanceMap distance, PredecessorMap pred,
                                  FwdDistanceMap fwd_distance,
                                  SuccessorMap succ, WeightMap weight,
                                  NodeGenerator node_generator_func,
                                  NcSelector select_neighborhood) {
  rrt_conn_visitor<Graph, RRGVisitor, PositionMap, WeightMap, DistanceMap,
                   PredecessorMap, FwdDistanceMap, SuccessorMap>
      conn_vis(vis, position, weight, distance, pred, fwd_distance, succ);

  generate_rrt_star_loop(g, super_space, conn_vis, lazy_node_connector(),
                         position, node_generator_func, select_neighborhood);
}

}  // namespace rrt_detail

template <typename Graph, pp::MetricSpace Space, typename RRTStarVisitor,
          typename NcSelector, typename PositionMap, typename WeightMap,
          typename DistanceMap, typename PredecessorMap,
          typename FwdDistanceMap = infinite_double_value_prop_map,
          typename SuccessorMap = null_vertex_prop_map<Graph>>
struct rrtstar_bundle {
  using Vertex = bagl::graph_vertex_descriptor_t<Graph>;

  Graph* g_;
  Vertex start_vertex_ = bagl::graph_traits<Graph>::null_vertex();
  Vertex goal_vertex_ = bagl::graph_traits<Graph>::null_vertex();
  const Space* super_space_;
  RRTStarVisitor vis_;
  NcSelector select_neighborhood_;
  PositionMap position_;
  WeightMap weight_;
  DistanceMap distance_;
  PredecessorMap predecessor_;
  FwdDistanceMap fwd_distance_;
  SuccessorMap successor_;

  rrtstar_bundle(Graph& g, Vertex start_vertex, const Space& super_space,
                 RRTStarVisitor vis, NcSelector select_neighborhood,
                 PositionMap position, WeightMap weight, DistanceMap distance,
                 PredecessorMap predecessor,
                 FwdDistanceMap fwd_distance = FwdDistanceMap(),
                 SuccessorMap successor = SuccessorMap())
      : g_(&g),
        start_vertex_(start_vertex),
        goal_vertex_(bagl::graph_traits<Graph>::null_vertex()),
        super_space_(&super_space),
        vis_(vis),
        select_neighborhood_(select_neighborhood),
        position_(position),
        weight_(weight),
        distance_(distance),
        predecessor_(predecessor),
        fwd_distance_(fwd_distance),
        successor_(successor) {}

  rrtstar_bundle(Graph& g, Vertex start_vertex, Vertex goal_vertex,
                 const Space& super_space, RRTStarVisitor vis,
                 NcSelector select_neighborhood, PositionMap position,
                 WeightMap weight, DistanceMap distance,
                 PredecessorMap predecessor,
                 FwdDistanceMap fwd_distance = FwdDistanceMap(),
                 SuccessorMap successor = SuccessorMap())
      : g_(&g),
        start_vertex_(start_vertex),
        goal_vertex_(goal_vertex),
        super_space_(&super_space),
        vis_(vis),
        select_neighborhood_(select_neighborhood),
        position_(position),
        weight_(weight),
        distance_(distance),
        predecessor_(predecessor),
        fwd_distance_(fwd_distance),
        successor_(successor) {}
};

/**
  * This function template creates a bundle of parameters to be fed to any of the
  * RRT* algorithms. This is mainly to simply the interface and the code of all these
  * different variants of the RRT* algorithm.
  * \tparam Graph The graph type that can store the generated roadmap, should model
  *         bagl::concepts::BidirectionalGraph and bagl::concepts::MutableGraph.
  * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
  * \tparam Space The topology type that represents the free-space, should model BGL's Topology concept.
  * \tparam Visitor The type of the RRT* visitor to be used, should model the RRGVisitorConcept.
  * \tparam PositionMap A property-map type that can store the position of each vertex-property object.
  * \tparam WeightMap This property-map type is used to store the weights of the edge-properties of the
  *         graph (cost of travel along an edge).
  * \tparam DistanceMap This property-map type is used to store the estimated distance of each vertex-property object
  *         to the goal.
  * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting
  *         vertex-property object together with its optimal predecessor.
  * \tparam NcSelector A functor type that can select a list of vertices of the graph that are
  *         the nearest-neighbors of a given vertex (or some other heuristic to select the neighbors).
  *         See classes in the topological_search.hpp header-file.
  *
  * \param g A mutable graph that should initially store the starting
  *        vertex (if not it will be randomly generated) and will store
  *        the generated graph once the algorithm has finished.
  * \param start_vertex The starting point of the algorithm, on the graph.
  * \param super_space A topology (as defined by the Born Again Graph Library). This topology
  *        should not include collision checking in its distance metric.
  * \param vis A RRT* visitor implementing the RRGVisitorConcept. This is the
  *        main point of customization and recording of results that the
  *        user can implement.
  * \param select_neighborhood A callable object (functor) that can select a list of
  *        vertices of the graph that ought to be connected to a new
  *        vertex. The list should be sorted in order of increasing "distance".
  * \param position A mapping that implements the MutablePropertyMap Concept. Also,
  *        the value_type of this map should be the same type as the topology's
  *        value_type.
  * \param weight The property-map which stores the weight of each edge-property object (the cost of travel
  *        along the edge).
  * \param distance The property-map which stores the estimated distance of each vertex to the goal.
  * \param predecessor The property-map which will store the resulting path by connecting
  *        vertices together with their optimal predecessor (follow in reverse to discover the
  *        complete path).
  */
template <bagl::concepts::MutableGraph Graph, pp::MetricSpace Space,
          RRGVisitor<Graph, Space> Visitor, typename NcSelector,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap>
rrtstar_bundle<Graph, Space, Visitor, NcSelector, PositionMap, WeightMap,
               DistanceMap, PredecessorMap>
make_rrtstar_bundle(Graph& g,
                    bagl::graph_vertex_descriptor_t<Graph> start_vertex,
                    const Space& super_space, Visitor vis,
                    NcSelector select_neighborhood, PositionMap position,
                    WeightMap weight, DistanceMap distance,
                    PredecessorMap predecessor) {
  return rrtstar_bundle<Graph, Space, Visitor, NcSelector, PositionMap,
                        WeightMap, DistanceMap, PredecessorMap>(
      g, start_vertex, super_space, vis, select_neighborhood, position, weight,
      distance, predecessor);
}

template <bagl::concepts::MutableGraph Graph, pp::MetricSpace Space,
          RRGVisitor<Graph, Space> Visitor, typename NcSelector,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap>
rrtstar_bundle<Graph, Space, Visitor, NcSelector, PositionMap, WeightMap,
               DistanceMap, PredecessorMap>
make_rrtstar_bundle(Graph& g,
                    bagl::graph_vertex_descriptor_t<Graph> start_vertex,
                    bagl::graph_vertex_descriptor_t<Graph> goal_vertex,
                    const Space& super_space, Visitor vis,
                    NcSelector select_neighborhood, PositionMap position,
                    WeightMap weight, DistanceMap distance,
                    PredecessorMap predecessor) {
  return rrtstar_bundle<Graph, Space, Visitor, NcSelector, PositionMap,
                        WeightMap, DistanceMap, PredecessorMap>(
      g, start_vertex, goal_vertex, super_space, vis, select_neighborhood,
      position, weight, distance, predecessor);
}

/**
  * This function template creates a bundle of parameters to be fed to any of the
  * RRT* algorithms. This is mainly to simply the interface and the code of all these
  * different variants of the RRT* algorithm.
  * \tparam Graph The graph type that can store the generated roadmap, should model
  *         bagl::concepts::BidirectionalGraph and bagl::concepts::MutableGraph.
  * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
  * \tparam Space The topology type that represents the free-space, should model BGL's Topology concept.
  * \tparam Visitor The type of the RRT* visitor to be used.
  * \tparam NcSelector A functor type that can select a list of vertices of the graph that are
  *         the nearest-neighbors of a given vertex (or some other heuristic to select the neighbors).
  *         See classes in the topological_search.hpp header-file.
  * \tparam PositionMap A property-map type that can store the position of each vertex-property object.
  * \tparam WeightMap This property-map type is used to store the weights of the edge-properties of the
  *         graph (cost of travel along an edge).
  * \tparam DistanceMap This property-map type is used to store the current best distance of each vertex-property object
  *         from the start.
  * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting
  *         vertex-property object together with its optimal predecessor.
  * \tparam FwdDistanceMap This property-map type is used to store the current best distance of each vertex-property
  *object
  *         to the goal.
  * \tparam SuccessorMap This property-map type is used to store the resulting path by connecting
  *         vertex-property object together with its optimal successor.
  *
  * \param g A mutable graph that should initially store the starting
  *        vertex (if not it will be randomly generated) and will store
  *        the generated graph once the algorithm has finished.
  * \param start_vertex The starting point of the algorithm, on the graph.
  * \param goal_vertex The goal point of the algorithm, on the graph.
  * \param super_space A topology (as defined by the Born Again Graph Library). This topology
  *        should not include collision checking in its distance metric.
  * \param vis A RRT* visitor implementing the RRGVisitorConcept. This is the
  *        main point of customization and recording of results that the
  *        user can implement.
  * \param select_neighborhood A callable object (functor) that can select a list of
  *        vertices of the graph that ought to be connected to a new
  *        vertex. The list should be sorted in order of increasing "distance".
  * \param position A mapping that implements the MutablePropertyMap Concept. Also,
  *        the value_type of this map should be the same type as the topology's
  *        value_type.
  * \param weight The property-map which stores the weight of each edge-property object (the cost of travel
  *        along the edge).
  * \param distance The property-map which stores the estimated distance of each vertex from the start.
  * \param predecessor The property-map which will store the resulting path by connecting
  *        vertices together with their optimal predecessor (follow in reverse to discover the
  *        complete path).
  * \param fwd_distance The property-map which stores the estimated distance of each vertex to the goal.
  * \param successor The property-map which will store the resulting path by connecting
  *        vertices together with their optimal successor (follow in order to discover the
  *        remaining path to the goal).
  */

template <bagl::concepts::MutableGraph Graph, pp::MetricSpace Space,
          RRGBidirVisitor<Graph, Space> Visitor, typename NcSelector,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> FwdDistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> SuccessorMap>
rrtstar_bundle<Graph, Space, Visitor, NcSelector, PositionMap, WeightMap,
               DistanceMap, PredecessorMap, FwdDistanceMap, SuccessorMap>
make_rrtstar_bundle(Graph& g,
                    bagl::graph_vertex_descriptor_t<Graph> start_vertex,
                    bagl::graph_vertex_descriptor_t<Graph> goal_vertex,
                    const Space& super_space, Visitor vis,
                    NcSelector select_neighborhood, PositionMap position,
                    WeightMap weight, DistanceMap distance,
                    PredecessorMap predecessor, FwdDistanceMap fwd_distance,
                    SuccessorMap successor) {
  return rrtstar_bundle<Graph, Space, Visitor, NcSelector, PositionMap,
                        WeightMap, DistanceMap, PredecessorMap, FwdDistanceMap,
                        SuccessorMap>(
      g, start_vertex, goal_vertex, super_space, vis, select_neighborhood,
      position, weight, distance, predecessor, fwd_distance, successor);
}

/**
  * This function template is the RRT* algorithm (refer to rrt_star.hpp dox).
  * \tparam Graph A mutable graph type that will represent the generated tree, should model
  *bagl::concepts::VertexListGraph and bagl::concepts::MutableGraph
  * \tparam Space A topology type that will represent the space in which the configurations (or positions) exist,
  *should model BGL's Topology concept
  * \tparam Visitor An RRT* visitor type that implements the customizations to this RRT* algorithm.
  * \tparam PositionMap A property-map type that can store the configurations (or positions) of the vertices.
  * \tparam DistanceMap This property-map type is used to store the estimated cost-to-go of each vertex to the start (or
  *goal).
  * \tparam PredecessorMap This property-map type is used to store the predecessor of each vertex.
  * \tparam WeightMap This property-map type is used to store the weights of the edges of the graph (cost of travel
  *along an edge).
  * \tparam Sampler This is a random-sampler over the topology.
  * \tparam NcSelector A functor type which can perform a neighborhood search of a point to a graph in the topology (see
  *topological_search.hpp).
  * \param g A mutable graph that should initially store the starting and goal
  *        vertex and will store the generated graph once the algorithm has finished.
  * \param super_space A topology (as defined by the Born Again Graph Library). Note
  *        that it should represent the entire configuration space (not collision-free space).
  * \param vis A RRT* visitor implementing the RRGVisitor. This is the
  *        main point of customization and recording of results that the
  *        user can implement.
  * \param position A mapping that implements the MutablePropertyMap Concept. Also,
  *        the value_type of this map should be the same type as the topology's
  *        value_type.
  * \param distance The property-map which stores the estimated cost-to-go of each vertex to the start (or goal).
  * \param pred The property-map which stores the predecessor of each vertex.
  * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
  *        along the edge).
  * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
  * \param select_neighborhood A callable object (functor) which can perform a
  *        nearest neighbor search of a point to a graph in the topology. (see star_neighborhood)
  *
  */
template <bagl::concepts::MutableGraph Graph, pp::MetricSpace Space,
          RRGVisitor<Graph, Space> Visitor,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap,
          bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
          pp::RandomSampler<Space> Sampler, typename NcSelector>
void generate_rrt_star(Graph& g, const Space& super_space, Visitor vis,
                       PositionMap position, DistanceMap distance,
                       PredecessorMap pred, WeightMap weight,
                       Sampler get_sample, NcSelector select_neighborhood) {
  if (num_vertices(g) == 0) {
    throw optim::infeasible_problem(
        "Cannot solve a RRT* problem without start position!");
  }

  rrt_detail::rrt_conn_visitor<Graph, Visitor, PositionMap, WeightMap,
                               DistanceMap, PredecessorMap>
      conn_vis(vis, position, weight, distance, pred);

  rrt_detail::generate_rrt_star_loop(
      g, super_space, conn_vis, lazy_node_connector(), position,
      rrg_node_generator<Space, Sampler, NcSelector>(&super_space, get_sample,
                                                     select_neighborhood),
      select_neighborhood);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the RRT* algorithm, without initialization of the existing graph.
  * \tparam RRTStarBundle A RRT* bundle type (see make_rrtstar_bundle()).
  * \tparam Sampler This is a random-sampler over the topology.
  * \param bdl A const-reference to a RRT* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
  */
template <typename RRTStarBundle, typename Sampler>
void generate_rrt_star(const RRTStarBundle& bdl, Sampler get_sample) {
  put(bdl.distance_, get_property(*bdl.g_, bdl.start_vertex_), 0.0);
  put(bdl.predecessor_, get_property(*bdl.g_, bdl.start_vertex_),
      bdl.start_vertex_);

  generate_rrt_star(*bdl.g_, *bdl.super_space_, bdl.vis_, bdl.position_,
                    bdl.distance_, bdl.predecessor_, bdl.weight_, get_sample,
                    bdl.select_neighborhood_);
}

/**
  * This function template is the Bi-directional RRT* algorithm (refer to rrt_star_bidir.hpp dox).
  * \tparam Graph A mutable graph type that will represent the generated tree, should model
  *bagl::concepts::VertexListGraph and bagl::concepts::MutableGraph
  * \tparam Space A topology type that will represent the space in which the configurations (or positions) exist,
  *should model BGL's Topology concept
  * \tparam Visitor An RRT* visitor type that implements the customizations to this RRT* algorithm.
  * \tparam PositionMap A property-map type that can store the configurations (or positions) of the vertices.
  * \tparam DistanceMap This property-map type is used to store the estimated cost-to-go of each vertex from the start.
  * \tparam PredecessorMap This property-map type is used to store the predecessor of each vertex.
  * \tparam FwdDistanceMap This property-map type is used to store the estimated cost-to-go of each vertex to the goal.
  * \tparam SuccessorMap This property-map type is used to store the successor of each vertex.
  * \tparam WeightMap This property-map type is used to store the weights of the edges of the graph (cost of travel
  *along an edge).
  * \tparam Sampler This is a random-sampler over the topology.
  * \tparam NcSelector A functor type which can perform a neighborhood search of a point to a graph in the topology (see
  *topological_search.hpp).
  * \param g A mutable graph that should initially store the starting and goal
  *        vertex and will store the generated graph once the algorithm has finished.
  * \param super_space A topology (as defined by the Born Again Graph Library). Note
  *        that it should represent the entire configuration space (not collision-free space).
  * \param vis A RRT* visitor. This is the
  *        main point of customization and recording of results that the
  *        user can implement.
  * \param position A mapping that implements the MutablePropertyMap Concept. Also,
  *        the value_type of this map should be the same type as the topology's
  *        value_type.
  * \param distance The property-map which stores the estimated cost-to-go of each vertex from the start.
  * \param pred The property-map which stores the predecessor of each vertex.
  * \param fwd_distance The property-map which stores the estimated cost-to-go of each vertex to the goal.
  * \param succ The property-map which stores the successor of each vertex.
  * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
  *        along the edge).
  * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
  * \param select_neighborhood A callable object (functor) which can perform a
  *        nearest neighbor search of a point to a graph in the topology. (see star_neighborhood)
  *
  */
template <bagl::concepts::MutableGraph Graph, pp::MetricSpace Space,
          RRGBidirVisitor<Graph, Space> Visitor,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> FwdDistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> SuccessorMap,
          bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
          pp::RandomSampler<Space> Sampler, typename NcSelector>
void generate_rrt_star_bidir(Graph& g, const Space& super_space, Visitor vis,
                             PositionMap position, DistanceMap distance,
                             PredecessorMap pred, FwdDistanceMap fwd_distance,
                             SuccessorMap succ, WeightMap weight,
                             Sampler get_sample,
                             NcSelector select_neighborhood) {
  if (num_vertices(g) == 0) {
    throw optim::infeasible_problem(
        "Cannot solve a bi-directional RRT* problem without start and goal "
        "positions!");
  }

  rrt_detail::rrt_conn_visitor<Graph, Visitor, PositionMap, WeightMap,
                               DistanceMap, PredecessorMap, FwdDistanceMap,
                               SuccessorMap>
      conn_vis(vis, position, weight, distance, pred, fwd_distance, succ);

  rrt_detail::generate_rrt_star_bidir_loop(
      g, super_space, conn_vis, lazy_node_connector(), position,
      rrg_bidir_generator(&super_space, get_sample, select_neighborhood, pred,
                          succ),
      select_neighborhood);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional RRT* algorithm, without initialization of the existing graph.
  * \tparam RRTStarBundle A RRT* bundle type (see make_rrtstar_bundle()).
  * \tparam Sampler This is a random-sampler over the topology.
  * \param bdl A const-reference to a RRT* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
  */
template <typename RRTStarBundle, typename Sampler>
void generate_rrt_star_bidir(const RRTStarBundle& bdl, Sampler get_sample) {
  using Graph = std::decay_t<decltype(*(bdl.g_))>;

  put(bdl.distance_, get_property(*bdl.g_, bdl.start_vertex_), 0.0);
  put(bdl.predecessor_, get_property(*bdl.g_, bdl.start_vertex_),
      bdl.start_vertex_);
  if (bdl.goal_vertex_ != bagl::graph_traits<Graph>::null_vertex()) {
    put(bdl.fwd_distance_, get_property(*bdl.g_, bdl.goal_vertex_), 0.0);
    put(bdl.successor_, get_property(*bdl.g_, bdl.goal_vertex_),
        bdl.goal_vertex_);
  }

  generate_rrt_star_bidir(*bdl.g_, *bdl.super_space_, bdl.vis_, bdl.position_,
                          bdl.distance_, bdl.predecessor_, bdl.fwd_distance_,
                          bdl.successor_, bdl.weight_, get_sample,
                          bdl.select_neighborhood_);
}

/**
  * This function template is the RRT* algorithm (refer to rrt_star.hpp dox).
  * This function uses a branch-and-bound heuristic to limit the number of nodes.
  * \tparam Graph A mutable graph type that will represent the generated tree, should model
  *bagl::concepts::VertexListGraph and bagl::concepts::MutableGraph
  * \tparam Space A topology type that will represent the space in which the configurations (or positions) exist,
  *should model BGL's Space concept
  * \tparam RRGVisitor An RRT* visitor type that implements the customizations to this RRT* algorithm.
  * \tparam PositionMap A property-map type that can store the configurations (or positions) of the vertices.
  * \tparam DistanceMap This property-map type is used to store the estimated cost-to-go of each vertex to the start (or
  *goal).
  * \tparam PredecessorMap This property-map type is used to store the predecessor of each vertex.
  * \tparam WeightMap This property-map type is used to store the weights of the edges of the graph (cost of travel
  *along an edge).
  * \tparam Sampler This is a random-sampler over the topology.
  * \tparam NcSelector A functor type which can perform a neighborhood search of a point to a graph in the topology (see
  *topological_search.hpp).
  * \param g A mutable graph that should initially store the starting and goal
  *        vertex and will store the generated graph once the algorithm has finished.
  * \param start_vertex The vertex from which the motion-graph is grown.
  * \param goal_vertex The vertex which we want to connect to the motion-graph.
  * \param super_space A topology (as defined by the Born Again Graph Library). Note
  *        that it should represent the entire configuration space (not collision-free space).
  * \param vis A RRT* visitor. This is the
  *        main point of customization and recording of results that the
  *        user can implement.
  * \param position A mapping that implements the MutablePropertyMap Concept. Also,
  *        the value_type of this map should be the same type as the topology's
  *        value_type.
  * \param distance The property-map which stores the estimated cost-to-go of each vertex to the start (or goal).
  * \param pred The property-map which stores the predecessor of each vertex.
  * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
  *        along the edge).
  * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
  * \param select_neighborhood A callable object (functor) which can perform a
  *        nearest neighbor search of a point to a graph in the topology. (see star_neighborhood)
  *
  */
template <bagl::concepts::MutableGraph Graph, pp::MetricSpace Space,
          RRGVisitor<Graph, Space> Visitor,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap,
          bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
          pp::RandomSampler<Space> Sampler, typename NcSelector>
void generate_bnb_rrt_star(Graph& g,
                           bagl::graph_vertex_descriptor_t<Graph> start_vertex,
                           bagl::graph_vertex_descriptor_t<Graph> goal_vertex,
                           const Space& super_space, Visitor vis,
                           PositionMap position, DistanceMap distance,
                           PredecessorMap pred, WeightMap weight,
                           Sampler get_sample, NcSelector select_neighborhood) {
  if ((num_vertices(g) == 0) ||
      (start_vertex == bagl::graph_traits<Graph>::null_vertex()) ||
      (goal_vertex == bagl::graph_traits<Graph>::null_vertex())) {
    generate_rrt_star(g, super_space, vis, position, distance, pred, weight,
                      get_sample, select_neighborhood);
    return;
  }

  rrt_detail::rrt_conn_visitor<Graph, Visitor, PositionMap, WeightMap,
                               DistanceMap, PredecessorMap>
      conn_vis(vis, position, weight, distance, pred);

  bnb_ordering_data<Graph> bnb_data(g, start_vertex, goal_vertex);

  rrt_detail::generate_rrt_star_loop(
      g, super_space, conn_vis, bnb_connector<Graph>(bnb_data), position,
      rrg_node_generator<Space, Sampler, NcSelector>(&super_space, get_sample,
                                                     select_neighborhood),
      select_neighborhood);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the RRT* algorithm, without initialization of the existing graph.
  * This function uses a branch-and-bound heuristic to limit the number of nodes.
  * \tparam RRTStarBundle A RRT* bundle type (see make_rrtstar_bundle()).
  * \tparam Sampler This is a random-sampler over the topology.
  * \param bdl A const-reference to a RRT* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
  */
template <typename RRTStarBundle, typename Sampler>
void generate_bnb_rrt_star(const RRTStarBundle& bdl, Sampler get_sample) {
  put(bdl.distance_, get_property(*bdl.g_, bdl.start_vertex_), 0.0);
  put(bdl.predecessor_, get_property(*bdl.g_, bdl.start_vertex_),
      bdl.start_vertex_);

  generate_bnb_rrt_star(*bdl.g_, bdl.start_vertex_, bdl.goal_vertex_,
                        *bdl.super_space_, bdl.vis_, bdl.position_,
                        bdl.distance_, bdl.predecessor_, bdl.weight_,
                        get_sample, bdl.select_neighborhood_);
}

/**
  * This function template is the Bi-directional RRT* algorithm (refer to rrt_star_bidir.hpp dox).
  * This function uses a branch-and-bound heuristic to limit the number of nodes.
  * \tparam Graph A mutable graph type that will represent the generated tree, should model
  *bagl::concepts::VertexListGraph and bagl::concepts::MutableGraph
  * \tparam Space A topology type that will represent the space in which the configurations (or positions) exist,
  *should model BGL's Topology concept
  * \tparam Visitor An RRT* visitor type that implements the customizations to this RRT* algorithm.
  * \tparam PositionMap A property-map type that can store the configurations (or positions) of the vertices.
  * \tparam DistanceMap This property-map type is used to store the estimated cost-to-go of each vertex from the start.
  * \tparam PredecessorMap This property-map type is used to store the predecessor of each vertex.
  * \tparam FwdDistanceMap This property-map type is used to store the estimated cost-to-go of each vertex to the goal.
  * \tparam SuccessorMap This property-map type is used to store the successor of each vertex.
  * \tparam WeightMap This property-map type is used to store the weights of the edges of the graph (cost of travel
  *along an edge).
  * \tparam Sampler This is a random-sampler over the topology.
  * \tparam NcSelector A functor type which can perform a neighborhood search of a point to a graph in the topology (see
  *topological_search.hpp).
  * \param g A mutable graph that should initially store the starting and goal
  *        vertex and will store the generated graph once the algorithm has finished.
  * \param start_vertex The vertex from which the motion-graph is grown.
  * \param goal_vertex The vertex which we want to connect to the motion-graph.
  * \param super_space A topology (as defined by the Born Again Graph Library). Note
  *        that it should represent the entire configuration space (not collision-free space).
  * \param vis A RRT* visitor. This is the
  *        main point of customization and recording of results that the
  *        user can implement.
  * \param position A mapping that implements the MutablePropertyMap Concept. Also,
  *        the value_type of this map should be the same type as the topology's
  *        value_type.
  * \param distance The property-map which stores the estimated cost-to-go of each vertex from the start.
  * \param pred The property-map which stores the predecessor of each vertex.
  * \param fwd_distance The property-map which stores the estimated cost-to-go of each vertex to the goal.
  * \param succ The property-map which stores the successor of each vertex.
  * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
  *        along the edge).
  * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
  * \param select_neighborhood A callable object (functor) which can perform a
  *        nearest neighbor search of a point to a graph in the topology. (see star_neighborhood)
  *
  */
template <bagl::concepts::MutableGraph Graph, pp::MetricSpace Space,
          RRGBidirVisitor<Graph, Space> Visitor,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PositionMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> DistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> PredecessorMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> FwdDistanceMap,
          bagl::concepts::ReadWriteVPropMemberMap<Graph> SuccessorMap,
          bagl::concepts::ReadableEPropMemberMap<Graph> WeightMap,
          pp::RandomSampler<Space> Sampler, typename NcSelector>
void generate_bnb_rrt_star_bidir(
    Graph& g, bagl::graph_vertex_descriptor_t<Graph> start_vertex,
    bagl::graph_vertex_descriptor_t<Graph> goal_vertex,
    const Space& super_space, Visitor vis, PositionMap position,
    DistanceMap distance, PredecessorMap pred, FwdDistanceMap fwd_distance,
    SuccessorMap succ, WeightMap weight, Sampler get_sample,
    NcSelector select_neighborhood) {
  if ((num_vertices(g) == 0) ||
      (start_vertex == bagl::graph_traits<Graph>::null_vertex()) ||
      (goal_vertex == bagl::graph_traits<Graph>::null_vertex())) {
    generate_rrt_star_bidir(g, super_space, vis, position, distance, pred,
                            fwd_distance, succ, weight, get_sample,
                            select_neighborhood);
    return;
  }

  rrt_detail::rrt_conn_visitor<Graph, Visitor, PositionMap, WeightMap,
                               DistanceMap, PredecessorMap, FwdDistanceMap,
                               SuccessorMap>
      conn_vis(vis, position, weight, distance, pred, fwd_distance, succ);

  bnb_ordering_data<Graph> bnb_data(g, start_vertex, goal_vertex);

  rrt_detail::generate_rrt_star_bidir_loop(
      g, super_space, conn_vis, bnb_connector<Graph>(bnb_data), position,
      rrg_bidir_generator(&super_space, get_sample, select_neighborhood, pred,
                          succ),
      select_neighborhood);
}

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the RRT* algorithm, without initialization of the existing graph.
  * This function uses a branch-and-bound heuristic to limit the number of nodes.
  * \tparam RRTStarBundle A RRT* bundle type (see make_rrtstar_bundle()).
  * \tparam Sampler This is a random-sampler over the topology.
  * \param bdl A const-reference to a RRT* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
  */
template <typename RRTStarBundle, typename Sampler>
void generate_bnb_rrt_star_bidir(const RRTStarBundle& bdl, Sampler get_sample) {
  using Graph = std::decay_t<decltype(*(bdl.g_))>;

  put(bdl.distance_, get_property(*bdl.g_, bdl.start_vertex_), 0.0);
  put(bdl.predecessor_, get_property(*bdl.g_, bdl.start_vertex_),
      bdl.start_vertex_);
  if (bdl.goal_vertex_ != bagl::graph_traits<Graph>::null_vertex()) {
    put(bdl.fwd_distance_, get_property(*bdl.g_, bdl.goal_vertex_), 0.0);
    put(bdl.successor_, get_property(*bdl.g_, bdl.goal_vertex_),
        bdl.goal_vertex_);
  }

  generate_bnb_rrt_star_bidir(*(bdl.g_), bdl.start_vertex_, bdl.goal_vertex_,
                              *(bdl.super_space_), bdl.vis_, bdl.position_,
                              bdl.distance_, bdl.predecessor_,
                              bdl.fwd_distance_, bdl.successor_, bdl.weight_,
                              get_sample, bdl.select_neighborhood_);
}

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_RRT_STAR_H_

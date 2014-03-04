/**
 * \file rrt_star.hpp
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

#ifndef REAK_RRT_STAR_HPP
#define REAK_RRT_STAR_HPP

#include <utility>
#include <limits>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>

#include "path_planning/metric_space_concept.hpp"
#include "path_planning/random_sampler_concept.hpp"

#include "sbmp_visitor_concepts.hpp"
#include "neighborhood_functors.hpp"
#include "node_generators.hpp"

#include "pruned_connector.hpp"
#include "lazy_connector.hpp"
#include "branch_and_bound_connector.hpp"

#include "optimization/optim_exceptions.hpp"

namespace ReaK {
  
namespace graph {

  
namespace detail {
  
  
  template <typename Graph, typename RRGVisitor, typename PositionMap, 
            typename WeightMap, typename DistanceMap, typename PredecessorMap,
            typename FwdDistanceMap = infinite_double_value_prop_map,  
            typename SuccessorMap = null_vertex_prop_map<Graph> >
  struct rrt_conn_visitor
  {
    
    rrt_conn_visitor(RRGVisitor vis, PositionMap pos, WeightMap weight, 
                     DistanceMap dist, PredecessorMap pred, 
                     FwdDistanceMap fwd_dist = FwdDistanceMap(),
                     SuccessorMap succ = SuccessorMap()) : 
                     m_vis(vis), m_position(pos), m_weight(weight), 
                     m_distance(dist), m_predecessor(pred),
                     m_fwd_distance(fwd_dist), m_successor(succ) { };
    
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::vertex_bundled VertexProp;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef typename Graph::edge_bundled EdgeProp;
    
    Vertex create_vertex(const PositionValue& p, Graph& g) const {
      
      VertexProp up;
      put(m_position, up, p);
      put(m_distance, up, std::numeric_limits<double>::infinity());
      put(m_predecessor, up, boost::graph_traits<Graph>::null_vertex());
      put(m_fwd_distance, up, std::numeric_limits<double>::infinity());
      put(m_successor, up, boost::graph_traits<Graph>::null_vertex());
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      Vertex u = add_vertex(std::move(up), g);
#else
      Vertex u = add_vertex(up, g);
#endif
      m_vis.vertex_added(u,g);
      
      return u;
    };
    
    void vertex_to_be_removed(Vertex u, Graph& g) const {
      m_vis.vertex_to_be_removed(u,g);
    };
    
    void vertex_added(Vertex v, Graph& g) const { m_vis.vertex_added(v,g); };
    void edge_added(Edge e, Graph& g) const { m_vis.edge_added(e,g); };
    
    void travel_explored(Vertex u, Vertex v, Graph& g) const { };
    void travel_succeeded(Vertex u, Vertex v, Graph& g) const { };
    void travel_failed(Vertex u, Vertex v, Graph& g) const { };
    void affected_vertex(Vertex, Graph&) const { };
    
    bool keep_going() const { return m_vis.keep_going(); };
    
    std::pair<bool, EdgeProp> can_be_connected(Vertex u, Vertex v, Graph& g) const { 
      return m_vis.can_be_connected(u,v,g);
    };
    
    boost::tuple<PositionValue, bool, EdgeProp> steer_towards_position(const PositionValue& p, Vertex u, Graph& g) const { 
      return m_vis.steer_towards_position(p, u, g);
    };
    
    boost::tuple<PositionValue, bool, EdgeProp> steer_back_to_position(const PositionValue& p, Vertex u, Graph& g) const { 
      return this->m_vis.steer_back_to_position(p, u, g);
    };
    
    RRGVisitor m_vis;
    PositionMap m_position;
    WeightMap m_weight;  // needed by generate_rrt_star_loop (given to connector call)
    DistanceMap m_distance;
    PredecessorMap m_predecessor;
    FwdDistanceMap m_fwd_distance;
    SuccessorMap m_successor;
  };
  
  
  template <typename Graph, typename Topology, typename RRTStarConnVisitor,
            typename MotionGraphConnector, typename PositionMap,
            typename NodeGenerator, typename NcSelector>
  void generate_rrt_star_loop(
    Graph& g, const Topology& super_space, RRTStarConnVisitor conn_vis,
    MotionGraphConnector connect_vertex, PositionMap position,
    NodeGenerator node_generator_func, NcSelector select_neighborhood) {
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    
    while(conn_vis.keep_going()) {
      
      PositionValue p_new; Vertex x_near; EdgeProp eprop;
      boost::tie(x_near, p_new, eprop) = node_generator_func(g, conn_vis, boost::bundle_prop_to_vertex_prop(position, g));
      
      if((x_near != boost::graph_traits<Graph>::null_vertex()) && 
          (get(conn_vis.m_distance, g[x_near]) != std::numeric_limits<double>::infinity())) {
        connect_vertex(p_new, x_near, eprop, g, 
                       super_space, conn_vis, conn_vis.m_position, 
                       conn_vis.m_distance, conn_vis.m_predecessor, 
                       conn_vis.m_weight, select_neighborhood);
      };
      
    };

  };
  
  template <typename Graph, typename Topology,
            typename RRTStarConnVisitor, typename MotionGraphConnector,
            typename PositionMap, typename NodeGenerator, typename NcSelector>
  void generate_rrt_star_bidir_loop(
      Graph& g, const Topology& super_space, RRTStarConnVisitor conn_vis,
      MotionGraphConnector connect_vertex, PositionMap position,
      NodeGenerator node_generator_func, NcSelector select_neighborhood) {
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    
    while(conn_vis.keep_going()) {
      
      Vertex x_near_pred, x_near_succ; 
      PositionValue p_new_pred, p_new_succ; 
      EdgeProp ep_pred, ep_succ;
      boost::tie(x_near_pred, p_new_pred, ep_pred, x_near_succ, p_new_succ, ep_succ) = 
        node_generator_func(g, conn_vis, boost::bundle_prop_to_vertex_prop(position, g));
      
      if( x_near_pred != boost::graph_traits<Graph>::null_vertex() ) {
        Vertex x_near_other = boost::graph_traits<Graph>::null_vertex();
        EdgeProp ep_other;
        connect_vertex(p_new_pred, x_near_pred, ep_pred, x_near_other, ep_other, g, 
                       super_space, conn_vis, conn_vis.m_position, 
                       conn_vis.m_distance, conn_vis.m_predecessor, 
                       conn_vis.m_fwd_distance, conn_vis.m_successor, 
                       conn_vis.m_weight, select_neighborhood);
      };
      if( x_near_succ != boost::graph_traits<Graph>::null_vertex() ) {
        Vertex x_near_other = boost::graph_traits<Graph>::null_vertex();
        EdgeProp ep_other;
        connect_vertex(p_new_succ, x_near_other, ep_other, x_near_succ, ep_succ, g, 
                       super_space, conn_vis, conn_vis.m_position, 
                       conn_vis.m_distance, conn_vis.m_predecessor, 
                       conn_vis.m_fwd_distance, conn_vis.m_successor, 
                       conn_vis.m_weight, select_neighborhood);
      };
      
    };

  };
  
  
  template <typename Graph, typename Topology, typename RRGVisitor,
            typename PositionMap, typename DistanceMap, typename PredecessorMap,
            typename WeightMap, typename NodeGenerator, typename NcSelector>
  void generate_rrt_star_loop(
    Graph& g, const Topology& super_space, RRGVisitor vis, PositionMap position, 
    DistanceMap distance, PredecessorMap pred,
    WeightMap weight, NodeGenerator node_generator_func, NcSelector select_neighborhood) {
    
    rrt_conn_visitor<Graph, RRGVisitor, PositionMap, WeightMap, DistanceMap, PredecessorMap> 
      conn_vis(vis, position, weight, distance, pred);
    
    generate_rrt_star_loop(g, super_space, conn_vis, lazy_node_connector(),
                           position, node_generator_func, select_neighborhood);
    
  };
  
  template <typename Graph, typename Topology, typename RRGVisitor, typename PositionMap,
            typename DistanceMap, typename PredecessorMap,
            typename FwdDistanceMap, typename SuccessorMap,
            typename WeightMap, typename NodeGenerator, typename NcSelector>
  void generate_rrt_star_bidir_loop(
      Graph& g, const Topology& super_space, RRGVisitor vis, PositionMap position,
      DistanceMap distance, PredecessorMap pred, FwdDistanceMap fwd_distance, SuccessorMap succ,
      WeightMap weight, NodeGenerator node_generator_func, NcSelector select_neighborhood) {
    
    rrt_conn_visitor<Graph, RRGVisitor, PositionMap, WeightMap, 
                     DistanceMap, PredecessorMap, FwdDistanceMap, SuccessorMap> 
      conn_vis(vis, position, weight, distance, pred, fwd_distance, succ);
    
    generate_rrt_star_loop(g, super_space, conn_vis, lazy_node_connector(),
                           position, node_generator_func, select_neighborhood);
    
  };
  
};




template <typename Graph, typename Topology, typename RRTStarVisitor, typename NcSelector,
          typename PositionMap, typename WeightMap,
          typename DistanceMap, typename PredecessorMap, 
          typename FwdDistanceMap = detail::infinite_double_value_prop_map, 
          typename SuccessorMap = detail::null_vertex_prop_map<Graph> >
struct rrtstar_bundle {
  typedef Graph             graph_type;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_type;
  typedef Topology          topology_type;
  typedef RRTStarVisitor    visitor_type;
  typedef NcSelector        nc_selector_type;
  typedef PositionMap       position_map_type;
  typedef WeightMap         weight_map_type;
  typedef DistanceMap       distance_map_type;
  typedef PredecessorMap    predecessor_map_type;
  typedef FwdDistanceMap    fwd_distance_map_type;
  typedef SuccessorMap      successor_map_type;
  
  graph_type* m_g;
  vertex_type m_start_vertex;
  vertex_type m_goal_vertex;
  const topology_type* m_super_space;
  visitor_type m_vis;
  nc_selector_type m_select_neighborhood;
  position_map_type m_position;
  weight_map_type m_weight;
  distance_map_type m_distance;
  predecessor_map_type m_predecessor;
  fwd_distance_map_type m_fwd_distance;
  successor_map_type m_successor;
  
  rrtstar_bundle(
    Graph &g, vertex_type start_vertex, const Topology& super_space, 
    RRTStarVisitor vis, NcSelector select_neighborhood,
    PositionMap position, WeightMap weight,
    DistanceMap distance, PredecessorMap predecessor,
    FwdDistanceMap fwd_distance = FwdDistanceMap(), SuccessorMap successor = SuccessorMap()) :
    m_g(&g), m_start_vertex(start_vertex), m_goal_vertex(boost::graph_traits<Graph>::null_vertex()), 
    m_super_space(&super_space), m_vis(vis), m_select_neighborhood(select_neighborhood),
    m_position(position), m_weight(weight), 
    m_distance(distance), m_predecessor(predecessor), 
    m_fwd_distance(fwd_distance), m_successor(successor) { };
  
  rrtstar_bundle(
    Graph &g, vertex_type start_vertex, vertex_type goal_vertex, const Topology& super_space, 
    RRTStarVisitor vis, NcSelector select_neighborhood,
    PositionMap position, WeightMap weight,
    DistanceMap distance, PredecessorMap predecessor,
    FwdDistanceMap fwd_distance = FwdDistanceMap(), SuccessorMap successor = SuccessorMap()) :
    m_g(&g), m_start_vertex(start_vertex), m_goal_vertex(goal_vertex), 
    m_super_space(&super_space), m_vis(vis), m_select_neighborhood(select_neighborhood),
    m_position(position), m_weight(weight), 
    m_distance(distance), m_predecessor(predecessor), 
    m_fwd_distance(fwd_distance), m_successor(successor) { };
  
  
};





/**
  * This function template creates a bundle of parameters to be fed to any of the
  * RRT* algorithms. This is mainly to simply the interface and the code of all these 
  * different variants of the RRT* algorithm.
  * \tparam Graph The graph type that can store the generated roadmap, should model 
  *         BidirectionalGraphConcept and MutableGraphConcept.
  * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
  * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
  * \tparam RRTStarVisitor The type of the RRT* visitor to be used, should model the RRGVisitorConcept.
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
  * \param super_space A topology (as defined by the Boost Graph Library). This topology 
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
template <typename Graph, typename Topology, typename RRTStarVisitor, typename NcSelector,
          typename PositionMap, typename WeightMap, typename DistanceMap, typename PredecessorMap>
rrtstar_bundle<Graph, Topology, RRTStarVisitor, NcSelector, PositionMap, WeightMap,
               DistanceMap, PredecessorMap> 
  make_rrtstar_bundle(Graph &g, 
    typename boost::graph_traits<Graph>::vertex_descriptor start_vertex, 
    const Topology& super_space, RRTStarVisitor vis, NcSelector select_neighborhood, 
    PositionMap position, WeightMap weight, DistanceMap distance, PredecessorMap predecessor) {
  
  BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((RRGVisitorConcept<RRTStarVisitor,Graph,Topology>));
  
  return rrtstar_bundle<Graph, Topology, RRTStarVisitor, NcSelector, PositionMap, WeightMap,
                        DistanceMap, PredecessorMap>(
    g, start_vertex, super_space, vis, select_neighborhood, position, weight, distance, predecessor);
};

template <typename Graph, typename Topology, typename RRTStarVisitor, typename NcSelector,
          typename PositionMap, typename WeightMap, typename DistanceMap, typename PredecessorMap>
rrtstar_bundle<Graph, Topology, RRTStarVisitor, NcSelector, PositionMap, WeightMap,
               DistanceMap, PredecessorMap> 
  make_rrtstar_bundle(Graph &g, 
    typename boost::graph_traits<Graph>::vertex_descriptor start_vertex, 
    typename boost::graph_traits<Graph>::vertex_descriptor goal_vertex, 
    const Topology& super_space, RRTStarVisitor vis, NcSelector select_neighborhood, 
    PositionMap position, WeightMap weight, DistanceMap distance, PredecessorMap predecessor) {
  
  BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((RRGVisitorConcept<RRTStarVisitor,Graph,Topology>));
  
  return rrtstar_bundle<Graph, Topology, RRTStarVisitor, NcSelector, PositionMap, WeightMap,
                        DistanceMap, PredecessorMap>(
    g, start_vertex, goal_vertex, super_space, vis, select_neighborhood, position, weight, distance, predecessor);
};


/**
  * This function template creates a bundle of parameters to be fed to any of the
  * RRT* algorithms. This is mainly to simply the interface and the code of all these 
  * different variants of the RRT* algorithm.
  * \tparam Graph The graph type that can store the generated roadmap, should model 
  *         BidirectionalGraphConcept and MutableGraphConcept.
  * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
  * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
  * \tparam RRTStarVisitor The type of the RRT* visitor to be used, should model the RRGVisitorConcept.
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
  * \tparam FwdDistanceMap This property-map type is used to store the current best distance of each vertex-property object
  *         to the goal.
  * \tparam SuccessorMap This property-map type is used to store the resulting path by connecting 
  *         vertex-property object together with its optimal successor.
  * 
  * \param g A mutable graph that should initially store the starting 
  *        vertex (if not it will be randomly generated) and will store 
  *        the generated graph once the algorithm has finished.
  * \param start_vertex The starting point of the algorithm, on the graph.
  * \param goal_vertex The goal point of the algorithm, on the graph.
  * \param super_space A topology (as defined by the Boost Graph Library). This topology 
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

template <typename Graph, typename Topology, typename RRTStarVisitor, typename NcSelector,
          typename PositionMap, typename WeightMap,
          typename DistanceMap, typename PredecessorMap, typename FwdDistanceMap, typename SuccessorMap >
rrtstar_bundle<Graph, Topology, RRTStarVisitor, NcSelector, PositionMap, WeightMap,
               DistanceMap, PredecessorMap, FwdDistanceMap, SuccessorMap>
  make_rrtstar_bundle(Graph &g, 
    typename boost::graph_traits<Graph>::vertex_descriptor start_vertex, 
    typename boost::graph_traits<Graph>::vertex_descriptor goal_vertex, 
    const Topology& super_space, RRTStarVisitor vis, NcSelector select_neighborhood, 
    PositionMap position, WeightMap weight, 
    DistanceMap distance, PredecessorMap predecessor, FwdDistanceMap fwd_distance, SuccessorMap successor) {
  
  BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((RRGBidirVisitorConcept<RRTStarVisitor,Graph,Topology>));
  
  return rrtstar_bundle<Graph, Topology, RRTStarVisitor, NcSelector, PositionMap, WeightMap,
                        DistanceMap, PredecessorMap, FwdDistanceMap, SuccessorMap>(
    g, start_vertex, goal_vertex, super_space, vis, select_neighborhood,
    position, weight, distance, predecessor, fwd_distance, successor);
};






/**
  * This function template is the RRT* algorithm (refer to rrt_star.hpp dox).
  * \tparam Graph A mutable graph type that will represent the generated tree, should model boost::VertexListGraphConcept and boost::MutableGraphConcept
  * \tparam Topology A topology type that will represent the space in which the configurations (or positions) exist, should model BGL's Topology concept
  * \tparam RRGVisitor An RRT* visitor type that implements the customizations to this RRT* algorithm, should model the RRGVisitorConcept.
  * \tparam PositionMap A property-map type that can store the configurations (or positions) of the vertices.
  * \tparam DistanceMap This property-map type is used to store the estimated cost-to-go of each vertex to the start (or goal).
  * \tparam PredecessorMap This property-map type is used to store the predecessor of each vertex.
  * \tparam WeightMap This property-map type is used to store the weights of the edges of the graph (cost of travel along an edge).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \tparam NcSelector A functor type which can perform a neighborhood search of a point to a graph in the topology (see topological_search.hpp).
  * \param g A mutable graph that should initially store the starting and goal 
  *        vertex and will store the generated graph once the algorithm has finished.
  * \param super_space A topology (as defined by the Boost Graph Library). Note 
  *        that it should represent the entire configuration space (not collision-free space).
  * \param vis A RRT* visitor implementing the RRGVisitorConcept. This is the 
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
template <typename Graph, typename Topology, typename RRGVisitor,
          typename PositionMap, typename DistanceMap, typename PredecessorMap,
          typename WeightMap, typename RandomSampler, typename NcSelector>
void generate_rrt_star(
    Graph& g, const Topology& super_space, RRGVisitor vis,
    PositionMap position, DistanceMap distance, PredecessorMap pred,
    WeightMap weight, RandomSampler get_sample, NcSelector select_neighborhood) {
  BOOST_CONCEPT_ASSERT((RRGVisitorConcept<RRGVisitor,Graph,Topology>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::RandomSamplerConcept<RandomSampler,Topology>));
  
  if(num_vertices(g) == 0) 
    throw optim::infeasible_problem("Cannot solve a RRT* problem without start position!");
  
  detail::rrt_conn_visitor<Graph, RRGVisitor, PositionMap, WeightMap, DistanceMap, PredecessorMap> 
    conn_vis(vis, position, weight, distance, pred);
  
  detail::generate_rrt_star_loop(g, super_space, conn_vis,
    lazy_node_connector(), position,
    rrg_node_generator<Topology, RandomSampler, NcSelector>(&super_space, get_sample, select_neighborhood),
    select_neighborhood);
  
};


/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the RRT* algorithm, without initialization of the existing graph.
  * \tparam RRTStarBundle A RRT* bundle type (see make_rrtstar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a RRT* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
  */
template <typename RRTStarBundle, typename RandomSampler>
void generate_rrt_star(const RRTStarBundle& bdl, RandomSampler get_sample) {
  
  put(bdl.m_distance, (*(bdl.m_g))[bdl.m_start_vertex], 0.0);
  put(bdl.m_predecessor, (*(bdl.m_g))[bdl.m_start_vertex], bdl.m_start_vertex);
  
  generate_rrt_star(
    *(bdl.m_g), *(bdl.m_super_space), bdl.m_vis, 
    bdl.m_position, bdl.m_distance, bdl.m_predecessor, 
    bdl.m_weight, get_sample, bdl.m_select_neighborhood);
};



/**
  * This function template is the Bi-directional RRT* algorithm (refer to rrt_star_bidir.hpp dox).
  * \tparam Graph A mutable graph type that will represent the generated tree, should model boost::VertexListGraphConcept and boost::MutableGraphConcept
  * \tparam Topology A topology type that will represent the space in which the configurations (or positions) exist, should model BGL's Topology concept
  * \tparam RRGVisitor An RRT* visitor type that implements the customizations to this RRT* algorithm, should model the RRGVisitorConcept.
  * \tparam PositionMap A property-map type that can store the configurations (or positions) of the vertices.
  * \tparam DistanceMap This property-map type is used to store the estimated cost-to-go of each vertex from the start.
  * \tparam PredecessorMap This property-map type is used to store the predecessor of each vertex.
  * \tparam FwdDistanceMap This property-map type is used to store the estimated cost-to-go of each vertex to the goal.
  * \tparam SuccessorMap This property-map type is used to store the successor of each vertex.
  * \tparam WeightMap This property-map type is used to store the weights of the edges of the graph (cost of travel along an edge).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \tparam NcSelector A functor type which can perform a neighborhood search of a point to a graph in the topology (see topological_search.hpp).
  * \param g A mutable graph that should initially store the starting and goal 
  *        vertex and will store the generated graph once the algorithm has finished.
  * \param super_space A topology (as defined by the Boost Graph Library). Note 
  *        that it should represent the entire configuration space (not collision-free space).
  * \param vis A RRT* visitor implementing the RRGVisitorConcept. This is the 
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
template <typename Graph, typename Topology, typename RRGVisitor, typename PositionMap,
          typename DistanceMap, typename PredecessorMap,
          typename FwdDistanceMap, typename SuccessorMap,
          typename WeightMap, typename RandomSampler, typename NcSelector>
void generate_rrt_star_bidir(
    Graph& g, const Topology& super_space, RRGVisitor vis, PositionMap position,
    DistanceMap distance, PredecessorMap pred, FwdDistanceMap fwd_distance, SuccessorMap succ,
    WeightMap weight, RandomSampler get_sample, NcSelector select_neighborhood) {
  BOOST_CONCEPT_ASSERT((RRGBidirVisitorConcept<RRGVisitor,Graph,Topology>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::RandomSamplerConcept<RandomSampler,Topology>));
  
  if(num_vertices(g) == 0) 
    throw optim::infeasible_problem("Cannot solve a bi-directional RRT* problem without start and goal positions!");
  
  detail::rrt_conn_visitor<Graph, RRGVisitor, PositionMap, WeightMap, 
                           DistanceMap, PredecessorMap, FwdDistanceMap, SuccessorMap> 
    conn_vis(vis, position, weight, distance, pred, fwd_distance, succ);
  
  detail::generate_rrt_star_bidir_loop(g, super_space, conn_vis,
    lazy_node_connector(), position,
    rrg_bidir_generator<Topology, RandomSampler, NcSelector, PredecessorMap, SuccessorMap>(&super_space, get_sample, select_neighborhood, pred, succ),
    select_neighborhood);
  
};


/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional RRT* algorithm, without initialization of the existing graph.
  * \tparam RRTStarBundle A RRT* bundle type (see make_rrtstar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a RRT* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
  */
template <typename RRTStarBundle, typename RandomSampler>
void generate_rrt_star_bidir(const RRTStarBundle& bdl, RandomSampler get_sample) {
  typedef typename RRTStarBundle::graph_type Graph;
  
  put(bdl.m_distance, (*(bdl.m_g))[bdl.m_start_vertex], 0.0);
  put(bdl.m_predecessor, (*(bdl.m_g))[bdl.m_start_vertex], bdl.m_start_vertex);
  if( bdl.m_goal_vertex != boost::graph_traits<Graph>::null_vertex() ) {
    put(bdl.m_fwd_distance, (*(bdl.m_g))[bdl.m_goal_vertex], 0.0);
    put(bdl.m_successor, (*(bdl.m_g))[bdl.m_goal_vertex], bdl.m_goal_vertex);
  };
  
  generate_rrt_star_bidir(
    *(bdl.m_g), *(bdl.m_super_space), bdl.m_vis, 
    bdl.m_position, bdl.m_distance, bdl.m_predecessor, 
    bdl.m_fwd_distance, bdl.m_successor, 
    bdl.m_weight, get_sample, bdl.m_select_neighborhood);
};



/**
  * This function template is the RRT* algorithm (refer to rrt_star.hpp dox).
  * This function uses a branch-and-bound heuristic to limit the number of nodes.
  * \tparam Graph A mutable graph type that will represent the generated tree, should model boost::VertexListGraphConcept and boost::MutableGraphConcept
  * \tparam Topology A topology type that will represent the space in which the configurations (or positions) exist, should model BGL's Topology concept
  * \tparam RRGVisitor An RRT* visitor type that implements the customizations to this RRT* algorithm, should model the RRGVisitorConcept.
  * \tparam PositionMap A property-map type that can store the configurations (or positions) of the vertices.
  * \tparam DistanceMap This property-map type is used to store the estimated cost-to-go of each vertex to the start (or goal).
  * \tparam PredecessorMap This property-map type is used to store the predecessor of each vertex.
  * \tparam WeightMap This property-map type is used to store the weights of the edges of the graph (cost of travel along an edge).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \tparam NcSelector A functor type which can perform a neighborhood search of a point to a graph in the topology (see topological_search.hpp).
  * \param g A mutable graph that should initially store the starting and goal 
  *        vertex and will store the generated graph once the algorithm has finished.
  * \param start_vertex The vertex from which the motion-graph is grown.
  * \param goal_vertex The vertex which we want to connect to the motion-graph.
  * \param super_space A topology (as defined by the Boost Graph Library). Note 
  *        that it should represent the entire configuration space (not collision-free space).
  * \param vis A RRT* visitor implementing the RRGVisitorConcept. This is the 
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
template <typename Graph, typename Topology, typename RRGVisitor, typename PositionMap,
          typename DistanceMap, typename PredecessorMap,
          typename WeightMap, typename RandomSampler, typename NcSelector>
void generate_bnb_rrt_star(Graph& g,
    typename boost::graph_traits<Graph>::vertex_descriptor start_vertex,
    typename boost::graph_traits<Graph>::vertex_descriptor goal_vertex,
    const Topology& super_space, RRGVisitor vis, PositionMap position,
    DistanceMap distance, PredecessorMap pred, WeightMap weight,
    RandomSampler get_sample, NcSelector select_neighborhood) {
  BOOST_CONCEPT_ASSERT((RRGVisitorConcept<RRGVisitor,Graph,Topology>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::RandomSamplerConcept<RandomSampler,Topology>));
  
  if( (num_vertices(g) == 0) ||
      ( start_vertex == boost::graph_traits<Graph>::null_vertex() ) ||
      ( goal_vertex  == boost::graph_traits<Graph>::null_vertex() ) ) {
    generate_rrt_star(g,super_space,vis,position,distance,pred,weight,get_sample,select_neighborhood);
    return;
  };
  
  detail::rrt_conn_visitor<Graph, RRGVisitor, PositionMap, WeightMap, DistanceMap, PredecessorMap> 
    conn_vis(vis, position, weight, distance, pred);
  
  detail::generate_rrt_star_loop(g, super_space, conn_vis,
    branch_and_bound_connector<Graph>(g, start_vertex, goal_vertex), position,
    rrg_node_generator<Topology, RandomSampler, NcSelector>(&super_space, get_sample, select_neighborhood),
    select_neighborhood);
  
};

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the RRT* algorithm, without initialization of the existing graph.
  * This function uses a branch-and-bound heuristic to limit the number of nodes.
  * \tparam RRTStarBundle A RRT* bundle type (see make_rrtstar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a RRT* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
  */
template <typename RRTStarBundle, typename RandomSampler>
void generate_bnb_rrt_star(const RRTStarBundle& bdl, RandomSampler get_sample) {
  
  put(bdl.m_distance, (*(bdl.m_g))[bdl.m_start_vertex], 0.0);
  put(bdl.m_predecessor, (*(bdl.m_g))[bdl.m_start_vertex], bdl.m_start_vertex);
  
  generate_bnb_rrt_star(
    *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, 
    *(bdl.m_super_space), bdl.m_vis, bdl.m_position, 
    bdl.m_distance, bdl.m_predecessor, bdl.m_weight, get_sample, bdl.m_select_neighborhood);
};


/**
  * This function template is the Bi-directional RRT* algorithm (refer to rrt_star_bidir.hpp dox).
  * This function uses a branch-and-bound heuristic to limit the number of nodes.
  * \tparam Graph A mutable graph type that will represent the generated tree, should model boost::VertexListGraphConcept and boost::MutableGraphConcept
  * \tparam Topology A topology type that will represent the space in which the configurations (or positions) exist, should model BGL's Topology concept
  * \tparam RRGVisitor An RRT* visitor type that implements the customizations to this RRT* algorithm, should model the RRGVisitorConcept.
  * \tparam PositionMap A property-map type that can store the configurations (or positions) of the vertices.
  * \tparam DistanceMap This property-map type is used to store the estimated cost-to-go of each vertex from the start.
  * \tparam PredecessorMap This property-map type is used to store the predecessor of each vertex.
  * \tparam FwdDistanceMap This property-map type is used to store the estimated cost-to-go of each vertex to the goal.
  * \tparam SuccessorMap This property-map type is used to store the successor of each vertex.
  * \tparam WeightMap This property-map type is used to store the weights of the edges of the graph (cost of travel along an edge).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \tparam NcSelector A functor type which can perform a neighborhood search of a point to a graph in the topology (see topological_search.hpp).
  * \param g A mutable graph that should initially store the starting and goal 
  *        vertex and will store the generated graph once the algorithm has finished.
  * \param start_vertex The vertex from which the motion-graph is grown.
  * \param goal_vertex The vertex which we want to connect to the motion-graph.
  * \param super_space A topology (as defined by the Boost Graph Library). Note 
  *        that it should represent the entire configuration space (not collision-free space).
  * \param vis A RRT* visitor implementing the RRGVisitorConcept. This is the 
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
template <typename Graph, typename Topology, typename RRGVisitor, typename PositionMap,
          typename DistanceMap, typename PredecessorMap,
          typename FwdDistanceMap, typename SuccessorMap,
          typename WeightMap, typename RandomSampler, typename NcSelector>
void generate_bnb_rrt_star_bidir( Graph& g,
    typename boost::graph_traits<Graph>::vertex_descriptor start_vertex,
    typename boost::graph_traits<Graph>::vertex_descriptor goal_vertex,
    const Topology& super_space, RRGVisitor vis, PositionMap position,
    DistanceMap distance, PredecessorMap pred, FwdDistanceMap fwd_distance, SuccessorMap succ,
    WeightMap weight, RandomSampler get_sample, NcSelector select_neighborhood) {
  BOOST_CONCEPT_ASSERT((RRGBidirVisitorConcept<RRGVisitor,Graph,Topology>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::RandomSamplerConcept<RandomSampler,Topology>));
  
  if( (num_vertices(g) == 0) ||
      ( start_vertex == boost::graph_traits<Graph>::null_vertex() ) ||
      ( goal_vertex  == boost::graph_traits<Graph>::null_vertex() ) ) {
    generate_rrt_star_bidir(g,super_space,vis,position,distance,pred,fwd_distance,succ,weight,get_sample,select_neighborhood);
    return;
  };
  
  detail::rrt_conn_visitor<Graph, RRGVisitor, PositionMap, WeightMap, 
                           DistanceMap, PredecessorMap, FwdDistanceMap, SuccessorMap> 
    conn_vis(vis, position, weight, distance, pred, fwd_distance, succ);
  
  detail::generate_rrt_star_bidir_loop(g, super_space, conn_vis,
    branch_and_bound_connector<Graph>(g, start_vertex, goal_vertex), position,
    rrg_bidir_generator<Topology, RandomSampler, NcSelector, PredecessorMap, SuccessorMap>(&super_space, get_sample, select_neighborhood, pred, succ),
    select_neighborhood);
  
};

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the RRT* algorithm, without initialization of the existing graph.
  * This function uses a branch-and-bound heuristic to limit the number of nodes.
  * \tparam RRTStarBundle A RRT* bundle type (see make_rrtstar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a RRT* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the free-space (obstacle-free sub-set of the topology).
  */
template <typename RRTStarBundle, typename RandomSampler>
void generate_bnb_rrt_star_bidir(const RRTStarBundle& bdl, RandomSampler get_sample) {
  typedef typename RRTStarBundle::graph_type Graph;
  
  put(bdl.m_distance, (*(bdl.m_g))[bdl.m_start_vertex], 0.0);
  put(bdl.m_predecessor, (*(bdl.m_g))[bdl.m_start_vertex], bdl.m_start_vertex);
  if( bdl.m_goal_vertex != boost::graph_traits<Graph>::null_vertex() ) {
    put(bdl.m_fwd_distance, (*(bdl.m_g))[bdl.m_goal_vertex], 0.0);
    put(bdl.m_successor, (*(bdl.m_g))[bdl.m_goal_vertex], bdl.m_goal_vertex);
  };
  
  generate_bnb_rrt_star_bidir(
    *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, 
    *(bdl.m_super_space), bdl.m_vis, bdl.m_position, 
    bdl.m_distance, bdl.m_predecessor, bdl.m_fwd_distance, bdl.m_successor, 
    bdl.m_weight, get_sample, bdl.m_select_neighborhood);
};



};

};


#endif


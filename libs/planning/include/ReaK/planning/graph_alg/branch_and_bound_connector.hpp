/**
 * \file branch_and_bound_connector.hpp
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

#ifndef REAK_BRANCH_AND_BOUND_CONNECTOR_HPP
#define REAK_BRANCH_AND_BOUND_CONNECTOR_HPP

#include <functional>
#include <boost/utility/enable_if.hpp>

#include <ReaK/topologies/spaces/metric_space_concept.hpp>

#include "lazy_connector.hpp"

#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/detail/d_ary_heap.hpp>

// BGL-Extra includes:
#include <boost/graph/more_property_tags.hpp>
#include <boost/graph/more_property_maps.hpp>

#include <stack>


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Graph */
namespace graph {
  
  
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
  * conn_vis.vertex_to_be_removed(u, g);  This function is called just before a vertex is removed from the motion-graph (i.e., pruned by the branch-and-bound heuristic).
  * 
  * \tparam Visitor The visitor class to be tested for modeling an AD* visitor concept.
  * \tparam Graph The graph type on which the visitor should be able to act.
  * \tparam Topology The topology type on which the visitor class is required to work with.
  */
template <typename ConnectorVisitor, typename Graph, typename Topology>
struct BNBConnectorVisitorConcept : MotionGraphConnectorVisitorConcept<ConnectorVisitor, Graph, Topology> {
  ConnectorVisitor conn_vis;
  Graph g;
  typename boost::graph_traits<Graph>::vertex_descriptor u;
  
  BOOST_CONCEPT_USAGE(BNBConnectorVisitorConcept)
  {
    conn_vis.vertex_to_be_removed(u, g);
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
struct branch_and_bound_connector {
  
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename Graph::edge_bundled EdgeProp;
  typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
  typedef typename boost::graph_traits<Graph>::in_edge_iterator InEdgeIter;
  
  typedef boost::vector_property_map<std::size_t> IndexInHeapMap;
  typedef boost::vector_property_map<double> KeyMap;
  typedef std::greater<double> KeyCompareType;  // <---- this is a max-heap.
  typedef boost::d_ary_heap_indirect<Vertex, 4, IndexInHeapMap, KeyMap, KeyCompareType> MutableQueue;
  
  Vertex start_vertex;
  Vertex goal_vertex;
  
  mutable IndexInHeapMap index_in_heap;
  mutable KeyMap key;
  mutable MutableQueue Q; //priority queue holding the OPEN set.
  
  
  
  branch_and_bound_connector(Graph& g,
                             Vertex aStartVertex, 
                             Vertex aGoalVertex) : 
                             start_vertex(aStartVertex), 
                             goal_vertex(aGoalVertex),
                             index_in_heap(), key(),
                             Q(key, index_in_heap, KeyCompareType()) { 
    
    typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
    for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui)
      put(index_in_heap,*ui, static_cast<std::size_t>(-1));
    
  };
  
  
  template <typename Topology, typename ConnectorVisitor,
            typename PositionMap, typename DistanceMap, typename PredecessorMap,
            typename WeightMap>
  void update_successors(
      Vertex v, Graph& g,
      const Topology& super_space, const ConnectorVisitor& conn_vis,
      PositionMap position, DistanceMap distance, PredecessorMap predecessor,
      WeightMap weight) const {
    
    // need to update all the children of the v node:
    std::stack< Vertex > incons;
    incons.push(v);
    while(!incons.empty()) {
      Vertex s = incons.top(); incons.pop();
      OutEdgeIter eo, eo_end;
      for(boost::tie(eo,eo_end) = out_edges(s, g); eo != eo_end; ++eo) {
        Vertex t = target(*eo, g);
        if(t == s)
          t = source(*eo, g);
        if(s != get(predecessor, g[t]))
          continue;
        put(distance, g[t], get(distance, g[s]) + get(weight, g[*eo]));
        
        conn_vis.affected_vertex(t,g);  // affected by changed distance value.
        
        put(key, t, get(distance, g[t]) + get(ReaK::pp::distance_metric, super_space)(get(position, g[t]), get(position, g[goal_vertex]), super_space) );
        Q.push_or_update(t);
        
        incons.push(t);
      };
    };
    
    if( get(predecessor, g[goal_vertex]) != boost::graph_traits<Graph>::null_vertex() ) {
      // prune all the worst nodes:
      double top_value = get(key, Q.top());
      while(top_value > get(distance, g[goal_vertex])) {
        conn_vis.vertex_to_be_removed(Q.top(), g);
        clear_vertex(Q.top(), g);
        remove_vertex(Q.top(), g);
        Q.pop();
        top_value = get(key, Q.top());
      };
    };
  };
  
  template <typename Topology, typename ConnectorVisitor,
            typename PositionMap, typename FwdDistanceMap, typename SuccessorMap,
            typename WeightMap>
  void update_predecessors(
      Vertex v, Graph& g,
      const Topology& super_space, const ConnectorVisitor& conn_vis,
      PositionMap position, FwdDistanceMap fwd_distance, SuccessorMap successor,
      WeightMap weight) const {
    
    // need to update all the children of the v node:
    std::stack< Vertex > incons;
    incons.push(v);
    while(!incons.empty()) {
      Vertex t = incons.top(); incons.pop();
      InEdgeIter ei, ei_end;
      for(boost::tie(ei,ei_end) = in_edges(t, g); ei != ei_end; ++ei) {
        Vertex s = source(*ei, g);
        if(t == s)
          s = target(*ei, g);
        if(t != get(successor, g[s]))
          continue;
        put(fwd_distance, g[s], get(fwd_distance, g[t]) + get(weight, g[*ei]));
        
        conn_vis.affected_vertex(s,g);  // affected by changed distance value.
        
        put(key, s, get(fwd_distance, g[s]) + get(ReaK::pp::distance_metric, super_space)(get(position, g[start_vertex]), get(position, g[s]), super_space) );
        Q.push_or_update(s);
        
        incons.push(s);
      };
    };
    
    if( get(successor, g[start_vertex]) != boost::graph_traits<Graph>::null_vertex() ) {
      // prune all the worst nodes:
      double top_value = get(key, Q.top());
      while(top_value > get(fwd_distance, g[start_vertex])) {
        conn_vis.vertex_to_be_removed(Q.top(), g);
        clear_vertex(Q.top(), g);
        remove_vertex(Q.top(), g);
        Q.pop();
        top_value = get(key, Q.top());
      };
    };
  };
  
  
  
  
  /**
   * This call operator takes a position value, the predecessor from which the new position was obtained,
   * the travel-record (as an edge property) that can do the travel from the predecessor to the new position,
   * and the other objects needed for motion planning, and it creates a new vertex for the new position and 
   * connects that new vertex to the motion-graph using a lazy and pruned strategy.
   * \note This version applies to a undirected graph (and undirected / symmetric distance metric).
   * 
   * \tparam Graph2 The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam SBAStarVisitor The type of the node-connector visitor to be used, should model the BNBConnectorVisitorConcept.
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
  template <typename Graph2, typename Topology, typename ConnectorVisitor,
            typename PositionMap, typename DistanceMap, typename PredecessorMap,
            typename WeightMap, typename NcSelector>
  typename boost::enable_if< boost::is_undirected_graph<Graph2> >::type operator()(
      const typename boost::property_traits<PositionMap>::value_type& p, 
      Vertex& x_near, EdgeProp& eprop, Graph2& g,
      const Topology& super_space, const ConnectorVisitor& conn_vis,
      PositionMap position, DistanceMap distance, PredecessorMap predecessor,
      WeightMap weight, NcSelector select_neighborhood) const {
    
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((BNBConnectorVisitorConcept<ConnectorVisitor,Graph2,Topology>));
    
    double dist_from_start = get(ReaK::pp::distance_metric, super_space)(get(position, g[start_vertex]), p, super_space);
    double dist_to_goal = get(ReaK::pp::distance_metric, super_space)(p, get(position, g[goal_vertex]), super_space);
    
    if( ( get(predecessor, g[goal_vertex]) != boost::graph_traits<Graph>::null_vertex() ) && 
        ( dist_from_start + dist_to_goal > get(distance, g[goal_vertex]) ) )
      return;
    
    std::vector<Vertex> Nc;
    select_neighborhood(p, std::back_inserter(Nc), g, super_space, boost::bundle_prop_to_vertex_prop(position, g)); 
    
    Vertex v = conn_vis.create_vertex(p, g);
    put(index_in_heap,v, static_cast<std::size_t>(-1));
    
    if( x_near != boost::graph_traits<Graph>::null_vertex() ) {
      conn_vis.travel_explored(x_near, v, g);
      conn_vis.travel_succeeded(x_near, v, g);
      conn_vis.affected_vertex(x_near, g);
    };
    
    lazy_node_connector::connect_best_predecessor(v, x_near, eprop, g, super_space, conn_vis, position, distance, predecessor, weight, Nc);
    pruned_node_connector::create_pred_edge(v, x_near, eprop, g, conn_vis, distance, predecessor, weight);
    
    if( ( get(predecessor, g[goal_vertex]) != boost::graph_traits<Graph>::null_vertex() ) && 
        ( get(distance, g[v]) + dist_to_goal > get(distance, g[goal_vertex]) ) ) {
      conn_vis.vertex_to_be_removed(v, g);
      clear_vertex(v, g);
      remove_vertex(v, g);
      return;
    };
    
    put(key, v, get(distance, g[v]) + dist_to_goal );
    Q.push(v);
    
    lazy_node_connector::connect_successors(v, x_near, g, super_space, conn_vis, position, distance, predecessor, weight, Nc);
    update_successors(v, g, super_space, conn_vis, position, distance, predecessor, weight);
    
  };
  
  
  /**
   * This call operator takes a position value, the predecessor from which the new position was obtained,
   * the travel-record (as an edge property) that can do the travel from the predecessor to the new position,
   * and the other objects needed for motion planning, and it creates a new vertex for the new position and 
   * connects that new vertex to the motion-graph using a lazy and pruned strategy.
   * \note This version applies to a directed graph (and directed / asymmetric distance metric).
   * 
   * \tparam Graph2 The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam SBAStarVisitor The type of the node-connector visitor to be used, should model the BNBConnectorVisitorConcept.
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
  template <typename Graph2, typename Topology, typename ConnectorVisitor,
            typename PositionMap, typename DistanceMap, typename PredecessorMap,
            typename WeightMap, typename NcSelector>
  typename boost::enable_if< boost::is_directed_graph<Graph2> >::type operator()(
      const typename boost::property_traits<PositionMap>::value_type& p, 
      Vertex& x_near, EdgeProp& eprop, Graph2& g,
      const Topology& super_space, const ConnectorVisitor& conn_vis,
      PositionMap position, DistanceMap distance, PredecessorMap predecessor,
      WeightMap weight, NcSelector select_neighborhood) const {
    
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((BNBConnectorVisitorConcept<ConnectorVisitor,Graph2,Topology>));
    
    double dist_from_start = get(ReaK::pp::distance_metric, super_space)(get(position, g[start_vertex]), p, super_space);
    double dist_to_goal = get(ReaK::pp::distance_metric, super_space)(p, get(position, g[goal_vertex]), super_space);
    
    if( ( get(predecessor, g[goal_vertex]) != boost::graph_traits<Graph>::null_vertex() ) && 
        ( dist_from_start + dist_to_goal > get(distance, g[goal_vertex]) ) )
      return;
    
    std::vector<Vertex> Pred, Succ;
    select_neighborhood(p, std::back_inserter(Pred), std::back_inserter(Succ), g, super_space, boost::bundle_prop_to_vertex_prop(position, g)); 
    
    Vertex v = conn_vis.create_vertex(p, g);
    put(index_in_heap,v, static_cast<std::size_t>(-1));
    
    if( x_near != boost::graph_traits<Graph>::null_vertex() ) {
      conn_vis.travel_explored(x_near, v, g);
      conn_vis.travel_succeeded(x_near, v, g);
      conn_vis.affected_vertex(x_near, g);
    };
    
    lazy_node_connector::connect_best_predecessor(v, x_near, eprop, g, super_space, conn_vis, position, distance, predecessor, weight, Pred);
    pruned_node_connector::create_pred_edge(v, x_near, eprop, g, conn_vis, distance, predecessor, weight);
    
    if( ( get(predecessor, g[goal_vertex]) != boost::graph_traits<Graph>::null_vertex() ) && 
        ( get(distance, g[v]) + dist_to_goal > get(distance, g[goal_vertex]) ) ) {
      conn_vis.vertex_to_be_removed(v, g);
      clear_vertex(v, g);
      remove_vertex(v, g);
      return;
    };
    
    put(key, v, get(distance, g[v]) + dist_to_goal );
    Q.push(v);
    
    lazy_node_connector::connect_successors(v, x_near, g, super_space, conn_vis, position, distance, predecessor, weight, Succ);
    update_successors(v, g, super_space, conn_vis, position, distance, predecessor, weight);
    
  };
  
  
  
  
  template <typename Graph2, typename Topology, typename ConnectorVisitor,
            typename PositionMap, typename DistanceMap, typename PredecessorMap,
            typename FwdDistanceMap, typename SuccessorMap, typename WeightMap, typename NcSelector>
  typename boost::enable_if< boost::is_undirected_graph<Graph2> >::type operator()(
      const typename boost::property_traits<PositionMap>::value_type& p, 
      Vertex& x_pred, 
      EdgeProp& eprop_pred, 
      Vertex& x_succ, 
      EdgeProp& eprop_succ, 
      Graph2& g, const Topology& super_space, const ConnectorVisitor& conn_vis,
      PositionMap position, DistanceMap distance, PredecessorMap predecessor,
      FwdDistanceMap fwd_distance, SuccessorMap successor,
      WeightMap weight, NcSelector select_neighborhood) const {
    
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((MotionGraphConnectorVisitorConcept<ConnectorVisitor,Graph,Topology>));
    
    using std::back_inserter;
    
    double dist_from_start = get(ReaK::pp::distance_metric, super_space)(get(position, g[start_vertex]), p, super_space);
    double dist_to_goal = get(ReaK::pp::distance_metric, super_space)(p, get(position, g[goal_vertex]), super_space);
    
    if( ( get(predecessor, g[goal_vertex]) != boost::graph_traits<Graph>::null_vertex() ) && 
        ( dist_from_start + dist_to_goal > get(distance, g[goal_vertex]) ) )
      return;
    
    std::vector<Vertex> Nc;
    select_neighborhood(p, back_inserter(Nc), g, super_space, boost::bundle_prop_to_vertex_prop(position, g)); 
    
    Vertex v = conn_vis.create_vertex(p, g);
    put(index_in_heap,v, static_cast<std::size_t>(-1));
    
    if( x_pred != boost::graph_traits<Graph>::null_vertex() ) {
      conn_vis.travel_explored(x_pred, v, g);
      conn_vis.travel_succeeded(x_pred, v, g);
      conn_vis.affected_vertex(x_pred, g);
    };
    if( x_succ != boost::graph_traits<Graph>::null_vertex() ) {
      conn_vis.travel_explored(v, x_succ, g);
      conn_vis.travel_succeeded(v, x_succ, g);
      conn_vis.affected_vertex(x_succ, g);
    };
    
    lazy_node_connector::connect_best_predecessor(v, x_pred, eprop_pred, g, super_space, conn_vis, position, distance, predecessor, weight, Nc);
    lazy_node_connector::connect_best_successor(v, x_succ, eprop_succ, g, super_space, conn_vis, position, fwd_distance, successor, weight, Nc);
    
    if( ( x_pred == boost::graph_traits<Graph>::null_vertex() ) &&
        ( x_succ == boost::graph_traits<Graph>::null_vertex() ) ){
      conn_vis.vertex_to_be_removed(v, g);
      clear_vertex(v, g);
      remove_vertex(v, g);
      return;
    };
    
    if( x_pred != boost::graph_traits<Graph>::null_vertex() )
      pruned_node_connector::create_pred_edge(v, x_pred, eprop_pred, g, conn_vis, distance, predecessor, weight);
    if( x_succ != boost::graph_traits<Graph>::null_vertex() )
      pruned_node_connector::create_succ_edge(v, x_succ, eprop_succ, g, conn_vis, fwd_distance, successor, weight);
    
    if( ( get(predecessor, g[goal_vertex]) != boost::graph_traits<Graph>::null_vertex() ) && 
        ( get(distance, g[v]) + dist_to_goal > get(distance, g[goal_vertex]) ) ) {
      conn_vis.vertex_to_be_removed(v, g);
      clear_vertex(v, g);
      remove_vertex(v, g);
      return;
    };
    put(key, v, get(distance, g[v]) + dist_to_goal );
    Q.push(v);
    
    lazy_node_connector::connect_successors(v, x_pred, g, super_space, conn_vis, position, distance, predecessor, weight, Nc, successor);
    update_successors(v, g, super_space, conn_vis, position, distance, predecessor, weight);
    lazy_node_connector::connect_predecessors(v, x_succ, g, super_space, conn_vis, position, fwd_distance, successor, weight, Nc, predecessor);
    update_predecessors(v, g, super_space, conn_vis, position, fwd_distance, successor, weight);
  };
  
  template <typename Graph2, typename Topology, typename ConnectorVisitor,
            typename PositionMap, typename DistanceMap, typename PredecessorMap,
            typename FwdDistanceMap, typename SuccessorMap, typename WeightMap, typename NcSelector>
  typename boost::enable_if< boost::is_directed_graph<Graph2> >::type operator()(
      const typename boost::property_traits<PositionMap>::value_type& p, 
      Vertex& x_pred, 
      EdgeProp& eprop_pred, 
      Vertex& x_succ, 
      EdgeProp& eprop_succ, 
      Graph2& g, const Topology& super_space, const ConnectorVisitor& conn_vis,
      PositionMap position, DistanceMap distance, PredecessorMap predecessor,
      FwdDistanceMap fwd_distance, SuccessorMap successor,
      WeightMap weight, NcSelector select_neighborhood) const {
    
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((MotionGraphConnectorVisitorConcept<ConnectorVisitor,Graph,Topology>));
    
    using std::back_inserter;
    
    double dist_from_start = get(ReaK::pp::distance_metric, super_space)(get(position, g[start_vertex]), p, super_space);
    double dist_to_goal = get(ReaK::pp::distance_metric, super_space)(p, get(position, g[goal_vertex]), super_space);
    
    if( ( get(predecessor, g[goal_vertex]) != boost::graph_traits<Graph>::null_vertex() ) && 
        ( dist_from_start + dist_to_goal > get(distance, g[goal_vertex]) ) )
      return;
    
    std::vector<Vertex> Pred, Succ;
    select_neighborhood(p, back_inserter(Pred), back_inserter(Succ), g, super_space, boost::bundle_prop_to_vertex_prop(position, g)); 
    
    Vertex v = conn_vis.create_vertex(p, g);
    put(index_in_heap,v, static_cast<std::size_t>(-1));
    
    if( x_pred != boost::graph_traits<Graph>::null_vertex() ) {
      conn_vis.travel_explored(x_pred, v, g);
      conn_vis.travel_succeeded(x_pred, v, g);
      conn_vis.affected_vertex(x_pred, g);
    };
    if( x_succ != boost::graph_traits<Graph>::null_vertex() ) {
      conn_vis.travel_explored(v, x_succ, g);
      conn_vis.travel_succeeded(v, x_succ, g);
      conn_vis.affected_vertex(x_succ, g);
    };
    
    lazy_node_connector::connect_best_predecessor(v, x_pred, eprop_pred, g, super_space, conn_vis, position, distance, predecessor, weight, Pred);
    lazy_node_connector::connect_best_successor(v, x_succ, eprop_succ, g, super_space, conn_vis, position, fwd_distance, successor, weight, Succ);
    
    if( ( x_pred == boost::graph_traits<Graph>::null_vertex() ) &&
        ( x_succ == boost::graph_traits<Graph>::null_vertex() ) ){
      conn_vis.vertex_to_be_removed(v, g);
      clear_vertex(v, g);
      remove_vertex(v, g);
      return;
    };
    
    if( x_pred != boost::graph_traits<Graph>::null_vertex() )
      pruned_node_connector::create_pred_edge(v, x_pred, eprop_pred, g, conn_vis, distance, predecessor, weight);
    if( x_succ != boost::graph_traits<Graph>::null_vertex() )
      pruned_node_connector::create_succ_edge(v, x_succ, eprop_succ, g, conn_vis, fwd_distance, successor, weight);
    
    if( ( get(predecessor, g[goal_vertex]) != boost::graph_traits<Graph>::null_vertex() ) && 
        ( get(distance, g[v]) + dist_to_goal > get(distance, g[goal_vertex]) ) ) {
      conn_vis.vertex_to_be_removed(v, g);
      clear_vertex(v, g);
      remove_vertex(v, g);
      return;
    };
    put(key, v, get(distance, g[v]) + dist_to_goal );
    Q.push(v);
    
    lazy_node_connector::connect_successors(v, x_pred, g, super_space, conn_vis, position, distance, predecessor, weight, Succ, successor);
    update_successors(v, g, super_space, conn_vis, position, distance, predecessor, weight);
    lazy_node_connector::connect_predecessors(v, x_succ, g, super_space, conn_vis, position, fwd_distance, successor, weight, Pred, predecessor);
    update_predecessors(v, g, super_space, conn_vis, position, fwd_distance, successor, weight);
  };
  
  
  
};


};

};

#endif

















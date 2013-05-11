/**
 * \file sbastar_search.hpp
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

#ifndef REAK_SBASTAR_SEARCH_HPP
#define REAK_SBASTAR_SEARCH_HPP

#include <functional>
#include <boost/utility/enable_if.hpp>

#include "path_planning/metric_space_concept.hpp"
#include "path_planning/prob_distribution_concept.hpp"

#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/exception.hpp>
#include <boost/graph/detail/d_ary_heap.hpp>

#include "bgl_more_property_maps.hpp"
#include "bgl_more_property_tags.hpp"
#include "bgl_raw_property_graph.hpp"

#include "motion_graph_connector.hpp"

#include <stack>


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Graph */
namespace graph {

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
   * The visitor class should model the boost::CopyConstructibleConcept.
   * 
   * Valid expressions:
   * 
   * vis.initialize_vertex(u,g);  A function that gets called whenever a vertex (u) is first initialized before the search.
   * 
   * vis.discover_vertex(u,g);  A function that gets called whenever a vertex (u) is added to the OPEN set (or updated in the OPEN set).
   * 
   * vis.examine_vertex(u,g);  A function that gets called whenever a vertex (u) is taken out of the OPEN set to be examined, this is called before it gets expanded.
   * 
   * vis.edge_relaxed(e,g);  A function that gets called whenever an edge (e) has been newly declared as a better alternative than the current surrounding edges (i.e. the edge is added to the optimal path).
   * 
   * vis.publish_path(g);  A function to notify the visitor that at least one A* round has completed and its resulting path (partial or complete) can be published (the path is encoded in the predecessor property-map).
   * 
   * b = vis.keep_going();  A function to check to see whether the task is finished (return false) or needs to keep going (true).
   * 
   * vis.vertex_added(u, g);  This function is called whenever a new vertex (u) has been added to the graph (g), but not yet connected.
   * 
   * vis.edge_added(e, g);  This function is called whenever a new edge (e) has been created between the last created vertex and its neighbor in the graph (g).
   * 
   * vis.travel_explored(u, v, g);  This function is called whenever a source-destination pair of vertices have been matched for a potential edge between them.
   * 
   * vis.travel_succeeded(u, v, g);  This function is called whenever a source-destination pair of vertices led to a successful travel.
   * 
   * vis.travel_failed(u, v, g);  This function is called whenever a source-destination pair of vertices led to a failure to travel (e.g., a collision occurred).
   * 
   * tie(pt,b,ep) = vis.random_walk(u, g);  This function is called to perform the expansion of the roadmap from a given vertex (u) in the graph (g). This function returns a newly generated position value that is a candidate to be added to the graph, and with an edge-property object associated with a new edge.
   * 
   * tie(b, ep) = vis.can_be_connected(u,v,g);  This function is called to attempt to steer from vertex u to vertex v, it returns true if a local path exists and is collision-free, and it also returns the edge-property of the edge that could connect those two vertices.
   * 
   * b = vis.has_search_potential(u, g);  This function is called to check if the key associated to the vertex is deemed minimally useful by some measure (e.g., a threshold value).
   * 
   * b = vis.should_close(u, g);  This function is called to check if a vertex should be closed according to some measure (e.g., a threshold value on key or density or constriction or anything else).
   * 
   * \tparam Visitor The visitor class to be tested for modeling an AD* visitor concept.
   * \tparam Graph The graph type on which the visitor should be able to act.
   * \tparam Topology The topology type on which the visitor class is required to work with.
   */
  template <typename Visitor, typename Graph, typename Topology>
  struct SBAStarVisitorConcept {
    BOOST_CONCEPT_USAGE(SBAStarVisitorConcept)
    {
      BOOST_CONCEPT_ASSERT((boost::CopyConstructibleConcept<Visitor>));
      vis.initialize_vertex(u, g);   //whenever the vertex is first initialized.
      vis.discover_vertex(u, g);     //whenever a vertex is added to the OPEN set (or updated in OPEN).
      vis.examine_vertex(u, g);      //whenever a vertex is taken out of OPEN, before it gets "expanded".
      vis.edge_relaxed(e, g);        //whenever it is newly decided that an edge is relaxed (has improved the distance for its target)
      vis.publish_path(g);           // notify the visitor that at least one A* round has completed and its resulting path (partial or complete) can be published (the path is encoded in the predecessor property-map).
      b = vis.keep_going();          // check to see whether the task is finished (return false) or needs to keep going (true).
      vis.vertex_added(u, g); 
      vis.edge_added(e, g);
      vis.travel_explored(u, u, g);
      vis.travel_succeeded(u, u, g);
      vis.travel_failed(u, u, g);
      boost::tie(pt, b, ep) = vis.random_walk(u, g);
      boost::tie(b, ep) = vis.can_be_connected(u, u, g);
      b = vis.has_search_potential(u, g); 
      b = vis.should_close(u, g); 
    }
    Visitor vis;
    Graph g;
    typename boost::graph_traits<Graph>::vertex_descriptor u;
    typename boost::graph_traits<Graph>::edge_descriptor e;
    typename Graph::edge_bundled ep;
    typename ReaK::pp::topology_traits<Topology>::point_type pt;
    bool b;
  };
  
  /**
   * This class is simply a "null" visitor for the SBA* algorithm. It is null in the sense that it
   * will do nothing on all accounts.
   */
  template <typename Topology>
  class default_sbastar_visitor {
    public:
      typedef typename ReaK::pp::topology_traits<Topology>::point_type PointType;
      
      default_sbastar_visitor() {};
      
      template <typename Vertex, typename Graph>
      void initialize_vertex(Vertex, const Graph&) const { };
      
      template <typename Vertex, typename Graph>
      void discover_vertex(Vertex, const Graph&) const { };
      
      template <typename Vertex, typename Graph>
      void examine_vertex(Vertex, const Graph&) const { };
      
      template <typename Edge, typename Graph>
      void edge_relaxed(Edge, const Graph&) const { };
      
      template <typename Graph>
      void publish_path(const Graph&) const { };
      
      template <typename Vertex, typename Graph>
      void vertex_added(Vertex, const Graph&) const { };
      
      template <typename Vertex, typename Graph>
      void vertex_to_be_removed(Vertex, Graph&) const;
      
      template <typename Edge, typename Graph>
      void edge_added(Edge, const Graph&) const { };
      
      template <typename Vertex, typename Graph>
      void travel_explored(Vertex, Vertex, Graph&) const { };
      
      template <typename Vertex, typename Graph>
      void travel_succeeded(Vertex, Vertex, Graph&) const { };
      
      template <typename Vertex, typename Graph>
      void travel_failed(Vertex, Vertex, Graph&) const { };
      
      template <typename Vertex, typename Graph>
      boost::tuple<PointType, bool, typename Graph::edge_bundled> random_walk(Vertex, const Graph&) const { 
        typedef typename Graph::edge_bundled EdgeProp;
        typedef boost::tuple<PointType, bool, EdgeProp> ResultType;
        return ResultType(PointType(), false, EdgeProp()); 
      };
      
      template <typename Vertex, typename Graph>
      std::pair<bool, typename Graph::edge_bundled> can_be_connected(Vertex, Vertex, const Graph&) const { 
        typedef typename Graph::edge_bundled EdgeProp;
        typedef std::pair<bool, EdgeProp> ResultType;
        return ResultType(false, EdgeProp()); 
      };
      
      bool keep_going() const { return true; };
      
      template <typename Vertex, typename Graph>
      bool has_search_potential(Vertex, const Graph&) const { return true; };
      
      template <typename Vertex, typename Graph>
      bool should_close(Vertex, const Graph&) const { return false; };
  };
  
  
  
  namespace detail {
  
#if 0
    // Old version of the connector:
    
    template <typename Topology,
              typename UniformCostVisitor,
              typename SBANodeConnector,
              typename UpdatableQueue, 
              typename IndexInHeapMap,
              typename AStarHeuristicMap, 
              typename PositionMap, 
              typename WeightMap,
              typename DensityMap,
              typename ConstrictionMap, 
              typename DistanceMap,  
              typename PredecessorMap,
              typename KeyMap,
              typename NcSelector>
    struct sbastar_bfs_visitor
    {

      sbastar_bfs_visitor(const Topology& super_space, UniformCostVisitor vis, SBANodeConnector connector,
                          UpdatableQueue& Q, IndexInHeapMap index_in_heap,  
                          AStarHeuristicMap heuristic, PositionMap pos, WeightMap weight, 
                          DensityMap density, ConstrictionMap constriction, DistanceMap dist, 
                          PredecessorMap pred, KeyMap key, NcSelector select_neighborhood) : 
                          m_vis(vis), m_Q(Q), m_index_in_heap(index_in_heap), 
                          m_heuristic(heuristic), m_position(pos), m_weight(weight),
                          m_density(density), m_constriction(constriction), m_distance(dist),
                          m_predecessor(pred), m_key(key), 
                          m_super_space(super_space), m_connect_vertex(connector), m_select_neighborhood(select_neighborhood) { };
#endif
                          
    template <typename UniformCostVisitor,
              typename UpdatableQueue, 
              typename IndexInHeapMap,
              typename AStarHeuristicMap, 
              typename PositionMap, 
              typename WeightMap,
              typename DensityMap,
              typename ConstrictionMap, 
              typename DistanceMap,  
              typename PredecessorMap,
              typename KeyMap>
    struct sbastar_bfs_visitor
    {

      sbastar_bfs_visitor(UniformCostVisitor vis, UpdatableQueue& Q, IndexInHeapMap index_in_heap,  
                          AStarHeuristicMap heuristic, PositionMap pos, WeightMap weight, 
                          DensityMap density, ConstrictionMap constriction, 
                          DistanceMap dist, PredecessorMap pred, KeyMap key) : 
                          m_vis(vis), m_Q(Q), m_index_in_heap(index_in_heap), 
                          m_heuristic(heuristic), m_position(pos), m_weight(weight),
                          m_density(density), m_constriction(constriction), 
                          m_distance(dist), m_predecessor(pred), m_key(key) { };
                          
      
      typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
      
      template <class Graph>
      typename boost::graph_traits<Graph>::vertex_descriptor create_vertex(const PositionValue& p, Graph& g) const {
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename Graph::vertex_bundled VertexProp;
        
        VertexProp up;
        put(m_position, up, p);
        put(m_distance, up, std::numeric_limits<double>::infinity());
#ifdef RK_ENABLE_CXX0X_FEATURES
        Vertex u = add_vertex(std::move(up), g);
#else
        Vertex u = add_vertex(up, g);
#endif
        put(m_predecessor, g[u], u);
        m_vis.vertex_added(u,g);
        put(m_index_in_heap, u, static_cast<std::size_t>(-1));
        put(m_key, u, 0.0);
        
        return u;
      };
      
      template <typename Vertex, typename Graph>
      std::pair< bool, typename Graph::edge_bundled > can_be_connected(Vertex u, Vertex v, Graph& g) const {
        return m_vis.can_be_connected(u, v, g);
      };
      
      template <typename Vertex, typename Graph>
      boost::tuple< PositionValue, bool, typename Graph::edge_bundled > random_walk(Vertex u, Graph& g) const {
        return m_vis.random_walk(u, g);
      };
      
      template <typename Vertex, typename Graph>
      void vertex_to_be_removed(Vertex u, Graph& g) const {
        m_vis.vertex_to_be_removed(u,g);
      };
      
      template <typename Edge, typename Graph>
      void edge_added(Edge e, Graph& g) const {
        m_vis.edge_added(e,g);
      };
      
      template <typename Vertex, typename Graph>
      void travel_explored(Vertex u, Vertex v, Graph& g) const {
        m_vis.travel_explored(u, v, g);
      };
      
      template <typename Vertex, typename Graph>
      void travel_succeeded(Vertex u, Vertex v, Graph& g) const {
        m_vis.travel_succeeded(u, v, g);
      };
      
      template <typename Vertex, typename Graph>
      void travel_failed(Vertex u, Vertex v, Graph& g) const {
        m_vis.travel_failed(u, v, g);
      };
      
      template <typename Vertex, typename Graph>
      void requeue_vertex(Vertex u, Graph& g) const {
        update_key(u,g); 
        if( ! m_vis.should_close(u, g) ) {
          m_Q.push_or_update(u);
          m_vis.discover_vertex(u, g);
        };
      };
      template <typename Vertex, typename Graph>
      void affected_vertex(Vertex u, Graph& g) const { requeue_vertex(u,g); }; // same function, different name.
      
    
    // Old version of the connector:
#if 0
      template <typename Vertex, typename EdgeProp, typename Graph>
      std::pair<typename boost::graph_traits<Graph>::edge_descriptor, bool> 
          create_edge(Vertex u, Vertex v, EdgeProp& ep, Graph& g) const {
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        
        m_vis.travel_explored(u, v, g);
        m_vis.travel_succeeded(u, v, g);
#ifdef RK_ENABLE_CXX0X_FEATURES
        std::pair<Edge, bool> e_new = add_edge(u, v, std::move(ep), g);
#else
        std::pair<Edge, bool> e_new = add_edge(u, v, ep, g);
#endif
        if(e_new.second)
          m_vis.edge_added(e_new.first, g);
        return e_new;
      };
      
      template <typename Vertex, typename Graph>
      std::pair<typename boost::graph_traits<Graph>::edge_descriptor, bool>
          attempt_connecting_edge(Vertex u, Vertex v, Graph& g) const {
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename Graph::edge_bundled EdgeProp;
        
        std::pair<bool, EdgeProp> connect_result = m_vis.can_be_connected(u, v, g);
        if(connect_result.first) {
          return create_edge(u, v, connect_result.second, g);
        } else {
          m_vis.travel_explored(u, v, g);
          m_vis.travel_failed(u, v, g);
          return std::pair<Edge, bool>(Edge(), false);
        };
      };
      
      template <typename Vertex, typename Graph>
      void requeue_vertex(Vertex u, const Graph& g) const { 
        update_key(u,g); 
        if( ! m_vis.should_close(u, g) ) {
          m_Q.push_or_update(u);
          m_vis.discover_vertex(u, g);
        };
      };

      template <class Vertex, class Graph>
      void update_vertex(Vertex u, Graph& g) const {
        boost::function_requires< boost::BidirectionalGraphConcept<Graph> >();
        typedef typename boost::graph_traits<Graph>::in_edge_iterator InEdgeIter;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        
        double g_u = get(m_distance, u); 
        if(g_u != 0.0) {
          g_u = std::numeric_limits<double>::infinity(); 
          Vertex pred_u = get(m_predecessor, u);
          Edge pred_e; 
          InEdgeIter ei, ei_end;
          for(boost::tie(ei,ei_end) = in_edges(u,g); ei != ei_end; ++ei) {
            double g_tmp = get(m_weight, g[*ei]) + get(m_distance, source(*ei,g)); 
            if(g_tmp < g_u) {
              g_u = g_tmp; 
              put(m_distance, u, g_u);
              put(m_predecessor, u, source(*ei,g)); 
              pred_e = *ei;
            };
          };
          g_u = get(m_distance, u); 
          if(pred_u != get(m_predecessor, u))  
            m_vis.edge_relaxed(pred_e, g);
        };
        
      };
      
      template <class Vertex, class Graph>
      void examine_vertex(Vertex u, Graph& g) const {
        typedef typename Graph::edge_bundled EdgeProp;
        
        m_vis.examine_vertex(u, g);
        
        PositionValue p_new; bool walk_succeeded; EdgeProp ep_new;
        boost::tie(p_new, walk_succeeded, ep_new) = m_vis.random_walk(u, g);
        if(walk_succeeded)
          m_connect_vertex(p_new, u, ep_new, g, m_super_space, *this, 
                           m_position, m_distance, m_predecessor, m_weight,
                           m_select_neighborhood);
        
      };
      
#endif
      
      template <class Vertex, class Graph>
      void examine_vertex(Vertex u, Graph& g) const {
        m_vis.examine_vertex(u, g);
      };
      
      template <typename Graph>
      void publish_path(const Graph& g) const { 
        m_vis.publish_path(g);
      };
      bool keep_going() const { 
        return m_vis.keep_going(); 
      };
      template <typename Vertex, typename Graph>
      bool has_search_potential(Vertex u, const Graph& g) const { 
        return m_vis.has_search_potential(u,g);
      };
      template <typename Vertex, typename Graph>
      bool should_close(Vertex u, const Graph& g) const { 
        return m_vis.should_close(u,g);
      };
      
      template <class Vertex, typename Graph>
      void update_key(Vertex u, Graph& g) const {
        double g_u = get(m_distance, g[u]);
        double h_u = get(m_heuristic, g[u]);
        // Key-value for the min-heap (priority-queue):
        put(m_key, u, (g_u + h_u) / (1.0 - get(m_constriction, g[u])) / (1.0 - get(m_density, g[u])));
      };

      UniformCostVisitor m_vis;
      UpdatableQueue& m_Q; 
      IndexInHeapMap m_index_in_heap;
      
      AStarHeuristicMap m_heuristic;
      PositionMap m_position;
      WeightMap m_weight;
      DensityMap m_density;
      ConstrictionMap m_constriction; 
      DistanceMap m_distance;
      PredecessorMap m_predecessor;
      KeyMap m_key;
      
#if 0
      // Old version of the connector:
      const Topology& m_super_space;
      SBANodeConnector m_connect_vertex;
      NcSelector m_select_neighborhood;
      
#endif
    };
    
    
    
    
    
    
#if 0
    // Old version of the connector:
    struct sbastar_node_connector {
      
      template <typename Graph,
                typename Topology,
                typename SBAVisitor,
                typename PositionMap,
                typename DistanceMap,
                typename PredecessorMap,
                typename WeightMap,
                typename NcSelector>
      typename boost::enable_if< boost::is_undirected_graph<Graph> >::type operator()(
          const typename boost::property_traits<PositionMap>::value_type& p, 
          typename boost::graph_traits<Graph>::vertex_descriptor u, 
          typename Graph::edge_bundled& ep, 
          Graph& g,
          const Topology& super_space,
          const SBAVisitor& sba_vis,
          PositionMap position,
          DistanceMap distance,
          PredecessorMap predecessor,
          WeightMap weight,
          NcSelector select_neighborhood) const {
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
        
        typedef boost::composite_property_map< 
          PositionMap, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t > > GraphPositionMap;
        GraphPositionMap g_position = GraphPositionMap(position, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t >(&g));
        
        std::vector<Vertex> Nc;
        select_neighborhood(p, std::back_inserter(Nc), g, super_space, g_position); 
        
        Vertex v = sba_vis.create_vertex(p, g);
        
        std::pair<Edge, bool> e_new = sba_vis.create_edge(u, v, ep, g);
        if( e_new.second ) {
          put(distance, v, get(distance, u) + get(weight, g[e_new.first]));
          put(predecessor, v, u);
          sba_vis.requeue_vertex(u,g);
        };
        
        for(typename std::vector<Vertex>::iterator it = Nc.begin(); it != Nc.end(); ++it) {
          if(*it == u)
            continue;
          e_new = sba_vis.attempt_connecting_edge(*it, v, g);
          if( e_new.second ) {
            sba_vis.update_vertex(*it,g);
            sba_vis.requeue_vertex(*it,g);
          };
        }; 
        sba_vis.update_vertex(v,g);
        sba_vis.requeue_vertex(v,g);
        
        // need to update all the children of the v node:
        std::stack<Vertex> incons;
        incons.push(v);
        while(!incons.empty()) {
          Vertex s = incons.top(); incons.pop();
          OutEdgeIter eo, eo_end;
          for(boost::tie(eo,eo_end) = out_edges(s,g); eo != eo_end; ++eo) {
            Vertex t = target(*eo, g);
            if(t == s)
              t = source(*eo, g);
            if(s != get(predecessor, t))
              continue;
            put(distance, t, get(distance, s) + get(weight, g[*eo]));
            sba_vis.requeue_vertex(t,g);
            incons.push(t);
          };
        };
      };
      
      template <typename Graph,
                typename Topology,
                typename SBAVisitor,
                typename PositionMap,
                typename DistanceMap,
                typename PredecessorMap,
                typename WeightMap,
                typename NcSelector>
      typename boost::enable_if< boost::is_directed_graph<Graph> >::type operator()(
          const typename boost::property_traits<PositionMap>::value_type& p, 
          typename boost::graph_traits<Graph>::vertex_descriptor u, 
          typename Graph::edge_bundled& ep, 
          Graph& g,
          const Topology& super_space,
          const SBAVisitor& sba_vis,
          PositionMap position,
          DistanceMap distance,
          PredecessorMap predecessor,
          WeightMap weight,
          NcSelector select_neighborhood) const {
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
        
        typedef boost::composite_property_map< 
          PositionMap, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t > > GraphPositionMap;
        GraphPositionMap g_position = GraphPositionMap(position, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t >(&g));
        
        std::vector<Vertex> Pred, Succ;
        select_neighborhood(p, std::back_inserter(Pred), std::back_inserter(Succ), g, super_space, g_position); 
        
        Vertex v = sba_vis.create_vertex(p, g);
        
        std::pair<Edge, bool> e_new = sba_vis.create_edge(u, v, ep, g);
        if( e_new.second ) {
          put(distance, v, get(distance, u) + get(weight, g[e_new.first]));
          put(predecessor, v, u);
          sba_vis.requeue_vertex(u,g);
        };
        
        for(typename std::vector<Vertex>::iterator it = Pred.begin(); it != Pred.end(); ++it) {
          if(*it == u)
            continue;
          sba_vis.attempt_connecting_edge(*it, v, g);
        };
        
        sba_vis.update_vertex(v, g);
        sba_vis.requeue_vertex(v,g);
        
        for(typename std::vector<Vertex>::iterator it = Succ.begin(); it != Succ.end(); ++it) {
          if(sba_vis.attempt_connecting_edge(v, *it, g).second) {
            sba_vis.update_vertex(*it,g);
            sba_vis.requeue_vertex(*it,g);
          };
        }; 
        
        // need to update all the children of the v node:
        std::stack<Vertex> incons;
        incons.push(v);
        while(!incons.empty()) {
          Vertex s = incons.top(); incons.pop();
          OutEdgeIter eo, eo_end;
          for(boost::tie(eo,eo_end) = out_edges(s,g); eo != eo_end; ++eo) {
            Vertex t = target(*eo, g);
            if(t == s)
              t = source(*eo, g);
            if(s != get(predecessor, t))
              continue;
            put(distance, t, get(distance, s) + get(weight, g[*eo]));
            sba_vis.requeue_vertex(t,g);
            incons.push(t);
          };
        };
      };
      
    };
    
    
    template <typename Graph, //this is the actual graph, should comply to BidirectionalGraphConcept.
              typename Vertex, //this is the type to describe a vertex in the graph.
              typename SBAStarBFSVisitor, //this is a visitor class that can perform special operations at event points.
              typename MutableQueue>
    inline void
    sbastar_search_loop(Graph &g, Vertex start_vertex, SBAStarBFSVisitor& sba_vis, MutableQueue& Q)
    { 
      while (sba_vis.keep_going()) {
        
        sba_vis.requeue_vertex(start_vertex,g);
        
        while (!Q.empty() && sba_vis.keep_going()) { 
          Vertex u = Q.top(); Q.pop();
          
          // stop if the best node does not meet the potential threshold.
          if( ! sba_vis.has_search_potential(u, g) )
            break;
          
          sba_vis.examine_vertex(u, g);
          // then push it back on the OPEN queue.
          sba_vis.requeue_vertex(u,g);
          
        }; // end while  (the queue is either empty or it contains vertices that still have low key values.
        
        sba_vis.publish_path(g);
        
      };
    };
    
    
    template <typename NodeConnector,
              typename Graph,
              typename Vertex,
              typename Topology,
              typename SBAStarVisitor,
              typename AStarHeuristicMap,
              typename PositionMap,
              typename WeightMap,
              typename DensityMap,
              typename ConstrictionMap,
              typename DistanceMap,
              typename PredecessorMap,
              typename KeyMap,
              typename NcSelector>
    inline void
    generate_sbastar_no_init_impl
      (Graph &g, Vertex start_vertex, const Topology& super_space, SBAStarVisitor vis,  // basic parameters
      AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
      DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
      PredecessorMap predecessor, KeyMap key,                          // properties resulting from the algorithm
      NcSelector select_neighborhood)
    {
      typedef typename boost::property_traits<KeyMap>::value_type KeyValue;
      typedef std::less<double> KeyCompareType;  // <---- this is a min-heap.
      typedef boost::vector_property_map<std::size_t> IndexInHeapMap;
      IndexInHeapMap index_in_heap;
      {
        typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
        for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
          put(index_in_heap,*ui, static_cast<std::size_t>(-1)); 
        };
      };
      
      typedef boost::d_ary_heap_indirect<Vertex, 4, IndexInHeapMap, KeyMap, KeyCompareType> MutableQueue;
      MutableQueue Q(key, index_in_heap, KeyCompareType()); //priority queue holding the OPEN set.
      
      sbastar_bfs_visitor<
        Topology, 
        SBAStarVisitor,
        NodeConnector,
        MutableQueue, 
        IndexInHeapMap,
        AStarHeuristicMap, 
        PositionMap, 
        WeightMap,
        DensityMap,
        ConstrictionMap, 
        DistanceMap,  
        PredecessorMap,
        KeyMap,
        NcSelector> bfs_vis(super_space, vis, NodeConnector(), Q, index_in_heap, 
                            hval, position, weight, 
                            density, constriction, distance,
                            predecessor, key, select_neighborhood);
      
      put(distance, start_vertex, 0.0);
      
      sbastar_search_loop(g, start_vertex, bfs_vis, Q);
      
    };
    
    
    
#endif
    
    
    
    struct sba_node_generator {
      
      template <typename Graph,
                typename SBAVisitor,
                typename PositionMap>
      inline 
      boost::tuple< typename boost::graph_traits<Graph>::vertex_descriptor,
                    typename boost::property_traits<PositionMap>::value_type,
                    typename Graph::edge_bundled > 
        operator()(typename boost::graph_traits<Graph>::vertex_descriptor u, Graph& g, 
                   const SBAVisitor& sba_vis, PositionMap) const {
        
        typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename Graph::edge_bundled EdgeProp;
        typedef boost::tuple< Vertex, PositionValue, EdgeProp > ResultType;
        
        PositionValue p_new; bool walk_succeeded; EdgeProp ep_new;
        boost::tie(p_new, walk_succeeded, ep_new) = sba_vis.random_walk(u, g);
        if( walk_succeeded )
          return ResultType(u, p_new, ep_new);
        
        return ResultType(boost::graph_traits<Graph>::null_vertex(), PositionValue(), EdgeProp());
      };
      
    };
  
    
    
    template <typename Graph,
              typename Vertex,
              typename Topology,
              typename SBAStarVisitor,
              typename MotionGraphConnector,
              typename SBANodeGenerator,
              typename MutableQueue,
              typename NcSelector>
    inline void
    sbastar_search_loop(Graph &g, Vertex start_vertex, const Topology& super_space, SBAStarVisitor& sba_vis, 
                        MotionGraphConnector connect_vertex, SBANodeGenerator sba_generate_node,
                        MutableQueue& Q, NcSelector select_neighborhood)
    { 
      typedef typename ReaK::pp::topology_traits<Topology>::point_type PositionValue;
      typedef typename Graph::edge_bundled EdgeProp;
      
      while (sba_vis.keep_going()) {
        
        sba_vis.requeue_vertex(start_vertex,g);
        
        while (!Q.empty() && sba_vis.keep_going()) { 
          Vertex u = Q.top(); Q.pop();
          
          // stop if the best node does not meet the potential threshold.
          if( ! sba_vis.has_search_potential(u, g) )
            break;
          
          sba_vis.examine_vertex(u, g);
          
          PositionValue p_new; Vertex x_near; EdgeProp eprop;
          boost::tie(x_near, p_new, eprop) = sba_generate_node(u, g, sba_vis, sba_vis.m_position);
          
          if(x_near != boost::graph_traits<Graph>::null_vertex()) {
            connect_vertex(p_new, x_near, eprop, g, 
                           super_space, sba_vis, sba_vis.m_position, 
                           sba_vis.m_distance, sba_vis.m_predecessor, 
                           sba_vis.m_weight, select_neighborhood);
          };
          
          // then push it back on the OPEN queue.
          sba_vis.requeue_vertex(u,g);
          
        }; // end while  (the queue is either empty or it contains vertices that still have low key values.
        
        sba_vis.publish_path(g);
        
      };
    };
    
    template <typename Graph,
              typename Vertex,
              typename Topology,
              typename SBAStarVisitor,
              typename NodeConnector,
              typename AStarHeuristicMap,
              typename PositionMap,
              typename WeightMap,
              typename DensityMap,
              typename ConstrictionMap,
              typename DistanceMap,
              typename PredecessorMap,
              typename KeyMap,
              typename NcSelector>
    inline void
    generate_sbastar_no_init_impl
      (Graph &g, Vertex start_vertex, const Topology& super_space, SBAStarVisitor vis,  // basic parameters
       NodeConnector connect_vertex, AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
       DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
       PredecessorMap predecessor, KeyMap key,                          // properties resulting from the algorithm
       NcSelector select_neighborhood)
    {
      typedef typename boost::property_traits<KeyMap>::value_type KeyValue;
      typedef std::less<double> KeyCompareType;  // <---- this is a min-heap.
      typedef boost::vector_property_map<std::size_t> IndexInHeapMap;
      IndexInHeapMap index_in_heap;
      {
        typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
        for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
          put(index_in_heap,*ui, static_cast<std::size_t>(-1)); 
        };
      };
      
      typedef boost::d_ary_heap_indirect<Vertex, 4, IndexInHeapMap, KeyMap, KeyCompareType> MutableQueue;
      MutableQueue Q(key, index_in_heap, KeyCompareType()); //priority queue holding the OPEN set.
      
      sbastar_bfs_visitor<
        SBAStarVisitor,
        MutableQueue, 
        IndexInHeapMap,
        AStarHeuristicMap, 
        PositionMap, 
        WeightMap,
        DensityMap,
        ConstrictionMap, 
        DistanceMap,  
        PredecessorMap,
        KeyMap> sba_bfs_vis(vis, Q, index_in_heap, hval, position, weight, 
                            density, constriction, distance, predecessor, key);
      
      put(distance, g[start_vertex], 0.0);
      
      sbastar_search_loop(g, start_vertex, super_space, sba_bfs_vis, 
                          connect_vertex, sba_node_generator(), 
                          Q, select_neighborhood);
      
    };
    
    
    
    
    
    template <typename Graph,
              typename SBAStarVisitor,
              typename DistanceMap,
              typename PredecessorMap,
              typename KeyMap>
    inline void initialize_sbastar_nodes(Graph &g, SBAStarVisitor vis, DistanceMap distance, 
                                         PredecessorMap predecessor, KeyMap key) {
      typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
      for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
        put(distance, g[*ui], std::numeric_limits<double>::infinity());
        put(key, *ui, 0.0);
        put(predecessor, g[*ui], *ui);
        vis.initialize_vertex(*ui, g);
      };
    };
  
  }; //end of detail namespace.
  
  
  
  
  template <typename Graph,
            typename Vertex,
            typename Topology,
            typename SBAStarVisitor,
            typename AStarHeuristicMap,
            typename PositionMap,
            typename WeightMap,
            typename DensityMap,
            typename ConstrictionMap,
            typename DistanceMap,
            typename PredecessorMap,
            typename KeyMap,
            typename NcSelector>
  struct sbastar_bundle {
    typedef Graph             graph_type;
    typedef Vertex            vertex_type;
    typedef Topology          topology_type;
    typedef SBAStarVisitor    visitor_type;
    typedef AStarHeuristicMap heuristic_map_type;
    typedef PositionMap       position_map_type;
    typedef WeightMap         weight_map_type;
    typedef DensityMap        density_map_type;
    typedef ConstrictionMap   constriction_map_type;
    typedef DistanceMap       distance_map_type;
    typedef PredecessorMap    predecessor_map_type;
    typedef KeyMap            key_map_type;
    typedef NcSelector        nc_selector_type;
    
    graph_type* m_g;
    vertex_type m_start_vertex;
    const topology_type* m_super_space;
    visitor_type m_vis;
    heuristic_map_type m_hval;
    position_map_type m_position;
    weight_map_type m_weight;
    density_map_type m_density;
    constriction_map_type m_constriction;
    distance_map_type m_distance;
    predecessor_map_type m_predecessor;
    key_map_type m_key;
    nc_selector_type m_select_neighborhood;
    
    sbastar_bundle(
      Graph &g, Vertex start_vertex, const Topology& super_space, SBAStarVisitor vis,
      AStarHeuristicMap hval, PositionMap position, WeightMap weight,
      DensityMap density, ConstrictionMap constriction, DistanceMap distance,
      PredecessorMap predecessor, KeyMap key, NcSelector select_neighborhood) :
      m_g(&g), m_start_vertex(start_vertex), m_super_space(&super_space), m_vis(vis),
      m_hval(hval), m_position(position), m_weight(weight), 
      m_density(density), m_constriction(constriction), m_distance(distance), 
      m_predecessor(predecessor), m_key(key), m_select_neighborhood(select_neighborhood) { };
    
    
  };
  
  
  
  
  /**
   * This function template creates a bundle of parameters to be fed to any of the
   * SBA* algorithms. This is mainly to simply the interface and the code of all these 
   * different variants of the SBA* algorithm.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam SBAStarVisitor The type of the SBA* visitor to be used, should model the SBAStarVisitorConcept.
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
   * \param hval The property-map of A* heuristic function values for each vertex.
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
   * \param key The property-map which stores the AD* key-values associated to each vertex.
   * \param select_neighborhood A callable object (functor) that can select a list of 
   *        vertices of the graph that ought to be connected to a new 
   *        vertex. The list should be sorted in order of increasing "distance".
   */
  template <typename Graph,
            typename Vertex,
            typename Topology,
            typename SBAStarVisitor,
            typename AStarHeuristicMap,
            typename PositionMap,
            typename WeightMap,
            typename DensityMap,
            typename ConstrictionMap,
            typename DistanceMap,
            typename PredecessorMap,
            typename KeyMap,
            typename NcSelector>
  sbastar_bundle<Graph, Vertex, Topology, SBAStarVisitor, 
                 AStarHeuristicMap, PositionMap, WeightMap, 
                 DensityMap, ConstrictionMap, DistanceMap, 
                 PredecessorMap, KeyMap, NcSelector>
    make_sbastar_bundle(
      Graph &g, Vertex start_vertex, const Topology& super_space, SBAStarVisitor vis, 
      AStarHeuristicMap hval, PositionMap position, WeightMap weight, 
      DensityMap density, ConstrictionMap constriction, DistanceMap distance, 
      PredecessorMap predecessor, KeyMap key, NcSelector select_neighborhood) {
    
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((SBAStarVisitorConcept<SBAStarVisitor,Graph,Topology>));
    
    return sbastar_bundle<Graph, Vertex, Topology, SBAStarVisitor,      
                          AStarHeuristicMap, PositionMap, WeightMap, 
                          DensityMap, ConstrictionMap, DistanceMap, 
                          PredecessorMap, KeyMap, NcSelector>(
                            g, start_vertex, super_space, vis, 
                            hval, position, weight, 
                            density, constriction, distance, 
                            predecessor, key, select_neighborhood);
  };
  
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the SBA* algorithm, without initialization of the existing graph.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   */
  template <typename SBAStarBundle>
  inline void generate_sbastar_no_init(const SBAStarBundle& bdl) {
    
    detail::generate_sbastar_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, motion_graph_connector(),
      bdl.m_hval, bdl.m_position, bdl.m_weight, bdl.m_density, bdl.m_constriction, 
      bdl.m_distance, bdl.m_predecessor, bdl.m_key, bdl.m_select_neighborhood);
    
#if 0
    // Old version of the connector:
    detail::generate_sbastar_no_init_impl< detail::sbastar_node_connector >(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
      bdl.m_hval, bdl.m_position, bdl.m_weight, bdl.m_density, bdl.m_constriction, 
      bdl.m_distance, bdl.m_predecessor, bdl.m_key, bdl.m_select_neighborhood);
#endif
  };
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the SBA* algorithm, with initialization of the existing graph to (re)start the search.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   */
  template <typename SBAStarBundle>
  inline void generate_sbastar(const SBAStarBundle& bdl) {
    
    detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_distance, bdl.m_predecessor, bdl.m_key);
    
    generate_sbastar_no_init(bdl);
    
  };
  
  

};

};

#endif

















/**
 * \file lazy_sbastar.hpp
 *
 * This library provides function templates and concepts that implement a Lazy Sampling-based A* search
 * algorithm. A Lazy-SBA* uses the A* search algorithm to drive the expansion of a roadmap into the free-space 
 * in order to connect a start and goal location. This algorithm has many customization points because there 
 * are many choices to be made in the method, such as how to find nearest neighbors for attempting to 
 * connect them through free-space, how to expand vertices, when to stop the algorithm, etc. 
 * All these customization points are left to the user to implement, some are defined by the 
 * SBAStarVisitorConcept (random-walk, edge-added, etc.).
 *
 * The Lazy-SBA* algorithm is a generalization of the A* algorithm where the neighborhood of a given node of 
 * the motion graph is not defined as a fixed set of neighbors (as in a classic A* over a fixed graph),
 * but rather as a region from which samples can be drawn (biased or not). In an ordinary A* algorithm,
 * vertices are closed when their entire neighborhood has been explored. In an SBA* algorithm, the same 
 * criteria cannot apply since samples could be drawn ad infinitum, so, instead, this concept of the 
 * neighborhood being fully explored is derived from the expected information gained (or conversely, the 
 * "surprisal") from drawing a new sample in the neighborhood. In this lazy version, the computation of the 
 * edge weights as a cost-to-go through the free-space is tentatively replaced by the cost-to-go in the 
 * configuration space (without obstacles), and collision along the path is only performed once the edge 
 * has relaxed (identified as a segment of the local optimal path).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2013
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

#ifndef REAK_LAZY_SBASTAR_HPP
#define REAK_LAZY_SBASTAR_HPP

#include "sbastar_search.hpp"

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

#include <stack>


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Graph */
namespace graph {
  
  namespace detail {
    
    struct lazy_sbastar_node_connector {
      
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
          double tentative_weight = get(ReaK::pp::distance_metric, super_space)(get(position,g[*it]), p, super_space);
          double g_in  = tentative_weight + get(distance, v);
          double g_out = tentative_weight + get(distance, *it);
          if(g_in < get(distance, *it)) {
            // edge is useful as an in-edge to (*it).
            e_new = sba_vis.attempt_connecting_edge(v, *it, g);
            if( e_new.second ) {
              put(distance, *it, g_in);
              Vertex old_pred = get(predecessor, *it);
              put(predecessor, *it, v); 
              sba_vis.edge_relaxed(e_new.first, g);
              remove_edge(old_pred, *it, g);
            };
          } else if(g_out < get(distance, v)) {
            // edge is useful as an in-edge to v.
            e_new = sba_vis.attempt_connecting_edge(*it, v, g);
            if( e_new.second ) {
              put(distance, v, g_out);
              Vertex old_pred = get(predecessor, v);
              put(predecessor, v, *it); 
              sba_vis.edge_relaxed(e_new.first, g);
              remove_edge(old_pred, v, g);
            };
          };
          sba_vis.requeue_vertex(*it,g);
        }; 
        
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
          
          double tentative_weight = get(ReaK::pp::distance_metric, super_space)(get(position,g[*it]), p, super_space);
          double g_out = tentative_weight + get(distance, *it);
          if(g_out < get(distance, v)) {
            // edge is useful as an in-edge to v.
            e_new = sba_vis.attempt_connecting_edge(*it, v, g);
            if( e_new.second ) {
              put(distance, v, g_out);
              Vertex old_pred = get(predecessor, v);
              put(predecessor, v, *it); 
              sba_vis.edge_relaxed(e_new.first, g);
              remove_edge(old_pred, v, g);
            };
          };
          sba_vis.requeue_vertex(*it,g);
        };
        
        sba_vis.requeue_vertex(v,g);
        
        for(typename std::vector<Vertex>::iterator it = Succ.begin(); it != Succ.end(); ++it) {
          
          double tentative_weight = get(ReaK::pp::distance_metric, super_space)(p, get(position,g[*it]), super_space);
          double g_in  = tentative_weight + get(distance, v);
          if(g_in < get(distance, *it)) {
            // edge is useful as an in-edge to (*it).
            e_new = sba_vis.attempt_connecting_edge(v, *it, g);
            if( e_new.second ) {
              put(distance, *it, g_in);
              Vertex old_pred = get(predecessor, *it);
              put(predecessor, *it, v); 
              sba_vis.edge_relaxed(e_new.first, g);
              remove_edge(old_pred, *it, g);
            };
          };
          sba_vis.requeue_vertex(*it,g);
        }; 
        
        sba_vis.requeue_vertex(v,g);
        
        // need to update all the children of the v node:
        std::stack<Vertex> incons;
        incons.push(v);
        while(!incons.empty()) {
          Vertex s = incons.top(); incons.pop();
          OutEdgeIter eo, eo_end;
          for(boost::tie(eo,eo_end) = out_edges(s,g); eo != eo_end; ++eo) {
            Vertex t = target(*eo, g);
            put(distance, t, get(distance, s) + get(weight, g[*eo]));
            
            sba_vis.requeue_vertex(t,g);
            
            incons.push(t);
          };
        };
      };
      
    };
    
  }; //end of detail namespace.
  
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the SBA* algorithm, without initialization of the existing graph.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam SBAStarVisitor The type of the SBA* visitor to be used, should model the SBAStarVisitorConcept.
   * \tparam AStarHeuristicMap This property-map type is used to obtain the heuristic-function values 
   *         for each vertex in the graph.
   * \tparam PositionMap A property-map type that can store the position of each vertex. 
   * \tparam WeightMap This property-map type is used to store the weights of the edge-properties of the 
   *         graph (cost of travel along an edge).
   * \tparam DensityMap A property-map type that can store the probability-measure of the expected common information 
   *         between a new sample and the current neighborhood for each vertex.
   * \tparam ConstrictionMap A property-map type that can store the probability-measure of sampling a colliding point 
   *         for each vertex.
   * \tparam DistanceMap This property-map type is used to store the estimated distance of each vertex 
   *         to the goal.
   * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting 
   *         vertex together with its optimal predecessor.
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
  inline void
  generate_lazy_sbastar_no_init
    (Graph &g, Vertex start_vertex, const Topology& super_space, SBAStarVisitor vis,  // basic parameters
     AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
     DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
     PredecessorMap predecessor, KeyMap key,                         // properties resulting from the algorithm
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
    
    detail::sbastar_bfs_visitor<
      Topology, 
      SBAStarVisitor,
      detail::lazy_sbastar_node_connector,
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
      NcSelector> bfs_vis(super_space, vis, detail::lazy_sbastar_node_connector(), Q, index_in_heap, 
                          hval, position, weight, 
                          density, constriction, distance,
                          predecessor, key, select_neighborhood);
    
    put(distance, start_vertex, 0.0);
    
    detail::sbastar_search_loop(g, start_vertex, bfs_vis, Q);
    
  };


   /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the SBA* algorithm, with initialization of the existing graph to (re)start the search.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam SBAStarVisitor The type of the SBA* visitor to be used, should model the SBAStarVisitorConcept.
   * \tparam AStarHeuristicMap This property-map type is used to obtain the heuristic-function values 
   *         for each vertex in the graph.
   * \tparam PositionMap A property-map type that can store the position of each vertex. 
   * \tparam WeightMap This property-map type is used to store the weights of the edge-properties of the 
   *         graph (cost of travel along an edge).
   * \tparam DensityMap A property-map type that can store the probability-measure of the expected common information 
   *         between a new sample and the current neighborhood for each vertex.
   * \tparam ConstrictionMap A property-map type that can store the probability-measure of sampling a colliding point for each vertex.
   * \tparam DistanceMap This property-map type is used to store the estimated distance of each vertex 
   *         to the goal.
   * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting 
   *         vertex together with its optimal predecessor.
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
  inline void
  generate_lazy_sbastar
    (Graph &g, Vertex start_vertex, const Topology& super_space, SBAStarVisitor vis,  // basic parameters
     AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
     DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
     PredecessorMap predecessor, KeyMap key,                          // properties resulting from the algorithm
     NcSelector select_neighborhood)
  {
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    //BOOST_CONCEPT_ASSERT((boost::MutablePropertyGraphConcept<Graph>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::PointDistributionConcept<Topology>));
    BOOST_CONCEPT_ASSERT((SBAStarVisitorConcept<SBAStarVisitor,Graph,Topology>));
    
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
    typedef typename Graph::vertex_bundled VertexProp;
    
    for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
      put(distance, *ui, std::numeric_limits<double>::infinity());
      put(key, *ui, 0.0);
      put(predecessor, *ui, *ui);
      vis.initialize_vertex(*ui, g);
    };

    generate_lazy_sbastar_no_init(
      g, start_vertex, super_space, vis, 
      hval, position, weight, density, constriction, distance,
      predecessor, key, select_neighborhood);

  };
  
  
  

};

};

#endif

















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

#include "bgl_more_property_maps.hpp"
#include "bgl_more_property_tags.hpp"
#include "bgl_raw_property_graph.hpp"


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
   * vis.examine_edge(e,g);  A function that gets called whenever an edge (e) is being looked at, as it comes out of the vertex that is currently being examined in the search.
   * 
   * vis.finish_vertex(u,g);  A function that gets called whenever a vertex (u) is put into the CLOSED set (after being explored by the current search pass).
   * 
   * vis.edge_relaxed(e,g);  A function that gets called whenever an edge (e) has been newly declared as a better alternative than the current surrounding edges (i.e. the edge is added to the optimal path).
   * 
   * vis.publish_path(g);  A function to notify the visitor that at least one A* round has completed and its resulting path (partial or complete) can be published (the path is encoded in the predecessor property-map).
   * 
   * new_threshold = vis.adjust_threshold(old_threshold, g);  A function to adjust the value of the key-value threshold for a given old-value and graph. The threshold controls what stays on the OPEN set (what still has potential).
   * 
   * b = vis.keep_going();  A function to check to see whether the task is finished (return false) or needs to keep going (true).
   * 
   * vis.vertex_added(u, g);  This function is called whenever a new vertex (u) has been added to the graph (g), but not yet connected.
   * 
   * vis.edge_added(e, g);  This function is called whenever a new edge (e) has been created between the last created vertex and its neighbor in the graph (g).
   * 
   * boost::tie(pt,b) = vis.random_walk(u, g);  This function is called to perform the expansion of the roadmap from a given vertex (u) in the graph (g). This function returns a newly generated position value that is a candidate to be added to the graph.
   * 
   * vis.examine_neighborhood(u, g);  This function is called to evaluate the probability-measures of the graph (g) around the given vertex (u). This value is used to prioritize the generation of new vertices by their potential.
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
      vis.examine_edge(e, g);        //whenever an edge is being looked at (an out_edge of the vertex under examination).
      vis.finish_vertex(u, g);       //whenever a vertex is added to the CLOSED set.
      vis.edge_relaxed(e, g);        //whenever it is newly decided that an edge is relaxed (has improved the distance for its target)
      vis.publish_path(g);           // notify the visitor that at least one A* round has completed and its resulting path (partial or complete) can be published (the path is encoded in the predecessor property-map).
      thresh = vis.adjust_threshold(thresh, g); // adjust the value of epsilon for a given old-value and last cummulative weight-change.
      b = vis.keep_going();          // check to see whether the task is finished (return false) or needs to keep going (true).
      vis.vertex_added(u, g); 
      vis.edge_added(e, g);
      boost::tie(pt,b) = vis.random_walk(u, g);
      vis.examine_neighborhood(u, g);
    }
    Visitor vis;
    Graph g;
    typename boost::graph_traits<Graph>::vertex_descriptor u;
    typename boost::graph_traits<Graph>::edge_descriptor e;
    typename ReaK::pp::topology_traits<Topology>::point_type pt;
    bool b;
    double thresh;
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
      void examine_edge(Edge, const Graph&) const { };
      
      template <typename Vertex, typename Graph>
      void finish_vertex(Vertex, const Graph&) const { };
      
      template <typename Edge, typename Graph>
      void edge_relaxed(Edge, const Graph&) const { };
      
      template <typename Graph>
      void publish_path(const Graph&) const { };
      
      template <typename Graph>
      double adjust_threshold(double old_thr, const Graph&) const { return old_thr * 0.5; };
      
      template <typename Vertex, typename Graph>
      void vertex_added(Vertex, const Graph&) const { };
      
      template <typename Edge, typename Graph>
      void edge_added(Edge, const Graph&) const { };
      
      template <typename Vertex, typename Graph>
      std::pair<PointType, bool> random_walk(Vertex, const Graph&) const { return std::make_pair(PointType(), false); };
      
      bool keep_going() const { return true; };
      
      template <typename Vertex, typename Graph>
      void examine_neighborhood(Vertex, const Graph&) const { };
  };
  
  
  
  namespace detail {
  
    template <typename Topology,
              typename UniformCostVisitor,
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
              typename ColorMap, 
              typename NcSelector>
    struct sbastar_bfs_visitor
    {
      typedef typename boost::property_traits<ColorMap>::value_type ColorValue;
      typedef boost::color_traits<ColorValue> Color;
      typedef typename boost::property_traits<PositionMap>::value_type PositionValue;

      sbastar_bfs_visitor(const Topology& free_space, UniformCostVisitor vis,
                          UpdatableQueue& Q, IndexInHeapMap index_in_heap,  
                          AStarHeuristicMap heuristic, PositionMap pos, WeightMap weight, 
                          DensityMap density, ConstrictionMap constriction, DistanceMap dist, 
                          PredecessorMap pred, KeyMap key, ColorMap col, NcSelector select_neighborhood) : 
                          m_free_space(free_space), m_vis(vis), 
                          m_Q(Q), m_index_in_heap(index_in_heap), 
                          m_heuristic(heuristic), m_position(pos), m_weight(weight),
                          m_density(density), m_constriction(constriction), m_distance(dist),
                          m_predecessor(pred), m_key(key), m_color(col), m_select_neighborhood(select_neighborhood) { };
      
      template <class Vertex, class Graph>
      void initialize_vertex(Vertex u, Graph& g) const {
        m_vis.initialize_vertex(u, g);
      };
      template <class Vertex, class Graph>
      void discover_vertex(Vertex u, Graph& g) const {
        m_vis.discover_vertex(u, g);
      };
      template <class Graph>
      typename boost::enable_if< boost::is_undirected_graph<Graph> >::type connect_vertex(const PositionValue& p, Graph& g) {
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename Graph::vertex_bundled VertexProp;
        
        typedef boost::composite_property_map< 
          PositionMap, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t > > GraphPositionMap;
        GraphPositionMap g_position = GraphPositionMap(m_position, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t >(&g));
        
        std::vector<Vertex> Nc;
        m_select_neighborhood(p, std::back_inserter(Nc), g, m_free_space, g_position); 
        
        VertexProp up;
        put(m_position, up, p);
#ifdef RK_ENABLE_CXX0X_FEATURES
        Vertex u = add_vertex(std::move(up), g);
#else
        Vertex u = add_vertex(up, g);
#endif
        m_vis.vertex_added(u,g);
        put(m_color, u, Color::white());
        put(m_index_in_heap, u, static_cast<std::size_t>(-1));
        put(m_distance, u, std::numeric_limits<double>::infinity());
        put(m_key, u, 0.0);
        put(m_predecessor, u, u);
        
        for(typename std::vector<Vertex>::iterator it = Nc.begin(); it != Nc.end(); ++it) {
          if((u != *it) && (get(ReaK::pp::distance_metric, m_free_space)(get(m_position,g[*it]), p, m_free_space) != std::numeric_limits<double>::infinity())) {
            //this means that u is reachable from *it.
            std::pair<Edge, bool> ep = add_edge(*it,u,g); 
            if(ep.second) { 
              m_vis.edge_added(ep.first, g); 
              update_vertex(*it,g);
            };
          };
        }; 
        update_vertex(u,g);
      };
      template <class Graph>
      typename boost::enable_if< boost::is_directed_graph<Graph> >::type connect_vertex(const PositionValue& p, Graph& g) {
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename Graph::vertex_bundled VertexProp;
        
        typedef boost::composite_property_map< 
          PositionMap, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t > > GraphPositionMap;
        GraphPositionMap g_position = GraphPositionMap(m_position, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t >(&g));
        
        std::vector<Vertex> Pred, Succ;
        m_select_neighborhood(p, std::back_inserter(Pred), std::back_inserter(Succ), g, m_free_space, g_position); 
        
        VertexProp up;
        put(m_position, up, p);
#ifdef RK_ENABLE_CXX0X_FEATURES
        Vertex u = add_vertex(std::move(up), g);
#else
        Vertex u = add_vertex(up, g);
#endif
        m_vis.vertex_added(u,g); 
        put(m_color, u, Color::white());
        put(m_index_in_heap, u, static_cast<std::size_t>(-1));
        put(m_distance, u, std::numeric_limits<double>::infinity());
        put(m_key, u, 0.0);
        put(m_predecessor, u, u);
        
        for(typename std::vector<Vertex>::iterator it = Pred.begin(); it != Pred.end(); ++it) {
          if((u != *it) && (get(ReaK::pp::distance_metric, m_free_space)(get(m_position,g[*it]), p, m_free_space) != std::numeric_limits<double>::infinity())) {
            //this means that u is reachable from *it.
            std::pair<Edge, bool> ep = add_edge(*it, u, g); 
            if(ep.second) {
              m_vis.edge_added(ep.first, g); 
              update_vertex(*it, g);
            };
          };
        };
        
        update_vertex(u, g);
        
        for(typename std::vector<Vertex>::iterator it = Succ.begin(); it != Succ.end(); ++it) {
          if((u != *it) && (get(ReaK::pp::distance_metric, m_free_space)(p, get(m_position,g[*it]), m_free_space) != std::numeric_limits<double>::infinity())) {
            //this means that u is reachable from *it.
            std::pair<Edge, bool> ep = add_edge(u, *it, g); 
            if(ep.second) {
              m_vis.edge_added(ep.first, g); 
              update_vertex(*it, g);
            };
          };
        }; 
      };
      
      template <class Vertex, class Graph>
      void examine_vertex(Vertex u, Graph& g) {
        m_vis.examine_vertex(u, g);
        
        std::pair< PositionValue, bool > p_new = m_vis.random_walk(u, g);
        if(p_new.second)
          connect_vertex(p_new.first, g);
        
        update_key(u,g);
      };
      
      template <class Edge, class Graph>
      void examine_edge(Edge e, Graph& g) const {
        if (get(m_weight, e) < 0.0)
          throw boost::negative_edge();
        m_vis.examine_edge(e, g);
      };
      template <class Vertex, class Graph>
      void finish_vertex(Vertex u, Graph& g) const {
        m_vis.finish_vertex(u, g);
      };
      template <typename Graph>
      void publish_path(const Graph& g) const { 
        m_vis.publish_path(g);
      };
      bool keep_going() const { 
        return m_vis.keep_going(); 
      };
      template <typename Graph>
      double adjust_threshold(double old_threshold, const Graph& g) const {
        return m_vis.adjust_threshold(old_threshold, g);
      };
      
      template <class Vertex, typename Graph>
      void update_key(Vertex u, Graph& g) const {
        m_vis.examine_neighborhood(u, g);  // <--- update the constriction / surprise probabilities from examining the neighborhood of u.
        double g_u = get(m_distance, u);
        double h_u = get(m_heuristic, u);
        double f_u = g_u + h_u;   // <--- no relaxation.
//         double f_u = g_u + 5.0 * h_u;   // <--- with relaxation.
        // Key-value for the min-heap (priority-queue):
        // key[u]  =  P( collision | N(u) ) * (1 - P( surprise | N(u) ) ) * total-distance 
        put(m_key, u, (1.0 - get(m_constriction, u)) * (1.0 - get(m_density, u)) / f_u);
      };

      template <class Vertex, class Graph>
      void update_vertex(Vertex u, Graph& g) {
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
            double g_tmp = get(m_weight, *ei) + get(m_distance, source(*ei,g)); 
            if(g_tmp < g_u) {
              g_u = g_tmp; 
              put(m_distance, u, g_u);
              put(m_predecessor, u, source(*ei,g)); 
              pred_e = *ei;
            };
          };
          g_u = get(m_distance, u); 
          if(pred_u != get(m_predecessor, u))  m_vis.edge_relaxed(pred_e, g);
        };
        
        update_key(u,g); 
        put(m_color, u, Color::gray()); 
        m_Q.push_or_update(u);                 m_vis.discover_vertex(u, g);
      };

      const Topology& m_free_space;
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
      ColorMap m_color;
      NcSelector m_select_neighborhood;
    };
    
    
    
    
    
    
    template <typename VertexListGraph, //this is the actual graph, should comply to BidirectionalGraphConcept.
              typename Vertex, //this is the type to describe a vertex in the graph.
              typename SBAStarBFSVisitor, //this is a visitor class that can perform special operations at event points.
              typename AStarHeuristicMap, //this the map of heuristic function value for each vertex.
              typename DistanceMap,
              typename PredecessorMap,
              typename KeyMap, //this is the map of key values associated to each vertex.
              typename ColorMap, //this is a color map for each vertex, i.e. white=not visited, gray=discovered, black=expanded.
              typename IndexInHeapMap,
              typename MutableQueue>
    inline void
    sbastar_search_loop
      (VertexListGraph &g, Vertex start_vertex, SBAStarBFSVisitor& bfs_vis, 
       AStarHeuristicMap heuristic, DistanceMap distance,
       PredecessorMap predecessor, KeyMap key, ColorMap color,
       IndexInHeapMap index_in_heap, MutableQueue& Q, double potential_threshold)
    {
      typedef typename boost::graph_traits<VertexListGraph>::edge_descriptor Edge;
      typedef typename boost::graph_traits<VertexListGraph>::out_edge_iterator OutEdgeIter;
      typedef typename boost::property_traits<ColorMap>::value_type ColorValue;
      typedef boost::color_traits<ColorValue> Color;
      
      double f_min = get(heuristic, start_vertex);
      
      while (bfs_vis.keep_going()) {
        
        Vertex s = start_vertex;
        put(distance, s, 0.0);
        bfs_vis.update_key(s,g);
        put(color, s, Color::gray());
        Q.push_or_update(s);                    bfs_vis.discover_vertex(s, g);
        
        while (!Q.empty() && bfs_vis.keep_going()) { 
          Vertex u = Q.top(); Q.pop();
          
          bfs_vis.examine_vertex(u, g);
          
          // stop if the best node does not meet the potential threshold.
          if( (get(key, u) < potential_threshold / f_min) ) {
            while(!Q.empty())
              Q.pop();
            break;
          };
          
          // stop if we have a node at the goal
          if( (get(heuristic, u) == 0.0) )
            break;
          
          OutEdgeIter eig, eig_end;
          for (boost::tie(eig, eig_end) = out_edges(u, g); eig != eig_end; ++eig) {
            bfs_vis.examine_edge(*eig, g);  
            bfs_vis.update_vertex(target(*eig, g), g);
          };
          
          // if the node still has a minimally good potential, then push it back on the OPEN queue.
          if( get(key, u) < potential_threshold / f_min ) {
            put(color, u, Color::gray());
            Q.push(u);                          bfs_vis.discover_vertex(u, g);
          } else {
            put(color, u, Color::black());      bfs_vis.finish_vertex(u, g);
          };
          
        }; // end while  (the queue is either empty or it contains vertices that still have low key values.
        
        bfs_vis.publish_path(g);
        
        potential_threshold = bfs_vis.adjust_threshold(potential_threshold, g);
        
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
   * \tparam WeightMap This property-map type is used to store the weights of the edges of the 
   *         graph (cost of travel along an edge).
   * \tparam DensityMap A property-map type that can store the probability-measure of the expected common information 
   *         between a new sample and the current neighborhood for each vertex.
   * \tparam ConstrictionMap A property-map type that can store the probability-measure of sampling a colliding point 
   *         for each vertex.
   * \tparam DistanceMap This property-map type is used to store the estimated distance of each vertex 
   *         to the goal.
   * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting 
   *         vertex together with its optimal predecessor.
   * \tparam KeyMap This property-map type is used to store the weights of the edges of the 
   *         graph (cost of travel along an edge).
   * \tparam ColorMap This property-map type is used to store the color-value of the vertices, colors 
   *         are used to mark vertices by their status in the A* algorithm (white = not visited, 
   *         gray = discovered (in OPEN), black = finished (in CLOSED)).
   * \tparam NcSelector A functor type that can select a list of vertices of the graph that are 
   *         the nearest-neighbors of a given vertex (or some other heuristic to select the neighbors). 
   *         See classes in the topological_search.hpp header-file.
   * 
   * \param g A mutable graph that should initially store the starting 
   *        vertex (if not it will be randomly generated) and will store 
   *        the generated graph once the algorithm has finished.
   * \param start_vertex The starting point of the algorithm, on the graph.
   * \param free_space A topology (as defined by the Boost Graph Library). Note 
   *        that it is required to generate only random points in 
   *        the free-space and to only allow interpolation within the free-space.
   * \param vis A SBA* visitor implementing the FADPRMVisitorConcept. This is the 
   *        main point of customization and recording of results that the 
   *        user can implement.
   * \param hval The property-map of A* heuristic function values for each vertex.
   * \param position A mapping that implements the MutablePropertyMap Concept. Also,
   *        the value_type of this map should be the same type as the topology's 
   *        value_type.
   * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
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
   * \param color The property-map which stores the color-value of the vertices, colors are used to mark
   *        vertices by their status in the AD* algorithm (white = not visited, gray = discovered (in OPEN), 
   *        black = finished (in CLOSED), green = recycled (not CLOSED, not OPEN), red = inconsistent (in INCONS)).
   * \param select_neighborhood A callable object (functor) that can select a list of 
   *        vertices of the graph that ought to be connected to a new 
   *        vertex. The list should be sorted in order of increasing "distance".
   * \param initial_threshold The initial threshold value that determines if vertices should still be in the OPEN
   *        set given their key value. Vertices with key values higher than the threshold are taken off the OPEN set. 
   *        The inner loop of the algorithm terminates when the vertex with the lowest key value is higher than the 
   *        threshold or the priority-queue is empty.
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
            typename ColorMap,
            typename NcSelector>
  inline void
  generate_sbastar_no_init
    (Graph &g, Vertex start_vertex, const Topology& free_space, SBAStarVisitor vis,  // basic parameters
     AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
     DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
     PredecessorMap predecessor, KeyMap key, ColorMap color,                         // properties resulting from the algorithm
     NcSelector select_neighborhood, 
     double initial_threshold)
  {
    typedef typename boost::property_traits<KeyMap>::value_type KeyValue;
    typedef std::greater<double> KeyCompareType;  // <---- this is a max-heap.
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
      ColorMap, 
      NcSelector> bfs_vis(free_space, vis, Q, index_in_heap, 
                          hval, position, weight, 
                          density, constriction, distance,
                          predecessor, key, color, select_neighborhood);
    
    detail::sbastar_search_loop(g, start_vertex, bfs_vis, 
                                hval, distance, predecessor, key, color,
                                index_in_heap, Q, initial_threshold);
    
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
   * \tparam WeightMap This property-map type is used to store the weights of the edges of the 
   *         graph (cost of travel along an edge).
   * \tparam DensityMap A property-map type that can store the probability-measure of the expected common information 
   *         between a new sample and the current neighborhood for each vertex.
   * \tparam ConstrictionMap A property-map type that can store the probability-measure of sampling a colliding point for each vertex.
   * \tparam DistanceMap This property-map type is used to store the estimated distance of each vertex 
   *         to the goal.
   * \tparam PredecessorMap This property-map type is used to store the resulting path by connecting 
   *         vertex together with its optimal predecessor.
   * \tparam KeyMap This property-map type is used to store the weights of the edges of the 
   *         graph (cost of travel along an edge).
   * \tparam ColorMap This property-map type is used to store the color-value of the vertices, colors 
   *         are used to mark vertices by their status in the A* algorithm (white = not visited, 
   *         gray = discovered (in OPEN), black = finished (in CLOSED)).
   * \tparam NcSelector A functor type that can select a list of vertices of the graph that are 
   *         the nearest-neighbors of a given vertex (or some other heuristic to select the neighbors). 
   *         See classes in the topological_search.hpp header-file.
   * 
   * \param g A mutable graph that should initially store the starting 
   *        vertex (if not it will be randomly generated) and will store 
   *        the generated graph once the algorithm has finished.
   * \param start_vertex The starting point of the algorithm, on the graph.
   * \param free_space A topology (as defined by the Boost Graph Library). Note 
   *        that it is required to generate only random points in 
   *        the free-space and to only allow interpolation within the free-space.
   * \param vis A SBA* visitor implementing the FADPRMVisitorConcept. This is the 
   *        main point of customization and recording of results that the 
   *        user can implement.
   * \param hval The property-map of A* heuristic function values for each vertex.
   * \param position A mapping that implements the MutablePropertyMap Concept. Also,
   *        the value_type of this map should be the same type as the topology's 
   *        value_type.
   * \param weight The property-map which stores the weight of each edge of the graph (the cost of travel
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
   * \param color The property-map which stores the color-value of the vertices, colors are used to mark
   *        vertices by their status in the AD* algorithm (white = not visited, gray = discovered (in OPEN), 
   *        black = finished (in CLOSED), green = recycled (not CLOSED, not OPEN), red = inconsistent (in INCONS)).
   * \param select_neighborhood A callable object (functor) that can select a list of 
   *        vertices of the graph that ought to be connected to a new 
   *        vertex. The list should be sorted in order of increasing "distance".
   * \param initial_threshold The initial threshold value that determines if vertices should still be in the OPEN
   *        set given their key value. Vertices with key values higher than the threshold are taken off the OPEN set. 
   *        The inner loop of the algorithm terminates when the vertex with the lowest key value is higher than the 
   *        threshold or the priority-queue is empty.
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
            typename ColorMap,
            typename NcSelector>
  inline void
  generate_sbastar
    (Graph &g, Vertex start_vertex, const Topology& free_space, SBAStarVisitor vis,  // basic parameters
     AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
     DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
     PredecessorMap predecessor, KeyMap key, ColorMap color,                         // properties resulting from the algorithm
     NcSelector select_neighborhood, double initial_threshold)
  {
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    //BOOST_CONCEPT_ASSERT((boost::MutablePropertyGraphConcept<Graph>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::PointDistributionConcept<Topology>));
    BOOST_CONCEPT_ASSERT((SBAStarVisitorConcept<SBAStarVisitor,Graph,Topology>));
    
    typedef typename boost::property_traits<ColorMap>::value_type ColorValue;
    typedef boost::color_traits<ColorValue> Color;
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
    typedef typename Graph::vertex_bundled VertexProp;
    
    if(num_vertices(g) == 0) {
      VertexProp up;
      PositionValue p = get(ReaK::pp::random_sampler, free_space)(free_space);
      put(position, up, p);
#ifdef RK_ENABLE_CXX0X_FEATURES
      Vertex u = add_vertex(std::move(up), g);
#else
      Vertex u = add_vertex(up, g);
#endif
      vis.vertex_added(u, g);
      start_vertex = u;
    };
    
    for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
      put(color, *ui, Color::white());
      put(distance, *ui, std::numeric_limits<double>::infinity());
      put(key, *ui, 0.0);
      put(predecessor, *ui, *ui);
      vis.initialize_vertex(*ui, g);
    };

    generate_sbastar_no_init(
      g, start_vertex, free_space, vis, 
      hval, position, weight, density, constriction, distance,
      predecessor, key, color, select_neighborhood, initial_threshold);

  };
  
  
  

};

};

#endif

















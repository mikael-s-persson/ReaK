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


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Graph */
namespace graph {
  
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
    struct lazy_sbastar_bfs_visitor : 
      sbastar_bfs_visitor< Topology, UniformCostVisitor, UpdatableQueue, IndexInHeapMap,
                           AStarHeuristicMap, PositionMap, WeightMap, DensityMap, ConstrictionMap, 
                           DistanceMap, PredecessorMap, KeyMap, ColorMap, NcSelector >
    {
      typedef sbastar_bfs_visitor< Topology, UniformCostVisitor, UpdatableQueue, IndexInHeapMap,        
                                   AStarHeuristicMap, PositionMap, WeightMap, DensityMap, ConstrictionMap, 
                                   DistanceMap, PredecessorMap, KeyMap, ColorMap, NcSelector > base_type;
      
      typedef typename boost::property_traits<ColorMap>::value_type ColorValue;
      typedef boost::color_traits<ColorValue> Color;
      typedef typename boost::property_traits<PositionMap>::value_type PositionValue;

      lazy_sbastar_bfs_visitor(const Topology& free_space, UniformCostVisitor vis,
                               UpdatableQueue& Q, IndexInHeapMap index_in_heap,  
                               AStarHeuristicMap heuristic, PositionMap pos, WeightMap weight, 
                               DensityMap density, ConstrictionMap constriction, DistanceMap dist, 
                               PredecessorMap pred, KeyMap key, ColorMap col, NcSelector select_neighborhood) : 
                               base_type(free_space, vis, Q, index_in_heap, heuristic, pos, weight, 
                                         density, constriction, dist, pred, key, col, select_neighborhood) { };
      
      template <class Graph>
      typename boost::enable_if< boost::is_undirected_graph<Graph> >::type connect_vertex(const PositionValue& p, Graph& g) {
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename Graph::vertex_bundled VertexProp;                //TODO: Add laziness here.
        
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
        typedef typename Graph::vertex_bundled VertexProp;                //TODO: Add laziness here.
        
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
      void examine_edge(Edge e, Graph& g) const {                //TODO: Add laziness here. (NOT SURE)
        if (get(m_weight, e) < 0.0)
          throw boost::negative_edge();
        m_vis.examine_edge(e, g);
      };
      
      template <class Vertex, class Graph>
      void update_vertex(Vertex u, Graph& g) {                //TODO: Add laziness here.
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
  generate_lazy_sbastar_no_init
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
    
    detail::lazy_sbastar_bfs_visitor<
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
  generate_lazy_sbastar
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

    generate_lazy_sbastar_no_init(
      g, start_vertex, free_space, vis, 
      hval, position, weight, density, constriction, distance,
      predecessor, key, color, select_neighborhood, initial_threshold);

  };
  
  
  

};

};

#endif

















/**
 * \file anytime_sbastar.hpp
 *
 * This library provides function templates and concepts that implement an Anytime Sampling-based A* search
 * algorithm. A ASBA* uses the Anytime A* search algorithm to drive the expansion of a roadmap into the free-space 
 * in order to connect a start and goal location. This algorithm has many customization points because there 
 * are many choices to be made in the method, such as how to find nearest neighbors for attempting to 
 * connect them through free-space, how to expand vertices, when to stop the algorithm, etc. 
 * All these customization points are left to the user to implement, some are defined by the 
 * ASBAStarVisitorConcept (random-walk, edge-added, etc.).
 *
 * The ASBA* algorithm is a generalization of the Anytime A* algorithm where the neighborhood of a given node of 
 * the motion graph is not defined as a fixed set of neighbors (as in a classic A* over a fixed graph),
 * but rather as a region from which samples can be drawn (biased or not). In an ordinary A* algorithm,
 * vertices are closed when their entire neighborhood has been explored. In an ASBA* algorithm, the same 
 * criteria cannot apply since samples could be drawn ad infinitum, so, instead, this concept of the 
 * neighborhood being fully explored is derived from the expected information gained (or conversely, the 
 * "surprisal") from drawing a new sample in the neighborhood.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2013
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

#ifndef REAK_ANYTIME_SBASTAR_HPP
#define REAK_ANYTIME_SBASTAR_HPP

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

#include "sbastar_search.hpp"
#include "lazy_sbastar.hpp"
#include "sbastar_rrtstar.hpp"

#include <stack>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Graph */
namespace graph {

  /**
   * This concept class defines the valid expressions required of a class to be used as a visitor 
   * class for the ASBA* algorithm. A visitor class is essentially a class that regroups a number of 
   * callback functions that can be used to inject customization into the ASBA* algorithm. In other 
   * words, the visitor pattern in generic programming is an implementation of IoC 
   * (Inversion of Control), since the ASBA* algorithm is in control of execution, but custom behavior can
   * be injected in several places, even blocking the algorithm if needed.
   * 
   * Required concepts:
   * 
   * The visitor class should model the SBAStarVisitorConcept.
   * 
   * Valid expressions:
   * 
   * d = vis.adjust_relaxation(d, g);  This function should return a new value for the relaxation factor used in the ASBA* algorithm.
   * 
   * \tparam Visitor The visitor class to be tested for modeling an ASBA* visitor concept.
   * \tparam Graph The graph type on which the visitor should be able to act.
   * \tparam Topology The topology type on which the visitor class is required to work with.
   */
  template <typename Visitor, typename Graph, typename Topology>
  struct ASBAStarVisitorConcept : SBAStarVisitorConcept<Visitor, Graph, Topology> {
    BOOST_CONCEPT_USAGE(ASBAStarVisitorConcept)
    {
      d = this->vis.adjust_relaxation(d, this->g); 
    }
    double d;
  };
  
  
  /**
   * This concept class defines the valid expressions required of a class to be used as a visitor 
   * class for the ASBA*-RRT* algorithm. A visitor class is essentially a class that regroups a number of 
   * callback functions that can be used to inject customization into the ASBA*-RRT* algorithm. In other 
   * words, the visitor pattern in generic programming is an implementation of IoC 
   * (Inversion of Control), since the ASBA*-RRT* algorithm is in control of execution, but custom behavior can
   * be injected in several places, even blocking the algorithm if needed.
   * 
   * Required concepts:
   * 
   * The visitor class should model the SBARRTStarVisitorConcept.
   * 
   * Valid expressions:
   * 
   * d = vis.adjust_relaxation(d, g);  This function should return a new value for the relaxation factor used in the ASBA* algorithm.
   * 
   * \tparam Visitor The visitor class to be tested for modeling an ASBA*-RRT* visitor concept.
   * \tparam Graph The graph type on which the visitor should be able to act.
   * \tparam Topology The topology type on which the visitor class is required to work with.
   */
  template <typename Visitor, typename Graph, typename Topology>
  struct ASBARRTStarVisitorConcept : SBARRTStarVisitorConcept<Visitor, Graph, Topology> {
    BOOST_CONCEPT_USAGE(ASBARRTStarVisitorConcept)
    {
      d = this->vis.adjust_relaxation(d, this->g); 
    }
    double d;
  };
  
  
  /**
   * This class is simply a "null" visitor for the ASBA* algorithm. It is null in the sense that it
   * will do nothing on all accounts.
   */
  template <typename Topology>
  class default_asbastar_visitor : public default_sbastar_visitor<Topology> {
    public:
      typedef typename default_sbastar_visitor<Topology>::PointType PointType;
      
      default_asbastar_visitor() : default_sbastar_visitor<Topology>() {};
      
      template <typename Graph>
      double initialize_vertex(double d, const Graph&) const { return d; };
      
      template <typename Vertex, typename Graph>
      boost::tuple<PointType, bool, typename Graph::edge_bundled> steer_towards_position(const PointType&, Vertex, const Graph&) const { 
        typedef typename Graph::edge_bundled EdgeProp;
        typedef boost::tuple<PointType, bool, EdgeProp> ResultType;
        return ResultType(PointType(), false, EdgeProp()); 
      };
  };
  
  
  
  
  
  namespace detail {
  
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
    struct anytime_sbastar_bfs_visitor : 
      sbastar_bfs_visitor<Topology, UniformCostVisitor, SBANodeConnector, 
                          UpdatableQueue, IndexInHeapMap, AStarHeuristicMap, 
                          PositionMap, WeightMap, DensityMap, ConstrictionMap, 
                          DistanceMap, PredecessorMap, KeyMap, NcSelector>
    {
      typedef sbastar_bfs_visitor<Topology, UniformCostVisitor, SBANodeConnector, 
                                  UpdatableQueue, IndexInHeapMap, AStarHeuristicMap, 
                                  PositionMap, WeightMap, DensityMap, ConstrictionMap, 
                                  DistanceMap, PredecessorMap, KeyMap, NcSelector> base_type;
      
      typedef typename base_type::PositionValue PositionValue;

      anytime_sbastar_bfs_visitor(const Topology& super_space, UniformCostVisitor vis, SBANodeConnector connector,
                                  UpdatableQueue& Q, IndexInHeapMap index_in_heap,  
                                  AStarHeuristicMap heuristic, PositionMap pos, WeightMap weight, 
                                  DensityMap density, ConstrictionMap constriction, DistanceMap dist, 
                                  PredecessorMap pred, KeyMap key, NcSelector select_neighborhood,
                                  double current_relaxation) : 
                                  base_type(super_space, vis, connector, Q, index_in_heap, heuristic, 
                                            pos, weight, density, constriction, dist, pred, key, select_neighborhood),
                                  m_current_relaxation(current_relaxation) { };
      
      template <typename Vertex, typename Graph>
      void requeue_vertex(Vertex u, const Graph& g) const { 
        update_key(u,g); 
        if( ! this->m_vis.should_close(u, g) ) {
          this->m_Q.push_or_update(u);
          this->m_vis.discover_vertex(u, g);
        };
      };
      
      template <class Vertex, class Graph>
      void examine_vertex(Vertex u, Graph& g) const {
        typedef typename Graph::edge_bundled EdgeProp;
        
        this->m_vis.examine_vertex(u, g);
        
        PositionValue p_new; bool walk_succeeded; EdgeProp ep_new;
        boost::tie(p_new, walk_succeeded, ep_new) = this->m_vis.random_walk(u, g);
        if(walk_succeeded)
          this->m_connect_vertex(p_new, u, ep_new, g, this->m_super_space, *this, 
                                 this->m_position, this->m_distance, this->m_predecessor, this->m_weight,
                                 this->m_select_neighborhood);
        
      };
      
      template <class Vertex, typename Graph>
      void update_key(Vertex u, Graph& g) const {
        double g_u = get(this->m_distance, u);
        double h_u = get(this->m_heuristic, u);
        // Key-value for the min-heap (priority-queue):
//         put(this->m_key, u, ((g_u + h_u) / (1.0 - get(this->m_constriction, u)) + m_current_relaxation * h_u) / (1.0 - get(this->m_density, u)));
        put(this->m_key, u, ((g_u + h_u) / (1.0 - get(this->m_constriction, u))) / (1.0 - get(this->m_density, u)) + m_current_relaxation * h_u);
      };
      
      template <typename Graph>
      void update_relaxation(const Graph& g) { 
        m_current_relaxation = this->m_vis.adjust_relaxation(m_current_relaxation, g);
        
        typedef typename boost::graph_traits<Graph>::vertex_iterator VIter;
        VIter vi, vi_end;
        for(boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
          requeue_vertex(*vi, g);
      };
      
      template <typename Graph>
      void publish_path(const Graph& g) { 
        this->m_vis.publish_path(g);
        update_relaxation(g);
      };
      
      template <typename Vertex, typename Graph>
      bool has_search_potential(Vertex u, const Graph& g) const { 
        return this->m_vis.has_search_potential(u,g);
      };
      
      double m_current_relaxation;
    };
    
    
    template <typename Topology,
              typename UniformCostVisitor,
              typename SBANodeConnector,
              typename RRTNodeGenerator,
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
    struct anytime_sbarrtstar_bfs_visitor :
      anytime_sbastar_bfs_visitor<Topology, UniformCostVisitor, SBANodeConnector, UpdatableQueue, IndexInHeapMap, 
                                  AStarHeuristicMap, PositionMap, WeightMap, DensityMap, ConstrictionMap, 
                                  DistanceMap, PredecessorMap, KeyMap, NcSelector>
    {
      typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
      typedef anytime_sbastar_bfs_visitor<Topology, UniformCostVisitor, SBANodeConnector, UpdatableQueue, IndexInHeapMap, 
                                          AStarHeuristicMap, PositionMap, WeightMap, DensityMap, ConstrictionMap, 
                                          DistanceMap, PredecessorMap, KeyMap, NcSelector> base_type;

      anytime_sbarrtstar_bfs_visitor(const Topology& super_space, UniformCostVisitor vis, 
                                     SBANodeConnector connector, RRTNodeGenerator generator,
                                     UpdatableQueue& Q, IndexInHeapMap index_in_heap,  
                                     AStarHeuristicMap heuristic, PositionMap pos, WeightMap weight, 
                                     DensityMap density, ConstrictionMap constriction, DistanceMap dist, 
                                     PredecessorMap pred, KeyMap key, NcSelector select_neighborhood,
                                     double current_relaxation) : 
                                     base_type(super_space, vis, connector, Q, index_in_heap, 
                                               heuristic, pos, weight, density, constriction, 
                                               dist, pred, key, select_neighborhood, current_relaxation),
                                     m_node_generator(generator) { };
      
      template <typename Graph>
      void add_exploratory_node(Graph& g) const {
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename Graph::edge_bundled EdgeProp;
        
        typedef boost::composite_property_map< 
          PositionMap, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t > > GraphPositionMap;
        GraphPositionMap g_position = GraphPositionMap(this->m_position, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t >(&g));
        
        while (true) {
          boost::tuple< Vertex, PositionValue, EdgeProp > gen_result = m_node_generator(g, this->m_vis, g_position);
          
          if(get(this->m_distance, get<0>(gen_result)) != std::numeric_limits<double>::infinity()) {
            this->m_connect_vertex(get<1>(gen_result), get<0>(gen_result), get<2>(gen_result), 
                                   g, this->m_super_space, *this, 
                                   this->m_position, this->m_distance, this->m_predecessor, this->m_weight,
                                   this->m_select_neighborhood);
            return;
          };
        };
      };
      
      RRTNodeGenerator m_node_generator;
    };
    
    
    
    
    template <typename NodeConnector,
              typename Graph,
              typename Vertex,
              typename Topology,
              typename ASBAStarVisitor,
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
    generate_anytime_sbastar_no_init_impl(
      Graph &g, Vertex start_vertex, const Topology& super_space, ASBAStarVisitor vis,  // basic parameters
      AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
      DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
      PredecessorMap predecessor, KeyMap key,                          // properties resulting from the algorithm
      NcSelector select_neighborhood, double init_relaxation)
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
      
      detail::anytime_sbastar_bfs_visitor<
        Topology, 
        ASBAStarVisitor,
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
                            predecessor, key, select_neighborhood, init_relaxation);
      
      put(distance, start_vertex, 0.0);
      
      detail::sbastar_search_loop(g, start_vertex, bfs_vis, Q);
      
    };
    
    
    
    template <typename NodeConnector,
              typename Graph,
              typename Vertex,
              typename Topology,
              typename SBARRTStarVisitor,
              typename AStarHeuristicMap,
              typename PositionMap,
              typename WeightMap,
              typename DensityMap,
              typename ConstrictionMap,
              typename DistanceMap,
              typename PredecessorMap,
              typename KeyMap,
              typename RandomSampler,
              typename NcSelector>
    inline void
    generate_anytime_sbarrtstar_no_init_impl(
      Graph &g, Vertex start_vertex, const Topology& super_space, SBARRTStarVisitor vis,  // basic parameters
      AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
      DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
      PredecessorMap predecessor, KeyMap key,                          // properties resulting from the algorithm
      RandomSampler get_sample,
      NcSelector select_neighborhood, 
      double init_relaxation,
      double SA_init_temperature)
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
      
      detail::anytime_sbarrtstar_bfs_visitor<
        Topology, 
        SBARRTStarVisitor,
        NodeConnector,
        rrg_node_generator<Topology, RandomSampler, NcSelector>,
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
        NcSelector> bfs_vis(super_space, vis, NodeConnector(), 
                            rrg_node_generator<Topology, RandomSampler, NcSelector>(&super_space, get_sample, select_neighborhood), 
                            Q, index_in_heap, 
                            hval, position, weight, 
                            density, constriction, distance,
                            predecessor, key, select_neighborhood, init_relaxation);
      
      put(distance, start_vertex, 0.0);
      
      detail::sbarrtstar_search_loop(g, start_vertex, bfs_vis, Q, SA_init_temperature);
      
    };
    
    
  }; //end of detail namespace.
  
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime SBA* algorithm, without initialization of the existing graph.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam ASBAStarVisitor The type of the ASBA* visitor to be used, should model the ASBAStarVisitorConcept.
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
   * \param vis A ASBA* visitor implementing the ASBAStarVisitorConcept. This is the 
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
   * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
   *        Should be greater than 0, the recommeded value is 10.
   */
  template <typename Graph,
            typename Vertex,
            typename Topology,
            typename ASBAStarVisitor,
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
  generate_anytime_sbastar_no_init(
    Graph &g, Vertex start_vertex, const Topology& super_space, ASBAStarVisitor vis,  // basic parameters
    AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
    DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
    PredecessorMap predecessor, KeyMap key,                          // properties resulting from the algorithm
    NcSelector select_neighborhood, double init_relaxation)
  {
    detail::generate_anytime_sbastar_no_init_impl< detail::sbastar_node_connector >(
      g, start_vertex, super_space, vis, 
      hval, position, weight, density, constriction, 
      distance, predecessor, key, select_neighborhood, init_relaxation);
  };
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime-SBA* algorithm, without initialization of the existing graph.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
   *        Should be greater than 0, the recommeded value is 10.
   */
  template <typename SBAStarBundle>
  inline void generate_anytime_sbastar_no_init(const SBAStarBundle& bdl, double init_relaxation) {
    detail::generate_anytime_sbastar_no_init_impl< detail::sbastar_node_connector >(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
      bdl.m_hval, bdl.m_position, bdl.m_weight, bdl.m_density, bdl.m_constriction, 
      bdl.m_distance, bdl.m_predecessor, bdl.m_key, 
      bdl.m_select_neighborhood, init_relaxation);
  };


   /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime SBA* algorithm, with initialization of the existing graph to (re)start the search.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam ASBAStarVisitor The type of the ASBA* visitor to be used, should model the ASBAStarVisitorConcept.
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
   * \param vis A ASBA* visitor implementing the ASBAStarVisitorConcept. This is the 
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
   * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
   *        Should be greater than 0, the recommeded value is 10.
   */
  template <typename Graph,
            typename Vertex,
            typename Topology,
            typename ASBAStarVisitor,
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
  generate_anytime_sbastar(
    Graph &g, Vertex start_vertex, const Topology& super_space, ASBAStarVisitor vis,  // basic parameters
    AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
    DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
    PredecessorMap predecessor, KeyMap key,                         // properties resulting from the algorithm
    NcSelector select_neighborhood, double init_relaxation)
  {
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((ASBAStarVisitorConcept<ASBAStarVisitor,Graph,Topology>));
    
    detail::initialize_sbastar_nodes(g, vis, distance, predecessor, key);
    
    generate_anytime_sbastar_no_init(
      g, start_vertex, super_space, vis, 
      hval, position, weight, density, constriction, distance,
      predecessor, key, select_neighborhood, init_relaxation);
    
  };
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime-SBA* algorithm, with initialization of the existing graph to (re)start the search.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   * \param init_relaxation The initial temperature of the Simulated Annealing when used 
   *        as the deciding factor between using RRT* or SBA* samples.
   */
  template <typename SBAStarBundle>
  inline void generate_anytime_sbastar(const SBAStarBundle& bdl, double init_relaxation) {
    BOOST_CONCEPT_ASSERT((ASBAStarVisitorConcept<typename SBAStarBundle::visitor_type,typename SBAStarBundle::graph_type,typename SBAStarBundle::topology_type>));
    
    detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_distance, bdl.m_predecessor, bdl.m_key);
    
    generate_anytime_sbastar_no_init(bdl, init_relaxation);
    
  };
  
  
  
  
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime Lazy-SBA* algorithm, without initialization of the existing graph.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam ASBAStarVisitor The type of the ASBA* visitor to be used, should model the ASBAStarVisitorConcept.
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
   * \param vis A ASBA* visitor implementing the ASBAStarVisitorConcept. This is the 
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
   * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
   *        Should be greater than 0, the recommeded value is 10.
   */
  template <typename Graph,
            typename Vertex,
            typename Topology,
            typename ASBAStarVisitor,
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
  generate_anytime_lazy_sbastar_no_init(
    Graph &g, Vertex start_vertex, const Topology& super_space, ASBAStarVisitor vis,  // basic parameters
    AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
    DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
    PredecessorMap predecessor, KeyMap key,                          // properties resulting from the algorithm
    NcSelector select_neighborhood, double init_relaxation)
  {
    detail::generate_anytime_sbastar_no_init_impl< detail::lazy_sbastar_node_connector >(
      g, start_vertex, super_space, vis, 
      hval, position, weight, density, constriction, 
      distance, predecessor, key, select_neighborhood, init_relaxation);
  };
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime-Lazy-SBA* algorithm, without initialization of the existing graph.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
   *        Should be greater than 0, the recommeded value is 10.
   */
  template <typename SBAStarBundle>
  inline void generate_anytime_lazy_sbastar_no_init(const SBAStarBundle& bdl, double init_relaxation) {
    detail::generate_anytime_sbastar_no_init_impl< detail::lazy_sbastar_node_connector >(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
      bdl.m_hval, bdl.m_position, bdl.m_weight, bdl.m_density, bdl.m_constriction, 
      bdl.m_distance, bdl.m_predecessor, bdl.m_key, 
      bdl.m_select_neighborhood, init_relaxation);
  };


   /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the search.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam ASBAStarVisitor The type of the ASBA* visitor to be used, should model the ASBAStarVisitorConcept.
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
   * \param vis A ASBA* visitor implementing the ASBAStarVisitorConcept. This is the 
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
   * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
   *        Should be greater than 0, the recommeded value is 10.
   */
  template <typename Graph,
            typename Vertex,
            typename Topology,
            typename ASBAStarVisitor,
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
  generate_anytime_lazy_sbastar(
    Graph &g, Vertex start_vertex, const Topology& super_space, ASBAStarVisitor vis,  // basic parameters
    AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
    DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
    PredecessorMap predecessor, KeyMap key,                         // properties resulting from the algorithm
    NcSelector select_neighborhood, double init_relaxation)
  {
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((ASBAStarVisitorConcept<ASBAStarVisitor,Graph,Topology>));
    
    detail::initialize_sbastar_nodes(g, vis, distance, predecessor, key);
    
    generate_anytime_lazy_sbastar_no_init(
      g, start_vertex, super_space, vis, 
      hval, position, weight, density, constriction, distance,
      predecessor, key, select_neighborhood, init_relaxation);
    
  };
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime-Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the search.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   * \param init_relaxation The initial temperature of the Simulated Annealing when used 
   *        as the deciding factor between using RRT* or SBA* samples.
   */
  template <typename SBAStarBundle>
  inline void generate_anytime_lazy_sbastar(const SBAStarBundle& bdl, double init_relaxation) {
    BOOST_CONCEPT_ASSERT((ASBAStarVisitorConcept<typename SBAStarBundle::visitor_type,typename SBAStarBundle::graph_type,typename SBAStarBundle::topology_type>));
    
    detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_distance, bdl.m_predecessor, bdl.m_key);
    
    generate_anytime_lazy_sbastar_no_init(bdl, init_relaxation);
    
  };
  
  
  
  
  
  
  
  
  
  
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime SBA*-RRT* algorithm, without initialization of the existing graph.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam SBARRTStarVisitor The type of the Anytime SBA*-RRT* visitor to be used, should model the ASBARRTStarVisitorConcept.
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
   * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
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
   * \param vis A Anytime SBA*-RRT* visitor implementing the SBARRTStarVisitorConcept. This is the 
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
   * \param get_sample A random sampler of positions in the space.
   * \param select_neighborhood A callable object (functor) that can select a list of 
   *        vertices of the graph that ought to be connected to a new 
   *        vertex. The list should be sorted in order of increasing "distance".
   * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
   *        Should be greater than 0, the recommeded value is 10.
   * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
   *        as the deciding factor between using RRT* or SBA* samples.
   */
  template <typename Graph,
            typename Vertex,
            typename Topology,
            typename SBARRTStarVisitor,
            typename AStarHeuristicMap,
            typename PositionMap,
            typename WeightMap,
            typename DensityMap,
            typename ConstrictionMap,
            typename DistanceMap,
            typename PredecessorMap,
            typename KeyMap,
            typename RandomSampler,
            typename NcSelector>
  inline void
  generate_anytime_sbarrtstar_no_init(
    Graph &g, Vertex start_vertex, const Topology& super_space, SBARRTStarVisitor vis,  // basic parameters
    AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
    DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
    PredecessorMap predecessor, KeyMap key,                          // properties resulting from the algorithm
    RandomSampler get_sample,
    NcSelector select_neighborhood, 
    double init_relaxation,
    double SA_init_temperature = 1.0)
  {
    detail::generate_anytime_sbarrtstar_no_init_impl< detail::sbastar_node_connector >(
      g, start_vertex, super_space, vis, hval, 
      position, weight, density, constriction, 
      distance, predecessor, key, get_sample, 
      select_neighborhood, init_relaxation, SA_init_temperature);
  };
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime-SBA*-RRT* algorithm, without initialization of the existing graph.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   * \param get_sample A random sampler of positions in the space.
   * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
   *        Should be greater than 0, the recommeded value is 10.
   * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
   *        as the deciding factor between using RRT* or SBA* samples.
   */
  template <typename SBAStarBundle, 
            typename RandomSampler>
  inline void generate_anytime_sbarrtstar_no_init(const SBAStarBundle& bdl, 
                                                  RandomSampler get_sample, 
                                                  double init_relaxation,
                                                  double SA_init_temperature = 1.0) {
    detail::generate_anytime_sbarrtstar_no_init_impl< detail::sbastar_node_connector >(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
      bdl.m_hval, bdl.m_position, bdl.m_weight, bdl.m_density, bdl.m_constriction, 
      bdl.m_distance, bdl.m_predecessor, bdl.m_key, get_sample, 
      bdl.m_select_neighborhood, init_relaxation, SA_init_temperature);
  };


   /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam SBARRTStarVisitor The type of the Anytime SBA*-RRT* visitor to be used, should model the ASBARRTStarVisitorConcept.
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
   * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
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
   * \param vis A Anytime SBA*-RRT* visitor implementing the ASBARRTStarVisitorConcept. This is the 
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
   * \param get_sample A random sampler of positions in the space.
   * \param select_neighborhood A callable object (functor) that can select a list of 
   *        vertices of the graph that ought to be connected to a new 
   *        vertex. The list should be sorted in order of increasing "distance".
   * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
   *        Should be greater than 0, the recommeded value is 10.
   * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
   *        as the deciding factor between using RRT* or SBA* samples.
   */
  template <typename Graph,
            typename Vertex,
            typename Topology,
            typename SBARRTStarVisitor,
            typename AStarHeuristicMap,
            typename PositionMap,
            typename WeightMap,
            typename DensityMap,
            typename ConstrictionMap,
            typename DistanceMap,
            typename PredecessorMap,
            typename KeyMap,
            typename RandomSampler,
            typename NcSelector>
  inline void
  generate_anytime_sbarrtstar(
    Graph &g, Vertex start_vertex, const Topology& super_space, SBARRTStarVisitor vis,  // basic parameters
    AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
    DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
    PredecessorMap predecessor, KeyMap key,                         // properties resulting from the algorithm
    RandomSampler get_sample,
    NcSelector select_neighborhood, 
    double init_relaxation,
    double SA_init_temperature = 1.0)
  {
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<SBARRTStarVisitor,Graph,Topology>));
    
    detail::initialize_sbastar_nodes(g, vis, distance, predecessor, key);
    
    generate_anytime_sbarrtstar_no_init(
      g, start_vertex, super_space, vis, 
      hval, position, weight, density, constriction, distance,
      predecessor, key, get_sample, select_neighborhood, 
      init_relaxation, SA_init_temperature);
    
  };
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   * \param get_sample A random sampler of positions in the space.
   * \param init_relaxation The initial temperature of the Simulated Annealing when used 
   *        as the deciding factor between using RRT* or SBA* samples.
   * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
   *        as the deciding factor between using RRT* or SBA* samples.
   */
  template <typename SBAStarBundle,
            typename RandomSampler>
  inline void generate_anytime_sbarrtstar(const SBAStarBundle& bdl, 
                                          RandomSampler get_sample, 
                                          double init_relaxation,
                                          double SA_init_temperature = 1.0) {
    BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<typename SBAStarBundle::visitor_type,typename SBAStarBundle::graph_type,typename SBAStarBundle::topology_type>));
    
    detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_distance, bdl.m_predecessor, bdl.m_key);
    
    generate_anytime_sbarrtstar_no_init(bdl, get_sample, init_relaxation, SA_init_temperature);
    
  };
  
  
  
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime Lazy-SBA*-RRT* algorithm, without initialization of the existing graph.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam SBARRTStarVisitor The type of the Anytime Lazy-SBA*-RRT* visitor to be used, should model the ASBARRTStarVisitorConcept.
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
   * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
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
   * \param vis A Anytime Lazy-SBA*-RRT* visitor implementing the ASBARRTStarVisitorConcept. This is the 
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
   * \param get_sample A random sampler of positions in the space.
   * \param select_neighborhood A callable object (functor) that can select a list of 
   *        vertices of the graph that ought to be connected to a new 
   *        vertex. The list should be sorted in order of increasing "distance".
   * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
   *        Should be greater than 0, the recommeded value is 10.
   * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
   *        as the deciding factor between using RRT* or SBA* samples.
   */
  template <typename Graph,
            typename Vertex,
            typename Topology,
            typename SBARRTStarVisitor,
            typename AStarHeuristicMap,
            typename PositionMap,
            typename WeightMap,
            typename DensityMap,
            typename ConstrictionMap,
            typename DistanceMap,
            typename PredecessorMap,
            typename KeyMap,
            typename RandomSampler,
            typename NcSelector>
  inline void
  generate_anytime_lazy_sbarrtstar_no_init(
    Graph &g, Vertex start_vertex, const Topology& super_space, SBARRTStarVisitor vis,  // basic parameters
    AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
    DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
    PredecessorMap predecessor, KeyMap key,                          // properties resulting from the algorithm
    RandomSampler get_sample,
    NcSelector select_neighborhood, 
    double init_relaxation,
    double SA_init_temperature = 1.0)
  {
    detail::generate_anytime_sbarrtstar_no_init_impl< detail::lazy_sbastar_node_connector >(
      g, start_vertex, super_space, vis, hval, 
      position, weight, density, constriction, 
      distance, predecessor, key, get_sample, 
      select_neighborhood, init_relaxation, SA_init_temperature);
  };
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime-SBA*-RRT* algorithm, without initialization of the existing graph.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   * \param get_sample A random sampler of positions in the space.
   * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
   *        Should be greater than 0, the recommeded value is 10.
   * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
   *        as the deciding factor between using RRT* or SBA* samples.
   */
  template <typename SBAStarBundle, 
            typename RandomSampler>
  inline void generate_anytime_lazy_sbarrtstar_no_init(const SBAStarBundle& bdl, 
                                                       RandomSampler get_sample, 
                                                       double init_relaxation,
                                                       double SA_init_temperature = 1.0) {
    detail::generate_anytime_sbarrtstar_no_init_impl< detail::lazy_sbastar_node_connector >(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
      bdl.m_hval, bdl.m_position, bdl.m_weight, bdl.m_density, bdl.m_constriction, 
      bdl.m_distance, bdl.m_predecessor, bdl.m_key, get_sample, 
      bdl.m_select_neighborhood, init_relaxation, SA_init_temperature);
  };


   /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime Lazy-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam SBARRTStarVisitor The type of the Anytime Lazy-SBA*-RRT* visitor to be used, should model the ASBARRTStarVisitorConcept.
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
   * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
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
   * \param vis A Anytime Lazy-SBA*-RRT* visitor implementing the ASBARRTStarVisitorConcept. This is the 
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
   * \param get_sample A random sampler of positions in the space.
   * \param select_neighborhood A callable object (functor) that can select a list of 
   *        vertices of the graph that ought to be connected to a new 
   *        vertex. The list should be sorted in order of increasing "distance".
   * \param init_relaxation The initial relaxation factor to use when computing the ASBA* key values.
   *        Should be greater than 0, the recommeded value is 10.
   * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
   *        as the deciding factor between using RRT* or SBA* samples.
   */
  template <typename Graph,
            typename Vertex,
            typename Topology,
            typename SBARRTStarVisitor,
            typename AStarHeuristicMap,
            typename PositionMap,
            typename WeightMap,
            typename DensityMap,
            typename ConstrictionMap,
            typename DistanceMap,
            typename PredecessorMap,
            typename KeyMap,
            typename RandomSampler,
            typename NcSelector>
  inline void
  generate_anytime_lazy_sbarrtstar(
    Graph &g, Vertex start_vertex, const Topology& super_space, SBARRTStarVisitor vis,  // basic parameters
    AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
    DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
    PredecessorMap predecessor, KeyMap key,                         // properties resulting from the algorithm
    RandomSampler get_sample,
    NcSelector select_neighborhood, 
    double init_relaxation,
    double SA_init_temperature = 1.0)
  {
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<SBARRTStarVisitor,Graph,Topology>));
    
    detail::initialize_sbastar_nodes(g, vis, distance, predecessor, key);
    
    generate_anytime_lazy_sbarrtstar_no_init(
      g, start_vertex, super_space, vis, 
      hval, position, weight, density, constriction, distance,
      predecessor, key, get_sample, select_neighborhood, 
      init_relaxation, SA_init_temperature);
    
  };
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the Anytime-Lazy-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
   * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
   * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
   * \param get_sample A random sampler of positions in the space.
   * \param init_relaxation The initial temperature of the Simulated Annealing when used 
   *        as the deciding factor between using RRT* or SBA* samples.
   * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
   *        as the deciding factor between using RRT* or SBA* samples.
   */
  template <typename SBAStarBundle,
            typename RandomSampler>
  inline void generate_anytime_lazy_sbarrtstar(const SBAStarBundle& bdl, 
                                               RandomSampler get_sample, 
                                               double init_relaxation,
                                               double SA_init_temperature = 1.0) {
    BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<typename SBAStarBundle::visitor_type,typename SBAStarBundle::graph_type,typename SBAStarBundle::topology_type>));
    
    detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_distance, bdl.m_predecessor, bdl.m_key);
    
    generate_anytime_lazy_sbarrtstar_no_init(bdl, get_sample, init_relaxation, SA_init_temperature);
    
  };
  
  
  

};

};

#endif

















/**
 * \file sbastar_rrtstar.hpp
 *
 * This library provides function templates and concepts that implement a Sampling-based A* search
 * algorithm with an RRT* exploratory phase. A SBA* uses the A* search algorithm to drive the expansion 
 * of a roadmap into the free-space in order to connect a start and goal location. When useful nodes
 * are exhausted, an number of RRT* iterations are performed to use the Voronoi bias to generate 
 * more useful nodes before continuing the SBA* iterations. This algorithm has many customization points 
 * because there are many choices to be made in the method, such as how to find nearest neighbors for 
 * attempting to connect them through free-space, how to expand vertices, when to stop the algorithm, etc. 
 * All these customization points are left to the user to implement, some are defined by the 
 * SBARRTStarVisitorConcept (random-walk, vertex-added, etc.).
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

#ifndef REAK_SBASTAR_RRTSTAR_HPP
#define REAK_SBASTAR_RRTSTAR_HPP

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
#include "node_generators.hpp"

#include "lazy_sbastar.hpp"

#include <stack>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Graph */
namespace graph {

  /**
   * This concept class defines the valid expressions required of a class to be used as a visitor 
   * class for the SBA*-RRT* algorithm. A visitor class is essentially a class that regroups a number of 
   * callback functions that can be used to inject customization into the SBA*-RRT* algorithm. In other 
   * words, the visitor pattern in generic programming is an implementation of IoC 
   * (Inversion of Control), since the SBA*-RRT* algorithm is in control of execution, but custom behavior can
   * be injected in several places, even blocking the algorithm if needed.
   * 
   * Required concepts:
   * 
   * the visitor should model the SBAStarVisitorConcept.
   * 
   * Valid expressions:
   * 
   * tie(p,b,ep) = vis.steer_towards_position(p,u,g);  This function is called to attempt to steer from vertex u to position p, it returns a std::pair with the position that could be reached and a boolean value to indicate whether any significant motion occurred (collision-free).
   * 
   * \tparam Visitor The visitor class to be tested for modeling an AD* visitor concept.
   * \tparam Graph The graph type on which the visitor should be able to act.
   * \tparam Topology The topology type on which the visitor class is required to work with.
   */
  template <typename Visitor, typename Graph, typename Topology>
  struct SBARRTStarVisitorConcept : SBAStarVisitorConcept<Visitor, Graph, Topology> {
    BOOST_CONCEPT_USAGE(SBARRTStarVisitorConcept)
    {
      boost::tie(this->pt,this->b,this->ep) = this->vis.steer_towards_position(this->pt,this->u,this->g);
    };
  };
  
  /**
   * This class is simply a "null" visitor for the SBA*-RRT* algorithm. It is null in the sense that it
   * will do nothing on all accounts.
   */
  template <typename Topology>
  class default_sbarrtstar_visitor : public default_sbastar_visitor<Topology> {
    public:
      typedef typename ReaK::pp::topology_traits<Topology>::point_type PointType;
      
      default_sbarrtstar_visitor() {};
      
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
    struct sbarrtstar_bfs_visitor :
      sbastar_bfs_visitor<Topology, UniformCostVisitor, SBANodeConnector, UpdatableQueue, IndexInHeapMap, 
                          AStarHeuristicMap, PositionMap, WeightMap, DensityMap, ConstrictionMap, 
                          DistanceMap, PredecessorMap, KeyMap, NcSelector>
    {
      typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
      typedef sbastar_bfs_visitor<Topology, UniformCostVisitor, SBANodeConnector, UpdatableQueue, IndexInHeapMap, 
                                  AStarHeuristicMap, PositionMap, WeightMap, DensityMap, ConstrictionMap, 
                                  DistanceMap, PredecessorMap, KeyMap, NcSelector> base_type;

      sbarrtstar_bfs_visitor(const Topology& super_space, UniformCostVisitor vis, 
                             SBANodeConnector connector, RRTNodeGenerator generator,
                             UpdatableQueue& Q, IndexInHeapMap index_in_heap,  
                             AStarHeuristicMap heuristic, PositionMap pos, WeightMap weight, 
                             DensityMap density, ConstrictionMap constriction, DistanceMap dist, 
                             PredecessorMap pred, KeyMap key, NcSelector select_neighborhood) : 
                             base_type(super_space, vis, connector, Q, index_in_heap, 
                                       heuristic, pos, weight, density, constriction, 
                                       dist, pred, key, select_neighborhood),
                             m_node_generator(generator) { };
      
      template <typename Graph>
      void add_exploratory_node(Graph& g) const {
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename Graph::edge_bundled EdgeProp;
        
        typedef boost::composite_property_map< 
          PositionMap, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t > > GraphPositionMap;
        GraphPositionMap g_position = GraphPositionMap(this->m_position, boost::whole_bundle_property_map< Graph, boost::vertex_bundle_t >(&g));
        
        boost::tuple< Vertex, PositionValue, EdgeProp > gen_result = m_node_generator(g, this->m_vis, g_position);
        
        this->m_connect_vertex(get<1>(gen_result), get<0>(gen_result), get<2>(gen_result), 
                               g, this->m_super_space, *this, 
                               this->m_position, this->m_distance, this->m_predecessor, this->m_weight,
                               this->m_select_neighborhood);
        
      };
      
      RRTNodeGenerator m_node_generator;
    };
    
    
    template <typename Graph, //this is the actual graph, should comply to BidirectionalGraphConcept.
              typename Vertex, //this is the type to describe a vertex in the graph.
              typename SBARRTStarBFSVisitor, //this is a visitor class that can perform special operations at event points.
              typename MutableQueue>
    inline void
    sbarrtstar_search_loop
      (Graph &g, Vertex start_vertex, SBARRTStarBFSVisitor& sba_vis, MutableQueue& Q)
    { 
      while (sba_vis.keep_going()) {
        
        sba_vis.requeue_vertex(start_vertex,g);
        
        while (!Q.empty() && sba_vis.keep_going()) { 
          Vertex u = Q.top(); Q.pop();
          
          sba_vis.examine_vertex(u, g);
          
          // stop if the best node does not meet the potential threshold.
          if( ! sba_vis.has_search_potential(u, g) )
            break;
          
          // if the node still has a minimally good potential, then push it back on the OPEN queue.
          sba_vis.requeue_vertex(u,g);
          
        }; // end while  (the queue is either empty or it contains vertices that still have low key values.
        
        sba_vis.publish_path(g);
        
        std::size_t max_vertex_count = num_vertices(g);
        max_vertex_count += 4 * (math::highest_set_bit(max_vertex_count) + 1);
        while((num_vertices(g) < max_vertex_count) && (sba_vis.keep_going()))
          sba_vis.add_exploratory_node(g);
        
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
   * \tparam SBARRTStarVisitor The type of the SBA*-RRT* visitor to be used, should model the SBARRTStarVisitorConcept.
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
   * \param vis A SBA*-RRT* visitor implementing the SBARRTStarVisitorConcept. This is the 
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
  generate_sbarrtstar_no_init
    (Graph &g, Vertex start_vertex, const Topology& super_space, SBARRTStarVisitor vis,  // basic parameters
     AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
     DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
     PredecessorMap predecessor, KeyMap key,                          // properties resulting from the algorithm
     RandomSampler get_sample,
     NcSelector select_neighborhood)
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
    
    detail::sbarrtstar_bfs_visitor<
      Topology, 
      SBARRTStarVisitor,
      detail::sbastar_node_connector,
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
      NcSelector> bfs_vis(super_space, vis, detail::sbastar_node_connector(), 
                          rrg_node_generator<Topology, RandomSampler, NcSelector>(&super_space, get_sample, select_neighborhood), 
                          Q, index_in_heap, 
                          hval, position, weight, 
                          density, constriction, distance,
                          predecessor, key, select_neighborhood);
    
    detail::sbastar_search_loop(g, start_vertex, bfs_vis, Q);
    
  };


   /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam SBARRTStarVisitor The type of the SBA*-RRT* visitor to be used, should model the SBARRTStarVisitorConcept.
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
   * \param vis A SBA*-RRT* visitor implementing the SBARRTStarVisitorConcept. This is the 
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
  generate_sbarrtstar
    (Graph &g, Vertex start_vertex, const Topology& super_space, SBARRTStarVisitor vis,  // basic parameters
     AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
     DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
     PredecessorMap predecessor, KeyMap key,                         // properties resulting from the algorithm
     RandomSampler get_sample,
     NcSelector select_neighborhood)
  {
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    //BOOST_CONCEPT_ASSERT((boost::MutablePropertyGraphConcept<Graph>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::PointDistributionConcept<Topology>));
    BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<SBARRTStarVisitor,Graph,Topology>));
    
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
    typedef typename Graph::vertex_bundled VertexProp;
    
    for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
      put(distance, *ui, std::numeric_limits<double>::infinity());
      put(key, *ui, 0.0);
      put(predecessor, *ui, *ui);
      vis.initialize_vertex(*ui, g);
    };

    generate_sbarrtstar_no_init(
      g, start_vertex, super_space, vis, 
      hval, position, weight, density, constriction, distance,
      predecessor, key, get_sample, select_neighborhood);

  };
  
  
  
  
  
  
  /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the SBA* algorithm, without initialization of the existing graph.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam SBARRTStarVisitor The type of the SBA*-RRT* visitor to be used, should model the SBARRTStarVisitorConcept.
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
   * \param vis A SBA*-RRT* visitor implementing the SBARRTStarVisitorConcept. This is the 
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
  generate_lazy_sbarrtstar_no_init
    (Graph &g, Vertex start_vertex, const Topology& super_space, SBARRTStarVisitor vis,  // basic parameters
     AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
     DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
     PredecessorMap predecessor, KeyMap key,                          // properties resulting from the algorithm
     RandomSampler get_sample,
     NcSelector select_neighborhood)
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
    
    detail::sbarrtstar_bfs_visitor<
      Topology, 
      SBARRTStarVisitor,
      detail::lazy_sbastar_node_connector,
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
      NcSelector> bfs_vis(super_space, vis, detail::lazy_sbastar_node_connector(), 
                          rrg_node_generator<Topology, RandomSampler, NcSelector>(&super_space, get_sample, select_neighborhood), 
                          Q, index_in_heap, 
                          hval, position, weight, 
                          density, constriction, distance,
                          predecessor, key, select_neighborhood);
    
    detail::sbastar_search_loop(g, start_vertex, bfs_vis, Q);
    
  };


   /**
   * This function template generates a roadmap to connect a goal location to a start location
   * using the SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
   * \tparam Graph The graph type that can store the generated roadmap, should model 
   *         BidirectionalGraphConcept and MutableGraphConcept.
   * \tparam Vertex The type to describe a vertex of the graph on which the search is performed.
   * \tparam Topology The topology type that represents the free-space, should model BGL's Topology concept.
   * \tparam SBARRTStarVisitor The type of the SBA*-RRT* visitor to be used, should model the SBARRTStarVisitorConcept.
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
   * \param vis A SBA*-RRT* visitor implementing the SBARRTStarVisitorConcept. This is the 
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
  generate_lazy_sbarrtstar
    (Graph &g, Vertex start_vertex, const Topology& super_space, SBARRTStarVisitor vis,  // basic parameters
     AStarHeuristicMap hval, PositionMap position, WeightMap weight,                 // properties provided by the caller.
     DensityMap density, ConstrictionMap constriction, DistanceMap distance,       // properties needed by the algorithm, filled by the visitor.
     PredecessorMap predecessor, KeyMap key,                         // properties resulting from the algorithm
     RandomSampler get_sample,
     NcSelector select_neighborhood)
  {
    BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
    //BOOST_CONCEPT_ASSERT((boost::MutablePropertyGraphConcept<Graph>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((ReaK::pp::PointDistributionConcept<Topology>));
    BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<SBARRTStarVisitor,Graph,Topology>));
    
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
    typedef typename Graph::vertex_bundled VertexProp;
    
    for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
      put(distance, *ui, std::numeric_limits<double>::infinity());
      put(key, *ui, 0.0);
      put(predecessor, *ui, *ui);
      vis.initialize_vertex(*ui, g);
    };

    generate_lazy_sbarrtstar_no_init(
      g, start_vertex, super_space, vis, 
      hval, position, weight, density, constriction, distance,
      predecessor, key, get_sample, select_neighborhood);

  };
  
  
  
  
  

};

};

#endif

















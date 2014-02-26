/**
 * \file sbastar_bidir.hpp
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
 * \date February 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_SBASTAR_BIDIR_HPP
#define REAK_SBASTAR_BIDIR_HPP

#include <utility>
#include <boost/tuple/tuple.hpp>

#include "path_planning/metric_space_concept.hpp"

#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/detail/d_ary_heap.hpp>

// BGL-Extra includes:
#include <boost/graph/more_property_tags.hpp>
#include <boost/graph/more_property_maps.hpp>

#include "sbastar_search.hpp"
#include "pruned_connector.hpp"
#include "sbmp_visitor_concepts.hpp"


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Graph */
namespace graph {

template <typename Visitor, typename Graph, typename Topology>
struct SBAStarBidirVisitorConcept : SBAStarVisitorConcept<Visitor,Graph,Topology> {
  
  BOOST_CONCEPT_ASSERT((NodeBackPushingVisitorConcept<Visitor,Graph,Topology>));
  BOOST_CONCEPT_ASSERT((NodeBackPullingVisitorConcept<Visitor,Graph,Topology>));
  
  BOOST_CONCEPT_USAGE(SBAStarBidirVisitorConcept)
  { }
};

template <typename Topology>
struct sbastar_bidir_visitor_archetype : 
  sbastar_visitor_archetype, 
  node_back_pulling_visitor_archetype,
  node_back_pushing_visitor_archetype { };


namespace detail {
                  
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
  struct sbastar_bidir_bfs_visitor
  {

    sbastar_bidir_bfs_visitor(
      UniformCostVisitor vis, UpdatableQueue& Q, IndexInHeapMap index_in_heap,  
      PositionMap pos, WeightMap weight, 
      DensityMap density, ConstrictionMap constriction, 
      DistanceMap dist, PredecessorMap pred, 
      FwdDistanceMap fwd_dist, SuccessorMap succ, KeyMap key) : 
      m_vis(vis), m_Q(Q), m_index_in_heap(index_in_heap), 
      m_position(pos), m_weight(weight),
      m_density(density), m_constriction(constriction), 
      m_distance(dist), m_predecessor(pred), 
      m_fwd_distance(fwd_dist), m_successor(succ), m_key(key) { };
    
    typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
    
    template <class Graph>
    typename boost::graph_traits<Graph>::vertex_descriptor create_vertex(const PositionValue& p, Graph& g) const {
      typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
      typedef typename Graph::vertex_bundled VertexProp;
      
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
      put(m_index_in_heap, u, static_cast<std::size_t>(-1));
      put(m_key, u, 0.0);
      
      return u;
    };
    
    template <typename Vertex, typename Graph>
    boost::tuple<PositionValue, bool, typename Graph::edge_bundled> steer_towards_position(const PositionValue& p, Vertex u, const Graph& g) const { 
      return m_vis.steer_towards_position(p, u, g);
    };
    
    template <typename Vertex, typename Graph>
    boost::tuple<PositionValue, bool, typename Graph::edge_bundled> steer_back_to_position(const PositionValue& p, Vertex u, const Graph& g) const { 
      return m_vis.steer_back_to_position(p, u, g);
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
    boost::tuple< PositionValue, bool, typename Graph::edge_bundled > random_back_walk(Vertex u, Graph& g) const {
      return m_vis.random_back_walk(u, g);
    };
    
    template <typename Vertex, typename Graph>
    void vertex_to_be_removed(Vertex u, Graph& g) const {
      put(m_key, u, -std::numeric_limits<double>::infinity());
      m_Q.push_or_update(u);
      m_Q.pop();
      m_vis.vertex_to_be_removed(u,g);
    };
    
    template <typename Vertex, typename Graph>
    void vertex_added(Vertex v, Graph& g) const { m_vis.vertex_added(v,g); };
    
    template <typename Edge, typename Graph>
    void edge_added(Edge e, Graph& g) const {
      m_vis.edge_added(e,g);
      m_vis.examine_edge(e,g);
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
    
    template <class Vertex, class Graph>
    void examine_vertex(Vertex u, Graph& g) const {
      m_vis.examine_vertex(u, g);
    };
    
    template <class Edge, class Graph>
    void examine_edge(Edge e, Graph& g) const {
      m_vis.examine_edge(e, g);
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
      m_vis.affected_vertex(u,g);
      double g_u = get(m_distance, g[u]);
      double h_u = get(m_fwd_distance, g[u]);
      // Key-value for the min-heap (priority-queue):
      put(m_key, u, (g_u + h_u) / (1.0 - get(m_constriction, g[u])) / (1.0 - get(m_density, g[u])));
    };

    UniformCostVisitor m_vis;
    UpdatableQueue& m_Q; 
    IndexInHeapMap m_index_in_heap;
    
    PositionMap m_position;
    WeightMap m_weight;
    DensityMap m_density;
    ConstrictionMap m_constriction; 
    DistanceMap m_distance;
    PredecessorMap m_predecessor; 
    DistanceMap m_fwd_distance;
    PredecessorMap m_successor;
    KeyMap m_key;
    
  };
  
  
  struct sba_bidir_node_generator {
    
    template <typename Graph,
              typename SBAVisitor,
              typename PositionMap>
    boost::tuple< typename boost::graph_traits<Graph>::vertex_descriptor,
                  typename boost::property_traits<PositionMap>::value_type,
                  typename Graph::edge_bundled,
                  typename boost::graph_traits<Graph>::vertex_descriptor,
                  typename boost::property_traits<PositionMap>::value_type,
                  typename Graph::edge_bundled > 
      operator()(typename boost::graph_traits<Graph>::vertex_descriptor u, Graph& g,
                 const SBAVisitor& sba_vis, PositionMap) const {
      typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
      typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
      typedef typename Graph::edge_bundled EdgeProp;
      typedef boost::tuple< Vertex, PositionValue, EdgeProp, Vertex, PositionValue, EdgeProp > ResultType;
      
      PositionValue p_exp, p_ret; 
      bool was_expanded, was_retracted; 
      EdgeProp ep_exp, ep_ret;
      if( get(sba_vis.m_predecessor, g[u]) != boost::graph_traits<Graph>::null_vertex() )
        boost::tie(p_exp, was_expanded, ep_exp) = sba_vis.random_walk(u, g);
      if( get(sba_vis.m_successor, g[u]) != boost::graph_traits<Graph>::null_vertex() )
        boost::tie(p_ret, was_retracted, ep_ret) = sba_vis.random_back_walk(u, g);
      
      return ResultType(( was_expanded ? u : boost::graph_traits<Graph>::null_vertex()), p_exp, ep_exp, 
                        ( was_retracted ? u : boost::graph_traits<Graph>::null_vertex()), p_ret, ep_ret);
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
  void sbastar_bidir_loop(
      Graph &g, const Topology& super_space, SBAStarVisitor& sba_vis, 
      MotionGraphConnector connect_vertex, SBANodeGenerator sba_generate_node,
      MutableQueue& Q, NcSelector select_neighborhood)
  { 
    typedef typename ReaK::pp::topology_traits<Topology>::point_type PositionValue;
    typedef typename Graph::edge_bundled EdgeProp;
    
    while (!Q.empty() && sba_vis.keep_going()) { 
      Vertex u = Q.top(); Q.pop();
    
      // stop if the best node does not meet the potential threshold.
      if( ! sba_vis.has_search_potential(u, g) )
        break;
      
      sba_vis.examine_vertex(u, g);
      
      Vertex x_near_pred, x_near_succ; 
      PositionValue p_new_pred, p_new_succ; 
      EdgeProp ep_pred, ep_succ;
      boost::tie(x_near_pred, p_new_pred, ep_pred, x_near_succ, p_new_succ, ep_succ) = 
        sba_generate_node(u, g, sba_vis, sba_vis.m_position);
      
      // then push it back on the OPEN queue.
      sba_vis.requeue_vertex(u,g);
      
      if( x_near_pred != boost::graph_traits<Graph>::null_vertex() ) {
        Vertex x_near_other = boost::graph_traits<Graph>::null_vertex();
        EdgeProp ep_other;
        connect_vertex(p_new_pred, x_near_pred, ep_pred, x_near_other, ep_other, g, 
                       super_space, sba_vis, sba_vis.m_position, 
                       sba_vis.m_distance, sba_vis.m_predecessor, 
                       sba_vis.m_fwd_distance, sba_vis.m_successor, 
                       sba_vis.m_weight, select_neighborhood);
      };
      if( x_near_succ != boost::graph_traits<Graph>::null_vertex() ) {
        Vertex x_near_other = boost::graph_traits<Graph>::null_vertex();
        EdgeProp ep_other;
        connect_vertex(p_new_succ, x_near_other, ep_other, x_near_succ, ep_succ, g, 
                       super_space, sba_vis, sba_vis.m_position, 
                       sba_vis.m_distance, sba_vis.m_predecessor, 
                       sba_vis.m_fwd_distance, sba_vis.m_successor, 
                       sba_vis.m_weight, select_neighborhood);
      };
      
    }; // end while 
    
  };
  
  template <typename Graph,
            typename Vertex,
            typename Topology,
            typename SBAStarVisitor,
            typename NodeConnector,
            typename PositionMap,
            typename WeightMap,
            typename DensityMap,
            typename ConstrictionMap,
            typename DistanceMap,
            typename PredecessorMap,
            typename FwdDistanceMap,
            typename SuccessorMap,
            typename KeyMap,
            typename NcSelector>
  void generate_sbastar_bidir_no_init_impl(
      Graph &g, Vertex start_vertex, Vertex goal_vertex, 
      const Topology& super_space, SBAStarVisitor vis, 
      NodeConnector connect_vertex, PositionMap position, WeightMap weight, 
      DensityMap density, ConstrictionMap constriction, 
      DistanceMap distance, PredecessorMap predecessor, 
      FwdDistanceMap fwd_distance, SuccessorMap successor, 
      KeyMap key, NcSelector select_neighborhood)
  {
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
      SBAStarVisitor, MutableQueue, IndexInHeapMap, 
      PositionMap, WeightMap, DensityMap, ConstrictionMap, 
      DistanceMap, PredecessorMap, FwdDistanceMap, SuccessorMap,
      KeyMap> sba_bfs_vis(vis, Q, index_in_heap, position, weight, 
                          density, constriction, distance, predecessor, 
                          fwd_distance, successor, key);
    
    put(distance, g[start_vertex], 0.0);
    put(predecessor, g[start_vertex], start_vertex);
    put(fwd_distance, g[goal_vertex], 0.0);
    put(successor, g[goal_vertex], goal_vertex);
    
    sbastar_bidir_loop(g, super_space, sba_bfs_vis, 
                       connect_vertex, sba_node_generator(), 
                       Q, select_neighborhood);
  };
  
  
  template <typename Graph,
            typename SBAStarVisitor,
            typename DistanceMap,
            typename PredecessorMap,
            typename KeyMap>
  void initialize_sbastar_bidir_nodes(
      Graph &g, SBAStarVisitor vis, 
      DistanceMap distance, PredecessorMap predecessor, 
      FwdDistanceMap fwd_distance, SuccessorMap successor, KeyMap key) {
    typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
    for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
      put(distance, g[*ui], std::numeric_limits<double>::infinity());
      put(fwd_distance, g[*ui], std::numeric_limits<double>::infinity());
      put(predecessor, g[*ui], boost::graph_traits<Graph>::null_vertex());
      put(successor, g[*ui], boost::graph_traits<Graph>::null_vertex());
      put(key, *ui, 0.0);
      vis.initialize_vertex(*ui, g);
    };
  };

}; //end of detail namespace.


template <typename Graph, typename Vertex, typename Topology, typename SBAStarVisitor,
          typename PositionMap, typename WeightMap, typename DensityMap, typename ConstrictionMap,
          typename DistanceMap, typename PredecessorMap, typename FwdDistanceMap, typename SuccessorMap,
          typename KeyMap, typename NcSelector>
struct sbastar_bidir_bundle {
  typedef Graph             graph_type;
  typedef Vertex            vertex_type;
  typedef Topology          topology_type;
  typedef SBAStarVisitor    visitor_type;
  typedef PositionMap       position_map_type;
  typedef WeightMap         weight_map_type;
  typedef DensityMap        density_map_type;
  typedef ConstrictionMap   constriction_map_type;
  typedef DistanceMap       distance_map_type;
  typedef PredecessorMap    predecessor_map_type;
  typedef FwdDistanceMap    fwd_distance_map_type;
  typedef SuccessorMap      successor_map_type;
  typedef KeyMap            key_map_type;
  typedef NcSelector        nc_selector_type;
  
  graph_type* m_g;
  vertex_type m_start_vertex;
  vertex_type m_goal_vertex;
  const topology_type* m_super_space;
  visitor_type m_vis;
  position_map_type m_position;
  weight_map_type m_weight;
  density_map_type m_density;
  constriction_map_type m_constriction;
  distance_map_type m_distance;
  predecessor_map_type m_predecessor;
  fwd_distance_map_type m_fwd_distance;
  successor_map_type m_successor;
  key_map_type m_key;
  nc_selector_type m_select_neighborhood;
  
  sbastar_bidir_bundle(
    Graph &g, Vertex start_vertex, Vertex goal_vertex, const Topology& super_space, SBAStarVisitor vis,
    PositionMap position, WeightMap weight,
    DensityMap density, ConstrictionMap constriction, DistanceMap distance,
    PredecessorMap predecessor, FwdDistanceMap fwd_distance, SuccessorMap successor, 
    KeyMap key, NcSelector select_neighborhood) :
    m_g(&g), m_start_vertex(start_vertex), m_goal_vertex(goal_vertex), 
    m_super_space(&super_space), m_vis(vis),
    m_position(position), m_weight(weight), 
    m_density(density), m_constriction(constriction), m_distance(distance), 
    m_predecessor(predecessor), m_fwd_distance(fwd_distance), m_successor(successor), 
    m_key(key), m_select_neighborhood(select_neighborhood) { };
  
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
template <typename Graph, typename Vertex, typename Topology, typename SBAStarVisitor,
          typename PositionMap, typename WeightMap, typename DensityMap, typename ConstrictionMap,
          typename DistanceMap, typename PredecessorMap, typename FwdDistanceMap, typename SuccessorMap,
          typename KeyMap, typename NcSelector>
sbastar_bidir_bundle<Graph, Vertex, Topology, SBAStarVisitor, 
                     PositionMap, WeightMap, 
                     DensityMap, ConstrictionMap, DistanceMap, 
                     PredecessorMap, FwdDistanceMap, SuccessorMap, KeyMap, NcSelector>
  make_sbastar_bidir_bundle(
    Graph &g, Vertex start_vertex, Vertex goal_vertex, const Topology& super_space, SBAStarVisitor vis, 
    PositionMap position, WeightMap weight, 
    DensityMap density, ConstrictionMap constriction, DistanceMap distance, 
    PredecessorMap predecessor, FwdDistanceMap fwd_distance, SuccessorMap successor, 
    KeyMap key, NcSelector select_neighborhood) {
  
  BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((SBAStarBidirVisitorConcept<SBAStarVisitor,Graph,Topology>));
  
  return sbastar_bidir_bundle<Graph, Vertex, Topology, SBAStarVisitor, 
                              PositionMap, WeightMap, DensityMap, ConstrictionMap, 
                              DistanceMap, PredecessorMap, FwdDistanceMap, SuccessorMap, 
                              KeyMap, NcSelector>(
                                g, start_vertex, goal_vertex, super_space, vis, 
                                position, weight, density, constriction, 
                                distance, predecessor, fwd_distance, successor, 
                                key, select_neighborhood);
};



/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the SBA* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bidir_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bidir_bundle().
  */
template <typename SBAStarBundle>
void generate_sbastar_bidir_no_init(const SBAStarBundle& bdl) {
  
  detail::generate_sbastar_bidir_no_init_impl(
    *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, 
    *(bdl.m_super_space), bdl.m_vis, pruned_node_connector(),
    bdl.m_position, bdl.m_weight, bdl.m_density, bdl.m_constriction, 
    bdl.m_distance, bdl.m_predecessor, bdl.m_fwd_distance, bdl.m_successor, 
    bdl.m_key, bdl.m_select_neighborhood);
  
};

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the SBA* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  */
template <typename SBAStarBundle>
void generate_sbastar_bidir(const SBAStarBundle& bdl) {
  
  detail::initialize_sbastar_bidir_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_distance, bdl.m_predecessor, 
                                         bdl.m_fwd_distance, bdl.m_successor, bdl.m_key);
  
  generate_sbastar_bidir_no_init(bdl);
  
};



};

};

#endif




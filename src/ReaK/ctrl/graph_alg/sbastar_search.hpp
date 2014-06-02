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

#include <utility>
#include <boost/tuple/tuple.hpp>

#include <ReaK/ctrl/path_planning/metric_space_concept.hpp>

#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/detail/d_ary_heap.hpp>

// BGL-Extra includes:
#include <boost/graph/more_property_tags.hpp>
#include <boost/graph/more_property_maps.hpp>

#include "pruned_connector.hpp"
#include "lazy_connector.hpp"
#include "branch_and_bound_connector.hpp"
#include "sbmp_visitor_concepts.hpp"


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
  * The visitor class should model SBMPVisitorConcept, NodePushingVisitorConcept, NodeReConnectVisitorConcept,
  * NeighborhoodTrackingVisitorConcept, and NodeExploringVisitorConcept.
  * 
  * Valid expressions:
  * 
  * vis.publish_path(g);  A function to notify the visitor that at least one A* round has completed and its resulting path (partial or complete) can be published (the path is encoded in the predecessor property-map).
  * 
  * \tparam Visitor The visitor class to be tested for modeling an AD* visitor concept.
  * \tparam Graph The graph type on which the visitor should be able to act.
  * \tparam Topology The topology type on which the visitor class is required to work with.
  */
template <typename Visitor, typename Graph, typename Topology>
struct SBAStarVisitorConcept : SBMPVisitorConcept<Visitor,Graph> {
  
  BOOST_CONCEPT_ASSERT((NodePushingVisitorConcept<Visitor,Graph,Topology>));
  BOOST_CONCEPT_ASSERT((NodeReConnectVisitorConcept<Visitor,Graph>));
  BOOST_CONCEPT_ASSERT((NeighborhoodTrackingVisitorConcept<Visitor,Graph>));
  BOOST_CONCEPT_ASSERT((NodeExploringVisitorConcept<Visitor,Graph>));
  
  BOOST_CONCEPT_USAGE(SBAStarVisitorConcept)
  {
    vis.publish_path(g);
  };
  Visitor vis;
  Graph g;
};


/**
 * This class is simply an archetype visitor for the SBA* algorithm.
 */
template <typename Topology>
struct sbastar_visitor_archetype : 
  sbmp_visitor_archetype, 
  node_pushing_visitor_archetype<Topology>, 
  node_reconnect_visitor_archetype, 
  neighborhood_tracking_visitor_archetype,
  node_exploring_visitor_archetype { 
  
  template <typename Graph>
  void publish_path(const Graph&) const { };
};


/**
  * This concept class defines the valid expressions required of a class to be used as a visitor 
  * class for the bi-directional SBA* algorithm. A visitor class is essentially a class that regroups a number of 
  * callback functions that can be used to inject customization into the bi-directional SBA* algorithm. 
  * In other words, the visitor pattern in generic programming is an implementation of IoC 
  * (Inversion of Control), since the bi-directional SBA* algorithm is in control of execution, 
  * but custom behavior can be injected in several places, even blocking the algorithm if needed.
  * 
  * Required concepts:
  * 
  * The visitor class should model SBAStarVisitorConcept and NodeBackPushingVisitorConcept.
  * 
  * \tparam Visitor The visitor class to be tested for modeling an AD* visitor concept.
  * \tparam Graph The graph type on which the visitor should be able to act.
  * \tparam Topology The topology type on which the visitor class is required to work with.
  */
template <typename Visitor, typename Graph, typename Topology>
struct SBAStarBidirVisitorConcept : SBAStarVisitorConcept<Visitor,Graph,Topology> {
  
  BOOST_CONCEPT_ASSERT((NodeBackPushingVisitorConcept<Visitor,Graph,Topology>));
  
  BOOST_CONCEPT_USAGE(SBAStarBidirVisitorConcept)
  { };
};

/**
 * This class is simply an archetype visitor for the bi-directional SBA* algorithm.
 */
template <typename Topology>
struct sbastar_bidir_visitor_archetype : 
  sbastar_visitor_archetype<Topology>, 
  node_back_pushing_visitor_archetype<Topology> { };


namespace detail {
                  
  template <typename Graph, typename UniformCostVisitor,
            typename UpdatableQueue, typename IndexInHeapMap,
            typename KeyMap, typename PositionMap, typename WeightMap,
            typename DensityMap, typename ConstrictionMap, 
            typename DistanceMap, typename PredecessorMap,
            typename FwdDistanceMap, typename SuccessorMap = null_vertex_prop_map<Graph> >
  struct sbastar_bfs_visitor
  {
    
    sbastar_bfs_visitor(UniformCostVisitor vis, UpdatableQueue& Q, IndexInHeapMap index_in_heap,  
                        KeyMap key, PositionMap pos, WeightMap weight, 
                        DensityMap density, ConstrictionMap constriction, 
                        DistanceMap dist, PredecessorMap pred, 
                        FwdDistanceMap fwd_dist, SuccessorMap succ = SuccessorMap()) : 
                        m_vis(vis), m_Q(Q), m_index_in_heap(index_in_heap), 
                        m_key(key), m_position(pos), m_weight(weight),
                        m_density(density), m_constriction(constriction), 
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
      put(m_index_in_heap, u, static_cast<std::size_t>(-1));
      put(m_key, u, 0.0);
      
      return u;
    };
    
    boost::tuple<PositionValue, bool, EdgeProp> steer_towards_position(const PositionValue& p, Vertex u, const Graph& g) const { 
      return m_vis.steer_towards_position(p, u, g);
    };
    
    boost::tuple<PositionValue, bool, EdgeProp> steer_back_to_position(const PositionValue& p, Vertex u, const Graph& g) const { 
      return this->m_vis.steer_back_to_position(p, u, g);
    };
    
    std::pair< bool, EdgeProp > can_be_connected(Vertex u, Vertex v, Graph& g) const {
      return m_vis.can_be_connected(u, v, g);
    };
    
    boost::tuple< PositionValue, bool, EdgeProp > random_walk(Vertex u, Graph& g) const {
      return m_vis.random_walk(u, g);
    };
    
    boost::tuple< PositionValue, bool, EdgeProp > random_back_walk(Vertex u, Graph& g) const {
      return this->m_vis.random_back_walk(u, g);
    };
    
    void vertex_to_be_removed(Vertex u, Graph& g) const {
      put(m_key, u, -std::numeric_limits<double>::infinity());
      m_Q.push_or_update(u);
      m_Q.pop();
      m_vis.vertex_to_be_removed(u,g);
    };
    
    void vertex_added(Vertex v, Graph& g) const { m_vis.vertex_added(v,g); };
    
    void edge_added(Edge e, Graph& g) const {
      m_vis.edge_added(e,g);
      m_vis.examine_edge(e,g);
    };
    
    void travel_explored(Vertex u, Vertex v, Graph& g) const {
      m_vis.travel_explored(u, v, g);
    };
    
    void travel_succeeded(Vertex u, Vertex v, Graph& g) const {
      m_vis.travel_succeeded(u, v, g);
    };
    
    void travel_failed(Vertex u, Vertex v, Graph& g) const {
      m_vis.travel_failed(u, v, g);
    };
    
    void requeue_vertex(Vertex u, Graph& g) const {
      update_key(u,g);
      if( ! m_vis.should_close(u, g) ) {
        m_Q.push_or_update(u);
        m_vis.discover_vertex(u, g);
      };
    };
    void affected_vertex(Vertex u, Graph& g) const { requeue_vertex(u,g); }; // same function, different name.
    
    void examine_vertex(Vertex u, Graph& g) const {
      m_vis.examine_vertex(u, g);
    };
    
    void examine_edge(Edge e, Graph& g) const {
      m_vis.examine_edge(e, g);
    };
    
    void publish_path(const Graph& g) const { 
      m_vis.publish_path(g);
    };
    bool keep_going() const { 
      return m_vis.keep_going(); 
    };
    bool has_search_potential(Vertex u, const Graph& g) const { 
      return m_vis.has_search_potential(u,g);
    };
    bool should_close(Vertex u, const Graph& g) const { 
      return m_vis.should_close(u,g);
    };
    
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
    
    KeyMap m_key;
    PositionMap m_position;
    WeightMap m_weight;
    DensityMap m_density;
    ConstrictionMap m_constriction; 
    DistanceMap m_distance;
    PredecessorMap m_predecessor;
    FwdDistanceMap m_fwd_distance;
    SuccessorMap m_successor;
    
  };
  
  
  
  
  
  struct sba_node_generator {
    
    template <typename Graph, typename SBAVisitor, typename PositionMap>
    boost::tuple< typename boost::graph_traits<Graph>::vertex_descriptor,
                  typename boost::property_traits<PositionMap>::value_type,
                  typename Graph::edge_bundled > 
      operator()(typename boost::graph_traits<Graph>::vertex_descriptor u, Graph& g, 
                  const SBAVisitor& sba_vis, PositionMap) const {
      typedef typename boost::property_traits<PositionMap>::value_type PositionValue;
      typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
      typedef typename Graph::edge_bundled EdgeProp;
      typedef boost::tuple< Vertex, PositionValue, EdgeProp > ResultType;
      
      PositionValue p_new; bool walk_succeeded = false; EdgeProp ep_new;
      boost::tie(p_new, walk_succeeded, ep_new) = sba_vis.random_walk(u, g);
      if( walk_succeeded )
        return ResultType(u, p_new, ep_new);
      
      return ResultType(boost::graph_traits<Graph>::null_vertex(), PositionValue(), EdgeProp());
    };
    
  };
  
  struct sba_bidir_node_generator {
    
    template <typename Graph, typename SBAVisitor, typename PositionMap>
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
      bool was_expanded = false;
      bool was_retracted = false; 
      EdgeProp ep_exp, ep_ret;
      if( get(sba_vis.m_predecessor, g[u]) != boost::graph_traits<Graph>::null_vertex() )
        boost::tie(p_exp, was_expanded, ep_exp) = sba_vis.random_walk(u, g);
      if( get(sba_vis.m_successor, g[u]) != boost::graph_traits<Graph>::null_vertex() )
        boost::tie(p_ret, was_retracted, ep_ret) = sba_vis.random_back_walk(u, g);
      
      return ResultType(( was_expanded ? u : boost::graph_traits<Graph>::null_vertex()), p_exp, ep_exp, 
                        ( was_retracted ? u : boost::graph_traits<Graph>::null_vertex()), p_ret, ep_ret);
    };
    
  };
  
  
  
  template <typename Graph, typename Topology, typename SBAStarVisitor, typename MotionGraphConnector, 
            typename SBANodeGenerator, typename MutableQueue, typename NcSelector>
  void sbastar_search_loop(
      Graph &g, const Topology& super_space, SBAStarVisitor& sba_vis, 
      MotionGraphConnector connect_vertex, SBANodeGenerator sba_generate_node,
      MutableQueue& Q, NcSelector select_neighborhood)
  { 
    typedef typename ReaK::pp::topology_traits<Topology>::point_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    
    while (!Q.empty() && sba_vis.keep_going()) { 
      Vertex u = Q.top(); Q.pop();
      
      // stop if the best node does not meet the potential threshold.
      if( ! sba_vis.has_search_potential(u, g) )
        break;
      
      sba_vis.examine_vertex(u, g);
      
      PositionValue p_new; Vertex x_near; EdgeProp eprop;
      boost::tie(x_near, p_new, eprop) = sba_generate_node(u, g, sba_vis, sba_vis.m_position);
      
      // then push it back on the OPEN queue.
      sba_vis.requeue_vertex(u,g);
      
      if(x_near != boost::graph_traits<Graph>::null_vertex()) {
        connect_vertex(p_new, x_near, eprop, g, 
                        super_space, sba_vis, sba_vis.m_position, 
                        sba_vis.m_distance, sba_vis.m_predecessor, 
                        sba_vis.m_weight, select_neighborhood);
      };
      
    }; // end while
    
  };
  
  template <typename Graph, typename Topology, typename SBAStarVisitor, typename MotionGraphConnector, 
            typename SBANodeGenerator, typename MutableQueue, typename NcSelector>
  void sbastar_bidir_loop(
      Graph &g, const Topology& super_space, SBAStarVisitor& sba_vis, 
      MotionGraphConnector connect_vertex, SBANodeGenerator sba_generate_node,
      MutableQueue& Q, NcSelector select_neighborhood)
  { 
    typedef typename ReaK::pp::topology_traits<Topology>::point_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
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
  
  
  
  template <typename Graph, typename Topology, typename SBAStarVisitor,
            typename NodeConnector, typename KeyMap, typename PositionMap, typename WeightMap,
            typename DensityMap, typename ConstrictionMap, typename DistanceMap, typename PredecessorMap,
            typename FwdDistanceMap, typename NcSelector>
  void generate_sbastar_no_init_impl(Graph &g, 
      typename boost::graph_traits<Graph>::vertex_descriptor start_vertex, 
      const Topology& super_space, SBAStarVisitor vis, 
      NodeConnector connect_vertex, KeyMap key, PositionMap position, WeightMap weight, 
      DensityMap density, ConstrictionMap constriction, 
      DistanceMap distance, PredecessorMap predecessor, 
      FwdDistanceMap fwd_distance, NcSelector select_neighborhood)
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    
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
    
    sbastar_bfs_visitor<Graph, SBAStarVisitor, MutableQueue, IndexInHeapMap,
                        KeyMap, PositionMap, WeightMap, DensityMap, ConstrictionMap, 
                        DistanceMap, PredecessorMap, FwdDistanceMap> sba_bfs_vis(
      vis, Q, index_in_heap, key, position, weight, 
      density, constriction, distance, predecessor, fwd_distance);
    
    put(distance, g[start_vertex], 0.0);
    put(predecessor, g[start_vertex], start_vertex);
    sba_bfs_vis.requeue_vertex(start_vertex,g);
    
    sbastar_search_loop(g, super_space, sba_bfs_vis, connect_vertex, sba_node_generator(), Q, select_neighborhood);
  };
  
  template <typename Graph, typename Topology, typename SBAStarVisitor,
            typename NodeConnector, typename KeyMap, typename PositionMap, typename WeightMap,
            typename DensityMap, typename ConstrictionMap, typename DistanceMap, typename PredecessorMap, 
            typename FwdDistanceMap, typename SuccessorMap, typename NcSelector>
  void generate_sbastar_bidir_no_init_impl(Graph &g, 
      typename boost::graph_traits<Graph>::vertex_descriptor start_vertex, 
      typename boost::graph_traits<Graph>::vertex_descriptor goal_vertex, 
      const Topology& super_space, SBAStarVisitor vis, 
      NodeConnector connect_vertex, KeyMap key, PositionMap position, WeightMap weight, 
      DensityMap density, ConstrictionMap constriction, 
      DistanceMap distance, PredecessorMap predecessor, 
      FwdDistanceMap fwd_distance, SuccessorMap successor, 
      NcSelector select_neighborhood)
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    
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
    
    sbastar_bfs_visitor<Graph, SBAStarVisitor, MutableQueue, IndexInHeapMap,
                        KeyMap, PositionMap, WeightMap, DensityMap, ConstrictionMap, 
                        DistanceMap, PredecessorMap, FwdDistanceMap, SuccessorMap> sba_bfs_vis(
      vis, Q, index_in_heap, key, position, weight, 
      density, constriction, distance, predecessor, fwd_distance, successor);
    
    put(distance, g[start_vertex], 0.0);
    put(predecessor, g[start_vertex], start_vertex);
    sba_bfs_vis.requeue_vertex(start_vertex,g);
    if( goal_vertex != boost::graph_traits<Graph>::null_vertex() ) {
      put(fwd_distance, g[goal_vertex], 0.0);
      put(successor, g[goal_vertex], goal_vertex);
      sba_bfs_vis.requeue_vertex(goal_vertex,g);
    };
    
    sbastar_bidir_loop(g, super_space, sba_bfs_vis, connect_vertex, sba_bidir_node_generator(), Q, select_neighborhood);
  };
  
  
  template <typename Graph, typename SBAStarVisitor, typename KeyMap, 
            typename DistanceMap, typename PredecessorMap>
  void initialize_sbastar_nodes(Graph &g, SBAStarVisitor vis, KeyMap key, 
                                DistanceMap distance, PredecessorMap predecessor) {
    typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
    for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
      put(key, *ui, 0.0);
      put(distance, g[*ui], std::numeric_limits<double>::infinity());
      put(predecessor, g[*ui], boost::graph_traits<Graph>::null_vertex());
      vis.initialize_vertex(*ui, g);
    };
  };
  
  template <typename Graph, typename SBAStarVisitor, typename KeyMap, 
            typename DistanceMap, typename PredecessorMap, typename FwdDistanceMap, typename SuccessorMap>
  void initialize_sbastar_nodes(Graph &g, SBAStarVisitor vis, KeyMap key, 
                                DistanceMap distance, PredecessorMap predecessor, 
                                FwdDistanceMap fwd_distance, SuccessorMap successor) {
    typename boost::graph_traits<Graph>::vertex_iterator ui, ui_end;
    for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
      put(key, *ui, 0.0);
      put(distance, g[*ui], std::numeric_limits<double>::infinity());
      put(predecessor, g[*ui], boost::graph_traits<Graph>::null_vertex());
      put(fwd_distance, g[*ui], std::numeric_limits<double>::infinity());
      put(successor, g[*ui], boost::graph_traits<Graph>::null_vertex());
      vis.initialize_vertex(*ui, g);
    };
  };
  
}; //end of detail namespace.




template <typename Graph, typename Topology, typename SBAStarVisitor, typename NcSelector,
          typename KeyMap, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, 
          typename DistanceMap, typename PredecessorMap, 
          typename FwdDistanceMap, typename SuccessorMap = detail::null_vertex_prop_map<Graph> >
struct sbastar_bundle {
  typedef Graph             graph_type;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_type;
  typedef Topology          topology_type;
  typedef SBAStarVisitor    visitor_type;
  typedef NcSelector        nc_selector_type;
  typedef KeyMap            key_map_type;
  typedef PositionMap       position_map_type;
  typedef WeightMap         weight_map_type;
  typedef DensityMap        density_map_type;
  typedef ConstrictionMap   constriction_map_type;
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
  key_map_type m_key;
  position_map_type m_position;
  weight_map_type m_weight;
  density_map_type m_density;
  constriction_map_type m_constriction;
  distance_map_type m_distance;
  predecessor_map_type m_predecessor;
  fwd_distance_map_type m_fwd_distance;
  successor_map_type m_successor;
  
  sbastar_bundle(
    Graph &g, vertex_type start_vertex, const Topology& super_space, 
    SBAStarVisitor vis, NcSelector select_neighborhood,
    KeyMap key, PositionMap position, WeightMap weight,
    DensityMap density, ConstrictionMap constriction, 
    DistanceMap distance, PredecessorMap predecessor,
    FwdDistanceMap fwd_distance, SuccessorMap successor = SuccessorMap()) :
    m_g(&g), m_start_vertex(start_vertex), m_goal_vertex(boost::graph_traits<Graph>::null_vertex()), 
    m_super_space(&super_space), m_vis(vis), m_select_neighborhood(select_neighborhood),
    m_key(key), m_position(position), m_weight(weight), 
    m_density(density), m_constriction(constriction), 
    m_distance(distance), m_predecessor(predecessor), 
    m_fwd_distance(fwd_distance), m_successor(successor) { };
  
  sbastar_bundle(
    Graph &g, vertex_type start_vertex, vertex_type goal_vertex, const Topology& super_space, 
    SBAStarVisitor vis, NcSelector select_neighborhood,
    KeyMap key, PositionMap position, WeightMap weight,
    DensityMap density, ConstrictionMap constriction, 
    DistanceMap distance, PredecessorMap predecessor,
    FwdDistanceMap fwd_distance, SuccessorMap successor = SuccessorMap()) :
    m_g(&g), m_start_vertex(start_vertex), m_goal_vertex(goal_vertex), 
    m_super_space(&super_space), m_vis(vis), m_select_neighborhood(select_neighborhood),
    m_key(key), m_position(position), m_weight(weight), 
    m_density(density), m_constriction(constriction), 
    m_distance(distance), m_predecessor(predecessor), 
    m_fwd_distance(fwd_distance), m_successor(successor) { };
  
  
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
  * \param select_neighborhood A callable object (functor) that can select a list of 
  *        vertices of the graph that ought to be connected to a new 
  *        vertex. The list should be sorted in order of increasing "distance".
  * \param key The property-map which stores the AD* key-values associated to each vertex.
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
  * \param fwd_distance The property-map which stores the estimated distance of each vertex to the goal.
  */
template <typename Graph, typename Topology, typename SBAStarVisitor, typename NcSelector,
          typename KeyMap, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, 
          typename DistanceMap, typename PredecessorMap, typename FwdDistanceMap >
sbastar_bundle<Graph, Topology, SBAStarVisitor, NcSelector, KeyMap, PositionMap, WeightMap,
               DensityMap, ConstrictionMap, DistanceMap, PredecessorMap, FwdDistanceMap> 
  make_sbastar_bundle(Graph &g, 
    typename boost::graph_traits<Graph>::vertex_descriptor start_vertex, 
    const Topology& super_space, SBAStarVisitor vis, NcSelector select_neighborhood, 
    KeyMap key, PositionMap position, WeightMap weight, 
    DensityMap density, ConstrictionMap constriction, DistanceMap distance, 
    PredecessorMap predecessor, FwdDistanceMap fwd_distance) {
  
  BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((SBAStarVisitorConcept<SBAStarVisitor,Graph,Topology>));
  
  return sbastar_bundle<Graph, Topology, SBAStarVisitor, NcSelector, KeyMap, PositionMap, WeightMap,
                        DensityMap, ConstrictionMap, DistanceMap, PredecessorMap, FwdDistanceMap>(
    g, start_vertex, super_space, vis, select_neighborhood, key, position, weight,
    density, constriction, distance, predecessor, fwd_distance);
};

template <typename Graph, typename Topology, typename SBAStarVisitor, typename NcSelector,
          typename KeyMap, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, 
          typename DistanceMap, typename PredecessorMap, typename FwdDistanceMap >
sbastar_bundle<Graph, Topology, SBAStarVisitor, NcSelector, KeyMap, PositionMap, WeightMap,
               DensityMap, ConstrictionMap, DistanceMap, PredecessorMap, FwdDistanceMap> 
  make_sbastar_bundle(Graph &g, 
    typename boost::graph_traits<Graph>::vertex_descriptor start_vertex, 
    typename boost::graph_traits<Graph>::vertex_descriptor goal_vertex, 
    const Topology& super_space, SBAStarVisitor vis, NcSelector select_neighborhood, 
    KeyMap key, PositionMap position, WeightMap weight, 
    DensityMap density, ConstrictionMap constriction, DistanceMap distance, 
    PredecessorMap predecessor, FwdDistanceMap fwd_distance) {
  
  BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((SBAStarVisitorConcept<SBAStarVisitor,Graph,Topology>));
  
  return sbastar_bundle<Graph, Topology, SBAStarVisitor, NcSelector, KeyMap, PositionMap, WeightMap,
                        DensityMap, ConstrictionMap, DistanceMap, PredecessorMap, FwdDistanceMap>(
    g, start_vertex, goal_vertex, super_space, vis, select_neighborhood, key, position, weight,
    density, constriction, distance, predecessor, fwd_distance);
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
  * \param goal_vertex The goal point of the algorithm, on the graph.
  * \param super_space A topology (as defined by the Boost Graph Library). This topology 
  *        should not include collision checking in its distance metric.
  * \param vis A SBA* visitor implementing the FADPRMVisitorConcept. This is the 
  *        main point of customization and recording of results that the 
  *        user can implement.
  * \param select_neighborhood A callable object (functor) that can select a list of 
  *        vertices of the graph that ought to be connected to a new 
  *        vertex. The list should be sorted in order of increasing "distance".
  * \param key The property-map which stores the AD* key-values associated to each vertex.
  * \param position A mapping that implements the MutablePropertyMap Concept. Also,
  *        the value_type of this map should be the same type as the topology's 
  *        value_type.
  * \param weight The property-map which stores the weight of each edge-property object (the cost of travel
  *        along the edge).
  * \param density A property-map that provides the expected common information associated with a sample drawn near 
  *        to a vertex w.r.t. the current neighborhood of that vertex.
  * \param constriction A property-map that provides the probability of a collision when a sample is drawn near to a 
  *        vertex (i.e., that a sample near this vertex will not be in the free-space).
  * \param distance The property-map which stores the estimated distance of each vertex from the start.
  * \param predecessor The property-map which will store the resulting path by connecting 
  *        vertices together with their optimal predecessor (follow in reverse to discover the 
  *        complete path).
  * \param fwd_distance The property-map which stores the estimated distance of each vertex to the goal.
  * \param successor The property-map which will store the resulting path by connecting 
  *        vertices together with their optimal successor (follow in order to discover the 
  *        remaining path to the goal).
  */

template <typename Graph, typename Topology, typename SBAStarVisitor, typename NcSelector,
          typename KeyMap, typename PositionMap, typename WeightMap,
          typename DensityMap, typename ConstrictionMap, 
          typename DistanceMap, typename PredecessorMap, typename FwdDistanceMap, typename SuccessorMap >
sbastar_bundle<Graph, Topology, SBAStarVisitor, NcSelector, KeyMap, PositionMap, WeightMap,
               DensityMap, ConstrictionMap, DistanceMap, PredecessorMap, FwdDistanceMap, SuccessorMap>
  make_sbastar_bundle(Graph &g, 
    typename boost::graph_traits<Graph>::vertex_descriptor start_vertex, 
    typename boost::graph_traits<Graph>::vertex_descriptor goal_vertex, 
    const Topology& super_space, SBAStarVisitor vis, NcSelector select_neighborhood, 
    KeyMap key, PositionMap position, WeightMap weight, 
    DensityMap density, ConstrictionMap constriction, DistanceMap distance, 
    PredecessorMap predecessor, FwdDistanceMap fwd_distance, SuccessorMap successor) {
  
  BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<Graph>));
  BOOST_CONCEPT_ASSERT((ReaK::pp::MetricSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((SBAStarBidirVisitorConcept<SBAStarVisitor,Graph,Topology>));
  
  return sbastar_bundle<Graph, Topology, SBAStarVisitor, NcSelector, KeyMap, PositionMap, WeightMap,
                        DensityMap, ConstrictionMap, DistanceMap, PredecessorMap, FwdDistanceMap, SuccessorMap>(
    g, start_vertex, goal_vertex, super_space, vis, select_neighborhood,
    key, position, weight, density, constriction, distance, predecessor, fwd_distance, successor);
};




/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the SBA* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  */
template <typename SBAStarBundle>
void generate_sbastar_no_init(const SBAStarBundle& bdl) {
  detail::generate_sbastar_no_init_impl(
    *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
    pruned_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight, 
    bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
    bdl.m_fwd_distance, bdl.m_select_neighborhood);
};

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the SBA* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  */
template <typename SBAStarBundle>
void generate_sbastar(const SBAStarBundle& bdl) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key, bdl.m_distance, bdl.m_predecessor);
  generate_sbastar_no_init(bdl);
};


/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional SBA* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  */
template <typename SBAStarBundle>
void generate_sbastar_bidir_no_init(const SBAStarBundle& bdl) {
  detail::generate_sbastar_bidir_no_init_impl(
    *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space), bdl.m_vis, 
    pruned_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight, 
    bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
    bdl.m_fwd_distance, bdl.m_successor, bdl.m_select_neighborhood);
};

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Bi-directional SBA* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  */
template <typename SBAStarBundle>
void generate_sbastar_bidir(const SBAStarBundle& bdl) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key, 
                                   bdl.m_distance, bdl.m_predecessor, 
                                   bdl.m_fwd_distance, bdl.m_successor);
  
  generate_sbastar_bidir_no_init(bdl);
};





/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Lazy-SBA* algorithm, without initialization of the existing graph.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_sbastar_no_init(const SBAStarBundle& bdl) {
  detail::generate_sbastar_no_init_impl(
    *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
    lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight, 
    bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
    bdl.m_fwd_distance, bdl.m_select_neighborhood);
};

/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the search.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_sbastar(const SBAStarBundle& bdl) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key, bdl.m_distance, bdl.m_predecessor);
  generate_lazy_sbastar_no_init(bdl);
};


/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Bi-directional Lazy-SBA* algorithm, without initialization of the existing graph.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_sbastar_bidir_no_init(const SBAStarBundle& bdl) {
  detail::generate_sbastar_bidir_no_init_impl(
    *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space), bdl.m_vis, 
    lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight, 
    bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
    bdl.m_fwd_distance, bdl.m_successor, bdl.m_select_neighborhood);
};

/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Bi-directional Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the search.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_sbastar_bidir(const SBAStarBundle& bdl) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key, 
                                   bdl.m_distance, bdl.m_predecessor, 
                                   bdl.m_fwd_distance, bdl.m_successor);
  
  generate_lazy_sbastar_bidir_no_init(bdl);
};




/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Lazy-SBA* algorithm, without initialization of the existing graph.
 * This function uses a branch-and-bound heuristic to limit the number of nodes.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_bnb_sbastar_no_init(const SBAStarBundle& bdl) {
  if( bdl.m_goal_vertex == boost::graph_traits<typename SBAStarBundle::graph_type>::null_vertex() ) {
    detail::generate_sbastar_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
      lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight, 
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
      bdl.m_fwd_distance, bdl.m_select_neighborhood);
  } else {
    detail::generate_sbastar_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
      branch_and_bound_connector<typename SBAStarBundle::graph_type>(
        *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex), 
      bdl.m_key, bdl.m_position, bdl.m_weight, 
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
      bdl.m_fwd_distance, bdl.m_select_neighborhood);
  };
};

/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the search.
 * This function uses a branch-and-bound heuristic to limit the number of nodes.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_bnb_sbastar(const SBAStarBundle& bdl) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key, bdl.m_distance, bdl.m_predecessor);
  generate_lazy_bnb_sbastar_no_init(bdl);
};



/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Bi-directional Lazy-SBA* algorithm, without initialization of the existing graph.
 * This function uses a branch-and-bound heuristic to limit the number of nodes.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_bnb_sbastar_bidir_no_init(const SBAStarBundle& bdl) {
  detail::generate_sbastar_bidir_no_init_impl(
    *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space), bdl.m_vis, 
    branch_and_bound_connector<typename SBAStarBundle::graph_type>(
      *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex), 
    bdl.m_key, bdl.m_position, bdl.m_weight, 
    bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
    bdl.m_fwd_distance, bdl.m_successor, bdl.m_select_neighborhood);
};

/**
 * This function template generates a roadmap to connect a goal location to a start location
 * using the Bi-directional Lazy-SBA* algorithm, with initialization of the existing graph to (re)start the search.
 * This function uses a branch-and-bound heuristic to limit the number of nodes.
 * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
 * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
 */
template <typename SBAStarBundle>
void generate_lazy_bnb_sbastar_bidir(const SBAStarBundle& bdl) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key, 
                                   bdl.m_distance, bdl.m_predecessor, 
                                   bdl.m_fwd_distance, bdl.m_successor);
  generate_lazy_bnb_sbastar_bidir_no_init(bdl);
};




};

};

#endif


















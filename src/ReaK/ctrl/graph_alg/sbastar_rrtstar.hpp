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

#include <utility>
#include <cmath>

#include <ReaK/ctrl/topologies/metric_space_concept.hpp>
#include <ReaK/ctrl/topologies/random_sampler_concept.hpp>
#include <ReaK/core/base/global_rng.hpp>

#include <boost/tuple/tuple.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/detail/d_ary_heap.hpp>

// BGL-Extra includes:
#include <boost/graph/more_property_tags.hpp>
#include <boost/graph/more_property_maps.hpp>

#include "sbastar_search.hpp"
#include "node_generators.hpp"
#include "lazy_connector.hpp"
#include "branch_and_bound_connector.hpp"

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
  * the visitor should model SBAStarVisitorConcept and NodePullingVisitorConcept.
  * 
  * \tparam Visitor The visitor class to be tested for modeling this visitor concept.
  * \tparam Graph The graph type on which the visitor should be able to act.
  * \tparam Topology The topology type on which the visitor class is required to work with.
  */
template <typename Visitor, typename Graph, typename Topology>
struct SBARRTStarVisitorConcept : SBAStarVisitorConcept<Visitor, Graph, Topology> {
  
  BOOST_CONCEPT_ASSERT((NodePullingVisitorConcept<Visitor, Graph, Topology>));
  
  BOOST_CONCEPT_USAGE(SBARRTStarVisitorConcept) { };
};

/**
  * This class is simply an archetype visitor for the SBA*-RRT* algorithm.
  */
template <typename Topology>
struct sbarrtstar_visitor_archetype : 
  sbastar_visitor_archetype<Topology>, 
  node_pulling_visitor_archetype { };

/**
  * This concept class defines the valid expressions required of a class to be used as a visitor 
  * class for the bi-directional SBA*-RRT* algorithm. A visitor class is essentially a class that regroups a number of 
  * callback functions that can be used to inject customization into the bi-directional SBA*-RRT* algorithm. 
  * In other words, the visitor pattern in generic programming is an implementation of IoC 
  * (Inversion of Control), since the bi-directional SBA*-RRT* algorithm is in control of execution, 
  * but custom behavior can be injected in several places, even blocking the algorithm if needed.
  * 
  * Required concepts:
  * 
  * The visitor class should model SBARRTStarVisitorConcept, NodeBackPullingVisitorConcept and NodeBackPushingVisitorConcept.
  * 
  * \tparam Visitor The visitor class to be tested for modeling an AD* visitor concept.
  * \tparam Graph The graph type on which the visitor should be able to act.
  * \tparam Topology The topology type on which the visitor class is required to work with.
  */
template <typename Visitor, typename Graph, typename Topology>
struct SBARRTStarBidirVisitorConcept : SBARRTStarVisitorConcept<Visitor,Graph,Topology> {
  
  BOOST_CONCEPT_ASSERT((NodeBackPushingVisitorConcept<Visitor,Graph,Topology>));
  BOOST_CONCEPT_ASSERT((NodeBackPullingVisitorConcept<Visitor,Graph,Topology>));
  
  BOOST_CONCEPT_USAGE(SBARRTStarBidirVisitorConcept)
  { }
};

/**
 * This class is simply an archetype visitor for the bi-directional SBA*-RRT* algorithm.
 */
template <typename Topology>
struct sbarrtstar_bidir_visitor_archetype : 
  sbarrtstar_visitor_archetype<Topology>, 
  node_back_pulling_visitor_archetype,
  node_back_pushing_visitor_archetype<Topology> { };



namespace detail {
  
  template <typename Graph, typename Topology, typename SBARRTStarVisitor, 
            typename MotionGraphConnector, typename SBANodeGenerator, typename RRTNodeGenerator,
            typename MutableQueue, typename NcSelector>
  void sbarrtstar_search_loop(
    Graph &g, const Topology& super_space, SBARRTStarVisitor& sba_vis, 
    MotionGraphConnector connect_vertex, SBANodeGenerator sba_generate_node, RRTNodeGenerator rrt_generate_node, 
    MutableQueue& Q, NcSelector select_neighborhood, double initial_temperature)
  {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename ReaK::pp::topology_traits<Topology>::point_type PositionValue;
    typedef typename Graph::edge_bundled EdgeProp;
    using std::exp; using std::log;
    std::size_t num_rrt_vertices = 0;
    std::size_t num_sba_vertices = 0;
    
    while (!Q.empty() && sba_vis.keep_going()) { 
      double entropy = 1.0 - exp( -initial_temperature / log( double(num_vertices(g)) ) );
      double rand_value = boost::uniform_01<ReaK::global_rng_type&,double>(ReaK::get_global_rng())(); // generate random-number between 0 and 1.
      bool use_sba_sampling = (rand_value > entropy);
      
      PositionValue p_new; Vertex x_near; EdgeProp eprop;
      if(use_sba_sampling) {
        Vertex u = Q.top(); Q.pop();
        
        // stop if the best nodes do not meet the potential threshold.
        while( ! sba_vis.has_search_potential(u, g) ) {
          if( Q.empty() ) {
            u = boost::graph_traits<Graph>::null_vertex();
            break;
          };
          u = Q.top(); Q.pop();
        };
        if( u == boost::graph_traits<Graph>::null_vertex() )
          break; // no more nodes with search potential.
        
        sba_vis.examine_vertex(u, g);
        
        boost::tie(x_near, p_new, eprop) = sba_generate_node(u, g, sba_vis, sba_vis.m_position);
        
        // then push it back on the OPEN queue.
        if( ( x_near != boost::graph_traits<Graph>::null_vertex() ) ||
            ( Q.empty() ) )
          sba_vis.requeue_vertex(u,g);
        
        ++num_sba_vertices;
      } else {
        boost::tie(x_near, p_new, eprop) = rrt_generate_node(g, sba_vis, boost::bundle_prop_to_vertex_prop(sba_vis.m_position, g));
        
        ++num_rrt_vertices;
      };
      
      if( ( x_near != boost::graph_traits<Graph>::null_vertex() ) && 
          ( get(sba_vis.m_predecessor, g[x_near]) != boost::graph_traits<Graph>::null_vertex() ) ) {
        connect_vertex(p_new, x_near, eprop, g, 
                       super_space, sba_vis, sba_vis.m_position, 
                       sba_vis.m_distance, sba_vis.m_predecessor, 
                       sba_vis.m_weight, select_neighborhood);
      };
      
    }; // end while
    
//     std::cout << " SBA* vertices generated = " << num_sba_vertices << std::endl;
//     std::cout << " RRT* vertices generated = " << num_rrt_vertices << std::endl;
  };
  
  template <typename Graph, typename Topology, typename SBARRTStarVisitor,
            typename MotionGraphConnector, typename SBANodeGenerator, typename RRTNodeGenerator, 
            typename MutableQueue, typename NcSelector>
  void sbarrtstar_bidir_loop(
      Graph &g, const Topology& super_space, SBARRTStarVisitor& sba_vis, 
      MotionGraphConnector connect_vertex, SBANodeGenerator sba_generate_node, RRTNodeGenerator rrt_generate_node, 
      MutableQueue& Q, NcSelector select_neighborhood, double initial_temperature)
  { 
    typedef typename ReaK::pp::topology_traits<Topology>::point_type PositionValue;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename Graph::edge_bundled EdgeProp;
    using std::exp; using std::log;
    std::size_t num_rrt_vertices = 0;
    std::size_t num_sba_vertices = 0;
    
    while (!Q.empty() && sba_vis.keep_going()) { 
      double entropy = 1.0 - exp( -initial_temperature / log( double(num_vertices(g)) ) );
      double rand_value = boost::uniform_01<ReaK::global_rng_type&,double>(ReaK::get_global_rng())(); // generate random-number between 0 and 1.
      bool use_sba_sampling = (rand_value > entropy);
      
      Vertex x_near_pred, x_near_succ; 
      PositionValue p_new_pred, p_new_succ; 
      EdgeProp ep_pred, ep_succ;
      if(use_sba_sampling) {
        Vertex u = Q.top(); Q.pop();
        
        // stop if the best nodes do not meet the potential threshold.
        while( ! sba_vis.has_search_potential(u, g) ) {
          if( Q.empty() ) {
            u = boost::graph_traits<Graph>::null_vertex();
            break;
          };
          u = Q.top(); Q.pop();
        };
        if( u == boost::graph_traits<Graph>::null_vertex() )
          break; // no more nodes with search potential.
        
        sba_vis.examine_vertex(u, g);
        
        boost::tie(x_near_pred, p_new_pred, ep_pred, x_near_succ, p_new_succ, ep_succ) = 
          sba_generate_node(u, g, sba_vis, sba_vis.m_position);
        
        // then push it back on the OPEN queue.
        if( ( x_near_pred != boost::graph_traits<Graph>::null_vertex() ) ||
            ( x_near_succ != boost::graph_traits<Graph>::null_vertex() ) || 
            ( Q.empty() ) )
          sba_vis.requeue_vertex(u,g);
        
        ++num_sba_vertices;
      } else {
        
        boost::tie(x_near_pred, p_new_pred, ep_pred, x_near_succ, p_new_succ, ep_succ) = 
          rrt_generate_node(g, sba_vis, boost::bundle_prop_to_vertex_prop(sba_vis.m_position, g));
        
        ++num_rrt_vertices;
      };
      
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
    
//     std::cout << " SBA* vertices generated = " << num_sba_vertices << std::endl;
//     std::cout << " RRT* vertices generated = " << num_rrt_vertices << std::endl;
  };
  
  
  
  
  template <typename Graph, typename Topology, typename SBARRTStarVisitor,
            typename NodeConnector, typename KeyMap, typename PositionMap, typename WeightMap,
            typename DensityMap, typename ConstrictionMap, typename DistanceMap, typename PredecessorMap,
            typename FwdDistanceMap, typename RandomSampler, typename NcSelector>
  void generate_sbarrtstar_no_init_impl(Graph &g, 
      typename boost::graph_traits<Graph>::vertex_descriptor start_vertex, 
      const Topology& super_space, SBARRTStarVisitor vis, 
      NodeConnector connect_vertex, KeyMap key, PositionMap position, WeightMap weight, 
      DensityMap density, ConstrictionMap constriction, 
      DistanceMap distance, PredecessorMap predecessor, 
      FwdDistanceMap fwd_distance, RandomSampler get_sample, 
      NcSelector select_neighborhood, double SA_init_temperature = 0.0)
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
    
    sbastar_bfs_visitor<Graph, SBARRTStarVisitor, MutableQueue, IndexInHeapMap,
                        KeyMap, PositionMap, WeightMap, DensityMap, ConstrictionMap, 
                        DistanceMap, PredecessorMap, FwdDistanceMap> sba_bfs_vis(
      vis, Q, index_in_heap, key, position, weight, 
      density, constriction, distance, predecessor, fwd_distance);
    
    put(distance, g[start_vertex], 0.0);
    put(predecessor, g[start_vertex], start_vertex);
    sba_bfs_vis.requeue_vertex(start_vertex,g);
    
    sbarrtstar_search_loop(g, super_space, sba_bfs_vis, connect_vertex, sba_node_generator(), 
                           rrg_node_generator<Topology, RandomSampler, NcSelector>(&super_space, get_sample, select_neighborhood), 
                           Q, select_neighborhood, SA_init_temperature);
  };
  
  template <typename Graph, typename Topology, typename SBARRTStarVisitor,
            typename NodeConnector, typename KeyMap, typename PositionMap, typename WeightMap,
            typename DensityMap, typename ConstrictionMap, typename DistanceMap, typename PredecessorMap, 
            typename FwdDistanceMap, typename SuccessorMap, typename RandomSampler, typename NcSelector>
  void generate_sbarrtstar_bidir_no_init_impl(Graph &g, 
      typename boost::graph_traits<Graph>::vertex_descriptor start_vertex, 
      typename boost::graph_traits<Graph>::vertex_descriptor goal_vertex, 
      const Topology& super_space, SBARRTStarVisitor vis, 
      NodeConnector connect_vertex, KeyMap key, PositionMap position, WeightMap weight, 
      DensityMap density, ConstrictionMap constriction, 
      DistanceMap distance, PredecessorMap predecessor, 
      FwdDistanceMap fwd_distance, SuccessorMap successor, 
      RandomSampler get_sample, NcSelector select_neighborhood, double SA_init_temperature = 0.0)
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
    
    sbastar_bfs_visitor<Graph, SBARRTStarVisitor, MutableQueue, IndexInHeapMap,
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
    
    sbarrtstar_bidir_loop(
      g, super_space, sba_bfs_vis, connect_vertex, sba_bidir_node_generator(), 
      rrg_bidir_generator<Topology, RandomSampler, NcSelector, PredecessorMap, SuccessorMap>(&super_space, get_sample, select_neighborhood, predecessor, successor), 
      Q, select_neighborhood, SA_init_temperature);
  };
  
  
  
  
}; //end of detail namespace.


/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches 
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_sbarrtstar_no_init(const SBAStarBundle& bdl, 
                                 RandomSampler get_sample, 
                                 double SA_init_temperature = 0.0) {
  BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<typename SBAStarBundle::visitor_type,typename SBAStarBundle::graph_type,typename SBAStarBundle::topology_type>));
  
  detail::generate_sbarrtstar_no_init_impl(
    *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
    pruned_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight, 
    bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
    bdl.m_fwd_distance, get_sample, bdl.m_select_neighborhood, SA_init_temperature);
};

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches 
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_sbarrtstar(const SBAStarBundle& bdl, 
                         RandomSampler get_sample, 
                         double SA_init_temperature = 0.0) {
  BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<typename SBAStarBundle::visitor_type,typename SBAStarBundle::graph_type,typename SBAStarBundle::topology_type>));
  
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key, bdl.m_distance, bdl.m_predecessor);
  
  generate_sbarrtstar_no_init(bdl, get_sample, SA_init_temperature);
  
};


/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches 
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_sbarrtstar_bidir_no_init(const SBAStarBundle& bdl, 
                                       RandomSampler get_sample, 
                                       double SA_init_temperature = 0.0) {
  detail::generate_sbarrtstar_bidir_no_init_impl(
    *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space), bdl.m_vis, 
    pruned_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight, 
    bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
    bdl.m_fwd_distance, bdl.m_successor, get_sample, bdl.m_select_neighborhood, SA_init_temperature);
};

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches 
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_sbarrtstar_bidir(const SBAStarBundle& bdl, 
                               RandomSampler get_sample, 
                               double SA_init_temperature = 0.0) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key, 
                                   bdl.m_distance, bdl.m_predecessor, 
                                   bdl.m_fwd_distance, bdl.m_successor);
  
  generate_sbarrtstar_bidir_no_init(bdl, get_sample, SA_init_temperature);
};



/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches 
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_lazy_sbarrtstar_no_init(const SBAStarBundle& bdl, 
                                      RandomSampler get_sample, 
                                      double SA_init_temperature = 0.0) {
  BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<typename SBAStarBundle::visitor_type,typename SBAStarBundle::graph_type,typename SBAStarBundle::topology_type>));
  
  detail::generate_sbarrtstar_no_init_impl(
    *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
    lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight, 
    bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
    bdl.m_fwd_distance, get_sample, bdl.m_select_neighborhood, SA_init_temperature);
};

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches 
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_lazy_sbarrtstar(const SBAStarBundle& bdl, 
                              RandomSampler get_sample, 
                              double SA_init_temperature = 0.0) {
  BOOST_CONCEPT_ASSERT((SBARRTStarVisitorConcept<typename SBAStarBundle::visitor_type,typename SBAStarBundle::graph_type,typename SBAStarBundle::topology_type>));
  
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key, bdl.m_distance, bdl.m_predecessor);
  
  generate_lazy_sbarrtstar_no_init(bdl, get_sample, SA_init_temperature);
  
};


/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches 
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_lazy_sbarrtstar_bidir_no_init(const SBAStarBundle& bdl, 
                                            RandomSampler get_sample, 
                                            double SA_init_temperature = 0.0) {
  detail::generate_sbarrtstar_bidir_no_init_impl(
    *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space), bdl.m_vis, 
    lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight, 
    bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
    bdl.m_fwd_distance, bdl.m_successor, get_sample, bdl.m_select_neighborhood, SA_init_temperature);
};

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches 
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_lazy_sbarrtstar_bidir(const SBAStarBundle& bdl, 
                                    RandomSampler get_sample, 
                                    double SA_init_temperature = 0.0) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key, 
                                   bdl.m_distance, bdl.m_predecessor, 
                                   bdl.m_fwd_distance, bdl.m_successor);
  
  generate_lazy_sbarrtstar_bidir_no_init(bdl, get_sample, SA_init_temperature);
  
};


/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches 
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_lazy_bnb_sbarrtstar_no_init(const SBAStarBundle& bdl, 
                                          RandomSampler get_sample, 
                                          double SA_init_temperature = 0.0) {
  if( bdl.m_goal_vertex == boost::graph_traits<typename SBAStarBundle::graph_type>::null_vertex() ) {
    detail::generate_sbarrtstar_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
      lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight, 
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
      bdl.m_fwd_distance, get_sample, bdl.m_select_neighborhood, SA_init_temperature);
  } else {
    detail::generate_sbarrtstar_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, *(bdl.m_super_space), bdl.m_vis, 
      branch_and_bound_connector<typename SBAStarBundle::graph_type>(
        *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex), 
      bdl.m_key, bdl.m_position, bdl.m_weight, 
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
      bdl.m_fwd_distance, get_sample, bdl.m_select_neighborhood, SA_init_temperature);
  };
};

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches 
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_lazy_bnb_sbarrtstar(const SBAStarBundle& bdl, 
                                  RandomSampler get_sample, 
                                  double SA_init_temperature = 0.0) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key, bdl.m_distance, bdl.m_predecessor);
  
  generate_lazy_bnb_sbarrtstar_no_init(bdl, get_sample, SA_init_temperature);
  
};


/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, without initialization of the existing graph.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches 
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_lazy_bnb_sbarrtstar_bidir_no_init(const SBAStarBundle& bdl, 
                                                RandomSampler get_sample, 
                                                double SA_init_temperature = 0.0) {
  if( bdl.m_goal_vertex == boost::graph_traits<typename SBAStarBundle::graph_type>::null_vertex() ) {
    detail::generate_sbarrtstar_bidir_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space), bdl.m_vis, 
      lazy_node_connector(), bdl.m_key, bdl.m_position, bdl.m_weight, 
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
      bdl.m_fwd_distance, bdl.m_successor, get_sample, bdl.m_select_neighborhood, SA_init_temperature);
  } else {
    detail::generate_sbarrtstar_bidir_no_init_impl(
      *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex, *(bdl.m_super_space), bdl.m_vis, 
      branch_and_bound_connector<typename SBAStarBundle::graph_type>(
        *(bdl.m_g), bdl.m_start_vertex, bdl.m_goal_vertex), 
      bdl.m_key, bdl.m_position, bdl.m_weight, 
      bdl.m_density, bdl.m_constriction, bdl.m_distance, bdl.m_predecessor, 
      bdl.m_fwd_distance, bdl.m_successor, get_sample, bdl.m_select_neighborhood, SA_init_temperature);
  };
};

/**
  * This function template generates a roadmap to connect a goal location to a start location
  * using the Lazy-SBA*-RRT* algorithm, with initialization of the existing graph to (re)start the search.
  * \tparam SBAStarBundle A SBA* bundle type (see make_sbastar_bundle()).
  * \tparam RandomSampler This is a random-sampler over the topology (see pp::RandomSamplerConcept).
  * \param bdl A const-reference to a SBA* bundle of parameters, see make_sbastar_bundle().
  * \param get_sample A random sampler of positions in the space.
  * \param SA_init_temperature The initial temperature of the Simulated Annealing when used 
  *        as the deciding factor between using RRT* or SBA* samples. When the value is negative,
  *        the RRT* algorithm is used until a solution is found, and then the algorithm switches 
  *        to SBA*.
  */
template <typename SBAStarBundle, typename RandomSampler>
void generate_lazy_bnb_sbarrtstar_bidir(const SBAStarBundle& bdl, 
                                        RandomSampler get_sample, 
                                        double SA_init_temperature = 0.0) {
  detail::initialize_sbastar_nodes(*(bdl.m_g), bdl.m_vis, bdl.m_key, 
                                   bdl.m_distance, bdl.m_predecessor, 
                                   bdl.m_fwd_distance, bdl.m_successor);
  
  generate_lazy_bnb_sbarrtstar_bidir_no_init(bdl, get_sample, SA_init_temperature);
  
};



};

};

#endif

















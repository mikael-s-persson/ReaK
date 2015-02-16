/**
 * \file fadprm_path_planner.tpp
 * 
 * This library contains template definitions of a class to solve path planning problems using the 
 * Flexible Anytime Dynamic Probabilistic Road-map (FADPRM) algorithm (or one of its variants). 
 * Given a C_free (configuration space restricted to non-colliding points) and a 
 * result reporting policy, this class will probabilistically construct a motion-graph 
 * that will connect a starting point and a goal point with a path through C-free 
 * that is as close as possible to the optimal path in terms of distance.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_FADPRM_PATH_PLANNER_TPP
#define REAK_FADPRM_PATH_PLANNER_TPP

#include "fadprm_path_planner.hpp"

#include <ReaK/planning/graph_alg/fadprm.hpp>
#include <ReaK/planning/graph_alg/neighborhood_functors.hpp>

#include "motion_graph_structures.hpp"

// BGL-Extra includes:
#include <boost/graph/more_property_maps.hpp>

#include "metric_space_search.hpp"
#include "topological_search.hpp"

#include "p2p_planning_query.hpp"
#include "path_planner_options.hpp"
#include "any_motion_graphs.hpp"
#include "planning_visitors.hpp"
#include "any_knn_synchro.hpp"
#include "density_plan_visitors.hpp"

#include <ReaK/core/base/misc_math.hpp>

#include <stack>

namespace ReaK {
  
namespace pp {
  

/**
 * This class template is used by the FADPRM path-planner as the visitor object needed to 
 * collaborate with the FADPRM algorithms to generate the motion-graph and path-planning solutions.
 * This class template models the FADPRMVisitorConcept.
 */
template <typename FreeSpaceType>
struct fadprm_planner_visitor : density_plan_visitor<FreeSpaceType, prm_density_calculator> {
  typedef density_plan_visitor<FreeSpaceType, prm_density_calculator> base_type;
  typedef fadprm_planner_visitor<FreeSpaceType> self;
  
  typedef typename base_type::planner_base_type planner_base_type;
  typedef typename base_type::query_type query_type;
  
  fadprm_planner_visitor(planner_base_type* aPlanner,
                         query_type* aQuery = NULL,
                         any_knn_synchro* aNNSynchro = NULL,
                         boost::any aStartNode = boost::any(),
                         boost::any aGoalNode = boost::any(),
                         double aDensityCutoff = 0.0) : 
                         base_type(aPlanner, aQuery, aNNSynchro, aStartNode, aGoalNode, aDensityCutoff) { };
  
  // AD* visitor functions:
  
  template <typename EdgeIter, typename Graph>
  std::pair<double, EdgeIter> detect_edge_change(EdgeIter ei, const Graph&) const {
    return std::pair<double, EdgeIter>(0.0, ei);
  };
  
  
};


template <typename FreeSpaceType>
void fadprm_planner<FreeSpaceType>::solve_planning_query(planning_query<FreeSpaceType>& aQuery) {
  
  this->reset_internal_state();
  
  typedef typename subspace_traits<FreeSpaceType>::super_space_type SuperSpace;
  typedef typename topology_traits<SuperSpace>::point_type PointType;
  
  typedef dense_mg_vertex< astar_mg_vertex<FreeSpaceType> > VertexProp;
  typedef optimal_mg_edge<FreeSpaceType> EdgeProp;
  
  typedef typename motion_segment_directionality<FreeSpaceType>::type DirectionalityTag;
  
  typedef mg_vertex_data<FreeSpaceType> BasicVertexProp;
  
  typedef boost::data_member_property_map<PointType, VertexProp > PositionMap;
  PositionMap pos_map = PositionMap(&VertexProp::position);
  
  //typedef boost::data_member_property_map<double, VertexProp > DensityMap;
  //DensityMap dens_map = DensityMap(&VertexProp::density);
  
  double space_dim = double( this->get_space_dimensionality() );
  double space_Lc = aQuery.get_heuristic_to_goal( aQuery.get_start_position() );
  
  shared_ptr<const SuperSpace> sup_space_ptr(&(this->m_space->get_super_space()),null_deleter());
  
  fadprm_planner_visitor<FreeSpaceType> vis(this, &aQuery);
  
  path_planning_p2p_query<FreeSpaceType>* p2p_query_ptr = reinterpret_cast< path_planning_p2p_query<FreeSpaceType>* >(aQuery.castTo(path_planning_p2p_query<FreeSpaceType>::getStaticObjectType()));
  
  
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES

#define RK_FADPRM_PLANNER_INITIALIZE_START_AND_GOAL \
  VertexProp vp_start; \
  vp_start.position = aQuery.get_start_position(); \
  Vertex start_node = add_vertex(std::move(vp_start), motion_graph); \
  motion_graph[start_node].density = 0.0; \
  motion_graph[start_node].heuristic_value = aQuery.get_heuristic_to_goal(motion_graph[start_node].position); \
  motion_graph[start_node].distance_accum = 0.0; \
  motion_graph[start_node].predecessor = start_node; \
  vis.m_start_node = boost::any( start_node ); \
  if( p2p_query_ptr ) { \
    VertexProp vp_goal; \
    vp_goal.position = p2p_query_ptr->goal_pos; \
    Vertex goal_node  = add_vertex(std::move(vp_goal),  motion_graph); \
    motion_graph[goal_node].density = 0.0; \
    motion_graph[goal_node].heuristic_value = 0.0; \
    motion_graph[goal_node].distance_accum = std::numeric_limits<double>::infinity(); \
    motion_graph[goal_node].predecessor = goal_node; \
    vis.m_goal_node = boost::any( goal_node ); \
  };
    
#else
    
#define RK_FADPRM_PLANNER_INITIALIZE_START_AND_GOAL \
  VertexProp vp_start; \
  vp_start.position = aQuery.get_start_position(); \
  Vertex start_node = add_vertex(vp_start, motion_graph); \
  motion_graph[start_node].density = 0.0; \
  motion_graph[start_node].heuristic_value = aQuery.get_heuristic_to_goal(motion_graph[start_node].position); \
  motion_graph[start_node].distance_accum = 0.0; \
  motion_graph[start_node].predecessor = start_node; \
  vis.m_start_node = boost::any( start_node ); \
  if( p2p_query_ptr ) { \
    VertexProp vp_goal; \
    vp_goal.position = p2p_query_ptr->goal_pos; \
    Vertex goal_node  = add_vertex(vp_goal,  motion_graph); \
    motion_graph[goal_node].density = 0.0; \
    motion_graph[goal_node].heuristic_value = 0.0; \
    motion_graph[goal_node].distance_accum = std::numeric_limits<double>::infinity(); \
    motion_graph[goal_node].predecessor = goal_node; \
    vis.m_goal_node = boost::any( goal_node ); \
  };
  
#endif
  
  
  
#define RK_FADPRM_PLANNER_SETUP_DVP_TREE_SYNCHRO(ARITY, TREE_STORAGE) \
  typedef typename boost::property_map< MotionGraphType, PointType BasicVertexProp::* >::type GraphPositionMap; \
  typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, ARITY, random_vp_chooser, TREE_STORAGE > SpacePartType; \
  SpacePartType space_part(motion_graph, sup_space_ptr, get(&BasicVertexProp::position, motion_graph)); \
    \
  typedef multi_dvp_tree_search<MotionGraphType, SpacePartType> NNFinderType; \
  NNFinderType nn_finder; \
  nn_finder.graph_tree_map[&motion_graph] = &space_part; \
    \
  type_erased_knn_synchro< MotionGraphType, NNFinderType > NN_synchro(nn_finder); \
  vis.m_nn_synchro = &NN_synchro; \
    \
  ReaK::graph::star_neighborhood< NNFinderType > nc_selector(nn_finder, space_dim, 3.0 * space_Lc);

//   ReaK::graph::fixed_neighborhood< NNFinderType > nc_selector(nn_finder, 10, this->get_sampling_radius());
  
  
#define RK_FADPRM_PLANNER_SETUP_ALT_TREE_SYNCHRO(ARITY, TREE_STORAGE) \
  typedef dvp_adjacency_list< \
    VertexProp, EdgeProp, SuperSpace, PositionMap, \
    ARITY, random_vp_chooser, TREE_STORAGE, \
    boost::vecBC, DirectionalityTag, boost::listBC > ALTGraph; \
  typedef typename ALTGraph::adj_list_type MotionGraphType; \
  typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex; \
   \
  ALTGraph space_part(sup_space_ptr, pos_map); \
  MotionGraphType motion_graph = space_part.get_adjacency_list(); \
   \
  typedef multi_dvp_tree_search<MotionGraphType, ALTGraph> NNFinderType; \
  NNFinderType nn_finder; \
  nn_finder.graph_tree_map[&motion_graph] = &space_part; \
   \
  any_knn_synchro NN_synchro; \
  vis.m_nn_synchro = &NN_synchro; \
   \
  ReaK::graph::star_neighborhood< NNFinderType > nc_selector(nn_finder, space_dim, 3.0 * space_Lc);
  
//   ReaK::graph::fixed_neighborhood< NNFinderType > nc_selector(nn_finder, 10, this->get_sampling_radius());
  
  
#define RK_FADPRM_PLANNER_MAKE_GENERATE_CALL_FADPRM_FUNCTION \
  ReaK::graph::generate_fadprm( \
    motion_graph, start_node, *(this->m_space), \
    get(&VertexProp::heuristic_value, motion_graph),  \
    vis, \
    get(&VertexProp::predecessor, motion_graph),  \
    get(&VertexProp::distance_accum, motion_graph),  \
    get(&VertexProp::key_value, motion_graph),  \
    get(&EdgeProp::weight, motion_graph),  \
    get(&VertexProp::density, motion_graph),  \
    pos_map,  \
    nc_selector, \
    get(&VertexProp::astar_color, motion_graph),  \
    this->m_initial_relaxation);
  
  
  
  if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == ADJ_LIST_MOTION_GRAPH) {
    
    typedef boost::adjacency_list_BC< boost::vecBC, boost::vecBC, DirectionalityTag, VertexProp, EdgeProp> MotionGraphType;
    typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex;
    
    MotionGraphType motion_graph;
    
    RK_FADPRM_PLANNER_INITIALIZE_START_AND_GOAL
    
    if((this->m_data_structure_flags & KNN_METHOD_MASK) == LINEAR_SEARCH_KNN) {
      
      typedef linear_neighbor_search<MotionGraphType> NNFinderType;
      NNFinderType nn_finder;
      
      ReaK::graph::star_neighborhood< NNFinderType > nc_selector(nn_finder, space_dim, 3.0 * space_Lc);
//       ReaK::graph::fixed_neighborhood< NNFinderType > nc_selector(nn_finder, 10, this->get_sampling_radius());
      
      any_knn_synchro NN_synchro;
      vis.m_nn_synchro = &NN_synchro;
      
      RK_FADPRM_PLANNER_MAKE_GENERATE_CALL_FADPRM_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
      
      RK_FADPRM_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, boost::bfl_d_ary_tree_storage<2>)
      
      RK_FADPRM_PLANNER_MAKE_GENERATE_CALL_FADPRM_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
      
      RK_FADPRM_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, boost::bfl_d_ary_tree_storage<4>)
      
      RK_FADPRM_PLANNER_MAKE_GENERATE_CALL_FADPRM_FUNCTION
      
#ifdef RK_PLANNERS_ENABLE_VEBL_TREE
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
      
      RK_FADPRM_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, boost::vebl_d_ary_tree_storage<2>)
      
      RK_FADPRM_PLANNER_MAKE_GENERATE_CALL_FADPRM_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
      
      RK_FADPRM_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, boost::vebl_d_ary_tree_storage<4>)
      
      RK_FADPRM_PLANNER_MAKE_GENERATE_CALL_FADPRM_FUNCTION
        
#endif
      
    };
    
#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT
    
  } else if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == DVP_ADJ_LIST_MOTION_GRAPH) {
    
    if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
      
      RK_FADPRM_PLANNER_SETUP_ALT_TREE_SYNCHRO(2, boost::bfl_d_ary_tree_storage<2>)
      
      RK_FADPRM_PLANNER_INITIALIZE_START_AND_GOAL
      
      RK_FADPRM_PLANNER_MAKE_GENERATE_CALL_FADPRM_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
      
      RK_FADPRM_PLANNER_SETUP_ALT_TREE_SYNCHRO(4, boost::bfl_d_ary_tree_storage<4>)
      
      RK_FADPRM_PLANNER_INITIALIZE_START_AND_GOAL
      
      RK_FADPRM_PLANNER_MAKE_GENERATE_CALL_FADPRM_FUNCTION
      
#ifdef RK_PLANNERS_ENABLE_VEBL_TREE
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
      
      RK_FADPRM_PLANNER_SETUP_ALT_TREE_SYNCHRO(2, boost::vebl_d_ary_tree_storage<2>)
      
      RK_FADPRM_PLANNER_INITIALIZE_START_AND_GOAL
      
      RK_FADPRM_PLANNER_MAKE_GENERATE_CALL_FADPRM_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
      
      RK_FADPRM_PLANNER_SETUP_ALT_TREE_SYNCHRO(4, boost::vebl_d_ary_tree_storage<4>)
      
      RK_FADPRM_PLANNER_INITIALIZE_START_AND_GOAL
      
      RK_FADPRM_PLANNER_MAKE_GENERATE_CALL_FADPRM_FUNCTION
        
#endif
      
    };
      
#endif
    
  };
  
#undef RK_FADPRM_PLANNER_INITIALIZE_START_AND_GOAL
#undef RK_FADPRM_PLANNER_SETUP_DVP_TREE_SYNCHRO
#undef RK_FADPRM_PLANNER_SETUP_ALT_TREE_SYNCHRO
#undef RK_FADPRM_PLANNER_MAKE_GENERATE_CALL_FADPRM_FUNCTION
  
};


};

};

#endif


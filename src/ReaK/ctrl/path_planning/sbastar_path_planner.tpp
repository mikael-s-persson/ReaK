/**
 * \file sbastar_path_planner.tpp
 * 
 * This library contains template definitions of a class to solve path planning problems using the 
 * Sampling-based A* algorithm (or one of its variants). Given a C_free (configuration space
 * restricted to non-colliding points) and a result reporting policy, this class 
 * will probabilistically construct a motion-graph that will connect a starting point 
 * and a goal point with a path through C-free that is as close as possible to the 
 * optimal path in terms of distance. The planner uses a selectable variant of the 
 * Sampling-based A* (SBA*) algorithm, including the basic version, the SBA*-RRT*
 * alternating algorithm, and the Anytime SBA* algorithm. In all cases, collision 
 * checking and connectivity can be either full or lazy (and pruned) to either construct
 * a full-connectivity graph containing only collision-free edges, or a single-query motion-tree
 * that includes only optimal edges (whose collisions are checked lazily).
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

#ifndef REAK_SBASTAR_PATH_PLANNER_TPP
#define REAK_SBASTAR_PATH_PLANNER_TPP

#include "sbastar_path_planner.hpp"

#include "graph_alg/lazy_sbastar.hpp"
#include "graph_alg/sbastar_rrtstar.hpp"
#include "graph_alg/anytime_sbastar.hpp"

#include "motion_graph_structures.hpp"


// BGL-Extra includes:
#include <boost/graph/more_property_tags.hpp>
#include <boost/graph/more_property_maps.hpp>


#include "metric_space_search.hpp"
#include "topological_search.hpp"

#include "p2p_planning_query.hpp"
#include "path_planner_options.hpp"
#include "graph_alg/neighborhood_functors.hpp"
#include "any_motion_graphs.hpp"
#include "density_plan_visitors.hpp"

namespace ReaK {
  
namespace pp {


template <typename FreeSpaceType, typename IsBidirPlanner>
struct sbastar_planner_bundle {
  
  typedef boost::mpl::and_< IsBidirPlanner, is_reversible_space<FreeSpaceType> > is_bidir;
  
  typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
  typedef typename topology_traits<super_space_type>::point_type point_type;
  
  typedef typename boost::mpl::if_< is_bidir,
    recursive_dense_mg_vertex< bidir_astar_mg_vertex<FreeSpaceType> >,
    recursive_dense_mg_vertex< astar_mg_vertex<FreeSpaceType> > >::type vertex_prop;
  typedef mg_vertex_data<FreeSpaceType> basic_vertex_prop;
  typedef optimal_mg_edge<FreeSpaceType> edge_prop;
  
  typedef typename motion_segment_directionality<FreeSpaceType>::type directionality_tag;
  
  typedef density_plan_visitor<FreeSpaceType, sbastar_density_calculator> visitor_type;
  
  typedef boost::data_member_property_map<point_type, vertex_prop > position_map;
  typedef boost::data_member_property_map<double, vertex_prop> density_map;
  typedef boost::data_member_property_map<double, vertex_prop> constriction_map;
  typedef boost::data_member_property_map<double, vertex_prop> distance_map;
  typedef boost::data_member_property_map<double, vertex_prop> fwd_distance_map;
  typedef boost::data_member_property_map<std::size_t, vertex_prop> predecessor_map;
  typedef typename boost::mpl::if_< is_bidir,
    boost::data_member_property_map<std::size_t, vertex_prop>,
    ReaK::graph::detail::null_vertex_prop_map<Graph> >::type successor_map;
  typedef boost::data_member_property_map<double, edge_prop > weight_map;
  
  struct ls_motion_graph {
    typedef boost::adjacency_list_BC< boost::vecBC, boost::poolBC,
      directionality_tag, vertex_prop, edge_prop> type;
    typedef typename boost::graph_traits<type>::vertex_descriptor vertex_type;
    
    typedef linear_neighbor_search<MotionGraphType> nn_finder_type;
    static nn_finder_type get_nn_finder() { return nn_finder_type(); };
    
    typedef any_knn_synchro nn_synchro_type;
    static nn_synchro_type get_nn_synchro() {
      return nn_synchro_type();
    };
    
    static type get_motion_graph() {
      return type();
    };
  };
  
  template <unsigned int Arity, typename TreeStorageTag>
  struct dvp_motion_graph {
    typedef boost::adjacency_list_BC< boost::vecBC, boost::poolBC,
      directionality_tag, vertex_prop, edge_prop> type;
    typedef typename boost::graph_traits<type>::vertex_descriptor vertex_type;
    
    typedef linear_neighbor_search<MotionGraphType> nn_finder_type;
    static nn_finder_type get_nn_finder() { return nn_finder_type; };
    
    typedef typename boost::property_map< type, point_type basic_vertex_prop::* >::type graph_position_map;
    typedef dvp_tree<vertex_type, super_space_type, graph_position_map, Arity, random_vp_chooser, TreeStorageTag > space_part_type;
    static space_part_type get_space_part(type& mg, shared_ptr<const super_space_type> s_ptr) {
      return space_part_type(mg, s_ptr, get(&basic_vertex_prop::position, mg));
    };
    
    typedef multi_dvp_tree_search<type, space_part_type> nn_finder_type;
    static nn_finder_type get_nn_finder(type& mg, space_part_type& space_part) { 
      nn_finder_type nn_finder;
      nn_finder.graph_tree_map[&mg] = &space_part;
      return nn_finder;
    };
    
    typedef type_erased_knn_synchro< type, nn_finder_type > nn_synchro_type;
    static nn_synchro_type get_nn_synchro(nn_finder_type& nn_finder) {
      return nn_synchro_type(nn_finder);
    };
    
    static type get_motion_graph() {
      return type();
    };
  };
  
  template <unsigned int Arity, typename TreeStorageTag>
  struct alt_motion_graph {
    typedef dvp_adjacency_list< vertex_prop, edge_prop, super_space_type, 
      position_map, Arity, random_vp_chooser, TreeStorageTag,
      boost::vecBC, directionality_tag, boost::listBC > alt_graph_type;
    typedef typename alt_graph_type::adj_list_type type;
    typedef alt_graph_type space_part_type;
    
    static space_part_type get_space_part(shared_ptr<const super_space_type> s_ptr) {
      return space_part_type(s_ptr, position_map(&vertex_prop::position));
    };
    
    typedef multi_dvp_tree_search<type, space_part_type> nn_finder_type;
    static nn_finder_type get_nn_finder(type& mg, space_part_type& space_part) { 
      nn_finder_type nn_finder;
      nn_finder.graph_tree_map[&mg] = &space_part;
      return nn_finder;
    };
    
    typedef any_knn_synchro nn_synchro_type;
    static nn_synchro_type get_nn_synchro() {
      return nn_synchro_type();
    };
    
    static type get_motion_graph(space_part_type& space_part) {
      return space_part.get_adjacency_list();
    };
  };
  
  typedef ReaK::graph::sbastar_bundle<Graph, super_space_type, visitor_type, NcSelector,
          typename boost::property_map<Graph, double vertex_prop::*>::type, 
          position_map, weight_map, density_map, constriction_map, 
          distance_map, predecessor_map, fwd_distance_map, successor_map > type;
};



template <typename FreeSpaceType>
void sbastar_planner<FreeSpaceType>::solve_planning_query(planning_query<FreeSpaceType>& aQuery) {
  
  this->reset_internal_state();
  
  double space_dim = double( this->get_space_dimensionality() );
  double space_Lc = aQuery.get_heuristic_to_goal( aQuery.get_start_position() );
  
  density_plan_visitor<FreeSpaceType, sbastar_density_calculator> vis(
    this, &aQuery, NULL, boost::any(), boost::any(), 
    this->m_init_dens_threshold);
  
  path_planning_p2p_query<FreeSpaceType>* p2p_query_ptr = reinterpret_cast< path_planning_p2p_query<FreeSpaceType>* >(aQuery.castTo(path_planning_p2p_query<FreeSpaceType>::getStaticObjectType()));
  
  typedef typename subspace_traits<FreeSpaceType>::super_space_type SuperSpace;
  shared_ptr<const SuperSpace> sup_space_ptr(&(this->m_space->get_super_space()),null_deleter());
  
  typedef typename topology_traits<SuperSpace>::point_type PointType;
  
  typedef recursive_dense_mg_vertex< astar_mg_vertex<FreeSpaceType> > VertexProp;
  typedef optimal_mg_edge<FreeSpaceType> EdgeProp;
  
  typedef typename motion_segment_directionality<FreeSpaceType>::type DirectionalityTag;
  
  typedef mg_vertex_data<FreeSpaceType> BasicVertexProp;
  
  typedef boost::data_member_property_map<PointType, VertexProp > PositionMap;
  PositionMap pos_map = PositionMap(&VertexProp::position);
  
  typedef boost::data_member_property_map<double, VertexProp> DensityMap;
  DensityMap dens_map = DensityMap(&VertexProp::density);
  
  typedef boost::data_member_property_map<double, VertexProp> ConstrictionMap;
  ConstrictionMap cons_map = ConstrictionMap(&VertexProp::constriction);
  
  typedef boost::data_member_property_map<double, VertexProp> DistanceMap;
  DistanceMap dist_map = DistanceMap(&VertexProp::distance_accum);
  
  typedef boost::data_member_property_map<double, VertexProp> HeuristicMap;
  HeuristicMap heuristic_map = HeuristicMap(&VertexProp::heuristic_value);
  
  typedef boost::data_member_property_map<std::size_t, VertexProp> PredecessorMap;
  PredecessorMap pred_map = PredecessorMap(&VertexProp::predecessor);
  
  typedef boost::data_member_property_map<double, EdgeProp > WeightMap;
  WeightMap weight_map = WeightMap(&EdgeProp::weight);
  
  // Some MACROs to reduce the size of the code below.
  
#define RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE \
  VertexProp vp_start; \
  vp_start.position = aQuery.get_start_position(); \
  Vertex start_node = add_vertex(vp_start, motion_graph); \
  vis.m_start_node = boost::any( start_node ); \
  if( p2p_query_ptr ) { \
    VertexProp vp_goal; \
    vp_goal.position = p2p_query_ptr->goal_pos; \
    Vertex goal_node  = add_vertex(vp_goal,  motion_graph); \
    vis.m_goal_node = boost::any( goal_node ); \
    vis.initialize_vertex(goal_node, motion_graph); \
  };\
  vis.initialize_vertex(start_node, motion_graph);
  
  
#define RK_SBASTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(ARITY, TREE_STORAGE) \
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


#define RK_SBASTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO(ARITY, TREE_STORAGE) \
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
  
  
  
#define RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, boost::any_cast<Vertex>( vis.m_start_node ), *sup_space_ptr, vis, \
        nc_selector, get(&VertexProp::key_value, motion_graph), \
        pos_map, weight_map, dens_map, cons_map, dist_map, pred_map, \
        heuristic_map)

#define RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE_WITH_GOAL \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, boost::any_cast<Vertex>( vis.m_start_node ), \
        ( vis.m_goal_node.empty() ? boost::graph_traits<MotionGraphType>::null_vertex() : boost::any_cast<Vertex>( vis.m_goal_node ) ), \
        *sup_space_ptr, vis, \
        nc_selector, get(&VertexProp::key_value, motion_graph), \
        pos_map, weight_map, dens_map, cons_map, dist_map, pred_map, \
        heuristic_map)


#define RK_SBASTAR_PLANNER_CALL_SBASTAR_FUNCTION \
    ReaK::graph::generate_sbastar( \
      RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE );
  
  
#define RK_SBASTAR_PLANNER_CALL_LAZY_SBASTAR_FUNCTION \
    ReaK::graph::generate_lazy_sbastar( \
      RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE );
  
#define RK_SBASTAR_PLANNER_CALL_LAZY_BNB_SBASTAR_FUNCTION \
    ReaK::graph::generate_lazy_bnb_sbastar( \
      RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE_WITH_GOAL );
  
#define RK_SBASTAR_PLANNER_CALL_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_sbarrtstar( \
      RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE_WITH_GOAL, \
      get(random_sampler, *sup_space_ptr), \
      this->m_SA_init_temperature);
  
  
#define RK_SBASTAR_PLANNER_CALL_LAZY_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_lazy_sbarrtstar( \
      RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE_WITH_GOAL, \
      get(random_sampler, *sup_space_ptr), \
      this->m_SA_init_temperature);
  
  
#define RK_SBASTAR_PLANNER_CALL_LAZY_BNB_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_lazy_bnb_sbarrtstar( \
      RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE_WITH_GOAL, \
      get(random_sampler, *sup_space_ptr), \
      this->m_SA_init_temperature);
   
  
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_SBASTAR_FUNCTION \
    ReaK::graph::generate_anytime_sbastar( \
      RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE, \
      this->m_init_relaxation);
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBASTAR_FUNCTION \
    ReaK::graph::generate_anytime_lazy_sbastar( \
      RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE, \
      this->m_init_relaxation);
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_BNB_SBASTAR_FUNCTION \
    ReaK::graph::generate_anytime_lazy_bnb_sbastar( \
      RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE_WITH_GOAL, \
      this->m_init_relaxation);
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_anytime_sbarrtstar( \
      RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE_WITH_GOAL, \
      get(random_sampler, *sup_space_ptr), \
      this->m_init_relaxation, this->m_SA_init_temperature);
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_anytime_lazy_sbarrtstar( \
      RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE_WITH_GOAL, \
      get(random_sampler, *sup_space_ptr), \
      this->m_init_relaxation, this->m_SA_init_temperature);
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_BNB_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_anytime_lazy_bnb_sbarrtstar( \
      RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE_WITH_GOAL, \
      get(random_sampler, *sup_space_ptr), \
      this->m_init_relaxation, this->m_SA_init_temperature);
  
  
#define RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION \
      if(((this->m_planning_method_flags & ADDITIONAL_PLANNING_BIAS_MASK) & PLAN_WITH_ANYTIME_HEURISTIC) && (this->m_init_relaxation > 1e-6)) { \
        if((this->m_planning_method_flags & COLLISION_CHECKING_POLICY_MASK) == EAGER_COLLISION_CHECKING) { \
          if((this->m_planning_method_flags & ADDITIONAL_PLANNING_BIAS_MASK) & PLAN_WITH_VORONOI_PULL) { \
            RK_SBASTAR_PLANNER_CALL_ANYTIME_SBARRTSTAR_FUNCTION \
          } else { /* assume nominal method only. */ \
            RK_SBASTAR_PLANNER_CALL_ANYTIME_SBASTAR_FUNCTION \
          }; \
        } else { /* assume lazy collision checking */ \
          if((this->m_planning_method_flags & ADDITIONAL_PLANNING_BIAS_MASK) & PLAN_WITH_VORONOI_PULL) { \
            if(this->m_planning_method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) { \
              RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_BNB_SBARRTSTAR_FUNCTION \
            } else { /* assume nominal method only. */ \
              RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBARRTSTAR_FUNCTION \
            }; \
          } else { /* assume nominal method only. */ \
            if(this->m_planning_method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) { \
              RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_BNB_SBASTAR_FUNCTION \
            } else { /* assume nominal method only. */ \
              RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBASTAR_FUNCTION \
            }; \
          }; \
        }; \
      } else { \
        if((this->m_planning_method_flags & COLLISION_CHECKING_POLICY_MASK) == EAGER_COLLISION_CHECKING) { \
          if((this->m_planning_method_flags & ADDITIONAL_PLANNING_BIAS_MASK) & PLAN_WITH_VORONOI_PULL) { \
            RK_SBASTAR_PLANNER_CALL_SBARRTSTAR_FUNCTION \
          } else { /* assume nominal method only. */ \
            RK_SBASTAR_PLANNER_CALL_SBASTAR_FUNCTION \
          }; \
        } else { /* assume lazy collision checking */ \
          if((this->m_planning_method_flags & ADDITIONAL_PLANNING_BIAS_MASK) & PLAN_WITH_VORONOI_PULL) { \
            if(this->m_planning_method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) { \
              RK_SBASTAR_PLANNER_CALL_LAZY_BNB_SBARRTSTAR_FUNCTION \
            } else { /* assume nominal method only. */ \
              RK_SBASTAR_PLANNER_CALL_LAZY_SBARRTSTAR_FUNCTION \
            }; \
          } else { /* assume nominal method only. */ \
            if(this->m_planning_method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) { \
              RK_SBASTAR_PLANNER_CALL_LAZY_BNB_SBASTAR_FUNCTION \
            } else { /* assume nominal method only. */ \
              RK_SBASTAR_PLANNER_CALL_LAZY_SBASTAR_FUNCTION \
            }; \
          }; \
        }; \
      };
  
  
  if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == ADJ_LIST_MOTION_GRAPH) {
    
    typedef boost::adjacency_list_BC< boost::vecBC, boost::poolBC,
      DirectionalityTag, VertexProp, EdgeProp> MotionGraphType;
    
    typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex;
    
    MotionGraphType motion_graph;
    
    RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
    
    if((this->m_data_structure_flags & KNN_METHOD_MASK) == LINEAR_SEARCH_KNN) {
      
      typedef linear_neighbor_search<MotionGraphType> NNFinderType;
      NNFinderType nn_finder;
      
      ReaK::graph::star_neighborhood< NNFinderType > nc_selector(nn_finder, space_dim, 3.0 * space_Lc);
//       ReaK::graph::fixed_neighborhood< NNFinderType > nc_selector(nn_finder, 10, this->get_sampling_radius());
      
      any_knn_synchro NN_synchro;
      vis.m_nn_synchro = &NN_synchro;
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, boost::bfl_d_ary_tree_storage<2>)
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, boost::bfl_d_ary_tree_storage<4>)
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
#ifdef RK_PLANNERS_ENABLE_VEBL_TREE
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, boost::vebl_d_ary_tree_storage<2>)
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, boost::vebl_d_ary_tree_storage<4>)
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
        
#endif
      
    };
    
#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT
    
  } else if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == DVP_ADJ_LIST_MOTION_GRAPH) {
    
    if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO(2, boost::bfl_d_ary_tree_storage<2>)
      
      RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO(4, boost::bfl_d_ary_tree_storage<4>)
      
      RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
#ifdef RK_PLANNERS_ENABLE_VEBL_TREE
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO(2, boost::vebl_d_ary_tree_storage<2>)
      
      RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO(4, boost::vebl_d_ary_tree_storage<4>)
      
      RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
        
#endif
      
    };
      
#endif
    
  };
  
#undef RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
#undef RK_SBASTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO
#undef RK_SBASTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO
#undef RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE
#undef RK_SBASTAR_PLANNER_MAKE_SBASTAR_BUNDLE_WITH_GOAL
#undef RK_SBASTAR_PLANNER_CALL_SBASTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_LAZY_SBASTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_LAZY_BNB_SBASTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_SBARRTSTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_LAZY_SBARRTSTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_LAZY_BNB_SBARRTSTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_ANYTIME_SBASTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBASTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_BNB_SBASTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_ANYTIME_SBARRTSTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBARRTSTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_BNB_SBARRTSTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
  
};


};

};

#endif


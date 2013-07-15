/**
 * \file rrtstar_path_planner.hpp
 * 
 * This library defines a class to solve path planning problems using the 
 * Rapidly-exploring Random Tree Star (RRT*) algorithm (or one of its variants). 
 * Given a C_free (configuration space restricted to non-colliding points) and a 
 * result reporting policy, this class will probabilistically construct a motion-graph 
 * that will connect a starting point and a goal point with a path through C-free 
 * that is as close as possible to the optimal path in terms of distance.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date August 2012
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

#ifndef REAK_RRTSTAR_PATH_PLANNER_HPP
#define REAK_RRTSTAR_PATH_PLANNER_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "motion_planner_base.hpp"

#include "metric_space_concept.hpp"
#include "any_sbmp_reporter.hpp"

#include "graph_alg/rrt_star.hpp"

#include "graph_alg/d_ary_bf_tree.hpp"
#include "graph_alg/d_ary_cob_tree.hpp"
#include "graph_alg/bgl_tree_adaptor.hpp"
#include "graph_alg/bgl_more_property_maps.hpp"
#include "graph_alg/pooled_adjacency_list.hpp"
#include "dvp_layout_adjacency_list.hpp"
#include "metric_space_search.hpp"
#include "topological_search.hpp"

#include "path_planner_options.hpp"
#include "graph_alg/neighborhood_functors.hpp"
#include "any_motion_graphs.hpp"
#include "planning_visitors.hpp"

namespace ReaK {
  
namespace pp {
  



/**
 * This class solves path planning problems using the 
 * Rapidly-exploring Random Tree Star (RRT*) algorithm (or one of its variants). 
 * Given a C_free (configuration space restricted to non-colliding points) and a 
 * result reporting policy, this class will probabilistically construct a motion-graph 
 * that will connect a starting point and a goal point with a path through C-free 
 * that is as close as possible to the optimal path in terms of distance.
 * \tparam FreeSpaceType The topology type on which to perform the planning, should be the C-free sub-space of a larger configuration space.
 * \tparam SBPPReporter The reporter type to use to report the progress of the path-planning.
 */
template <typename FreeSpaceType>
class rrtstar_planner : public sample_based_planner<FreeSpaceType> {
  public:
    typedef sample_based_planner<FreeSpaceType> base_type;
    typedef rrtstar_planner<FreeSpaceType> self;
    
    typedef FreeSpaceType space_type;
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    
    BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));
    
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename topology_traits< super_space_type >::point_difference_type point_difference_type;
    
  public:
    
    /**
     * This function computes a valid path in the C-free. If it cannot 
     * achieve a valid path, an exception will be thrown. This algorithmic
     * path solver class is such that any settings that ought to be set for the 
     * path planning algorithm should be set before calling this function, otherwise
     * the function is likely to fail.
     * \param aQuery The query object that defines as input the parameters of the query, 
     *               and as output, the recorded solutions.
     */
    virtual void solve_planning_query(planning_query<FreeSpaceType>& aQuery);
    
    /**
     * Parametrized constructor.
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     * \param aMaxVertexCount The maximum number of samples to generate during the motion planning.
     * \param aProgressInterval The number of new samples between each "progress report".
     * \param aDataStructureFlags An integer flags representing the kind of motion graph data-structure to use in the 
     *                            planning algorithm. Can be ADJ_LIST_MOTION_GRAPH or DVP_ADJ_LIST_MOTION_GRAPH.
     *                            Any combination of those two and of KNN method flags to use for nearest
     *                            neighbor queries in the graph. KNN method flags can be LINEAR_SEARCH_KNN, 
     *                            DVP_BF2_TREE_KNN, DVP_BF4_TREE_KNN, DVP_COB2_TREE_KNN, or DVP_COB4_TREE_KNN.
     *                            See path_planner_options.hpp documentation.
     * \param aPlanningMethodFlags The integer flags that identify various options to use with this planner.
     *                             The options available include EAGER_COLLISION_CHECKING or LAZY_COLLISION_CHECKING, 
     *                             NOMINAL_PLANNER_ONLY or any combination of PLAN_WITH_VORONOI_PULL, 
     *                             PLAN_WITH_NARROW_PASSAGE_PUSH and PLAN_WITH_ANYTIME_HEURISTIC, UNIDIRECTIONAL_PLANNING 
     *                             or BIDIRECTIONAL_PLANNING, and USE_BRANCH_AND_BOUND_PRUNING_FLAG. 
     * \param aSteerProgressTolerance The steer progress tolerance to be used by this planner when making connections.
     * \param aConnectionTolerance The connection tolerance to be used by this planner when making connections.
     * \param aSpaceDimensionality The dimensionality of the space used by this planner.
     * \param aReporter The path-planning reporter to be used by this planner.
     */
    rrtstar_planner(const shared_ptr< space_type >& aWorld = shared_ptr< space_type >(), 
                    std::size_t aMaxVertexCount = 5000, 
                    std::size_t aProgressInterval = 100,
                    std::size_t aDataStructureFlags = ADJ_LIST_MOTION_GRAPH | DVP_BF2_TREE_KNN,
                    std::size_t aPlanningMethodFlags = BIDIRECTIONAL_PLANNING,
                    double aSteerProgressTolerance = 0.1,
                    double aConnectionTolerance = 0.1,
                    std::size_t aSpaceDimensionality = 1,
                    const any_sbmp_reporter_chain<space_type>& aReporter = any_sbmp_reporter_chain<space_type>()) :
                    base_type("rrtstar_planner", aWorld, aMaxVertexCount, aProgressInterval,
                              aDataStructureFlags, aPlanningMethodFlags,
                              aSteerProgressTolerance, aConnectionTolerance, 1.0, 
                              aSpaceDimensionality, aReporter) { };
    
    virtual ~rrtstar_planner() { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2460009,1,"rrtstar_planner",base_type)
};


template <typename FreeSpaceType>
void rrtstar_planner<FreeSpaceType>::solve_planning_query(planning_query<FreeSpaceType>& aQuery) {
  
  this->reset_internal_state();
  
  typedef typename subspace_traits<FreeSpaceType>::super_space_type SuperSpace;
  typedef typename topology_traits<SuperSpace>::point_type PointType;
  
  typedef optimal_mg_vertex<FreeSpaceType> VertexProp;
  typedef optimal_mg_edge<FreeSpaceType> EdgeProp;
  
  typedef boost::data_member_property_map<PointType, VertexProp > PositionMap;
  PositionMap pos_map = PositionMap(&VertexProp::position);
  
  typedef boost::data_member_property_map<double, VertexProp > CostMap;
  CostMap cost_map = CostMap(&VertexProp::distance_accum);
  
  typedef boost::data_member_property_map<std::size_t, VertexProp > PredMap;
  PredMap pred_map = PredMap(&VertexProp::predecessor);
  
  typedef boost::data_member_property_map<double, EdgeProp > WeightMap;
  WeightMap weight_map = WeightMap(&EdgeProp::weight);
  
  double space_dim = double( this->get_space_dimensionality() );
  double space_Lc = aQuery.get_heuristic_to_goal( aQuery.get_start_position() );
  
  shared_ptr<const SuperSpace> sup_space_ptr(&(this->m_space->get_super_space()),null_deleter());
  
  planning_visitor<FreeSpaceType> vis(this, &aQuery);
  
  VertexProp vp_start;
  vp_start.position = aQuery.get_start_position();
  
  path_planning_p2p_query<FreeSpaceType>* p2p_query_ptr = reinterpret_cast< path_planning_p2p_query<FreeSpaceType>* >(aQuery.castTo(path_planning_p2p_query<FreeSpaceType>::getStaticObjectType()));
  
  
  
#define RK_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE \
  Vertex vs = add_vertex(vp_start, motion_graph); \
  motion_graph[vs].distance_accum = 0.0; \
  motion_graph[vs].predecessor = vs; \
  vis.m_start_node = boost::any( vs ); \
  if( p2p_query_ptr ) { \
    VertexProp vp_goal; \
    vp_goal.position = p2p_query_ptr->goal_pos; \
    Vertex vg = add_vertex(vp_goal, motion_graph); \
    motion_graph[vg].distance_accum = std::numeric_limits<double>::infinity(); \
    motion_graph[vg].predecessor = vg; \
    vis.m_goal_node = boost::any( vg ); \
  };
  
#define RK_RRTSTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(ARITY, TREE_STORAGE) \
  typedef boost::property_map< MotionGraphType, PointType VertexProp::* >::type GraphPositionMap; \
  typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, ARITY, random_vp_chooser, TREE_STORAGE > SpacePartType; \
  SpacePartType space_part(motion_graph, sup_space_ptr, get(&VertexProp::position, motion_graph)); \
   \
  typedef multi_dvp_tree_search<MotionGraphType, SpacePartType> NNFinderType; \
  NNFinderType nn_finder; \
  nn_finder.graph_tree_map[&motion_graph] = &space_part; \
   \
  type_erased_knn_synchro< MotionGraphType, NNFinderType > NN_synchro(nn_finder); \
  vis.m_nn_synchro = &NN_synchro;
  
#define RK_RRTSTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO(ARITY, TREE_STORAGE) \
  typedef dvp_adjacency_list< \
    VertexProp, EdgeProp, SuperSpace, PositionMap, \
    ARITY, random_vp_chooser, TREE_STORAGE, \
    boost::vecS, boost::bidirectionalS, boost::listS > ALTGraph; \
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
  vis.m_nn_synchro = &NN_synchro;
  
  
#define RK_RRTSTAR_PLANNER_CALL_RRTSTAR_FUNCTION \
  ReaK::graph::generate_rrt_star( \
    motion_graph, \
    *sup_space_ptr, \
    vis, \
    pos_map, \
    cost_map, \
    pred_map, \
    weight_map, \
    get(random_sampler, *sup_space_ptr), \
    ReaK::graph::star_neighborhood< NNFinderType >( \
      nn_finder, \
      space_dim, 3.0 * space_Lc), \
    this->get_max_vertex_count());

#define RK_RRTSTAR_PLANNER_CALL_RRTSTAR_BNB_FUNCTION \
  ReaK::graph::generate_bnb_rrt_star( \
    motion_graph, \
    boost::any_cast<Vertex>( vis.m_start_node ), \
    boost::any_cast<Vertex>( vis.m_goal_node ), \
    *sup_space_ptr, \
    vis, \
    pos_map, \
    cost_map, \
    pred_map, \
    weight_map, \
    get(random_sampler, *sup_space_ptr), \
    ReaK::graph::star_neighborhood< NNFinderType >( \
      nn_finder, \
      space_dim, 3.0 * space_Lc), \
    this->get_max_vertex_count());
  
  
#define RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION \
  if(this->m_planning_method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) { \
    RK_RRTSTAR_PLANNER_CALL_RRTSTAR_BNB_FUNCTION \
  } else { /* assume nominal method only. */ \
    RK_RRTSTAR_PLANNER_CALL_RRTSTAR_FUNCTION \
  };
  
  
  if((this->m_planning_method_flags & PLANNING_DIRECTIONALITY_MASK) == UNIDIRECTIONAL_PLANNING) {
    
    if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == ADJ_LIST_MOTION_GRAPH) {
      
      typedef boost::pooled_adjacency_list< 
        boost::undirectedS, VertexProp, EdgeProp,
        boost::no_property, boost::listS> MotionGraphType;
      
      typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex;
      
      MotionGraphType motion_graph;
      
      RK_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE
      
      if((this->m_data_structure_flags & KNN_METHOD_MASK) == LINEAR_SEARCH_KNN) {
        
        typedef linear_neighbor_search<> NNFinderType;
        NNFinderType nn_finder;
        
        any_knn_synchro NN_synchro;
        vis.m_nn_synchro = &NN_synchro;
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, graph::d_ary_bf_tree_storage<2>)
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, graph::d_ary_bf_tree_storage<4>)
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, graph::d_ary_cob_tree_storage<2>)
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, graph::d_ary_cob_tree_storage<4>)
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      };
      
    } else if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == DVP_ADJ_LIST_MOTION_GRAPH) {
      
      if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO(2, graph::d_ary_bf_tree_storage<2>)
        
        RK_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO(4, graph::d_ary_bf_tree_storage<4>)
        
        RK_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO(2, graph::d_ary_cob_tree_storage<2>)
        
        RK_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO(4, graph::d_ary_cob_tree_storage<4>)
        
        RK_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      };
      
    };
    
    
#undef RK_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE
#undef RK_RRTSTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO
#undef RK_RRTSTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO
#undef RK_RRTSTAR_PLANNER_CALL_RRTSTAR_FUNCTION
#undef RK_RRTSTAR_PLANNER_CALL_RRTSTAR_BNB_FUNCTION
#undef RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
    
    
  } else {
#if 0    
    
    if(p2p_query_ptr == NULL)
      return;
    
    
#define RK_RRTSTAR_PLANNER_INIT_BIDIR_START_AND_GOAL_NODE \
  Vertex vs = add_vertex(vp_start, motion_graph1); \
  motion_graph1[vs].distance_accum = 0.0; \
  motion_graph1[vs].predecessor = vs; \
  vis.m_start_node = boost::any( vs ); \
  VertexProp vp_goal; \
  vp_goal.position = p2p_query_ptr->goal_pos; \
  Vertex vg = add_vertex(vp_goal, motion_graph2); \
  motion_graph2[vg].distance_accum = std::numeric_limits<double>::infinity(); \
  motion_graph2[vg].predecessor = vg; \
  vis.m_goal_node = boost::any( vg );
    
    
#define RK_RRTSTAR_PLANNER_SETUP_BIDIR_DVP_TREE_SYNCHRO(ARITY, TREE_STORAGE) \
  typedef boost::property_map< MotionGraphType, PointType VertexProp::* >::type GraphPositionMap; \
  typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, ARITY, random_vp_chooser, TREE_STORAGE > SpacePartType; \
   \
  SpacePartType space_part1(motion_graph1, sup_space_ptr, get(&VertexProp::position, motion_graph1)); \
  SpacePartType space_part2(motion_graph2, sup_space_ptr, get(&VertexProp::position, motion_graph2)); \
   \
  typedef multi_dvp_tree_search<MotionGraphType, SpacePartType> NNFinderType; \
  NNFinderType nn_finder; \
  nn_finder.graph_tree_map[&motion_graph1] = &space_part1; \
  nn_finder.graph_tree_map[&motion_graph2] = &space_part2; \
   \
  type_erased_knn_synchro< MotionGraphType, NNFinderType > NN_synchro(nn_finder); \
  vis.m_nn_synchro = &NN_synchro;
  
  
#define RK_RRTSTAR_PLANNER_SETUP_BIDIR_ALT_TREE_SYNCHRO(ARITY, TREE_STORAGE) \
  typedef dvp_adjacency_list< \
    VertexProp, EdgeProp, SuperSpace, PositionMap, \
    ARITY, random_vp_chooser, TREE_STORAGE, \
    boost::vecS, boost::bidirectionalS, boost::listS > ALTGraph; \
  typedef typename ALTGraph::adj_list_type MotionGraphType; \
  typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex; \
   \
  ALTGraph space_part1(sup_space_ptr, pos_map); \
  ALTGraph space_part2(sup_space_ptr, pos_map); \
   \
  MotionGraphType motion_graph1 = space_part1.get_adjacency_list(); \
  MotionGraphType motion_graph2 = space_part2.get_adjacency_list(); \
   \
  typedef multi_dvp_tree_search<MotionGraphType, ALTGraph> NNFinderType; \
  NNFinderType nn_finder; \
  nn_finder.graph_tree_map[&motion_graph1] = &space_part1; \
  nn_finder.graph_tree_map[&motion_graph2] = &space_part2; \
   \
  any_knn_synchro NN_synchro; \
  vis.m_nn_synchro = &NN_synchro;
  
  
#define RK_RRTSTAR_PLANNER_CALL_BIDIR_RRTSTAR_FUNCTION \
  ReaK::graph::generate_bidir_rrt_star( \
    motion_graph1, motion_graph2, \
    *sup_space_ptr, \
    vis, \
    pos_map, \
    cost_map, \
    pred_map, \
    weight_map, \
    get(random_sampler, *sup_space_ptr), \
    ReaK::graph::star_neighborhood< NNFinderType >( \
      nn_finder, \
      space_dim, 3.0 * space_Lc), \
    this->get_max_vertex_count());
  
  
    
    if((m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == ADJ_LIST_MOTION_GRAPH) {
      
      typedef boost::adjacency_list< boost::vecS, boost::listS, boost::bidirectionalS,
                             VertexProp,
                             EdgeProp,
                             boost::vecS> MotionGraphType;
      typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex;
      
      typedef boost::composite_property_map< 
        PositionMap, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t > > GraphPositionMap;
      
      MotionGraphType motion_graph1;
      MotionGraphType motion_graph2;
      
      RK_RRTSTAR_PLANNER_INIT_BIDIR_START_AND_GOAL_NODE
      
      if(m_knn_flag == LINEAR_SEARCH_KNN) {
        
        typedef linear_neighbor_search<> NNFinderType;
        NNFinderType nn_finder;
        
        any_knn_synchro NN_synchro;
        vis.m_nn_synchro = &NN_synchro;
        
        RK_RRTSTAR_PLANNER_CALL_BIDIR_RRTSTAR_FUNCTION
        
      } else if(m_knn_flag == DVP_BF2_TREE_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_BIDIR_DVP_TREE_SYNCHRO(2, graph::d_ary_bf_tree_storage<2>)
        
        RK_RRTSTAR_PLANNER_CALL_BIDIR_RRTSTAR_FUNCTION
        
      } else if(m_knn_flag == DVP_BF4_TREE_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_BIDIR_DVP_TREE_SYNCHRO(4, graph::d_ary_bf_tree_storage<4>)
        
        RK_RRTSTAR_PLANNER_CALL_BIDIR_RRTSTAR_FUNCTION
        
      } else if(m_knn_flag == DVP_COB2_TREE_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_BIDIR_DVP_TREE_SYNCHRO(2, graph::d_ary_cob_tree_storage<2>)
        
        RK_RRTSTAR_PLANNER_CALL_BIDIR_RRTSTAR_FUNCTION
        
      } else if(m_knn_flag == DVP_COB4_TREE_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_BIDIR_DVP_TREE_SYNCHRO(4, graph::d_ary_cob_tree_storage<4>)
        
        RK_RRTSTAR_PLANNER_CALL_BIDIR_RRTSTAR_FUNCTION
        
      };
      
    } else if((m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == DVP_ADJ_LIST_MOTION_GRAPH) {
      
      if(m_knn_flag == DVP_ALT_BF2_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_BIDIR_ALT_TREE_SYNCHRO(2, graph::d_ary_bf_tree_storage<2>)
        
        RK_RRTSTAR_PLANNER_INIT_BIDIR_START_AND_GOAL_NODE
        
        RK_RRTSTAR_PLANNER_CALL_BIDIR_RRTSTAR_FUNCTION
        
      } else if(m_knn_flag == DVP_ALT_BF4_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_BIDIR_ALT_TREE_SYNCHRO(4, graph::d_ary_bf_tree_storage<4>)
        
        RK_RRTSTAR_PLANNER_INIT_BIDIR_START_AND_GOAL_NODE
        
        RK_RRTSTAR_PLANNER_CALL_BIDIR_RRTSTAR_FUNCTION
        
      } else if(m_knn_flag == DVP_ALT_COB2_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_BIDIR_ALT_TREE_SYNCHRO(2, graph::d_ary_cob_tree_storage<2>)
        
        RK_RRTSTAR_PLANNER_INIT_BIDIR_START_AND_GOAL_NODE
        
        RK_RRTSTAR_PLANNER_CALL_BIDIR_RRTSTAR_FUNCTION
        
      } else if(m_knn_flag == DVP_ALT_COB4_KNN) {
        
        RK_RRTSTAR_PLANNER_SETUP_BIDIR_ALT_TREE_SYNCHRO(4, graph::d_ary_cob_tree_storage<4>)
        
        RK_RRTSTAR_PLANNER_INIT_BIDIR_START_AND_GOAL_NODE
        
        RK_RRTSTAR_PLANNER_CALL_BIDIR_RRTSTAR_FUNCTION
        
      };
      
    };
    
#undef RK_RRTSTAR_PLANNER_INIT_BIDIR_START_AND_GOAL_NODE
#undef RK_RRTSTAR_PLANNER_SETUP_BIDIR_DVP_TREE_SYNCHRO
#undef RK_RRTSTAR_PLANNER_SETUP_BIDIR_ALT_TREE_SYNCHRO
#undef RK_RRTSTAR_PLANNER_CALL_BIDIR_RRTSTAR_FUNCTION
    
    
#endif
  };
  
  
};


};

};

#endif


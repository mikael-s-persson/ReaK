/**
 * \file rrt_path_planner.hpp
 * 
 * This library defines a class to solve path planning problems using the 
 * Rapidly-exploring Random Tree (RRT) algorithm (or one of its variants). 
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

#ifndef REAK_RRT_PATH_PLANNER_HPP
#define REAK_RRT_PATH_PLANNER_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "motion_planner_base.hpp"

#include "metric_space_concept.hpp"
#include "any_sbmp_reporter.hpp"

#include "graph_alg/rr_tree.hpp"

#include "graph_alg/d_ary_bf_tree.hpp"
#include "graph_alg/d_ary_cob_tree.hpp"
#include "graph_alg/bgl_tree_adaptor.hpp"
#include "graph_alg/bgl_more_property_maps.hpp"
#include "dvp_layout_adjacency_list.hpp"
#include "metric_space_search.hpp"
#include "topological_search.hpp"

#include "path_planner_options.hpp"
#include "any_motion_graphs.hpp"
#include "planning_visitors.hpp"


namespace ReaK {
  
namespace pp {


/**
 * This class solves path planning problems using the 
 * Rapidly-exploring Random Tree (RRT) algorithm (or one of its variants). 
 * Given a C_free (configuration space restricted to non-colliding points) and a 
 * result reporting policy, this class will probabilistically construct a motion-graph 
 * that will connect a starting point and a goal point with a path through C-free 
 * that is as close as possible to the optimal path in terms of distance.
 * \tparam FreeSpaceType The topology type on which to perform the planning, should be the C-free sub-space of a larger configuration space.
 */
template <typename FreeSpaceType>
class rrt_planner : public sample_based_planner<FreeSpaceType> {
  public:
    typedef sample_based_planner<FreeSpaceType> base_type;
    typedef rrt_planner<FreeSpaceType> self;
    
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
     * \param aReporter The path-planning reporter to be used by this planner.
     */
    rrt_planner(const shared_ptr< space_type >& aWorld = shared_ptr< space_type >(), 
                std::size_t aMaxVertexCount = 5000, 
                std::size_t aProgressInterval = 100,
                std::size_t aDataStructureFlags = ADJ_LIST_MOTION_GRAPH | DVP_BF2_TREE_KNN,
                std::size_t aPlanningMethodFlags = BIDIRECTIONAL_PLANNING,
                double aSteerProgressTolerance = 0.1,
                double aConnectionTolerance = 0.1,
                const any_sbmp_reporter_chain<space_type>& aReporter = any_sbmp_reporter_chain<space_type>()) :
                base_type("rrt_planner", aWorld, aMaxVertexCount, aProgressInterval,
                          aDataStructureFlags, aPlanningMethodFlags,
                          aSteerProgressTolerance, aConnectionTolerance, 1.0, 1, aReporter) { };
    
    virtual ~rrt_planner() { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2460007,1,"rrt_planner",base_type)
};



template <typename FreeSpaceType>
void rrt_planner<FreeSpaceType>::solve_planning_query(planning_query<FreeSpaceType>& aQuery) {
  
  this->reset_internal_state();
  
  typedef typename subspace_traits<FreeSpaceType>::super_space_type SuperSpace;
  typedef typename topology_traits<SuperSpace>::point_type PointType;
  
  typedef mg_vertex_data<FreeSpaceType> VertexProp;
  typedef mg_edge_data<FreeSpaceType> EdgeProp;
  
  typedef boost::data_member_property_map<PointType, VertexProp > PositionMap;
  PositionMap pos_map = PositionMap(&VertexProp::position);
  
  shared_ptr<const SuperSpace> sup_space_ptr(&(this->m_space->get_super_space()),null_deleter());
  
  planning_visitor<FreeSpaceType> vis(this, &aQuery);
  
  VertexProp vp_start;
  vp_start.position = aQuery.get_start_position();
  
  if((this->m_planning_method_flags & PLANNING_DIRECTIONALITY_MASK) == UNIDIRECTIONAL_PLANNING) {
    
    
#define RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(ARITY, TREE_STORAGE) \
  typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex; \
  typedef typename boost::property_map< MotionGraphType, PointType VertexProp::* >::type GraphPositionMap; \
  typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, ARITY, random_vp_chooser, TREE_STORAGE > SpacePartType; \
  SpacePartType space_part(motion_graph, sup_space_ptr, get(&VertexProp::position, motion_graph)); \
   \
  typedef multi_dvp_tree_search<MotionGraphType, SpacePartType> NNFinderType; \
  NNFinderType nn_finder; \
  nn_finder.graph_tree_map[&motion_graph] = &space_part; \
   \
  type_erased_knn_synchro< MotionGraphType, NNFinderType > NN_synchro(nn_finder); \
  vis.m_nn_synchro = &NN_synchro;
  
#define RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(ARITY, TREE_STORAGE) \
  typedef dvp_adjacency_list< \
    VertexProp, EdgeProp, SuperSpace, PositionMap, \
    ARITY, random_vp_chooser, TREE_STORAGE, \
    boost::vecS, boost::bidirectionalS, boost::listS > ALTGraph; \
  typedef typename ALTGraph::adj_list_type MotionGraphType; \
   \
  ALTGraph space_part(sup_space_ptr, pos_map); \
  MotionGraphType motion_graph = space_part.get_adjacency_list(); \
  vis.m_start_node = boost::any( create_root(vp_start, motion_graph) ); \
   \
  typedef multi_dvp_tree_search<MotionGraphType, ALTGraph> NNFinderType; \
  NNFinderType nn_finder; \
  nn_finder.graph_tree_map[&motion_graph] = &space_part; \
   \
  any_knn_synchro NN_synchro; \
  vis.m_nn_synchro = &NN_synchro;
  
  
#define RK_RRT_PLANNER_CALL_RRT_FUNCTION \
  ReaK::graph::generate_rrt( \
    motion_graph, *sup_space_ptr, \
    vis, pos_map, get(random_sampler, *sup_space_ptr), \
    nn_finder);
    
    
    if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == ADJ_LIST_MOTION_GRAPH) {
      
      typedef boost::adjacency_list< 
        boost::vecS, boost::listS, boost::bidirectionalS,
        VertexProp, EdgeProp, boost::vecS> MotionGraphType;
      
      MotionGraphType motion_graph;
      vis.m_start_node = boost::any( create_root(vp_start, motion_graph) );
      
      if((this->m_data_structure_flags & KNN_METHOD_MASK) == LINEAR_SEARCH_KNN) {
        
        any_knn_synchro NN_synchro;
        vis.m_nn_synchro = &NN_synchro;
        linear_neighbor_search<> nn_finder;
        
        RK_RRT_PLANNER_CALL_RRT_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, graph::d_ary_bf_tree_storage<2>)
        
        RK_RRT_PLANNER_CALL_RRT_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, graph::d_ary_bf_tree_storage<4>)
        
        RK_RRT_PLANNER_CALL_RRT_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, graph::d_ary_cob_tree_storage<2>)
        
        RK_RRT_PLANNER_CALL_RRT_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, graph::d_ary_cob_tree_storage<4>)
        
        RK_RRT_PLANNER_CALL_RRT_FUNCTION
        
      };
      
    } else if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == DVP_ADJ_LIST_MOTION_GRAPH) {
      
      if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(2, graph::d_ary_bf_tree_storage<2>)
        
        RK_RRT_PLANNER_CALL_RRT_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(4, graph::d_ary_bf_tree_storage<4>)
        
        RK_RRT_PLANNER_CALL_RRT_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(2, graph::d_ary_cob_tree_storage<2>)
        
        RK_RRT_PLANNER_CALL_RRT_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(4, graph::d_ary_cob_tree_storage<4>)
        
        RK_RRT_PLANNER_CALL_RRT_FUNCTION
        
      };
      
    };
    
#undef RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO
#undef RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO
#undef RK_RRT_PLANNER_CALL_RRT_FUNCTION
    
  } else {
    path_planning_p2p_query<FreeSpaceType>* p2p_query_ptr = reinterpret_cast< path_planning_p2p_query<FreeSpaceType>* >(aQuery.castTo(path_planning_p2p_query<FreeSpaceType>::getStaticObjectType()));
    if(p2p_query_ptr == NULL)
      return;
    
    VertexProp vp_goal;
    vp_goal.position = p2p_query_ptr->goal_pos;
    
    
#define RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(ARITY, TREE_STORAGE) \
  typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex; \
  typedef typename boost::property_map< MotionGraphType, PointType VertexProp::* >::type GraphPositionMap; \
  typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, ARITY, random_vp_chooser, TREE_STORAGE > SpacePartType; \
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
  
  
#define RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(ARITY, TREE_STORAGE) \
  typedef dvp_adjacency_list< \
    VertexProp, EdgeProp, SuperSpace, PositionMap, \
    ARITY, random_vp_chooser, TREE_STORAGE, \
    boost::vecS, boost::bidirectionalS, boost::listS > ALTGraph; \
  typedef typename ALTGraph::adj_list_type MotionGraphType; \
   \
  ALTGraph space_part1(sup_space_ptr, pos_map); \
  ALTGraph space_part2(sup_space_ptr, pos_map); \
   \
  MotionGraphType motion_graph1 = space_part1.get_adjacency_list(); \
  MotionGraphType motion_graph2 = space_part2.get_adjacency_list(); \
   \
  vis.m_start_node = boost::any( create_root(vp_start, motion_graph1) ); \
  vis.m_goal_node  = boost::any( create_root(vp_goal,  motion_graph2) ); \
   \
  typedef multi_dvp_tree_search<MotionGraphType, ALTGraph> NNFinderType; \
  NNFinderType nn_finder; \
  nn_finder.graph_tree_map[&motion_graph1] = &space_part1; \
  nn_finder.graph_tree_map[&motion_graph2] = &space_part2; \
   \
  any_knn_synchro NN_synchro; \
  vis.m_nn_synchro = &NN_synchro;
  
  
#define RK_RRT_PLANNER_CALL_BIRRT_FUNCTION \
  ReaK::graph::generate_bidirectional_rrt( \
    motion_graph1, motion_graph2, *sup_space_ptr, \
    vis, pos_map, get(random_sampler, *sup_space_ptr), \
    nn_finder);
    
    
    if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == ADJ_LIST_MOTION_GRAPH) {
      
      typedef boost::adjacency_list< 
        boost::vecS, boost::listS, boost::bidirectionalS,
        VertexProp, EdgeProp, boost::vecS> MotionGraphType;
      
      MotionGraphType motion_graph1;
      MotionGraphType motion_graph2;
      
      vis.m_start_node = boost::any( create_root(vp_start, motion_graph1) );
      vis.m_goal_node  = boost::any( create_root(vp_goal,  motion_graph2) );
      
      if((this->m_data_structure_flags & KNN_METHOD_MASK) == LINEAR_SEARCH_KNN) {
        
        any_knn_synchro NN_synchro;
        vis.m_nn_synchro = &NN_synchro;
        linear_neighbor_search<> nn_finder;
        
        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(2, graph::d_ary_bf_tree_storage<2>)
        
        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(4, graph::d_ary_bf_tree_storage<4>)
        
        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(2, graph::d_ary_cob_tree_storage<2>)
        
        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(4, graph::d_ary_cob_tree_storage<4>)
        
        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION
        
      };
      
    } else if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == DVP_ADJ_LIST_MOTION_GRAPH) {
      
      if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(2, graph::d_ary_bf_tree_storage<2>)
        
        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(4, graph::d_ary_bf_tree_storage<4>)
        
        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(2, graph::d_ary_cob_tree_storage<2>)
        
        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
        
        RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(4, graph::d_ary_cob_tree_storage<4>)
        
        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION
        
      };
      
    };
    
#undef RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO
#undef RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO
#undef RK_RRT_PLANNER_CALL_BIRRT_FUNCTION
    
  };
  
  
  
};


};

};

#endif


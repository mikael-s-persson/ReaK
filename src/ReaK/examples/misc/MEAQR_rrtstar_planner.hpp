/**
 * \file MEAQR_rrtstar_planner.hpp
 * 
 * This library defines a class
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

#ifndef REAK_MEAQR_RRTSTAR_PLANNER_HPP
#define REAK_MEAQR_RRTSTAR_PLANNER_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "path_planning/motion_planner_base.hpp"

#include "path_planning/metric_space_concept.hpp"
#include "path_planning/any_sbmp_reporter.hpp"

#include "graph_alg/rrt_star.hpp"

#include "path_planning/motion_graph_structures.hpp"


// BGL-Extra includes:
#include <boost/graph/tree_adaptor.hpp>
#include <boost/graph/more_property_tags.hpp>
#include <boost/graph/more_property_maps.hpp>


#include "path_planning/metric_space_search.hpp"
#include "path_planning/topological_search.hpp"

#include "path_planning/p2p_planning_query.hpp"
#include "path_planning/path_planner_options.hpp"
#include "graph_alg/neighborhood_functors.hpp"
#include "path_planning/any_motion_graphs.hpp"
#include "path_planning/planning_visitors.hpp"

#include "MEAQR_topology.hpp"
#include "topologies/fixed_topology_random_sampler.hpp"

namespace ReaK {
  
namespace pp {


/**
 * This class is a RRT*-based path-planner over the a MEAQR-controlled system-topology.
 * \tparam StateSpace The topology type of state-space of the dynamic system under MEAQR control.
 * \tparam StateSpaceSystem The type of the dynamic system under MEAQR control.
 * \tparam SBPPReporter The reporter type to use to report the progress of the path-planning.
 */
template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
class MEAQR_rrtstar_planner : public sample_based_planner< MEAQR_topology_with_CD<StateSpace, StateSpaceSystem, StateSpaceSampler> > {
  public:
    typedef MEAQR_topology_with_CD<StateSpace, StateSpaceSystem, StateSpaceSampler> space_type;
    typedef typename subspace_traits<space_type>::super_space_type super_space_type;
    
    typedef sample_based_planner< space_type > base_type;
    typedef MEAQR_rrtstar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler> self;
    
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename topology_traits< super_space_type >::point_difference_type point_difference_type;
    
    typedef typename steerable_space_traits< super_space_type >::steer_record_type steer_record_type;
    
  protected:
    
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
    virtual void solve_planning_query(planning_query<space_type>& aQuery);
    
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
    MEAQR_rrtstar_planner(const shared_ptr< space_type >& aWorld = shared_ptr< space_type >(), 
                          std::size_t aMaxVertexCount = 5000, 
                          std::size_t aProgressInterval = 100,
                          std::size_t aDataStructureFlags = ADJ_LIST_MOTION_GRAPH | DVP_BF2_TREE_KNN,
                          std::size_t aPlanningMethodFlags = BIDIRECTIONAL_PLANNING,
                          double aSteerProgressTolerance = 0.1,
                          double aConnectionTolerance = 0.1,
                          std::size_t aSpaceDimensionality = 1,
                          const any_sbmp_reporter_chain<space_type>& aReporter = any_sbmp_reporter_chain<space_type>()) :
                          base_type("MEAQR_rrtstar_planner", aWorld, aMaxVertexCount, aProgressInterval,
                                    aDataStructureFlags, aPlanningMethodFlags,
                                    aSteerProgressTolerance, aConnectionTolerance, 1.0, 
                                    aSpaceDimensionality, aReporter) { };
    
    virtual ~MEAQR_rrtstar_planner() { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
    };
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC246000D,1,"MEAQR_rrtstar_planner",base_type)
};




template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
struct MEAQR_rrtstar_visitor : planning_visitor< MEAQR_topology_with_CD<StateSpace, StateSpaceSystem, StateSpaceSampler> > {
  typedef MEAQR_topology_with_CD<StateSpace, StateSpaceSystem, StateSpaceSampler> space_type;
  typedef planning_visitor< space_type > base_type;
  
  typedef typename base_type::planner_base_type planner_base_type;
  typedef typename base_type::query_type query_type;
  
  MEAQR_rrtstar_visitor(planner_base_type* aPlanner,
                        query_type* aQuery = NULL,
                        any_knn_synchro* aNNSynchro = NULL,
                        boost::any aStartNode = boost::any(),
                        boost::any aGoalNode = boost::any()) :
                        base_type(aPlanner, aQuery, aNNSynchro, aStartNode, aGoalNode) { };
  
  typedef typename topology_traits< space_type >::point_type point_type;
  typedef optimal_mg_vertex<space_type> VertexProp;
  typedef optimal_mg_edge<space_type> EdgeProp;
  
  template <typename Vertex, typename Graph>
  boost::tuple<point_type, bool, EdgeProp> steer_towards_position(const point_type& p, Vertex u, Graph& g) const {
    typedef typename EdgeProp::steer_record_type SteerRec;
    typedef boost::tuple<point_type, bool, EdgeProp> ResultType;
    
    // First, try to bring the state-space point within the time-horizon:
    double total_dist = get(distance_metric, this->m_query->space->get_super_space())(g[u].position, p, this->m_query->space->get_super_space());
    double max_cost_to_go = 0.75 * this->m_query->space->get_max_time_horizon() * this->m_query->space->get_idle_power_cost(g[u].position);
    point_type p_dest = p;
    while(total_dist > max_cost_to_go) {
      p_dest = point_type( this->m_query->space->get_state_space().move_position_toward(g[u].position.x, 0.5, p_dest.x) );
      total_dist = get(distance_metric, this->m_query->space->get_super_space())(g[u].position, p_dest, this->m_query->space->get_super_space());
    };
    
    // Then, steer to that point, recording the path of the steering function.
    std::pair<point_type, SteerRec> steer_result = this->m_query->space->steer_position_toward(g[u].position, 0.8, p_dest);
    
    // Check if the progress in the state-space was significant (at least 0.1 of the best-case).
    double best_case_dist = get(distance_metric, this->m_query->space->get_state_space())(g[u].position.x, p_dest.x, this->m_query->space->get_state_space());
    double actual_dist = get(distance_metric, this->m_query->space->get_state_space())(g[u].position.x, steer_result.first.x, this->m_query->space->get_state_space());
    
    if(actual_dist > this->m_planner->get_steer_progress_tolerance() * best_case_dist) {
//       std::cout << "Steered successfully!" << std::endl;
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      return ResultType(steer_result.first, true, EdgeProp(0.8 * total_dist, std::move(steer_result.second)));
#else
      return ResultType(steer_result.first, true, EdgeProp(0.8 * total_dist, steer_result.second));
#endif
    } else {
      return ResultType(steer_result.first, false, EdgeProp());
    };
  };
  
  template <typename Vertex, typename Graph>
  std::pair<bool,EdgeProp> can_be_connected(Vertex u, Vertex v, const Graph& g) const {
    typedef typename EdgeProp::steer_record_type SteerRec;
    typedef std::pair<bool,EdgeProp> ResultType;
    
    std::pair<point_type, SteerRec> steer_result = this->m_query->space->steer_position_toward(g[u].position, 1.0, g[v].position);
    
    // NOTE Differs from rrtstar_path_planner HERE:
    double best_case_dist = get(distance_metric, this->m_query->space->get_state_space())(g[u].position.x, g[v].position.x, this->m_query->space->get_state_space());
    double diff_dist = get(distance_metric, this->m_query->space->get_state_space())(steer_result.first.x, g[v].position.x, this->m_query->space->get_state_space());
    
    if(diff_dist < this->m_planner->get_connection_tolerance() * best_case_dist) {
//       std::cout << "Connected successfully!" << std::endl;
      return ResultType(true, EdgeProp(
        get(distance_metric, this->m_query->space->get_super_space())(g[u].position, g[v].position, this->m_query->space->get_super_space()),
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        std::move(steer_result.second)
#else
        steer_result.second
#endif
      ));
    } else {
      return ResultType(false,EdgeProp());
    };
  };
  
};





template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
void MEAQR_rrtstar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler>::solve_planning_query(planning_query< typename MEAQR_rrtstar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler>::space_type >& aQuery) {
  
  this->reset_internal_state();
  
  typedef typename MEAQR_rrtstar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler>::space_type FreeSpaceType;
  typedef typename subspace_traits<FreeSpaceType>::super_space_type SuperSpace;
  typedef typename topology_traits<SuperSpace>::point_type PointType;
  
  typedef optimal_mg_vertex<FreeSpaceType> VertexProp;
  typedef optimal_mg_edge<FreeSpaceType> EdgeProp;
  
  typedef mg_vertex_data<FreeSpaceType> BasicVertexProp;
  
  
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
  
  
  typedef fixed_topology_random_sampler< SuperSpace > SamplerType;
  SamplerType get_sample = SamplerType( sup_space_ptr.get() );
  
  
  MEAQR_rrtstar_visitor<StateSpace, StateSpaceSystem, StateSpaceSampler> vis(this, &aQuery);
  
  path_planning_p2p_query<FreeSpaceType>* p2p_query_ptr = reinterpret_cast< path_planning_p2p_query<FreeSpaceType>* >(aQuery.castTo(path_planning_p2p_query<FreeSpaceType>::getStaticObjectType()));
  
  typedef boost::adjacency_list_BC< boost::vecBC, boost::vecBC, boost::bidirectionalS, VertexProp, EdgeProp> MotionGraphType;
  typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex;
  
  MotionGraphType motion_graph;
  
  
#define RK_MEAQR_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE \
  VertexProp vp_start; \
  vp_start.position = aQuery.get_start_position(); \
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
  
  
#define RK_MEAQR_RRTSTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(ARITY, TREE_STORAGE) \
    typedef typename boost::property_map< MotionGraphType, PointType BasicVertexProp::* >::type GraphPositionMap; \
    typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, ARITY, \
                     random_vp_chooser, TREE_STORAGE, \
                     no_position_caching_policy > SpacePartType; \
    SpacePartType space_part(motion_graph, sup_space_ptr, get(&BasicVertexProp::position, motion_graph)); \
     \
    typedef multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> NNFinderType; \
    NNFinderType nn_finder; \
    nn_finder.graph_tree_map[&motion_graph] = &space_part; \
     \
    type_erased_knn_synchro< MotionGraphType, NNFinderType > NN_synchro(nn_finder); \
    vis.m_nn_synchro = &NN_synchro;
  
  
#define RK_MEAQR_RRTSTAR_PLANNER_CALL_RRTSTAR_FUNCTION \
    graph::detail::generate_rrt_star_loop( \
      motion_graph, *sup_space_ptr, vis, pos_map, cost_map, pred_map, weight_map, \
      graph::rrg_node_generator< SuperSpace, SamplerType, ReaK::graph::fixed_neighborhood< NNFinderType > >( \
        sup_space_ptr.get(), get_sample,  \
        graph::fixed_neighborhood< NNFinderType >(nn_finder, 5, std::numeric_limits<double>::infinity())), \
      graph::star_neighborhood< NNFinderType >(nn_finder, space_dim, 3.0 * space_Lc));
//       graph::fixed_neighborhood< NNFinderType >(nn_finder, 10, std::numeric_limits<double>::infinity()));
  
  
  
  RK_MEAQR_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE
  
  if((this->m_data_structure_flags & KNN_METHOD_MASK) == LINEAR_SEARCH_KNN) {
    
    typedef linear_pred_succ_search<MotionGraphType> NNFinderType;
    NNFinderType nn_finder;
    
    any_knn_synchro NN_synchro;
    vis.m_nn_synchro = &NN_synchro;
    
    RK_MEAQR_RRTSTAR_PLANNER_CALL_RRTSTAR_FUNCTION
    
  } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
    
    RK_MEAQR_RRTSTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, boost::bfl_d_ary_tree_storage<2>)
    
    RK_MEAQR_RRTSTAR_PLANNER_CALL_RRTSTAR_FUNCTION
    
  } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
    
    RK_MEAQR_RRTSTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, boost::bfl_d_ary_tree_storage<4>)
    
    RK_MEAQR_RRTSTAR_PLANNER_CALL_RRTSTAR_FUNCTION
        
#ifdef RK_PLANNERS_ENABLE_VEBL_TREE
    
  } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
    
    RK_MEAQR_RRTSTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, boost::vebl_d_ary_tree_storage<2>)
    
    RK_MEAQR_RRTSTAR_PLANNER_CALL_RRTSTAR_FUNCTION
    
  } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
    
    RK_MEAQR_RRTSTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, boost::vebl_d_ary_tree_storage<4>)
    
    RK_MEAQR_RRTSTAR_PLANNER_CALL_RRTSTAR_FUNCTION
        
#endif
    
  };
  
#undef RK_MEAQR_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE
#undef RK_MEAQR_RRTSTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO
#undef RK_MEAQR_RRTSTAR_PLANNER_CALL_RRTSTAR_FUNCTION
  
};


};

};

#endif


/**
 * \file sbastar_path_planner.hpp
 * 
 * This library defines a class to solve path planning problems using the 
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

#ifndef REAK_SBASTAR_PATH_PLANNER_HPP
#define REAK_SBASTAR_PATH_PLANNER_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "motion_planner_base.hpp"

#include "metric_space_concept.hpp"
#include "any_sbmp_reporter.hpp"

#include "graph_alg/lazy_sbastar.hpp"
#include "graph_alg/sbastar_rrtstar.hpp"
#include "graph_alg/anytime_sbastar.hpp"

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
#include "density_plan_visitors.hpp"

namespace ReaK {
  
namespace pp {


/**
 * This class solves path planning problems using the 
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
 * \tparam FreeSpaceType The topology type on which to perform the planning, should be the C-free sub-space of a larger configuration space.
 */
template <typename FreeSpaceType>
class sbastar_planner : public sample_based_planner<FreeSpaceType> {
  public:
    typedef sample_based_planner<FreeSpaceType> base_type;
    typedef sbastar_planner<FreeSpaceType> self;
    
    typedef FreeSpaceType space_type;
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    
    BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));
    
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename topology_traits< super_space_type >::point_difference_type point_difference_type;
    
  protected:
    double m_init_dens_threshold;
    double m_init_relaxation;
    double m_SA_init_temperature;
    
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
     * Returns the initial density-value threshold used by this planner.
     * \return The initial density-value threshold used by this planner.
     */
    double get_initial_density_threshold() const { return m_init_dens_threshold; };
    /**
     * Sets the initial density-value threshold to be used by this planner.
     * \param aInitialThreshold The initial density-value threshold to be used by this planner.
     */
    void set_initial_density_threshold(double aInitialThreshold) { m_init_dens_threshold = aInitialThreshold; };
    
    /**
     * Returns the initial relaxation factor used by this planner.
     * \return The initial relaxation factor used by this planner.
     */
    double get_initial_relaxation() const { return m_init_relaxation; };
    /**
     * Sets the initial relaxation factor to be used by this planner.
     * \param aInitialRelaxation The initial relaxation factor to be used by this planner.
     */
    void set_initial_relaxation(double aInitialRelaxation) { m_init_relaxation = aInitialRelaxation; };
    
    /**
     * Returns the initial Simulated Annealing temperature use by this planner, if the 
     * added-bias is set to an exploratory bias (e.g., PLAN_WITH_VORONOI_PULL). If negative,
     * then simulated annealing is not used, and the exploratory bias (if any) is applied 
     * only when SBA* seaching stalls (isn't progressing anymore).
     * \return The initial Simulated Annealing temperature used by this planner.
     */
    double get_initial_SA_temperature() const { return m_SA_init_temperature; };
    /**
     * Sets the initial Simulated Annealing temperature use by this planner, if the 
     * added-bias is set to an exploratory bias (e.g., PLAN_WITH_VORONOI_PULL). If negative,
     * then simulated annealing is not used, and the exploratory bias (if any) is applied 
     * only when SBA* seaching stalls (isn't progressing anymore).
     * \param aInitialSATemperature The initial Simulated Annealing temperature to be used by this planner.
     */
    void set_initial_SA_temperature(double aInitialSATemperature) { m_SA_init_temperature = aInitialSATemperature; };
    
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
     * \param aSamplingRadius The sampling radius to be used by this planner when doing random walks.
     * \param aSpaceDimensionality The dimensionality of the space used by this planner.
     * \param aReporter The path-planning reporter to be used by this planner.
     */
    sbastar_planner(const shared_ptr< space_type >& aWorld = shared_ptr< space_type >(), 
                    std::size_t aMaxVertexCount = 5000, 
                    std::size_t aProgressInterval = 100,
                    std::size_t aDataStructureFlags = ADJ_LIST_MOTION_GRAPH | DVP_BF2_TREE_KNN,
                    std::size_t aPlanningMethodFlags = LAZY_COLLISION_CHECKING | NOMINAL_PLANNER_ONLY,
                    double aSteerProgressTolerance = 0.1,
                    double aConnectionTolerance = 0.1,
                    double aSamplingRadius = 1.0,
                    std::size_t aSpaceDimensionality = 1,
                    const any_sbmp_reporter_chain<space_type>& aReporter = any_sbmp_reporter_chain<space_type>()) :
                    base_type("sbastar_planner", aWorld, aMaxVertexCount, aProgressInterval,
                              aDataStructureFlags, aPlanningMethodFlags,
                              aSteerProgressTolerance, aConnectionTolerance, 
                              aSamplingRadius, aSpaceDimensionality, aReporter),
                    m_init_dens_threshold(0.8),
                    m_init_relaxation(0.0),
                    m_SA_init_temperature(-1.0) { };
    
    virtual ~sbastar_planner() { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_init_dens_threshold)
        & RK_SERIAL_SAVE_WITH_NAME(m_init_relaxation)
        & RK_SERIAL_SAVE_WITH_NAME(m_SA_init_temperature);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_init_dens_threshold)
        & RK_SERIAL_LOAD_WITH_NAME(m_init_relaxation)
        & RK_SERIAL_LOAD_WITH_NAME(m_SA_init_temperature);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC246000C,1,"sbastar_planner",base_type)
};


template <typename FreeSpaceType>
void sbastar_planner<FreeSpaceType>::solve_planning_query(planning_query<FreeSpaceType>& aQuery) {
  
  this->reset_internal_state();
  
  typedef typename subspace_traits<FreeSpaceType>::super_space_type SuperSpace;
  typedef typename topology_traits<SuperSpace>::point_type PointType;
  
  typedef recursive_dense_mg_vertex< astar_mg_vertex<FreeSpaceType> > VertexProp;
  typedef optimal_mg_edge<FreeSpaceType> EdgeProp;
  
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
  
  double space_dim = double( this->get_space_dimensionality() );
  double space_Lc = aQuery.get_heuristic_to_goal( aQuery.get_start_position() );
  
  shared_ptr<const SuperSpace> sup_space_ptr(&(this->m_space->get_super_space()),null_deleter());
  
  density_plan_visitor<FreeSpaceType, sbastar_density_calculator> vis(
    this, &aQuery, NULL, boost::any(), boost::any(), 
    this->m_init_dens_threshold);
  
  path_planning_p2p_query<FreeSpaceType>* p2p_query_ptr = reinterpret_cast< path_planning_p2p_query<FreeSpaceType>* >(aQuery.castTo(path_planning_p2p_query<FreeSpaceType>::getStaticObjectType()));
  
  // Some MACROs to reduce the size of the code below.
  
#ifdef RK_ENABLE_CXX0X_FEATURES

#define RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE \
  VertexProp vp_start; \
  vp_start.position = aQuery.get_start_position(); \
  Vertex start_node = add_vertex(std::move(vp_start), motion_graph); \
  motion_graph[start_node].constriction = 0.0; \
  motion_graph[start_node].collision_count = 0; \
  motion_graph[start_node].density = 0.0; \
  motion_graph[start_node].expansion_trials = 0; \
  motion_graph[start_node].heuristic_value = aQuery.get_heuristic_to_goal(motion_graph[start_node].position); \
  motion_graph[start_node].distance_accum = 0.0; \
  motion_graph[start_node].predecessor = start_node; \
  vis.m_start_node = boost::any( start_node ); \
  if( p2p_query_ptr ) { \
    VertexProp vp_goal; \
    vp_goal.position = p2p_query_ptr->goal_pos; \
    Vertex goal_node  = add_vertex(std::move(vp_goal),  motion_graph); \
    motion_graph[goal_node].constriction = 0.0; \
    motion_graph[goal_node].collision_count = 0; \
    motion_graph[goal_node].density = 0.0; \
    motion_graph[goal_node].expansion_trials = 0; \
    motion_graph[goal_node].heuristic_value = 0.0; \
    motion_graph[goal_node].distance_accum = std::numeric_limits<double>::infinity(); \
    motion_graph[goal_node].predecessor = goal_node; \
    vis.m_goal_node = boost::any( goal_node ); \
  };
    
#else
    
#define RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE \
  VertexProp vp_start; \
  vp_start.position = aQuery.get_start_position(); \
  Vertex start_node = add_vertex(vp_start, motion_graph); \
  motion_graph[start_node].constriction = 0.0; \
  motion_graph[start_node].collision_count = 0; \
  motion_graph[start_node].density = 0.0; \
  motion_graph[start_node].expansion_trials = 0; \
  motion_graph[start_node].heuristic_value = aQuery.get_heuristic_to_goal(motion_graph[start_node].position); \
  motion_graph[start_node].distance_accum = 0.0; \
  motion_graph[start_node].predecessor = start_node; \
  vis.m_start_node = boost::any( start_node ); \
  if( p2p_query_ptr ) { \
    VertexProp vp_goal; \
    vp_goal.position = p2p_query_ptr->goal_pos; \
    Vertex goal_node  = add_vertex(vp_goal,  motion_graph); \
    motion_graph[goal_node].constriction = 0.0; \
    motion_graph[goal_node].collision_count = 0; \
    motion_graph[goal_node].density = 0.0; \
    motion_graph[goal_node].expansion_trials = 0; \
    motion_graph[goal_node].heuristic_value = 0.0; \
    motion_graph[goal_node].distance_accum = std::numeric_limits<double>::infinity(); \
    motion_graph[goal_node].predecessor = goal_node; \
    vis.m_goal_node = boost::any( goal_node ); \
  };
  
#endif
  
  
#define RK_SBASTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(ARITY, TREE_STORAGE) \
  typedef typename boost::property_map< MotionGraphType, PointType VertexProp::* >::type GraphPositionMap; \
  typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, ARITY, random_vp_chooser, TREE_STORAGE > SpacePartType; \
  SpacePartType space_part(motion_graph, sup_space_ptr, get(&VertexProp::position, motion_graph)); \
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
    boost::vecS, boost::undirectedS, boost::listS > ALTGraph; \
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
  
  
  
  
    
    
#define RK_SBASTAR_PLANNER_CALL_SBASTAR_FUNCTION \
    ReaK::graph::generate_sbastar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, boost::any_cast<Vertex>( vis.m_start_node ), *sup_space_ptr, vis, \
        heuristic_map, pos_map, weight_map, dens_map, cons_map, dist_map, pred_map, \
        get(&VertexProp::key_value, motion_graph), nc_selector)\
      );
  
  
#define RK_SBASTAR_PLANNER_CALL_LAZY_SBASTAR_FUNCTION \
    ReaK::graph::generate_lazy_sbastar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, boost::any_cast<Vertex>( vis.m_start_node ), *sup_space_ptr, vis, \
        heuristic_map, pos_map, weight_map, dens_map, cons_map, dist_map, pred_map, \
        get(&VertexProp::key_value, motion_graph), nc_selector) );
  
#define RK_SBASTAR_PLANNER_CALL_LAZY_BNB_SBASTAR_FUNCTION \
    ReaK::graph::generate_lazy_bnb_sbastar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, boost::any_cast<Vertex>( vis.m_start_node ), *sup_space_ptr, vis, \
        heuristic_map, pos_map, weight_map, dens_map, cons_map, dist_map, pred_map, \
        get(&VertexProp::key_value, motion_graph), nc_selector), \
      boost::any_cast<Vertex>( vis.m_goal_node ) );
  
#define RK_SBASTAR_PLANNER_CALL_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_sbarrtstar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, boost::any_cast<Vertex>( vis.m_start_node ), *sup_space_ptr, vis, \
        heuristic_map, pos_map, weight_map, dens_map, cons_map, dist_map, pred_map, \
        get(&VertexProp::key_value, motion_graph), nc_selector), \
      boost::any_cast<Vertex>( vis.m_goal_node ), \
      get(random_sampler, *sup_space_ptr), \
      this->m_SA_init_temperature);
  
  
#define RK_SBASTAR_PLANNER_CALL_LAZY_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_lazy_sbarrtstar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, boost::any_cast<Vertex>( vis.m_start_node ), *sup_space_ptr, vis, \
        heuristic_map, pos_map, weight_map, dens_map, cons_map, dist_map, pred_map, \
        get(&VertexProp::key_value, motion_graph), nc_selector), \
      boost::any_cast<Vertex>( vis.m_goal_node ), \
      get(random_sampler, *sup_space_ptr), \
      this->m_SA_init_temperature);
  
  
#define RK_SBASTAR_PLANNER_CALL_LAZY_BNB_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_lazy_bnb_sbarrtstar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, boost::any_cast<Vertex>( vis.m_start_node ), *sup_space_ptr, vis, \
        heuristic_map, pos_map, weight_map, dens_map, cons_map, dist_map, pred_map, \
        get(&VertexProp::key_value, motion_graph), nc_selector), \
      boost::any_cast<Vertex>( vis.m_goal_node ), \
      get(random_sampler, *sup_space_ptr), \
      this->m_SA_init_temperature);
   
  
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_SBASTAR_FUNCTION \
    ReaK::graph::generate_anytime_sbastar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, boost::any_cast<Vertex>( vis.m_start_node ), *sup_space_ptr, vis, \
        heuristic_map, pos_map, weight_map, dens_map, cons_map, dist_map, pred_map, \
        get(&VertexProp::key_value, motion_graph), nc_selector), \
      this->m_init_relaxation);
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBASTAR_FUNCTION \
    ReaK::graph::generate_anytime_lazy_sbastar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, boost::any_cast<Vertex>( vis.m_start_node ), *sup_space_ptr, vis, \
        heuristic_map, pos_map, weight_map, dens_map, cons_map, dist_map, pred_map, \
        get(&VertexProp::key_value, motion_graph), nc_selector), \
      this->m_init_relaxation);
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_BNB_SBASTAR_FUNCTION \
    ReaK::graph::generate_anytime_lazy_bnb_sbastar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, boost::any_cast<Vertex>( vis.m_start_node ), *sup_space_ptr, vis, \
        heuristic_map, pos_map, weight_map, dens_map, cons_map, dist_map, pred_map, \
        get(&VertexProp::key_value, motion_graph), nc_selector), \
      boost::any_cast<Vertex>( vis.m_goal_node ), \
      this->m_init_relaxation);
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_anytime_sbarrtstar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, boost::any_cast<Vertex>( vis.m_start_node ), *sup_space_ptr, vis, \
        heuristic_map, pos_map, weight_map, dens_map, cons_map, dist_map, pred_map, \
        get(&VertexProp::key_value, motion_graph), nc_selector), \
      boost::any_cast<Vertex>( vis.m_goal_node ), \
      get(random_sampler, *sup_space_ptr), \
      this->m_init_relaxation, this->m_SA_init_temperature);
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_anytime_lazy_sbarrtstar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, boost::any_cast<Vertex>( vis.m_start_node ), *sup_space_ptr, vis, \
        heuristic_map, pos_map, weight_map, dens_map, cons_map, dist_map, pred_map, \
        get(&VertexProp::key_value, motion_graph), nc_selector), \
      boost::any_cast<Vertex>( vis.m_goal_node ), \
      get(random_sampler, *sup_space_ptr), \
      this->m_init_relaxation, this->m_SA_init_temperature);
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_BNB_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_anytime_lazy_bnb_sbarrtstar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, boost::any_cast<Vertex>( vis.m_start_node ), *sup_space_ptr, vis, \
        heuristic_map, pos_map, weight_map, dens_map, cons_map, dist_map, pred_map, \
        get(&VertexProp::key_value, motion_graph), nc_selector), \
      boost::any_cast<Vertex>( vis.m_goal_node ), \
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
    
    typedef boost::pooled_adjacency_list< 
      boost::undirectedS, VertexProp, EdgeProp, 
      boost::no_property, boost::listS> MotionGraphType;
    
    typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex;
    
    MotionGraphType motion_graph;
    
    RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
    
    if((this->m_data_structure_flags & KNN_METHOD_MASK) == LINEAR_SEARCH_KNN) {
      
      typedef linear_neighbor_search<> NNFinderType;
      NNFinderType nn_finder;
      
      ReaK::graph::star_neighborhood< NNFinderType > nc_selector(nn_finder, space_dim, 3.0 * space_Lc);
//       ReaK::graph::fixed_neighborhood< NNFinderType > nc_selector(nn_finder, 10, this->get_sampling_radius());
      
      any_knn_synchro NN_synchro;
      vis.m_nn_synchro = &NN_synchro;
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, graph::d_ary_bf_tree_storage<2>)
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, graph::d_ary_bf_tree_storage<4>)
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, graph::d_ary_cob_tree_storage<2>)
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, graph::d_ary_cob_tree_storage<4>)
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    };
    
  } else if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == DVP_ADJ_LIST_MOTION_GRAPH) {
    
    if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO(2, graph::d_ary_bf_tree_storage<2>)
      
      RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO(4, graph::d_ary_bf_tree_storage<4>)
      
      RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO(2, graph::d_ary_cob_tree_storage<2>)
      
      RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
      
      RK_SBASTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO(4, graph::d_ary_cob_tree_storage<4>)
      
      RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    };
    
  };
  
#undef RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
#undef RK_SBASTAR_PLANNER_SETUP_DVP_TREE_SYNCHRO
#undef RK_SBASTAR_PLANNER_SETUP_ALT_TREE_SYNCHRO
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


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
#include "sbmp_reporter_concept.hpp"

#include "metric_space_concept.hpp"
#include "seq_path_wrapper.hpp"
#include "interpolation/point_to_point_path.hpp"
#include "basic_sbmp_reporters.hpp"

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
#include "rrt_path_planner.hpp"
#include "lin_alg/arithmetic_tuple.hpp"
#include "any_motion_graphs.hpp"

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
template <typename FreeSpaceType, 
          typename SBPPReporter = no_sbmp_report>
class rrtstar_path_planner : public sample_based_planner< path_planner_base<FreeSpaceType> > {
  public:
    typedef sample_based_planner< path_planner_base<FreeSpaceType> > base_type;
    typedef rrtstar_path_planner<FreeSpaceType, SBPPReporter> self;
    
    typedef FreeSpaceType space_type;
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    
    BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));
    
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename topology_traits< super_space_type >::point_difference_type point_difference_type;
    
  protected:
    SBPPReporter m_reporter;
    point_type m_start_pos;
    point_type m_goal_pos;
    std::size_t max_num_results;
    
    std::map<double, shared_ptr< seq_path_base< super_space_type > > > m_solutions;
    
  public:
    
    /**
     * Returns the best solution distance obtained after solve_path() has been called.
     * \return The best solution distance obtained after solve_path() has been called.
     */
    double get_best_solution_distance() const {
      if(m_solutions.size() == 0)
        return std::numeric_limits<double>::infinity();
      else
        return m_solutions.begin()->first;
    };
    
    /**
     * This function constructs a solution path (if one is found) and invokes the path-planning 
     * reporter to report on that solution path.
     * \note This function is for internal use by the path-planning algorithm (a visitor callback).
     * \param start_node The start node in the motion-graph.
     * \param goal_node The goal node in the motion-graph.
     * \param g The current motion-graph.
     */
    template <typename Vertex, typename Graph>
    void create_solution_path(Vertex start_node, Vertex goal_node, Graph& g) {
      
      double goal_distance = g[goal_node].distance_accum;
      
      if(goal_distance < std::numeric_limits<double>::infinity()) {
        //Draw the edges of the current best solution:
        
        shared_ptr< super_space_type > sup_space_ptr(&(this->m_space->get_super_space()),null_deleter());
        shared_ptr< seq_path_wrapper< point_to_point_path<super_space_type> > > new_sol(new seq_path_wrapper< point_to_point_path<super_space_type> >("rrtstar_solution", point_to_point_path<super_space_type>(sup_space_ptr,get(distance_metric, this->m_space->get_super_space()))));
        point_to_point_path<super_space_type>& waypoints = new_sol->get_underlying_path();
        std::set<Vertex> path;
        
        Vertex v = goal_node;
        point_type p_v = g[v].position;
        Vertex u = g[v].predecessor;
        point_type p_u = g[u].position;
        
        waypoints.push_front(p_v);
        waypoints.push_front(p_u);
        path.insert(v);
      
        while((u != start_node) && (path.insert(u).second)) {
          v = u; p_v = p_u;
          u = g[v].predecessor; 
          p_u = g[u].position; 
          waypoints.push_front(p_u);
        };
        
        if(u == start_node) {
          m_solutions[goal_distance] = new_sol;
          m_reporter.draw_solution(*(this->m_space), m_solutions[goal_distance]);
        };
      };
    };
    
    /**
     * This function constructs a solution path (if one is found) and invokes the path-planning 
     * reporter to report on that solution path. This is the bi-directional version which is 
     * called when a joining vertex is found between the two motion-graphs.
     * \note This function is for internal use by the path-planning algorithm (a visitor callback).
     * \param u1 The latest added node in the first motion-graph.
     * \param u2 The latest added node in the second motion-graph.
     * \param g1 The first motion-graph.
     * \param g2 The second motion-graph.
     */
    template <typename Vertex, typename Graph>
    void joining_vertex_found(Vertex u1, Vertex u2, Graph& g1, Graph& g2) {
      double total_dist = g1[u1].distance_accum + g2[u2].distance_accum
        + get(distance_metric, this->m_space->get_super_space())(g1[u1].position, g2[u2].position, this->m_space->get_super_space());
      
      if((m_solutions.size()) && (total_dist >= m_solutions.begin()->first))
        return;
      
      shared_ptr< super_space_type > sup_space_ptr(&(this->m_space->get_super_space()),null_deleter());
      shared_ptr< seq_path_wrapper< point_to_point_path<super_space_type> > > new_sol(new seq_path_wrapper< point_to_point_path<super_space_type> >("birrt_solution", point_to_point_path<super_space_type>(sup_space_ptr,get(distance_metric, this->m_space->get_super_space()))));
      point_to_point_path<super_space_type>& waypoints = new_sol->get_underlying_path();
      
      waypoints.push_front(g1[u1].position);
      while( ( g1[u1].predecessor != Graph::null_vertex() ) && ( g1[u1].predecessor != u1 ) ) {
        u1 = g1[u1].predecessor;
        waypoints.push_front(g1[u1].position);
      };
      
      waypoints.push_back(g2[u2].position);
      while( ( g2[u2].predecessor != Graph::null_vertex() ) && ( g2[u2].predecessor != u2 ) ) {
        u2 = g2[u2].predecessor;
        waypoints.push_back(g2[u2].position);
      };
      
      m_solutions[total_dist] = new_sol;
      m_reporter.draw_solution(*(this->m_space), m_solutions[total_dist]);
    };
    
    /**
     * Returns true if the solver should keep on going trying to solve the path-planning problem.
     * \note This function is for internal use by the path-planning algorithm (a visitor callback).
     * \return True if the solver should keep on going trying to solve the path-planning problem.
     */
    bool keep_going() const {
      return (max_num_results > m_solutions.size()) && !(this->has_reached_max_iterations());
    };
    
    /**
     * This function invokes the path-planning reporter to report on the progress of the path-planning
     * solver.
     * \note This function is for internal use by the path-planning algorithm (a visitor callback).
     * \param g The current motion-graph.
     */
    template <typename Graph>
    void register_progress(Graph& g) {
      this->report_progress(g, m_reporter);
    };
    
    virtual void reset_internal_state() {
      base_type::reset_internal_state();
      m_solutions.clear();
    };
    
    /**
     * This function computes a valid path in the C-free. If it cannot 
     * achieve a valid path, an exception will be thrown. This algorithmic
     * path solver class is such that any settings that ought to be set for the 
     * path planning algorithm should be set before calling this function, otherwise
     * the function is likely to fail.
     * \return The path object that can be used to map out the path.
     */
    virtual shared_ptr< seq_path_base< super_space_type > > solve_path();
    
    /**
     * Returns a const-reference to the path-planning reporter used by this planner.
     * \return A const-reference to the path-planning reporter used by this planner.
     */
    const SBPPReporter& get_reporter() const { return m_reporter; };
    /**
     * Sets the path-planning reporter to be used by this planner.
     * \param aNewReporter The path-planning reporter to be used by this planner.
     */
    void set_reporter(const SBPPReporter& aNewReporter) { m_reporter = aNewReporter; };
    
    /**
     * Returns a const-reference to the start position used by this planner.
     * \return A const-reference to the start position used by this planner.
     */
    const point_type& get_start_pos() const { return m_start_pos; };
    /**
     * Sets the start position to be used by this planner.
     * \param aNewReporter The start position to be used by this planner.
     */
    void set_start_pos(const point_type& aStartPos) { m_start_pos = aStartPos; };
    
    /**
     * Returns a const-reference to the goal position used by this planner.
     * \return A const-reference to the goal position used by this planner.
     */
    const point_type& get_goal_pos() const { return m_goal_pos; };
    /**
     * Sets the goal position to be used by this planner.
     * \param aNewReporter The goal position to be used by this planner.
     */
    void set_goal_pos(const point_type& aGoalPos) { m_goal_pos = aGoalPos; };
    
    /**
     * Returns the maximum number of solutions that this planner should register.
     * \note Most probabilistic path-planners produce an initial solution that isn't perfect, and so, 
     *       more solutions should be sought to eventually have a more optimal one (shorter distance).
     * \return The maximum number of solutions that this planner should register.
     */
    std::size_t get_max_result_count() const { return max_num_results; };
    /**
     * Sets the maximum number of solutions that this planner should register.
     * \note Most probabilistic path-planners produce an initial solution that isn't perfect, and so, 
     *       more solutions should be sought to eventually have a more optimal one (shorter distance).
     * \param aMaxResultCount The maximum number of solutions that this planner should register.
     */
    void set_max_result_count(std::size_t aMaxResultCount) { max_num_results = aMaxResultCount; };
    
    
    /**
     * Parametrized constructor.
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     * \param aStartPos The position value of the starting location.
     * \param aGoalPos The position value of the goal location.
     * \param aMaxVertexCount The maximum number of samples to generate during the motion planning.
     * \param aProgressInterval The number of new samples between each "progress report".
     * \param aDataStructureFlags An integer flags representing the kind of motion graph data-structure to use in the 
     *                            planning algorithm. Can be ADJ_LIST_MOTION_GRAPH or DVP_ADJ_LIST_MOTION_GRAPH.
     *                            Any combination of those two and of KNN method flags to use for nearest
     *                            neighbor queries in the graph. KNN method flags can be LINEAR_SEARCH_KNN, 
     *                            DVP_BF2_TREE_KNN, DVP_BF4_TREE_KNN, DVP_COB2_TREE_KNN, or DVP_COB4_TREE_KNN.
     *                            See path_planner_options.hpp documentation.
     * \param aPlanningMethodFlags The integer flags that identify various options to use with this planner.
     *                             The options available include only USE_BRANCH_AND_BOUND_PRUNING_FLAG. 
     *                             See path_planner_options.hpp documentation.
     * \param aReporter The SBPP reporter object to use to report results and progress.
     * \param aMaxResultCount The maximum number of successful start-goal connections to make before 
     *                        stopping the path planner (the higher the number the more likely that a 
     *                        good path will be found, however, running time can become much longer).
     */
    rrtstar_path_planner(const shared_ptr< space_type >& aWorld = shared_ptr< space_type >(), 
                         const point_type& aStartPos = point_type(),
                         const point_type& aGoalPos = point_type(),
                         std::size_t aMaxVertexCount = 5000, 
                         std::size_t aProgressInterval = 100,
                         std::size_t aDataStructureFlags = ADJ_LIST_MOTION_GRAPH | DVP_BF2_TREE_KNN,
                         std::size_t aPlanningMethodFlags = UNIDIRECTIONAL_PLANNING,
                         SBPPReporter aReporter = SBPPReporter(),
                         std::size_t aMaxResultCount = 50) :
                         base_type("rrtstar_planner", aWorld, aMaxVertexCount, aProgressInterval, aDataStructureFlags, aPlanningMethodFlags),
                         m_reporter(aReporter),
                         m_start_pos(aStartPos),
                         m_goal_pos(aGoalPos),
                         max_num_results(aMaxResultCount),
                         m_solutions() { };
    
    virtual ~rrtstar_path_planner() { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_reporter)
        & RK_SERIAL_SAVE_WITH_NAME(m_start_pos)
        & RK_SERIAL_SAVE_WITH_NAME(m_goal_pos)
        & RK_SERIAL_SAVE_WITH_NAME(max_num_results);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_reporter)
        & RK_SERIAL_LOAD_WITH_NAME(m_start_pos)
        & RK_SERIAL_LOAD_WITH_NAME(m_goal_pos)
        & RK_SERIAL_LOAD_WITH_NAME(max_num_results);
      m_solutions.clear();
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2460009,1,"rrtstar_path_planner",base_type)
};




/**
 * This class template is used by the RRT* path-planner as the visitor object needed to 
 * collaborate with the RRT* algorithms to generate the motion-graph and path-planning solutions.
 * This class template models the RRTStarVisitorConcept.
 * As with most planning algorithms in ReaK, the algorithm is really made up of a high-level 
 * algorithmic logic in the form of function templates, and a number of customization points 
 * collected as member functions of an algorithm visitor class that implement the problem-specific 
 * behaviors (random-walks / local-planning, progress reporting, completion criteria, etc.). 
 */
template <typename FreeSpaceType, typename MotionGraph, typename NNFinderSynchro, typename SBPPReporter = no_sbmp_report>
struct rrtstar_planner_visitor {
  typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
  
  shared_ptr< FreeSpaceType > m_space;
  rrtstar_path_planner<FreeSpaceType,SBPPReporter>* m_planner;
  NNFinderSynchro m_nn_synchro;
  Vertex m_start_node;
  Vertex m_goal_node;
  
  rrtstar_planner_visitor(const shared_ptr< FreeSpaceType >& aSpace, 
                          rrtstar_path_planner<FreeSpaceType,SBPPReporter>* aPlanner,
                          NNFinderSynchro aNNSynchro,
                          Vertex aStartNode, Vertex aGoalNode) : 
                          m_space(aSpace), m_planner(aPlanner), m_nn_synchro(aNNSynchro),
                          m_start_node(aStartNode), m_goal_node(aGoalNode) { };
  
  typedef typename topology_traits<FreeSpaceType>::point_type PointType;
  typedef optimal_mg_edge<FreeSpaceType> EdgeProp;
  
  template <typename Vertex, typename Graph>
  void vertex_added(Vertex u, Graph& g) const {
    m_nn_synchro.added_vertex(u,g);
    
    // Call progress reporter...
    m_planner->register_progress(g);
    
    if((in_degree(m_goal_node,g)) && (g[m_goal_node].distance_accum < m_planner->get_best_solution_distance()))
      m_planner->create_solution_path(m_start_node, m_goal_node, g);
  };
  
  template <typename Vertex, typename Graph>
  void vertex_to_be_removed(Vertex u, Graph& g) const {
    m_nn_synchro.removed_vertex(u,g);
  };
  
  template <typename EdgeType, typename Graph>
  void edge_added(EdgeType e, Graph& g) const { };
  
  bool is_position_free(PointType p) const {
    return m_space->is_free(p);
  };
  
  bool keep_going() const {
    return m_planner->keep_going();
  };
  
  template <typename Vertex, typename Graph>
  void joining_vertex_found(Vertex u1, Vertex u2, Graph& g1, Graph& g2) const {
    m_planner->joining_vertex_found(u1, u2, g1, g2);
  };
  
  template <typename Vertex, typename Graph>
  boost::tuple<PointType, bool, EdgeProp> steer_towards_position(const PointType& p, Vertex u, Graph& g) const {
    typedef boost::tuple<PointType, bool, EdgeProp> ResultType;
    PointType result_p = m_space->move_position_toward(g[u].position, 1.0, p);
    double best_case_dist = get(distance_metric, m_space->get_super_space())(g[u].position, p, m_space->get_super_space());
    double actual_dist = get(distance_metric, m_space->get_super_space())(g[u].position, result_p, m_space->get_super_space());
    if(actual_dist > 0.1 * best_case_dist)
      return ResultType(result_p, true, EdgeProp(actual_dist));
    else
      return ResultType(result_p, false, EdgeProp(actual_dist));
  };
  
  template <typename Vertex, typename Graph>
  std::pair<bool, EdgeProp> can_be_connected(Vertex u, Vertex v, const Graph& g) const {
    double dist = get(distance_metric, *m_space)(g[u].position, g[v].position, *m_space);
    return std::pair<bool, EdgeProp>((dist < std::numeric_limits<double>::infinity()), EdgeProp(dist));
  };
  
};






template <typename FreeSpaceType, 
          typename SBPPReporter>
shared_ptr< seq_path_base< typename rrtstar_path_planner<FreeSpaceType,SBPPReporter>::super_space_type > > 
  rrtstar_path_planner<FreeSpaceType,SBPPReporter>::solve_path() {
  using ReaK::to_vect;
  
  this->m_solutions.clear();
  
  typedef typename subspace_traits<FreeSpaceType>::super_space_type SuperSpace;
  typedef typename topology_traits<SuperSpace>::point_type PointType;
  typedef boost::data_member_property_map<PointType, optimal_mg_vertex<FreeSpaceType> > PositionMap;
  PositionMap pos_map = PositionMap(&mg_vertex_data<FreeSpaceType>::position);
  typedef boost::data_member_property_map<double, optimal_mg_vertex<FreeSpaceType> > CostMap;
  CostMap cost_map = CostMap(&optimal_mg_vertex<FreeSpaceType>::distance_accum);
  typedef boost::data_member_property_map<std::size_t, optimal_mg_vertex<FreeSpaceType> > PredMap;
  PredMap pred_map = PredMap(&optimal_mg_vertex<FreeSpaceType>::predecessor);
  typedef boost::data_member_property_map<double, optimal_mg_edge<FreeSpaceType> > WeightMap;
  WeightMap weight_map = WeightMap(&optimal_mg_edge<FreeSpaceType>::weight);
  
  double space_dim = double((to_vect<double>(this->m_space->get_super_space().difference(this->m_goal_pos,this->m_start_pos))).size()); 
  double space_Lc = get(distance_metric,this->m_space->get_super_space())(this->m_start_pos, this->m_goal_pos, this->m_space->get_super_space());
  
  
#define RK_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE \
      optimal_mg_vertex<FreeSpaceType> vs_p; \
      vs_p.position = this->m_start_pos; \
      Vertex vs = add_vertex(vs_p, motion_graph); \
      motion_graph[vs].distance_accum = 0.0; \
      motion_graph[vs].predecessor = vs; \
      optimal_mg_vertex<FreeSpaceType> vg_p; \
      vg_p.position = this->m_goal_pos; \
      Vertex vg = add_vertex(vg_p, motion_graph); \
      motion_graph[vg].distance_accum = std::numeric_limits<double>::infinity(); \
      motion_graph[vg].predecessor = vg;
  
  
#define RK_RRTSTAR_PLANNER_CALL_RRTSTAR_FUNCTION \
        ReaK::graph::generate_rrt_star( \
          motion_graph, \
          this->m_space->get_super_space(), \
          vis, \
          pos_map, \
          cost_map, \
          pred_map, \
          weight_map, \
          get(random_sampler, this->m_space->get_super_space()), \
          ReaK::graph::star_neighborhood< NNFinderType >( \
            nn_finder, \
            space_dim, 3.0 * space_Lc), \
          this->m_max_vertex_count);

#define RK_RRTSTAR_PLANNER_CALL_RRTSTAR_BNB_FUNCTION \
        ReaK::graph::generate_bnb_rrt_star( \
          motion_graph, \
          vs, vg, \
          this->m_space->get_super_space(), \
          vis, \
          pos_map, \
          cost_map, \
          pred_map, \
          weight_map, \
          get(random_sampler, this->m_space->get_super_space()), \
          ReaK::graph::star_neighborhood< NNFinderType >( \
            nn_finder, \
            space_dim, 3.0 * space_Lc), \
          this->m_max_vertex_count);
  
  
#define RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION \
      if(this->m_planning_method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) { \
        RK_RRTSTAR_PLANNER_CALL_RRTSTAR_BNB_FUNCTION \
      } else { /* assume nominal method only. */ \
        RK_RRTSTAR_PLANNER_CALL_RRTSTAR_FUNCTION \
      };
  
  
  if((this->m_planning_method_flags & PLANNING_DIRECTIONALITY_MASK) == UNIDIRECTIONAL_PLANNING) {
    
    if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == ADJ_LIST_MOTION_GRAPH) {
      
      typedef boost::pooled_adjacency_list< 
        boost::undirectedS,
        optimal_mg_vertex<FreeSpaceType>,
        optimal_mg_edge<FreeSpaceType>,
        boost::no_property,
        boost::listS> MotionGraphType;
      
      typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex;
      typedef typename MotionGraphType::vertex_property_type VertexProp;
      typedef boost::composite_property_map< 
        PositionMap, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t > > GraphPositionMap;
      
      MotionGraphType motion_graph;
      GraphPositionMap g_pos_map = GraphPositionMap(pos_map, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t >(&motion_graph));
      
      RK_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE
      
      if((this->m_data_structure_flags & KNN_METHOD_MASK) == LINEAR_SEARCH_KNN) {
        
        typedef linear_neighbor_search<> NNFinderType;
        NNFinderType nn_finder;
        
        rrtstar_planner_visitor<FreeSpaceType, MotionGraphType, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), vs, vg);
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 2, 
                         random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2> > SpacePartType;
        SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
        
        typedef multi_dvp_tree_search<MotionGraphType, SpacePartType> NNFinderType;
        NNFinderType nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrtstar_planner_visitor<FreeSpaceType, MotionGraphType, NNFinderType, SBPPReporter> vis(this->m_space, this, nn_finder, vs, vg);
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 4, 
                         random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4> > SpacePartType;
        SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
        
        typedef multi_dvp_tree_search<MotionGraphType, SpacePartType> NNFinderType;
        NNFinderType nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrtstar_planner_visitor<FreeSpaceType, MotionGraphType, NNFinderType, SBPPReporter> vis(this->m_space, this, nn_finder, vs, vg);
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 2, 
                         random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2> > SpacePartType;
        SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
        
        typedef multi_dvp_tree_search<MotionGraphType, SpacePartType> NNFinderType;
        NNFinderType nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrtstar_planner_visitor<FreeSpaceType, MotionGraphType, NNFinderType, SBPPReporter> vis(this->m_space, this, nn_finder, vs, vg);
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 4, 
                         random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4> > SpacePartType;
        SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
        
        typedef multi_dvp_tree_search<MotionGraphType, SpacePartType> NNFinderType;
        NNFinderType nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrtstar_planner_visitor<FreeSpaceType, MotionGraphType, NNFinderType, SBPPReporter> vis(this->m_space, this, nn_finder, vs, vg);
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      };
      
    } else if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == DVP_ADJ_LIST_MOTION_GRAPH) {
      
      if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
        
        typedef dvp_adjacency_list<
          optimal_mg_vertex<FreeSpaceType>,
          optimal_mg_edge<FreeSpaceType>,
          SuperSpace,
          PositionMap,
          2, random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2>,
          boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
        
        ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
        
        MotionGraph motion_graph = space_part.get_adjacency_list();
        
        RK_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE
        
        typedef multi_dvp_tree_search<MotionGraph, ALTGraph> NNFinderType;
        NNFinderType nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrtstar_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), vs, vg);
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
        
        typedef dvp_adjacency_list<
          optimal_mg_vertex<FreeSpaceType>,
          optimal_mg_edge<FreeSpaceType>,
          SuperSpace,
          PositionMap,
          4, random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4>,
          boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
        
        ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
        
        MotionGraph motion_graph = space_part.get_adjacency_list();
        
        RK_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE
        
        typedef multi_dvp_tree_search<MotionGraph, ALTGraph> NNFinderType;
        NNFinderType nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrtstar_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), vs, vg);
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
        
        typedef dvp_adjacency_list<
          optimal_mg_vertex<FreeSpaceType>,
          optimal_mg_edge<FreeSpaceType>,
          SuperSpace,
          PositionMap,
          2, random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2>,
          boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
        
        ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
        
        MotionGraph motion_graph = space_part.get_adjacency_list();
        
        RK_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE
        
        typedef multi_dvp_tree_search<MotionGraph, ALTGraph> NNFinderType;
        NNFinderType nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrtstar_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), vs, vg);
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
        
        typedef dvp_adjacency_list<
          optimal_mg_vertex<FreeSpaceType>,
          optimal_mg_edge<FreeSpaceType>,
          SuperSpace,
          PositionMap,
          4, random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4>,
          boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
        
        ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
        
        MotionGraph motion_graph = space_part.get_adjacency_list();
        
        RK_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE
        
        typedef multi_dvp_tree_search<MotionGraph, ALTGraph> NNFinderType;
        NNFinderType nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrtstar_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), vs, vg);
        
        RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
        
      };
      
    };
    
  } else {
#if 0    
    if((m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == ADJ_LIST_MOTION_GRAPH) {
      
      typedef boost::adjacency_list< boost::vecS, boost::listS, boost::bidirectionalS,
                             optimal_mg_vertex<FreeSpaceType>,
                             optimal_mg_edge<FreeSpaceType>,
                             boost::vecS> MotionGraphType;
      typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex;
      
      typedef boost::composite_property_map< 
        PositionMap, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t > > GraphPositionMap;
      
      MotionGraphType motion_graph1;
      GraphPositionMap g1_pos_map = GraphPositionMap(pos_map, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t >(&motion_graph1));
      Vertex v1 = add_vertex(motion_graph1);
      motion_graph1[v1].position = this->m_start_pos;
      motion_graph1[v1].distance_accum = 0.0;
      
      MotionGraphType motion_graph2;
      GraphPositionMap g2_pos_map = GraphPositionMap(pos_map, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t >(&motion_graph2));
      Vertex v2 = add_vertex(motion_graph2);
      motion_graph2[v2].position = this->m_goal_pos;
      motion_graph2[v2].distance_accum = 0.0;
      
      if(m_knn_flag == LINEAR_SEARCH_KNN) {
        rrtstar_planner_visitor<FreeSpaceType, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro());
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, pos_map, get(random_sampler, this->m_space->get_super_space()), 
          linear_neighbor_search<>(), this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_BF2_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 2, 
                         random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2> > SpacePartType;
        SpacePartType space_part1(motion_graph1, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  g1_pos_map);
        SpacePartType space_part2(motion_graph2, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  g2_pos_map);
        
        multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrtstar_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, pos_map, get(random_sampler, this->m_space->get_super_space()), 
          nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_BF4_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 4, 
                         random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4> > SpacePartType;
        SpacePartType space_part1(motion_graph1, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  g1_pos_map);
        SpacePartType space_part2(motion_graph2, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  g2_pos_map);
        
        multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrtstar_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, pos_map, get(random_sampler, this->m_space->get_super_space()), 
          nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_COB2_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 2, 
                         random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2> > SpacePartType;
        SpacePartType space_part1(motion_graph1, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  g1_pos_map);
        SpacePartType space_part2(motion_graph2, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  g2_pos_map);
        
        multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrtstar_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, pos_map, get(random_sampler, this->m_space->get_super_space()), 
          nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_COB4_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 4, 
                         random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4> > SpacePartType;
        SpacePartType space_part1(motion_graph1, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  g1_pos_map);
        SpacePartType space_part2(motion_graph2, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  g2_pos_map);
        
        multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrtstar_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, pos_map, get(random_sampler, this->m_space->get_super_space()), 
          nn_finder, this->m_max_vertex_count);
        
      };
      
    } else if((m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == DVP_ADJ_LIST_MOTION_GRAPH) {
      
      if(m_knn_flag == DVP_ALT_BF2_KNN) {
        
        typedef dvp_adjacency_list<
          optimal_mg_vertex<FreeSpaceType>,
          optimal_mg_edge<FreeSpaceType>,
          SuperSpace,
          PositionMap,
          2, random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2>,
          boost::vecS, boost::bidirectionalS, boost::listS > ALTGraph;
        
        ALTGraph space_part1(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                             pos_map);
        
        ALTGraph space_part2(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                             pos_map);
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        
        MotionGraph motion_graph1 = space_part1.get_adjacency_list();
        optimal_mg_vertex<FreeSpaceType> v1_p;
        v1_p.position = this->m_start_pos;
        v1_p.distance_accum = 0.0;
        create_root(v1_p, motion_graph1);
        
        MotionGraph motion_graph2 = space_part2.get_adjacency_list();
        optimal_mg_vertex<FreeSpaceType> v2_p;
        v2_p.position = this->m_goal_pos;
        v2_p.distance_accum = 0.0;
        create_root(v2_p, motion_graph2);
        
        multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrtstar_planner_visitor<FreeSpaceType, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro());
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, pos_map, get(random_sampler, this->m_space->get_super_space()), 
          nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_ALT_BF4_KNN) {
        
        typedef dvp_adjacency_list<
          optimal_mg_vertex<FreeSpaceType>,
          optimal_mg_edge<FreeSpaceType>,
          SuperSpace,
          PositionMap,
          4, random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4>,
          boost::vecS, boost::bidirectionalS, boost::listS > ALTGraph;
        
        ALTGraph space_part1(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
        
        ALTGraph space_part2(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        
        MotionGraph motion_graph1 = space_part1.get_adjacency_list();
        optimal_mg_vertex<FreeSpaceType> v1_p;
        v1_p.position = this->m_start_pos;
        v1_p.distance_accum = 0.0;
        create_root(v1_p, motion_graph1);
        
        MotionGraph motion_graph2 = space_part2.get_adjacency_list();
        optimal_mg_vertex<FreeSpaceType> v2_p;
        v2_p.position = this->m_goal_pos;
        v2_p.distance_accum = 0.0;
        create_root(v2_p, motion_graph2);
        
        multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrtstar_planner_visitor<FreeSpaceType, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro());
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, pos_map, get(random_sampler, this->m_space->get_super_space()), 
          nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_ALT_COB2_KNN) {
        
        typedef dvp_adjacency_list<
          optimal_mg_vertex<FreeSpaceType>,
          optimal_mg_edge<FreeSpaceType>,
          SuperSpace,
          PositionMap,
          2, random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2>,
          boost::vecS, boost::bidirectionalS, boost::listS > ALTGraph;
        
        ALTGraph space_part1(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
        
        ALTGraph space_part2(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        
        MotionGraph motion_graph1 = space_part1.get_adjacency_list();
        optimal_mg_vertex<FreeSpaceType> v1_p;
        v1_p.position = this->m_start_pos;
        v1_p.distance_accum = 0.0;
        create_root(v1_p, motion_graph1);
        
        MotionGraph motion_graph2 = space_part2.get_adjacency_list();
        optimal_mg_vertex<FreeSpaceType> v2_p;
        v2_p.position = this->m_goal_pos;
        v2_p.distance_accum = 0.0;
        create_root(v2_p, motion_graph2);
        
        multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrtstar_planner_visitor<FreeSpaceType, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro());
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, pos_map, get(random_sampler, this->m_space->get_super_space()), 
          nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_ALT_COB4_KNN) {
        
        typedef dvp_adjacency_list<
          optimal_mg_vertex<FreeSpaceType>,
          optimal_mg_edge<FreeSpaceType>,
          SuperSpace,
          PositionMap,
          4, random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4>,
          boost::vecS, boost::bidirectionalS, boost::listS > ALTGraph;
        
        ALTGraph space_part1(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
        
        ALTGraph space_part2(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        
        MotionGraph motion_graph1 = space_part1.get_adjacency_list();
        optimal_mg_vertex<FreeSpaceType> v1_p;
        v1_p.position = this->m_start_pos;
        v1_p.distance_accum = 0.0;
        create_root(v1_p, motion_graph1);
        
        MotionGraph motion_graph2 = space_part2.get_adjacency_list();
        optimal_mg_vertex<FreeSpaceType> v2_p;
        v2_p.position = this->m_goal_pos;
        v2_p.distance_accum = 0.0;
        create_root(v2_p, motion_graph2);
        
        multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrtstar_planner_visitor<FreeSpaceType, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro());
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, pos_map, get(random_sampler, this->m_space->get_super_space()), 
          nn_finder, this->m_max_vertex_count);
        
      };
      
    };
#endif
  };
  
#undef RK_RRTSTAR_PLANNER_INIT_START_AND_GOAL_NODE
#undef RK_RRTSTAR_PLANNER_CALL_RRTSTAR_FUNCTION
#undef RK_RRTSTAR_PLANNER_CALL_RRTSTAR_BNB_FUNCTION
#undef RK_RRTSTAR_PLANNER_CALL_APPROPRIATE_RRTSTAR_PLANNER_FUNCTION
  
  if(m_solutions.size())
    return m_solutions.begin()->second;
  else
    return shared_ptr< seq_path_base< SuperSpace > >();
};


};

};

#endif


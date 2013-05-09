/**
 * \file MEAQR_sbastar_planner.hpp
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

#ifndef REAK_MEAQR_SBASTAR_PLANNER_HPP
#define REAK_MEAQR_SBASTAR_PLANNER_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "path_planning/motion_planner_base.hpp"
#include "path_planning/sbmp_reporter_concept.hpp"

#include "path_planning/metric_space_concept.hpp"
#include "path_planning/seq_path_wrapper.hpp"
#include "interpolation/point_to_point_path.hpp"
#include "path_planning/basic_sbmp_reporters.hpp"

#include "graph_alg/lazy_sbastar.hpp"
#include "graph_alg/sbastar_rrtstar.hpp"
#include "graph_alg/anytime_sbastar.hpp"

#include "graph_alg/d_ary_bf_tree.hpp"
#include "graph_alg/d_ary_cob_tree.hpp"
#include "graph_alg/bgl_tree_adaptor.hpp"
#include "graph_alg/bgl_more_property_maps.hpp"
#include "path_planning/metric_space_search.hpp"
#include "path_planning/topological_search.hpp"
#include "path_planning/path_planner_options.hpp"
#include "graph_alg/neighborhood_functors.hpp"
#include "lin_alg/arithmetic_tuple.hpp"

#include "MEAQR_topology.hpp"
#include "topologies/fixed_topology_random_sampler.hpp"

#include "base/misc_math.hpp"

namespace ReaK {
  
namespace pp {
  
  
  

/**
 * This POD type contains the data required on a per-vertex basis for the 
 * Sampling-based A* path-planning algorithms that rely on a MEAQR topology.
 * \tparam StateSpace The topology type that represents the states of the underlying dynamic system.
 * \tparam StateSpaceSystem The type of the underlying dynamic system.
 * \tparam StateSpaceSampler The type of the sampler used on the state-space.
 */
template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
struct MEAQR_sbastar_vdata {
  typedef typename topology_traits< MEAQR_topology<StateSpace, StateSpaceSystem, StateSpaceSampler> >::point_type PointType;
  /// The position associated to the vertex.
  PointType position;
  /// The constriction associated to the vertex, which represents the probability that a sample drawn from the neighborhood of the vertex will be colliding (or unreachable by a collision-free path).
  double constriction;
  /// Keeps track of the number of neighbors of the vertex that could not be connected to it due to a collision.
  std::size_t collision_count;
  /// The density associated to the vertex, which represents the probability that a sample drawn from the neighborhood of the vertex will not yield any information gain.
  double density;
  /// Keeps track of the number of neighbors of the vertex.
  std::size_t expansion_trials;
  /// The heuristic-value associated to the vertex, i.e., the bird-flight distance to the goal.
  double heuristic_value;
  /// The travel-distance accumulated in the vertex, i.e., the travel-distance from the start vertex to this vertex.
  double distance_accum;
  /// The key-value associated to the vertex, computed by the SBA* algorithm.
  double key_value;
  /// The predecessor associated to the vertex, i.e., following the predecessor links starting at the goal node yields a backward trace of the optimal path.
  std::size_t predecessor;
  
  /**
   * Default constructor.
   */
  MEAQR_sbastar_vdata() : position(PointType()),
                          constriction(0.2), collision_count(1), density(0.0), expansion_trials(0),
                          heuristic_value(0.0), distance_accum(0.0), key_value(0.0), predecessor(0) { };
};

/**
 * This POD type contains the data required on a per-edge basis for the 
 * Sampling-based A* path-planning algorithms that rely on a MEAQR topology.
 * \tparam StateSpace The topology type that represents the states of the underlying dynamic system.
 * \tparam StateSpaceSystem The type of the underlying dynamic system.
 * \tparam StateSpaceSampler The type of the sampler used on the state-space.
 */
template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
struct MEAQR_sbastar_edata {
  typedef typename MEAQR_topology<StateSpace, StateSpaceSystem, StateSpaceSampler>::steer_record_type steer_record_type;
  
  /// The cost-to-go associated to the edge (from source to target).
  double astar_weight;
  /// The steering-record associated to the edge, because the MEAQR topology is a steerable-space, and thus, records the path between two points (which is non-trivial).
  steer_record_type steer_rec;
  
  /**
   * Default constructor.
   * \param aWeight The cost-to-go to be associated to this edge.
   * \param aSteerRec The steering-record to be associated to this edge, because the MEAQR topology is a steerable-space, and thus, records the path between two points (which is non-trivial).
   */
  MEAQR_sbastar_edata(double aWeight = 0.0, const steer_record_type& aSteerRec = steer_record_type()) : astar_weight(aWeight), steer_rec(aSteerRec) { };
#ifdef RK_ENABLE_CXX11_FEATURES
  /**
   * Default constructor with move-semantics.
   * \param aWeight The cost-to-go to be associated to this edge.
   * \param aSteerRec The steering-record to be moved into to this edge, because the MEAQR topology is a steerable-space, and thus, records the path between two points (which is non-trivial).
   */
  MEAQR_sbastar_edata(double aWeight, steer_record_type&& aSteerRec) : astar_weight(aWeight), steer_rec(std::move(aSteerRec)) { };
#endif
};




/**
 * This class is a SBA*-based path-planner over the a MEAQR-controlled system-topology.
 * \tparam StateSpace The topology type of state-space of the dynamic system under MEAQR control.
 * \tparam StateSpaceSystem The type of the dynamic system under MEAQR control.
 * \tparam SBPPReporter The reporter type to use to report the progress of the path-planning.
 */
template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler,
          typename SBPPReporter = no_sbmp_report>
class MEAQR_sbastar_planner : public sample_based_planner< path_planner_base< MEAQR_topology_with_CD<StateSpace, StateSpaceSystem, StateSpaceSampler> > > {
  public:
    typedef MEAQR_topology_with_CD<StateSpace, StateSpaceSystem, StateSpaceSampler> space_type;
    typedef typename subspace_traits<space_type>::super_space_type super_space_type;
    
    typedef sample_based_planner< path_planner_base< space_type > > base_type;
    typedef MEAQR_sbastar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler, SBPPReporter> self;
    
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename topology_traits< super_space_type >::point_difference_type point_difference_type;
    
  protected:
    SBPPReporter m_reporter;
    point_type m_start_pos;
    point_type m_goal_pos;
    double m_init_key_threshold;
    double m_init_dens_threshold;
    double m_init_relaxation;
    double m_SA_init_temperature;
    double m_sampling_radius;
    std::size_t m_current_num_results;
    std::size_t max_num_results;
    bool has_reached_max_vertices;
    std::size_t m_knn_flag;
    std::size_t m_collision_check_flag;
    std::size_t m_added_bias_flags;
    
    double m_current_key_threshold;
    double m_current_dens_threshold;
    double m_current_relaxation;
    
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
     * Returns true if the solver should keep on going trying to solve the path-planning problem.
     * \note This function is for internal use by the path-planning algorithm (a visitor callback).
     * \return True if the solver should keep on going trying to solve the path-planning problem.
     */
    bool keep_going() const {
      return (max_num_results > m_current_num_results) && !has_reached_max_vertices;
    };
    
    /**
     * This function invokes the path-planning reporter to report on the progress of the path-planning
     * solver.
     * \note This function is for internal use by the path-planning algorithm (a visitor callback).
     * \param g The current motion-graph.
     */
    template <typename Graph, typename PositionMap>
    void report_progress(Graph& g, PositionMap g_pos) {
      if(num_vertices(g) % this->m_progress_interval == 0)
        m_reporter.draw_motion_graph(*(this->m_space), g, g_pos);
      has_reached_max_vertices = (num_vertices(g) >= this->m_max_vertex_count);
    };
    
    /**
     * This function computes the heuristic distance value of a given node, i.e., the bird-flight 
     * distance to the goal position.
     * \note This function is used internally by the path-planning algorithm (a visitor callback).
     * \param u The node for which the heuristic value is sought.
     * \param g The current motion-graph.
     * \return The heuristic distance value of the given node.
     */
    template <typename Graph>
    double heuristic(typename boost::graph_traits<Graph>::vertex_descriptor u, const Graph& g) const {
      return get(distance_metric, this->m_space->get_super_space())(g[u].position, this->m_goal_pos, this->m_space->get_super_space());
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
        shared_ptr< seq_path_wrapper< point_to_point_path<super_space_type> > > new_sol(new seq_path_wrapper< point_to_point_path<super_space_type> >("MEAQR_sbastar_solution", point_to_point_path<super_space_type>(sup_space_ptr,get(distance_metric, this->m_space->get_super_space()))));
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
          ++m_current_num_results;
          m_solutions[goal_distance] = new_sol;
          m_reporter.draw_solution(*(this->m_space), m_solutions[goal_distance]);
        };
      };
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
     * Returns the initial key-value threshold used by this planner.
     * \return The initial key-value threshold used by this planner.
     */
    double get_initial_key_threshold() const { return m_init_key_threshold; };
    /**
     * Sets the initial key-value threshold to be used by this planner.
     * \param aInitialThreshold The initial key-value threshold to be used by this planner.
     */
    void set_initial_key_threshold(double aInitialThreshold) { m_init_key_threshold = aInitialThreshold; };
    
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
     * Returns the initial Simulated Annealing temperature used by this planner.
     * \return The initial Simulated Annealing temperature used by this planner.
     */
    double get_initial_SA_temperature() const { return m_SA_init_temperature; };
    /**
     * Sets the initial Simulated Annealing temperature to be used by this planner.
     * \param aInitialSATemperature The initial Simulated Annealing temperature to be used by this planner.
     */
    void set_initial_SA_temperature(double aInitialSATemperature) { m_SA_init_temperature = aInitialSATemperature; };
    
    /**
     * Returns the current key-value threshold used by this planner.
     * \return The current key-value threshold used by this planner.
     */
    double get_current_key_threshold() const { return m_current_key_threshold; };
    /**
     * Sets the current key-value threshold to be used by this planner.
     * \param aCurrentThreshold The current key-value threshold to be used by this planner.
     */
    void set_current_key_threshold(double aCurrentThreshold) { m_current_key_threshold = aCurrentThreshold; };
    
    /**
     * Returns the current density-value threshold used by this planner.
     * \return The current density-value threshold used by this planner.
     */
    double get_current_density_threshold() const { return m_current_dens_threshold; };
    /**
     * Sets the current density-value threshold to be used by this planner.
     * \param aCurrentThreshold The current density-value threshold to be used by this planner.
     */
    void set_current_density_threshold(double aCurrentThreshold) { m_current_dens_threshold = aCurrentThreshold; };
    
    /**
     * Returns the current relaxation factor used by this planner.
     * \return The current relaxation factor used by this planner.
     */
    double get_current_relaxation() const { return m_current_relaxation; };
    /**
     * Sets the current relaxation factor to be used by this planner.
     * \param aCurrentRelaxation The current relaxation factor to be used by this planner.
     */
    void set_current_relaxation(double aCurrentRelaxation) { m_current_relaxation = aCurrentRelaxation; };
    
    
    /**
     * Returns the sampling radius (in the topology's distance metric) used by this planner.
     * \return The sampling radius used by this planner.
     */
    double get_sampling_radius() const { return m_sampling_radius; };
    /**
     * Sets the sampling radius (in the topology's distance metric) to be used by this planner.
     * \param aSamplingRadius The sampling radius to be used by this planner.
     */
    void set_sampling_radius(double aSamplingRadius) { m_sampling_radius = aSamplingRadius; };
    
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
     * Returns the integer flag that identifies the kind of K-nearest-neighbor method to use (see path_planner_options.hpp).
     * \return The integer flag that identifies the kind of K-nearest-neighbor method to use (see path_planner_options.hpp).
     */
    std::size_t get_knn_flag() const { return m_knn_flag; };
    /**
     * Sets the integer flag that identifies the kind of K-nearest-neighbor method to use (see path_planner_options.hpp).
     * \param aKNNMethodFlag The integer flag that identifies the kind of K-nearest-neighbor method to use (see path_planner_options.hpp).
     */
    void set_knn_flag(std::size_t aKNNMethodFlag) { m_knn_flag = aKNNMethodFlag; };
    
    /**
     * Returns the integer flag that identifies the kind of collision-checking to use (eager or lazy) (see path_planner_options.hpp).
     * \return The integer flag that identifies the kind of collision-checking to use (eager or lazy) (see path_planner_options.hpp).
     */
    std::size_t get_collision_check_flag() const { return m_collision_check_flag; };
    /**
     * Sets the integer flag that identifies the kind of collision-checking to use (eager or lazy) (see path_planner_options.hpp).
     * \param aCollisionCheckFlag The integer flag that identifies the kind of collision-checking to use (eager or lazy) (see path_planner_options.hpp).
     */
    void set_collision_check_flag(std::size_t aCollisionCheckFlag) { m_collision_check_flag = aCollisionCheckFlag; };
    
    /**
     * Returns the integer flag that identifies the kind of added bias to use (no-bias, RRT*, etc.) (see path_planner_options.hpp).
     * \return The integer flag that identifies the kind of added bias to use (no-bias, RRT*, etc.) (see path_planner_options.hpp).
     */
    std::size_t get_added_bias_flags() const { return m_added_bias_flags; };
    /**
     * Sets the integer flag that identifies the kind of added bias to use (no-bias, RRT*, etc.) (see path_planner_options.hpp).
     * \param aAddedBiasFlags The integer flag that identifies the kind of added bias to use (no-bias, RRT*, etc.) (see path_planner_options.hpp).
     */
    void set_added_bias_flags(std::size_t aAddedBiasFlags) { m_added_bias_flags = aAddedBiasFlags; };
    
    
    /**
     * Parametrized constructor.
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     * \param aStartPos The position value of the starting location.
     * \param aGoalPos The position value of the goal location.
     * \param aInitialKeyThreshold The initial key-value threshold for exploring nodes, should be somewhere between 0 and 1.
     * \param aInitialDensityThreshold The initial density-value threshold for exploring nodes, should be somewhere between 0 (reject no points based on density) and 1 (reject all points, requires zero density).
     * \param aInitialRelaxation The initial relaxation factor for exploring nodes, should be greater than 0 (if equal to or less than 0, then the plain SBA* algorithms are used instead of the Anytime SBA* algorithms).
     * \param aSamplingRadius The radius of the sampled space around a given point when doing random walks.
     * \param aMaxVertexCount The maximum number of samples to generate during the motion planning.
     * \param aProgressInterval The number of new samples between each "progress report".
     * \param aBiDirFlag An integer flag representing the directionality of the RRT algorithm 
     *                   used (either UNIDIRECTIONAL_RRT or BIDIRECTIONAL_RRT).
     * \param aGraphKindFlag An integer flag representing the kind of motion graph to use in the 
     *                       RRT algorithm. Can be ADJ_LIST_MOTION_GRAPH or DVP_ADJ_LIST_MOTION_GRAPH.
     * \param aKNNMethodFlag An integer flag representing the kind of KNN method to use for nearest
     *                       neighbor queries in the graph. Can be LINEAR_SEARCH_KNN, DVP_BF2_TREE_KNN,
     *                       DVP_BF4_TREE_KNN, DVP_COB2_TREE_KNN, or DVP_COB4_TREE_KNN when the 
     *                       motion graph is of kind ADJ_LIST_MOTION_GRAPH. Can be DVP_ALT_BF2_TREE_KNN,
     *                       DVP_ALT_BF4_TREE_KNN, DVP_ALT_COB2_TREE_KNN, or DVP_ALT_COB4_TREE_KNN when 
     *                       the motion graph is of kind DVP_ADJ_LIST_MOTION_GRAPH.
     * \param aCollisionCheckFlag The integer flag that identifies the kind of collision-checking to use (eager or lazy) (see path_planner_options.hpp).
     * \param aAddedBiasFlags The integer flag that identifies the kind of added bias to use (no-bias, RRT*, etc.) (see path_planner_options.hpp).
     * \param aReporter The SBPP reporter object to use to report results and progress.
     * \param aMaxResultCount The maximum number of successful start-goal connections to make before 
     *                        stopping the path planner (the higher the number the more likely that a 
     *                        good path will be found, however, running time can become much longer).
     * \param aInitialSATemperature The initial Simulated Annealing temperature use by this planner, if the 
     *                              added-bias is set to an exploratory bias (e.g., PLAN_WITH_VORONOI_PULL). If negative,
     *                              then simulated annealing is not used, and the exploratory bias (if any) is applied 
     *                              only when SBA* seaching stalls (isn't progressing anymore).
     */
    MEAQR_sbastar_planner(const shared_ptr< space_type >& aWorld = shared_ptr< space_type >(), 
                          const point_type& aStartPos = point_type(),
                          const point_type& aGoalPos = point_type(),
                          double aInitialKeyThreshold = 0.8,
                          double aInitialDensityThreshold = 0.8,
                          double aInitialRelaxation = 0.0,
                          double aSamplingRadius = 1.0,
                          std::size_t aMaxVertexCount = 5000, 
                          std::size_t aProgressInterval = 100,
                          std::size_t aKNNMethodFlag = DVP_BF2_TREE_KNN,
                          std::size_t aCollisionCheckFlag = LAZY_COLLISION_CHECKING,
                          std::size_t aAddedBiasFlags = NOMINAL_PLANNER_ONLY,
                          SBPPReporter aReporter = SBPPReporter(),
                          std::size_t aMaxResultCount = 50,
                          double aInitialSATemperature = -1.0) :
                          base_type("MEAQR_sbastar_planner", aWorld, aMaxVertexCount, aProgressInterval),
                          m_reporter(aReporter),
                          m_start_pos(aStartPos),
                          m_goal_pos(aGoalPos),
                          m_init_key_threshold(aInitialKeyThreshold),
                          m_init_dens_threshold(aInitialDensityThreshold),
                          m_init_relaxation(aInitialRelaxation),
                          m_SA_init_temperature(aInitialSATemperature),
                          m_sampling_radius(aSamplingRadius),
                          m_current_num_results(0),
                          max_num_results(aMaxResultCount),
                          has_reached_max_vertices(false),
                          m_knn_flag(aKNNMethodFlag),
                          m_collision_check_flag(aCollisionCheckFlag), 
                          m_added_bias_flags(aAddedBiasFlags),
                          m_current_key_threshold(m_init_key_threshold),
                          m_current_dens_threshold(m_init_dens_threshold),
                          m_current_relaxation(m_init_relaxation) { };
    
    virtual ~MEAQR_sbastar_planner() { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_reporter)
        & RK_SERIAL_SAVE_WITH_NAME(m_start_pos)
        & RK_SERIAL_SAVE_WITH_NAME(m_goal_pos)
        & RK_SERIAL_SAVE_WITH_NAME(m_init_key_threshold)
        & RK_SERIAL_SAVE_WITH_NAME(m_init_dens_threshold)
        & RK_SERIAL_SAVE_WITH_NAME(m_init_relaxation)
        & RK_SERIAL_SAVE_WITH_NAME(m_SA_init_temperature)
        & RK_SERIAL_SAVE_WITH_NAME(m_sampling_radius)
        & RK_SERIAL_SAVE_WITH_NAME(max_num_results)
        & RK_SERIAL_SAVE_WITH_NAME(m_knn_flag)
        & RK_SERIAL_SAVE_WITH_NAME(m_collision_check_flag)
        & RK_SERIAL_SAVE_WITH_NAME(m_added_bias_flags);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_reporter)
        & RK_SERIAL_LOAD_WITH_NAME(m_start_pos)
        & RK_SERIAL_LOAD_WITH_NAME(m_goal_pos)
        & RK_SERIAL_LOAD_WITH_NAME(m_init_key_threshold)
        & RK_SERIAL_LOAD_WITH_NAME(m_init_dens_threshold)
        & RK_SERIAL_LOAD_WITH_NAME(m_init_relaxation)
        & RK_SERIAL_LOAD_WITH_NAME(m_SA_init_temperature)
        & RK_SERIAL_LOAD_WITH_NAME(m_sampling_radius)
        & RK_SERIAL_LOAD_WITH_NAME(max_num_results)
        & RK_SERIAL_LOAD_WITH_NAME(m_knn_flag)
        & RK_SERIAL_LOAD_WITH_NAME(m_collision_check_flag)
        & RK_SERIAL_LOAD_WITH_NAME(m_added_bias_flags);
      has_reached_max_vertices = false;
      m_solutions.clear();
      m_current_num_results = 0;
      m_current_key_threshold = m_init_key_threshold;
      m_current_dens_threshold = m_init_dens_threshold;
      m_current_relaxation = m_init_relaxation;
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC246000E,1,"MEAQR_sbastar_planner",base_type)
};




/**
 * This class template is used by the MEAQR SBA* path-planner as the visitor object needed to 
 * collaborate with the SBA* algorithms to generate the motion-graph and path-planning solutions.
 * This class template models the SBAStarVisitorConcept and SBARRTStarVisitorConcept.
 * As with most planning algorithms in ReaK, the algorithm is really made up of a high-level 
 * algorithmic logic in the form of function templates, and a number of customization points 
 * collected as member functions of an algorithm visitor class that implement the problem-specific 
 * behaviors (random-walks / local-planning, heuristic computation, progress reporting, 
 * completion criteria, etc.). 
 */
template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler, typename MotionGraph, typename NNFinderSynchro, typename SBPPReporter = no_sbmp_report>
struct MEAQR_sbastar_visitor {
  typedef MEAQR_topology_with_CD<StateSpace, StateSpaceSystem, StateSpaceSampler> space_type;
  typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
  
  shared_ptr< space_type > m_space;
  MEAQR_sbastar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler, SBPPReporter>* m_planner;
  NNFinderSynchro m_nn_synchro;
  Vertex m_start_node;
  Vertex m_goal_node;
  double m_space_dim;
  double m_space_Lc;
  
  MEAQR_sbastar_visitor(const shared_ptr< space_type >& aSpace, 
                        MEAQR_sbastar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler, SBPPReporter>* aPlanner,
                        NNFinderSynchro aNNSynchro,
                        Vertex aStartNode, Vertex aGoalNode, double aSpaceDim, double aSpaceLc) : 
                        m_space(aSpace), m_planner(aPlanner), m_nn_synchro(aNNSynchro),
                        m_start_node(aStartNode), m_goal_node(aGoalNode), 
                        m_space_dim(aSpaceDim), m_space_Lc(aSpaceLc) { };
  
  typedef typename topology_traits< space_type >::point_type PointType;
  typedef MEAQR_sbastar_edata<StateSpace, StateSpaceSystem, StateSpaceSampler> EdgeProp;
  
  
  
  
  template <typename Vertex, typename Graph>
  void init_vertex_properties(Vertex u, Graph& g) const {
    g[u].heuristic_value = m_planner->heuristic(u,g);
    g[u].constriction = 0.0;
    g[u].collision_count = 0;  // r
    g[u].density = 0.0;
    g[u].expansion_trials = 0;  // m
  };
  
  // compute sample similarity directly as KL-divergence between Gaussians:
  
  double compute_sample_similarity(double travel_dist, double event_radius) const {
    using std::exp; using std::log;
    
    double sig2_n = 0.25 * m_planner->get_sampling_radius() * m_planner->get_sampling_radius();
    double sig2_x = 0.25 * event_radius * event_radius;
    return exp(-travel_dist * travel_dist / (sig2_x * 2.0) - 0.5 * m_space_dim * ( sig2_n / sig2_x - 1.0 - log(sig2_n / sig2_x) ) );
  };
  
  double compute_sample_similarity(double travel_dist) const {
    using std::exp;
    
    return exp(-travel_dist * travel_dist / (0.25 * m_planner->get_sampling_radius() * m_planner->get_sampling_radius() * 2.0));
  };
  
  // keep the sample similarity weighted by the sample probability (and its binomial converse).
  // that is, assume the existing density to reflect the overall density and the newly computed 
  // sample similarity to reflect the density in its relatively probable region (binomial).
  
  template <typename Vertex, typename Graph>
  void register_explored_sample(Vertex u, Graph& g, double samp_sim) const {
    g[u].density = g[u].density * (1.0 - samp_sim) + samp_sim * samp_sim;
  };
  template <typename Vertex, typename Graph>
  void register_failed_sample(Vertex u, Graph& g, double samp_sim) const {
    g[u].constriction = g[u].constriction * (1.0 - samp_sim) + samp_sim * samp_sim;
  };
  
  
  
  
  
  
  
  template <typename Vertex, typename Graph>
  void vertex_added(Vertex u, Graph& g) const {
    m_nn_synchro.added_vertex(u,g);
    
    init_vertex_properties(u,g);
    
    // Call progress reporter...
    m_planner->report_progress(g, get(&EdgeProp::steer_rec,g));
    
    if((in_degree(m_goal_node,g)) && (g[m_goal_node].distance_accum < m_planner->get_best_solution_distance()))
      m_planner->create_solution_path(m_start_node, m_goal_node, g);
    
  };
  
  template <typename EdgeType, typename Graph>
  void edge_added(EdgeType e, Graph& g) const { };
  
  template <typename Vertex, typename Graph>
  void travel_explored(Vertex u, Vertex v, Graph& g) const { 
    double dist = get(distance_metric, m_space->get_super_space())(g[u].position, g[v].position, m_space->get_super_space());
    double samp_sim = compute_sample_similarity(dist);
    if(samp_sim < std::numeric_limits<double>::epsilon())
      return;
    register_explored_sample(u,g,samp_sim);
    register_explored_sample(v,g,samp_sim);
  };
  
  template <typename Vertex, typename Graph>
  void travel_succeeded(Vertex, Vertex, Graph&) const { 
    // nothing to do (record explorations and failures only).
  };
  
  template <typename Vertex, typename Graph>
  void travel_failed(Vertex u, Vertex v, Graph& g) const { 
    
    using std::log; using std::exp;
    double dist = get(distance_metric, m_space->get_super_space())(g[u].position, g[v].position, m_space->get_super_space());
    
    // assume collision occured half-way.  
    // assume a blob of collision half-way and occupying half of the interval between u and v (radius 1/4 of travel_dist).
    double samp_sim = compute_sample_similarity(0.5 * dist, 0.25 * dist);
    if(samp_sim < std::numeric_limits<double>::epsilon())
      return;
    register_failed_sample(u,g,samp_sim);
    register_failed_sample(v,g,samp_sim);
  };
      
  
  bool keep_going() const {
    if( ( m_planner->get_initial_relaxation() > 1e-6 ) && ( m_planner->get_current_relaxation() < 1e-6 ) )
      return false;
    return m_planner->keep_going();
  };
  
  // NOTE: This is the main thing that must be revised:
  template <typename Vertex, typename Graph>
  boost::tuple<PointType, bool, EdgeProp > random_walk(Vertex u, Graph& g) const {
    typedef boost::tuple<PointType, bool, EdgeProp > ResultType;
    typedef typename EdgeProp::steer_record_type SteerRec;
    typedef typename topology_traits<StateSpace>::point_difference_type state_difference_type;
    using std::exp;
    using std::log;
    using ReaK::to_vect;
    using ReaK::from_vect;
    
    typename point_distribution_traits< 
      typename subspace_traits<space_type>::super_space_type >::random_sampler_type 
      get_sample = get(random_sampler, m_space->get_super_space());
    typename metric_space_traits< 
      typename subspace_traits<space_type>::super_space_type >::distance_metric_type 
      get_distance = get(distance_metric, m_space->get_super_space());
    
    boost::variate_generator< pp::global_rng_type&, boost::normal_distribution<double> > var_rnd(pp::get_global_rng(), boost::normal_distribution<double>());
    
    unsigned int i = 0;
    do {
      
      PointType p_rnd = get_sample(m_space->get_super_space());
      double dist = get_distance(g[u].position, p_rnd, m_space->get_super_space());
      double target_dist = var_rnd() * m_planner->get_sampling_radius();
      
      PointType p_inter( m_space->get_state_space().move_position_toward(g[u].position.x, target_dist / dist, p_rnd.x) );
      target_dist = get_distance(g[u].position, p_inter, m_space->get_super_space());
      
      std::pair<PointType, SteerRec> steer_result = m_space->steer_position_toward(g[u].position, 0.8, p_inter);
      dist = get_distance(g[u].position, steer_result.first, m_space->get_super_space());
      if( dist < 0.7 * target_dist ) {
        // this means that we had a collision before reaching the target distance, 
        // must record that to the constriction statistic:
        double samp_sim = compute_sample_similarity(target_dist, (target_dist - dist));
        if(samp_sim > std::numeric_limits<double>::epsilon())
          register_failed_sample(u,g,samp_sim);
      } else
#ifdef RK_ENABLE_CXX11_FEATURES
        return ResultType(steer_result.first, true, EdgeProp(dist, std::move(steer_result.second)));
#else
        return ResultType(steer_result.first, true, EdgeProp(dist, steer_result.second));
#endif
    } while(++i <= 10);
    return ResultType(g[u].position, false, EdgeProp(std::numeric_limits<double>::infinity()));
  };
  
  
  template <typename Vertex, typename Graph>
  boost::tuple<PointType, bool, EdgeProp> steer_towards_position(const PointType& p, Vertex u, Graph& g) const {
    typedef typename EdgeProp::steer_record_type SteerRec;
    typedef boost::tuple<PointType, bool, EdgeProp> ResultType;
    
    // First, try to bring the state-space point within the time-horizon:
    double total_dist = get(distance_metric, m_space->get_super_space())(g[u].position, p, m_space->get_super_space());
    double max_cost_to_go = 0.75 * m_space->get_max_time_horizon() * m_space->get_idle_power_cost(g[u].position);
    PointType p_dest = p;
    while(total_dist > max_cost_to_go) {
      p_dest = PointType( m_space->get_state_space().move_position_toward(g[u].position.x, 0.5, p_dest.x) );
      total_dist = get(distance_metric, m_space->get_super_space())(g[u].position, p_dest, m_space->get_super_space());
    };
    
    // Then, steer to that point, recording the path of the steering function.
    std::pair<PointType, SteerRec> steer_result = m_space->steer_position_toward(g[u].position, 0.8, p_dest);
    
    // Check if the progress in the state-space was significant (at least 0.1 of the best-case).
    double best_case_dist = get(distance_metric, m_space->get_state_space())(g[u].position.x, p_dest.x, m_space->get_state_space());
    double actual_dist = get(distance_metric, m_space->get_state_space())(g[u].position.x, steer_result.first.x, m_space->get_state_space());
    
    if(actual_dist > 0.1 * best_case_dist) {
//       std::cout << "Steered successfully!" << std::endl;
#ifdef RK_ENABLE_CXX11_FEATURES
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
    
    std::pair<PointType, SteerRec> steer_result = m_space->steer_position_toward(g[u].position, 1.0, g[v].position);
    
    // NOTE Differs from rrtstar_path_planner HERE:
    double best_case_dist = get(distance_metric, m_space->get_state_space())(g[u].position.x, g[v].position.x, m_space->get_state_space());
    double diff_dist = get(distance_metric, m_space->get_state_space())(steer_result.first.x, g[v].position.x, m_space->get_state_space());
    
    if(diff_dist < 0.05 * best_case_dist) {
//       std::cout << "Connected successfully!" << std::endl;
      return ResultType(true, EdgeProp(
        get(distance_metric, m_space->get_super_space())(g[u].position, g[v].position, m_space->get_super_space()),
#ifdef RK_ENABLE_CXX11_FEATURES
        std::move(steer_result.second)
#else
        steer_result.second
#endif
      ));
    } else {
      return ResultType(false,EdgeProp());
    };
  };
  
  
  template <typename Vertex, typename Graph>
  void initialize_vertex(Vertex u, Graph& g) const {
    typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
    
    init_vertex_properties(u,g);
    
    OutEdgeIter ei, ei_end;
    for(boost::tie(ei,ei_end) = out_edges(u,g); ei != ei_end; ++ei) {
      double samp_sim = compute_sample_similarity(g[*ei].astar_weight);
      if(samp_sim > std::numeric_limits<double>::epsilon())
        register_explored_sample(u,g,samp_sim);
    };
    
  };
  
  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex, const Graph&) const { };
  template <typename Vertex, typename Graph>
  void examine_vertex(Vertex, const Graph&) const { };
  template <typename Edge, typename Graph>
  void edge_relaxed(Edge, const Graph&) const { };
  
  template <typename Graph>
  void publish_path(const Graph& g) const {
//     std::cout << "Publishing path..." << std::endl;
    // try to create a goal connection path
    if((in_degree(m_goal_node,g)) && (g[m_goal_node].distance_accum < m_planner->get_best_solution_distance()))
      m_planner->create_solution_path(m_start_node, m_goal_node, g); 
    
//     m_planner->set_current_key_threshold( 0.75 * m_planner->get_current_key_threshold() );
//     m_planner->set_current_density_threshold( 0.95 * m_planner->get_current_density_threshold() );
    
//     std::cout << " new key-value threshold =\t" << m_planner->get_current_key_threshold() << std::endl;
//     std::cout << " new density-value threshold =\t" << m_planner->get_current_density_threshold() << std::endl;
  };
  
  template <typename Vertex, typename Graph>
  bool has_search_potential(Vertex u, const Graph& g) const { 
    if( m_planner->get_initial_relaxation() > 1e-6 ) {
      // assume we are running a Anytime SBA* algorithm.
      if(in_degree(m_goal_node,g) && ( g[u].heuristic_value < m_planner->get_current_key_threshold() * (g[u].distance_accum + g[u].heuristic_value) / (1.0 - g[u].constriction) / (1.0 - g[u].density) ))
        return false;
      else 
        return true;
    };
    if( (m_planner->get_current_key_threshold() < 1e-6) || (g[u].key_value < m_space_Lc / m_planner->get_current_key_threshold()) )
      return true;
    else 
      return false;
  };
  
  template <typename Vertex, typename Graph>
  bool should_close(Vertex u, const Graph& g) const {
    if(u == m_goal_node)  // never enqueue the goal node.
      return true;
    if( m_planner->get_initial_relaxation() > 1e-6 ) {
      // assume we are running a Anytime SBA* algorithm.
      if((1.0 - g[u].constriction) * (1.0 - g[u].density) < m_planner->get_current_density_threshold())
        return true;
      else 
        return false;
    };
    if((1.0 - g[u].constriction) * (1.0 - g[u].density) < m_planner->get_current_density_threshold())
      return true;
    else 
      return false;
  };
  
  template <typename Graph>
  double adjust_relaxation(double, const Graph& g) const {
    m_planner->set_current_relaxation(m_planner->get_current_relaxation() * 0.5);
//     std::cout << " new relaxation =\t" << m_planner->get_current_relaxation() << std::endl;
    return m_planner->get_current_relaxation();
  };
  
};






template <typename StateSpace, 
          typename StateSpaceSystem, typename StateSpaceSampler, 
          typename SBPPReporter>
shared_ptr< seq_path_base< typename MEAQR_sbastar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler, SBPPReporter>::super_space_type > > 
  MEAQR_sbastar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler, SBPPReporter>::solve_path() {
  using ReaK::to_vect;
  
  this->has_reached_max_vertices = false;
  this->m_current_num_results = 0;
  this->m_solutions.clear();
  this->m_current_key_threshold = this->m_init_key_threshold;
  this->m_current_dens_threshold = this->m_init_dens_threshold;
  this->m_current_relaxation = this->m_init_relaxation;
  
  
  typedef typename MEAQR_sbastar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler, SBPPReporter>::space_type SpaceType;
  typedef typename MEAQR_sbastar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler, SBPPReporter>::super_space_type SuperSpace;
  typedef typename SuperSpace::IHAQR_space_type IHAQRSpace;
  typedef typename topology_traits<SuperSpace>::point_type PointType;
  typedef MEAQR_sbastar_vdata<StateSpace, StateSpaceSystem, StateSpaceSampler> VertexProp;
  typedef MEAQR_sbastar_edata<StateSpace, StateSpaceSystem, StateSpaceSampler> EdgeProp;
  
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
  WeightMap weight_map = WeightMap(&EdgeProp::astar_weight);
  
  double space_dim = double((to_vect<double>(this->m_space->get_state_space().difference(this->m_goal_pos.x,this->m_start_pos.x))).size()); 
  double space_Lc = get(distance_metric,this->m_space->get_super_space())(this->m_start_pos, this->m_goal_pos, this->m_space->get_super_space());
  
  double max_radius = 2.0 * m_sampling_radius;
    
  typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::bidirectionalS,
                                 VertexProp,
                                 EdgeProp, 
                                 boost::listS> MotionGraphType;
  typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex;
  typedef boost::composite_property_map< 
    PositionMap, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t > > GraphPositionMap;
  
  MotionGraphType motion_graph;
  GraphPositionMap g_pos_map = GraphPositionMap(pos_map, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t >(&motion_graph));
  
  VertexProp vs_p, vg_p;
  vs_p.position = this->m_start_pos;
  vg_p.position = this->m_goal_pos;
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  Vertex start_node = add_vertex(std::move(vs_p), motion_graph);
  Vertex goal_node  = add_vertex(std::move(vg_p), motion_graph);
#else
  Vertex start_node = add_vertex(vs_p, motion_graph);
  Vertex goal_node  = add_vertex(vg_p, motion_graph);
#endif
  motion_graph[start_node].constriction     = 0.0;
  motion_graph[start_node].collision_count  = 0;
  motion_graph[start_node].density          = 0.0;
  motion_graph[start_node].expansion_trials = 0;
  motion_graph[start_node].heuristic_value  = heuristic(start_node,motion_graph);  // distance to goal node.
  motion_graph[start_node].distance_accum   = 0.0;
  motion_graph[start_node].key_value        = space_Lc + this->m_init_relaxation * space_Lc;
  motion_graph[start_node].predecessor      = start_node;
  
  motion_graph[goal_node].constriction      = 0.0;
  motion_graph[goal_node].collision_count   = 0;
  motion_graph[goal_node].density           = 0.0;
  motion_graph[goal_node].expansion_trials  = 0;
  motion_graph[goal_node].heuristic_value   = 0.0;
  motion_graph[goal_node].distance_accum    = std::numeric_limits<double>::infinity();
  motion_graph[goal_node].key_value         = std::numeric_limits<double>::infinity();
  motion_graph[goal_node].predecessor       = goal_node;
  
  
#define RK_MEAQR_SBASTAR_PLANNER_CALL_SBASTAR_FUNCTION \
    ReaK::graph::generate_sbastar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, start_node, this->m_space->get_super_space(), vis, \
        heuristic_map,  \
        pos_map,  \
        weight_map, \
        dens_map,  \
        cons_map,  \
        dist_map, \
        pred_map,  \
        get(&VertexProp::key_value, motion_graph),  \
        nc_selector)\
      );
  
  
#define RK_MEAQR_SBASTAR_PLANNER_CALL_LAZY_SBASTAR_FUNCTION \
    ReaK::graph::generate_lazy_sbastar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, start_node, this->m_space->get_super_space(), vis, \
        heuristic_map,  \
        pos_map,  \
        weight_map, \
        dens_map,  \
        cons_map,  \
        dist_map, \
        pred_map,  \
        get(&VertexProp::key_value, motion_graph),  \
        nc_selector)\
      );
  
  
#define RK_MEAQR_SBASTAR_PLANNER_CALL_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_sbarrtstar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, start_node, this->m_space->get_super_space(), vis, \
        heuristic_map,  \
        pos_map,  \
        weight_map, \
        dens_map,  \
        cons_map,  \
        dist_map, \
        pred_map,  \
        get(&VertexProp::key_value, motion_graph),  \
        nc_selector), \
      get(random_sampler, this->m_space->get_super_space()), \
      this->m_SA_init_temperature);
  
  
#define RK_MEAQR_SBASTAR_PLANNER_CALL_LAZY_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_lazy_sbarrtstar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, start_node, this->m_space->get_super_space(), vis, \
        heuristic_map,  \
        pos_map,  \
        weight_map, \
        dens_map,  \
        cons_map,  \
        dist_map, \
        pred_map,  \
        get(&VertexProp::key_value, motion_graph),  \
        nc_selector), \
      get(random_sampler, this->m_space->get_super_space()), \
      this->m_SA_init_temperature);
   
   
   
   
#define RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_SBASTAR_FUNCTION \
    ReaK::graph::generate_anytime_sbastar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, start_node, this->m_space->get_super_space(), vis, \
        heuristic_map,  \
        pos_map,  \
        weight_map, \
        dens_map,  \
        cons_map,  \
        dist_map, \
        pred_map,  \
        get(&VertexProp::key_value, motion_graph),  \
        nc_selector), \
      this->m_init_relaxation);
  
  
#define RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBASTAR_FUNCTION \
    ReaK::graph::generate_anytime_lazy_sbastar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, start_node, this->m_space->get_super_space(), vis, \
        heuristic_map,  \
        pos_map,  \
        weight_map, \
        dens_map,  \
        cons_map,  \
        dist_map, \
        pred_map,  \
        get(&VertexProp::key_value, motion_graph),  \
        nc_selector), \
      this->m_init_relaxation);
  
  
#define RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_anytime_sbarrtstar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, start_node, this->m_space->get_super_space(), vis, \
        heuristic_map,  \
        pos_map,  \
        weight_map, \
        dens_map,  \
        cons_map,  \
        dist_map, \
        pred_map,  \
        get(&VertexProp::key_value, motion_graph),  \
        nc_selector), \
      get(random_sampler, this->m_space->get_super_space()), \
      this->m_init_relaxation, this->m_SA_init_temperature);
  
  
#define RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_anytime_lazy_sbarrtstar( \
      ReaK::graph::make_sbastar_bundle( \
        motion_graph, start_node, this->m_space->get_super_space(), vis, \
        heuristic_map,  \
        pos_map,  \
        weight_map, \
        dens_map,  \
        cons_map,  \
        dist_map, \
        pred_map,  \
        get(&VertexProp::key_value, motion_graph),  \
        nc_selector), \
      get(random_sampler, this->m_space->get_super_space()), \
      this->m_init_relaxation, this->m_SA_init_temperature);
  
  
  
  if(m_knn_flag == LINEAR_SEARCH_KNN) {
    MEAQR_sbastar_visitor<StateSpace, StateSpaceSystem, StateSpaceSampler, MotionGraphType, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node, space_dim, space_Lc);
    
    ReaK::graph::fixed_neighborhood< linear_pred_succ_search<> > nc_selector(
      linear_pred_succ_search<>(), 
      10, max_radius);
//   ReaK::graph::star_neighborhood< linear_pred_succ_search<> > nc_selector(
//     linear_pred_succ_search<>(), 
//     space_dim, 3.0 * space_Lc);
    
    
    if(this->m_init_relaxation < 1e-6) {
      if(m_collision_check_flag == EAGER_COLLISION_CHECKING) {
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_SBASTAR_FUNCTION
        };
      } else { // assume lazy collision checking
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_LAZY_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_LAZY_SBASTAR_FUNCTION
        };
      };
    } else {
      if(m_collision_check_flag == EAGER_COLLISION_CHECKING) {
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_SBASTAR_FUNCTION
        };
      } else { // assume lazy collision checking
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBASTAR_FUNCTION
        };
      };
    };
    
  } else if(m_knn_flag == DVP_BF2_TREE_KNN) {
    
    typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 2, 
                     random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2>, 
                     no_position_caching_policy > SpacePartType;
    SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
    
    multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> nn_finder;
    nn_finder.graph_tree_map[&motion_graph] = &space_part;
    
    MEAQR_sbastar_visitor<StateSpace, StateSpaceSystem, StateSpaceSampler, MotionGraphType, multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node, space_dim, space_Lc);
    
    ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> > nc_selector(
      nn_finder, 
      10, max_radius);
//   ReaK::graph::star_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> > nc_selector(
//     nn_finder, 
//     space_dim, 3.0 * space_Lc);
    
    
    if(this->m_init_relaxation < 1e-6) {
      if(m_collision_check_flag == EAGER_COLLISION_CHECKING) {
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_SBASTAR_FUNCTION
        };
      } else { // assume lazy collision checking
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_LAZY_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_LAZY_SBASTAR_FUNCTION
        };
      };
    } else {
      if(m_collision_check_flag == EAGER_COLLISION_CHECKING) {
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_SBASTAR_FUNCTION
        };
      } else { // assume lazy collision checking
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBASTAR_FUNCTION
        };
      };
    };
    
  } else if(m_knn_flag == DVP_BF4_TREE_KNN) {
    
    typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 4, 
                      random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4>, 
                     no_position_caching_policy > SpacePartType;
    SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
    
    multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> nn_finder;
    nn_finder.graph_tree_map[&motion_graph] = &space_part;
    
    MEAQR_sbastar_visitor<StateSpace, StateSpaceSystem, StateSpaceSampler, MotionGraphType, multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node, space_dim, space_Lc);
    
    ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> > nc_selector(
      nn_finder, 
      10, max_radius);
//   ReaK::graph::star_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> > nc_selector(
//     nn_finder, 
//     space_dim, 3.0 * space_Lc);
    
    if(this->m_init_relaxation < 1e-6) {
      if(m_collision_check_flag == EAGER_COLLISION_CHECKING) {
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_SBASTAR_FUNCTION
        };
      } else { // assume lazy collision checking
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_LAZY_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_LAZY_SBASTAR_FUNCTION
        };
      };
    } else {
      if(m_collision_check_flag == EAGER_COLLISION_CHECKING) {
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_SBASTAR_FUNCTION
        };
      } else { // assume lazy collision checking
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBASTAR_FUNCTION
        };
      };
    };
    
  } else if(m_knn_flag == DVP_COB2_TREE_KNN) {
    
    typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 2, 
                      random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2>, 
                     no_position_caching_policy > SpacePartType;
    SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
    
    multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> nn_finder;
    nn_finder.graph_tree_map[&motion_graph] = &space_part;
    
    MEAQR_sbastar_visitor<StateSpace, StateSpaceSystem, StateSpaceSampler, MotionGraphType, multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node, space_dim, space_Lc);
    
    ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> > nc_selector(
      nn_finder, 
      10, max_radius);
//   ReaK::graph::star_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> > nc_selector(
//     nn_finder, 
//     space_dim, 3.0 * space_Lc);
    
    if(this->m_init_relaxation < 1e-6) {
      if(m_collision_check_flag == EAGER_COLLISION_CHECKING) {
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_SBASTAR_FUNCTION
        };
      } else { // assume lazy collision checking
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_LAZY_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_LAZY_SBASTAR_FUNCTION
        };
      };
    } else {
      if(m_collision_check_flag == EAGER_COLLISION_CHECKING) {
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_SBASTAR_FUNCTION
        };
      } else { // assume lazy collision checking
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBASTAR_FUNCTION
        };
      };
    };
    
  } else if(m_knn_flag == DVP_COB4_TREE_KNN) {
    
    typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 4, 
                      random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4>, 
                     no_position_caching_policy > SpacePartType;
    SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
    
    multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> nn_finder;
    nn_finder.graph_tree_map[&motion_graph] = &space_part;
    
    MEAQR_sbastar_visitor<StateSpace, StateSpaceSystem, StateSpaceSampler, MotionGraphType, multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node, space_dim, space_Lc);
    
    ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> > nc_selector(
      nn_finder, 
      10, max_radius);
//   ReaK::graph::star_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> > nc_selector(
//     nn_finder, 
//     space_dim, 3.0 * space_Lc);
    
    if(this->m_init_relaxation < 1e-6) {
      if(m_collision_check_flag == EAGER_COLLISION_CHECKING) {
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_SBASTAR_FUNCTION
        };
      } else { // assume lazy collision checking
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_LAZY_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_LAZY_SBASTAR_FUNCTION
        };
      };
    } else {
      if(m_collision_check_flag == EAGER_COLLISION_CHECKING) {
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_SBASTAR_FUNCTION
        };
      } else { // assume lazy collision checking
        if(m_added_bias_flags == PLAN_WITH_VORONOI_PULL) {
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBARRTSTAR_FUNCTION
        } else { // assume nominal method only.
          RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBASTAR_FUNCTION
        };
      };
    };
    
  };
  
#undef RK_MEAQR_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
#undef RK_MEAQR_SBASTAR_PLANNER_CALL_SBARRTSTAR_FUNCTION
#undef RK_MEAQR_SBASTAR_PLANNER_CALL_SBASTAR_FUNCTION
#undef RK_MEAQR_SBASTAR_PLANNER_CALL_LAZY_SBARRTSTAR_FUNCTION
#undef RK_MEAQR_SBASTAR_PLANNER_CALL_LAZY_SBASTAR_FUNCTION
#undef RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_SBARRTSTAR_FUNCTION
#undef RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_SBASTAR_FUNCTION
#undef RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBARRTSTAR_FUNCTION
#undef RK_MEAQR_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBASTAR_FUNCTION
  
  if(m_solutions.size())
    return m_solutions.begin()->second;
  else
    return shared_ptr< seq_path_base< SuperSpace > >();
};


};

};

#endif


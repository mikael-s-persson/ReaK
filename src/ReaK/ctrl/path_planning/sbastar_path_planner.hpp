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
#include "sbmp_reporter_concept.hpp"

#include "metric_space_concept.hpp"
#include "seq_path_wrapper.hpp"
#include "interpolation/point_to_point_path.hpp"
#include "basic_sbmp_reporters.hpp"

#include "graph_alg/lazy_sbastar.hpp"
#include "graph_alg/sbastar_rrtstar.hpp"
#include "graph_alg/anytime_sbastar.hpp"

#include "graph_alg/d_ary_bf_tree.hpp"
#include "graph_alg/d_ary_cob_tree.hpp"
#include "graph_alg/bgl_tree_adaptor.hpp"
#include "graph_alg/bgl_more_property_maps.hpp"
#include "graph_alg/pooled_adjacency_list.hpp"
#include "graph_alg/bgl_tree_adaptor.hpp"
#include "dvp_layout_adjacency_list.hpp"
#include "metric_space_search.hpp"
#include "topological_search.hpp"
#include "path_planner_options.hpp"
#include "graph_alg/neighborhood_functors.hpp"
#include "lin_alg/arithmetic_tuple.hpp"

#include "base/misc_math.hpp"

#include <stack>

namespace ReaK {
  
namespace pp {
  
  
  
/**
 * This POD type contains the data required on a per-vertex basis for the 
 * Sampling-based A* path-planning algorithms.
 * \tparam FreeSpaceType The topology type on which to perform the planning, should be the C-free sub-space of a larger configuration space.
 */
template <typename FreeSpaceType>
struct sbastar_vertex_data {
  /// The position associated to the vertex.
  typename topology_traits<FreeSpaceType>::point_type position;
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
  sbastar_vertex_data() : position(typename topology_traits<FreeSpaceType>::point_type()),
                          constriction(0.0), collision_count(0), density(0.0), expansion_trials(0),
                          heuristic_value(0.0), distance_accum(0.0), key_value(0.0),
                          predecessor(0) { };
};

/**
 * This POD type contains the data required on a per-edge basis for the 
 * Sampling-based A* path-planning algorithms.
 * \tparam FreeSpaceType The topology type on which to perform the planning, should be the C-free sub-space of a larger configuration space.
 */
template <typename FreeSpaceType>
struct sbastar_edge_data { 
  /// The travel-distance associated to the edge (from source to target).
  double astar_weight;
  
  /**
   * Default constructor.
   * \param aWeight The travel-distance to be associated to this edge.
   */
  sbastar_edge_data(double aWeight = 0.0) : astar_weight(aWeight) { };
};



/**
 * This stateless functor type can be used to print out the information about a given SBA* vertex.
 * This functor type can be used as a printing policy type for the vlist_sbmp_report class 
 * template that prints the list of vertices to a file.
 * \note This is mostly useful for debugging purposes (recording all information about the 
 *       motion-graph), it should not be used as the "output" of the path-planner.
 */
struct sbastar_vprinter : serialization::serializable {
  
  /**
   * This call operator prints all the SBA* information about a given vertex 
   * to a given output-stream.
   * \tparam Vertex The vertex-descriptor type for the motion-graph.
   * \tparam Graph The motion-graph type used by the SBA* planning algorithm.
   * \param out The output-stream to which to print the SBA* information about the vertex.
   * \param u The vertex whose information is to be printed.
   * \param g The motion-graph to which the vertex belongs.
   */
  template <typename Vertex, typename Graph>
  void operator()(std::ostream& out, Vertex u, const Graph& g) const {
    using ReaK::to_vect;
    vect_n<double> v_pos = to_vect<double>(g[u].position);
    for(std::size_t i = 0; i < v_pos.size(); ++i)
      out << " " << std::setw(10) << v_pos[i];
    out << " " << std::setw(10) << g[u].constriction 
        << " " << std::setw(10) << g[u].collision_count 
        << " " << std::setw(10) << g[u].density 
        << " " << std::setw(10) << g[u].expansion_trials 
        << " " << std::setw(10) << g[u].heuristic_value 
        << " " << std::setw(10) << g[u].distance_accum 
        << " " << std::setw(10) << g[u].key_value << std::endl;
  };
  
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const { };
  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { };
  
  RK_RTTI_MAKE_ABSTRACT_1BASE(sbastar_vprinter,0xC2460012,1,"sbastar_vprinter",serialization::serializable)
};





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
 * \tparam SBPPReporter The reporter type to use to report the progress of the path-planning.
 */
template <typename FreeSpaceType, 
          typename SBPPReporter = no_sbmp_report>
class sbastar_path_planner : public sample_based_planner< path_planner_base<FreeSpaceType> > {
  public:
    typedef sample_based_planner< path_planner_base<FreeSpaceType> > base_type;
    typedef sbastar_path_planner<FreeSpaceType, SBPPReporter> self;
    
    typedef FreeSpaceType space_type;
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    
    BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));
    
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
    
    double m_current_key_threshold;
    double m_current_dens_threshold;
    double m_current_relaxation;
    
    std::map<double, shared_ptr< seq_path_base< super_space_type > > > m_solutions;
    
    std::size_t m_iteration_count;
    
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
    template <typename Graph>
    void report_progress(Graph& g) {
      m_iteration_count++;
      if(m_iteration_count % this->m_progress_interval == 0)
        m_reporter.draw_motion_graph(*(this->m_space), g, get(&sbastar_vertex_data<FreeSpaceType>::position,g));
      has_reached_max_vertices = (m_iteration_count >= this->m_max_vertex_count);
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
        shared_ptr< seq_path_wrapper< point_to_point_path<super_space_type> > > new_sol(new seq_path_wrapper< point_to_point_path<super_space_type> >("sbastar_solution", point_to_point_path<super_space_type>(sup_space_ptr,get(distance_metric, this->m_space->get_super_space()))));
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
     * Returns the sampling radius (in the topology's distance metric) used by this planner when doing random walks.
     * \return The sampling radius used by this planner.
     */
    double get_sampling_radius() const { return m_sampling_radius; };
    /**
     * Sets the sampling radius (in the topology's distance metric) to be used by this planner when doing random walks.
     * \param aSamplingRadius The sampling radius to be used by this planner when doing random walks.
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
     *                             The options available include EAGER_COLLISION_CHECKING or LAZY_COLLISION_CHECKING,
     *                             NOMINAL_PLANNER_ONLY or any combination of PLAN_WITH_VORONOI_PULL and PLAN_WITH_ANYTIME_HEURISTIC,
     *                             and USE_BRANCH_AND_BOUND_PRUNING_FLAG. See path_planner_options.hpp documentation.
     * \param aReporter The SBPP reporter object to use to report results and progress.
     * \param aMaxResultCount The maximum number of successful start-goal connections to make before 
     *                        stopping the path planner (the higher the number the more likely that a 
     *                        good path will be found, however, running time can become much longer).
     */
    sbastar_path_planner(const shared_ptr< space_type >& aWorld = shared_ptr< space_type >(), 
                         const point_type& aStartPos = point_type(),
                         const point_type& aGoalPos = point_type(),
                         std::size_t aMaxVertexCount = 5000, 
                         std::size_t aProgressInterval = 100,
                         std::size_t aDataStructureFlags = ADJ_LIST_MOTION_GRAPH | DVP_BF2_TREE_KNN,
                         std::size_t aPlanningMethodFlags = LAZY_COLLISION_CHECKING | NOMINAL_PLANNER_ONLY,
                         SBPPReporter aReporter = SBPPReporter(),
                         std::size_t aMaxResultCount = 50) :
                         base_type("sbastar_planner", aWorld, aMaxVertexCount, aProgressInterval, aDataStructureFlags, aPlanningMethodFlags),
                         m_reporter(aReporter),
                         m_start_pos(aStartPos),
                         m_goal_pos(aGoalPos),
                         m_init_key_threshold(0.8),
                         m_init_dens_threshold(0.8),
                         m_init_relaxation(0.0),
                         m_SA_init_temperature(-1.0),
                         m_sampling_radius(1.0),
                         m_current_num_results(0),
                         max_num_results(aMaxResultCount),
                         has_reached_max_vertices(false),
                         m_current_key_threshold(m_init_key_threshold),
                         m_current_dens_threshold(m_init_dens_threshold),
                         m_current_relaxation(m_init_relaxation),
                         m_solutions(),
                         m_iteration_count(0) { };
    
    virtual ~sbastar_path_planner() { };
    
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
        & RK_SERIAL_SAVE_WITH_NAME(max_num_results);
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
        & RK_SERIAL_LOAD_WITH_NAME(max_num_results);
      has_reached_max_vertices = false;
      m_solutions.clear();
      m_current_num_results = 0;
      m_current_key_threshold = m_init_key_threshold;
      m_current_dens_threshold = m_init_dens_threshold;
      m_current_relaxation = m_init_relaxation;
      m_iteration_count = 0;
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC246000C,1,"sbastar_path_planner",base_type)
};



/**
 * This class template is used by the SBA* path-planner as the visitor object needed to 
 * collaborate with the SBA* algorithms to generate the motion-graph and path-planning solutions.
 * This class template models the SBAStarVisitorConcept and SBARRTStarVisitorConcept.
 * As with most planning algorithms in ReaK, the algorithm is really made up of a high-level 
 * algorithmic logic in the form of function templates, and a number of customization points 
 * collected as member functions of an algorithm visitor class that implement the problem-specific 
 * behaviors (random-walks / local-planning, heuristic computation, progress reporting, 
 * completion criteria, etc.). 
 */
template <typename FreeSpaceType, typename MotionGraph, typename NNFinderSynchro, typename SBPPReporter = no_sbmp_report>
struct sbastar_planner_visitor {
  typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
  
  shared_ptr< FreeSpaceType > m_space;
  sbastar_path_planner<FreeSpaceType,SBPPReporter>* m_planner;
  NNFinderSynchro m_nn_synchro;
  Vertex m_start_node;
  Vertex m_goal_node;
  double m_space_dim;
  double m_space_Lc;
  
  sbastar_planner_visitor(const shared_ptr< FreeSpaceType >& aSpace, 
                          sbastar_path_planner<FreeSpaceType,SBPPReporter>* aPlanner,
                          NNFinderSynchro aNNSynchro,
                          Vertex aStartNode, Vertex aGoalNode, double aSpaceDim, double aSpaceLc) : 
                          m_space(aSpace), m_planner(aPlanner), m_nn_synchro(aNNSynchro),
                          m_start_node(aStartNode), m_goal_node(aGoalNode), 
                          m_space_dim(aSpaceDim), m_space_Lc(aSpaceLc) { };
  
  typedef typename topology_traits<FreeSpaceType>::point_type PointType;
  typedef sbastar_edge_data<FreeSpaceType> EdgeProp;
  
  
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
    
//     if(travel_dist > m_planner->get_sampling_radius())
//       return 0.0;
    double sig2_n = 0.25 * m_planner->get_sampling_radius() * m_planner->get_sampling_radius();
    double sig2_x = 0.25 * event_radius * event_radius;
    return exp(-travel_dist * travel_dist / (sig2_x * 2.0) - 0.5 * m_space_dim * ( sig2_n / sig2_x - 1.0 - log(sig2_n / sig2_x) ) );
  };
  
  double compute_sample_similarity(double travel_dist) const {
    using std::exp;
    
//     if(travel_dist > m_planner->get_sampling_radius())
//       return 0.0;
    return exp(-travel_dist * travel_dist / (0.25 * m_planner->get_sampling_radius() * m_planner->get_sampling_radius() * 2.0));
  };
  
#if 0
  // keep the average sample similarity:
  
  template <typename Vertex, typename Graph>
  void register_explored_sample(Vertex u, Graph& g, double samp_sim) const {
    g[u].density = (g[u].expansion_trials * g[u].density + samp_sim) / (g[u].expansion_trials + 1);
    ++(g[u].expansion_trials);
//     std::cout << " Vertex " << u << " has density: " << g[u].density << std::endl;
  };
  
  template <typename Vertex, typename Graph>
  void register_failed_sample(Vertex u, Graph& g, double samp_sim) const {
    g[u].constriction = (g[u].collision_count * g[u].constriction + samp_sim) / (g[u].collision_count + 1);
    ++(g[u].collision_count);
//     std::cout << " Vertex " << u << " has constriction: " << g[u].constriction << std::endl;
  };
#endif
  
  // NOTE: This seems to be the best one.
#if 0
  // keep the average sample similarity weighted by the sample probability (and its binomial converse).
  
  template <typename Vertex, typename Graph>
  void register_explored_sample(Vertex u, Graph& g, double samp_sim) const {
    double tmp_density = g[u].density * (1.0 - samp_sim) + samp_sim * samp_sim;
    g[u].density = (g[u].expansion_trials * g[u].density + tmp_density) / (g[u].expansion_trials + 1);
    ++(g[u].expansion_trials);
//     std::cout << " Vertex " << u << " has density: " << g[u].density << std::endl;
  };
  template <typename Vertex, typename Graph>
  void register_failed_sample(Vertex u, Graph& g, double samp_sim) const {
    double tmp_density = g[u].constriction * (1.0 - samp_sim) + samp_sim * samp_sim;
    g[u].constriction = (g[u].collision_count * g[u].constriction + tmp_density) / (g[u].collision_count + 1);
    ++(g[u].collision_count);
//     std::cout << " Vertex " << u << " has constriction: " << g[u].constriction << std::endl;
  };
#endif
  
  // NOTE: This one seems to work quite well too.
#if 1
  // keep the sample similarity weighted by the sample probability (and its binomial converse).
  // that is, assume the existing density to reflect the overall density and the newly computed 
  // sample similarity to reflect the density in its relatively probable region (binomial).
  
  template <typename Vertex, typename Graph>
  void register_explored_sample(Vertex u, Graph& g, double samp_sim) const {
    g[u].density = g[u].density * (1.0 - samp_sim) + samp_sim * samp_sim;
//     std::cout << " Vertex " << u << " has density: " << g[u].density << std::endl;
  };
  template <typename Vertex, typename Graph>
  void register_failed_sample(Vertex u, Graph& g, double samp_sim) const {
    g[u].constriction = g[u].constriction * (1.0 - samp_sim) + samp_sim * samp_sim;
//     std::cout << " Vertex " << u << " has constriction: " << g[u].constriction << std::endl;
  };
#endif
  
#if 0
  // keep track only of the maximum sample similarity in the region of a node:
  
  template <typename Vertex, typename Graph>
  void register_explored_sample(Vertex u, Graph& g, double samp_sim) const {
    if(samp_sim > g[u].density)
      g[u].density = samp_sim;
//     std::cout << " Vertex " << u << " has density: " << g[u].density << std::endl;
  };
  template <typename Vertex, typename Graph>
  void register_failed_sample(Vertex u, Graph& g, double samp_sim) const {
    if(samp_sim > g[u].constriction)
      g[u].constriction = samp_sim;
//     std::cout << " Vertex " << u << " has constriction: " << g[u].constriction << std::endl;
  };
#endif
  
  
  
  template <typename Vertex, typename Graph>
  void vertex_added(Vertex u, Graph& g) const {
    m_nn_synchro.added_vertex(u,g);
    
    init_vertex_properties(u,g);
    
    // Call progress reporter...
    m_planner->report_progress(g);
    
    if((in_degree(m_goal_node,g)) && (g[m_goal_node].distance_accum < m_planner->get_best_solution_distance()))
      m_planner->create_solution_path(m_start_node, m_goal_node, g);
  };
  
  template <typename Vertex, typename Graph>
  void vertex_to_be_removed(Vertex u, Graph& g) const {
    m_nn_synchro.removed_vertex(u,g);
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
  
  template <typename Vertex, typename Graph>
  void affected_vertex(Vertex, Graph&) const { 
    // nothing to do.
  };
  
      
  
  bool keep_going() const {
    if( ( m_planner->get_initial_relaxation() > 1e-6 ) && ( m_planner->get_current_relaxation() < 1e-6 ) )
      return false;
    return m_planner->keep_going();
  };
  
  template <typename Vertex, typename Graph>
  boost::tuple<PointType, bool, EdgeProp > random_walk(Vertex u, Graph& g) const {
    typedef boost::tuple<PointType, bool, EdgeProp > ResultType;
    using std::exp;
    using std::log;
    using std::fabs;
    
    typename point_distribution_traits< typename subspace_traits<FreeSpaceType>::super_space_type >::random_sampler_type get_sample = get(random_sampler, m_space->get_super_space());
    typename metric_space_traits< typename subspace_traits<FreeSpaceType>::super_space_type >::distance_metric_type get_distance = get(distance_metric, m_space->get_super_space());
    
    boost::variate_generator< pp::global_rng_type&, boost::normal_distribution<double> > var_rnd(pp::get_global_rng(), boost::normal_distribution<double>());
    
    unsigned int i = 0;
//     PointType p_rnd = g[m_goal_node].position;
    PointType p_rnd = get_sample(m_space->get_super_space());
    do {
//       PointType p_rnd = get_sample(m_space->get_super_space());
      double dist = get_distance(g[u].position, p_rnd, m_space->get_super_space());
      double target_dist = boost::uniform_01<global_rng_type&,double>(get_global_rng())() * m_planner->get_sampling_radius();
//       double target_dist = fabs(var_rnd()) * m_planner->get_sampling_radius();
      PointType p_v = m_space->move_position_toward(g[u].position, target_dist / dist, p_rnd);
      dist = get_distance(g[u].position, p_v, m_space->get_super_space());
      if( dist < 0.9 * target_dist ) {
//       if(( dist < 0.9 * target_dist ) && ( dist < 0.95 * m_planner->get_sampling_radius() )) {
        // this means that we had a collision before reaching the target distance, 
        // must record that to the constriction statistic:
        double samp_sim = compute_sample_similarity(target_dist, (target_dist - dist));
        if(samp_sim > std::numeric_limits<double>::epsilon())
          register_failed_sample(u,g,samp_sim);
        
        p_rnd = get_sample(m_space->get_super_space());
      } else
        return ResultType(p_v, true, EdgeProp(dist));
    } while(++i <= 10);
    return ResultType(g[u].position, false, EdgeProp(std::numeric_limits<double>::infinity()));
  };
  
  // for the SBA*-RRT* variant:
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
  template <class Edge, class Graph>
  void examine_edge(Edge e, Graph& g) const { };
  template <typename Edge, typename Graph>
  void edge_relaxed(Edge, const Graph&) const { };
  
  template <typename Graph>
  void publish_path(const Graph& g) const {
    // try to create a goal connection path
    if((in_degree(m_goal_node,g)) && (g[m_goal_node].distance_accum < m_planner->get_best_solution_distance()))
      m_planner->create_solution_path(m_start_node, m_goal_node, g); 
    
//     m_planner->set_current_key_threshold( 0.5 * m_planner->get_current_key_threshold() );
//     m_planner->set_current_density_threshold( 0.95 * m_planner->get_current_density_threshold() );
    
//     std::cout << " new key-value threshold =\t" << m_planner->get_current_key_threshold() << std::endl;
//     std::cout << " new density-value threshold =\t" << m_planner->get_current_density_threshold() << std::endl;
  };
  
  template <typename Vertex, typename Graph>
  bool has_search_potential(Vertex u, const Graph& g) const { 
    if( m_planner->get_initial_relaxation() > 1e-6 ) {
      // assume we are running a Anytime SBA* algorithm.
//       std::cout << " A* key = " << std::setw(10) << (g[u].distance_accum + g[u].heuristic_value) 
//                 << " weighted by density = " << std::setw(10) << ((g[u].distance_accum + g[u].heuristic_value) / (1.0 - g[u].constriction) / (1.0 - g[u].density))
//                 << " compared to key-value = " << std::setw(10) << (g[u].key_value)
//                 << " and compared to relaxed heuristic: " << std::setw(10) << (g[u].heuristic_value * m_planner->get_current_relaxation()) << std::endl;
//       std::cout << " sampling potential = " << std::setw(10) << ((1.0 - g[u].constriction) * (1.0 - g[u].density)) 
//                 << " compared to relaxation potential: " << std::setw(10) << (1.0 / (1.0 + m_planner->get_current_relaxation())) << std::endl;
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






template <typename FreeSpaceType, 
          typename SBPPReporter>
shared_ptr< seq_path_base< typename sbastar_path_planner<FreeSpaceType,SBPPReporter>::super_space_type > > 
  sbastar_path_planner<FreeSpaceType,SBPPReporter>::solve_path() {
  using ReaK::to_vect;
  
  this->has_reached_max_vertices = false;
  this->m_current_num_results = 0;
  this->m_solutions.clear();
  this->m_current_key_threshold = this->m_init_key_threshold;
  this->m_current_dens_threshold = this->m_init_dens_threshold;
  this->m_current_relaxation = this->m_init_relaxation;
  this->m_iteration_count = 0;
  
  typedef typename subspace_traits<FreeSpaceType>::super_space_type SuperSpace;
  typedef typename topology_traits<SuperSpace>::point_type PointType;
  typedef sbastar_vertex_data<FreeSpaceType> VertexProp;
  typedef sbastar_edge_data<FreeSpaceType> EdgeProp;
  
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
  
  double space_dim = double((to_vect<double>(this->m_space->get_super_space().difference(this->m_goal_pos,this->m_start_pos))).size()); 
  double space_Lc = get(distance_metric,this->m_space->get_super_space())(this->m_start_pos, this->m_goal_pos, this->m_space->get_super_space());
  
  
  // Some MACROs to reduce the size of the code below.
  
#ifdef RK_ENABLE_CXX0X_FEATURES

#define RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE \
    VertexProp vs_p, vg_p; \
    vs_p.position = this->m_start_pos; \
    vg_p.position = this->m_goal_pos; \
    Vertex start_node = add_vertex(std::move(vs_p), motion_graph); \
    Vertex goal_node = add_vertex(std::move(vg_p), motion_graph); \
    motion_graph[start_node].constriction = 0.0; \
    motion_graph[start_node].collision_count = 0; \
    motion_graph[start_node].density = 0.0; \
    motion_graph[start_node].expansion_trials = 0; \
    motion_graph[start_node].heuristic_value = heuristic(start_node,motion_graph); \
    motion_graph[start_node].distance_accum = 0.0; \
    motion_graph[start_node].key_value = space_Lc + this->m_init_relaxation * space_Lc; \
    motion_graph[start_node].predecessor = start_node; \
     \
    motion_graph[goal_node].constriction = 0.0; \
    motion_graph[goal_node].collision_count = 0; \
    motion_graph[goal_node].density = 0.0; \
    motion_graph[goal_node].expansion_trials = 0; \
    motion_graph[goal_node].heuristic_value = 0.0; \
    motion_graph[goal_node].distance_accum = std::numeric_limits<double>::infinity(); \
    motion_graph[goal_node].key_value = std::numeric_limits<double>::infinity(); \
    motion_graph[goal_node].predecessor = goal_node;
    
#else
    
#define RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE \
    VertexProp vs_p, vg_p; \
    vs_p.position = this->m_start_pos; \
    vg_p.position = this->m_goal_pos; \
    Vertex start_node = add_vertex(vs_p, motion_graph); \
    Vertex goal_node = add_vertex(vg_p, motion_graph); \
    motion_graph[start_node].constriction = 0.0; \
    motion_graph[start_node].collision_count = 0; \
    motion_graph[start_node].density = 0.0; \
    motion_graph[start_node].expansion_trials = 0; \
    motion_graph[start_node].heuristic_value = heuristic(start_node,motion_graph); \
    motion_graph[start_node].distance_accum = 0.0; \
    motion_graph[start_node].key_value = space_Lc + this->m_init_relaxation * space_Lc; \
    motion_graph[start_node].predecessor = start_node; \
     \
    motion_graph[goal_node].constriction = 0.0; \
    motion_graph[goal_node].collision_count = 0; \
    motion_graph[goal_node].density = 0.0; \
    motion_graph[goal_node].expansion_trials = 0; \
    motion_graph[goal_node].heuristic_value = 0.0; \
    motion_graph[goal_node].distance_accum = std::numeric_limits<double>::infinity(); \
    motion_graph[goal_node].key_value = std::numeric_limits<double>::infinity(); \
    motion_graph[goal_node].predecessor = goal_node;
  
#endif
    
    
#define RK_SBASTAR_PLANNER_CALL_SBASTAR_FUNCTION \
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
  
  
#define RK_SBASTAR_PLANNER_CALL_LAZY_SBASTAR_FUNCTION \
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
  
#define RK_SBASTAR_PLANNER_CALL_LAZY_BNB_SBASTAR_FUNCTION \
    ReaK::graph::generate_lazy_bnb_sbastar( \
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
      goal_node \
      );
  
#define RK_SBASTAR_PLANNER_CALL_SBARRTSTAR_FUNCTION \
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
  
  
#define RK_SBASTAR_PLANNER_CALL_LAZY_SBARRTSTAR_FUNCTION \
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
  
  
#define RK_SBASTAR_PLANNER_CALL_LAZY_BNB_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_lazy_bnb_sbarrtstar( \
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
      goal_node, \
      get(random_sampler, this->m_space->get_super_space()), \
      this->m_SA_init_temperature);
   
  
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_SBASTAR_FUNCTION \
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
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBASTAR_FUNCTION \
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
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_BNB_SBASTAR_FUNCTION \
    ReaK::graph::generate_anytime_lazy_bnb_sbastar( \
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
      goal_node, \
      this->m_init_relaxation);
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_SBARRTSTAR_FUNCTION \
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
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBARRTSTAR_FUNCTION \
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
  
  
#define RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_BNB_SBARRTSTAR_FUNCTION \
    ReaK::graph::generate_anytime_lazy_bnb_sbarrtstar( \
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
      goal_node, \
      get(random_sampler, this->m_space->get_super_space()), \
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
    typedef boost::composite_property_map< 
      PositionMap, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t > > GraphPositionMap;
    
    MotionGraphType motion_graph;
    GraphPositionMap g_pos_map = GraphPositionMap(pos_map, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t >(&motion_graph));
    
    RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
    
    if((this->m_data_structure_flags & KNN_METHOD_MASK) == LINEAR_SEARCH_KNN) {
      sbastar_planner_visitor<FreeSpaceType, MotionGraphType, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node, space_dim, space_Lc);
      
//       ReaK::graph::fixed_neighborhood< linear_neighbor_search<> > nc_selector(
//         linear_neighbor_search<>(), 
//         10, m_sampling_radius);
      ReaK::graph::star_neighborhood< linear_neighbor_search<> > nc_selector(
        linear_neighbor_search<>(), 
        space_dim, 3.0 * space_Lc);
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
      
      typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 2, 
                       random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2> > SpacePartType;
      SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
      
      multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraphType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node, space_dim, space_Lc);
      
//       ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> > nc_selector(
//         nn_finder, 
//         10, max_radius);
      
      ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> > nc_selector(
        nn_finder, 
        space_dim, 3.0 * space_Lc);
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
      
      typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 4, 
                       random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4> > SpacePartType;
      SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
      
      multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraphType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node, space_dim, space_Lc);
      
//       ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> > nc_selector(
//         nn_finder, 
//         10, m_sampling_radius);
      
      ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> > nc_selector(
        nn_finder, 
        space_dim, 3.0 * space_Lc);
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
      
      typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 2, 
                       random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2> > SpacePartType;
      SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
      
      multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraphType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node, space_dim, space_Lc);
      
//       ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> > nc_selector(
//         nn_finder, 
//         10, m_sampling_radius);
      
      ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> > nc_selector(
        nn_finder, 
        space_dim, 3.0 * space_Lc);
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
      
      typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 4, 
                       random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4> > SpacePartType;
      SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
      
      multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraphType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node, space_dim, space_Lc);
      
//       ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> > nc_selector(
//         nn_finder, 
//         10, m_sampling_radius);
      
      ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> > nc_selector(
        nn_finder, 
        space_dim, 3.0 * space_Lc);
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    };
    
  } else if((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) == DVP_ADJ_LIST_MOTION_GRAPH) {
    
    if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {
      
      typedef dvp_adjacency_list<
        VertexProp,
        EdgeProp,
        SuperSpace,
        PositionMap,
        2, random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2>,
        boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
      
      ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
      
      typedef typename ALTGraph::adj_list_type MotionGraph;
      typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
      
      MotionGraph motion_graph = space_part.get_adjacency_list();
      
      RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
      
      multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node, space_dim, space_Lc);
      
//       ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> > nc_selector(
//         nn_finder, 
//         10, m_sampling_radius);
      
      ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> > nc_selector(
        nn_finder, 
        space_dim, 3.0 * space_Lc);
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF4_TREE_KNN) {
      
      typedef dvp_adjacency_list<
        VertexProp,
        EdgeProp,
        SuperSpace,
        PositionMap,
        4, random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4>,
        boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
      
      ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
      
      typedef typename ALTGraph::adj_list_type MotionGraph;
      typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
      
      MotionGraph motion_graph = space_part.get_adjacency_list();
      
      RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
      
      multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node, space_dim, space_Lc);
      
//       ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> > nc_selector(
//         nn_finder, 
//         10, m_sampling_radius);
      
      ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> > nc_selector(
        nn_finder, 
        space_dim, 3.0 * space_Lc);
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB2_TREE_KNN) {
      
      typedef dvp_adjacency_list<
        VertexProp,
        EdgeProp,
        SuperSpace,
        PositionMap,
        2, random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2>,
        boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
      
      ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
      
      typedef typename ALTGraph::adj_list_type MotionGraph;
      typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
      
      MotionGraph motion_graph = space_part.get_adjacency_list();
      
      RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
      
      multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node, space_dim, space_Lc);
      
//       ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> > nc_selector(
//         nn_finder, 
//         10, m_sampling_radius);
      
      ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> > nc_selector(
        nn_finder, 
        space_dim, 3.0 * space_Lc);
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    } else if((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_COB4_TREE_KNN) {
      
      typedef dvp_adjacency_list<
        VertexProp,
        EdgeProp,
        SuperSpace,
        PositionMap,
        4, random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4>,
        boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
      
      ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
      
      typedef typename ALTGraph::adj_list_type MotionGraph;
      typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
      
      MotionGraph motion_graph = space_part.get_adjacency_list();
      
      RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
      
      multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node, space_dim, space_Lc);
      
//       ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> > nc_selector(
//         nn_finder, 
//         10, m_sampling_radius);
      
      ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> > nc_selector(
        nn_finder, 
        space_dim, 3.0 * space_Lc);
      
      RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
      
    };
    
  };
  
#undef RK_SBASTAR_PLANNER_INIT_START_AND_GOAL_NODE
#undef RK_SBASTAR_PLANNER_CALL_SBARRTSTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_SBASTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_LAZY_SBARRTSTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_LAZY_SBASTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_ANYTIME_SBARRTSTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_ANYTIME_SBASTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBARRTSTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_ANYTIME_LAZY_SBASTAR_FUNCTION
#undef RK_SBASTAR_PLANNER_CALL_APPROPRIATE_SBASTAR_PLANNER_FUNCTION
  
  if(m_solutions.size())
    return m_solutions.begin()->second;
  else
    return shared_ptr< seq_path_base< SuperSpace > >();
};


};

};

#endif


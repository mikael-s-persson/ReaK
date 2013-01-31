/**
 * \file sbastar_path_planner.hpp
 * 
 * This library defines a class
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2013
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

#ifndef REAK_SBASTAR_PATH_PLANNER_HPP
#define REAK_SBASTAR_PATH_PLANNER_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "motion_planner_base.hpp"
#include "sbmp_reporter_concept.hpp"

#include "metric_space_concept.hpp"
#include "path_base.hpp"
#include "interpolation/point_to_point_path.hpp"
#include "basic_sbmp_reporters.hpp"

#include "graph_alg/sbastar_search.hpp"

#include "graph_alg/d_ary_bf_tree.hpp"
#include "graph_alg/d_ary_cob_tree.hpp"
#include "graph_alg/bgl_tree_adaptor.hpp"
#include "graph_alg/bgl_more_property_maps.hpp"
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
  
  
  

template <typename FreeSpaceType>
struct sbastar_vertex_data {
  typename topology_traits<FreeSpaceType>::point_type position;
  double constriction;
  std::size_t collision_count;  // r
  double density;
  std::size_t expansion_trials;  // m
  double heuristic_value;
  double distance_accum;
  double key_value;
  boost::default_color_type astar_color;
  std::size_t predecessor;
  
  sbastar_vertex_data() : position(typename topology_traits<FreeSpaceType>::point_type()),
                          constriction(0.2), collision_count(1), density(0.0), expansion_trials(0),
                          heuristic_value(0.0), distance_accum(0.0), key_value(0.0),
                          astar_color(), predecessor(0) { };
};

template <typename FreeSpaceType>
struct sbastar_edge_data { 
  double astar_weight; //for A*
  
  sbastar_edge_data() : astar_weight(0.0) { };
};





/**
 * This class is a FADPRM-based path-planner over the given topology.
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
    double m_initial_threshold;
    double m_sampling_radius;
    std::size_t max_num_results;
    bool has_reached_max_vertices;
    std::size_t m_graph_kind_flag;
    std::size_t m_knn_flag;
    
    std::map<double, shared_ptr< path_base< super_space_type > > > m_solutions;
    
  public:
    
    bool keep_going() const {
      return (max_num_results > m_solutions.size()) && !has_reached_max_vertices;
    };
    
    template <typename Graph>
    double adjust_threshold(double old_thr, const Graph&) const { 
      return old_thr * 2.0;  // geometrically progress towards infinity.
    };
    
    template <typename Graph>
    void report_progress(Graph& g) {
      if(num_vertices(g) % this->m_progress_interval == 0)
        m_reporter.draw_motion_graph(*(this->m_space), g, get(&sbastar_vertex_data<FreeSpaceType>::position,g));
      has_reached_max_vertices = (num_vertices(g) >= this->m_max_vertex_count);
    };
    
    template <typename Graph>
    double heuristic(typename boost::graph_traits<Graph>::vertex_descriptor u, const Graph& g) const {
      return get(distance_metric, this->m_space->get_super_space())(g[u].position, this->m_goal_pos, this->m_space->get_super_space());
    };
    
    template <typename Vertex, typename Graph>
    void create_solution_path(Vertex start_node, Vertex goal_node, Graph& g) {
      
      double goal_distance = g[goal_node].distance_accum;
      
      if(goal_distance < std::numeric_limits<double>::infinity()) {
        //Draw the edges of the current best solution:
        
        shared_ptr< super_space_type > sup_space_ptr(&(this->m_space->get_super_space()),null_deleter());
        shared_ptr< path_wrapper< point_to_point_path<super_space_type> > > new_sol(new path_wrapper< point_to_point_path<super_space_type> >("sbastar_solution", point_to_point_path<super_space_type>(sup_space_ptr,get(distance_metric, this->m_space->get_super_space()))));
        point_to_point_path<super_space_type>& waypoints = new_sol->get_underlying_path();
        std::set<Vertex> path;
        
        Vertex v = goal_node;
        point_type p_v = this->m_goal_pos;
        Vertex u = g[v].predecessor;
        point_type p_u = g[u].position;
        
        waypoints.push_front(m_goal_pos);
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
     * This function computes a valid path in the C-free. If it cannot 
     * achieve a valid path, an exception will be thrown. This algorithmic
     * path solver class is such that any settings that ought to be set for the 
     * path planning algorithm should be set before calling this function, otherwise
     * the function is likely to fail.
     * \return The path object that can be used to map out the path.
     */
    virtual shared_ptr< path_base< super_space_type > > solve_path();
    
    const SBPPReporter& get_reporter() const { return m_reporter; };
    void set_reporter(const SBPPReporter& aNewReporter) { m_reporter = aNewReporter; };
    
    const point_type& get_start_pos() const { return m_start_pos; };
    void set_start_pos(const point_type& aStartPos) { m_start_pos = aStartPos; };
    
    const point_type& get_goal_pos() const { return m_goal_pos; };
    void set_goal_pos(const point_type& aGoalPos) { m_goal_pos = aGoalPos; };
    
    double get_initial_threshold() const { return m_initial_threshold; };
    void set_initial_threshold(double aInitialThreshold) { m_initial_threshold = aInitialThreshold; };
    
    double get_sampling_radius() const { return m_sampling_radius; };
    void set_sampling_radius(double aSamplingRadius) { m_sampling_radius = aSamplingRadius; };
    
    std::size_t get_max_result_count() const { return max_num_results; };
    void set_max_result_count(std::size_t aMaxResultCount) { max_num_results = aMaxResultCount; };
    
    std::size_t get_graph_kind_flag() const { return m_graph_kind_flag; };
    void set_graph_kind_flag(std::size_t aGraphKindFlag) { m_graph_kind_flag = aGraphKindFlag; };
    
    std::size_t get_knn_flag() const { return m_knn_flag; };
    void set_knn_flag(std::size_t aKNNMethodFlag) { m_knn_flag = aKNNMethodFlag; };
    
    
    /**
     * Parametrized constructor.
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     * \param aStartPos The position value of the starting location.
     * \param aGoalPos The position value of the goal location.
     * \param aInitialThreshold The initial threshold for exploring nodes, should be somewhere between 0 and 1.
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
     * \param aReporter The SBPP reporter object to use to report results and progress.
     * \param aMaxResultCount The maximum number of successful start-goal connections to make before 
     *                        stopping the path planner (the higher the number the more likely that a 
     *                        good path will be found, however, running time can become much longer).
     */
    sbastar_path_planner(const shared_ptr< space_type >& aWorld = shared_ptr< space_type >(), 
                         const point_type& aStartPos = point_type(),
                         const point_type& aGoalPos = point_type(),
                         double aInitialThreshold = 0.5,
                         double aSamplingRadius = 1.0,
                         std::size_t aMaxVertexCount = 5000, 
                         std::size_t aProgressInterval = 100,
                         std::size_t aGraphKindFlag = ADJ_LIST_MOTION_GRAPH,
                         std::size_t aKNNMethodFlag = DVP_BF2_TREE_KNN,
                         SBPPReporter aReporter = SBPPReporter(),
                         std::size_t aMaxResultCount = 50) :
                         base_type("sbastar_planner", aWorld, aMaxVertexCount, aProgressInterval),
                         m_reporter(aReporter),
                         m_start_pos(aStartPos),
                         m_goal_pos(aGoalPos),
                         m_initial_threshold(aInitialThreshold),
                         m_sampling_radius(aSamplingRadius),
                         max_num_results(aMaxResultCount),
                         has_reached_max_vertices(false),
                         m_graph_kind_flag(aGraphKindFlag),
                         m_knn_flag(aKNNMethodFlag) { };
    
    virtual ~sbastar_path_planner() { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_reporter)
        & RK_SERIAL_SAVE_WITH_NAME(m_start_pos)
        & RK_SERIAL_SAVE_WITH_NAME(m_goal_pos)
        & RK_SERIAL_SAVE_WITH_NAME(m_initial_threshold)
        & RK_SERIAL_SAVE_WITH_NAME(m_sampling_radius)
        & RK_SERIAL_SAVE_WITH_NAME(max_num_results)
        & RK_SERIAL_SAVE_WITH_NAME(m_graph_kind_flag)
        & RK_SERIAL_SAVE_WITH_NAME(m_knn_flag);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_reporter)
        & RK_SERIAL_LOAD_WITH_NAME(m_start_pos)
        & RK_SERIAL_LOAD_WITH_NAME(m_goal_pos)
        & RK_SERIAL_LOAD_WITH_NAME(m_initial_threshold)
        & RK_SERIAL_LOAD_WITH_NAME(m_sampling_radius)
        & RK_SERIAL_LOAD_WITH_NAME(max_num_results)
        & RK_SERIAL_LOAD_WITH_NAME(m_graph_kind_flag)
        & RK_SERIAL_LOAD_WITH_NAME(m_knn_flag);
      has_reached_max_vertices = false;
      m_solutions.clear();
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC246000C,1,"sbastar_path_planner",base_type)
};




template <typename FreeSpaceType, typename MotionGraph, typename NNFinderSynchro, typename SBPPReporter = no_sbmp_report>
struct sbastar_planner_visitor {
  typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
  
  shared_ptr< FreeSpaceType > m_space;
  sbastar_path_planner<FreeSpaceType,SBPPReporter>* m_planner;
  NNFinderSynchro m_nn_synchro;
  Vertex m_start_node;
  Vertex m_goal_node;
  double m_space_dim;
  
  sbastar_planner_visitor(const shared_ptr< FreeSpaceType >& aSpace, 
                          sbastar_path_planner<FreeSpaceType,SBPPReporter>* aPlanner,
                          NNFinderSynchro aNNSynchro,
                          Vertex aStartNode, Vertex aGoalNode, double aSpaceDim) : 
                          m_space(aSpace), m_planner(aPlanner), m_nn_synchro(aNNSynchro),
                          m_start_node(aStartNode), m_goal_node(aGoalNode), m_space_dim(aSpaceDim) { };
  
  typedef typename topology_traits<FreeSpaceType>::point_type PointType;
  
  template <typename Vertex, typename Graph>
  void vertex_added(Vertex u, Graph& g) const {
    m_nn_synchro.added_vertex(u,g);
    
    g[u].heuristic_value = get(distance_metric, m_space->get_super_space())(
      g[m_start_node].position,
      g[u].position,
      m_space->get_super_space());
    
    g[u].constriction = 0.0;
    g[u].collision_count = 0;  // r
    g[u].density = 0.0;
    g[u].expansion_trials = 0;  // m
    
    // Call progress reporter...
    m_planner->report_progress(g);
    
  };
  
  template <typename EdgeType, typename Graph>
  void edge_added(EdgeType e, Graph& g) const {
    using std::exp;
    
    g[e].astar_weight = get(distance_metric, m_space->get_super_space())(
      g[source(e,g)].position,
      g[target(e,g)].position,
      m_space->get_super_space());
    
    double exp_value = exp(-g[e].astar_weight * g[e].astar_weight / (m_planner->get_sampling_radius() * m_planner->get_sampling_radius() * 2.0));
    
    g[source(e,g)].density = (g[source(e,g)].expansion_trials * g[source(e,g)].density + exp_value) / (g[source(e,g)].expansion_trials + 1);
    ++(g[source(e,g)].expansion_trials);
    g[target(e,g)].density = (g[target(e,g)].expansion_trials * g[target(e,g)].density + exp_value) / (g[target(e,g)].expansion_trials + 1);
    ++(g[target(e,g)].expansion_trials);
    
  };
  
  bool keep_going() const {
    return m_planner->keep_going();
  };
  
  template <typename Vertex, typename Graph>
  std::pair<PointType, bool> random_walk(Vertex u, Graph& g) const {
    using std::exp;
    using std::log;
    
    typename point_distribution_traits< typename subspace_traits<FreeSpaceType>::super_space_type >::random_sampler_type get_sample = get(random_sampler, m_space->get_super_space());
    typename metric_space_traits< typename subspace_traits<FreeSpaceType>::super_space_type >::distance_metric_type get_distance = get(distance_metric, m_space->get_super_space());
    
    unsigned int i = 0;
    do {
      PointType p_rnd = get_sample(m_space->get_super_space());
      double dist = get_distance(g[u].position, p_rnd, m_space->get_super_space());
      double target_dist = boost::uniform_01<global_rng_type&,double>(get_global_rng())() * m_planner->get_sampling_radius();
      PointType p_v = m_space->move_position_toward(g[u].position, target_dist / dist, p_rnd);
      dist = get_distance(g[u].position, p_v, m_space->get_super_space());
      if( dist < 0.9 * target_dist ) {
        // this means that we had a collision before reaching the target distance, 
        // must record that to the constriction statistic:
        double sig2_n = m_planner->get_sampling_radius() * m_planner->get_sampling_radius();
        double sig2_x = (target_dist - dist) * (target_dist - dist);
        double exp_D_KL = exp(-target_dist * target_dist / (sig2_x * 2.0) - 0.5 * m_space_dim * ( sig2_n / sig2_x - 1.0 - log(sig2_n / sig2_x) ) );
        g[u].constriction = ( g[u].collision_count * g[u].constriction + exp_D_KL ) / (g[u].collision_count + 1);
        ++(g[u].collision_count);
        // and to the expansion attempts statistic:
        g[u].density = ( g[u].expansion_trials * g[u].density + exp_D_KL ) / (g[u].expansion_trials + 1);
        ++(g[u].expansion_trials);
      } else
        return std::make_pair(p_v, true);
    } while(++i <= 10);
    return std::make_pair(g[u].position, false);
  };
    
  template <typename Vertex, typename Graph>
  void examine_neighborhood(Vertex, Graph&) const {
    
    // NOTE : don't know what is needed here exactly, maybe nothing at all 
    // since density and contriction are computed in added-edge and random-walk.
    
  };
  
  
  template <typename Vertex, typename Graph>
  void initialize_vertex(Vertex u, Graph& g) const {
    typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
    using std::exp;
    
    g[u].heuristic_value = get(ReaK::pp::distance_metric, m_space->get_super_space())(
      g[m_start_node].position,
      g[u].position,
      m_space->get_super_space());
    
    g[u].constriction = 0.0;
    g[u].collision_count = 0;  // r
    
    g[u].density = 0.0;
    g[u].expansion_trials = 0;  // m
    OutEdgeIter ei, ei_end;
    const double denom = m_planner->get_sampling_radius() * m_planner->get_sampling_radius() * 2.0;
    for(boost::tie(ei,ei_end) = out_edges(u,g); ei != ei_end; ++ei) {
      g[u].density = (g[u].expansion_trials * g[u].density + exp(-g[*ei].astar_weight * g[*ei].astar_weight / denom)) / (g[u].expansion_trials + 1);
      ++(g[u].expansion_trials);
    };
    
  };
  
  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex, const Graph&) const { };
  template <typename Vertex, typename Graph>
  void examine_vertex(Vertex, const Graph&) const { };
  template <typename Edge, typename Graph>
  void examine_edge(Edge, const Graph&) const { };
  template <typename Edge, typename Graph>
  void edge_relaxed(Edge, const Graph&) const { };
  template <typename Vertex, typename Graph>
  void finish_vertex(Vertex, const Graph&) const { };
  
  template <typename Graph>
  void publish_path(const Graph& g) const {
    // try to create a goal connection path
    m_planner->create_solution_path(m_start_node, m_goal_node, g); 
  };
  
  template <typename Graph>
  double adjust_threshold(double old_thr, const Graph&) const { 
    return old_thr * 0.5;  // geometrically progress towards 0, from above.
  };
  
  
};






template <typename FreeSpaceType, 
          typename SBPPReporter>
shared_ptr< path_base< typename sbastar_path_planner<FreeSpaceType,SBPPReporter>::super_space_type > > 
  sbastar_path_planner<FreeSpaceType,SBPPReporter>::solve_path() {
  using ReaK::to_vect;
  
  this->has_reached_max_vertices = false;
  this->m_solutions.clear();
  
  typedef typename subspace_traits<FreeSpaceType>::super_space_type SuperSpace;
  typedef typename topology_traits<SuperSpace>::point_type PointType;
  
  typedef boost::data_member_property_map<PointType, sbastar_vertex_data<FreeSpaceType> > PositionMap;
  PositionMap pos_map = PositionMap(&sbastar_vertex_data<FreeSpaceType>::position);
  
  double space_dim = double((to_vect<double>(this->m_space->get_super_space().difference(this->m_goal_pos,this->m_start_pos))).size()); 
  double space_Lc = get(distance_metric,this->m_space->get_super_space())(this->m_start_pos, this->m_goal_pos, this->m_space->get_super_space());
  
//   double max_radius = 2.0 * m_sampling_radius;
  
  if(m_graph_kind_flag == ADJ_LIST_MOTION_GRAPH) {
    
    typedef boost::adjacency_list< 
      boost::vecS, boost::vecS, boost::undirectedS,
      sbastar_vertex_data<FreeSpaceType>,
      sbastar_edge_data<FreeSpaceType>, boost::listS> MotionGraphType;
    typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex;
    typedef typename MotionGraphType::vertex_property_type VertexProp;
    typedef boost::composite_property_map< 
      PositionMap, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t > > GraphPositionMap;
    
    MotionGraphType motion_graph;
    GraphPositionMap g_pos_map = GraphPositionMap(pos_map, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t >(&motion_graph));
    
    sbastar_vertex_data<FreeSpaceType> vs_p, vg_p;
    vs_p.position = this->m_start_pos;
    vg_p.position = this->m_goal_pos;
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    Vertex start_node = add_vertex(std::move(vs_p), motion_graph);
    Vertex goal_node = add_vertex(std::move(vg_p), motion_graph);
#else
    Vertex start_node = add_vertex(vs_p, motion_graph);
    Vertex goal_node = add_vertex(vg_p, motion_graph);
#endif
    motion_graph[start_node].constriction = 0.0;
    motion_graph[start_node].collision_count = 0;
    motion_graph[start_node].density = 0.0;
    motion_graph[start_node].expansion_trials = 0;
    motion_graph[start_node].heuristic_value = 0.0;  // distance to start node.
    motion_graph[start_node].distance_accum = std::numeric_limits<double>::infinity();
    motion_graph[start_node].key_value = 0.0;
    motion_graph[start_node].astar_color = boost::color_traits<boost::default_color_type>::white();
    motion_graph[start_node].predecessor = start_node;
    
    motion_graph[goal_node].constriction = 0.0;
    motion_graph[goal_node].collision_count = 0;
    motion_graph[goal_node].density = 0.0;
    motion_graph[goal_node].expansion_trials = 0;
    motion_graph[goal_node].heuristic_value = space_Lc;
    motion_graph[goal_node].distance_accum = 0.0;
    motion_graph[goal_node].key_value = 1.0 / space_Lc;
    motion_graph[goal_node].astar_color = boost::color_traits<boost::default_color_type>::white();
    motion_graph[goal_node].predecessor = goal_node;
    
    
    if(m_knn_flag == LINEAR_SEARCH_KNN) {
      sbastar_planner_visitor<FreeSpaceType, MotionGraphType, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node, space_dim);
      
      ReaK::graph::generate_sbastar(
        motion_graph, goal_node, *(this->m_space), vis,
        get(&sbastar_vertex_data<FreeSpaceType>::heuristic_value, motion_graph), 
        pos_map, 
        get(&sbastar_edge_data<FreeSpaceType>::astar_weight, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::density, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::constriction, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::distance_accum, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::predecessor, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::key_value, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::astar_color, motion_graph),
//         ReaK::graph::fixed_neighborhood< linear_neighbor_search<> >(
//           linear_neighbor_search<>(), 
//           10, max_radius),
        ReaK::graph::star_neighborhood< linear_neighbor_search<> >(
          linear_neighbor_search<>(), 
          space_dim, 3.0 * space_Lc),
        this->m_initial_threshold);
      
    } else if(m_knn_flag == DVP_BF2_TREE_KNN) {
      
      typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 2, 
                       random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2> > SpacePartType;
      SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
      
      multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraphType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node, space_dim);
      
      ReaK::graph::generate_sbastar(
        motion_graph, goal_node, *(this->m_space), vis,
        get(&sbastar_vertex_data<FreeSpaceType>::heuristic_value, motion_graph), 
        pos_map, 
        get(&sbastar_edge_data<FreeSpaceType>::astar_weight, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::density, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::constriction, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::distance_accum, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::predecessor, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::key_value, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::astar_color, motion_graph),
//         ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> >(
//           nn_finder, 
//           10, max_radius),
        ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> >(
          nn_finder, 
          space_dim, 3.0 * space_Lc),
        this->m_initial_threshold);
      
    } else if(m_knn_flag == DVP_BF4_TREE_KNN) {
      
      typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 4, 
                       random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4> > SpacePartType;
      SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
      
      multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraphType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node, space_dim);
      
      ReaK::graph::generate_sbastar(
        motion_graph, goal_node, *(this->m_space), vis,
        get(&sbastar_vertex_data<FreeSpaceType>::heuristic_value, motion_graph), 
        pos_map, 
        get(&sbastar_edge_data<FreeSpaceType>::astar_weight, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::density, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::constriction, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::distance_accum, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::predecessor, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::key_value, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::astar_color, motion_graph),
//         ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> >(
//           nn_finder, 
//           10, max_radius),
        ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> >(
          nn_finder, 
          space_dim, 3.0 * space_Lc),
        this->m_initial_threshold);
      
    } else if(m_knn_flag == DVP_COB2_TREE_KNN) {
      
      typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 2, 
                       random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2> > SpacePartType;
      SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
      
      multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraphType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node, space_dim);
      
      ReaK::graph::generate_sbastar(
        motion_graph, goal_node, *(this->m_space), vis,
        get(&sbastar_vertex_data<FreeSpaceType>::heuristic_value, motion_graph), 
        pos_map, 
        get(&sbastar_edge_data<FreeSpaceType>::astar_weight, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::density, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::constriction, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::distance_accum, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::predecessor, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::key_value, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::astar_color, motion_graph),
//         ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> >(
//           nn_finder, 
//           10, max_radius),
        ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> >(
          nn_finder, 
          space_dim, 3.0 * space_Lc),
        this->m_initial_threshold);
      
    } else if(m_knn_flag == DVP_COB4_TREE_KNN) {
      
      typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 4, 
                       random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4> > SpacePartType;
      SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
      
      multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraphType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder, start_node, goal_node, space_dim);
      
      ReaK::graph::generate_sbastar(
        motion_graph, goal_node, *(this->m_space), vis,
        get(&sbastar_vertex_data<FreeSpaceType>::heuristic_value, motion_graph), 
        pos_map, 
        get(&sbastar_edge_data<FreeSpaceType>::astar_weight, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::density, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::constriction, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::distance_accum, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::predecessor, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::key_value, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::astar_color, motion_graph),
//         ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> >(
//           nn_finder, 
//           10, max_radius),
        ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraphType, SpacePartType> >(
          nn_finder, 
          space_dim, 3.0 * space_Lc),
        this->m_initial_threshold);
      
    };
    
  } else if(m_graph_kind_flag == DVP_ADJ_LIST_MOTION_GRAPH) {
    
    if(m_knn_flag == DVP_ALT_BF2_KNN) {
      
      typedef dvp_adjacency_list<
        sbastar_vertex_data<FreeSpaceType>,
        sbastar_edge_data<FreeSpaceType>,
        SuperSpace,
        PositionMap,
        2, random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2>,
        boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
      
      ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
      
      typedef typename ALTGraph::adj_list_type MotionGraph;
      typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
      
      MotionGraph motion_graph = space_part.get_adjacency_list();
      
      sbastar_vertex_data<FreeSpaceType> vs_p, vg_p;
      vs_p.position = this->m_start_pos;
      vg_p.position = this->m_goal_pos;
    
#ifdef RK_ENABLE_CXX0X_FEATURES
      Vertex start_node = add_vertex(std::move(vs_p), motion_graph);
      Vertex goal_node = add_vertex(std::move(vg_p), motion_graph);
#else
      Vertex start_node = add_vertex(vs_p, motion_graph);
      Vertex goal_node = add_vertex(vg_p, motion_graph);
#endif
      motion_graph[start_node].constriction = 0.0;
      motion_graph[start_node].collision_count = 0;
      motion_graph[start_node].density = 0.0;
      motion_graph[start_node].expansion_trials = 0;
      motion_graph[start_node].heuristic_value = 0.0;  // distance to start node.
      motion_graph[start_node].distance_accum = std::numeric_limits<double>::infinity();
      motion_graph[start_node].key_value = 0.0;
      motion_graph[start_node].astar_color = boost::color_traits<boost::default_color_type>::white();
      motion_graph[start_node].predecessor = start_node;
      
      motion_graph[goal_node].constriction = 0.0;
      motion_graph[goal_node].collision_count = 0;
      motion_graph[goal_node].density = 0.0;
      motion_graph[goal_node].expansion_trials = 0;
      motion_graph[goal_node].heuristic_value = space_Lc;
      motion_graph[goal_node].distance_accum = 0.0;
      motion_graph[goal_node].key_value = 1.0 / space_Lc;
      motion_graph[goal_node].astar_color = boost::color_traits<boost::default_color_type>::white();
      motion_graph[goal_node].predecessor = goal_node;
      
      multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node, space_dim);
      
      ReaK::graph::generate_sbastar(
        motion_graph, goal_node, *(this->m_space), vis,
        get(&sbastar_vertex_data<FreeSpaceType>::heuristic_value, motion_graph), 
        pos_map, 
        get(&sbastar_edge_data<FreeSpaceType>::astar_weight, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::density, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::constriction, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::distance_accum, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::predecessor, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::key_value, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::astar_color, motion_graph),
//         ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> >(
//           nn_finder, 
//           10, max_radius),
        ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> >(
          nn_finder, 
          space_dim, 3.0 * space_Lc),
        this->m_initial_threshold);
      
    } else if(m_knn_flag == DVP_ALT_BF4_KNN) {
      
      typedef dvp_adjacency_list<
        sbastar_vertex_data<FreeSpaceType>,
        sbastar_edge_data<FreeSpaceType>,
        SuperSpace,
        PositionMap,
        4, random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4>,
        boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
      
      ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
      
      typedef typename ALTGraph::adj_list_type MotionGraph;
      typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
      
      MotionGraph motion_graph = space_part.get_adjacency_list();
      
      sbastar_vertex_data<FreeSpaceType> vs_p, vg_p;
      vs_p.position = this->m_start_pos;
      vg_p.position = this->m_goal_pos;
    
#ifdef RK_ENABLE_CXX0X_FEATURES
      Vertex start_node = add_vertex(std::move(vs_p), motion_graph);
      Vertex goal_node = add_vertex(std::move(vg_p), motion_graph);
#else
      Vertex start_node = add_vertex(vs_p, motion_graph);
      Vertex goal_node = add_vertex(vg_p, motion_graph);
#endif
      motion_graph[start_node].constriction = 0.0;
      motion_graph[start_node].collision_count = 0;
      motion_graph[start_node].density = 0.0;
      motion_graph[start_node].expansion_trials = 0;
      motion_graph[start_node].heuristic_value = 0.0;  // distance to start node.
      motion_graph[start_node].distance_accum = std::numeric_limits<double>::infinity();
      motion_graph[start_node].key_value = 0.0;
      motion_graph[start_node].astar_color = boost::color_traits<boost::default_color_type>::white();
      motion_graph[start_node].predecessor = start_node;
      
      motion_graph[goal_node].constriction = 0.0;
      motion_graph[goal_node].collision_count = 0;
      motion_graph[goal_node].density = 0.0;
      motion_graph[goal_node].expansion_trials = 0;
      motion_graph[goal_node].heuristic_value = space_Lc;
      motion_graph[goal_node].distance_accum = 0.0;
      motion_graph[goal_node].key_value = 1.0 / space_Lc;
      motion_graph[goal_node].astar_color = boost::color_traits<boost::default_color_type>::white();
      motion_graph[goal_node].predecessor = goal_node;
      
      multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node, space_dim);
      
      ReaK::graph::generate_sbastar(
        motion_graph, goal_node, *(this->m_space), vis,
        get(&sbastar_vertex_data<FreeSpaceType>::heuristic_value, motion_graph), 
        pos_map, 
        get(&sbastar_edge_data<FreeSpaceType>::astar_weight, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::density, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::constriction, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::distance_accum, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::predecessor, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::key_value, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::astar_color, motion_graph),
//         ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> >(
//           nn_finder, 
//           10, max_radius),
        ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> >(
          nn_finder, 
          space_dim, 3.0 * space_Lc),
        this->m_initial_threshold);
      
    } else if(m_knn_flag == DVP_ALT_COB2_KNN) {
      
      typedef dvp_adjacency_list<
        sbastar_vertex_data<FreeSpaceType>,
        sbastar_edge_data<FreeSpaceType>,
        SuperSpace,
        PositionMap,
        2, random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2>,
        boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
      
      ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
      
      typedef typename ALTGraph::adj_list_type MotionGraph;
      typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
      
      MotionGraph motion_graph = space_part.get_adjacency_list();
      
      sbastar_vertex_data<FreeSpaceType> vs_p, vg_p;
      vs_p.position = this->m_start_pos;
      vg_p.position = this->m_goal_pos;
    
#ifdef RK_ENABLE_CXX0X_FEATURES
      Vertex start_node = add_vertex(std::move(vs_p), motion_graph);
      Vertex goal_node = add_vertex(std::move(vg_p), motion_graph);
#else
      Vertex start_node = add_vertex(vs_p, motion_graph);
      Vertex goal_node = add_vertex(vg_p, motion_graph);
#endif
      motion_graph[start_node].constriction = 0.0;
      motion_graph[start_node].collision_count = 0;
      motion_graph[start_node].density = 0.0;
      motion_graph[start_node].expansion_trials = 0;
      motion_graph[start_node].heuristic_value = 0.0;  // distance to start node.
      motion_graph[start_node].distance_accum = std::numeric_limits<double>::infinity();
      motion_graph[start_node].key_value = 0.0;
      motion_graph[start_node].astar_color = boost::color_traits<boost::default_color_type>::white();
      motion_graph[start_node].predecessor = start_node;
      
      motion_graph[goal_node].constriction = 0.0;
      motion_graph[goal_node].collision_count = 0;
      motion_graph[goal_node].density = 0.0;
      motion_graph[goal_node].expansion_trials = 0;
      motion_graph[goal_node].heuristic_value = space_Lc;
      motion_graph[goal_node].distance_accum = 0.0;
      motion_graph[goal_node].key_value = 1.0 / space_Lc;
      motion_graph[goal_node].astar_color = boost::color_traits<boost::default_color_type>::white();
      motion_graph[goal_node].predecessor = goal_node;
      
      multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node, space_dim);
      
      ReaK::graph::generate_sbastar(
        motion_graph, goal_node, *(this->m_space), vis,
        get(&sbastar_vertex_data<FreeSpaceType>::heuristic_value, motion_graph), 
        pos_map, 
        get(&sbastar_edge_data<FreeSpaceType>::astar_weight, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::density, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::constriction, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::distance_accum, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::predecessor, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::key_value, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::astar_color, motion_graph),
//         ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> >(
//           nn_finder, 
//           10, max_radius),
        ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> >(
          nn_finder, 
          space_dim, 3.0 * space_Lc),
        this->m_initial_threshold);
      
    } else if(m_knn_flag == DVP_ALT_COB4_KNN) {
      
      typedef dvp_adjacency_list<
        sbastar_vertex_data<FreeSpaceType>,
        sbastar_edge_data<FreeSpaceType>,
        SuperSpace,
        PositionMap,
        4, random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4>,
        boost::vecS, boost::undirectedS, boost::listS > ALTGraph;
      
      ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
      
      typedef typename ALTGraph::adj_list_type MotionGraph;
      typedef typename boost::graph_traits<MotionGraph>::vertex_descriptor Vertex;
      
      MotionGraph motion_graph = space_part.get_adjacency_list();
      
      sbastar_vertex_data<FreeSpaceType> vs_p, vg_p;
      vs_p.position = this->m_start_pos;
      vg_p.position = this->m_goal_pos;
    
#ifdef RK_ENABLE_CXX0X_FEATURES
      Vertex start_node = add_vertex(std::move(vs_p), motion_graph);
      Vertex goal_node = add_vertex(std::move(vg_p), motion_graph);
#else
      Vertex start_node = add_vertex(vs_p, motion_graph);
      Vertex goal_node = add_vertex(vg_p, motion_graph);
#endif
      motion_graph[start_node].constriction = 0.0;
      motion_graph[start_node].collision_count = 0;
      motion_graph[start_node].density = 0.0;
      motion_graph[start_node].expansion_trials = 0;
      motion_graph[start_node].heuristic_value = 0.0;  // distance to start node.
      motion_graph[start_node].distance_accum = std::numeric_limits<double>::infinity();
      motion_graph[start_node].key_value = 0.0;
      motion_graph[start_node].astar_color = boost::color_traits<boost::default_color_type>::white();
      motion_graph[start_node].predecessor = start_node;
      
      motion_graph[goal_node].constriction = 0.0;
      motion_graph[goal_node].collision_count = 0;
      motion_graph[goal_node].density = 0.0;
      motion_graph[goal_node].expansion_trials = 0;
      motion_graph[goal_node].heuristic_value = space_Lc;
      motion_graph[goal_node].distance_accum = 0.0;
      motion_graph[goal_node].key_value = 1.0 / space_Lc;
      motion_graph[goal_node].astar_color = boost::color_traits<boost::default_color_type>::white();
      motion_graph[goal_node].predecessor = goal_node;
      
      multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
      nn_finder.graph_tree_map[&motion_graph] = &space_part;
      
      sbastar_planner_visitor<FreeSpaceType, MotionGraph, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro(), start_node, goal_node, space_dim);
      
      ReaK::graph::generate_sbastar(
        motion_graph, goal_node, *(this->m_space), vis,
        get(&sbastar_vertex_data<FreeSpaceType>::heuristic_value, motion_graph), 
        pos_map, 
        get(&sbastar_edge_data<FreeSpaceType>::astar_weight, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::density, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::constriction, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::distance_accum, motion_graph),
        get(&sbastar_vertex_data<FreeSpaceType>::predecessor, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::key_value, motion_graph), 
        get(&sbastar_vertex_data<FreeSpaceType>::astar_color, motion_graph),
//         ReaK::graph::fixed_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> >(
//           nn_finder, 
//           10, max_radius),
        ReaK::graph::star_neighborhood< multi_dvp_tree_search<MotionGraph, ALTGraph> >(
          nn_finder, 
          space_dim, 3.0 * space_Lc),
        this->m_initial_threshold);
      
    };
    
  };
  
  if(m_solutions.size())
    return m_solutions.begin()->second;
  else
    return shared_ptr< path_base< SuperSpace > >();
};


};

};

#endif


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
#include "path_planning/sbmp_reporter_concept.hpp"

#include "path_planning/metric_space_concept.hpp"
#include "path_planning/seq_path_wrapper.hpp"
#include "interpolation/point_to_point_path.hpp"
#include "path_planning/basic_sbmp_reporters.hpp"

#include "graph_alg/rrt_star.hpp"

#include "graph_alg/d_ary_bf_tree.hpp"
#include "graph_alg/d_ary_cob_tree.hpp"
#include "graph_alg/bgl_tree_adaptor.hpp"
#include "graph_alg/bgl_more_property_maps.hpp"
#include "path_planning/dvp_layout_adjacency_list.hpp"
#include "path_planning/metric_space_search.hpp"
#include "path_planning/topological_search.hpp"
#include "path_planning/path_planner_options.hpp"
#include "path_planning/rrt_path_planner.hpp"
#include "path_planning/rrtstar_path_planner.hpp"
#include "lin_alg/arithmetic_tuple.hpp"

#include "MEAQR_topology.hpp"
#include "topologies/fixed_topology_random_sampler.hpp"

namespace ReaK {
  
namespace pp {
  


template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
struct MEAQR_rrtstar_vdata {
  typedef typename topology_traits< MEAQR_topology<StateSpace, StateSpaceSystem, StateSpaceSampler> >::point_type PointType;
  
  PointType position;
  double distance_accum;
  /// The predecessor associated to the vertex, i.e., following the predecessor links starting at the goal node yields a backward trace of the optimal path.
  std::size_t predecessor;
  
  MEAQR_rrtstar_vdata() : position(PointType()), distance_accum(0.0) { };
};

template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler>
struct MEAQR_rrtstar_edata {
  typedef typename MEAQR_topology<StateSpace, StateSpaceSystem, StateSpaceSampler>::steer_record_type steer_record_type;
  
  double weight;
  steer_record_type steer_rec;
  
  MEAQR_rrtstar_edata(double aWeight = 0.0, const steer_record_type& aSteerRec = steer_record_type()) : weight(aWeight), steer_rec(aSteerRec) { };
#ifdef RK_ENABLE_CXX11_FEATURES
  MEAQR_rrtstar_edata(double aWeight, steer_record_type&& aSteerRec) : weight(aWeight), steer_rec(std::move(aSteerRec)) { };
#endif
};





/**
 * This class is a RRT*-based path-planner over the a MEAQR-controlled system-topology.
 * \tparam StateSpace The topology type of state-space of the dynamic system under MEAQR control.
 * \tparam StateSpaceSystem The type of the dynamic system under MEAQR control.
 * \tparam SBPPReporter The reporter type to use to report the progress of the path-planning.
 */
template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler,
          typename SBPPReporter = no_sbmp_report>
class MEAQR_rrtstar_planner : public sample_based_planner< path_planner_base< MEAQR_topology_with_CD<StateSpace, StateSpaceSystem, StateSpaceSampler> > > {
  public:
    typedef MEAQR_topology_with_CD<StateSpace, StateSpaceSystem, StateSpaceSampler> space_type;
    typedef typename subspace_traits<space_type>::super_space_type super_space_type;
    
    typedef sample_based_planner< path_planner_base< space_type > > base_type;
    typedef MEAQR_rrtstar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler, SBPPReporter> self;
    
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename topology_traits< super_space_type >::point_difference_type point_difference_type;
    
    typedef typename steerable_space_traits< super_space_type >::steer_record_type steer_record_type;
    
  protected:
    SBPPReporter m_reporter;
    point_type m_start_pos;
    point_type m_goal_pos;
    std::size_t max_num_results;
    bool has_reached_max_vertices;
    std::size_t m_knn_flag;
    
    std::map<double, shared_ptr< seq_path_base< super_space_type > > > m_solutions;
    
    
  public:
    
    template <typename Vertex, typename Graph>
    void check_goal_connection(const point_type& p_u, Vertex u, Graph& g) {
      typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
      typedef typename steer_record_type::point_fraction_iterator SteerIter;
      
      std::pair< point_type, steer_record_type > steer_result = this->m_space->steer_position_toward(p_u, 1.0, this->m_goal_pos);
      // NOTE Differs from rrtstar_path_planner HERE:
      double diff_dist = get(distance_metric, this->m_space->get_state_space())(steer_result.first.x, this->m_goal_pos.x, this->m_space->get_state_space());
      double actual_dist = get(distance_metric, this->m_space->get_state_space())(p_u.x, this->m_goal_pos.x, this->m_space->get_state_space());
      
      if(diff_dist > 0.05 * actual_dist)
        return;
      
      std::cout << " starting to print out solution... " << std::endl;
      
      double solutions_total_dist = g[u].distance_accum + get(distance_metric, this->m_space->get_super_space())(p_u, this->m_goal_pos, this->m_space->get_super_space());
//       if(solutions_total_dist >= this->m_solutions.begin()->first)
//         return;
      
      shared_ptr< super_space_type > sup_space_ptr(&(this->m_space->get_super_space()),null_deleter());
      shared_ptr< seq_path_wrapper< discrete_point_path<super_space_type> > > new_sol(new seq_path_wrapper< discrete_point_path<super_space_type> >("MEAQR_rrtstar_solution", discrete_point_path<super_space_type>(sup_space_ptr,get(distance_metric, this->m_space->get_super_space()))));
      discrete_point_path<super_space_type>& waypoints = new_sol->get_underlying_path();
      
      for(SteerIter it = steer_result.second.end_fraction_travel(); it != steer_result.second.begin_fraction_travel(); it -= 10.0)
        waypoints.push_front( point_type(*it) );
      
      while(in_degree(u,g)) {
        Edge e = *(in_edges(u,g).first);
        u = source(e,g);
        for(SteerIter it = g[e].steer_rec.end_fraction_travel(); it != g[e].steer_rec.begin_fraction_travel(); it -= 10.0) {
          waypoints.push_front( point_type(*it) );
        };
      };
      
      this->m_solutions[solutions_total_dist] = new_sol;
      this->m_reporter.draw_solution(*(this->m_space), this->m_solutions[solutions_total_dist]);
      
      std::cout << "Done!" << std::endl;
    };
    
    // NOTE: This is the same as rrtstar_path_planner
    template <typename Vertex, typename Graph>
    void joining_vertex_found(Vertex u1, Vertex u2, Graph& g1, Graph& g2) {
      typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
      typedef typename steer_record_type::point_fraction_iterator SteerIter;
      
      double total_dist = g1[u1].distance_accum + g2[u2].distance_accum
        + get(distance_metric, this->m_space->get_super_space())(g1[u1].position, g2[u2].position, this->m_space->get_super_space());
      
//       if(total_dist >= m_solutions.begin()->first)
//         return;
      
      shared_ptr< super_space_type > sup_space_ptr(&(this->m_space->get_super_space()),null_deleter());
      shared_ptr< seq_path_wrapper< discrete_point_path<super_space_type> > > new_sol(new seq_path_wrapper< discrete_point_path<super_space_type> >("birrt_solution", point_to_point_path<super_space_type>(sup_space_ptr,get(distance_metric, this->m_space->get_super_space()))));
      discrete_point_path<super_space_type>& waypoints = new_sol->get_underlying_path();
      
      while(in_degree(u1,g1)) {
        Edge e = *(in_edges(u1,g1).first);
        u1 = source(e,g1);
        for(SteerIter it = g1[e].steer_rec.end_fraction_travel(); it != g1[e].steer_rec.begin_fraction_travel(); it -= 10.0)
          waypoints.push_front( point_type(*it) );
      };
      
      while(in_degree(u2,g2)) {
        Edge e = *(in_edges(u2,g2).first);
        u2 = source(e,g2);
        for(SteerIter it = g2[e].steer_rec.begin_fraction_travel(); it != g2[e].steer_rec.end_fraction_travel(); it += 10.0)
          waypoints.push_back( point_type(*it) );
      };
      
      m_solutions[total_dist] = new_sol;
      m_reporter.draw_solution(*(this->m_space), m_solutions[total_dist]);
    };
    
    // NOTE: This is the same as rrtstar_path_planner
    bool keep_going() const {
      return (max_num_results > m_solutions.size()) && !has_reached_max_vertices;
    };
    
    // NOTE: This is the same as rrtstar_path_planner
    template <typename Graph, typename PositionMap>
    void report_progress(Graph& g, PositionMap g_pos) {
      if(num_vertices(g) % this->m_progress_interval == 0)
        m_reporter.draw_motion_graph(*(this->m_space), g, g_pos);
      has_reached_max_vertices = (num_vertices(g) >= this->m_max_vertex_count);
    };
    
    // NOTE: This is the same as rrtstar_path_planner
    template <typename Graph, typename PositionMap>
    void report_progress(Graph& g1, Graph& g2, PositionMap g1_pos, PositionMap g2_pos) {
      if((num_vertices(g1) + num_vertices(g2)) % this->m_progress_interval == 0) {
        m_reporter.draw_motion_graph(*(this->m_space), g1, g1_pos);
        m_reporter.draw_motion_graph(*(this->m_space), g2, g2_pos);
      };
      has_reached_max_vertices = (num_vertices(g1) + num_vertices(g2) >= this->m_max_vertex_count);
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
    
    // NOTE: This is the same as rrtstar_path_planner
    /**
     * Parametrized constructor.
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     * \param aStartPos The position value of the starting location.
     * \param aGoalPos The position value of the goal location.
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
    MEAQR_rrtstar_planner(const shared_ptr< space_type >& aWorld = shared_ptr< space_type >(), 
                          const point_type& aStartPos = point_type(),
                          const point_type& aGoalPos = point_type(),
                          std::size_t aMaxVertexCount = 5000, 
                          std::size_t aProgressInterval = 100,
                          std::size_t aKNNMethodFlag = DVP_BF2_TREE_KNN,
                          SBPPReporter aReporter = SBPPReporter(),
                          std::size_t aMaxResultCount = 50) :
                          base_type("MEAQR_rrtstar_planner", aWorld, aMaxVertexCount, aProgressInterval),
                          m_reporter(aReporter),
                          m_start_pos(aStartPos),
                          m_goal_pos(aGoalPos),
                          max_num_results(aMaxResultCount),
                          has_reached_max_vertices(false),
                          m_knn_flag(aKNNMethodFlag),
                          m_solutions() { };
    
    virtual ~MEAQR_rrtstar_planner() { };
    
    const SBPPReporter& get_reporter() const { return m_reporter; };
    void set_reporter(const SBPPReporter& aNewReporter) { m_reporter = aNewReporter; };
    
    const point_type& get_start_pos() const { return m_start_pos; };
    void set_start_pos(const point_type& aStartPos) { m_start_pos = aStartPos; };
    
    const point_type& get_goal_pos() const { return m_goal_pos; };
    void set_goal_pos(const point_type& aGoalPos) { m_goal_pos = aGoalPos; };
    
    std::size_t get_max_result_count() const { return max_num_results; };
    void set_max_result_count(std::size_t aMaxResultCount) { max_num_results = aMaxResultCount; };
    
    std::size_t get_knn_flag() const { return m_knn_flag; };
    void set_knn_flag(std::size_t aKNNMethodFlag) { m_knn_flag = aKNNMethodFlag; };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    // NOTE: This is the same as rrtstar_path_planner
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_reporter)
        & RK_SERIAL_SAVE_WITH_NAME(m_start_pos)
        & RK_SERIAL_SAVE_WITH_NAME(m_goal_pos)
        & RK_SERIAL_SAVE_WITH_NAME(max_num_results)
        & RK_SERIAL_SAVE_WITH_NAME(m_knn_flag);
    };
    
    // NOTE: This is the same as rrtstar_path_planner
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_reporter)
        & RK_SERIAL_LOAD_WITH_NAME(m_start_pos)
        & RK_SERIAL_LOAD_WITH_NAME(m_goal_pos)
        & RK_SERIAL_LOAD_WITH_NAME(max_num_results)
        & RK_SERIAL_LOAD_WITH_NAME(m_knn_flag);
      has_reached_max_vertices = false;
      m_solutions.clear();
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC246000D,1,"MEAQR_rrtstar_planner",base_type)
};




template <typename StateSpace, typename StateSpaceSystem, typename StateSpaceSampler, typename NNFinderSynchro, typename SBPPReporter = no_sbmp_report>
struct MEAQR_rrtstar_visitor {
  typedef MEAQR_topology_with_CD<StateSpace, StateSpaceSystem, StateSpaceSampler> space_type;
  shared_ptr< space_type > m_space;
  MEAQR_rrtstar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler, SBPPReporter>* m_planner;
  NNFinderSynchro m_nn_synchro;
  
  MEAQR_rrtstar_visitor(const shared_ptr< space_type >& aSpace, 
                        MEAQR_rrtstar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler, SBPPReporter>* aPlanner,
                        NNFinderSynchro aNNSynchro) : 
                        m_space(aSpace), m_planner(aPlanner), m_nn_synchro(aNNSynchro) { };
  
  typedef typename topology_traits< space_type >::point_type PointType;
  typedef MEAQR_rrtstar_edata<StateSpace, StateSpaceSystem, StateSpaceSampler> EdgeProp;
  
  template <typename Vertex, typename Graph>
  void vertex_added(Vertex u, Graph& g) const {
    m_nn_synchro.added_vertex(u,g);
    std::cout << "\r" << std::setw(10) << num_vertices(g) << std::flush;
  };
  
  template <typename Vertex, typename Graph>
  void vertex_added(Vertex u, Graph& g, Graph&) const {
    m_nn_synchro.added_vertex(u,g);
  };
  
  template <typename EdgeType, typename Graph>
  void edge_added(EdgeType e, Graph& g) const {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor VertexType;
    VertexType v = target(e,g);
    
    // Call progress reporter...
//     m_planner->report_progress(g, get(&MEAQR_rrtstar_vdata<StateSpace,StateSpaceSystem, StateSpaceSampler>::position,g));
    m_planner->report_progress(g, get(&MEAQR_rrtstar_edata<StateSpace,StateSpaceSystem, StateSpaceSampler>::steer_rec,g));
    
    // Check if a straight path to goal is possible...
    m_planner->check_goal_connection(g[v].position,v,g);
    
  };
  
  template <typename EdgeType, typename Graph>
  void edge_added(EdgeType e, Graph& g1, Graph& g2) const {
    
    // Call progress reporter...
//     m_planner->report_progress(g1, g2, 
//                                get(&MEAQR_rrtstar_vdata<StateSpace,StateSpaceSystem, StateSpaceSampler>::position,g1), 
//                                get(&MEAQR_rrtstar_vdata<StateSpace,StateSpaceSystem, StateSpaceSampler>::position,g2));
    m_planner->report_progress(g1, g2, 
                               get(&MEAQR_rrtstar_edata<StateSpace,StateSpaceSystem, StateSpaceSampler>::steer_rec,g1), 
                               get(&MEAQR_rrtstar_edata<StateSpace,StateSpaceSystem, StateSpaceSampler>::steer_rec,g2));
  };
  
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
  
};






template <typename StateSpace, 
          typename StateSpaceSystem, typename StateSpaceSampler, 
          typename SBPPReporter>
shared_ptr< seq_path_base< typename MEAQR_rrtstar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler, SBPPReporter>::super_space_type > > 
  MEAQR_rrtstar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler, SBPPReporter>::solve_path() {
  using ReaK::to_vect;
  
  typedef typename MEAQR_rrtstar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler, SBPPReporter>::space_type SpaceType;
  typedef typename MEAQR_rrtstar_planner<StateSpace, StateSpaceSystem, StateSpaceSampler, SBPPReporter>::super_space_type SuperSpace;
  typedef typename SuperSpace::IHAQR_space_type IHAQRSpace;
  typedef typename topology_traits<SuperSpace>::point_type PointType;
  
  typedef boost::data_member_property_map<PointType, MEAQR_rrtstar_vdata<StateSpace, StateSpaceSystem, StateSpaceSampler> > PositionMap;
  PositionMap pos_map = PositionMap(&MEAQR_rrtstar_vdata<StateSpace, StateSpaceSystem, StateSpaceSampler>::position);
  
  typedef boost::data_member_property_map<std::size_t, MEAQR_rrtstar_vdata<StateSpace, StateSpaceSystem, StateSpaceSampler> > PredecessorMap;
  PredecessorMap pred_map = PredecessorMap(&MEAQR_rrtstar_vdata<StateSpace, StateSpaceSystem, StateSpaceSampler>::predecessor);
  
  typedef boost::data_member_property_map<double, MEAQR_rrtstar_vdata<StateSpace, StateSpaceSystem, StateSpaceSampler> > CostMap;
  CostMap cost_map = CostMap(&MEAQR_rrtstar_vdata<StateSpace, StateSpaceSystem, StateSpaceSampler>::distance_accum);
  
  typedef boost::data_member_property_map<double, MEAQR_rrtstar_edata<StateSpace, StateSpaceSystem, StateSpaceSampler> > WeightMap;
  WeightMap weight_map = WeightMap(&MEAQR_rrtstar_edata<StateSpace, StateSpaceSystem, StateSpaceSampler>::weight);
  
  typedef fixed_topology_random_sampler< SuperSpace > SamplerType;
  SamplerType get_sample = SamplerType( &(this->m_space->get_super_space()) );
  
//   double space_dim = double((to_vect<double>(this->m_space->get_state_space().difference(this->m_goal_pos.x,this->m_start_pos.x))).size()); 
  double space_Lc = get(distance_metric,this->m_space->get_super_space())(this->m_start_pos, this->m_goal_pos, this->m_space->get_super_space());
  
  typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::bidirectionalS,
                                 MEAQR_rrtstar_vdata<StateSpace, StateSpaceSystem, StateSpaceSampler>,
                                 MEAQR_rrtstar_edata<StateSpace, StateSpaceSystem, StateSpaceSampler>,
                                 boost::vecS> MotionGraphType;
  typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex;
  typedef typename MotionGraphType::vertex_property_type VertexProp;
  typedef boost::composite_property_map< 
    PositionMap, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t > > GraphPositionMap;
  
  MotionGraphType motion_graph;
  GraphPositionMap g_pos_map = GraphPositionMap(pos_map, boost::whole_bundle_property_map< MotionGraphType, boost::vertex_bundle_t >(&motion_graph));
  
  MEAQR_rrtstar_vdata<StateSpace, StateSpaceSystem, StateSpaceSampler> v_p;
  v_p.position = this->m_start_pos;
  v_p.distance_accum = 0.0;
  add_vertex(v_p, motion_graph);
  
  if(this->m_knn_flag == LINEAR_SEARCH_KNN) {
    MEAQR_rrtstar_visitor<StateSpace, StateSpaceSystem, StateSpaceSampler, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro());
    
    ReaK::graph::detail::generate_rrt_star_loop(
      motion_graph, 
      this->m_space->get_super_space(),
      vis, 
      pos_map, 
      cost_map,
      pred_map,
      weight_map,
      ReaK::graph::rrg_node_generator<
        SuperSpace, 
        SamplerType, 
        ReaK::graph::fixed_neighborhood< linear_pred_succ_search<> > >(
          &(this->m_space->get_super_space()), 
          get_sample, 
          ReaK::graph::fixed_neighborhood< linear_pred_succ_search<> >(
            linear_pred_succ_search<>(), 
            5, std::numeric_limits<double>::infinity())),
//       get(distance_metric, this->m_space->get_super_space()),
      ReaK::graph::star_neighborhood< linear_pred_succ_search<> >(
        linear_pred_succ_search<>(), 
        1, 3.0 * space_Lc), 
//       ReaK::graph::fixed_neighborhood< linear_pred_succ_search<> >(
//         linear_pred_succ_search<>(), 
//         5, std::numeric_limits<double>::infinity()),
      this->m_max_vertex_count);
    
  } else if(this->m_knn_flag == DVP_BF2_TREE_KNN) {
    
    typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 2, 
                     random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2>, 
                     no_position_caching_policy > SpacePartType;
    SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
    
    multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> nn_finder;
    nn_finder.graph_tree_map[&motion_graph] = &space_part;
    
    MEAQR_rrtstar_visitor<StateSpace, StateSpaceSystem, StateSpaceSampler, multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
    
    
    ReaK::graph::detail::generate_rrt_star_loop(
      motion_graph, 
      this->m_space->get_super_space(),
      vis, 
      pos_map, 
      cost_map,
      pred_map,
      weight_map,
      ReaK::graph::rrg_node_generator<
        SuperSpace, 
        SamplerType, 
        ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> > >(
          &(this->m_space->get_super_space()), 
          get_sample, 
          ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> >(
            nn_finder, 
            5, std::numeric_limits<double>::infinity())),
//       get(distance_metric, this->m_space->get_super_space()),
      ReaK::graph::star_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> >(
        nn_finder, 
        1, 3.0 * space_Lc),
//       ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> >(
//         nn_finder, 
//         5, std::numeric_limits<double>::infinity()),
      this->m_max_vertex_count);
    
  } else if(this->m_knn_flag == DVP_BF4_TREE_KNN) {
    
    typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 4, 
                     random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4>, 
                     no_position_caching_policy > SpacePartType;
    SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
    
    multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> nn_finder;
    nn_finder.graph_tree_map[&motion_graph] = &space_part;
    
    MEAQR_rrtstar_visitor<StateSpace, StateSpaceSystem, StateSpaceSampler, multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
    
    ReaK::graph::detail::generate_rrt_star_loop(
      motion_graph, 
      this->m_space->get_super_space(),
      vis, 
      pos_map, 
      cost_map,
      pred_map,
      weight_map,
      ReaK::graph::rrg_node_generator<
        SuperSpace, 
        SamplerType, 
        ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> > >(
          &(this->m_space->get_super_space()), 
          get_sample, 
          ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> >(
            nn_finder, 
            5, std::numeric_limits<double>::infinity())),
//       get(distance_metric, this->m_space->get_super_space()),
      ReaK::graph::star_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> >(
        nn_finder, 
        1, 3.0 * space_Lc), 
//       ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> >(
//         nn_finder, 
//         5, std::numeric_limits<double>::infinity()),
      this->m_max_vertex_count);
    
  } else if(this->m_knn_flag == DVP_COB2_TREE_KNN) {
    
    typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 2, 
                     random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2>, 
                     no_position_caching_policy > SpacePartType;
    SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
    
    multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> nn_finder;
    nn_finder.graph_tree_map[&motion_graph] = &space_part;
    
    MEAQR_rrtstar_visitor<StateSpace, StateSpaceSystem, StateSpaceSampler, multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
    
    
    ReaK::graph::detail::generate_rrt_star_loop(
      motion_graph, 
      this->m_space->get_super_space(),
      vis, 
      pos_map, 
      cost_map,
      pred_map,
      weight_map,
      ReaK::graph::rrg_node_generator<
        SuperSpace, 
        SamplerType, 
        ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> > >(
          &(this->m_space->get_super_space()), 
          get_sample, 
          ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> >(
            nn_finder, 
            5, std::numeric_limits<double>::infinity())),
//       get(distance_metric, this->m_space->get_super_space()),
      ReaK::graph::star_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> >(
        nn_finder, 
        1, 3.0 * space_Lc), 
//       ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> >(
//         nn_finder, 
//         5, std::numeric_limits<double>::infinity()),
      this->m_max_vertex_count);
    
  } else if(this->m_knn_flag == DVP_COB4_TREE_KNN) {
    
    
    typedef dvp_tree<Vertex, SuperSpace, GraphPositionMap, 4, 
                     random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4>, 
                     no_position_caching_policy > SpacePartType;
    SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), g_pos_map);
    
    multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> nn_finder;
    nn_finder.graph_tree_map[&motion_graph] = &space_part;
    
    MEAQR_rrtstar_visitor<StateSpace, StateSpaceSystem, StateSpaceSampler, multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
    
    
    ReaK::graph::detail::generate_rrt_star_loop(
      motion_graph, 
      this->m_space->get_super_space(),
      vis, 
      pos_map, 
      cost_map,
      pred_map,
      weight_map,
      ReaK::graph::rrg_node_generator<
        SuperSpace, 
        SamplerType, 
        ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> > >(
          &(this->m_space->get_super_space()), 
          get_sample, 
          ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> >(
            nn_finder, 
            5, std::numeric_limits<double>::infinity())),
//       get(distance_metric, this->m_space->get_super_space()),
      ReaK::graph::star_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> >(
        nn_finder, 
        1, 3.0 * space_Lc), 
//       ReaK::graph::fixed_neighborhood< multi_dvp_tree_pred_succ_search<MotionGraphType, SpacePartType> >(
//         nn_finder, 
//         5, std::numeric_limits<double>::infinity()),
      this->m_max_vertex_count);
    
  };
  
  if(this->m_solutions.size())
    return this->m_solutions.begin()->second;
  else
    return shared_ptr< seq_path_base< SuperSpace > >();
};


};

};

#endif


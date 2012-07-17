/**
 * \file rrt_path_planner.hpp
 * 
 * This library defines a class
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
#include "sbmp_reporter_concept.hpp"

#include "metric_space_concept.hpp"
#include "path_base.hpp"
#include "interpolation/point_to_point_path.hpp"
#include "basic_sbmp_reporters.hpp"

#include "graph_alg/rr_tree.hpp"

namespace ReaK {
  
namespace pp {
  

const std::size_t UNIDIRECTIONAL_RRT = 0;
const std::size_t BIDIRECTIONAL_RRT = 1;


template <typename FreeSpaceType>
struct rrt_vertex_data {
  typename topology_traits<FreeSpaceType>::point_type position;
  double distance_accum;
};

template <typename FreeSpaceType>
struct rrt_edge_data { };





/**
 * This class is a RRT-based path-planner over the given topology.
 * \tparam FreeSpaceType The topology type on which to perform the planning, should be the C-free sub-space of a larger configuration space.
 * \tparam SBPPReporter The reporter type to use to report the progress of the path-planning.
 */
template <typename FreeSpaceType, 
          typename SBPPReporter = no_sbmp_report>
class rrt_path_planner : public sample_based_planner< path_planner_base<FreeSpaceType> > {
  public:
    typedef sample_based_planner< path_planner_base<FreeSpaceType> > base_type;
    typedef rrt_path_planner<FreeSpaceType, SBPPReporter> self;
    
    typedef FreeSpaceType space_type;
    typedef subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    
    BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));
    
    typedef topology_traits< super_space_type >::point_type point_type;
    typedef topology_traits< super_space_type >::point_difference_type point_difference_type;
    
  protected:
    SBPPReporter m_reporter;
    point_type m_start_pos;
    point_type m_goal_pos;
    std::size_t max_num_results;
    bool has_reached_max_vertices;
    std::size_t m_bidir_flag;
    std::size_t m_graph_kind_flag;
    std::size_t m_knn_flag;
    
    std::map<double, shared_ptr< path_base< super_space_type > > > m_solutions;
    
  public:
    
    template <typename Vertex, typename Graph>
    void check_goal_connection(const point_type& p_u, Vertex u, Graph& g) {
      point_type result_p = this->m_space->move_position_toward(p_u, 1.0, m_goal_pos);
      double best_case_dist = get(distance_metric, this->m_space->get_super_space())(p_u, m_goal_pos, this->m_space->get_super_space());
      double actual_dist = get(distance_metric, this->m_space->get_super_space())(p_u, result_p, this->m_space->get_super_space());
      
      if(actual_dist < 0.99 * best_case_dist)
        return;
      
      double solutions_total_dist = actual_dist + g[u].distance_accum;
      if(solutions_total_dist >= m_solutions.begin()->first)
        return;
      
      shared_ptr< path_wrapper< point_to_point_path<FreeSpaceType> > > new_sol(new path_wrapper< point_to_point_path<FreeSpaceType> >("rrt_solution"));
      point_to_point_path<FreeSpaceType>& waypoints = new_sol->get_underlying_path();
      
      waypoints.push_front(goal_pos);
      waypoints.push_front(p_u);
      
      while(in_degree(u, g)) {
        u = source(*(in_edges(u,g).first),g);
        waypoints.push_front(g[u].position);
      };
      
      m_solutions[solutions_total_dist] = new_sol;
      m_reporter.draw_solution(*(this->m_space), m_solutions[solutions_total_dist]);
    };
    
    template <typename Vertex, typename Graph>
    void joining_vertex_found(Vertex u1, Vertex u2, Graph& g1, Graph& g2) {
      double total_dist = g1[u1].distance_accum + g2[u2].distance_accum
        + get(distance_metric, this->m_space->get_super_space())(g1[u1].position, g2[u2].position, this->m_space->get_super_space());
      
      if(total_dist >= m_solutions.begin()->first)
        return;
      
      shared_ptr< path_wrapper< point_to_point_path<FreeSpaceType> > > new_sol(new path_wrapper< point_to_point_path<FreeSpaceType> >("birrt_solution"));
      point_to_point_path<FreeSpaceType>& waypoints = new_sol->get_underlying_path();
      
      waypoints.push_front(g1[u1].position);
      while(in_degree(u1, g1)) {
        u1 = source(*(in_edges(u1,g1).first),g1);
        waypoints.push_front(g1[u1].position);
      };
      
      waypoints.push_back(g2[u2].position);
      while(in_degree(u2, g2)) {
        u2 = source(*(in_edges(u2,g2).first),g2);
        waypoints.push_back(g2[u2].position);
      };
      
      m_solutions[total_dist] = new_sol;
      m_reporter.draw_solution(*(this->m_space), m_solutions[total_dist]);
    };
    
    bool keep_going() const {
      return (max_num_results > m_solutions.size()) && !has_reached_max_vertices;
    };
    
    template <typename Graph>
    void report_progress(Graph& g) {
      if(num_vertices(g) % this->m_progress_interval == 0)
        m_reporter.draw_motion_graph(*(this->m_space), g, get(&rrt_vertex_data<FreeSpaceType>::position,g));
      has_reached_max_vertices = (num_vertices(g) >= this->m_max_vertex_count);
    };
    
    template <typename Graph>
    void report_progress(Graph& g1, Graph& g2) {
      if(num_vertices(g) % this->m_progress_interval == 0) {
        m_reporter.draw_motion_graph(*(this->m_space), g1, get(&rrt_vertex_data<FreeSpaceType>::position,g1));
        m_reporter.draw_motion_graph(*(this->m_space), g2, get(&rrt_vertex_data<FreeSpaceType>::position,g2));
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
    virtual shared_ptr< path_base< super_space_type > > solve_motion();
    
    /**
     * Parametrized constructor.
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     * \param aMaxVertexCount The maximum number of samples to generate during the motion planning.
     * \param aProgressInterval The number of new samples between each "progress report".
     */
    rrt_path_planner(const shared_ptr< space_type >& aWorld, 
                     std::size_t aMaxVertexCount = 5000, 
                     std::size_t aProgressInterval = 100) :
                     base_type("rrt_planner", aWorld, aMaxVertexCount, aProgressInterval) { };
    
    virtual ~rrt_path_planner() { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
//       A & RK_SERIAL_SAVE_WITH_NAME(m_max_vertex_count);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
//       A & RK_SERIAL_LOAD_WITH_NAME(m_max_vertex_count);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2460007,1,"rrt_path_planner",base_type)
};




template <typename FreeSpaceType, typename NNFinderSynchro, typename SBPPReporter = no_sbmp_report>
struct rrt_planner_visitor {
  shared_ptr< FreeSpaceType > m_space;
  rrt_path_planner<FreeSpaceType,SBPPReporter>* m_planner;
  NNFinderSynchro m_nn_synchro;
  
  rrt_planner_visitor(const shared_ptr< FreeSpaceType >& aSpace, 
                      rrt_path_planner<FreeSpaceType,SBPPReporter>* aPlanner,
                      NNFinderSynchro aNNSynchro) : 
                      m_space(aSpace), m_planner(aPlanner), m_nn_synchro(aNNSynchro) { };
  
  typedef typename topology_traits<FreeSpaceType>::point_type PointType;
  
  template <typename Vertex, typename Graph>
  void vertex_added(Vertex u, Graph& g) const {
    m_nn_synchro.added_vertex(u,g);
  };
  
  template <typename Vertex, typename Graph>
  void vertex_added(Vertex u, Graph& g, Graph&) const {
    m_nn_synchro.added_vertex(u,g);
  };
  
  template <typename EdgeType, typename Graph>
  void edge_added(EdgeType e, Graph& g) const {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor VertexType;
    VertexType u = source(e,g);
    VertexType v = target(e,g);
    
    g[v].distance_accum = g[u].distance_accum + get(distance_metric, m_space->get_super_space())(g[u].position, g[v].position, m_space->get_super_space())
    
    // Call progress reporter...
    m_planner->report_progress(g);
    
    // Check if a straight path to goal is possible...
    m_planner->check_goal_connection(g[v].position);
    
  };
  
  template <typename EdgeType, typename Graph>
  void edge_added(EdgeType e, Graph& g1, Graph& g2) const {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor VertexType;
    VertexType u = source(e,g1);
    VertexType v = target(e,g1);
    
    g1[v].distance_accum = g1[u].distance_accum + get(distance_metric, m_space->get_super_space())(g1[u].position, g1[v].position, m_space->get_super_space())
    
    // Call progress reporter...
    m_planner->report_progress(g1,g2);
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
  std::pair<PointType, bool> steer_towards_position(const PointType& p, Vertex u, Graph& g) const {
    PointType result_p = m_space->move_position_toward(g[u].position, 1.0, p);
    double best_case_dist = get(distance_metric, m_space->get_super_space())(g[u].position, p, m_space->get_super_space());
    double actual_dist = get(distance_metric, m_space->get_super_space())(g[u].position, result_p, m_space->get_super_space());
    if(actual_dist > 0.1 * best_case_dist)
      return std::make_pair(result_p, true);
    else
      return std::make_pair(result_p, false);
  };
  
  
  
};






template <typename FreeSpaceType, 
          typename SBPPReporter>
shared_ptr< path_base< rrt_path_planner<FreeSpaceType,SBPPReporter>::super_space_type > > 
  rrt_path_planner<FreeSpaceType,SBPPReporter>::solve_motion() {
  
  typedef typename subspace_traits<FreeSpaceType>::super_space_type SuperSpace;
  typedef typename topology_traits<SuperSpace>::point_type PointType;
  
  if(m_bidir_flag == UNIDIRECTIONAL_RRT) {
    
    if(m_graph_kind_flag == ADJ_LIST_MOTION_GRAPH) {
      
      typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::bidirectionalS,
                             rrt_vertex_data<FreeSpaceType>,
                             rrt_edge_data<FreeSpaceType>,
                             boost::vecS> MotionGraphType;
      typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex;
      typedef typename boost::property_map<MotionGraphType, &rrt_vertex_data<FreeSpaceType>::position >::type PositionMap;
      
      MotionGraphType motion_graph;
      Vertex v = add_vertex(motion_graph);
      motion_graph[v].position = this->start_pos;
      motion_graph[v].distance_accum = 0.0;
      PositionMap pos_map = get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph);
      
      if(m_knn_flag == LINEAR_SEARCH_KNN) {
        rrt_planner_visitor<FreeSpaceType, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro());
        
        ReaK::graph::generate_rrt(motion_graph, this->m_space->get_super_space(),
                                  vis, get(random_sampler, this->m_space->get_super_space()), 
                                  linear_neighbor_search<>(), this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_BF2_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, PositionMap, 2, 
                         random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2> > SpacePartType;
        SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
        
        multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_rrt(motion_graph, this->m_space->get_super_space(),
                                  vis, get(random_sampler, this->m_space->get_super_space()), 
                                  pos_map, nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_BF4_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, PositionMap, 4, 
                         random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4> > SpacePartType;
        SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
        
        multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_rrt(motion_graph, this->m_space->get_super_space(),
                                  vis, get(random_sampler, this->m_space->get_super_space()), 
                                  pos_map, nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_COB2_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, PositionMap, 2, 
                         random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2> > SpacePartType;
        SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
        
        multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_rrt(motion_graph, this->m_space->get_super_space(),
                                  vis, get(random_sampler, this->m_space->get_super_space()), 
                                  pos_map, nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_COB4_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, PositionMap, 4, 
                         random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4> > SpacePartType;
        SpacePartType space_part(motion_graph, ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), pos_map);
        
        multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_rrt(motion_graph, this->m_space->get_super_space(),
                                  vis, get(random_sampler, this->m_space->get_super_space()), 
                                  pos_map, nn_finder, this->m_max_vertex_count);
        
      };
      
    } else if(m_graph_kind_flag == DVP_ADJ_LIST_MOTION_GRAPH) {
      
      if(m_knn_flag == DVP_ALT_BF2_KNN) {
        
        typedef dvp_adjacency_list<
          rrt_vertex_data<FreeSpaceType>,
          rrt_edge_data<FreeSpaceType>,
          SuperSpace,
          data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >,
          2, random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2> > ALTGraph;
        
        ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                            data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >(&rrt_vertex_data<FreeSpaceType>::position));
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        
        MotionGraph motion_graph = space_part.get_adjacency_list();
        
        rrt_vertex_data<FreeSpaceType> v_p;
        v_p.position = this->start_pos;
        v_p.distance_accum = 0.0;
        add_vertex(v_p, motion_graph);
        
        multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraph, ALTGraph>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_rrt(motion_graph, this->m_space->get_super_space(),
                                  vis, get(random_sampler, this->m_space->get_super_space()), 
                                  get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph), 
                                  nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_ALT_BF4_KNN) {
        
        typedef dvp_adjacency_list<
          rrt_vertex_data<FreeSpaceType>,
          rrt_edge_data<FreeSpaceType>,
          SuperSpace,
          data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >,
          4, random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4> > ALTGraph;
        
        ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                            data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >(&rrt_vertex_data<FreeSpaceType>::position));
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        
        MotionGraph motion_graph = space_part.get_adjacency_list();
        
        rrt_vertex_data<FreeSpaceType> v_p;
        v_p.position = this->start_pos;
        v_p.distance_accum = 0.0;
        add_vertex(v_p, motion_graph);
        
        multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraph, ALTGraph>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_rrt(motion_graph, this->m_space->get_super_space(),
                                  vis, get(random_sampler, this->m_space->get_super_space()), 
                                  get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph), 
                                  nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_ALT_COB2_KNN) {
        
        typedef dvp_adjacency_list<
          rrt_vertex_data<FreeSpaceType>,
          rrt_edge_data<FreeSpaceType>,
          SuperSpace,
          data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >,
          2, random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2> > ALTGraph;
        
        ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                            data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >(&rrt_vertex_data<FreeSpaceType>::position));
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        
        MotionGraph motion_graph = space_part.get_adjacency_list();
        
        rrt_vertex_data<FreeSpaceType> v_p;
        v_p.position = this->start_pos;
        v_p.distance_accum = 0.0;
        add_vertex(v_p, motion_graph);
        
        multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraph, ALTGraph>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_rrt(motion_graph, this->m_space->get_super_space(),
                                  vis, get(random_sampler, this->m_space->get_super_space()), 
                                  get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph), 
                                  nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_ALT_COB4_KNN) {
        
        typedef dvp_adjacency_list<
          rrt_vertex_data<FreeSpaceType>,
          rrt_edge_data<FreeSpaceType>,
          SuperSpace,
          data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >,
          4, random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4> > ALTGraph;
        
        ALTGraph space_part(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                            data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >(&rrt_vertex_data<FreeSpaceType>::position));
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        
        MotionGraph motion_graph = space_part.get_adjacency_list();
        
        rrt_vertex_data<FreeSpaceType> v_p;
        v_p.position = this->start_pos;
        v_p.distance_accum = 0.0;
        add_vertex(v_p, motion_graph);
        
        multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
        nn_finder.graph_tree_map[&motion_graph] = &space_part;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraph, ALTGraph>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_rrt(motion_graph, this->m_space->get_super_space(),
                                  vis, get(random_sampler, this->m_space->get_super_space()), 
                                  get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph), 
                                  nn_finder, this->m_max_vertex_count);
        
      };
      
    };
    
  } else {
    
    if(m_graph_kind_flag == ADJ_LIST_MOTION_GRAPH) {
      
      typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::bidirectionalS,
                             rrt_vertex_data<FreeSpaceType>,
                             rrt_edge_data<FreeSpaceType>,
                             boost::vecS> MotionGraphType;
      typedef typename boost::graph_traits<MotionGraphType>::vertex_descriptor Vertex;
      typedef typename boost::property_map<MotionGraphType, &rrt_vertex_data<FreeSpaceType>::position >::type PositionMap;
      
      MotionGraphType motion_graph1;
      Vertex v1 = add_vertex(motion_graph1);
      motion_graph1[v1].position = this->start_pos;
      motion_graph1[v1].distance_accum = 0.0;
      PositionMap pos_map1 = get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph1);
      
      MotionGraphType motion_graph2;
      Vertex v2 = add_vertex(motion_graph2);
      motion_graph2[v2].position = this->goal_pos;
      motion_graph2[v2].distance_accum = 0.0;
      PositionMap pos_map2 = get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph2);
      
      if(m_knn_flag == LINEAR_SEARCH_KNN) {
        rrt_planner_visitor<FreeSpaceType, no_NNfinder_synchro, SBPPReporter> vis(this->m_space, this, no_NNfinder_synchro());
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, get(random_sampler, this->m_space->get_super_space()), 
          pos_map1, pos_map2, linear_neighbor_search<>(), this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_BF2_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, PositionMap, 2, 
                         random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2> > SpacePartType;
        SpacePartType space_part1(motion_graph1, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  pos_map1);
        SpacePartType space_part2(motion_graph2, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  pos_map2);
        
        multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, get(random_sampler, this->m_space->get_super_space()), 
          pos_map1, pos_map2, nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_BF4_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, PositionMap, 4, 
                         random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4> > SpacePartType;
        SpacePartType space_part1(motion_graph1, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  pos_map1);
        SpacePartType space_part2(motion_graph2, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  pos_map2);
        
        multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, get(random_sampler, this->m_space->get_super_space()), 
          pos_map1, pos_map2, nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_COB2_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, PositionMap, 2, 
                         random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2> > SpacePartType;
        SpacePartType space_part1(motion_graph1, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  pos_map1);
        SpacePartType space_part2(motion_graph2, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  pos_map2);
        
        multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, get(random_sampler, this->m_space->get_super_space()), 
          pos_map1, pos_map2, nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_COB4_TREE_KNN) {
        
        typedef dvp_tree<Vertex, SuperSpace, PositionMap, 4, 
                         random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4> > SpacePartType;
        SpacePartType space_part1(motion_graph1, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  pos_map1);
        SpacePartType space_part2(motion_graph2, 
                                  ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                                  pos_map2);
        
        multi_dvp_tree_search<MotionGraphType, SpacePartType> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraphType, SpacePartType>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, get(random_sampler, this->m_space->get_super_space()), 
          pos_map1, pos_map2, nn_finder, this->m_max_vertex_count);
        
      };
      
    } else if(m_graph_kind_flag == DVP_ADJ_LIST_MOTION_GRAPH) {
      
      if(m_knn_flag == DVP_ALT_BF2_KNN) {
        
        typedef dvp_adjacency_list<
          rrt_vertex_data<FreeSpaceType>,
          rrt_edge_data<FreeSpaceType>,
          SuperSpace,
          data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >,
          2, random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<2> > ALTGraph;
        
        ALTGraph space_part1(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                             data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >(&rrt_vertex_data<FreeSpaceType>::position));
        
        ALTGraph space_part2(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                             data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >(&rrt_vertex_data<FreeSpaceType>::position));
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        
        MotionGraph motion_graph1 = space_part1.get_adjacency_list();
        rrt_vertex_data<FreeSpaceType> v1_p;
        v1_p.position = this->start_pos;
        v1_p.distance_accum = 0.0;
        add_vertex(v1_p, motion_graph1);
        
        MotionGraph motion_graph2 = space_part2.get_adjacency_list();
        rrt_vertex_data<FreeSpaceType> v2_p;
        v2_p.position = this->goal_pos;
        v2_p.distance_accum = 0.0;
        add_vertex(v2_p, motion_graph2);
        
        multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraph, ALTGraph>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, get(random_sampler, this->m_space->get_super_space()), 
          get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph1), 
          get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph2), 
          nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_ALT_BF4_KNN) {
        
        typedef dvp_adjacency_list<
          rrt_vertex_data<FreeSpaceType>,
          rrt_edge_data<FreeSpaceType>,
          SuperSpace,
          data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >,
          4, random_vp_chooser, ReaK::graph::d_ary_bf_tree_storage<4> > ALTGraph;
        
        ALTGraph space_part1(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                             data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >(&rrt_vertex_data<FreeSpaceType>::position));
        
        ALTGraph space_part2(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                             data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >(&rrt_vertex_data<FreeSpaceType>::position));
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        
        MotionGraph motion_graph1 = space_part1.get_adjacency_list();
        rrt_vertex_data<FreeSpaceType> v1_p;
        v1_p.position = this->start_pos;
        v1_p.distance_accum = 0.0;
        add_vertex(v1_p, motion_graph1);
        
        MotionGraph motion_graph2 = space_part2.get_adjacency_list();
        rrt_vertex_data<FreeSpaceType> v2_p;
        v2_p.position = this->goal_pos;
        v2_p.distance_accum = 0.0;
        add_vertex(v2_p, motion_graph2);
        
        multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraph, ALTGraph>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, get(random_sampler, this->m_space->get_super_space()), 
          get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph1), 
          get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph2), 
          nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_ALT_COB2_KNN) {
        
        typedef dvp_adjacency_list<
          rrt_vertex_data<FreeSpaceType>,
          rrt_edge_data<FreeSpaceType>,
          SuperSpace,
          data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >,
          2, random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<2> > ALTGraph;
        
        ALTGraph space_part1(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                             data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >(&rrt_vertex_data<FreeSpaceType>::position));
        
        ALTGraph space_part2(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                             data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >(&rrt_vertex_data<FreeSpaceType>::position));
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        
        MotionGraph motion_graph1 = space_part1.get_adjacency_list();
        rrt_vertex_data<FreeSpaceType> v1_p;
        v1_p.position = this->start_pos;
        v1_p.distance_accum = 0.0;
        add_vertex(v1_p, motion_graph1);
        
        MotionGraph motion_graph2 = space_part2.get_adjacency_list();
        rrt_vertex_data<FreeSpaceType> v2_p;
        v2_p.position = this->goal_pos;
        v2_p.distance_accum = 0.0;
        add_vertex(v2_p, motion_graph2);
        
        multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraph, ALTGraph>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, get(random_sampler, this->m_space->get_super_space()), 
          get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph1), 
          get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph2), 
          nn_finder, this->m_max_vertex_count);
        
      } else if(m_knn_flag == DVP_ALT_COB4_KNN) {
        
        typedef dvp_adjacency_list<
          rrt_vertex_data<FreeSpaceType>,
          rrt_edge_data<FreeSpaceType>,
          SuperSpace,
          data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >,
          4, random_vp_chooser, ReaK::graph::d_ary_cob_tree_storage<4> > ALTGraph;
        
        ALTGraph space_part1(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                             data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >(&rrt_vertex_data<FreeSpaceType>::position));
        
        ALTGraph space_part2(ReaK::shared_ptr<const SuperSpace>(&(this->m_space->get_super_space()),null_deleter()), 
                             data_member_property_map<PointType, rrt_vertex_data<FreeSpaceType> >(&rrt_vertex_data<FreeSpaceType>::position));
        
        typedef typename ALTGraph::adj_list_type MotionGraph;
        
        MotionGraph motion_graph1 = space_part1.get_adjacency_list();
        rrt_vertex_data<FreeSpaceType> v1_p;
        v1_p.position = this->start_pos;
        v1_p.distance_accum = 0.0;
        add_vertex(v1_p, motion_graph1);
        
        MotionGraph motion_graph2 = space_part2.get_adjacency_list();
        rrt_vertex_data<FreeSpaceType> v2_p;
        v2_p.position = this->goal_pos;
        v2_p.distance_accum = 0.0;
        add_vertex(v2_p, motion_graph2);
        
        multi_dvp_tree_search<MotionGraph, ALTGraph> nn_finder;
        nn_finder.graph_tree_map[&motion_graph1] = &space_part1;
        nn_finder.graph_tree_map[&motion_graph2] = &space_part2;
        
        rrt_planner_visitor<FreeSpaceType, multi_dvp_tree_search<MotionGraph, ALTGraph>, SBPPReporter> vis(this->m_space, this, nn_finder);
        
        ReaK::graph::generate_bidirectional_rrt(
          motion_graph1, motion_graph2, this->m_space->get_super_space(),
          vis, get(random_sampler, this->m_space->get_super_space()), 
          get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph1), 
          get(&rrt_vertex_data<FreeSpaceType>::position, motion_graph2), 
          nn_finder, this->m_max_vertex_count);
        
      };
      
    };
  };
  
  if(m_solutions.size())
    return m_solutions.begin()->second;
  else
    return shared_ptr< path_base< rrt_path_planner<FreeSpaceType,SBPPReporter>::super_space_type > >();
};


};

};

#endif


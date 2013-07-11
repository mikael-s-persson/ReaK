/**
 * \file planning_queries.hpp
 * 
 * This library defines class templates to encode a path-planning or motion-planning query as the 
 * contract to be fulfilled by path-planners used in ReaK. 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2013
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

#ifndef REAK_PLANNING_QUERIES_HPP
#define REAK_PLANNING_QUERIES_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "metric_space_concept.hpp"
#include "steerable_space_concept.hpp"
#include "random_sampler_concept.hpp"
#include "subspace_concept.hpp"

#include "trajectory_base.hpp"
#include "seq_path_base.hpp"

#include "any_motion_graphs.hpp"

#include <map>

namespace ReaK {
  
namespace pp {

  
#if 0
/**
 * This class is the basic OOP interface for a motion planner. 
 * OOP-style planners are useful to hide away
 * the cumbersome details of calling the underlying planning algorithms which are 
 * generic programming (GP) style and thus provide a lot more flexibility but are difficult
 * to deal with in the user-space. The OOP planners are meant to offer a much simpler interface,
 * i.e., a member function that "solves the problem" and returns the solution path or trajectory.
 */
template <typename FreeSpaceType>
class motion_planner_base : public named_object {
  public:
    typedef motion_planner_base<FreeSpaceType> self;
    typedef FreeSpaceType space_type;
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    
    BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));
    
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename topology_traits< super_space_type >::point_difference_type point_difference_type;
    
    
  public:
    
    shared_ptr< space_type > space;
    
    
    /**
     * This function computes a valid trajectory in the C-free. If it cannot 
     * achieve a valid trajectory, an exception will be thrown. This algorithmic
     * motion solver class is such that any settings that ought to be set for the 
     * motion planning algorithm should be set before calling this function, otherwise
     * the function is likely to fail.
     * \return The trajectory object that can be used to map out the motion in time.
     */
    virtual shared_ptr< trajectory_base< super_space_type > > solve_motion() = 0;
    
    /**
     * This function is called to reset the internal state of the planner.
     */
    virtual void reset_internal_state() = 0;
    
    /**
     * Parametrized constructor.
     * \param aName The name for this object.
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     */
    motion_planner_base(const std::string& aName,
                        const shared_ptr< space_type >& aWorld) :
                        named_object(),
                        m_space(aWorld) { setName(aName); };
    
    virtual ~motion_planner_base() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_space);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_space);
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2460000,1,"motion_planner_base",named_object)
};
#endif



namespace detail {
  
  template <typename FreeSpaceType, typename Graph>
  typename boost::enable_if< is_steerable_space< FreeSpaceType >,
  bool >::type register_solution_path_impl(const FreeSpaceType& space, const Graph& g, 
                                           typename boost::graph_traits< Graph >::vertex_descriptor start_node, 
                                           typename boost::graph_traits< Graph >::vertex_descriptor goal_node,
                                           const mg_vertex_data< FreeSpaceType >& goal_data,
                                           const typename topology_traits< FreeSpaceType >::point_type& goal_pos,
                                           std::map<double, shared_ptr< seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type > > >& solutions) {
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename boost::graph_traits< Graph >::edge_descriptor Edge;
    typedef typename steerable_space_traits< super_space_type >::steer_record_type SteerRecordType;
    typedef typename SteerRecordType::point_fraction_iterator SteerIter;
    
    shared_ptr< super_space_type > sup_space_ptr(&(space.get_super_space()),null_deleter());
    
    double solutions_total_dist = get(distance_metric, *sup_space_ptr)(goal_data.position, goal_pos, *sup_space_ptr);
    
    shared_ptr< seq_path_wrapper< discrete_point_path<super_space_type> > > new_sol(new seq_path_wrapper< discrete_point_path<super_space_type> >("planning_solution", discrete_point_path<super_space_type>(sup_space_ptr,get(distance_metric, *sup_space_ptr))));
    discrete_point_path<super_space_type>& waypoints = new_sol->get_underlying_path();
    std::set<Vertex> path;
    
    Vertex v = goal_node;
    waypoints.push_front(g[v].position);
    
    while((in_degree(v, g)) && (v != start_node) && (path.insert(v).second)) {
      Edge e = *(in_edges(v, g).first);
      Vertex u = source(e, g);
      for(SteerIter it = g[e].steer_record.end_fraction_travel(); it != g[e].steer_record.begin_fraction_travel(); it -= 1.0)
        waypoints.push_front( point_type(*it) );
      solutions_total_dist += get(distance_metric, *sup_space_ptr)(g[u].position, g[v].position, *sup_space_ptr);
      v = u;
      waypoints.push_front(g[v].position);
    };
    
    if((v == start_node) && ((solutions.empty()) || (solutions_total_dist < solutions.begin()->first))) {
      solutions[solutions_total_dist] = new_sol;
//       reporter.draw_solution(space, solutions[solutions_total_dist]);
      return true;
    };
    
    return false;
  };
  
  template <typename FreeSpaceType, typename Graph>
  typename boost::disable_if< is_steerable_space< FreeSpaceType >,
  bool >::type register_solution_path_impl(const FreeSpaceType& space, const Graph& g, 
                                           typename boost::graph_traits< Graph >::vertex_descriptor start_node, 
                                           typename boost::graph_traits< Graph >::vertex_descriptor goal_node,
                                           const mg_vertex_data< FreeSpaceType >& goal_data,
                                           const typename topology_traits< FreeSpaceType >::point_type& goal_pos,
                                           std::map<double, shared_ptr< seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type > > >& solutions) {
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    
    shared_ptr< super_space_type > sup_space_ptr(&(space.get_super_space()),null_deleter());
    
    double solutions_total_dist = get(distance_metric, *sup_space_ptr)(goal_data.position, goal_pos, *sup_space_ptr);
    
    shared_ptr< seq_path_wrapper< point_to_point_path<super_space_type> > > new_sol(new seq_path_wrapper< point_to_point_path<super_space_type> >("planning_solution", point_to_point_path<super_space_type>(sup_space_ptr,get(distance_metric, *sup_space_ptr))));
    point_to_point_path<super_space_type>& waypoints = new_sol->get_underlying_path();
    
    waypoints.push_front(goal_pos);
    waypoints.push_front(goal_data.position);
    
    while(in_degree(goal_node, g) && (goal_node != start_node)) {
      Vertex v = source(*(in_edges(goal_node, g).first), g);
      solutions_total_dist += get(distance_metric, *sup_space_ptr)(g[v].position, g[goal_node].position, *sup_space_ptr);
      waypoints.push_front(g[v].position);
      goal_node = v;
    };
    
    if((goal_node == start_node) && ((solutions.empty()) || (solutions_total_dist < solutions.begin()->first))) {
      solutions[solutions_total_dist] = new_sol;
//       reporter.draw_solution(space, solutions[solutions_total_dist]);
      return true;
    };
    
    return false;
  };
  
  
  
  template <typename FreeSpaceType, typename Graph>
  typename boost::enable_if< is_steerable_space< FreeSpaceType >,
  bool >::type register_solution_path_impl(const FreeSpaceType& space, const Graph& g, 
                                           typename boost::graph_traits< Graph >::vertex_descriptor start_node, 
                                           typename boost::graph_traits< Graph >::vertex_descriptor goal_node,
                                           const optimal_mg_vertex< FreeSpaceType >& goal_data,
                                           const typename topology_traits< FreeSpaceType >::point_type& goal_pos,
                                           std::map<double, shared_ptr< seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type > > >& solutions) {
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename boost::graph_traits< Graph >::in_edge_iterator InEdgeIter;
    typedef typename steerable_space_traits< super_space_type >::steer_record_type SteerRecordType;
    typedef typename SteerRecordType::point_fraction_iterator SteerIter;
    
    shared_ptr< super_space_type > sup_space_ptr(&(space.get_super_space()),null_deleter());
    
    if( ! (goal_data.distance_accum < std::numeric_limits<double>::infinity()) ||
        ( (!solutions.empty()) && (goal_data.distance_accum >= solutions.begin()->first) ) )
      return false;
    
    shared_ptr< seq_path_wrapper< discrete_point_path<super_space_type> > > new_sol(new seq_path_wrapper< discrete_point_path<super_space_type> >("planning_solution", discrete_point_path<super_space_type>(sup_space_ptr,get(distance_metric, *sup_space_ptr))));
    discrete_point_path<super_space_type>& waypoints = new_sol->get_underlying_path();
    std::set<Vertex> path;
    
    Vertex v = goal_node;
    waypoints.push_front(g[v].position);
    
    while((v != start_node) && (path.insert(v).second)) {
      Vertex u = g[v].predecessor;
      std::pair<InEdgeIter,InEdgeIter> er = in_edges(v, g);
      while( ( er.first != er.second ) && ( source(*(er.first), g) != u ) )
        ++(er.first);
      if(er.first == er.second)
        break;
      for(SteerIter it = g[*(er.first)].steer_record.end_fraction_travel(); it != g[*(er.first)].steer_record.begin_fraction_travel(); it -= 1.0)
        waypoints.push_front( point_type(*it) );
      v = u;
      waypoints.push_front(g[v].position);
    };
    
    if(v == start_node) {
      solutions[goal_data.distance_accum] = new_sol;
//       reporter.draw_solution(space, solutions[goal_data.distance_accum]);
      return true;
    };
    
    return false;
  };
  
  
  
  template <typename FreeSpaceType, typename Graph>
  typename boost::disable_if< is_steerable_space< FreeSpaceType >,
  bool >::type register_solution_path_impl(const FreeSpaceType& space, const Graph& g, 
                                           typename boost::graph_traits< Graph >::vertex_descriptor start_node, 
                                           typename boost::graph_traits< Graph >::vertex_descriptor goal_node,
                                           const optimal_mg_vertex< FreeSpaceType >& goal_data,
                                           const typename topology_traits< FreeSpaceType >::point_type& goal_pos,
                                           std::map<double, shared_ptr< seq_path_base< typename subspace_traits<FreeSpaceType>::super_space_type > > >& solutions) {
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    
    shared_ptr< super_space_type > sup_space_ptr(&(space.get_super_space()),null_deleter());
    
    if( ! (goal_data.distance_accum < std::numeric_limits<double>::infinity()) ||
        ( (!solutions.empty()) && (goal_data.distance_accum >= solutions.begin()->first) ) )
      return false;
    
    shared_ptr< seq_path_wrapper< point_to_point_path<super_space_type> > > new_sol(new seq_path_wrapper< point_to_point_path<super_space_type> >("planning_solution", point_to_point_path<super_space_type>(sup_space_ptr,get(distance_metric, *sup_space_ptr))));
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
      solutions[goal_data.distance_accum] = new_sol;
//       reporter.draw_solution(space, solutions[goal_data.distance_accum]);
      return true;
    };
    
    return false;
  };
  
};





/**
 * This class is the basic OOP interface for a path planner. 
 * OOP-style planners are useful to hide away
 * the cumbersome details of calling the underlying planning algorithms which are 
 * generic programming (GP) style and thus provide a lot more flexibility but are difficult
 * to deal with in the user-space. The OOP planners are meant to offer a much simpler interface,
 * i.e., a member function that "solves the problem" and returns the solution path or trajectory.
 */
template <typename FreeSpaceType>
class path_planning_query : public named_object {
  public:
    typedef path_planner_base<FreeSpaceType> self;
    typedef FreeSpaceType space_type;
    typedef typename subspace_traits<FreeSpaceType>::super_space_type super_space_type;
    
    BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));
    
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename topology_traits< super_space_type >::point_difference_type point_difference_type;
    
    
  public:
    
    shared_ptr< space_type > space;
    point_type start_pos;
    point_type goal_pos;
    std::size_t max_num_results;
    
    std::map<double, shared_ptr< seq_path_base< super_space_type > > > solutions;
    
    
    /**
     * Returns the best solution distance registered in this query object.
     * \return The best solution distance registered in this query object.
     */
    double get_best_solution_distance() const {
      if(solutions.size() == 0)
        return std::numeric_limits<double>::infinity();
      else
        return solutions.begin()->first;
    };
    
    /**
     * Returns true if the solver should keep on going trying to solve the path-planning problem.
     * \return True if the solver should keep on going trying to solve the path-planning problem.
     */
    bool keep_going() const {
      return (max_num_results > solutions.size());
    };
    
    
    /**
     * This function registers a solution path (if one is found) and invokes the path-planning 
     * reporter to report on that solution path.
     * \note This function is for internal use by the path-planning algorithm (a visitor callback).
     * \param start_node The start node in the motion-graph.
     * \param goal_node The goal node in the motion-graph.
     * \param g The current motion-graph.
     * \return True if a new solution was registered.
     */
    template <typename Vertex, typename Graph>
    bool register_solution_path(Vertex start_node, Vertex goal_node, const Graph& g) {
      return detail::register_solution_path_impl(*space, g, start_node, goal_node, g[goal_node], goal_pos, solutions);
    };
    
    
    /**
     * This function is called to reset the internal state of the planner.
     */
    void reset_solution_records() {
      solutions.clear();
    };
    
    /**
     * Parametrized constructor.
     * \param aName The name for this object.
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     */
    path_planner_base(const std::string& aName,
                        const shared_ptr< space_type >& aWorld) :
                        named_object(),
                        m_space(aWorld) { setName(aName); };
    
    virtual ~path_planner_base() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_space);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_space);
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2460001,1,"path_planner_base",named_object)
};


};

};

#endif


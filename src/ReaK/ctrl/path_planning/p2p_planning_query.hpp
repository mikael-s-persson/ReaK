/**
 * \file p2p_planning_query.hpp
 * 
 * This library defines class templates for a path-planning query to travel from a start point 
 * to a fixed goal point, i.e., a point-to-point query.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2013
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

#ifndef REAK_P2P_PLANNING_QUERY_HPP
#define REAK_P2P_PLANNING_QUERY_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/named_object.hpp>

#include <ReaK/ctrl/topologies/metric_space_concept.hpp>
#include <ReaK/ctrl/topologies/subspace_concept.hpp>
#include <ReaK/ctrl/topologies/steerable_space_concept.hpp>
#include <ReaK/ctrl/topologies/random_sampler_concept.hpp>

#include "planning_queries.hpp"
#include <ReaK/ctrl/interpolation/seq_path_wrapper.hpp>
#include <ReaK/ctrl/interpolation/point_to_point_path.hpp>
#include <ReaK/ctrl/interpolation/discrete_point_path.hpp>
#include <ReaK/ctrl/interpolation/seq_trajectory_wrapper.hpp>
#include <ReaK/ctrl/interpolation/point_to_point_trajectory.hpp>
#include <ReaK/ctrl/interpolation/discrete_point_trajectory.hpp>

#include "any_motion_graphs.hpp"

#include "solution_path_factories.hpp"

#include <boost/mpl/if.hpp>

#include <map>

namespace ReaK {
  
namespace pp {


/**
 * This class is the basic OOP interface for a path planner. 
 * OOP-style planners are useful to hide away
 * the cumbersome details of calling the underlying planning algorithms which are 
 * generic programming (GP) style and thus provide a lot more flexibility but are difficult
 * to deal with in the user-space. The OOP planners are meant to offer a much simpler interface,
 * i.e., a member function that "solves the problem" and returns the solution path or trajectory.
 */
template <typename FreeSpaceType>
class path_planning_p2p_query : public planning_query<FreeSpaceType> {
  public:
    typedef path_planning_p2p_query<FreeSpaceType> self;
    typedef planning_query<FreeSpaceType> base_type;
    typedef typename base_type::space_type space_type;
    typedef typename base_type::super_space_type super_space_type;
    
    typedef typename base_type::point_type point_type;
    typedef typename base_type::point_difference_type point_difference_type;
    
    typedef typename base_type::solution_record_ptr solution_record_ptr;
    
    typedef typename boost::mpl::if_< 
      is_steerable_space< space_type >,
      typename boost::mpl::if_< 
        is_temporal_space<space_type>,
        seq_trajectory_wrapper< discrete_point_trajectory<super_space_type> >,
        seq_path_wrapper< discrete_point_path<super_space_type> > >::type,
      typename boost::mpl::if_< 
        is_temporal_space<space_type>,
        seq_trajectory_wrapper< point_to_point_trajectory<super_space_type> >,
        seq_path_wrapper< point_to_point_path<super_space_type> > >::type >::type solution_path_wrapper;
    
  public:
    
    point_type start_pos;
    point_type goal_pos;
    std::size_t max_num_results;
    
    std::map<double, solution_record_ptr > solutions;
    
    
    /**
     * Returns the best solution distance registered in this query object.
     * \return The best solution distance registered in this query object.
     */
    virtual double get_best_solution_distance() const {
      if(solutions.size() == 0)
        return std::numeric_limits<double>::infinity();
      else
        return solutions.begin()->first;
    };
    
    /**
     * Returns true if the solver should keep on going trying to solve the path-planning problem.
     * \return True if the solver should keep on going trying to solve the path-planning problem.
     */
    virtual bool keep_going() const {
      return (max_num_results > solutions.size());
    };
    
    /**
     * This function is called to reset the internal state of the planner.
     */
    virtual void reset_solution_records() {
      solutions.clear();
    };
    
    virtual const point_type& get_start_position() const { return start_pos; };
    
    virtual double get_distance_to_goal(const point_type& pos) { 
      return get(distance_metric, *(this->space))(pos, goal_pos, *(this->space));
    };
    
    virtual double get_heuristic_to_goal(const point_type& pos) { 
      return get(distance_metric, this->space->get_super_space())(pos, goal_pos, this->space->get_super_space());
    };
    
  protected:
    
    virtual solution_record_ptr 
      register_solution_from_optimal_mg(graph::any_graph::vertex_descriptor start_node, 
                                        graph::any_graph::vertex_descriptor goal_node, 
                                        double goal_distance,
                                        graph::any_graph& g) {
      return detail::register_optimal_solution_path_impl<solution_path_wrapper>(*(this->space), g, start_node, goal_node, goal_pos, goal_distance, solutions);
    };
    
    virtual solution_record_ptr 
      register_solution_from_basic_mg(graph::any_graph::vertex_descriptor start_node, 
                                      graph::any_graph::vertex_descriptor goal_node, 
                                      double goal_distance,
                                      graph::any_graph& g) {
      return detail::register_basic_solution_path_impl<solution_path_wrapper>(*(this->space), g, start_node, goal_node, goal_pos, goal_distance, solutions);
    };
    
    virtual solution_record_ptr 
      register_joining_point_from_optimal_mg(graph::any_graph::vertex_descriptor start_node, 
                                             graph::any_graph::vertex_descriptor goal_node, 
                                             graph::any_graph::vertex_descriptor join1_node, 
                                             graph::any_graph::vertex_descriptor join2_node, 
                                             double joining_distance,
                                             graph::any_graph& g1, 
                                             graph::any_graph& g2) {
      return detail::register_optimal_solution_path_impl<solution_path_wrapper>(*(this->space), g1, g2, start_node, goal_node, join1_node, join2_node, joining_distance, solutions);
    };
    
    virtual solution_record_ptr 
      register_joining_point_from_basic_mg(graph::any_graph::vertex_descriptor start_node, 
                                           graph::any_graph::vertex_descriptor goal_node, 
                                           graph::any_graph::vertex_descriptor join1_node, 
                                           graph::any_graph::vertex_descriptor join2_node, 
                                           double joining_distance,
                                           graph::any_graph& g1, 
                                           graph::any_graph& g2) {
      return detail::register_basic_solution_path_impl<solution_path_wrapper>(*(this->space), g1, g2, start_node, goal_node, join1_node, join2_node, joining_distance, solutions);
    };
    
    
  public:
    
    
    /**
     * Parametrized constructor.
     * \param aName The name for this object.
     * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
     */
    path_planning_p2p_query(const std::string& aName,
                            const shared_ptr< space_type >& aWorld,
                            const point_type& aStartPos,
                            const point_type& aGoalPos,
                            std::size_t aMaxNumResults = 1) :
                            base_type(aName, aWorld), start_pos(aStartPos), goal_pos(aGoalPos), max_num_results(aMaxNumResults) { };
    
    virtual ~path_planning_p2p_query() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(start_pos)
        & RK_SERIAL_SAVE_WITH_NAME(goal_pos)
        & RK_SERIAL_SAVE_WITH_NAME(max_num_results);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(start_pos)
        & RK_SERIAL_LOAD_WITH_NAME(goal_pos)
        & RK_SERIAL_LOAD_WITH_NAME(max_num_results);
      solutions.clear();
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2460017,1,"path_planning_p2p_query",base_type)
};


};

};

#endif


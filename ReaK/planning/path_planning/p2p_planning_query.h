/**
 * \file p2p_planning_query.h
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

#ifndef REAK_PLANNING_PATH_PLANNING_P2P_PLANNING_QUERY_H_
#define REAK_PLANNING_PATH_PLANNING_P2P_PLANNING_QUERY_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/core/base/named_object.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/random_sampler_concept.h"
#include "ReaK/topologies/spaces/steerable_space_concept.h"
#include "ReaK/topologies/spaces/subspace_concept.h"

#include "ReaK/planning/path_planning/planning_queries.h"
#include "ReaK/topologies/interpolation/discrete_point_path.h"
#include "ReaK/topologies/interpolation/discrete_point_trajectory.h"
#include "ReaK/topologies/interpolation/point_to_point_path.h"
#include "ReaK/topologies/interpolation/point_to_point_trajectory.h"
#include "ReaK/topologies/interpolation/seq_path_wrapper.h"
#include "ReaK/topologies/interpolation/seq_trajectory_wrapper.h"

#include "ReaK/planning/path_planning/any_motion_graphs.h"

#include "ReaK/planning/path_planning/solution_path_factories.h"

#include <map>
#include <type_traits>

namespace ReaK::pp {

namespace details {

template <SteerableSpace Space>
requires TemporalSpace<Space> auto p2p_solution_representation_impl() {
  return seq_trajectory_wrapper<discrete_point_trajectory<Space>>();
}
template <SteerableSpace Space>
auto p2p_solution_representation_impl() {
  return seq_path_wrapper<discrete_point_path<Space>>();
}
template <TemporalSpace Space>
auto p2p_solution_representation_impl() {
  return seq_trajectory_wrapper<point_to_point_trajectory<Space>>();
}
template <Topology Space>
auto p2p_solution_representation_impl() {
  return seq_path_wrapper<point_to_point_path<Space>>();
}

}  // namespace details

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
  using self = path_planning_p2p_query<FreeSpaceType>;
  using base_type = planning_query<FreeSpaceType>;
  using space_type = typename base_type::space_type;
  using super_space_type = typename base_type::super_space_type;

  using point_type = typename base_type::point_type;
  using point_difference_type = typename base_type::point_difference_type;

  using solution_record_ptr = typename base_type::solution_record_ptr;

  using solution_path_wrapper =
      decltype(details::p2p_solution_representation_impl<super_space_type>());

  point_type start_pos;
  point_type goal_pos;
  std::size_t max_num_results;

  std::map<double, solution_record_ptr> solutions;

  /**
   * Returns the best solution distance registered in this query object.
   * \return The best solution distance registered in this query object.
   */
  double get_best_solution_distance() const override {
    if (solutions.size() == 0) {
      return std::numeric_limits<double>::infinity();
    }
    return solutions.begin()->first;
  }

  /**
   * Returns true if the solver should keep on going trying to solve the path-planning problem.
   * \return True if the solver should keep on going trying to solve the path-planning problem.
   */
  bool keep_going() const override {
    return (max_num_results > solutions.size());
  }

  /**
   * This function is called to reset the internal state of the planner.
   */
  void reset_solution_records() override { solutions.clear(); }

  const point_type& get_start_position() const override { return start_pos; }

  double get_distance_to_goal(const point_type& pos) override {
    return get(distance_metric, *(this->space))(pos, goal_pos, *(this->space));
  }

  double get_heuristic_to_goal(const point_type& pos) override {
    return get(distance_metric, this->space->get_super_space())(
        pos, goal_pos, this->space->get_super_space());
  }

 protected:
  solution_record_ptr register_solution_from_optimal_mg(
      bagl::dynamic_graph_observer::vertex_descriptor start_node,
      bagl::dynamic_graph_observer::vertex_descriptor goal_node,
      double goal_distance, bagl::dynamic_graph_observer& g) override {
    return detail::register_optimal_solution_path_impl<solution_path_wrapper>(
        *(this->space), g, start_node, goal_node, goal_pos, goal_distance,
        solutions);
  }

  solution_record_ptr register_solution_from_basic_mg(
      bagl::dynamic_graph_observer::vertex_descriptor start_node,
      bagl::dynamic_graph_observer::vertex_descriptor goal_node,
      double goal_distance, bagl::dynamic_graph_observer& g) override {
    return detail::register_basic_solution_path_impl<solution_path_wrapper>(
        *(this->space), g, start_node, goal_node, goal_pos, goal_distance,
        solutions);
  }

  solution_record_ptr register_joining_point_from_optimal_mg(
      bagl::dynamic_graph_observer::vertex_descriptor start_node,
      bagl::dynamic_graph_observer::vertex_descriptor goal_node,
      bagl::dynamic_graph_observer::vertex_descriptor join1_node,
      bagl::dynamic_graph_observer::vertex_descriptor join2_node,
      double joining_distance, bagl::dynamic_graph_observer& g1,
      bagl::dynamic_graph_observer& g2) override {
    return detail::register_optimal_solution_path_impl<solution_path_wrapper>(
        *(this->space), g1, g2, start_node, goal_node, join1_node, join2_node,
        joining_distance, solutions);
  }

  solution_record_ptr register_joining_point_from_basic_mg(
      bagl::dynamic_graph_observer::vertex_descriptor start_node,
      bagl::dynamic_graph_observer::vertex_descriptor goal_node,
      bagl::dynamic_graph_observer::vertex_descriptor join1_node,
      bagl::dynamic_graph_observer::vertex_descriptor join2_node,
      double joining_distance, bagl::dynamic_graph_observer& g1,
      bagl::dynamic_graph_observer& g2) override {
    return detail::register_basic_solution_path_impl<solution_path_wrapper>(
        *(this->space), g1, g2, start_node, goal_node, join1_node, join2_node,
        joining_distance, solutions);
  }

 public:
  /**
   * Parametrized constructor.
   * \param aName The name for this object.
   * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
   */
  path_planning_p2p_query(const std::string& aName,
                          const std::shared_ptr<space_type>& aWorld,
                          const point_type& aStartPos,
                          const point_type& aGoalPos,
                          std::size_t aMaxNumResults = 1)
      : base_type(aName, aWorld),
        start_pos(aStartPos),
        goal_pos(aGoalPos),
        max_num_results(aMaxNumResults) {}

  ~path_planning_p2p_query() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_type::save(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(start_pos) &
        RK_SERIAL_SAVE_WITH_NAME(goal_pos) &
        RK_SERIAL_SAVE_WITH_NAME(max_num_results);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base_type::load(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(start_pos) &
        RK_SERIAL_LOAD_WITH_NAME(goal_pos) &
        RK_SERIAL_LOAD_WITH_NAME(max_num_results);
    solutions.clear();
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2460017, 1, "path_planning_p2p_query",
                              base_type)
};

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_P2P_PLANNING_QUERY_H_

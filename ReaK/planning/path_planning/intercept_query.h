/**
 * \file intercept_query.h
 *
 * This library defines class templates to encode a motion-planning query to intercept a
 * given trajectory.
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

#ifndef REAK_PLANNING_PATH_PLANNING_INTERCEPT_QUERY_H_
#define REAK_PLANNING_PATH_PLANNING_INTERCEPT_QUERY_H_

#include "ReaK/core/base/named_object.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/random_sampler_concept.h"
#include "ReaK/topologies/spaces/steerable_space_concept.h"
#include "ReaK/topologies/spaces/subspace_concept.h"

#include "ReaK/planning/path_planning/planning_queries.h"

#include "ReaK/topologies/interpolation/discrete_point_trajectory.h"
#include "ReaK/topologies/interpolation/point_to_point_trajectory.h"
#include "ReaK/topologies/interpolation/seq_trajectory_wrapper.h"

#include "ReaK/planning/path_planning/any_motion_graphs.h"

#include "ReaK/planning/path_planning/solution_path_factories.h"

#include "ReaK/math/optimization/optim_exceptions.h"

#include <map>
#include <set>
#include <type_traits>

namespace ReaK::pp {

/**
 * This class template is used to act as a query object for a motion planner, in this class,
 * the planning problem is a target interception problem. This class implements a linear
 * search through the target trajectory to compute the distance between any point and the
 * interception point.
 * \tparam FreeSpaceType The topology of the planning problem, including obstacles (i.e., the free-space).
 * \tparam TargetTrajectory The trajectory that the target follows (or will follow).
 */
template <typename FreeSpaceType, typename TargetTrajectory>
class motion_plan_intercept_query : public planning_query<FreeSpaceType> {
 public:
  using self = motion_plan_intercept_query<FreeSpaceType, TargetTrajectory>;
  using base_type = planning_query<FreeSpaceType>;
  using space_type = typename base_type::space_type;
  using super_space_type = typename base_type::super_space_type;

  using point_type = typename base_type::point_type;
  using point_difference_type = typename base_type::point_difference_type;

  using solution_record_ptr = typename base_type::solution_record_ptr;

  using solution_trajectory_wrapper = std::conditional_t<
      is_steerable_space_v<space_type>,
      seq_trajectory_wrapper<discrete_point_trajectory<super_space_type>>,
      seq_trajectory_wrapper<point_to_point_trajectory<super_space_type>>>;

  using target_trajectory_type = TargetTrajectory;

  point_type start_pos;
  std::shared_ptr<target_trajectory_type> goal_traj;
  std::set<point_type, temporal_point_time_ordering> goal_knots;
  std::size_t max_num_results;
  double maximum_horizon;
  double horizon_decimation;

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
   * Returns true if the solver should keep on going trying to solve the motion-planning problem.
   * \return True if the solver should keep on going trying to solve the motion-planning problem.
   */
  bool keep_going() const override {
    return (max_num_results > solutions.size());
  }

  /**
   * This function is called to reset the internal state of the planner.
   */
  void reset_solution_records() override { solutions.clear(); }

  const point_type& get_start_position() const override { return start_pos; }

  void set_start_position(const point_type& pt) { start_pos = pt; }

  void clear_cached_intercepts() { goal_knots.clear(); }

  void scan_for_goal_knots() {
    int init_knot_count = goal_knots.size();
    for (double t = goal_traj->get_start_time();
         t < goal_traj->get_end_time() + horizon_decimation * 0.5;
         t += horizon_decimation) {
      try {
        point_type goal_pos = goal_traj->get_point_at_time(t);
        if (this->space->is_free(goal_pos)) {
          goal_knots.insert(goal_pos);
        }
      } catch ([[maybe_unused]] std::exception& e) {}
    }

    if (init_knot_count == goal_knots.size()) {
      throw optim::infeasible_problem(
          "Interception problem is infeasible! Not a single point on the "
          "target "
          "trajectory is feasible and collision-free!");
    }
  }

  std::pair<point_type, double> get_distance_position_to_goal(
      const point_type& pos) {
    for (const auto& knot : goal_knots) {
      if (knot.time <= pos.time) {
        continue;
      }
      double tmp =
          get(distance_metric, *(this->space))(pos, knot, *(this->space));
      if (tmp != std::numeric_limits<double>::infinity()) {
        return {knot, tmp};
      }
    }

    for (double t = pos.time;
         t < maximum_horizon + pos.time + horizon_decimation * 0.5;
         t += horizon_decimation) {
      try {
        point_type goal_pos = goal_traj->get_point_at_time(t);
        double tmp =
            get(distance_metric, *(this->space))(pos, goal_pos, *(this->space));
        if (tmp != std::numeric_limits<double>::infinity()) {
          goal_knots.insert(goal_pos);
          return {goal_pos, tmp};
        }
      } catch ([[maybe_unused]] std::exception& e) {}
    }

    return {pos, std::numeric_limits<double>::infinity()};
  }

  double get_distance_to_goal(const point_type& pos) override {
    return get_distance_position_to_goal(pos).second;
  }

  double get_heuristic_to_goal(const point_type& pos) override {
    for (const auto& knot : goal_knots) {
      if (knot.time <= pos.time) {
        continue;
      }
      double tmp = get(distance_metric, this->space->get_super_space())(
          pos, knot, this->space->get_super_space());
      if (tmp != std::numeric_limits<double>::infinity()) {
        return tmp;
      }
    }

    for (double t = pos.time;
         t < maximum_horizon + pos.time + horizon_decimation * 0.5;
         t += horizon_decimation) {
      try {
        point_type goal_pos = goal_traj->get_point_at_time(t);
        double tmp = get(distance_metric, this->space->get_super_space())(
            pos, goal_pos, this->space->get_super_space());
        if (tmp != std::numeric_limits<double>::infinity()) {
          goal_knots.insert(goal_pos);
          return tmp;
        }
      } catch ([[maybe_unused]] std::exception& e) {}
    }

    return std::numeric_limits<double>::infinity();
  }

 protected:
  solution_record_ptr register_solution_from_optimal_mg(
      bagl::dynamic_graph_observer::vertex_descriptor start_node,
      bagl::dynamic_graph_observer::vertex_descriptor goal_node,
      double goal_distance, bagl::dynamic_graph_observer& g) override {
    auto position = bagl::get_dynamic_property_map<const point_type&>(
        "vertex_position", g.get_properties());
    std::pair<point_type, double> goal_pos_dist =
        get_distance_position_to_goal(position[goal_node]);
    return detail::register_optimal_solution_path_impl<
        solution_trajectory_wrapper>(*(this->space), g, start_node, goal_node,
                                     goal_pos_dist.first, goal_pos_dist.second,
                                     solutions);
  }

  solution_record_ptr register_solution_from_basic_mg(
      bagl::dynamic_graph_observer::vertex_descriptor start_node,
      bagl::dynamic_graph_observer::vertex_descriptor goal_node,
      double goal_distance, bagl::dynamic_graph_observer& g) override {
    auto position = bagl::get_dynamic_property_map<const point_type&>(
        "vertex_position", g.get_properties());
    std::pair<point_type, double> goal_pos_dist =
        get_distance_position_to_goal(position[goal_node]);
    return detail::register_basic_solution_path_impl<
        solution_trajectory_wrapper>(*(this->space), g, start_node, goal_node,
                                     goal_pos_dist.first, goal_pos_dist.second,
                                     solutions);
  }

  solution_record_ptr register_joining_point_from_optimal_mg(
      bagl::dynamic_graph_observer::vertex_descriptor start_node,
      bagl::dynamic_graph_observer::vertex_descriptor goal_node,
      bagl::dynamic_graph_observer::vertex_descriptor join1_node,
      bagl::dynamic_graph_observer::vertex_descriptor join2_node,
      double joining_distance, bagl::dynamic_graph_observer& g1,
      bagl::dynamic_graph_observer& g2) override {
    return detail::register_optimal_solution_path_impl<
        solution_trajectory_wrapper>(*(this->space), g1, g2, start_node,
                                     goal_node, join1_node, join2_node,
                                     joining_distance, solutions);
  }

  solution_record_ptr register_joining_point_from_basic_mg(
      bagl::dynamic_graph_observer::vertex_descriptor start_node,
      bagl::dynamic_graph_observer::vertex_descriptor goal_node,
      bagl::dynamic_graph_observer::vertex_descriptor join1_node,
      bagl::dynamic_graph_observer::vertex_descriptor join2_node,
      double joining_distance, bagl::dynamic_graph_observer& g1,
      bagl::dynamic_graph_observer& g2) override {
    return detail::register_basic_solution_path_impl<
        solution_trajectory_wrapper>(*(this->space), g1, g2, start_node,
                                     goal_node, join1_node, join2_node,
                                     joining_distance, solutions);
  }

 public:
  /**
   * Parametrized constructor.
   * \param aName The name for this object.
   * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
   * \param aStartPos The starting position (current position) of the planning problem (current config of robot).
   * \param aGoalTraj The goal trajectory, i.e., the trajectory of the interception target.
   * \param aMaximumHorizon The maximum time forward to try to check a successful interception when computing the
   * "distance to goal".
   * \param aHorizonDecimation The time interval between checks of a successful interception when computing the
   * "distance to goal".
   * \param aMaxNumResults The maximum number of solutions to record and consider the planning successful once reached.
   */
  motion_plan_intercept_query(
      const std::string& aName, const std::shared_ptr<space_type>& aWorld,
      const point_type& aStartPos,
      const std::shared_ptr<target_trajectory_type>& aGoalTraj,
      double aMaximumHorizon, double aHorizonDecimation,
      std::size_t aMaxNumResults = 1)
      : base_type(aName, aWorld),
        start_pos(aStartPos),
        goal_traj(aGoalTraj),
        max_num_results(aMaxNumResults),
        maximum_horizon(aMaximumHorizon),
        horizon_decimation(aHorizonDecimation) {
    scan_for_goal_knots();
  }

  ~motion_plan_intercept_query() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_type::save(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(start_pos) &
        RK_SERIAL_SAVE_WITH_NAME(goal_traj) &
        RK_SERIAL_SAVE_WITH_NAME(max_num_results) &
        RK_SERIAL_SAVE_WITH_NAME(maximum_horizon) &
        RK_SERIAL_SAVE_WITH_NAME(horizon_decimation);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base_type::load(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(start_pos) &
        RK_SERIAL_LOAD_WITH_NAME(goal_traj) &
        RK_SERIAL_LOAD_WITH_NAME(max_num_results) &
        RK_SERIAL_LOAD_WITH_NAME(maximum_horizon) &
        RK_SERIAL_LOAD_WITH_NAME(horizon_decimation);
    solutions.clear();
    goal_knots.clear();
    scan_for_goal_knots();
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2460018, 1,
                              "motion_plan_intercept_query", base_type)
};

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_INTERCEPT_QUERY_H_

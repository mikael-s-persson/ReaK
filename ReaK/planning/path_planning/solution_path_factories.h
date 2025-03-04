/**
 * \file solution_path_factories.h
 *
 * This library contains implementations details for path-planning queries to use for constructing
 * solution paths or trajectories.
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

#ifndef REAK_PLANNING_PATH_PLANNING_SOLUTION_PATH_FACTORIES_H_
#define REAK_PLANNING_PATH_PLANNING_SOLUTION_PATH_FACTORIES_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/planning/path_planning/any_motion_graphs.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/steerable_space_concept.h"
#include "ReaK/topologies/spaces/subspace_concept.h"

#include <type_traits>

#include <map>

namespace ReaK::pp::detail {

template <typename SolutionWrapperType, typename FreeSpaceType,
          typename SolutionRecPtr>
SolutionRecPtr register_basic_solution_path_impl(
    FreeSpaceType& space, const bagl::dynamic_graph_observer& g,
    bagl::dynamic_graph_observer::vertex_descriptor start_node,
    bagl::dynamic_graph_observer::vertex_descriptor goal_node,
    const topology_point_type_t<FreeSpaceType>& goal_pos, double goal_distance,
    std::map<double, SolutionRecPtr>& solutions) {
  using super_space_type =
      typename subspace_traits<FreeSpaceType>::super_space_type;
  using point_type = topology_point_type_t<super_space_type>;
  using Vertex = bagl::dynamic_graph_observer::vertex_descriptor;
  using Edge = bagl::dynamic_graph_observer::edge_descriptor;

  using SolutionPathType = typename SolutionWrapperType::wrapped_type;

  std::shared_ptr<super_space_type> sup_space_ptr(&(space.get_super_space()),
                                                  null_deleter());

  auto position = bagl::get_dynamic_property_map<const point_type&>(
      "vertex_position", g.get_properties());

  double solutions_total_dist = goal_distance;

  auto new_sol = std::make_shared<SolutionWrapperType>(
      "planning_solution",
      SolutionPathType(sup_space_ptr, get(distance_metric, *sup_space_ptr)));
  SolutionPathType& waypoints = new_sol->get_wrapped_object();

  if (goal_distance > 0.0) {
    if constexpr (is_steerable_space_v<FreeSpaceType>) {
      auto [goal_steer_pt, goal_steer_record] =
          space.steer_position_toward(position[goal_node], 1.0, goal_pos);
      RK_UNUSED(goal_steer_pt);
      waypoints.push_front(goal_pos);
      for (auto it = goal_steer_record.end_fraction_travel();
           it != goal_steer_record.begin_fraction_travel(); it -= 1.0) {
        waypoints.push_front(point_type(*it));
      }
    } else {
      waypoints.push_front(goal_pos);
    }
  }

  waypoints.push_front(position[goal_node]);

  while ((in_degree(goal_node, g)) &&
         (!g.are_vertices_equal(goal_node, start_node))) {
    Edge e = *(in_edges(goal_node, g).begin());
    Vertex u = source(e, g);
    if constexpr (is_steerable_space_v<FreeSpaceType>) {
      using SteerRecordType = steerable_space_steer_record_t<super_space_type>;
      auto steer_record =
          bagl::get_dynamic_property_map<const SteerRecordType&>(
              "edge_steer_record", g.get_properties());
      const auto& sr = steer_record[e];
      for (auto it = sr.end_fraction_travel(); it != sr.begin_fraction_travel();
           it -= 1.0) {
        waypoints.push_front(point_type(*it));
      }
    }
    solutions_total_dist += get(distance_metric, *sup_space_ptr)(
        position[u], position[goal_node], *sup_space_ptr);
    waypoints.push_front(position[u]);
    goal_node = u;
  }

  if (g.are_vertices_equal(goal_node, start_node) &&
      ((solutions.empty()) ||
       (solutions_total_dist < solutions.begin()->first))) {
    solutions[solutions_total_dist] = new_sol;
    return new_sol;
  }

  return SolutionRecPtr();
};

template <typename SolutionWrapperType, typename FreeSpaceType,
          typename SolutionRecPtr>
SolutionRecPtr register_optimal_solution_path_impl(
    FreeSpaceType& space, const bagl::dynamic_graph_observer& g,
    bagl::dynamic_graph_observer::vertex_descriptor start_node,
    bagl::dynamic_graph_observer::vertex_descriptor goal_node,
    const topology_point_type_t<FreeSpaceType>& goal_pos, double goal_distance,
    std::map<double, SolutionRecPtr>& solutions) {
  using super_space_type =
      typename subspace_traits<FreeSpaceType>::super_space_type;
  using point_type = topology_point_type_t<super_space_type>;
  using Vertex = bagl::dynamic_graph_observer::vertex_descriptor;

  using SolutionPathType = typename SolutionWrapperType::wrapped_type;

  std::shared_ptr<super_space_type> sup_space_ptr(&(space.get_super_space()),
                                                  null_deleter());

  auto position = bagl::get_dynamic_property_map<const point_type&>(
      "vertex_position", g.get_properties());
  auto predecessor = bagl::get_dynamic_property_map<const std::size_t&>(
      "vertex_predecessor", g.get_properties());
  auto distance_accum = bagl::get_dynamic_property_map<const double&>(
      "vertex_distance_accum", g.get_properties());

  double solutions_total_dist = distance_accum[goal_node] + goal_distance;

  if (!(solutions_total_dist < std::numeric_limits<double>::infinity()) ||
      ((!solutions.empty()) &&
       (solutions_total_dist >= solutions.begin()->first))) {
    return SolutionRecPtr();
  }

  auto new_sol = std::make_shared<SolutionWrapperType>(
      "planning_solution",
      SolutionPathType(sup_space_ptr, get(distance_metric, *sup_space_ptr)));
  SolutionPathType& waypoints = new_sol->get_wrapped_object();

  if (goal_distance > 0.0) {
    if constexpr (is_steerable_space_v<FreeSpaceType>) {
      auto [goal_steer_pt, goal_steer_record] =
          space.steer_position_toward(position[goal_node], 1.0, goal_pos);
      RK_UNUSED(goal_steer_pt);
      waypoints.push_front(goal_pos);
      for (auto it = goal_steer_record.end_fraction_travel();
           it != goal_steer_record.begin_fraction_travel(); it -= 1.0) {
        waypoints.push_front(point_type(*it));
      }
    } else {
      waypoints.push_front(goal_pos);
    }
  }

  waypoints.push_front(position[goal_node]);

  while (!g.are_vertices_equal(goal_node, start_node)) {
    goal_node = Vertex(std::any(predecessor[goal_node]));
    if constexpr (is_steerable_space_v<FreeSpaceType>) {
      auto ei_rg = in_edges(goal_node, g);
      auto ei = ei_rg.begin();
      auto ei_end = ei_rg.end();
      while ((ei != ei_end) &&
             (!g.are_vertices_equal(goal_node, source(*ei, g)))) {
        ++ei;
      }
      if (ei == ei_end) {
        break;
      }
      using SteerRecordType = steerable_space_steer_record_t<super_space_type>;
      auto steer_record =
          bagl::get_dynamic_property_map<const SteerRecordType&>(
              "edge_steer_record", g.get_properties());
      const auto& sr = steer_record[*ei];
      for (auto it = sr.end_fraction_travel(); it != sr.begin_fraction_travel();
           it -= 1.0) {
        waypoints.push_front(point_type(*it));
      }
    }
    waypoints.push_front(position[goal_node]);
  }

  if (g.are_vertices_equal(goal_node, start_node)) {
    solutions[solutions_total_dist] = new_sol;
    return new_sol;
  }

  return SolutionRecPtr();
}

template <typename SolutionWrapperType, typename FreeSpaceType,
          typename SolutionRecPtr>
SolutionRecPtr register_basic_solution_path_impl(
    FreeSpaceType& space, const bagl::dynamic_graph_observer& g1,
    const bagl::dynamic_graph_observer& g2,
    bagl::dynamic_graph_observer::vertex_descriptor start_node,
    bagl::dynamic_graph_observer::vertex_descriptor goal_node,
    bagl::dynamic_graph_observer::vertex_descriptor join1_node,
    bagl::dynamic_graph_observer::vertex_descriptor join2_node,
    double joining_distance, std::map<double, SolutionRecPtr>& solutions) {
  using super_space_type =
      typename subspace_traits<FreeSpaceType>::super_space_type;
  using point_type = topology_point_type_t<super_space_type>;
  using Vertex = bagl::dynamic_graph_observer::vertex_descriptor;
  using Edge = bagl::dynamic_graph_observer::edge_descriptor;

  using SolutionPathType = typename SolutionWrapperType::wrapped_type;

  std::shared_ptr<super_space_type> sup_space_ptr(&(space.get_super_space()),
                                                  null_deleter());

  auto position1 = bagl::get_dynamic_property_map<const point_type&>(
      "vertex_position", g1.get_properties());
  auto position2 = bagl::get_dynamic_property_map<const point_type&>(
      "vertex_position", g2.get_properties());

  double solutions_total_dist = joining_distance;

  auto new_sol = std::make_shared<SolutionWrapperType>(
      "planning_solution",
      SolutionPathType(sup_space_ptr, get(distance_metric, *sup_space_ptr)));
  SolutionPathType& waypoints = new_sol->get_wrapped_object();

  if constexpr (is_steerable_space_v<FreeSpaceType>) {
    if (joining_distance > 0.0) {
      auto [join_steer_pt, join_steer_record] = space.steer_position_toward(
          position1[join1_node], 1.0, position2[join2_node]);
      for (auto it = join_steer_record.end_fraction_travel();
           it != join_steer_record.begin_fraction_travel(); it -= 1.0) {
        waypoints.push_front(point_type(*it));
      }
    }
  }

  waypoints.push_front(position1[join1_node]);

  while ((in_degree(join1_node, g1)) &&
         (!g1.are_vertices_equal(join1_node, start_node))) {
    Edge e = *(in_edges(join1_node, g1).begin());
    Vertex u = source(e, g1);
    if constexpr (is_steerable_space_v<FreeSpaceType>) {
      using SteerRecordType = steerable_space_steer_record_t<super_space_type>;
      auto steer_record1 =
          bagl::get_dynamic_property_map<const SteerRecordType&>(
              "edge_steer_record", g1.get_properties());
      const auto& sr = steer_record1[e];
      for (auto it = sr.end_fraction_travel(); it != sr.begin_fraction_travel();
           it -= 1.0) {
        waypoints.push_front(point_type(*it));
      }
    }
    solutions_total_dist += get(distance_metric, *sup_space_ptr)(
        position1[u], position1[join1_node], *sup_space_ptr);
    join1_node = u;
    waypoints.push_front(position1[join1_node]);
  }

  waypoints.push_back(position2[join2_node]);

  while ((in_degree(join2_node, g2)) &&
         (!g2.are_vertices_equal(join2_node, goal_node))) {
    Edge e = *(in_edges(join2_node, g2).begin());
    Vertex u = source(e, g2);
    if constexpr (is_steerable_space_v<FreeSpaceType>) {
      using SteerRecordType = steerable_space_steer_record_t<super_space_type>;
      auto steer_record2 =
          bagl::get_dynamic_property_map<const SteerRecordType&>(
              "edge_steer_record", g2.get_properties());
      const auto& sr = steer_record2[e];
      for (auto it = sr.end_fraction_travel(); it != sr.begin_fraction_travel();
           it -= 1.0) {
        waypoints.push_back(point_type(*it));
      }
    }
    solutions_total_dist += get(distance_metric, *sup_space_ptr)(
        position2[u], position2[join2_node], *sup_space_ptr);
    join2_node = u;
    waypoints.push_back(position2[join2_node]);
  }

  if (g1.are_vertices_equal(join1_node, start_node) &&
      g2.are_vertices_equal(join2_node, goal_node) &&
      ((solutions.empty()) ||
       (solutions_total_dist < solutions.begin()->first))) {
    solutions[solutions_total_dist] = new_sol;
    return new_sol;
  }

  return SolutionRecPtr();
}

template <typename SolutionWrapperType, typename FreeSpaceType,
          typename SolutionRecPtr>
SolutionRecPtr register_optimal_solution_path_impl(
    FreeSpaceType& space, const bagl::dynamic_graph_observer& g1,
    const bagl::dynamic_graph_observer& g2,
    bagl::dynamic_graph_observer::vertex_descriptor start_node,
    bagl::dynamic_graph_observer::vertex_descriptor goal_node,
    bagl::dynamic_graph_observer::vertex_descriptor join1_node,
    bagl::dynamic_graph_observer::vertex_descriptor join2_node,
    double joining_distance, std::map<double, SolutionRecPtr>& solutions) {
  using super_space_type =
      typename subspace_traits<FreeSpaceType>::super_space_type;
  using point_type = topology_point_type_t<super_space_type>;
  using Vertex = bagl::dynamic_graph_observer::vertex_descriptor;

  using SolutionPathType = typename SolutionWrapperType::wrapped_type;

  std::shared_ptr<super_space_type> sup_space_ptr(&(space.get_super_space()),
                                                  null_deleter());

  auto position1 = bagl::get_dynamic_property_map<const point_type&>(
      "vertex_position", g1.get_properties());
  auto predecessor1 = bagl::get_dynamic_property_map<const std::size_t&>(
      "vertex_predecessor", g1.get_properties());
  auto distance_accum1 = bagl::get_dynamic_property_map<const double&>(
      "vertex_distance_accum", g1.get_properties());
  auto position2 = bagl::get_dynamic_property_map<const point_type&>(
      "vertex_position", g2.get_properties());
  auto predecessor2 = bagl::get_dynamic_property_map<const std::size_t&>(
      "vertex_predecessor", g2.get_properties());
  auto distance_accum2 = bagl::get_dynamic_property_map<const double&>(
      "vertex_distance_accum", g2.get_properties());

  double solutions_total_dist = distance_accum1[join1_node] +
                                distance_accum2[join2_node] + joining_distance;

  if (!(solutions_total_dist < std::numeric_limits<double>::infinity()) ||
      ((!solutions.empty()) &&
       (solutions_total_dist >= solutions.begin()->first))) {
    return SolutionRecPtr();
  }

  auto new_sol = std::make_shared<SolutionWrapperType>(
      "planning_solution",
      SolutionPathType(sup_space_ptr, get(distance_metric, *sup_space_ptr)));
  SolutionPathType& waypoints = new_sol->get_wrapped_object();

  if constexpr (is_steerable_space_v<FreeSpaceType>) {
    if (joining_distance > 0.0) {
      auto [join_steer_pt, join_steer_record] = space.steer_position_toward(
          position1[join1_node], 1.0, position2[join2_node]);
      for (auto it = join_steer_record.end_fraction_travel();
           it != join_steer_record.begin_fraction_travel(); it -= 1.0) {
        waypoints.push_front(point_type(*it));
      }
    }
  }

  waypoints.push_front(position1[join1_node]);

  while (!g1.are_vertices_equal(join1_node, start_node)) {
    join1_node = Vertex(std::any(predecessor1[join1_node]));
    if constexpr (is_steerable_space_v<FreeSpaceType>) {
      auto ei_rg = in_edges(join1_node, g1);
      auto ei = ei_rg.begin();
      auto ei_end = ei_rg.end();
      while ((ei != ei_end) &&
             (!g1.are_vertices_equal(join1_node, source(*ei, g1)))) {
        ++ei;
      }
      if (ei == ei_end) {
        break;
      }
      using SteerRecordType = steerable_space_steer_record_t<super_space_type>;
      auto steer_record1 =
          bagl::get_dynamic_property_map<const SteerRecordType&>(
              "edge_steer_record", g1.get_properties());
      const auto& sr = steer_record1[*ei];
      for (auto it = sr.end_fraction_travel(); it != sr.begin_fraction_travel();
           it -= 1.0) {
        waypoints.push_front(point_type(*it));
      }
    }
    waypoints.push_front(position1[join1_node]);
  }

  waypoints.push_back(position2[join2_node]);

  while (!g2.are_vertices_equal(join2_node, goal_node)) {
    join2_node = Vertex(std::any(predecessor2[join2_node]));
    if constexpr (is_steerable_space_v<FreeSpaceType>) {
      auto ei_rg = in_edges(join2_node, g2);
      auto ei = ei_rg.begin();
      auto ei_end = ei_rg.end();
      while ((ei != ei_end) &&
             (!g2.are_vertices_equal(join2_node, source(*ei, g2)))) {
        ++ei;
      }
      if (ei == ei_end) {
        break;
      }
      using SteerRecordType = steerable_space_steer_record_t<super_space_type>;
      auto steer_record2 =
          bagl::get_dynamic_property_map<const SteerRecordType&>(
              "edge_steer_record", g2.get_properties());
      const auto& sr = steer_record2[*ei];
      for (auto it = sr.end_fraction_travel(); it != sr.begin_fraction_travel();
           it -= 1.0) {
        waypoints.push_back(point_type(*it));
      }
    }
    waypoints.push_back(position2[join2_node]);
  }

  if (g1.are_vertices_equal(join1_node, start_node) &&
      g2.are_vertices_equal(join2_node, goal_node)) {
    solutions[solutions_total_dist] = new_sol;
    return new_sol;
  }

  return SolutionRecPtr();
}

}  // namespace ReaK::pp::detail

#endif  // REAK_PLANNING_PATH_PLANNING_SOLUTION_PATH_FACTORIES_H_

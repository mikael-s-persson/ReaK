/**
 * \file planning_queries.h
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

#ifndef REAK_PLANNING_PATH_PLANNING_PLANNING_QUERIES_H_
#define REAK_PLANNING_PATH_PLANNING_PLANNING_QUERIES_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/core/base/named_object.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "ReaK/topologies/spaces/subspace_concept.h"

#include "ReaK/topologies/interpolation/seq_path_base.h"
#include "ReaK/topologies/interpolation/seq_trajectory_base.h"

#include "ReaK/planning/path_planning/any_motion_graphs.h"

#include <type_traits>

#include <map>

namespace ReaK::pp {

/**
 * This class is the basic OOP interface for a path planner.
 * OOP-style planners are useful to hide away
 * the cumbersome details of calling the underlying planning algorithms which are
 * generic programming (GP) style and thus provide a lot more flexibility but are difficult
 * to deal with in the user-space. The OOP planners are meant to offer a much simpler interface,
 * i.e., a member function that "solves the problem" and returns the solution path or trajectory.
 */
template <typename FreeSpaceType>
class planning_query : public named_object {
 public:
  using self = planning_query<FreeSpaceType>;
  using space_type = FreeSpaceType;
  using super_space_type =
      typename subspace_traits<FreeSpaceType>::super_space_type;

  BOOST_CONCEPT_ASSERT((SubSpaceConcept<FreeSpaceType>));

  using point_type = topology_point_type_t<super_space_type>;
  using point_difference_type =
      topology_point_difference_type_t<super_space_type>;

  using solution_base_type =
      std::conditional_t<is_temporal_space_v<space_type>,
                         seq_trajectory_base<super_space_type>,
                         seq_path_base<super_space_type>>;

  using solution_record_ptr = std::shared_ptr<solution_base_type>;

  std::shared_ptr<space_type> space;

  /**
   * Returns the best solution distance registered in this query object.
   * \return The best solution distance registered in this query object.
   */
  virtual double get_best_solution_distance() const {
    return std::numeric_limits<double>::infinity();
  }

  /**
   * Returns true if the solver should keep on going trying to solve the path-planning problem.
   * \return True if the solver should keep on going trying to solve the path-planning problem.
   */
  virtual bool keep_going() const { return true; }

  /**
   * This function is called to reset the internal state of the planner.
   */
  virtual void reset_solution_records() {}

  virtual const point_type& get_start_position() const = 0;

  /**
   * This function returns the distance of the collision-free travel from the given point to the
   * goal region.
   * \param pos The position from which to try and reach the goal.
   * \return The distance of the collision-free travel from the given point to the goal region.
   */
  virtual double get_distance_to_goal(const point_type& pos) {
    return std::numeric_limits<double>::infinity();
  }

  /**
   * This function returns the heuristic distance of the bird-flight travel from the given
   * point to the goal region.
   * \param pos The position from which to reach the goal.
   * \return The heuristic distance of the bird-flight travel from the given
   * point to the goal region.
   */
  virtual double get_heuristic_to_goal(const point_type& pos) {
    return std::numeric_limits<double>::infinity();
  }

 protected:
  virtual solution_record_ptr register_solution_from_optimal_mg(
      graph::any_graph::vertex_descriptor start_node,
      graph::any_graph::vertex_descriptor goal_node, double goal_distance,
      graph::any_graph& g) = 0;

  virtual solution_record_ptr register_solution_from_basic_mg(
      graph::any_graph::vertex_descriptor start_node,
      graph::any_graph::vertex_descriptor goal_node, double goal_distance,
      graph::any_graph& g) = 0;

  virtual solution_record_ptr register_joining_point_from_optimal_mg(
      graph::any_graph::vertex_descriptor start_node,
      graph::any_graph::vertex_descriptor goal_node,
      graph::any_graph::vertex_descriptor join1_node,
      graph::any_graph::vertex_descriptor join2_node, double goal_distance,
      graph::any_graph& g1, graph::any_graph& g2) = 0;

  virtual solution_record_ptr register_joining_point_from_basic_mg(
      graph::any_graph::vertex_descriptor start_node,
      graph::any_graph::vertex_descriptor goal_node,
      graph::any_graph::vertex_descriptor join1_node,
      graph::any_graph::vertex_descriptor join2_node, double goal_distance,
      graph::any_graph& g1, graph::any_graph& g2) = 0;

 public:
  /**
   * This function registers a solution path (if one is found) and invokes the path-planning
   * reporter to report on that solution path.
   * \note This function works for optimal motion graphs (optimal: uses a shortest-distance rewiring strategy).
   * \param start_node The start node in the motion-graph.
   * \param goal_node The goal node in the motion-graph.
   * \param g The current motion-graph.
   * \return True if a new solution was registered.
   */
  template <typename Vertex, typename Graph>
  solution_record_ptr register_solution(Vertex start_node, Vertex goal_node,
                                        double goal_distance, Graph& g) {
    if constexpr (std::is_convertible_v<graph::graph_vertex_bundle_t<Graph>*,
                                        optimal_mg_vertex<FreeSpaceType>*>) {
      using TEGraph = any_optimal_motion_graph<FreeSpaceType, Graph>;
      using TEVertex = graph::graph_vertex_t<TEGraph>;

      TEGraph te_g(&g);
      return register_solution_from_optimal_mg(TEVertex(std::any(start_node)),
                                               TEVertex(std::any(goal_node)),
                                               goal_distance, te_g);
    } else {
      using TEGraph = any_motion_graph<FreeSpaceType, Graph>;
      using TEVertex = graph::graph_vertex_t<TEGraph>;

      TEGraph te_g(&g);
      return register_solution_from_basic_mg(TEVertex(std::any(start_node)),
                                             TEVertex(std::any(goal_node)),
                                             goal_distance, te_g);
    }
  }

  /**
   * This function registers a solution path (if one is found) and invokes the path-planning
   * reporter to report on that solution path.
   * \note This function works for optimal motion graphs (optimal: uses a shortest-distance rewiring strategy).
   * \param start_node The start node in the motion-graph.
   * \param goal_node The goal node in the motion-graph.
   * \param g The current motion-graph.
   * \return True if a new solution was registered.
   */
  template <typename Vertex, typename Graph>
  solution_record_ptr register_joining_point(
      Vertex start_node, Vertex goal_node, Vertex join1_node, Vertex join2_node,
      double joining_distance, Graph& g1, Graph& g2) {
    if constexpr (std::is_convertible_v<graph::graph_vertex_bundle_t<Graph>*,
                                        optimal_mg_vertex<FreeSpaceType>*>) {
      using TEGraph = any_optimal_motion_graph<FreeSpaceType, Graph>;
      using TEVertex = graph::graph_vertex_t<TEGraph>;

      TEGraph te_g1(&g1);
      TEGraph te_g2(&g2);
      return register_joining_point_from_optimal_mg(
          TEVertex(std::any(start_node)), TEVertex(std::any(goal_node)),
          TEVertex(std::any(join1_node)), TEVertex(std::any(join2_node)),
          joining_distance, te_g1, te_g2);
    } else {
      using TEGraph = any_motion_graph<FreeSpaceType, Graph>;
      using TEVertex = graph::graph_vertex_t<TEGraph>;

      TEGraph te_g1(&g1);
      TEGraph te_g2(&g2);
      return register_joining_point_from_basic_mg(
          TEVertex(std::any(start_node)), TEVertex(std::any(goal_node)),
          TEVertex(std::any(join1_node)), TEVertex(std::any(join2_node)),
          joining_distance, te_g1, te_g2);
    }
  }

  /**
   * Parametrized constructor.
   * \param aName The name for this object.
   * \param aWorld A topology which represents the C-free (obstacle-free configuration space).
   */
  planning_query(const std::string& aName,
                 const std::shared_ptr<space_type>& aWorld)
      : named_object(), space(aWorld) {
    setName(aName);
  }

  ~planning_query() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(space);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(space);
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2460016, 1, "planning_query",
                              named_object)
};

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_PLANNING_QUERIES_H_

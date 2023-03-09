/**
 * \file density_plan_visitors.h
 *
 * This library defines
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

#ifndef REAK_PLANNING_PATH_PLANNING_DENSITY_PLAN_VISITORS_H_
#define REAK_PLANNING_PATH_PLANNING_DENSITY_PLAN_VISITORS_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/core/base/named_object.h"

#include "ReaK/planning/path_planning/density_calculators.h"
#include "ReaK/planning/path_planning/planning_visitors.h"

namespace ReaK::pp {

/**
 * This class template
 */
template <typename FreeSpaceType,
          typename DensityCalc = sbastar_density_calculator>
struct density_plan_visitor
    : planning_visitor_base<density_plan_visitor<FreeSpaceType, DensityCalc>,
                            FreeSpaceType> {

  using self = density_plan_visitor<FreeSpaceType, DensityCalc>;
  using base_type = planning_visitor_base<self, FreeSpaceType>;

  using space_type = typename base_type::space_type;
  using planner_base_type = typename base_type::planner_base_type;
  using query_type = typename base_type::query_type;

  double m_density_cutoff;
  DensityCalc m_density_calc;

  explicit density_plan_visitor(planner_base_type* aPlanner,
                                query_type* aQuery = nullptr,
                                any_knn_synchro* aNNSynchro = nullptr,
                                std::any aStartNode = std::any(),
                                std::any aGoalNode = std::any(),
                                double aDensityCutoff = 0.0,
                                DensityCalc aDensityCalc = DensityCalc())
      : base_type(aPlanner, aQuery, aNNSynchro, aStartNode, aGoalNode),
        m_density_cutoff(aDensityCutoff),
        m_density_calc(aDensityCalc) {}

  /***************************************************
                  NodeExploringVisitorConcept
  ***************************************************/

  template <typename Vertex, typename Graph>
  void init_nonrecursive_density(Vertex u, Graph& g) const {
    m_density_calc.update_density(u, g, *(this->m_query->space),
                                  this->m_planner->get_sampling_radius(),
                                  this->m_planner->get_space_dimensionality());
  }

  template <typename Vertex, typename Graph>
  void init_recursive_density(Vertex u, Graph& g) const {
    g[u].constriction = 0.0;
    g[u].collision_count = 0;
    g[u].density = 0.0;
    g[u].expansion_trials = 0;

    m_density_calc.update_density(u, g, *(this->m_query->space),
                                  this->m_planner->get_sampling_radius(),
                                  this->m_planner->get_space_dimensionality());
  }

  template <typename BaseType, typename Vertex, typename Graph>
  void dispatched_initialize_vertex(dense_mg_vertex<BaseType>& vp, Vertex u,
                                    Graph& g) const {
    base_type::dispatched_initialize_vertex(static_cast<BaseType&>(vp), u, g);
    init_nonrecursive_density(u, g);
  }

  template <typename BaseType, typename Vertex, typename Graph>
  void dispatched_initialize_vertex(recursive_dense_mg_vertex<BaseType>& vp,
                                    Vertex u, Graph& g) const {
    base_type::dispatched_initialize_vertex(static_cast<BaseType&>(vp), u, g);
    init_recursive_density(u, g);
  }

  template <typename Vertex, typename Graph>
  void initialize_vertex(Vertex u, Graph& g) const {
    dispatched_initialize_vertex(g[u], u, g);
  }
  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex /*unused*/, const Graph& /*unused*/) const {}
  template <typename Vertex, typename Graph>
  void examine_vertex(Vertex /*unused*/, const Graph& /*unused*/) const {}
  template <typename Edge, typename Graph>
  void examine_edge(Edge /*unused*/, const Graph& /*unused*/) const {}

  bool dispatched_heuristic_potential(
      const mg_vertex_data<space_type>& /*unused*/) const {
    return true;
  }
  bool dispatched_heuristic_potential(
      const astar_mg_vertex<space_type>& vp,
      const astar_mg_vertex<space_type>& sp) const {
    return (vp.heuristic_value >
            std::numeric_limits<double>::epsilon() * sp.heuristic_value);
  }
  bool dispatched_heuristic_potential(
      const bidir_astar_mg_vertex<space_type>& vp,
      const bidir_astar_mg_vertex<space_type>& sp) const {
    return (vp.fwd_distance_accum >
            std::numeric_limits<double>::epsilon() * sp.fwd_distance_accum);
  }

  template <typename BaseType>
  bool dispatched_density_cutoff_test(
      const dense_mg_vertex<BaseType>& vp) const {
    return ((1.0 - vp.density) > m_density_cutoff);
  }
  template <typename BaseType>
  bool dispatched_density_cutoff_test(
      const recursive_dense_mg_vertex<BaseType>& vp) const {
    return ((1.0 - vp.constriction) * (1.0 - vp.density) > m_density_cutoff);
  }

  template <typename Vertex, typename Graph>
  bool has_search_potential(Vertex u, const Graph& g) const {
    if (!this->m_goal_node.has_value()) {
      return dispatched_heuristic_potential(
                 g[u], g[std::any_cast<Vertex>(this->m_start_node)]) &&
             dispatched_density_cutoff_test(g[u]);
    }
    return (u != std::any_cast<Vertex>(this->m_goal_node)) &&
           dispatched_density_cutoff_test(g[u]);
  }
  template <typename Vertex, typename Graph>
  bool should_close(Vertex u, const Graph& g) const {
    return !has_search_potential(u, g);
  }

  /***************************************************
                  AnytimeHeuristicVisitorConcept  (Anytime A* search)
  ***************************************************/

  template <typename Graph>
  double adjust_relaxation(double old_relaxation, const Graph& g) const {
    return old_relaxation * 0.5;
  }

  /***************************************************
                  NeighborhoodTrackingVisitorConcept
  ***************************************************/

  template <typename BaseType, typename Vertex, typename Graph>
  void dispatched_travel_succeeded(dense_mg_vertex<BaseType>& /*unused*/,
                                   Vertex /*unused*/, Vertex /*unused*/,
                                   Graph& /*unused*/) const {}

  template <typename BaseType, typename Vertex, typename Graph>
  void dispatched_travel_succeeded(
      recursive_dense_mg_vertex<BaseType>& /*unused*/, Vertex u, Vertex v,
      Graph& g) const {
    m_density_calc.travel_succeeded(
        u, v, g, *(this->m_query->space),
        this->m_planner->get_sampling_radius(),
        this->m_planner->get_space_dimensionality());
  }

  template <typename Vertex, typename Graph>
  void travel_succeeded(Vertex u, Vertex v, Graph& g) const {
    dispatched_travel_succeeded(g[u], u, v, g);
  }

  template <typename BaseType, typename Vertex, typename Graph>
  void dispatched_travel_explored(dense_mg_vertex<BaseType>& /*unused*/,
                                  Vertex /*unused*/, Vertex /*unused*/,
                                  Graph& /*unused*/) const {}

  template <typename BaseType, typename Vertex, typename Graph>
  void dispatched_travel_explored(
      recursive_dense_mg_vertex<BaseType>& /*unused*/, Vertex u, Vertex v,
      Graph& g) const {
    m_density_calc.travel_explored(u, v, g, *(this->m_query->space),
                                   this->m_planner->get_sampling_radius(),
                                   this->m_planner->get_space_dimensionality());
  }

  template <typename Vertex, typename Graph>
  void travel_explored(Vertex u, Vertex v, Graph& g) const {
    dispatched_travel_explored(g[u], u, v, g);
  }

  template <typename BaseType, typename Vertex, typename Graph>
  void dispatched_travel_failed(dense_mg_vertex<BaseType>& /*unused*/,
                                Vertex /*unused*/, Vertex /*unused*/,
                                Graph& /*unused*/) const {}

  template <typename BaseType, typename Vertex, typename Graph>
  void dispatched_travel_failed(recursive_dense_mg_vertex<BaseType>& /*unused*/,
                                Vertex u, Vertex v, Graph& g) const {
    m_density_calc.travel_failed(u, v, g, *(this->m_query->space),
                                 this->m_planner->get_sampling_radius(),
                                 this->m_planner->get_space_dimensionality());
  }

  template <typename Vertex, typename Graph>
  void travel_failed(Vertex u, Vertex v, Graph& g) const {
    dispatched_travel_failed(g[u], u, v, g);
  }

  template <typename BaseType, typename Vertex, typename Graph>
  void dispatched_affected_vertex(dense_mg_vertex<BaseType>& /*unused*/,
                                  Vertex u, Graph& g) const {
    init_nonrecursive_density(u, g);
  }

  template <typename BaseType, typename Vertex, typename Graph>
  void dispatched_affected_vertex(
      recursive_dense_mg_vertex<BaseType>& /*unused*/, Vertex /*unused*/,
      Graph& /*unused*/) const {}

  template <typename Vertex, typename Graph>
  void affected_vertex(Vertex u, Graph& g) const {
    dispatched_affected_vertex(g[u], u, g);
  }
};

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_DENSITY_PLAN_VISITORS_H_

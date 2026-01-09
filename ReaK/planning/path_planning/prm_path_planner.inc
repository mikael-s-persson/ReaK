/**
 * \file prm_path_planner.inc
 *
 * This library contains template definitions of a class to solve path planning problems using the
 * Probabilistic Road-map (PRM) algorithm (or one of its variants).
 * Given a C_free (configuration space restricted to non-colliding points) and a
 * result reporting policy, this class will probabilistically construct a motion-graph
 * that will connect a starting point and a goal point with a path through C-free
 * that is as close as possible to the optimal path in terms of distance.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_PLANNING_PATH_PLANNING_PRM_PATH_PLANNER_INC_
#define REAK_PLANNING_PATH_PLANNING_PRM_PATH_PLANNER_INC_

#include "ReaK/planning/path_planning/prm_path_planner.h"

#include "ReaK/planning/graph_alg/neighborhood_functors.h"
#include "ReaK/planning/graph_alg/probabilistic_roadmap.h"

#include "ReaK/planning/path_planning/motion_graph_structures.h"

#include "bagl/astar_search.h"
#include "bagl/more_property_maps.h"

#include "ReaK/planning/path_planning/metric_space_search.h"
#include "ReaK/planning/path_planning/topological_search.h"

#include "ReaK/planning/path_planning/any_motion_graphs.h"
#include "ReaK/planning/path_planning/density_plan_visitors.h"
#include "ReaK/planning/path_planning/p2p_planning_query.h"
#include "ReaK/planning/path_planning/path_planner_options.h"

namespace ReaK::pp {

/**
 * This class template is used by the FADPRM path-planner as the visitor object needed to
 * collaborate with the FADPRM algorithms to generate the motion-graph and path-planning solutions.
 * This class template models the FADPRMVisitorConcept.
 */
template <SubSpace FreeSpaceType>
struct prm_planner_visitor
    : density_plan_visitor<FreeSpaceType, prm_density_calculator> {
  using base_type = density_plan_visitor<FreeSpaceType, prm_density_calculator>;
  using self = prm_planner_visitor<FreeSpaceType>;

  using planner_base_type = typename base_type::planner_base_type;
  using query_type = typename base_type::query_type;

  explicit prm_planner_visitor(planner_base_type* aPlanner,
                               query_type* aQuery = nullptr,
                               any_knn_synchro* aNNSynchro = nullptr,
                               std::any aStartNode = std::any(),
                               std::any aGoalNode = std::any(),
                               double aDensityCutoff = 0.0)
      : base_type(aPlanner, aQuery, aNNSynchro, aStartNode, aGoalNode,
                  aDensityCutoff) {}

  template <typename Graph>
  struct astar_heuristic_getter {
    using Vertex = bagl::graph_vertex_descriptor_t<Graph>;
    const Graph* p_g;
    explicit astar_heuristic_getter(const Graph* pG) : p_g(pG) {}
    double operator()(Vertex u) const { return (*p_g)[u].heuristic_value; }
  };

  template <typename Graph>
  void publish_path(Graph& g) const {

    using Vertex = bagl::graph_vertex_descriptor_t<Graph>;
    using VertexProp = dense_mg_vertex<astar_mg_vertex<FreeSpaceType>>;
    using EdgeProp = optimal_mg_edge<FreeSpaceType>;

    Vertex start_node = std::any_cast<Vertex>(this->m_start_node);
    Vertex goal_node = std::any_cast<Vertex>(this->m_goal_node);

    bagl::astar_search(
        g, start_node, astar_heuristic_getter<Graph>(&g),
        bagl::default_astar_visitor(), get(&VertexProp::predecessor, g),
        get(&VertexProp::key_value, g), get(&VertexProp::distance_accum, g),
        get(&EdgeProp::weight, g), bagl::identity_property_map(),
        get(&VertexProp::astar_color, g), std::less<>(), std::plus<>(),
        std::numeric_limits<double>::infinity(), double(0.0));

    this->dispatched_register_solution(start_node, goal_node, goal_node, g,
                                       g[goal_node]);
  }
};

template <SubSpace FreeSpaceType>
void prm_planner<FreeSpaceType>::solve_planning_query(
    planning_query<FreeSpaceType>& aQuery) {

  this->reset_internal_state();

  using SuperSpace = typename subspace_traits<FreeSpaceType>::super_space_type;
  using PointType = topology_point_type_t<SuperSpace>;

  using VertexProp = dense_mg_vertex<astar_mg_vertex<FreeSpaceType>>;
  using EdgeProp = optimal_mg_edge<FreeSpaceType>;

  using DirectionalityTag =
      typename motion_segment_directionality<FreeSpaceType>::type;

  using BasicVertexProp = mg_vertex_data<FreeSpaceType>;

  using PositionMap = bagl::data_member_property_map<PointType, VertexProp>;
  auto pos_map = PositionMap(&VertexProp::position);

  using DensityMap = bagl::data_member_property_map<double, VertexProp>;
  auto dens_map = DensityMap(&VertexProp::density);

  auto space_dim = double(this->get_space_dimensionality());
  double space_Lc = aQuery.get_heuristic_to_goal(aQuery.get_start_position());

  std::shared_ptr<const SuperSpace> sup_space_ptr(
      &(this->m_space->get_super_space()), null_deleter());

  density_plan_visitor<FreeSpaceType, prm_density_calculator> vis(this,
                                                                  &aQuery);

  auto* p2p_query_ptr =
      reinterpret_cast<path_planning_p2p_query<FreeSpaceType>*>(aQuery.cast_to(
          path_planning_p2p_query<FreeSpaceType>::get_static_object_type()));

#define RK_PRM_PLANNER_INITIALIZE_START_AND_GOAL                       \
  VertexProp vp_start;                                                 \
  vp_start.position = aQuery.get_start_position();                     \
  Vertex start_node = add_vertex(motion_graph, std::move(vp_start));   \
  motion_graph[start_node].density = 0.0;                              \
  motion_graph[start_node].heuristic_value =                           \
      aQuery.get_heuristic_to_goal(motion_graph[start_node].position); \
  motion_graph[start_node].distance_accum = 0.0;                       \
  motion_graph[start_node].predecessor = start_node;                   \
  vis.m_start_node = std::any(start_node);                             \
  if (p2p_query_ptr) {                                                 \
    VertexProp vp_goal;                                                \
    vp_goal.position = p2p_query_ptr->goal_pos;                        \
    Vertex goal_node = add_vertex(motion_graph, std::move(vp_goal));   \
    motion_graph[goal_node].density = 0.0;                             \
    motion_graph[goal_node].heuristic_value = 0.0;                     \
    motion_graph[goal_node].distance_accum =                           \
        std::numeric_limits<double>::infinity();                       \
    motion_graph[goal_node].predecessor = goal_node;                   \
    vis.m_goal_node = std::any(goal_node);                             \
  }

#define RK_PRM_PLANNER_SETUP_DVP_TREE_SYNCHRO(ARITY, TREE_STORAGE)            \
  using GraphPositionMap =                                                    \
      bagl::property_map_t<MotionGraphType, PointType BasicVertexProp::*>;    \
  using SpacePartType = dvp_tree<Vertex, SuperSpace, GraphPositionMap, ARITY, \
                                 random_vp_chooser, TREE_STORAGE>;            \
  SpacePartType space_part(motion_graph, sup_space_ptr,                       \
                           get(&BasicVertexProp::position, motion_graph));    \
                                                                              \
  using NNFinderType = multi_dvp_tree_search<MotionGraphType, SpacePartType>; \
  NNFinderType nn_finder;                                                     \
  nn_finder.graph_tree_map[&motion_graph] = &space_part;                      \
                                                                              \
  type_erased_knn_synchro<MotionGraphType, NNFinderType> NN_synchro(          \
      nn_finder);                                                             \
  vis.m_nn_synchro = &NN_synchro;                                             \
                                                                              \
  ReaK::graph::star_neighborhood<NNFinderType> nc_selector(                   \
      nn_finder, space_dim, 3.0 * space_Lc);

#define RK_PRM_PLANNER_SETUP_ALT_TREE_SYNCHRO(ARITY, TREE_STORAGE)             \
  using ALTGraph =                                                             \
      dvp_adjacency_list<VertexProp, EdgeProp, SuperSpace, PositionMap, ARITY, \
                         random_vp_chooser, TREE_STORAGE, bagl::vec_s,         \
                         DirectionalityTag, bagl::list_s>;                     \
  using MotionGraphType = typename ALTGraph::adj_list_type;                    \
  using Vertex = bagl::graph_vertex_descriptor_t<MotionGraphType>;             \
                                                                               \
  ALTGraph space_part(sup_space_ptr, pos_map);                                 \
  MotionGraphType motion_graph = space_part.get_adjacency_list();              \
                                                                               \
  using NNFinderType = multi_dvp_tree_search<MotionGraphType, ALTGraph>;       \
  NNFinderType nn_finder;                                                      \
  nn_finder.graph_tree_map[&motion_graph] = &space_part;                       \
                                                                               \
  any_knn_synchro NN_synchro;                                                  \
  vis.m_nn_synchro = &NN_synchro;                                              \
                                                                               \
  ReaK::graph::star_neighborhood<NNFinderType> nc_selector(                    \
      nn_finder, space_dim, 3.0 * space_Lc);

#define RK_PRM_PLANNER_MAKE_GENERATE_PRM_CALL                              \
  ReaK::graph::generate_prm(motion_graph, *sup_space_ptr, vis, pos_map,    \
                            get(random_sampler, *sup_space_ptr), dens_map, \
                            nc_selector, 0.2);

  if ((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) ==
      ADJ_LIST_MOTION_GRAPH) {

    using MotionGraphType =
        bagl::adjacency_list<bagl::vec_s, bagl::vec_s, DirectionalityTag,
                             VertexProp, EdgeProp>;
    using Vertex = bagl::graph_vertex_descriptor_t<MotionGraphType>;

    MotionGraphType motion_graph;

    RK_PRM_PLANNER_INITIALIZE_START_AND_GOAL

    if ((this->m_data_structure_flags & KNN_METHOD_MASK) == LINEAR_SEARCH_KNN) {

      using NNFinderType = linear_neighbor_search<MotionGraphType>;
      NNFinderType nn_finder;

      ReaK::graph::star_neighborhood<NNFinderType> nc_selector(
          nn_finder, space_dim, 3.0 * space_Lc);

      any_knn_synchro NN_synchro;
      vis.m_nn_synchro = &NN_synchro;

      RK_PRM_PLANNER_MAKE_GENERATE_PRM_CALL

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_BF2_TREE_KNN) {

      RK_PRM_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, bagl::bfl_d_ary_tree_storage<2>)

      RK_PRM_PLANNER_MAKE_GENERATE_PRM_CALL

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_BF4_TREE_KNN) {

      RK_PRM_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, bagl::bfl_d_ary_tree_storage<4>)

      RK_PRM_PLANNER_MAKE_GENERATE_PRM_CALL

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_COB2_TREE_KNN) {

      RK_PRM_PLANNER_SETUP_DVP_TREE_SYNCHRO(2,
                                            bagl::vebl_d_ary_tree_storage<2>)

      RK_PRM_PLANNER_MAKE_GENERATE_PRM_CALL

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_COB4_TREE_KNN) {

      RK_PRM_PLANNER_SETUP_DVP_TREE_SYNCHRO(4,
                                            bagl::vebl_d_ary_tree_storage<4>)

      RK_PRM_PLANNER_MAKE_GENERATE_PRM_CALL

#endif
    }

#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT

  } else if ((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) ==
             DVP_ADJ_LIST_MOTION_GRAPH) {

    if ((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {

      RK_PRM_PLANNER_SETUP_ALT_TREE_SYNCHRO(2, bagl::bfl_d_ary_tree_storage<2>)

      RK_PRM_PLANNER_INITIALIZE_START_AND_GOAL

      RK_PRM_PLANNER_MAKE_GENERATE_PRM_CALL

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_BF4_TREE_KNN) {

      RK_PRM_PLANNER_SETUP_ALT_TREE_SYNCHRO(4, bagl::bfl_d_ary_tree_storage<4>)

      RK_PRM_PLANNER_INITIALIZE_START_AND_GOAL

      RK_PRM_PLANNER_MAKE_GENERATE_PRM_CALL

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_COB2_TREE_KNN) {

      RK_PRM_PLANNER_SETUP_ALT_TREE_SYNCHRO(2,
                                            bagl::vebl_d_ary_tree_storage<2>)

      RK_PRM_PLANNER_INITIALIZE_START_AND_GOAL

      RK_PRM_PLANNER_MAKE_GENERATE_PRM_CALL

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_COB4_TREE_KNN) {

      RK_PRM_PLANNER_SETUP_ALT_TREE_SYNCHRO(4,
                                            bagl::vebl_d_ary_tree_storage<4>)

      RK_PRM_PLANNER_INITIALIZE_START_AND_GOAL

      RK_PRM_PLANNER_MAKE_GENERATE_PRM_CALL

#endif
    }

#endif
  }

#undef RK_PRM_PLANNER_INITIALIZE_START_AND_GOAL
#undef RK_PRM_PLANNER_SETUP_DVP_TREE_SYNCHRO
#undef RK_PRM_PLANNER_SETUP_ALT_TREE_SYNCHRO
#undef RK_PRM_PLANNER_MAKE_GENERATE_PRM_CALL
}

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_PRM_PATH_PLANNER_INC_

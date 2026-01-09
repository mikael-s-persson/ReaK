/**
 * \file rrt_path_planner.inc
 *
 * This library contains template definitions of a class to solve path planning problems using the
 * Rapidly-exploring Random Tree (RRT) algorithm (or one of its variants).
 * Given a C_free (configuration space restricted to non-colliding points) and a
 * result reporting policy, this class will probabilistically construct a motion-graph
 * that will connect a starting point and a goal point with a path through C-free
 * that is as close as possible to the optimal path in terms of distance.
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

#ifndef REAK_PLANNING_PATH_PLANNING_RRT_PATH_PLANNER_INC_
#define REAK_PLANNING_PATH_PLANNING_RRT_PATH_PLANNER_INC_

#include "ReaK/planning/path_planning/rrt_path_planner.h"

#include "ReaK/planning/graph_alg/rr_tree.h"
#include "ReaK/planning/path_planning/motion_graph_structures.h"

#include "ReaK/planning/path_planning/metric_space_search.h"
#include "ReaK/planning/path_planning/topological_search.h"

#include "bagl/more_property_maps.h"
#include "bagl/tree_adaptor.h"

#include "ReaK/planning/path_planning/any_motion_graphs.h"
#include "ReaK/planning/path_planning/p2p_planning_query.h"
#include "ReaK/planning/path_planning/path_planner_options.h"
#include "ReaK/planning/path_planning/planning_visitors.h"

namespace ReaK::pp {

template <SubSpace FreeSpaceType>
void rrt_planner<FreeSpaceType>::solve_planning_query(
    planning_query<FreeSpaceType>& aQuery) {

  this->reset_internal_state();

  using SuperSpace = typename subspace_traits<FreeSpaceType>::super_space_type;
  using PointType = topology_point_type_t<SuperSpace>;

  using VertexProp = mg_vertex_data<FreeSpaceType>;
  using EdgeProp = mg_edge_data<FreeSpaceType>;

  using DirectionalityTag =
      typename motion_segment_directionality<FreeSpaceType>::type;

  using BasicVertexProp = mg_vertex_data<FreeSpaceType>;

  using PositionMap = bagl::data_member_property_map<PointType, VertexProp>;
  auto pos_map = PositionMap(&VertexProp::position);

  std::shared_ptr<const SuperSpace> sup_space_ptr(
      &(this->m_space->get_super_space()), null_deleter());

  planning_visitor<FreeSpaceType> vis(this, &aQuery);

  VertexProp vp_start;
  vp_start.position = aQuery.get_start_position();

  if ((this->m_planning_method_flags & PLANNING_DIRECTIONALITY_MASK) ==
      UNIDIRECTIONAL_PLANNING) {

#define RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(ARITY, TREE_STORAGE)            \
  using Vertex = bagl::graph_vertex_descriptor_t<MotionGraphType>;            \
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
  vis.m_nn_synchro = &NN_synchro;

#define RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(ARITY, TREE_STORAGE)             \
  using ALTGraph =                                                             \
      dvp_adjacency_list<VertexProp, EdgeProp, SuperSpace, PositionMap, ARITY, \
                         random_vp_chooser, TREE_STORAGE, bagl::vec_s,         \
                         DirectionalityTag, bagl::list_s>;                     \
  using MotionGraphType = typename ALTGraph::adj_list_type;                    \
                                                                               \
  ALTGraph space_part(sup_space_ptr, pos_map);                                 \
  MotionGraphType motion_graph = space_part.get_adjacency_list();              \
  vis.m_start_node = std::any(create_root(motion_graph, vp_start));            \
                                                                               \
  using NNFinderType = multi_dvp_tree_search<MotionGraphType, ALTGraph>;       \
  NNFinderType nn_finder;                                                      \
  nn_finder.graph_tree_map[&motion_graph] = &space_part;                       \
                                                                               \
  any_knn_synchro NN_synchro;                                                  \
  vis.m_nn_synchro = &NN_synchro;

#define RK_RRT_PLANNER_CALL_RRT_FUNCTION                                \
  ReaK::graph::generate_rrt(motion_graph, *sup_space_ptr, vis, pos_map, \
                            get(random_sampler, *sup_space_ptr), nn_finder);

    if ((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) ==
        ADJ_LIST_MOTION_GRAPH) {

      using MotionGraphType =
          bagl::adjacency_list<bagl::vec_s, bagl::list_s,
                               DirectionalityTag, VertexProp, EdgeProp>;

      MotionGraphType motion_graph;
      vis.m_start_node = std::any(create_root(motion_graph, vp_start));

      if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
          LINEAR_SEARCH_KNN) {

        any_knn_synchro NN_synchro;
        vis.m_nn_synchro = &NN_synchro;
        linear_neighbor_search<MotionGraphType> nn_finder;

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_BF2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, bagl::bfl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_BF4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, bagl::bfl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(2, bagl::vebl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(4, bagl::vebl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

#endif
      }

#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT

    } else if ((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) ==
               DVP_ADJ_LIST_MOTION_GRAPH) {

      if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
          DVP_BF2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(2, bagl::bfl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_BF4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(4, bagl::bfl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(2, bagl::vebl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(4, bagl::vebl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

#endif
      }

#endif
    }

#undef RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO
#undef RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO
#undef RK_RRT_PLANNER_CALL_RRT_FUNCTION

  } else {
    auto* p2p_query_ptr = reinterpret_cast<
        path_planning_p2p_query<FreeSpaceType>*>(aQuery.cast_to(
        path_planning_p2p_query<FreeSpaceType>::get_static_object_type()));
    if (p2p_query_ptr == nullptr) {
      return;
    }

    VertexProp vp_goal;
    vp_goal.position = p2p_query_ptr->goal_pos;

#define RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(ARITY, TREE_STORAGE)        \
  using Vertex = bagl::graph_vertex_descriptor_t<MotionGraphType>;            \
  using GraphPositionMap =                                                    \
      bagl::property_map_t<MotionGraphType, PointType BasicVertexProp::*>;    \
  using SpacePartType = dvp_tree<Vertex, SuperSpace, GraphPositionMap, ARITY, \
                                 random_vp_chooser, TREE_STORAGE>;            \
  SpacePartType space_part1(motion_graph1, sup_space_ptr,                     \
                            get(&BasicVertexProp::position, motion_graph1));  \
  SpacePartType space_part2(motion_graph2, sup_space_ptr,                     \
                            get(&BasicVertexProp::position, motion_graph2));  \
                                                                              \
  using NNFinderType = multi_dvp_tree_search<MotionGraphType, SpacePartType>; \
  NNFinderType nn_finder;                                                     \
  nn_finder.graph_tree_map[&motion_graph1] = &space_part1;                    \
  nn_finder.graph_tree_map[&motion_graph2] = &space_part2;                    \
                                                                              \
  type_erased_knn_synchro<MotionGraphType, NNFinderType> NN_synchro(          \
      nn_finder);                                                             \
  vis.m_nn_synchro = &NN_synchro;

#define RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(ARITY, TREE_STORAGE)         \
  using ALTGraph =                                                             \
      dvp_adjacency_list<VertexProp, EdgeProp, SuperSpace, PositionMap, ARITY, \
                         random_vp_chooser, TREE_STORAGE, bagl::vec_s,         \
                         DirectionalityTag, bagl::list_s>;                     \
  using MotionGraphType = typename ALTGraph::adj_list_type;                    \
                                                                               \
  ALTGraph space_part1(sup_space_ptr, pos_map);                                \
  ALTGraph space_part2(sup_space_ptr, pos_map);                                \
                                                                               \
  MotionGraphType motion_graph1 = space_part1.get_adjacency_list();            \
  MotionGraphType motion_graph2 = space_part2.get_adjacency_list();            \
                                                                               \
  vis.m_start_node = std::any(create_root(motion_graph1, vp_start));           \
  vis.m_goal_node = std::any(create_root(motion_graph2, vp_goal));             \
                                                                               \
  using NNFinderType = multi_dvp_tree_search<MotionGraphType, ALTGraph>;       \
  NNFinderType nn_finder;                                                      \
  nn_finder.graph_tree_map[&motion_graph1] = &space_part1;                     \
  nn_finder.graph_tree_map[&motion_graph2] = &space_part2;                     \
                                                                               \
  any_knn_synchro NN_synchro;                                                  \
  vis.m_nn_synchro = &NN_synchro;

#define RK_RRT_PLANNER_CALL_BIRRT_FUNCTION                        \
  ReaK::graph::generate_bidirectional_rrt(                        \
      motion_graph1, motion_graph2, *sup_space_ptr, vis, pos_map, \
      get(random_sampler, *sup_space_ptr), nn_finder);

    if ((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) ==
        ADJ_LIST_MOTION_GRAPH) {

      using MotionGraphType =
          bagl::adjacency_list<bagl::vec_s, bagl::list_s,
                               DirectionalityTag, VertexProp, EdgeProp>;

      MotionGraphType motion_graph1;
      MotionGraphType motion_graph2;

      vis.m_start_node = std::any(create_root(motion_graph1, vp_start));
      vis.m_goal_node = std::any(create_root(motion_graph2, vp_goal));

      if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
          LINEAR_SEARCH_KNN) {

        any_knn_synchro NN_synchro;
        vis.m_nn_synchro = &NN_synchro;
        linear_neighbor_search<MotionGraphType> nn_finder;

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_BF2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(
            2, bagl::bfl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_BF4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(
            4, bagl::bfl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(
            2, bagl::vebl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(
            4, bagl::vebl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

#endif
      }

#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT

    } else if ((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) ==
               DVP_ADJ_LIST_MOTION_GRAPH) {

      if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
          DVP_BF2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(
            2, bagl::bfl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_BF4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(
            4, bagl::bfl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(
            2, bagl::vebl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(
            4, bagl::vebl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

#endif
      }

#endif
    }

#undef RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO
#undef RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO
#undef RK_RRT_PLANNER_CALL_BIRRT_FUNCTION
  }
}

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_RRT_PATH_PLANNER_INC_

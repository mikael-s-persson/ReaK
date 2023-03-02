/**
 * \file rrt_path_planner.tpp
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

#ifndef REAK_RRT_PATH_PLANNER_TPP
#define REAK_RRT_PATH_PLANNER_TPP

#include "rrt_path_planner.hpp"

#include <ReaK/planning/graph_alg/rr_tree.hpp>
#include "motion_graph_structures.hpp"

#include "metric_space_search.hpp"
#include "topological_search.hpp"

// BGL-Extra includes:
#include <boost/graph/more_property_maps.hpp>
#include <boost/graph/more_property_tags.hpp>
#include <boost/graph/tree_adaptor.hpp>

#include "any_motion_graphs.hpp"
#include "p2p_planning_query.hpp"
#include "path_planner_options.hpp"
#include "planning_visitors.hpp"

namespace ReaK::pp {

template <typename FreeSpaceType>
void rrt_planner<FreeSpaceType>::solve_planning_query(
    planning_query<FreeSpaceType>& aQuery) {

  this->reset_internal_state();

  using SuperSpace = typename subspace_traits<FreeSpaceType>::super_space_type;
  using PointType = topology_point_type_t<SuperSpace>;

  using VertexProp = mg_vertex_data<FreeSpaceType>;
  using EdgeProp = mg_edge_data<FreeSpaceType>;

  using DirectionalityTag = typename motion_segment_directionality<FreeSpaceType>::type;

  using BasicVertexProp = mg_vertex_data<FreeSpaceType>;

  using PositionMap = boost::data_member_property_map<PointType, VertexProp>;
  PositionMap pos_map = PositionMap(&VertexProp::position);

  std::shared_ptr<const SuperSpace> sup_space_ptr(
      &(this->m_space->get_super_space()), null_deleter());

  planning_visitor<FreeSpaceType> vis(this, &aQuery);

  VertexProp vp_start;
  vp_start.position = aQuery.get_start_position();

  if ((this->m_planning_method_flags & PLANNING_DIRECTIONALITY_MASK) ==
      UNIDIRECTIONAL_PLANNING) {

#define RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(ARITY, TREE_STORAGE)             \
  using Vertex = graph::graph_vertex_t<MotionGraphType>;                       \
  using GraphPositionMap = typename boost::property_map<                       \
      MotionGraphType, PointType BasicVertexProp::*>::type;                    \
  using SpacePartType = dvp_tree<Vertex, SuperSpace, GraphPositionMap, ARITY,  \
      random_vp_chooser, TREE_STORAGE>;                                        \
  SpacePartType space_part(motion_graph, sup_space_ptr,                        \
                           get(&BasicVertexProp::position, motion_graph));     \
                                                                               \
  using NNFinderType = multi_dvp_tree_search<MotionGraphType, SpacePartType>;  \
  NNFinderType nn_finder;                                                      \
  nn_finder.graph_tree_map[&motion_graph] = &space_part;                       \
                                                                               \
  type_erased_knn_synchro<MotionGraphType, NNFinderType> NN_synchro(           \
      nn_finder);                                                              \
  vis.m_nn_synchro = &NN_synchro;

#define RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(ARITY, TREE_STORAGE)           \
  using ALTGraph = dvp_adjacency_list<VertexProp, EdgeProp, SuperSpace,      \
      PositionMap, ARITY, random_vp_chooser, TREE_STORAGE, boost::vecBC,     \
      DirectionalityTag, boost::listBC>;                                     \
  using MotionGraphType = typename ALTGraph::adj_list_type;                  \
                                                                             \
  ALTGraph space_part(sup_space_ptr, pos_map);                               \
  MotionGraphType motion_graph = space_part.get_adjacency_list();            \
  vis.m_start_node = std::any(create_root(vp_start, motion_graph));          \
                                                                             \
  using NNFinderType = multi_dvp_tree_search<MotionGraphType, ALTGraph>;     \
  NNFinderType nn_finder;                                                    \
  nn_finder.graph_tree_map[&motion_graph] = &space_part;                     \
                                                                             \
  any_knn_synchro NN_synchro;                                                \
  vis.m_nn_synchro = &NN_synchro;

#define RK_RRT_PLANNER_CALL_RRT_FUNCTION                                \
  ReaK::graph::generate_rrt(motion_graph, *sup_space_ptr, vis, pos_map, \
                            get(random_sampler, *sup_space_ptr), nn_finder);

    if ((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) ==
        ADJ_LIST_MOTION_GRAPH) {

      using MotionGraphType = boost::adjacency_list_BC<
          boost::vecBC, boost::listBC, DirectionalityTag, VertexProp, EdgeProp>;

      MotionGraphType motion_graph;
      vis.m_start_node = std::any(create_root(vp_start, motion_graph));

      if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
          LINEAR_SEARCH_KNN) {

        any_knn_synchro NN_synchro;
        vis.m_nn_synchro = &NN_synchro;
        linear_neighbor_search<MotionGraphType> nn_finder;

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_BF2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(2,
                                              boost::bfl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_BF4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(4,
                                              boost::bfl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(2,
                                              boost::vebl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO(4,
                                              boost::vebl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

#endif
      }

#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT

    } else if ((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) ==
               DVP_ADJ_LIST_MOTION_GRAPH) {

      if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
          DVP_BF2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(2,
                                              boost::bfl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_BF4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(4,
                                              boost::bfl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(2,
                                              boost::vebl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO(4,
                                              boost::vebl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_RRT_FUNCTION

#endif
      }

#endif
    }

#undef RK_RRT_PLANNER_SETUP_DVP_TREE_SYNCHRO
#undef RK_RRT_PLANNER_SETUP_ALT_TREE_SYNCHRO
#undef RK_RRT_PLANNER_CALL_RRT_FUNCTION

  } else {
    path_planning_p2p_query<FreeSpaceType>* p2p_query_ptr =
        reinterpret_cast<path_planning_p2p_query<FreeSpaceType>*>(aQuery.castTo(
            path_planning_p2p_query<FreeSpaceType>::getStaticObjectType()));
    if (p2p_query_ptr == nullptr) {
      return;
    }

    VertexProp vp_goal;
    vp_goal.position = p2p_query_ptr->goal_pos;

#define RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(ARITY, TREE_STORAGE)         \
  using Vertex = graph::graph_vertex_t<MotionGraphType>;                       \
  using GraphPositionMap = typename boost::property_map<MotionGraphType,       \
      PointType BasicVertexProp::*>::type;                                     \
  using SpacePartType = dvp_tree<Vertex, SuperSpace, GraphPositionMap, ARITY,  \
      random_vp_chooser, TREE_STORAGE>;                                        \
  SpacePartType space_part1(motion_graph1, sup_space_ptr,                      \
                            get(&BasicVertexProp::position, motion_graph1));   \
  SpacePartType space_part2(motion_graph2, sup_space_ptr,                      \
                            get(&BasicVertexProp::position, motion_graph2));   \
                                                                               \
  using NNFinderType = multi_dvp_tree_search<MotionGraphType, SpacePartType>;  \
  NNFinderType nn_finder;                                                      \
  nn_finder.graph_tree_map[&motion_graph1] = &space_part1;                     \
  nn_finder.graph_tree_map[&motion_graph2] = &space_part2;                     \
                                                                               \
  type_erased_knn_synchro<MotionGraphType, NNFinderType> NN_synchro(           \
      nn_finder);                                                              \
  vis.m_nn_synchro = &NN_synchro;

#define RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(ARITY, TREE_STORAGE)       \
  using ALTGraph = dvp_adjacency_list<VertexProp, EdgeProp, SuperSpace,      \
      PositionMap, ARITY, random_vp_chooser, TREE_STORAGE, boost::vecBC,     \
      DirectionalityTag, boost::listBC>;                                     \
  using MotionGraphType = typename ALTGraph::adj_list_type;                  \
                                                                             \
  ALTGraph space_part1(sup_space_ptr, pos_map);                              \
  ALTGraph space_part2(sup_space_ptr, pos_map);                              \
                                                                             \
  MotionGraphType motion_graph1 = space_part1.get_adjacency_list();          \
  MotionGraphType motion_graph2 = space_part2.get_adjacency_list();          \
                                                                             \
  vis.m_start_node = std::any(create_root(vp_start, motion_graph1));         \
  vis.m_goal_node = std::any(create_root(vp_goal, motion_graph2));           \
                                                                             \
  using NNFinderType = multi_dvp_tree_search<MotionGraphType, ALTGraph>;     \
  NNFinderType nn_finder;                                                    \
  nn_finder.graph_tree_map[&motion_graph1] = &space_part1;                   \
  nn_finder.graph_tree_map[&motion_graph2] = &space_part2;                   \
                                                                             \
  any_knn_synchro NN_synchro;                                                \
  vis.m_nn_synchro = &NN_synchro;

#define RK_RRT_PLANNER_CALL_BIRRT_FUNCTION                        \
  ReaK::graph::generate_bidirectional_rrt(                        \
      motion_graph1, motion_graph2, *sup_space_ptr, vis, pos_map, \
      get(random_sampler, *sup_space_ptr), nn_finder);

    if ((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) ==
        ADJ_LIST_MOTION_GRAPH) {

      using MotionGraphType = boost::adjacency_list_BC<
          boost::vecBC, boost::listBC, DirectionalityTag, VertexProp, EdgeProp>;

      MotionGraphType motion_graph1;
      MotionGraphType motion_graph2;

      vis.m_start_node = std::any(create_root(vp_start, motion_graph1));
      vis.m_goal_node = std::any(create_root(vp_goal, motion_graph2));

      if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
          LINEAR_SEARCH_KNN) {

        any_knn_synchro NN_synchro;
        vis.m_nn_synchro = &NN_synchro;
        linear_neighbor_search<MotionGraphType> nn_finder;

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_BF2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(
            2, boost::bfl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_BF4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(
            4, boost::bfl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(
            2, boost::vebl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_DVP_TREE_SYNCHRO(
            4, boost::vebl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

#endif
      }

#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT

    } else if ((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) ==
               DVP_ADJ_LIST_MOTION_GRAPH) {

      if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
          DVP_BF2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(
            2, boost::bfl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_BF4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(
            4, boost::bfl_d_ary_tree_storage<4>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB2_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(
            2, boost::vebl_d_ary_tree_storage<2>)

        RK_RRT_PLANNER_CALL_BIRRT_FUNCTION

      } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
                 DVP_COB4_TREE_KNN) {

        RK_RRT_PLANNER_SETUP_TWO_ALT_TREE_SYNCHRO(
            4, boost::vebl_d_ary_tree_storage<4>)

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

#endif

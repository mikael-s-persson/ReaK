/**
 * \file rrtstar_path_planner.inc
 *
 * This library contains template definitions of a class to solve path planning problems using the
 * Rapidly-exploring Random Tree Star (RRT*) algorithm (or one of its variants).
 * Given a C_free (configuration space restricted to non-colliding points) and a
 * result reporting policy, this class will probabilistically construct a motion-graph
 * that will connect a starting point and a goal point with a path through C-free
 * that is as close as possible to the optimal path in terms of distance.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date August 2012
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

#ifndef REAK_PLANNING_PATH_PLANNING_RRTSTAR_PATH_PLANNER_INC_
#define REAK_PLANNING_PATH_PLANNING_RRTSTAR_PATH_PLANNER_INC_

#include "ReaK/planning/path_planning/rrtstar_path_planner.h"

#include "ReaK/planning/graph_alg/rrt_star.h"

#include "ReaK/planning/path_planning/motion_graph_structures.h"

#include "bagl/more_property_maps.h"
#include "bagl/tree_adaptor.h"

#include <type_traits>

#include "ReaK/planning/path_planning/metric_space_search.h"
#include "ReaK/planning/path_planning/topological_search.h"

#include "ReaK/planning/graph_alg/neighborhood_functors.h"
#include "ReaK/planning/path_planning/any_motion_graphs.h"
#include "ReaK/planning/path_planning/p2p_planning_query.h"
#include "ReaK/planning/path_planning/path_planner_options.h"
#include "ReaK/planning/path_planning/planning_visitors.h"

namespace ReaK::pp {

template <SubSpace FreeSpaceType, bool IsBidirPlanner = false>
struct rrtstar_bundle_factory {

  static constexpr bool is_bidir =
      IsBidirPlanner && is_reversible_space_v<FreeSpaceType>;

  using super_space_type =
      typename subspace_traits<FreeSpaceType>::super_space_type;
  using point_type = topology_point_type_t<super_space_type>;

  using vertex_prop =
      std::conditional_t<is_bidir, bidir_optimal_mg_vertex<FreeSpaceType>,
                         optimal_mg_vertex<FreeSpaceType>>;
  using basic_vertex_prop = mg_vertex_data<FreeSpaceType>;
  using edge_prop = optimal_mg_edge<FreeSpaceType>;

  using directionality_tag =
      typename motion_segment_directionality<FreeSpaceType>::type;

  using visitor_type = planning_visitor<FreeSpaceType>;

  using position_map = bagl::data_member_property_map<point_type, vertex_prop>;
  using distance_map = bagl::data_member_property_map<double, vertex_prop>;
  using predecessor_map =
      bagl::data_member_property_map<std::size_t, vertex_prop>;
  using fwd_distance_map =
      std::conditional_t<is_bidir,
                         bagl::data_member_property_map<double, vertex_prop>,
                         ReaK::graph::infinite_double_value_prop_map>;

  template <typename Graph>
  using successor_map_t = std::conditional_t<
      is_bidir, bagl::data_member_property_map<std::size_t, vertex_prop>,
      ReaK::graph::null_vertex_prop_map<Graph>>;

  using weight_map = bagl::data_member_property_map<double, edge_prop>;

  struct ls_motion_graph {
    using type =
        bagl::adjacency_list<bagl::vec_s, bagl::pool_s,
                             directionality_tag, vertex_prop, edge_prop>;
    using vertex_type = bagl::graph_vertex_descriptor_t<type>;

    struct space_part_type {
      template <typename GraphPositionMap>
      space_part_type(const type& /*unused*/,
                      const std::shared_ptr<const super_space_type>& /*unused*/,
                      GraphPositionMap /*unused*/) {}
    };

    using nn_finder_type = linear_neighbor_search<type>;
    static nn_finder_type get_nn_finder(type& /*unused*/,
                                        space_part_type& /*unused*/) {
      return nn_finder_type();
    }

    using nn_synchro_type = any_knn_synchro;
    static nn_synchro_type get_nn_synchro(nn_finder_type& /*unused*/) {
      return {};
    }

    static type get_motion_graph() { return type(); }
  };

  template <unsigned int Arity, typename TreeStorageTag>
  struct dvp_motion_graph {
    using type =
        bagl::adjacency_list<bagl::vec_s, bagl::pool_s,
                                 directionality_tag, vertex_prop, edge_prop>;
    using vertex_type = bagl::graph_vertex_descriptor_t<type>;

    using graph_position_map =
        bagl::property_map_t<type, point_type basic_vertex_prop::*>;
    using space_part_type =
        dvp_tree<vertex_type, super_space_type, graph_position_map, Arity,
                 random_vp_chooser, TreeStorageTag>;
    static space_part_type get_space_part(
        type& mg, std::shared_ptr<const super_space_type> s_ptr) {
      return space_part_type(mg, s_ptr, get(&basic_vertex_prop::position, mg));
    }

    using nn_finder_type = multi_dvp_tree_search<type, space_part_type>;
    static nn_finder_type get_nn_finder(type& mg, space_part_type& space_part) {
      nn_finder_type nn_finder;
      nn_finder.graph_tree_map[&mg] = &space_part;
      return nn_finder;
    }

    using nn_synchro_type = type_erased_knn_synchro<type, nn_finder_type>;
    static nn_synchro_type get_nn_synchro(nn_finder_type& nn_finder) {
      return nn_synchro_type(nn_finder);
    }

    static type get_motion_graph() { return type(); }
  };

#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT

  template <unsigned int Arity, typename TreeStorageTag>
  struct alt_motion_graph {
    using alt_graph_type =
        dvp_adjacency_list<vertex_prop, edge_prop, super_space_type,
                           position_map, Arity, random_vp_chooser,
                           TreeStorageTag, bagl::vec_s, directionality_tag,
                           bagl::list_s>;
    using type = typename alt_graph_type::adj_list_type;
    using space_part_type = alt_graph_type;
    using vertex_type = bagl::graph_vertex_descriptor_t<type>;

    static space_part_type get_space_part(
        std::shared_ptr<const super_space_type> s_ptr) {
      return space_part_type(s_ptr, position_map(&vertex_prop::position));
    }

    using nn_finder_type = multi_dvp_tree_search<type, space_part_type>;
    static nn_finder_type get_nn_finder(type& mg, space_part_type& space_part) {
      nn_finder_type nn_finder;
      nn_finder.graph_tree_map[&mg] = &space_part;
      return nn_finder;
    }

    using nn_synchro_type = any_knn_synchro;
    static nn_synchro_type get_nn_synchro(nn_finder_type&) {
      return nn_synchro_type();
    }

    static type get_motion_graph(space_part_type& space_part) {
      return space_part.get_adjacency_list();
    }
  };

#endif

  template <typename MotionGraphType>
  static void init_motion_graph(MotionGraphType& motion_graph,
                                visitor_type& vis,
                                planning_query<FreeSpaceType>& query) {
    using Vertex = bagl::graph_vertex_descriptor_t<MotionGraphType>;
    vertex_prop vp_start;
    vp_start.position = query.get_start_position();
    Vertex start_node = add_vertex(motion_graph, std::move(vp_start));
    vis.m_start_node = std::any(start_node);
    auto* p2p_query_ptr =
        reinterpret_cast<path_planning_p2p_query<FreeSpaceType>*>(query.cast_to(
            path_planning_p2p_query<FreeSpaceType>::get_static_object_type()));
    if (p2p_query_ptr) {
      vertex_prop vp_goal;
      vp_goal.position = p2p_query_ptr->goal_pos;
      Vertex goal_node = add_vertex(motion_graph, std::move(vp_goal));
      vis.m_goal_node = std::any(goal_node);
      vis.initialize_vertex(goal_node, motion_graph);
    }
    vis.initialize_vertex(start_node, motion_graph);
  }

  template <typename VertexProp>
  static fwd_distance_map dispatched_make_fwd_dist_map() {
    if constexpr (std::is_convertible_v<
                      VertexProp*, bidir_optimal_mg_vertex<FreeSpaceType>*>) {
      return fwd_distance_map(&VertexProp::fwd_distance_accum);
    } else {
      return fwd_distance_map();
    }
  }

  template <typename MotionGraphType, typename VertexProp>
  static successor_map_t<MotionGraphType> dispatched_make_succ_map() {
    if constexpr (std::is_convertible_v<
                      VertexProp*, bidir_optimal_mg_vertex<FreeSpaceType>*>) {
      return successor_map_t<MotionGraphType>(&VertexProp::successor);
    } else {
      return successor_map_t<MotionGraphType>();
    }
  }

  template <typename MotionGraphType, typename NcSelector>
  static ReaK::graph::rrtstar_bundle<
      MotionGraphType, super_space_type, visitor_type, NcSelector, position_map,
      weight_map, distance_map, predecessor_map, fwd_distance_map,
      successor_map_t<MotionGraphType>>
  make_bundle(MotionGraphType& motion_graph, visitor_type& vis,
              NcSelector nc_selector,
              std::shared_ptr<const super_space_type> s_ptr) {
    using Vertex = bagl::graph_vertex_descriptor_t<MotionGraphType>;
    return ReaK::graph::make_rrtstar_bundle(
        motion_graph, std::any_cast<Vertex>(vis.m_start_node),
        (!vis.m_goal_node.has_value()
             ? bagl::graph_traits<MotionGraphType>::null_vertex()
             : std::any_cast<Vertex>(vis.m_goal_node)),
        *s_ptr, vis, nc_selector, position_map(&vertex_prop::position),
        weight_map(&edge_prop::weight),
        distance_map(&vertex_prop::distance_accum),
        predecessor_map(&vertex_prop::predecessor),
        dispatched_make_fwd_dist_map<vertex_prop>(),
        dispatched_make_succ_map<MotionGraphType, vertex_prop>());
  }

  template <typename MotionGraphType, typename NcSelector>
  static void make_call_to_planner(
      MotionGraphType& motion_graph, visitor_type& vis, NcSelector nc_selector,
      std::shared_ptr<const super_space_type> s_ptr, std::size_t method_flags) {
    if constexpr (is_bidir) {
      if (method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) {
        ReaK::graph::generate_bnb_rrt_star_bidir(
            make_bundle(motion_graph, vis, nc_selector, s_ptr),
            get(random_sampler, *s_ptr));
      } else { /* assume nominal method only. */
        ReaK::graph::generate_rrt_star_bidir(
            make_bundle(motion_graph, vis, nc_selector, s_ptr),
            get(random_sampler, *s_ptr));
      }
    } else {
      if (method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) {
        ReaK::graph::generate_bnb_rrt_star(
            make_bundle(motion_graph, vis, nc_selector, s_ptr),
            get(random_sampler, *s_ptr));
      } else { /* assume nominal method only. */
        ReaK::graph::generate_rrt_star(
            make_bundle(motion_graph, vis, nc_selector, s_ptr),
            get(random_sampler, *s_ptr));
      }
    }
  }
};

template <SubSpace FreeSpaceType>
template <typename RRTStarFactory>
void rrtstar_planner<FreeSpaceType>::solve_planning_query_impl(
    planning_query<FreeSpaceType>& aQuery) {

  this->reset_internal_state();

  auto space_dim = double(this->get_space_dimensionality());
  double space_Lc = aQuery.get_heuristic_to_goal(aQuery.get_start_position());

  using SuperSpace = typename subspace_traits<FreeSpaceType>::super_space_type;
  std::shared_ptr<const SuperSpace> sup_space_ptr(
      &(this->m_space->get_super_space()), null_deleter());

  // Some MACROs to reduce the size of the code below.

#define RK_RRTSTAR_PLANNER_SETUP_LS_OR_DVP_SUPPORT_STRUCTURES                  \
  using MotionGraphType = typename MGFactory::type;                            \
  MotionGraphType motion_graph = MGFactory::get_motion_graph();                \
  RRTStarFactory::init_motion_graph(motion_graph, vis, aQuery);                \
                                                                               \
  using SpacePartType = typename MGFactory::space_part_type;                   \
  using BasicVProp = typename RRTStarFactory::basic_vertex_prop;               \
  SpacePartType space_part(motion_graph, sup_space_ptr,                        \
                           get(&BasicVProp::position, motion_graph));          \
                                                                               \
  using NNFinderType = typename MGFactory::nn_finder_type;                     \
  NNFinderType nn_finder = MGFactory::get_nn_finder(motion_graph, space_part); \
                                                                               \
  ReaK::graph::star_neighborhood<NNFinderType> nc_selector(                    \
      nn_finder, space_dim, 3.0 * space_Lc);                                   \
                                                                               \
  using NNSynchroType = typename MGFactory::nn_synchro_type;                   \
  NNSynchroType NN_synchro = MGFactory::get_nn_synchro(nn_finder);             \
  vis.m_nn_synchro = &NN_synchro;

#define RK_RRTSTAR_PLANNER_SETUP_ALT_SUPPORT_STRUCTURES                        \
  using MotionGraphType = typename MGFactory::type;                            \
  using SpacePartType = typename MGFactory::space_part_type;                   \
  using PosMap = typename RRTStarFactory::position_map;                        \
  using VProp = typename RRTStarFactory::vertex_prop;                          \
  SpacePartType space_part(sup_space_ptr, PosMap(&VProp::position));           \
                                                                               \
  MotionGraphType motion_graph = MGFactory::get_motion_graph(space_part);      \
                                                                               \
  using NNFinderType = typename MGFactory::nn_finder_type;                     \
  NNFinderType nn_finder = MGFactory::get_nn_finder(motion_graph, space_part); \
                                                                               \
  ReaK::graph::star_neighborhood<NNFinderType> nc_selector(                    \
      nn_finder, space_dim, 3.0 * space_Lc);                                   \
                                                                               \
  using NNSynchroType = typename MGFactory::nn_synchro_type;                   \
  NNSynchroType NN_synchro = MGFactory::get_nn_synchro(nn_finder);             \
  vis.m_nn_synchro = &NN_synchro;                                              \
                                                                               \
  RRTStarFactory::init_motion_graph(motion_graph, vis, aQuery);

  using VisitorType = typename RRTStarFactory::visitor_type;
  VisitorType vis(this, &aQuery);

  if ((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) ==
      ADJ_LIST_MOTION_GRAPH) {

    if ((this->m_data_structure_flags & KNN_METHOD_MASK) == LINEAR_SEARCH_KNN) {

      using MGFactory = typename RRTStarFactory::ls_motion_graph;

      RK_RRTSTAR_PLANNER_SETUP_LS_OR_DVP_SUPPORT_STRUCTURES

      RRTStarFactory::make_call_to_planner(motion_graph, vis, nc_selector,
                                           sup_space_ptr,
                                           this->m_planning_method_flags);

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_BF2_TREE_KNN) {

      using MGFactory = typename RRTStarFactory::template dvp_motion_graph<
          2, bagl::bfl_d_ary_tree_storage<2>>;

      RK_RRTSTAR_PLANNER_SETUP_LS_OR_DVP_SUPPORT_STRUCTURES

      RRTStarFactory::make_call_to_planner(motion_graph, vis, nc_selector,
                                           sup_space_ptr,
                                           this->m_planning_method_flags);

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_BF4_TREE_KNN) {

      using MGFactory = typename RRTStarFactory::template dvp_motion_graph<
          4, bagl::bfl_d_ary_tree_storage<4>>;

      RK_RRTSTAR_PLANNER_SETUP_LS_OR_DVP_SUPPORT_STRUCTURES

      RRTStarFactory::make_call_to_planner(motion_graph, vis, nc_selector,
                                           sup_space_ptr,
                                           this->m_planning_method_flags);

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_COB2_TREE_KNN) {

      using MGFactory = typename RRTStarFactory::template dvp_motion_graph<
          2, bagl::vebl_d_ary_tree_storage<2>>;

      RK_RRTSTAR_PLANNER_SETUP_LS_OR_DVP_SUPPORT_STRUCTURES

      RRTStarFactory::make_call_to_planner(motion_graph, vis, nc_selector,
                                           sup_space_ptr,
                                           this->m_planning_method_flags);

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_COB4_TREE_KNN) {

      using MGFactory = typename RRTStarFactory::template dvp_motion_graph<
          4, bagl::vebl_d_ary_tree_storage<4>>;

      RK_RRTSTAR_PLANNER_SETUP_LS_OR_DVP_SUPPORT_STRUCTURES

      RRTStarFactory::make_call_to_planner(motion_graph, vis, nc_selector,
                                           sup_space_ptr,
                                           this->m_planning_method_flags);

#endif
    }

#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT

  } else if ((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) ==
             DVP_ADJ_LIST_MOTION_GRAPH) {

    if ((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {

      using MGFactory = typename RRTStarFactory::template alt_motion_graph<
          2, bagl::bfl_d_ary_tree_storage<2>>;

      RK_RRTSTAR_PLANNER_SETUP_ALT_SUPPORT_STRUCTURES

      RRTStarFactory::make_call_to_planner(motion_graph, vis, nc_selector,
                                           sup_space_ptr,
                                           this->m_planning_method_flags);

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_BF4_TREE_KNN) {

      using MGFactory = typename RRTStarFactory::template alt_motion_graph<
          4, bagl::bfl_d_ary_tree_storage<4>>;

      RK_RRTSTAR_PLANNER_SETUP_ALT_SUPPORT_STRUCTURES

      RRTStarFactory::make_call_to_planner(motion_graph, vis, nc_selector,
                                           sup_space_ptr,
                                           this->m_planning_method_flags);

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_COB2_TREE_KNN) {

      using MGFactory = typename RRTStarFactory::template alt_motion_graph<
          2, bagl::vebl_d_ary_tree_storage<2>>;

      RK_RRTSTAR_PLANNER_SETUP_ALT_SUPPORT_STRUCTURES

      RRTStarFactory::make_call_to_planner(motion_graph, vis, nc_selector,
                                           sup_space_ptr,
                                           this->m_planning_method_flags);

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_COB4_TREE_KNN) {

      using MGFactory = typename RRTStarFactory::template alt_motion_graph<
          4, bagl::vebl_d_ary_tree_storage<4>>;

      RK_RRTSTAR_PLANNER_SETUP_ALT_SUPPORT_STRUCTURES

      RRTStarFactory::make_call_to_planner(motion_graph, vis, nc_selector,
                                           sup_space_ptr,
                                           this->m_planning_method_flags);

#endif
    }

#endif
  }

#undef RK_RRTSTAR_PLANNER_SETUP_LS_OR_DVP_SUPPORT_STRUCTURES
#undef RK_RRTSTAR_PLANNER_SETUP_ALT_SUPPORT_STRUCTURES
}

template <SubSpace FreeSpaceType>
void rrtstar_planner<FreeSpaceType>::solve_planning_query(
    planning_query<FreeSpaceType>& aQuery) {
  if ((this->m_planning_method_flags & PLANNING_DIRECTIONALITY_MASK) ==
      UNIDIRECTIONAL_PLANNING) {
    this->solve_planning_query_impl<
        rrtstar_bundle_factory<FreeSpaceType, false>>(aQuery);
  } else { /* NOTE: Bi-directional version: */
    this->solve_planning_query_impl<
        rrtstar_bundle_factory<FreeSpaceType, true>>(aQuery);
  }
}

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_RRTSTAR_PATH_PLANNER_INC_

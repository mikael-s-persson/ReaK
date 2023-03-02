/**
 * \file sbastar_path_planner.tpp
 *
 * This library contains template definitions of a class to solve path planning problems using the
 * Sampling-based A* algorithm (or one of its variants). Given a C_free (configuration space
 * restricted to non-colliding points) and a result reporting policy, this class
 * will probabilistically construct a motion-graph that will connect a starting point
 * and a goal point with a path through C-free that is as close as possible to the
 * optimal path in terms of distance. The planner uses a selectable variant of the
 * Sampling-based A* (SBA*) algorithm, including the basic version, the SBA*-RRT*
 * alternating algorithm, and the Anytime SBA* algorithm. In all cases, collision
 * checking and connectivity can be either full or lazy (and pruned) to either construct
 * a full-connectivity graph containing only collision-free edges, or a single-query motion-tree
 * that includes only optimal edges (whose collisions are checked lazily).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2013
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

#ifndef REAK_SBASTAR_PATH_PLANNER_TPP
#define REAK_SBASTAR_PATH_PLANNER_TPP

#include "sbastar_path_planner.hpp"

#include <ReaK/planning/graph_alg/anytime_sbastar.hpp>
#include <ReaK/planning/graph_alg/lazy_sbastar.hpp>
#include <ReaK/planning/graph_alg/sbastar_rrtstar.hpp>

#include "motion_graph_structures.hpp"

// BGL-Extra includes:
#include <boost/graph/more_property_maps.hpp>
#include <boost/graph/more_property_tags.hpp>

#include <type_traits>

#include "metric_space_search.hpp"
#include "topological_search.hpp"

#include <ReaK/planning/graph_alg/neighborhood_functors.hpp>
#include "any_motion_graphs.hpp"
#include "density_plan_visitors.hpp"
#include "p2p_planning_query.hpp"
#include "path_planner_options.hpp"

namespace ReaK::pp {

template <typename FreeSpaceType, bool IsBidirPlanner = false>
struct sbastar_bundle_factory {

  static constexpr bool is_bidir = IsBidirPlanner && is_reversible_space_v<FreeSpaceType>;

  using super_space_type = typename subspace_traits<FreeSpaceType>::super_space_type;
  using point_type = topology_point_type_t<super_space_type>;

  using vertex_prop = std::conditional_t<is_bidir,
                                         recursive_dense_mg_vertex<bidir_astar_mg_vertex<FreeSpaceType>>,
                                         recursive_dense_mg_vertex<astar_mg_vertex<FreeSpaceType>>>;
  using basic_vertex_prop = mg_vertex_data<FreeSpaceType>;
  using edge_prop = optimal_mg_edge<FreeSpaceType>;

  using directionality_tag = typename motion_segment_directionality<FreeSpaceType>::type;

  using visitor_type = density_plan_visitor<FreeSpaceType, sbastar_density_calculator>;

  using position_map = boost::data_member_property_map<point_type, vertex_prop>;
  using density_map = boost::data_member_property_map<double, vertex_prop>;
  using constriction_map = boost::data_member_property_map<double, vertex_prop>;
  using distance_map = boost::data_member_property_map<double, vertex_prop>;
  using fwd_distance_map = boost::data_member_property_map<double, vertex_prop>;
  using predecessor_map = boost::data_member_property_map<std::size_t, vertex_prop>;

  template <typename Graph>
  using successor_map_t = std::conditional_t<is_bidir,
                                             boost::data_member_property_map<std::size_t, vertex_prop>,
                                             ReaK::graph::detail::null_vertex_prop_map<Graph>>;

  using weight_map = boost::data_member_property_map<double, edge_prop>;

  struct ls_motion_graph {
    using type = boost::adjacency_list_BC<boost::vecBC, boost::poolBC, directionality_tag,
        vertex_prop, edge_prop>;
    using vertex_type = graph::graph_vertex_t<type>;

    struct space_part_type {
      template <typename GraphPositionMap>
      space_part_type(const type&,
                      const std::shared_ptr<const super_space_type>&,
                      GraphPositionMap){}
    };

    using nn_finder_type = linear_neighbor_search<type>;
    static nn_finder_type get_nn_finder(type&, space_part_type&) {
      return nn_finder_type();
    }

    using nn_synchro_type = any_knn_synchro;
    static nn_synchro_type get_nn_synchro(nn_finder_type&) {
      return nn_synchro_type();
    }

    static type get_motion_graph() { return type(); }
  };

  template <unsigned int Arity, typename TreeStorageTag>
  struct dvp_motion_graph {
    using type = boost::adjacency_list_BC<boost::vecBC, boost::poolBC, directionality_tag,
        vertex_prop, edge_prop>;
    using vertex_type = graph::graph_vertex_t<type>;

    using graph_position_map = typename boost::property_map<type,
        point_type basic_vertex_prop::*>::type;
    using space_part_type = dvp_tree<vertex_type, super_space_type, graph_position_map,
        Arity, random_vp_chooser, TreeStorageTag>;
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
    using alt_graph_type = dvp_adjacency_list<vertex_prop, edge_prop, super_space_type,
        position_map, Arity, random_vp_chooser, TreeStorageTag, boost::vecBC,
        directionality_tag, boost::listBC>;
    using type = typename alt_graph_type::adj_list_type;
    using space_part_type = alt_graph_type;
    using vertex_type = graph::graph_vertex_t<type>;

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
    using Vertex = graph::graph_vertex_t<MotionGraphType>;
    vertex_prop vp_start;
    vp_start.position = query.get_start_position();
    Vertex start_node = add_vertex(vp_start, motion_graph);
    vis.m_start_node = std::any(start_node);
    path_planning_p2p_query<FreeSpaceType>* p2p_query_ptr =
        reinterpret_cast<path_planning_p2p_query<FreeSpaceType>*>(query.castTo(
            path_planning_p2p_query<FreeSpaceType>::getStaticObjectType()));
    if (p2p_query_ptr) {
      vertex_prop vp_goal;
      vp_goal.position = p2p_query_ptr->goal_pos;
      Vertex goal_node = add_vertex(vp_goal, motion_graph);
      vis.m_goal_node = std::any(goal_node);
      vis.initialize_vertex(goal_node, motion_graph);
    }
    vis.initialize_vertex(start_node, motion_graph);
  }

  template <typename VertexProp>
  static auto dispatched_make_fwd_dist_map() {
    if constexpr (std::is_convertible_v<VertexProp*, bidir_optimal_mg_vertex<FreeSpaceType>*>) {
      return fwd_distance_map(&VertexProp::fwd_distance_accum);
    } else {
      return fwd_distance_map(&VertexProp::heuristic_value);
    }
  }

  template <typename MotionGraphType, typename VertexProp>
  static auto dispatched_make_succ_map() {
    if constexpr (std::is_convertible_v<VertexProp*, bidir_optimal_mg_vertex<FreeSpaceType>*>) {
      return successor_map_t<MotionGraphType>(&VertexProp::successor);
    } else {
      return successor_map_t<MotionGraphType>();
    }
  }

  template <typename MotionGraphType, typename NcSelector>
  static ReaK::graph::sbastar_bundle<
      MotionGraphType, super_space_type, visitor_type, NcSelector,
      typename boost::property_map<MotionGraphType,
                                   double vertex_prop::*>::type,
      position_map, weight_map, density_map, constriction_map, distance_map,
      predecessor_map, fwd_distance_map,
      successor_map_t<MotionGraphType>>
  make_bundle(MotionGraphType& motion_graph, visitor_type& vis,
              NcSelector nc_selector,
              std::shared_ptr<const super_space_type> s_ptr) {
    using Vertex = graph::graph_vertex_t<MotionGraphType>;
    return ReaK::graph::make_sbastar_bundle(
        motion_graph, std::any_cast<Vertex>(vis.m_start_node),
        (!vis.m_goal_node.has_value()
             ? boost::graph_traits<MotionGraphType>::null_vertex()
             : std::any_cast<Vertex>(vis.m_goal_node)),
        *s_ptr, vis, nc_selector, get(&vertex_prop::key_value, motion_graph),
        position_map(&vertex_prop::position), weight_map(&edge_prop::weight),
        density_map(&vertex_prop::density),
        constriction_map(&vertex_prop::constriction),
        distance_map(&vertex_prop::distance_accum),
        predecessor_map(&vertex_prop::predecessor),
        dispatched_make_fwd_dist_map<vertex_prop>(),
        dispatched_make_succ_map<MotionGraphType, vertex_prop>());
  }

  template <typename MotionGraphType, typename NcSelector>
  static void make_call_to_planner_unidir(MotionGraphType& motion_graph,
                                  visitor_type& vis, NcSelector nc_selector,
                                  std::shared_ptr<const super_space_type> s_ptr,
                                  std::size_t method_flags, double init_relax,
                                  double init_temp) {
    using namespace graph;

    if (((method_flags & ADDITIONAL_PLANNING_BIAS_MASK) &
         PLAN_WITH_ANYTIME_HEURISTIC) &&
        (init_relax > 1e-6)) {
      if ((method_flags & COLLISION_CHECKING_POLICY_MASK) ==
          EAGER_COLLISION_CHECKING) {
        if ((method_flags & ADDITIONAL_PLANNING_BIAS_MASK) &
            PLAN_WITH_VORONOI_PULL) {
          generate_anytime_sbarrtstar(
              make_bundle(motion_graph, vis, nc_selector, s_ptr),
              get(random_sampler, *s_ptr), init_relax, init_temp);
        } else { /* assume nominal method only. */
          generate_anytime_sbastar(
              make_bundle(motion_graph, vis, nc_selector, s_ptr), init_relax);
        }
      } else { /* assume lazy collision checking */
        if ((method_flags & ADDITIONAL_PLANNING_BIAS_MASK) &
            PLAN_WITH_VORONOI_PULL) {
          if (method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) {
            generate_anytime_lazy_bnb_sbarrtstar(
                make_bundle(motion_graph, vis, nc_selector, s_ptr),
                get(random_sampler, *s_ptr), init_relax, init_temp);
          } else { /* assume nominal method only. */
            generate_anytime_lazy_sbarrtstar(
                make_bundle(motion_graph, vis, nc_selector, s_ptr),
                get(random_sampler, *s_ptr), init_relax, init_temp);
          }
        } else { /* assume nominal method only. */
          if (method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) {
            generate_anytime_lazy_bnb_sbastar(
                make_bundle(motion_graph, vis, nc_selector, s_ptr), init_relax);
          } else { /* assume nominal method only. */
            generate_anytime_lazy_sbastar(
                make_bundle(motion_graph, vis, nc_selector, s_ptr), init_relax);
          }
        }
      }
    } else {
      if ((method_flags & COLLISION_CHECKING_POLICY_MASK) ==
          EAGER_COLLISION_CHECKING) {
        if ((method_flags & ADDITIONAL_PLANNING_BIAS_MASK) &
            PLAN_WITH_VORONOI_PULL) {
          generate_sbarrtstar(
              make_bundle(motion_graph, vis, nc_selector, s_ptr),
              get(random_sampler, *s_ptr), init_temp);
        } else { /* assume nominal method only. */
          generate_sbastar(make_bundle(motion_graph, vis, nc_selector, s_ptr));
        }
      } else { /* assume lazy collision checking */
        if ((method_flags & ADDITIONAL_PLANNING_BIAS_MASK) &
            PLAN_WITH_VORONOI_PULL) {
          if (method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) {
            generate_lazy_bnb_sbarrtstar(
                make_bundle(motion_graph, vis, nc_selector, s_ptr),
                get(random_sampler, *s_ptr), init_temp);
          } else { /* assume nominal method only. */
            generate_lazy_sbarrtstar(
                make_bundle(motion_graph, vis, nc_selector, s_ptr),
                get(random_sampler, *s_ptr), init_temp);
          }
        } else { /* assume nominal method only. */
          if (method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) {
            generate_lazy_bnb_sbastar(
                make_bundle(motion_graph, vis, nc_selector, s_ptr));
          } else { /* assume nominal method only. */
            generate_lazy_sbastar(
                make_bundle(motion_graph, vis, nc_selector, s_ptr));
          }
        }
      }
    }
  }

  template <typename MotionGraphType, typename NcSelector>
  static void make_call_to_planner_bidir(MotionGraphType& motion_graph,
                                  visitor_type& vis, NcSelector nc_selector,
                                  std::shared_ptr<const super_space_type> s_ptr,
                                  std::size_t method_flags, double init_relax,
                                  double init_temp) {
    using namespace graph;

    if (((method_flags & ADDITIONAL_PLANNING_BIAS_MASK) &
         PLAN_WITH_ANYTIME_HEURISTIC) &&
        (init_relax > 1e-6)) {
      if ((method_flags & COLLISION_CHECKING_POLICY_MASK) ==
          EAGER_COLLISION_CHECKING) {
        if ((method_flags & ADDITIONAL_PLANNING_BIAS_MASK) &
            PLAN_WITH_VORONOI_PULL) {
          generate_anytime_sbarrtstar_bidir(
              make_bundle(motion_graph, vis, nc_selector, s_ptr),
              get(random_sampler, *s_ptr), init_relax, init_temp);
        } else { /* assume nominal method only. */
          generate_anytime_sbastar_bidir(
              make_bundle(motion_graph, vis, nc_selector, s_ptr), init_relax);
        }
      } else { /* assume lazy collision checking */
        if ((method_flags & ADDITIONAL_PLANNING_BIAS_MASK) &
            PLAN_WITH_VORONOI_PULL) {
          if (method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) {
            generate_anytime_lazy_bnb_sbarrtstar_bidir(
                make_bundle(motion_graph, vis, nc_selector, s_ptr),
                get(random_sampler, *s_ptr), init_relax, init_temp);
          } else { /* assume nominal method only. */
            generate_anytime_lazy_sbarrtstar_bidir(
                make_bundle(motion_graph, vis, nc_selector, s_ptr),
                get(random_sampler, *s_ptr), init_relax, init_temp);
          }
        } else { /* assume nominal method only. */
          if (method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) {
            generate_anytime_lazy_bnb_sbastar_bidir(
                make_bundle(motion_graph, vis, nc_selector, s_ptr), init_relax);
          } else { /* assume nominal method only. */
            generate_anytime_lazy_sbastar_bidir(
                make_bundle(motion_graph, vis, nc_selector, s_ptr), init_relax);
          }
        }
      }
    } else {
      if ((method_flags & COLLISION_CHECKING_POLICY_MASK) ==
          EAGER_COLLISION_CHECKING) {
        if ((method_flags & ADDITIONAL_PLANNING_BIAS_MASK) &
            PLAN_WITH_VORONOI_PULL) {
          generate_sbarrtstar_bidir(
              make_bundle(motion_graph, vis, nc_selector, s_ptr),
              get(random_sampler, *s_ptr), init_temp);
        } else { /* assume nominal method only. */
          generate_sbastar_bidir(
              make_bundle(motion_graph, vis, nc_selector, s_ptr));
        }
      } else { /* assume lazy collision checking */
        if ((method_flags & ADDITIONAL_PLANNING_BIAS_MASK) &
            PLAN_WITH_VORONOI_PULL) {
          if (method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) {
            generate_lazy_bnb_sbarrtstar_bidir(
                make_bundle(motion_graph, vis, nc_selector, s_ptr),
                get(random_sampler, *s_ptr), init_temp);
          } else { /* assume nominal method only. */
            generate_lazy_sbarrtstar_bidir(
                make_bundle(motion_graph, vis, nc_selector, s_ptr),
                get(random_sampler, *s_ptr), init_temp);
          }
        } else { /* assume nominal method only. */
          if (method_flags & USE_BRANCH_AND_BOUND_PRUNING_FLAG) {
            generate_lazy_bnb_sbastar_bidir(
                make_bundle(motion_graph, vis, nc_selector, s_ptr));
          } else { /* assume nominal method only. */
            generate_lazy_sbastar_bidir(
                make_bundle(motion_graph, vis, nc_selector, s_ptr));
          }
        }
      }
    }
  }

  template <typename MotionGraphType, typename NcSelector>
  static void make_call_to_planner(
      MotionGraphType& motion_graph, visitor_type& vis, NcSelector nc_selector,
      std::shared_ptr<const super_space_type> s_ptr, std::size_t method_flags,
      double init_relax, double init_temp) {
    if constexpr (is_bidir) {
      make_call_to_planner_bidir(motion_graph, vis, nc_selector,
                                 s_ptr, method_flags, init_relax, init_temp);
    } else {
      make_call_to_planner_unidir(motion_graph, vis, nc_selector,
                                  s_ptr, method_flags, init_relax, init_temp);
    }
  }
};

template <typename FreeSpaceType>
template <typename SBAStarFactory>
void sbastar_planner<FreeSpaceType>::solve_planning_query_impl(
    planning_query<FreeSpaceType>& aQuery) {

  this->reset_internal_state();

  double space_dim = double(this->get_space_dimensionality());
  double space_Lc = aQuery.get_heuristic_to_goal(aQuery.get_start_position());

  using SuperSpace = typename subspace_traits<FreeSpaceType>::super_space_type;
  std::shared_ptr<const SuperSpace> sup_space_ptr(
      &(this->m_space->get_super_space()), null_deleter());

  // Some MACROs to reduce the size of the code below.

#define RK_SBASTAR_PLANNER_SETUP_LS_OR_DVP_SUPPORT_STRUCTURES                  \
  using MotionGraphType = typename MGFactory::type;                            \
  MotionGraphType motion_graph = MGFactory::get_motion_graph();                \
  SBAStarFactory::init_motion_graph(motion_graph, vis, aQuery);                \
                                                                               \
  using SpacePartType = typename MGFactory::space_part_type;                   \
  using BasicVProp = typename SBAStarFactory::basic_vertex_prop;               \
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

#define RK_SBASTAR_PLANNER_SETUP_ALT_SUPPORT_STRUCTURES                        \
  using MotionGraphType = typename MGFactory::type;                            \
  using SpacePartType = typename MGFactory::space_part_type;                   \
  using PosMap = typename SBAStarFactory::position_map;                        \
  using VProp = typename SBAStarFactory::vertex_prop;                          \
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
  SBAStarFactory::init_motion_graph(motion_graph, vis, aQuery);

  using VisitorType = typename SBAStarFactory::visitor_type;
  VisitorType vis(this, &aQuery, NULL, std::any(), std::any(),
                  this->m_init_dens_threshold);

  if ((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) ==
      ADJ_LIST_MOTION_GRAPH) {

    if ((this->m_data_structure_flags & KNN_METHOD_MASK) == LINEAR_SEARCH_KNN) {

      using MGFactory = typename SBAStarFactory::ls_motion_graph;

      RK_SBASTAR_PLANNER_SETUP_LS_OR_DVP_SUPPORT_STRUCTURES

      SBAStarFactory::make_call_to_planner(
          motion_graph, vis, nc_selector, sup_space_ptr,
          this->m_planning_method_flags, this->m_init_relaxation,
          this->m_SA_init_temperature);

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_BF2_TREE_KNN) {

      using MGFactory = typename SBAStarFactory::template dvp_motion_graph<2, boost::bfl_d_ary_tree_storage<2>>;

      RK_SBASTAR_PLANNER_SETUP_LS_OR_DVP_SUPPORT_STRUCTURES

      SBAStarFactory::make_call_to_planner(
          motion_graph, vis, nc_selector, sup_space_ptr,
          this->m_planning_method_flags, this->m_init_relaxation,
          this->m_SA_init_temperature);

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_BF4_TREE_KNN) {

      using MGFactory = typename SBAStarFactory::template dvp_motion_graph<4, boost::bfl_d_ary_tree_storage<4>>;

      RK_SBASTAR_PLANNER_SETUP_LS_OR_DVP_SUPPORT_STRUCTURES

      SBAStarFactory::make_call_to_planner(
          motion_graph, vis, nc_selector, sup_space_ptr,
          this->m_planning_method_flags, this->m_init_relaxation,
          this->m_SA_init_temperature);

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_COB2_TREE_KNN) {

      using MGFactory = typename SBAStarFactory::template dvp_motion_graph<2, boost::vebl_d_ary_tree_storage<2>>;

      RK_SBASTAR_PLANNER_SETUP_LS_OR_DVP_SUPPORT_STRUCTURES

      SBAStarFactory::make_call_to_planner(
          motion_graph, vis, nc_selector, sup_space_ptr,
          this->m_planning_method_flags, this->m_init_relaxation,
          this->m_SA_init_temperature);

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_COB4_TREE_KNN) {

      using MGFactory = typename SBAStarFactory::template dvp_motion_graph<4, boost::vebl_d_ary_tree_storage<4>>;

      RK_SBASTAR_PLANNER_SETUP_LS_OR_DVP_SUPPORT_STRUCTURES

      SBAStarFactory::make_call_to_planner(
          motion_graph, vis, nc_selector, sup_space_ptr,
          this->m_planning_method_flags, this->m_init_relaxation,
          this->m_SA_init_temperature);

#endif
    }

#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT

  } else if ((this->m_data_structure_flags & MOTION_GRAPH_STORAGE_MASK) ==
             DVP_ADJ_LIST_MOTION_GRAPH) {

    if ((this->m_data_structure_flags & KNN_METHOD_MASK) == DVP_BF2_TREE_KNN) {

      using MGFactory = typename SBAStarFactory::template alt_motion_graph<2, boost::bfl_d_ary_tree_storage<2>>;

      RK_SBASTAR_PLANNER_SETUP_ALT_SUPPORT_STRUCTURES

      SBAStarFactory::make_call_to_planner(
          motion_graph, vis, nc_selector, sup_space_ptr,
          this->m_planning_method_flags, this->m_init_relaxation,
          this->m_SA_init_temperature);

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_BF4_TREE_KNN) {

      using MGFactory = typename SBAStarFactory::template alt_motion_graph<4, boost::bfl_d_ary_tree_storage<4>>;

      RK_SBASTAR_PLANNER_SETUP_ALT_SUPPORT_STRUCTURES

      SBAStarFactory::make_call_to_planner(
          motion_graph, vis, nc_selector, sup_space_ptr,
          this->m_planning_method_flags, this->m_init_relaxation,
          this->m_SA_init_temperature);

#ifdef RK_PLANNERS_ENABLE_VEBL_TREE

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_COB2_TREE_KNN) {

      using MGFactory = typename SBAStarFactory::template alt_motion_graph<2, boost::vebl_d_ary_tree_storage<2>>;

      RK_SBASTAR_PLANNER_SETUP_ALT_SUPPORT_STRUCTURES

      SBAStarFactory::make_call_to_planner(
          motion_graph, vis, nc_selector, sup_space_ptr,
          this->m_planning_method_flags, this->m_init_relaxation,
          this->m_SA_init_temperature);

    } else if ((this->m_data_structure_flags & KNN_METHOD_MASK) ==
               DVP_COB4_TREE_KNN) {

      using MGFactory = typename SBAStarFactory::template alt_motion_graph<4, boost::vebl_d_ary_tree_storage<4>>;

      RK_SBASTAR_PLANNER_SETUP_ALT_SUPPORT_STRUCTURES

      SBAStarFactory::make_call_to_planner(
          motion_graph, vis, nc_selector, sup_space_ptr,
          this->m_planning_method_flags, this->m_init_relaxation,
          this->m_SA_init_temperature);

#endif
    }

#endif
  }

#undef RK_SBASTAR_PLANNER_SETUP_LS_OR_DVP_SUPPORT_STRUCTURES
#undef RK_SBASTAR_PLANNER_SETUP_ALT_SUPPORT_STRUCTURES
}

template <typename FreeSpaceType>
void sbastar_planner<FreeSpaceType>::solve_planning_query(
    planning_query<FreeSpaceType>& aQuery) {
  if ((this->m_planning_method_flags & PLANNING_DIRECTIONALITY_MASK) ==
      UNIDIRECTIONAL_PLANNING) {
    this->solve_planning_query_impl<
        sbastar_bundle_factory<FreeSpaceType, false>>(aQuery);
  } else { /* NOTE: Bi-directional version: */
    this->solve_planning_query_impl<
        sbastar_bundle_factory<FreeSpaceType, true>>(aQuery);
  }
}

}  // namespace ReaK::pp

#endif

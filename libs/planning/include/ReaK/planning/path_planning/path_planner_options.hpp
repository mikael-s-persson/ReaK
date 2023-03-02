/**
 * \file path_planner_options.hpp
 *
 * This library defines the options available when creating a path-planner.
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

#ifndef REAK_PATH_PLANNER_OPTIONS_HPP
#define REAK_PATH_PLANNER_OPTIONS_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/shared_object.hpp>

namespace ReaK::pp {

/// This mask indicates the storage strategy used for the motion-graph.
static constexpr std::size_t MOTION_GRAPH_STORAGE_MASK = 0x03;

/// This flag indicates that the motion-graph should be stored as an adjacency-list graph.
static constexpr std::size_t ADJ_LIST_MOTION_GRAPH = 0;
/// This flag indicates that the motion-graph should be stored as an adjacency-list graph that is overlaid on a dynamic
/// vantage-point tree.
static constexpr std::size_t DVP_ADJ_LIST_MOTION_GRAPH = 1;
/// This flag indicates that the motion-graph should be stored as a linked-tree (a tree with links, i.e., like a
/// linked-list).
static constexpr std::size_t LINKED_TREE_MOTION_GRAPH = 2;

/// This mask indicates the nearest-neighbor method used during the motion-planning.
static constexpr std::size_t KNN_METHOD_MASK = 0x0F << 2;

/// This flag indicates that the nearest-neighbor queries should be done via a linear search method.
static constexpr std::size_t LINEAR_SEARCH_KNN = 0;
/// This flag indicates that the nearest-neighbor queries should be done via an approximate linear search method (i.e.,
/// some partial and randomized method that looks at only a fraction of the points for an approximate KNN set).
static constexpr std::size_t APPROX_LINEAR_SEARCH_KNN = 1 << 2;
/// This flag indicates that the nearest-neighbor queries should be done via a best-first search through a tree method
/// (this method will be approximate in most cases).
static constexpr std::size_t DVP_LINKED_TREE_KNN = 2 << 2;
/// This flag indicates that the nearest-neighbor queries should be done via a DVP-tree laid out on a contiguous storage
/// in breadth-first layout of arity 2 (binary).
static constexpr std::size_t DVP_BF2_TREE_KNN = 3 << 2;
/// This flag indicates that the nearest-neighbor queries should be done via a DVP-tree laid out on a contiguous storage
/// in breadth-first layout of arity 4.
static constexpr std::size_t DVP_BF4_TREE_KNN = 4 << 2;
/// This flag indicates that the nearest-neighbor queries should be done via a DVP-tree laid out on a contiguous storage
/// in cache-oblivious breadth-first layout of arity 2 (binary).
static constexpr std::size_t DVP_COB2_TREE_KNN = 5 << 2;
/// This flag indicates that the nearest-neighbor queries should be done via a DVP-tree laid out on a contiguous storage
/// in cache-oblivious breadth-first layout of arity 4.
static constexpr std::size_t DVP_COB4_TREE_KNN = 6 << 2;

/// This mask indicates the collision checking policy that should be preferred (generally should be lazy, unless the
/// bird-flight distance metric is not good as a reflection of the collision-free travel cost).
static constexpr std::size_t COLLISION_CHECKING_POLICY_MASK = 0x03;

/// This flag indicates that an eager collision checking policy should be preferred (generally only useful if the
/// bird-flight distance metric is not good as a reflection of the collision-free travel cost).
static constexpr std::size_t EAGER_COLLISION_CHECKING = 0;
/// This flag indicates that a lazy collision checking policy should be preferred, meaning that the bird-flight distance
/// will be used tentatively for edge-costs until the edges are considered for inclusion in the optimal path, at which
/// point collision is checked.
static constexpr std::size_t LAZY_COLLISION_CHECKING = 1;

/// This mask indicates the kinds of additional biases should be applied during the motion planning.
static constexpr std::size_t ADDITIONAL_PLANNING_BIAS_MASK = 0x3F << 2;

/// This flag indicates that the plain version (no bias) of an algorithm is to be used.
static constexpr std::size_t NOMINAL_PLANNER_ONLY = 0;
/// This flag indicates that a Voronoi pull (i.e., RRT-style expansion) should be used if the main planning algorithm
/// (usually a more greedy algorithm) gets "stuck". Note, this can be combined with other biasing flags.
static constexpr std::size_t PLAN_WITH_VORONOI_PULL = 0x01 << 2;
/// This flag indicates that a narrow-passage push (i.e., PRM-style expansion) should be used if the main planning
/// algorithm (usually a more greedy algorithm) gets "stuck". Note, this can be combined with other biasing flags.
static constexpr std::size_t PLAN_WITH_NARROW_PASSAGE_PUSH = 0x02 << 2;
/// This flag indicates that an anytime heuristic should be used. This is for optimizing planning algorithms and causes
/// them to more eagerly seek a connection (complete path) before optimizing the solution. Note, this can be combined
/// with other biasing flags.
static constexpr std::size_t PLAN_WITH_ANYTIME_HEURISTIC = 0x04 << 2;

/// This mask indicates the kinds of additional biases should be applied during the motion planning.
static constexpr std::size_t PLANNING_DIRECTIONALITY_MASK = 0x3 << 8;

/// This flag indicates that the motion-graph should be grown uni-directionally (from start).
static constexpr std::size_t UNIDIRECTIONAL_PLANNING = 0;
/// This flag indicates that the motion-graph should be grown bi-directionally (from start AND goal).
static constexpr std::size_t BIDIRECTIONAL_PLANNING = 1 << 8;

/// This flag indicates a branch-and-bound pruning strategy should be use during the motion planning.
static constexpr std::size_t USE_BRANCH_AND_BOUND_PRUNING_FLAG = 0x01 << 10;

class planning_option_collection : public shared_object {
 public:
  std::size_t planning_algo;
  std::size_t max_vertices;
  std::size_t prog_interval;
  std::size_t max_results;
  std::size_t planning_options;
  std::size_t store_policy;
  std::size_t knn_method;
  double init_SA_temp;
  double init_relax;
  double max_random_walk;
  double start_delay;

  planning_option_collection()
      : planning_algo(0),
        max_vertices(5000),
        prog_interval(10),
        max_results(50),
        planning_options(0),
        store_policy(0),
        knn_method(0),
        init_SA_temp(0.0),
        init_relax(0.0),
        max_random_walk(1.0),
        start_delay(20.0) {}

  std::string get_planning_algo_str() const {
    switch (planning_algo) {
      case 0:
        return "rrt";
      case 1:
        return "rrt_star";
      case 2:
        return "prm";
      case 3:
        return "sba_star";
      case 4:
        return "fadprm";
      default:
        return "rrt";
    }
  }

  std::string get_knn_method_str() const {
    switch (knn_method & KNN_METHOD_MASK) {
      case LINEAR_SEARCH_KNN:
        return "linear";
      case DVP_BF2_TREE_KNN:
        return "bf2";
      case DVP_BF4_TREE_KNN:
        return "bf4";
      case DVP_COB2_TREE_KNN:
        return "cob2";
      case DVP_COB4_TREE_KNN:
        return "cob4";
      default:
        return "";
    }
  }

  std::string get_mg_storage_str() const {
    switch (store_policy & MOTION_GRAPH_STORAGE_MASK) {
      case ADJ_LIST_MOTION_GRAPH:
        return "adj_list";
      case DVP_ADJ_LIST_MOTION_GRAPH:
        return "dvp_adj_list";
      case LINKED_TREE_MOTION_GRAPH:
        return "linked_tree";
      default:
        return "";
    }
  }

  std::string get_planner_qualifier_str() const {
    std::string result = "";
    if (planning_options & BIDIRECTIONAL_PLANNING) {
      result += "_bidir";
    }
    if (planning_options & USE_BRANCH_AND_BOUND_PRUNING_FLAG) {
      result += "_bnb";
    }
    if (planning_options & PLAN_WITH_ANYTIME_HEURISTIC) {
      result += "_any";
    }
    if (planning_options & PLAN_WITH_VORONOI_PULL) {
      result += "_sa";
    }
    return result;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A, unsigned int) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(planning_algo) &
        RK_SERIAL_SAVE_WITH_NAME(max_vertices) &
        RK_SERIAL_SAVE_WITH_NAME(prog_interval) &
        RK_SERIAL_SAVE_WITH_NAME(max_results) &
        RK_SERIAL_SAVE_WITH_NAME(planning_options) &
        RK_SERIAL_SAVE_WITH_NAME(store_policy) &
        RK_SERIAL_SAVE_WITH_NAME(knn_method) &
        RK_SERIAL_SAVE_WITH_NAME(init_SA_temp) &
        RK_SERIAL_SAVE_WITH_NAME(init_relax) &
        RK_SERIAL_SAVE_WITH_NAME(max_random_walk) &
        RK_SERIAL_SAVE_WITH_NAME(start_delay);
  }

  void load(serialization::iarchive& A, unsigned int) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(planning_algo) &
        RK_SERIAL_LOAD_WITH_NAME(max_vertices) &
        RK_SERIAL_LOAD_WITH_NAME(prog_interval) &
        RK_SERIAL_LOAD_WITH_NAME(max_results) &
        RK_SERIAL_LOAD_WITH_NAME(planning_options) &
        RK_SERIAL_LOAD_WITH_NAME(store_policy) &
        RK_SERIAL_LOAD_WITH_NAME(knn_method) &
        RK_SERIAL_LOAD_WITH_NAME(init_SA_temp) &
        RK_SERIAL_LOAD_WITH_NAME(init_relax) &
        RK_SERIAL_LOAD_WITH_NAME(max_random_walk) &
        RK_SERIAL_LOAD_WITH_NAME(start_delay);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(planning_option_collection, 0xC2460019, 1,
                              "planning_option_collection", shared_object)
};

}  // namespace ReaK::pp

#endif

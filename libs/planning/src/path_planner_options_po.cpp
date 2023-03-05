
/*
 *    Copyright 2023 Sven Mikael Persson
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
#include "ReaK/planning/path_planning/path_planner_options_po.hpp"
#include "ReaK/planning/path_planning/path_planner_options.hpp"

#include "ReaK/core/serialization/archiver_factory.hpp"

#include "absl/flags/flag.h"

// Planning algorithm options
ABSL_FLAG(std::string, pp_planner_options, "",
          "Specify the file containing the planner-options data.");
ABSL_FLAG(std::string, pp_planner_alg, "",
          "Specify the planner algorithm to use, can be any of (rrt, rrt_star, "
          "prm, sba_star, fadprm).");
ABSL_FLAG(int, pp_max_vertices, 5000,
          "Maximum number of vertices during runs (default is 5000).");
ABSL_FLAG(int, pp_max_results, 50,
          "Maximum number of result-paths during runs (default is 50).");
ABSL_FLAG(
    int, pp_prog_interval, 10,
    "Number of vertices between progress reports during runs (default is 10).");
ABSL_FLAG(double, pp_max_random_walk, 1.0,
          "Specify the maximum random-walk distance allowed in the algorithm "
          "(default: 1.0). Only meaningful for expanding algorithms (SBA*, "
          "FADPRM, PRM).");
ABSL_FLAG(bool, pp_no_lazy_connect, false,
          "If set, disable lazy connection strategy during planning.");
ABSL_FLAG(bool, pp_bi_directional, false,
          "Specify whether to use a bi-directional algorithm or not during "
          "planning. Only supported for some algorithms (RRT, RRT*, SBA*).");
ABSL_FLAG(bool, pp_with_bnb, false,
          "Specify whether to use a Branch-and-bound or not during planning to "
          "prune useless nodes from the motion-graph. Only supported for "
          "optimizing algorithms (RRT*, SBA*).");
ABSL_FLAG(double, pp_relaxation_factor, 0.0,
          "Specify the initial relaxation factor for the algorithm (default: "
          "0.0). Only supported for heuristic-driven algorithms.");
ABSL_FLAG(double, pp_start_delay, 0.0,
          "Specify the starting time delay for the algorithm (default: 0.0). "
          "Only for dynamic problems (temporal space).");
ABSL_FLAG(bool, pp_with_voronoi_pull, false,
          "Specify whether to use a Voronoi pull or not to add an exploratory "
          "bias to the search (default: not).");
ABSL_FLAG(double, pp_sa_temperature, -1.0,
          "Specify the initial Simulated Annealing temperature for algorithms "
          "that work on a exploration-exploitation schedule (e.g., SA-SBA*).");
#ifdef RK_PLANNERS_ENABLE_VEBL_TREE
ABSL_FLAG(std::string, pp_knn_method, "",
          "Specify the KNN method to use (supported options: linear, bf2, bf4, "
          "cob2, cob4) (default: bf2).");
#else
ABSL_FLAG(std::string, pp_knn_method, "",
          "Specify the KNN method to use (supported options: linear, bf2, bf4) "
          "(default: bf2).");
#endif
#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT
ABSL_FLAG(std::string, pp_mg_storage, "",
          "Specify the KNN method to use (supported options: adj-list, "
          "dvp-adj-list) (default: adj-list).");
#else
ABSL_FLAG(std::string, pp_mg_storage, "",
          "Specify the KNN method to use (supported options: adj-list) "
          "(default: adj-list).");
#endif

namespace ReaK::pp {

planning_option_collection get_planning_option_from_flags() {
  planning_option_collection plan_options;

  if (!absl::GetFlag(FLAGS_pp_planner_options).empty()) {
    try {
      (*serialization::open_iarchive(
          absl::GetFlag(FLAGS_pp_planner_options))) >>
          plan_options;
    } catch (std::exception& e) {
      RK_UNUSED(e);
    }
  }

  if (!absl::GetFlag(FLAGS_pp_planner_alg).empty()) {
    plan_options.planning_algo = 0;  // RRT (default)
    if (absl::GetFlag(FLAGS_pp_planner_alg) == "rrt_star") {
      plan_options.planning_algo = 1;
    } else if (absl::GetFlag(FLAGS_pp_planner_alg) == "prm") {
      plan_options.planning_algo = 2;
    } else if (absl::GetFlag(FLAGS_pp_planner_alg) == "sba_star") {
      plan_options.planning_algo = 3;
    } else if (absl::GetFlag(FLAGS_pp_planner_alg) == "fadprm") {
      plan_options.planning_algo = 4;
    }
  }

  if (absl::GetFlag(FLAGS_pp_max_vertices) != 5000) {
    plan_options.max_vertices = absl::GetFlag(FLAGS_pp_max_vertices);
  }
  if (absl::GetFlag(FLAGS_pp_max_results) != 50) {
    plan_options.max_results = absl::GetFlag(FLAGS_pp_max_results);
  }
  if (absl::GetFlag(FLAGS_pp_prog_interval) != 10) {
    plan_options.prog_interval = absl::GetFlag(FLAGS_pp_prog_interval);
  }

  if (!absl::GetFlag(FLAGS_pp_knn_method).empty()) {
    plan_options.knn_method = 0;
    if ((absl::GetFlag(FLAGS_pp_knn_method) == "linear") &&
        (absl::GetFlag(FLAGS_pp_mg_storage) == "adj-list")) {
      plan_options.knn_method |= LINEAR_SEARCH_KNN;
    } else if (absl::GetFlag(FLAGS_pp_knn_method) == "bf4") {
      plan_options.knn_method |= DVP_BF4_TREE_KNN;
    }
#ifdef RK_PLANNERS_ENABLE_VEBL_TREE
    else if (absl::GetFlag(FLAGS_pp_knn_method) == "cob2") {
      plan_options.knn_method |= DVP_COB2_TREE_KNN;
    } else if (absl::GetFlag(FLAGS_pp_knn_method) == "cob4") {
      plan_options.knn_method |= DVP_COB4_TREE_KNN;
    }
#endif
    else {
      plan_options.knn_method |= DVP_BF2_TREE_KNN;
    }
  }

  if (!absl::GetFlag(FLAGS_pp_mg_storage).empty()) {
    plan_options.store_policy = 0;
#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT
    if (absl::GetFlag(FLAGS_pp_mg_storage) == "dvp-adj-list") {
      plan_options.store_policy |= DVP_ADJ_LIST_MOTION_GRAPH;
    } else
#endif
    {
      plan_options.store_policy |= ADJ_LIST_MOTION_GRAPH;
    }
  }

  if (!absl::GetFlag(FLAGS_pp_no_lazy_connect)) {
    plan_options.planning_options |=
        LAZY_COLLISION_CHECKING;  // use lazy, if supported.
  }

  if (absl::GetFlag(FLAGS_pp_bi_directional)) {
    plan_options.planning_options |= BIDIRECTIONAL_PLANNING;
  }

  if (absl::GetFlag(FLAGS_pp_with_bnb)) {
    plan_options.planning_options |= USE_BRANCH_AND_BOUND_PRUNING_FLAG;
  }

  if (absl::GetFlag(FLAGS_pp_relaxation_factor) > 1e-6) {
    plan_options.init_relax = absl::GetFlag(FLAGS_pp_relaxation_factor);
    plan_options.planning_options |= PLAN_WITH_ANYTIME_HEURISTIC;
  }

  if (absl::GetFlag(FLAGS_pp_sa_temperature) > -1.0) {
    plan_options.init_SA_temp = absl::GetFlag(FLAGS_pp_sa_temperature);
  }

  if (absl::GetFlag(FLAGS_pp_max_random_walk) != 1.0) {
    plan_options.max_random_walk = absl::GetFlag(FLAGS_pp_max_random_walk);
  }

  if (absl::GetFlag(FLAGS_pp_with_voronoi_pull)) {
    plan_options.planning_options |= PLAN_WITH_VORONOI_PULL;
  }

  return plan_options;
}

}  // namespace ReaK::pp

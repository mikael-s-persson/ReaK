
/*
 *    Copyright 2011 Sven Mikael Persson
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

#include <fstream>
#include <iostream>

#include "ReaK/topologies/spaces/ptrobot2D_test_world.h"

// #define RK_DISABLE_RRT_PLANNER
// #define RK_DISABLE_RRTSTAR_PLANNER
// #define RK_DISABLE_PRM_PLANNER
// #define RK_DISABLE_FADPRM_PLANNER
// #define RK_DISABLE_SBASTAR_PLANNER

#include "ReaK/planning/path_planning/planner_exec_engines.h"

#include "ReaK/planning/path_planning/path_planner_options_po.h"

#include "ReaK/math/optimization/optim_exceptions.h"

#include "ReaK/core/serialization/archiver_factory.h"

#include <filesystem>
#include <memory>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

// I/O options
ABSL_FLAG(std::string, input_file, "",
          "Specify the input file (default is stdin).");
ABSL_FLAG(std::string, output_path, "pp_results",
          "Specify the output path (default is pp_results).");

// Monte-Carlo options
ABSL_FLAG(
    bool, monte_carlo, false,
    "Specify that monte-carlo runs should be performed (default is not).");
ABSL_FLAG(int, mc_runs, 100,
          "Number of monte-carlo runs to average out (default is 100).");

// Single-run options
ABSL_FLAG(bool, single_run, false,
          "Specify that single runs should be performed (default is not).");

// More planning options
ABSL_FLAG(
    double, max_edge_length, 20.0,
    "Maximum length (in pixels) of edges of the motion-graph (default is 20).");

// File generation options
ABSL_FLAG(std::string, generate_all_files, "",
          "Specify that all configuration files should be generated with the "
          "given file-name prefix (file-name without suffix and extension).");
ABSL_FLAG(std::string, generate_planner_options, "",
          "Specify that the planner options file should be generated with the "
          "given file-name prefix (file-name without extension).");
ABSL_FLAG(bool, generate_xml, false,
          "If set, output results in XML format (rkx) (default).");
ABSL_FLAG(bool, generate_protobuf, false,
          "If set, output results in protobuf format (pb).");
ABSL_FLAG(bool, generate_binary, false,
          "If set, output results in binary format (rkb).");

namespace fs = std::filesystem;

int main(int argc, char** argv) {

  using namespace ReaK;
  using namespace pp;

  absl::ParseCommandLine(argc, argv);

  if (static_cast<int>(!absl::GetFlag(FLAGS_input_file).empty()) *
              (static_cast<int>(absl::GetFlag(FLAGS_monte_carlo)) +
               static_cast<int>(absl::GetFlag(FLAGS_single_run))) +
          static_cast<int>(!absl::GetFlag(FLAGS_generate_all_files).empty()) +
          static_cast<int>(
              absl::GetFlag(FLAGS_generate_planner_options).empty()) <
      1) {
    std::cout << "Error: There was no action specified! This program is "
                 "designed to perform Monte-Carlo runs, single "
                 "runs (with output), or generate the configuration files to "
                 "construct scenarios. You must specify at "
                 "least one of these actions to be performed!"
              << std::endl;
    return 1;
  }

  std::string output_path_name = absl::GetFlag(FLAGS_output_path);
  while (output_path_name.back() == '/') {
    output_path_name.erase(output_path_name.length() - 1, 1);
  }

  fs::create_directory(output_path_name.c_str());

  planning_option_collection plan_options = get_planning_option_from_flags();

  std::string knn_method_str = plan_options.get_knn_method_str();
  std::string mg_storage_str = plan_options.get_mg_storage_str();
  std::string planner_qualifier_str = plan_options.get_planner_qualifier_str();
  std::string planner_name_str = plan_options.get_planning_algo_str() + "_" +
                                 planner_qualifier_str + "_" + mg_storage_str +
                                 "_" + knn_method_str;

  double max_radius = absl::GetFlag(FLAGS_max_edge_length);
  plan_options.max_random_walk = max_radius;

  // Do the generations if required:

  if (!absl::GetFlag(FLAGS_generate_all_files).empty() ||
      !absl::GetFlag(FLAGS_generate_planner_options).empty()) {
    std::string file_name;
    if (absl::GetFlag(FLAGS_generate_planner_options).empty()) {
      file_name = absl::GetFlag(FLAGS_generate_all_files) + "_planner";
    } else {
      file_name = absl::GetFlag(FLAGS_generate_planner_options);
    }
    if (absl::GetFlag(FLAGS_generate_protobuf)) {
      file_name += ".pb";
    } else if (absl::GetFlag(FLAGS_generate_binary)) {
      file_name += ".rkb";
    } else {
      file_name += ".rkx";
    }

    try {
      (*serialization::open_oarchive(file_name)) << plan_options;
    } catch ([[maybe_unused]] std::exception& e) {
      std::cerr << "Error: Could not generate the planner options file!"
                << std::endl;
    }
    if (absl::GetFlag(FLAGS_input_file).empty()) {
      return 0;
    }
  }

  std::string world_file_name = absl::GetFlag(FLAGS_input_file);
  std::string world_file_name_only(
      std::find(world_file_name.rbegin(), world_file_name.rend(), '/').base(),
      std::find(world_file_name.rbegin(), world_file_name.rend(), '.').base() -
          1);

  std::shared_ptr<ptrobot2D_test_world> world_map =
      std::make_shared<ptrobot2D_test_world>(world_file_name, max_radius, 1.0);

  if (absl::GetFlag(FLAGS_monte_carlo)) {
    monte_carlo_mp_engine mc_eng(absl::GetFlag(FLAGS_mc_runs), planner_name_str,
                                 output_path_name + "/" + world_file_name_only);
    try {
      execute_p2p_planner(world_map, plan_options, 2, mc_eng,
                          world_map->get_start_pos(),
                          world_map->get_goal_pos());
    } catch (std::exception& e) {
      std::cerr
          << "Error: An exception was raised during the planning:\nwhat(): "
          << e.what() << std::endl;
      return 2;
    }
  }

  if (absl::GetFlag(FLAGS_single_run)) {
    differ_report_mp_engine sr_eng(
        0.25 * plan_options.max_random_walk, planner_name_str,
        output_path_name + "/" + world_file_name_only);
    try {
      execute_p2p_planner(world_map, plan_options, 2, sr_eng,
                          world_map->get_start_pos(),
                          world_map->get_goal_pos());
    } catch (std::exception& e) {
      std::cerr
          << "Error: An exception was raised during the planning:\nwhat(): "
          << e.what() << std::endl;
      return 3;
    }
  }

  return 0;
}

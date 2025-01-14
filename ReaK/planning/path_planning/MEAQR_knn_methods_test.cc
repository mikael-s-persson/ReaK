
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

#include "bagl/adjacency_list.h"
#include "bagl/properties.h"
#include "bagl/topology.h"

#include "ReaK/math/lin_alg/vect_alg.h"
#include "ReaK/planning/path_planning/metric_space_search.h"
#include "ReaK/planning/path_planning/topological_search.h"

#include "ReaK/control/controllers/IHAQR_topology.h"
#include "ReaK/control/controllers/MEAQR_topology.h"
#include "ReaK/control/systems/quadrotor_system.h"
#include "ReaK/topologies/spaces/se3_random_samplers.h"

#include "ReaK/mbd/kte/kte_map_chain.h"

#include "ReaK/core/serialization/xml_archiver.h"
#include "ReaK/math/optimization/optim_exceptions.h"

#include <chrono>
#include <filesystem>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

ABSL_FLAG(std::string, space_def_file, "models/quadrotor_spaces.xml",
          "Specify the space-definition file (default is "
          "models/quadrotor_spaces.xml).");
ABSL_FLAG(std::string, output_path, "vp_results",
          "Specify the output path (default is vp_results).");

namespace fs = std::filesystem;

using namespace ReaK;
using namespace pp;
using namespace ctrl;

using X8_IHAQR_space_type =
    IHAQR_topology<quadrotor_system::state_space_type, quadrotor_system,
                   position_only_sampler>;
using X8_MEAQR_space_type =
    MEAQR_topology<quadrotor_system::state_space_type, quadrotor_system,
                   position_only_sampler>;

using MEAQR_PointType = X8_MEAQR_space_type::point_type;

struct MEAQR_vprop {
  MEAQR_PointType position;
};

int main(int argc, char** argv) {

  using namespace std::chrono;

  absl::ParseCommandLine(argc, argv);

  std::string output_path_name = absl::GetFlag(FLAGS_output_path);
  while (output_path_name[output_path_name.length() - 1] == '/') {
    output_path_name.erase(output_path_name.length() - 1, 1);
  }

  fs::create_directory(output_path_name.c_str());

  std::string space_def_file_name = absl::GetFlag(FLAGS_space_def_file);
  std::string space_def_file_name_only(
      std::find(space_def_file_name.rbegin(), space_def_file_name.rend(), '/')
          .base(),
      std::find(space_def_file_name.rbegin(), space_def_file_name.rend(), '.')
              .base() -
          1);
  if (!fs::exists(space_def_file_name)) {
    std::cout << "Error: input file '" << space_def_file_name
              << "' does not exist!" << std::endl
              << "Correct usage of this program is:" << std::endl;
    return 1;
  }

  std::shared_ptr<quadrotor_system> X8_sys;
  std::shared_ptr<X8_IHAQR_space_type> X8_IHAQR_space;
  std::shared_ptr<X8_MEAQR_space_type> X8_MEAQR_space;
  {
    serialization::xml_iarchive in(space_def_file_name);
    in >> X8_sys >> X8_IHAQR_space >> X8_MEAQR_space;
  }

  using WorldGridType =
      bagl::adjacency_list<bagl::vec_s, bagl::vec_s, bagl::undirected_s,
                           MEAQR_vprop, bagl::no_property>;
  using VertexType = bagl::graph_traits<WorldGridType>::vertex_descriptor;
  using PositionMapType =
      bagl::property_map<WorldGridType, MEAQR_PointType MEAQR_vprop::*>::type;

  using WorldPartition2 =
      ReaK::pp::dvp_tree<VertexType, X8_MEAQR_space_type, PositionMapType, 2>;
  using WorldPartition4 =
      ReaK::pp::dvp_tree<VertexType, X8_MEAQR_space_type, PositionMapType, 4>;

  const std::array<unsigned int, 16> grid_sizes = {
      100,  200,  300,  400,  500,  600,  800,  1000,
      1500, 2000, 2500, 3000, 4000, 5000, 7500, 10000};

  std::string output_file_name = output_path_name + "/X8_knn_times.dat";
  std::ofstream outFile(output_file_name.c_str());
  outFile << "N\tVP2\tVP4\tLS\t (all times in micro-seconds per query)"
          << std::endl;

  for (unsigned int sz : grid_sizes) {
    WorldGridType grid;

    PositionMapType m_position(get(&MEAQR_vprop::position, grid));

    std::cout << "Generating " << sz << " nodes..." << std::endl;
    for (unsigned int j = 0; j < sz; ++j) {
      std::cout << "\r" << std::setw(10) << j << std::flush;
      VertexType v = add_vertex(grid);
      put(m_position, v, X8_MEAQR_space->random_point());
    }
    std::cout << std::endl << "Done!" << std::endl;

    outFile << sz;
    std::cout << "N = " << sz << std::endl;

    {
      std::cout << "VP2 ..." << std::endl;
      WorldPartition2 part2_fresh(grid, X8_MEAQR_space, m_position);

      ReaK::pp::multi_dvp_tree_search<WorldGridType, WorldPartition2>
          nn_finder2_fresh;
      nn_finder2_fresh.graph_tree_map[&grid] = &part2_fresh;
      high_resolution_clock::time_point t_start = high_resolution_clock::now();
      for (unsigned int j = 0; j < 1000; ++j) {
        nn_finder2_fresh(X8_MEAQR_space->random_point(), grid, *X8_MEAQR_space,
                         m_position);
      };
      high_resolution_clock::duration dt =
          high_resolution_clock::now() - t_start;
      outFile << "\t" << duration_cast<microseconds>(dt).count() * 0.001;
      std::cout << "Done!" << std::endl;
    }

    {
      std::cout << "VP4 ..." << std::endl;
      WorldPartition4 part4_fresh(grid, X8_MEAQR_space, m_position);

      ReaK::pp::multi_dvp_tree_search<WorldGridType, WorldPartition4>
          nn_finder4_fresh;
      nn_finder4_fresh.graph_tree_map[&grid] = &part4_fresh;
      high_resolution_clock::time_point t_start = high_resolution_clock::now();
      for (unsigned int j = 0; j < 1000; ++j) {
        nn_finder4_fresh(X8_MEAQR_space->random_point(), grid, *X8_MEAQR_space,
                         m_position);
      };
      high_resolution_clock::duration dt =
          high_resolution_clock::now() - t_start;
      outFile << "\t" << duration_cast<microseconds>(dt).count() * 0.001;
      std::cout << "Done!" << std::endl;
    }

    {
      std::cout << "Linear Search ..." << std::endl;
      ReaK::pp::linear_neighbor_search<WorldGridType> lnn_finder;
      high_resolution_clock::time_point t_start = high_resolution_clock::now();
      for (unsigned int j = 0; j < 1000; ++j) {
        lnn_finder(X8_MEAQR_space->random_point(), grid, *X8_MEAQR_space,
                   m_position);
      };
      high_resolution_clock::duration dt =
          high_resolution_clock::now() - t_start;
      outFile << "\t" << duration_cast<microseconds>(dt).count() * 0.001;
      std::cout << "Done!" << std::endl;
    }

    outFile << std::endl;
  }
}

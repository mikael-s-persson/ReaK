
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
#include "ReaK/topologies/spaces/hyperbox_topology.h"

#include <chrono>

int main() {

  using namespace std::chrono;

  using TopologyType = ReaK::pp::hyperbox_topology<ReaK::vect<double, 6>>;

  using PointType = TopologyType::point_type;

  using WorldGridVertexProperties =
      bagl::property<bagl::vertex_position_t, PointType, bagl::no_property>;

  using WorldGridEdgeProperties = bagl::no_property;

  using WorldGridType =
      bagl::adjacency_list<bagl::vec_s, bagl::vec_s, bagl::undirected_s,
                           WorldGridVertexProperties, WorldGridEdgeProperties>;

  using VertexType = bagl::graph_traits<WorldGridType>::vertex_descriptor;
  using WorldPartition4 = ReaK::pp::dvp_tree<
      VertexType, TopologyType,
      bagl::property_map_t<WorldGridType, bagl::vertex_position_t>, 4>;

  using WorldPartition2 = ReaK::pp::dvp_tree<
      VertexType, TopologyType,
      bagl::property_map_t<WorldGridType, bagl::vertex_position_t>, 2>;

  const std::vector<unsigned int> grid_sizes = {
      100,    200,    300,     400,     500,    800,   1000,  1100,
      1300,   1500,   1700,    1900,    2000,   2200,  2500,  3000,
      3500,   4000,   4500,    5000,    6000,   7000,  8000,  9000,
      10000,  12000,  15000,   20000,   25000,  30000, 50000, 100000,
      200000, 500000, 1000000, 2000000, 5000000};

  std::ofstream outFile("test_vp_results/dvp_umap_vecS_6.dat");
  outFile << "N\tVP2\tVP4\t (all times in micro-seconds per query per vertex)"
          << std::endl;

  for (unsigned int grid_size : grid_sizes) {
    WorldGridType grid;
    TopologyType m_space("",
                         ReaK::vect<double, 6>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                         ReaK::vect<double, 6>(1.0, 1.0, 1.0, 1.0, 1.0, 1.0));
    auto m_position = get(bagl::vertex_position, grid);

    for (unsigned int j = 0; j < grid_size; ++j) {
      VertexType v = add_vertex(grid);
      put(m_position, v, m_space.random_point());
    }

    outFile << grid_size;
    std::cout << "N = " << grid_size << std::endl;

    {
      WorldPartition2 part2_fresh(
          grid,
          std::shared_ptr<const TopologyType>(&m_space, ReaK::null_deleter()),
          m_position);

      ReaK::pp::multi_dvp_tree_search<WorldGridType, WorldPartition2>
          nn_finder2_fresh;
      nn_finder2_fresh.graph_tree_map[&grid] = &part2_fresh;
      high_resolution_clock::time_point t_start = high_resolution_clock::now();
      for (unsigned int j = 0; j < 1000; ++j) {
        nn_finder2_fresh(m_space.random_point(), grid, m_space, m_position);
      }
      high_resolution_clock::duration dt =
          high_resolution_clock::now() - t_start;
      outFile << "\t"
              << duration_cast<microseconds>(dt).count() * 0.001 /
                     double(grid_size);
      std::cout << "VP2-fresh" << std::endl;
    }

    {
      WorldPartition4 part4_fresh(
          grid,
          std::shared_ptr<const TopologyType>(&m_space, ReaK::null_deleter()),
          m_position);

      ReaK::pp::multi_dvp_tree_search<WorldGridType, WorldPartition4>
          nn_finder4_fresh;
      nn_finder4_fresh.graph_tree_map[&grid] = &part4_fresh;
      high_resolution_clock::time_point t_start = high_resolution_clock::now();
      for (unsigned int j = 0; j < 1000; ++j) {
        nn_finder4_fresh(m_space.random_point(), grid, m_space, m_position);
      }
      high_resolution_clock::duration dt =
          high_resolution_clock::now() - t_start;
      outFile << "\t"
              << duration_cast<microseconds>(dt).count() * 0.001 /
                     double(grid_size);
      std::cout << "VP4-fresh" << std::endl;
    }

    /*
    {
      ReaK::pp::linear_neighbor_search<WorldGridType> lnn_finder;
      high_resolution_clock::time_point t_start = high_resolution_clock::now();
      for(unsigned int j=0;j<1000;++j) {
        lnn_finder(m_space.random_point(),grid,m_space,m_position);
      };
      high_resolution_clock::duration dt = high_resolution_clock::now() - t_start;
      outFile << "\t" << duration_cast<microseconds>(dt).count() * 0.001 / double(grid_size);
      std::cout << "Linear Search" << std::endl;
    }
    */
    outFile << std::endl;
  }
}

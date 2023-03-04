
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

#include <boost/graph/adjacency_list_BC.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/topology.hpp>

#include <ReaK/math/lin_alg/vect_alg.hpp>
#include <ReaK/planning/path_planning/dvp_layout_adjacency_list.hpp>
#include <ReaK/planning/path_planning/topological_search.hpp>
#include <ReaK/topologies/spaces/hyperbox_topology.hpp>

#include <chrono>
#include <memory>

using TopologyType = ReaK::pp::hyperbox_topology<ReaK::vect<double, 6>>;

using PointType = TopologyType::point_type;

struct WorldGridVertexProperties {
  PointType pos;

  explicit WorldGridVertexProperties(PointType aPos = PointType())
      : pos(aPos){};
};

struct WorldGridEdgeProperties {};

int main() {

  using namespace std::chrono;

  using PositionMap =
      boost::data_member_property_map<PointType, WorldGridVertexProperties>;

  using WorldPartition2BF = ReaK::pp::dvp_adjacency_list<
      WorldGridVertexProperties, WorldGridEdgeProperties, TopologyType,
      PositionMap, 2, ReaK::pp::random_vp_chooser,
      boost::bfl_d_ary_tree_storage<2>, boost::vecBC, boost::undirectedS,
      boost::vecBC>;

  using WorldGrid2BF = WorldPartition2BF::adj_list_type;

  const std::vector<unsigned int> grid_sizes = {
      100,  200,  300,  400,  500,   800,   1000,  1100,  1300,  1500,
      1700, 1900, 2000, 2200, 2500,  3000,  3500,  4000,  4500,  5000,
      6000, 7000, 8000, 9000, 10000, 12000, 15000, 20000, 25000, 30000};

  std::ofstream outFile("test_vp_results/dvp_adj_list.dat");
  outFile << "N\tVP2\t (all times in micro-seconds per query per vertex)"
          << std::endl;

  for (unsigned int grid_size : grid_sizes) {
    std::shared_ptr<TopologyType> m_space = std::make_shared<TopologyType>(
        "", ReaK::vect<double, 6>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        ReaK::vect<double, 6>(1.0, 1.0, 1.0, 1.0, 1.0, 1.0));

    WorldPartition2BF dvp2(m_space,
                           PositionMap(&WorldGridVertexProperties::pos));
    WorldGrid2BF g2 = dvp2.get_adjacency_list();

    using RandSampler =
        ReaK::pp::point_distribution_traits<TopologyType>::random_sampler_type;
    RandSampler get_sample = get(ReaK::pp::random_sampler, *m_space);

    for (unsigned int j = 0; j < grid_size; ++j) {
      auto vp = WorldGridVertexProperties(get_sample(*m_space));
      add_vertex(vp, g2);
    };

    outFile << grid_size;
    std::cout << "N = " << grid_size << std::endl;

    {
      ReaK::pp::multi_dvp_tree_search<WorldGrid2BF, WorldPartition2BF>
          nn_finder2;
      nn_finder2.graph_tree_map[&g2] = &dvp2;
      high_resolution_clock::time_point t_start = high_resolution_clock::now();
      for (unsigned int j = 0; j < 1000; ++j) {
        nn_finder2(get_sample(*m_space), g2, *m_space,
                   get(&WorldGridVertexProperties::pos, g2));
      };
      high_resolution_clock::duration dt =
          high_resolution_clock::now() - t_start;
      outFile << "\t"
              << duration_cast<microseconds>(dt).count() * 0.001 /
                     double(grid_size);
      std::cout << "VP2-fresh" << std::endl;
    };

    outFile << std::endl;
  };
};

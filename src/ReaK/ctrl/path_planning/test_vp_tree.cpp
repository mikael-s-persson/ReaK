
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

#include <iostream>
#include <fstream>

#include <boost/graph/adjacency_list_BC.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/properties.hpp>

#include <ReaK/ctrl/path_planning/topological_search.hpp>
#include <ReaK/ctrl/path_planning/metric_space_search.hpp>
#include <ReaK/ctrl/topologies/hyperbox_topology.hpp>
#include <ReaK/core/lin_alg/vect_alg.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>


int main() {
  typedef ReaK::pp::hyperbox_topology< ReaK::vect<double,6> > TopologyType;
  
  typedef TopologyType::point_type PointType;

  typedef boost::property< boost::vertex_position_t, PointType, boost::no_property > WorldGridVertexProperties;

  typedef boost::no_property WorldGridEdgeProperties;

  typedef boost::adjacency_list_BC< boost::vecBC, boost::vecBC, boost::undirectedS,
                                    WorldGridVertexProperties, WorldGridEdgeProperties > WorldGridType;

  typedef boost::graph_traits<WorldGridType>::vertex_descriptor VertexType;
  typedef ReaK::pp::dvp_tree<VertexType, 
                             TopologyType, 
                             boost::property_map<WorldGridType, boost::vertex_position_t>::type, 
                             4> WorldPartition4;

  typedef ReaK::pp::dvp_tree<VertexType, 
                             TopologyType, 
                             boost::property_map<WorldGridType, boost::vertex_position_t>::type, 
                             2> WorldPartition2;
  
   const unsigned int grid_sizes[] = {100, 200, 300, 400, 500, 800, 1000, 1100, 1300, 1500, 1700, 
                                      1900, 2000, 2200, 2500, 3000, 3500, 4000, 4500, 5000, 6000,
                                      7000, 8000, 9000, 10000, 12000, 15000, 20000, 25000, 30000,
                                      50000, 100000, 200000, 500000, 1000000, 2000000, 5000000};
                                     
//  const unsigned int grid_sizes[] = {50000, 100000, 200000, 500000, 1000000,
//                                     2000000, 5000000, 10000000, 20000000};
  
  std::ofstream outFile("test_vp_results/dvp_umap_vecS_6.dat");
  outFile << "N\tVP2\tVP4\t (all times in micro-seconds per query per vertex)" << std::endl;
  
  for(int i=0;i<37;++i) {
    WorldGridType grid;
    TopologyType m_space("",ReaK::vect<double,6>(0.0,0.0,0.0,0.0,0.0,0.0),ReaK::vect<double,6>(1.0,1.0,1.0,1.0,1.0,1.0));
    boost::property_map<WorldGridType, boost::vertex_position_t>::type m_position(get(boost::vertex_position, grid));
    
    for(unsigned int j=0;j < grid_sizes[i];++j) {
      VertexType v = add_vertex(grid);
      put(m_position,v,m_space.random_point()); 
    };
    
    outFile << grid_sizes[i];
    std::cout << "N = " << grid_sizes[i] << std::endl;
    
    {
    WorldPartition2 part2_fresh(grid,ReaK::shared_ptr<const TopologyType>(&m_space,ReaK::null_deleter()),m_position);
    
    ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition2> nn_finder2_fresh;
    nn_finder2_fresh.graph_tree_map[&grid] = &part2_fresh;
    boost::posix_time::ptime t_start = boost::posix_time::microsec_clock::local_time();
    for(unsigned int j=0;j<1000;++j) {
      nn_finder2_fresh(m_space.random_point(),grid,m_space,m_position);
    };
    boost::posix_time::time_duration dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.001 / double(grid_sizes[i]);
    std::cout << "VP2-fresh" << std::endl;
    };
    
    {
    WorldPartition4 part4_fresh(grid,ReaK::shared_ptr<const TopologyType>(&m_space,ReaK::null_deleter()),m_position);
    
    ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition4> nn_finder4_fresh;
    nn_finder4_fresh.graph_tree_map[&grid] = &part4_fresh;
    boost::posix_time::ptime t_start = boost::posix_time::microsec_clock::local_time();
    for(unsigned int j=0;j<1000;++j) {
      nn_finder4_fresh(m_space.random_point(),grid,m_space,m_position);
    };
    boost::posix_time::time_duration dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.001 / double(grid_sizes[i]);
    std::cout << "VP4-fresh" << std::endl;
    };
    
    /*
    ReaK::pp::linear_neighbor_search<WorldGridType> lnn_finder;
    t_start = boost::posix_time::microsec_clock::local_time();
    for(unsigned int j=0;j<1000;++j) {
      lnn_finder(m_space.random_point(),grid,m_space,m_position);
    };
    dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.001 / double(grid_sizes[i]);
    std::cout << "Linear Search" << std::endl;
    */    
    outFile << std::endl;
  };

};














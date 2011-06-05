
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

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/properties.hpp>

#include "topological_search.hpp"
#include "metric_space_search.hpp"

#include <boost/date_time/posix_time/posix_time.hpp>


namespace boost {

  enum vertex_position_t { vertex_position };

  BOOST_INSTALL_PROPERTY(vertex, position);

};

int main() {
  typedef boost::hypercube_topology<6,boost::minstd_rand> TopologyType;
  
  typedef TopologyType::point_type PointType;

  typedef boost::property< boost::vertex_position_t, PointType, boost::no_property > WorldGridVertexProperties;

  typedef boost::no_property WorldGridEdgeProperties;

  typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::bidirectionalS,
                                 WorldGridVertexProperties,
	  	                 WorldGridEdgeProperties,
		                 boost::vecS> WorldGridType;

  typedef boost::adjacency_list_traits<boost::vecS,boost::vecS,boost::bidirectionalS,boost::vecS>::vertex_descriptor VertexType;
  typedef ReaK::pp::dvp_tree<VertexType, 
                             TopologyType, 
			     boost::property_map<WorldGridType, boost::vertex_position_t>::type, 
			     4> WorldPartition4;

  typedef ReaK::pp::dvp_tree<VertexType, 
                             TopologyType, 
			     boost::property_map<WorldGridType, boost::vertex_position_t>::type, 
			     2> WorldPartition2;

  const unsigned int grid_sizes[] = {100, 200, 500, 1000, 2000, 5000, 10000};
  
  std::ofstream outFile("test_vp_tree_results_6.dat");
  outFile << "N\tVP2\tVP2-f\tVP4\tVP4-f\tLS" << std::endl;
  
  for(int i=0;i<7;++i) {
    WorldGridType grid;
    boost::minstd_rand m_rng(std::time(0));
    TopologyType m_space(m_rng);
    boost::property_map<WorldGridType, boost::vertex_position_t>::type m_position(boost::get(boost::vertex_position, grid));
    
    WorldPartition2 part2(grid,m_space,m_position);
    WorldPartition4 part4(grid,m_space,m_position);
    
    for(int j=0;j<grid_sizes[i];++j) {
      VertexType v = boost::add_vertex(grid);
      boost::put(m_position,v,m_space.random_point());
      part2.insert(v);
      part4.insert(v);
    };
    
    WorldPartition2 part2_fresh(grid,m_space,m_position);
    WorldPartition4 part4_fresh(grid,m_space,m_position);
    
    outFile << grid_sizes[i];
    std::cout << "N = " << grid_sizes[i] << std::endl;
    
    ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition2> nn_finder2;
    nn_finder2.graph_tree_map[&grid] = &part2;
    boost::posix_time::ptime t_start = boost::posix_time::microsec_clock::local_time();
    for(unsigned int i=0;i<10000;++i) {
      nn_finder2(m_space.random_point(),grid,m_space,m_position);
    };
    boost::posix_time::time_duration dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.0001;
    std::cout << "VP2" << std::endl;
    
    ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition2> nn_finder2_fresh;
    nn_finder2_fresh.graph_tree_map[&grid] = &part2_fresh;
    t_start = boost::posix_time::microsec_clock::local_time();
    for(unsigned int i=0;i<10000;++i) {
      nn_finder2_fresh(m_space.random_point(),grid,m_space,m_position);
    };
    dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.0001;
    std::cout << "VP2-fresh" << std::endl;
    
    ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition4> nn_finder4;
    nn_finder4.graph_tree_map[&grid] = &part4;
    t_start = boost::posix_time::microsec_clock::local_time();
    for(unsigned int i=0;i<10000;++i) {
      nn_finder4(m_space.random_point(),grid,m_space,m_position);
    };
    dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.0001;
    std::cout << "VP4" << std::endl;
    
    ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition4> nn_finder4_fresh;
    nn_finder4_fresh.graph_tree_map[&grid] = &part4_fresh;
    t_start = boost::posix_time::microsec_clock::local_time();
    for(unsigned int i=0;i<10000;++i) {
      nn_finder4_fresh(m_space.random_point(),grid,m_space,m_position);
    };
    dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.0001;
    std::cout << "VP4-fresh" << std::endl;
    
    ReaK::pp::linear_neighbor_search<> lnn_finder;
    t_start = boost::posix_time::microsec_clock::local_time();
    for(unsigned int i=0;i<10000;++i) {
      lnn_finder(m_space.random_point(),grid,m_space,m_position);
    };
    dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.0001;
    std::cout << "Linear Search" << std::endl;
        
    outFile << std::endl;
  };

};














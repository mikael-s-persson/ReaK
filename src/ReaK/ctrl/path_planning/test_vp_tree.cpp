
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
#include "topologies/hyperbox_topology.hpp"
#include "lin_alg/vect_alg.hpp"

#include <boost/date_time/posix_time/posix_time.hpp>


namespace boost {

  enum vertex_position_t { vertex_position };

  BOOST_INSTALL_PROPERTY(vertex, position);

};

int main() {
  typedef ReaK::pp::hyperbox_topology< ReaK::vect<double,6> > TopologyType;
  
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
  
  const unsigned int grid_sizes[] = {100, 200, 300, 400, 500, 800, 1000, 1100, 1300, 1500, 1700, 
                                     1900, 2000, 2200, 2500, 3000, 3500, 4000, 4500, 5000, 6000,
				     7000, 8000, 9000, 10000, 12000, 15000, 20000, 25000};
  
  std::ofstream outFile("test_vp_results/dvp_umap_vecS_6.dat");
  outFile << "N\tVP2\tVP2d\tVP2e\tVP2-fg\tVP2-fd\tVP2-f\tVP2-fe\tVP4\tVP4d\tVP4e\tVP4-fg\tVP4-fd\tVP4-f\tVP4-fe\tLS\t (all times in micro-seconds per query per vertex)" << std::endl;
  
  for(int i=0;i<29;++i) {
    WorldGridType grid;
    std::size_t eval_count;
    TopologyType m_space("",ReaK::vect<double,6>(0.0,0.0,0.0,0.0,0.0,0.0),ReaK::vect<double,6>(1.0,1.0,1.0,1.0,1.0,1.0));
    boost::property_map<WorldGridType, boost::vertex_position_t>::type m_position(get(boost::vertex_position, grid));
    
    WorldPartition2 part2(grid,m_space,m_position);
    WorldPartition4 part4(grid,m_space,m_position);
    
    for(unsigned int j=0;j < grid_sizes[i];++j) {
      VertexType v = add_vertex(grid);
      put(m_position,v,m_space.random_point()); 
      part2.insert(v); 
      part4.insert(v); 
    };
    
    outFile << grid_sizes[i];
    std::cout << "N = " << grid_sizes[i] << std::endl;
    
    ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition2> nn_finder2;
    nn_finder2.graph_tree_map[&grid] = &part2;
    boost::posix_time::ptime t_start = boost::posix_time::microsec_clock::local_time();
    eval_count = 0;
    for(unsigned int j=0;j<10000;++j) {
      nn_finder2(m_space.random_point(),grid,m_space,m_position);
      eval_count += part2.get_dist_eval();
    };
    boost::posix_time::time_duration dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.0001 / double(grid_sizes[i]);
    std::cout << "VP2" << std::endl;
    
    outFile << "\t" << part2.depth();
    
    outFile << "\t" << eval_count * 0.0001;
    
    t_start = boost::posix_time::microsec_clock::local_time();
    WorldPartition2 part2_fresh(grid,m_space,m_position);
    dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.0001 / double(grid_sizes[i]);
    
    outFile << "\t" << part2_fresh.depth();
    
    ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition2> nn_finder2_fresh;
    nn_finder2_fresh.graph_tree_map[&grid] = &part2_fresh;
    t_start = boost::posix_time::microsec_clock::local_time();
    eval_count = 0;
    for(unsigned int j=0;j<10000;++j) {
      nn_finder2_fresh(m_space.random_point(),grid,m_space,m_position);
      eval_count += part2_fresh.get_dist_eval();
    };
    dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.0001 / double(grid_sizes[i]);
    std::cout << "VP2-fresh" << std::endl;
    
    outFile << "\t" << eval_count * 0.0001;
    
    ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition4> nn_finder4;
    nn_finder4.graph_tree_map[&grid] = &part4;
    t_start = boost::posix_time::microsec_clock::local_time();
    eval_count = 0;
    for(unsigned int j=0;j<10000;++j) {
      nn_finder4(m_space.random_point(),grid,m_space,m_position);
      eval_count += part4.get_dist_eval();
    };
    dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.0001 / double(grid_sizes[i]);
    std::cout << "VP4" << std::endl;
    
    outFile << "\t" << part4.depth();
    
    outFile << "\t" << eval_count * 0.0001;
    
    t_start = boost::posix_time::microsec_clock::local_time();
    WorldPartition4 part4_fresh(grid,m_space,m_position);
    dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.0001 / double(grid_sizes[i]);
    
    outFile << "\t" << part4_fresh.depth();
    
    ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition4> nn_finder4_fresh;
    nn_finder4_fresh.graph_tree_map[&grid] = &part4_fresh;
    t_start = boost::posix_time::microsec_clock::local_time();
    eval_count = 0;
    for(unsigned int j=0;j<10000;++j) {
      nn_finder4_fresh(m_space.random_point(),grid,m_space,m_position);
      eval_count += part4_fresh.get_dist_eval();
    };
    dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.0001 / double(grid_sizes[i]);
    std::cout << "VP4-fresh" << std::endl;
    
    outFile << "\t" << eval_count * 0.0001;
    
    ReaK::pp::linear_neighbor_search<> lnn_finder;
    t_start = boost::posix_time::microsec_clock::local_time();
    eval_count = 0;
    for(unsigned int j=0;j<10000;++j) {
      lnn_finder(m_space.random_point(),grid,m_space,m_position);
    };
    dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.0001 / double(grid_sizes[i]);
    std::cout << "Linear Search" << std::endl;
        
    outFile << std::endl;
  };

};














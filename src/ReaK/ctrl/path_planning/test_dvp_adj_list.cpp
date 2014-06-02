
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
#include <ReaK/ctrl/path_planning/dvp_layout_adjacency_list.hpp>
#include <ReaK/ctrl/topologies/hyperbox_topology.hpp>
#include <ReaK/core/lin_alg/vect_alg.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>


typedef ReaK::pp::hyperbox_topology< ReaK::vect<double,6> > TopologyType;
  
typedef TopologyType::point_type PointType;

struct WorldGridVertexProperties {
  PointType pos;
  
  WorldGridVertexProperties(PointType aPos = PointType()) : pos(aPos) { };
};

struct WorldGridEdgeProperties { };


int main() {
  
  typedef boost::data_member_property_map< PointType, WorldGridVertexProperties > PositionMap;
  
  typedef ReaK::pp::dvp_adjacency_list< 
    WorldGridVertexProperties,
    WorldGridEdgeProperties,
    TopologyType, PositionMap,
    2, ReaK::pp::random_vp_chooser,
    boost::bfl_d_ary_tree_storage<2>,
    boost::vecBC, boost::undirectedS, boost::vecBC > WorldPartition2BF;
  
  typedef WorldPartition2BF::adj_list_type WorldGrid2BF;
  
  const unsigned int grid_sizes[] = {100, 200, 300, 400, 500, 800, 1000, 1100, 1300, 1500, 1700, 
                                     1900, 2000, 2200, 2500, 3000, 3500, 4000, 4500, 5000, 6000,
                                     7000, 8000, 9000, 10000, 12000, 15000, 20000, 25000, 30000};
                                     
//  const unsigned int grid_sizes[] = {50000, 100000, 200000, 500000, 1000000,
//                                     2000000, 5000000, 10000000, 20000000};
  
  std::ofstream outFile("test_vp_results/dvp_adj_list.dat");
  outFile << "N\tVP2\t (all times in micro-seconds per query per vertex)" << std::endl;
  
  for(int i=0;i<30;++i) {
    ReaK::shared_ptr<TopologyType> m_space = ReaK::shared_ptr<TopologyType>(new TopologyType("",ReaK::vect<double,6>(0.0,0.0,0.0,0.0,0.0,0.0),ReaK::vect<double,6>(1.0,1.0,1.0,1.0,1.0,1.0)));
    
    WorldPartition2BF dvp2(m_space, PositionMap(&WorldGridVertexProperties::pos));
    WorldGrid2BF g2 = dvp2.get_adjacency_list();
    
    typedef ReaK::pp::point_distribution_traits< TopologyType >::random_sampler_type RandSampler;
    RandSampler get_sample = get(ReaK::pp::random_sampler, *m_space);
    
    for(unsigned int j=0;j < grid_sizes[i];++j) {
      WorldGridVertexProperties vp = WorldGridVertexProperties(get_sample(*m_space));
      add_vertex(vp,g2);
    };
    
    outFile << grid_sizes[i];
    std::cout << "N = " << grid_sizes[i] << std::endl;
    
    {
    ReaK::pp::multi_dvp_tree_search<WorldGrid2BF, WorldPartition2BF> nn_finder2;
    nn_finder2.graph_tree_map[&g2] = &dvp2;
    boost::posix_time::ptime t_start = boost::posix_time::microsec_clock::local_time();
    for(unsigned int j=0;j<1000;++j) {
      nn_finder2(get_sample(*m_space),g2,*m_space,get(&WorldGridVertexProperties::pos,g2));
    };
    boost::posix_time::time_duration dt = boost::posix_time::microsec_clock::local_time() - t_start;
    outFile << "\t" << dt.total_microseconds() * 0.001 / double(grid_sizes[i]);
    std::cout << "VP2-fresh" << std::endl;
    };
    
    outFile << std::endl;
  };

};














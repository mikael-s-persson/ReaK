
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

#include "path_planning/topological_search.hpp"
#include "path_planning/metric_space_search.hpp"
#include "lin_alg/vect_alg.hpp"

#include "IHAQR_topology.hpp"
#include "MEAQR_topology.hpp"
#include "quadrotor_system.hpp"
#include "topologies/se3_random_samplers.hpp"

#include "mbd_kte/kte_map_chain.hpp"

#include "serialization/xml_archiver.hpp"
#include "optimization/optim_exceptions.hpp"

#include <boost/date_time/posix_time/posix_time.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;


using namespace ReaK;
using namespace pp;
using namespace ctrl;


typedef IHAQR_topology< quadrotor_system::state_space_type, quadrotor_system, position_only_sampler > X8_IHAQR_space_type;
typedef MEAQR_topology< quadrotor_system::state_space_type, quadrotor_system, position_only_sampler > X8_MEAQR_space_type;

typedef X8_MEAQR_space_type::point_type MEAQR_PointType;

struct MEAQR_vprop {
  MEAQR_PointType position;
};


int main(int argc, char ** argv) {
  
  po::options_description generic_options("Generic options");
  generic_options.add_options()
    ("help,h", "produce this help message.")
  ;
  
  po::options_description io_options("I/O options");
  io_options.add_options()
    ("space-def-file,s", po::value< std::string >()->default_value("models/quadrotor_spaces.xml"), "specify the space-definition file (default is models/quadrotor_spaces.xml)")
    ("output-path,o", po::value< std::string >()->default_value("vp_results"), "specify the output path (default is vp_results)")
  ;
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options).add(io_options);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);
  
  if(vm.count("help")) {
    std::cout << cmdline_options << std::endl;
    return 1;
  };
  
  std::string output_path_name = vm["output-path"].as<std::string>();
  while(output_path_name[output_path_name.length()-1] == '/') 
    output_path_name.erase(output_path_name.length()-1, 1);
  
  fs::create_directory(output_path_name.c_str());
  
  
  std::string space_def_file_name = vm["space-def-file"].as<std::string>();
  std::string space_def_file_name_only(std::find(space_def_file_name.rbegin(),space_def_file_name.rend(),'/').base(), std::find(space_def_file_name.rbegin(),space_def_file_name.rend(),'.').base()-1);
  if( ! fs::exists(space_def_file_name) ) {
    std::cout << "Error: input file '" << space_def_file_name << "' does not exist!" << std::endl
              << "Correct usage of this program is:" << std::endl;
    std::cout << cmdline_options << std::endl;
    return 1;
  };
  
  shared_ptr< quadrotor_system >    X8_sys;
  shared_ptr< X8_IHAQR_space_type > X8_IHAQR_space;
  shared_ptr< X8_MEAQR_space_type > X8_MEAQR_space;
  {
    serialization::xml_iarchive in(space_def_file_name);
    in >> X8_sys >> X8_IHAQR_space >> X8_MEAQR_space;
  };
  
  
  typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
                                 MEAQR_vprop, boost::no_property, boost::no_property, boost::listS> WorldGridType;
  typedef boost::graph_traits< WorldGridType >::vertex_descriptor VertexType;
  typedef boost::property_map<WorldGridType, MEAQR_PointType MEAQR_vprop::* >::type PositionMapType;
  
  typedef ReaK::pp::dvp_tree<VertexType, X8_MEAQR_space_type, PositionMapType, 2> WorldPartition2;
  typedef ReaK::pp::dvp_tree<VertexType, X8_MEAQR_space_type, PositionMapType, 4> WorldPartition4;
  
  const unsigned int grid_sizes[] = {100, 200, 300, 400, 500, 600, 800, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 7500, 10000};
  
  std::string output_file_name = output_path_name + "/X8_knn_times.dat";
  std::ofstream outFile(output_file_name.c_str());
  outFile << "N\tVP2\tVP4\tLS\t (all times in micro-seconds per query)" << std::endl;
  
  for(unsigned int i = 0; i < sizeof(grid_sizes) / sizeof(const unsigned int); ++i) {
    WorldGridType grid;
    
    PositionMapType m_position(get(&MEAQR_vprop::position, grid));
    
    std::cout << "Generating " << grid_sizes[i] << " nodes..." << std::endl;
    for(unsigned int j = 0; j < grid_sizes[i]; ++j) {
      std::cout << "\r" << std::setw(10) << j << std::flush;
      VertexType v = add_vertex(grid);
      put(m_position, v, X8_MEAQR_space->random_point()); 
    };
    std::cout << std::endl << "Done!" << std::endl;
    
    outFile << grid_sizes[i];
    std::cout << "N = " << grid_sizes[i] << std::endl;
    
    {
      std::cout << "VP2 ..." << std::endl;
      WorldPartition2 part2_fresh(grid, X8_MEAQR_space, m_position);
      
      ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition2> nn_finder2_fresh;
      nn_finder2_fresh.graph_tree_map[&grid] = &part2_fresh;
      boost::posix_time::ptime t_start = boost::posix_time::microsec_clock::local_time();
      for(unsigned int j = 0; j < 1000; ++j) {
        nn_finder2_fresh(X8_MEAQR_space->random_point(), grid, *X8_MEAQR_space, m_position);
      };
      boost::posix_time::time_duration dt = boost::posix_time::microsec_clock::local_time() - t_start;
      outFile << "\t" << dt.total_microseconds() * 0.001;
      std::cout << "Done!" << std::endl;
    };
    
    {
      std::cout << "VP4 ..." << std::endl;
      WorldPartition4 part4_fresh(grid, X8_MEAQR_space, m_position);
      
      ReaK::pp::multi_dvp_tree_search<WorldGridType,WorldPartition4> nn_finder4_fresh;
      nn_finder4_fresh.graph_tree_map[&grid] = &part4_fresh;
      boost::posix_time::ptime t_start = boost::posix_time::microsec_clock::local_time();
      for(unsigned int j = 0; j < 1000; ++j) {
        nn_finder4_fresh(X8_MEAQR_space->random_point(), grid, *X8_MEAQR_space, m_position);
      };
      boost::posix_time::time_duration dt = boost::posix_time::microsec_clock::local_time() - t_start;
      outFile << "\t" << dt.total_microseconds() * 0.001;
      std::cout << "Done!" << std::endl;
    };
    
    
    {
      std::cout << "Linear Search ..." << std::endl;
      ReaK::pp::linear_neighbor_search<WorldGridType> lnn_finder;
      boost::posix_time::ptime t_start = boost::posix_time::microsec_clock::local_time();
      for(unsigned int j = 0; j < 1000; ++j) {
        lnn_finder(X8_MEAQR_space->random_point(), grid, *X8_MEAQR_space, m_position);
      };
      boost::posix_time::time_duration dt = boost::posix_time::microsec_clock::local_time() - t_start;
      outFile << "\t" << dt.total_microseconds() * 0.001;
      std::cout << "Done!" << std::endl;
    };
    
    outFile << std::endl;
  };

};














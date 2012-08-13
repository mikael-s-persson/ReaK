
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

#include "topologies/ptrobot2D_test_world.hpp"
#include "rrt_path_planner.hpp"
#include "prm_path_planner.hpp"
#include "rrtstar_path_planner.hpp"

#include "basic_sbmp_reporters.hpp"


int main(int argc, char** argv) {
  
  if(argc < 4) {
    std::cout << "Usage: ./test_planners [world_map.bmp] [number_of_runs] [max_vertices]" << std::endl
    << "\tworld_map.bmp:\t Image representing the pt-robot world." << std::endl
    << "\tnumber_of_runs:\t The number of runs to average out." << std::endl
    << "\tmax_vertices:\t The vertex limit count on the planner." << std::endl;
    return 0;
  };
  
  std::string world_file_name = argv[1];
  std::string world_file_name_only(std::find(world_file_name.rbegin(),world_file_name.rend(),'/').base(), std::find(world_file_name.rbegin(),world_file_name.rend(),'.').base()-1);
  
  std::size_t run_count = 0;
  std::stringstream(argv[2]) >> run_count;
  std::size_t max_vertices = 0;
  std::stringstream(argv[3]) >> max_vertices;
  std::size_t max_vertices_100 = max_vertices / 100;
  
  ReaK::shared_ptr< ReaK::pp::ptrobot2D_test_world > world_map =
    ReaK::shared_ptr< ReaK::pp::ptrobot2D_test_world >(new ReaK::pp::ptrobot2D_test_world(world_file_name, 10, 1.0));
  
//   ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> > > 
//     rrt_plan(world_map, 
// 	     world_map->get_start_pos(), 
// 	     world_map->get_goal_pos(),
// 	     10000, 
// 	     100,
// 	     ReaK::pp::BIDIRECTIONAL_RRT,
// 	     ReaK::pp::ADJ_LIST_MOTION_GRAPH,
// 	     ReaK::pp::DVP_BF2_TREE_KNN,
// 	     ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> >("pp_results/rrt/adstar_test_world_", 5),
// 	     50);
  
  
  std::ofstream timing_output("pp_results/" + world_file_name_only + "_times.txt");
  
  
  /**********************************************************************************
   * 
   * 
   *     Unidirectional Rapidly-exploring Random Tree Path-Planners
   * 
   * 
   * *******************************************************************************/
  
  std::cout << "Running RRT with Uni-dir, adj-list, dvp-bf2..." << std::endl;
  timing_output << "RRT, Uni-dir, adj-list, dvp-bf2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::UNIDIRECTIONAL_RRT,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_BF2_TREE_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT with Uni-dir, adj-list, dvp-bf4..." << std::endl;
  timing_output << "RRT, Uni-dir, adj-list, dvp-bf4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::UNIDIRECTIONAL_RRT,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_BF4_TREE_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT with Uni-dir, adj-list, dvp-cob2..." << std::endl;
  timing_output << "RRT, Uni-dir, adj-list, dvp-cob2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::UNIDIRECTIONAL_RRT,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_COB2_TREE_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT with Uni-dir, adj-list, dvp-cob4..." << std::endl;
  timing_output << "RRT, Uni-dir, adj-list, dvp-cob4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::UNIDIRECTIONAL_RRT,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_COB4_TREE_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT with Uni-dir, adj-list, linear-search..." << std::endl;
  timing_output << "RRT, Uni-dir, adj-list, linear-search" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::UNIDIRECTIONAL_RRT,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::LINEAR_SEARCH_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  
  std::cout << "Running RRT with Uni-dir, dvp-adj-list-bf2..." << std::endl;
  timing_output << "RRT, Uni-dir, dvp-adj-list-bf2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::UNIDIRECTIONAL_RRT,
               ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_ALT_BF2_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT with Uni-dir, dvp-adj-list-bf4..." << std::endl;
  timing_output << "RRT, Uni-dir, dvp-adj-list-bf4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::UNIDIRECTIONAL_RRT,
               ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_ALT_BF4_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT with Uni-dir, dvp-adj-list-cob2..." << std::endl;
  timing_output << "RRT, Uni-dir, dvp-adj-list-cob2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::UNIDIRECTIONAL_RRT,
               ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_ALT_COB2_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT with Uni-dir, dvp-adj-list-cob4..." << std::endl;
  timing_output << "RRT, Uni-dir, dvp-adj-list-cob4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::UNIDIRECTIONAL_RRT,
               ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_ALT_COB4_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Outputting RRT with Uni-dir, adj-list, dvp-bf2..." << std::endl;
  {
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> > > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::UNIDIRECTIONAL_RRT,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_BF2_TREE_KNN,
               ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> >("pp_results/rrt/" + world_file_name_only + "_", 5),
               10);
    
    rrt_plan.solve_path();
    
  };
  std::cout << "Done!" << std::endl;
  
  
  
  
  
  /**********************************************************************************
   * 
   * 
   *     Bidirectional Rapidly-exploring Random Tree Path-Planners
   * 
   * 
   * *******************************************************************************/
  
  
  std::cout << "Running RRT with Bi-dir, adj-list, dvp-bf2..." << std::endl;
  timing_output << "RRT, Bi-dir, adj-list, dvp-bf2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::BIDIRECTIONAL_RRT,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_BF2_TREE_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
      ss >> v_count >> t_val;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT with Bi-dir, adj-list, dvp-bf4..." << std::endl;
  timing_output << "RRT, Bi-dir, adj-list, dvp-bf4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::BIDIRECTIONAL_RRT,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_BF4_TREE_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
      ss >> v_count >> t_val;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT with Bi-dir, adj-list, dvp-cob2..." << std::endl;
  timing_output << "RRT, Bi-dir, adj-list, dvp-cob2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::BIDIRECTIONAL_RRT,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_COB2_TREE_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT with Bi-dir, adj-list, dvp-cob4..." << std::endl;
  timing_output << "RRT, Bi-dir, adj-list, dvp-cob4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::BIDIRECTIONAL_RRT,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_COB4_TREE_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
      ss >> v_count >> t_val;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT with Bi-dir, adj-list, linear-search..." << std::endl;
  timing_output << "RRT, Bi-dir, adj-list, linear-search" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::BIDIRECTIONAL_RRT,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::LINEAR_SEARCH_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1;
      ss >> v_count >> t_val;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  
  std::cout << "Running RRT with Bi-dir, dvp-adj-list-bf2..." << std::endl;
  timing_output << "RRT, Bi-dir, dvp-adj-list-bf2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::BIDIRECTIONAL_RRT,
               ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_ALT_BF2_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
      ss >> v_count >> t_val;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT with Bi-dir, dvp-adj-list-bf4..." << std::endl;
  timing_output << "RRT, Bi-dir, dvp-adj-list-bf4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::BIDIRECTIONAL_RRT,
               ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_ALT_BF4_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
      ss >> v_count >> t_val;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT with Bi-dir, dvp-adj-list-cob2..." << std::endl;
  timing_output << "RRT, Bi-dir, dvp-adj-list-cob2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::BIDIRECTIONAL_RRT,
               ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_ALT_COB2_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
      ss >> v_count >> t_val;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT with Bi-dir, dvp-adj-list-cob4..." << std::endl;
  timing_output << "RRT, Bi-dir, dvp-adj-list-cob4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::BIDIRECTIONAL_RRT,
               ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_ALT_COB4_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    rrt_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
      ss >> v_count >> t_val;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
    
  
  std::cout << "Outputting RRT with Bi-dir, adj-list, dvp-bf2..." << std::endl;
  {
    
    ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> > > 
      rrt_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::BIDIRECTIONAL_RRT,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_BF2_TREE_KNN,
               ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> >("pp_results/birrt/" + world_file_name_only + "_", 5),
               10);
    
    rrt_plan.solve_path();
    
  };
  std::cout << "Done!" << std::endl;
  
  
  
  
  
  
  
  /**********************************************************************************
   * 
   * 
   *     Probabilistic Roadmap Path-Planners
   * 
   * 
   * *******************************************************************************/
  
  
  std::cout << "Running PRM with adj-list, dvp-bf2..." << std::endl;
  timing_output << "PRM, adj-list, dvp-bf2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      prm_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_BF2_TREE_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    prm_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running PRM with adj-list, dvp-bf4..." << std::endl;
  timing_output << "PRM, adj-list, dvp-bf4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      prm_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_BF4_TREE_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    prm_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running PRM with adj-list, dvp-cob2..." << std::endl;
  timing_output << "PRM, adj-list, dvp-cob2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      prm_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_COB2_TREE_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    prm_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running PRM with adj-list, dvp-cob4..." << std::endl;
  timing_output << "PRM, adj-list, dvp-cob4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      prm_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_COB4_TREE_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    prm_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running PRM with adj-list, linear-search..." << std::endl;
  timing_output << "PRM, adj-list, linear-search" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      prm_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::LINEAR_SEARCH_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    prm_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  
  
  std::cout << "Running PRM with dvp-adj-list-bf2..." << std::endl;
  timing_output << "PRM, dvp-adj-list-bf2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      prm_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_ALT_BF2_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    prm_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running PRM with dvp-adj-list-bf4..." << std::endl;
  timing_output << "PRM, dvp-adj-list-bf4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      prm_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_ALT_BF4_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    prm_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running PRM with dvp-adj-list-cob2..." << std::endl;
  timing_output << "PRM, dvp-adj-list-cob2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      prm_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_ALT_COB2_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    prm_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running PRM with dvp-adj-list-cob4..." << std::endl;
  timing_output << "PRM, dvp-adj-list-cob4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      prm_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_ALT_COB4_KNN,
               ReaK::pp::timing_sbmp_report<>(ss),
               1);
    
    prm_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Outputting PRM with dvp-adj-list-cob4..." << std::endl;
  {
    
    ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> > > 
      prm_plan(world_map, 
               world_map->get_start_pos(), 
               world_map->get_goal_pos(),
               max_vertices, 
               100,
               ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
               ReaK::pp::DVP_ALT_COB4_KNN,
               ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> >("pp_results/prm/" + world_file_name_only + "_", 5),
               10);
    
    prm_plan.solve_path();
    
  };
  std::cout << "Done!" << std::endl;
  
  
  
  
  
  
  /**********************************************************************************
   * 
   * 
   *     Unidirectional Rapidly-exploring Random Tree Star Path-Planners
   * 
   * 
   * *******************************************************************************/
  
  
  std::cout << "Running RRT* with Uni-dir, adj-list, dvp-bf2..." << std::endl;
  timing_output << "RRT*, Uni-dir, adj-list, dvp-bf2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrtstar_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   max_vertices, 
                   100,
                   ReaK::pp::UNIDIRECTIONAL_RRT,
                   ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                   ReaK::pp::DVP_BF2_TREE_KNN,
                   ReaK::pp::timing_sbmp_report<>(ss),
                   1);
    
    rrtstar_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT* with Uni-dir, adj-list, dvp-bf4..." << std::endl;
  timing_output << "RRT*, Uni-dir, adj-list, dvp-bf4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrtstar_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   max_vertices, 
                   100,
                   ReaK::pp::UNIDIRECTIONAL_RRT,
                   ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                   ReaK::pp::DVP_BF4_TREE_KNN,
                   ReaK::pp::timing_sbmp_report<>(ss),
                   1);
    
    rrtstar_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT* with Uni-dir, adj-list, dvp-cob2..." << std::endl;
  timing_output << "RRT*, Uni-dir, adj-list, dvp-cob2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrtstar_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   max_vertices, 
                   100,
                   ReaK::pp::UNIDIRECTIONAL_RRT,
                   ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                   ReaK::pp::DVP_COB2_TREE_KNN,
                   ReaK::pp::timing_sbmp_report<>(ss),
                   1);
    
    rrtstar_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT* with Uni-dir, adj-list, dvp-cob4..." << std::endl;
  timing_output << "RRT*, Uni-dir, adj-list, dvp-cob4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrtstar_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   max_vertices, 
                   100,
                   ReaK::pp::UNIDIRECTIONAL_RRT,
                   ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                   ReaK::pp::DVP_COB4_TREE_KNN,
                   ReaK::pp::timing_sbmp_report<>(ss),
                   1);
    
    rrtstar_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT* with Uni-dir, adj-list, linear-search..." << std::endl;
  timing_output << "RRT*, Uni-dir, adj-list, linear-search" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrtstar_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   max_vertices, 
                   100,
                   ReaK::pp::UNIDIRECTIONAL_RRT,
                   ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                   ReaK::pp::LINEAR_SEARCH_KNN,
                   ReaK::pp::timing_sbmp_report<>(ss),
                   1);
    
    rrtstar_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  
  std::cout << "Running RRT* with Uni-dir, dvp-adj-list-bf2..." << std::endl;
  timing_output << "RRT*, Uni-dir, dvp-adj-list-bf2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrtstar_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   max_vertices, 
                   100,
                   ReaK::pp::UNIDIRECTIONAL_RRT,
                   ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                   ReaK::pp::DVP_ALT_BF2_KNN,
                   ReaK::pp::timing_sbmp_report<>(ss),
                   1);
    
    rrtstar_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT* with Uni-dir, dvp-adj-list-bf4..." << std::endl;
  timing_output << "RRT*, Uni-dir, dvp-adj-list-bf4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrtstar_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   max_vertices, 
                   100,
                   ReaK::pp::UNIDIRECTIONAL_RRT,
                   ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                   ReaK::pp::DVP_ALT_BF4_KNN,
                   ReaK::pp::timing_sbmp_report<>(ss),
                   1);
    
    rrtstar_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT* with Uni-dir, dvp-adj-list-cob2..." << std::endl;
  timing_output << "RRT*, Uni-dir, dvp-adj-list-cob2" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrtstar_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   max_vertices, 
                   100,
                   ReaK::pp::UNIDIRECTIONAL_RRT,
                   ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                   ReaK::pp::DVP_ALT_COB2_KNN,
                   ReaK::pp::timing_sbmp_report<>(ss),
                   1);
    
    rrtstar_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  std::cout << "Running RRT* with Uni-dir, dvp-adj-list-cob4..." << std::endl;
  timing_output << "RRT*, Uni-dir, dvp-adj-list-cob4" << std::endl;
  {
  std::vector< std::pair<int, double> > avg_times(max_vertices_100, std::pair<int, double>(0,0.0));
  for(std::size_t i = 0; i < run_count; ++i) {
    std::stringstream ss;
    
    ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::timing_sbmp_report<> > 
      rrtstar_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   max_vertices, 
                   100,
                   ReaK::pp::UNIDIRECTIONAL_RRT,
                   ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                   ReaK::pp::DVP_ALT_COB4_KNN,
                   ReaK::pp::timing_sbmp_report<>(ss),
                   1);
    
    rrtstar_plan.solve_path();
    
    int v_count, t_val;
    int j = 0;
    while(ss >> v_count) {
      ss >> t_val;
      avg_times[j].second = (double(t_val) + double(avg_times[j].first) * avg_times[j].second) / double(avg_times[j].first + 1);
      avg_times[j].first += 1; ++j;
    };
  };
  for(std::size_t i = 0; i < max_vertices_100; ++i) {
    if(avg_times[i].first)
      timing_output << std::setw(6) << (i+1)*100 << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second << std::endl; 
  };
  };
  std::cout << "Done!" << std::endl;
  
  
  
  std::cout << "Outputting RRT* with Uni-dir, adj-list, dvp-bf4..." << std::endl;
  {
    
    ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> > > 
      rrtstar_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   max_vertices, 
                   100,
                   ReaK::pp::UNIDIRECTIONAL_RRT,
                   ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                   ReaK::pp::DVP_BF4_TREE_KNN,
                   ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> >("pp_results/rrt_star/" + world_file_name_only + "_", 5),
                   10);
    
    rrtstar_plan.solve_path();
    
  };
  std::cout << "Done!" << std::endl;
  
  
  
  return 0;
};














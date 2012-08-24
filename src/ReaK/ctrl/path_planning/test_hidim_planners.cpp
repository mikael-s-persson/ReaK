
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

#include "topologies/hyperbox_topology.hpp"
#include "topologies/se3_topologies.hpp"
#include "topologies/no_obstacle_space.hpp"
#include "interpolation/sustained_velocity_pulse.hpp"
#include "rrt_path_planner.hpp"
#include "prm_path_planner.hpp"
#include "rrtstar_path_planner.hpp"

#include "basic_sbmp_reporters.hpp"




template <typename SpaceType>
void test_planners_on_space(ReaK::shared_ptr< SpaceType > world_map, 
                            std::ostream& timing_output,
                            std::size_t run_count, std::size_t max_vertices) {
  
  std::size_t max_vertices_100 = max_vertices / 100;
  
  std::cout << "*****************************************************************" << std::endl
            << "*            Running tests on '" << world_map->getName() << "'" << std::endl
            << "*****************************************************************" << std::endl;
  
#if 0
            
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrt_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::prm_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::prm_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::prm_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::prm_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::prm_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::prm_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::prm_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::prm_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::prm_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
  
#endif
  
  
  
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrtstar_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrtstar_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrtstar_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrtstar_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrtstar_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrtstar_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrtstar_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrtstar_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
    world_map->set_start_pos(world_map->random_point());
    world_map->set_goal_pos(world_map->random_point());
    
    ReaK::pp::rrtstar_path_planner< SpaceType, ReaK::pp::timing_sbmp_report<> > 
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
  
  
  
  
  
};










int main(int argc, char** argv) {
  
  if(argc < 3) {
    std::cout << "Usage: ./test_hidim_planners [number_of_runs] [max_vertices]" << std::endl
    << "\tnumber_of_runs:\t The number of runs to average out." << std::endl
    << "\tmax_vertices:\t The vertex limit count on the planner." << std::endl;
    return 0;
  };
  
  std::size_t run_count = 0;
  std::stringstream(argv[1]) >> run_count;
  std::size_t max_vertices = 0;
  std::stringstream(argv[2]) >> max_vertices;
  
  
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double, 6> > > World6DType;
  ReaK::shared_ptr< World6DType > world_6D =
    ReaK::shared_ptr< World6DType >(
      new World6DType(
        "world_6D_no_obstacles",
        ReaK::pp::hyperbox_topology< ReaK::vect<double, 6> >("world_6D",
                                                             ReaK::vect<double,6>(0.0,0.0,0.0,0.0,0.0,0.0),
                                                             ReaK::vect<double,6>(1.0,1.0,1.0,1.0,1.0,1.0)),
        0.1));
  
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double, 12> > > World12DType;
  ReaK::shared_ptr< World12DType > world_12D =
    ReaK::shared_ptr< World12DType >(
      new World12DType(
        "world_12D_no_obstacles",
        ReaK::pp::hyperbox_topology< ReaK::vect<double, 12> >("world_12D",
                                                              ReaK::vect<double,12>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                              ReaK::vect<double,12>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
        0.5));
  
  typedef ReaK::pp::no_obstacle_space< typename ReaK::pp::rl_se3_1st_order_topology<double>::type > WorldRLSE3Type;
  ReaK::shared_ptr< WorldRLSE3Type > world_RLSE3 =
    ReaK::shared_ptr< WorldRLSE3Type >(
      new WorldRLSE3Type(
        "world_RLSE3_no_obstacles",
        ReaK::pp::make_rl_se3_space("world_RLSE3", 
                                    ReaK::vect<double,3>(-10.0, -10.0, -10.0),
                                    ReaK::vect<double,3>( 10.0,  10.0,  10.0),
                                    2.0,
                                    1.5,
                                    0.2,
                                    0.1),
        0.5));
          
  {
    std::ofstream timing_output("pp_results/" + world_6D->getName() + "_times.txt");
    
    test_planners_on_space(world_6D, timing_output, run_count, max_vertices);
  };
  
  {
    std::ofstream timing_output("pp_results/" + world_12D->getName() + "_times.txt");
    
    test_planners_on_space(world_12D, timing_output, run_count, max_vertices);
  };
  
  {
    std::ofstream timing_output("pp_results/" + world_RLSE3->getName() + "_times.txt");
    
    test_planners_on_space(world_RLSE3, timing_output, run_count, max_vertices);
  };
  
  
  
  return 0;
};














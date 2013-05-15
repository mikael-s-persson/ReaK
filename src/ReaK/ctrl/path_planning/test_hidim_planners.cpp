
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


// #define RK_ENABLE_TEST_URRT_PLANNER
// #define RK_ENABLE_TEST_BRRT_PLANNER
#define RK_ENABLE_TEST_RRTSTAR_PLANNER
// #define RK_ENABLE_TEST_PRM_PLANNER
// #define RK_ENABLE_TEST_FADPRM_PLANNER
#define RK_ENABLE_TEST_SBASTAR_PLANNER


#if defined(RK_ENABLE_TEST_URRT_PLANNER) || defined(RK_ENABLE_TEST_BRRT_PLANNER)
#include "rrt_path_planner.hpp"
#endif

#if defined(RK_ENABLE_TEST_PRM_PLANNER)
#include "prm_path_planner.hpp"
#endif

#if defined(RK_ENABLE_TEST_RRTSTAR_PLANNER)
#include "rrtstar_path_planner.hpp"
#endif

#if defined(RK_ENABLE_TEST_FADPRM_PLANNER)
#include "fadprm_path_planner.hpp"
#endif

#if defined(RK_ENABLE_TEST_SBASTAR_PLANNER)
#include "sbastar_path_planner.hpp"
#endif


#include "basic_sbmp_reporters.hpp"


#include <boost/program_options.hpp>

namespace po = boost::program_options;






#define RK_HIDIM_PLANNER_DO_RRT     0x00000001
#define RK_HIDIM_PLANNER_DO_BIRRT   0x00000002
#define RK_HIDIM_PLANNER_DO_PRM     0x00000004
#define RK_HIDIM_PLANNER_DO_FADPRM  0x00000008
#define RK_HIDIM_PLANNER_DO_SBASTAR 0x00000010
#define RK_HIDIM_PLANNER_DO_RRTSTAR 0x00000020

#define RK_HIDIM_PLANNER_DO_COB 0x00000100
#define RK_HIDIM_PLANNER_DO_ALT 0x00000200



double sba_potential_cutoff;
double sba_density_cutoff;
double sba_relaxation;
double sba_sa_temperature;
bool sba_use_voronoi_pull;


template <typename SpaceType>
void run_monte_carlo_tests(
    std::size_t mc_run_count,
    std::size_t mc_num_records,
    ReaK::pp::sample_based_planner< ReaK::pp::path_planner_base<SpaceType> >& planner,
    std::stringstream& time_rec_ss,
    std::stringstream& cost_rec_ss,
    std::ostream& result_output) {
  std::vector< std::size_t > vertex_counts(mc_num_records, 0);
  std::vector< std::size_t > num_remaining_planners(mc_num_records, 0);
  std::vector< std::size_t > num_successful_planners(mc_num_records, 0);
  
  std::vector< double > time_values(mc_num_records, 0.0);
  std::vector< double > best_costs(mc_num_records, 0.0);
  std::vector< double > worst_costs(mc_num_records, 1.0e10);
  std::vector< double > avg_costs(mc_num_records, 0.0);
  
  for(std::size_t i = 0; i < mc_run_count; ++i) {
    
    planner.solve_path();
    
    std::size_t v_count, t_val; 
    std::size_t j = 0;
    while(time_rec_ss >> v_count) {
      time_rec_ss >> t_val;
      vertex_counts[j] = v_count;
      time_values[j] = (double(t_val) + double(num_remaining_planners[j]) * time_values[j]) / double(num_remaining_planners[j] + 1);
      num_remaining_planners[j] += 1; 
      ++j;
    };
    
    double c_val;
    j = 0;
    while(cost_rec_ss >> v_count) {
      cost_rec_ss >> c_val;
      if(c_val < best_costs[j])
        best_costs[j] = c_val;
      if(c_val > worst_costs[j])
        worst_costs[j] = c_val;
      if(c_val < 1.0e9) {
        avg_costs[j] = (double(c_val) + double(num_successful_planners[j]) * avg_costs[j]) / double(num_successful_planners[j] + 1);
        num_successful_planners[j] += 1;
      };
      ++j;
    };
    
    while(j < mc_num_records) {
      if(c_val < best_costs[j])
        best_costs[j] = c_val;
      if(c_val > worst_costs[j])
        worst_costs[j] = c_val;
      if(c_val < 1.0e9) {
        avg_costs[j] = (double(c_val) + double(num_successful_planners[j]) * avg_costs[j]) / double(num_successful_planners[j] + 1);
        num_successful_planners[j] += 1;
      };
      ++j;
    };
  };
  for(std::size_t i = 0; i < mc_num_records; ++i) {
    result_output << std::setw(9) << vertex_counts[i] 
           << " " << std::setw(9) << num_remaining_planners[i] 
           << " " << std::setw(9) << num_successful_planners[i] 
           << " " << std::setw(9) << time_values[i] 
           << " " << std::setw(9) << best_costs[i] 
           << " " << std::setw(9) << worst_costs[i] 
           << " " << std::setw(9) << avg_costs[i] << std::endl; 
  };
};






template <typename SpaceType>
void test_planners_on_space(ReaK::shared_ptr< SpaceType > world_map, 
                            std::ostream& timing_output,
                            std::size_t mc_run_count, 
                            std::size_t mc_max_vertices, 
                            std::size_t mc_prog_interval,
                            std::size_t mc_results,
                            std::size_t mc_flags) {
  
  std::size_t mc_max_vertices_100 = mc_max_vertices / mc_prog_interval;
  
  std::cout << "*****************************************************************" << std::endl
            << "*            Running tests on '" << world_map->getName() << "'" << std::endl
            << "*****************************************************************" << std::endl;
  
  
  typedef ReaK::pp::timing_sbmp_report< ReaK::pp::least_cost_sbmp_report<> > ReporterType;
  
#ifdef RK_ENABLE_TEST_URRT_PLANNER
  
  if(mc_flags & RK_HIDIM_PLANNER_DO_RRT) {
    
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
      std::stringstream ss, ss2;
      
      ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
        rrt_plan(world_map, 
                 world_map->get_start_pos(), 
                 world_map->get_goal_pos(),
                 mc_max_vertices, 
                 mc_prog_interval,
                 ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF2_TREE_KNN,
                 ReaK::pp::UNIDIRECTIONAL_PLANNING,
                 ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                 mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    std::cout << "Running RRT with Uni-dir, adj-list, dvp-bf4..." << std::endl;
    timing_output << "RRT, Uni-dir, adj-list, dvp-bf4" << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
        rrt_plan(world_map, 
                 world_map->get_start_pos(), 
                 world_map->get_goal_pos(),
                 mc_max_vertices, 
                 mc_prog_interval,
                 ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF4_TREE_KNN,
                 ReaK::pp::UNIDIRECTIONAL_PLANNING,
                 ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                 mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    if(mc_flags & RK_HIDIM_PLANNER_DO_COB) {
      std::cout << "Running RRT with Uni-dir, adj-list, dvp-cob2..." << std::endl;
      timing_output << "RRT, Uni-dir, adj-list, dvp-cob2" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
          rrt_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   mc_max_vertices, 
                   mc_prog_interval,
                   ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB2_TREE_KNN,
                   ReaK::pp::UNIDIRECTIONAL_PLANNING,
                   ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                   mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      std::cout << "Running RRT with Uni-dir, adj-list, dvp-cob4..." << std::endl;
      timing_output << "RRT, Uni-dir, adj-list, dvp-cob4" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
          rrt_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   mc_max_vertices, 
                   mc_prog_interval,
                   ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB4_TREE_KNN,
                   ReaK::pp::UNIDIRECTIONAL_PLANNING,
                   ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                   mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
    };
    
    
    std::cout << "Running RRT with Uni-dir, adj-list, linear-search..." << std::endl;
    timing_output << "RRT, Uni-dir, adj-list, linear-search" << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
        rrt_plan(world_map, 
                 world_map->get_start_pos(), 
                 world_map->get_goal_pos(),
                 mc_max_vertices, 
                 mc_prog_interval,
                 ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::LINEAR_SEARCH_KNN,
                 ReaK::pp::UNIDIRECTIONAL_PLANNING,
                 ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                 mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    if(mc_flags & RK_HIDIM_PLANNER_DO_ALT) {
      std::cout << "Running RRT with Uni-dir, dvp-adj-list-bf2..." << std::endl;
      timing_output << "RRT, Uni-dir, dvp-adj-list-bf2" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
          rrt_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   mc_max_vertices, 
                   mc_prog_interval,
                   ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF2_TREE_KNN,
                   ReaK::pp::UNIDIRECTIONAL_PLANNING,
                   ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                   mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT with Uni-dir, dvp-adj-list-bf4..." << std::endl;
      timing_output << "RRT, Uni-dir, dvp-adj-list-bf4" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
          rrt_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   mc_max_vertices, 
                   mc_prog_interval,
                   ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF4_TREE_KNN,
                   ReaK::pp::UNIDIRECTIONAL_PLANNING,
                   ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                   mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      if(mc_flags & RK_HIDIM_PLANNER_DO_COB) {
        std::cout << "Running RRT with Uni-dir, dvp-adj-list-cob2..." << std::endl;
        timing_output << "RRT, Uni-dir, dvp-adj-list-cob2" << std::endl;
        {
          std::stringstream ss, ss2;
          
          ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
            rrt_plan(world_map, 
                     world_map->get_start_pos(), 
                     world_map->get_goal_pos(),
                     mc_max_vertices, 
                     mc_prog_interval,
                     ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB2_TREE_KNN,
                     ReaK::pp::UNIDIRECTIONAL_PLANNING,
                     ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                     mc_results);
          
          run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
        };
        std::cout << "Done!" << std::endl;
        
        
        std::cout << "Running RRT with Uni-dir, dvp-adj-list-cob4..." << std::endl;
        timing_output << "RRT, Uni-dir, dvp-adj-list-cob4" << std::endl;
        {
          std::stringstream ss, ss2;
          
          ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
            rrt_plan(world_map, 
                     world_map->get_start_pos(), 
                     world_map->get_goal_pos(),
                     mc_max_vertices, 
                     mc_prog_interval,
                     ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB4_TREE_KNN,
                     ReaK::pp::UNIDIRECTIONAL_PLANNING,
                     ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                     mc_results);
          
          run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
        };
        std::cout << "Done!" << std::endl;
      };
    };
  };
  
#endif
  
  
#ifdef RK_ENABLE_TEST_BRRT_PLANNER
  
  
  if(mc_flags & RK_HIDIM_PLANNER_DO_BIRRT) {
    
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
      std::stringstream ss, ss2;
      
      ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
        rrt_plan(world_map, 
                 world_map->get_start_pos(), 
                 world_map->get_goal_pos(),
                 mc_max_vertices, 
                 mc_prog_interval,
                 ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF2_TREE_KNN,
                 ReaK::pp::BIDIRECTIONAL_PLANNING,
                 ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                 mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    std::cout << "Running RRT with Bi-dir, adj-list, dvp-bf4..." << std::endl;
    timing_output << "RRT, Bi-dir, adj-list, dvp-bf4" << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
        rrt_plan(world_map, 
                 world_map->get_start_pos(), 
                 world_map->get_goal_pos(),
                 mc_max_vertices, 
                 mc_prog_interval,
                 ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF4_TREE_KNN,
                 ReaK::pp::BIDIRECTIONAL_PLANNING,
                 ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                 mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    if(mc_flags & RK_HIDIM_PLANNER_DO_COB) {
      std::cout << "Running RRT with Bi-dir, adj-list, dvp-cob2..." << std::endl;
      timing_output << "RRT, Bi-dir, adj-list, dvp-cob2" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
          rrt_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   mc_max_vertices, 
                   mc_prog_interval,
                   ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB2_TREE_KNN,
                   ReaK::pp::BIDIRECTIONAL_PLANNING,
                   ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                   mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT with Bi-dir, adj-list, dvp-cob4..." << std::endl;
      timing_output << "RRT, Bi-dir, adj-list, dvp-cob4" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
          rrt_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   mc_max_vertices, 
                   mc_prog_interval,
                   ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB4_TREE_KNN,
                   ReaK::pp::BIDIRECTIONAL_PLANNING,
                   ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                   mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
    };
    
    
    std::cout << "Running RRT with Bi-dir, adj-list, linear-search..." << std::endl;
    timing_output << "RRT, Bi-dir, adj-list, linear-search" << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
        rrt_plan(world_map, 
                 world_map->get_start_pos(), 
                 world_map->get_goal_pos(),
                 mc_max_vertices, 
                 mc_prog_interval,
                 ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::LINEAR_SEARCH_KNN,
                 ReaK::pp::BIDIRECTIONAL_PLANNING,
                 ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                 mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    
    if(mc_flags & RK_HIDIM_PLANNER_DO_ALT) {
      std::cout << "Running RRT with Bi-dir, dvp-adj-list-bf2..." << std::endl;
      timing_output << "RRT, Bi-dir, dvp-adj-list-bf2" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
          rrt_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   mc_max_vertices, 
                   mc_prog_interval,
                   ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF2_TREE_KNN,
                   ReaK::pp::BIDIRECTIONAL_PLANNING,
                   ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                   mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT with Bi-dir, dvp-adj-list-bf4..." << std::endl;
      timing_output << "RRT, Bi-dir, dvp-adj-list-bf4" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
          rrt_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   mc_max_vertices, 
                   mc_prog_interval,
                   ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF4_TREE_KNN,
                   ReaK::pp::BIDIRECTIONAL_PLANNING,
                   ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                   mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      if(mc_flags & RK_HIDIM_PLANNER_DO_COB) {
        std::cout << "Running RRT with Bi-dir, dvp-adj-list-cob2..." << std::endl;
        timing_output << "RRT, Bi-dir, dvp-adj-list-cob2" << std::endl;
        {
          std::stringstream ss, ss2;
          
          ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
            rrt_plan(world_map, 
                     world_map->get_start_pos(), 
                     world_map->get_goal_pos(),
                     mc_max_vertices, 
                     mc_prog_interval,
                     ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB2_TREE_KNN,
                     ReaK::pp::BIDIRECTIONAL_PLANNING,
                     ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                     mc_results);
          
          run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
        };
        std::cout << "Done!" << std::endl;
        
        
        std::cout << "Running RRT with Bi-dir, dvp-adj-list-cob4..." << std::endl;
        timing_output << "RRT, Bi-dir, dvp-adj-list-cob4" << std::endl;
        {
          std::stringstream ss, ss2;
          
          ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
            rrt_plan(world_map, 
                     world_map->get_start_pos(), 
                     world_map->get_goal_pos(),
                     mc_max_vertices, 
                     mc_prog_interval,
                     ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB4_TREE_KNN,
                     ReaK::pp::BIDIRECTIONAL_PLANNING,
                     ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                     mc_results);
          
          run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrt_plan,ss,ss2,timing_output);
        };
        std::cout << "Done!" << std::endl;
      };
    };
  };
  
#endif
  
  
  
#ifdef RK_ENABLE_TEST_PRM_PLANNER
  
  
  if(mc_flags & RK_HIDIM_PLANNER_DO_PRM) {
    
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
      std::stringstream ss, ss2;
      
      ReaK::pp::prm_path_planner< SpaceType, ReporterType > 
        prm_plan(world_map, 
                 world_map->get_start_pos(), 
                 world_map->get_goal_pos(),
                 mc_max_vertices, 
                 mc_prog_interval,
                 ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF2_TREE_KNN,
                 ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                 mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,prm_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    std::cout << "Running PRM with adj-list, dvp-bf4..." << std::endl;
    timing_output << "PRM, adj-list, dvp-bf4" << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::prm_path_planner< SpaceType, ReporterType > 
        prm_plan(world_map, 
                 world_map->get_start_pos(), 
                 world_map->get_goal_pos(),
                 mc_max_vertices, 
                 mc_prog_interval,
                 ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF4_TREE_KNN,
                 ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                 mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,prm_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    if(mc_flags & RK_HIDIM_PLANNER_DO_COB) {
      std::cout << "Running PRM with adj-list, dvp-cob2..." << std::endl;
      timing_output << "PRM, adj-list, dvp-cob2" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::prm_path_planner< SpaceType, ReporterType > 
          prm_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   mc_max_vertices, 
                   mc_prog_interval,
                   ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB2_TREE_KNN,
                   ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                   mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,prm_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running PRM with adj-list, dvp-cob4..." << std::endl;
      timing_output << "PRM, adj-list, dvp-cob4" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::prm_path_planner< SpaceType, ReporterType > 
          prm_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   mc_max_vertices, 
                   mc_prog_interval,
                   ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB4_TREE_KNN,
                   ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                   mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,prm_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
    };
    
    
    std::cout << "Running PRM with adj-list, linear-search..." << std::endl;
    timing_output << "PRM, adj-list, linear-search" << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::prm_path_planner< SpaceType, ReporterType > 
        prm_plan(world_map, 
                 world_map->get_start_pos(), 
                 world_map->get_goal_pos(),
                 mc_max_vertices, 
                 mc_prog_interval,
                 ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::LINEAR_SEARCH_KNN,
                 ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                 mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,prm_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    
    
    if(mc_flags & RK_HIDIM_PLANNER_DO_ALT) {
      std::cout << "Running PRM with dvp-adj-list-bf2..." << std::endl;
      timing_output << "PRM, dvp-adj-list-bf2" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::prm_path_planner< SpaceType, ReporterType > 
          prm_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   mc_max_vertices, 
                   mc_prog_interval,
                   ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF2_TREE_KNN,
                   ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                   mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,prm_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running PRM with dvp-adj-list-bf4..." << std::endl;
      timing_output << "PRM, dvp-adj-list-bf4" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::prm_path_planner< SpaceType, ReporterType > 
          prm_plan(world_map, 
                   world_map->get_start_pos(), 
                   world_map->get_goal_pos(),
                   mc_max_vertices, 
                   mc_prog_interval,
                   ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF4_TREE_KNN,
                   ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                   mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,prm_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      if(mc_flags & RK_HIDIM_PLANNER_DO_COB) {
        std::cout << "Running PRM with dvp-adj-list-cob2..." << std::endl;
        timing_output << "PRM, dvp-adj-list-cob2" << std::endl;
        {
          std::stringstream ss, ss2;
          
          ReaK::pp::prm_path_planner< SpaceType, ReporterType > 
            prm_plan(world_map, 
                     world_map->get_start_pos(), 
                     world_map->get_goal_pos(),
                     mc_max_vertices, 
                     mc_prog_interval,
                     ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB2_TREE_KNN,
                     ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                     mc_results);
          
          run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,prm_plan,ss,ss2,timing_output);
        };
        std::cout << "Done!" << std::endl;
        
        
        std::cout << "Running PRM with dvp-adj-list-cob4..." << std::endl;
        timing_output << "PRM, dvp-adj-list-cob4" << std::endl;
        {
          std::stringstream ss, ss2;
          
          ReaK::pp::prm_path_planner< SpaceType, ReporterType > 
            prm_plan(world_map, 
                     world_map->get_start_pos(), 
                     world_map->get_goal_pos(),
                     mc_max_vertices, 
                     mc_prog_interval,
                     ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB4_TREE_KNN,
                     ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                     mc_results);
          
          run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,prm_plan,ss,ss2,timing_output);
        };
        std::cout << "Done!" << std::endl;
      };
    };
  };
  
#endif
  
  
  
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
  
  if(mc_flags & RK_HIDIM_PLANNER_DO_FADPRM) {
    
    /**********************************************************************************
    * 
    * 
    *     Flexible Anytime-Dynamic Probabilistic Roadmap Path-Planners
    * 
    * 
    * *******************************************************************************/

    
    std::cout << "Running FADPRM with adj-list, dvp-bf2..." << std::endl;
    timing_output << "FADPRM, adj-list, dvp-bf2" << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::fadprm_path_planner< SpaceType, ReporterType > 
        fadprm_plan(
          world_map, 
          world_map->get_start_pos(),
          world_map->get_goal_pos(),
          10.0,
          mc_max_vertices, 
          mc_prog_interval,
          ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF2_TREE_KNN,
          ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
          mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,fadprm_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    
    std::cout << "Running FADPRM with adj-list, dvp-bf4..." << std::endl;
    timing_output << "FADPRM, adj-list, dvp-bf4" << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::fadprm_path_planner< SpaceType, ReporterType > 
        fadprm_plan(
          world_map, 
          world_map->get_start_pos(), 
          world_map->get_goal_pos(),
          10.0,
          mc_max_vertices, 
          mc_prog_interval,
          ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF4_TREE_KNN,
          ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
          mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,fadprm_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    if(mc_flags & RK_HIDIM_PLANNER_DO_COB) {
      std::cout << "Running FADPRM with adj-list, dvp-cob2..." << std::endl;
      timing_output << "FADPRM, adj-list, dvp-cob2" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::fadprm_path_planner< SpaceType, ReporterType > 
          fadprm_plan(
            world_map, 
            world_map->get_start_pos(), 
            world_map->get_goal_pos(),
            10.0,
            mc_max_vertices, 
            mc_prog_interval,
            ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB2_TREE_KNN,
            ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
            mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,fadprm_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running FADPRM with adj-list, dvp-cob4..." << std::endl;
      timing_output << "FADPRM, adj-list, dvp-cob4" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::fadprm_path_planner< SpaceType, ReporterType > 
          fadprm_plan(
            world_map, 
            world_map->get_start_pos(), 
            world_map->get_goal_pos(),
            10.0,
            mc_max_vertices, 
            mc_prog_interval,
            ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB4_TREE_KNN,
            ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
            mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,fadprm_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
    };
    
    
    std::cout << "Running FADPRM with adj-list, linear-search..." << std::endl;
    timing_output << "FADPRM, adj-list, linear-search" << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::fadprm_path_planner< SpaceType, ReporterType > 
        fadprm_plan(
          world_map, 
          world_map->get_start_pos(), 
          world_map->get_goal_pos(),
          10.0,
          mc_max_vertices, 
          mc_prog_interval,
          ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::LINEAR_SEARCH_KNN,
          ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
          mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,fadprm_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    
    
    if(mc_flags & RK_HIDIM_PLANNER_DO_ALT) {
      std::cout << "Running FADPRM with dvp-adj-list-bf2..." << std::endl;
      timing_output << "FADPRM, dvp-adj-list-bf2" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::fadprm_path_planner< SpaceType, ReporterType > 
          fadprm_plan(
            world_map, 
            world_map->get_start_pos(), 
            world_map->get_goal_pos(),
            10.0,
            mc_max_vertices, 
            mc_prog_interval,
            ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF2_TREE_KNN,
            ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
            mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,fadprm_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running FADPRM with dvp-adj-list-bf4..." << std::endl;
      timing_output << "FADPRM, dvp-adj-list-bf4" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::fadprm_path_planner< SpaceType, ReporterType > 
          fadprm_plan(
            world_map, 
            world_map->get_start_pos(), 
            world_map->get_goal_pos(),
            10.0,
            mc_max_vertices, 
            mc_prog_interval,
            ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF4_TREE_KNN,
            ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
            mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,fadprm_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      if(mc_flags & RK_HIDIM_PLANNER_DO_COB) {
        std::cout << "Running FADPRM with dvp-adj-list-cob2..." << std::endl;
        timing_output << "FADPRM, dvp-adj-list-cob2" << std::endl;
        {
          std::stringstream ss, ss2;
          
          ReaK::pp::fadprm_path_planner< SpaceType, ReporterType > 
            fadprm_plan(
              world_map, 
              world_map->get_start_pos(), 
              world_map->get_goal_pos(),
              10.0,
              mc_max_vertices, 
              mc_prog_interval,
              ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB2_TREE_KNN,
              ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
              mc_results);
          
          run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,fadprm_plan,ss,ss2,timing_output);
        };
        std::cout << "Done!" << std::endl;
        
        
        std::cout << "Running FADPRM with dvp-adj-list-cob4..." << std::endl;
        timing_output << "FADPRM, dvp-adj-list-cob4" << std::endl;
        {
          std::stringstream ss, ss2;
          
          ReaK::pp::fadprm_path_planner< SpaceType, ReporterType > 
            fadprm_plan(
              world_map, 
              world_map->get_start_pos(), 
              world_map->get_goal_pos(),
              10.0,
              mc_max_vertices, 
              mc_prog_interval,
              ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB4_TREE_KNN,
              ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
              mc_results);
          
          run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,fadprm_plan,ss,ss2,timing_output);
        };
        std::cout << "Done!" << std::endl;
      };
    };
  };
  
#endif
  
  
  
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
  
  if(mc_flags & RK_HIDIM_PLANNER_DO_SBASTAR) {
    
    /**********************************************************************************
    * 
    * 
    *    Sampling-based A-Star Path-Planners
    * 
    * 
    * *******************************************************************************/
    
    std::cout << "Running SBA* with adj-list, dvp-bf2..." << std::endl;
    timing_output << "SBA*, adj-list, dvp-bf2" << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::sbastar_path_planner< SpaceType, ReporterType > 
        sbastar_plan(world_map, 
                      world_map->get_start_pos(), 
                      world_map->get_goal_pos(),
                      mc_max_vertices, 
                      mc_prog_interval,
                      ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF2_TREE_KNN,
                      ReaK::pp::LAZY_COLLISION_CHECKING | ReaK::pp::PLAN_WITH_ANYTIME_HEURISTIC | ( sba_use_voronoi_pull ? ReaK::pp::PLAN_WITH_VORONOI_PULL : ReaK::pp::NOMINAL_PLANNER_ONLY ),
                      ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                      mc_results);
      
      sbastar_plan.set_initial_key_threshold(sba_potential_cutoff);
      sbastar_plan.set_initial_density_threshold(sba_density_cutoff);
      sbastar_plan.set_initial_relaxation(sba_relaxation);
      sbastar_plan.set_initial_SA_temperature(sba_sa_temperature);
      sbastar_plan.set_sampling_radius( world_map->get_max_edge_length() );
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,sbastar_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    std::cout << "Running SBA* with adj-list, dvp-bf4..." << std::endl;
    timing_output << "SBA*, adj-list, dvp-bf4" << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::sbastar_path_planner< SpaceType, ReporterType > 
        sbastar_plan(world_map, 
                      world_map->get_start_pos(), 
                      world_map->get_goal_pos(),
                      mc_max_vertices, 
                      mc_prog_interval,
                      ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF4_TREE_KNN,
                      ReaK::pp::LAZY_COLLISION_CHECKING | ReaK::pp::PLAN_WITH_ANYTIME_HEURISTIC | ( sba_use_voronoi_pull ? ReaK::pp::PLAN_WITH_VORONOI_PULL : ReaK::pp::NOMINAL_PLANNER_ONLY ),
                      ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                      mc_results);
      
      sbastar_plan.set_initial_key_threshold(sba_potential_cutoff);
      sbastar_plan.set_initial_density_threshold(sba_density_cutoff);
      sbastar_plan.set_initial_relaxation(sba_relaxation);
      sbastar_plan.set_initial_SA_temperature(sba_sa_temperature);
      sbastar_plan.set_sampling_radius( world_map->get_max_edge_length() );
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,sbastar_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    if(mc_flags & RK_HIDIM_PLANNER_DO_COB) {
      std::cout << "Running SBA* with adj-list, dvp-cob2..." << std::endl;
      timing_output << "SBA*, adj-list, dvp-cob2" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::sbastar_path_planner< SpaceType, ReporterType > 
          sbastar_plan(world_map, 
                        world_map->get_start_pos(), 
                        world_map->get_goal_pos(),
                        mc_max_vertices, 
                        mc_prog_interval,
                        ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB2_TREE_KNN,
                        ReaK::pp::LAZY_COLLISION_CHECKING | ReaK::pp::PLAN_WITH_ANYTIME_HEURISTIC | ( sba_use_voronoi_pull ? ReaK::pp::PLAN_WITH_VORONOI_PULL : ReaK::pp::NOMINAL_PLANNER_ONLY ),
                        ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                        mc_results);
        
        sbastar_plan.set_initial_key_threshold(sba_potential_cutoff);
        sbastar_plan.set_initial_density_threshold(sba_density_cutoff);
        sbastar_plan.set_initial_relaxation(sba_relaxation);
        sbastar_plan.set_initial_SA_temperature(sba_sa_temperature);
        sbastar_plan.set_sampling_radius( world_map->get_max_edge_length() );
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,sbastar_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running SBA* with adj-list, dvp-cob4..." << std::endl;
      timing_output << "SBA*, adj-list, dvp-cob4" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::sbastar_path_planner< SpaceType, ReporterType > 
          sbastar_plan(world_map, 
                        world_map->get_start_pos(), 
                        world_map->get_goal_pos(),
                        mc_max_vertices, 
                        mc_prog_interval,
                        ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB4_TREE_KNN,
                        ReaK::pp::LAZY_COLLISION_CHECKING | ReaK::pp::PLAN_WITH_ANYTIME_HEURISTIC | ( sba_use_voronoi_pull ? ReaK::pp::PLAN_WITH_VORONOI_PULL : ReaK::pp::NOMINAL_PLANNER_ONLY ),
                        ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                        mc_results);
        
        sbastar_plan.set_initial_key_threshold(sba_potential_cutoff);
        sbastar_plan.set_initial_density_threshold(sba_density_cutoff);
        sbastar_plan.set_initial_relaxation(sba_relaxation);
        sbastar_plan.set_initial_SA_temperature(sba_sa_temperature);
        sbastar_plan.set_sampling_radius( world_map->get_max_edge_length() );
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,sbastar_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
    };
    
    
    std::cout << "Running SBA* with adj-list, linear-search..." << std::endl;
    timing_output << "SBA*, adj-list, linear-search" << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::sbastar_path_planner< SpaceType, ReporterType > 
        sbastar_plan(world_map, 
                      world_map->get_start_pos(), 
                      world_map->get_goal_pos(),
                      mc_max_vertices, 
                      mc_prog_interval,
                      ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::LINEAR_SEARCH_KNN,
                      ReaK::pp::LAZY_COLLISION_CHECKING | ReaK::pp::PLAN_WITH_ANYTIME_HEURISTIC | ( sba_use_voronoi_pull ? ReaK::pp::PLAN_WITH_VORONOI_PULL : ReaK::pp::NOMINAL_PLANNER_ONLY ),
                      ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                      mc_results);
      
      sbastar_plan.set_initial_key_threshold(sba_potential_cutoff);
      sbastar_plan.set_initial_density_threshold(sba_density_cutoff);
      sbastar_plan.set_initial_relaxation(sba_relaxation);
      sbastar_plan.set_initial_SA_temperature(sba_sa_temperature);
      sbastar_plan.set_sampling_radius( world_map->get_max_edge_length() );
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,sbastar_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    
    if(mc_flags & RK_HIDIM_PLANNER_DO_ALT) {
      std::cout << "Running SBA* with dvp-adj-list-bf2..." << std::endl;
      timing_output << "SBA*, dvp-adj-list-bf2" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::sbastar_path_planner< SpaceType, ReporterType > 
          sbastar_plan(world_map, 
                        world_map->get_start_pos(), 
                        world_map->get_goal_pos(),
                        mc_max_vertices, 
                        mc_prog_interval,
                        ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF2_TREE_KNN,
                        ReaK::pp::LAZY_COLLISION_CHECKING | ReaK::pp::PLAN_WITH_ANYTIME_HEURISTIC | ( sba_use_voronoi_pull ? ReaK::pp::PLAN_WITH_VORONOI_PULL : ReaK::pp::NOMINAL_PLANNER_ONLY ),
                        ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                        mc_results);
        
        sbastar_plan.set_initial_key_threshold(sba_potential_cutoff);
        sbastar_plan.set_initial_density_threshold(sba_density_cutoff);
        sbastar_plan.set_initial_relaxation(sba_relaxation);
        sbastar_plan.set_initial_SA_temperature(sba_sa_temperature);
        sbastar_plan.set_sampling_radius( world_map->get_max_edge_length() );
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,sbastar_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running SBA* with dvp-adj-list-bf4..." << std::endl;
      timing_output << "SBA*, dvp-adj-list-bf4" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::sbastar_path_planner< SpaceType, ReporterType > 
          sbastar_plan(world_map, 
                        world_map->get_start_pos(), 
                        world_map->get_goal_pos(),
                        mc_max_vertices, 
                        mc_prog_interval,
                        ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF4_TREE_KNN,
                        ReaK::pp::LAZY_COLLISION_CHECKING | ReaK::pp::PLAN_WITH_ANYTIME_HEURISTIC | ( sba_use_voronoi_pull ? ReaK::pp::PLAN_WITH_VORONOI_PULL : ReaK::pp::NOMINAL_PLANNER_ONLY ),
                        ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                        mc_results);
        
        sbastar_plan.set_initial_key_threshold(sba_potential_cutoff);
        sbastar_plan.set_initial_density_threshold(sba_density_cutoff);
        sbastar_plan.set_initial_relaxation(sba_relaxation);
        sbastar_plan.set_initial_SA_temperature(sba_sa_temperature);
        sbastar_plan.set_sampling_radius( world_map->get_max_edge_length() );
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,sbastar_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      if(mc_flags & RK_HIDIM_PLANNER_DO_COB) {
        std::cout << "Running SBA* with dvp-adj-list-cob2..." << std::endl;
        timing_output << "SBA*, dvp-adj-list-cob2" << std::endl;
        {
          std::stringstream ss, ss2;
          
          ReaK::pp::sbastar_path_planner< SpaceType, ReporterType > 
            sbastar_plan(world_map, 
                          world_map->get_start_pos(), 
                          world_map->get_goal_pos(),
                          mc_max_vertices, 
                          mc_prog_interval,
                          ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB2_TREE_KNN,
                          ReaK::pp::LAZY_COLLISION_CHECKING | ReaK::pp::PLAN_WITH_ANYTIME_HEURISTIC | ( sba_use_voronoi_pull ? ReaK::pp::PLAN_WITH_VORONOI_PULL : ReaK::pp::NOMINAL_PLANNER_ONLY ),
                          ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                          mc_results);
          
          sbastar_plan.set_initial_key_threshold(sba_potential_cutoff);
          sbastar_plan.set_initial_density_threshold(sba_density_cutoff);
          sbastar_plan.set_initial_relaxation(sba_relaxation);
          sbastar_plan.set_initial_SA_temperature(sba_sa_temperature);
          sbastar_plan.set_sampling_radius( world_map->get_max_edge_length() );
          
          run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,sbastar_plan,ss,ss2,timing_output);
        };
        std::cout << "Done!" << std::endl;
        
        
        std::cout << "Running SBA* with dvp-adj-list-cob4..." << std::endl;
        timing_output << "SBA*, dvp-adj-list-cob4" << std::endl;
        {
          std::stringstream ss, ss2;
          
          ReaK::pp::sbastar_path_planner< SpaceType, ReporterType > 
            sbastar_plan(world_map, 
                          world_map->get_start_pos(), 
                          world_map->get_goal_pos(),
                          mc_max_vertices, 
                          mc_prog_interval,
                          ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB4_TREE_KNN,
                          ReaK::pp::LAZY_COLLISION_CHECKING | ReaK::pp::PLAN_WITH_ANYTIME_HEURISTIC | ( sba_use_voronoi_pull ? ReaK::pp::PLAN_WITH_VORONOI_PULL : ReaK::pp::NOMINAL_PLANNER_ONLY ),
                          ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                          mc_results);
          
          sbastar_plan.set_initial_key_threshold(sba_potential_cutoff);
          sbastar_plan.set_initial_density_threshold(sba_density_cutoff);
          sbastar_plan.set_initial_relaxation(sba_relaxation);
          sbastar_plan.set_initial_SA_temperature(sba_sa_temperature);
          sbastar_plan.set_sampling_radius( world_map->get_max_edge_length() );
          
          run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,sbastar_plan,ss,ss2,timing_output);
        };
        std::cout << "Done!" << std::endl;
      };
    };
  };
    
#endif
  
  
  
  
  
#ifdef RK_ENABLE_TEST_RRTSTAR_PLANNER
  
  if(mc_flags & RK_HIDIM_PLANNER_DO_RRTSTAR) {
    
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
      std::stringstream ss, ss2;
      
      ReaK::pp::rrtstar_path_planner< SpaceType, ReporterType > 
        rrtstar_plan(world_map, 
                     world_map->get_start_pos(), 
                     world_map->get_goal_pos(),
                     mc_max_vertices, 
                     mc_prog_interval,
                     ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF2_TREE_KNN,
                     ReaK::pp::UNIDIRECTIONAL_PLANNING | ReaK::pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG,
                     ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                     mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrtstar_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    std::cout << "Running RRT* with Uni-dir, adj-list, dvp-bf4..." << std::endl;
    timing_output << "RRT*, Uni-dir, adj-list, dvp-bf4" << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::rrtstar_path_planner< SpaceType, ReporterType > 
        rrtstar_plan(world_map, 
                     world_map->get_start_pos(), 
                     world_map->get_goal_pos(),
                     mc_max_vertices, 
                     mc_prog_interval,
                     ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF4_TREE_KNN,
                     ReaK::pp::UNIDIRECTIONAL_PLANNING | ReaK::pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG,
                     ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                     mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrtstar_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    if(mc_flags & RK_HIDIM_PLANNER_DO_COB) {
      std::cout << "Running RRT* with Uni-dir, adj-list, dvp-cob2..." << std::endl;
      timing_output << "RRT*, Uni-dir, adj-list, dvp-cob2" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrtstar_path_planner< SpaceType, ReporterType > 
          rrtstar_plan(world_map, 
                       world_map->get_start_pos(), 
                       world_map->get_goal_pos(),
                       mc_max_vertices, 
                       mc_prog_interval,
                       ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB2_TREE_KNN,
                       ReaK::pp::UNIDIRECTIONAL_PLANNING | ReaK::pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG,
                       ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                       mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrtstar_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT* with Uni-dir, adj-list, dvp-cob4..." << std::endl;
      timing_output << "RRT*, Uni-dir, adj-list, dvp-cob4" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrtstar_path_planner< SpaceType, ReporterType > 
          rrtstar_plan(world_map, 
                       world_map->get_start_pos(), 
                       world_map->get_goal_pos(),
                       mc_max_vertices, 
                       mc_prog_interval,
                       ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB4_TREE_KNN,
                       ReaK::pp::UNIDIRECTIONAL_PLANNING | ReaK::pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG,
                       ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                       mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrtstar_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
    };
    
    
    std::cout << "Running RRT* with Uni-dir, adj-list, linear-search..." << std::endl;
    timing_output << "RRT*, Uni-dir, adj-list, linear-search" << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::rrtstar_path_planner< SpaceType, ReporterType > 
        rrtstar_plan(world_map, 
                     world_map->get_start_pos(), 
                     world_map->get_goal_pos(),
                     mc_max_vertices, 
                     mc_prog_interval,
                     ReaK::pp::ADJ_LIST_MOTION_GRAPH | ReaK::pp::LINEAR_SEARCH_KNN,
                     ReaK::pp::UNIDIRECTIONAL_PLANNING | ReaK::pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG,
                     ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                     mc_results);
      
      run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrtstar_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
    
    
    if(mc_flags & RK_HIDIM_PLANNER_DO_ALT) {
      std::cout << "Running RRT* with Uni-dir, dvp-adj-list-bf2..." << std::endl;
      timing_output << "RRT*, Uni-dir, dvp-adj-list-bf2" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrtstar_path_planner< SpaceType, ReporterType > 
          rrtstar_plan(world_map, 
                       world_map->get_start_pos(), 
                       world_map->get_goal_pos(),
                       mc_max_vertices, 
                       mc_prog_interval,
                       ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF2_TREE_KNN,
                       ReaK::pp::UNIDIRECTIONAL_PLANNING | ReaK::pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG,
                       ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                       mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrtstar_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT* with Uni-dir, dvp-adj-list-bf4..." << std::endl;
      timing_output << "RRT*, Uni-dir, dvp-adj-list-bf4" << std::endl;
      {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrtstar_path_planner< SpaceType, ReporterType > 
          rrtstar_plan(world_map, 
                       world_map->get_start_pos(), 
                       world_map->get_goal_pos(),
                       mc_max_vertices, 
                       mc_prog_interval,
                       ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_BF4_TREE_KNN,
                       ReaK::pp::UNIDIRECTIONAL_PLANNING | ReaK::pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG,
                       ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                       mc_results);
        
        run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrtstar_plan,ss,ss2,timing_output);
      };
      std::cout << "Done!" << std::endl;
      
      
      if(mc_flags & RK_HIDIM_PLANNER_DO_COB) {
        std::cout << "Running RRT* with Uni-dir, dvp-adj-list-cob2..." << std::endl;
        timing_output << "RRT*, Uni-dir, dvp-adj-list-cob2" << std::endl;
        {
          std::stringstream ss, ss2;
          
          ReaK::pp::rrtstar_path_planner< SpaceType, ReporterType > 
            rrtstar_plan(world_map, 
                         world_map->get_start_pos(), 
                         world_map->get_goal_pos(),
                         mc_max_vertices, 
                         mc_prog_interval,
                         ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB2_TREE_KNN,
                         ReaK::pp::UNIDIRECTIONAL_PLANNING | ReaK::pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG,
                         ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                         mc_results);
          
          run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrtstar_plan,ss,ss2,timing_output);
        };
        std::cout << "Done!" << std::endl;
        
        
        std::cout << "Running RRT* with Uni-dir, dvp-adj-list-cob4..." << std::endl;
        timing_output << "RRT*, Uni-dir, dvp-adj-list-cob4" << std::endl;
        {
          std::stringstream ss, ss2;
          
          ReaK::pp::rrtstar_path_planner< SpaceType, ReporterType > 
            rrtstar_plan(world_map, 
                         world_map->get_start_pos(), 
                         world_map->get_goal_pos(),
                         mc_max_vertices, 
                         mc_prog_interval,
                         ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH | ReaK::pp::DVP_COB4_TREE_KNN,
                         ReaK::pp::UNIDIRECTIONAL_PLANNING | ReaK::pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG,
                         ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                         mc_results);
          
          run_monte_carlo_tests(mc_run_count,mc_max_vertices_100,rrtstar_plan,ss,ss2,timing_output);
        };
        std::cout << "Done!" << std::endl;
      };
    };
  };
    
#endif
  
  
};










int main(int argc, char** argv) {
  
  
  po::options_description generic_options("Generic options");
  generic_options.add_options()
    ("help,h", "produce this help message.")
  ;
  
  po::options_description io_options("I/O options");
  io_options.add_options()
    ("output-path,o", po::value< std::string >()->default_value("pp_results"), "specify the output path (default is pp_results)")
  ;
  
  po::options_description mc_options("Monte-Carlo options");
  mc_options.add_options()
    ("mc-runs", po::value< std::size_t >()->default_value(500), "number of monte-carlo runs to average out")
    ("mc-vertices", po::value< std::size_t >()->default_value(20000), "maximum number of vertices during monte-carlo runs")
    ("mc-prog-interval", po::value< std::size_t >()->default_value(100), "number of vertices between progress reports during monte-carlo runs")
    ("mc-results", po::value< std::size_t >()->default_value(5), "maximum number of result-paths during monte-carlo runs")
    ("mc-dvp-alt", "do monte-carlo runs with the DVP-adjacency-list-tree layout (default is not)")
    ("mc-cob-tree", "do monte-carlo runs with cache-oblivious b-trees (default is not)")
  ;
  
  po::options_description planner_select_options("Planner selection options");
  planner_select_options.add_options()
#ifdef RK_ENABLE_TEST_URRT_PLANNER
    ("rrt", "specify that the uni-directional RRT algorithm should be run")
#endif
#ifdef RK_ENABLE_TEST_BRRT_PLANNER
    ("bi-rrt", "specify that the bi-directional RRT algorithm should be run")
#endif
#ifdef RK_ENABLE_TEST_RRTSTAR_PLANNER
    ("rrt-star", "specify that the RRT* algorithm should be run")
#endif
#ifdef RK_ENABLE_TEST_PRM_PLANNER
    ("prm", "specify that the PRM algorithm should be run")
#endif
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
    ("fadprm", "specify that the FADPRM algorithm should be run")
#endif
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
    ("sba-star", "specify that the SBA* algorithm should be run")
    ("sba-potential-cutoff", po::value< double >()->default_value(0.0), "specify the potential cutoff for the SBA* algorithm")
    ("sba-density-cutoff", po::value< double >()->default_value(0.5), "specify the density cutoff for the SBA* algorithm")
    ("sba-relaxation", po::value< double >()->default_value(0.0), "specify the initial relaxation factor for the Anytime SBA* algorithm")
    ("sba-with-voronoi-pull", "specify whether to use a Voronoi pull or not as a method to add an exploratory bias to the search")
    ("sba-sa-temperature", po::value< double >()->default_value(-1.0), "specify the initial Simulated Annealing temperature for the SBA*-RRT* algorithms")
#endif
    ("all-planners,a", "specify that all supported planners should be run (default if no particular planner is specified)")
  ;
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options).add(io_options).add(mc_options).add(planner_select_options);
  
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
  
  bool run_all_planners = false;
  if(vm.count("all-planners") || (vm.count("rrt") + vm.count("bi-rrt") + vm.count("rrt-star") + vm.count("prm") + vm.count("fadprm") + vm.count("sba-star") == 0)) 
    run_all_planners = true;
  
  std::size_t mc_run_count        = vm["mc-runs"].as<std::size_t>();
  std::size_t mc_max_vertices     = vm["mc-vertices"].as<std::size_t>();
  std::size_t mc_prog_interval    = vm["mc-prog-interval"].as<std::size_t>();
  std::size_t mc_results          = vm["mc-results"].as<std::size_t>();
  
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
  sba_potential_cutoff = vm["sba-potential-cutoff"].as<double>();
  sba_density_cutoff   = vm["sba-density-cutoff"].as<double>();
  sba_relaxation       = vm["sba-relaxation"].as<double>();
  sba_sa_temperature   = vm["sba-sa-temperature"].as<double>();
  sba_use_voronoi_pull = vm.count("sba-with-voronoi-pull");
#endif
  
  std::size_t mc_flags = 0;
#ifdef RK_ENABLE_TEST_URRT_PLANNER
  if(run_all_planners || vm.count("rrt"))
    mc_flags |= RK_HIDIM_PLANNER_DO_RRT;
#endif
#ifdef RK_ENABLE_TEST_BRRT_PLANNER
  if(run_all_planners || vm.count("bi-rrt"))
    mc_flags |= RK_HIDIM_PLANNER_DO_BIRRT;
#endif
#ifdef RK_ENABLE_TEST_RRTSTAR_PLANNER
  if(run_all_planners || vm.count("rrt-star"))
    mc_flags |= RK_HIDIM_PLANNER_DO_RRTSTAR;
#endif
#ifdef RK_ENABLE_TEST_PRM_PLANNER
  if(run_all_planners || vm.count("prm"))
    mc_flags |= RK_HIDIM_PLANNER_DO_PRM;
#endif
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
  if(run_all_planners || vm.count("fadprm"))
    mc_flags |= RK_HIDIM_PLANNER_DO_FADPRM;
#endif
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
  if(run_all_planners || vm.count("sba-star"))
    mc_flags |= RK_HIDIM_PLANNER_DO_SBASTAR;
#endif
  if(vm.count("mc-cob-tree"))
    mc_flags |= RK_HIDIM_PLANNER_DO_COB;
  if(vm.count("mc-dvp-alt"))
    mc_flags |= RK_HIDIM_PLANNER_DO_ALT;
  
  
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double, 6> > > World6DType;
  ReaK::shared_ptr< World6DType > world_6D =
    ReaK::shared_ptr< World6DType >(
      new World6DType(
        "world_6D_no_obstacles",
        ReaK::pp::hyperbox_topology< ReaK::vect<double, 6> >("world_6D",
                                                             ReaK::vect<double,6>(0.0,0.0,0.0,0.0,0.0,0.0),
                                                             ReaK::vect<double,6>(1.0,1.0,1.0,1.0,1.0,1.0)),
        0.1));
  world_6D->set_start_pos(ReaK::vect<double,6>(0.05,0.05,0.05,0.05,0.05,0.05));
  world_6D->set_goal_pos(ReaK::vect<double,6>(0.95,0.95,0.95,0.95,0.95,0.95));
  
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double, 12> > > World12DType;
  ReaK::shared_ptr< World12DType > world_12D =
    ReaK::shared_ptr< World12DType >(
      new World12DType(
        "world_12D_no_obstacles",
        ReaK::pp::hyperbox_topology< ReaK::vect<double, 12> >("world_12D",
                                                              ReaK::vect<double,12>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                              ReaK::vect<double,12>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
        0.5));
  world_12D->set_start_pos(ReaK::vect<double,12>(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
  world_12D->set_goal_pos( ReaK::vect<double,12>(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95));
  
  typedef ReaK::pp::no_obstacle_space< typename ReaK::pp::se3_1st_order_rl_topology<double>::type > WorldRLSE3Type;
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
  typedef ReaK::arithmetic_tuple< ReaK::arithmetic_tuple< ReaK::vect<double,3>,    ReaK::vect<double,3> >,
                                  ReaK::arithmetic_tuple< ReaK::unit_quat<double>, ReaK::vect<double,3> > > rlse3_point;
  typedef ReaK::arithmetic_tuple< ReaK::vect<double,3>,    ReaK::vect<double,3> > rlp3_point;
  typedef ReaK::arithmetic_tuple< ReaK::unit_quat<double>, ReaK::vect<double,3> > rlq3_point;
  world_RLSE3->set_start_pos(
    rlse3_point(rlp3_point(ReaK::vect<double,3>(-9.5,-9.5,-9.5),ReaK::vect<double,3>()),rlq3_point())
  );
  world_RLSE3->set_goal_pos(
    rlse3_point(rlp3_point(ReaK::vect<double,3>( 9.5, 9.5, 9.5),ReaK::vect<double,3>()),rlq3_point())
  );
  
  
  {
    std::ofstream timing_output(output_path_name + "/" + world_6D->getName() + "_times.txt");
    
    test_planners_on_space(world_6D, timing_output, mc_run_count, mc_max_vertices, mc_prog_interval, mc_results, mc_flags);
  };
  
  {
    std::ofstream timing_output(output_path_name + "/" + world_12D->getName() + "_times.txt");
    
    test_planners_on_space(world_12D, timing_output, mc_run_count, mc_max_vertices, mc_prog_interval, mc_results, mc_flags);
  };
  
  {
    std::ofstream timing_output(output_path_name + "/" + world_RLSE3->getName() + "_times.txt");
    
    test_planners_on_space(world_RLSE3, timing_output, mc_run_count, mc_max_vertices, mc_prog_interval, mc_results, mc_flags);
  };
  
  
  
  return 0;
};












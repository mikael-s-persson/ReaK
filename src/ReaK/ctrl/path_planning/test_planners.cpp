
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


int main(int argc, char** argv) {
  
  po::options_description generic_options("Generic options");
  generic_options.add_options()
    ("help,h", "produce this help message.")
  ;
  
  po::options_description io_options("I/O options");
  io_options.add_options()
    ("input-file,i", po::value< std::string >(), "specify the input file (default is stdin)")
    ("output-path,o", po::value< std::string >()->default_value("pp_results"), "specify the output path (default is pp_results)")
  ;
  
  po::options_description mc_options("Monte-Carlo options");
  mc_options.add_options()
    ("monte-carlo,m", "specify that monte-carlo runs should be performed (default is not)")
    ("mc-runs", po::value< std::size_t >()->default_value(100), "number of monte-carlo runs to average out (default is 100)")
    ("mc-vertices", po::value< std::size_t >()->default_value(5000), "maximum number of vertices during monte-carlo runs (default is 5000)")
    ("mc-prog-interval", po::value< std::size_t >()->default_value(100), "number of vertices between progress reports during monte-carlo runs (default is 5000)")
    ("mc-results", po::value< std::size_t >()->default_value(50), "maximum number of result-paths during monte-carlo runs (default is 50)")
    ("mc-dvp-alt", "do monte-carlo runs with the DVP-adjacency-list-tree layout (default is not)")
    ("mc-cob-tree", "do monte-carlo runs with cache-oblivious b-trees (default is not)")
  ;
  
  po::options_description single_options("Single-run options");
  single_options.add_options()
    ("single-run,s", "specify that single runs should be performed (default is not)")
    ("max-vertices", po::value< std::size_t >()->default_value(5000), "maximum number of vertices during single runs (default is 5000)")
    ("max-results", po::value< std::size_t >()->default_value(50), "maximum number of result-paths during single runs (default is 50)")
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
    ("sba-potential-cutoff", po::value< double >()->default_value(0.9), "specify the potential cutoff for the SBA* algorithm")
    ("sba-density-cutoff", po::value< double >()->default_value(0.9), "specify the density cutoff for the SBA* algorithm")
#endif
    ("all-planners,a", "specify that all supported planners should be run (default if no particular planner is specified)")
  ;
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options).add(io_options).add(mc_options).add(single_options).add(planner_select_options);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);
  
  if(vm.count("help") || (vm.count("input-file") == 0) || (vm.count("monte-carlo") + vm.count("single-run") < 1)) {
    std::cout << cmdline_options << std::endl;
    return 1;
  };
  
  std::string world_file_name = vm["input-file"].as<std::string>();
  std::string world_file_name_only(std::find(world_file_name.rbegin(),world_file_name.rend(),'/').base(), std::find(world_file_name.rbegin(),world_file_name.rend(),'.').base()-1);
  
  std::string output_path_name = vm["output-path"].as<std::string>();
  while(output_path_name[output_path_name.length()-1] == '/') 
    output_path_name.erase(output_path_name.length()-1, 1);
  
  bool run_all_planners = false;
  if(vm.count("all-planners") || (vm.count("rrt") + vm.count("bi-rrt") + vm.count("rrt-star") + vm.count("prm") + vm.count("fadprm") + vm.count("sba-star") == 0)) 
    run_all_planners = true;
  
  ReaK::shared_ptr< ReaK::pp::ptrobot2D_test_world > world_map =
    ReaK::shared_ptr< ReaK::pp::ptrobot2D_test_world >(new ReaK::pp::ptrobot2D_test_world(world_file_name, 20, 1.0));
  
  if(vm.count("monte-carlo")) {
    
    std::size_t mc_run_count        = vm["mc-runs"].as<std::size_t>();
    std::size_t mc_max_vertices     = vm["mc-vertices"].as<std::size_t>();
    std::size_t mc_prog_interval    = vm["mc-prog-interval"].as<std::size_t>();
    std::size_t mc_max_vertices_100 = mc_max_vertices / mc_prog_interval;
    std::size_t mc_results          = vm["mc-results"].as<std::size_t>();
    
    std::ofstream timing_output(output_path_name + "/" + world_file_name_only + "_times.txt");
    
    typedef ReaK::pp::timing_sbmp_report< ReaK::pp::least_cost_sbmp_report<> > ReporterType;
    
#ifdef RK_ENABLE_TEST_URRT_PLANNER
    
    if(run_all_planners || vm.count("rrt")) {
      
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
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10, 0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::UNIDIRECTIONAL_RRT,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_BF2_TREE_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1])
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT with Uni-dir, adj-list, dvp-bf4..." << std::endl;
      timing_output << "RRT, Uni-dir, adj-list, dvp-bf4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::UNIDIRECTIONAL_RRT,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_BF4_TREE_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      if(vm.count("mc-cob-tree")) {
      std::cout << "Running RRT with Uni-dir, adj-list, dvp-cob2..." << std::endl;
      timing_output << "RRT, Uni-dir, adj-list, dvp-cob2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::UNIDIRECTIONAL_RRT,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_COB2_TREE_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      std::cout << "Running RRT with Uni-dir, adj-list, dvp-cob4..." << std::endl;
      timing_output << "RRT, Uni-dir, adj-list, dvp-cob4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::UNIDIRECTIONAL_RRT,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_COB4_TREE_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      };
      
      
      std::cout << "Running RRT with Uni-dir, adj-list, linear-search..." << std::endl;
      timing_output << "RRT, Uni-dir, adj-list, linear-search" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::UNIDIRECTIONAL_RRT,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::LINEAR_SEARCH_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      if(vm.count("mc-dvp-alt")) {
      std::cout << "Running RRT with Uni-dir, dvp-adj-list-bf2..." << std::endl;
      timing_output << "RRT, Uni-dir, dvp-adj-list-bf2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::UNIDIRECTIONAL_RRT,
                  ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_ALT_BF2_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT with Uni-dir, dvp-adj-list-bf4..." << std::endl;
      timing_output << "RRT, Uni-dir, dvp-adj-list-bf4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::UNIDIRECTIONAL_RRT,
                  ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_ALT_BF4_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      if(vm.count("mc-cob-tree")) {
      std::cout << "Running RRT with Uni-dir, dvp-adj-list-cob2..." << std::endl;
      timing_output << "RRT, Uni-dir, dvp-adj-list-cob2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::UNIDIRECTIONAL_RRT,
                  ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_ALT_COB2_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT with Uni-dir, dvp-adj-list-cob4..." << std::endl;
      timing_output << "RRT, Uni-dir, dvp-adj-list-cob4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::UNIDIRECTIONAL_RRT,
                  ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_ALT_COB4_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      };
      };
      
    };
    
#endif
    
    
#ifdef RK_ENABLE_TEST_BRRT_PLANNER
    
    
    if(run_all_planners || vm.count("bi-rrt")) {
      
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
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::BIDIRECTIONAL_RRT,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_BF2_TREE_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
          ss >> v_count >> t_val;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT with Bi-dir, adj-list, dvp-bf4..." << std::endl;
      timing_output << "RRT, Bi-dir, adj-list, dvp-bf4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::BIDIRECTIONAL_RRT,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_BF4_TREE_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
          ss >> v_count >> t_val;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      if(vm.count("mc-cob-tree")) {
      std::cout << "Running RRT with Bi-dir, adj-list, dvp-cob2..." << std::endl;
      timing_output << "RRT, Bi-dir, adj-list, dvp-cob2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::BIDIRECTIONAL_RRT,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_COB2_TREE_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT with Bi-dir, adj-list, dvp-cob4..." << std::endl;
      timing_output << "RRT, Bi-dir, adj-list, dvp-cob4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::BIDIRECTIONAL_RRT,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_COB4_TREE_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
          ss >> v_count >> t_val;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      };
      
      
      std::cout << "Running RRT with Bi-dir, adj-list, linear-search..." << std::endl;
      timing_output << "RRT, Bi-dir, adj-list, linear-search" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::BIDIRECTIONAL_RRT,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::LINEAR_SEARCH_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
          ss >> v_count >> t_val;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      
      if(vm.count("mc-dvp-alt")) {
      std::cout << "Running RRT with Bi-dir, dvp-adj-list-bf2..." << std::endl;
      timing_output << "RRT, Bi-dir, dvp-adj-list-bf2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::BIDIRECTIONAL_RRT,
                  ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_ALT_BF2_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
          ss >> v_count >> t_val;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT with Bi-dir, dvp-adj-list-bf4..." << std::endl;
      timing_output << "RRT, Bi-dir, dvp-adj-list-bf4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::BIDIRECTIONAL_RRT,
                  ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_ALT_BF4_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
          ss >> v_count >> t_val;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      if(vm.count("mc-cob-tree")) {
      std::cout << "Running RRT with Bi-dir, dvp-adj-list-cob2..." << std::endl;
      timing_output << "RRT, Bi-dir, dvp-adj-list-cob2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::BIDIRECTIONAL_RRT,
                  ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_ALT_COB2_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
          ss >> v_count >> t_val;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT with Bi-dir, dvp-adj-list-cob4..." << std::endl;
      timing_output << "RRT, Bi-dir, dvp-adj-list-cob4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrt_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::BIDIRECTIONAL_RRT,
                  ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_ALT_COB4_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        rrt_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
          ss >> v_count >> t_val;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      };
      };
      
    };
    
  #endif
    
    
    
  #ifdef RK_ENABLE_TEST_PRM_PLANNER
    
    
    if(run_all_planners || vm.count("prm")) {
      
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
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          prm_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_BF2_TREE_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        prm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running PRM with adj-list, dvp-bf4..." << std::endl;
      timing_output << "PRM, adj-list, dvp-bf4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          prm_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_BF4_TREE_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        prm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      if(vm.count("mc-cob-tree")) {
      std::cout << "Running PRM with adj-list, dvp-cob2..." << std::endl;
      timing_output << "PRM, adj-list, dvp-cob2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          prm_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_COB2_TREE_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        prm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running PRM with adj-list, dvp-cob4..." << std::endl;
      timing_output << "PRM, adj-list, dvp-cob4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          prm_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_COB4_TREE_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        prm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      };
      
      
      std::cout << "Running PRM with adj-list, linear-search..." << std::endl;
      timing_output << "PRM, adj-list, linear-search" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          prm_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::LINEAR_SEARCH_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        prm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      
      
      if(vm.count("mc-dvp-alt")) {
      std::cout << "Running PRM with dvp-adj-list-bf2..." << std::endl;
      timing_output << "PRM, dvp-adj-list-bf2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          prm_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_ALT_BF2_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        prm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running PRM with dvp-adj-list-bf4..." << std::endl;
      timing_output << "PRM, dvp-adj-list-bf4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          prm_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_ALT_BF4_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        prm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      if(vm.count("mc-cob-tree")) {
      std::cout << "Running PRM with dvp-adj-list-cob2..." << std::endl;
      timing_output << "PRM, dvp-adj-list-cob2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          prm_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_ALT_COB2_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        prm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running PRM with dvp-adj-list-cob4..." << std::endl;
      timing_output << "PRM, dvp-adj-list-cob4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          prm_plan(world_map, 
                  world_map->get_start_pos(), 
                  world_map->get_goal_pos(),
                  mc_max_vertices, 
                  mc_prog_interval,
                  ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                  ReaK::pp::DVP_ALT_COB4_KNN,
                  ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                  mc_results);
        
        prm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      };
      };
      
    };
    
#endif
    
    
    
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
    
    if(run_all_planners || vm.count("fadprm")) {
      
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
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::fadprm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          fadprm_plan(
            world_map, 
            world_map->get_start_pos(),
            world_map->get_goal_pos(),
            10.0,
            mc_max_vertices, 
            mc_prog_interval,
            ReaK::pp::ADJ_LIST_MOTION_GRAPH,
            ReaK::pp::DVP_BF2_TREE_KNN,
            ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
            mc_results);
        
        fadprm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      
      std::cout << "Running FADPRM with adj-list, dvp-bf4..." << std::endl;
      timing_output << "FADPRM, adj-list, dvp-bf4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::fadprm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          fadprm_plan(
            world_map, 
            world_map->get_start_pos(), 
            world_map->get_goal_pos(),
            10.0,
            mc_max_vertices, 
            mc_prog_interval,
            ReaK::pp::ADJ_LIST_MOTION_GRAPH,
            ReaK::pp::DVP_BF4_TREE_KNN,
            ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
            mc_results);
        
        fadprm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      if(vm.count("mc-cob-tree")) {
      std::cout << "Running FADPRM with adj-list, dvp-cob2..." << std::endl;
      timing_output << "FADPRM, adj-list, dvp-cob2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::fadprm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          fadprm_plan(
            world_map, 
            world_map->get_start_pos(), 
            world_map->get_goal_pos(),
            10.0,
            mc_max_vertices, 
            mc_prog_interval,
            ReaK::pp::ADJ_LIST_MOTION_GRAPH,
            ReaK::pp::DVP_COB2_TREE_KNN,
            ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
            mc_results);
        
        fadprm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running FADPRM with adj-list, dvp-cob4..." << std::endl;
      timing_output << "FADPRM, adj-list, dvp-cob4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::fadprm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          fadprm_plan(
            world_map, 
            world_map->get_start_pos(), 
            world_map->get_goal_pos(),
            10.0,
            mc_max_vertices, 
            mc_prog_interval,
            ReaK::pp::ADJ_LIST_MOTION_GRAPH,
            ReaK::pp::DVP_COB4_TREE_KNN,
            ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
            mc_results);
        
        fadprm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      };
      
      
      std::cout << "Running FADPRM with adj-list, linear-search..." << std::endl;
      timing_output << "FADPRM, adj-list, linear-search" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::fadprm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          fadprm_plan(
            world_map, 
            world_map->get_start_pos(), 
            world_map->get_goal_pos(),
            10.0,
            mc_max_vertices, 
            mc_prog_interval,
            ReaK::pp::ADJ_LIST_MOTION_GRAPH,
            ReaK::pp::LINEAR_SEARCH_KNN,
            ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
            mc_results);
        
        fadprm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      
      
      if(vm.count("mc-dvp-alt")) {
      std::cout << "Running FADPRM with dvp-adj-list-bf2..." << std::endl;
      timing_output << "FADPRM, dvp-adj-list-bf2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::fadprm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          fadprm_plan(
            world_map, 
            world_map->get_start_pos(), 
            world_map->get_goal_pos(),
            10.0,
            mc_max_vertices, 
            mc_prog_interval,
            ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
            ReaK::pp::DVP_ALT_BF2_KNN,
            ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
            mc_results);
        
        fadprm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running FADPRM with dvp-adj-list-bf4..." << std::endl;
      timing_output << "FADPRM, dvp-adj-list-bf4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::fadprm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          fadprm_plan(
            world_map, 
            world_map->get_start_pos(), 
            world_map->get_goal_pos(),
            10.0,
            mc_max_vertices, 
            mc_prog_interval,
            ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
            ReaK::pp::DVP_ALT_BF4_KNN,
            ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
            mc_results);
        
        fadprm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      if(vm.count("mc-cob-tree")) {
      std::cout << "Running FADPRM with dvp-adj-list-cob2..." << std::endl;
      timing_output << "FADPRM, dvp-adj-list-cob2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::fadprm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          fadprm_plan(
            world_map, 
            world_map->get_start_pos(), 
            world_map->get_goal_pos(),
            10.0,
            mc_max_vertices, 
            mc_prog_interval,
            ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
            ReaK::pp::DVP_ALT_COB2_KNN,
            ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
            mc_results);
        
        fadprm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running FADPRM with dvp-adj-list-cob4..." << std::endl;
      timing_output << "FADPRM, dvp-adj-list-cob4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::fadprm_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          fadprm_plan(
            world_map, 
            world_map->get_start_pos(), 
            world_map->get_goal_pos(),
            10.0,
            mc_max_vertices, 
            mc_prog_interval,
            ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
            ReaK::pp::DVP_ALT_COB4_KNN,
            ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
            mc_results);
        
        fadprm_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      };
      };
      
    };
    
#endif
    
    
    
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
    
    if(run_all_planners || vm.count("sba-star")) {
      
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
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::sbastar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          sbastar_plan(world_map, 
                       world_map->get_start_pos(), 
                       world_map->get_goal_pos(),
                       vm["sba-potential-cutoff"].as<double>(),
                       vm["sba-density-cutoff"].as<double>(),
                       world_map->get_max_edge_length(),
                       mc_max_vertices, 
                       mc_prog_interval,
                       ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                       ReaK::pp::DVP_BF2_TREE_KNN,
                       ReaK::pp::LAZY_COLLISION_CHECKING,
                       ReaK::pp::PLAN_WITH_VORONOI_PULL,
                       ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                       mc_results);
        
        sbastar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running SBA* with adj-list, dvp-bf4..." << std::endl;
      timing_output << "SBA*, adj-list, dvp-bf4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::sbastar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          sbastar_plan(world_map, 
                       world_map->get_start_pos(), 
                       world_map->get_goal_pos(),
                       vm["sba-potential-cutoff"].as<double>(),
                       vm["sba-density-cutoff"].as<double>(),
                       world_map->get_max_edge_length(),
                       mc_max_vertices, 
                       mc_prog_interval,
                       ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                       ReaK::pp::DVP_BF4_TREE_KNN,
                       ReaK::pp::LAZY_COLLISION_CHECKING,
                       ReaK::pp::PLAN_WITH_VORONOI_PULL,
                       ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                       mc_results);
        
        sbastar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      if(vm.count("mc-cob-tree")) {
      std::cout << "Running SBA* with adj-list, dvp-cob2..." << std::endl;
      timing_output << "SBA*, adj-list, dvp-cob2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::sbastar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          sbastar_plan(world_map, 
                       world_map->get_start_pos(), 
                       world_map->get_goal_pos(),
                       vm["sba-potential-cutoff"].as<double>(),
                       vm["sba-density-cutoff"].as<double>(),
                       world_map->get_max_edge_length(),
                       mc_max_vertices, 
                       mc_prog_interval,
                       ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                       ReaK::pp::DVP_COB2_TREE_KNN,
                       ReaK::pp::LAZY_COLLISION_CHECKING,
                       ReaK::pp::PLAN_WITH_VORONOI_PULL,
                       ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                       mc_results);
        
        sbastar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running SBA* with adj-list, dvp-cob4..." << std::endl;
      timing_output << "SBA*, adj-list, dvp-cob4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::sbastar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          sbastar_plan(world_map, 
                       world_map->get_start_pos(), 
                       world_map->get_goal_pos(),
                       vm["sba-potential-cutoff"].as<double>(),
                       vm["sba-density-cutoff"].as<double>(),
                       world_map->get_max_edge_length(),
                       mc_max_vertices, 
                       mc_prog_interval,
                       ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                       ReaK::pp::DVP_COB4_TREE_KNN,
                       ReaK::pp::LAZY_COLLISION_CHECKING,
                       ReaK::pp::PLAN_WITH_VORONOI_PULL,
                       ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                       mc_results);
        
        sbastar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      };
      
      
      std::cout << "Running SBA* with adj-list, linear-search..." << std::endl;
      timing_output << "SBA*, adj-list, linear-search" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::sbastar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          sbastar_plan(world_map, 
                       world_map->get_start_pos(), 
                       world_map->get_goal_pos(),
                       vm["sba-potential-cutoff"].as<double>(),
                       vm["sba-density-cutoff"].as<double>(),
                       world_map->get_max_edge_length(),
                       mc_max_vertices, 
                       mc_prog_interval,
                       ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                       ReaK::pp::LINEAR_SEARCH_KNN,
                       ReaK::pp::LAZY_COLLISION_CHECKING,
                       ReaK::pp::PLAN_WITH_VORONOI_PULL,
                       ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                       mc_results);
        
        sbastar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      
      if(vm.count("mc-dvp-alt")) {
      std::cout << "Running SBA* with dvp-adj-list-bf2..." << std::endl;
      timing_output << "SBA*, dvp-adj-list-bf2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::sbastar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          sbastar_plan(world_map, 
                       world_map->get_start_pos(), 
                       world_map->get_goal_pos(),
                       vm["sba-potential-cutoff"].as<double>(),
                       vm["sba-density-cutoff"].as<double>(),
                       world_map->get_max_edge_length(),
                       mc_max_vertices, 
                       mc_prog_interval,
                       ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                       ReaK::pp::DVP_ALT_BF2_KNN,
                       ReaK::pp::LAZY_COLLISION_CHECKING,
                       ReaK::pp::PLAN_WITH_VORONOI_PULL,
                       ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                       mc_results);
        
        sbastar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running SBA* with dvp-adj-list-bf4..." << std::endl;
      timing_output << "SBA*, dvp-adj-list-bf4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::sbastar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          sbastar_plan(world_map, 
                       world_map->get_start_pos(), 
                       world_map->get_goal_pos(),
                       vm["sba-potential-cutoff"].as<double>(),
                       vm["sba-density-cutoff"].as<double>(),
                       world_map->get_max_edge_length(),
                       mc_max_vertices, 
                       mc_prog_interval,
                       ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                       ReaK::pp::DVP_ALT_BF4_KNN,
                       ReaK::pp::LAZY_COLLISION_CHECKING,
                       ReaK::pp::PLAN_WITH_VORONOI_PULL,
                       ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                       mc_results);
        
        sbastar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      if(vm.count("mc-cob-tree")) {
      std::cout << "Running SBA* with dvp-adj-list-cob2..." << std::endl;
      timing_output << "SBA*, dvp-adj-list-cob2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::sbastar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          sbastar_plan(world_map, 
                       world_map->get_start_pos(), 
                       world_map->get_goal_pos(),
                       vm["sba-potential-cutoff"].as<double>(),
                       vm["sba-density-cutoff"].as<double>(),
                       world_map->get_max_edge_length(),
                       mc_max_vertices, 
                       mc_prog_interval,
                       ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                       ReaK::pp::DVP_ALT_COB2_KNN,
                       ReaK::pp::LAZY_COLLISION_CHECKING,
                       ReaK::pp::PLAN_WITH_VORONOI_PULL,
                       ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                       mc_results);
        
        sbastar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running SBA* with dvp-adj-list-cob4..." << std::endl;
      timing_output << "SBA*, dvp-adj-list-cob4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::sbastar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          sbastar_plan(world_map, 
                       world_map->get_start_pos(), 
                       world_map->get_goal_pos(),
                       vm["sba-potential-cutoff"].as<double>(),
                       vm["sba-density-cutoff"].as<double>(),
                       world_map->get_max_edge_length(),
                       mc_max_vertices, 
                       mc_prog_interval,
                       ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                       ReaK::pp::DVP_ALT_COB4_KNN,
                       ReaK::pp::LAZY_COLLISION_CHECKING,
                       ReaK::pp::PLAN_WITH_VORONOI_PULL,
                       ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                       mc_results);
        
        sbastar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      };
      };
      
    };
      
#endif
    
    
    
    
    
#ifdef RK_ENABLE_TEST_RRTSTAR_PLANNER
    
    if(run_all_planners || vm.count("rrt-star")) {
      
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
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrtstar_plan(world_map, 
                      world_map->get_start_pos(), 
                      world_map->get_goal_pos(),
                      mc_max_vertices, 
                      mc_prog_interval,
                      ReaK::pp::UNIDIRECTIONAL_RRT,
                      ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                      ReaK::pp::DVP_BF2_TREE_KNN,
                      ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                      mc_results);
        
        rrtstar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT* with Uni-dir, adj-list, dvp-bf4..." << std::endl;
      timing_output << "RRT*, Uni-dir, adj-list, dvp-bf4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrtstar_plan(world_map, 
                      world_map->get_start_pos(), 
                      world_map->get_goal_pos(),
                      mc_max_vertices, 
                      mc_prog_interval,
                      ReaK::pp::UNIDIRECTIONAL_RRT,
                      ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                      ReaK::pp::DVP_BF4_TREE_KNN,
                      ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                      mc_results);
        
        rrtstar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      if(vm.count("mc-cob-tree")) {
      std::cout << "Running RRT* with Uni-dir, adj-list, dvp-cob2..." << std::endl;
      timing_output << "RRT*, Uni-dir, adj-list, dvp-cob2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrtstar_plan(world_map, 
                      world_map->get_start_pos(), 
                      world_map->get_goal_pos(),
                      mc_max_vertices, 
                      mc_prog_interval,
                      ReaK::pp::UNIDIRECTIONAL_RRT,
                      ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                      ReaK::pp::DVP_COB2_TREE_KNN,
                      ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                      mc_results);
        
        rrtstar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT* with Uni-dir, adj-list, dvp-cob4..." << std::endl;
      timing_output << "RRT*, Uni-dir, adj-list, dvp-cob4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrtstar_plan(world_map, 
                      world_map->get_start_pos(), 
                      world_map->get_goal_pos(),
                      mc_max_vertices, 
                      mc_prog_interval,
                      ReaK::pp::UNIDIRECTIONAL_RRT,
                      ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                      ReaK::pp::DVP_COB4_TREE_KNN,
                      ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                      mc_results);
        
        rrtstar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      };
      
      
      std::cout << "Running RRT* with Uni-dir, adj-list, linear-search..." << std::endl;
      timing_output << "RRT*, Uni-dir, adj-list, linear-search" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrtstar_plan(world_map, 
                      world_map->get_start_pos(), 
                      world_map->get_goal_pos(),
                      mc_max_vertices, 
                      mc_prog_interval,
                      ReaK::pp::UNIDIRECTIONAL_RRT,
                      ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                      ReaK::pp::LINEAR_SEARCH_KNN,
                      ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                      mc_results);
        
        rrtstar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      
      if(vm.count("mc-dvp-alt")) {
      std::cout << "Running RRT* with Uni-dir, dvp-adj-list-bf2..." << std::endl;
      timing_output << "RRT*, Uni-dir, dvp-adj-list-bf2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrtstar_plan(world_map, 
                      world_map->get_start_pos(), 
                      world_map->get_goal_pos(),
                      mc_max_vertices, 
                      mc_prog_interval,
                      ReaK::pp::UNIDIRECTIONAL_RRT,
                      ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                      ReaK::pp::DVP_ALT_BF2_KNN,
                      ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                      mc_results);
        
        rrtstar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT* with Uni-dir, dvp-adj-list-bf4..." << std::endl;
      timing_output << "RRT*, Uni-dir, dvp-adj-list-bf4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrtstar_plan(world_map, 
                      world_map->get_start_pos(), 
                      world_map->get_goal_pos(),
                      mc_max_vertices, 
                      mc_prog_interval,
                      ReaK::pp::UNIDIRECTIONAL_RRT,
                      ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                      ReaK::pp::DVP_ALT_BF4_KNN,
                      ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                      mc_results);
        
        rrtstar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      if(vm.count("mc-cob-tree")) {
      std::cout << "Running RRT* with Uni-dir, dvp-adj-list-cob2..." << std::endl;
      timing_output << "RRT*, Uni-dir, dvp-adj-list-cob2" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrtstar_plan(world_map, 
                      world_map->get_start_pos(), 
                      world_map->get_goal_pos(),
                      mc_max_vertices, 
                      mc_prog_interval,
                      ReaK::pp::UNIDIRECTIONAL_RRT,
                      ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                      ReaK::pp::DVP_ALT_COB2_KNN,
                      ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                      mc_results);
        
        rrtstar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      
      
      std::cout << "Running RRT* with Uni-dir, dvp-adj-list-cob4..." << std::endl;
      timing_output << "RRT*, Uni-dir, dvp-adj-list-cob4" << std::endl;
      {
      std::vector< std::pair<int, ReaK::vect<double,3> > > avg_times(mc_max_vertices_100, std::pair<int, ReaK::vect<double,3> >(0,ReaK::vect<double,3>(0.0,1.0e10,0.0)));
      for(std::size_t i = 0; i < mc_run_count; ++i) {
        std::stringstream ss, ss2;
        
        ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReporterType > 
          rrtstar_plan(world_map, 
                      world_map->get_start_pos(), 
                      world_map->get_goal_pos(),
                      mc_max_vertices, 
                      mc_prog_interval,
                      ReaK::pp::UNIDIRECTIONAL_RRT,
                      ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH,
                      ReaK::pp::DVP_ALT_COB4_KNN,
                      ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                      mc_results);
        
        rrtstar_plan.solve_path();
        
        int v_count, t_val; double c_val;
        int j = 0;
        while(ss >> v_count) {
          ss >> t_val; ss2 >> v_count >> c_val;
          avg_times[j].second[0] = (double(t_val) + double(avg_times[j].first) * avg_times[j].second[0]) / double(avg_times[j].first + 1);
          if(c_val < avg_times[j].second[1]) 
            avg_times[j].second[1] = c_val;
          if(c_val > avg_times[j].second[2])
            avg_times[j].second[2] = c_val;
          avg_times[j].first += 1; ++j;
        };
      };
      for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
        if(avg_times[i].first)
          timing_output << std::setw(6) << (i+1)*mc_prog_interval << " " << std::setw(6) << avg_times[i].first << " " << std::setw(10) << avg_times[i].second[0] << " " << std::setw(10) << avg_times[i].second[1] << " " << std::setw(10) << avg_times[i].second[2] << std::endl; 
      };
      };
      std::cout << "Done!" << std::endl;
      };
      };
      
    };
      
#endif
      
  };
  
  
  
  
  if(vm.count("single-run")) {
      
    std::size_t sr_max_vertices     = vm["max-vertices"].as<std::size_t>();
    std::size_t sr_results          = vm["max-results"].as<std::size_t>();
    
#ifdef RK_ENABLE_TEST_URRT_PLANNER
    
    if(run_all_planners || vm.count("rrt")) {
      std::cout << "Outputting RRT with Uni-dir, adj-list, dvp-bf4..." << std::endl;
      
      ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> > > 
        rrt_plan(world_map, 
                world_map->get_start_pos(), 
                world_map->get_goal_pos(),
                sr_max_vertices, 
                10,
                ReaK::pp::UNIDIRECTIONAL_RRT,
                ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                ReaK::pp::DVP_BF4_TREE_KNN,
                ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> >(output_path_name + "/rrt/" + world_file_name_only + "_", 5),
                sr_results);
      
      rrt_plan.solve_path();
      
      std::cout << "Done!" << std::endl;
    };
      
#endif
    
#ifdef RK_ENABLE_TEST_BRRT_PLANNER
    
    if(run_all_planners || vm.count("bi-rrt")) {
      std::cout << "Outputting RRT with Bi-dir, adj-list, dvp-bf4..." << std::endl;
      
      ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> > > 
        rrt_plan(world_map, 
                world_map->get_start_pos(), 
                world_map->get_goal_pos(),
                sr_max_vertices, 
                100,
                ReaK::pp::BIDIRECTIONAL_RRT,
                ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                ReaK::pp::DVP_BF4_TREE_KNN,
                ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> >(output_path_name + "/birrt/" + world_file_name_only + "_", 5),
                sr_results);
      
      rrt_plan.solve_path();
      
      std::cout << "Done!" << std::endl;
    };
    
#endif
    
#ifdef RK_ENABLE_TEST_PRM_PLANNER
    
    if(run_all_planners || vm.count("prm")) {
      std::cout << "Outputting PRM with adj-list dvp-bf4..." << std::endl;
      
      ReaK::pp::prm_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> > > 
        prm_plan(world_map, 
                world_map->get_start_pos(), 
                world_map->get_goal_pos(),
                sr_max_vertices, 
                10,
                ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                ReaK::pp::DVP_BF4_TREE_KNN,
                ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> >(output_path_name + "/prm/" + world_file_name_only + "_", 5),
                sr_results);
      
      prm_plan.solve_path();
      
      std::cout << "Done!" << std::endl;
    };
    
#endif
    
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
    
    if(run_all_planners || vm.count("fadprm")) {
      std::cout << "Outputting FADPRM with adj-list dvp-bf4..." << std::endl;
      
      ReaK::pp::fadprm_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> > > 
        fadprm_plan(
          world_map, 
          world_map->get_start_pos(), 
          world_map->get_goal_pos(),
          0.1,
          sr_max_vertices, 
          10,
          ReaK::pp::ADJ_LIST_MOTION_GRAPH,
          ReaK::pp::DVP_BF4_TREE_KNN,
          ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> >(output_path_name + "/fadprm/" + world_file_name_only + "_", 5),
          sr_results);
      
      fadprm_plan.solve_path();
      
      std::cout << "Done!" << std::endl;
    };
    
#endif
    
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
    
    if(run_all_planners || vm.count("sba-star")) {
      std::cout << "Outputting SBA* with adj-list dvp-bf4..." << std::endl;
      
      ReaK::pp::sbastar_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> > > 
        sbastar_plan(
          world_map, 
          world_map->get_start_pos(), 
          world_map->get_goal_pos(),
          vm["sba-potential-cutoff"].as<double>(),
          vm["sba-density-cutoff"].as<double>(),
          world_map->get_max_edge_length(),
          sr_max_vertices, 
          10,
          ReaK::pp::ADJ_LIST_MOTION_GRAPH,
          ReaK::pp::DVP_BF4_TREE_KNN,
          ReaK::pp::LAZY_COLLISION_CHECKING,
          ReaK::pp::PLAN_WITH_VORONOI_PULL,
          ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> >(output_path_name + "/sbastar/" + world_file_name_only + "_", 5),
          sr_results);
      
      sbastar_plan.solve_path();
      
      std::cout << "Done!" << std::endl;
    };
    
#endif
    
#ifdef RK_ENABLE_TEST_RRTSTAR_PLANNER
    
    if(run_all_planners || vm.count("rrt-star")) {
      std::cout << "Outputting RRT* with Uni-dir, adj-list, dvp-bf4..." << std::endl;
      
      ReaK::pp::rrtstar_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> > > 
        rrtstar_plan(world_map, 
                    world_map->get_start_pos(), 
                    world_map->get_goal_pos(),
                    sr_max_vertices, 
                    10,
                    ReaK::pp::UNIDIRECTIONAL_RRT,
                    ReaK::pp::ADJ_LIST_MOTION_GRAPH,
                    ReaK::pp::DVP_BF4_TREE_KNN,
                    ReaK::pp::differ_sbmp_report_to_space< ReaK::pp::print_sbmp_progress<> >(output_path_name + "/rrt_star/" + world_file_name_only + "_", 5),
                    sr_results);
      
      rrtstar_plan.solve_path();
      
      std::cout << "Done!" << std::endl;
    };
    
#endif
    
  };
  
  return 0;
};















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
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;






#define RK_HIDIM_PLANNER_DO_RRT     0x00000001
#define RK_HIDIM_PLANNER_DO_BIRRT   0x00000002
#define RK_HIDIM_PLANNER_DO_PRM     0x00000004
#define RK_HIDIM_PLANNER_DO_FADPRM  0x00000008
#define RK_HIDIM_PLANNER_DO_SBASTAR 0x00000010
#define RK_HIDIM_PLANNER_DO_RRTSTAR 0x00000020

#define RK_HIDIM_PLANNER_DO_COB 0x00000100
#define RK_HIDIM_PLANNER_DO_ALT 0x00000200




std::size_t mc_run_count = 0;
std::size_t mc_max_vertices = 0;
std::size_t mc_prog_interval = 0;
std::size_t mc_max_vertices_100 = 0;
std::size_t mc_results = 0;
std::size_t mc_flags = 0;

std::size_t data_struct_flags = 0;
std::string data_struct_str = "";


#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
double fadprm_relaxation = 0.0;
#endif

#ifdef RK_ENABLE_TEST_RRTSTAR_PLANNER
std::size_t rrtstar_opt_flags = ReaK::pp::UNIDIRECTIONAL_PLANNING;
std::string rrtstar_qualifier = "";
#endif

#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
double sba_potential_cutoff = 0.0;
double sba_density_cutoff = 0.0;
double sba_relaxation = 0.0;
double sba_sa_temperature = 0.0;
bool sba_use_voronoi_pull = false;
std::size_t sba_opt_flags = ReaK::pp::UNIDIRECTIONAL_PLANNING | ReaK::pp::LAZY_COLLISION_CHECKING;
std::string sba_qualifier = "_lazy";
#endif


template <typename PlannerType>
void run_monte_carlo_tests(
    const PlannerType& planner,
    std::stringstream& time_rec_ss,
    std::stringstream& cost_rec_ss,
    std::ostream& result_output) {
  std::vector< double > vertex_counts(mc_max_vertices_100, 0.0);
  std::vector< std::size_t > num_remaining_planners(mc_max_vertices_100, 0);
  std::vector< std::size_t > num_successful_planners(mc_max_vertices_100, 0);
  
  std::vector< double > time_values(mc_max_vertices_100, 0.0);
  std::vector< double > best_costs(mc_max_vertices_100, 1.0e10);
  std::vector< double > worst_costs(mc_max_vertices_100, 0.0);
  std::vector< double > avg_costs(mc_max_vertices_100, 0.0);
  
  cost_rec_ss << std::fixed;
  
  for(std::size_t i = 0; i < mc_run_count; ++i) {
    time_rec_ss.clear();
    time_rec_ss.seekg(0, time_rec_ss.end);
    cost_rec_ss.clear();
    cost_rec_ss.seekg(0, cost_rec_ss.end);
    
    PlannerType planner_tmp = planner;
    
    planner_tmp.solve_path();
    
    std::size_t v_count = 0, t_val = 0; 
    std::string tmp;
    std::size_t j = 0;
    while( std::getline(time_rec_ss, tmp) && (tmp.size()) ) {
      std::stringstream ss_tmp(tmp);
      ss_tmp >> v_count >> t_val;
      vertex_counts[j] = (double(v_count) + double(num_remaining_planners[j]) * vertex_counts[j]) / double(num_remaining_planners[j] + 1);
      time_values[j] = (double(t_val) + double(num_remaining_planners[j]) * time_values[j]) / double(num_remaining_planners[j] + 1);
      num_remaining_planners[j] += 1; 
      ++j;
    };
    
    double c_val = 1e10;
    j = 0;
    while( std::getline(cost_rec_ss, tmp) && (tmp.size()) ) {
      std::stringstream ss_tmp(tmp);
      ss_tmp >> v_count >> c_val;
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
    
    while(j < mc_max_vertices_100) {
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
  for(std::size_t i = 0; i < mc_max_vertices_100; ++i) {
    result_output << std::setw(9) << i 
           << " " << std::setw(9) << vertex_counts[i] 
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
                            std::ostream& timing_output) {
  
  std::cout << "*****************************************************************" << std::endl
            << "*            Running tests on '" << world_map->getName() << "'" << std::endl
            << "*****************************************************************" << std::endl;
  
  
  typedef ReaK::pp::timing_sbmp_report< ReaK::pp::least_cost_sbmp_report<> > ReporterType;
  
#ifdef RK_ENABLE_TEST_URRT_PLANNER
  
  if(mc_flags & RK_HIDIM_PLANNER_DO_RRT) {
    
    std::cout << "Running RRT with Uni-dir, " << data_struct_str << std::endl;
    timing_output << "RRT, Uni-dir, " << data_struct_str << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
        rrt_plan(world_map, 
                world_map->get_start_pos(), 
                world_map->get_goal_pos(),
                mc_max_vertices, 
                mc_prog_interval,
                data_struct_flags,
                ReaK::pp::UNIDIRECTIONAL_PLANNING,
                ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                mc_results);
      
      run_monte_carlo_tests(rrt_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
  };
  
#endif
  
  
#ifdef RK_ENABLE_TEST_BRRT_PLANNER
  
  if(mc_flags & RK_HIDIM_PLANNER_DO_BIRRT) {
    
    std::cout << "Running RRT with Bi-dir, " << data_struct_str << std::endl;
    timing_output << "RRT, Bi-dir, " << data_struct_str << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::rrt_path_planner< SpaceType, ReporterType > 
        rrt_plan(world_map, 
                world_map->get_start_pos(), 
                world_map->get_goal_pos(),
                mc_max_vertices, 
                mc_prog_interval,
                data_struct_flags,
                ReaK::pp::BIDIRECTIONAL_PLANNING,
                ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                mc_results);
      
      run_monte_carlo_tests(rrt_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
  };
  
#endif
  
  
  
#ifdef RK_ENABLE_TEST_PRM_PLANNER
  
  if(mc_flags & RK_HIDIM_PLANNER_DO_PRM) {
    
    std::cout << "Running PRM with " << data_struct_str << std::endl;
    timing_output << "PRM, " << data_struct_str << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::prm_path_planner< SpaceType, ReporterType > 
        prm_plan(world_map, 
                world_map->get_start_pos(), 
                world_map->get_goal_pos(),
                mc_max_vertices, 
                mc_prog_interval,
                data_struct_flags,
                ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                mc_results);
      
      run_monte_carlo_tests(prm_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
  };
  
#endif
  
  
  
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
  
  if(mc_flags & RK_HIDIM_PLANNER_DO_FADPRM) {
    
    std::cout << "Running FADPRM with " << data_struct_str << std::endl;
    timing_output << "FADPRM, " << data_struct_str << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::fadprm_path_planner< SpaceType, ReporterType > 
        fadprm_plan(
          world_map, 
          world_map->get_start_pos(),
          world_map->get_goal_pos(),
          fadprm_relaxation,
          mc_max_vertices, 
          mc_prog_interval,
          data_struct_flags,
          ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
          mc_results);
      
      run_monte_carlo_tests(fadprm_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
  };
  
#endif
  
  
  
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
  
  if(mc_flags & RK_HIDIM_PLANNER_DO_SBASTAR) {
    
    std::cout << "Running SBA* with " << data_struct_str << std::endl;
    timing_output << "SBA*, " << data_struct_str << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::sbastar_path_planner< SpaceType, ReporterType > 
        sbastar_plan(world_map, 
                    world_map->get_start_pos(), 
                    world_map->get_goal_pos(),
                    mc_max_vertices, 
                    mc_prog_interval,
                    data_struct_flags,
                    sba_opt_flags,
                    ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                    mc_results);
      
      sbastar_plan.set_initial_key_threshold(sba_potential_cutoff);
      sbastar_plan.set_initial_density_threshold(sba_density_cutoff);
      sbastar_plan.set_initial_relaxation(sba_relaxation);
      sbastar_plan.set_initial_SA_temperature(sba_sa_temperature);
      sbastar_plan.set_sampling_radius( world_map->get_max_edge_length() );
      
      run_monte_carlo_tests(sbastar_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
  };
    
#endif
  
  
#ifdef RK_ENABLE_TEST_RRTSTAR_PLANNER
  
  if(mc_flags & RK_HIDIM_PLANNER_DO_RRTSTAR) {
    
    std::cout << "Running RRT* with Uni-dir, " << data_struct_str << std::endl;
    timing_output << "RRT*, Uni-dir, " << data_struct_str << std::endl;
    {
      std::stringstream ss, ss2;
      
      ReaK::pp::rrtstar_path_planner< SpaceType, ReporterType > 
        rrtstar_plan(world_map, 
                    world_map->get_start_pos(), 
                    world_map->get_goal_pos(),
                    mc_max_vertices, 
                    mc_prog_interval,
                    data_struct_flags,
                    rrtstar_opt_flags,
                    ReporterType(ss, ReaK::pp::least_cost_sbmp_report<>(ss2)),
                    mc_results);
      
      run_monte_carlo_tests(rrtstar_plan,ss,ss2,timing_output);
    };
    std::cout << "Done!" << std::endl;
    
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
    ("mc-space", po::value< std::string >()->default_value("e3"), "the type of configuration space to use for the monte-carlo runs (default: e3 (euclidean-3)), can range from e3 to e20.")
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
    ("rrt-star-with-bnb", "specify whether to use a Branch-and-bound or not during RRT* as a method to prune useless nodes from the motion-graph")
#endif
#ifdef RK_ENABLE_TEST_PRM_PLANNER
    ("prm", "specify that the PRM algorithm should be run")
#endif
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
    ("fadprm", "specify that the FADPRM algorithm should be run")
    ("fadprm-relaxation", po::value< double >()->default_value(10.0), "specify the initial relaxation factor for the FADPRM algorithm (default: 10.0)")
#endif
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
    ("sba-star", "specify that the SBA* algorithm should be run")
    ("sba-potential-cutoff", po::value< double >()->default_value(0.0), "specify the potential cutoff for the SBA* algorithm")
    ("sba-density-cutoff", po::value< double >()->default_value(0.5), "specify the density cutoff for the SBA* algorithm")
    ("sba-relaxation", po::value< double >()->default_value(0.0), "specify the initial relaxation factor for the Anytime SBA* algorithm")
    ("sba-with-voronoi-pull", "specify whether to use a Voronoi pull or not as a method to add an exploratory bias to the search")
    ("sba-sa-temperature", po::value< double >()->default_value(-1.0), "specify the initial Simulated Annealing temperature for the SBA*-RRT* algorithms")
    ("sba-with-bnb", "specify whether to use a Branch-and-bound or not during SBA* as a method to prune useless nodes from the motion-graph")
#endif
    ("all-planners,a", "specify that all supported planners should be run (default if no particular planner is specified)")
    ("knn-method", po::value< std::string >()->default_value("bf2"), "specify the KNN method to use (options: linear, bf2, bf4, cob2, cob4) (default: bf2)")
    ("mg-storage", po::value< std::string >()->default_value("adj-list"), "specify the KNN method to use (options: adj-list, dvp-adj-list) (default: adj-list)")
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
  
  fs::create_directory(output_path_name.c_str());
  
  bool run_all_planners = false;
  if(vm.count("all-planners") || (vm.count("rrt") + vm.count("bi-rrt") + vm.count("rrt-star") + vm.count("prm") + vm.count("fadprm") + vm.count("sba-star") == 0)) 
    run_all_planners = true;
  
  
  data_struct_flags = 0;
  std::string knn_method_str = "bf2";
  if((vm["knn-method"].as<std::string>() == "linear") && (vm["mg-storage"].as<std::string>() == "adj-list")) {
    data_struct_flags |= ReaK::pp::LINEAR_SEARCH_KNN;
    knn_method_str = "linear";
  } else if(vm["knn-method"].as<std::string>() == "bf4") {
    data_struct_flags |= ReaK::pp::DVP_BF4_TREE_KNN;
    knn_method_str = "bf4";
  } else if(vm["knn-method"].as<std::string>() == "cob2") {
    data_struct_flags |= ReaK::pp::DVP_COB2_TREE_KNN;
    knn_method_str = "cob2";
  } else if(vm["knn-method"].as<std::string>() == "cob4") {
    data_struct_flags |= ReaK::pp::DVP_COB4_TREE_KNN;
    knn_method_str = "cob4";
  } else {
    data_struct_flags |= ReaK::pp::DVP_BF2_TREE_KNN;
  };
  
  std::string mg_storage_str = "adj-list";
  if(vm["mg-storage"].as<std::string>() == "dvp-adj-list") {
    data_struct_flags |= ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH;
    mg_storage_str = "dvp-adj-list";
  } else {
    data_struct_flags |= ReaK::pp::ADJ_LIST_MOTION_GRAPH;
  };
  
  data_struct_str = mg_storage_str + ", " + knn_method_str;
  
  
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
  fadprm_relaxation       = vm["fadprm-relaxation"].as<double>();
#endif
  
#ifdef RK_ENABLE_TEST_RRTSTAR_PLANNER
  rrtstar_opt_flags = ReaK::pp::UNIDIRECTIONAL_PLANNING;
  rrtstar_qualifier = "";
  
  if( vm.count("rrt-star-with-bnb") ) {
    rrtstar_opt_flags |= ReaK::pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG;
    rrtstar_qualifier += "_bnb";
  };
#endif
  
  
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
  sba_opt_flags = ReaK::pp::UNIDIRECTIONAL_PLANNING | ReaK::pp::LAZY_COLLISION_CHECKING;
  sba_qualifier = "_lazy";
  
  if(vm["sba-relaxation"].as<double>() > 1e-6) {
    sba_opt_flags |= ReaK::pp::PLAN_WITH_ANYTIME_HEURISTIC;
    sba_qualifier += "_any";
  };
  
  if( vm.count("sba-with-voronoi-pull") ) {
    sba_opt_flags |= ReaK::pp::PLAN_WITH_VORONOI_PULL;
    sba_qualifier += "_sa";
  };
  
  if( vm.count("sba-with-bnb") ) {
    sba_opt_flags |= ReaK::pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG;
    sba_qualifier += "_bnb";
  };
  
  sba_potential_cutoff = vm["sba-potential-cutoff"].as<double>();
  sba_density_cutoff   = vm["sba-density-cutoff"].as<double>();
  sba_relaxation       = vm["sba-relaxation"].as<double>();
  sba_sa_temperature   = vm["sba-sa-temperature"].as<double>();
  sba_use_voronoi_pull = vm.count("sba-with-voronoi-pull");
#endif
  
  
  
  mc_run_count        = vm["mc-runs"].as<std::size_t>();
  mc_max_vertices     = vm["mc-vertices"].as<std::size_t>();
  mc_prog_interval    = vm["mc-prog-interval"].as<std::size_t>();
  mc_max_vertices_100 = mc_max_vertices / mc_prog_interval;
  mc_results          = vm["mc-results"].as<std::size_t>();
  
  
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
  
  
  
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double, 3> > > World3DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double, 4> > > World4DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double, 5> > > World5DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double, 6> > > World6DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double, 7> > > World7DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double, 8> > > World8DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double, 9> > > World9DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double,10> > > World10DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double,11> > > World11DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double,12> > > World12DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double,13> > > World13DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double,14> > > World14DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double,15> > > World15DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double,16> > > World16DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double,17> > > World17DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double,18> > > World18DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double,19> > > World19DType;
  typedef ReaK::pp::no_obstacle_space< ReaK::pp::hyperbox_topology< ReaK::vect<double,20> > > World20DType;
  
  
  if(vm["mc-space"].as<std::string>() == "e3") {
    ReaK::shared_ptr< World3DType > world_3D =
      ReaK::shared_ptr< World3DType >(
        new World3DType(
          "world_3D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 3> >("world_3D",
                                                              ReaK::vect<double,3>(0.0,0.0,0.0),
                                                              ReaK::vect<double,3>(1.0,1.0,1.0)),
          0.1 * std::sqrt(3.0)));
    world_3D->set_start_pos(ReaK::vect<double,3>(0.05,0.05,0.05));
    world_3D->set_goal_pos(ReaK::vect<double,3>(0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e3_times.txt");
    
    test_planners_on_space(world_3D, timing_output);
    
  } else if(vm["mc-space"].as<std::string>() == "e4") {
    ReaK::shared_ptr< World4DType > world_4D =
      ReaK::shared_ptr< World4DType >(
        new World4DType(
          "world_4D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 4> >("world_4D",
                                                              ReaK::vect<double,4>(0.0,0.0,0.0,0.0),
                                                              ReaK::vect<double,4>(1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(4.0)));
    world_4D->set_start_pos(ReaK::vect<double,4>(0.05,0.05,0.05,0.05));
    world_4D->set_goal_pos(ReaK::vect<double,4>(0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e4_times.txt");
    
    test_planners_on_space(world_4D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e5") {
    ReaK::shared_ptr< World5DType > world_5D =
      ReaK::shared_ptr< World5DType >(
        new World5DType(
          "world_5D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 5> >("world_5D",
                                                              ReaK::vect<double,5>(0.0,0.0,0.0,0.0,0.0),
                                                              ReaK::vect<double,5>(1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(5.0)));
    world_5D->set_start_pos(ReaK::vect<double,5>(0.05,0.05,0.05,0.05,0.05));
    world_5D->set_goal_pos(ReaK::vect<double,5>(0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e5_times.txt");
    
    test_planners_on_space(world_5D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e6") {
    ReaK::shared_ptr< World6DType > world_6D =
      ReaK::shared_ptr< World6DType >(
        new World6DType(
          "world_6D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 6> >("world_6D",
                                                              ReaK::vect<double,6>(0.0,0.0,0.0,0.0,0.0,0.0),
                                                              ReaK::vect<double,6>(1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(6.0)));
    world_6D->set_start_pos(ReaK::vect<double,6>(0.05,0.05,0.05,0.05,0.05,0.05));
    world_6D->set_goal_pos(ReaK::vect<double,6>(0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e6_times.txt");
    
    test_planners_on_space(world_6D, timing_output);
    
    
  }
#if 0
  else if(vm["mc-space"].as<std::string>() == "e7") {
    ReaK::shared_ptr< World7DType > world_7D =
      ReaK::shared_ptr< World7DType >(
        new World7DType(
          "world_7D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 7> >("world_7D",
                                                              ReaK::vect<double,7>(0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                              ReaK::vect<double,7>(1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(7.0)));
    world_7D->set_start_pos(ReaK::vect<double,7>(0.05,0.05,0.05,0.05,0.05,0.05,0.05));
    world_7D->set_goal_pos(ReaK::vect<double,7>(0.95,0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e7_times.txt");
    
    test_planners_on_space(world_7D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e8") {
    ReaK::shared_ptr< World8DType > world_8D =
      ReaK::shared_ptr< World8DType >(
        new World8DType(
          "world_8D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 8> >("world_8D",
                                                              ReaK::vect<double,8>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                              ReaK::vect<double,8>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(8.0)));
    world_8D->set_start_pos(ReaK::vect<double,8>(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
    world_8D->set_goal_pos(ReaK::vect<double,8>(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e8_times.txt");
    
    test_planners_on_space(world_8D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e9") {
    ReaK::shared_ptr< World9DType > world_9D =
      ReaK::shared_ptr< World9DType >(
        new World9DType(
          "world_9D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 9> >("world_9D",
                                                              ReaK::vect<double,9>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                              ReaK::vect<double,9>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(9.0)));
    world_9D->set_start_pos(ReaK::vect<double,9>(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
    world_9D->set_goal_pos(ReaK::vect<double,9>(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e9_times.txt");
    
    test_planners_on_space(world_9D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e10") {
    ReaK::shared_ptr< World10DType > world_10D =
      ReaK::shared_ptr< World10DType >(
        new World10DType(
          "world_10D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 10> >("world_10D",
                                                              ReaK::vect<double,10>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                              ReaK::vect<double,10>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(10.0)));
    world_10D->set_start_pos(ReaK::vect<double,10>(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
    world_10D->set_goal_pos(ReaK::vect<double,10>(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e10_times.txt");
    
    test_planners_on_space(world_10D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e11") {
    ReaK::shared_ptr< World11DType > world_11D =
      ReaK::shared_ptr< World11DType >(
        new World11DType(
          "world_11D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 11> >("world_11D",
                                                                ReaK::vect<double,11>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                                ReaK::vect<double,11>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(11.0)));
    world_11D->set_start_pos(ReaK::vect<double,11>(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
    world_11D->set_goal_pos( ReaK::vect<double,11>(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e11_times.txt");
    
    test_planners_on_space(world_11D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e12") {
    ReaK::shared_ptr< World12DType > world_12D =
      ReaK::shared_ptr< World12DType >(
        new World12DType(
          "world_12D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 12> >("world_12D",
                                                                ReaK::vect<double,12>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                                ReaK::vect<double,12>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(12.0)));
    world_12D->set_start_pos(ReaK::vect<double,12>(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
    world_12D->set_goal_pos( ReaK::vect<double,12>(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e12_times.txt");
    
    test_planners_on_space(world_12D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e13") {
    ReaK::shared_ptr< World13DType > world_13D =
      ReaK::shared_ptr< World13DType >(
        new World13DType(
          "world_13D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 13> >("world_13D",
                                                                ReaK::vect<double,13>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                                ReaK::vect<double,13>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(13.0)));
    world_13D->set_start_pos(ReaK::vect<double,13>(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
    world_13D->set_goal_pos( ReaK::vect<double,13>(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e13_times.txt");
    
    test_planners_on_space(world_13D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e14") {
    ReaK::shared_ptr< World14DType > world_14D =
      ReaK::shared_ptr< World14DType >(
        new World14DType(
          "world_14D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 14> >("world_14D",
                                                                ReaK::vect<double,14>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                                ReaK::vect<double,14>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(14.0)));
    world_14D->set_start_pos(ReaK::vect<double,14>(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
    world_14D->set_goal_pos( ReaK::vect<double,14>(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e14_times.txt");
    
    test_planners_on_space(world_14D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e15") {
    ReaK::shared_ptr< World15DType > world_15D =
      ReaK::shared_ptr< World15DType >(
        new World15DType(
          "world_15D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 15> >("world_15D",
                                                                ReaK::vect<double,15>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                                ReaK::vect<double,15>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(15.0)));
    world_15D->set_start_pos(ReaK::vect<double,15>(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
    world_15D->set_goal_pos( ReaK::vect<double,15>(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e15_times.txt");
    
    test_planners_on_space(world_15D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e16") {
    ReaK::shared_ptr< World16DType > world_16D =
      ReaK::shared_ptr< World16DType >(
        new World16DType(
          "world_16D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 16> >("world_16D",
                                                                ReaK::vect<double,16>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                                ReaK::vect<double,16>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(16.0)));
    world_16D->set_start_pos(ReaK::vect<double,16>(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
    world_16D->set_goal_pos( ReaK::vect<double,16>(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e16_times.txt");
    
    test_planners_on_space(world_16D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e17") {
    ReaK::shared_ptr< World17DType > world_17D =
      ReaK::shared_ptr< World17DType >(
        new World17DType(
          "world_17D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 17> >("world_17D",
                                                                ReaK::vect<double,17>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                                ReaK::vect<double,17>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(17.0)));
    world_17D->set_start_pos(ReaK::vect<double,17>(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
    world_17D->set_goal_pos( ReaK::vect<double,17>(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e17_times.txt");
    
    test_planners_on_space(world_17D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e18") {
    ReaK::shared_ptr< World18DType > world_18D =
      ReaK::shared_ptr< World18DType >(
        new World18DType(
          "world_18D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 18> >("world_18D",
                                                                ReaK::vect<double,18>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                                ReaK::vect<double,18>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(18.0)));
    world_18D->set_start_pos(ReaK::vect<double,18>(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
    world_18D->set_goal_pos( ReaK::vect<double,18>(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e18_times.txt");
    
    test_planners_on_space(world_18D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e19") {
    ReaK::shared_ptr< World19DType > world_19D =
      ReaK::shared_ptr< World19DType >(
        new World19DType(
          "world_19D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 19> >("world_19D",
                                                                ReaK::vect<double,19>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                                ReaK::vect<double,19>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(19.0)));
    world_19D->set_start_pos(ReaK::vect<double,19>(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
    world_19D->set_goal_pos( ReaK::vect<double,19>(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e19_times.txt");
    
    test_planners_on_space(world_19D, timing_output);
    
    
  } else if(vm["mc-space"].as<std::string>() == "e20") {
    ReaK::shared_ptr< World20DType > world_20D =
      ReaK::shared_ptr< World20DType >(
        new World20DType(
          "world_20D_no_obstacles",
          ReaK::pp::hyperbox_topology< ReaK::vect<double, 20> >("world_20D",
                                                                ReaK::vect<double,20>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                                                                ReaK::vect<double,20>(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)),
          0.1 * std::sqrt(20.0)));
    world_20D->set_start_pos(ReaK::vect<double,20>(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05));
    world_20D->set_goal_pos( ReaK::vect<double,20>(0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95));
    
    std::ofstream timing_output(output_path_name + "/e20_times.txt");
    
    test_planners_on_space(world_20D, timing_output);
    
    
  };
#endif
  
  
#if 0
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
#endif
  
  
  return 0;
};












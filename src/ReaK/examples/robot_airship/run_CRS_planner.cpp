
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

// #define RK_ENABLE_TEST_RRT_PLANNER
#define RK_ENABLE_TEST_RRTSTAR_PLANNER
// #define RK_ENABLE_TEST_PRM_PLANNER
// #define RK_ENABLE_TEST_FADPRM_PLANNER
#define RK_ENABLE_TEST_SBASTAR_PLANNER


#if defined(RK_ENABLE_TEST_RRT_PLANNER)
#include "path_planning/rrt_path_planner.hpp"
#endif

#if defined(RK_ENABLE_TEST_PRM_PLANNER)
#include "path_planning/prm_path_planner.hpp"
#endif

#if defined(RK_ENABLE_TEST_RRTSTAR_PLANNER)
#include "path_planning/rrtstar_manip_planners.hpp"
#endif

#if defined(RK_ENABLE_TEST_FADPRM_PLANNER)
#include "path_planning/fadprm_path_planner.hpp"
#endif

#if defined(RK_ENABLE_TEST_SBASTAR_PLANNER)
#include "path_planning/sbastar_manip_planners.hpp"
#endif

#include "path_planning/p2p_planning_query.hpp"
#include "path_planning/path_planner_options.hpp"

#include "optimization/optim_exceptions.hpp"

#include "mbd_kte/kte_map_chain.hpp"
#include "kte_models/manip_dynamics_model.hpp"

#include "topologies/manip_planning_traits.hpp"
#include "topologies/manip_P3R3R_workspaces.hpp"
#include "topologies/Ndof_linear_spaces.hpp"
#include "topologies/Ndof_cubic_spaces.hpp"
#include "topologies/Ndof_quintic_spaces.hpp"
#include "topologies/Ndof_svp_spaces.hpp"
#include "topologies/Ndof_sap_spaces.hpp"


#include "basic_sbmp_reporters.hpp"
#include "vlist_sbmp_report.hpp"



#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;



typedef ReaK::pp::ptrobot2D_test_world TestTopology;
typedef ReaK::pp::timing_sbmp_report< ReaK::pp::least_cost_sbmp_report<> > MCReporterType;

void run_monte_carlo_tests(
    std::size_t mc_run_count,
    std::size_t mc_num_records,
    ReaK::pp::sample_based_planner< TestTopology >& planner,
    ReaK::pp::planning_query< TestTopology >& mc_query,
    std::stringstream& time_rec_ss,
    std::stringstream& cost_rec_ss,
    std::stringstream& sol_rec_ss,
    std::ostream& result_output,
    std::ostream& first_sol_event_output) {
  std::vector< double > vertex_counts(mc_num_records, 0.0);
  std::vector< std::size_t > num_remaining_planners(mc_num_records, 0);
  std::vector< std::size_t > num_successful_planners(mc_num_records, 0);
  
  std::vector< double > time_values(mc_num_records, 0.0);
  std::vector< double > best_costs(mc_num_records, 1.0e10);
  std::vector< double > worst_costs(mc_num_records, 0.0);
  std::vector< double > avg_costs(mc_num_records, 0.0);
  
  cost_rec_ss << std::fixed;
  sol_rec_ss << std::fixed;
  
  for(std::size_t i = 0; i < mc_run_count; ++i) {
    time_rec_ss.clear();
    time_rec_ss.seekg(0, time_rec_ss.end);
    cost_rec_ss.clear();
    cost_rec_ss.seekg(0, cost_rec_ss.end);
    sol_rec_ss.clear();
    sol_rec_ss.seekg(0, sol_rec_ss.end);
    
    mc_query.reset_solution_records();
    planner.reset_internal_state();
    planner.solve_planning_query(mc_query);
    
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
    
    std::string first_sol_event;
    std::getline(sol_rec_ss, first_sol_event);
    if(first_sol_event != "")
      first_sol_event_output << first_sol_event << std::endl;
  };
  for(std::size_t i = 0; i < mc_num_records; ++i) {
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






/****************************** Notes on the CRS path-planner code **********************************

Input files:

 - planner-config : instance of the planning_option_collection class. Contains most of the planning options 
                    like method options, storage, max-results, progress interval, max-vertices, relax factor,
                    initial SA temperature, etc... Many of these could be overridden with CL options or for 
                    Monte-Carlo needs.
 - chaser-target data : instance of the chaser_target_data class. Contains the KTE / geom / proxy 
                        models for the chaser, target and environment for the planning scenario.
 - start joint position  (could be embedded in the chaser model)
 - target pose           (could be embedded in the target model)
 - space configuration : specifies the order, interpolator, min/max travel, temporality, 
                         rate-limited'ness, and output space-order.

*****************************************************************************************************/



int main(int argc, char** argv) {
  
  std::string config_file;
  
  po::options_description generic_options("Generic options");
  generic_options.add_options()
    ("help,h", "produce this help message.")
    ("config,c", po::value< std::string >(&config_file)->default_value("run_CRS_planner.cfg"),
                  "configuration file-name (can contain any or all options, will be overriden by command-line options).")
  ;
  
  po::options_description io_options("I/O options");
  io_options.add_options()
    ("planner-options", po::value< std::string >(), "specify the file containing the planner-options data.")
    ("chaser-target-env", po::value< std::string >(), "specify the file containing the chaser-target-env models.")
    ("space-definition", po::value< std::string >(), "specify the file containing the space settings (order, interp, etc.).")
    ("start-configuration", po::value< std::string >(), "specify the file containing the start configuration of the chaser (P3R3R-manipulator), if not specified, the chaser configuration of the chaser model will be used.")
    ("target-pose", po::value< std::string >(), "specify the file containing the target pose of the capture target (satellite / airship), if not specified, the pose of the target model will be used.")
    
    ("output-path,o", po::value< std::string >()->default_value("pp_results"), "specify the output path (default is pp_results)")
  ;
  
  po::options_description mc_options("Monte-Carlo options");
  mc_options.add_options()
    ("monte-carlo,m", "specify that monte-carlo runs should be performed (default is not)")
    ("mc-runs", po::value< std::size_t >()->default_value(100), "number of monte-carlo runs to average out (default is 100)")
  ;
  
  po::options_description single_options("Single-run options");
  single_options.add_options()
    ("single-run,s", "specify that single runs should be performed (default is not)")
  ;
  
  std::string available_algs;
  {
    std::stringstream available_alg_names;
    std::ostream_iterator< const char* > alg_names_iter(available_alg_names, ", ");
#ifdef RK_ENABLE_TEST_RRT_PLANNER
    *(alg_names_iter++) = "rrt";
#endif
#ifdef RK_ENABLE_TEST_RRTSTAR_PLANNER
    *(alg_names_iter++) = "rrt_star";
#endif
#ifdef RK_ENABLE_TEST_PRM_PLANNER
    *(alg_names_iter++) = "prm";
#endif
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
    *(alg_names_iter++) = "fadprm";
#endif
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
    *(alg_names_iter++) = "sba_star";
#endif
    available_algs = available_alg_names.str();
  };
  
  po::options_description planner_select_options("Planner selection options");
  planner_select_options.add_options()
    ("bi-directional", "specify whether to use a bi-directional algorithm or not during planning. Only supported for some algorithms (RRT).")
    ("with-bnb", "specify whether to use a Branch-and-bound or not during planning to prune useless nodes from the motion-graph. Only supported for optimizing algorithms.")
    ("relaxation-factor", po::value< double >()->default_value(0.0), "specify the initial relaxation factor for the algorithm (default: 0.0). Only supported for heuristic-driven algorithms.")
    ("density-cutoff", po::value< double >()->default_value(0.0), "specify the density cutoff (default: 0.0). Only supported for density-driven algorithms.")
    ("with-voronoi-pull", "specify whether to use a Voronoi pull or not to add an exploratory bias to the search (default: not).")
    ("sa-temperature", po::value< double >()->default_value(-1.0), "specify the initial Simulated Annealing temperature for algorithms that work on a exploration-exploitation schedule (e.g., SA-SBA*).")
    ("planner-alg", po::value< std::string >(), "specify the planner algorithm to use, can be any of (" + available_algs + ").")

    ("knn-method", po::value< std::string >(), 
#ifdef RK_PLANNERS_ENABLE_VEBL_TREE
     "specify the KNN method to use (supported options: linear, bf2, bf4, cob2, cob4) (default: bf2)"
#else
     "specify the KNN method to use (supported options: linear, bf2, bf4) (default: bf2)"
#endif
    )

    ("mg-storage", po::value< std::string >(), 
#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT
     "specify the KNN method to use (supported options: adj-list, dvp-adj-list) (default: adj-list)"
#else
     "specify the KNN method to use (supported options: adj-list) (default: adj-list)"
#endif
    )
    
    ("max-vertices", po::value< std::size_t >()->default_value(5000), "maximum number of vertices during runs (default is 5000)")
    ("max-results", po::value< std::size_t >()->default_value(50), "maximum number of result-paths during runs (default is 50)")
    ("prog-interval", po::value< std::size_t >()->default_value(10), "number of vertices between progress reports during runs (default is 10)")
  ;
  
  po::options_description generate_options("File generation options");
  generate_options.add_options()
    ("generate-all-files", po::value< std::string >(), "specify that all configuration files should be generated with the given file-name prefix (file-name without suffix and extension).")
    ("generate-planner-options", po::value< std::string >(), "specify that the planner options file should be generated with the given file-name prefix (file-name without extension).")
    ("generate-chaser-target-env", po::value< std::string >(), "specify that the chaser-target-env file should be generated with the given file-name prefix (file-name without extension).")
    ("generate-space-definition", po::value< std::string >(), "specify that the space-definition file should be generated with the given file-name prefix (file-name without extension).")
    ("generate-start-config", po::value< std::string >(), "specify that the start-configuration file should be generated with the given file-name prefix (file-name without extension).")
    ("generate-target-pose", po::value< std::string >(), "specify that the target-pose file should be generated with the given file-name prefix (file-name without extension).")
    
    ("generate-xml",      "if set, output results in XML format (rkx) (default).")
    ("generate-protobuf", "if set, output results in protobuf format (pbuf).")
    ("generate-binary",   "if set, output results in binary format (rkb).")
  ;
  
  
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options).add(io_options).add(mc_options).add(single_options).add(planner_select_options).add(generate_options);
  
  po::options_description config_file_options;
  config_file_options.add(io_options).add(mc_options).add(single_options).add(planner_select_options);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);
  
  {
    std::ifstream ifs(config_file.c_str());
    if(ifs) {
      po::store(po::parse_config_file(ifs, config_file_options), vm);
      po::notify(vm);
    };
  };
  
  
  
  if(vm.count("help")) {
    std::cout << cmdline_options << std::endl;
    return 0;
  };
  
  if( vm.count("monte-carlo") + vm.count("single-run")
       + vm.count("generate-all-files") + vm.count("generate-planner-options") + vm.count("generate-chaser-target-env") 
       + vm.count("generate-space-definition") + vm.count("generate-start-config") + vm.count("generate-target-pose") < 1 ) {
    std::cout << "Error: There was no action specified! This program is designed to perform Monte-Carlo runs, single runs (with output), or generate the configuration files to construct scenarios. You must specify at least one of these actions to be performed!" << std::endl;
    std::cout << cmdline_options << std::endl;
    return 1;
  };
  
  if( !vm.count("planner-alg") || ( available_algs.find(vm["planner-alg"].as<std::string>()) == std::string::npos ) ) {
    std::cout << "Error: Invalid planning algorithm selected! The planner '" << vm["planner-alg"].as<std::string>() << "' is not supported, the list of supported algorithms is (" << available_algs << ")." << std::endl;
    std::cout << cmdline_options << std::endl;
    return 2;
  };
  
  std::string world_file_name = vm["input-file"].as<std::string>();
  std::string world_file_name_only(std::find(world_file_name.rbegin(),world_file_name.rend(),'/').base(), std::find(world_file_name.rbegin(),world_file_name.rend(),'.').base()-1);
  
  std::string output_path_name = vm["output-path"].as<std::string>();
  while(output_path_name[output_path_name.length()-1] == '/') 
    output_path_name.erase(output_path_name.length()-1, 1);
  
  fs::create_directory(output_path_name.c_str());
  
  
  ReaK::pp::planning_option_collection plan_options;
  /* default planner options: */
  plan_options.planning_algo = 0;
  plan_options.planning_options = 0;
  plan_options.max_vertices = 5000;
  plan_options.prog_interval = 10;
  plan_options.max_results = 50;
  plan_options.knn_method = 0;
  plan_options.store_policy = 0;
  plan_options.init_SA_temp = 0.0;
  plan_options.init_relax = 0.0;
  plan_options.max_random_walk = 1.0;
  plan_options.start_delay = 20.0;
  
  
  /* NOTE: Insert input-file loading code here. */
  
  
  if( vm["max-vertices"].as< std::size_t >() != 5000 )
    plan_options.max_vertices = vm["max-vertices"].as< std::size_t >();
  if( vm["max-results"].as< std::size_t >() != 50 )
    plan_options.max_results = vm["max-results"].as< std::size_t >();
  if( vm["prog-interval"].as< std::size_t >() != 10 )
    plan_options.prog_interval = vm["prog-interval"].as< std::size_t >();
  
  if(vm.count("knn-method")) {
    plan_options.knn_method = 0;
    if((vm["knn-method"].as<std::string>() == "linear") && (vm["mg-storage"].as<std::string>() == "adj-list"))
      plan_options.knn_method |= ReaK::pp::LINEAR_SEARCH_KNN;
    else if(vm["knn-method"].as<std::string>() == "bf4")
      plan_options.knn_method |= ReaK::pp::DVP_BF4_TREE_KNN;
#ifdef RK_PLANNERS_ENABLE_VEBL_TREE
    else if(vm["knn-method"].as<std::string>() == "cob2")
      plan_options.knn_method |= ReaK::pp::DVP_COB2_TREE_KNN;
    else if(vm["knn-method"].as<std::string>() == "cob4")
      plan_options.knn_method |= ReaK::pp::DVP_COB4_TREE_KNN;
#endif
    else
      plan_options.knn_method |= ReaK::pp::DVP_BF2_TREE_KNN;
  };
  
  std::string knn_method_str = "bf2";
  if( plan_options.knn_method & ReaK::pp::LINEAR_SEARCH_KNN )
    knn_method_str = "linear";
  else if( plan_options.knn_method & ReaK::pp::DVP_BF4_TREE_KNN )
    knn_method_str = "bf4";
  else if( plan_options.knn_method & ReaK::pp::DVP_COB2_TREE_KNN )
    knn_method_str = "cob2";
  else if( plan_options.knn_method & ReaK::pp::DVP_COB4_TREE_KNN )
    knn_method_str = "cob4";
  
  if( vm.count("mg-storage") ) {
    plan_options.store_policy = 0;
#ifdef RK_PLANNERS_ENABLE_DVP_ADJ_LIST_LAYOUT
    if(vm["mg-storage"].as<std::string>() == "dvp-adj-list")
      plan_options.store_policy |= ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH;
    else 
#endif
      plan_options.store_policy |= ReaK::pp::ADJ_LIST_MOTION_GRAPH;
  };
  
  std::string mg_storage_str = "adj-list";
  if( plan_options.store_policy & ReaK::pp::DVP_ADJ_LIST_MOTION_GRAPH )
    mg_storage_str = "dvp-adj-list";
  
  
  plan_options.planning_options |= ReaK::pp::LAZY_COLLISION_CHECKING;  // never use eager, always lazy, if supported.
  std::string planner_qualifier_str = "";
  
  if( vm.count("bi-directional") )  {
    plan_options.planning_options |= ReaK::pp::BIDIRECTIONAL_PLANNING;
    planner_qualifier_str += "_bidir";
  };
  
  if( vm.count("with-bnb") ) {
    plan_options.planning_options |= ReaK::pp::USE_BRANCH_AND_BOUND_PRUNING_FLAG;
    planner_qualifier_str += "_bnb";
  };
  
  if(vm["relaxation-factor"].as<double>() > 1e-6) {
    plan_options.init_relax = vm["relaxation-factor"].as<double>();
    plan_options.planning_options |= ReaK::pp::PLAN_WITH_ANYTIME_HEURISTIC;
    planner_qualifier_str += "_any";
  };
  
  if(vm["sa-temperature"].as<double>() > -1.0)  // default has been overriden
    plan_options.init_relax = vm["sa-temperature"].as<double>();
  
  if( vm.count("with-voronoi-pull") ) {
    plan_options.planning_options |= ReaK::pp::PLAN_WITH_VORONOI_PULL;
    planner_qualifier_str += "_sa";
  };
  
//   ("density-cutoff", po::value< double >()->default_value(0.0), "specify the density cutoff (default: 0.0). Only supported for density-driven algorithms.")
  
  
  ReaK::shared_ptr< ReaK::pp::ptrobot2D_test_world > world_map =
    ReaK::shared_ptr< ReaK::pp::ptrobot2D_test_world >(new ReaK::pp::ptrobot2D_test_world(world_file_name, 20, 1.0));
  
  if(vm.count("monte-carlo")) {
    
    std::size_t mc_run_count        = vm["mc-runs"].as<std::size_t>();
    std::size_t mc_max_vertices_100 = plan_options.max_vertices / plan_options.prog_interval;
    
    std::ofstream timing_output(output_path_name + "/" + world_file_name_only + "_times.txt");
    std::ofstream sol_events_output(output_path_name + "/" + world_file_name_only + "_solutions.txt");
    
    std::stringstream time_ss, cost_ss, sol_ss;
    
    //typedef ReaK::pp::timing_sbmp_report< ReaK::pp::least_cost_sbmp_report<> > ReporterType;
    ReaK::pp::any_sbmp_reporter_chain< ReaK::pp::ptrobot2D_test_world > report_chain;
    report_chain.add_reporter( ReaK::pp::timing_sbmp_report<>(time_ss) );
    report_chain.add_reporter( ReaK::pp::least_cost_sbmp_report<>(cost_ss, &sol_ss) );
    
    ReaK::pp::path_planning_p2p_query< ReaK::pp::ptrobot2D_test_world > mc_query(
      "mc_planning_query",
      world_map,
      world_map->get_start_pos(),
      world_map->get_goal_pos(),
      plan_options.max_results);
    
    std::cout << "Running " << vm["planner-alg"].as< std::string >() << " with " << planner_qualifier_str << ", " << mg_storage_str << ", " << knn_method_str << std::endl;
    timing_output << vm["planner-alg"].as< std::string >() << ", " << planner_qualifier_str << ", " << mg_storage_str << ", " << knn_method_str << std::endl;
    sol_events_output << vm["planner-alg"].as< std::string >() << ", Solutions" << std::endl;
    
#ifdef RK_ENABLE_TEST_RRT_PLANNER
    if(vm["planner-alg"].as< std::string >() == "rrt") {
      ReaK::pp::rrt_planner< ReaK::pp::ptrobot2D_test_world > rrt_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        plan_options.planning_options, 0.1, 0.05, report_chain);
      
      run_monte_carlo_tests(mc_run_count, mc_max_vertices_100, rrt_plan, mc_query, time_ss, cost_ss, sol_ss, timing_output, sol_events_output);
    };
#endif
    
#ifdef RK_ENABLE_TEST_PRM_PLANNER
    if(vm["planner-alg"].as< std::string >() == "prm") {
      ReaK::pp::prm_planner< ReaK::pp::ptrobot2D_test_world > prm_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        0.1, 0.05, world_map->get_max_edge_length(), 2, report_chain);
      
      run_monte_carlo_tests(mc_run_count, mc_max_vertices_100, prm_plan, mc_query, time_ss, cost_ss, sol_ss, timing_output, sol_events_output);
    };
#endif
    
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
    if(vm["planner-alg"].as< std::string >() == "fadprm") {
      ReaK::pp::fadprm_planner< ReaK::pp::ptrobot2D_test_world > fadprm_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        0.1, 0.05, world_map->get_max_edge_length(), 2, report_chain);
      
      fadprm_plan.set_initial_relaxation(plan_options.init_relax);
      
      run_monte_carlo_tests(mc_run_count, mc_max_vertices_100, fadprm_plan, mc_query, time_ss, cost_ss, sol_ss, timing_output, sol_events_output);
    };
#endif
    
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
    if(vm["planner-alg"].as< std::string >() == "sba_star") {
      ReaK::pp::sbastar_planner< ReaK::pp::ptrobot2D_test_world > sbastar_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        plan_options.planning_options, 0.1, 0.05, world_map->get_max_edge_length(), 2, report_chain);
      
      sbastar_plan.set_initial_density_threshold(0.0);
      sbastar_plan.set_initial_relaxation(plan_options.init_relax);
      sbastar_plan.set_initial_SA_temperature(plan_options.init_SA_temp);
      
      run_monte_carlo_tests(mc_run_count, mc_max_vertices_100, sbastar_plan, mc_query, time_ss, cost_ss, sol_ss, timing_output, sol_events_output);
    };
#endif
    
#ifdef RK_ENABLE_TEST_RRTSTAR_PLANNER
    if(vm["planner-alg"].as< std::string >() == "rrt_star") {
      ReaK::pp::rrtstar_planner< ReaK::pp::ptrobot2D_test_world > rrtstar_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        plan_options.planning_options, 0.1, 0.05, 2, report_chain);
      
      run_monte_carlo_tests(mc_run_count, mc_max_vertices_100, rrtstar_plan, mc_query, time_ss, cost_ss, sol_ss, timing_output, sol_events_output);
    };
#endif
    std::cout << "Done!" << std::endl;
    
  };
  
  
  
  
  if(vm.count("single-run")) {
    
    ReaK::pp::path_planning_p2p_query< ReaK::pp::ptrobot2D_test_world > sr_query(
      "sr_planning_query",
      world_map,
      world_map->get_start_pos(),
      world_map->get_goal_pos(),
      plan_options.max_results);
    
    ReaK::pp::differ_sbmp_report_to_space<> image_report("", 0.25 * world_map->get_max_edge_length());
    
    std::cout << "Outputting " << vm["planner-alg"].as< std::string >() << " with " << planner_qualifier_str << ", " << mg_storage_str << ", " << knn_method_str << std::endl;
    
    std::string qualified_output_path = output_path_name + "/" + vm["planner-alg"].as< std::string >() + planner_qualifier_str;
    fs::create_directory(qualified_output_path.c_str());
    
    ReaK::pp::any_sbmp_reporter_chain< ReaK::pp::ptrobot2D_test_world > report_chain;
    image_report.file_path = qualified_output_path + "/" + world_file_name_only + "_";
    report_chain.add_reporter( image_report );
    report_chain.add_reporter( ReaK::pp::print_sbmp_progress<>() );
    
    ReaK::shared_ptr< ReaK::pp::sample_based_planner< ReaK::pp::ptrobot2D_test_world > > p_planner;
    
#ifdef RK_ENABLE_TEST_RRT_PLANNER
    if(vm["planner-alg"].as< std::string >() == "rrt") {
      p_planner = ReaK::shared_ptr< ReaK::pp::sample_based_planner< ReaK::pp::ptrobot2D_test_world > >(
        new ReaK::pp::rrt_planner< ReaK::pp::ptrobot2D_test_world >(
          world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
          plan_options.planning_options, 0.1, 0.05, report_chain)
      );
    };
#endif
    
#ifdef RK_ENABLE_TEST_PRM_PLANNER
    if(vm["planner-alg"].as< std::string >() == "prm") {
      p_planner = ReaK::shared_ptr< ReaK::pp::sample_based_planner< ReaK::pp::ptrobot2D_test_world > >(
        new ReaK::pp::prm_planner< ReaK::pp::ptrobot2D_test_world > prm_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        0.1, 0.05, world_map->get_max_edge_length(), 2, report_chain)
      );
    };
#endif
    
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
    if(vm["planner-alg"].as< std::string >() == "fadprm") {
      ReaK::shared_ptr< ReaK::pp::fadprm_planner< ReaK::pp::ptrobot2D_test_world > > tmp(
        new ReaK::pp::fadprm_planner< ReaK::pp::ptrobot2D_test_world > fadprm_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        0.1, 0.05, world_map->get_max_edge_length(), 2, report_chain)
      );
      tmp->set_initial_relaxation(plan_options.init_relax);
      
      p_planner = tmp;
    };
#endif
    
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
    if(vm["planner-alg"].as< std::string >() == "sba_star") {
      ReaK::shared_ptr< ReaK::pp::sbastar_planner< ReaK::pp::ptrobot2D_test_world > > tmp(
        new ReaK::pp::sbastar_planner< ReaK::pp::ptrobot2D_test_world > sbastar_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        plan_options.planning_options, 0.1, 0.05, world_map->get_max_edge_length(), 2, report_chain)
      );
      tmp->set_initial_density_threshold(0.0);
      tmp->set_initial_relaxation(plan_options.init_relax);
      tmp->set_initial_SA_temperature(plan_options.init_SA_temp);
      
      p_planner = tmp;
    };
#endif
    
#ifdef RK_ENABLE_TEST_RRTSTAR_PLANNER
    if(vm["planner-alg"].as< std::string >() == "rrt_star") {
      p_planner = ReaK::shared_ptr< ReaK::pp::sample_based_planner< ReaK::pp::ptrobot2D_test_world > >(
        new ReaK::pp::rrtstar_planner< ReaK::pp::ptrobot2D_test_world > rrtstar_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        plan_options.planning_options, 0.1, 0.05, 2, report_chain)
      );
    };
#endif
    
    if(!p_planner) {
      std::cout << "Error: Failed to construct a suitable planner! The planner was selected as '" << vm["planner-alg"].as<std::string>() << "'." << std::endl;
      return 4;
    };
    
    sr_query.reset_solution_records();
    p_planner->solve_planning_query(sr_query);
    
    std::cout << "Done!" << std::endl;
    
  };
  
  return 0;
};















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

#include "path_planning/path_planner_options_po.hpp"
#include "path_planning/planning_space_options_po.hpp"
#include "kte_models/chaser_target_model_data_po.hpp"

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







template <typename ManipMdlType, typename InterpTag, int Order>
void CRS_execute_static_planner_impl(const ReaK::kte::chaser_target_data& scene_data, 
                                     const ReaK::pp::planning_option_collection& plan_options,
                                     const ReaK::pp::planning_space_options& space_def,
                                     const ReaK::vect_n<double>& jt_start, 
                                     const ReaK::vect_n<double>& jt_desired) {
  using namespace ReaK;
  using namespace pp;
  
  shared_ptr< ManipMdlType > chaser_concrete_model = rtti::rk_dynamic_ptr_cast<ManipMdlType>(scene_data.chaser_kin_model);
  if( !chaser_concrete_model )
    return;
  
  typedef typename manip_static_workspace< ManipMdlType, Order >::rl_workspace_type static_workspace_type;
  typedef typename manip_pp_traits< ManipMdlType, Order >::rl_jt_space_type rl_jt_space_type;
  typedef typename manip_pp_traits< ManipMdlType, Order >::jt_space_type jt_space_type;
  typedef typename manip_DK_map< ManipMdlType, Order >::rl_map_type rlDK_map_type;
  
  typedef typename topology_traits< rl_jt_space_type >::point_type rl_point_type;
  typedef typename topology_traits< jt_space_type >::point_type point_type;
  
  typedef typename subspace_traits<static_workspace_type>::super_space_type static_super_space_type;  // SuperSpaceType
  
  std::size_t workspace_dims = (Order + 1) * manip_pp_traits< ManipMdlType, Order >::degrees_of_freedom;
  
  // Create the planning spaces:
  shared_ptr< static_workspace_type > workspace = 
    make_manip_static_workspace<Order>(InterpTag(),
      chaser_concrete_model, scene_data.chaser_jt_limits, space_def.min_travel);
  
  shared_ptr< rl_jt_space_type > jt_space = 
    make_manip_rl_jt_space<Order>(chaser_concrete_model, scene_data.chaser_jt_limits);
  
  shared_ptr< jt_space_type > normal_jt_space = 
    make_manip_jt_space<Order>(chaser_concrete_model, scene_data.chaser_jt_limits);
  
  (*workspace) << scene_data.chaser_target_proxy;
  for(std::size_t i = 0; i < scene_data.chaser_env_proxies.size(); ++i)
    (*workspace) << scene_data.chaser_env_proxies[i];
  
  
  // Create the start and goal points:
  rl_point_type start_point, goal_point;
  point_type start_inter, goal_inter;
  
  start_inter = normal_jt_space->origin();
  get<0>(start_inter) = jt_start;
  start_point = scene_data.chaser_jt_limits->map_to_space(start_inter, *normal_jt_space, *jt_space);
  
  goal_inter = normal_jt_space->origin();
  get<0>(goal_inter) = jt_desired;
  goal_point = scene_data.chaser_jt_limits->map_to_space(goal_inter, *normal_jt_space, *jt_space);
  
  
  // Create the reporter chain.
  any_sbmp_reporter_chain< static_workspace_type > report_chain;
  
  report_chain.add_reporter( print_sbmp_progress<>() );
  report_chain.add_reporter( timing_sbmp_report<>() );
  
  
  // Create the point-to-point query:
  path_planning_p2p_query< static_workspace_type > pp_query("pp_query", workspace,
    start_point, goal_point, plan_options.max_results);
  
  
  // Create the planner:
  shared_ptr< sample_based_planner< static_workspace_type > > workspace_planner;
  
#ifdef RK_ENABLE_TEST_RRT_PLANNER
  if( plan_options.planning_algo == 0 ) { // RRT
    
    workspace_planner = shared_ptr< sample_based_planner< static_workspace_type > >(
      new rrt_planner< static_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, report_chain));
    
  } else 
#endif
#ifdef RK_ENABLE_TEST_RRTSTAR_PLANNER
  if( plan_options.planning_algo == 1 ) { // RRT*
    
    workspace_planner = shared_ptr< sample_based_planner< static_workspace_type > >(
      new rrtstar_planner< static_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, workspace_dims, report_chain));
    
  } else 
#endif
#ifdef RK_ENABLE_TEST_PRM_PLANNER
  if( plan_options.planning_algo == 2 ) { // PRM
    
    workspace_planner = shared_ptr< sample_based_planner< static_workspace_type > >(
      new prm_planner< static_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, plan_options.max_random_walk, workspace_dims, report_chain));
    
  } else 
#endif
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
  if( plan_options.planning_algo == 4 ) { // FADPRM
    
    shared_ptr< fadprm_planner< static_workspace_type > > tmp(
      new fadprm_planner< static_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        0.1, 0.05, plan_options.max_random_walk, workspace_dims, report_chain));
    
    tmp->set_initial_relaxation(plan_options.init_relax);
    
    workspace_planner = tmp;
    
  } else 
#endif
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
  if( plan_options.planning_algo == 3 ) { // SBA*
    
    shared_ptr< sbastar_planner< static_workspace_type > > tmp(
      new sbastar_planner< static_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, plan_options.max_random_walk, workspace_dims, report_chain));
    
    tmp->set_initial_density_threshold(0.0);
    tmp->set_initial_relaxation(plan_options.init_relax);
    tmp->set_initial_SA_temperature(plan_options.init_SA_temp);
    
    workspace_planner = tmp;
    
  } else 
#endif
  { };
  
  if(!workspace_planner)
    return;
  
  
  // Solve the planning problem:
  pp_query.reset_solution_records();
  workspace_planner->solve_planning_query(pp_query);
  
  
  // Report the results:
  shared_ptr< seq_path_base< static_super_space_type > > bestsol_rlpath;
  if(pp_query.solutions.size())
    bestsol_rlpath = pp_query.solutions.begin()->second;
  std::cout << "The shortest distance is: " << pp_query.get_best_solution_distance() << std::endl;
  
  
  // Restore model's state:
  chaser_concrete_model->setJointPositions( jt_start );
  chaser_concrete_model->doDirectMotion();
  
};




template <typename Topology>
void run_monte_carlo_tests(
    std::size_t mc_run_count,
    std::size_t mc_num_records,
    ReaK::pp::sample_based_planner< Topology >& planner,
    ReaK::pp::planning_query< Topology >& mc_query,
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
 - space configuration : specifies the order, interpolator, min/max travel, temporality, 
                         rate-limited'ness, and output space-order.
 - start joint position  (could be embedded in the chaser model)
 - target pose           (could be embedded in the target model)

*****************************************************************************************************/



int main(int argc, char** argv) {
  
  using namespace ReaK;
  using namespace pp;
  using namespace kte;
  
  std::string config_file;
  
  po::options_description generic_options("Generic options");
  generic_options.add_options()
    ("help,h", "produce this help message.")
    ("config,c", po::value< std::string >(&config_file)->default_value("run_CRS_planner.cfg"),
                  "configuration file-name (can contain any or all options, will be overriden by command-line options).")
  ;
  
  po::options_description io_options("I/O options");
  io_options.add_options()
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
  
  po::options_description planner_select_options = get_planning_option_po_desc();
  planner_select_options.add_options()
    ("planner-alg", po::value< std::string >(), "specify the planner algorithm to use, can be any of (" + available_algs + ").");
  
  po::options_description space_def_options = get_planning_space_options_po_desc();
  
  po::options_description scene_data_options = get_chaser_target_data_po_desc();
  
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
  cmdline_options.add(generic_options).add(io_options).add(mc_options).add(single_options)
                 .add(planner_select_options).add(scene_data_options).add(space_def_options).add(generate_options);
  
  po::options_description config_file_options;
  config_file_options.add(io_options).add(mc_options).add(single_options)
                     .add(planner_select_options).add(scene_data_options).add(space_def_options);
  
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
  
  std::string output_path_name = vm["output-path"].as<std::string>();
  while(output_path_name[output_path_name.length()-1] == '/') 
    output_path_name.erase(output_path_name.length()-1, 1);
  
  fs::create_directory(output_path_name.c_str());
  
  
  planning_option_collection plan_options = get_planning_option_from_po(vm);
  
  std::string knn_method_str = plan_options.get_knn_method_str();
  std::string mg_storage_str = plan_options.get_mg_storage_str();
  std::string planner_qualifier_str = plan_options.get_planner_qualifier_str();
  
  
  planning_space_options space_def = get_planning_space_options_from_po(vm);
  
  chaser_target_data scene_data = get_chaser_target_data_from_po(vm);
  
  vect_n<double> jt_start = scene_data.chaser_kin_model->getJointPositions(); 
  if( vm.count("start-configuration") ) {
    try {
      vect_n<double> jt_start_tmp = jt_start; 
      (*serialization::open_iarchive(vm["start-configuration"].as< std::string >()))
        >> jt_start_tmp;
      jt_start = jt_start_tmp;
    } catch(std::exception& e) { 
      std::cerr << "Error: Could not load the start-configuration file!" << std::endl;
    };
  };
  
  frame_3D<double> target_frame = scene_data.target_frame->getGlobalFrame();
  if( vm.count("target-pose") ) {
    try {
      frame_3D<double> target_frame_tmp = target_frame;
      (*serialization::open_iarchive(vm["target-pose"].as< std::string >()))
        >> target_frame_tmp;
      target_frame = target_frame_tmp;
    } catch(std::exception& e) { 
      std::cerr << "Error: Could not load the target-pose file!" << std::endl;
    };
  };
  
  shared_ptr< frame_3D<double> > dep_EE_frame = scene_data.chaser_kin_model->getDependentFrame3D(0)->mFrame;
  
  vect_n<double> jt_desired(7,0.0);
  try {
    *dep_EE_frame = target_frame;
    scene_data.chaser_kin_model->doInverseMotion();
    jt_desired = scene_data.chaser_kin_model->getJointPositions();
  } catch( optim::infeasible_problem& e ) { RK_UNUSED(e);
    std::cerr << "Error: The target frame cannot be reached! No inverse kinematics solution possible!" << std::endl;
    return 10;
  };
  scene_data.chaser_kin_model->setJointPositions(jt_start);
  scene_data.chaser_kin_model->doDirectMotion();
  
  
  
  // Do the generations if required:
  
  if( vm.count("generate-all-files") + vm.count("generate-planner-options") > 0 ) {
    std::string file_name;
    if( vm.count("generate-planner-options") == 0 ) {
      file_name = vm["generate-all-files"].as< std::string >() + "_planner";
    } else {
      file_name = vm["generate-planner-options"].as< std::string >();
    };
    if( vm.count("generate-protobuf") ) 
      file_name += ".pbuf";
    else if( vm.count("generate-binary") )
      file_name += ".rkb";
    else 
      file_name += ".rkx";
    
    try {
      (*serialization::open_oarchive(file_name)) << plan_options;
    } catch( std::exception& e ) { 
      std::cerr << "Error: Could not generate the planner options file!" << std::endl;
    };
  };
  
  if( vm.count("generate-all-files") + vm.count("generate-chaser-target-env") > 0 ) {
    std::string file_name;
    if( vm.count("generate-chaser-target-env") == 0 ) {
      file_name = vm["generate-all-files"].as< std::string >() + "_models";
    } else {
      file_name = vm["generate-chaser-target-env"].as< std::string >();
    };
    if( vm.count("generate-protobuf") ) 
      file_name += ".pbuf";
    else if( vm.count("generate-binary") )
      file_name += ".rkb";
    else 
      file_name += ".rkx";
    
    try {
      (*serialization::open_oarchive(file_name)) << scene_data;
    } catch( std::exception& e ) { 
      std::cerr << "Error: Could not generate the chaser-target-env model file!" << std::endl;
    };
  };
  
  if( vm.count("generate-all-files") + vm.count("generate-space-definition") > 0 ) {
    std::string file_name;
    if( vm.count("generate-space-definition") == 0 ) {
      file_name = vm["generate-all-files"].as< std::string >() + "_space";
    } else {
      file_name = vm["generate-space-definition"].as< std::string >();
    };
    if( vm.count("generate-protobuf") ) 
      file_name += ".pbuf";
    else if( vm.count("generate-binary") )
      file_name += ".rkb";
    else 
      file_name += ".rkx";
    
    try {
      (*serialization::open_oarchive(file_name)) << space_def;
    } catch( std::exception& e ) { 
      std::cerr << "Error: Could not generate the space-definition file!" << std::endl;
    };
  };
  
  if( vm.count("generate-all-files") + vm.count("generate-start-config") > 0 ) {
    std::string file_name;
    if( vm.count("generate-start-config") == 0 ) {
      file_name = vm["generate-all-files"].as< std::string >() + "_start_config";
    } else {
      file_name = vm["generate-start-config"].as< std::string >();
    };
    if( vm.count("generate-protobuf") ) 
      file_name += ".pbuf";
    else if( vm.count("generate-binary") )
      file_name += ".rkb";
    else 
      file_name += ".rkx";
    
    try {
      (*serialization::open_oarchive(file_name)) << jt_start;
    } catch( std::exception& e ) { 
      std::cerr << "Error: Could not generate the start-configuration file!" << std::endl;
    };
  };
  
  if( vm.count("generate-all-files") + vm.count("generate-target-pose") > 0 ) {
    std::string file_name;
    if( vm.count("generate-target-pose") == 0 ) {
      file_name = vm["generate-all-files"].as< std::string >() + "_target_pose";
    } else {
      file_name = vm["generate-target-pose"].as< std::string >();
    };
    if( vm.count("generate-protobuf") ) 
      file_name += ".pbuf";
    else if( vm.count("generate-binary") )
      file_name += ".rkb";
    else 
      file_name += ".rkx";
    
    try {
      (*serialization::open_oarchive(file_name)) << target_frame;
    } catch( std::exception& e ) { 
      std::cerr << "Error: Could not generate the target-pose file!" << std::endl;
    };
  };
  
  
  
  shared_ptr< ptrobot2D_test_world > world_map =
    shared_ptr< ptrobot2D_test_world >(new ptrobot2D_test_world(world_file_name, 20, 1.0));
  
  path_planning_p2p_query< ptrobot2D_test_world > pp_query(
    "planning_query",
    world_map,
    world_map->get_start_pos(),
    world_map->get_goal_pos(),
    plan_options.max_results);
  
  any_sbmp_reporter_chain< ptrobot2D_test_world > report_chain;
  
  if(vm.count("monte-carlo")) {
    
    std::size_t mc_run_count        = vm["mc-runs"].as<std::size_t>();
    std::size_t mc_max_vertices_100 = plan_options.max_vertices / plan_options.prog_interval;
    
    std::ofstream timing_output(output_path_name + "/" + world_file_name_only + "_times.txt");
    std::ofstream sol_events_output(output_path_name + "/" + world_file_name_only + "_solutions.txt");
    
    std::stringstream time_ss, cost_ss, sol_ss;
    
    report_chain.add_reporter( timing_sbmp_report<>(time_ss) );
    report_chain.add_reporter( least_cost_sbmp_report<>(cost_ss, &sol_ss) );
    
    std::cout << "Running " << vm["planner-alg"].as< std::string >() << " with " << planner_qualifier_str << ", " << mg_storage_str << ", " << knn_method_str << std::endl;
    timing_output << vm["planner-alg"].as< std::string >() << ", " << planner_qualifier_str << ", " << mg_storage_str << ", " << knn_method_str << std::endl;
    sol_events_output << vm["planner-alg"].as< std::string >() << ", Solutions" << std::endl;
    
    shared_ptr< sample_based_planner< ptrobot2D_test_world > > p_planner;
    
#ifdef RK_ENABLE_TEST_RRT_PLANNER
    if(vm["planner-alg"].as< std::string >() == "rrt") {
      p_planner = shared_ptr< sample_based_planner< ptrobot2D_test_world > >(
        new rrt_planner< ptrobot2D_test_world >(
          world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
          plan_options.planning_options, 0.1, 0.05, report_chain)
      );
    };
#endif
    
#ifdef RK_ENABLE_TEST_PRM_PLANNER
    if(vm["planner-alg"].as< std::string >() == "prm") {
      p_planner = shared_ptr< sample_based_planner< ptrobot2D_test_world > >(
        new prm_planner< ptrobot2D_test_world > prm_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        0.1, 0.05, world_map->get_max_edge_length(), 2, report_chain)
      );
    };
#endif
    
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
    if(vm["planner-alg"].as< std::string >() == "fadprm") {
      shared_ptr< fadprm_planner< ptrobot2D_test_world > > tmp(
        new fadprm_planner< ptrobot2D_test_world > fadprm_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        0.1, 0.05, world_map->get_max_edge_length(), 2, report_chain)
      );
      tmp->set_initial_relaxation(plan_options.init_relax);
      
      p_planner = tmp;
    };
#endif
    
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
    if(vm["planner-alg"].as< std::string >() == "sba_star") {
      shared_ptr< sbastar_planner< ptrobot2D_test_world > > tmp(
        new sbastar_planner< ptrobot2D_test_world > sbastar_plan(
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
      p_planner = shared_ptr< sample_based_planner< ptrobot2D_test_world > >(
        new rrtstar_planner< ptrobot2D_test_world > rrtstar_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        plan_options.planning_options, 0.1, 0.05, 2, report_chain)
      );
    };
#endif
    
    if(!p_planner) {
      std::cout << "Error: Failed to construct a suitable planner! The planner was selected as '" << vm["planner-alg"].as<std::string>() << "'." << std::endl;
      return 4;
    };
    
    run_monte_carlo_tests(mc_run_count, mc_max_vertices_100, *p_planner, pp_query, time_ss, cost_ss, sol_ss, timing_output, sol_events_output);
    
    std::cout << "Done!" << std::endl;
    
  };
  
  
  
  
  if(vm.count("single-run")) {
    
    
    std::cout << "Outputting " << vm["planner-alg"].as< std::string >() << " with " << planner_qualifier_str << ", " << mg_storage_str << ", " << knn_method_str << std::endl;
    
    std::string qualified_output_path = output_path_name + "/" + vm["planner-alg"].as< std::string >() + planner_qualifier_str;
    fs::create_directory(qualified_output_path.c_str());
    
    differ_sbmp_report_to_space<> image_report("", 0.25 * world_map->get_max_edge_length());
    image_report.file_path = qualified_output_path + "/" + world_file_name_only + "_";
    report_chain.add_reporter( image_report );
    report_chain.add_reporter( print_sbmp_progress<>() );
    
    shared_ptr< sample_based_planner< ptrobot2D_test_world > > p_planner;
    
#ifdef RK_ENABLE_TEST_RRT_PLANNER
    if(vm["planner-alg"].as< std::string >() == "rrt") {
      p_planner = shared_ptr< sample_based_planner< ptrobot2D_test_world > >(
        new rrt_planner< ptrobot2D_test_world >(
          world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
          plan_options.planning_options, 0.1, 0.05, report_chain)
      );
    };
#endif
    
#ifdef RK_ENABLE_TEST_PRM_PLANNER
    if(vm["planner-alg"].as< std::string >() == "prm") {
      p_planner = shared_ptr< sample_based_planner< ptrobot2D_test_world > >(
        new prm_planner< ptrobot2D_test_world > prm_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        0.1, 0.05, world_map->get_max_edge_length(), 2, report_chain)
      );
    };
#endif
    
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
    if(vm["planner-alg"].as< std::string >() == "fadprm") {
      shared_ptr< fadprm_planner< ptrobot2D_test_world > > tmp(
        new fadprm_planner< ptrobot2D_test_world > fadprm_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        0.1, 0.05, world_map->get_max_edge_length(), 2, report_chain)
      );
      tmp->set_initial_relaxation(plan_options.init_relax);
      
      p_planner = tmp;
    };
#endif
    
#ifdef RK_ENABLE_TEST_SBASTAR_PLANNER
    if(vm["planner-alg"].as< std::string >() == "sba_star") {
      shared_ptr< sbastar_planner< ptrobot2D_test_world > > tmp(
        new sbastar_planner< ptrobot2D_test_world > sbastar_plan(
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
      p_planner = shared_ptr< sample_based_planner< ptrobot2D_test_world > >(
        new rrtstar_planner< ptrobot2D_test_world > rrtstar_plan(
        world_map, plan_options.max_vertices, plan_options.prog_interval, plan_options.store_policy | plan_options.knn_method, 
        plan_options.planning_options, 0.1, 0.05, 2, report_chain)
      );
    };
#endif
    
    if(!p_planner) {
      std::cout << "Error: Failed to construct a suitable planner! The planner was selected as '" << vm["planner-alg"].as<std::string>() << "'." << std::endl;
      return 4;
    };
    
    pp_query.reset_solution_records();
    p_planner->solve_planning_query(pp_query);
    
    std::cout << "Done!" << std::endl;
    
  };
  
  return 0;
};















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

#include "mbd_kte/inertia.hpp"
#include "mbd_kte/driving_actuator.hpp"
#include "mbd_kte/state_measures.hpp"
#include "mbd_kte/free_joints.hpp"
#include "mbd_kte/kte_map_chain.hpp"

#include "kte_models/manip_dynamics_model.hpp"

#include "shapes/colored_model.hpp"
#include "shapes/sphere.hpp"
#include "shapes/box.hpp"
#include "shapes/coord_arrows_3D.hpp"
#include "proximity/proxy_query_model.hpp"

#include "topologies/manip_planning_traits.hpp"
#include "topologies/manip_P3R3R_workspaces.hpp"
#include "topologies/Ndof_linear_spaces.hpp"
#include "topologies/Ndof_cubic_spaces.hpp"
#include "topologies/Ndof_quintic_spaces.hpp"
#include "topologies/Ndof_svp_spaces.hpp"
#include "topologies/Ndof_sap_spaces.hpp"


#include "path_planning/basic_sbmp_reporters.hpp"
#include "path_planning/vlist_sbmp_report.hpp"



#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;


namespace {

struct monte_carlo_engine {
  
  std::size_t mc_run_count;
  
  std::ofstream timing_output;
  std::ofstream sol_events_output;
  
  std::stringstream time_ss;
  std::stringstream cost_ss;
  std::stringstream sol_ss;
  
  monte_carlo_engine(std::size_t aMCRuns, 
                     const std::string& aPlannerName,
                     const std::string& aOutputPathStem) : 
                     mc_run_count(aMCRuns) {
    
    fs::create_directory(aOutputPathStem.c_str());
    
    timing_output.open(aOutputPathStem + "/" + aPlannerName + "_times.txt");
    sol_events_output.open(aOutputPathStem + "/" + aPlannerName + "_solutions.txt");
    
    timing_output << aPlannerName << std::endl;
    sol_events_output << aPlannerName << ", Solutions" << std::endl;
    std::cout << "Running " << aPlannerName << std::endl;
  };
  
  template <typename Topology>
  ReaK::shared_ptr< ReaK::pp::any_sbmp_reporter_chain< Topology > > create_reporter(ReaK::shared_ptr< Topology >) {
    
    ReaK::shared_ptr< ReaK::pp::any_sbmp_reporter_chain< Topology > > report_chain(
      new ReaK::pp::any_sbmp_reporter_chain< Topology >());
    
    report_chain->add_reporter( ReaK::pp::timing_sbmp_report<>(time_ss) );
    report_chain->add_reporter( ReaK::pp::least_cost_sbmp_report<>(cost_ss, &sol_ss) );
    
    return report_chain;
  };
  
  template <typename Topology>
  void operator()(const ReaK::pp::planning_option_collection& plan_options,
                  ReaK::shared_ptr< ReaK::pp::sample_based_planner< Topology > > planner,
                  ReaK::pp::planning_query< Topology >& mc_query) {
    std::size_t mc_num_records = plan_options.max_vertices / plan_options.prog_interval;
    
    std::vector< double > vertex_counts(mc_num_records, 0.0);
    std::vector< std::size_t > num_remaining_planners(mc_num_records, 0);
    std::vector< std::size_t > num_successful_planners(mc_num_records, 0);
    std::vector< double > time_values(mc_num_records, 0.0);
    std::vector< double > best_costs(mc_num_records, 1.0e10);
    std::vector< double > worst_costs(mc_num_records, 0.0);
    std::vector< double > avg_costs(mc_num_records, 0.0);
    
    cost_ss << std::fixed;
    sol_ss << std::fixed;
    
    for(std::size_t i = 0; i < mc_run_count; ++i) {
      time_ss.clear();
      cost_ss.clear();
      sol_ss.clear();
      
      mc_query.reset_solution_records();
      planner->reset_internal_state();
      planner->solve_planning_query(mc_query);
      
      std::size_t v_count = 0, t_val = 0; 
      std::string tmp;
      std::size_t j = 0;
      while( std::getline(time_ss, tmp) && (tmp.size()) ) {
        std::stringstream ss_tmp(tmp);
        ss_tmp >> v_count >> t_val;
        vertex_counts[j] = (double(v_count) + double(num_remaining_planners[j]) * vertex_counts[j]) / double(num_remaining_planners[j] + 1);
        time_values[j] = (double(t_val) + double(num_remaining_planners[j]) * time_values[j]) / double(num_remaining_planners[j] + 1);
        num_remaining_planners[j] += 1; 
        ++j;
      };
      
      double c_val = 1e10;
      j = 0;
      while( std::getline(cost_ss, tmp) && (tmp.size()) ) {
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
      std::getline(sol_ss, first_sol_event);
      if(first_sol_event != "")
        sol_events_output << first_sol_event << std::endl;
    };
    for(std::size_t i = 0; i < mc_num_records; ++i) {
      timing_output << std::setw(9) << i 
            << " " << std::setw(9) << vertex_counts[i] 
            << " " << std::setw(9) << num_remaining_planners[i] 
            << " " << std::setw(9) << num_successful_planners[i] 
            << " " << std::setw(9) << time_values[i] 
            << " " << std::setw(9) << best_costs[i] 
            << " " << std::setw(9) << worst_costs[i] 
            << " " << std::setw(9) << avg_costs[i] << std::endl; 
    };
    
    std::cout << "Done!" << std::endl;
    
    
  };
  
  
};




struct single_run_engine {
  
  single_run_engine(const std::string& aPlannerName,
                    const std::string& aOutputPathStem) {
    
    fs::create_directory(aOutputPathStem.c_str());
    
    std::string qualified_output_path = aOutputPathStem + "/" + aPlannerName;
    fs::create_directory(qualified_output_path.c_str());
    
    // TODO: open some output file to output the qualified_output_path folder.
    
    std::cout << "Outputting " << aPlannerName << std::endl;
    
  };
  
  template <typename Topology>
  ReaK::shared_ptr< ReaK::pp::any_sbmp_reporter_chain< Topology > > create_reporter(ReaK::shared_ptr< Topology >) {
    
    ReaK::shared_ptr< ReaK::pp::any_sbmp_reporter_chain< Topology > > report_chain(
      new ReaK::pp::any_sbmp_reporter_chain< Topology >());
    
    report_chain->add_reporter( ReaK::pp::print_sbmp_progress<>() );
    
    return report_chain;
  };
  
  template <typename Topology>
  void operator()(const ReaK::pp::planning_option_collection& plan_options,
                  ReaK::shared_ptr< ReaK::pp::sample_based_planner< Topology > > planner,
                  ReaK::pp::planning_query< Topology >& pp_query) {
    
    pp_query.reset_solution_records();
    planner->reset_internal_state();
    planner->solve_planning_query(pp_query);
    
    std::cout << "Done!" << std::endl;
    
    // Report the results:
    // TODO: output something... (I don't know exactly, maybe vlist printer).
//     shared_ptr< seq_path_base< static_super_space_type > > bestsol_rlpath;
//     if(pp_query.solutions.size())
//       bestsol_rlpath = pp_query.solutions.begin()->second;
    
    std::cout << "The shortest distance is: " << pp_query.get_best_solution_distance() << std::endl;
    
    
  };
  
  
};



};




template <typename ManipMdlType, typename InterpTag, int Order, typename PlanEngine>
void CRS_execute_static_planner(const ReaK::kte::chaser_target_data& scene_data, 
                                const ReaK::pp::planning_option_collection& plan_options,
                                const ReaK::pp::planning_space_options& space_def,
                                PlanEngine& engine,
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
  
  typedef typename topology_traits< rl_jt_space_type >::point_type rl_point_type;
  typedef typename topology_traits< jt_space_type >::point_type point_type;
  
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
  shared_ptr< any_sbmp_reporter_chain< static_workspace_type > > p_report_chain = engine.create_reporter(workspace);
  
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
        0.1, 0.05, *p_report_chain));
    
  } else 
#endif
#ifdef RK_ENABLE_TEST_RRTSTAR_PLANNER
  if( plan_options.planning_algo == 1 ) { // RRT*
    
    workspace_planner = shared_ptr< sample_based_planner< static_workspace_type > >(
      new rrtstar_planner< static_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, workspace_dims, *p_report_chain));
    
  } else 
#endif
#ifdef RK_ENABLE_TEST_PRM_PLANNER
  if( plan_options.planning_algo == 2 ) { // PRM
    
    workspace_planner = shared_ptr< sample_based_planner< static_workspace_type > >(
      new prm_planner< static_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, plan_options.max_random_walk, workspace_dims, *p_report_chain));
    
  } else 
#endif
#ifdef RK_ENABLE_TEST_FADPRM_PLANNER
  if( plan_options.planning_algo == 4 ) { // FADPRM
    
    shared_ptr< fadprm_planner< static_workspace_type > > tmp(
      new fadprm_planner< static_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        0.1, 0.05, plan_options.max_random_walk, workspace_dims, *p_report_chain));
    
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
        0.1, 0.05, plan_options.max_random_walk, workspace_dims, *p_report_chain));
    
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
  engine(plan_options, workspace_planner, pp_query);
  
  
  // Restore model's state:
  chaser_concrete_model->setJointPositions( jt_start );
  chaser_concrete_model->doDirectMotion();
  
};




template <typename PlanEngine>
void CRS_launch_static_planner(const ReaK::kte::chaser_target_data& scene_data, 
                               const ReaK::pp::planning_option_collection& plan_options,
                               const ReaK::pp::planning_space_options& space_def,
                               PlanEngine& engine,
                               const ReaK::vect_n<double>& jt_start, 
                               const ReaK::vect_n<double>& jt_desired) {
  
  try {
    if((space_def.get_space_order() == 0) && (space_def.get_interp_id() == 0)) { 
      CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics, ReaK::pp::linear_interpolation_tag, 0>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
    } else 
#if 0
    if((space_def.get_space_order() == 1) && (space_def.get_interp_id() == 0)) {
      CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics, ReaK::pp::linear_interpolation_tag, 1>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
    } else 
    if((space_def.get_space_order() == 2) && (space_def.get_interp_id() == 0)) {
      CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics, ReaK::pp::linear_interpolation_tag, 2>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
    } else 
#endif
    if((space_def.get_space_order() == 1) && (space_def.get_interp_id() == 1)) {
      CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics, ReaK::pp::cubic_hermite_interpolation_tag, 1>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
    } else 
#if 0
    if((space_def.get_space_order() == 2) && (space_def.get_interp_id() == 1)) {
      CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics, ReaK::pp::cubic_hermite_interpolation_tag, 2>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
    } else 
#endif
    if((space_def.get_space_order() == 2) && (space_def.get_interp_id() == 2)) {
      CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics, ReaK::pp::quintic_hermite_interpolation_tag, 2>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
    } else 
    if((space_def.get_space_order() == 1) && (space_def.get_interp_id() == 3)) {
      CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics, ReaK::pp::svp_Ndof_interpolation_tag, 1>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
    } else 
#if 0
    if((space_def.get_space_order() == 2) && (space_def.get_interp_id() == 3)) {
      CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics, ReaK::pp::svp_Ndof_interpolation_tag, 2>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
    } else 
#endif
    if((space_def.get_space_order() == 2) && (space_def.get_interp_id() == 4)) {
      CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics, ReaK::pp::sap_Ndof_interpolation_tag, 2>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
    };
  } catch(std::exception& e) {
    std::cerr << "Error: An exception was raised during the planning:\nwhat(): " << e.what() << std::endl;
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
    
    ("output-path,o", po::value< std::string >()->default_value("pp_results"), "specify the output path (default is pp_results).")
    ("result-file-prefix", po::value< std::string >(), "specify the prefix to apply to the result output files.")
  ;
  
  po::options_description mc_options("Monte-Carlo options");
  mc_options.add_options()
    ("monte-carlo,m", "specify that monte-carlo runs should be performed (default is not).")
    ("mc-runs", po::value< std::size_t >()->default_value(100), "number of monte-carlo runs to average out (default is 100).")
  ;
  
  po::options_description single_options("Single-run options");
  single_options.add_options()
    ("single-run,s", "specify that single runs should be performed (default is not).")
  ;
  
  po::options_description planner_select_options = get_planning_option_po_desc();
  
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
  
  std::string output_path_name = vm["output-path"].as<std::string>();
  while(output_path_name[output_path_name.length()-1] == '/') 
    output_path_name.erase(output_path_name.length()-1, 1);
  
  fs::create_directory(output_path_name.c_str());
  
  std::string result_file_prefix = "CRS_static_scene";
  if( vm.count("result-file-prefix") )
    result_file_prefix = vm["result-file-prefix"].as<std::string>();
  
  
  
  planning_option_collection plan_options = get_planning_option_from_po(vm);
  
  std::string knn_method_str = plan_options.get_knn_method_str();
  std::string mg_storage_str = plan_options.get_mg_storage_str();
  std::string planner_qualifier_str = plan_options.get_planner_qualifier_str();
  std::string planner_name_str = plan_options.get_planning_algo_str() + "_" + planner_qualifier_str + "_" + mg_storage_str + "_" + knn_method_str;
  
  
  planning_space_options space_def = get_planning_space_options_from_po(vm);
  
  chaser_target_data scene_data = get_chaser_target_data_from_po(vm);
  
  vect_n<double> jt_start(7,0.0);
  if( scene_data.chaser_kin_model ) {
    jt_start = scene_data.chaser_kin_model->getJointPositions(); 
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
  };
  
  frame_3D<double> target_frame;
  if( scene_data.target_kin_model ) {
    target_frame = scene_data.target_kin_model->getFrame3D(0)->getGlobalFrame();
    if( vm.count("target-pose") ) {
      try {
        frame_3D<double> target_frame_tmp = target_frame;
        (*serialization::open_iarchive(vm["target-pose"].as< std::string >()))
          >> target_frame_tmp;
        target_frame = target_frame_tmp;
        *(scene_data.target_kin_model->getFrame3D(0)) = target_frame;
        scene_data.target_kin_model->doDirectMotion();
      } catch(std::exception& e) { 
        std::cerr << "Error: Could not load the target-pose file!" << std::endl;
      };
    };
  };
  
  
  vect_n<double> jt_desired(7,0.0);
  if( scene_data.chaser_kin_model ) {
    shared_ptr< frame_3D<double> > dep_EE_frame = scene_data.chaser_kin_model->getDependentFrame3D(0)->mFrame;
    if( vm.count("monte-carlo") + vm.count("single-run") > 0 ) {
      try {
        *dep_EE_frame = scene_data.target_frame->getGlobalFrame();
        scene_data.chaser_kin_model->doInverseMotion();
        jt_desired = scene_data.chaser_kin_model->getJointPositions();
      } catch( optim::infeasible_problem& e ) { RK_UNUSED(e);
        std::cerr << "Error: The target frame cannot be reached! No inverse kinematics solution possible!" << std::endl;
        return 10;
      };
      scene_data.chaser_kin_model->setJointPositions(jt_start);
      scene_data.chaser_kin_model->doDirectMotion();
    };
  };
  
  
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
  
  
  // Do the Monte-Carlo runs if required:
  
  if(vm.count("monte-carlo")) {
    
    monte_carlo_engine mc_eng(vm["mc-runs"].as<std::size_t>(), planner_name_str, output_path_name + "/" + result_file_prefix);
    
    CRS_launch_static_planner(scene_data, plan_options, space_def, mc_eng, jt_start, jt_desired);
    
  };
  
  
  // Do a single run if required:
  
  if(vm.count("single-run")) {
    
    single_run_engine sr_eng(planner_name_str, output_path_name + "/" + result_file_prefix);
    
    CRS_launch_static_planner(scene_data, plan_options, space_def, sr_eng, jt_start, jt_desired);
    
  };
  
  return 0;
};














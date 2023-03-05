
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

#include <fstream>
#include <iostream>

#define RK_DISABLE_RRT_PLANNER
// #define RK_DISABLE_RRTSTAR_PLANNER
#define RK_DISABLE_PRM_PLANNER
#define RK_DISABLE_FADPRM_PLANNER
// #define RK_DISABLE_SBASTAR_PLANNER

// disable template definitions for the planners, because this program uses extern templates for the planners:
#define RK_DISABLE_PLANNER_DEFINITIONS

#include "ReaK/planning/path_planning/planner_exec_engines.hpp"
#include "ReaK/planning/path_planning/planner_exec_intercept.hpp"

#include "ReaK/mbd/models/chaser_target_model_data_po.hpp"
#include "ReaK/planning/path_planning/path_planner_options_po.hpp"
#include "ReaK/planning/path_planning/planning_space_options_po.hpp"

#include "ReaK/math/optimization/optim_exceptions.hpp"

#include "ReaK/mbd/kte/driving_actuator.hpp"
#include "ReaK/mbd/kte/free_joints.hpp"
#include "ReaK/mbd/kte/inertia.hpp"
#include "ReaK/mbd/kte/kte_map_chain.hpp"
#include "ReaK/mbd/kte/state_measures.hpp"

#include "ReaK/mbd/models/manip_dynamics_model.hpp"

#include "ReaK/geometry/proximity/proxy_query_model.hpp"
#include "ReaK/geometry/shapes/box.hpp"
#include "ReaK/geometry/shapes/colored_model.hpp"
#include "ReaK/geometry/shapes/coord_arrows_3D.hpp"
#include "ReaK/geometry/shapes/sphere.hpp"

#include "ReaK/topologies/interpolation/discrete_point_trajectory.hpp"
#include "ReaK/topologies/interpolation/trajectory_base.hpp"
#include "ReaK/topologies/spaces/manip_P3R3R_workspaces.hpp"
#include "ReaK/topologies/spaces/manip_planning_traits.hpp"

#include "ReaK/topologies/interpolation/transformed_trajectory.hpp"
#include "ReaK/topologies/spaces/direct_inverse_kin_topomap.hpp"
#include "ReaK/topologies/spaces/proxy_traj_applicator.hpp"

#include "ReaK/topologies/spaces/Ndof_cubic_spaces.hpp"
#include "ReaK/topologies/spaces/Ndof_linear_spaces.hpp"
#include "ReaK/topologies/spaces/Ndof_quintic_spaces.hpp"
#include "ReaK/topologies/spaces/Ndof_sap_spaces.hpp"
#include "ReaK/topologies/spaces/Ndof_svp_spaces.hpp"

#include <filesystem>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

// I/O options
ABSL_FLAG(std::string, start_configuration, "",
          "Specify the file containing the start configuration of the chaser "
          "(P3R3R-manipulator), if not specified, the chaser configuration of "
          "the chaser model will be used.");
ABSL_FLAG(std::string, target_pose, "",
          "Specify the file containing the target pose of the capture target "
          "(satellite / airship), if not specified, the pose of the target "
          "model will be used.");
ABSL_FLAG(std::string, target_trajectory, "",
          "Specify the file containing the target trajectory of the capture "
          "target (satellite / airship), if not specified, the pose of the "
          "target model will be used.");
ABSL_FLAG(std::string, output_path, "pp_results",
          "Specify the output path (default is pp_results).");
ABSL_FLAG(std::string, result_file_prefix, "",
          "Specify the prefix to apply to the result output files.");

// Monte-Carlo options
ABSL_FLAG(
    bool, monte_carlo, false,
    "Specify that monte-carlo runs should be performed (default is not).");
ABSL_FLAG(int, mc_runs, 100,
          "Number of monte-carlo runs to average out (default is 100).");

// Single-run options
ABSL_FLAG(bool, single_run, false,
          "Specify that single runs should be performed (default is not).");

// File generation options
ABSL_FLAG(std::string, generate_all_files, "",
          "Specify that all configuration files should be generated with the "
          "given file-name prefix (file-name without suffix and extension).");
ABSL_FLAG(std::string, generate_planner_options, "",
          "Specify that the planner options file should be generated with the "
          "given file-name prefix (file-name without extension).");
ABSL_FLAG(std::string, generate_chaser_target_env, "",
          "Specify that the chaser-target-env file should be generated with "
          "the given file-name prefix (file-name without extension).");
ABSL_FLAG(std::string, generate_space_definition, "",
          "Specify that the space-definition file should be generated with the "
          "given file-name prefix (file-name without extension).");
ABSL_FLAG(std::string, generate_start_config, "",
          "Specify that the start-configuration file should be generated with "
          "the given file-name prefix (file-name without extension).");
ABSL_FLAG(std::string, generate_target_pose, "",
          "Specify that the target-pose file should be generated with the "
          "given file-name prefix (file-name without extension).");
ABSL_FLAG(bool, generate_xml, false,
          "If set, output results in XML format (rkx) (default).");
ABSL_FLAG(bool, generate_protobuf, false,
          "If set, output results in protobuf format (pb).");
ABSL_FLAG(bool, generate_binary, false,
          "If set, output results in binary format (rkb).");

namespace fs = std::filesystem;

using sat_state_space_type = ReaK::pp::se3_1st_order_topology<double>::type;
using sat_temporal_space_type =
    ReaK::pp::temporal_space<sat_state_space_type,
                             ReaK::pp::time_poisson_topology,
                             ReaK::pp::time_distance_only>;
using sat_trajectory_type = ReaK::pp::trajectory_base<sat_temporal_space_type>;

template <typename ManipMdlType, typename InterpTag, int Order,
          typename PlanEngine>
void CRS_execute_static_planner(
    const ReaK::kte::chaser_target_data& scene_data,
    const ReaK::pp::planning_option_collection& plan_options,
    const ReaK::pp::planning_space_options& space_def, PlanEngine& engine,
    const ReaK::vect_n<double>& jt_start,
    const ReaK::vect_n<double>& jt_desired) {
  using namespace ReaK;
  using namespace pp;

  std::shared_ptr<ManipMdlType> chaser_concrete_model =
      rtti::rk_dynamic_ptr_cast<ManipMdlType>(scene_data.chaser_kin_model);
  if (!chaser_concrete_model) {
    return;
  }

  using static_workspace_type =
      typename manip_static_workspace<ManipMdlType, Order>::rl_workspace_type;
  using rl_jt_space_type =
      typename manip_pp_traits<ManipMdlType, Order>::rl_jt_space_type;
  using jt_space_type =
      typename manip_pp_traits<ManipMdlType, Order>::jt_space_type;

  using rl_point_type = typename topology_traits<rl_jt_space_type>::point_type;
  using point_type = typename topology_traits<jt_space_type>::point_type;

  std::size_t workspace_dims =
      (Order + 1) * manip_pp_traits<ManipMdlType, Order>::degrees_of_freedom;

  // Create the planning spaces:
  std::shared_ptr<static_workspace_type> workspace =
      make_manip_static_workspace<Order>(InterpTag(), chaser_concrete_model,
                                         scene_data.chaser_jt_limits,
                                         space_def.min_travel);

  std::shared_ptr<rl_jt_space_type> jt_space = make_manip_rl_jt_space<Order>(
      chaser_concrete_model, scene_data.chaser_jt_limits);

  std::shared_ptr<jt_space_type> normal_jt_space = make_manip_jt_space<Order>(
      chaser_concrete_model, scene_data.chaser_jt_limits);

  (*workspace) << scene_data.chaser_target_proxy;
  for (const auto& chaser_env_proxie : scene_data.chaser_env_proxies) {
    (*workspace) << chaser_env_proxie;
  }

  // Create the start and goal points:
  rl_point_type start_point;
  rl_point_type goal_point;
  point_type start_inter;
  point_type goal_inter;

  start_inter = normal_jt_space->origin();
  get<0>(start_inter) = jt_start;
  start_point = joint_limits_mapping<double>(scene_data.chaser_jt_limits)
                    .map_to_space(start_inter, *normal_jt_space, *jt_space);

  goal_inter = normal_jt_space->origin();
  get<0>(goal_inter) = jt_desired;
  goal_point = joint_limits_mapping<double>(scene_data.chaser_jt_limits)
                   .map_to_space(goal_inter, *normal_jt_space, *jt_space);

  execute_p2p_planner(workspace, plan_options, workspace_dims, engine,
                      start_point, goal_point);

  // Restore model's state:
  chaser_concrete_model->setJointPositions(jt_start);
  chaser_concrete_model->doDirectMotion();
}

template <typename PlanEngine>
void CRS_launch_static_planner(
    const ReaK::kte::chaser_target_data& scene_data,
    const ReaK::pp::planning_option_collection& plan_options,
    const ReaK::pp::planning_space_options& space_def, PlanEngine& engine,
    const ReaK::vect_n<double>& jt_start,
    const ReaK::vect_n<double>& jt_desired) {

  if ((space_def.get_space_order() == 0) && (space_def.get_interp_id() == 0)) {
    CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics,
                               ReaK::pp::linear_interpolation_tag, 0>(
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
      if ((space_def.get_space_order() == 1) &&
          (space_def.get_interp_id() == 1)) {
    CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics,
                               ReaK::pp::cubic_hermite_interpolation_tag, 1>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
  } else
#if 0
    if((space_def.get_space_order() == 2) && (space_def.get_interp_id() == 1)) {
      CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics, ReaK::pp::cubic_hermite_interpolation_tag, 2>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
    } else
#endif
      if ((space_def.get_space_order() == 2) &&
          (space_def.get_interp_id() == 2)) {
    CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics,
                               ReaK::pp::quintic_hermite_interpolation_tag, 2>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
  } else if ((space_def.get_space_order() == 1) &&
             (space_def.get_interp_id() == 3)) {
    CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics,
                               ReaK::pp::svp_Ndof_interpolation_tag, 1>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
  } else
#if 0
    if((space_def.get_space_order() == 2) && (space_def.get_interp_id() == 3)) {
      CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics, ReaK::pp::svp_Ndof_interpolation_tag, 2>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
    } else
#endif
      if ((space_def.get_space_order() == 2) &&
          (space_def.get_interp_id() == 4)) {
    CRS_execute_static_planner<ReaK::kte::manip_P3R3R_kinematics,
                               ReaK::pp::sap_Ndof_interpolation_tag, 2>(
        scene_data, plan_options, space_def, engine, jt_start, jt_desired);
  }
}

template <typename ManipMdlType, typename InterpTag, int Order,
          typename PlanEngine>
void CRS_execute_dynamic_planner(
    const ReaK::kte::chaser_target_data& scene_data,
    const ReaK::pp::planning_option_collection& plan_options,
    const ReaK::pp::planning_space_options& space_def, PlanEngine& engine,
    const ReaK::vect_n<double>& jt_start,
    const std::shared_ptr<sat_trajectory_type>& target_state_traj) {
  using namespace ReaK;
  using namespace pp;

  std::shared_ptr<ManipMdlType> chaser_concrete_model =
      rtti::rk_dynamic_ptr_cast<ManipMdlType>(scene_data.chaser_kin_model);

  if (!chaser_concrete_model) {
    std::cout << "Could not cast the chaser model to the given type!"
              << std::endl;
    if (scene_data.chaser_kin_model) {
      std::cout << "Because the chaser model is a null pointer!" << std::endl;
    }
    if (dynamic_cast<ManipMdlType*>(scene_data.chaser_kin_model.get())) {
      std::cout << "But the C++RTTI could perform the dynamic-cast correctly!"
                << std::endl;
    }
    return;
  }

  vect_n<double> target_jt_original =
      scene_data.target_kin_model->getJointPositions();

  using dynamic_workspace_type =
      typename manip_dynamic_workspace<ManipMdlType, Order>::rl_workspace_type;
  using rl_jt_space_type =
      typename manip_pp_traits<ManipMdlType, Order>::rl_jt_space_type;
  using jt_space_type =
      typename manip_pp_traits<ManipMdlType, Order>::jt_space_type;
  using rlIK_map_type = typename manip_IK_map<ManipMdlType, Order>::rl_map_type;

  using rl_point_type = typename topology_traits<rl_jt_space_type>::point_type;
  using point_type = typename topology_traits<jt_space_type>::point_type;

  using dynamic_super_space_type = typename subspace_traits<
      dynamic_workspace_type>::super_space_type;  // SuperSpaceType

  std::size_t workspace_dims =
      (Order + 1) * manip_pp_traits<ManipMdlType, Order>::degrees_of_freedom;

  std::shared_ptr<frame_3D<double>> EE_frame =
      chaser_concrete_model->getDependentFrame3D(0)->mFrame;

  // Create the workspace for the chaser manipulator:
  std::shared_ptr<dynamic_workspace_type> workspace =
      make_manip_dynamic_workspace<Order>(
          InterpTag(), chaser_concrete_model, scene_data.chaser_jt_limits,
          space_def.min_travel,
          target_state_traj->get_end_time() -
              target_state_traj->get_start_time(),
          plan_options.start_delay + target_state_traj->get_start_time());

  std::shared_ptr<rl_jt_space_type> jt_space = make_manip_rl_jt_space<Order>(
      chaser_concrete_model, scene_data.chaser_jt_limits);

  std::shared_ptr<jt_space_type> normal_jt_space = make_manip_jt_space<Order>(
      chaser_concrete_model, scene_data.chaser_jt_limits);

  // Add the proximity-query models for the target and the environment to the workspace's collision-checking:
  (*workspace) << scene_data.chaser_target_proxy;

  for (const auto& chaser_env_proxie : scene_data.chaser_env_proxies) {
    (*workspace) << chaser_env_proxie;
  }

  // Add the proximity-query model updater for the target's trajectory:
  using target_proxy_applicator_type =
      proxy_traj_applicator<sat_trajectory_type>;
  using target_jt_space_type =
      typename target_proxy_applicator_type::joint_space_type;

  std::shared_ptr<target_proxy_applicator_type> target_proxy_applicator(
      new target_proxy_applicator_type(
          make_any_model_applicator<target_jt_space_type>(
              manip_direct_kin_map(scene_data.target_kin_model)),
          target_state_traj));

  workspace->add_proxy_model_updater(target_proxy_applicator);

  // Create the starting configuration as a joint-space configuration.
  using rl_temporal_point_type =
      typename topology_traits<dynamic_workspace_type>::point_type;

  point_type start_inter = normal_jt_space->origin();
  get<0>(start_inter) = jt_start;
  rl_point_type start_point =
      joint_limits_mapping<double>(scene_data.chaser_jt_limits)
          .map_to_space(start_inter, *normal_jt_space, *jt_space);
  rl_temporal_point_type temporal_start_point(
      plan_options.start_delay + target_state_traj->get_start_time(),
      start_point);

  // Create the interception query object:

  using dkik_map_type = manip_dk_ik_map<manip_direct_kin_map, rlIK_map_type>;
  using temporal_dkik_map_type = temporal_topo_map<dkik_map_type>;
  using mapped_traj_type =
      transformed_trajectory<dynamic_super_space_type, sat_trajectory_type,
                             temporal_dkik_map_type>;

  std::shared_ptr<mapped_traj_type> mapped_goal_traj(new mapped_traj_type(
      std::shared_ptr<dynamic_super_space_type>(&(workspace->get_super_space()),
                                                null_deleter()),
      target_state_traj,
      temporal_dkik_map_type(dkik_map_type(
          manip_direct_kin_map(scene_data.target_kin_model),
          rlIK_map_type(
              chaser_concrete_model,
              joint_limits_mapping<double>(scene_data.chaser_jt_limits))))));

  execute_intercept_planner(workspace, plan_options, workspace_dims + 1,
                            target_state_traj->get_end_time(),
                            space_def.min_travel, engine, temporal_start_point,
                            mapped_goal_traj);

  // Difference with the GUI's dynexec "manual" code (now done within the 'execute_intercept_planner' function and the
  // engine:
  //   plan_stopper = std::bind(&planner_base< dynamic_workspace_type >::stop, workspace_planner);
  //   pp_query.reset_solution_records();
  //   workspace_planner->solve_planning_query(pp_query);
  //   plan_stopper = std::function<void()>(); // clear the stopper function-pointer.

#if 0
  std::shared_ptr< seq_trajectory_base< dynamic_super_space_type > > bestsol_rltraj;
  if(pp_query.solutions.size())
    bestsol_rltraj = pp_query.solutions.begin()->second;
  std::cout << "The shortest distance is: " << pp_query.get_best_solution_distance() << std::endl;
#endif

#if 0
  sol_trace.reset();
  if(bestsol_rltraj) {
    
    typedef hyperbox_topology< vect<double,7> > manip_cspace_type;
    typedef temporal_space<manip_cspace_type, time_poisson_topology, time_distance_only> temporal_space_type;
    typedef discrete_point_trajectory< temporal_space_type > trajectory_type;
    typedef trajectory_wrapper<trajectory_type> wrapped_traj_type;
    
    std::shared_ptr< wrapped_traj_type > tmp_sol_traj(new wrapped_traj_type("chaser_sol_trajectory", trajectory_type()));
    
    trajectory_type& sol_traj = tmp_sol_traj->get_underlying_trajectory();
    
    typedef typename seq_trajectory_base< dynamic_super_space_type >::point_time_iterator PtIter;
    typedef typename spatial_trajectory_traits<trajectory_type>::point_type TCSpacePointType;
    
    for(PtIter it = bestsol_rltraj->begin_time_travel(); it != bestsol_rltraj->end_time_travel(); it += 0.001) {
      rl_temporal_point_type cur_pt = *it;
      sol_traj.push_back( TCSpacePointType(cur_pt.time, get<0>(joint_limits_mapping<double>(scene_data.chaser_jt_limits).map_to_space(cur_pt.pt, *jt_space, *normal_jt_space))) );
    };
    sol_trace = tmp_sol_traj;
  };
#endif

  chaser_concrete_model->setJointPositions(jt_start);
  chaser_concrete_model->doDirectMotion();

  scene_data.target_kin_model->setJointPositions(target_jt_original);
  scene_data.target_kin_model->doDirectMotion();
}

template <typename PlanEngine>
void CRS_launch_dynamic_planner(
    const ReaK::kte::chaser_target_data& scene_data,
    const ReaK::pp::planning_option_collection& plan_options,
    const ReaK::pp::planning_space_options& space_def, PlanEngine& engine,
    const ReaK::vect_n<double>& jt_start,
    const std::shared_ptr<sat_trajectory_type>& sat_trajectory) {
  using namespace ReaK;
  using namespace pp;

  if ((space_def.get_space_order() == 0) && (space_def.get_interp_id() == 0)) {
    CRS_execute_dynamic_planner<kte::manip_P3R3R_kinematics,
                                linear_interpolation_tag, 0>(
        scene_data, plan_options, space_def, engine, jt_start, sat_trajectory);
  } else
#if 0
  if((space_def.get_space_order() == 1) && (space_def.get_interp_id() == 0)) {
    CRS_execute_dynamic_planner<kte::manip_P3R3R_kinematics, linear_interpolation_tag, 1>(
      scene_data, plan_options, space_def, engine, jt_start, sat_trajectory);
  } else 
  if((space_def.get_space_order() == 2) && (space_def.get_interp_id() == 0)) {
    CRS_execute_dynamic_planner<kte::manip_P3R3R_kinematics, linear_interpolation_tag, 2>(
      scene_data, plan_options, space_def, engine, jt_start, sat_trajectory);
  } else
#endif
      if ((space_def.get_space_order() == 1) &&
          (space_def.get_interp_id() == 1)) {
    CRS_execute_dynamic_planner<kte::manip_P3R3R_kinematics,
                                cubic_hermite_interpolation_tag, 1>(
        scene_data, plan_options, space_def, engine, jt_start, sat_trajectory);
  } else
#if 0
  if((space_def.get_space_order() == 2) && (space_def.get_interp_id() == 1)) {
    CRS_execute_dynamic_planner<kte::manip_P3R3R_kinematics, cubic_hermite_interpolation_tag, 2>(
      scene_data, plan_options, space_def, engine, jt_start, sat_trajectory);
  } else
#endif
      if ((space_def.get_space_order() == 2) &&
          (space_def.get_interp_id() == 2)) {
    CRS_execute_dynamic_planner<kte::manip_P3R3R_kinematics,
                                quintic_hermite_interpolation_tag, 2>(
        scene_data, plan_options, space_def, engine, jt_start, sat_trajectory);
  } else if ((space_def.get_space_order() == 1) &&
             (space_def.get_interp_id() == 3)) {
    CRS_execute_dynamic_planner<kte::manip_P3R3R_kinematics,
                                svp_Ndof_interpolation_tag, 1>(
        scene_data, plan_options, space_def, engine, jt_start, sat_trajectory);
  } else
#if 0
  if((space_def.get_space_order() == 2) && (space_def.get_interp_id() == 3)) {
    CRS_execute_dynamic_planner<kte::manip_P3R3R_kinematics, svp_Ndof_interpolation_tag, 2>(
      scene_data, plan_options, space_def, engine, jt_start, sat_trajectory);
  } else
#endif
      if ((space_def.get_space_order() == 2) &&
          (space_def.get_interp_id() == 4)) {
    CRS_execute_dynamic_planner<kte::manip_P3R3R_kinematics,
                                sap_Ndof_interpolation_tag, 2>(
        scene_data, plan_options, space_def, engine, jt_start, sat_trajectory);
  }
}

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

  absl::ParseCommandLine(argc, argv);

  if (static_cast<int>(absl::GetFlag(FLAGS_monte_carlo)) +
          static_cast<int>(absl::GetFlag(FLAGS_single_run)) +
          static_cast<int>(!absl::GetFlag(FLAGS_generate_all_files).empty()) +
          static_cast<int>(
              !absl::GetFlag(FLAGS_generate_planner_options).empty()) +
          static_cast<int>(
              !absl::GetFlag(FLAGS_generate_chaser_target_env).empty()) +
          static_cast<int>(
              !absl::GetFlag(FLAGS_generate_space_definition).empty()) +
          static_cast<int>(
              !absl::GetFlag(FLAGS_generate_start_config).empty()) +
          static_cast<int>(!absl::GetFlag(FLAGS_generate_target_pose).empty()) <
      1) {
    std::cout << "Error: There was no action specified! This program is "
                 "designed to perform Monte-Carlo runs, single "
                 "runs (with output), or generate the configuration files to "
                 "construct scenarios. You must specify at "
                 "least one of these actions to be performed!"
              << std::endl;
    return 1;
  }

  std::string output_path_name = absl::GetFlag(FLAGS_output_path);
  while (output_path_name.back() == '/') {
    output_path_name.erase(output_path_name.length() - 1, 1);
  }

  fs::create_directory(output_path_name.c_str());

  std::string result_file_prefix = "CRS_static_scene";
  if (!absl::GetFlag(FLAGS_result_file_prefix).empty()) {
    result_file_prefix = absl::GetFlag(FLAGS_result_file_prefix);
  }

  planning_option_collection plan_options = get_planning_option_from_flags();

  std::string knn_method_str = plan_options.get_knn_method_str();
  std::string mg_storage_str = plan_options.get_mg_storage_str();
  std::string planner_qualifier_str = plan_options.get_planner_qualifier_str();
  std::string planner_name_str = plan_options.get_planning_algo_str() + "_" +
                                 planner_qualifier_str + "_" + mg_storage_str +
                                 "_" + knn_method_str;

  planning_space_options space_def = get_planning_space_options_from_flags();

  chaser_target_data scene_data = get_chaser_target_data_from_flags();

  vect_n<double> jt_start(7, 0.0);
  if (scene_data.chaser_kin_model) {
    jt_start = scene_data.chaser_kin_model->getJointPositions();
    if (!absl::GetFlag(FLAGS_start_configuration).empty()) {
      try {
        vect_n<double> jt_start_tmp = jt_start;
        (*serialization::open_iarchive(
            absl::GetFlag(FLAGS_start_configuration))) >>
            jt_start_tmp;
        jt_start = jt_start_tmp;
      } catch (std::exception& e) {
        RK_UNUSED(e);
        std::cerr << "Error: Could not load the start-configuration file!"
                  << std::endl;
      }
    }
  }

  frame_3D<double> target_frame;
  if (scene_data.target_kin_model) {
    target_frame = scene_data.target_kin_model->getFrame3D(0)->getGlobalFrame();
    if (!absl::GetFlag(FLAGS_target_pose).empty()) {
      try {
        frame_3D<double> target_frame_tmp = target_frame;
        (*serialization::open_iarchive(absl::GetFlag(FLAGS_target_pose))) >>
            target_frame_tmp;
        target_frame = target_frame_tmp;
        *(scene_data.target_kin_model->getFrame3D(0)) = target_frame;
        scene_data.target_kin_model->doDirectMotion();
      } catch (std::exception& e) {
        RK_UNUSED(e);
        std::cerr << "Error: Could not load the target-pose file!" << std::endl;
      }
    }
  }

  vect_n<double> jt_desired(7, 0.0);
  if (absl::GetFlag(FLAGS_target_trajectory).empty() &&
      scene_data.chaser_kin_model) {
    std::shared_ptr<frame_3D<double>> dep_EE_frame =
        scene_data.chaser_kin_model->getDependentFrame3D(0)->mFrame;
    if (absl::GetFlag(FLAGS_monte_carlo) ||
        static_cast<int>(absl::GetFlag(FLAGS_single_run)) > 0) {
      try {
        *dep_EE_frame = scene_data.target_frame->getGlobalFrame();
        scene_data.chaser_kin_model->doInverseMotion();
        jt_desired = scene_data.chaser_kin_model->getJointPositions();
      } catch (optim::infeasible_problem& e) {
        RK_UNUSED(e);
        std::cerr << "Error: The target frame cannot be reached! No inverse "
                     "kinematics solution possible!"
                  << std::endl;
        return 10;
      }
      scene_data.chaser_kin_model->setJointPositions(jt_start);
      scene_data.chaser_kin_model->doDirectMotion();
    }
  }

  std::shared_ptr<sat_trajectory_type> target_state_traj;
  if (scene_data.target_kin_model) {
    if (!absl::GetFlag(FLAGS_target_trajectory).empty()) {
      try {
        using wrapped_traj_type = trajectory_wrapper<
            discrete_point_trajectory<sat_temporal_space_type>>;
        std::shared_ptr<wrapped_traj_type> tmp_traj(new wrapped_traj_type());

        (*serialization::open_iarchive(
            absl::GetFlag(FLAGS_target_trajectory))) &
            RK_SERIAL_LOAD_WITH_ALIAS("se3_trajectory",
                                      tmp_traj->get_underlying_trajectory());

        target_state_traj = tmp_traj;
      } catch (std::exception& e) {
        RK_UNUSED(e);
        std::cerr << "Error: Could not load the target-trajectory file!"
                  << std::endl;
        return 11;
      }
    }
  }

  // Do the generations if required:

  if (!absl::GetFlag(FLAGS_generate_all_files).empty() ||
      !absl::GetFlag(FLAGS_generate_planner_options).empty()) {
    std::string file_name;
    if (absl::GetFlag(FLAGS_generate_planner_options).empty()) {
      file_name = absl::GetFlag(FLAGS_generate_all_files) + "_planner";
    } else {
      file_name = absl::GetFlag(FLAGS_generate_planner_options);
    }
    if (absl::GetFlag(FLAGS_generate_protobuf)) {
      file_name += ".pbuf";
    } else if (absl::GetFlag(FLAGS_generate_binary)) {
      file_name += ".rkb";
    } else {
      file_name += ".rkx";
    }

    try {
      (*serialization::open_oarchive(file_name)) << plan_options;
    } catch (std::exception& e) {
      RK_UNUSED(e);
      std::cerr << "Error: Could not generate the planner options file!"
                << std::endl;
    }
  }

  if (!absl::GetFlag(FLAGS_generate_all_files).empty() ||
      !absl::GetFlag(FLAGS_generate_chaser_target_env).empty()) {
    std::string file_name;
    if (absl::GetFlag(FLAGS_generate_chaser_target_env).empty()) {
      file_name = absl::GetFlag(FLAGS_generate_all_files) + "_models";
    } else {
      file_name = absl::GetFlag(FLAGS_generate_chaser_target_env);
    }
    if (absl::GetFlag(FLAGS_generate_protobuf)) {
      file_name += ".pbuf";
    } else if (absl::GetFlag(FLAGS_generate_binary)) {
      file_name += ".rkb";
    } else {
      file_name += ".rkx";
    }

    try {
      (*serialization::open_oarchive(file_name)) << scene_data;
    } catch (std::exception& e) {
      RK_UNUSED(e);
      std::cerr << "Error: Could not generate the chaser-target-env model file!"
                << std::endl;
    }
  }

  if (!absl::GetFlag(FLAGS_generate_all_files).empty() ||
      !absl::GetFlag(FLAGS_generate_space_definition).empty()) {
    std::string file_name;
    if (absl::GetFlag(FLAGS_generate_space_definition).empty()) {
      file_name = absl::GetFlag(FLAGS_generate_all_files) + "_space";
    } else {
      file_name = absl::GetFlag(FLAGS_generate_space_definition);
    }
    if (absl::GetFlag(FLAGS_generate_protobuf)) {
      file_name += ".pbuf";
    } else if (absl::GetFlag(FLAGS_generate_binary)) {
      file_name += ".rkb";
    } else {
      file_name += ".rkx";
    }

    try {
      (*serialization::open_oarchive(file_name)) << space_def;
    } catch (std::exception& e) {
      RK_UNUSED(e);
      std::cerr << "Error: Could not generate the space-definition file!"
                << std::endl;
    }
  }

  if (!absl::GetFlag(FLAGS_generate_all_files).empty() ||
      !absl::GetFlag(FLAGS_generate_start_config).empty()) {
    std::string file_name;
    if (absl::GetFlag(FLAGS_generate_start_config).empty()) {
      file_name = absl::GetFlag(FLAGS_generate_all_files) + "_start_config";
    } else {
      file_name = absl::GetFlag(FLAGS_generate_start_config);
    }
    if (absl::GetFlag(FLAGS_generate_protobuf)) {
      file_name += ".pbuf";
    } else if (absl::GetFlag(FLAGS_generate_binary)) {
      file_name += ".rkb";
    } else {
      file_name += ".rkx";
    }

    try {
      (*serialization::open_oarchive(file_name)) << jt_start;
    } catch (std::exception& e) {
      RK_UNUSED(e);
      std::cerr << "Error: Could not generate the start-configuration file!"
                << std::endl;
    }
  }

  if (!absl::GetFlag(FLAGS_generate_all_files).empty() ||
      !absl::GetFlag(FLAGS_generate_target_pose).empty()) {
    std::string file_name;
    if (absl::GetFlag(FLAGS_generate_target_pose).empty()) {
      file_name = absl::GetFlag(FLAGS_generate_all_files) + "_target_pose";
    } else {
      file_name = absl::GetFlag(FLAGS_generate_target_pose);
    }
    if (absl::GetFlag(FLAGS_generate_protobuf)) {
      file_name += ".pbuf";
    } else if (absl::GetFlag(FLAGS_generate_binary)) {
      file_name += ".rkb";
    } else {
      file_name += ".rkx";
    }

    try {
      (*serialization::open_oarchive(file_name)) << target_frame;
    } catch (std::exception& e) {
      RK_UNUSED(e);
      std::cerr << "Error: Could not generate the target-pose file!"
                << std::endl;
    }
  }

  if (space_def.is_temporal_space()) {

    // Do the Monte-Carlo runs if required:
    if (absl::GetFlag(FLAGS_monte_carlo)) {
      monte_carlo_mp_engine mc_eng(absl::GetFlag(FLAGS_mc_runs),
                                   planner_name_str,
                                   output_path_name + "/" + result_file_prefix);
      try {
        CRS_launch_dynamic_planner(scene_data, plan_options, space_def, mc_eng,
                                   jt_start, target_state_traj);
      } catch (std::exception& e) {
        std::cerr
            << "Error: An exception was raised during the planning:\nwhat(): "
            << e.what() << std::endl;
        return 2;
      }
    }

    // Do a single run if required:
    if (absl::GetFlag(FLAGS_single_run)) {
      vlist_print_mp_engine sr_eng(planner_name_str,
                                   output_path_name + "/" + result_file_prefix);
      try {
        CRS_launch_dynamic_planner(scene_data, plan_options, space_def, sr_eng,
                                   jt_start, target_state_traj);
      } catch (std::exception& e) {
        std::cerr
            << "Error: An exception was raised during the planning:\nwhat(): "
            << e.what() << std::endl;
        return 3;
      }
    }

  } else {

    // Do the Monte-Carlo runs if required:
    if (absl::GetFlag(FLAGS_monte_carlo)) {
      monte_carlo_mp_engine mc_eng(absl::GetFlag(FLAGS_mc_runs),
                                   planner_name_str,
                                   output_path_name + "/" + result_file_prefix);
      try {
        CRS_launch_static_planner(scene_data, plan_options, space_def, mc_eng,
                                  jt_start, jt_desired);
      } catch (std::exception& e) {
        std::cerr
            << "Error: An exception was raised during the planning:\nwhat(): "
            << e.what() << std::endl;
        return 2;
      }
    }

    // Do a single run if required:
    if (absl::GetFlag(FLAGS_single_run)) {
      vlist_print_mp_engine sr_eng(planner_name_str,
                                   output_path_name + "/" + result_file_prefix);
      try {
        CRS_launch_static_planner(scene_data, plan_options, space_def, sr_eng,
                                  jt_start, jt_desired);
      } catch (std::exception& e) {
        std::cerr
            << "Error: An exception was raised during the planning:\nwhat(): "
            << e.what() << std::endl;
        return 3;
      }
    }
  }

  return 0;
}

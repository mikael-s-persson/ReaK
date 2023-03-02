
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

#include "CRS_planner_impl.hpp"

#include <QMessageBox>

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>

#include <ReaK/geometry/proximity/proxy_query_model.hpp>
#include <ReaK/mbd/coin3D/oi_scene_graph.hpp>

#include <ReaK/mbd/kte/kte_map_chain.hpp>
#include <ReaK/mbd/models/manip_dynamics_model.hpp>

#include <ReaK/topologies/spaces/manip_P3R3R_workspaces.hpp>

#include <ReaK/planning/path_planning/p2p_planning_query.hpp>

#include <ReaK/planning/path_planning/rrtstar_manip_planners.hpp>
#include <ReaK/planning/path_planning/sbastar_manip_planners.hpp>

#if 0
#include <ReaK/planning/path_planning/fadprm_path_planner.hpp>
#include <ReaK/planning/path_planning/prm_path_planner.hpp>
#include <ReaK/planning/path_planning/rrt_path_planner.hpp>
#endif

#include "CRS_planner_data.hpp"

#include <ReaK/planning/path_planning/frame_tracer_coin3d.hpp>

#include <ReaK/math/optimization/optim_exceptions.hpp>

#include <ReaK/topologies/interpolation/discrete_point_trajectory.hpp>
#include <ReaK/topologies/interpolation/trajectory_base.hpp>
#include <ReaK/topologies/spaces/manip_planning_traits.hpp>

#include <ReaK/topologies/spaces/Ndof_cubic_spaces.hpp>
#include <ReaK/topologies/spaces/Ndof_linear_spaces.hpp>
#include <ReaK/topologies/spaces/Ndof_quintic_spaces.hpp>
#include <ReaK/topologies/spaces/Ndof_sap_spaces.hpp>
#include <ReaK/topologies/spaces/Ndof_svp_spaces.hpp>

template <typename ManipMdlType, typename InterpTag, int Order,
          typename ManipCSpaceTrajectory>
void CRS_execute_static_planner_impl(
    const ReaK::kte::chaser_target_data& scene_data,
    const ReaK::pp::planning_option_collection& plan_options,
    SoSwitch* sw_motion_graph, SoSwitch* sw_solutions, double min_travel,
    bool print_timing, bool print_counter, const ReaK::vect_n<double>& jt_start,
    const ReaK::vect_n<double>& jt_desired,
    std::shared_ptr<ManipCSpaceTrajectory>& sol_trace,
    std::function<void()>& plan_stopper) {
  using namespace ReaK;
  using namespace pp;

  std::shared_ptr<ManipMdlType> chaser_concrete_model =
      rtti::rk_dynamic_ptr_cast<ManipMdlType>(scene_data.chaser_kin_model);
  if (!chaser_concrete_model)
    return;

  typedef
      typename manip_static_workspace<ManipMdlType, Order>::rl_workspace_type
          static_workspace_type;
  typedef typename manip_pp_traits<ManipMdlType, Order>::rl_jt_space_type
      rl_jt_space_type;
  typedef typename manip_pp_traits<ManipMdlType, Order>::jt_space_type
      jt_space_type;
  typedef typename manip_DK_map<ManipMdlType, Order>::rl_map_type rlDK_map_type;

  typedef typename topology_traits<rl_jt_space_type>::point_type rl_point_type;
  typedef typename topology_traits<jt_space_type>::point_type point_type;

  typedef typename subspace_traits<static_workspace_type>::super_space_type
      static_super_space_type;  // SuperSpaceType

  std::size_t workspace_dims =
      (Order + 1) * manip_pp_traits<ManipMdlType, Order>::degrees_of_freedom;

  std::shared_ptr<frame_3D<double>> EE_frame =
      chaser_concrete_model->getDependentFrame3D(0)->mFrame;

  std::shared_ptr<static_workspace_type> workspace =
      make_manip_static_workspace<Order>(InterpTag(), chaser_concrete_model,
                                         scene_data.chaser_jt_limits,
                                         min_travel);

  std::shared_ptr<rl_jt_space_type> jt_space = make_manip_rl_jt_space<Order>(
      chaser_concrete_model, scene_data.chaser_jt_limits);

  std::shared_ptr<jt_space_type> normal_jt_space = make_manip_jt_space<Order>(
      chaser_concrete_model, scene_data.chaser_jt_limits);

  (*workspace) << scene_data.chaser_target_proxy;
  for (std::size_t i = 0; i < scene_data.chaser_env_proxies.size(); ++i)
    (*workspace) << scene_data.chaser_env_proxies[i];

  rl_point_type start_point, goal_point;
  point_type start_inter, goal_inter;

  start_inter = normal_jt_space->origin();
  get<0>(start_inter) = jt_start;
  start_point = joint_limits_mapping<double>(scene_data.chaser_jt_limits)
                    .map_to_space(start_inter, *normal_jt_space, *jt_space);

  goal_inter = normal_jt_space->origin();
  get<0>(goal_inter) = jt_desired;
  goal_point = joint_limits_mapping<double>(scene_data.chaser_jt_limits)
                   .map_to_space(goal_inter, *normal_jt_space, *jt_space);

  // Create the reporter chain.
  any_sbmp_reporter_chain<static_workspace_type> report_chain;

  // Create the frame tracing reporter.
  typedef frame_tracer_3D<rl_jt_space_type> frame_reporter_type;

  frame_reporter_type temp_reporter(
      make_any_model_applicator<rl_jt_space_type>(rlDK_map_type(
          chaser_concrete_model, scene_data.chaser_jt_limits, normal_jt_space)),
      0.5 * min_travel, (sw_motion_graph != nullptr));

  if ((sw_motion_graph) || (sw_solutions)) {
    temp_reporter.add_traced_frame(EE_frame);
    report_chain.add_reporter(std::ref(temp_reporter));
  };

  if (print_counter)
    report_chain.add_reporter(print_sbmp_progress<>());

  if (print_timing)
    report_chain.add_reporter(timing_sbmp_report<>());

  path_planning_p2p_query<static_workspace_type> pp_query(
      "pp_query", workspace, start_point, goal_point, plan_options.max_results);

  std::shared_ptr<sample_based_planner<static_workspace_type>>
      workspace_planner;

#if 0
  if( plan_options.planning_algo == 0 ) { // RRT
    
    workspace_planner = std::shared_ptr< sample_based_planner< static_workspace_type > >(
      new rrt_planner< static_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, report_chain));
    
  } else
#endif
  if (plan_options.planning_algo == 1) {  // RRT*

    workspace_planner =
        std::shared_ptr<sample_based_planner<static_workspace_type>>(
            new rrtstar_planner<static_workspace_type>(
                workspace, plan_options.max_vertices,
                plan_options.prog_interval,
                plan_options.store_policy | plan_options.knn_method,
                plan_options.planning_options, 0.1, 0.05, workspace_dims,
                report_chain));

  } else
#if 0
  if( plan_options.planning_algo == 2 ) { // PRM
    
    workspace_planner = std::shared_ptr< sample_based_planner< static_workspace_type > >(
      new prm_planner< static_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, plan_options.max_random_walk, workspace_dims, report_chain));
    
  } else
#endif
      if (plan_options.planning_algo == 3) {  // SBA*

    std::shared_ptr<sbastar_planner<static_workspace_type>> tmp(
        new sbastar_planner<static_workspace_type>(
            workspace, plan_options.max_vertices, plan_options.prog_interval,
            plan_options.store_policy | plan_options.knn_method,
            plan_options.planning_options, 0.1, 0.05,
            plan_options.max_random_walk, workspace_dims, report_chain));

    tmp->set_initial_density_threshold(0.0);
    tmp->set_initial_relaxation(plan_options.init_relax);
    tmp->set_initial_SA_temperature(plan_options.init_SA_temp);

    workspace_planner = tmp;

  }
#if 0
  else if( plan_options.planning_algo == 4 ) { // FADPRM
    
    std::shared_ptr< fadprm_planner< static_workspace_type > > tmp(
      new fadprm_planner< static_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        0.1, 0.05, plan_options.max_random_walk, workspace_dims, report_chain));
    
    tmp->set_initial_relaxation(plan_options.init_relax);
    
    workspace_planner = tmp;
    
  }
#endif
  else {
    std::cout << "The desired planner algorithm is not supported!" << std::endl;
  };

  if (!workspace_planner) {
    std::cout << "The workspace planner was not created successfully!"
              << std::endl;
    return;
  };

  plan_stopper =
      std::bind(&planner_base<static_workspace_type>::stop, workspace_planner);

  pp_query.reset_solution_records();
  workspace_planner->solve_planning_query(pp_query);

  plan_stopper = std::function<void()>();

  std::shared_ptr<seq_path_base<static_super_space_type>> bestsol_rlpath;
  if (pp_query.solutions.size())
    bestsol_rlpath = pp_query.solutions.begin()->second;
  std::cout << "The shortest distance is: "
            << pp_query.get_best_solution_distance() << std::endl;

  sol_trace.reset();
  if (bestsol_rlpath) {
    typedef hyperbox_topology<vect<double, 7>> manip_cspace_type;
    typedef temporal_space<manip_cspace_type, time_poisson_topology,
                           time_distance_only>
        temporal_space_type;
    typedef discrete_point_trajectory<temporal_space_type> trajectory_type;
    typedef trajectory_wrapper<trajectory_type> wrapped_traj_type;

    std::shared_ptr<wrapped_traj_type> tmp_sol_traj(
        new wrapped_traj_type("chaser_sol_trajectory", trajectory_type()));
    trajectory_type& sol_traj = tmp_sol_traj->get_underlying_trajectory();

    typedef
        typename seq_path_base<static_super_space_type>::point_fraction_iterator
            PtIter;
    typedef typename spatial_trajectory_traits<trajectory_type>::point_type
        TCSpacePointType;
    double t = 0.0;
    for (PtIter it = bestsol_rlpath->begin_fraction_travel();
         it != bestsol_rlpath->end_fraction_travel(); it += 0.1, t += 0.1) {
      sol_traj.push_back(TCSpacePointType(
          t, get<0>(joint_limits_mapping<double>(scene_data.chaser_jt_limits)
                        .map_to_space(*it, *jt_space, *normal_jt_space))));
    };
    sol_trace = tmp_sol_traj;
  };

  // Check the motion-graph separator and solution separators
  //  add them to the switches.

  if (sw_motion_graph) {
    SoSeparator* mg_sep =
        temp_reporter.get_motion_graph_tracer(EE_frame).get_separator();
    if (mg_sep)
      mg_sep->ref();

    SoDB::writelock();
    sw_motion_graph->removeAllChildren();
    if (mg_sep) {
      sw_motion_graph->addChild(mg_sep);
      mg_sep->unref();
    };
    SoDB::writeunlock();
  };

  if (sw_solutions) {
    SoSeparator* sol_sep = nullptr;
    if (temp_reporter.get_solution_count()) {
      sol_sep = temp_reporter.get_solution_tracer(EE_frame, 0).get_separator();
      if (sol_sep)
        sol_sep->ref();
    };

    SoDB::writelock();
    sw_solutions->removeAllChildren();
    if (sol_sep) {
      sw_solutions->addChild(sol_sep);
      sol_sep->unref();
    };
    SoDB::writeunlock();
  };

  chaser_concrete_model->setJointPositions(jt_start);
  chaser_concrete_model->doDirectMotion();
};

void CRSPlannerGUI::executePlanner() {
  using namespace ReaK;
  using namespace pp;

  std::shared_ptr<frame_3D<double>> EE_frame =
      ct_config.sceneData.chaser_kin_model->getDependentFrame3D(0)->mFrame;

  vect_n<double> jt_desired(7, 0.0);

  vect_n<double> jt_previous =
      ct_config.sceneData.chaser_kin_model->getJointPositions();

  *EE_frame = *(ct_config.sceneData.target_frame);
  ct_config.sceneData.chaser_kin_model->doInverseMotion();
  jt_desired = ct_config.sceneData.chaser_kin_model->getJointPositions();

  ct_config.sceneData.chaser_kin_model->setJointPositions(jt_previous);
  ct_config.sceneData.chaser_kin_model->doDirectMotion();

  vect_n<double> jt_start =
      ct_config.sceneData.chaser_kin_model->getJointPositions();

  SoSwitch* sw_motion_graph = nullptr;
  if (plan_alg_config.outputMotionGraph())
    sw_motion_graph = view3d_menu.getDisplayGroup("Motion-Graph", true);

  SoSwitch* sw_solutions = nullptr;
  if (plan_alg_config.outputSolution())
    sw_solutions = view3d_menu.getDisplayGroup("Solution(s)", true);

  bool print_timing = plan_alg_config.outputTiming();
  bool print_counter = plan_alg_config.outputNodeCounter();

  if ((space_config.space_order == 0) && (space_config.interp_id == 0)) {
    CRS_execute_static_planner_impl<kte::manip_P3R3R_kinematics,
                                    linear_interpolation_tag, 0>(
        ct_config.sceneData, plan_alg_config.planOptions, sw_motion_graph,
        sw_solutions, space_config.min_travel, print_timing, print_counter,
        jt_start, jt_desired, sol_anim.trajectory, stop_planner);
  } else
#if 0
  if((space_config.space_order == 1) && (space_config.interp_id == 0)) {
    CRS_execute_static_planner_impl<kte::manip_P3R3R_kinematics, linear_interpolation_tag, 1>(ct_config.sceneData, plan_alg_config.planOptions,
      sw_motion_graph, sw_solutions, space_config.min_travel, print_timing, print_counter, 
      jt_start, jt_desired, sol_anim.trajectory, stop_planner);
  } else 
  if((space_config.space_order == 2) && (space_config.interp_id == 0)) {
    CRS_execute_static_planner_impl<kte::manip_P3R3R_kinematics, linear_interpolation_tag, 2>(ct_config.sceneData, plan_alg_config.planOptions,
      sw_motion_graph, sw_solutions, space_config.min_travel, print_timing, print_counter, 
      jt_start, jt_desired, sol_anim.trajectory, stop_planner);
  } else
#endif
      if ((space_config.space_order == 1) && (space_config.interp_id == 1)) {
    CRS_execute_static_planner_impl<kte::manip_P3R3R_kinematics,
                                    cubic_hermite_interpolation_tag, 1>(
        ct_config.sceneData, plan_alg_config.planOptions, sw_motion_graph,
        sw_solutions, space_config.min_travel, print_timing, print_counter,
        jt_start, jt_desired, sol_anim.trajectory, stop_planner);
  } else
#if 0
  if((space_config.space_order == 2) && (space_config.interp_id == 1)) {
    CRS_execute_static_planner_impl<kte::manip_P3R3R_kinematics, cubic_hermite_interpolation_tag, 2>(ct_config.sceneData, plan_alg_config.planOptions,
      sw_motion_graph, sw_solutions, space_config.min_travel, print_timing, print_counter, 
      jt_start, jt_desired, sol_anim.trajectory, stop_planner);
  } else
#endif
      if ((space_config.space_order == 2) && (space_config.interp_id == 2)) {
    CRS_execute_static_planner_impl<kte::manip_P3R3R_kinematics,
                                    quintic_hermite_interpolation_tag, 2>(
        ct_config.sceneData, plan_alg_config.planOptions, sw_motion_graph,
        sw_solutions, space_config.min_travel, print_timing, print_counter,
        jt_start, jt_desired, sol_anim.trajectory, stop_planner);
  } else if ((space_config.space_order == 1) && (space_config.interp_id == 3)) {
    CRS_execute_static_planner_impl<kte::manip_P3R3R_kinematics,
                                    svp_Ndof_interpolation_tag, 1>(
        ct_config.sceneData, plan_alg_config.planOptions, sw_motion_graph,
        sw_solutions, space_config.min_travel, print_timing, print_counter,
        jt_start, jt_desired, sol_anim.trajectory, stop_planner);
  } else
#if 0
  if((space_config.space_order == 2) && (space_config.interp_id == 3)) {
    CRS_execute_static_planner_impl<kte::manip_P3R3R_kinematics, svp_Ndof_interpolation_tag, 2>(ct_config.sceneData, plan_alg_config.planOptions,
      sw_motion_graph, sw_solutions, space_config.min_travel, print_timing, print_counter, 
      jt_start, jt_desired, sol_anim.trajectory, stop_planner);
  } else
#endif
      if ((space_config.space_order == 2) && (space_config.interp_id == 4)) {
    CRS_execute_static_planner_impl<kte::manip_P3R3R_kinematics,
                                    sap_Ndof_interpolation_tag, 2>(
        ct_config.sceneData, plan_alg_config.planOptions, sw_motion_graph,
        sw_solutions, space_config.min_travel, print_timing, print_counter,
        jt_start, jt_desired, sol_anim.trajectory, stop_planner);
  };
};


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

#include <ReaK/ctrl/topologies/manip_P3R3R_workspaces.hpp>

#include "CRS_planner_data.hpp"

#include <ReaK/ctrl/path_planning/path_planner_options.hpp>
#include <ReaK/ctrl/kte_models/chaser_target_model_data.hpp>

#include <ReaK/ctrl/interpolation/discrete_point_trajectory.hpp>
#include <ReaK/ctrl/interpolation/trajectory_base.hpp>
#include <ReaK/ctrl/topologies/manip_planning_traits.hpp>

#include <ReaK/ctrl/topologies/proxy_traj_applicator.hpp>
#include <ReaK/ctrl/topologies/direct_inverse_kin_topomap.hpp>
#include <ReaK/ctrl/interpolation/transformed_trajectory.hpp>





template <typename ManipMdlType, typename InterpTag, int Order, typename TargetStateTrajectory, typename ManipCSpaceTrajectory>
void CRS_execute_ee_control_impl(const ReaK::kte::chaser_target_data& scene_data, 
                                 const ReaK::shared_ptr< TargetStateTrajectory >& target_state_traj) {
  using namespace ReaK;
  using namespace pp;
  
  shared_ptr< ManipMdlType > chaser_concrete_model = rtti::rk_dynamic_ptr_cast<ManipMdlType>(scene_data.chaser_kin_model);
  if( !chaser_concrete_model ) {
    std::cout << "Could not cast the chaser model to the given type!" << std::endl;
    if( scene_data.chaser_kin_model )
      std::cout << "Because the chaser model is a null pointer!" << std::endl;
    else
      std::cout << "Because the RTTI could not perform the dynamic-cast!" << std::endl;
    return;
  };
  
  double current_time = 0.0; // TODO
  
  typedef typename spatial_trajectory_traits<TargetStateTrajectory>::space_topology TargetSpaceType;
  typedef typename spatial_trajectory_traits<TargetStateTrajectory>::point_type TargetTempPointType;
  
  TargetTempPointType target_tp = target_state_traj->get_point_at_time(current_time);
  target_state_traj->get_temporal_space().
  
  detail::write_joint_coordinates_impl(target_tp.pt, space_in, scene_data.target_kin_model);
  scene_data.target_kin_model->doDirectMotion();
  
  vect_n<double> target_jt_original = scene_data.target_kin_model->getJointPositions(); 
  
  typedef typename subspace_traits<dynamic_workspace_type>::super_space_type dynamic_super_space_type;  // SuperSpaceType
  
  std::size_t workspace_dims = (Order + 1) * manip_pp_traits< ManipMdlType, Order >::degrees_of_freedom;
  
  shared_ptr< frame_3D<double> > EE_frame = chaser_concrete_model->getDependentFrame3D(0)->mFrame;
  
  // Create the workspace for the chaser manipulator:
  shared_ptr< dynamic_workspace_type > workspace = 
    make_manip_dynamic_workspace<Order>(InterpTag(),
      chaser_concrete_model, scene_data.chaser_jt_limits, 
      min_travel, target_state_traj->get_end_time() - target_state_traj->get_start_time(), 
      plan_options.start_delay + target_state_traj->get_start_time());
  
  shared_ptr< rl_jt_space_type > jt_space = 
    make_manip_rl_jt_space<Order>(chaser_concrete_model, scene_data.chaser_jt_limits);
  
  shared_ptr< jt_space_type > normal_jt_space = 
    make_manip_jt_space<Order>(chaser_concrete_model, scene_data.chaser_jt_limits);
  
  
  // Add the proximity-query models for the target and the environment to the workspace's collision-checking:
  (*workspace) << scene_data.chaser_target_proxy;
  for(std::size_t i = 0; i < scene_data.chaser_env_proxies.size(); ++i)
    (*workspace) << scene_data.chaser_env_proxies[i];
  
  
  // Add the proximity-query model updater for the target's trajectory:
  typedef proxy_traj_applicator<TargetStateTrajectory> target_proxy_applicator_type;
  typedef typename target_proxy_applicator_type::joint_space_type target_jt_space_type;
  
  shared_ptr< target_proxy_applicator_type > target_proxy_applicator( 
    new target_proxy_applicator_type(
      make_any_model_applicator< target_jt_space_type >(manip_direct_kin_map(scene_data.target_kin_model)),
      target_state_traj));
  
  workspace->add_proxy_model_updater(target_proxy_applicator);
  
  
  // Create the starting configuration as a joint-space configuration.
  typedef typename topology_traits< dynamic_workspace_type >::point_type rl_temporal_point_type;
  
  point_type start_inter = normal_jt_space->origin();
  get<0>(start_inter) = jt_start;
  rl_point_type start_point = scene_data.chaser_jt_limits->map_to_space(start_inter, *normal_jt_space, *jt_space);
  
  rl_temporal_point_type temporal_start_point(plan_options.start_delay + target_state_traj->get_start_time(), start_point);
  
  
  // Create the reporter chain.
  any_sbmp_reporter_chain< dynamic_workspace_type > report_chain;
  
  // Create the frame tracing reporter.
  typedef frame_tracer_3D< dynamic_super_space_type > frame_reporter_type;
  
  frame_reporter_type temp_reporter(
    make_any_model_applicator< dynamic_super_space_type >( 
      rlDK_map_type(chaser_concrete_model, scene_data.chaser_jt_limits, normal_jt_space), 
      extract_spatial_component(), jt_space),
    0.5 * min_travel, (sw_motion_graph != NULL));
  
  if((sw_motion_graph) || (sw_solutions)) {
    temp_reporter.add_traced_frame(EE_frame);
    report_chain.add_reporter( boost::ref(temp_reporter) );
  };
  
  if( print_counter )
    report_chain.add_reporter( print_sbmp_progress<>() );
  
  if( print_timing )
    report_chain.add_reporter( timing_sbmp_report<>() );
  
  
  // Create the interception query object:
  
  typedef manip_dk_ik_map< manip_direct_kin_map, rlIK_map_type > dkik_map_type;
  typedef temporal_topo_map< dkik_map_type > temporal_dkik_map_type;
  typedef transformed_trajectory<dynamic_super_space_type, TargetStateTrajectory, temporal_dkik_map_type> mapped_traj_type;
  
  shared_ptr<mapped_traj_type> mapped_goal_traj( new mapped_traj_type(
    shared_ptr<dynamic_super_space_type>(&(workspace->get_super_space()), null_deleter()), 
    target_state_traj, 
    temporal_dkik_map_type(
      dkik_map_type(
        manip_direct_kin_map(scene_data.target_kin_model),
        rlIK_map_type(chaser_concrete_model, scene_data.chaser_jt_limits)
      ))));
  
  motion_plan_intercept_query< dynamic_workspace_type, mapped_traj_type > pp_query(
    "intercept_query", workspace, temporal_start_point, mapped_goal_traj,
    target_state_traj->get_end_time(), min_travel, plan_options.max_results);
  
  shared_ptr< sample_based_planner< dynamic_workspace_type > > workspace_planner;
  
  std::cout << "The selected planning algorithm is of index = " << plan_options.planning_algo << std::endl;
  switch(plan_options.planning_algo) {
    case 0:
      std::cout << "RRT" << std::endl;
      break;
    case 1:
      std::cout << "RRT*" << std::endl;
      break;
    case 2:
      std::cout << "PRM" << std::endl;
      break;
    case 3:
      std::cout << "SBA*" << std::endl;
      break;
    case 4:
      std::cout << "FADPRM" << std::endl;
      break;
    default:
      std::cout << "Unkown algorithm type" << std::endl;
      break;
  };
  
#if 0
  if( plan_options.planning_algo == 0 ) { // RRT
    
    workspace_planner = shared_ptr< sample_based_planner< dynamic_workspace_type > >(
      new rrt_planner< dynamic_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, report_chain));
    
  } else 
#endif
  if( plan_options.planning_algo == 1 ) { // RRT*
    
    workspace_planner = shared_ptr< sample_based_planner< dynamic_workspace_type > >(
      new rrtstar_planner< dynamic_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, workspace_dims + 1, report_chain));
    
  } else 
#if 0
  if( plan_options.planning_algo == 2 ) { // PRM
    
    workspace_planner = shared_ptr< sample_based_planner< dynamic_workspace_type > >(
      new prm_planner< dynamic_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, plan_options.max_random_walk, workspace_dims + 1, report_chain));
    
  } else 
#endif
  if( plan_options.planning_algo == 3 ) { // SBA*
    
    shared_ptr< sbastar_planner< dynamic_workspace_type > > tmp(
      new sbastar_planner< dynamic_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, plan_options.max_random_walk, workspace_dims + 1, report_chain));
    
    tmp->set_initial_density_threshold(0.0);
    tmp->set_initial_relaxation(plan_options.init_relax);
    tmp->set_initial_SA_temperature(plan_options.init_SA_temp);
    
    workspace_planner = tmp;
    
  } 
#if 0
  else if( plan_options.planning_algo == 4 ) { // FADPRM
    
    shared_ptr< fadprm_planner< dynamic_workspace_type > > tmp(
      new fadprm_planner< dynamic_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        0.1, 0.05, plan_options.max_random_walk, workspace_dims + 1, report_chain));
    
    tmp->set_initial_relaxation(plan_options.init_relax);
    
    workspace_planner = tmp;
    
  }
#endif
  else {
    std::cout << "The desired planner algorithm is not supported!" << std::endl;
  };
  
  
  if(!workspace_planner) {
    std::cout << "The workspace planner was not created successfully!" << std::endl;
    return;
  };
  
  pp_query.reset_solution_records();
  workspace_planner->solve_planning_query(pp_query);
  
  shared_ptr< seq_trajectory_base< dynamic_super_space_type > > bestsol_rltraj;
  if(pp_query.solutions.size())
    bestsol_rltraj = pp_query.solutions.begin()->second;
  std::cout << "The shortest distance is: " << pp_query.get_best_solution_distance() << std::endl;
  
  
  
  sol_trace.reset();
  if(bestsol_rltraj) {
    typedef hyperbox_topology< vect<double,7> > manip_cspace_type;
    typedef temporal_space<manip_cspace_type, time_poisson_topology, time_distance_only> temporal_space_type;
    typedef discrete_point_trajectory< temporal_space_type > trajectory_type;
    typedef trajectory_wrapper<trajectory_type> wrapped_traj_type;
    
    shared_ptr< wrapped_traj_type > tmp_sol_traj(new wrapped_traj_type("chaser_sol_trajectory", trajectory_type()));
    trajectory_type& sol_traj = tmp_sol_traj->get_underlying_trajectory();
    
    typedef typename seq_trajectory_base< dynamic_super_space_type >::point_time_iterator PtIter;
    typedef typename spatial_trajectory_traits<trajectory_type>::point_type TCSpacePointType;
    for(PtIter it = bestsol_rltraj->begin_time_travel(); it != bestsol_rltraj->end_time_travel(); it += 0.001) {
      rl_temporal_point_type cur_pt = *it;
      sol_traj.push_back( TCSpacePointType(cur_pt.time, get<0>(scene_data.chaser_jt_limits->map_to_space(cur_pt.pt, *jt_space, *normal_jt_space))) );
    };
    sol_trace = tmp_sol_traj;
  };
  
  
  
  // Check the motion-graph separator and solution separators
  //  add them to the switches.
  if(sw_motion_graph) {
    SoSeparator* mg_sep = temp_reporter.get_motion_graph_tracer(EE_frame).get_separator();
    if(mg_sep)
      mg_sep->ref();
    
    sw_motion_graph->removeAllChildren();
    if(mg_sep) {
      sw_motion_graph->addChild(mg_sep);
      mg_sep->unref();
    };
  };
  
  if(sw_solutions) {
    SoSeparator* sol_sep = NULL;
    if( temp_reporter.get_solution_count() ) {
      sol_sep = temp_reporter.get_solution_tracer(EE_frame, 0).get_separator();
      if(sol_sep)
        sol_sep->ref();
    };
    
    sw_solutions->removeAllChildren();
    if(sol_sep) {
      sw_solutions->addChild(sol_sep);
      sol_sep->unref();
    };
  };
  
  chaser_concrete_model->setJointPositions( jt_start );
  chaser_concrete_model->doDirectMotion();
  
  scene_data.target_kin_model->setJointPositions(target_jt_original);
  scene_data.target_kin_model->doDirectMotion();
  
};












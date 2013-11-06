
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

#include "CRS_planner2_impl.hpp"

#include "path_planning/path_planner_options.hpp"
#include "kte_models/chaser_target_model_data.hpp"


#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/sensors/SoTimerSensor.h>  // for SoTimerSensor

#include "shapes/oi_scene_graph.hpp"
#include "proximity/proxy_query_model.hpp"

#include "mbd_kte/kte_map_chain.hpp"
#include "kte_models/manip_dynamics_model.hpp"

#include "topologies/manip_P3R3R_workspaces.hpp"


#include "path_planning/intercept_query.hpp"

#include "path_planning/rrtstar_manip_planners.hpp"
#include "path_planning/sbastar_manip_planners.hpp"

#if 0
#include "path_planning/rrt_path_planner.hpp"
#include "path_planning/prm_path_planner.hpp"
#include "path_planning/fadprm_path_planner.hpp"
#endif

#include "CRS_planners_utility.hpp"

#include "CRS_planner_data.hpp"

#include "path_planning/frame_tracer_coin3d.hpp"

#include "optimization/optim_exceptions.hpp"

#include "topologies/manip_planning_traits.hpp"

#include "topologies/proxy_traj_applicator.hpp"
#include "topologies/direct_inverse_kin_topomap.hpp"

#include "topologies/Ndof_linear_spaces.hpp"
#include "topologies/Ndof_cubic_spaces.hpp"
#include "topologies/Ndof_quintic_spaces.hpp"
#include "topologies/Ndof_svp_spaces.hpp"
#include "topologies/Ndof_sap_spaces.hpp"

#include "path_planning/transformed_trajectory.hpp"

#include <chrono>


using namespace ReaK;





template <typename InterpTag, int Order, typename TargetStateTrajectory>
void CRS_execute_dynamic_planner_impl(const kte::chaser_target_data& scene_data, 
                                      const pp::planning_option_collection& plan_options,
                                      CRS_coin_nodes* draw_data,
                                      const vect_n<double>& jt_start, double start_delay,
                                      const shared_ptr< TargetStateTrajectory >& target_state_traj,
                                      std::vector< vect<double,7> >& sol_trace) {
  
  shared_ptr< kte::manip_P3R3R_kinematics > chaser_P3R3R_model = rtti::rk_dynamic_ptr_cast<kte::manip_P3R3R_kinematics>(scene_data.chaser_kin_model);
  if( !chaser_P3R3R_model )
    return;
  
  typedef typename pp::manip_dynamic_workspace< kte::manip_P3R3R_kinematics, Order >::rl_workspace_type dynamic_workspace_type;
  typedef typename pp::manip_pp_traits< kte::manip_P3R3R_kinematics, Order >::rl_jt_space_type rl_jt_space_type;
  typedef typename pp::manip_pp_traits< kte::manip_P3R3R_kinematics, Order >::jt_space_type jt_space_type;
  typedef typename pp::manip_pp_traits< kte::manip_P3R3R_kinematics, Order >::ee_space_type ee_space_type;
  typedef typename pp::manip_DK_map< kte::manip_P3R3R_kinematics, Order >::rl_map_type rlDK_map_type;
  typedef typename pp::manip_IK_map< kte::manip_P3R3R_kinematics, Order >::rl_map_type rlIK_map_type;
  
  typedef typename pp::topology_traits< rl_jt_space_type >::point_type rl_point_type;
  typedef typename pp::topology_traits< jt_space_type >::point_type point_type;
  
  typedef typename pp::subspace_traits<dynamic_workspace_type>::super_space_type dynamic_super_space_type;  // SuperSpaceType
  
  std::size_t workspace_dims = Order * pp::manip_pp_traits< kte::manip_P3R3R_kinematics, Order >::degrees_of_freedom;
  
  shared_ptr< frame_3D<double> > EE_frame = chaser_P3R3R_model->getDependentFrame3D(0)->mFrame;
  
  // Create the workspace for the chaser manipulator:
  shared_ptr< dynamic_workspace_type > workspace = 
    pp::make_manip_dynamic_workspace<Order>(InterpTag(),
      chaser_P3R3R_model, scene_data.chaser_jt_limits, 
      plan_options.min_travel, plan_options.max_travel);
  
  shared_ptr< rl_jt_space_type > jt_space = 
    pp::make_manip_rl_jt_space<Order>(chaser_P3R3R_model, scene_data.chaser_jt_limits);
  
  shared_ptr< jt_space_type > normal_jt_space = 
    pp::make_manip_jt_space<Order>(chaser_P3R3R_model, scene_data.chaser_jt_limits);
  
  
  // Add the proximity-query models for the target and the environment to the workspace's collision-checking:
  (*workspace) << scene_data.chaser_target_proxy;
  for(std::size_t i = 0; i < scene_data.chaser_env_proxies.size(); ++i)
    (*workspace) << scene_data.chaser_env_proxies[i];
  
  
  // Add the proximity-query model updater for the target's trajectory:
  typedef pp::proxy_traj_applicator<TargetStateTrajectory> target_proxy_applicator_type;
  typedef typename target_proxy_applicator_type::joint_space_type target_jt_space_type;
  
  shared_ptr< target_proxy_applicator_type > target_proxy_applicator( 
    new target_proxy_applicator_type(
      pp::make_any_model_applicator< target_jt_space_type >(pp::manip_direct_kin_map(scene_data.target_kin_model)),
      target_state_traj));
  
  workspace->add_proxy_model_updater(target_proxy_applicator);
  
  
  // Create the starting configuration as a joint-space configuration.
  typedef typename pp::topology_traits< dynamic_workspace_type >::point_type rl_temporal_point_type;
  
  point_type start_inter = normal_jt_space->origin();
  get<0>(start_inter) = jt_start;
  rl_point_type start_point = scene_data.chaser_jt_limits->map_to_space(start_inter, *normal_jt_space, *jt_space);
  
  rl_temporal_point_type temporal_start_point(start_delay, start_point);
  
  
  // Create the frame tracing reporter.
  typedef pp::frame_tracer_3D< dynamic_super_space_type, pp::print_sbmp_progress<> > frame_reporter_type;
  
  frame_reporter_type temp_reporter(
    pp::make_any_model_applicator< dynamic_super_space_type >( 
      rlDK_map_type(chaser_P3R3R_model, scene_data.chaser_jt_limits, normal_jt_space), 
      pp::extract_spatial_component(), jt_space),
    0.5 * plan_options.min_travel);
  temp_reporter.add_traced_frame(EE_frame);
  
  pp::any_sbmp_reporter_chain< dynamic_workspace_type > report_chain;
  report_chain.add_reporter( boost::ref(temp_reporter) );
  
  
  // Create the interception query object:
  
  typedef pp::manip_dk_ik_map< pp::manip_direct_kin_map, rlIK_map_type > dkik_map_type;
  typedef pp::temporal_topo_map< dkik_map_type > temporal_dkik_map_type;
  typedef pp::transformed_trajectory<dynamic_super_space_type, TargetStateTrajectory, temporal_dkik_map_type> mapped_traj_type;
  
  shared_ptr<mapped_traj_type> mapped_goal_traj( new mapped_traj_type(
    shared_ptr<dynamic_super_space_type>(&(workspace->get_super_space()), null_deleter()), 
    target_state_traj, 
    temporal_dkik_map_type(
      dkik_map_type(
        pp::manip_direct_kin_map(scene_data.target_kin_model),
        rlIK_map_type(chaser_P3R3R_model, scene_data.chaser_jt_limits)
      ))));
  
  pp::motion_plan_intercept_query< dynamic_workspace_type, mapped_traj_type > pp_query(
    "intercept_query", workspace, temporal_start_point, mapped_goal_traj,
    target_state_traj->get_end_time(), plan_options.min_travel, plan_options.max_results);
  
  shared_ptr< pp::sample_based_planner< dynamic_workspace_type > > workspace_planner;
  
#if 0
  if( plan_options.planning_algo == 0 ) { // RRT
    
    workspace_planner = shared_ptr< pp::sample_based_planner< dynamic_workspace_type > >(
      new pp::rrt_planner< dynamic_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, report_chain));
    
  } else 
#endif
  if( plan_options.planning_algo == 1 ) { // RRT*
    
    workspace_planner = shared_ptr< pp::sample_based_planner< dynamic_workspace_type > >(
      new pp::rrtstar_planner< dynamic_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, workspace_dims + 1, report_chain));
    
  } else 
#if 0
  if( plan_options.planning_algo == 2 ) { // PRM
    
    workspace_planner = shared_ptr< pp::sample_based_planner< dynamic_workspace_type > >(
      new pp::prm_planner< dynamic_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, plan_options.max_travel, workspace_dims + 1, report_chain));
    
  } else 
#endif
  if( plan_options.planning_algo == 3 ) { // SBA*
    
    shared_ptr< pp::sbastar_planner< dynamic_workspace_type > > tmp(
      new pp::sbastar_planner< dynamic_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, plan_options.max_travel, workspace_dims + 1, report_chain));
    
    tmp->set_initial_density_threshold(0.0);
    tmp->set_initial_relaxation(plan_options.init_relax);
    tmp->set_initial_SA_temperature(plan_options.init_SA_temp);
    
    workspace_planner = tmp;
    
  } 
#if 0
  else if( plan_options.planning_algo == 4 ) { // FADPRM
    
    shared_ptr< pp::fadprm_planner< dynamic_workspace_type > > tmp(
      new pp::fadprm_planner< dynamic_workspace_type >(
        workspace, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        0.1, 0.05, plan_options.max_travel, workspace_dims + 1, report_chain));
    
    tmp->set_initial_relaxation(plan_options.init_relax);
    
    workspace_planner = tmp;
    
  }
#endif
  ;
  
  
  if(!workspace_planner)
    return;
  
  pp_query.reset_solution_records();
  workspace_planner->solve_planning_query(pp_query);
  
  shared_ptr< pp::seq_trajectory_base< dynamic_super_space_type > > bestsol_rltraj;
  if(pp_query.solutions.size())
    bestsol_rltraj = pp_query.solutions.begin()->second;
  std::cout << "The shortest distance is: " << pp_query.get_best_solution_distance() << std::endl;
  
  sol_trace.clear();
  if(bestsol_rltraj) {
    typedef typename pp::seq_trajectory_base< dynamic_super_space_type >::point_time_iterator PtIter;
    for(PtIter it = bestsol_rltraj->begin_time_travel(); it != bestsol_rltraj->end_time_travel(); it += 0.05)
      sol_trace.push_back( vect<double,7>(get<0>(scene_data.chaser_jt_limits->map_to_space((*it).pt, *jt_space, *normal_jt_space))) );
  };
  
  SoSeparator* mg_sep = temp_reporter.get_motion_graph_tracer(EE_frame).get_separator();
  mg_sep->ref();
  SoSeparator* sol_sep = NULL;
  if( temp_reporter.get_solution_count() ) {
    sol_sep = temp_reporter.get_solution_tracer(EE_frame, 0).get_separator();
    sol_sep->ref();
  };
  
  chaser_P3R3R_model->setJointPositions( vect_n<double>( get<0>(start_inter) ) );
  chaser_P3R3R_model->doDirectMotion();
  
  
  // Check the motion-graph separator and solution separators
  //  add them to the switches.
  if(mg_sep) {
    draw_data->sw_motion_graph->removeAllChildren();
    draw_data->sw_motion_graph->addChild(mg_sep);
    mg_sep->unref();
  };
  
  draw_data->sw_solutions->removeAllChildren();
  if(sol_sep) {
    draw_data->sw_solutions->addChild(sol_sep);
    sol_sep->unref();
  };
  
};






void CRSPlannerGUI::executeDynamicPlanner() {
  
  vect_n<double> jt_start;
  if(configs.check_current_start->isChecked()) {
    jt_start = scene_data.chaser_kin_model->getJointPositions(); 
  } else {
    std::stringstream ss(configs.custom_start_edit->text().toStdString());
    ss >> jt_start;
  };
  
  // update the planning options record:
  onConfigsChanged();
  
  if((plan_options.space_order == 0) && (plan_options.interp_id == 0)) { 
    CRS_execute_dynamic_planner_impl<pp::linear_interpolation_tag, 0>(
      scene_data, plan_options, draw_data, jt_start, 10.0, target_anim->target_trajectory, sol_anim->bestsol_trajectory);
  } else 
#if 0
  if((plan_options.space_order == 1) && (plan_options.interp_id == 0)) {
    CRS_execute_dynamic_planner_impl<pp::linear_interpolation_tag, 1>(
      scene_data, plan_options, draw_data, jt_start, 10.0, target_anim->target_trajectory, sol_anim->bestsol_trajectory);
  } else 
  if((plan_options.space_order == 2) && (plan_options.interp_id == 0)) {
    CRS_execute_dynamic_planner_impl<pp::linear_interpolation_tag, 2>(
      scene_data, plan_options, draw_data, jt_start, 10.0, target_anim->target_trajectory, sol_anim->bestsol_trajectory);
  } else 
#endif
  if((plan_options.space_order == 1) && (plan_options.interp_id == 1)) {
    CRS_execute_dynamic_planner_impl<pp::cubic_hermite_interpolation_tag, 1>(
      scene_data, plan_options, draw_data, jt_start, 10.0, target_anim->target_trajectory, sol_anim->bestsol_trajectory);
  } else 
#if 0
  if((plan_options.space_order == 2) && (plan_options.interp_id == 1)) {
    CRS_execute_dynamic_planner_impl<pp::cubic_hermite_interpolation_tag, 2>(
      scene_data, plan_options, draw_data, jt_start, 10.0, target_anim->target_trajectory, sol_anim->bestsol_trajectory);
  } else 
#endif
  if((plan_options.space_order == 2) && (plan_options.interp_id == 2)) {
    CRS_execute_dynamic_planner_impl<pp::quintic_hermite_interpolation_tag, 2>(
      scene_data, plan_options, draw_data, jt_start, 10.0, target_anim->target_trajectory, sol_anim->bestsol_trajectory);
  } else 
  if((plan_options.space_order == 1) && (plan_options.interp_id == 3)) {
    CRS_execute_dynamic_planner_impl<pp::svp_Ndof_interpolation_tag, 1>(
      scene_data, plan_options, draw_data, jt_start, 10.0, target_anim->target_trajectory, sol_anim->bestsol_trajectory);
  } else 
#if 0
  if((plan_options.space_order == 2) && (plan_options.interp_id == 3)) {
    CRS_execute_static_planner_impl<pp::svp_Ndof_interpolation_tag, 2>(scene_data, plan_options,draw_data, jt_start, jt_desired, sol_anim->bestsol_trajectory);
  } else 
#endif
  if((plan_options.space_order == 2) && (plan_options.interp_id == 4)) {
    CRS_execute_dynamic_planner_impl<pp::sap_Ndof_interpolation_tag, 2>(
      scene_data, plan_options, draw_data, jt_start, 10.0, target_anim->target_trajectory, sol_anim->bestsol_trajectory);
  };
  
  
};











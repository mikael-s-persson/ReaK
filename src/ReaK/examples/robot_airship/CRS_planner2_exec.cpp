
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


#include <QApplication>
#include <QMessageBox>
#include <QFileDialog>
#include <QMainWindow>
#include <QDir>


#include "CRS_A465_geom_model.hpp"

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

#include "CRS_workspaces.hpp"
#include "CRS_rrtstar_planners.hpp"
#include "CRS_sbastar_planners.hpp"

#if 0
#include "CRS_rrt_planners.hpp"
#include "CRS_prm_planners.hpp"
#include "CRS_fadprm_planners.hpp"
#endif

#include "CRS_planners_utility.hpp"

#include "CRS_planner_data.hpp"

#include "path_planning/frame_tracer_coin3d.hpp"

#include "optimization/optim_exceptions.hpp"

#include "topologies/manip_planning_traits.hpp"


#include <chrono>


using namespace ReaK;





template <typename InterpTag, int Order>
void CRS_execute_static_planner_impl(robot_airship::scenario_data* scene_data, 
                                     CRS_planning_options* plan_options,
                                     CRS_coin_nodes* draw_data,
                                     const vect_n<double>& jt_start, 
                                     const vect_n<double>& jt_desired,
                                     std::vector< vect<double,7> >& sol_trace
                                    ) {
  
  typedef typename pp::manip_static_workspace< kte::manip_P3R3R_kinematics, Order, InterpTag >::rl_workspace_type static_workspace_type;
  typedef typename pp::manip_pp_traits< kte::manip_P3R3R_kinematics, Order >::rl_jt_space_type rl_jt_space_type;
  typedef typename pp::manip_pp_traits< kte::manip_P3R3R_kinematics, Order >::jt_space_type jt_space_type;
  typedef typename pp::manip_pp_traits< kte::manip_P3R3R_kinematics, Order >::ee_space_type ee_space_type;
  typedef typename pp::manip_DK_map< kte::manip_P3R3R_kinematics, Order >::rl_map_type rlDK_map_type;
  
  typedef typename pp::topology_traits< rl_jt_space_type >::point_type rl_point_type;
  typedef typename pp::topology_traits< jt_space_type >::point_type point_type;
  
  typedef typename pp::subspace_traits<static_workspace_type>::super_space_type static_super_space_type;  // SuperSpaceType
  
  typedef pp::frame_tracer_3D< rlDK_map_type, rl_jt_space_type, pp::identity_topo_map, pp::print_sbmp_progress<> > frame_reporter_type;
  
  std::size_t workspace_dims = Order * pp::manip_pp_traits< kte::manip_P3R3R_kinematics, Order >::degrees_of_freedom;
  
  shared_ptr< frame_3D<double> > EE_frame = scene_data->chaser_kin_model->getDependentFrame3D(0)->mFrame;
  
  shared_ptr< static_workspace_type > workspace = 
    pp::make_manip_static_workspace<Order, InterpTag>(
      scene_data->chaser_kin_model, scene_data->chaser_jt_limits, 
      plan_options->min_travel, plan_options->max_travel);
  
  shared_ptr< rl_jt_space_type > jt_space = 
    pp::make_manip_rl_jt_space<Order>(scene_data->chaser_kin_model, scene_data->chaser_jt_limits);
  
  shared_ptr< jt_space_type > normal_jt_space = 
    pp::make_manip_jt_space<Order>(scene_data->chaser_kin_model, scene_data->chaser_jt_limits);
  
  
  (*workspace) << scene_data->chaser_target_proxy;
  for(std::size_t i = 0; i < scene_data->chaser_env_proxies.size(); ++i)
    (*workspace) << scene_data->chaser_env_proxies[i];
  
  
  rl_point_type start_point, goal_point;
  point_type start_inter, goal_inter;
  
  start_inter = normal_jt_space->origin();
  get<0>(start_inter) = jt_start;
  start_point = scene_data->chaser_jt_limits->map_to_space(start_inter, *normal_jt_space, *jt_space);
  
  goal_inter = normal_jt_space->origin();
  get<0>(goal_inter) = jt_desired;
  goal_point = scene_data->chaser_jt_limits->map_to_space(goal_inter, *normal_jt_space, *jt_space);
  
  
  frame_reporter_type temp_reporter(
    rlDK_map_type(scene_data->chaser_kin_model, scene_data->chaser_jt_limits, normal_jt_space),
    jt_space, 0.5 * plan_options->min_travel);
  temp_reporter.add_traced_frame(EE_frame);
  
  pp::any_sbmp_reporter_chain< static_workspace_type > report_chain;
  report_chain.add_reporter( boost::ref(temp_reporter) );
  
  pp::path_planning_p2p_query< static_workspace_type > pp_query("pp_query", workspace,
    start_point, goal_point, plan_options->max_results);
  
  
  shared_ptr< pp::sample_based_planner< static_workspace_type > > workspace_planner;
  
#if 0
  if( plan_options->planning_algo == 0 ) { // RRT
    
    workspace_planner = shared_ptr< pp::sample_based_planner< static_workspace_type > >(
      new pp::rrt_planner< static_workspace_type >(
        workspace, plan_options->max_vertices, plan_options->prog_interval,
        plan_options->store_policy | plan_options->knn_method,
        plan_options->planning_options,
        0.1, 0.05, report_chain));
    
  } else 
#endif
  if( plan_options->planning_algo == 1 ) { // RRT*
    
    workspace_planner = shared_ptr< pp::sample_based_planner< static_workspace_type > >(
      new pp::rrtstar_planner< static_workspace_type >(
        workspace, plan_options->max_vertices, plan_options->prog_interval,
        plan_options->store_policy | plan_options->knn_method,
        plan_options->planning_options,
        0.1, 0.05, workspace_dims, report_chain));
    
  } else 
#if 0
  if( plan_options->planning_algo == 2 ) { // PRM
    
    workspace_planner = shared_ptr< pp::sample_based_planner< static_workspace_type > >(
      new pp::prm_planner< static_workspace_type >(
        workspace, plan_options->max_vertices, plan_options->prog_interval,
        plan_options->store_policy | plan_options->knn_method,
        plan_options->planning_options,
        0.1, 0.05, plan_options->max_travel, workspace_dims, report_chain));
    
  } else 
#endif
  if( plan_options->planning_algo == 3 ) { // SBA*
    
    shared_ptr< pp::sbastar_planner< static_workspace_type > > tmp(
      new pp::sbastar_planner< static_workspace_type >(
        workspace, plan_options->max_vertices, plan_options->prog_interval,
        plan_options->store_policy | plan_options->knn_method,
        plan_options->planning_options,
        0.1, 0.05, plan_options->max_travel, workspace_dims, report_chain));
    
    tmp->set_initial_density_threshold(0.0);
    tmp->set_initial_relaxation(plan_options->init_relax);
    tmp->set_initial_SA_temperature(plan_options->init_SA_temp);
    
    workspace_planner = tmp;
    
  } 
#if 0
  else if( plan_options->planning_algo == 4 ) { // FADPRM
    
    shared_ptr< pp::fadprm_planner< static_workspace_type > > tmp(
      new pp::fadprm_planner< static_workspace_type >(
        workspace, plan_options->max_vertices, plan_options->prog_interval,
        plan_options->store_policy | plan_options->knn_method,
        0.1, 0.05, plan_options->max_travel, workspace_dims, report_chain));
    
    tmp->set_initial_relaxation(plan_options->init_relax);
    
    workspace_planner = tmp;
    
  }
#endif
  ;
  
  
  if(!workspace_planner)
    return;
  
  pp_query.reset_solution_records();
  workspace_planner->solve_planning_query(pp_query);
  
  shared_ptr< pp::seq_path_base< static_super_space_type > > bestsol_rlpath;
  if(pp_query.solutions.size())
    bestsol_rlpath = pp_query.solutions.begin()->second;
  std::cout << "The shortest distance is: " << pp_query.get_best_solution_distance() << std::endl;
  
  sol_trace.clear();
  if(bestsol_rlpath) {
    typedef typename pp::seq_path_base< static_super_space_type >::point_fraction_iterator PtIter;
    for(PtIter it = bestsol_rlpath->begin_fraction_travel(); it != bestsol_rlpath->end_fraction_travel(); it += 0.1)
      sol_trace.push_back( vect<double,7>(get<0>(scene_data->chaser_jt_limits->map_to_space(*it, *jt_space, *normal_jt_space))) );
  };
  
  SoSeparator* mg_sep = temp_reporter.get_motion_graph_tracer(EE_frame).get_separator();
  mg_sep->ref();
  SoSeparator* sol_sep = NULL;
  if( temp_reporter.get_solution_count() ) {
    sol_sep = temp_reporter.get_solution_tracer(EE_frame, 0).get_separator();
    sol_sep->ref();
  };
  
  scene_data->chaser_kin_model->setJointPositions( vect_n<double>( get<0>(start_inter) ) );
  scene_data->chaser_kin_model->doDirectMotion();
  
  
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









void CRSPlannerGUI::executePlanner() {
  
  shared_ptr< frame_3D<double> > EE_frame = scene_data->chaser_kin_model->getDependentFrame3D(0)->mFrame;
  
  vect_n<double> jt_desired(7,0.0);
  if(configs.check_ik_goal->isChecked()) {
    vect_n<double> jt_previous = scene_data->chaser_kin_model->getJointPositions();
    try {
      frame_3D<double> tf = scene_data->target_frame->getFrameRelativeTo(EE_frame);
      EE_frame->addBefore(tf);
      scene_data->chaser_kin_model->doInverseMotion();
      jt_desired = scene_data->chaser_kin_model->getJointPositions();
    } catch( optim::infeasible_problem& e ) { RK_UNUSED(e);
      QMessageBox::critical(this,
                    "Inverse Kinematics Error!",
                    "The target frame cannot be reached! No inverse kinematics solution possible!",
                    QMessageBox::Ok);
      return;
    };
    scene_data->chaser_kin_model->setJointPositions(jt_previous);
    scene_data->chaser_kin_model->doDirectMotion();
  } else {
    std::stringstream ss(configs.custom_goal_edit->text().toStdString());
    ss >> jt_desired;
  };
  
  vect_n<double> jt_start;
  if(configs.check_current_start->isChecked()) {
    jt_start = scene_data->chaser_kin_model->getJointPositions(); 
  } else {
    std::stringstream ss(configs.custom_start_edit->text().toStdString());
    ss >> jt_start;
  };
  
  
  // update the planning options record:
  onConfigsChanged();
  
  
  if((plan_options->space_order == 0) && (plan_options->interp_id == 0)) { 
    CRS_execute_static_planner_impl<pp::linear_interpolation_tag, 0>(scene_data, plan_options,draw_data, jt_start, jt_desired, sol_anim->bestsol_trajectory);
  } else 
#if 0
  if((plan_options->space_order == 1) && (plan_options->interp_id == 0)) {
    CRS_execute_static_planner_impl<pp::linear_interpolation_tag, 1>(scene_data, plan_options,draw_data, jt_start, jt_desired, sol_anim->bestsol_trajectory);
  } else 
  if((plan_options->space_order == 2) && (plan_options->interp_id == 0)) {
    CRS_execute_static_planner_impl<pp::linear_interpolation_tag, 2>(scene_data, plan_options,draw_data, jt_start, jt_desired, sol_anim->bestsol_trajectory);
  } else 
#endif
  if((plan_options->space_order == 1) && (plan_options->interp_id == 1)) {
    CRS_execute_static_planner_impl<pp::cubic_hermite_interpolation_tag, 1>(scene_data, plan_options,draw_data, jt_start, jt_desired, sol_anim->bestsol_trajectory);
  } else 
#if 0
  if((plan_options->space_order == 2) && (plan_options->interp_id == 1)) {
    CRS_execute_static_planner_impl<pp::cubic_hermite_interpolation_tag, 2>(scene_data, plan_options,draw_data, jt_start, jt_desired, sol_anim->bestsol_trajectory);
  } else 
#endif
  if((plan_options->space_order == 2) && (plan_options->interp_id == 2)) {
    CRS_execute_static_planner_impl<pp::quintic_hermite_interpolation_tag, 2>(scene_data, plan_options,draw_data, jt_start, jt_desired, sol_anim->bestsol_trajectory);
  } else 
  if((plan_options->space_order == 1) && (plan_options->interp_id == 3)) {
    CRS_execute_static_planner_impl<pp::svp_Ndof_interpolation_tag, 1>(scene_data, plan_options,draw_data, jt_start, jt_desired, sol_anim->bestsol_trajectory);
  } else 
#if 0
  if((plan_options->space_order == 2) && (plan_options->interp_id == 3)) {
    CRS_execute_static_planner_impl<pp::svp_Ndof_interpolation_tag, 2>(scene_data, plan_options,draw_data, jt_start, jt_desired, sol_anim->bestsol_trajectory);
  } else 
#endif
  if((plan_options->space_order == 2) && (plan_options->interp_id == 4)) {
    CRS_execute_static_planner_impl<pp::sap_Ndof_interpolation_tag, 2>(scene_data, plan_options,draw_data, jt_start, jt_desired, sol_anim->bestsol_trajectory);
  };
  
  
};





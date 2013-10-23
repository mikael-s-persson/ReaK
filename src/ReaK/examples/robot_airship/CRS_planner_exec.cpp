
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

#include "CRS_planner_impl.h"


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

#include "CRS_workspaces.hpp"
#include "CRS_rrtstar_planners.hpp"
#include "CRS_sbastar_planners.hpp"

#if 0
#include "CRS_rrt_planners.hpp"
#include "CRS_prm_planners.hpp"
#include "CRS_fadprm_planners.hpp"
#endif

#include "CRS_planner_data.hpp"

#include "path_planning/frame_tracer_coin3d.hpp"

#include "optimization/optim_exceptions.hpp"


#include <chrono>


using namespace ReaK;



void CRSPlannerGUI::executePlanner() {
  
  vect_n<double> jt_desired(7,0.0);
  if(configs.check_ik_goal->isChecked()) {
    try {
      
      jt_desired = mdl_data->builder.compute_inverse_kinematics(mdl_data->target_frame.getGlobalPose());
      
    } catch( optim::infeasible_problem& e ) { RK_UNUSED(e);
      QMessageBox::critical(this,
                    "Inverse Kinematics Error!",
                    "The target frame cannot be reached! No inverse kinematics solution possible!",
                    QMessageBox::Ok);
      return;
    };
  } else {
    std::stringstream ss(configs.custom_goal_edit->text().toStdString());
    ss >> jt_desired;
  };
  
  vect_n<double> jt_start;
  if(configs.check_current_start->isChecked()) {
    jt_start.resize(7);
    jt_start[0] = mdl_data->builder.track_joint_coord->q;
    jt_start[1] = mdl_data->builder.arm_joint_1_coord->q;
    jt_start[2] = mdl_data->builder.arm_joint_2_coord->q;
    jt_start[3] = mdl_data->builder.arm_joint_3_coord->q;
    jt_start[4] = mdl_data->builder.arm_joint_4_coord->q; 
    jt_start[5] = mdl_data->builder.arm_joint_5_coord->q; 
    jt_start[6] = mdl_data->builder.arm_joint_6_coord->q; 
  } else {
    std::stringstream ss(configs.custom_start_edit->text().toStdString());
    ss >> jt_start;
  };
  
  
  // update the planning options record:
  onConfigsChanged();
  
  
  SoSeparator* mg_sep = NULL;
  std::vector< SoSeparator* > sol_seps;
  
  
  /*
   * Below are a few rather large MACROs that are used to generate all the code for the different planner-space
   * combinations.
   */
  
#define RK_CRS_PLANNER_GENERATE_PLANNER_IC_0(WORKSPACE) \
        typedef pp::subspace_traits<WORKSPACE>::super_space_type SuperSpaceType; \
        shared_ptr< robot_airship::CRS3D_jspace_rl_o0_type > jt_space(new robot_airship::CRS3D_jspace_rl_o0_type(mdl_data->builder.get_rl_joint_space_0th())); \
        shared_ptr< robot_airship::CRS3D_jspace_o0_type > normal_jt_space(new robot_airship::CRS3D_jspace_o0_type(mdl_data->builder.get_joint_space_0th())); \
         \
        shared_ptr<WORKSPACE>  \
          workspace(new WORKSPACE( \
            *jt_space, \
            mdl_data->manip_kin_mdl, \
            mdl_data->manip_jt_limits, \
            plan_options->min_travel, \
            plan_options->max_travel)); \
        std::size_t workspace_dims = 7; RK_UNUSED(workspace_dims); \
         \
        (*workspace) << mdl_data->robot_lab_proxy << mdl_data->robot_airship_proxy; \
         \
        typedef pp::topology_traits< robot_airship::CRS3D_jspace_rl_o0_type >::point_type RLPointType; \
        typedef pp::topology_traits< robot_airship::CRS3D_jspace_o0_type >::point_type PointType; \
        RLPointType start_point, goal_point; \
        PointType start_inter, goal_inter; \
        start_inter = normal_jt_space->origin(); \
        get<0>(start_inter) = vect<double,7>(jt_start); \
        start_point = mdl_data->manip_jt_limits->map_to_space(start_inter, *normal_jt_space, *jt_space); \
         \
        goal_inter = normal_jt_space->origin(); \
        get<0>(goal_inter) = vect<double,7>(jt_desired); \
        goal_point = mdl_data->manip_jt_limits->map_to_space(goal_inter, *normal_jt_space, *jt_space); \
         \
        typedef pp::frame_tracer_3D< robot_airship::CRS3D_rlDK_o0_type, robot_airship::CRS3D_jspace_rl_o0_type, pp::identity_topo_map, pp::print_sbmp_progress<> > frame_reporter_type; \
        frame_reporter_type temp_reporter( \
          robot_airship::CRS3D_rlDK_o0_type(mdl_data->manip_kin_mdl, mdl_data->manip_jt_limits, normal_jt_space), \
          jt_space, 0.5 * plan_options->min_travel); \
        temp_reporter.add_traced_frame(mdl_data->builder.arm_joint_6_end); \
         \
        pp::any_sbmp_reporter_chain< WORKSPACE > report_chain; \
        report_chain.add_reporter( boost::ref(temp_reporter) ); \
         \
        pp::path_planning_p2p_query< WORKSPACE > pp_query("pp_query", workspace, \
          start_point, goal_point, plan_options->max_results);
        
#define RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(WORKSPACE) \
        typedef pp::subspace_traits<WORKSPACE>::super_space_type SuperSpaceType; \
        shared_ptr< robot_airship::CRS3D_jspace_rl_o1_type > jt_space(new robot_airship::CRS3D_jspace_rl_o1_type(mdl_data->builder.get_rl_joint_space_1st())); \
        shared_ptr< robot_airship::CRS3D_jspace_o1_type > normal_jt_space(new robot_airship::CRS3D_jspace_o1_type(mdl_data->builder.get_joint_space_1st())); \
         \
        shared_ptr<WORKSPACE>  \
          workspace(new WORKSPACE( \
            *jt_space, \
            mdl_data->manip_kin_mdl, \
            mdl_data->manip_jt_limits, \
            plan_options->min_travel, \
            plan_options->max_travel)); \
        std::size_t workspace_dims = 14; RK_UNUSED(workspace_dims); \
         \
        (*workspace) << mdl_data->robot_lab_proxy << mdl_data->robot_airship_proxy; \
         \
        typedef pp::topology_traits< robot_airship::CRS3D_jspace_rl_o1_type >::point_type RLPointType; \
        typedef pp::topology_traits< robot_airship::CRS3D_jspace_o1_type >::point_type PointType; \
        RLPointType start_point, goal_point; \
        PointType start_inter, goal_inter; \
        start_inter = normal_jt_space->origin(); \
        get<0>(start_inter) = vect<double,7>(jt_start); \
        start_point = mdl_data->manip_jt_limits->map_to_space(start_inter, *normal_jt_space, *jt_space); \
         \
        goal_inter = normal_jt_space->origin(); \
        get<0>(goal_inter) = vect<double,7>(jt_desired); \
        goal_point = mdl_data->manip_jt_limits->map_to_space(goal_inter, *normal_jt_space, *jt_space); \
         \
        typedef pp::frame_tracer_3D< robot_airship::CRS3D_rlDK_o1_type, robot_airship::CRS3D_jspace_rl_o1_type, pp::identity_topo_map, pp::print_sbmp_progress<> > frame_reporter_type; \
        frame_reporter_type temp_reporter( \
          robot_airship::CRS3D_rlDK_o1_type(mdl_data->manip_kin_mdl, mdl_data->manip_jt_limits, normal_jt_space), \
          jt_space, 0.5 * plan_options->min_travel); \
        temp_reporter.add_traced_frame(mdl_data->builder.arm_joint_6_end); \
         \
        pp::any_sbmp_reporter_chain< WORKSPACE > report_chain; \
        report_chain.add_reporter( boost::ref(temp_reporter) ); \
         \
        pp::path_planning_p2p_query< WORKSPACE > pp_query("pp_query", workspace, \
          start_point, goal_point, plan_options->max_results);
        
#define RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(WORKSPACE) \
        typedef pp::subspace_traits<WORKSPACE>::super_space_type SuperSpaceType; \
        shared_ptr< robot_airship::CRS3D_jspace_rl_o2_type > jt_space(new robot_airship::CRS3D_jspace_rl_o2_type(mdl_data->builder.get_rl_joint_space())); \
        shared_ptr< robot_airship::CRS3D_jspace_o2_type > normal_jt_space(new robot_airship::CRS3D_jspace_o2_type(mdl_data->builder.get_joint_space())); \
         \
        shared_ptr<WORKSPACE>  \
          workspace(new WORKSPACE( \
            *jt_space, \
            mdl_data->manip_kin_mdl, \
            mdl_data->manip_jt_limits, \
            plan_options->min_travel, \
            plan_options->max_travel)); \
        std::size_t workspace_dims = 21; RK_UNUSED(workspace_dims); \
         \
        (*workspace) << mdl_data->robot_lab_proxy << mdl_data->robot_airship_proxy; \
         \
        typedef pp::topology_traits< robot_airship::CRS3D_jspace_rl_o2_type >::point_type RLPointType; \
        typedef pp::topology_traits< robot_airship::CRS3D_jspace_o2_type >::point_type PointType; \
        RLPointType start_point, goal_point; \
        PointType start_inter, goal_inter; \
        start_inter = normal_jt_space->origin(); \
        get<0>(start_inter) = vect<double,7>(jt_start); \
        start_point = mdl_data->manip_jt_limits->map_to_space(start_inter, *normal_jt_space, *jt_space); \
         \
        goal_inter = normal_jt_space->origin(); \
        get<0>(goal_inter) = vect<double,7>(jt_desired); \
        goal_point = mdl_data->manip_jt_limits->map_to_space(goal_inter, *normal_jt_space, *jt_space); \
         \
        typedef pp::frame_tracer_3D< robot_airship::CRS3D_rlDK_o2_type, robot_airship::CRS3D_jspace_rl_o2_type, pp::identity_topo_map, pp::print_sbmp_progress<> > frame_reporter_type; \
        frame_reporter_type temp_reporter( \
          robot_airship::CRS3D_rlDK_o2_type(mdl_data->manip_kin_mdl, mdl_data->manip_jt_limits, normal_jt_space), \
          jt_space, 0.5 * plan_options->min_travel); \
        temp_reporter.add_traced_frame(mdl_data->builder.arm_joint_6_end); \
         \
        pp::any_sbmp_reporter_chain< WORKSPACE > report_chain; \
        report_chain.add_reporter( boost::ref(temp_reporter) ); \
         \
        pp::path_planning_p2p_query< WORKSPACE > pp_query("pp_query", workspace, \
          start_point, goal_point, plan_options->max_results);
        
        
#if 0
        
#define RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(WORKSPACE) \
        pp::rrt_planner< WORKSPACE > workspace_planner( \
            workspace, plan_options->max_vertices, plan_options->prog_interval, \
            plan_options->store_policy | plan_options->knn_method, \
            plan_options->planning_options, 0.1, 0.05, report_chain); \
         \
        pp_query.reset_solution_records(); \
        workspace_planner.solve_planning_query(pp_query); \
         \
        shared_ptr< pp::seq_path_base< SuperSpaceType > > bestsol_rlpath; \
        if(pp_query.solutions.size()) \
          bestsol_rlpath = pp_query.solutions.begin()->second; \
        std::cout << "The shortest distance is: " << pp_query.get_best_solution_distance() << std::endl; \
         \
        sol_anim->bestsol_trajectory.clear(); \
        if(bestsol_rlpath) { \
          typedef pp::seq_path_base< SuperSpaceType >::point_fraction_iterator PtIter; \
          for(PtIter it = bestsol_rlpath->begin_fraction_travel(); it != bestsol_rlpath->end_fraction_travel(); it += 0.1) \
            sol_anim->bestsol_trajectory.push_back( get<0>(mdl_data->manip_jt_limits->map_to_space(*it, *jt_space, *normal_jt_space)) ); \
        }; \
         \
        mg_sep = temp_reporter.get_motion_graph_tracer(mdl_data->builder.arm_joint_6_end).get_separator(); \
        mg_sep->ref(); \
        for(std::size_t i = 0; i < temp_reporter.get_solution_count(); ++i) { \
          sol_seps.push_back(temp_reporter.get_solution_tracer(mdl_data->builder.arm_joint_6_end, i).get_separator()); \
          sol_seps.back()->ref(); \
        }; \
         \
        mdl_data->builder.track_joint_coord->q = get<0>(start_inter)[0]; \
        mdl_data->builder.arm_joint_1_coord->q = get<0>(start_inter)[1]; \
        mdl_data->builder.arm_joint_2_coord->q = get<0>(start_inter)[2]; \
        mdl_data->builder.arm_joint_3_coord->q = get<0>(start_inter)[3]; \
        mdl_data->builder.arm_joint_4_coord->q = get<0>(start_inter)[4]; \
        mdl_data->builder.arm_joint_5_coord->q = get<0>(start_inter)[5]; \
        mdl_data->builder.arm_joint_6_coord->q = get<0>(start_inter)[6]; \
        mdl_data->kin_chain->doMotion();
  
#define RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(WORKSPACE) \
        pp::prm_planner< WORKSPACE > workspace_planner( \
          workspace, plan_options->max_vertices, plan_options->prog_interval, \
          plan_options->store_policy | plan_options->knn_method, \
          0.1, 0.05, plan_options->max_travel, workspace_dims, report_chain); \
         \
        pp_query.reset_solution_records(); \
        workspace_planner.solve_planning_query(pp_query); \
         \
        shared_ptr< pp::seq_path_base< SuperSpaceType > > bestsol_rlpath; \
        if(pp_query.solutions.size()) \
          bestsol_rlpath = pp_query.solutions.begin()->second; \
        std::cout << "The shortest distance is: " << pp_query.get_best_solution_distance() << std::endl; \
         \
        sol_anim->bestsol_trajectory.clear(); \
        if(bestsol_rlpath) { \
          typedef pp::seq_path_base< SuperSpaceType >::point_fraction_iterator PtIter; \
          for(PtIter it = bestsol_rlpath->begin_fraction_travel(); it != bestsol_rlpath->end_fraction_travel(); it += 0.1) \
            sol_anim->bestsol_trajectory.push_back( get<0>(mdl_data->manip_jt_limits->map_to_space(*it, *jt_space, *normal_jt_space)) ); \
        }; \
         \
        mg_sep = temp_reporter.get_motion_graph_tracer(mdl_data->builder.arm_joint_6_end).get_separator(); \
        mg_sep->ref(); \
        for(std::size_t i = 0; i < temp_reporter.get_solution_count(); ++i) { \
          sol_seps.push_back(temp_reporter.get_solution_tracer(mdl_data->builder.arm_joint_6_end, i).get_separator()); \
          sol_seps.back()->ref(); \
        }; \
         \
        mdl_data->builder.track_joint_coord->q = get<0>(start_inter)[0]; \
        mdl_data->builder.arm_joint_1_coord->q = get<0>(start_inter)[1]; \
        mdl_data->builder.arm_joint_2_coord->q = get<0>(start_inter)[2]; \
        mdl_data->builder.arm_joint_3_coord->q = get<0>(start_inter)[3]; \
        mdl_data->builder.arm_joint_4_coord->q = get<0>(start_inter)[4]; \
        mdl_data->builder.arm_joint_5_coord->q = get<0>(start_inter)[5]; \
        mdl_data->builder.arm_joint_6_coord->q = get<0>(start_inter)[6]; \
        mdl_data->kin_chain->doMotion();
  
#define RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(WORKSPACE) \
        pp::fadprm_planner< WORKSPACE > workspace_planner( \
          workspace, plan_options->max_vertices, plan_options->prog_interval, \
          plan_options->store_policy | plan_options->knn_method, \
          0.1, 0.05, plan_options->max_travel, workspace_dims, report_chain); \
         \
        workspace_planner.set_initial_relaxation(plan_options->init_relax); \
         \
        pp_query.reset_solution_records(); \
        workspace_planner.solve_planning_query(pp_query); \
         \
        shared_ptr< pp::seq_path_base< SuperSpaceType > > bestsol_rlpath; \
        if(pp_query.solutions.size()) \
          bestsol_rlpath = pp_query.solutions.begin()->second; \
        std::cout << "The shortest distance is: " << pp_query.get_best_solution_distance() << std::endl; \
         \
        sol_anim->bestsol_trajectory.clear(); \
        if(bestsol_rlpath) { \
          typedef pp::seq_path_base< SuperSpaceType >::point_fraction_iterator PtIter; \
          for(PtIter it = bestsol_rlpath->begin_fraction_travel(); it != bestsol_rlpath->end_fraction_travel(); it += 0.1) \
            sol_anim->bestsol_trajectory.push_back( get<0>(mdl_data->manip_jt_limits->map_to_space(*it, *jt_space, *normal_jt_space)) ); \
        }; \
         \
        mg_sep = temp_reporter.get_motion_graph_tracer(mdl_data->builder.arm_joint_6_end).get_separator(); \
        mg_sep->ref(); \
        for(std::size_t i = 0; i < temp_reporter.get_solution_count(); ++i) { \
          sol_seps.push_back(temp_reporter.get_solution_tracer(mdl_data->builder.arm_joint_6_end, i).get_separator()); \
          sol_seps.back()->ref(); \
        }; \
         \
        mdl_data->builder.track_joint_coord->q = get<0>(start_inter)[0]; \
        mdl_data->builder.arm_joint_1_coord->q = get<0>(start_inter)[1]; \
        mdl_data->builder.arm_joint_2_coord->q = get<0>(start_inter)[2]; \
        mdl_data->builder.arm_joint_3_coord->q = get<0>(start_inter)[3]; \
        mdl_data->builder.arm_joint_4_coord->q = get<0>(start_inter)[4]; \
        mdl_data->builder.arm_joint_5_coord->q = get<0>(start_inter)[5]; \
        mdl_data->builder.arm_joint_6_coord->q = get<0>(start_inter)[6]; \
        mdl_data->kin_chain->doMotion();
  
#endif
  
#define RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(WORKSPACE) \
        pp::rrtstar_planner< WORKSPACE > workspace_planner( \
          workspace, plan_options->max_vertices, plan_options->prog_interval, \
          plan_options->store_policy | plan_options->knn_method, \
          plan_options->planning_options, \
          0.1, 0.05, workspace_dims, report_chain); \
         \
        pp_query.reset_solution_records(); \
        workspace_planner.solve_planning_query(pp_query); \
         \
        shared_ptr< pp::seq_path_base< SuperSpaceType > > bestsol_rlpath; \
        if(pp_query.solutions.size()) \
          bestsol_rlpath = pp_query.solutions.begin()->second; \
        std::cout << "The shortest distance is: " << pp_query.get_best_solution_distance() << std::endl; \
         \
        sol_anim->bestsol_trajectory.clear(); \
        if(bestsol_rlpath) { \
          typedef pp::seq_path_base< SuperSpaceType >::point_fraction_iterator PtIter; \
          for(PtIter it = bestsol_rlpath->begin_fraction_travel(); it != bestsol_rlpath->end_fraction_travel(); it += 0.1) \
            sol_anim->bestsol_trajectory.push_back( get<0>(mdl_data->manip_jt_limits->map_to_space(*it, *jt_space, *normal_jt_space)) ); \
        }; \
         \
        mg_sep = temp_reporter.get_motion_graph_tracer(mdl_data->builder.arm_joint_6_end).get_separator(); \
        mg_sep->ref(); \
        for(std::size_t i = 0; i < temp_reporter.get_solution_count(); ++i) { \
          sol_seps.push_back(temp_reporter.get_solution_tracer(mdl_data->builder.arm_joint_6_end, i).get_separator()); \
          sol_seps.back()->ref(); \
        }; \
         \
        mdl_data->builder.track_joint_coord->q = get<0>(start_inter)[0]; \
        mdl_data->builder.arm_joint_1_coord->q = get<0>(start_inter)[1]; \
        mdl_data->builder.arm_joint_2_coord->q = get<0>(start_inter)[2]; \
        mdl_data->builder.arm_joint_3_coord->q = get<0>(start_inter)[3]; \
        mdl_data->builder.arm_joint_4_coord->q = get<0>(start_inter)[4]; \
        mdl_data->builder.arm_joint_5_coord->q = get<0>(start_inter)[5]; \
        mdl_data->builder.arm_joint_6_coord->q = get<0>(start_inter)[6]; \
        mdl_data->kin_chain->doMotion();
  
#define RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(WORKSPACE) \
        pp::sbastar_planner< WORKSPACE > workspace_planner( \
          workspace, plan_options->max_vertices, plan_options->prog_interval, \
          plan_options->store_policy | plan_options->knn_method, \
          plan_options->planning_options, \
          0.1, 0.05, plan_options->max_travel, workspace_dims, report_chain); \
         \
        workspace_planner.set_initial_density_threshold(0.0); \
        workspace_planner.set_initial_relaxation(plan_options->init_relax); \
        workspace_planner.set_initial_SA_temperature(plan_options->init_SA_temp); \
         \
        pp_query.reset_solution_records(); \
        workspace_planner.solve_planning_query(pp_query); \
         \
        shared_ptr< pp::seq_path_base< SuperSpaceType > > bestsol_rlpath; \
        if(pp_query.solutions.size()) \
          bestsol_rlpath = pp_query.solutions.begin()->second; \
        std::cout << "The shortest distance is: " << pp_query.get_best_solution_distance() << std::endl; \
         \
        sol_anim->bestsol_trajectory.clear(); \
        if(bestsol_rlpath) { \
          typedef pp::seq_path_base< SuperSpaceType >::point_fraction_iterator PtIter; \
          for(PtIter it = bestsol_rlpath->begin_fraction_travel(); it != bestsol_rlpath->end_fraction_travel(); it += 0.1) \
            sol_anim->bestsol_trajectory.push_back( get<0>(mdl_data->manip_jt_limits->map_to_space(*it, *jt_space, *normal_jt_space)) ); \
        }; \
         \
        mg_sep = temp_reporter.get_motion_graph_tracer(mdl_data->builder.arm_joint_6_end).get_separator(); \
        mg_sep->ref(); \
        for(std::size_t i = 0; i < temp_reporter.get_solution_count(); ++i) { \
          sol_seps.push_back(temp_reporter.get_solution_tracer(mdl_data->builder.arm_joint_6_end, i).get_separator()); \
          sol_seps.back()->ref(); \
        }; \
         \
        mdl_data->builder.track_joint_coord->q = get<0>(start_inter)[0]; \
        mdl_data->builder.arm_joint_1_coord->q = get<0>(start_inter)[1]; \
        mdl_data->builder.arm_joint_2_coord->q = get<0>(start_inter)[2]; \
        mdl_data->builder.arm_joint_3_coord->q = get<0>(start_inter)[3]; \
        mdl_data->builder.arm_joint_4_coord->q = get<0>(start_inter)[4]; \
        mdl_data->builder.arm_joint_5_coord->q = get<0>(start_inter)[5]; \
        mdl_data->builder.arm_joint_6_coord->q = get<0>(start_inter)[6]; \
        mdl_data->kin_chain->doMotion();
  
  
  
  
  switch(plan_options->planning_algo) {
#if 0
    case 0:  // RRT
    {
      
      if((plan_options->space_order == 0) && (plan_options->interp_id == 0)) { 
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_0(robot_airship::CRS3D_workspace_o0_i1_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(robot_airship::CRS3D_workspace_o0_i1_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_i1_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_i1_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i1_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i1_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_i3_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_i3_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i3_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i3_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 2)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i5_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i5_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_svp_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_svp_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_svp_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_svp_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 4)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_sap_type)
        RK_CRS_PLANNER_GENERATE_RRT_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_sap_type)
      };
      
    }; break;
#endif
    case 1:  // RRT*
    {
      
      if((plan_options->space_order == 0) && (plan_options->interp_id == 0)) { 
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_0(robot_airship::CRS3D_workspace_o0_i1_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o0_i1_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_i1_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_i1_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i1_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i1_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_i3_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_i3_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i3_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i3_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 2)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i5_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i5_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_svp_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_svp_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_svp_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_svp_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 4)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_sap_type)
        RK_CRS_PLANNER_GENERATE_RRTSTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_sap_type)
      };
      
    }; break;
#if 0
    case 2:  // PRM
    {
      
      if((plan_options->space_order == 0) && (plan_options->interp_id == 0)) { 
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_0(robot_airship::CRS3D_workspace_o0_i1_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o0_i1_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_i1_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_i1_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i1_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i1_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_i3_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_i3_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i3_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i3_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 2)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i5_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i5_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_svp_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_svp_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_svp_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_svp_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 4)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_sap_type)
        RK_CRS_PLANNER_GENERATE_PRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_sap_type)
      };
      
    }; break;
#endif
    case 3:  // SBA*
    { 
      
      if((plan_options->space_order == 0) && (plan_options->interp_id == 0)) { 
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_0(robot_airship::CRS3D_workspace_o0_i1_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o0_i1_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_i1_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_i1_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i1_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i1_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_i3_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_i3_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i3_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i3_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 2)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i5_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i5_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_svp_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_svp_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_svp_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_svp_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 4)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_sap_type)
        RK_CRS_PLANNER_GENERATE_SBASTAR_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_sap_type)
      };
      
    }; break;
#if 0
    case 4:  // FADPRM
    { 
      
      if((plan_options->space_order == 0) && (plan_options->interp_id == 0)) { 
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_0(robot_airship::CRS3D_workspace_o0_i1_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o0_i1_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_i1_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_i1_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 0)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i1_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i1_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_i3_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_i3_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 1)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i3_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i3_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 2)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_i5_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_i5_type)
      } else if((plan_options->space_order == 1) && (plan_options->interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_1(robot_airship::CRS3D_workspace_o1_svp_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o1_svp_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 3)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_svp_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_svp_type)
      } else if((plan_options->space_order == 2) && (plan_options->interp_id == 4)) {
        RK_CRS_PLANNER_GENERATE_PLANNER_IC_2(robot_airship::CRS3D_workspace_o2_sap_type)
        RK_CRS_PLANNER_GENERATE_FADPRM_PLANNER_CALL(robot_airship::CRS3D_workspace_o2_sap_type)
      };
      
    }; break;
#endif
    case 5:  // LSBA*
    { 
      
    };// break;
    case 6:  // ???
    { 
      
    };// break;
    default:
      QMessageBox::information(this,
                  "Planner Not Supported!",
                  "Sorry, the planning method you selected is not yet supported!",
                  QMessageBox::Ok);
  };
  
  // Check the motion-graph separator and solution separators
  //  add them to the switches.
  if(mg_sep) {
    draw_data->sw_motion_graph->removeAllChildren();
    draw_data->sw_motion_graph->addChild(mg_sep);
    mg_sep->unref();
  };
  
  draw_data->sw_solutions->removeAllChildren();
  if((configs.check_print_best->isChecked()) && (sol_seps.size())) {
    draw_data->sw_solutions->addChild(sol_seps[0]);
    sol_seps[0]->unref();
  };
  
};





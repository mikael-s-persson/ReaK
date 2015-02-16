/**
 * \file planner_exec_intercept.hpp
 * 
 * This library defines functions and classes useful to execute path-planners for interception problems.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2015
 */

/*
 *    Copyright 2015 Sven Mikael Persson
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

#ifndef REAK_PLANNER_EXEC_INTERCEPT_HPP
#define REAK_PLANNER_EXEC_INTERCEPT_HPP

#include "planner_exec_engines.hpp"

#include "intercept_query.hpp"

#include <ReaK/topologies/interpolation/trajectory_base.hpp>

namespace ReaK {

namespace pp {


template <typename Topology, typename PlanEngine, typename GoalTrajType>
void execute_intercept_planner(const shared_ptr< Topology >& world_topo,
                               const planning_option_collection& plan_options,
                               std::size_t world_dimensionality,
                               double time_horizon, double time_step,
                               PlanEngine& engine,
                               const typename topology_traits<Topology>::point_type& p_start,
                               const shared_ptr<GoalTrajType>& goal_traj) {
  
  // NOTE: From the original CRS dynexec functions, this should have:
  //        (Topology == dynamic_workspace_type)
  //        (GoalTrajType == mapped_traj_type)
  //        (p_start == temporal_start_point)
  //        (world_topo == workspace)
  //        (time_horizon == target_state_traj->get_end_time())
  //        (time_step == min_travel)
  
  // Create the reporter chain.
  shared_ptr< any_sbmp_reporter_chain< Topology > > p_report_chain = engine.create_reporter(world_topo);
  
  // Create the point-to-point query:
  motion_plan_intercept_query< Topology, GoalTrajType > pp_query(
    "intercept_query", world_topo, p_start, goal_traj, time_horizon, time_step, plan_options.max_results);
  
  // Create the planner:
  shared_ptr< sample_based_planner< Topology > > world_planner;
  
#ifndef RK_DISABLE_RRT_PLANNER
  if( plan_options.planning_algo == 0 ) { // RRT
    
    world_planner = shared_ptr< sample_based_planner< Topology > >(
      new rrt_planner< Topology >(
        world_topo, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, *p_report_chain));
    
  } else 
#endif
#ifndef RK_DISABLE_RRTSTAR_PLANNER
  if( plan_options.planning_algo == 1 ) { // RRT*
    
    world_planner = shared_ptr< sample_based_planner< Topology > >(
      new rrtstar_planner< Topology >(
        world_topo, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, world_dimensionality, *p_report_chain));
    
  } else 
#endif
#ifndef RK_DISABLE_PRM_PLANNER
  if( plan_options.planning_algo == 2 ) { // PRM
    
    world_planner = shared_ptr< sample_based_planner< Topology > >(
      new prm_planner< Topology >(
        world_topo, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        0.1, 0.05, plan_options.max_random_walk, world_dimensionality, *p_report_chain));
    
  } else 
#endif
#ifndef RK_DISABLE_FADPRM_PLANNER
  if( plan_options.planning_algo == 4 ) { // FADPRM
    
    shared_ptr< fadprm_planner< Topology > > tmp(
      new fadprm_planner< Topology >(
        world_topo, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        0.1, 0.05, plan_options.max_random_walk, world_dimensionality, *p_report_chain));
    
    tmp->set_initial_relaxation(plan_options.init_relax);
    
    world_planner = tmp;
    
  } else 
#endif
#ifndef RK_DISABLE_SBASTAR_PLANNER
  if( plan_options.planning_algo == 3 ) { // SBA*
    
    shared_ptr< sbastar_planner< Topology > > tmp(
      new sbastar_planner< Topology >(
        world_topo, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options,
        0.1, 0.05, plan_options.max_random_walk, world_dimensionality, *p_report_chain));
    
    tmp->set_initial_density_threshold(0.0);
    tmp->set_initial_relaxation(plan_options.init_relax);
    tmp->set_initial_SA_temperature(plan_options.init_SA_temp);
    
    world_planner = tmp;
    
  } else 
#endif
  { };
  
  if(!world_planner)
    return;
  
  // Solve the planning problem:
  engine(plan_options, world_planner, pp_query);
  
};





};


};


#endif



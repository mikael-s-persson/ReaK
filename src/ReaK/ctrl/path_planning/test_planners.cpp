
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

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/properties.hpp>

#include "topologies/ptrobot2D_test_world.hpp"
#include "rrt_path_planner.hpp"

#include "basic_sbmp_reporters.hpp"


int main() {
  
  ReaK::shared_ptr< ReaK::pp::ptrobot2D_test_world > world_map =
    ReaK::shared_ptr< ReaK::pp::ptrobot2D_test_world >(new ReaK::pp::ptrobot2D_test_world("test_world.bmp", 10, 1.0));
  
  ReaK::pp::rrt_path_planner< ReaK::pp::ptrobot2D_test_world, ReaK::pp::differ_sbmp_report_to_space<> > 
    rrt_plan(world_map, 
	     world_map->get_start_pos(), 
	     world_map->get_goal_pos(),
	     10000, 
	     100,
	     ReaK::pp::BIDIRECTIONAL_RRT,
	     ReaK::pp::ADJ_LIST_MOTION_GRAPH,
	     ReaK::pp::DVP_BF2_TREE_KNN,
	     ReaK::pp::differ_sbmp_report_to_space<>("pp_results/rrt/test_world_", 5),
	     50);
  
  
  
  
};














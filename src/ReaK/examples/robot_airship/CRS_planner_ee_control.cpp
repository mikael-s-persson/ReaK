
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
  
  detail::write_joint_coordinates_impl(target_tp.pt, 
                                       target_state_traj->get_temporal_space().get_space_topology(), 
                                       scene_data.target_kin_model);
  scene_data.target_kin_model->doDirectMotion();
  
  shared_ptr< frame_3D<double> > capture_frame = scene_data.target_kin_model->getDependentFrame3D(0)->mFrame;
  
  
  vect_n<double> q_current = chaser_concrete_model->getJointPositions();
  vect_n<double> qd_current = chaser_concrete_model->getJointVelocities();
  chaser_concrete_model->doDirectMotion();
  
  shared_ptr< frame_3D<double> > EE_current = chaser_concrete_model->getDependentFrame3D(0)->mFrame;
  
  mat<double, mat_structure::rectangular> Jac_current;
  chaser_concrete_model->getJacobianMatrix(Jac_current);
  
  
  
  
  
};












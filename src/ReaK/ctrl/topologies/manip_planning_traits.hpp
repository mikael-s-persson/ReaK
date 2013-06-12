/**
 * \file manip_planning_traits.hpp
 * 
 * This library defines a class
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_MANIP_PLANNING_TRAITS_HPP
#define REAK_MANIP_PLANNING_TRAITS_HPP

#include "base/defs.hpp"

#include "manip_free_workspace.hpp"
#include "manip_free_svp_workspace.hpp"
#include "manip_free_sap_workspace.hpp"
#include "manip_free_dynamic_workspace.hpp"
#include "manip_free_svp_dynamic_workspace.hpp"
#include "manip_free_sap_dynamic_workspace.hpp"

#include "direct_kinematics_topomap.hpp"
#include "inverse_kinematics_topomap.hpp"

namespace ReaK {

namespace pp {


template <typename ManipMdlType>
struct manip_pp_traits {
  BOOST_STATIC_CONSTANT(std::size_t, degrees_of_freedom = ManipMdlType::degrees_of_freedom);
  
//   Ndof_rl_space<double, ManipMdlType::degrees_of_freedom, 0>::type
  typedef ManipMdlType::rl_o0_jt_space_type rl_o0_jt_space_type;
//   Ndof_space<double, ManipMdlType::degrees_of_freedom, 0>::type
  typedef ManipMdlType::o0_jt_space_type o0_jt_space_type;
//   se3_0th_order_topology<double>::type
  typedef ManipMdlType::o0_ee_space_type o0_ee_space_type;
  
//   Ndof_rl_space<double, ManipMdlType::degrees_of_freedom, 1>::type
  typedef ManipMdlType::rl_o1_jt_space_type rl_o1_jt_space_type;
//   Ndof_space<double, ManipMdlType::degrees_of_freedom, 1>::type
  typedef ManipMdlType::o1_jt_space_type o1_jt_space_type;
//   se3_1st_order_topology<double>::type
  typedef ManipMdlType::o1_ee_space_type o1_ee_space_type;
  
//   Ndof_rl_space<double, ManipMdlType::degrees_of_freedom, 2>::type
  typedef ManipMdlType::rl_o2_jt_space_type rl_o2_jt_space_type;
//   Ndof_space<double, ManipMdlType::degrees_of_freedom, 2>::type
  typedef ManipMdlType::o2_jt_space_type o2_jt_space_type;
//   se3_2nd_order_topology<double>::type
  typedef ManipMdlType::o2_ee_space_type o2_ee_space_type;
};


template <typename ManipMdlType, typename InterpTag = linear_interpolation_tag>
struct manip_static_workspace {
  typedef manip_pp_traits<ManipMdlType> Traits;
  
  typedef manip_quasi_static_env<typename Traits::rl_o0_jt_space_type, InterpTag>  rl_o0_workspace_type;
  typedef manip_quasi_static_env<typename Traits::rl_o1_jt_space_type, InterpTag>  rl_o1_workspace_type;
  typedef manip_quasi_static_env<typename Traits::rl_o2_jt_space_type, InterpTag>  rl_o2_workspace_type;
  
  typedef manip_quasi_static_env<typename Traits::o0_jt_space_type, InterpTag>  o0_workspace_type;
  typedef manip_quasi_static_env<typename Traits::o1_jt_space_type, InterpTag>  o1_workspace_type;
  typedef manip_quasi_static_env<typename Traits::o2_jt_space_type, InterpTag>  o2_workspace_type;
};


template <typename ManipMdlType, typename InterpTag = linear_interpolation_tag>
struct manip_dynamic_workspace {
  typedef manip_pp_traits<ManipMdlType> Traits;
  
  typedef manip_dynamic_env<typename Traits::rl_o0_jt_space_type, InterpTag>  rl_o0_workspace_type;
  typedef manip_dynamic_env<typename Traits::rl_o1_jt_space_type, InterpTag>  rl_o1_workspace_type;
  typedef manip_dynamic_env<typename Traits::rl_o2_jt_space_type, InterpTag>  rl_o2_workspace_type;
};


template <typename ManipMdlType>
struct manip_DK_map {
  typedef manip_pp_traits<ManipMdlType> Traits;
  
  typedef manip_rl_direct_kin_map< joint_limits_collection<double>, typename Traits::o0_jt_space_type > rl_o0_map_type;
  typedef manip_rl_direct_kin_map< joint_limits_collection<double>, typename Traits::o1_jt_space_type > rl_o1_map_type;
  typedef manip_rl_direct_kin_map< joint_limits_collection<double>, typename Traits::o2_jt_space_type > rl_o2_map_type;
  
  typedef manip_direct_kin_map o0_map_type;
  typedef manip_direct_kin_map o1_map_type;
  typedef manip_direct_kin_map o2_map_type;
};


template <typename ManipMdlType>
struct manip_IK_map {
  typedef manip_rl_inverse_kin_map< joint_limits_collection<double> > rl_o0_map_type;
  typedef manip_rl_inverse_kin_map< joint_limits_collection<double> > rl_o1_map_type;
  typedef manip_rl_inverse_kin_map< joint_limits_collection<double> > rl_o2_map_type;
  
  typedef manip_inverse_kin_map o0_map_type;
  typedef manip_inverse_kin_map o1_map_type;
  typedef manip_inverse_kin_map o2_map_type;
};




template <typename InterpTag, typename ManipMdlType>
shared_ptr< typename manip_static_workspace< ManipMdlType, InterpTag >::rl_o0_workspace_type >
  make_manip_o0_static_workspace(
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits,
    double min_travel, double max_travel) {
  typedef typename manip_static_workspace< ManipMdlType, InterpTag >::rl_o0_workspace_type WorkspaceType;
  
  return shared_ptr<WorkspaceType>(new WorkspaceType(
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits),
    manip_kin_mdl,
    manip_jt_limits,
    min_travel,
    max_travel));
};

template <typename InterpTag, typename ManipMdlType>
shared_ptr< typename manip_static_workspace< ManipMdlType, InterpTag >::rl_o1_workspace_type >
  make_manip_o1_static_workspace(
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits,
    double min_travel, double max_travel) {
  typedef typename manip_static_workspace< ManipMdlType, InterpTag >::rl_o1_workspace_type WorkspaceType;
  
  return shared_ptr<WorkspaceType>(new WorkspaceType(
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits, 
      manip_jt_limits->gen_accel_limits),
    manip_kin_mdl,
    manip_jt_limits,
    min_travel,
    max_travel));
};

template <typename InterpTag, typename ManipMdlType>
shared_ptr< typename manip_static_workspace< ManipMdlType, InterpTag >::rl_o2_workspace_type >
  make_manip_o2_static_workspace(
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits,
    double min_travel, double max_travel) {
  typedef typename manip_static_workspace< ManipMdlType, InterpTag >::rl_o2_workspace_type WorkspaceType;
  
  return shared_ptr<WorkspaceType>(new WorkspaceType(
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits, 
      manip_jt_limits->gen_accel_limits, 
      manip_jt_limits->gen_jerk_limits),
    manip_kin_mdl,
    manip_jt_limits,
    min_travel,
    max_travel));
};


template <typename InterpTag, typename ManipMdlType>
shared_ptr< typename manip_dynamic_workspace< ManipMdlType, InterpTag >::rl_o0_workspace_type >
  make_manip_o0_dynamic_workspace(
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits,
    double min_travel, double max_travel) {
  typedef typename manip_dynamic_workspace< ManipMdlType, InterpTag >::rl_o0_workspace_type WorkspaceType;
  
  return shared_ptr<WorkspaceType>(new WorkspaceType(
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits),
    manip_kin_mdl,
    manip_jt_limits,
    min_travel,
    max_travel));
};

template <typename InterpTag, typename ManipMdlType>
shared_ptr< typename manip_dynamic_workspace< ManipMdlType, InterpTag >::rl_o1_workspace_type >
  make_manip_o1_dynamic_workspace(
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits,
    double min_travel, double max_travel) {
  typedef typename manip_dynamic_workspace< ManipMdlType, InterpTag >::rl_o1_workspace_type WorkspaceType;
  
  return shared_ptr<WorkspaceType>(new WorkspaceType(
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits, 
      manip_jt_limits->gen_accel_limits),
    manip_kin_mdl,
    manip_jt_limits,
    min_travel,
    max_travel));
};

template <typename InterpTag, typename ManipMdlType>
shared_ptr< typename manip_dynamic_workspace< ManipMdlType, InterpTag >::rl_o2_workspace_type >
  make_manip_o2_dynamic_workspace(
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits,
    double min_travel, double max_travel) {
  typedef typename manip_dynamic_workspace< ManipMdlType, InterpTag >::rl_o2_workspace_type WorkspaceType;
  
  return shared_ptr<WorkspaceType>(new WorkspaceType(
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits, 
      manip_jt_limits->gen_accel_limits, 
      manip_jt_limits->gen_jerk_limits),
    manip_kin_mdl,
    manip_jt_limits,
    min_travel,
    max_travel));
};





};

};

#endif


/**
 * \file CRS_workspaces.hpp
 * 
 * This library defines a class
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_CRS_WORKSPACES_HPP
#define REAK_CRS_WORKSPACES_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "CRS_A465_models.hpp"

#include "interpolation/linear_interp.hpp"
#include "interpolation/cubic_hermite_interp.hpp"
#include "interpolation/quintic_hermite_interp.hpp"
#include "interpolation/sustained_velocity_pulse.hpp"
#include "interpolation/sustained_acceleration_pulse.hpp"
#include "topologies/manip_free_workspace.hpp"
#include "topologies/manip_free_svp_workspace.hpp"
#include "topologies/manip_free_sap_workspace.hpp"

namespace ReaK {
  
namespace robot_airship {
  
typedef CRS_A465_model_builder::rate_limited_joint_space_0th_type CRS3D_jspace_rl_o0_type;
typedef CRS_A465_model_builder::rate_limited_joint_space_1st_type CRS3D_jspace_rl_o1_type;
typedef CRS_A465_model_builder::rate_limited_joint_space_type     CRS3D_jspace_rl_o2_type;

typedef CRS_A465_model_builder::joint_space_0th_type CRS3D_jspace_o0_type;
typedef CRS_A465_model_builder::joint_space_1st_type CRS3D_jspace_o1_type;
typedef CRS_A465_model_builder::joint_space_type     CRS3D_jspace_o2_type;

typedef pp::manip_quasi_static_env<CRS3D_jspace_rl_o0_type, pp::linear_interpolation_tag> CRS3D_workspace_o0_i1_type;
typedef pp::manip_quasi_static_env<CRS3D_jspace_rl_o1_type, pp::linear_interpolation_tag> CRS3D_workspace_o1_i1_type;
typedef pp::manip_quasi_static_env<CRS3D_jspace_rl_o2_type, pp::linear_interpolation_tag> CRS3D_workspace_o2_i1_type;

typedef pp::manip_quasi_static_env<CRS3D_jspace_rl_o1_type, pp::cubic_hermite_interpolation_tag> CRS3D_workspace_o1_i3_type;
typedef pp::manip_quasi_static_env<CRS3D_jspace_rl_o2_type, pp::cubic_hermite_interpolation_tag> CRS3D_workspace_o2_i3_type;

typedef pp::manip_quasi_static_env<CRS3D_jspace_rl_o2_type, pp::quintic_hermite_interpolation_tag> CRS3D_workspace_o2_i5_type;

typedef pp::manip_quasi_static_env<CRS3D_jspace_rl_o1_type, pp::svp_Ndof_interpolation_tag> CRS3D_workspace_o1_svp_type;
typedef pp::manip_quasi_static_env<CRS3D_jspace_rl_o2_type, pp::svp_Ndof_interpolation_tag> CRS3D_workspace_o2_svp_type;

typedef pp::manip_quasi_static_env<CRS3D_jspace_rl_o2_type, pp::sap_Ndof_interpolation_tag> CRS3D_workspace_o2_sap_type;


typedef pp::manip_rl_direct_kin_map< pp::joint_limits_collection<double>, CRS3D_jspace_o0_type > CRS3D_rlDK_o0_type;
typedef pp::manip_rl_direct_kin_map< pp::joint_limits_collection<double>, CRS3D_jspace_o1_type > CRS3D_rlDK_o1_type;
typedef pp::manip_rl_direct_kin_map< pp::joint_limits_collection<double>, CRS3D_jspace_o2_type > CRS3D_rlDK_o2_type;


// typedef pp::frame_tracer_3D< CRS3D_rlDK_o0_type, CRS3D_jspace_rl_o0_type, pp::identity_topo_map, pp::print_sbmp_progress<> > CRS3D_rl_o0_tracer;
// typedef pp::frame_tracer_3D< CRS3D_rlDK_o1_type, CRS3D_jspace_rl_o1_type, pp::identity_topo_map, pp::print_sbmp_progress<> > CRS3D_rl_o1_tracer;
// typedef pp::frame_tracer_3D< CRS3D_rlDK_o2_type, CRS3D_jspace_rl_o2_type, pp::identity_topo_map, pp::print_sbmp_progress<> > CRS3D_rl_o2_tracer;


};

};

#endif


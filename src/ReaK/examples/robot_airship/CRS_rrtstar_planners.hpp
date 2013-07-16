/**
 * \file CRS_rrtstar_planners.hpp
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

#ifndef REAK_CRS_RRTSTAR_PLANNERS_HPP
#define REAK_CRS_RRTSTAR_PLANNERS_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "CRS_workspaces.hpp"

#include "path_planning/rrtstar_path_planner.hpp"


namespace ReaK {
  
namespace robot_airship {
  
//typedef pp::rrtstar_planner<WORKSPACE> CRS3D_rrtstar_o0_i1_type;

typedef pp::rrtstar_planner<CRS3D_workspace_o0_i1_type> CRS3D_rrtstar_o0_i1_type;
typedef pp::rrtstar_planner<CRS3D_workspace_o1_i1_type> CRS3D_rrtstar_o1_i1_type;
typedef pp::rrtstar_planner<CRS3D_workspace_o2_i1_type> CRS3D_rrtstar_o2_i1_type;

typedef pp::rrtstar_planner<CRS3D_workspace_o1_i3_type> CRS3D_rrtstar_o1_i3_type;
typedef pp::rrtstar_planner<CRS3D_workspace_o2_i3_type> CRS3D_rrtstar_o2_i3_type;

typedef pp::rrtstar_planner<CRS3D_workspace_o2_i5_type> CRS3D_rrtstar_o2_i5_type;

typedef pp::rrtstar_planner<CRS3D_workspace_o1_svp_type> CRS3D_rrtstar_o1_svp_type;
typedef pp::rrtstar_planner<CRS3D_workspace_o2_svp_type> CRS3D_rrtstar_o2_svp_type;

typedef pp::rrtstar_planner<CRS3D_workspace_o2_sap_type> CRS3D_rrtstar_o2_sap_type;


};


namespace pp {


#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

//extern template class rrtstar_planner<WORKSPACE>;

extern template class rrtstar_planner<robot_airship::CRS3D_workspace_o0_i1_type>;
extern template class rrtstar_planner<robot_airship::CRS3D_workspace_o1_i1_type>;
extern template class rrtstar_planner<robot_airship::CRS3D_workspace_o2_i1_type>;

extern template class rrtstar_planner<robot_airship::CRS3D_workspace_o1_i3_type>;
extern template class rrtstar_planner<robot_airship::CRS3D_workspace_o2_i3_type>;

extern template class rrtstar_planner<robot_airship::CRS3D_workspace_o2_i5_type>;

extern template class rrtstar_planner<robot_airship::CRS3D_workspace_o1_svp_type>;
extern template class rrtstar_planner<robot_airship::CRS3D_workspace_o2_svp_type>;

extern template class rrtstar_planner<robot_airship::CRS3D_workspace_o2_sap_type>;

#endif


};

};

#endif


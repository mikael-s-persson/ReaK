
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


#include "base/defs.hpp"

#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

#include "CRS_rrt_planners.hpp"

namespace ReaK {

namespace pp {

template class rrt_planner<robot_airship::CRS3D_workspace_o1_svp_type>;
template class rrt_planner<robot_airship::CRS3D_workspace_o2_svp_type>;


};

};

#else

namespace ReaK {

namespace pp {

void dummy_CRS_rrt_planners_svp_externs_symbol() { };

};

};

#endif



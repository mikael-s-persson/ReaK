
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

#include "topologies/manip_free_workspace.hpp"

namespace ReaK {

namespace pp {


#define RK_MANIP_FREE_WORKSPACE_MAKE_QUINTIC_QSTAT_ENV_FOR_JOINTS_DEFS(NDOF) \
template class manip_quasi_static_env< Ndof_2nd_order_rl_space<double, NDOF, euclidean_tuple_distance>::type, quintic_hermite_interpolation_tag>;\
\
template class manip_quasi_static_env< Ndof_2nd_order_rl_space<double, NDOF, inf_norm_tuple_distance>::type, quintic_hermite_interpolation_tag>; \
\
template class manip_quasi_static_env< Ndof_rl_space<double, NDOF, 2>::type, quintic_hermite_interpolation_tag>;


RK_MANIP_FREE_WORKSPACE_MAKE_QUINTIC_QSTAT_ENV_FOR_JOINTS_DEFS(1)
RK_MANIP_FREE_WORKSPACE_MAKE_QUINTIC_QSTAT_ENV_FOR_JOINTS_DEFS(2)
RK_MANIP_FREE_WORKSPACE_MAKE_QUINTIC_QSTAT_ENV_FOR_JOINTS_DEFS(3)
RK_MANIP_FREE_WORKSPACE_MAKE_QUINTIC_QSTAT_ENV_FOR_JOINTS_DEFS(4)
RK_MANIP_FREE_WORKSPACE_MAKE_QUINTIC_QSTAT_ENV_FOR_JOINTS_DEFS(5)
RK_MANIP_FREE_WORKSPACE_MAKE_QUINTIC_QSTAT_ENV_FOR_JOINTS_DEFS(6)
RK_MANIP_FREE_WORKSPACE_MAKE_QUINTIC_QSTAT_ENV_FOR_JOINTS_DEFS(7)
RK_MANIP_FREE_WORKSPACE_MAKE_QUINTIC_QSTAT_ENV_FOR_JOINTS_DEFS(8)
RK_MANIP_FREE_WORKSPACE_MAKE_QUINTIC_QSTAT_ENV_FOR_JOINTS_DEFS(9)
RK_MANIP_FREE_WORKSPACE_MAKE_QUINTIC_QSTAT_ENV_FOR_JOINTS_DEFS(10)


template class manip_quasi_static_env< metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type, quintic_hermite_interpolation_tag>;
template class manip_quasi_static_env< metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type, quintic_hermite_interpolation_tag>;

};

};

#else

namespace ReaK {

namespace pp {

void dummy_manip_free_workspace_externs_3_symbol() { };

};

};

#endif



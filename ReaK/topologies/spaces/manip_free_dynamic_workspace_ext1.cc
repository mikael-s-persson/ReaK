
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


#include "ReaK/topologies/spaces/manip_free_dynamic_workspace.h"

namespace ReaK::pp {

#if 0

#define RK_MANIP_FREE_WORKSPACE_MAKE_LIN_DYN_ENV_FOR_JOINTS_DEFS(NDOF)    \
  template class manip_dynamic_env<                                       \
      Ndof_0th_order_rl_space_t<double, NDOF, euclidean_tuple_distance>>; \
  template class manip_dynamic_env<                                       \
      Ndof_1st_order_rl_space_t<double, NDOF, euclidean_tuple_distance>>; \
  template class manip_dynamic_env<                                       \
      Ndof_2nd_order_rl_space_t<double, NDOF, euclidean_tuple_distance>>; \
                                                                          \
  template class manip_dynamic_env<                                       \
      Ndof_0th_order_rl_space_t<double, NDOF, inf_norm_tuple_distance>>;  \
  template class manip_dynamic_env<                                       \
      Ndof_1st_order_rl_space_t<double, NDOF, inf_norm_tuple_distance>>;  \
  template class manip_dynamic_env<                                       \
      Ndof_2nd_order_rl_space_t<double, NDOF, inf_norm_tuple_distance>>;  \
                                                                          \
  template class manip_dynamic_env<Ndof_rl_space_t<double, NDOF, 0>>;     \
  template class manip_dynamic_env<Ndof_rl_space_t<double, NDOF, 1>>;     \
  template class manip_dynamic_env<Ndof_rl_space_t<double, NDOF, 2>>;

#else

#define RK_MANIP_FREE_WORKSPACE_MAKE_LIN_DYN_ENV_FOR_JOINTS_DEFS(NDOF) \
  template class manip_dynamic_env<Ndof_rl_space_t<double, NDOF, 0>>;  \
  template class manip_dynamic_env<Ndof_rl_space_t<double, NDOF, 1>>;  \
  template class manip_dynamic_env<Ndof_rl_space_t<double, NDOF, 2>>;

#endif

RK_MANIP_FREE_WORKSPACE_MAKE_LIN_DYN_ENV_FOR_JOINTS_DEFS(1)
RK_MANIP_FREE_WORKSPACE_MAKE_LIN_DYN_ENV_FOR_JOINTS_DEFS(2)
RK_MANIP_FREE_WORKSPACE_MAKE_LIN_DYN_ENV_FOR_JOINTS_DEFS(3)
RK_MANIP_FREE_WORKSPACE_MAKE_LIN_DYN_ENV_FOR_JOINTS_DEFS(4)
RK_MANIP_FREE_WORKSPACE_MAKE_LIN_DYN_ENV_FOR_JOINTS_DEFS(5)
RK_MANIP_FREE_WORKSPACE_MAKE_LIN_DYN_ENV_FOR_JOINTS_DEFS(6)
RK_MANIP_FREE_WORKSPACE_MAKE_LIN_DYN_ENV_FOR_JOINTS_DEFS(7)
RK_MANIP_FREE_WORKSPACE_MAKE_LIN_DYN_ENV_FOR_JOINTS_DEFS(8)
RK_MANIP_FREE_WORKSPACE_MAKE_LIN_DYN_ENV_FOR_JOINTS_DEFS(9)
RK_MANIP_FREE_WORKSPACE_MAKE_LIN_DYN_ENV_FOR_JOINTS_DEFS(10)

template class manip_dynamic_env<
    metric_space_array_t<se2_0th_order_rl_topology_t<double>, 1>>;
template class manip_dynamic_env<
    metric_space_array_t<se2_1st_order_rl_topology_t<double>, 1>>;
template class manip_dynamic_env<
    metric_space_array_t<se2_2nd_order_rl_topology_t<double>, 1>>;

template class manip_dynamic_env<
    metric_space_array_t<se3_0th_order_rl_topology_t<double>, 1>>;
template class manip_dynamic_env<
    metric_space_array_t<se3_1st_order_rl_topology_t<double>, 1>>;
template class manip_dynamic_env<
    metric_space_array_t<se3_2nd_order_rl_topology_t<double>, 1>>;

}  // namespace ReaK::pp

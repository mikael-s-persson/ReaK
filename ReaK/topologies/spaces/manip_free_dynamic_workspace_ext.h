/**
 * \file manip_free_dynamic_workspace_ext.h
 *
 *
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date December 2012
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

#ifndef REAK_TOPOLOGIES_SPACES_MANIP_FREE_DYNAMIC_WORKSPACE_EXT_H_
#define REAK_TOPOLOGIES_SPACES_MANIP_FREE_DYNAMIC_WORKSPACE_EXT_H_


#include "ReaK/topologies/spaces/differentiable_space.h"
#include "ReaK/topologies/spaces/metric_space_tuple.h"
#include "ReaK/topologies/spaces/rate_limited_spaces.h"

#include "ReaK/topologies/spaces/ndof_spaces.h"
#include "ReaK/topologies/spaces/se2_topologies.h"
#include "ReaK/topologies/spaces/se3_topologies.h"

namespace ReaK::pp {

#define RK_MANIP_FREE_WORKSPACE_MAKE_DYN_ENV_FOR_JOINTS(NDOF)                \
  extern template class manip_dynamic_env<Ndof_rl_space_t<double, NDOF, 0>>; \
  extern template class manip_dynamic_env<Ndof_rl_space_t<double, NDOF, 1>>; \
  extern template class manip_dynamic_env<Ndof_rl_space_t<double, NDOF, 2>>;

RK_MANIP_FREE_WORKSPACE_MAKE_DYN_ENV_FOR_JOINTS(1)
RK_MANIP_FREE_WORKSPACE_MAKE_DYN_ENV_FOR_JOINTS(2)
RK_MANIP_FREE_WORKSPACE_MAKE_DYN_ENV_FOR_JOINTS(3)
RK_MANIP_FREE_WORKSPACE_MAKE_DYN_ENV_FOR_JOINTS(4)
RK_MANIP_FREE_WORKSPACE_MAKE_DYN_ENV_FOR_JOINTS(5)
RK_MANIP_FREE_WORKSPACE_MAKE_DYN_ENV_FOR_JOINTS(6)
RK_MANIP_FREE_WORKSPACE_MAKE_DYN_ENV_FOR_JOINTS(7)
RK_MANIP_FREE_WORKSPACE_MAKE_DYN_ENV_FOR_JOINTS(8)
RK_MANIP_FREE_WORKSPACE_MAKE_DYN_ENV_FOR_JOINTS(9)
RK_MANIP_FREE_WORKSPACE_MAKE_DYN_ENV_FOR_JOINTS(10)

extern template class manip_dynamic_env<
    metric_space_array_t<se2_0th_order_rl_topology_t<double>, 1>>;
extern template class manip_dynamic_env<
    metric_space_array_t<se2_1st_order_rl_topology_t<double>, 1>>;
extern template class manip_dynamic_env<
    metric_space_array_t<se2_2nd_order_rl_topology_t<double>, 1>>;

extern template class manip_dynamic_env<
    metric_space_array_t<se3_0th_order_rl_topology_t<double>, 1>>;
extern template class manip_dynamic_env<
    metric_space_array_t<se3_1st_order_rl_topology_t<double>, 1>>;
extern template class manip_dynamic_env<
    metric_space_array_t<se3_2nd_order_rl_topology_t<double>, 1>>;

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_MANIP_FREE_DYNAMIC_WORKSPACE_EXT_H_

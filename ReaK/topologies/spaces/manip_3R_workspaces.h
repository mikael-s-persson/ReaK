/**
 * \file manip_3R_workspaces.h
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

#ifndef REAK_TOPOLOGIES_SPACES_MANIP_3R_WORKSPACES_H_
#define REAK_TOPOLOGIES_SPACES_MANIP_3R_WORKSPACES_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/mbd/models/manip_3R_arm.h"
#include "ReaK/topologies/spaces/manip_planning_traits.h"
#include "ReaK/topologies/spaces/se2_topologies.h"
#include "ReaK/topologies/spaces/se3_topologies.h"

namespace ReaK::pp {

template <int Order>
struct manip_pp_traits<kte::manip_3R_2D_kinematics, Order> {
  static constexpr std::size_t degrees_of_freedom = 3;

  using rl_jt_space_type = Ndof_rl_space_t<double, 3, Order>;
  using jt_space_type = Ndof_space_t<double, 3, Order>;
  using ee_space_type = se2_topology_t<double, Order>;
};

template <int Order>
struct manip_pp_traits<kte::manip_3R_3D_kinematics, Order> {
  static constexpr std::size_t degrees_of_freedom = 3;

  using rl_jt_space_type = Ndof_rl_space_t<double, 3, Order>;
  using jt_space_type = Ndof_space_t<double, 3, Order>;
  using ee_space_type = se3_topology_t<double, Order>;
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_MANIP_3R_WORKSPACES_H_

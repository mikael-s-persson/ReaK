
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

#include "ReaK/core/base/defs.hpp"

#include "ReaK/topologies/spaces/joint_space_limits.tpp"
#include "ReaK/topologies/spaces/joint_space_topologies.hpp"

namespace ReaK::pp {

#define RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(NDOF)                 \
                                                                             \
  template Ndof_0th_order_rl_space_t<double, NDOF>                           \
  joint_limits_mapping<double>::make_rl_joint_space(                         \
      const Ndof_0th_order_space_t<double, NDOF>&) const;                    \
  template Ndof_1st_order_rl_space_t<double, NDOF>                           \
  joint_limits_mapping<double>::make_rl_joint_space(                         \
      const Ndof_1st_order_space_t<double, NDOF>&) const;                    \
  template Ndof_2nd_order_rl_space_t<double, NDOF>                           \
  joint_limits_mapping<double>::make_rl_joint_space(                         \
      const Ndof_2nd_order_space_t<double, NDOF>&) const;                    \
                                                                             \
  template Ndof_0th_order_space_t<double, NDOF>                              \
  joint_limits_mapping<double>::make_normal_joint_space(                     \
      const Ndof_0th_order_rl_space_t<double, NDOF>&) const;                 \
  template Ndof_1st_order_space_t<double, NDOF>                              \
  joint_limits_mapping<double>::make_normal_joint_space(                     \
      const Ndof_1st_order_rl_space_t<double, NDOF>&) const;                 \
  template Ndof_2nd_order_space_t<double, NDOF>                              \
  joint_limits_mapping<double>::make_normal_joint_space(                     \
      const Ndof_2nd_order_rl_space_t<double, NDOF>&) const;                 \
                                                                             \
  template topology_point_type_t<Ndof_0th_order_rl_space_t<double, NDOF>>    \
  joint_limits_mapping<double>::map_to_space(                                \
      const topology_point_type_t<Ndof_0th_order_space_t<double, NDOF>>& pt, \
      const Ndof_0th_order_space_t<double, NDOF>&,                           \
      const Ndof_0th_order_rl_space_t<double, NDOF>&) const;                 \
  template topology_point_type_t<Ndof_1st_order_rl_space_t<double, NDOF>>    \
  joint_limits_mapping<double>::map_to_space(                                \
      const topology_point_type_t<Ndof_1st_order_space_t<double, NDOF>>& pt, \
      const Ndof_1st_order_space_t<double, NDOF>&,                           \
      const Ndof_1st_order_rl_space_t<double, NDOF>&) const;                 \
  template topology_point_type_t<Ndof_2nd_order_rl_space_t<double, NDOF>>    \
  joint_limits_mapping<double>::map_to_space(                                \
      const topology_point_type_t<Ndof_2nd_order_space_t<double, NDOF>>& pt, \
      const Ndof_2nd_order_space_t<double, NDOF>&,                           \
      const Ndof_2nd_order_rl_space_t<double, NDOF>&) const;                 \
                                                                             \
  template topology_point_type_t<Ndof_0th_order_space_t<double, NDOF>>       \
  joint_limits_mapping<double>::map_to_space(                                \
      const topology_point_type_t<Ndof_0th_order_rl_space_t<double, NDOF>>&  \
          pt,                                                                \
      const Ndof_0th_order_rl_space_t<double, NDOF>&,                        \
      const Ndof_0th_order_space_t<double, NDOF>&) const;                    \
  template topology_point_type_t<Ndof_1st_order_space_t<double, NDOF>>       \
  joint_limits_mapping<double>::map_to_space(                                \
      const topology_point_type_t<Ndof_1st_order_rl_space_t<double, NDOF>>&  \
          pt,                                                                \
      const Ndof_1st_order_rl_space_t<double, NDOF>&,                        \
      const Ndof_1st_order_space_t<double, NDOF>&) const;                    \
  template topology_point_type_t<Ndof_2nd_order_space_t<double, NDOF>>       \
  joint_limits_mapping<double>::map_to_space(                                \
      const topology_point_type_t<Ndof_2nd_order_rl_space_t<double, NDOF>>&  \
          pt,                                                                \
      const Ndof_2nd_order_rl_space_t<double, NDOF>&,                        \
      const Ndof_2nd_order_space_t<double, NDOF>&) const;

RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(1)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(2)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(3)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(4)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(5)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(6)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(7)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(8)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(9)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(10)

}  // namespace ReaK::pp

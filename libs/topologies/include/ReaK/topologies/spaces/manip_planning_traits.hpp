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

#include "ReaK/core/base/defs.hpp"

#include "ReaK/topologies/spaces/manip_free_dynamic_workspace.hpp"
#include "ReaK/topologies/spaces/manip_free_workspace.hpp"

#include "ReaK/topologies/spaces/direct_kinematics_topomap.hpp"
#include "ReaK/topologies/spaces/inverse_kinematics_topomap.hpp"

#include <type_traits>

namespace ReaK::pp {

template <typename ManipMdlType, int Order>
struct manip_pp_traits {
  static constexpr std::size_t degrees_of_freedom =
      ManipMdlType::degrees_of_freedom;

  //   Ndof_rl_space<double, ManipMdlType::degrees_of_freedom, Order>::type
  using rl_jt_space_type = std::conditional_t<
      (Order == 0), typename ManipMdlType::rl_o0_jt_space_type,
      std::conditional_t<(Order == 1),
                         typename ManipMdlType::rl_o1_jt_space_type,
                         typename ManipMdlType::rl_o2_jt_space_type>>;

  //   Ndof_space<double, ManipMdlType::degrees_of_freedom, Order>::type
  using jt_space_type = std::conditional_t<
      (Order == 0), typename ManipMdlType::o0_jt_space_type,
      std::conditional_t<(Order == 1), typename ManipMdlType::o1_jt_space_type,
                         typename ManipMdlType::o2_jt_space_type>>;

  //   se3_0th_order_topology<double>::type
  //   se3_1st_order_topology<double>::type
  //   se3_2nd_order_topology<double>::type
  using ee_space_type = std::conditional_t<
      (Order == 0), typename ManipMdlType::o0_ee_space_type,
      std::conditional_t<(Order == 1), typename ManipMdlType::o1_ee_space_type,
                         typename ManipMdlType::o2_ee_space_type>>;
};

template <typename ManipMdlType, int Order>
struct manip_static_workspace {
  using Traits = manip_pp_traits<ManipMdlType, Order>;

  using rl_workspace_type =
      manip_quasi_static_env<typename Traits::rl_jt_space_type>;
  using workspace_type = manip_quasi_static_env<typename Traits::jt_space_type>;
};

template <typename ManipMdlType, int Order>
struct manip_dynamic_workspace {
  using Traits = manip_pp_traits<ManipMdlType, Order>;

  using rl_workspace_type =
      manip_dynamic_env<typename Traits::rl_jt_space_type>;
};

template <typename ManipMdlType, int Order>
struct manip_DK_map {
  using Traits = manip_pp_traits<ManipMdlType, Order>;

  using rl_map_type = manip_rl_direct_kin_map<joint_limits_mapping<double>,
                                              typename Traits::jt_space_type>;
  using map_type = manip_direct_kin_map;
};

template <typename ManipMdlType, int Order>
struct manip_IK_map {
  using rl_map_type = manip_rl_inverse_kin_map<joint_limits_mapping<double>>;
  using map_type = manip_inverse_kin_map;
};

template <int Order, typename ManipMdlType>
auto make_Ndof_space(
    const std::shared_ptr<ManipMdlType>& manip_kin_mdl,
    const std::shared_ptr<kte::joint_limits_collection<double>>&
        manip_jt_limits) {
  if constexpr (Order == 0) {
    return make_Ndof_space<
        manip_pp_traits<ManipMdlType, Order>::degrees_of_freedom>(
        manip_kin_mdl->getJointPositionLowerBounds(),
        manip_kin_mdl->getJointPositionUpperBounds());
  } else if constexpr (Order == 1) {
    return make_Ndof_space<
        manip_pp_traits<ManipMdlType, Order>::degrees_of_freedom>(
        manip_kin_mdl->getJointPositionLowerBounds(),
        manip_kin_mdl->getJointPositionUpperBounds(),
        manip_jt_limits->gen_speed_limits);
  } else {
    return make_Ndof_space<
        manip_pp_traits<ManipMdlType, Order>::degrees_of_freedom>(
        manip_kin_mdl->getJointPositionLowerBounds(),
        manip_kin_mdl->getJointPositionUpperBounds(),
        manip_jt_limits->gen_speed_limits, manip_jt_limits->gen_accel_limits);
  }
}

template <int Order, typename ManipMdlType>
auto make_manip_jt_space(
    const std::shared_ptr<ManipMdlType>& manip_kin_mdl,
    const std::shared_ptr<kte::joint_limits_collection<double>>&
        manip_jt_limits) {
  using JtspaceType =
      typename manip_pp_traits<ManipMdlType, Order>::jt_space_type;
  return std::make_shared<JtspaceType>(
      make_Ndof_space<Order>(manip_kin_mdl, manip_jt_limits));
}

template <int Order, typename ManipMdlType>
auto make_Ndof_rl_space(
    const std::shared_ptr<ManipMdlType>& manip_kin_mdl,
    const std::shared_ptr<kte::joint_limits_collection<double>>&
        manip_jt_limits) {
  if constexpr (Order == 0) {
    return make_Ndof_rl_space<
        manip_pp_traits<ManipMdlType, Order>::degrees_of_freedom>(
        manip_kin_mdl->getJointPositionLowerBounds(),
        manip_kin_mdl->getJointPositionUpperBounds(),
        manip_jt_limits->gen_speed_limits);
  } else if constexpr (Order == 1) {
    return make_Ndof_rl_space<
        manip_pp_traits<ManipMdlType, Order>::degrees_of_freedom>(
        manip_kin_mdl->getJointPositionLowerBounds(),
        manip_kin_mdl->getJointPositionUpperBounds(),
        manip_jt_limits->gen_speed_limits, manip_jt_limits->gen_accel_limits);
  } else {
    return make_Ndof_rl_space<
        manip_pp_traits<ManipMdlType, Order>::degrees_of_freedom>(
        manip_kin_mdl->getJointPositionLowerBounds(),
        manip_kin_mdl->getJointPositionUpperBounds(),
        manip_jt_limits->gen_speed_limits, manip_jt_limits->gen_accel_limits,
        manip_jt_limits->gen_jerk_limits);
  }
}

template <int Order, typename ManipMdlType>
auto make_manip_rl_jt_space(
    const std::shared_ptr<ManipMdlType>& manip_kin_mdl,
    const std::shared_ptr<kte::joint_limits_collection<double>>&
        manip_jt_limits) {
  using JtspaceType =
      typename manip_pp_traits<ManipMdlType, Order>::rl_jt_space_type;
  return std::make_shared<JtspaceType>(
      make_Ndof_rl_space<Order>(manip_kin_mdl, manip_jt_limits));
}

template <int Order, typename InterpTag, typename ManipMdlType>
auto make_manip_static_workspace(
    InterpTag interp_tag, const std::shared_ptr<ManipMdlType>& manip_kin_mdl,
    const std::shared_ptr<kte::joint_limits_collection<double>>&
        manip_jt_limits,
    double min_travel) {
  using WorkspaceType =
      typename manip_static_workspace<ManipMdlType, Order>::rl_workspace_type;
  using JtspaceType =
      typename manip_pp_traits<ManipMdlType, Order>::jt_space_type;
  using RLJtspaceType =
      typename manip_pp_traits<ManipMdlType, Order>::rl_jt_space_type;

  return std::make_shared<WorkspaceType>(
      interp_tag, make_Ndof_rl_space<Order>(manip_kin_mdl, manip_jt_limits),
      make_any_model_applicator<RLJtspaceType>(
          manip_rl_direct_kin_map<joint_limits_mapping<double>, JtspaceType>(
              manip_kin_mdl, joint_limits_mapping<double>(manip_jt_limits),
              make_manip_jt_space<Order>(manip_kin_mdl, manip_jt_limits))),
      min_travel);
}

template <int Order, typename InterpTag, typename ManipMdlType>
auto make_manip_dynamic_workspace(
    InterpTag interp_tag, const std::shared_ptr<ManipMdlType>& manip_kin_mdl,
    const std::shared_ptr<kte::joint_limits_collection<double>>&
        manip_jt_limits,
    double min_travel, double max_horizon, double time_delay) {
  using WorkspaceType =
      typename manip_dynamic_workspace<ManipMdlType, Order>::rl_workspace_type;
  using JtspaceType =
      typename manip_pp_traits<ManipMdlType, Order>::jt_space_type;
  using RLJtspaceType =
      typename manip_pp_traits<ManipMdlType, Order>::rl_jt_space_type;

  return std::make_shared<WorkspaceType>(
      interp_tag, make_Ndof_rl_space<Order>(manip_kin_mdl, manip_jt_limits),
      make_any_model_applicator<RLJtspaceType>(
          manip_rl_direct_kin_map<joint_limits_mapping<double>, JtspaceType>(
              manip_kin_mdl, joint_limits_mapping<double>(manip_jt_limits),
              make_manip_jt_space<Order>(manip_kin_mdl, manip_jt_limits))),
      min_travel, max_horizon, time_delay);
}

}  // namespace ReaK::pp

#endif

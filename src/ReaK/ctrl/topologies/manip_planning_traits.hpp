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
#include "manip_free_dynamic_workspace.hpp"

#include "direct_kinematics_topomap.hpp"
#include "inverse_kinematics_topomap.hpp"

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/equal_to.hpp>

namespace ReaK {

namespace pp {


template <typename ManipMdlType, int Order>
struct manip_pp_traits {
  BOOST_STATIC_CONSTANT(std::size_t, degrees_of_freedom = ManipMdlType::degrees_of_freedom);
  
//   Ndof_rl_space<double, ManipMdlType::degrees_of_freedom, Order>::type
  typedef typename boost::mpl::if_<
    boost::mpl::equal_to< boost::mpl::int_<0>, boost::mpl::int_<Order> >,
    typename ManipMdlType::rl_o0_jt_space_type,
    typename boost::mpl::if_<
      boost::mpl::equal_to< boost::mpl::int_<1>, boost::mpl::int_<Order> >,
      typename ManipMdlType::rl_o1_jt_space_type,
      typename ManipMdlType::rl_o2_jt_space_type 
    >::type
  >::type rl_jt_space_type;
  
//   Ndof_space<double, ManipMdlType::degrees_of_freedom, Order>::type
  typedef typename boost::mpl::if_<
    boost::mpl::equal_to< boost::mpl::int_<0>, boost::mpl::int_<Order> >,
    typename ManipMdlType::o0_jt_space_type,
    typename boost::mpl::if_<
      boost::mpl::equal_to< boost::mpl::int_<1>, boost::mpl::int_<Order> >,
      typename ManipMdlType::o1_jt_space_type,
      typename ManipMdlType::o2_jt_space_type 
    >::type
  >::type jt_space_type;
  
//   se3_0th_order_topology<double>::type
//   se3_1st_order_topology<double>::type
//   se3_2nd_order_topology<double>::type
  typedef typename boost::mpl::if_<
    boost::mpl::equal_to< boost::mpl::int_<0>, boost::mpl::int_<Order> >,
    typename ManipMdlType::o0_ee_space_type,
    typename boost::mpl::if_<
      boost::mpl::equal_to< boost::mpl::int_<1>, boost::mpl::int_<Order> >,
      typename ManipMdlType::o1_ee_space_type,
      typename ManipMdlType::o2_ee_space_type 
    >::type
  >::type ee_space_type;
  
};


template <typename ManipMdlType, int Order>
struct manip_static_workspace {
  typedef manip_pp_traits<ManipMdlType, Order> Traits;
  
  typedef manip_quasi_static_env<typename Traits::rl_jt_space_type>  rl_workspace_type;
  typedef manip_quasi_static_env<typename Traits::jt_space_type>     workspace_type;
  
};


template <typename ManipMdlType, int Order>
struct manip_dynamic_workspace {
  typedef manip_pp_traits<ManipMdlType, Order> Traits;
  
  typedef manip_dynamic_env<typename Traits::rl_jt_space_type> rl_workspace_type;
};


template <typename ManipMdlType, int Order>
struct manip_DK_map {
  typedef manip_pp_traits<ManipMdlType, Order> Traits;
  
  typedef manip_rl_direct_kin_map< joint_limits_collection<double>, typename Traits::jt_space_type > rl_map_type;
  typedef manip_direct_kin_map map_type;
};


template <typename ManipMdlType, int Order>
struct manip_IK_map {
  typedef manip_rl_inverse_kin_map< joint_limits_collection<double> > rl_map_type;
  typedef manip_inverse_kin_map map_type;
};






template <int Order, typename ManipMdlType>
typename boost::enable_if< boost::mpl::equal_to< boost::mpl::int_<0>, boost::mpl::int_<Order> >,
  shared_ptr< typename manip_pp_traits< ManipMdlType, Order >::jt_space_type > 
>::type make_manip_jt_space(
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >&) {
  typedef typename manip_pp_traits< ManipMdlType, Order >::jt_space_type JtspaceType;
  
  return shared_ptr<JtspaceType>(new JtspaceType(
    make_Ndof_space< manip_pp_traits< ManipMdlType, Order >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds())));
};

template <int Order, typename ManipMdlType>
typename boost::enable_if< boost::mpl::equal_to< boost::mpl::int_<1>, boost::mpl::int_<Order> >,
  shared_ptr< typename manip_pp_traits< ManipMdlType, Order >::jt_space_type > 
>::type make_manip_jt_space(
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits) {
  typedef typename manip_pp_traits< ManipMdlType, Order >::jt_space_type JtspaceType;
  
  return shared_ptr<JtspaceType>(new JtspaceType(
    make_Ndof_space< manip_pp_traits< ManipMdlType, Order >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits)));
};

template <int Order, typename ManipMdlType>
typename boost::enable_if< boost::mpl::equal_to< boost::mpl::int_<2>, boost::mpl::int_<Order> >,
  shared_ptr< typename manip_pp_traits< ManipMdlType, Order >::jt_space_type > 
>::type make_manip_jt_space(
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits) {
  typedef typename manip_pp_traits< ManipMdlType, Order >::jt_space_type JtspaceType;
  
  return shared_ptr<JtspaceType>(new JtspaceType(
    make_Ndof_space< manip_pp_traits< ManipMdlType, Order >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits, 
      manip_jt_limits->gen_accel_limits)));
};






template <int Order, typename ManipMdlType>
typename boost::enable_if< boost::mpl::equal_to< boost::mpl::int_<0>, boost::mpl::int_<Order> >,
  shared_ptr< typename manip_pp_traits< ManipMdlType, Order >::rl_jt_space_type > 
>::type make_manip_rl_jt_space(
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits) {
  typedef typename manip_pp_traits< ManipMdlType, Order >::rl_jt_space_type JtspaceType;
  
  return shared_ptr<JtspaceType>(new JtspaceType(
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType, Order >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits)));
};

template <int Order, typename ManipMdlType>
typename boost::enable_if< boost::mpl::equal_to< boost::mpl::int_<1>, boost::mpl::int_<Order> >,
  shared_ptr< typename manip_pp_traits< ManipMdlType, Order >::rl_jt_space_type > 
>::type make_manip_rl_jt_space(
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits) {
  typedef typename manip_pp_traits< ManipMdlType, Order >::rl_jt_space_type JtspaceType;
  
  return shared_ptr<JtspaceType>(new JtspaceType(
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType, Order >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits, 
      manip_jt_limits->gen_accel_limits)));
};

template <int Order, typename ManipMdlType>
typename boost::enable_if< boost::mpl::equal_to< boost::mpl::int_<2>, boost::mpl::int_<Order> >,
  shared_ptr< typename manip_pp_traits< ManipMdlType, Order >::rl_jt_space_type > 
>::type make_manip_rl_jt_space(
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits) {
  typedef typename manip_pp_traits< ManipMdlType, Order >::rl_jt_space_type JtspaceType;
  
  return shared_ptr<JtspaceType>(new JtspaceType(
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType, Order >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits, 
      manip_jt_limits->gen_accel_limits, 
      manip_jt_limits->gen_jerk_limits)));
};








template <int Order, typename InterpTag, typename ManipMdlType>
typename boost::enable_if< boost::mpl::equal_to< boost::mpl::int_<0>, boost::mpl::int_<Order> >,
  shared_ptr< typename manip_static_workspace< ManipMdlType, Order >::rl_workspace_type > 
>::type make_manip_static_workspace(
    InterpTag interp_tag,
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits,
    double min_travel) {
  typedef typename manip_static_workspace< ManipMdlType, Order >::rl_workspace_type WorkspaceType;
  typedef typename manip_pp_traits< ManipMdlType, Order >::jt_space_type JtspaceType;
  typedef typename manip_pp_traits< ManipMdlType, Order >::rl_jt_space_type RLJtspaceType;
  
  return shared_ptr<WorkspaceType>(new WorkspaceType(
    interp_tag,
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType, Order >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits),
    make_any_model_applicator<RLJtspaceType>(manip_rl_direct_kin_map< joint_limits_collection<double>, JtspaceType >(
      manip_kin_mdl, manip_jt_limits, make_manip_jt_space<Order>(manip_kin_mdl, manip_jt_limits))),
    min_travel));
};

template <int Order, typename InterpTag, typename ManipMdlType>
typename boost::enable_if< boost::mpl::equal_to< boost::mpl::int_<1>, boost::mpl::int_<Order> >,
  shared_ptr< typename manip_static_workspace< ManipMdlType, Order >::rl_workspace_type > 
>::type make_manip_static_workspace(
    InterpTag interp_tag,
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits,
    double min_travel) {
  typedef typename manip_static_workspace< ManipMdlType, Order >::rl_workspace_type WorkspaceType;
  typedef typename manip_pp_traits< ManipMdlType, Order >::jt_space_type JtspaceType;
  typedef typename manip_pp_traits< ManipMdlType, Order >::rl_jt_space_type RLJtspaceType;
  
  return shared_ptr<WorkspaceType>(new WorkspaceType(
    interp_tag,
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType, Order >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits, 
      manip_jt_limits->gen_accel_limits),
    make_any_model_applicator<RLJtspaceType>(manip_rl_direct_kin_map< joint_limits_collection<double>, JtspaceType >(
      manip_kin_mdl, manip_jt_limits, make_manip_jt_space<Order>(manip_kin_mdl, manip_jt_limits))),
    min_travel));
};

template <int Order, typename InterpTag, typename ManipMdlType>
typename boost::enable_if< boost::mpl::equal_to< boost::mpl::int_<2>, boost::mpl::int_<Order> >,
  shared_ptr< typename manip_static_workspace< ManipMdlType, Order >::rl_workspace_type > 
>::type make_manip_static_workspace(
    InterpTag interp_tag,
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits,
    double min_travel) {
  typedef typename manip_static_workspace< ManipMdlType, Order >::rl_workspace_type WorkspaceType;
  typedef typename manip_pp_traits< ManipMdlType, Order >::jt_space_type JtspaceType;
  typedef typename manip_pp_traits< ManipMdlType, Order >::rl_jt_space_type RLJtspaceType;
  
  return shared_ptr<WorkspaceType>(new WorkspaceType(
    interp_tag,
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType, Order >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits, 
      manip_jt_limits->gen_accel_limits, 
      manip_jt_limits->gen_jerk_limits),
    make_any_model_applicator<RLJtspaceType>(manip_rl_direct_kin_map< joint_limits_collection<double>, JtspaceType >(
      manip_kin_mdl, manip_jt_limits, make_manip_jt_space<Order>(manip_kin_mdl, manip_jt_limits))),
    min_travel));
};


template <int Order, typename InterpTag, typename ManipMdlType>
typename boost::enable_if< boost::mpl::equal_to< boost::mpl::int_<0>, boost::mpl::int_<Order> >,
  shared_ptr< typename manip_dynamic_workspace< ManipMdlType, Order >::rl_workspace_type > 
>::type make_manip_dynamic_workspace(
    InterpTag interp_tag,
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits,
    double min_travel, double max_horizon, double time_delay) {
  typedef typename manip_dynamic_workspace< ManipMdlType, Order >::rl_workspace_type WorkspaceType;
  typedef typename manip_pp_traits< ManipMdlType, Order >::jt_space_type JtspaceType;
  typedef typename manip_pp_traits< ManipMdlType, Order >::rl_jt_space_type RLJtspaceType;
  
  return shared_ptr<WorkspaceType>(new WorkspaceType(
    interp_tag,
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType, Order >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits),
    make_any_model_applicator<RLJtspaceType>(manip_rl_direct_kin_map< joint_limits_collection<double>, JtspaceType >(
      manip_kin_mdl, manip_jt_limits, make_manip_jt_space<Order>(manip_kin_mdl, manip_jt_limits))),
    min_travel, max_horizon, time_delay));
};

template <int Order, typename InterpTag, typename ManipMdlType>
typename boost::enable_if< boost::mpl::equal_to< boost::mpl::int_<1>, boost::mpl::int_<Order> >,
  shared_ptr< typename manip_dynamic_workspace< ManipMdlType, Order >::rl_workspace_type > 
>::type make_manip_dynamic_workspace(
    InterpTag interp_tag,
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits,
    double min_travel, double max_horizon, double time_delay) {
  typedef typename manip_dynamic_workspace< ManipMdlType, Order >::rl_workspace_type WorkspaceType;
  typedef typename manip_pp_traits< ManipMdlType, Order >::jt_space_type JtspaceType;
  typedef typename manip_pp_traits< ManipMdlType, Order >::rl_jt_space_type RLJtspaceType;
  
  return shared_ptr<WorkspaceType>(new WorkspaceType(
    interp_tag,
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType, Order >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits, 
      manip_jt_limits->gen_accel_limits),
    make_any_model_applicator<RLJtspaceType>(manip_rl_direct_kin_map< joint_limits_collection<double>, JtspaceType >(
      manip_kin_mdl, manip_jt_limits, make_manip_jt_space<Order>(manip_kin_mdl, manip_jt_limits))),
    min_travel, max_horizon, time_delay));
};

template <int Order, typename InterpTag, typename ManipMdlType>
typename boost::enable_if< boost::mpl::equal_to< boost::mpl::int_<2>, boost::mpl::int_<Order> >,
  shared_ptr< typename manip_dynamic_workspace< ManipMdlType, Order >::rl_workspace_type > 
>::type make_manip_dynamic_workspace(
    InterpTag interp_tag,
    const shared_ptr< ManipMdlType >& manip_kin_mdl,
    const shared_ptr< joint_limits_collection<double> >& manip_jt_limits,
    double min_travel, double max_horizon, double time_delay) {
  typedef typename manip_dynamic_workspace< ManipMdlType, Order >::rl_workspace_type WorkspaceType;
  typedef typename manip_pp_traits< ManipMdlType, Order >::jt_space_type JtspaceType;
  typedef typename manip_pp_traits< ManipMdlType, Order >::rl_jt_space_type RLJtspaceType;
  
  return shared_ptr<WorkspaceType>(new WorkspaceType(
    interp_tag,
    make_Ndof_rl_space< manip_pp_traits< ManipMdlType, Order >::degrees_of_freedom >(
      manip_kin_mdl->getJointPositionLowerBounds(), 
      manip_kin_mdl->getJointPositionUpperBounds(), 
      manip_jt_limits->gen_speed_limits, 
      manip_jt_limits->gen_accel_limits, 
      manip_jt_limits->gen_jerk_limits),
    make_any_model_applicator<RLJtspaceType>(manip_rl_direct_kin_map< joint_limits_collection<double>, JtspaceType >(
      manip_kin_mdl, manip_jt_limits, make_manip_jt_space<Order>(manip_kin_mdl, manip_jt_limits))),
    min_travel, max_horizon, time_delay));
};





};

};

#endif


/**
 * \file manip_SSRMS_workspaces.hpp
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

#ifndef REAK_MANIP_SSRMS_WORKSPACES_HPP
#define REAK_MANIP_SSRMS_WORKSPACES_HPP

#include "base/defs.hpp"

#include "kte_models/manip_SSRMS_arm.hpp"
#include "topologies/manip_planning_traits.hpp"
#include "topologies/se3_topologies.hpp"

namespace ReaK {

namespace pp {
  


template <>
struct manip_pp_traits< kte::manip_SSRMS_kinematics > {
  BOOST_STATIC_CONSTANT(std::size_t, degrees_of_freedom = 7);
  
  typedef Ndof_rl_space<double, 7, 0>::type rl_o0_jt_space_type;
  typedef Ndof_space<double, 7, 0>::type o0_jt_space_type;
  typedef se3_0th_order_topology<double>::type o0_ee_space_type;
  
  typedef Ndof_rl_space<double, 7, 1>::type rl_o1_jt_space_type;
  typedef Ndof_space<double, 7, 1>::type o1_jt_space_type;
  typedef se3_1st_order_topology<double>::type o1_ee_space_type;
  
  typedef Ndof_rl_space<double, 7, 2>::type rl_o2_jt_space_type;
  typedef Ndof_space<double, 7, 2>::type o2_jt_space_type;
  typedef se3_2nd_order_topology<double>::type o2_ee_space_type;
};


};

};

#endif


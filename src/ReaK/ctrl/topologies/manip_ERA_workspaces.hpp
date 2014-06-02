/**
 * \file manip_ERA_workspaces.hpp
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

#ifndef REAK_MANIP_ERA_WORKSPACES_HPP
#define REAK_MANIP_ERA_WORKSPACES_HPP

#include <ReaK/core/base/defs.hpp>

#include <ReaK/ctrl/kte_models/manip_ERA_arm.hpp>
#include <ReaK/ctrl/topologies/manip_planning_traits.hpp>
#include <ReaK/ctrl/topologies/se3_topologies.hpp>

namespace ReaK {

namespace pp {
  


template <int Order>
struct manip_pp_traits< kte::manip_ERA_kinematics, Order > {
  BOOST_STATIC_CONSTANT(std::size_t, degrees_of_freedom = 7);
  
  typedef typename Ndof_rl_space<double, 7, Order>::type rl_jt_space_type;
  typedef typename Ndof_space<double, 7, Order>::type jt_space_type;
  typedef typename se3_topology<double, Order>::type ee_space_type;
};


};

};

#endif


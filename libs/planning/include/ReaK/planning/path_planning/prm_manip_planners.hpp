/**
 * \file prm_manip_planners.hpp
 *
 * This library defines extern template declarations for prm planner instantiations.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date January 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_PRM_MANIP_PLANNERS_HPP
#define REAK_PRM_MANIP_PLANNERS_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/named_object.hpp>

#include <ReaK/topologies/spaces/Ndof_spaces.hpp>

#include <ReaK/topologies/spaces/manip_free_dynamic_workspace.hpp>
#include <ReaK/topologies/spaces/manip_free_workspace.hpp>

#include <ReaK/planning/path_planning/prm_path_planner.hpp>

namespace ReaK::pp {

// extern template class prm_planner<WORKSPACE>;

#define RK_PRM_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(NDOF) \
  extern template class prm_planner<                              \
      manip_quasi_static_env<Ndof_rl_space_t<double, NDOF, 0>>>;  \
  extern template class prm_planner<                              \
      manip_quasi_static_env<Ndof_rl_space_t<double, NDOF, 1>>>;  \
  extern template class prm_planner<                              \
      manip_quasi_static_env<Ndof_rl_space_t<double, NDOF, 2>>>;

RK_PRM_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(0)
RK_PRM_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(1)
RK_PRM_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(2)
RK_PRM_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(3)
RK_PRM_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(4)
RK_PRM_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(5)
RK_PRM_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(6)
RK_PRM_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(7)
RK_PRM_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(8)
RK_PRM_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(9)
RK_PRM_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(10)

#define RK_PRM_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(NDOF) \
  extern template class prm_planner<                               \
      manip_dynamic_env<Ndof_rl_space_t<double, NDOF, 0>>>;        \
  extern template class prm_planner<                               \
      manip_dynamic_env<Ndof_rl_space_t<double, NDOF, 1>>>;        \
  extern template class prm_planner<                               \
      manip_dynamic_env<Ndof_rl_space_t<double, NDOF, 2>>>;

RK_PRM_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(0)
RK_PRM_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(1)
RK_PRM_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(2)
RK_PRM_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(3)
RK_PRM_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(4)
RK_PRM_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(5)
RK_PRM_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(6)
RK_PRM_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(7)
RK_PRM_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(8)
RK_PRM_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(9)
RK_PRM_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(10)

}  // namespace ReaK::pp

#endif

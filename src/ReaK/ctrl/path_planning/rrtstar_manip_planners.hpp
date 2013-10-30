/**
 * \file rrtstar_manip_planners.hpp
 * 
 * This library defines extern template declarations for rrtstar planner instantiations.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2013
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

#ifndef REAK_RRTSTAR_MANIP_PLANNERS_HPP
#define REAK_RRTSTAR_MANIP_PLANNERS_HPP

#include "base/defs.hpp"
#include "base/named_object.hpp"

#include "topologies/Ndof_spaces.hpp"

#include "topologies/manip_free_workspace.hpp"
#include "topologies/manip_free_svp_Ndof_workspace.hpp"
#include "topologies/manip_free_sap_Ndof_workspace.hpp"

#include "topologies/manip_free_dynamic_workspace.hpp"
#include "topologies/manip_free_svp_Ndof_dynamic_workspace.hpp"
#include "topologies/manip_free_sap_Ndof_dynamic_workspace.hpp"

#include "path_planning/rrtstar_path_planner.hpp"



#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

namespace ReaK {

namespace pp {

//extern template class rrtstar_planner<WORKSPACE>;


#define RK_RRTSTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(NDOF) \
extern template class rrtstar_planner< \
  manip_quasi_static_env< typename Ndof_rl_space<double, NDOF, 0>::type, linear_interpolation_tag> >; \
extern template class rrtstar_planner< \
  manip_quasi_static_env< typename Ndof_rl_space<double, NDOF, 1>::type, linear_interpolation_tag> >; \
extern template class rrtstar_planner< \
  manip_quasi_static_env< typename Ndof_rl_space<double, NDOF, 2>::type, linear_interpolation_tag> >; \
 \
extern template class rrtstar_planner< \
  manip_quasi_static_env< typename Ndof_rl_space<double, NDOF, 1>::type, cubic_hermite_interpolation_tag> >; \
extern template class rrtstar_planner< \
  manip_quasi_static_env< typename Ndof_rl_space<double, NDOF, 2>::type, cubic_hermite_interpolation_tag> >; \
 \
extern template class rrtstar_planner< \
  manip_quasi_static_env< typename Ndof_rl_space<double, NDOF, 2>::type, quintic_hermite_interpolation_tag> >; \
 \
extern template class rrtstar_planner< \
  manip_quasi_static_env< typename Ndof_rl_space<double, NDOF, 1>::type, svp_Ndof_interpolation_tag> >; \
extern template class rrtstar_planner< \
  manip_quasi_static_env< typename Ndof_rl_space<double, NDOF, 2>::type, svp_Ndof_interpolation_tag> >; \
 \
extern template class rrtstar_planner< \
  manip_quasi_static_env< typename Ndof_rl_space<double, NDOF, 2>::type, sap_Ndof_interpolation_tag> >;


// RK_RRTSTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(1)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(2)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(3)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(4)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(5)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(6)
RK_RRTSTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(7)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(8)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(9)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_DECL(10)



#define RK_RRTSTAR_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(NDOF) \
extern template class rrtstar_planner< \
  manip_dynamic_env< typename Ndof_rl_space<double, NDOF, 0>::type, linear_interpolation_tag> >; \
extern template class rrtstar_planner< \
  manip_dynamic_env< typename Ndof_rl_space<double, NDOF, 1>::type, linear_interpolation_tag> >; \
extern template class rrtstar_planner< \
  manip_dynamic_env< typename Ndof_rl_space<double, NDOF, 2>::type, linear_interpolation_tag> >; \
 \
extern template class rrtstar_planner< \
  manip_dynamic_env< typename Ndof_rl_space<double, NDOF, 1>::type, cubic_hermite_interpolation_tag> >; \
extern template class rrtstar_planner< \
  manip_dynamic_env< typename Ndof_rl_space<double, NDOF, 2>::type, cubic_hermite_interpolation_tag> >; \
 \
extern template class rrtstar_planner< \
  manip_dynamic_env< typename Ndof_rl_space<double, NDOF, 2>::type, quintic_hermite_interpolation_tag> >; \
 \
extern template class rrtstar_planner< \
  manip_dynamic_env< typename Ndof_rl_space<double, NDOF, 1>::type, svp_Ndof_interpolation_tag> >; \
extern template class rrtstar_planner< \
  manip_dynamic_env< typename Ndof_rl_space<double, NDOF, 2>::type, svp_Ndof_interpolation_tag> >; \
 \
extern template class rrtstar_planner< \
  manip_dynamic_env< typename Ndof_rl_space<double, NDOF, 2>::type, sap_Ndof_interpolation_tag> >;


// RK_RRTSTAR_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(1)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(2)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(3)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(4)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(5)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(6)
RK_RRTSTAR_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(7)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(8)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(9)
// RK_RRTSTAR_MANIP_PLANNERS_MAKE_DYNAMIC_MANIP_EXTERN_DECL(10)



};

};

#else

#include "path_planning/rrtstar_path_planner.tpp"

#endif


#endif


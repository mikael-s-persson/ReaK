
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


#include "base/defs.hpp"

#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

#include "path_planning/sbastar_manip_planners.hpp"

#include "path_planning/sbastar_path_planner.tpp"

namespace ReaK {

namespace pp {


#define RK_SBASTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_INSTANCES(NDOF) \
template class sbastar_planner< \
  manip_quasi_static_env< typename Ndof_rl_space<double, NDOF, 0>::type > >; \
template class sbastar_planner< \
  manip_quasi_static_env< typename Ndof_rl_space<double, NDOF, 1>::type > >; \
template class sbastar_planner< \
  manip_quasi_static_env< typename Ndof_rl_space<double, NDOF, 2>::type > >;

RK_SBASTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_INSTANCES(RK_PLANNERS_DIM)
  
//RK_SBASTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_INSTANCES(0)
// RK_SBASTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_INSTANCES(1)
// RK_SBASTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_INSTANCES(2)
// RK_SBASTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_INSTANCES(3)
// RK_SBASTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_INSTANCES(4)
// RK_SBASTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_INSTANCES(5)
// RK_SBASTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_INSTANCES(6)
// RK_SBASTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_INSTANCES(7)
// RK_SBASTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_INSTANCES(8)
// RK_SBASTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_INSTANCES(9)
// RK_SBASTAR_MANIP_PLANNERS_MAKE_STATIC_MANIP_EXTERN_INSTANCES(10)

};

};

#else

namespace ReaK {

namespace pp {

#define CONCAT_NAMES(X,Y) X##Y

void CONCAT_NAMES(dummy_sbastar_manip_planners_static_externs_symbol_, RK_PLANNERS_DIM)() { };

};

};

#endif



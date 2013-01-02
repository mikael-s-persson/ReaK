/**
 * \file manip_free_dynamic_workspace_ext.hpp
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



#ifndef REAK_MANIP_FREE_DYNAMIC_WORKSPACE_EXT_HPP
#define REAK_MANIP_FREE_DYNAMIC_WORKSPACE_EXT_HPP

#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

#include "base/defs.hpp"
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "differentiable_space.hpp"
#include "metric_space_tuple.hpp"
#include "rate_limited_spaces.hpp"

#include "joint_space_topologies.hpp"
#include "se2_topologies.hpp"
#include "se3_topologies.hpp"

#include "interpolation/linear_interp.hpp"
#include "interpolation/cubic_hermite_interp.hpp"
#include "interpolation/quintic_hermite_interp.hpp"
#include "interpolation/sustained_velocity_pulse.hpp"
#include "interpolation/sustained_acceleration_pulse.hpp"

namespace ReaK {

namespace pp {


#define RK_MANIP_FREE_WORKSPACE_MAKE_DYN_ENV_FOR_JOINTS(NDOF) \
extern template class manip_dynamic_env< Ndof_0th_order_rl_space<double, NDOF, euclidean_tuple_distance>::type, linear_interpolation_tag>;\
extern template class manip_dynamic_env< Ndof_1st_order_rl_space<double, NDOF, euclidean_tuple_distance>::type, linear_interpolation_tag>;\
extern template class manip_dynamic_env< Ndof_2nd_order_rl_space<double, NDOF, euclidean_tuple_distance>::type, linear_interpolation_tag>;\
\
extern template class manip_dynamic_env< Ndof_1st_order_rl_space<double, NDOF, euclidean_tuple_distance>::type, cubic_hermite_interpolation_tag>;\
extern template class manip_dynamic_env< Ndof_2nd_order_rl_space<double, NDOF, euclidean_tuple_distance>::type, cubic_hermite_interpolation_tag>;\
\
extern template class manip_dynamic_env< Ndof_2nd_order_rl_space<double, NDOF, euclidean_tuple_distance>::type, quintic_hermite_interpolation_tag>;\
\
extern template class manip_dynamic_env< Ndof_1st_order_rl_space<double, NDOF, euclidean_tuple_distance>::type, svp_interpolation_tag>;\
extern template class manip_dynamic_env< Ndof_2nd_order_rl_space<double, NDOF, euclidean_tuple_distance>::type, svp_interpolation_tag>;\
\
extern template class manip_dynamic_env< Ndof_2nd_order_rl_space<double, NDOF, euclidean_tuple_distance>::type, sap_interpolation_tag>;\
\
\
extern template class manip_dynamic_env< Ndof_0th_order_rl_space<double, NDOF, inf_norm_tuple_distance>::type, linear_interpolation_tag>;\
extern template class manip_dynamic_env< Ndof_1st_order_rl_space<double, NDOF, inf_norm_tuple_distance>::type, linear_interpolation_tag>;\
extern template class manip_dynamic_env< Ndof_2nd_order_rl_space<double, NDOF, inf_norm_tuple_distance>::type, linear_interpolation_tag>;\
\
extern template class manip_dynamic_env< Ndof_1st_order_rl_space<double, NDOF, inf_norm_tuple_distance>::type, cubic_hermite_interpolation_tag>;\
extern template class manip_dynamic_env< Ndof_2nd_order_rl_space<double, NDOF, inf_norm_tuple_distance>::type, cubic_hermite_interpolation_tag>;\
\
extern template class manip_dynamic_env< Ndof_2nd_order_rl_space<double, NDOF, inf_norm_tuple_distance>::type, quintic_hermite_interpolation_tag>;\
\
extern template class manip_dynamic_env< Ndof_1st_order_rl_space<double, NDOF, inf_norm_tuple_distance>::type, svp_interpolation_tag>;\
extern template class manip_dynamic_env< Ndof_2nd_order_rl_space<double, NDOF, inf_norm_tuple_distance>::type, svp_interpolation_tag>;\
\
extern template class manip_dynamic_env< Ndof_2nd_order_rl_space<double, NDOF, inf_norm_tuple_distance>::type, sap_interpolation_tag>;

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


extern template class manip_dynamic_env< metric_space_array< se2_0th_order_rl_topology<double>::type, 1>::type, linear_interpolation_tag>;
extern template class manip_dynamic_env< metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type, linear_interpolation_tag>;
extern template class manip_dynamic_env< metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type, linear_interpolation_tag>;

extern template class manip_dynamic_env< metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type, cubic_hermite_interpolation_tag>;
extern template class manip_dynamic_env< metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type, cubic_hermite_interpolation_tag>;

extern template class manip_dynamic_env< metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type, quintic_hermite_interpolation_tag>;

extern template class manip_dynamic_env< metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type, svp_interpolation_tag>;
extern template class manip_dynamic_env< metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type, svp_interpolation_tag>;

extern template class manip_dynamic_env< metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type, sap_interpolation_tag>;


extern template class manip_dynamic_env< metric_space_array< se3_0th_order_rl_topology<double>::type, 1>::type, linear_interpolation_tag>;
extern template class manip_dynamic_env< metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type, linear_interpolation_tag>;
extern template class manip_dynamic_env< metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type, linear_interpolation_tag>;

extern template class manip_dynamic_env< metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type, cubic_hermite_interpolation_tag>;
extern template class manip_dynamic_env< metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type, cubic_hermite_interpolation_tag>;

extern template class manip_dynamic_env< metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type, quintic_hermite_interpolation_tag>;

extern template class manip_dynamic_env< metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type, svp_interpolation_tag>;
extern template class manip_dynamic_env< metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type, svp_interpolation_tag>;

extern template class manip_dynamic_env< metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type, sap_interpolation_tag>;


};

};





#endif

#endif



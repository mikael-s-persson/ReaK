/**
 * \file manip_workspaces.hpp
 *
 * This library defines a class
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
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

#ifndef REAK_MANIP_WORKSPACES_HPP
#define REAK_MANIP_WORKSPACES_HPP

#include "ReaK/core/base/defs.hpp"
#include "ReaK/core/base/named_object.hpp"

#include "ReaK/topologies/spaces/joint_space_limits.hpp"
#include "ReaK/topologies/spaces/joint_space_topologies.hpp"
#include "ReaK/topologies/spaces/se3_topologies.hpp"

#include "ReaK/topologies/interpolation/cubic_hermite_interp.hpp"
#include "ReaK/topologies/interpolation/linear_interp.hpp"
#include "ReaK/topologies/interpolation/quintic_hermite_interp.hpp"
#include "ReaK/topologies/interpolation/sustained_acceleration_pulse.hpp"
#include "ReaK/topologies/interpolation/sustained_velocity_pulse.hpp"

#include "ReaK/topologies/spaces/manip_free_dynamic_workspace.hpp"
#include "ReaK/topologies/spaces/manip_free_workspace.hpp"

namespace ReaK::pp {

// Rate-limited N-dof Joint-space of order 0, with infinity-norm distance metric.
using rl_jspace_1d_o0_type =
    Ndof_0th_order_rl_space_t<double, 1, inf_norm_tuple_distance>;
using rl_jspace_2d_o0_type =
    Ndof_0th_order_rl_space_t<double, 2, inf_norm_tuple_distance>;
using rl_jspace_3d_o0_type =
    Ndof_0th_order_rl_space_t<double, 3, inf_norm_tuple_distance>;
using rl_jspace_4d_o0_type =
    Ndof_0th_order_rl_space_t<double, 4, inf_norm_tuple_distance>;
using rl_jspace_5d_o0_type =
    Ndof_0th_order_rl_space_t<double, 5, inf_norm_tuple_distance>;
using rl_jspace_6d_o0_type =
    Ndof_0th_order_rl_space_t<double, 6, inf_norm_tuple_distance>;
using rl_jspace_7d_o0_type =
    Ndof_0th_order_rl_space_t<double, 7, inf_norm_tuple_distance>;
using rl_jspace_8d_o0_type =
    Ndof_0th_order_rl_space_t<double, 8, inf_norm_tuple_distance>;
using rl_jspace_9d_o0_type =
    Ndof_0th_order_rl_space_t<double, 9, inf_norm_tuple_distance>;

// Rate-limited N-dof Joint-space of order 1, with infinity-norm distance metric.
using rl_jspace_1d_o1_type =
    Ndof_1st_order_rl_space_t<double, 1, inf_norm_tuple_distance>;
using rl_jspace_2d_o1_type =
    Ndof_1st_order_rl_space_t<double, 2, inf_norm_tuple_distance>;
using rl_jspace_3d_o1_type =
    Ndof_1st_order_rl_space_t<double, 3, inf_norm_tuple_distance>;
using rl_jspace_4d_o1_type =
    Ndof_1st_order_rl_space_t<double, 4, inf_norm_tuple_distance>;
using rl_jspace_5d_o1_type =
    Ndof_1st_order_rl_space_t<double, 5, inf_norm_tuple_distance>;
using rl_jspace_6d_o1_type =
    Ndof_1st_order_rl_space_t<double, 6, inf_norm_tuple_distance>;
using rl_jspace_7d_o1_type =
    Ndof_1st_order_rl_space_t<double, 7, inf_norm_tuple_distance>;
using rl_jspace_8d_o1_type =
    Ndof_1st_order_rl_space_t<double, 8, inf_norm_tuple_distance>;
using rl_jspace_9d_o1_type =
    Ndof_1st_order_rl_space_t<double, 9, inf_norm_tuple_distance>;

// Rate-limited N-dof Joint-space of order 2, with infinity-norm distance metric.
using rl_jspace_1d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 1, inf_norm_tuple_distance>;
using rl_jspace_2d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 2, inf_norm_tuple_distance>;
using rl_jspace_3d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 3, inf_norm_tuple_distance>;
using rl_jspace_4d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 4, inf_norm_tuple_distance>;
using rl_jspace_5d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 5, inf_norm_tuple_distance>;
using rl_jspace_6d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 6, inf_norm_tuple_distance>;
using rl_jspace_7d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 7, inf_norm_tuple_distance>;
using rl_jspace_8d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 8, inf_norm_tuple_distance>;
using rl_jspace_9d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 9, inf_norm_tuple_distance>;

// Rate-limited N-dof Joint-space of order 0, with euclidean distance metric.
using rle_jspace_1d_o0_type =
    Ndof_0th_order_rl_space_t<double, 1, euclidean_tuple_distance>;
using rle_jspace_2d_o0_type =
    Ndof_0th_order_rl_space_t<double, 2, euclidean_tuple_distance>;
using rle_jspace_3d_o0_type =
    Ndof_0th_order_rl_space_t<double, 3, euclidean_tuple_distance>;
using rle_jspace_4d_o0_type =
    Ndof_0th_order_rl_space_t<double, 4, euclidean_tuple_distance>;
using rle_jspace_5d_o0_type =
    Ndof_0th_order_rl_space_t<double, 5, euclidean_tuple_distance>;
using rle_jspace_6d_o0_type =
    Ndof_0th_order_rl_space_t<double, 6, euclidean_tuple_distance>;
using rle_jspace_7d_o0_type =
    Ndof_0th_order_rl_space_t<double, 7, euclidean_tuple_distance>;
using rle_jspace_8d_o0_type =
    Ndof_0th_order_rl_space_t<double, 8, euclidean_tuple_distance>;
using rle_jspace_9d_o0_type =
    Ndof_0th_order_rl_space_t<double, 9, euclidean_tuple_distance>;

// Rate-limited N-dof Joint-space of order 1, with euclidean distance metric.
using rle_jspace_1d_o1_type =
    Ndof_1st_order_rl_space_t<double, 1, euclidean_tuple_distance>;
using rle_jspace_2d_o1_type =
    Ndof_1st_order_rl_space_t<double, 2, euclidean_tuple_distance>;
using rle_jspace_3d_o1_type =
    Ndof_1st_order_rl_space_t<double, 3, euclidean_tuple_distance>;
using rle_jspace_4d_o1_type =
    Ndof_1st_order_rl_space_t<double, 4, euclidean_tuple_distance>;
using rle_jspace_5d_o1_type =
    Ndof_1st_order_rl_space_t<double, 5, euclidean_tuple_distance>;
using rle_jspace_6d_o1_type =
    Ndof_1st_order_rl_space_t<double, 6, euclidean_tuple_distance>;
using rle_jspace_7d_o1_type =
    Ndof_1st_order_rl_space_t<double, 7, euclidean_tuple_distance>;
using rle_jspace_8d_o1_type =
    Ndof_1st_order_rl_space_t<double, 8, euclidean_tuple_distance>;
using rle_jspace_9d_o1_type =
    Ndof_1st_order_rl_space_t<double, 9, euclidean_tuple_distance>;

// Rate-limited N-dof Joint-space of order 2, with euclidean distance metric.
using rle_jspace_1d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 1, euclidean_tuple_distance>;
using rle_jspace_2d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 2, euclidean_tuple_distance>;
using rle_jspace_3d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 3, euclidean_tuple_distance>;
using rle_jspace_4d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 4, euclidean_tuple_distance>;
using rle_jspace_5d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 5, euclidean_tuple_distance>;
using rle_jspace_6d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 6, euclidean_tuple_distance>;
using rle_jspace_7d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 7, euclidean_tuple_distance>;
using rle_jspace_8d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 8, euclidean_tuple_distance>;
using rle_jspace_9d_o2_type =
    Ndof_2nd_order_rl_space_t<double, 9, euclidean_tuple_distance>;

// N-dof Joint-space of order 0, with euclidean distance metric.
using jspace_1d_o0_type =
    Ndof_0th_order_space_t<double, 1, euclidean_tuple_distance>;
using jspace_2d_o0_type =
    Ndof_0th_order_space_t<double, 2, euclidean_tuple_distance>;
using jspace_3d_o0_type =
    Ndof_0th_order_space_t<double, 3, euclidean_tuple_distance>;
using jspace_4d_o0_type =
    Ndof_0th_order_space_t<double, 4, euclidean_tuple_distance>;
using jspace_5d_o0_type =
    Ndof_0th_order_space_t<double, 5, euclidean_tuple_distance>;
using jspace_6d_o0_type =
    Ndof_0th_order_space_t<double, 6, euclidean_tuple_distance>;
using jspace_7d_o0_type =
    Ndof_0th_order_space_t<double, 7, euclidean_tuple_distance>;
using jspace_8d_o0_type =
    Ndof_0th_order_space_t<double, 8, euclidean_tuple_distance>;
using jspace_9d_o0_type =
    Ndof_0th_order_space_t<double, 9, euclidean_tuple_distance>;

// N-dof Joint-space of order 1, with euclidean distance metric.
using jspace_1d_o1_type =
    Ndof_1st_order_space_t<double, 1, euclidean_tuple_distance>;
using jspace_2d_o1_type =
    Ndof_1st_order_space_t<double, 2, euclidean_tuple_distance>;
using jspace_3d_o1_type =
    Ndof_1st_order_space_t<double, 3, euclidean_tuple_distance>;
using jspace_4d_o1_type =
    Ndof_1st_order_space_t<double, 4, euclidean_tuple_distance>;
using jspace_5d_o1_type =
    Ndof_1st_order_space_t<double, 5, euclidean_tuple_distance>;
using jspace_6d_o1_type =
    Ndof_1st_order_space_t<double, 6, euclidean_tuple_distance>;
using jspace_7d_o1_type =
    Ndof_1st_order_space_t<double, 7, euclidean_tuple_distance>;
using jspace_8d_o1_type =
    Ndof_1st_order_space_t<double, 8, euclidean_tuple_distance>;
using jspace_9d_o1_type =
    Ndof_1st_order_space_t<double, 9, euclidean_tuple_distance>;

// N-dof Joint-space of order 2, with euclidean distance metric.
using jspace_1d_o2_type =
    Ndof_2nd_order_space_t<double, 1, euclidean_tuple_distance>;
using jspace_2d_o2_type =
    Ndof_2nd_order_space_t<double, 2, euclidean_tuple_distance>;
using jspace_3d_o2_type =
    Ndof_2nd_order_space_t<double, 3, euclidean_tuple_distance>;
using jspace_4d_o2_type =
    Ndof_2nd_order_space_t<double, 4, euclidean_tuple_distance>;
using jspace_5d_o2_type =
    Ndof_2nd_order_space_t<double, 5, euclidean_tuple_distance>;
using jspace_6d_o2_type =
    Ndof_2nd_order_space_t<double, 6, euclidean_tuple_distance>;
using jspace_7d_o2_type =
    Ndof_2nd_order_space_t<double, 7, euclidean_tuple_distance>;
using jspace_8d_o2_type =
    Ndof_2nd_order_space_t<double, 8, euclidean_tuple_distance>;
using jspace_9d_o2_type =
    Ndof_2nd_order_space_t<double, 9, euclidean_tuple_distance>;

// 2D End-effector-space (SE(2)) of order N.
using eespace_2D_o0_type = se2_0th_order_topology_t<double>;
using eespace_2D_o1_type = se2_1st_order_topology_t<double>;
using eespace_2D_o2_type = se2_2nd_order_topology_t<double>;

// Rate-limited 2D End-effector-space (SE(2)) of order N.
using rl_eespace_2D_o0_type = se2_0th_order_rl_topology_t<double>;
using rl_eespace_2D_o1_type = se2_1st_order_rl_topology_t<double>;
using rl_eespace_2D_o2_type = se2_2nd_order_rl_topology_t<double>;

// 3D End-effector-space (SE(3)) of order N.
using eespace_3D_o0_type = se3_0th_order_topology_t<double>;
using eespace_3D_o1_type = se3_1st_order_topology_t<double>;
using eespace_3D_o2_type = se3_2nd_order_topology_t<double>;

// Rate-limited 3D End-effector-space (SE(3)) of order N.
using rl_eespace_3D_o0_type = se3_0th_order_rl_topology_t<double>;
using rl_eespace_3D_o1_type = se3_1st_order_rl_topology_t<double>;
using rl_eespace_3D_o2_type = se3_2nd_order_rl_topology_t<double>;

// Rate-limited N-dof static-workspace of order 0, with linear interpolation.
using rl_wspace_1d_o0_i1_type =
    manip_quasi_static_env<rl_jspace_1d_o0_type, linear_interpolation_tag>;
using rl_wspace_2d_o0_i1_type =
    manip_quasi_static_env<rl_jspace_2d_o0_type, linear_interpolation_tag>;
using rl_wspace_3d_o0_i1_type =
    manip_quasi_static_env<rl_jspace_3d_o0_type, linear_interpolation_tag>;
using rl_wspace_4d_o0_i1_type =
    manip_quasi_static_env<rl_jspace_4d_o0_type, linear_interpolation_tag>;
using rl_wspace_5d_o0_i1_type =
    manip_quasi_static_env<rl_jspace_5d_o0_type, linear_interpolation_tag>;
using rl_wspace_6d_o0_i1_type =
    manip_quasi_static_env<rl_jspace_6d_o0_type, linear_interpolation_tag>;
using rl_wspace_7d_o0_i1_type =
    manip_quasi_static_env<rl_jspace_7d_o0_type, linear_interpolation_tag>;
using rl_wspace_8d_o0_i1_type =
    manip_quasi_static_env<rl_jspace_8d_o0_type, linear_interpolation_tag>;
using rl_wspace_9d_o0_i1_type =
    manip_quasi_static_env<rl_jspace_9d_o0_type, linear_interpolation_tag>;

// Rate-limited N-dof static-workspace of order 1, with linear interpolation.
using rl_wspace_1d_o1_i1_type =
    manip_quasi_static_env<rl_jspace_1d_o1_type, linear_interpolation_tag>;
using rl_wspace_2d_o1_i1_type =
    manip_quasi_static_env<rl_jspace_2d_o1_type, linear_interpolation_tag>;
using rl_wspace_3d_o1_i1_type =
    manip_quasi_static_env<rl_jspace_3d_o1_type, linear_interpolation_tag>;
using rl_wspace_4d_o1_i1_type =
    manip_quasi_static_env<rl_jspace_4d_o1_type, linear_interpolation_tag>;
using rl_wspace_5d_o1_i1_type =
    manip_quasi_static_env<rl_jspace_5d_o1_type, linear_interpolation_tag>;
using rl_wspace_6d_o1_i1_type =
    manip_quasi_static_env<rl_jspace_6d_o1_type, linear_interpolation_tag>;
using rl_wspace_7d_o1_i1_type =
    manip_quasi_static_env<rl_jspace_7d_o1_type, linear_interpolation_tag>;
using rl_wspace_8d_o1_i1_type =
    manip_quasi_static_env<rl_jspace_8d_o1_type, linear_interpolation_tag>;
using rl_wspace_9d_o1_i1_type =
    manip_quasi_static_env<rl_jspace_9d_o1_type, linear_interpolation_tag>;

// Rate-limited N-dof static-workspace of order 2, with linear interpolation.
using rl_wspace_1d_o2_i1_type =
    manip_quasi_static_env<rl_jspace_1d_o2_type, linear_interpolation_tag>;
using rl_wspace_2d_o2_i1_type =
    manip_quasi_static_env<rl_jspace_2d_o2_type, linear_interpolation_tag>;
using rl_wspace_3d_o2_i1_type =
    manip_quasi_static_env<rl_jspace_3d_o2_type, linear_interpolation_tag>;
using rl_wspace_4d_o2_i1_type =
    manip_quasi_static_env<rl_jspace_4d_o2_type, linear_interpolation_tag>;
using rl_wspace_5d_o2_i1_type =
    manip_quasi_static_env<rl_jspace_5d_o2_type, linear_interpolation_tag>;
using rl_wspace_6d_o2_i1_type =
    manip_quasi_static_env<rl_jspace_6d_o2_type, linear_interpolation_tag>;
using rl_wspace_7d_o2_i1_type =
    manip_quasi_static_env<rl_jspace_7d_o2_type, linear_interpolation_tag>;
using rl_wspace_8d_o2_i1_type =
    manip_quasi_static_env<rl_jspace_8d_o2_type, linear_interpolation_tag>;
using rl_wspace_9d_o2_i1_type =
    manip_quasi_static_env<rl_jspace_9d_o2_type, linear_interpolation_tag>;

// Rate-limited N-dof static-workspace of order 1, with cubic interpolation.
using rl_wspace_1d_o1_i3_type =
    manip_quasi_static_env<rl_jspace_1d_o1_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_2d_o1_i3_type =
    manip_quasi_static_env<rl_jspace_2d_o1_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_3d_o1_i3_type =
    manip_quasi_static_env<rl_jspace_3d_o1_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_4d_o1_i3_type =
    manip_quasi_static_env<rl_jspace_4d_o1_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_5d_o1_i3_type =
    manip_quasi_static_env<rl_jspace_5d_o1_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_6d_o1_i3_type =
    manip_quasi_static_env<rl_jspace_6d_o1_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_7d_o1_i3_type =
    manip_quasi_static_env<rl_jspace_7d_o1_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_8d_o1_i3_type =
    manip_quasi_static_env<rl_jspace_8d_o1_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_9d_o1_i3_type =
    manip_quasi_static_env<rl_jspace_9d_o1_type,
                           cubic_hermite_interpolation_tag>;

// Rate-limited N-dof static-workspace of order 2, with cubic interpolation.
using rl_wspace_1d_o2_i3_type =
    manip_quasi_static_env<rl_jspace_1d_o2_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_2d_o2_i3_type =
    manip_quasi_static_env<rl_jspace_2d_o2_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_3d_o2_i3_type =
    manip_quasi_static_env<rl_jspace_3d_o2_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_4d_o2_i3_type =
    manip_quasi_static_env<rl_jspace_4d_o2_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_5d_o2_i3_type =
    manip_quasi_static_env<rl_jspace_5d_o2_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_6d_o2_i3_type =
    manip_quasi_static_env<rl_jspace_6d_o2_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_7d_o2_i3_type =
    manip_quasi_static_env<rl_jspace_7d_o2_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_8d_o2_i3_type =
    manip_quasi_static_env<rl_jspace_8d_o2_type,
                           cubic_hermite_interpolation_tag>;
using rl_wspace_9d_o2_i3_type =
    manip_quasi_static_env<rl_jspace_9d_o2_type,
                           cubic_hermite_interpolation_tag>;

///< Rate-limited N-dof static-workspace of order 2, with quintic interpolation.
using rl_wspace_1d_o2_i5_type =
    manip_quasi_static_env<rl_jspace_1d_o2_type,
                           quintic_hermite_interpolation_tag>;
using rl_wspace_2d_o2_i5_type =
    manip_quasi_static_env<rl_jspace_2d_o2_type,
                           quintic_hermite_interpolation_tag>;
using rl_wspace_3d_o2_i5_type =
    manip_quasi_static_env<rl_jspace_3d_o2_type,
                           quintic_hermite_interpolation_tag>;
using rl_wspace_4d_o2_i5_type =
    manip_quasi_static_env<rl_jspace_4d_o2_type,
                           quintic_hermite_interpolation_tag>;
using rl_wspace_5d_o2_i5_type =
    manip_quasi_static_env<rl_jspace_5d_o2_type,
                           quintic_hermite_interpolation_tag>;
using rl_wspace_6d_o2_i5_type =
    manip_quasi_static_env<rl_jspace_6d_o2_type,
                           quintic_hermite_interpolation_tag>;
using rl_wspace_7d_o2_i5_type =
    manip_quasi_static_env<rl_jspace_7d_o2_type,
                           quintic_hermite_interpolation_tag>;
using rl_wspace_8d_o2_i5_type =
    manip_quasi_static_env<rl_jspace_8d_o2_type,
                           quintic_hermite_interpolation_tag>;
using rl_wspace_9d_o2_i5_type =
    manip_quasi_static_env<rl_jspace_9d_o2_type,
                           quintic_hermite_interpolation_tag>;

// Rate-limited N-dof static-workspace of order 1, with SVP interpolation.
using rl_wspace_1d_o1_svp_type =
    manip_quasi_static_env<rl_jspace_1d_o1_type, svp_interpolation_tag>;
using rl_wspace_2d_o1_svp_type =
    manip_quasi_static_env<rl_jspace_2d_o1_type, svp_interpolation_tag>;
using rl_wspace_3d_o1_svp_type =
    manip_quasi_static_env<rl_jspace_3d_o1_type, svp_interpolation_tag>;
using rl_wspace_4d_o1_svp_type =
    manip_quasi_static_env<rl_jspace_4d_o1_type, svp_interpolation_tag>;
using rl_wspace_5d_o1_svp_type =
    manip_quasi_static_env<rl_jspace_5d_o1_type, svp_interpolation_tag>;
using rl_wspace_6d_o1_svp_type =
    manip_quasi_static_env<rl_jspace_6d_o1_type, svp_interpolation_tag>;
using rl_wspace_7d_o1_svp_type =
    manip_quasi_static_env<rl_jspace_7d_o1_type, svp_interpolation_tag>;
using rl_wspace_8d_o1_svp_type =
    manip_quasi_static_env<rl_jspace_8d_o1_type, svp_interpolation_tag>;
using rl_wspace_9d_o1_svp_type =
    manip_quasi_static_env<rl_jspace_9d_o1_type, svp_interpolation_tag>;

// Rate-limited N-dof static-workspace of order 2, with SVP interpolation.
using rl_wspace_1d_o2_svp_type =
    manip_quasi_static_env<rl_jspace_1d_o2_type, svp_interpolation_tag>;
using rl_wspace_2d_o2_svp_type =
    manip_quasi_static_env<rl_jspace_2d_o2_type, svp_interpolation_tag>;
using rl_wspace_3d_o2_svp_type =
    manip_quasi_static_env<rl_jspace_3d_o2_type, svp_interpolation_tag>;
using rl_wspace_4d_o2_svp_type =
    manip_quasi_static_env<rl_jspace_4d_o2_type, svp_interpolation_tag>;
using rl_wspace_5d_o2_svp_type =
    manip_quasi_static_env<rl_jspace_5d_o2_type, svp_interpolation_tag>;
using rl_wspace_6d_o2_svp_type =
    manip_quasi_static_env<rl_jspace_6d_o2_type, svp_interpolation_tag>;
using rl_wspace_7d_o2_svp_type =
    manip_quasi_static_env<rl_jspace_7d_o2_type, svp_interpolation_tag>;
using rl_wspace_8d_o2_svp_type =
    manip_quasi_static_env<rl_jspace_8d_o2_type, svp_interpolation_tag>;
using rl_wspace_9d_o2_svp_type =
    manip_quasi_static_env<rl_jspace_9d_o2_type, svp_interpolation_tag>;

// Rate-limited N-dof static-workspace of order 2, with SAP interpolation.
using rl_wspace_1d_o2_sap_type =
    manip_quasi_static_env<rl_jspace_1d_o2_type, sap_interpolation_tag>;
using rl_wspace_2d_o2_sap_type =
    manip_quasi_static_env<rl_jspace_2d_o2_type, sap_interpolation_tag>;
using rl_wspace_3d_o2_sap_type =
    manip_quasi_static_env<rl_jspace_3d_o2_type, sap_interpolation_tag>;
using rl_wspace_4d_o2_sap_type =
    manip_quasi_static_env<rl_jspace_4d_o2_type, sap_interpolation_tag>;
using rl_wspace_5d_o2_sap_type =
    manip_quasi_static_env<rl_jspace_5d_o2_type, sap_interpolation_tag>;
using rl_wspace_6d_o2_sap_type =
    manip_quasi_static_env<rl_jspace_6d_o2_type, sap_interpolation_tag>;
using rl_wspace_7d_o2_sap_type =
    manip_quasi_static_env<rl_jspace_7d_o2_type, sap_interpolation_tag>;
using rl_wspace_8d_o2_sap_type =
    manip_quasi_static_env<rl_jspace_8d_o2_type, sap_interpolation_tag>;
using rl_wspace_9d_o2_sap_type =
    manip_quasi_static_env<rl_jspace_9d_o2_type, sap_interpolation_tag>;

// Rate-limited N-dof dynamic-workspace of order 0, with linear interpolation.
using rl_wspace_1d_o0_i1_type =
    manip_dynamic_env<rl_jspace_1d_o0_type, linear_interpolation_tag>;
using rl_wspace_2d_o0_i1_type =
    manip_dynamic_env<rl_jspace_2d_o0_type, linear_interpolation_tag>;
using rl_wspace_3d_o0_i1_type =
    manip_dynamic_env<rl_jspace_3d_o0_type, linear_interpolation_tag>;
using rl_wspace_4d_o0_i1_type =
    manip_dynamic_env<rl_jspace_4d_o0_type, linear_interpolation_tag>;
using rl_wspace_5d_o0_i1_type =
    manip_dynamic_env<rl_jspace_5d_o0_type, linear_interpolation_tag>;
using rl_wspace_6d_o0_i1_type =
    manip_dynamic_env<rl_jspace_6d_o0_type, linear_interpolation_tag>;
using rl_wspace_7d_o0_i1_type =
    manip_dynamic_env<rl_jspace_7d_o0_type, linear_interpolation_tag>;
using rl_wspace_8d_o0_i1_type =
    manip_dynamic_env<rl_jspace_8d_o0_type, linear_interpolation_tag>;
using rl_wspace_9d_o0_i1_type =
    manip_dynamic_env<rl_jspace_9d_o0_type, linear_interpolation_tag>;

// Rate-limited N-dof dynamic-workspace of order 1, with linear interpolation.
using rl_wspace_1d_o1_i1_type =
    manip_dynamic_env<rl_jspace_1d_o1_type, linear_interpolation_tag>;
using rl_wspace_2d_o1_i1_type =
    manip_dynamic_env<rl_jspace_2d_o1_type, linear_interpolation_tag>;
using rl_wspace_3d_o1_i1_type =
    manip_dynamic_env<rl_jspace_3d_o1_type, linear_interpolation_tag>;
using rl_wspace_4d_o1_i1_type =
    manip_dynamic_env<rl_jspace_4d_o1_type, linear_interpolation_tag>;
using rl_wspace_5d_o1_i1_type =
    manip_dynamic_env<rl_jspace_5d_o1_type, linear_interpolation_tag>;
using rl_wspace_6d_o1_i1_type =
    manip_dynamic_env<rl_jspace_6d_o1_type, linear_interpolation_tag>;
using rl_wspace_7d_o1_i1_type =
    manip_dynamic_env<rl_jspace_7d_o1_type, linear_interpolation_tag>;
using rl_wspace_8d_o1_i1_type =
    manip_dynamic_env<rl_jspace_8d_o1_type, linear_interpolation_tag>;
using rl_wspace_9d_o1_i1_type =
    manip_dynamic_env<rl_jspace_9d_o1_type, linear_interpolation_tag>;

// Rate-limited N-dof dynamic-workspace of order 2, with linear interpolation.
using rl_wspace_1d_o2_i1_type =
    manip_dynamic_env<rl_jspace_1d_o2_type, linear_interpolation_tag>;
using rl_wspace_2d_o2_i1_type =
    manip_dynamic_env<rl_jspace_2d_o2_type, linear_interpolation_tag>;
using rl_wspace_3d_o2_i1_type =
    manip_dynamic_env<rl_jspace_3d_o2_type, linear_interpolation_tag>;
using rl_wspace_4d_o2_i1_type =
    manip_dynamic_env<rl_jspace_4d_o2_type, linear_interpolation_tag>;
using rl_wspace_5d_o2_i1_type =
    manip_dynamic_env<rl_jspace_5d_o2_type, linear_interpolation_tag>;
using rl_wspace_6d_o2_i1_type =
    manip_dynamic_env<rl_jspace_6d_o2_type, linear_interpolation_tag>;
using rl_wspace_7d_o2_i1_type =
    manip_dynamic_env<rl_jspace_7d_o2_type, linear_interpolation_tag>;
using rl_wspace_8d_o2_i1_type =
    manip_dynamic_env<rl_jspace_8d_o2_type, linear_interpolation_tag>;
using rl_wspace_9d_o2_i1_type =
    manip_dynamic_env<rl_jspace_9d_o2_type, linear_interpolation_tag>;

// Rate-limited N-dof dynamic-workspace of order 1, with cubic interpolation.
using rl_wspace_1d_o1_i3_type =
    manip_dynamic_env<rl_jspace_1d_o1_type, cubic_hermite_interpolation_tag>;
using rl_wspace_2d_o1_i3_type =
    manip_dynamic_env<rl_jspace_2d_o1_type, cubic_hermite_interpolation_tag>;
using rl_wspace_3d_o1_i3_type =
    manip_dynamic_env<rl_jspace_3d_o1_type, cubic_hermite_interpolation_tag>;
using rl_wspace_4d_o1_i3_type =
    manip_dynamic_env<rl_jspace_4d_o1_type, cubic_hermite_interpolation_tag>;
using rl_wspace_5d_o1_i3_type =
    manip_dynamic_env<rl_jspace_5d_o1_type, cubic_hermite_interpolation_tag>;
using rl_wspace_6d_o1_i3_type =
    manip_dynamic_env<rl_jspace_6d_o1_type, cubic_hermite_interpolation_tag>;
using rl_wspace_7d_o1_i3_type =
    manip_dynamic_env<rl_jspace_7d_o1_type, cubic_hermite_interpolation_tag>;
using rl_wspace_8d_o1_i3_type =
    manip_dynamic_env<rl_jspace_8d_o1_type, cubic_hermite_interpolation_tag>;
using rl_wspace_9d_o1_i3_type =
    manip_dynamic_env<rl_jspace_9d_o1_type, cubic_hermite_interpolation_tag>;

// Rate-limited N-dof dynamic-workspace of order 2, with cubic interpolation.
using rl_wspace_1d_o2_i3_type =
    manip_dynamic_env<rl_jspace_1d_o2_type, cubic_hermite_interpolation_tag>;
using rl_wspace_2d_o2_i3_type =
    manip_dynamic_env<rl_jspace_2d_o2_type, cubic_hermite_interpolation_tag>;
using rl_wspace_3d_o2_i3_type =
    manip_dynamic_env<rl_jspace_3d_o2_type, cubic_hermite_interpolation_tag>;
using rl_wspace_4d_o2_i3_type =
    manip_dynamic_env<rl_jspace_4d_o2_type, cubic_hermite_interpolation_tag>;
using rl_wspace_5d_o2_i3_type =
    manip_dynamic_env<rl_jspace_5d_o2_type, cubic_hermite_interpolation_tag>;
using rl_wspace_6d_o2_i3_type =
    manip_dynamic_env<rl_jspace_6d_o2_type, cubic_hermite_interpolation_tag>;
using rl_wspace_7d_o2_i3_type =
    manip_dynamic_env<rl_jspace_7d_o2_type, cubic_hermite_interpolation_tag>;
using rl_wspace_8d_o2_i3_type =
    manip_dynamic_env<rl_jspace_8d_o2_type, cubic_hermite_interpolation_tag>;
using rl_wspace_9d_o2_i3_type =
    manip_dynamic_env<rl_jspace_9d_o2_type, cubic_hermite_interpolation_tag>;

// Rate-limited N-dof dynamic-workspace of order 2, with quintic interpolation.
using rl_wspace_1d_o2_i5_type =
    manip_dynamic_env<rl_jspace_1d_o2_type, quintic_hermite_interpolation_tag>;
using rl_wspace_2d_o2_i5_type =
    manip_dynamic_env<rl_jspace_2d_o2_type, quintic_hermite_interpolation_tag>;
using rl_wspace_3d_o2_i5_type =
    manip_dynamic_env<rl_jspace_3d_o2_type, quintic_hermite_interpolation_tag>;
using rl_wspace_4d_o2_i5_type =
    manip_dynamic_env<rl_jspace_4d_o2_type, quintic_hermite_interpolation_tag>;
using rl_wspace_5d_o2_i5_type =
    manip_dynamic_env<rl_jspace_5d_o2_type, quintic_hermite_interpolation_tag>;
using rl_wspace_6d_o2_i5_type =
    manip_dynamic_env<rl_jspace_6d_o2_type, quintic_hermite_interpolation_tag>;
using rl_wspace_7d_o2_i5_type =
    manip_dynamic_env<rl_jspace_7d_o2_type, quintic_hermite_interpolation_tag>;
using rl_wspace_8d_o2_i5_type =
    manip_dynamic_env<rl_jspace_8d_o2_type, quintic_hermite_interpolation_tag>;
using rl_wspace_9d_o2_i5_type =
    manip_dynamic_env<rl_jspace_9d_o2_type, quintic_hermite_interpolation_tag>;

// Rate-limited N-dof dynamic-workspace of order 1, with SVP interpolation.
using rl_wspace_1d_o1_svp_type =
    manip_dynamic_env<rl_jspace_1d_o1_type, svp_interpolation_tag>;
using rl_wspace_2d_o1_svp_type =
    manip_dynamic_env<rl_jspace_2d_o1_type, svp_interpolation_tag>;
using rl_wspace_3d_o1_svp_type =
    manip_dynamic_env<rl_jspace_3d_o1_type, svp_interpolation_tag>;
using rl_wspace_4d_o1_svp_type =
    manip_dynamic_env<rl_jspace_4d_o1_type, svp_interpolation_tag>;
using rl_wspace_5d_o1_svp_type =
    manip_dynamic_env<rl_jspace_5d_o1_type, svp_interpolation_tag>;
using rl_wspace_6d_o1_svp_type =
    manip_dynamic_env<rl_jspace_6d_o1_type, svp_interpolation_tag>;
using rl_wspace_7d_o1_svp_type =
    manip_dynamic_env<rl_jspace_7d_o1_type, svp_interpolation_tag>;
using rl_wspace_8d_o1_svp_type =
    manip_dynamic_env<rl_jspace_8d_o1_type, svp_interpolation_tag>;
using rl_wspace_9d_o1_svp_type =
    manip_dynamic_env<rl_jspace_9d_o1_type, svp_interpolation_tag>;

// Rate-limited N-dof dynamic-workspace of order 2, with SVP interpolation.
using rl_wspace_1d_o2_svp_type =
    manip_dynamic_env<rl_jspace_1d_o2_type, svp_interpolation_tag>;
using rl_wspace_2d_o2_svp_type =
    manip_dynamic_env<rl_jspace_2d_o2_type, svp_interpolation_tag>;
using rl_wspace_3d_o2_svp_type =
    manip_dynamic_env<rl_jspace_3d_o2_type, svp_interpolation_tag>;
using rl_wspace_4d_o2_svp_type =
    manip_dynamic_env<rl_jspace_4d_o2_type, svp_interpolation_tag>;
using rl_wspace_5d_o2_svp_type =
    manip_dynamic_env<rl_jspace_5d_o2_type, svp_interpolation_tag>;
using rl_wspace_6d_o2_svp_type =
    manip_dynamic_env<rl_jspace_6d_o2_type, svp_interpolation_tag>;
using rl_wspace_7d_o2_svp_type =
    manip_dynamic_env<rl_jspace_7d_o2_type, svp_interpolation_tag>;
using rl_wspace_8d_o2_svp_type =
    manip_dynamic_env<rl_jspace_8d_o2_type, svp_interpolation_tag>;
using rl_wspace_9d_o2_svp_type =
    manip_dynamic_env<rl_jspace_9d_o2_type, svp_interpolation_tag>;

// Rate-limited N-dof dynamic-workspace of order 2, with SAP interpolation.
using rl_wspace_1d_o2_sap_type =
    manip_dynamic_env<rl_jspace_1d_o2_type, sap_interpolation_tag>;
using rl_wspace_2d_o2_sap_type =
    manip_dynamic_env<rl_jspace_2d_o2_type, sap_interpolation_tag>;
using rl_wspace_3d_o2_sap_type =
    manip_dynamic_env<rl_jspace_3d_o2_type, sap_interpolation_tag>;
using rl_wspace_4d_o2_sap_type =
    manip_dynamic_env<rl_jspace_4d_o2_type, sap_interpolation_tag>;
using rl_wspace_5d_o2_sap_type =
    manip_dynamic_env<rl_jspace_5d_o2_type, sap_interpolation_tag>;
using rl_wspace_6d_o2_sap_type =
    manip_dynamic_env<rl_jspace_6d_o2_type, sap_interpolation_tag>;
using rl_wspace_7d_o2_sap_type =
    manip_dynamic_env<rl_jspace_7d_o2_type, sap_interpolation_tag>;
using rl_wspace_8d_o2_sap_type =
    manip_dynamic_env<rl_jspace_8d_o2_type, sap_interpolation_tag>;
using rl_wspace_9d_o2_sap_type =
    manip_dynamic_env<rl_jspace_9d_o2_type, sap_interpolation_tag>;

// Direct Kinematics mapping for a 0th-order N-dof system.
using DKmap_1d_o0_type = manip_direct_kin_map;
using DKmap_2d_o0_type = manip_direct_kin_map;
using DKmap_3d_o0_type = manip_direct_kin_map;
using DKmap_4d_o0_type = manip_direct_kin_map;
using DKmap_5d_o0_type = manip_direct_kin_map;
using DKmap_6d_o0_type = manip_direct_kin_map;
using DKmap_7d_o0_type = manip_direct_kin_map;
using DKmap_8d_o0_type = manip_direct_kin_map;
using DKmap_9d_o0_type = manip_direct_kin_map;

// Direct Kinematics mapping for a 1st-order N-dof system.
using DKmap_1d_o1_type = manip_direct_kin_map;
using DKmap_2d_o1_type = manip_direct_kin_map;
using DKmap_3d_o1_type = manip_direct_kin_map;
using DKmap_4d_o1_type = manip_direct_kin_map;
using DKmap_5d_o1_type = manip_direct_kin_map;
using DKmap_6d_o1_type = manip_direct_kin_map;
using DKmap_7d_o1_type = manip_direct_kin_map;
using DKmap_8d_o1_type = manip_direct_kin_map;
using DKmap_9d_o1_type = manip_direct_kin_map;

// Direct Kinematics mapping for a 2nd-order N-dof system.
using DKmap_1d_o2_type = manip_direct_kin_map;
using DKmap_2d_o2_type = manip_direct_kin_map;
using DKmap_3d_o2_type = manip_direct_kin_map;
using DKmap_4d_o2_type = manip_direct_kin_map;
using DKmap_5d_o2_type = manip_direct_kin_map;
using DKmap_6d_o2_type = manip_direct_kin_map;
using DKmap_7d_o2_type = manip_direct_kin_map;
using DKmap_8d_o2_type = manip_direct_kin_map;
using DKmap_9d_o2_type = manip_direct_kin_map;

// Direct Kinematics mapping for a rate-limited 0th-order N-dof system.
using rl_DKmap_1d_o0_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o0_type>;
using rl_DKmap_2d_o0_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_2d_o0_type>;
using rl_DKmap_3d_o0_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_3d_o0_type>;
using rl_DKmap_4d_o0_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_4d_o0_type>;
using rl_DKmap_5d_o0_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_5d_o0_type>;
using rl_DKmap_6d_o0_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_6d_o0_type>;
using rl_DKmap_7d_o0_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_7d_o0_type>;
using rl_DKmap_8d_o0_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_8d_o0_type>;
using rl_DKmap_9d_o0_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_9d_o0_type>;

// Direct Kinematics mapping for a rate-limited 1st-order N-dof system.
using rl_DKmap_1d_o1_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o1_type>;
using rl_DKmap_2d_o1_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o1_type>;
using rl_DKmap_3d_o1_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o1_type>;
using rl_DKmap_4d_o1_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o1_type>;
using rl_DKmap_5d_o1_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o1_type>;
using rl_DKmap_6d_o1_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o1_type>;
using rl_DKmap_7d_o1_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o1_type>;
using rl_DKmap_8d_o1_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o1_type>;
using rl_DKmap_9d_o1_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o1_type>;

// Direct Kinematics mapping for a rate-limited 2nd-order N-dof system.
using rl_DKmap_1d_o2_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o2_type>;
using rl_DKmap_2d_o2_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o2_type>;
using rl_DKmap_3d_o2_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o2_type>;
using rl_DKmap_4d_o2_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o2_type>;
using rl_DKmap_5d_o2_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o2_type>;
using rl_DKmap_6d_o2_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o2_type>;
using rl_DKmap_7d_o2_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o2_type>;
using rl_DKmap_8d_o2_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o2_type>;
using rl_DKmap_9d_o2_type =
    manip_rl_direct_kin_map<joint_limits_mapping<double>, jspace_1d_o2_type>;

}  // namespace ReaK::pp

#endif

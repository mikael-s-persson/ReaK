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

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/named_object.hpp>

#include "joint_space_topologies.hpp"
#include "se3_topologies.hpp"
#include "joint_space_limits.hpp"

#include <ReaK/topologies/interpolation/linear_interp.hpp>
#include <ReaK/topologies/interpolation/cubic_hermite_interp.hpp>
#include <ReaK/topologies/interpolation/quintic_hermite_interp.hpp>
#include <ReaK/topologies/interpolation/sustained_velocity_pulse.hpp>
#include <ReaK/topologies/interpolation/sustained_acceleration_pulse.hpp>

#include "manip_free_workspace.hpp"
#include "manip_free_dynamic_workspace.hpp"

namespace ReaK {

namespace pp {


typedef Ndof_0th_order_rl_space< double, 1, inf_norm_tuple_distance >::type
  rl_jspace_1d_o0_type; ///< Rate-limited 1-dof Joint-space of order 0, with infinity-norm distance metric.
typedef Ndof_0th_order_rl_space< double, 2, inf_norm_tuple_distance >::type
  rl_jspace_2d_o0_type; ///< Rate-limited 2-dof Joint-space of order 0, with infinity-norm distance metric.
typedef Ndof_0th_order_rl_space< double, 3, inf_norm_tuple_distance >::type
  rl_jspace_3d_o0_type; ///< Rate-limited 3-dof Joint-space of order 0, with infinity-norm distance metric.
typedef Ndof_0th_order_rl_space< double, 4, inf_norm_tuple_distance >::type
  rl_jspace_4d_o0_type; ///< Rate-limited 4-dof Joint-space of order 0, with infinity-norm distance metric.
typedef Ndof_0th_order_rl_space< double, 5, inf_norm_tuple_distance >::type
  rl_jspace_5d_o0_type; ///< Rate-limited 5-dof Joint-space of order 0, with infinity-norm distance metric.
typedef Ndof_0th_order_rl_space< double, 6, inf_norm_tuple_distance >::type
  rl_jspace_6d_o0_type; ///< Rate-limited 6-dof Joint-space of order 0, with infinity-norm distance metric.
typedef Ndof_0th_order_rl_space< double, 7, inf_norm_tuple_distance >::type
  rl_jspace_7d_o0_type; ///< Rate-limited 7-dof Joint-space of order 0, with infinity-norm distance metric.
typedef Ndof_0th_order_rl_space< double, 8, inf_norm_tuple_distance >::type
  rl_jspace_8d_o0_type; ///< Rate-limited 8-dof Joint-space of order 0, with infinity-norm distance metric.
typedef Ndof_0th_order_rl_space< double, 9, inf_norm_tuple_distance >::type
  rl_jspace_9d_o0_type; ///< Rate-limited 9-dof Joint-space of order 0, with infinity-norm distance metric.

typedef Ndof_1st_order_rl_space< double, 1, inf_norm_tuple_distance >::type
  rl_jspace_1d_o1_type; ///< Rate-limited 1-dof Joint-space of order 1, with infinity-norm distance metric.
typedef Ndof_1st_order_rl_space< double, 2, inf_norm_tuple_distance >::type
  rl_jspace_2d_o1_type; ///< Rate-limited 2-dof Joint-space of order 1, with infinity-norm distance metric.
typedef Ndof_1st_order_rl_space< double, 3, inf_norm_tuple_distance >::type
  rl_jspace_3d_o1_type; ///< Rate-limited 3-dof Joint-space of order 1, with infinity-norm distance metric.
typedef Ndof_1st_order_rl_space< double, 4, inf_norm_tuple_distance >::type
  rl_jspace_4d_o1_type; ///< Rate-limited 4-dof Joint-space of order 1, with infinity-norm distance metric.
typedef Ndof_1st_order_rl_space< double, 5, inf_norm_tuple_distance >::type
  rl_jspace_5d_o1_type; ///< Rate-limited 5-dof Joint-space of order 1, with infinity-norm distance metric.
typedef Ndof_1st_order_rl_space< double, 6, inf_norm_tuple_distance >::type
  rl_jspace_6d_o1_type; ///< Rate-limited 6-dof Joint-space of order 1, with infinity-norm distance metric.
typedef Ndof_1st_order_rl_space< double, 7, inf_norm_tuple_distance >::type
  rl_jspace_7d_o1_type; ///< Rate-limited 7-dof Joint-space of order 1, with infinity-norm distance metric.
typedef Ndof_1st_order_rl_space< double, 8, inf_norm_tuple_distance >::type
  rl_jspace_8d_o1_type; ///< Rate-limited 8-dof Joint-space of order 1, with infinity-norm distance metric.
typedef Ndof_1st_order_rl_space< double, 9, inf_norm_tuple_distance >::type
  rl_jspace_9d_o1_type; ///< Rate-limited 9-dof Joint-space of order 1, with infinity-norm distance metric.

typedef Ndof_2nd_order_rl_space< double, 1, inf_norm_tuple_distance >::type
  rl_jspace_1d_o2_type; ///< Rate-limited 1-dof Joint-space of order 2, with infinity-norm distance metric.
typedef Ndof_2nd_order_rl_space< double, 2, inf_norm_tuple_distance >::type
  rl_jspace_2d_o2_type; ///< Rate-limited 2-dof Joint-space of order 2, with infinity-norm distance metric.
typedef Ndof_2nd_order_rl_space< double, 3, inf_norm_tuple_distance >::type
  rl_jspace_3d_o2_type; ///< Rate-limited 3-dof Joint-space of order 2, with infinity-norm distance metric.
typedef Ndof_2nd_order_rl_space< double, 4, inf_norm_tuple_distance >::type
  rl_jspace_4d_o2_type; ///< Rate-limited 4-dof Joint-space of order 2, with infinity-norm distance metric.
typedef Ndof_2nd_order_rl_space< double, 5, inf_norm_tuple_distance >::type
  rl_jspace_5d_o2_type; ///< Rate-limited 5-dof Joint-space of order 2, with infinity-norm distance metric.
typedef Ndof_2nd_order_rl_space< double, 6, inf_norm_tuple_distance >::type
  rl_jspace_6d_o2_type; ///< Rate-limited 6-dof Joint-space of order 2, with infinity-norm distance metric.
typedef Ndof_2nd_order_rl_space< double, 7, inf_norm_tuple_distance >::type
  rl_jspace_7d_o2_type; ///< Rate-limited 7-dof Joint-space of order 2, with infinity-norm distance metric.
typedef Ndof_2nd_order_rl_space< double, 8, inf_norm_tuple_distance >::type
  rl_jspace_8d_o2_type; ///< Rate-limited 8-dof Joint-space of order 2, with infinity-norm distance metric.
typedef Ndof_2nd_order_rl_space< double, 9, inf_norm_tuple_distance >::type
  rl_jspace_9d_o2_type; ///< Rate-limited 9-dof Joint-space of order 2, with infinity-norm distance metric.


typedef Ndof_0th_order_rl_space< double, 1, euclidean_tuple_distance >::type
  rle_jspace_1d_o0_type; ///< Rate-limited 1-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_rl_space< double, 2, euclidean_tuple_distance >::type
  rle_jspace_2d_o0_type; ///< Rate-limited 2-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_rl_space< double, 3, euclidean_tuple_distance >::type
  rle_jspace_3d_o0_type; ///< Rate-limited 3-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_rl_space< double, 4, euclidean_tuple_distance >::type
  rle_jspace_4d_o0_type; ///< Rate-limited 4-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_rl_space< double, 5, euclidean_tuple_distance >::type
  rle_jspace_5d_o0_type; ///< Rate-limited 5-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_rl_space< double, 6, euclidean_tuple_distance >::type
  rle_jspace_6d_o0_type; ///< Rate-limited 6-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_rl_space< double, 7, euclidean_tuple_distance >::type
  rle_jspace_7d_o0_type; ///< Rate-limited 7-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_rl_space< double, 8, euclidean_tuple_distance >::type
  rle_jspace_8d_o0_type; ///< Rate-limited 8-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_rl_space< double, 9, euclidean_tuple_distance >::type
  rle_jspace_9d_o0_type; ///< Rate-limited 9-dof Joint-space of order 0, with euclidean distance metric.

typedef Ndof_1st_order_rl_space< double, 1, euclidean_tuple_distance >::type
  rle_jspace_1d_o1_type; ///< Rate-limited 1-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_rl_space< double, 2, euclidean_tuple_distance >::type
  rle_jspace_2d_o1_type; ///< Rate-limited 2-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_rl_space< double, 3, euclidean_tuple_distance >::type
  rle_jspace_3d_o1_type; ///< Rate-limited 3-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_rl_space< double, 4, euclidean_tuple_distance >::type
  rle_jspace_4d_o1_type; ///< Rate-limited 4-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_rl_space< double, 5, euclidean_tuple_distance >::type
  rle_jspace_5d_o1_type; ///< Rate-limited 5-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_rl_space< double, 6, euclidean_tuple_distance >::type
  rle_jspace_6d_o1_type; ///< Rate-limited 6-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_rl_space< double, 7, euclidean_tuple_distance >::type
  rle_jspace_7d_o1_type; ///< Rate-limited 7-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_rl_space< double, 8, euclidean_tuple_distance >::type
  rle_jspace_8d_o1_type; ///< Rate-limited 8-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_rl_space< double, 9, euclidean_tuple_distance >::type
  rle_jspace_9d_o1_type; ///< Rate-limited 9-dof Joint-space of order 1, with euclidean distance metric.

typedef Ndof_2nd_order_rl_space< double, 1, euclidean_tuple_distance >::type
  rle_jspace_1d_o2_type; ///< Rate-limited 1-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_rl_space< double, 2, euclidean_tuple_distance >::type
  rle_jspace_2d_o2_type; ///< Rate-limited 2-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_rl_space< double, 3, euclidean_tuple_distance >::type
  rle_jspace_3d_o2_type; ///< Rate-limited 3-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_rl_space< double, 4, euclidean_tuple_distance >::type
  rle_jspace_4d_o2_type; ///< Rate-limited 4-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_rl_space< double, 5, euclidean_tuple_distance >::type
  rle_jspace_5d_o2_type; ///< Rate-limited 5-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_rl_space< double, 6, euclidean_tuple_distance >::type
  rle_jspace_6d_o2_type; ///< Rate-limited 6-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_rl_space< double, 7, euclidean_tuple_distance >::type
  rle_jspace_7d_o2_type; ///< Rate-limited 7-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_rl_space< double, 8, euclidean_tuple_distance >::type
  rle_jspace_8d_o2_type; ///< Rate-limited 8-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_rl_space< double, 9, euclidean_tuple_distance >::type
  rle_jspace_9d_o2_type; ///< Rate-limited 9-dof Joint-space of order 2, with euclidean distance metric.


typedef Ndof_0th_order_space< double, 1, euclidean_tuple_distance >::type
  jspace_1d_o0_type; ///< 1-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_space< double, 2, euclidean_tuple_distance >::type
  jspace_2d_o0_type; ///< 2-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_space< double, 3, euclidean_tuple_distance >::type
  jspace_3d_o0_type; ///< 3-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_space< double, 4, euclidean_tuple_distance >::type
  jspace_4d_o0_type; ///< 4-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_space< double, 5, euclidean_tuple_distance >::type
  jspace_5d_o0_type; ///< 5-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_space< double, 6, euclidean_tuple_distance >::type
  jspace_6d_o0_type; ///< 6-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_space< double, 7, euclidean_tuple_distance >::type
  jspace_7d_o0_type; ///< 7-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_space< double, 8, euclidean_tuple_distance >::type
  jspace_8d_o0_type; ///< 8-dof Joint-space of order 0, with euclidean distance metric.
typedef Ndof_0th_order_space< double, 9, euclidean_tuple_distance >::type
  jspace_9d_o0_type; ///< 9-dof Joint-space of order 0, with euclidean distance metric.

typedef Ndof_1st_order_space< double, 1, euclidean_tuple_distance >::type
  jspace_1d_o1_type; ///< 1-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_space< double, 2, euclidean_tuple_distance >::type
  jspace_2d_o1_type; ///< 2-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_space< double, 3, euclidean_tuple_distance >::type
  jspace_3d_o1_type; ///< 3-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_space< double, 4, euclidean_tuple_distance >::type
  jspace_4d_o1_type; ///< 4-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_space< double, 5, euclidean_tuple_distance >::type
  jspace_5d_o1_type; ///< 5-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_space< double, 6, euclidean_tuple_distance >::type
  jspace_6d_o1_type; ///< 6-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_space< double, 7, euclidean_tuple_distance >::type
  jspace_7d_o1_type; ///< 7-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_space< double, 8, euclidean_tuple_distance >::type
  jspace_8d_o1_type; ///< 8-dof Joint-space of order 1, with euclidean distance metric.
typedef Ndof_1st_order_space< double, 9, euclidean_tuple_distance >::type
  jspace_9d_o1_type; ///< 9-dof Joint-space of order 1, with euclidean distance metric.

typedef Ndof_2nd_order_space< double, 1, euclidean_tuple_distance >::type
  jspace_1d_o2_type; ///< 1-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_space< double, 2, euclidean_tuple_distance >::type
  jspace_2d_o2_type; ///< 2-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_space< double, 3, euclidean_tuple_distance >::type
  jspace_3d_o2_type; ///< 3-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_space< double, 4, euclidean_tuple_distance >::type
  jspace_4d_o2_type; ///< 4-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_space< double, 5, euclidean_tuple_distance >::type
  jspace_5d_o2_type; ///< 5-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_space< double, 6, euclidean_tuple_distance >::type
  jspace_6d_o2_type; ///< 6-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_space< double, 7, euclidean_tuple_distance >::type
  jspace_7d_o2_type; ///< 7-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_space< double, 8, euclidean_tuple_distance >::type
  jspace_8d_o2_type; ///< 8-dof Joint-space of order 2, with euclidean distance metric.
typedef Ndof_2nd_order_space< double, 9, euclidean_tuple_distance >::type
  jspace_9d_o2_type; ///< 9-dof Joint-space of order 2, with euclidean distance metric.


typedef se2_0th_order_topology< double >::type eespace_2D_o0_type; ///< 2D End-effector-space (SE(2)) of order 0.
typedef se2_1st_order_topology< double >::type eespace_2D_o1_type; ///< 2D End-effector-space (SE(2)) of order 1.
typedef se2_2nd_order_topology< double >::type eespace_2D_o2_type; ///< 2D End-effector-space (SE(2)) of order 2.

typedef se2_0th_order_rl_topology< double >::type
  rl_eespace_2D_o0_type; ///< Rate-limited 2D End-effector-space (SE(2)) of order 0.
typedef se2_1st_order_rl_topology< double >::type
  rl_eespace_2D_o1_type; ///< Rate-limited 2D End-effector-space (SE(2)) of order 1.
typedef se2_2nd_order_rl_topology< double >::type
  rl_eespace_2D_o2_type; ///< Rate-limited 2D End-effector-space (SE(2)) of order 2.

typedef se3_0th_order_topology< double >::type eespace_3D_o0_type; ///< 3D End-effector-space (SE(3)) of order 0.
typedef se3_1st_order_topology< double >::type eespace_3D_o1_type; ///< 3D End-effector-space (SE(3)) of order 1.
typedef se3_2nd_order_topology< double >::type eespace_3D_o2_type; ///< 3D End-effector-space (SE(3)) of order 2.

typedef se3_0th_order_rl_topology< double >::type
  rl_eespace_3D_o0_type; ///< Rate-limited 3D End-effector-space (SE(3)) of order 0.
typedef se3_1st_order_rl_topology< double >::type
  rl_eespace_3D_o1_type; ///< Rate-limited 3D End-effector-space (SE(3)) of order 1.
typedef se3_2nd_order_rl_topology< double >::type
  rl_eespace_3D_o2_type; ///< Rate-limited 3D End-effector-space (SE(3)) of order 2.


typedef manip_quasi_static_env< rl_jspace_1d_o0_type, linear_interpolation_tag >
  rl_wspace_1d_o0_i1_type; ///< Rate-limited 1-dof static-workspace of order 0, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_2d_o0_type, linear_interpolation_tag >
  rl_wspace_2d_o0_i1_type; ///< Rate-limited 2-dof static-workspace of order 0, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_3d_o0_type, linear_interpolation_tag >
  rl_wspace_3d_o0_i1_type; ///< Rate-limited 3-dof static-workspace of order 0, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_4d_o0_type, linear_interpolation_tag >
  rl_wspace_4d_o0_i1_type; ///< Rate-limited 4-dof static-workspace of order 0, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_5d_o0_type, linear_interpolation_tag >
  rl_wspace_5d_o0_i1_type; ///< Rate-limited 5-dof static-workspace of order 0, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_6d_o0_type, linear_interpolation_tag >
  rl_wspace_6d_o0_i1_type; ///< Rate-limited 6-dof static-workspace of order 0, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_7d_o0_type, linear_interpolation_tag >
  rl_wspace_7d_o0_i1_type; ///< Rate-limited 7-dof static-workspace of order 0, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_8d_o0_type, linear_interpolation_tag >
  rl_wspace_8d_o0_i1_type; ///< Rate-limited 8-dof static-workspace of order 0, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_9d_o0_type, linear_interpolation_tag >
  rl_wspace_9d_o0_i1_type; ///< Rate-limited 9-dof static-workspace of order 0, with linear interpolation.

typedef manip_quasi_static_env< rl_jspace_1d_o1_type, linear_interpolation_tag >
  rl_wspace_1d_o1_i1_type; ///< Rate-limited 1-dof static-workspace of order 1, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_2d_o1_type, linear_interpolation_tag >
  rl_wspace_2d_o1_i1_type; ///< Rate-limited 2-dof static-workspace of order 1, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_3d_o1_type, linear_interpolation_tag >
  rl_wspace_3d_o1_i1_type; ///< Rate-limited 3-dof static-workspace of order 1, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_4d_o1_type, linear_interpolation_tag >
  rl_wspace_4d_o1_i1_type; ///< Rate-limited 4-dof static-workspace of order 1, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_5d_o1_type, linear_interpolation_tag >
  rl_wspace_5d_o1_i1_type; ///< Rate-limited 5-dof static-workspace of order 1, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_6d_o1_type, linear_interpolation_tag >
  rl_wspace_6d_o1_i1_type; ///< Rate-limited 6-dof static-workspace of order 1, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_7d_o1_type, linear_interpolation_tag >
  rl_wspace_7d_o1_i1_type; ///< Rate-limited 7-dof static-workspace of order 1, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_8d_o1_type, linear_interpolation_tag >
  rl_wspace_8d_o1_i1_type; ///< Rate-limited 8-dof static-workspace of order 1, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_9d_o1_type, linear_interpolation_tag >
  rl_wspace_9d_o1_i1_type; ///< Rate-limited 9-dof static-workspace of order 1, with linear interpolation.

typedef manip_quasi_static_env< rl_jspace_1d_o2_type, linear_interpolation_tag >
  rl_wspace_1d_o2_i1_type; ///< Rate-limited 1-dof static-workspace of order 2, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_2d_o2_type, linear_interpolation_tag >
  rl_wspace_2d_o2_i1_type; ///< Rate-limited 2-dof static-workspace of order 2, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_3d_o2_type, linear_interpolation_tag >
  rl_wspace_3d_o2_i1_type; ///< Rate-limited 3-dof static-workspace of order 2, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_4d_o2_type, linear_interpolation_tag >
  rl_wspace_4d_o2_i1_type; ///< Rate-limited 4-dof static-workspace of order 2, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_5d_o2_type, linear_interpolation_tag >
  rl_wspace_5d_o2_i1_type; ///< Rate-limited 5-dof static-workspace of order 2, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_6d_o2_type, linear_interpolation_tag >
  rl_wspace_6d_o2_i1_type; ///< Rate-limited 6-dof static-workspace of order 2, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_7d_o2_type, linear_interpolation_tag >
  rl_wspace_7d_o2_i1_type; ///< Rate-limited 7-dof static-workspace of order 2, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_8d_o2_type, linear_interpolation_tag >
  rl_wspace_8d_o2_i1_type; ///< Rate-limited 8-dof static-workspace of order 2, with linear interpolation.
typedef manip_quasi_static_env< rl_jspace_9d_o2_type, linear_interpolation_tag >
  rl_wspace_9d_o2_i1_type; ///< Rate-limited 9-dof static-workspace of order 2, with linear interpolation.

typedef manip_quasi_static_env< rl_jspace_1d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_1d_o1_i3_type; ///< Rate-limited 1-dof static-workspace of order 1, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_2d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_2d_o1_i3_type; ///< Rate-limited 2-dof static-workspace of order 1, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_3d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_3d_o1_i3_type; ///< Rate-limited 3-dof static-workspace of order 1, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_4d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_4d_o1_i3_type; ///< Rate-limited 4-dof static-workspace of order 1, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_5d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_5d_o1_i3_type; ///< Rate-limited 5-dof static-workspace of order 1, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_6d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_6d_o1_i3_type; ///< Rate-limited 6-dof static-workspace of order 1, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_7d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_7d_o1_i3_type; ///< Rate-limited 7-dof static-workspace of order 1, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_8d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_8d_o1_i3_type; ///< Rate-limited 8-dof static-workspace of order 1, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_9d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_9d_o1_i3_type; ///< Rate-limited 9-dof static-workspace of order 1, with cubic interpolation.

typedef manip_quasi_static_env< rl_jspace_1d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_1d_o2_i3_type; ///< Rate-limited 1-dof static-workspace of order 2, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_2d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_2d_o2_i3_type; ///< Rate-limited 2-dof static-workspace of order 2, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_3d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_3d_o2_i3_type; ///< Rate-limited 3-dof static-workspace of order 2, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_4d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_4d_o2_i3_type; ///< Rate-limited 4-dof static-workspace of order 2, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_5d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_5d_o2_i3_type; ///< Rate-limited 5-dof static-workspace of order 2, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_6d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_6d_o2_i3_type; ///< Rate-limited 6-dof static-workspace of order 2, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_7d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_7d_o2_i3_type; ///< Rate-limited 7-dof static-workspace of order 2, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_8d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_8d_o2_i3_type; ///< Rate-limited 8-dof static-workspace of order 2, with cubic interpolation.
typedef manip_quasi_static_env< rl_jspace_9d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_9d_o2_i3_type; ///< Rate-limited 9-dof static-workspace of order 2, with cubic interpolation.

typedef manip_quasi_static_env< rl_jspace_1d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_1d_o2_i5_type; ///< Rate-limited 1-dof static-workspace of order 2, with quintic interpolation.
typedef manip_quasi_static_env< rl_jspace_2d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_2d_o2_i5_type; ///< Rate-limited 2-dof static-workspace of order 2, with quintic interpolation.
typedef manip_quasi_static_env< rl_jspace_3d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_3d_o2_i5_type; ///< Rate-limited 3-dof static-workspace of order 2, with quintic interpolation.
typedef manip_quasi_static_env< rl_jspace_4d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_4d_o2_i5_type; ///< Rate-limited 4-dof static-workspace of order 2, with quintic interpolation.
typedef manip_quasi_static_env< rl_jspace_5d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_5d_o2_i5_type; ///< Rate-limited 5-dof static-workspace of order 2, with quintic interpolation.
typedef manip_quasi_static_env< rl_jspace_6d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_6d_o2_i5_type; ///< Rate-limited 6-dof static-workspace of order 2, with quintic interpolation.
typedef manip_quasi_static_env< rl_jspace_7d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_7d_o2_i5_type; ///< Rate-limited 7-dof static-workspace of order 2, with quintic interpolation.
typedef manip_quasi_static_env< rl_jspace_8d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_8d_o2_i5_type; ///< Rate-limited 8-dof static-workspace of order 2, with quintic interpolation.
typedef manip_quasi_static_env< rl_jspace_9d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_9d_o2_i5_type; ///< Rate-limited 9-dof static-workspace of order 2, with quintic interpolation.

typedef manip_quasi_static_env< rl_jspace_1d_o1_type, svp_interpolation_tag >
  rl_wspace_1d_o1_svp_type; ///< Rate-limited 1-dof static-workspace of order 1, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_2d_o1_type, svp_interpolation_tag >
  rl_wspace_2d_o1_svp_type; ///< Rate-limited 2-dof static-workspace of order 1, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_3d_o1_type, svp_interpolation_tag >
  rl_wspace_3d_o1_svp_type; ///< Rate-limited 3-dof static-workspace of order 1, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_4d_o1_type, svp_interpolation_tag >
  rl_wspace_4d_o1_svp_type; ///< Rate-limited 4-dof static-workspace of order 1, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_5d_o1_type, svp_interpolation_tag >
  rl_wspace_5d_o1_svp_type; ///< Rate-limited 5-dof static-workspace of order 1, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_6d_o1_type, svp_interpolation_tag >
  rl_wspace_6d_o1_svp_type; ///< Rate-limited 6-dof static-workspace of order 1, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_7d_o1_type, svp_interpolation_tag >
  rl_wspace_7d_o1_svp_type; ///< Rate-limited 7-dof static-workspace of order 1, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_8d_o1_type, svp_interpolation_tag >
  rl_wspace_8d_o1_svp_type; ///< Rate-limited 8-dof static-workspace of order 1, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_9d_o1_type, svp_interpolation_tag >
  rl_wspace_9d_o1_svp_type; ///< Rate-limited 9-dof static-workspace of order 1, with SVP interpolation.

typedef manip_quasi_static_env< rl_jspace_1d_o2_type, svp_interpolation_tag >
  rl_wspace_1d_o2_svp_type; ///< Rate-limited 1-dof static-workspace of order 2, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_2d_o2_type, svp_interpolation_tag >
  rl_wspace_2d_o2_svp_type; ///< Rate-limited 2-dof static-workspace of order 2, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_3d_o2_type, svp_interpolation_tag >
  rl_wspace_3d_o2_svp_type; ///< Rate-limited 3-dof static-workspace of order 2, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_4d_o2_type, svp_interpolation_tag >
  rl_wspace_4d_o2_svp_type; ///< Rate-limited 4-dof static-workspace of order 2, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_5d_o2_type, svp_interpolation_tag >
  rl_wspace_5d_o2_svp_type; ///< Rate-limited 5-dof static-workspace of order 2, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_6d_o2_type, svp_interpolation_tag >
  rl_wspace_6d_o2_svp_type; ///< Rate-limited 6-dof static-workspace of order 2, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_7d_o2_type, svp_interpolation_tag >
  rl_wspace_7d_o2_svp_type; ///< Rate-limited 7-dof static-workspace of order 2, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_8d_o2_type, svp_interpolation_tag >
  rl_wspace_8d_o2_svp_type; ///< Rate-limited 8-dof static-workspace of order 2, with SVP interpolation.
typedef manip_quasi_static_env< rl_jspace_9d_o2_type, svp_interpolation_tag >
  rl_wspace_9d_o2_svp_type; ///< Rate-limited 9-dof static-workspace of order 2, with SVP interpolation.

typedef manip_quasi_static_env< rl_jspace_1d_o2_type, sap_interpolation_tag >
  rl_wspace_1d_o2_sap_type; ///< Rate-limited 1-dof static-workspace of order 2, with SAP interpolation.
typedef manip_quasi_static_env< rl_jspace_2d_o2_type, sap_interpolation_tag >
  rl_wspace_2d_o2_sap_type; ///< Rate-limited 2-dof static-workspace of order 2, with SAP interpolation.
typedef manip_quasi_static_env< rl_jspace_3d_o2_type, sap_interpolation_tag >
  rl_wspace_3d_o2_sap_type; ///< Rate-limited 3-dof static-workspace of order 2, with SAP interpolation.
typedef manip_quasi_static_env< rl_jspace_4d_o2_type, sap_interpolation_tag >
  rl_wspace_4d_o2_sap_type; ///< Rate-limited 4-dof static-workspace of order 2, with SAP interpolation.
typedef manip_quasi_static_env< rl_jspace_5d_o2_type, sap_interpolation_tag >
  rl_wspace_5d_o2_sap_type; ///< Rate-limited 5-dof static-workspace of order 2, with SAP interpolation.
typedef manip_quasi_static_env< rl_jspace_6d_o2_type, sap_interpolation_tag >
  rl_wspace_6d_o2_sap_type; ///< Rate-limited 6-dof static-workspace of order 2, with SAP interpolation.
typedef manip_quasi_static_env< rl_jspace_7d_o2_type, sap_interpolation_tag >
  rl_wspace_7d_o2_sap_type; ///< Rate-limited 7-dof static-workspace of order 2, with SAP interpolation.
typedef manip_quasi_static_env< rl_jspace_8d_o2_type, sap_interpolation_tag >
  rl_wspace_8d_o2_sap_type; ///< Rate-limited 8-dof static-workspace of order 2, with SAP interpolation.
typedef manip_quasi_static_env< rl_jspace_9d_o2_type, sap_interpolation_tag >
  rl_wspace_9d_o2_sap_type; ///< Rate-limited 9-dof static-workspace of order 2, with SAP interpolation.


typedef manip_dynamic_env< rl_jspace_1d_o0_type, linear_interpolation_tag >
  rl_wspace_1d_o0_i1_type; ///< Rate-limited 1-dof dynamic-workspace of order 0, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_2d_o0_type, linear_interpolation_tag >
  rl_wspace_2d_o0_i1_type; ///< Rate-limited 2-dof dynamic-workspace of order 0, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_3d_o0_type, linear_interpolation_tag >
  rl_wspace_3d_o0_i1_type; ///< Rate-limited 3-dof dynamic-workspace of order 0, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_4d_o0_type, linear_interpolation_tag >
  rl_wspace_4d_o0_i1_type; ///< Rate-limited 4-dof dynamic-workspace of order 0, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_5d_o0_type, linear_interpolation_tag >
  rl_wspace_5d_o0_i1_type; ///< Rate-limited 5-dof dynamic-workspace of order 0, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_6d_o0_type, linear_interpolation_tag >
  rl_wspace_6d_o0_i1_type; ///< Rate-limited 6-dof dynamic-workspace of order 0, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_7d_o0_type, linear_interpolation_tag >
  rl_wspace_7d_o0_i1_type; ///< Rate-limited 7-dof dynamic-workspace of order 0, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_8d_o0_type, linear_interpolation_tag >
  rl_wspace_8d_o0_i1_type; ///< Rate-limited 8-dof dynamic-workspace of order 0, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_9d_o0_type, linear_interpolation_tag >
  rl_wspace_9d_o0_i1_type; ///< Rate-limited 9-dof dynamic-workspace of order 0, with linear interpolation.

typedef manip_dynamic_env< rl_jspace_1d_o1_type, linear_interpolation_tag >
  rl_wspace_1d_o1_i1_type; ///< Rate-limited 1-dof dynamic-workspace of order 1, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_2d_o1_type, linear_interpolation_tag >
  rl_wspace_2d_o1_i1_type; ///< Rate-limited 2-dof dynamic-workspace of order 1, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_3d_o1_type, linear_interpolation_tag >
  rl_wspace_3d_o1_i1_type; ///< Rate-limited 3-dof dynamic-workspace of order 1, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_4d_o1_type, linear_interpolation_tag >
  rl_wspace_4d_o1_i1_type; ///< Rate-limited 4-dof dynamic-workspace of order 1, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_5d_o1_type, linear_interpolation_tag >
  rl_wspace_5d_o1_i1_type; ///< Rate-limited 5-dof dynamic-workspace of order 1, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_6d_o1_type, linear_interpolation_tag >
  rl_wspace_6d_o1_i1_type; ///< Rate-limited 6-dof dynamic-workspace of order 1, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_7d_o1_type, linear_interpolation_tag >
  rl_wspace_7d_o1_i1_type; ///< Rate-limited 7-dof dynamic-workspace of order 1, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_8d_o1_type, linear_interpolation_tag >
  rl_wspace_8d_o1_i1_type; ///< Rate-limited 8-dof dynamic-workspace of order 1, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_9d_o1_type, linear_interpolation_tag >
  rl_wspace_9d_o1_i1_type; ///< Rate-limited 9-dof dynamic-workspace of order 1, with linear interpolation.

typedef manip_dynamic_env< rl_jspace_1d_o2_type, linear_interpolation_tag >
  rl_wspace_1d_o2_i1_type; ///< Rate-limited 1-dof dynamic-workspace of order 2, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_2d_o2_type, linear_interpolation_tag >
  rl_wspace_2d_o2_i1_type; ///< Rate-limited 2-dof dynamic-workspace of order 2, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_3d_o2_type, linear_interpolation_tag >
  rl_wspace_3d_o2_i1_type; ///< Rate-limited 3-dof dynamic-workspace of order 2, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_4d_o2_type, linear_interpolation_tag >
  rl_wspace_4d_o2_i1_type; ///< Rate-limited 4-dof dynamic-workspace of order 2, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_5d_o2_type, linear_interpolation_tag >
  rl_wspace_5d_o2_i1_type; ///< Rate-limited 5-dof dynamic-workspace of order 2, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_6d_o2_type, linear_interpolation_tag >
  rl_wspace_6d_o2_i1_type; ///< Rate-limited 6-dof dynamic-workspace of order 2, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_7d_o2_type, linear_interpolation_tag >
  rl_wspace_7d_o2_i1_type; ///< Rate-limited 7-dof dynamic-workspace of order 2, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_8d_o2_type, linear_interpolation_tag >
  rl_wspace_8d_o2_i1_type; ///< Rate-limited 8-dof dynamic-workspace of order 2, with linear interpolation.
typedef manip_dynamic_env< rl_jspace_9d_o2_type, linear_interpolation_tag >
  rl_wspace_9d_o2_i1_type; ///< Rate-limited 9-dof dynamic-workspace of order 2, with linear interpolation.

typedef manip_dynamic_env< rl_jspace_1d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_1d_o1_i3_type; ///< Rate-limited 1-dof dynamic-workspace of order 1, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_2d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_2d_o1_i3_type; ///< Rate-limited 2-dof dynamic-workspace of order 1, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_3d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_3d_o1_i3_type; ///< Rate-limited 3-dof dynamic-workspace of order 1, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_4d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_4d_o1_i3_type; ///< Rate-limited 4-dof dynamic-workspace of order 1, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_5d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_5d_o1_i3_type; ///< Rate-limited 5-dof dynamic-workspace of order 1, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_6d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_6d_o1_i3_type; ///< Rate-limited 6-dof dynamic-workspace of order 1, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_7d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_7d_o1_i3_type; ///< Rate-limited 7-dof dynamic-workspace of order 1, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_8d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_8d_o1_i3_type; ///< Rate-limited 8-dof dynamic-workspace of order 1, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_9d_o1_type, cubic_hermite_interpolation_tag >
  rl_wspace_9d_o1_i3_type; ///< Rate-limited 9-dof dynamic-workspace of order 1, with cubic interpolation.

typedef manip_dynamic_env< rl_jspace_1d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_1d_o2_i3_type; ///< Rate-limited 1-dof dynamic-workspace of order 2, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_2d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_2d_o2_i3_type; ///< Rate-limited 2-dof dynamic-workspace of order 2, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_3d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_3d_o2_i3_type; ///< Rate-limited 3-dof dynamic-workspace of order 2, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_4d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_4d_o2_i3_type; ///< Rate-limited 4-dof dynamic-workspace of order 2, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_5d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_5d_o2_i3_type; ///< Rate-limited 5-dof dynamic-workspace of order 2, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_6d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_6d_o2_i3_type; ///< Rate-limited 6-dof dynamic-workspace of order 2, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_7d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_7d_o2_i3_type; ///< Rate-limited 7-dof dynamic-workspace of order 2, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_8d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_8d_o2_i3_type; ///< Rate-limited 8-dof dynamic-workspace of order 2, with cubic interpolation.
typedef manip_dynamic_env< rl_jspace_9d_o2_type, cubic_hermite_interpolation_tag >
  rl_wspace_9d_o2_i3_type; ///< Rate-limited 9-dof dynamic-workspace of order 2, with cubic interpolation.

typedef manip_dynamic_env< rl_jspace_1d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_1d_o2_i5_type; ///< Rate-limited 1-dof dynamic-workspace of order 2, with quintic interpolation.
typedef manip_dynamic_env< rl_jspace_2d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_2d_o2_i5_type; ///< Rate-limited 2-dof dynamic-workspace of order 2, with quintic interpolation.
typedef manip_dynamic_env< rl_jspace_3d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_3d_o2_i5_type; ///< Rate-limited 3-dof dynamic-workspace of order 2, with quintic interpolation.
typedef manip_dynamic_env< rl_jspace_4d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_4d_o2_i5_type; ///< Rate-limited 4-dof dynamic-workspace of order 2, with quintic interpolation.
typedef manip_dynamic_env< rl_jspace_5d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_5d_o2_i5_type; ///< Rate-limited 5-dof dynamic-workspace of order 2, with quintic interpolation.
typedef manip_dynamic_env< rl_jspace_6d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_6d_o2_i5_type; ///< Rate-limited 6-dof dynamic-workspace of order 2, with quintic interpolation.
typedef manip_dynamic_env< rl_jspace_7d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_7d_o2_i5_type; ///< Rate-limited 7-dof dynamic-workspace of order 2, with quintic interpolation.
typedef manip_dynamic_env< rl_jspace_8d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_8d_o2_i5_type; ///< Rate-limited 8-dof dynamic-workspace of order 2, with quintic interpolation.
typedef manip_dynamic_env< rl_jspace_9d_o2_type, quintic_hermite_interpolation_tag >
  rl_wspace_9d_o2_i5_type; ///< Rate-limited 9-dof dynamic-workspace of order 2, with quintic interpolation.

typedef manip_dynamic_env< rl_jspace_1d_o1_type, svp_interpolation_tag >
  rl_wspace_1d_o1_svp_type; ///< Rate-limited 1-dof dynamic-workspace of order 1, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_2d_o1_type, svp_interpolation_tag >
  rl_wspace_2d_o1_svp_type; ///< Rate-limited 2-dof dynamic-workspace of order 1, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_3d_o1_type, svp_interpolation_tag >
  rl_wspace_3d_o1_svp_type; ///< Rate-limited 3-dof dynamic-workspace of order 1, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_4d_o1_type, svp_interpolation_tag >
  rl_wspace_4d_o1_svp_type; ///< Rate-limited 4-dof dynamic-workspace of order 1, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_5d_o1_type, svp_interpolation_tag >
  rl_wspace_5d_o1_svp_type; ///< Rate-limited 5-dof dynamic-workspace of order 1, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_6d_o1_type, svp_interpolation_tag >
  rl_wspace_6d_o1_svp_type; ///< Rate-limited 6-dof dynamic-workspace of order 1, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_7d_o1_type, svp_interpolation_tag >
  rl_wspace_7d_o1_svp_type; ///< Rate-limited 7-dof dynamic-workspace of order 1, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_8d_o1_type, svp_interpolation_tag >
  rl_wspace_8d_o1_svp_type; ///< Rate-limited 8-dof dynamic-workspace of order 1, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_9d_o1_type, svp_interpolation_tag >
  rl_wspace_9d_o1_svp_type; ///< Rate-limited 9-dof dynamic-workspace of order 1, with SVP interpolation.

typedef manip_dynamic_env< rl_jspace_1d_o2_type, svp_interpolation_tag >
  rl_wspace_1d_o2_svp_type; ///< Rate-limited 1-dof dynamic-workspace of order 2, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_2d_o2_type, svp_interpolation_tag >
  rl_wspace_2d_o2_svp_type; ///< Rate-limited 2-dof dynamic-workspace of order 2, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_3d_o2_type, svp_interpolation_tag >
  rl_wspace_3d_o2_svp_type; ///< Rate-limited 3-dof dynamic-workspace of order 2, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_4d_o2_type, svp_interpolation_tag >
  rl_wspace_4d_o2_svp_type; ///< Rate-limited 4-dof dynamic-workspace of order 2, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_5d_o2_type, svp_interpolation_tag >
  rl_wspace_5d_o2_svp_type; ///< Rate-limited 5-dof dynamic-workspace of order 2, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_6d_o2_type, svp_interpolation_tag >
  rl_wspace_6d_o2_svp_type; ///< Rate-limited 6-dof dynamic-workspace of order 2, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_7d_o2_type, svp_interpolation_tag >
  rl_wspace_7d_o2_svp_type; ///< Rate-limited 7-dof dynamic-workspace of order 2, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_8d_o2_type, svp_interpolation_tag >
  rl_wspace_8d_o2_svp_type; ///< Rate-limited 8-dof dynamic-workspace of order 2, with SVP interpolation.
typedef manip_dynamic_env< rl_jspace_9d_o2_type, svp_interpolation_tag >
  rl_wspace_9d_o2_svp_type; ///< Rate-limited 9-dof dynamic-workspace of order 2, with SVP interpolation.

typedef manip_dynamic_env< rl_jspace_1d_o2_type, sap_interpolation_tag >
  rl_wspace_1d_o2_sap_type; ///< Rate-limited 1-dof dynamic-workspace of order 2, with SAP interpolation.
typedef manip_dynamic_env< rl_jspace_2d_o2_type, sap_interpolation_tag >
  rl_wspace_2d_o2_sap_type; ///< Rate-limited 2-dof dynamic-workspace of order 2, with SAP interpolation.
typedef manip_dynamic_env< rl_jspace_3d_o2_type, sap_interpolation_tag >
  rl_wspace_3d_o2_sap_type; ///< Rate-limited 3-dof dynamic-workspace of order 2, with SAP interpolation.
typedef manip_dynamic_env< rl_jspace_4d_o2_type, sap_interpolation_tag >
  rl_wspace_4d_o2_sap_type; ///< Rate-limited 4-dof dynamic-workspace of order 2, with SAP interpolation.
typedef manip_dynamic_env< rl_jspace_5d_o2_type, sap_interpolation_tag >
  rl_wspace_5d_o2_sap_type; ///< Rate-limited 5-dof dynamic-workspace of order 2, with SAP interpolation.
typedef manip_dynamic_env< rl_jspace_6d_o2_type, sap_interpolation_tag >
  rl_wspace_6d_o2_sap_type; ///< Rate-limited 6-dof dynamic-workspace of order 2, with SAP interpolation.
typedef manip_dynamic_env< rl_jspace_7d_o2_type, sap_interpolation_tag >
  rl_wspace_7d_o2_sap_type; ///< Rate-limited 7-dof dynamic-workspace of order 2, with SAP interpolation.
typedef manip_dynamic_env< rl_jspace_8d_o2_type, sap_interpolation_tag >
  rl_wspace_8d_o2_sap_type; ///< Rate-limited 8-dof dynamic-workspace of order 2, with SAP interpolation.
typedef manip_dynamic_env< rl_jspace_9d_o2_type, sap_interpolation_tag >
  rl_wspace_9d_o2_sap_type; ///< Rate-limited 9-dof dynamic-workspace of order 2, with SAP interpolation.


typedef manip_direct_kin_map DKmap_1d_o0_type; ///< Direct Kinematics mapping for a 0th-order 1-dof system.
typedef manip_direct_kin_map DKmap_2d_o0_type; ///< Direct Kinematics mapping for a 0th-order 2-dof system.
typedef manip_direct_kin_map DKmap_3d_o0_type; ///< Direct Kinematics mapping for a 0th-order 3-dof system.
typedef manip_direct_kin_map DKmap_4d_o0_type; ///< Direct Kinematics mapping for a 0th-order 4-dof system.
typedef manip_direct_kin_map DKmap_5d_o0_type; ///< Direct Kinematics mapping for a 0th-order 5-dof system.
typedef manip_direct_kin_map DKmap_6d_o0_type; ///< Direct Kinematics mapping for a 0th-order 6-dof system.
typedef manip_direct_kin_map DKmap_7d_o0_type; ///< Direct Kinematics mapping for a 0th-order 7-dof system.
typedef manip_direct_kin_map DKmap_8d_o0_type; ///< Direct Kinematics mapping for a 0th-order 8-dof system.
typedef manip_direct_kin_map DKmap_9d_o0_type; ///< Direct Kinematics mapping for a 0th-order 9-dof system.

typedef manip_direct_kin_map DKmap_1d_o1_type; ///< Direct Kinematics mapping for a 1st-order 1-dof system.
typedef manip_direct_kin_map DKmap_2d_o1_type; ///< Direct Kinematics mapping for a 1st-order 2-dof system.
typedef manip_direct_kin_map DKmap_3d_o1_type; ///< Direct Kinematics mapping for a 1st-order 3-dof system.
typedef manip_direct_kin_map DKmap_4d_o1_type; ///< Direct Kinematics mapping for a 1st-order 4-dof system.
typedef manip_direct_kin_map DKmap_5d_o1_type; ///< Direct Kinematics mapping for a 1st-order 5-dof system.
typedef manip_direct_kin_map DKmap_6d_o1_type; ///< Direct Kinematics mapping for a 1st-order 6-dof system.
typedef manip_direct_kin_map DKmap_7d_o1_type; ///< Direct Kinematics mapping for a 1st-order 7-dof system.
typedef manip_direct_kin_map DKmap_8d_o1_type; ///< Direct Kinematics mapping for a 1st-order 8-dof system.
typedef manip_direct_kin_map DKmap_9d_o1_type; ///< Direct Kinematics mapping for a 1st-order 9-dof system.

typedef manip_direct_kin_map DKmap_1d_o2_type; ///< Direct Kinematics mapping for a 2nd-order 1-dof system.
typedef manip_direct_kin_map DKmap_2d_o2_type; ///< Direct Kinematics mapping for a 2nd-order 2-dof system.
typedef manip_direct_kin_map DKmap_3d_o2_type; ///< Direct Kinematics mapping for a 2nd-order 3-dof system.
typedef manip_direct_kin_map DKmap_4d_o2_type; ///< Direct Kinematics mapping for a 2nd-order 4-dof system.
typedef manip_direct_kin_map DKmap_5d_o2_type; ///< Direct Kinematics mapping for a 2nd-order 5-dof system.
typedef manip_direct_kin_map DKmap_6d_o2_type; ///< Direct Kinematics mapping for a 2nd-order 6-dof system.
typedef manip_direct_kin_map DKmap_7d_o2_type; ///< Direct Kinematics mapping for a 2nd-order 7-dof system.
typedef manip_direct_kin_map DKmap_8d_o2_type; ///< Direct Kinematics mapping for a 2nd-order 8-dof system.
typedef manip_direct_kin_map DKmap_9d_o2_type; ///< Direct Kinematics mapping for a 2nd-order 9-dof system.


typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o0_type >
  rl_DKmap_1d_o0_type; ///< Direct Kinematics mapping for a rate-limited 0th-order 1-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_2d_o0_type >
  rl_DKmap_2d_o0_type; ///< Direct Kinematics mapping for a rate-limited 0th-order 2-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_3d_o0_type >
  rl_DKmap_3d_o0_type; ///< Direct Kinematics mapping for a rate-limited 0th-order 3-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_4d_o0_type >
  rl_DKmap_4d_o0_type; ///< Direct Kinematics mapping for a rate-limited 0th-order 4-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_5d_o0_type >
  rl_DKmap_5d_o0_type; ///< Direct Kinematics mapping for a rate-limited 0th-order 5-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_6d_o0_type >
  rl_DKmap_6d_o0_type; ///< Direct Kinematics mapping for a rate-limited 0th-order 6-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_7d_o0_type >
  rl_DKmap_7d_o0_type; ///< Direct Kinematics mapping for a rate-limited 0th-order 7-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_8d_o0_type >
  rl_DKmap_8d_o0_type; ///< Direct Kinematics mapping for a rate-limited 0th-order 8-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_9d_o0_type >
  rl_DKmap_9d_o0_type; ///< Direct Kinematics mapping for a rate-limited 0th-order 9-dof system.

typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o1_type >
  rl_DKmap_1d_o1_type; ///< Direct Kinematics mapping for a rate-limited 1st-order 1-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o1_type >
  rl_DKmap_2d_o1_type; ///< Direct Kinematics mapping for a rate-limited 1st-order 2-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o1_type >
  rl_DKmap_3d_o1_type; ///< Direct Kinematics mapping for a rate-limited 1st-order 3-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o1_type >
  rl_DKmap_4d_o1_type; ///< Direct Kinematics mapping for a rate-limited 1st-order 4-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o1_type >
  rl_DKmap_5d_o1_type; ///< Direct Kinematics mapping for a rate-limited 1st-order 5-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o1_type >
  rl_DKmap_6d_o1_type; ///< Direct Kinematics mapping for a rate-limited 1st-order 6-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o1_type >
  rl_DKmap_7d_o1_type; ///< Direct Kinematics mapping for a rate-limited 1st-order 7-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o1_type >
  rl_DKmap_8d_o1_type; ///< Direct Kinematics mapping for a rate-limited 1st-order 8-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o1_type >
  rl_DKmap_9d_o1_type; ///< Direct Kinematics mapping for a rate-limited 1st-order 9-dof system.

typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o2_type >
  rl_DKmap_1d_o2_type; ///< Direct Kinematics mapping for a rate-limited 2nd-order 1-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o2_type >
  rl_DKmap_2d_o2_type; ///< Direct Kinematics mapping for a rate-limited 2nd-order 2-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o2_type >
  rl_DKmap_3d_o2_type; ///< Direct Kinematics mapping for a rate-limited 2nd-order 3-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o2_type >
  rl_DKmap_4d_o2_type; ///< Direct Kinematics mapping for a rate-limited 2nd-order 4-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o2_type >
  rl_DKmap_5d_o2_type; ///< Direct Kinematics mapping for a rate-limited 2nd-order 5-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o2_type >
  rl_DKmap_6d_o2_type; ///< Direct Kinematics mapping for a rate-limited 2nd-order 6-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o2_type >
  rl_DKmap_7d_o2_type; ///< Direct Kinematics mapping for a rate-limited 2nd-order 7-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o2_type >
  rl_DKmap_8d_o2_type; ///< Direct Kinematics mapping for a rate-limited 2nd-order 8-dof system.
typedef manip_rl_direct_kin_map< joint_limits_mapping< double >, jspace_1d_o2_type >
  rl_DKmap_9d_o2_type; ///< Direct Kinematics mapping for a rate-limited 2nd-order 9-dof system.
};
};

#endif

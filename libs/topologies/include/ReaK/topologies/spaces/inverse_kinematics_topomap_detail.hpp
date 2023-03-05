/*
 * \file inverse_kinematics_topomap_detail.hpp
 *
 * DETAILS OF:
 * This library provides classes that define topological mappings between a joint-space (generalized
 * coordinates) and the end-effector frame of a serial manipulator.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2012
 */

/*
 *    Copyright 2011 Sven Mikael Persson
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

#ifndef REAK_INVERSE_KINEMATICS_TOPOMAP_DETAIL_HPP
#define REAK_INVERSE_KINEMATICS_TOPOMAP_DETAIL_HPP

#include "ReaK/math/lin_alg/arithmetic_tuple.hpp"

#include "ReaK/topologies/spaces/Ndof_spaces.hpp"
#include "ReaK/topologies/spaces/joint_space_limits.hpp"
#include "ReaK/topologies/spaces/joint_space_topologies.hpp"
#include "ReaK/topologies/spaces/se2_topologies.hpp"
#include "ReaK/topologies/spaces/se3_topologies.hpp"

#include <type_traits>

namespace ReaK::pp::detail {
namespace {

template <typename PointType, typename InSpace>
void read_one_joint_coord_impl(
    PointType& pt, const InSpace& space_in, std::size_t& gen_i,
    std::size_t& f2d_i, std::size_t& f3d_i,
    const std::shared_ptr<kte::direct_kinematics_model>& model) {
  if constexpr (is_normal_joint_space_v<InSpace>) {
    set_gen_coord(pt, *(model->getCoord(gen_i++)));
  } else if constexpr (is_Ndof_space_v<InSpace>) {
    for (int i = 0; i < get<0>(pt).size(); ++i) {
      set_gen_coord(pt, i, *(model->getCoord(gen_i++)));
    }
  } else if constexpr (is_se2_space_v<InSpace>) {
    set_frame_2D(pt, *(model->getFrame2D(f2d_i++)));
  } else if constexpr (is_se3_space_v<InSpace>) {
    set_frame_3D(pt, *(model->getFrame3D(f3d_i++)));
  } else {
    tuple_for_each(pt, space_in, [&](auto& pt_elem, const auto& space_in_elem) {
      read_one_joint_coord_impl(pt_elem, space_in_elem, gen_i, f2d_i, f3d_i,
                                model);
    });
  }
}

template <typename PointType, typename InSpaceTuple>
void read_joint_coordinates_impl(
    PointType& pt, const InSpaceTuple& space_in,
    const std::shared_ptr<kte::direct_kinematics_model>& model) {
  std::size_t gen_i = 0;
  std::size_t f2d_i = 0;
  std::size_t f3d_i = 0;
  read_one_joint_coord_impl(pt, space_in, gen_i, f2d_i, f3d_i, model);
}

template <typename PointType, typename InSpace>
void write_one_dependent_coord_impl(
    const PointType& pt, const InSpace& space_in, std::size_t& gen_i,
    std::size_t& f2d_i, std::size_t& f3d_i,
    const std::shared_ptr<kte::direct_kinematics_model>& model) {
  if constexpr (is_normal_joint_space_v<InSpace>) {
    *(model->getDependentCoord(gen_i++)->mFrame) = get_gen_coord(pt);
  } else if constexpr (is_Ndof_space_v<InSpace>) {
    for (std::size_t i = 0; i < get<0>(pt).size(); ++i) {
      *(model->getDependentCoord(gen_i++)->mFrame) = get_gen_coord(pt, i);
    }
  } else if constexpr (is_se2_space_v<InSpace>) {
    *(model->getDependentFrame2D(f2d_i++)->mFrame) = get_frame_2D(pt);
  } else if constexpr (is_se3_space_v<InSpace>) {
    *(model->getDependentFrame3D(f3d_i++)->mFrame) = get_frame_3D(pt);
  } else {
    tuple_for_each(pt, space_in,
                   [&](const auto& pt_elem, const auto& space_in_elem) {
                     read_one_joint_coord_impl(pt_elem, space_in_elem, gen_i,
                                               f2d_i, f3d_i, model);
                   });
  }
}

template <typename PointType, typename InSpaceTuple>
void write_dependent_coordinates_impl(
    const PointType& pt, const InSpaceTuple& space_in,
    const std::shared_ptr<kte::direct_kinematics_model>& model) {
  std::size_t gen_i = 0;
  std::size_t f2d_i = 0;
  std::size_t f3d_i = 0;
  write_one_dependent_coord_impl(pt, space_in, gen_i, f2d_i, f3d_i, model);
}

}  // namespace

}  // namespace ReaK::pp::detail

#endif

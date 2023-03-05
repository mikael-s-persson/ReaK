/**
 * \file direct_kinematics_topomap.hpp
 *
 * This library provides classes that define topological mappings from the joint-space
 * to the end-effector frame of a serial manipulator (or any open kinematics KTE chain,
 * see ReaK::kte::manipulator_kinematics_model).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2012
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

#ifndef REAK_DIRECT_KINEMATICS_TOPOMAP_HPP
#define REAK_DIRECT_KINEMATICS_TOPOMAP_HPP

#include "ReaK/core/base/defs.hpp"
#include "ReaK/core/base/shared_object.hpp"

#include "ReaK/mbd/models/direct_kinematics_model.hpp"

#include "ReaK/topologies/spaces/direct_kinematics_topomap_detail.hpp"

namespace ReaK::pp {

/**
 * This class implements the forward kinematics mappings associated to a given manipulator kinematics
 * model. This class assumes that the manipulator model has a number of joint coordinates (both
 * generalized and frames), and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
class manip_direct_kin_map : public shared_object {
 public:
  using self = manip_direct_kin_map;

  /** This data member points to a manipulator kinematics model to use for the mappings performed. */
  std::shared_ptr<kte::direct_kinematics_model> model;

  explicit manip_direct_kin_map(
      const std::shared_ptr<kte::direct_kinematics_model>& aModel =
          std::shared_ptr<kte::direct_kinematics_model>())
      : model(aModel) {}

  /**
   * This function template applies a forward kinematics calculation on the
   * manipulator model with the given joint-space state.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space (joint-space).
   * \param pt The point in the input space, i.e. the joint coordinates.
   * \param space_in The input space, i.e. the joint-space.
   */
  template <typename PointType, typename InSpace>
  void apply_to_model(const PointType& pt, const InSpace& space_in) const {
    detail::write_joint_coordinates_impl(pt, space_in, model);
    model->doDirectMotion();
  }

  /**
   * This function template performs a forward kinematics calculation on the
   * manipulator model.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space (joint-space).
   * \tparam OutSpace The type of the output space (end-effector space).
   * \param pt The point in the input space, i.e. the joint coordinates.
   * \param space_in The input space, i.e. the joint-space.
   * \param space_out The output space, i.e. the end-effector space.
   * \return A point in the output space, i.e. the end-effector coordinates.
   */
  template <typename PointType, typename InSpace, typename OutSpace>
  typename topology_traits<OutSpace>::point_type map_to_space(
      const PointType& pt, const InSpace& space_in,
      const OutSpace& space_out) const {

    apply_to_model(pt, space_in);

    typename topology_traits<OutSpace>::point_type result;
    detail::read_dependent_coordinates_impl(result, space_out, model);

    return result;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(model);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(model);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400012, 1, "manip_direct_kin_map",
                              shared_object)
};

/**
 * This class implements the forward kinematics mappings associated to a given manipulator kinematics
 * model. This class assumes that the manipulator model has a number of rate-limited joint coordinates
 * (both generalized and frames), and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 * \tparam RateLimitMap The type of the mapping between rate-limited joint-spaces and normal joint-spaces.
 */
template <typename RateLimitMap, typename NormalJointSpace>
class manip_rl_direct_kin_map : public shared_object {
 public:
  using self = manip_rl_direct_kin_map<RateLimitMap, NormalJointSpace>;

  /** This data member points to a manipulator kinematics model to use for the mappings performed. */
  std::shared_ptr<kte::direct_kinematics_model> model;
  /** This data member holds a mapping between the rate-limited joint space and the normal joint-space. */
  RateLimitMap joint_limits_map;
  std::shared_ptr<NormalJointSpace> normal_jt_space;

  explicit manip_rl_direct_kin_map(
      const std::shared_ptr<kte::direct_kinematics_model>& aModel =
          std::shared_ptr<kte::direct_kinematics_model>(),
      const RateLimitMap& aJointLimitMap = RateLimitMap(),
      const std::shared_ptr<NormalJointSpace>& aNormalJtSpace =
          std::shared_ptr<NormalJointSpace>())
      : model(aModel),
        joint_limits_map(aJointLimitMap),
        normal_jt_space(aNormalJtSpace) {}

  /**
   * This function template applies a forward kinematics calculation on the
   * manipulator model with the given joint-space state.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space (joint-space).
   * \param pt The point in the input space, i.e. the joint coordinates.
   * \param space_in The input space, i.e. the joint-space.
   */
  template <typename PointType, typename InSpace>
  void apply_to_model(const PointType& pt, const InSpace& space_in) const {
    typename topology_traits<NormalJointSpace>::point_type pt_inter =
        joint_limits_map.map_to_space(pt, space_in, *normal_jt_space);
    detail::write_joint_coordinates_impl(pt_inter, *normal_jt_space, model);
    model->doDirectMotion();
  }

  /**
   * This function template performs a forward kinematics calculation on the
   * manipulator model.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space (rate-limited joint-space).
   * \tparam OutSpace The type of the output space (end-effector space).
   * \param pt The point in the input space, i.e. the rate-limited joint coordinates.
   * \param space_in The input space, i.e. the rate-limited joint-space.
   * \param space_out The output space, i.e. the end-effector space.
   * \return A point in the output space, i.e. the end-effector coordinates.
   */
  template <typename PointType, typename InSpace, typename OutSpace>
  typename topology_traits<OutSpace>::point_type map_to_space(
      const PointType& pt, const InSpace& space_in,
      const OutSpace& space_out) const {

    apply_to_model(pt, space_in);

    typename topology_traits<OutSpace>::point_type result;
    detail::read_dependent_coordinates_impl(result, space_out, model);

    return result;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(model) &
        RK_SERIAL_SAVE_WITH_NAME(joint_limits_map) &
        RK_SERIAL_SAVE_WITH_NAME(normal_jt_space);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(model) &
        RK_SERIAL_LOAD_WITH_NAME(joint_limits_map) &
        RK_SERIAL_LOAD_WITH_NAME(normal_jt_space);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400013, 1, "manip_rl_direct_kin_map",
                              shared_object)
};

}  // namespace ReaK::pp

#endif

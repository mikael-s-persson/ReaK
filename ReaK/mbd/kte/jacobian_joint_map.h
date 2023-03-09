/**
 * \file jacobian_joint_map.h
 *
 * This library declares types for jacobian mappings of generalized coordinates. These jacobian frames are
 * associated to the motion in a generalized coordinate by the joints that have these generalized coordinates
 * as input.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date April 2010
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

#ifndef REAK_MBD_KTE_JACOBIAN_JOINT_MAP_H_
#define REAK_MBD_KTE_JACOBIAN_JOINT_MAP_H_

#include <map>
#include <utility>
#include "ReaK/core/base/defs.h"
#include "ReaK/math/kinetostatics/kinetostatics.h"
#include "ReaK/math/kinetostatics/motion_jacobians.h"

namespace ReaK::kte {

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian generalized coordinate. */
using jacobian_joint_map_gen =
    std::map<std::shared_ptr<gen_coord<double>>,
             std::shared_ptr<jacobian_gen_gen<double>>>;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 2D frame. */
using jacobian_joint_map_2D =
    std::map<std::shared_ptr<gen_coord<double>>,
             std::shared_ptr<jacobian_gen_2D<double>>>;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 3D frame. */
using jacobian_joint_map_3D =
    std::map<std::shared_ptr<gen_coord<double>>,
             std::shared_ptr<jacobian_gen_3D<double>>>;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian generalized coordinate. */
using jacobian_joint2D_map_gen =
    std::map<std::shared_ptr<frame_2D<double>>,
             std::shared_ptr<jacobian_2D_gen<double>>>;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 2D frame. */
using jacobian_joint2D_map_2D =
    std::map<std::shared_ptr<frame_2D<double>>,
             std::shared_ptr<jacobian_2D_2D<double>>>;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 3D frame. */
using jacobian_joint2D_map_3D =
    std::map<std::shared_ptr<frame_2D<double>>,
             std::shared_ptr<jacobian_2D_3D<double>>>;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian generalized coordinate. */
using jacobian_joint3D_map_gen =
    std::map<std::shared_ptr<frame_3D<double>>,
             std::shared_ptr<jacobian_3D_gen<double>>>;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 2D frame. */
using jacobian_joint3D_map_2D =
    std::map<std::shared_ptr<frame_3D<double>>,
             std::shared_ptr<jacobian_3D_2D<double>>>;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 3D frame. */
using jacobian_joint3D_map_3D =
    std::map<std::shared_ptr<frame_3D<double>>,
             std::shared_ptr<jacobian_3D_3D<double>>>;

class joint_dependent_gen_coord : public shared_object {
 public:
  std::shared_ptr<gen_coord<double>>
      mFrame;  ///< Holds the generalized coordinate.

  jacobian_joint_map_gen
      mUpStreamJoints;  ///< Holds the jacobian mappings of up-stream joints.
  jacobian_joint2D_map_gen
      mUpStream2DJoints;  ///< Holds the jacobian mappings of up-stream 2D joints.
  jacobian_joint3D_map_gen
      mUpStream3DJoints;  ///< Holds the jacobian mappings of up-stream 3D joints.

  /**
   * Parametrized and Default constructor.
   * \param aFrame The generalized coordinate.
   * \param aUpStreamJoints The jacobian mappings of up-stream joints.
   * \param aUpStream2DJoints The jacobian mappings of up-stream 2D joints.
   * \param aUpStream3DJoints The jacobian mappings of up-stream 3D joints.
   */
  explicit joint_dependent_gen_coord(
      std::shared_ptr<gen_coord<double>> aFrame =
          std::shared_ptr<gen_coord<double>>(),
      jacobian_joint_map_gen aUpStreamJoints = jacobian_joint_map_gen(),
      jacobian_joint2D_map_gen aUpStream2DJoints = jacobian_joint2D_map_gen(),
      jacobian_joint3D_map_gen aUpStream3DJoints = jacobian_joint3D_map_gen())
      : mFrame(std::move(aFrame)),
        mUpStreamJoints(std::move(aUpStreamJoints)),
        mUpStream2DJoints(std::move(aUpStream2DJoints)),
        mUpStream3DJoints(std::move(aUpStream3DJoints)) {}

  /**
   * Default destructor.
   */
  ~joint_dependent_gen_coord() override = default;

  joint_dependent_gen_coord& add_joint(
      const std::shared_ptr<gen_coord<double>>& aJointFrame,
      const std::shared_ptr<jacobian_gen_gen<double>>& aJointJacobian) {
    mUpStreamJoints[aJointFrame] = aJointJacobian;
    return *this;
  }

  joint_dependent_gen_coord& add_joint(
      const std::shared_ptr<frame_2D<double>>& aJointFrame,
      const std::shared_ptr<jacobian_2D_gen<double>>& aJointJacobian) {
    mUpStream2DJoints[aJointFrame] = aJointJacobian;
    return *this;
  }

  joint_dependent_gen_coord& add_joint(
      const std::shared_ptr<frame_3D<double>>& aJointFrame,
      const std::shared_ptr<jacobian_3D_gen<double>>& aJointJacobian) {
    mUpStream3DJoints[aJointFrame] = aJointJacobian;
    return *this;
  }

  joint_dependent_gen_coord& remove_joint(
      const std::shared_ptr<gen_coord<double>>& aJointFrame) {
    mUpStreamJoints.erase(mUpStreamJoints.find(aJointFrame));
    return *this;
  }

  joint_dependent_gen_coord& remove_joint(
      const std::shared_ptr<frame_2D<double>>& aJointFrame) {
    mUpStream2DJoints.erase(mUpStream2DJoints.find(aJointFrame));
    return *this;
  }

  joint_dependent_gen_coord& remove_joint(
      const std::shared_ptr<frame_3D<double>>& aJointFrame) {
    mUpStream3DJoints.erase(mUpStream3DJoints.find(aJointFrame));
    return *this;
  }

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(mFrame) &
        RK_SERIAL_SAVE_WITH_NAME(mUpStreamJoints) &
        RK_SERIAL_SAVE_WITH_NAME(mUpStream2DJoints) &
        RK_SERIAL_SAVE_WITH_NAME(mUpStream3DJoints);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(mFrame) &
        RK_SERIAL_LOAD_WITH_NAME(mUpStreamJoints) &
        RK_SERIAL_LOAD_WITH_NAME(mUpStream2DJoints) &
        RK_SERIAL_LOAD_WITH_NAME(mUpStream3DJoints);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(joint_dependent_gen_coord, 0xC2000002, 1,
                              "joint_dependent_gen_coord", shared_object)
};

class joint_dependent_frame_2D : public shared_object {
 public:
  std::shared_ptr<frame_2D<double>>
      mFrame;  ///< Holds the generalized coordinate.

  jacobian_joint_map_2D
      mUpStreamJoints;  ///< Holds the jacobian mappings of up-stream joints.
  jacobian_joint2D_map_2D
      mUpStream2DJoints;  ///< Holds the jacobian mappings of up-stream 2D joints.
  jacobian_joint3D_map_2D
      mUpStream3DJoints;  ///< Holds the jacobian mappings of up-stream 3D joints.

  /**
   * Parametrized and Default constructor.
   * \param aFrame The generalized coordinate.
   * \param aUpStreamJoints The jacobian mappings of up-stream joints.
   * \param aUpStream2DJoints The jacobian mappings of up-stream 2D joints.
   * \param aUpStream3DJoints The jacobian mappings of up-stream 3D joints.
   */
  explicit joint_dependent_frame_2D(
      std::shared_ptr<frame_2D<double>> aFrame =
          std::shared_ptr<frame_2D<double>>(),
      jacobian_joint_map_2D aUpStreamJoints = jacobian_joint_map_2D(),
      jacobian_joint2D_map_2D aUpStream2DJoints = jacobian_joint2D_map_2D(),
      jacobian_joint3D_map_2D aUpStream3DJoints = jacobian_joint3D_map_2D())
      : mFrame(std::move(aFrame)),
        mUpStreamJoints(std::move(aUpStreamJoints)),
        mUpStream2DJoints(std::move(aUpStream2DJoints)),
        mUpStream3DJoints(std::move(aUpStream3DJoints)) {}

  /**
   * Default destructor.
   */
  ~joint_dependent_frame_2D() override = default;

  joint_dependent_frame_2D& add_joint(
      const std::shared_ptr<gen_coord<double>>& aJointFrame,
      const std::shared_ptr<jacobian_gen_2D<double>>& aJointJacobian) {
    mUpStreamJoints[aJointFrame] = aJointJacobian;
    return *this;
  }

  joint_dependent_frame_2D& add_joint(
      const std::shared_ptr<frame_2D<double>>& aJointFrame,
      const std::shared_ptr<jacobian_2D_2D<double>>& aJointJacobian) {
    mUpStream2DJoints[aJointFrame] = aJointJacobian;
    return *this;
  }

  joint_dependent_frame_2D& add_joint(
      const std::shared_ptr<frame_3D<double>>& aJointFrame,
      const std::shared_ptr<jacobian_3D_2D<double>>& aJointJacobian) {
    mUpStream3DJoints[aJointFrame] = aJointJacobian;
    return *this;
  }

  joint_dependent_frame_2D& remove_joint(
      const std::shared_ptr<gen_coord<double>>& aJointFrame) {
    mUpStreamJoints.erase(mUpStreamJoints.find(aJointFrame));
    return *this;
  }

  joint_dependent_frame_2D& remove_joint(
      const std::shared_ptr<frame_2D<double>>& aJointFrame) {
    mUpStream2DJoints.erase(mUpStream2DJoints.find(aJointFrame));
    return *this;
  }

  joint_dependent_frame_2D& remove_joint(
      const std::shared_ptr<frame_3D<double>>& aJointFrame) {
    mUpStream3DJoints.erase(mUpStream3DJoints.find(aJointFrame));
    return *this;
  }

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(mFrame) &
        RK_SERIAL_SAVE_WITH_NAME(mUpStreamJoints) &
        RK_SERIAL_SAVE_WITH_NAME(mUpStream2DJoints) &
        RK_SERIAL_SAVE_WITH_NAME(mUpStream3DJoints);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(mFrame) &
        RK_SERIAL_LOAD_WITH_NAME(mUpStreamJoints) &
        RK_SERIAL_LOAD_WITH_NAME(mUpStream2DJoints) &
        RK_SERIAL_LOAD_WITH_NAME(mUpStream3DJoints);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(joint_dependent_frame_2D, 0xC2000003, 1,
                              "joint_dependent_frame_2D", shared_object)
};

class joint_dependent_frame_3D : public shared_object {
 public:
  std::shared_ptr<frame_3D<double>>
      mFrame;  ///< Holds the generalized coordinate.

  jacobian_joint_map_3D
      mUpStreamJoints;  ///< Holds the jacobian mappings of up-stream joints.
  jacobian_joint2D_map_3D
      mUpStream2DJoints;  ///< Holds the jacobian mappings of up-stream 2D joints.
  jacobian_joint3D_map_3D
      mUpStream3DJoints;  ///< Holds the jacobian mappings of up-stream 3D joints.

  /**
   * Parametrized and Default constructor.
   * \param aFrame The 3D frame.
   * \param aUpStreamJoints The jacobian mappings of up-stream joints.
   * \param aUpStream2DJoints The jacobian mappings of up-stream 2D joints.
   * \param aUpStream3DJoints The jacobian mappings of up-stream 3D joints.
   */
  explicit joint_dependent_frame_3D(
      std::shared_ptr<frame_3D<double>> aFrame =
          std::shared_ptr<frame_3D<double>>(),
      jacobian_joint_map_3D aUpStreamJoints = jacobian_joint_map_3D(),
      jacobian_joint2D_map_3D aUpStream2DJoints = jacobian_joint2D_map_3D(),
      jacobian_joint3D_map_3D aUpStream3DJoints = jacobian_joint3D_map_3D())
      : mFrame(std::move(aFrame)),
        mUpStreamJoints(std::move(aUpStreamJoints)),
        mUpStream2DJoints(std::move(aUpStream2DJoints)),
        mUpStream3DJoints(std::move(aUpStream3DJoints)) {}

  /**
   * Default destructor.
   */
  ~joint_dependent_frame_3D() override = default;

  joint_dependent_frame_3D& add_joint(
      const std::shared_ptr<gen_coord<double>>& aJointFrame,
      const std::shared_ptr<jacobian_gen_3D<double>>& aJointJacobian) {
    mUpStreamJoints[aJointFrame] = aJointJacobian;
    return *this;
  }

  joint_dependent_frame_3D& add_joint(
      const std::shared_ptr<frame_2D<double>>& aJointFrame,
      const std::shared_ptr<jacobian_2D_3D<double>>& aJointJacobian) {
    mUpStream2DJoints[aJointFrame] = aJointJacobian;
    return *this;
  }

  joint_dependent_frame_3D& add_joint(
      const std::shared_ptr<frame_3D<double>>& aJointFrame,
      const std::shared_ptr<jacobian_3D_3D<double>>& aJointJacobian) {
    mUpStream3DJoints[aJointFrame] = aJointJacobian;
    return *this;
  }

  joint_dependent_frame_3D& remove_joint(
      const std::shared_ptr<gen_coord<double>>& aJointFrame) {
    mUpStreamJoints.erase(mUpStreamJoints.find(aJointFrame));
    return *this;
  }

  joint_dependent_frame_3D& remove_joint(
      const std::shared_ptr<frame_2D<double>>& aJointFrame) {
    mUpStream2DJoints.erase(mUpStream2DJoints.find(aJointFrame));
    return *this;
  }

  joint_dependent_frame_3D& remove_joint(
      const std::shared_ptr<frame_3D<double>>& aJointFrame) {
    mUpStream3DJoints.erase(mUpStream3DJoints.find(aJointFrame));
    return *this;
  }

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(mFrame) &
        RK_SERIAL_SAVE_WITH_NAME(mUpStreamJoints) &
        RK_SERIAL_SAVE_WITH_NAME(mUpStream2DJoints) &
        RK_SERIAL_SAVE_WITH_NAME(mUpStream3DJoints);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(mFrame) &
        RK_SERIAL_LOAD_WITH_NAME(mUpStreamJoints) &
        RK_SERIAL_LOAD_WITH_NAME(mUpStream2DJoints) &
        RK_SERIAL_LOAD_WITH_NAME(mUpStream3DJoints);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(joint_dependent_frame_3D, 0xC2000004, 1,
                              "joint_dependent_frame_3D", shared_object)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_JACOBIAN_JOINT_MAP_H_

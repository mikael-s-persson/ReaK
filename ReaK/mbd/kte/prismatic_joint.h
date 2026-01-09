/**
 * \file prismatic_joint.h
 *
 * This library declares the KTE models for prismatic joints in 2D and 3D space. These
 * models implement a model of a single degree-of-freedom linear joint (i.e. prismatic joint),
 * allowing no rotation and only displacement along its predefined, fixed axis.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2010
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

#ifndef REAK_MBD_KTE_PRISMATIC_JOINT_H_
#define REAK_MBD_KTE_PRISMATIC_JOINT_H_

#include <utility>
#include "ReaK/math/kinetostatics/motion_jacobians.h"

#include "ReaK/mbd/kte/reacting_kte.h"

namespace ReaK::kte {

/**
 * This class implements a prismatic joint in 2D space. A generalized coordinate is used to represent the
 * joint's linear displacement between a base coordinate frame to an end coordinate frame.
 */
class prismatic_joint_2D : public reacting_kte_gen {
 protected:
  std::shared_ptr<gen_coord<double>>
      mCoord;  ///< Generalized coordinate representing the joint's linear displacement.
  vect<double, 2>
      mAxis;  ///< Joint axis, as a fixed vector, in the coordinate system of the base frame.
  std::shared_ptr<frame_2D<double>>
      mBase;  ///< The coordinate frame at the base of the joint.
  std::shared_ptr<frame_2D<double>>
      mEnd;  ///< The coordinate frame just after the joint transformations are applied.

  std::shared_ptr<jacobian_gen_2D<double>>
      mJacobian;  ///< The Jacobian frame produced by this joint.

 public:
  /**
   * Sets the joint's coordinate.
   * \param aPtr The new joint's coordinate.
   */
  void setCoord(const std::shared_ptr<gen_coord<double>>& aPtr) {
    mCoord = aPtr;
  }
  /**
   * Returns the joint's coordinate.
   * \return The joint's coordinate.
   */
  std::shared_ptr<gen_coord<double>> Coord() const { return mCoord; }

  /**
   * Sets the joint's axis (relative to base frame).
   * \param aValue The new joint's axis (relative to base frame).
   */
  void setAxis(const vect<double, 2>& aValue) { mAxis = aValue; }
  /**
   * Returns the joint's axis (relative to base frame).
   * \return The joint's axis (relative to base frame).
   */
  vect<double, 2> Axis() const { return mAxis; }

  /**
   * Sets the joint's base frame.
   * \param aPtr The new joint's base frame.
   */
  void setBaseFrame(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mBase = aPtr;
  }
  /**
   * Returns the joint's base frame.
   * \return The joint's base frame.
   */
  std::shared_ptr<frame_2D<double>> BaseFrame() const { return mBase; }

  /**
   * Sets the joint's output frame.
   * \param aPtr The new joint's output frame.
   */
  void setEndFrame(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mEnd = aPtr;
  }
  /**
   * Returns the joint's output frame.
   * \return The joint's output frame.
   */
  std::shared_ptr<frame_2D<double>> EndFrame() const { return mEnd; }

  /**
   * Sets the joint's Jacobian.
   * \param aPtr The new joint's Jacobian.
   */
  void setJacobian(const std::shared_ptr<jacobian_gen_2D<double>>& aPtr) {
    mJacobian = aPtr;
  }
  /**
   * Returns the joint's Jacobian.
   * \return The joint's Jacobian.
   */
  std::shared_ptr<jacobian_gen_2D<double>> Jacobian() const {
    return mJacobian;
  }

  /**
   * Default constructor.
   */
  explicit prismatic_joint_2D(const std::string& aName = "")
      : reacting_kte_gen(aName) {}

  /**
   * Parametrized constructor.
   * \param aName the name of this KTE model.
   * \param aCoord the generalized coordinate associated with the displacement of this joint.
   * \param aAxis the joint axis along which linear displacement occurs, represented in base frame coordinates.
   * \param aBase the coordinate frame at the base of the joint.
   * \param aEnd the coordinate frame just after the joint transformations are applied.
   * \param aJacobian a pointer to contain the Jacobian frame produced by this joint, default value will disable the
   * Jacobian frame's calculation.
   */
  prismatic_joint_2D(const std::string& aName,
                     std::shared_ptr<gen_coord<double>> aCoord,
                     const vect<double, 2>& aAxis,
                     std::shared_ptr<frame_2D<double>> aBase,
                     std::shared_ptr<frame_2D<double>> aEnd,
                     std::shared_ptr<jacobian_gen_2D<double>> aJacobian =
                         std::shared_ptr<jacobian_gen_2D<double>>())
      : reacting_kte_gen(aName),
        mCoord(std::move(aCoord)),
        mAxis(aAxis),
        mBase(std::move(aBase)),
        mEnd(std::move(aEnd)),
        mJacobian(std::move(aJacobian)) {}

  /**
   * Default destructor.
   */
  ~prismatic_joint_2D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override;

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override;

  void clearForce() override;

  void applyReactionForce(double aForce) override;

  void save(serialization::oarchive& A, unsigned int Version) const override {
    reacting_kte_gen::save(
        A, reacting_kte_gen::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(mCoord) & RK_SERIAL_SAVE_WITH_NAME(mAxis) &
        RK_SERIAL_SAVE_WITH_NAME(mBase) & RK_SERIAL_SAVE_WITH_NAME(mEnd) &
        RK_SERIAL_SAVE_WITH_NAME(mJacobian);
  }

  void load(serialization::iarchive& A, unsigned int Version) override {
    reacting_kte_gen::load(
        A, reacting_kte_gen::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(mCoord) & RK_SERIAL_LOAD_WITH_NAME(mAxis) &
        RK_SERIAL_LOAD_WITH_NAME(mBase) & RK_SERIAL_LOAD_WITH_NAME(mEnd) &
        RK_SERIAL_LOAD_WITH_NAME(mJacobian);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(prismatic_joint_2D, 0xC2100005, 1,
                              "prismatic_joint_2D", reacting_kte_gen)
};

/**
 * This class implements a prismatic joint in 3D space. A generalized coordinate is used to represent the
 * joint's linear displacement between a base coordinate frame to an end coordinate frame.
 */
class prismatic_joint_3D : public reacting_kte_gen {
 protected:
  std::shared_ptr<gen_coord<double>>
      mCoord;  ///< Generalized coordinate representing the joint's linear displacement.
  vect<double, 3>
      mAxis;  ///< Joint axis, as a fixed vector, in the coordinate system of the base frame.
  std::shared_ptr<frame_3D<double>>
      mBase;  ///< The coordinate frame at the base of the joint.
  std::shared_ptr<frame_3D<double>>
      mEnd;  ///< The coordinate frame just after the joint transformations are applied.

  std::shared_ptr<jacobian_gen_3D<double>>
      mJacobian;  ///< The Jacobian frame produced by this joint.

 public:
  /**
   * Sets the joint's coordinate.
   * \param aPtr The new joint's coordinate.
   */
  void setCoord(const std::shared_ptr<gen_coord<double>>& aPtr) {
    mCoord = aPtr;
  }
  /**
   * Returns the joint's coordinate.
   * \return The joint's coordinate.
   */
  std::shared_ptr<gen_coord<double>> Coord() const { return mCoord; }

  /**
   * Sets the joint's axis (relative to base frame).
   * \param aValue The new joint's axis (relative to base frame).
   */
  void setAxis(const vect<double, 3>& aValue) { mAxis = aValue; }
  /**
   * Returns the joint's axis (relative to base frame).
   * \return The joint's axis (relative to base frame).
   */
  vect<double, 3> Axis() const { return mAxis; }

  /**
   * Sets the joint's base frame.
   * \param aPtr The new joint's base frame.
   */
  void setBaseFrame(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mBase = aPtr;
  }
  /**
   * Returns the joint's base frame.
   * \return The joint's base frame.
   */
  std::shared_ptr<frame_3D<double>> BaseFrame() const { return mBase; }

  /**
   * Sets the joint's output frame.
   * \param aPtr The new joint's output frame.
   */
  void setEndFrame(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mEnd = aPtr;
  }
  /**
   * Returns the joint's output frame.
   * \return The joint's output frame.
   */
  std::shared_ptr<frame_3D<double>> EndFrame() const { return mEnd; }

  /**
   * Sets the joint's Jacobian.
   * \param aPtr The new joint's Jacobian.
   */
  void setJacobian(const std::shared_ptr<jacobian_gen_3D<double>>& aPtr) {
    mJacobian = aPtr;
  }
  /**
   * Returns the joint's Jacobian.
   * \return The joint's Jacobian.
   */
  std::shared_ptr<jacobian_gen_3D<double>> Jacobian() const {
    return mJacobian;
  }

  /**
   * Default constructor.
   */
  explicit prismatic_joint_3D(const std::string& aName = "")
      : reacting_kte_gen(aName) {}

  /**
   * Parametrized constructor.
   * \param aName the name of this KTE model.
   * \param aCoord the generalized coordinate associated with the displacement of this joint.
   * \param aAxis the joint axis along which linear displacement occurs, represented in base frame coordinates.
   * \param aBase the coordinate frame at the base of the joint.
   * \param aEnd the coordinate frame just after the joint transformations are applied.
   * \param aJacobian a pointer to contain the Jacobian frame produced by this joint, default value will disable the
   * Jacobian frame's calculation.
   */
  prismatic_joint_3D(const std::string& aName,
                     std::shared_ptr<gen_coord<double>> aCoord,
                     const vect<double, 3>& aAxis,
                     std::shared_ptr<frame_3D<double>> aBase,
                     std::shared_ptr<frame_3D<double>> aEnd,
                     std::shared_ptr<jacobian_gen_3D<double>> aJacobian =
                         std::shared_ptr<jacobian_gen_3D<double>>())
      : reacting_kte_gen(aName),
        mCoord(std::move(aCoord)),
        mAxis(aAxis),
        mBase(std::move(aBase)),
        mEnd(std::move(aEnd)),
        mJacobian(std::move(aJacobian)) {}

  /**
   * Default destructor.
   */
  ~prismatic_joint_3D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override;

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override;

  void clearForce() override;

  void applyReactionForce(double aForce) override;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    reacting_kte_gen::save(
        A, reacting_kte_gen::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(mCoord) & RK_SERIAL_SAVE_WITH_NAME(mAxis) &
        RK_SERIAL_SAVE_WITH_NAME(mBase) & RK_SERIAL_SAVE_WITH_NAME(mEnd) &
        RK_SERIAL_SAVE_WITH_NAME(mJacobian);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    reacting_kte_gen::load(
        A, reacting_kte_gen::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(mCoord) & RK_SERIAL_LOAD_WITH_NAME(mAxis) &
        RK_SERIAL_LOAD_WITH_NAME(mBase) & RK_SERIAL_LOAD_WITH_NAME(mEnd) &
        RK_SERIAL_LOAD_WITH_NAME(mJacobian);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(prismatic_joint_3D, 0xC2100006, 1,
                              "prismatic_joint_3D", reacting_kte_gen)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_PRISMATIC_JOINT_H_

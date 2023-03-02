/**
 * \file free_joints.hpp
 *
 * This library provides classes to handle free-joints. A free joint is simple a
 * frame transformation that is applied to an input frame to obtain the output
 * frame. This way, one can represent a general motion of a body.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2011
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

#ifndef REAK_FREE_JOINTS_HPP
#define REAK_FREE_JOINTS_HPP

#include <ReaK/math/kinetostatics/motion_jacobians.hpp>
#include <utility>
#include "reacting_kte.hpp"

namespace ReaK::kte {

/**
 * This class implements a free joint in 2D space. A 2D frame is used to represent the
 * joint's motion between a base coordinate frame to an end coordinate frame.
 */
class free_joint_2D : public reacting_kte_2D {
 protected:
  std::shared_ptr<frame_2D<double>>
      mCoord;  ///< The coordinate frame representing the joint's motion.
  std::shared_ptr<frame_2D<double>>
      mBase;  ///< The coordinate frame at the base of the joint.
  std::shared_ptr<frame_2D<double>>
      mEnd;  ///< The coordinate frame just after the joint transformations are applied.

  std::shared_ptr<jacobian_2D_2D<double>>
      mJacobian;  ///< The Jacobian frame produced by this joint.

 public:
  /**
   * Sets the joint's space coordinate.
   * \param aPtr The new joint's space coordinate.
   */
  void setCoord(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mCoord = aPtr;
  }
  /**
   * Returns the joint's space coordinate.
   * \return The joint's space coordinate.
   */
  std::shared_ptr<frame_2D<double>> Coord() const { return mCoord; }

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
  void setJacobian(const std::shared_ptr<jacobian_2D_2D<double>>& aPtr) {
    mJacobian = aPtr;
  }
  /**
   * Returns the joint's Jacobian.
   * \return The joint's Jacobian.
   */
  std::shared_ptr<jacobian_2D_2D<double>> Jacobian() const { return mJacobian; }

  /**
   * Default constructor.
   */
  explicit free_joint_2D(const std::string& aName = "")
      : reacting_kte_2D(aName) {}

  /**
   * Parametrized constructor.
   * \param aName the name of this KTE model.
   * \param aCoord the coordinate frame associated with the displacement of this joint.
   * \param aBase the coordinate frame at the base of the joint.
   * \param aEnd the coordinate frame just after the joint transformations are applied.
   * \param aJacobian a pointer to contain the Jacobian frame produced by this joint, default value will disable the
   * Jacobian frame's calculation.
   */
  free_joint_2D(const std::string& aName,
                std::shared_ptr<frame_2D<double>> aCoord,
                std::shared_ptr<frame_2D<double>> aBase,
                std::shared_ptr<frame_2D<double>> aEnd,
                std::shared_ptr<jacobian_2D_2D<double>> aJacobian =
                    std::shared_ptr<jacobian_2D_2D<double>>())
      : reacting_kte_2D(aName),
        mCoord(std::move(aCoord)),
        mBase(std::move(aBase)),
        mEnd(std::move(aEnd)),
        mJacobian(std::move(aJacobian)) {}

  /**
   * Default destructor.
   */
  ~free_joint_2D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override;

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override;

  void clearForce() override;

  void applyReactionForce(vect<double, 2> aForce, double aTorque) override;

  void save(serialization::oarchive& A, unsigned int Version) const override {
    reacting_kte_2D::save(
        A, reacting_kte_2D::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mCoord) & RK_SERIAL_SAVE_WITH_NAME(mBase) &
        RK_SERIAL_SAVE_WITH_NAME(mEnd) & RK_SERIAL_SAVE_WITH_NAME(mJacobian);
  }

  void load(serialization::iarchive& A, unsigned int Version) override {
    reacting_kte_2D::load(
        A, reacting_kte_2D::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mCoord) & RK_SERIAL_LOAD_WITH_NAME(mBase) &
        RK_SERIAL_LOAD_WITH_NAME(mEnd) & RK_SERIAL_LOAD_WITH_NAME(mJacobian);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(free_joint_2D, 0xC2100041, 1, "free_joint_2D",
                              reacting_kte_2D)
};

/**
 * This class implements a free joint in 3D space. A 3D frame is used to represent the
 * joint's motion between a base coordinate frame to an end coordinate frame.
 */
class free_joint_3D : public reacting_kte_3D {
 protected:
  std::shared_ptr<frame_3D<double>>
      mCoord;  ///< The coordinate frame representing the joint's motion.
  std::shared_ptr<frame_3D<double>>
      mBase;  ///< The coordinate frame at the base of the joint.
  std::shared_ptr<frame_3D<double>>
      mEnd;  ///< The coordinate frame just after the joint transformations are applied.

  std::shared_ptr<jacobian_3D_3D<double>>
      mJacobian;  ///< The Jacobian frame produced by this joint.

 public:
  /**
   * Sets the joint's space coordinate.
   * \param aPtr The new joint's space coordinate.
   */
  void setCoord(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mCoord = aPtr;
  }
  /**
   * Returns the joint's space coordinate.
   * \return The joint's space coordinate.
   */
  std::shared_ptr<frame_3D<double>> Coord() const { return mCoord; }

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
  void setJacobian(const std::shared_ptr<jacobian_3D_3D<double>>& aPtr) {
    mJacobian = aPtr;
  }
  /**
   * Returns the joint's Jacobian.
   * \return The joint's Jacobian.
   */
  std::shared_ptr<jacobian_3D_3D<double>> Jacobian() const { return mJacobian; }

  /**
   * Default constructor.
   */
  explicit free_joint_3D(const std::string& aName = "")
      : reacting_kte_3D(aName) {}

  /**
   * Parametrized constructor.
   * \param aName the name of this KTE model.
   * \param aCoord the coordinate frame associated with the motion of this joint.
   * \param aBase the coordinate frame at the base of the joint.
   * \param aEnd the coordinate frame just after the joint transformations are applied.
   * \param aJacobian a pointer to contain the Jacobian frame produced by this joint, default value will disable the
   * Jacobian frame's calculation.
   */
  free_joint_3D(const std::string& aName,
                std::shared_ptr<frame_3D<double>> aCoord,
                std::shared_ptr<frame_3D<double>> aBase,
                std::shared_ptr<frame_3D<double>> aEnd,
                std::shared_ptr<jacobian_3D_3D<double>> aJacobian =
                    std::shared_ptr<jacobian_3D_3D<double>>())
      : reacting_kte_3D(aName),
        mCoord(std::move(aCoord)),
        mBase(std::move(aBase)),
        mEnd(std::move(aEnd)),
        mJacobian(std::move(aJacobian)) {}

  /**
   * Default destructor.
   */
  ~free_joint_3D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override;

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override;

  void clearForce() override;

  void applyReactionForce(vect<double, 3> aForce,
                          vect<double, 3> aTorque) override;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    reacting_kte_3D::save(
        A, reacting_kte_3D::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mCoord) & RK_SERIAL_SAVE_WITH_NAME(mBase) &
        RK_SERIAL_SAVE_WITH_NAME(mEnd) & RK_SERIAL_SAVE_WITH_NAME(mJacobian);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    reacting_kte_3D::load(
        A, reacting_kte_3D::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mCoord) & RK_SERIAL_LOAD_WITH_NAME(mBase) &
        RK_SERIAL_LOAD_WITH_NAME(mEnd) & RK_SERIAL_LOAD_WITH_NAME(mJacobian);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(free_joint_3D, 0xC2100042, 1, "free_joint_3D",
                              reacting_kte_3D)
};

}  // namespace ReaK::kte

#endif

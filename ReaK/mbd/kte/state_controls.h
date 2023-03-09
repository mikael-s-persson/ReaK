/**
 * \file state_controls.h
 *
 * This library provides a number of classes which can be used as KTE models which
 * control the motion-state variables (position, velocity, etc.) directly. Note that
 * these classes don't implement control systems, simply direct control. If a control
 * system is required, one can make a KTE for that controller and insert it between
 * a frame that is controlled via the classes provided here and the input frame of a
 * joint.
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

#ifndef REAK_MBD_KTE_STATE_CONTROLS_H_
#define REAK_MBD_KTE_STATE_CONTROLS_H_

#include <utility>
#include "ReaK/math/kinetostatics/kinetostatics.h"
#include "ReaK/mbd/kte/kte_map.h"
#include "ReaK/mbd/kte/kte_system_input.h"

namespace ReaK::kte {

/**
 * This class can be used as a system input to set the value of the position of a
 * generalized coordinate.
 */
class position_control_gen : public kte_map, public system_input {
 private:
  std::shared_ptr<gen_coord<double>> mAnchor;
  double mPosDesired;

 public:
  /**
   * Sets the frame on which position control is exerted.
   * \param aPtr The frame on which position control is exerted.
   */
  void setAnchor(const std::shared_ptr<gen_coord<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which position control is exerted.
   * \return The frame on which position control is exerted.
   */
  std::shared_ptr<gen_coord<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the desired position.
   * \param aValue The new desired position.
   */
  void setPosDesired(double aValue) { mPosDesired = aValue; }
  /**
   * Returns the desired position, for read-only access.
   * \return the desired position for read-only access.
   */
  double PosDesired() const { return mPosDesired; }

  unsigned int getInputCount() const override { return 1; }
  void setInput(unsigned int i, double aValue) override {
    mPosDesired = aValue;
  }
  double getInput(unsigned int i) const override { return mPosDesired; }

  /**
   * Default constructor.
   */
  explicit position_control_gen(const std::string& aName = "")
      : kte_map(aName), mPosDesired(0.0) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  position_control_gen(const std::string& aName,
                       std::shared_ptr<ReaK::gen_coord<double>> aAnchor)
      : kte_map(aName),
        system_input(aName),
        mAnchor(std::move(aAnchor)),
        mPosDesired(0.0) {}

  /**
   * Default destructor.
   */
  ~position_control_gen() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mAnchor->q = mPosDesired;
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::save(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::load(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(position_control_gen, 0xC2100043, 1,
                              "position_control_gen", kte_map, system_input)
};

/**
 * This class can be used as a system input to set the value of the position of a
 * 2D coordinate frame.
 */
class position_control_2D : public kte_map, public system_input {
 private:
  std::shared_ptr<frame_2D<double>> mAnchor;
  vect<double, 2> mPosDesired;

 public:
  /**
   * Sets the frame on which position control is exerted.
   * \param aPtr The frame on which position control is exerted.
   */
  void setAnchor(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which position control is exerted.
   * \return The frame on which position control is exerted.
   */
  std::shared_ptr<frame_2D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the desired position.
   * \param aValue The new desired position.
   */
  void setPosDesired(const vect<double, 2>& aValue) { mPosDesired = aValue; }
  /**
   * Returns the desired position.
   * \return the desired position.
   */
  vect<double, 2> PosDesired() const { return mPosDesired; }

  unsigned int getInputCount() const override { return 2; }
  void setInput(unsigned int i, double aValue) override {
    if (i < 2) {
      mPosDesired[i] = aValue;
    } else {
      mPosDesired[0] = aValue;
    }
  }
  double getInput(unsigned int i) const override {
    if (i < 2) {
      return mPosDesired[i];
    }
    return mPosDesired[0];
  }

  /**
   * Default constructor.
   */
  explicit position_control_2D(const std::string& aName = "")
      : kte_map(aName), system_input(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  position_control_2D(const std::string& aName,
                      std::shared_ptr<ReaK::frame_2D<double>> aAnchor)
      : kte_map(aName), system_input(aName), mAnchor(std::move(aAnchor)) {}

  /**
   * Default destructor.
   */
  ~position_control_2D() override = default;
  ;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mAnchor->Position = mPosDesired;
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::save(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::load(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(position_control_2D, 0xC2100044, 1,
                              "position_control_2D", kte_map, system_input)
};

/**
 * This class can be used as a system input to set the value of the position of a
 * 3D coordinate frame.
 */
class position_control_3D : public kte_map, public system_input {
 private:
  std::shared_ptr<frame_3D<double>> mAnchor;
  vect<double, 3> mPosDesired;

 public:
  /**
   * Sets the frame on which position control is exerted.
   * \param aPtr The frame on which position control is exerted.
   */
  void setAnchor(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which position control is exerted.
   * \return The frame on which position control is exerted.
   */
  std::shared_ptr<frame_3D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the desired position.
   * \param aValue The new desired position.
   */
  void setPosDesired(const vect<double, 3>& aValue) { mPosDesired = aValue; }
  /**
   * Returns the desired position.
   * \return the desired position.
   */
  vect<double, 3> PosDesired() const { return mPosDesired; }

  unsigned int getInputCount() const override { return 3; }
  void setInput(unsigned int i, double aValue) override {
    if (i < 3) {
      mPosDesired[i] = aValue;
    } else {
      mPosDesired[0] = aValue;
    }
  }
  double getInput(unsigned int i) const override {
    if (i < 3) {
      return mPosDesired[i];
    }
    return mPosDesired[0];
  }

  /**
   * Default constructor.
   */
  explicit position_control_3D(const std::string& aName = "")
      : kte_map(aName), system_input(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  position_control_3D(const std::string& aName,
                      std::shared_ptr<ReaK::frame_3D<double>> aAnchor)
      : kte_map(aName), system_input(aName), mAnchor(std::move(aAnchor)) {}

  /**
   * Default destructor.
   */
  ~position_control_3D() override = default;
  ;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mAnchor->Position = mPosDesired;
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::save(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::load(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(position_control_3D, 0xC2100045, 1,
                              "position_control_3D", kte_map, system_input)
};

/**
 * This class can be used as a system input to set the value of the rotation of a
 * 2D coordinate frame.
 */
class rotation_control_2D : public kte_map, public system_input {
 private:
  std::shared_ptr<frame_2D<double>> mAnchor;
  double mAngleDesired;

 public:
  /**
   * Sets the frame on which rotation control is exerted.
   * \param aPtr The frame on which rotation control is exerted.
   */
  void setAnchor(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which rotation control is exerted.
   * \return The frame on which rotation control is exerted.
   */
  std::shared_ptr<frame_2D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the desired rotation.
   * \param aValue The new desired rotation.
   */
  void setAngleDesired(double aValue) { mAngleDesired = aValue; }
  /**
   * Returns the desired rotation.
   * \return the desired rotation.
   */
  double AngleDesired() const { return mAngleDesired; }

  unsigned int getInputCount() const override { return 1; }
  void setInput(unsigned int /*i*/, double aValue) override {
    mAngleDesired = aValue;
  }
  double getInput(unsigned int /*i*/) const override { return mAngleDesired; }

  /**
   * Default constructor.
   */
  explicit rotation_control_2D(const std::string& aName = "")
      : kte_map(aName), system_input(aName), mAngleDesired() {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  rotation_control_2D(const std::string& aName,
                      std::shared_ptr<ReaK::frame_2D<double>> aAnchor)
      : kte_map(aName),
        system_input(aName),
        mAnchor(std::move(aAnchor)),
        mAngleDesired() {}

  /**
   * Default destructor.
   */
  ~rotation_control_2D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mAnchor->Rotation.setAngle(mAngleDesired);
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::save(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::load(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(rotation_control_2D, 0xC2100046, 1,
                              "rotation_control_2D", kte_map, system_input)
};

/**
 * This class can be used as a system input to set the value of the rotation of a
 * 3D coordinate frame.
 */
class rotation_control_3D : public kte_map, public system_input {
 private:
  std::shared_ptr<frame_3D<double>> mAnchor;
  vect<double, 4> mQuatDesired;

 public:
  /**
   * Sets the frame on which rotation control is exerted.
   * \param aPtr The frame on which rotation control is exerted.
   */
  void setAnchor(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which rotation control is exerted.
   * \return The frame on which rotation control is exerted.
   */
  std::shared_ptr<frame_3D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the desired rotation.
   * \param aValue The new desired rotation.
   */
  void setQuatDesired(const vect<double, 4>& aValue) { mQuatDesired = aValue; }
  /**
   * Returns the desired rotation.
   * \return the desired rotation.
   */
  vect<double, 4> QuatDesired() const { return mQuatDesired; }

  unsigned int getInputCount() const override { return 4; }
  void setInput(unsigned int i, double aValue) override {
    if (i < 4) {
      mQuatDesired[i] = aValue;
    } else {
      mQuatDesired[0] = aValue;
    }
  }
  double getInput(unsigned int i) const override {
    if (i < 4) {
      return mQuatDesired[i];
    }
    return mQuatDesired[0];
  }

  /**
   * Default constructor.
   */
  explicit rotation_control_3D(const std::string& aName = "")
      : kte_map(aName), system_input(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  rotation_control_3D(const std::string& aName,
                      std::shared_ptr<ReaK::frame_3D<double>> aAnchor)
      : kte_map(aName), system_input(aName), mAnchor(std::move(aAnchor)) {}

  /**
   * Default destructor.
   */
  ~rotation_control_3D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mAnchor->Quat = quaternion<double>(mQuatDesired);
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::save(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::load(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(rotation_control_3D, 0xC2100047, 1,
                              "rotation_control_3D", kte_map, system_input)
};

/**
 * This class can be used as a system input to set the value of the velocity of a
 * generalized coordinate.
 */
class velocity_control_gen : public kte_map, public system_input {
 private:
  std::shared_ptr<gen_coord<double>> mAnchor;
  double mVelDesired;

 public:
  /**
   * Sets the frame on which velocity control is exerted.
   * \param aPtr The frame on which velocity control is exerted.
   */
  void setAnchor(const std::shared_ptr<gen_coord<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which velocity control is exerted.
   * \return The frame on which velocity control is exerted.
   */
  std::shared_ptr<gen_coord<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the desired velocity.
   * \param aValue The new desired velocity.
   */
  void setVelDesired(double aValue) { mVelDesired = aValue; }
  /**
   * Returns the desired velocity.
   * \return the desired velocity.
   */
  double VelDesired() const { return mVelDesired; }

  unsigned int getInputCount() const override { return 1; }
  void setInput(unsigned int /*i*/, double aValue) override {
    mVelDesired = aValue;
  }
  double getInput(unsigned int /*i*/) const override { return mVelDesired; }

  /**
   * Default constructor.
   */
  explicit velocity_control_gen(const std::string& aName = "")
      : kte_map(aName), system_input(aName), mVelDesired(0.0) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  velocity_control_gen(const std::string& aName,
                       std::shared_ptr<ReaK::gen_coord<double>> aAnchor)
      : kte_map(aName),
        system_input(aName),
        mAnchor(std::move(aAnchor)),
        mVelDesired(0.0) {}

  /**
   * Default destructor.
   */
  ~velocity_control_gen() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mAnchor->q_dot = mVelDesired;
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::save(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::load(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(velocity_control_gen, 0xC2100048, 1,
                              "velocity_control_gen", kte_map, system_input)
};

/**
 * This class can be used as a system input to set the value of the velocity of a
 * 2D coordinate frame.
 */
class velocity_control_2D : public kte_map, public system_input {
 private:
  std::shared_ptr<frame_2D<double>> mAnchor;
  vect<double, 2> mVelDesired;

 public:
  /**
   * Sets the frame on which velocity control is exerted.
   * \param aPtr The frame on which velocity control is exerted.
   */
  void setAnchor(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which velocity control is exerted.
   * \return The frame on which velocity control is exerted.
   */
  std::shared_ptr<frame_2D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the desired velocity.
   * \param aValue The new desired velocity.
   */
  void setVelDesired(const vect<double, 2>& aValue) { mVelDesired = aValue; }
  /**
   * Returns the desired velocity.
   * \return the desired velocity.
   */
  vect<double, 2> VelDesired() const { return mVelDesired; }

  unsigned int getInputCount() const override { return 2; }
  void setInput(unsigned int i, double aValue) override {
    if (i < 2) {
      mVelDesired[i] = aValue;
    } else {
      mVelDesired[0] = aValue;
    }
  }
  double getInput(unsigned int i) const override {
    if (i < 2) {
      return mVelDesired[i];
    }
    return mVelDesired[0];
  }

  /**
   * Default constructor.
   */
  explicit velocity_control_2D(const std::string& aName = "")
      : kte_map(aName), system_input(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  velocity_control_2D(const std::string& aName,
                      std::shared_ptr<ReaK::frame_2D<double>> aAnchor)
      : kte_map(aName), system_input(aName), mAnchor(std::move(aAnchor)) {}

  /**
   * Default destructor.
   */
  ~velocity_control_2D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mAnchor->Velocity = mVelDesired;
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::save(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::load(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(velocity_control_2D, 0xC2100049, 1,
                              "velocity_control_2D", kte_map, system_input)
};

/**
 * This class can be used as a system input to set the value of the velocity of a
 * 3D coordinate frame.
 */
class velocity_control_3D : public kte_map, public system_input {
 private:
  std::shared_ptr<frame_3D<double>> mAnchor;
  vect<double, 3> mVelDesired;

 public:
  /**
   * Sets the frame on which velocity control is exerted.
   * \param aPtr The frame on which velocity control is exerted.
   */
  void setAnchor(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which velocity control is exerted.
   * \return The frame on which velocity control is exerted.
   */
  std::shared_ptr<frame_3D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the desired velocity.
   * \param aValue The new desired velocity.
   */
  void setVelDesired(const vect<double, 3>& aValue) { mVelDesired = aValue; }
  /**
   * Returns the desired velocity.
   * \return the desired velocity.
   */
  vect<double, 3> VelDesired() const { return mVelDesired; }

  unsigned int getInputCount() const override { return 3; }
  void setInput(unsigned int i, double aValue) override {
    if (i < 3) {
      mVelDesired[i] = aValue;
    } else {
      mVelDesired[0] = aValue;
    }
  }
  double getInput(unsigned int i) const override {
    if (i < 3) {
      return mVelDesired[i];
    }
    return mVelDesired[0];
  }

  /**
   * Default constructor.
   */
  explicit velocity_control_3D(const std::string& aName = "")
      : kte_map(aName), system_input(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  velocity_control_3D(const std::string& aName,
                      std::shared_ptr<ReaK::frame_3D<double>> aAnchor)
      : kte_map(aName), system_input(aName), mAnchor(std::move(aAnchor)) {}

  /**
   * Default destructor.
   */
  ~velocity_control_3D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mAnchor->Velocity = mVelDesired;
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::save(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::load(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(velocity_control_3D, 0xC210004A, 1,
                              "velocity_control_3D", kte_map, system_input)
};

/**
 * This class can be used as a system input to set the value of the angular velocity of a
 * 2D coordinate frame.
 */
class ang_velocity_control_2D : public kte_map, public system_input {
 private:
  std::shared_ptr<frame_2D<double>> mAnchor;
  double mAngVelDesired;

 public:
  /**
   * Sets the frame on which angular velocity control is exerted.
   * \param aPtr The frame on which angular velocity control is exerted.
   */
  void setAnchor(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which angular velocity control is exerted.
   * \return The frame on which angular velocity control is exerted.
   */
  std::shared_ptr<frame_2D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the desired angular velocity.
   * \param aValue The new desired angular velocity.
   */
  void setAngVelDesired(double aValue) { mAngVelDesired = aValue; }
  /**
   * Returns the desired angular velocity.
   * \return the desired angular velocity.
   */
  double AngVelDesired() const { return mAngVelDesired; }

  unsigned int getInputCount() const override { return 1; }
  void setInput(unsigned int i, double aValue) override {
    mAngVelDesired = aValue;
  }
  double getInput(unsigned int i) const override { return mAngVelDesired; }

  /**
   * Default constructor.
   */
  explicit ang_velocity_control_2D(const std::string& aName = "")
      : kte_map(aName), system_input(aName), mAngVelDesired() {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  ang_velocity_control_2D(const std::string& aName,
                          std::shared_ptr<ReaK::frame_2D<double>> aAnchor)
      : kte_map(aName),
        system_input(aName),
        mAnchor(std::move(aAnchor)),
        mAngVelDesired() {}

  /**
   * Default destructor.
   */
  ~ang_velocity_control_2D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mAnchor->AngVelocity = mAngVelDesired;
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::save(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::load(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(ang_velocity_control_2D, 0xC210004B, 1,
                              "ang_velocity_control_2D", kte_map, system_input)
};

/**
 * This class can be used as a system input to set the value of the angular velocity of a
 * 3D coordinate frame.
 */
class ang_velocity_control_3D : public kte_map, public system_input {
 private:
  std::shared_ptr<frame_3D<double>> mAnchor;
  vect<double, 3> mAngVelDesired;

 public:
  /**
   * Sets the frame on which angular velocity control is exerted.
   * \param aPtr The frame on which angular velocity control is exerted.
   */
  void setAnchor(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which angular velocity control is exerted.
   * \return The frame on which angular velocity control is exerted.
   */
  std::shared_ptr<frame_3D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the desired angular velocity.
   * \param aValue The new desired angular velocity.
   */
  void setAngVelDesired(const vect<double, 3>& aValue) {
    mAngVelDesired = aValue;
  }
  /**
   * Returns the desired angular velocity.
   * \return the desired angular velocity.
   */
  vect<double, 3> AngVelDesired() const { return mAngVelDesired; }

  unsigned int getInputCount() const override { return 3; }
  void setInput(unsigned int i, double aValue) override {
    if (i < 3) {
      mAngVelDesired[i] = aValue;
    } else {
      mAngVelDesired[0] = aValue;
    }
  }
  double getInput(unsigned int i) const override {
    if (i < 3) {
      return mAngVelDesired[i];
    }
    return mAngVelDesired[0];
  }

  /**
   * Default constructor.
   */
  explicit ang_velocity_control_3D(const std::string& aName = "")
      : kte_map(aName), system_input(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  ang_velocity_control_3D(const std::string& aName,
                          std::shared_ptr<ReaK::frame_3D<double>> aAnchor)
      : kte_map(aName), system_input(aName), mAnchor(std::move(aAnchor)) {}

  /**
   * Default destructor.
   */
  ~ang_velocity_control_3D() override = default;
  ;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mAnchor->AngVelocity = mAngVelDesired;
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::save(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    system_input::load(A, system_input::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(ang_velocity_control_3D, 0xC210004C, 1,
                              "ang_velocity_control_3D", kte_map, system_input)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_STATE_CONTROLS_H_

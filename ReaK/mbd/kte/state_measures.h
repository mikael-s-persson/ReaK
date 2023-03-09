/**
 * \file state_measures.h
 *
 * This library provides a number of classes which can be used as KTE models which
 * measures the motion-state variables (position, velocity, etc.) directly. Note that
 * these classes don't implement measurement systems, simply direct measure. If a measurement
 * system is required, one can make a KTE for that and insert it between
 * a frame that is measured via the classes provided here and the output of the measurement.
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

#ifndef REAK_MBD_KTE_STATE_MEASURES_H_
#define REAK_MBD_KTE_STATE_MEASURES_H_

#include <utility>
#include "ReaK/math/kinetostatics/kinetostatics.h"
#include "ReaK/mbd/kte/kte_map.h"
#include "ReaK/mbd/kte/kte_system_output.h"

namespace ReaK::kte {

/**
 * This class can be used as a system output to get the value of the position of a
 * generalized coordinate.
 */
class position_measure_gen : public kte_map, public system_output {
 private:
  std::shared_ptr<gen_coord<double>> mAnchor;
  double mPosMeasure;

 public:
  /**
   * Sets the frame on which position measure is taken.
   * \param aPtr The new frame on which position measure is taken.
   */
  void setAnchor(const std::shared_ptr<gen_coord<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which position measure is taken.
   * \return The frame on which position measure is taken.
   */
  std::shared_ptr<gen_coord<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the measured position.
   * \param aValue The new measured position.
   */
  void setPosMeasure(double aValue) { mPosMeasure = aValue; }
  /**
   * Returns the measured position.
   * \return the measured position.
   */
  double PosMeasure() const { return mPosMeasure; }

  unsigned int getOutputCount() const override { return 1; }
  double getOutput(unsigned int /*i*/) const override { return mPosMeasure; }

  /**
   * Default constructor.
   */
  explicit position_measure_gen(const std::string& aName = "")
      : kte_map(aName), mPosMeasure(0.0) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  position_measure_gen(const std::string& aName,
                       std::shared_ptr<ReaK::gen_coord<double>> aAnchor)
      : kte_map(aName),
        system_output(aName),
        mAnchor(std::move(aAnchor)),
        mPosMeasure(0.0) {}

  /**
   * Default destructor.
   */
  ~position_measure_gen() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mPosMeasure = mAnchor->q;
    } else {
      mPosMeasure = 0.0;
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
    if (mAnchor) {
      mPosMeasure = mAnchor->q;
    } else {
      mPosMeasure = 0.0;
    }
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(position_measure_gen, 0xC2100035, 1,
                              "position_measure_gen", kte_map, system_output)
};

/**
 * This class can be used as a system output to get the value of the position of a
 * 2D coordinate frame.
 */
class position_measure_2D : public kte_map, public system_output {
 private:
  std::shared_ptr<frame_2D<double>> mAnchor;
  vect<double, 2> mPosMeasure;

 public:
  /**
   * Sets the frame on which position measure is taken.
   * \param aPtr The new frame on which position measure is taken.
   */
  void setAnchor(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which position measure is taken.
   * \return The frame on which position measure is taken.
   */
  std::shared_ptr<frame_2D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the measured position.
   * \param aValue The new measured position.
   */
  void setPosMeasure(const vect<double, 2>& aValue) { mPosMeasure = aValue; }
  /**
   * Returns the measured position.
   * \return the measured position.
   */
  vect<double, 2> PosMeasure() const { return mPosMeasure; }

  unsigned int getOutputCount() const override { return 2; }
  double getOutput(unsigned int i) const override {
    if (i < 2) {
      return mPosMeasure[i];
    }
    return mPosMeasure[0];
  }

  /**
   * Default constructor.
   */
  explicit position_measure_2D(const std::string& aName = "")
      : kte_map(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  position_measure_2D(const std::string& aName,
                      std::shared_ptr<ReaK::frame_2D<double>> aAnchor)
      : kte_map(aName), system_output(aName), mAnchor(std::move(aAnchor)) {}

  /**
   * Default destructor.
   */
  ~position_measure_2D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mPosMeasure = mAnchor->Position;
    } else {
      mPosMeasure = vect<double, 2>();
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
    if (mAnchor) {
      mPosMeasure = mAnchor->Position;
    } else {
      mPosMeasure = vect<double, 2>();
    }
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(position_measure_2D, 0xC2100036, 1,
                              "position_measure_2D", kte_map, system_output)
};

/**
 * This class can be used as a system output to get the value of the position of a
 * 3D coordinate frame.
 */
class position_measure_3D : public kte_map, public system_output {
 private:
  std::shared_ptr<frame_3D<double>> mAnchor;
  vect<double, 3> mPosMeasure;

 public:
  /**
   * Sets the frame on which position measure is taken.
   * \param aPtr The new frame on which position measure is taken.
   */
  void setAnchor(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which position measure is taken.
   * \return The frame on which position measure is taken.
   */
  std::shared_ptr<frame_3D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the measured position.
   * \param aValue The new measured position.
   */
  void setPosMeasure(const vect<double, 3>& aValue) { mPosMeasure = aValue; }
  /**
   * Returns the measured position.
   * \return the measured position.
   */
  vect<double, 3> PosMeasure() const { return mPosMeasure; }

  unsigned int getOutputCount() const override { return 3; }
  double getOutput(unsigned int i) const override {
    if (i < 3) {
      return mPosMeasure[i];
    }
    return mPosMeasure[0];
  }

  /**
   * Default constructor.
   */
  explicit position_measure_3D(const std::string& aName = "")
      : kte_map(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  position_measure_3D(const std::string& aName,
                      std::shared_ptr<frame_3D<double>> aAnchor)
      : kte_map(aName), system_output(aName), mAnchor(std::move(aAnchor)) {}

  /**
   * Default destructor.
   */
  ~position_measure_3D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mPosMeasure = mAnchor->Position;
    } else {
      mPosMeasure = vect<double, 3>();
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
    if (mAnchor) {
      mPosMeasure = mAnchor->Position;
    } else {
      mPosMeasure = vect<double, 3>();
    }
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(position_measure_3D, 0xC2100037, 1,
                              "position_measure_3D", kte_map, system_output)
};

/**
 * This class can be used as a system output to get the value of the rotation of a
 * 2D coordinate frame.
 */
class rotation_measure_2D : public kte_map, public system_output {
 private:
  std::shared_ptr<frame_2D<double>> mAnchor;
  double mAngleMeasure;

 public:
  /**
   * Sets the frame on which rotation measure is taken.
   * \param aPtr The new frame on which rotation measure is taken.
   */
  void setAnchor(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which rotation measure is taken.
   * \return The frame on which rotation measure is taken.
   */
  std::shared_ptr<frame_2D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the measured rotation.
   * \param aValue The new measured rotation.
   */
  void setAngleMeasure(double aValue) { mAngleMeasure = aValue; }
  /**
   * Returns the measured rotation.
   * \return the measured rotation.
   */
  double AngleMeasure() const { return mAngleMeasure; }

  unsigned int getOutputCount() const override { return 1; }
  double getOutput(unsigned int /*i*/) const override { return mAngleMeasure; }

  /**
   * Default constructor.
   */
  explicit rotation_measure_2D(const std::string& aName = "")
      : kte_map(aName), mAngleMeasure() {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  rotation_measure_2D(const std::string& aName,
                      std::shared_ptr<frame_2D<double>> aAnchor)
      : kte_map(aName),
        system_output(aName),
        mAnchor(std::move(aAnchor)),
        mAngleMeasure() {}

  /**
   * Default destructor.
   */
  ~rotation_measure_2D() override = default;
  ;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mAngleMeasure = mAnchor->Rotation.getAngle();
    } else {
      mAngleMeasure = 0.0;
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
    if (mAnchor) {
      mAngleMeasure = mAnchor->Rotation.getAngle();
    } else {
      mAngleMeasure = 0.0;
    }
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(rotation_measure_2D, 0xC2100038, 1,
                              "rotation_measure_2D", kte_map, system_output)
};

/**
 * This class can be used as a system output to get the value of the rotation of a
 * 3D coordinate frame.
 */
class rotation_measure_3D : public kte_map, public system_output {
 private:
  std::shared_ptr<frame_3D<double>> mAnchor;
  quaternion<double> mQuatMeasure;

 public:
  /**
   * Sets the frame on which rotation measure is taken.
   * \param aPtr The new frame on which rotation measure is taken.
   */
  void setAnchor(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which rotation measure is taken.
   * \return The frame on which rotation measure is taken.
   */
  std::shared_ptr<frame_3D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the measured rotation.
   * \param aValue The new measured rotation.
   */
  void setQuatMeasure(const quaternion<double>& aValue) {
    mQuatMeasure = aValue;
  }
  /**
   * Returns the measured rotation.
   * \return the measured rotation.
   */
  quaternion<double> QuatMeasure() const { return mQuatMeasure; }

  unsigned int getOutputCount() const override { return 4; }
  double getOutput(unsigned int i) const override {
    if (i < 4) {
      return mQuatMeasure[i];
    }
    return mQuatMeasure[0];
  }

  /**
   * Default constructor.
   */
  explicit rotation_measure_3D(const std::string& aName = "")
      : kte_map(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  rotation_measure_3D(const std::string& aName,
                      std::shared_ptr<frame_3D<double>> aAnchor)
      : kte_map(aName), system_output(aName), mAnchor(std::move(aAnchor)) {}

  /**
   * Default destructor.
   */
  ~rotation_measure_3D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mQuatMeasure = mAnchor->Quat;
    } else {
      mQuatMeasure = quaternion<double>();
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
    if (mAnchor) {
      mQuatMeasure = mAnchor->Quat;
    } else {
      mQuatMeasure = quaternion<double>();
    }
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(rotation_measure_3D, 0xC2100039, 1,
                              "rotation_measure_3D", kte_map, system_output)
};

/**
 * This class can be used as a system output to get the value of the velocity of a
 * generalized coordinate.
 */
class velocity_measure_gen : public kte_map, public system_output {
 private:
  std::shared_ptr<gen_coord<double>> mAnchor;
  double mVelMeasure;

 public:
  /**
   * Sets the frame on which velocity measure is taken.
   * \param aPtr The new frame on which velocity measure is taken.
   */
  void setAnchor(const std::shared_ptr<gen_coord<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which velocity measure is taken.
   * \return The frame on which velocity measure is taken.
   */
  std::shared_ptr<gen_coord<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the measured velocity.
   * \param aValue The new measured velocity.
   */
  void setVelMeasure(double aValue) { mVelMeasure = aValue; }
  /**
   * Returns the measured velocity.
   * \return the measured velocity.
   */
  double VelMeasure() const { return mVelMeasure; }

  unsigned int getOutputCount() const override { return 1; }
  double getOutput(unsigned int /*i*/) const override { return mVelMeasure; }

  /**
   * Default constructor.
   */
  explicit velocity_measure_gen(const std::string& aName = "")
      : kte_map(aName), mVelMeasure(0.0) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  velocity_measure_gen(const std::string& aName,
                       std::shared_ptr<gen_coord<double>> aAnchor)
      : kte_map(aName),
        system_output(aName),
        mAnchor(std::move(aAnchor)),
        mVelMeasure(0.0) {}

  /**
   * Default destructor.
   */
  ~velocity_measure_gen() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mVelMeasure = mAnchor->q_dot;
    } else {
      mVelMeasure = 0.0;
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
    if (mAnchor) {
      mVelMeasure = mAnchor->q_dot;
    } else {
      mVelMeasure = 0.0;
    }
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(velocity_measure_gen, 0xC210003A, 1,
                              "velocity_measure_gen", kte_map, system_output)
};

/**
 * This class can be used as a system output to get the value of the velocity of a
 * 2D coordinate frame.
 */
class velocity_measure_2D : public kte_map, public system_output {
 private:
  std::shared_ptr<frame_2D<double>> mAnchor;
  vect<double, 2> mVelMeasure;

 public:
  /**
   * Sets the frame on which velocity measure is taken.
   * \param aPtr The new frame on which velocity measure is taken.
   */
  void setAnchor(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which velocity measure is taken.
   * \return The frame on which velocity measure is taken.
   */
  std::shared_ptr<frame_2D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the measured velocity.
   * \param aValue The new measured velocity.
   */
  void setVelMeasure(const vect<double, 2>& aValue) { mVelMeasure = aValue; }
  /**
   * Returns the measured velocity.
   * \return the measured velocity.
   */
  vect<double, 2> VelMeasure() const { return mVelMeasure; }

  unsigned int getOutputCount() const override { return 2; }
  double getOutput(unsigned int i) const override {
    if (i < 2) {
      return mVelMeasure[i];
    }
    return mVelMeasure[0];
  }

  /**
   * Default constructor.
   */
  explicit velocity_measure_2D(const std::string& aName = "")
      : kte_map(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  velocity_measure_2D(const std::string& aName,
                      std::shared_ptr<frame_2D<double>> aAnchor)
      : kte_map(aName), system_output(aName), mAnchor(std::move(aAnchor)) {}

  /**
   * Default destructor.
   */
  ~velocity_measure_2D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mVelMeasure = mAnchor->Velocity;
    } else {
      mVelMeasure = vect<double, 2>();
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
    if (mAnchor) {
      mVelMeasure = mAnchor->Velocity;
    } else {
      mVelMeasure = vect<double, 2>();
    }
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(velocity_measure_2D, 0xC210003B, 1,
                              "velocity_measure_2D", kte_map, system_output)
};

/**
 * This class can be used as a system output to get the value of the velocity of a
 * 3D coordinate frame.
 */
class velocity_measure_3D : public kte_map, public system_output {
 private:
  std::shared_ptr<frame_3D<double>> mAnchor;
  vect<double, 3> mVelMeasure;

 public:
  /**
   * Sets the frame on which velocity measure is taken.
   * \param aPtr The new frame on which velocity measure is taken.
   */
  void setAnchor(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which velocity measure is taken.
   * \return The frame on which velocity measure is taken.
   */
  std::shared_ptr<frame_3D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the measured velocity.
   * \param aValue The new measured velocity.
   */
  void setVelMeasure(const vect<double, 3>& aValue) { mVelMeasure = aValue; }
  /**
   * Returns the measured velocity.
   * \return the measured velocity.
   */
  vect<double, 3> VelMeasure() const { return mVelMeasure; }

  unsigned int getOutputCount() const override { return 3; }
  double getOutput(unsigned int i) const override {
    if (i < 3) {
      return mVelMeasure[i];
    }
    return mVelMeasure[0];
  }

  /**
   * Default constructor.
   */
  explicit velocity_measure_3D(const std::string& aName = "")
      : kte_map(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  velocity_measure_3D(const std::string& aName,
                      std::shared_ptr<frame_3D<double>> aAnchor)
      : kte_map(aName), system_output(aName), mAnchor(std::move(aAnchor)) {}

  /**
   * Default destructor.
   */
  ~velocity_measure_3D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mVelMeasure = mAnchor->Velocity;
    } else {
      mVelMeasure = vect<double, 3>();
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
    if (mAnchor) {
      mVelMeasure = mAnchor->Velocity;
    } else {
      mVelMeasure = vect<double, 3>();
    }
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(velocity_measure_3D, 0xC210003C, 1,
                              "velocity_measure_3D", kte_map, system_output)
};

/**
 * This class can be used as a system output to get the value of the angular velocity of a
 * 3D coordinate frame.
 */
class ang_velocity_measure_2D : public kte_map, public system_output {
 private:
  std::shared_ptr<frame_2D<double>> mAnchor;
  double mAngVelMeasure;

 public:
  /**
   * Sets the frame on which angular velocity measure is taken.
   * \param aPtr The new frame on which angular velocity measure is taken.
   */
  void setAnchor(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which angular velocity measure is taken.
   * \return The frame on which angular velocity measure is taken.
   */
  std::shared_ptr<frame_2D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the measured angular velocity.
   * \param aValue The new measured angular velocity.
   */
  void setAngVelMeasure(double aValue) { mAngVelMeasure = aValue; }
  /**
   * Returns the measured angular velocity.
   * \return the measured angular velocity.
   */
  double AngVelMeasure() const { return mAngVelMeasure; }

  unsigned int getOutputCount() const override { return 1; }
  double getOutput(unsigned int /*i*/) const override { return mAngVelMeasure; }

  /**
   * Default constructor.
   */
  explicit ang_velocity_measure_2D(const std::string& aName = "")
      : kte_map(aName), mAngVelMeasure() {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  ang_velocity_measure_2D(const std::string& aName,
                          std::shared_ptr<frame_2D<double>> aAnchor)
      : kte_map(aName),
        system_output(aName),
        mAnchor(std::move(aAnchor)),
        mAngVelMeasure() {}

  /**
   * Default destructor.
   */
  ~ang_velocity_measure_2D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mAngVelMeasure = mAnchor->AngVelocity;
    } else {
      mAngVelMeasure = 0.0;
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
    if (mAnchor) {
      mAngVelMeasure = mAnchor->AngVelocity;
    } else {
      mAngVelMeasure = 0.0;
    }
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(ang_velocity_measure_2D, 0xC210003B, 1,
                              "ang_velocity_measure_2D", kte_map, system_output)
};

/**
 * This class can be used as a system output to get the value of the angular velocity of a
 * 3D coordinate frame.
 */
class ang_velocity_measure_3D : public kte_map, public system_output {
 private:
  std::shared_ptr<frame_3D<double>> mAnchor;
  vect<double, 3> mAngVelMeasure;

 public:
  /**
   * Sets the frame on which angular velocity measure is taken.
   * \param aPtr The new frame on which angular velocity measure is taken.
   */
  void setAnchor(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor = aPtr;
  }
  /**
   * Returns the frame on which angular velocity measure is taken.
   * \return The frame on which angular velocity measure is taken.
   */
  std::shared_ptr<frame_3D<double>> Anchor() const { return mAnchor; }

  /**
   * Sets the measured angular velocity.
   * \param aValue The new measured angular velocity.
   */
  void setAngVelMeasure(const vect<double, 3>& aValue) {
    mAngVelMeasure = aValue;
  }
  /**
   * Returns the measured angular velocity.
   * \return the measured angular velocity.
   */
  vect<double, 3> AngVelMeasure() const { return mAngVelMeasure; }

  unsigned int getOutputCount() const override { return 3; }
  double getOutput(unsigned int i) const override {
    if (i < 3) {
      return mAngVelMeasure[i];
    }
    return mAngVelMeasure[0];
  }

  /**
   * Default constructor.
   */
  explicit ang_velocity_measure_3D(const std::string& aName = "")
      : kte_map(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor the coordinate from which position is measured.
   */
  ang_velocity_measure_3D(const std::string& aName,
                          std::shared_ptr<frame_3D<double>> aAnchor)
      : kte_map(aName), system_output(aName), mAnchor(std::move(aAnchor)) {}

  /**
   * Default destructor.
   */
  ~ang_velocity_measure_3D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    if (mAnchor) {
      mAngVelMeasure = mAnchor->AngVelocity;
    } else {
      mAngVelMeasure = vect<double, 3>();
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {}

  void clearForce() override {}

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor);
    if (mAnchor) {
      mAngVelMeasure = mAnchor->AngVelocity;
    } else {
      mAngVelMeasure = vect<double, 3>();
    }
  }

  RK_RTTI_MAKE_CONCRETE_2BASE(ang_velocity_measure_3D, 0xC210003C, 1,
                              "ang_velocity_measure_3D", kte_map, system_output)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_STATE_MEASURES_H_

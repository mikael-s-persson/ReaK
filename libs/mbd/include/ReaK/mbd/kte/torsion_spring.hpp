/**
 * \file torsion_spring.hpp
 *
 * This library declares the KTE models for linear torsion springs.
 * Here spring classes are available for 2D and 3D point-to-point angular difference.
 * The model of the spring is a basic linear, constant
 * torsion stiffness with optional staturation value based on Hooke's Law to obtain the torque that will oppose the
 *relative
 * angle between the anchors.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date May 2010
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

#ifndef REAK_TORSION_SPRING_HPP
#define REAK_TORSION_SPRING_HPP

#include "ReaK/mbd/kte/kte_map.hpp"

#include <utility>
#include "ReaK/math/kinetostatics/kinetostatics.hpp"

namespace ReaK::kte {

/**
 * This class implements a torsion spring model in 2D space. The model of the torsion spring is a basic linear, constant
 * torsion stiffness with optional staturation value based on Hooke's Law to obtain the torque that will oppose the
 * relative
 * angle between the anchors.
 */
class torsion_spring_2D : public kte_map {
 private:
  std::shared_ptr<frame_2D<double>> mAnchor1;  ///< Holds the first 2D frame.
  std::shared_ptr<frame_2D<double>> mAnchor2;  ///< Holds the second 2D frame.
  double mStiffness;  ///< Holds the torsion stiffness of the spring.
  double
      mSaturation;  ///< Holds the saturation torque, or maximum torque the spring can exert, if 0 there is no saturation.

 public:
  /**
   * Sets the first anchor frame of the torsion spring.
   * \param aPtr The new first anchor frame of the torsion spring.
   */
  void setAnchor1(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor1 = aPtr;
  }
  /**
   * Returns the first anchor frame of the torsion spring.
   * \return The first anchor frame of the torsion spring.
   */
  std::shared_ptr<frame_2D<double>> Anchor1() const { return mAnchor1; }

  /**
   * Sets the second anchor frame of the torsion spring.
   * \param aPtr The new second anchor frame of the torsion spring.
   */
  void setAnchor2(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor2 = aPtr;
  }
  /**
   * Returns the second anchor frame of the torsion spring.
   * \return The second anchor frame of the torsion spring.
   */
  std::shared_ptr<frame_2D<double>> Anchor2() const { return mAnchor2; }

  /**
   * Sets the stiffness value of the torsion spring.
   * \param aValue The new stiffness value of the torsion spring.
   */
  void setStiffness(double aValue) { mStiffness = aValue; }
  /**
   * Returns the stiffness value of the torsion spring.
   * \return The stiffness value of the torsion spring.
   */
  double Stiffness() const { return mStiffness; }

  /**
   * Sets the saturation torque value of the torsion spring (0 implies no saturation at all).
   * \param aValue The new saturation torque value of the torsion spring (0 implies no saturation at all).
   */
  void setSaturation(double aValue) { mSaturation = aValue; }
  /**
   * Returns the value of the saturation torque of the torsion spring (0 implies no saturation at all).
   * \return The value of the saturation torque of the torsion spring (0 implies no saturation at all).
   */
  double Saturation() const { return mSaturation; }

  /**
   * Default constructor.
   */
  explicit torsion_spring_2D(const std::string& aName = "")
      : kte_map(aName),

        mStiffness(0.0),
        mSaturation(0.0) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor1 first attach point of the spring.
   * \param aAnchor2 second attach point of the spring.
   * \param aStiffness torsion stiffness coefficient (in Nm/rad).
   * \param aSaturation saturation torque of the spring, default 0 will disable saturation.
   */
  torsion_spring_2D(const std::string& aName,
                    std::shared_ptr<frame_2D<double>> aAnchor1,
                    std::shared_ptr<frame_2D<double>> aAnchor2,
                    double aStiffness, double aSaturation = 0.0)
      : kte_map(aName),
        mAnchor1(std::move(aAnchor1)),
        mAnchor2(std::move(aAnchor2)),
        mStiffness(aStiffness),
        mSaturation(aSaturation) {}

  /**
   * Default destructor.
   */
  ~torsion_spring_2D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override;

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override;

  void clearForce() override;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor1) & RK_SERIAL_SAVE_WITH_NAME(mAnchor2) &
        RK_SERIAL_SAVE_WITH_NAME(mStiffness) &
        RK_SERIAL_SAVE_WITH_NAME(mSaturation);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor1) & RK_SERIAL_LOAD_WITH_NAME(mAnchor2) &
        RK_SERIAL_LOAD_WITH_NAME(mStiffness) &
        RK_SERIAL_LOAD_WITH_NAME(mSaturation);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(torsion_spring_2D, 0xC210002C, 1,
                              "torsion_spring_2D", kte_map)
};

/**
 * This class implements a torsion spring model in 3D space. The model of the torsion spring is a basic linear, constant
 * torsion stiffness with optional staturation value based on Hooke's Law to obtain the torque that will oppose the
 * relative
 * angle between the anchors.
 */
class torsion_spring_3D : public kte_map {
 private:
  std::shared_ptr<frame_3D<double>> mAnchor1;  ///< Holds the first 3D frame.
  std::shared_ptr<frame_3D<double>> mAnchor2;  ///< Holds the second 2D frame.
  double mStiffness;  ///< Holds the torsion stiffness of the spring.
  double
      mSaturation;  ///< Holds the saturation torque, or maximum torque the spring can exert, if 0 there is no saturation.

 public:
  /**
   * Sets the first anchor frame of the torsion spring.
   * \param aPtr The new first anchor frame of the torsion spring.
   */
  void setAnchor1(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor1 = aPtr;
  }
  /**
   * Returns the first anchor frame of the torsion spring.
   * \return The first anchor frame of the torsion spring.
   */
  std::shared_ptr<frame_3D<double>> Anchor1() const { return mAnchor1; }

  /**
   * Sets the second anchor frame of the torsion spring.
   * \param aPtr The new second anchor frame of the torsion spring.
   */
  void setAnchor2(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor2 = aPtr;
  }
  /**
   * Returns the second anchor frame of the torsion spring.
   * \return The second anchor frame of the torsion spring.
   */
  std::shared_ptr<frame_3D<double>> Anchor2() const { return mAnchor2; }

  /**
   * Sets the stiffness value of the torsion spring.
   * \param aValue The new stiffness value of the torsion spring.
   */
  void setStiffness(double aValue) { mStiffness = aValue; }
  /**
   * Returns the stiffness value of the torsion spring.
   * \return The stiffness value of the torsion spring.
   */
  double Stiffness() const { return mStiffness; }

  /**
   * Sets the saturation torque value of the torsion spring (0 implies no saturation at all).
   * \param aValue The new saturation torque value of the torsion spring (0 implies no saturation at all).
   */
  void setSaturation(double aValue) { mSaturation = aValue; }
  /**
   * Returns the value of the saturation torque of the torsion spring (0 implies no saturation at all).
   * \return The value of the saturation torque of the torsion spring (0 implies no saturation at all).
   */
  double Saturation() const { return mSaturation; }

  /**
   * Default constructor.
   */
  explicit torsion_spring_3D(const std::string& aName = "")
      : kte_map(aName),

        mStiffness(0.0),
        mSaturation(0.0) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor1 first attach point of the spring.
   * \param aAnchor2 second attach point of the spring.
   * \param aStiffness torsion stiffness coefficient (in Nm/rad).
   * \param aSaturation saturation torque of the spring, default 0 will disable saturation.
   */
  torsion_spring_3D(const std::string& aName,
                    std::shared_ptr<frame_3D<double>> aAnchor1,
                    std::shared_ptr<frame_3D<double>> aAnchor2,
                    double aStiffness, double aSaturation = 0.0)
      : kte_map(aName),
        mAnchor1(std::move(aAnchor1)),
        mAnchor2(std::move(aAnchor2)),
        mStiffness(aStiffness),
        mSaturation(aSaturation) {}

  /**
   * Default destructor.
   */
  ~torsion_spring_3D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override;

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override;

  void clearForce() override;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor1) & RK_SERIAL_SAVE_WITH_NAME(mAnchor2) &
        RK_SERIAL_SAVE_WITH_NAME(mStiffness) &
        RK_SERIAL_SAVE_WITH_NAME(mSaturation);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor1) & RK_SERIAL_LOAD_WITH_NAME(mAnchor2) &
        RK_SERIAL_LOAD_WITH_NAME(mStiffness) &
        RK_SERIAL_LOAD_WITH_NAME(mSaturation);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(torsion_spring_3D, 0xC210002D, 1,
                              "torsion_spring_3D", kte_map)
};

}  // namespace ReaK::kte

#endif

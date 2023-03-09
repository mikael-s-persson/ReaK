/**
 * \file flexible_beam.h
 *
 * This library declares KTE models for a flexible beam in 2D and 3D. The flexible beam is
 * implemented as a linear iso-tropic stiffness for both elongation and torsion. This
 * constitutes a simplified Euler-Bernoulli beam, akin to beam elements in finite element methods.
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

#ifndef REAK_MBD_KTE_FLEXIBLE_BEAM_H_
#define REAK_MBD_KTE_FLEXIBLE_BEAM_H_

#include <utility>
#include "ReaK/math/kinetostatics/kinetostatics.h"
#include "ReaK/mbd/kte/kte_map.h"

namespace ReaK::kte {

/**
 * This class implements the 2D model of a flexible beam. The beam is made of three 2D frames, two inputs and
 * one output (kinematically speaking). The differencial displacement between the two anchors is used to compute
 * restitution forces and torques at the anchors. The average or common-mode or net motion of the two anchors is
 * used to compute the kinematics of the center of the beam, the object-frame. Furthermore, the net forces applied
 * on the object-frame is rigidly transmitted and evenly split to the anchors, in addition to the restitution forces.
 */
class flexible_beam_2D : public kte_map {
 private:
  std::shared_ptr<frame_2D<double>>
      mAnchor1;  ///< Holds the first end of the beam.
  std::shared_ptr<frame_2D<double>>
      mAnchor2;  ///< Holds the second end of the beam.
  std::shared_ptr<frame_2D<double>>
      mObjectFrame;    ///< Holds the center frame of the beam (i.e. its bulk).
  double mRestLength;  ///< The undeformed length of the beam.
  double
      mStiffness;  ///< The linear stiffness of the beam, stress-strain relation, iso-tropically.
  double
      mTorsionStiffness;  ///< The angular or torsion stiffness of the beam, iso-tropically.

 public:
  /**
   * Sets the first anchor frame of the beam.
   * \param aPtr A pointer to the new first anchor frame of the beam.
   */
  void setAnchor1(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor1 = aPtr;
  }
  /**
   * Returns a const-reference to the first anchor frame of the beam.
   * \return A const-reference to the first anchor frame of the beam.
   */
  std::shared_ptr<frame_2D<double>> Anchor1() const { return mAnchor1; }

  /**
   * Sets the first anchor frame of the beam.
   * \param aPtr A pointer to the new first anchor frame of the beam.
   */
  void setAnchor2(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor2 = aPtr;
  }
  /**
   * Returns a const-reference to the second anchor frame of the beam.
   * \return A const-reference to the second anchor frame of the beam.
   */
  std::shared_ptr<frame_2D<double>> Anchor2() const { return mAnchor2; }

  /**
   * Sets the center frame of the beam (i.e. its bulk).
   * \param aPtr A pointer to the new center frame of the beam (i.e. its bulk).
   */
  void setCenterFrame(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mObjectFrame = aPtr;
  }
  /**
   * Returns the center frame of the beam (i.e. its bulk).
   * \return The center frame of the beam (i.e. its bulk).
   */
  std::shared_ptr<frame_2D<double>> CenterFrame() const { return mObjectFrame; }

  /**
   * Sets the rest-length of the beam.
   * \param aValue The new rest-length of the beam.
   */
  void setRestLength(double aValue) { mRestLength = aValue; }
  /**
   * Returns the rest-length of the beam.
   * \return The rest-length of the beam.
   */
  double RestLength() const { return mRestLength; }

  /**
   * Sets the stiffness value of the beam.
   * \param aValue The new stiffness value of the beam.
   */
  void setStiffness(double aValue) { mStiffness = aValue; }
  /**
   * Returns the stiffness value of the beam.
   * \return The stiffness value of the beam.
   */
  double Stiffness() const { return mStiffness; }

  /**
   * Sets the torsion stiffness value of the beam.
   * \param aValue The new torsion stiffness value of the beam.
   */
  void setTorsionStiffness(double aValue) { mTorsionStiffness = aValue; }
  /**
   * Returns the torsion stiffness value of the beam.
   * \return The torsion stiffness value of the beam.
   */
  double TorsionStiffness() const { return mTorsionStiffness; }

  /**
   * Default constructor.
   */
  explicit flexible_beam_2D(const std::string& aName = "")
      : kte_map(aName),

        mRestLength(0.0),
        mStiffness(0.0),
        mTorsionStiffness(0.0) {}

  /**
   * Parametrized constructor.
   * \param aName the name of the KTE model.
   * \param aAnchor1 the first end of the beam (kinematic input).
   * \param aAnchor2 the second end of the beam (kinematic input).
   * \param aObjectFrame the center frame of the beam (i.e. its bulk).
   * \param aRestLength the undeformed length of the beam.
   * \param aStiffness the linear stiffness of the beam, stress-strain relation, iso-tropically.
   * \param aTorsionStiffness the angular or torsion stiffness of the beam, iso-tropically.
   */
  flexible_beam_2D(const std::string& aName,
                   std::shared_ptr<frame_2D<double>> aAnchor1,
                   std::shared_ptr<frame_2D<double>> aAnchor2,
                   std::shared_ptr<frame_2D<double>> aObjectFrame,
                   double aRestLength, double aStiffness,
                   double aTorsionStiffness)
      : kte_map(aName),
        mAnchor1(std::move(aAnchor1)),
        mAnchor2(std::move(aAnchor2)),
        mObjectFrame(std::move(aObjectFrame)),
        mRestLength(aRestLength),
        mStiffness(aStiffness),
        mTorsionStiffness(aTorsionStiffness) {}

  /**
   * Default destructor.
   */
  ~flexible_beam_2D() override = default;

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
        RK_SERIAL_SAVE_WITH_NAME(mObjectFrame) &
        RK_SERIAL_SAVE_WITH_NAME(mRestLength) &
        RK_SERIAL_SAVE_WITH_NAME(mStiffness) &
        RK_SERIAL_SAVE_WITH_NAME(mTorsionStiffness);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor1) & RK_SERIAL_LOAD_WITH_NAME(mAnchor2) &
        RK_SERIAL_LOAD_WITH_NAME(mObjectFrame) &
        RK_SERIAL_LOAD_WITH_NAME(mRestLength) &
        RK_SERIAL_LOAD_WITH_NAME(mStiffness) &
        RK_SERIAL_LOAD_WITH_NAME(mTorsionStiffness);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(flexible_beam_2D, 0xC210001D, 1,
                              "flexible_beam_2D", kte_map)
};

/**
 * This class implements the 3D model of a flexible beam. The beam is made of three 3D frames, two inputs and
 * one output (kinematically speaking). The differencial displacement between the two anchors is used to compute
 * restitution forces and torques at the anchors. The average or common-mode or net motion of the two anchors is
 * used to compute the kinematics of the center of the beam, the object-frame. Furthermore, the net forces applied
 * on the object-frame is rigidly transmitted and evenly split to the anchors, in addition to the restitution forces.
 */
class flexible_beam_3D : public kte_map {
 private:
  std::shared_ptr<frame_3D<double>>
      mAnchor1;  ///< Holds the first end of the beam.
  std::shared_ptr<frame_3D<double>>
      mAnchor2;  ///< Holds the second end of the beam.
  std::shared_ptr<frame_3D<double>>
      mObjectFrame;    ///< Holds the center frame of the beam (i.e. its bulk).
  double mRestLength;  ///< The undeformed length of the beam.
  double
      mStiffness;  ///< The linear stiffness of the beam, stress-strain relation, iso-tropically.
  double
      mTorsionStiffness;  ///< The angular or torsion stiffness of the beam, iso-tropically.

 public:
  /**
   * Sets the first anchor frame of the beam.
   * \param aPtr A pointer to the new first anchor frame of the beam.
   */
  void setAnchor1(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor1 = aPtr;
  }
  /**
   * Returns a const-reference to the first anchor frame of the beam.
   * \return A const-reference to the first anchor frame of the beam.
   */
  std::shared_ptr<frame_3D<double>> Anchor1() const { return mAnchor1; }

  /**
   * Sets the first anchor frame of the beam.
   * \param aPtr A pointer to the new first anchor frame of the beam.
   */
  void setAnchor2(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor2 = aPtr;
  }
  /**
   * Returns a const-reference to the second anchor frame of the beam.
   * \return A const-reference to the second anchor frame of the beam.
   */
  std::shared_ptr<frame_3D<double>> Anchor2() const { return mAnchor2; }

  /**
   * Sets the center frame of the beam (i.e. its bulk).
   * \param aPtr A pointer to the new center frame of the beam (i.e. its bulk).
   */
  void setCenterFrame(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mObjectFrame = aPtr;
  }
  /**
   * Returns the center frame of the beam (i.e. its bulk).
   * \return The center frame of the beam (i.e. its bulk).
   */
  std::shared_ptr<frame_3D<double>> CenterFrame() const { return mObjectFrame; }

  /**
   * Sets the rest-length of the beam.
   * \param aValue The new rest-length of the beam.
   */
  void setRestLength(double aValue) { mRestLength = aValue; }
  /**
   * Returns the rest-length of the beam.
   * \return The rest-length of the beam.
   */
  double RestLength() const { return mRestLength; }

  /**
   * Sets the stiffness value of the beam.
   * \param aValue The new stiffness value of the beam.
   */
  void setStiffness(double aValue) { mStiffness = aValue; }
  /**
   * Returns the stiffness value of the beam.
   * \return The stiffness value of the beam.
   */
  double Stiffness() const { return mStiffness; }

  /**
   * Sets the torsion stiffness value of the beam.
   * \param aValue The new torsion stiffness value of the beam.
   */
  void setTorsionStiffness(double aValue) { mTorsionStiffness = aValue; }
  /**
   * Returns the torsion stiffness value of the beam.
   * \return The torsion stiffness value of the beam.
   */
  double TorsionStiffness() const { return mTorsionStiffness; }

  /**
   * Default constructor.
   */
  explicit flexible_beam_3D(const std::string& aName = "")
      : kte_map(aName),

        mRestLength(0.0),
        mStiffness(0.0),
        mTorsionStiffness(0.0) {}

  /**
   * Parametrized constructor.
   * \param aName the name of the KTE model.
   * \param aAnchor1 the first end of the beam (kinematic input).
   * \param aAnchor2 the second end of the beam (kinematic input).
   * \param aObjectFrame the center frame of the beam (i.e. its bulk).
   * \param aRestLength the undeformed length of the beam.
   * \param aStiffness the linear stiffness of the beam, stress-strain relation, iso-tropically.
   * \param aTorsionStiffness the angular or torsion stiffness of the beam, iso-tropically.
   */
  flexible_beam_3D(const std::string& aName,
                   std::shared_ptr<frame_3D<double>> aAnchor1,
                   std::shared_ptr<frame_3D<double>> aAnchor2,
                   std::shared_ptr<frame_3D<double>> aObjectFrame,
                   double aRestLength, double aStiffness,
                   double aTorsionStiffness)
      : kte_map(aName),
        mAnchor1(std::move(aAnchor1)),
        mAnchor2(std::move(aAnchor2)),
        mObjectFrame(std::move(aObjectFrame)),
        mRestLength(aRestLength),
        mStiffness(aStiffness),
        mTorsionStiffness(aTorsionStiffness) {}

  /**
   * Default destructor.
   */
  ~flexible_beam_3D() override = default;

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
        RK_SERIAL_SAVE_WITH_NAME(mObjectFrame) &
        RK_SERIAL_SAVE_WITH_NAME(mRestLength) &
        RK_SERIAL_SAVE_WITH_NAME(mStiffness) &
        RK_SERIAL_SAVE_WITH_NAME(mTorsionStiffness);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor1) & RK_SERIAL_LOAD_WITH_NAME(mAnchor2) &
        RK_SERIAL_LOAD_WITH_NAME(mObjectFrame) &
        RK_SERIAL_LOAD_WITH_NAME(mRestLength) &
        RK_SERIAL_LOAD_WITH_NAME(mStiffness) &
        RK_SERIAL_LOAD_WITH_NAME(mTorsionStiffness);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(flexible_beam_3D, 0xC210001D, 1,
                              "flexible_beam_3D", kte_map)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_FLEXIBLE_BEAM_H_

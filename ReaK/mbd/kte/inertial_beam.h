/**
 * \file inertial_beam.h
 *
 * This library declares KTE models for an inertial beam in 2D and 3D. This compliant
 * beam model does not apply any restitution forces but rather models only the inertia
 * of a beam with two, lumped end-masses, handling also the moment of inertia generated
 * by the relative rotational motion of the two end-points.
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

#ifndef REAK_MBD_KTE_INERTIAL_BEAM_H_
#define REAK_MBD_KTE_INERTIAL_BEAM_H_

#include <utility>
#include "ReaK/math/kinetostatics/kinetostatics.h"
#include "ReaK/mbd/kte/kte_map.h"

namespace ReaK::kte {

/**
 * This class implements the inertial forces of a 2D beam member between two end-points (anchors).
 * The inertia is conveniently modelled as two lumped end-masses.
 */
class inertial_beam_2D : public kte_map {
 private:
  std::shared_ptr<frame_2D<double>>
      mAnchor1;  ///< Holds the first end of the beam.
  std::shared_ptr<frame_2D<double>>
      mAnchor2;  ///< Holds the second end of the beam.
  double
      mMass;  ///< Holds the total mass of the beam, half of which is lumped at each end.

 public:
  /**
   * Sets the first anchor frame of the inertial beam.
   * \param aPtr A pointer to the new first anchor frame of the inertial beam.
   */
  void setAnchor1(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor1 = aPtr;
  }
  /**
   * Returns a const-reference to the first anchor frame of the inertial beam.
   * \return A const-reference to the first anchor frame of the inertial beam.
   */
  std::shared_ptr<frame_2D<double>> Anchor1() const { return mAnchor1; }

  /**
   * Sets the first anchor frame of the inertial beam.
   * \param aPtr A pointer to the new first anchor frame of the inertial beam.
   */
  void setAnchor2(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor2 = aPtr;
  }
  /**
   * Returns a const-reference to the second anchor frame of the inertial beam.
   * \return A const-reference to the second anchor frame of the inertial beam.
   */
  std::shared_ptr<frame_2D<double>> Anchor2() const { return mAnchor2; }

  /**
   * Sets the mass of the inertial beam.
   * \param aValue The new mass of the inertial beam.
   */
  void setMass(double aValue) { mMass = aValue; }
  /**
   * Returns the mass of the inertial beam.
   * \return The mass of the inertial beam.
   */
  double Mass() const { return mMass; }

  /**
   * Default constructor.
   */
  explicit inertial_beam_2D(const std::string& aName = "")
      : kte_map(aName), mMass(0.0) {}

  /**
   * Parametrized constructor.
   * \param aName the name of the KTE model.
   * \param aAnchor1 the first end of the beam.
   * \param aAnchor2 the second end of the beam.
   * \param aMass the total mass of the beam, half of which is lumped at each end.
   */
  inertial_beam_2D(const std::string& aName,
                   std::shared_ptr<frame_2D<double>> aAnchor1,
                   std::shared_ptr<frame_2D<double>> aAnchor2, double aMass)
      : kte_map(aName),
        mAnchor1(std::move(aAnchor1)),
        mAnchor2(std::move(aAnchor2)),
        mMass(aMass) {}

  /**
   * Default destructor.
   */
  ~inertial_beam_2D() override = default;

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
        RK_SERIAL_SAVE_WITH_NAME(mMass);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor1) & RK_SERIAL_LOAD_WITH_NAME(mAnchor2) &
        RK_SERIAL_LOAD_WITH_NAME(mMass);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(inertial_beam_2D, 0xC210001E, 1,
                              "inertial_beam_2D", kte_map)
};

/**
 * This class implements the inertial forces of a 3D beam member between two end-points (anchors).
 * The inertia is conveniently modelled as two lumped end-masses.
 */
class inertial_beam_3D : public kte_map {
 private:
  std::shared_ptr<frame_3D<double>>
      mAnchor1;  ///< Holds the first end of the beam.
  std::shared_ptr<frame_3D<double>>
      mAnchor2;  ///< Holds the second end of the beam.
  double
      mMass;  ///< Holds the total mass of the beam, half of which is lumped at each end.

 public:
  /**
   * Sets the first anchor frame of the inertial beam.
   * \param aPtr A pointer to the new first anchor frame of the inertial beam.
   */
  void setAnchor1(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor1 = aPtr;
  }
  /**
   * Returns a const-reference to the first anchor frame of the inertial beam.
   * \return A const-reference to the first anchor frame of the inertial beam.
   */
  std::shared_ptr<frame_3D<double>> Anchor1() const { return mAnchor1; }

  /**
   * Sets the first anchor frame of the inertial beam.
   * \param aPtr A pointer to the new first anchor frame of the inertial beam.
   */
  void setAnchor2(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor2 = aPtr;
  }
  /**
   * Returns a const-reference to the second anchor frame of the inertial beam.
   * \return A const-reference to the second anchor frame of the inertial beam.
   */
  std::shared_ptr<frame_3D<double>> Anchor2() const { return mAnchor2; }

  /**
   * Sets the mass of the inertial beam.
   * \param aValue The new mass of the inertial beam.
   */
  void setMass(double aValue) { mMass = aValue; }
  /**
   * Returns the mass of the inertial beam.
   * \return The mass of the inertial beam.
   */
  double Mass() const { return mMass; }

  /**
   * Default constructor.
   */
  explicit inertial_beam_3D(const std::string& aName = "")
      : kte_map(aName), mMass(0.0) {}

  /**
   * Parametrized constructor.
   * \param aName the name of the KTE model.
   * \param aAnchor1 the first end of the beam.
   * \param aAnchor2 the second end of the beam.
   * \param aMass the total mass of the beam, half of which is lumped at each end.
   */
  inertial_beam_3D(const std::string& aName,
                   std::shared_ptr<frame_3D<double>> aAnchor1,
                   std::shared_ptr<frame_3D<double>> aAnchor2, double aMass)
      : kte_map(aName),
        mAnchor1(std::move(aAnchor1)),
        mAnchor2(std::move(aAnchor2)),
        mMass(aMass) {}

  /**
   * Default destructor.
   */
  ~inertial_beam_3D() override = default;

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
        RK_SERIAL_SAVE_WITH_NAME(mMass);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor1) & RK_SERIAL_LOAD_WITH_NAME(mAnchor2) &
        RK_SERIAL_LOAD_WITH_NAME(mMass);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(inertial_beam_3D, 0xC210001F, 1,
                              "inertial_beam_3D", kte_map)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_INERTIAL_BEAM_H_

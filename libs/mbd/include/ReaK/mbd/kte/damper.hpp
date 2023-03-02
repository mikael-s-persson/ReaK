/**
 * \file damper.hpp
 *
 * This library declares the KTE models for linear dampers, often reference to as dashpots.
 * Here damper classes are available for 2D and 3D point-to-point dashpots as well as a
 * damper between two generalized coordinates. The model of the damper is a basic linear, constant
 * damping coefficient that multiplies the velocity to obtain the force that will oppose the relative
 * velocity between the anchors.
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

#ifndef REAK_DAMPER_HPP
#define REAK_DAMPER_HPP

#include "kte_map.hpp"

#include <ReaK/math/kinetostatics/kinetostatics.hpp>
#include <utility>

namespace ReaK::kte {

/** This class defines a damper acting between two generalized coordinates. */
class damper_gen : public kte_map {
 private:
  std::shared_ptr<gen_coord<double>>
      mAnchor1;  ///< Holds the first generalized coordinate.
  std::shared_ptr<gen_coord<double>>
      mAnchor2;     ///< Holds the second generalized coordinate.
  double mDamping;  ///< The damping coefficient (in Ns/m or Nms/rad).

 public:
  /**
   * Sets the first anchor frame of the damper.
   * \param aPtr A pointer to the new first anchor frame of the damper.
   */
  void setAnchor1(const std::shared_ptr<gen_coord<double>>& aPtr) {
    mAnchor1 = aPtr;
  }
  /**
   * Returns a const-reference to the first anchor frame of the damper.
   * \return A const-reference to the first anchor frame of the damper.
   */
  std::shared_ptr<gen_coord<double>> Anchor1() const { return mAnchor1; }

  /**
   * Sets the first anchor frame of the damper.
   * \param aPtr A pointer to the new first anchor frame of the damper.
   */
  void setAnchor2(const std::shared_ptr<gen_coord<double>>& aPtr) {
    mAnchor2 = aPtr;
  }
  /**
   * Returns a const-reference to the second anchor frame of the damper.
   * \return A const-reference to the second anchor frame of the damper.
   */
  std::shared_ptr<gen_coord<double>> Anchor2() const { return mAnchor2; }

  /**
   * Sets the damping factor of the damper.
   * \param aValue The new damping factor of the damper.
   */
  void setDamping(double aValue) { mDamping = aValue; }
  /**
   * Returns the damping factor of the damper.
   * \return The damping factor of the damper.
   */
  double Damping() const { return mDamping; }

  /**
   * Default constructor.
   */
  explicit damper_gen(const std::string& aName = "")
      : kte_map(aName), mDamping(0.0) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor1 first attach point of the damper.
   * \param aAnchor2 second attach point of the damper.
   * \param aDamping damping coefficient (in Ns/m for a linear generalized coord. or Nms/rad for an angular generalized
   * coord.).
   */
  damper_gen(const std::string& aName,
             std::shared_ptr<gen_coord<double>> aAnchor1,
             std::shared_ptr<gen_coord<double>> aAnchor2, double aDamping)
      : kte_map(aName),
        mAnchor1(std::move(aAnchor1)),
        mAnchor2(std::move(aAnchor2)),
        mDamping(aDamping) {}

  /**
   * Default destructor.
   */
  ~damper_gen() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override;

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override;

  void clearForce() override;

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor1) & RK_SERIAL_SAVE_WITH_NAME(mAnchor2) &
        RK_SERIAL_SAVE_WITH_NAME(mDamping);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor1) & RK_SERIAL_LOAD_WITH_NAME(mAnchor2) &
        RK_SERIAL_LOAD_WITH_NAME(mDamping);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(damper_gen, 0xC2100010, 1, "damper_gen", kte_map)
};

/** This class defines a damper acting between two 2D frames. */
class damper_2D : public kte_map {
 private:
  std::shared_ptr<frame_2D<double>> mAnchor1;  ///< Holds the first 2D frame.
  std::shared_ptr<frame_2D<double>> mAnchor2;  ///< Holds the second 2D frame.
  double mDamping;  ///< The damping coefficient (in Ns/m).

 public:
  /**
   * Sets the first anchor frame of the damper.
   * \param aPtr A pointer to the new first anchor frame of the damper.
   */
  void setAnchor1(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor1 = aPtr;
  }
  /**
   * Returns a const-reference to the first anchor frame of the damper.
   * \return A const-reference to the first anchor frame of the damper.
   */
  std::shared_ptr<frame_2D<double>> Anchor1() const { return mAnchor1; }

  /**
   * Sets the first anchor frame of the damper.
   * \param aPtr A pointer to the new first anchor frame of the damper.
   */
  void setAnchor2(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mAnchor2 = aPtr;
  }
  /**
   * Returns a const-reference to the second anchor frame of the damper.
   * \return A const-reference to the second anchor frame of the damper.
   */
  std::shared_ptr<frame_2D<double>> Anchor2() const { return mAnchor2; }

  /**
   * Sets the damping factor of the damper.
   * \param aValue The new damping factor of the damper.
   */
  void setDamping(double aValue) { mDamping = aValue; }
  /**
   * Returns the damping factor of the damper.
   * \return The damping factor of the damper.
   */
  double Damping() const { return mDamping; }

  /**
   * Default constructor.
   */
  explicit damper_2D(const std::string& aName = "")
      : kte_map(aName), mDamping(0.0) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor1 first attach point of the damper.
   * \param aAnchor2 second attach point of the damper.
   * \param aDamping damping coefficient (in Ns/m).
   */
  damper_2D(const std::string& aName,
            std::shared_ptr<frame_2D<double>> aAnchor1,
            std::shared_ptr<frame_2D<double>> aAnchor2, double aDamping)
      : kte_map(aName),
        mAnchor1(std::move(aAnchor1)),
        mAnchor2(std::move(aAnchor2)),
        mDamping(aDamping) {}

  /**
   * Default destructor.
   */
  ~damper_2D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override;

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override;

  void clearForce() override;

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor1) & RK_SERIAL_SAVE_WITH_NAME(mAnchor2) &
        RK_SERIAL_SAVE_WITH_NAME(mDamping);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor1) & RK_SERIAL_LOAD_WITH_NAME(mAnchor2) &
        RK_SERIAL_LOAD_WITH_NAME(mDamping);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(damper_2D, 0xC2100011, 1, "damper_2D", kte_map)
};

/** This class defines a damper acting between two 3D frames. */
class damper_3D : public kte_map {
 private:
  std::shared_ptr<frame_3D<double>> mAnchor1;  ///< Holds the first 3D frame.
  std::shared_ptr<frame_3D<double>> mAnchor2;  ///< Holds the second 3D frame.
  double mDamping;  ///< The damping coefficient (in Ns/m).

 public:
  /**
   * Sets the first anchor frame of the damper.
   * \param aPtr A pointer to the new first anchor frame of the damper.
   */
  void setAnchor1(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor1 = aPtr;
  }
  /**
   * Returns a const-reference to the first anchor frame of the damper.
   * \return A const-reference to the first anchor frame of the damper.
   */
  std::shared_ptr<frame_3D<double>> Anchor1() const { return mAnchor1; }

  /**
   * Sets the first anchor frame of the damper.
   * \param aPtr A pointer to the new first anchor frame of the damper.
   */
  void setAnchor2(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mAnchor2 = aPtr;
  }
  /**
   * Returns a const-reference to the second anchor frame of the damper.
   * \return A const-reference to the second anchor frame of the damper.
   */
  std::shared_ptr<frame_3D<double>> Anchor2() const { return mAnchor2; }

  /**
   * Sets the damping factor of the damper.
   * \param aValue The new damping factor of the damper.
   */
  void setDamping(double aValue) { mDamping = aValue; }
  /**
   * Returns the damping factor of the damper.
   * \return The damping factor of the damper.
   */
  double Damping() const { return mDamping; }

  /**
   * Default constructor.
   */
  explicit damper_3D(const std::string& aName = "")
      : kte_map(aName), mDamping(0.0) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aAnchor1 first attach point of the damper.
   * \param aAnchor2 second attach point of the damper.
   * \param aDamping damping coefficient (in Ns/m).
   */
  damper_3D(const std::string& aName,
            std::shared_ptr<frame_3D<double>> aAnchor1,
            std::shared_ptr<frame_3D<double>> aAnchor2, double aDamping)
      : kte_map(aName),
        mAnchor1(std::move(aAnchor1)),
        mAnchor2(std::move(aAnchor2)),
        mDamping(aDamping) {}

  /**
   * Default destructor.
   */
  ~damper_3D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override;

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override;

  void clearForce() override;

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mAnchor1) & RK_SERIAL_SAVE_WITH_NAME(mAnchor2) &
        RK_SERIAL_SAVE_WITH_NAME(mDamping);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mAnchor1) & RK_SERIAL_LOAD_WITH_NAME(mAnchor2) &
        RK_SERIAL_LOAD_WITH_NAME(mDamping);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(damper_3D, 0xC2100012, 1, "damper_3D", kte_map)
};

}  // namespace ReaK::kte

#endif

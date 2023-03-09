/**
 * \file line_point_mindist.h
 *
 * This library implements geometric models used to compute the kinematics of an end-frame which
 * should follow the minimum distance from the base-frame to a given (fixed) line in 2D or 3D space.
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

#ifndef REAK_MBD_KTE_LINE_POINT_MINDIST_H_
#define REAK_MBD_KTE_LINE_POINT_MINDIST_H_

#include "ReaK/mbd/kte/kte_map.h"

#include <utility>
#include "ReaK/math/kinetostatics/kinetostatics.h"

namespace ReaK::kte {

/**
 * This class implements the geometric calculation of a 2D end-frame which follows a (fixed) line while
 * maintain the minimum distance to a base-frame.
 */
class line_point_mindist_2D : public kte_map {
 private:
  std::shared_ptr<frame_2D<double>>
      mBase;  ///< Holds the base-frame, or kinematic input, or free-point.
  std::shared_ptr<frame_2D<double>>
      mEnd;  ///< Holds the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
  vect<double, 2>
      mTangent;  ///< Holds the tangent vector of the (fixed) line in global coordinates.
  vect<double, 2>
      mOriginMinDist;  ///< Holds the minimum-distance vector from origin to the line, in global coordinates.

 public:
  /**
   * Sets the base-frame, or kinematic input, or free-point.
   * \param aPtr The new base-frame, or kinematic input, or free-point.
   */
  void setBaseFrame(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mBase = aPtr;
  }
  /**
   * Returns the base-frame, or kinematic input, or free-point.
   * \return The base-frame, or kinematic input, or free-point.
   */
  std::shared_ptr<frame_2D<double>> BaseFrame() const { return mBase; }

  /**
   * Sets the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
   * \param aPtr The new end-frame, or kinematic output, or min-dist point to base-frame, on the line.
   */
  void setEndFrame(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mEnd = aPtr;
  }
  /**
   * Returns the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
   * \return The end-frame, or kinematic output, or min-dist point to base-frame, on the line.
   */
  std::shared_ptr<frame_2D<double>> EndFrame() const { return mEnd; }

  /**
   * Sets the global tangent vector of the line.
   * \param aValue The new global tangent vector of the line.
   */
  void setTangent(const vect<double, 2>& aValue) { mTangent = aValue; }
  /**
   * Returns the global tangent vector of the line.
   * \return The global tangent vector of the line.
   */
  vect<double, 2> Tangent() const { return mTangent; }

  /**
   * Sets the minimum-distance vector from origin to the line, in global coordinates.
   * \param aValue The new minimum-distance vector from origin to the line, in global coordinates.
   */
  void setOriginMinDist(const vect<double, 2>& aValue) {
    mOriginMinDist = aValue;
  }
  /**
   * Returns the minimum-distance vector from origin to the line, in global coordinates.
   * \return The minimum-distance vector from origin to the line, in global coordinates.
   */
  vect<double, 2> OriginMinDist() const { return mOriginMinDist; }

  /**
   * Default constructor.
   */
  explicit line_point_mindist_2D(const std::string& aName = "")
      : kte_map(aName) {}

  /**
   * Parametrized constructor.
   * \param aName the name for the KTE model.
   * \param aBase the base-frame, or kinematic input, or free-point.
   * \param aEnd the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
   * \param aTangent the tangent vector of the (fixed) line in global coordinates.
   * \param aOriginMinDist the minimum-distance vector from origin to the line, in global coordinates.
   */
  line_point_mindist_2D(const std::string& aName,
                        std::shared_ptr<frame_2D<double>> aBase,
                        std::shared_ptr<frame_2D<double>> aEnd,
                        const vect<double, 2>& aTangent,
                        const vect<double, 2>& aOriginMinDist)
      : kte_map(aName),
        mBase(std::move(aBase)),
        mEnd(std::move(aEnd)),
        mTangent(aTangent),
        mOriginMinDist(aOriginMinDist) {}

  /**
   * Default destructor.
   */
  ~line_point_mindist_2D() override = default;

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
    A& RK_SERIAL_SAVE_WITH_NAME(mBase) & RK_SERIAL_SAVE_WITH_NAME(mEnd) &
        RK_SERIAL_SAVE_WITH_NAME(mTangent) &
        RK_SERIAL_SAVE_WITH_NAME(mOriginMinDist);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mBase) & RK_SERIAL_LOAD_WITH_NAME(mEnd) &
        RK_SERIAL_LOAD_WITH_NAME(mTangent) &
        RK_SERIAL_LOAD_WITH_NAME(mOriginMinDist);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(line_point_mindist_2D, 0xC2100031, 1,
                              "line_point_mindist_2D", kte_map)
};

/**
 * This class implements the geometric calculation of a 3D end-frame which follows a (fixed) line while
 * maintain the minimum distance to a base-frame.
 */
class line_point_mindist_3D : public kte_map {
 private:
  std::shared_ptr<frame_3D<double>>
      mBase;  ///< Holds the base-frame, or kinematic input, or free-point.
  std::shared_ptr<frame_3D<double>>
      mEnd;  ///< Holds the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
  vect<double, 3>
      mTangent;  ///< Holds the tangent vector of the (fixed) line in global coordinates.
  vect<double, 3>
      mOriginMinDist;  ///< Holds the minimum-distance vector from origin to the line, in global coordinates.

 public:
  /**
   * Sets the base-frame, or kinematic input, or free-point.
   * \param aPtr The new base-frame, or kinematic input, or free-point.
   */
  void setBaseFrame(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mBase = aPtr;
  }
  /**
   * Returns the base-frame, or kinematic input, or free-point.
   * \return The base-frame, or kinematic input, or free-point.
   */
  std::shared_ptr<frame_3D<double>> BaseFrame() const { return mBase; }

  /**
   * Sets the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
   * \param aPtr The new end-frame, or kinematic output, or min-dist point to base-frame, on the line.
   */
  void setEndFrame(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mEnd = aPtr;
  }
  /**
   * Returns the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
   * \return The end-frame, or kinematic output, or min-dist point to base-frame, on the line.
   */
  std::shared_ptr<frame_3D<double>> EndFrame() const { return mEnd; }

  /**
   * Sets the global tangent vector of the line.
   * \param aValue The new global tangent vector of the line.
   */
  void setTangent(const vect<double, 3>& aValue) { mTangent = aValue; }
  /**
   * Returns the global tangent vector of the line.
   * \return The global tangent vector of the line.
   */
  vect<double, 3> Tangent() const { return mTangent; }

  /**
   * Sets the minimum-distance vector from origin to the line, in global coordinates.
   * \param aValue The new minimum-distance vector from origin to the line, in global coordinates.
   */
  void setOriginMinDist(const vect<double, 3>& aValue) {
    mOriginMinDist = aValue;
  }
  /**
   * Returns the minimum-distance vector from origin to the line, in global coordinates.
   * \return The minimum-distance vector from origin to the line, in global coordinates.
   */
  vect<double, 3> OriginMinDist() const { return mOriginMinDist; }

  /**
   * Default constructor.
   */
  explicit line_point_mindist_3D(const std::string& aName = "")
      : kte_map(aName) {}

  /**
   * Parametrized constructor.
   * \param aName the name for the KTE model.
   * \param aBase the base-frame, or kinematic input, or free-point.
   * \param aEnd the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
   * \param aTangent the tangent vector of the (fixed) line in global coordinates.
   * \param aOriginMinDist the minimum-distance vector from origin to the line, in global coordinates.
   */
  line_point_mindist_3D(const std::string& aName,
                        std::shared_ptr<frame_3D<double>> aBase,
                        std::shared_ptr<frame_3D<double>> aEnd,
                        const vect<double, 3>& aTangent,
                        const vect<double, 3>& aOriginMinDist)
      : kte_map(aName),
        mBase(std::move(aBase)),
        mEnd(std::move(aEnd)),
        mTangent(aTangent),
        mOriginMinDist(aOriginMinDist) {}

  /**
   * Default destructor.
   */
  ~line_point_mindist_3D() override = default;

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
    A& RK_SERIAL_SAVE_WITH_NAME(mBase) & RK_SERIAL_SAVE_WITH_NAME(mEnd) &
        RK_SERIAL_SAVE_WITH_NAME(mTangent) &
        RK_SERIAL_SAVE_WITH_NAME(mOriginMinDist);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mBase) & RK_SERIAL_LOAD_WITH_NAME(mEnd) &
        RK_SERIAL_LOAD_WITH_NAME(mTangent) &
        RK_SERIAL_LOAD_WITH_NAME(mOriginMinDist);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(line_point_mindist_3D, 0xC2100032, 1,
                              "line_point_mindist_3D", kte_map)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_LINE_POINT_MINDIST_H_

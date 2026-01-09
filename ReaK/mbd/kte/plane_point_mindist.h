/**
 * \file plane_point_mindist.h
 *
 * This library implements geometric models used to compute the kinematics of an end-frame which
 * should follow the minimum distance from the base-frame to a given (fixed) plance in 3D space.
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

#ifndef REAK_MBD_KTE_PLANE_POINT_MINDIST_H_
#define REAK_MBD_KTE_PLANE_POINT_MINDIST_H_

#include "ReaK/mbd/kte/kte_map.h"

#include <utility>
#include "ReaK/math/kinetostatics/kinetostatics.h"

namespace ReaK::kte {

/**
 * This class implements the geometric calculation of a 3D end-frame which follows a (fixed) plane while
 * maintain the minimum distance to a base-frame.
 */
class plane_point_mindist_3D : public kte_map {
 private:
  std::shared_ptr<frame_3D<double>>
      mBase;  ///< Holds the base-frame, or kinematic input, or free-point.
  std::shared_ptr<frame_3D<double>>
      mEnd;  ///< Holds the end-frame, or kinematic output, or min-dist point to base-frame, on the plane.
  vect<double, 3>
      mNormal;  ///< Holds the normal vector of the plane, in global coordinates.
  double
      mOrigin;  ///< Holds the distance to the origin, i.e., mOrigin * mNormal is the vector from the plane to a plane
                ///parallel and intersecting the origin.

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
   * Sets the end-frame, or kinematic output, or min-dist point to base-frame, on the plane.
   * \param aPtr The new end-frame, or kinematic output, or min-dist point to base-frame, on the plane.
   */
  void setEndFrame(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mEnd = aPtr;
  }
  /**
   * Returns the end-frame, or kinematic output, or min-dist point to base-frame, on the plane.
   * \return The end-frame, or kinematic output, or min-dist point to base-frame, on the plane.
   */
  std::shared_ptr<frame_3D<double>> EndFrame() const { return mEnd; }

  /**
   * Sets the normal vector of the plane, in global coordinates.
   * \param aValue The new normal vector of the plane, in global coordinates.
   */
  void setNormal(const vect<double, 3>& aValue) { mNormal = aValue; }
  /**
   * Returns the normal vector of the plane, in global coordinates.
   * \return The normal vector of the plane, in global coordinates.
   */
  vect<double, 3> Normal() const { return mNormal; }

  /**
   * Sets the distance to the origin, i.e., aOrigin * aNormal is the vector from the plane to a plane parallel and
   * intersecting the origin.
   * \param aValue The new distance to the origin, i.e., aOrigin * aNormal is the vector from the plane to a plane
   * parallel and intersecting the origin.
   */
  void setOrigin(double& aValue) { mOrigin = aValue; }
  /**
   * Returns the distance to the origin, i.e., aOrigin * aNormal is the vector from the plane to a plane parallel and
   * intersecting the origin.
   * \return The distance to the origin, i.e., aOrigin * aNormal is the vector from the plane to a plane parallel and
   * intersecting the origin.
   */
  double Origin() const { return mOrigin; }

  /**
   * Default constructor.
   */
  explicit plane_point_mindist_3D(const std::string& aName = "")
      : kte_map(aName), mOrigin(0.0) {}

  /**
   * Parametrized constructor.
   * \param aName the name of the KTE model.
   * \param aBase the base-frame, or kinematic input, or free-point.
   * \param aEnd the end-frame, or kinematic output, or min-dist point to base-frame, on the plane.
   * \param aNormal the normal vector of the plane, in global coordinates.
   * \param aOrigin the distance to the origin, i.e., aOrigin * aNormal is the vector from the plane to a plane parallel
   * and intersecting the origin.
   */
  plane_point_mindist_3D(const std::string& aName,
                         std::shared_ptr<frame_3D<double>> aBase,
                         std::shared_ptr<frame_3D<double>> aEnd,
                         const vect<double, 3>& aNormal, double aOrigin)
      : kte_map(aName),
        mBase(std::move(aBase)),
        mEnd(std::move(aEnd)),
        mNormal(aNormal),
        mOrigin(aOrigin) {}

  /**
   * Default destructor.
   */
  ~plane_point_mindist_3D() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override;

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override;

  void clearForce() override;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(mBase) & RK_SERIAL_SAVE_WITH_NAME(mEnd) &
        RK_SERIAL_SAVE_WITH_NAME(mNormal) & RK_SERIAL_SAVE_WITH_NAME(mOrigin);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(mBase) & RK_SERIAL_LOAD_WITH_NAME(mEnd) &
        RK_SERIAL_LOAD_WITH_NAME(mNormal) & RK_SERIAL_LOAD_WITH_NAME(mOrigin);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(plane_point_mindist_3D, 0xC2100030, 1,
                              "plane_point_mindist_3D", kte_map)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_PLANE_POINT_MINDIST_H_

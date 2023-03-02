/**
 * \file geometry_3D.hpp
 *
 * This library declares the base-class for all 3D geometric objects (renderable).
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date April 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_GEOMETRY_3D_HPP
#define REAK_GEOMETRY_3D_HPP

#include <ReaK/core/base/named_object.hpp>

#include <ReaK/math/kinetostatics/pose_3D.hpp>

namespace ReaK::geom {

/** This class is the base-class for all 3D geometric objects (renderable). */
class geometry_3D : public named_object {
 protected:
  std::shared_ptr<pose_3D<double>> mAnchor;
  pose_3D<double> mPose;

 public:
  /**
   * This function returns the anchor of the geometry.
   * \return A shared-pointer to the anchor pose object.
   */
  const std::shared_ptr<pose_3D<double>>& getAnchor() const { return mAnchor; }
  /**
   * This function sets the anchor of the geometry.
   * \param aAnchor A shared-pointer to the new anchor pose object.
   */
  void setAnchor(const std::shared_ptr<pose_3D<double>>& aAnchor);

  /**
   * This function returns the pose of the geometry.
   * \return The pose object.
   */
  const pose_3D<double>& getPose() const { return mPose; }
  /**
   * This function sets the pose of the geometry.
   * \param aPose The new pose object.
   */
  void setPose(const pose_3D<double>& aPose);

  /**
   * Default constructor.
   * \param aName The name of the object.
   * \param aAnchor The anchor object for the geometry.
   * \param aPose The pose of the geometry (relative to the anchor).
   */
  explicit geometry_3D(const std::string& aName = "",
                       std::shared_ptr<pose_3D<double>> aAnchor =
                           std::shared_ptr<pose_3D<double>>(),
                       const pose_3D<double>& aPose = pose_3D<double>());

  /**
   * Default destructor.
   */
  ~geometry_3D() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override;

  void load(ReaK::serialization::iarchive& A, unsigned int /*unused*/) override;

  RK_RTTI_MAKE_CONCRETE_1BASE(geometry_3D, 0xC3100001, 1, "geometry_3D",
                              named_object)
};

}  // namespace ReaK::geom

#endif

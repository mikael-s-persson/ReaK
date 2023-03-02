/**
 * \file capped_cylinder.hpp
 *
 * This library declares a class to represent capped cylinders.
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

#ifndef REAK_CAPPED_CYLINDER_HPP
#define REAK_CAPPED_CYLINDER_HPP

#include "shape_3D.hpp"

namespace ReaK::geom {

/** This class represents a cylinder aligned along the z-axis of its center pose. */
class capped_cylinder : public shape_3D {
 protected:
  double mLength;
  double mRadius;

 public:
  /**
   * This function returns the maximum radius of the shape (radius of the sphere that bounds the shape).
   * \return The maximum radius of the shape.
   */
  double getBoundingRadius() const override;

  /**
   * This function returns the length of the cylinder.
   * \return The length of the cylinder.
   */
  double getLength() const { return mLength; }
  /**
   * This function sets the length of the cylinder.
   * \param aLength The new length of the cylinder.
   */
  void setLength(double aLength) { mLength = aLength; }

  /**
   * This function returns the radius of the cylinder.
   * \return The radius of the cylinder.
   */
  double getRadius() const { return mRadius; }
  /**
   * This function sets the radius of the cylinder.
   * \param aRadius The new radius of the cylinder.
   */
  void setRadius(double aRadius) { mRadius = aRadius; }

  /**
   * Default constructor.
   * \param aName The name of the object.
   * \param aAnchor The anchor object for the geometry.
   * \param aPose The pose of the geometry (relative to the anchor).
   * \param aLength The length of the cylinder.
   * \param aRadius The radius of the cylinder.
   */
  explicit capped_cylinder(const std::string& aName = "",
                           const std::shared_ptr<pose_3D<double>>& aAnchor =
                               std::shared_ptr<pose_3D<double>>(),
                           const pose_3D<double>& aPose = pose_3D<double>(),
                           double aLength = 1.0, double aRadius = 1.0);

  /**
   * Default destructor.
   */
  ~capped_cylinder() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override;

  void load(ReaK::serialization::iarchive& A, unsigned int /*unused*/) override;

  RK_RTTI_MAKE_CONCRETE_1BASE(capped_cylinder, 0xC3100011, 1, "capped_cylinder",
                              shape_3D)
};

}  // namespace ReaK::geom

#endif

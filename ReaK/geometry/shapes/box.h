/**
 * \file box.h
 *
 * This library declares a class to represent boxes.
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

#ifndef REAK_GEOMETRY_SHAPES_BOX_H_
#define REAK_GEOMETRY_SHAPES_BOX_H_

#include "ReaK/geometry/shapes/shape_3D.h"

namespace ReaK::geom {

/** This class represents a box aligned along its center pose. */
class box : public shape_3D {
 protected:
  vect<double, 3> mDimensions;

 public:
  /**
   * This function returns the maximum radius of the shape (radius of the sphere that bounds the shape).
   * \return The maximum radius of the shape.
   */
  double getBoundingRadius() const override;

  /**
   * This function returns the dimensions of the box.
   * \return The dimensions of the box.
   */
  const vect<double, 3>& getDimensions() const { return mDimensions; }
  /**
   * This function sets the dimensions of the box.
   * \param aDimensions The new dimensions of the box.
   */
  void setDimensions(const vect<double, 3>& aDimensions) {
    mDimensions = aDimensions;
  }

  /**
   * Default constructor.
   * \param aName The name of the object.
   * \param aAnchor The anchor object for the geometry.
   * \param aPose The pose of the geometry (relative to the anchor).
   * \param aDimensions The dimensions of the box.
   */
  explicit box(const std::string& aName = "",
               const std::shared_ptr<pose_3D<double>>& aAnchor =
                   std::shared_ptr<pose_3D<double>>(),
               const pose_3D<double>& aPose = pose_3D<double>(),
               const vect<double, 3>& aDimensions = (vect<double, 3>(1.0, 1.0,
                                                                     1.0)));

  /**
   * Default destructor.
   */
  ~box() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override;

  void load(ReaK::serialization::iarchive& A, unsigned int /*unused*/) override;

  RK_RTTI_MAKE_CONCRETE_1BASE(box, 0xC3100013, 1, "box", shape_3D)
};

}  // namespace ReaK::geom

#endif  // REAK_GEOMETRY_SHAPES_BOX_H_

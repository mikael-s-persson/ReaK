/**
 * \file coord_arrows_3D.h
 *
 * This library declares a class to render 3D coordinate arrows.
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

#ifndef REAK_GEOMETRY_SHAPES_COORD_ARROWS_3D_H_
#define REAK_GEOMETRY_SHAPES_COORD_ARROWS_3D_H_

#include "ReaK/geometry/shapes/geometry_3D.h"

namespace ReaK::geom {

/** This class defines a class to render 3D coordinate arrows. */
class coord_arrows_3D : public geometry_3D {
 protected:
  double mArrowLength;

 public:
  /**
   * This function returns the length of the arrows representing the coordinate axes.
   * \return The length of the arrows representing the coordinate axes.
   */
  double getArrowLength() const { return mArrowLength; }
  /**
   * This function sets the new length of the arrows representing the coordinate axes.
   * \param aArrowLength The new length of the arrows representing the coordinate axes.
   */
  void setArrowLength(double aArrowLength) { mArrowLength = aArrowLength; }

  /**
   * Default constructor.
   * \param aName The name of the object.
   * \param aAnchor The anchor object for the geometry.
   * \param aPose The pose of the geometry (relative to the anchor).
   * \param aArrowLength The length of the arrows representing the coordinate axes.
   */
  explicit coord_arrows_3D(const std::string& aName = "",
                           const std::shared_ptr<pose_3D<double>>& aAnchor =
                               std::shared_ptr<pose_3D<double>>(),
                           const pose_3D<double>& aPose = pose_3D<double>(),
                           double aArrowLength = 1.0);

  /**
   * Default destructor.
   */
  ~coord_arrows_3D() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override;

  void load(ReaK::serialization::iarchive& A, unsigned int /*unused*/) override;

  RK_RTTI_MAKE_CONCRETE_1BASE(coord_arrows_3D, 0xC3100003, 1, "coord_arrows_3D",
                              geometry_3D)
};

}  // namespace ReaK::geom

#endif  // REAK_GEOMETRY_SHAPES_COORD_ARROWS_3D_H_

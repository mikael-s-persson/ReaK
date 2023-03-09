/**
 * \file prox_plane_cylinder.h
 *
 * This library declares a class for proximity queries between a plane and a cylinder.
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

#ifndef REAK_GEOMETRY_PROXIMITY_PROX_PLANE_CYLINDER_H_
#define REAK_GEOMETRY_PROXIMITY_PROX_PLANE_CYLINDER_H_

#include "ReaK/geometry/proximity/proximity_finder_3D.h"

#include "ReaK/geometry/shapes/cylinder.h"
#include "ReaK/geometry/shapes/plane.h"

namespace ReaK::geom {

proximity_record_3D compute_proximity(const plane& aPlane,
                                      const shape_3D_precompute_pack& aPack1,
                                      const cylinder& aCylinder,
                                      const shape_3D_precompute_pack& aPack2);

proximity_record_3D compute_proximity(const cylinder& aCylinder,
                                      const shape_3D_precompute_pack& aPack1,
                                      const plane& aPlane,
                                      const shape_3D_precompute_pack& aPack2);

/**
 * This class is for proximity queries between a plane and a cylinder.
 */
class prox_plane_cylinder : public proximity_finder_3D {
 protected:
  const plane* mPlane;
  const cylinder* mCylinder;

 public:
  /** This function performs the proximity query on its associated shapes. */
  proximity_record_3D computeProximity(
      const shape_3D_precompute_pack& aPack1,
      const shape_3D_precompute_pack& aPack2) override;

  /**
   * Default constructor.
   * \param aPlane The plane involved in the proximity query.
   * \param aCylinder The cylinder involved in the proximity query.
   */
  explicit prox_plane_cylinder(const plane* aPlane = nullptr,
                               const cylinder* aCylinder = nullptr);

  /** Destructor. */
  ~prox_plane_cylinder() override = default;
};

}  // namespace ReaK::geom

#endif  // REAK_GEOMETRY_PROXIMITY_PROX_PLANE_CYLINDER_H_

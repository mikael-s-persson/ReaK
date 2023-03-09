/**
 * \file prox_plane_plane.h
 *
 * This library declares a class for proximity queries between two planes.
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

#ifndef REAK_GEOMETRY_PROXIMITY_PROX_PLANE_PLANE_H_
#define REAK_GEOMETRY_PROXIMITY_PROX_PLANE_PLANE_H_

#include "ReaK/geometry/proximity/proximity_finder_3D.h"

#include "ReaK/geometry/shapes/plane.h"

namespace ReaK::geom {

void compute_proximity_of_point(const plane& aPlane,
                                const pose_3D<double>& aPlGblPose,
                                const vect<double, 3>& aPoint,
                                vect<double, 3>& aPointRec, double& aDistance);

proximity_record_3D compute_proximity(const plane& aPlane1,
                                      const shape_3D_precompute_pack& aPack1,
                                      const plane& aPlane2,
                                      const shape_3D_precompute_pack& aPack2);

/**
 * This class is for proximity queries between two planes.
 */
class prox_plane_plane : public proximity_finder_3D {
 protected:
  const plane* mPlane1;
  const plane* mPlane2;

 public:
  /** This function performs the proximity query on its associated shapes. */
  proximity_record_3D computeProximity(
      const shape_3D_precompute_pack& aPack1,
      const shape_3D_precompute_pack& aPack2) override;

  /**
   * Default constructor.
   * \param aPlane1 The first plane involved in the proximity query.
   * \param aPlane2 The second plane involved in the proximity query.
   */
  explicit prox_plane_plane(const plane* aPlane1 = nullptr,
                            const plane* aPlane2 = nullptr);

  /** Destructor. */
  ~prox_plane_plane() override = default;
};

}  // namespace ReaK::geom

#endif  // REAK_GEOMETRY_PROXIMITY_PROX_PLANE_PLANE_H_

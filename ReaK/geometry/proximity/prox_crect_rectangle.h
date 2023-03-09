/**
 * \file prox_crect_rectangle.h
 *
 * This library declares a class for proximity queries between a capped rectangle and a rectangle.
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

#ifndef REAK_GEOMETRY_PROXIMITY_PROX_CRECT_RECTANGLE_H_
#define REAK_GEOMETRY_PROXIMITY_PROX_CRECT_RECTANGLE_H_

#include "ReaK/geometry/proximity/proximity_finder_2D.h"

#include "ReaK/geometry/shapes/capped_rectangle.h"
#include "ReaK/geometry/shapes/rectangle.h"

namespace ReaK::geom {

proximity_record_2D compute_proximity_of_line(const rectangle& aRectangle,
                                              const pose_2D<double>& aGblPose,
                                              const vect<double, 2>& ln_c,
                                              const vect<double, 2>& ln_t,
                                              double half_length);

proximity_record_2D compute_proximity(const capped_rectangle& aCRect,
                                      const shape_2D_precompute_pack& aPack1,
                                      const rectangle& aRectangle,
                                      const shape_2D_precompute_pack& aPack2);

proximity_record_2D compute_proximity(const rectangle& aRectangle,
                                      const shape_2D_precompute_pack& aPack1,
                                      const capped_rectangle& aCRect,
                                      const shape_2D_precompute_pack& aPack2);

/**
 * This class is for proximity queries between a capped rectangle and a rectangle.
 */
class prox_crect_rectangle : public proximity_finder_2D {
 protected:
  const capped_rectangle* mCRect;
  const rectangle* mRectangle;

 public:
  /** This function performs the proximity query on its associated shapes. */
  proximity_record_2D computeProximity(
      const shape_2D_precompute_pack& aPack1,
      const shape_2D_precompute_pack& aPack2) override;

  /**
   * Default constructor.
   * \param aCRect The first capped rectangle involved in the proximity query.
   * \param aRectangle The second capped rectangle involved in the proximity query.
   */
  explicit prox_crect_rectangle(const capped_rectangle* aCRect = nullptr,
                                const rectangle* aRectangle = nullptr);

  /** Destructor. */
  ~prox_crect_rectangle() override = default;
};

}  // namespace ReaK::geom

#endif  // REAK_GEOMETRY_PROXIMITY_PROX_CRECT_RECTANGLE_H_

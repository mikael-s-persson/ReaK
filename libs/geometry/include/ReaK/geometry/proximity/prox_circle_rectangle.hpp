/**
 * \file prox_circle_rectangle.hpp
 *
 * This library declares a class for proximity queries between a circle and a rectangle.
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

#ifndef REAK_PROX_CIRCLE_RECTANGLE_HPP
#define REAK_PROX_CIRCLE_RECTANGLE_HPP

#include "ReaK/geometry/proximity/proximity_finder_2D.hpp"

#include "ReaK/geometry/shapes/circle.hpp"
#include "ReaK/geometry/shapes/rectangle.hpp"

namespace ReaK::geom {

proximity_record_2D compute_proximity(const circle& aCircle,
                                      const shape_2D_precompute_pack& aPack1,
                                      const rectangle& aRectangle,
                                      const shape_2D_precompute_pack& aPack2);

proximity_record_2D compute_proximity(const rectangle& aRectangle,
                                      const shape_2D_precompute_pack& aPack1,
                                      const circle& aCircle,
                                      const shape_2D_precompute_pack& aPack2);

/**
 * This class is for proximity queries between a circle and a rectangle.
 */
class prox_circle_rectangle : public proximity_finder_2D {
 protected:
  const circle* mCircle;
  const rectangle* mRectangle;

 public:
  /** This function performs the proximity query on its associated shapes. */
  proximity_record_2D computeProximity(
      const shape_2D_precompute_pack& aPack1,
      const shape_2D_precompute_pack& aPack2) override;

  /**
   * Default constructor.
   * \param aCircle The circle involved in the proximity query.
   * \param aRectangle The rectangle involved in the proximity query.
   */
  explicit prox_circle_rectangle(const circle* aCircle = nullptr,
                                 const rectangle* aRectangle = nullptr);

  /** Destructor. */
  ~prox_circle_rectangle() override = default;
};

}  // namespace ReaK::geom

#endif

/**
 * \file proximity_finder_2D.h
 *
 * This library declares the base-class for all proximity finders (that perform the proximity queries) between 2D
 *shapes.
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

#ifndef REAK_GEOMETRY_PROXIMITY_PROXIMITY_FINDER_2D_H_
#define REAK_GEOMETRY_PROXIMITY_PROXIMITY_FINDER_2D_H_

#include "ReaK/geometry/shapes/shape_2D.h"

#include "ReaK/geometry/proximity/proximity_record_2D.h"

namespace ReaK::geom {

/**
 * This class is the base-class for a proximity query with 2D shapes.
 */
class proximity_finder_2D {
 public:
  /** This function performs the proximity query on its associated shapes. */
  virtual proximity_record_2D computeProximity(
      const shape_2D_precompute_pack& aPack1,
      const shape_2D_precompute_pack& aPack2) = 0;

  /** Default constructor. */
  proximity_finder_2D() = default;

  /** Destructor. */
  virtual ~proximity_finder_2D() = default;
};

}  // namespace ReaK::geom

#endif  // REAK_GEOMETRY_PROXIMITY_PROXIMITY_FINDER_2D_H_

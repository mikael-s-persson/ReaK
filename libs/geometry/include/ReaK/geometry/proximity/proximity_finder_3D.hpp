/**
 * \file proximity_finder_3D.hpp
 *
 * This library declares the base-class for all proximity finders (that perform the proximity queries) between 3D
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

#ifndef REAK_PROXIMITY_FINDER_3D_HPP
#define REAK_PROXIMITY_FINDER_3D_HPP

#include "ReaK/geometry/shapes/shape_3D.hpp"

#include "ReaK/geometry/proximity/proximity_record_3D.hpp"

namespace ReaK::geom {

/**
 * This class is the base-class for a proximity query with 3D shapes.
 */
class proximity_finder_3D {
 public:
  /** This function performs the proximity query on its associated shapes. */
  virtual proximity_record_3D computeProximity(
      const shape_3D_precompute_pack& aPack1,
      const shape_3D_precompute_pack& aPack2) = 0;

  /** Default constructor. */
  proximity_finder_3D() = default;

  /** Destructor. */
  virtual ~proximity_finder_3D() = default;
};

}  // namespace ReaK::geom

#endif

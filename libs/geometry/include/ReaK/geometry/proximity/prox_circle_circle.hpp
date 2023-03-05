/**
 * \file prox_circle_circle.hpp
 *
 * This library declares a class for proximity queries between circles.
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

#ifndef REAK_PROX_CIRCLE_CIRCLE_HPP
#define REAK_PROX_CIRCLE_CIRCLE_HPP

#include "ReaK/geometry/proximity/proximity_finder_2D.hpp"

#include "ReaK/geometry/shapes/circle.hpp"

namespace ReaK::geom {

proximity_record_2D compute_proximity(const circle& aCircle1,
                                      const shape_2D_precompute_pack& aPack1,
                                      const circle& aCircle2,
                                      const shape_2D_precompute_pack& aPack2);

/**
 * This class is for proximity queries between circles.
 */
class prox_circle_circle : public proximity_finder_2D {
 protected:
  const circle* mCircle1;
  const circle* mCircle2;

 public:
  /** This function performs the proximity query on its associated shapes. */
  proximity_record_2D computeProximity(
      const shape_2D_precompute_pack& aPack1,
      const shape_2D_precompute_pack& aPack2) override;

  /**
   * Default constructor.
   * \param aCircle1 The first circle involved in the proximity query.
   * \param aCircle2 The second circle involved in the proximity query.
   */
  explicit prox_circle_circle(const circle* aCircle1 = nullptr,
                              const circle* aCircle2 = nullptr);

  /** Destructor. */
  ~prox_circle_circle() override = default;
};

}  // namespace ReaK::geom

#endif

/**
 * \file prox_crect_crect.hpp
 *
 * This library declares a class for proximity queries between two capped rectangles.
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

#ifndef REAK_PROX_CRECT_CRECT_HPP
#define REAK_PROX_CRECT_CRECT_HPP

#include "proximity_finder_2D.hpp"

#include <ReaK/geometry/shapes/capped_rectangle.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


proximity_record_2D compute_proximity( const capped_rectangle& aCRect1, const shape_2D_precompute_pack& aPack1,
                                       const capped_rectangle& aCRect2, const shape_2D_precompute_pack& aPack2 );

/**
 * This class is for proximity queries between two capped rectangles.
 */
class prox_crect_crect : public proximity_finder_2D {
protected:
  const capped_rectangle* mCRect1;
  const capped_rectangle* mCRect2;

public:
  /** This function performs the proximity query on its associated shapes. */
  virtual proximity_record_2D computeProximity( const shape_2D_precompute_pack& aPack1,
                                                const shape_2D_precompute_pack& aPack2 );

  /**
   * Default constructor.
   * \param aCRect1 The first capped rectangle involved in the proximity query.
   * \param aCRect2 The second capped rectangle involved in the proximity query.
   */
  prox_crect_crect( const capped_rectangle* aCRect1 = nullptr, const capped_rectangle* aCRect2 = nullptr );

  /** Destructor. */
  virtual ~prox_crect_crect(){};
};
};
};

#endif

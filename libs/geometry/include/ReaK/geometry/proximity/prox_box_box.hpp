/**
 * \file .hpp
 *
 * This library declares
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


#ifndef REAK_PROX_BOX_BOX_HPP
#define REAK_PROX_BOX_BOX_HPP

#include "proximity_finder_3D.hpp"

#include <ReaK/geometry/shapes/box.hpp>


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


proximity_record_3D compute_proximity( const box& aBox1, const shape_3D_precompute_pack& aPack1, const box& aBox2,
                                       const shape_3D_precompute_pack& aPack2 );


/**
 * This class is for proximity queries between a box and a box.
 */
class prox_box_box : public proximity_finder_3D {
protected:
  const box* mBox1;
  const box* mBox2;

public:
  /** This function performs the proximity query on its associated shapes. */
  virtual proximity_record_3D computeProximity( const shape_3D_precompute_pack& aPack1,
                                                const shape_3D_precompute_pack& aPack2 );

  /**
   * Default constructor.
   * \param aBox1 The first box involved in the proximity query.
   * \param aBox2 The second box involved in the proximity query.
   */
  prox_box_box( const box* aBox1 = nullptr, const box* aBox2 = nullptr );

  /** Destructor. */
  virtual ~prox_box_box(){};
};
};
};

#endif

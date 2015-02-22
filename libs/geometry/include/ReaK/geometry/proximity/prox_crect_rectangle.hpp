/**
 * \file prox_crect_rectangle.hpp
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

#ifndef REAK_PROX_CRECT_RECTANGLE_HPP
#define REAK_PROX_CRECT_RECTANGLE_HPP

#include "proximity_finder_2D.hpp"

#include <ReaK/geometry/shapes/capped_rectangle.hpp>
#include <ReaK/geometry/shapes/rectangle.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/**
 * This class is for proximity queries between a capped rectangle and a rectangle.
 */
class prox_crect_rectangle : public proximity_finder_2D {
  protected:
    
    const capped_rectangle* mCRect;
    const rectangle* mRectangle;
    
    static void computeProximityOfLine(const rectangle&, const pose_2D<double>&, 
                                       const vect<double,2>&, const vect<double,2>&, 
                                       double, proximity_record_2D&);
    
  public:
    
    /** This function performs the proximity query on its associated shapes. */
    virtual void computeProximity(const shape_2D_precompute_pack& aPack1, 
                                  const shape_2D_precompute_pack& aPack2);
    
    /** 
     * Default constructor.
     * \param aCRect The first capped rectangle involved in the proximity query.
     * \param aRectangle The second capped rectangle involved in the proximity query.
     */
    prox_crect_rectangle(const capped_rectangle* aCRect = NULL,
                         const rectangle* aRectangle = NULL);
    
    /** Destructor. */
    virtual ~prox_crect_rectangle() { };
    
};


};

};

#endif


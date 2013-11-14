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

#include "proximity_finder_2D.hpp"

#include "shapes/circle.hpp"
#include "shapes/rectangle.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/**
 * This class is for proximity queries between a circle and a rectangle.
 */
class prox_circle_rectangle : public proximity_finder_2D {
  protected:
    
    shared_ptr< circle > mCircle;
    shared_ptr< rectangle > mRectangle;
    
  public:
    
    /** Returns the first shape involved in the proximity query. */
    virtual shared_ptr< shape_2D > getShape1() const;
    /** Returns the second shape involved in the proximity query. */
    virtual shared_ptr< shape_2D > getShape2() const;
    
    /** This function performs the proximity query on its associated shapes. */
    virtual void computeProximity();
    
    /** 
     * Default constructor.
     * \param aCircle The circle involved in the proximity query.
     * \param aRectangle The rectangle involved in the proximity query.
     */
    prox_circle_rectangle(const shared_ptr< circle >& aCircle = shared_ptr< circle >(),
                          const shared_ptr< rectangle >& aRectangle = shared_ptr< rectangle >());
    
    /** Destructor. */
    virtual ~prox_circle_rectangle() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const;
    
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int);
    
    RK_RTTI_MAKE_CONCRETE_1BASE(prox_circle_rectangle,0xC3200007,1,"prox_circle_rectangle",proximity_finder_2D)
    
};


};

};

#endif











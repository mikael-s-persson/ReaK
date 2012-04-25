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

#include "proximity_finder_2D.hpp"

#include "shapes/circle.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/**
 * This class is for proximity queries between circles.
 */
class prox_circle_circle : public proximity_finder_2D {
  protected:
    
    shared_ptr< circle > mCircle1;
    shared_ptr< circle > mCircle2;
    
  public:
    
    /** Returns the first shape involved in the proximity query. */
    virtual shared_ptr< shape_2D > getShape1() const;
    /** Returns the second shape involved in the proximity query. */
    virtual shared_ptr< shape_2D > getShape2() const;
    
    /** This function performs the proximity query on its associated shapes. */
    virtual void computeProximity();
    
    /** 
     * Default constructor.
     * \param aCircle1 The first circle involved in the proximity query.
     * \param aCircle2 The second circle involved in the proximity query.
     */
    prox_circle_circle(const shared_ptr< circle >& aCircle1 = shared_ptr< circle >(),
                       const shared_ptr< circle >& aCircle2 = shared_ptr< circle >());
    
    /** Destructor. */
    virtual ~prox_circle_circle() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const;
    
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int);
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(prox_circle_circle,0xC3200005,1,"prox_circle_circle",proximity_finder_2D)
    
};


};

};

#endif











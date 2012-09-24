/**
 * \file prox_ccylinder_box.hpp
 *
 * This library declares a class for proximity queries between a capped cylinder and a box.
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

#ifndef REAK_PROX_CCYLINDER_BOX_HPP
#define REAK_PROX_CCYLINDER_BOX_HPP

#include "proximity_finder_3D.hpp"

#include "shapes/capped_cylinder.hpp"
#include "shapes/box.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/**
 * This class is for proximity queries between a capped cylinder and a box.
 */
class prox_ccylinder_box : public proximity_finder_3D {
  protected:
    
    shared_ptr< capped_cylinder > mCCylinder;
    shared_ptr< box > mBox;
    
    static void computeProximityOfLine(const shared_ptr< box >&, const vect<double,3>&, const vect<double,3>&, double, proximity_record_3D&);
    
  public:
    
    /** Returns the first shape involved in the proximity query. */
    virtual shared_ptr< shape_3D > getShape1() const;
    /** Returns the second shape involved in the proximity query. */
    virtual shared_ptr< shape_3D > getShape2() const;
    
    /** This function performs the proximity query on its associated shapes. */
    virtual void computeProximity();
    
    /** 
     * Default constructor. 
     * \param aCCylinder The capped cylinder involved in the proximity query.
     * \param aBox The box involved in the proximity query.
     */
    prox_ccylinder_box(const shared_ptr< capped_cylinder >& aCCylinder = shared_ptr< capped_cylinder >(),
                       const shared_ptr< box >& aBox = shared_ptr< box >());
    
    /** Destructor. */
    virtual ~prox_ccylinder_box() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const;
    
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int);
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(prox_ccylinder_box,0xC3200016,1,"prox_ccylinder_box",proximity_finder_3D)
    
};


};

};

#endif











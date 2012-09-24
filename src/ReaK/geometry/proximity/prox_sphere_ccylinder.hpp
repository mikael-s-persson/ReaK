/**
 * \file prox_sphere_ccylinder.hpp
 *
 * This library declares a class for proximity queries between a sphere and a capped cylinder.
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

#ifndef REAK_PROX_SPHERE_CCYLINDER_HPP
#define REAK_PROX_SPHERE_CCYLINDER_HPP

#include "proximity_finder_3D.hpp"

#include "shapes/sphere.hpp"
#include "shapes/capped_cylinder.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/**
 * This class is for proximity queries between a sphere and a cylinder.
 */
class prox_sphere_ccylinder : public proximity_finder_3D {
  protected:
    
    shared_ptr< sphere > mSphere;
    shared_ptr< capped_cylinder > mCCylinder;
    
  public:
    
    /** Returns the first shape involved in the proximity query. */
    virtual shared_ptr< shape_3D > getShape1() const;
    /** Returns the second shape involved in the proximity query. */
    virtual shared_ptr< shape_3D > getShape2() const;
    
    /** This function performs the proximity query on its associated shapes. */
    virtual void computeProximity();
    
    /** 
     * Default constructor. 
     * \param aSphere The sphere involved in the proximity query.
     * \param aCCylinder The capped cylinder involved in the proximity query.
     */
    prox_sphere_ccylinder(const shared_ptr< sphere >& aSphere = shared_ptr< sphere >(),
                          const shared_ptr< capped_cylinder >& aCCylinder = shared_ptr< capped_cylinder >());
    
    /** Destructor. */
    virtual ~prox_sphere_ccylinder() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const;
    
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int);
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(prox_sphere_ccylinder,0xC3200011,1,"prox_sphere_ccylinder",proximity_finder_3D)
    
};


};

};

#endif











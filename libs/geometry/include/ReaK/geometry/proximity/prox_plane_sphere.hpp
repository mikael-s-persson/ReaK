/**
 * \file prox_plane_sphere.hpp
 *
 * This library declares a class for proximity queries between a plane and a sphere.
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

#ifndef REAK_PROX_PLANE_SPHERE_HPP
#define REAK_PROX_PLANE_SPHERE_HPP

#include "proximity_finder_3D.hpp"

#include <ReaK/geometry/shapes/sphere.hpp>
#include <ReaK/geometry/shapes/plane.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


proximity_record_3D compute_proximity(const plane& aPlane, 
                                      const shape_3D_precompute_pack& aPack1,
                                      const sphere& aSphere, 
                                      const shape_3D_precompute_pack& aPack2);

proximity_record_3D compute_proximity(const sphere& aSphere, 
                                      const shape_3D_precompute_pack& aPack1,
                                      const plane& aPlane, 
                                      const shape_3D_precompute_pack& aPack2);

/**
 * This class is for proximity queries between a plane and a sphere.
 */
class prox_plane_sphere : public proximity_finder_3D {
  protected:
    
    const plane* mPlane;
    const sphere* mSphere;
    
  public:
    
    /** This function performs the proximity query on its associated shapes. */
    virtual void computeProximity(const shape_3D_precompute_pack& aPack1, 
                                  const shape_3D_precompute_pack& aPack2);
    
    /** 
     * Default constructor. 
     * \param aPlane The plane involved in the proximity query.
     * \param aSphere The sphere involved in the proximity query.
     */
    prox_plane_sphere(const plane* aPlane = NULL,
                      const sphere* aSphere = NULL);
    
    /** Destructor. */
    virtual ~prox_plane_sphere() { };
    
};


};

};

#endif


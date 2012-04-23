/**
 * \file ellipsoid.hpp
 *
 * This library declares a class to represent ellipsoids (3D ellipse).
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

#ifndef REAK_ELLIPSOID_HPP
#define REAK_ELLIPSOID_HPP

#include "shape_3D.hpp"

#include "color.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/** This class represents an ellipsoids (3D ellipse). */
class ellipsoid : public shape_3D {
  protected:
    
    vect<double,3> mRadii;
    
    RGBA_color mColor;
    
  public:
    
    /**
     * This function returns the radii of the ellipsoid (x, y and z radius, or principle axes).
     * \return The radii of the ellipsoid (x, y and z radius, or principle axes).
     */
    const vect<double,3>& getRadii() const { return mRadii; };
    /**
     * This function sets the radii of the ellipsoid (x, y and z radius, or principle axes).
     * \param aRadii The new radii of the ellipsoid (x, y and z radius, or principle axes).
     */
    void setRadii(const vect<double,3>& aRadii) { mRadii = aRadii; };
    
    /** 
     * This function returns the color of the ellipsoid.
     * \return The color.
     */
    const RGBA_color& getColor() const { return mColor; };
    /** 
     * This function sets the color of the ellipsoid.
     * \param aColor The new color.
     */
    void setColor(const RGBA_color& aColor) { mColor = aColor; };
    
    virtual void render() const;
    
    /**
     * Default constructor.
     * \param aName The name of the object.
     * \param aAnchor The anchor object for the geometry.
     * \param aPose The pose of the geometry (relative to the anchor).
     * \param aRadii The radii of the ellipsoid (x, y and z radius, or principle axes).
     */
    ellipsoid(const std::string& aName = "",
              const shared_ptr< pose_3D<double> >& aAnchor = shared_ptr< pose_3D<double> >(),
	      const pose_3D<double>& aPose = pose_3D<double>(),
	      const vect<double,3>& aRadii = vect<double,3>(1.0,1.0,1.0),
	      const RGBA_color& aColor = RGBA_color(255,255,255,255));
    
    /**
     * Default destructor.
     */
    virtual ~ellipsoid() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const;

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int);

    RK_RTTI_MAKE_CONCRETE_1BASE(ellipsoid,0xC3100011,1,"ellipsoid",shape_3D)

};


};

};

#endif











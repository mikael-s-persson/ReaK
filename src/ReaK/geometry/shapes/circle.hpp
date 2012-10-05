/**
 * \file circle.hpp
 *
 * This library declares a class to represent a circle.
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

#ifndef REAK_CIRCLE_HPP
#define REAK_CIRCLE_HPP

#include "shape_2D.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/** This class represents a circle. */
class circle : public shape_2D {
  protected:
    
    double mRadius;
    
  public:
    
    /**
     * This function returns the maximum radius of the shape (radius of the circle that bounds the shape).
     * \return The maximum radius of the shape.
     */
    virtual double getBoundingRadius() const;
    
    /**
     * This function returns the radius of the circle.
     * \return The radius of the circle.
     */
    double getRadius() const { return mRadius; };
    /**
     * This function sets the radius of the circle.
     * \param aRadius The new radius for the circle.
     */
    void setRadius(double aRadius) { mRadius = aRadius; };
    
    /**
     * Default constructor.
     * \param aName The name of the object.
     * \param aAnchor The anchor object for the geometry.
     * \param aPose The pose of the geometry (relative to the anchor).
     * \param aRadius The radius of the circle.
     */
    circle(const std::string& aName = "",
           const shared_ptr< pose_2D<double> >& aAnchor = shared_ptr< pose_2D<double> >(),
	   const pose_2D<double>& aPose = pose_2D<double>(),
	   double aRadius = 1.0);
    
    /**
     * Default destructor.
     */
    virtual ~circle() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const;

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int);

    RK_RTTI_MAKE_ABSTRACT_1BASE(circle,0xC310000C,1,"circle",shape_2D)

};


};

};

#endif











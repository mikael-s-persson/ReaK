/**
 * \file capped_rectangle.hpp
 *
 * This library declares a class to represent a capped rectangle in 2D.
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

#ifndef REAK_CAPPED_RECTANGLE_HPP
#define REAK_CAPPED_RECTANGLE_HPP

#include "shape_2D.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/** This class represents a capped rectangle in 2D (aligned about its center pose, with x-axis ends capped with circles). */
class capped_rectangle : public shape_2D {
  protected:
    
    vect<double,2> mDimensions;
    
  public:
    
    /** 
     * This function returns the dimensions of the capped rectangle.
     * \return The dimensions of the capped rectangle.
     */
    const vect<double,2>& getDimensions() const { return mDimensions; };
    /** 
     * This function sets the new dimensions of the capped rectangle.
     * \param aDimensions The new dimensions of the capped rectangle.
     */
    void setDimensions(const vect<double,2>& aDimensions) { mDimensions = aDimensions; };
    
    virtual void render() const;
    
    /**
     * Default constructor.
     * \param aName The name of the object.
     * \param aAnchor The anchor object for the geometry.
     * \param aPose The pose of the geometry (relative to the anchor).
     * \param aDimensions The dimensions.
     */
    capped_rectangle(const std::string& aName = "",
	             const shared_ptr< pose_2D<double> >& aAnchor = shared_ptr< pose_2D<double> >(),
	             const pose_2D<double>& aPose = pose_2D<double>(),
	             const vect<double,2>& aDimensions = vect<double,2>());
    
    /**
     * Default destructor.
     */
    virtual ~capped_rectangle() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const;

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int);

    RK_RTTI_MAKE_CONCRETE_1BASE(capped_rectangle,0xC310000D,1,"capped_rectangle",shape_2D)
    
};


};

};

#endif











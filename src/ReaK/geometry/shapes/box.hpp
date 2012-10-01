/**
 * \file box.hpp
 *
 * This library declares a class to represent boxes.
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

#ifndef REAK_BOX_HPP
#define REAK_BOX_HPP

#include "shape_3D.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/** This class represents a box aligned along its center pose. */
class box : public shape_3D {
  protected:
    
    vect<double,3> mDimensions;
    
  public:
    
    /**
     * This function returns the dimensions of the box.
     * \return The dimensions of the box.
     */
    const vect<double,3>& getDimensions() const { return mDimensions; };
    /**
     * This function sets the dimensions of the box.
     * \param aDimensions The new dimensions of the box.
     */
    void setDimensions(const vect<double,3>& aDimensions) { mDimensions = aDimensions; };
    
    /**
     * Default constructor.
     * \param aName The name of the object.
     * \param aAnchor The anchor object for the geometry.
     * \param aPose The pose of the geometry (relative to the anchor).
     * \param aDimensions The dimensions of the box.
     */
    box(const std::string& aName = "",
        const shared_ptr< pose_3D<double> >& aAnchor = shared_ptr< pose_3D<double> >(),
	const pose_3D<double>& aPose = pose_3D<double>(),
	const vect<double,3>& aDimensions = vect<double,3>(1.0,1.0,1.0));
    
    /**
     * Default destructor.
     */
    virtual ~box() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const;

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int);

    RK_RTTI_MAKE_CONCRETE_1BASE(box,0xC3100013,1,"box",shape_3D)

};


};

};

#endif











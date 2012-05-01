/**
 * \file grid_3D.hpp
 *
 * This library declares a class for a 3D grid to be rendered.
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

#ifndef REAK_GRID_3D_HPP
#define REAK_GRID_3D_HPP

#include "geometry_3D.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/** This class is a class for a 3D grid to be rendered. */
class grid_3D : public geometry_3D {
  protected:
    
    vect<double,3> mDimensions;
    vect<std::size_t,3> mSquareCounts;
    
  public:
    
    /** 
     * This function returns the dimensions of the grid.
     * \return The dimensions.
     */
    const vect<double,3>& getDimensions() const { return mDimensions; };
    /** 
     * This function sets the new dimensions of the grid.
     * \param aDimensions The new dimensions.
     */
    void setDimensions(const vect<double,3>& aDimensions) { mDimensions = aDimensions; };
    
    /** 
     * This function returns the square-counts of the grid.
     * \return The square-counts.
     */
    const vect<std::size_t,3>& getSquareCounts() const { return mSquareCounts; };
    /** 
     * This function sets the square-counts of the grid.
     * \param mSquareCounts The new square-counts.
     */
    void setSquareCounts(const vect<std::size_t,3>& aSquareCounts) { mSquareCounts = aSquareCounts; };
    
    virtual void render() const;
    
    /**
     * Default constructor.
     * \param aName The name of the object.
     * \param aAnchor The anchor object for the geometry.
     * \param aPose The pose of the geometry (relative to the anchor).
     * \param aDimensions The dimensions.
     * \param aSquareCounts The square-counts.
     */
    grid_3D(const std::string& aName = "",
	    const shared_ptr< pose_3D<double> >& aAnchor = shared_ptr< pose_3D<double> >(),
	    const pose_3D<double>& aPose = pose_3D<double>(),
	    const vect<double,3>& aDimensions = vect<double,3>(1.0,1.0,1.0),
	    const vect<std::size_t,3>& aSquareCounts = vect<std::size_t,3>(10,10,10));
    
    /**
     * Default destructor.
     */
    virtual ~grid_3D() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const;

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int);

    RK_RTTI_MAKE_CONCRETE_1BASE(grid_3D,0xC3100007,1,"grid_3D",geometry_3D)

};


};

};

#endif











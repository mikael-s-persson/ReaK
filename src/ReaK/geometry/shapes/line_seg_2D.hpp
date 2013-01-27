/**
 * \file line_seg_2D.hpp
 *
 * This library declares a 2D line-segment class to render a line in a 2D scene.
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

#ifndef REAK_LINE_SEG_2D_HPP
#define REAK_LINE_SEG_2D_HPP

#include "geometry_2D.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/** This class defines a 2D line-segment class to render a line in a 2D scene. */
class line_seg_2D : public geometry_2D {
  protected:
    vect<double,2> mStart;
    vect<double,2> mEnd;
    
  public:
    
    /** 
     * This function returns the start point of the line-segment.
     * \return The start point.
     */
    const vect<double,2>& getStart() const { return mStart; };
    /** 
     * This function sets the new start point of the line-segment.
     * \param aStart The new start point.
     */
    void setStart(const vect<double,2>& aStart) { mStart = aStart; };
    
    /** 
     * This function returns the end point of the line-segment.
     * \return The end point.
     */
    const vect<double,2>& getEnd() const { return mEnd; };
    /** 
     * This function sets the end point of the line-segment.
     * \param aEnd The new end point.
     */
    void setEnd(const vect<double,2>& aEnd) { mEnd = aEnd; };
    
    /**
     * Default constructor.
     * \param aName The name of the object.
     * \param aAnchor The anchor object for the geometry.
     * \param aPose The pose of the geometry (relative to the anchor).
     * \param aStart The start point of the line-segment.
     * \param aEnd The end point of the line-segment.
     */
    line_seg_2D(const std::string& aName = "",
                const shared_ptr< pose_2D<double> >& aAnchor = shared_ptr< pose_2D<double> >(),
		const pose_2D<double>& aPose = pose_2D<double>(),
		const vect<double,2>& aStart = (vect<double,2>()),
		const vect<double,2>& aEnd = (vect<double,2>()));
    
    /**
     * Default destructor.
     */
    virtual ~line_seg_2D() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const;

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int);

    RK_RTTI_MAKE_CONCRETE_1BASE(line_seg_2D,0xC3100004,1,"line_seg_2D",geometry_2D)

};


};

};

#endif











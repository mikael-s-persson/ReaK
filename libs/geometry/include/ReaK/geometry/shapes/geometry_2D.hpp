/**
 * \file geometry_2D.hpp
 *
 * This library declares the base-class for all 2D geometric objects (renderable).
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

#ifndef REAK_GEOMETRY_2D_HPP
#define REAK_GEOMETRY_2D_HPP

#include <ReaK/core/base/named_object.hpp>

#include <ReaK/math/kinetostatics/pose_2D.hpp>


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/** This class is the base-class for all 2D geometric objects (renderable). */
class geometry_2D : public named_object {
  protected:
    shared_ptr< pose_2D<double> > mAnchor;
    pose_2D<double> mPose;
    
  public:
    
    /** 
     * This function returns the anchor of the geometry.
     * \return A shared-pointer to the anchor pose object.
     */
    const shared_ptr< pose_2D<double> >& getAnchor() const { return mAnchor; };
    /** 
     * This function sets the anchor of the geometry.
     * \param aAnchor A shared-pointer to the new anchor pose object.
     */
    void setAnchor(const shared_ptr< pose_2D<double> >& aAnchor);
    
    /** 
     * This function returns the pose of the geometry.
     * \return The pose object.
     */
    const pose_2D<double>& getPose() const { return mPose; };
    /** 
     * This function sets the pose of the geometry.
     * \param aPose The new pose object.
     */
    void setPose(const pose_2D<double>& aPose);
    
    /**
     * Default constructor.
     * \param aName The name of the object.
     * \param aAnchor The anchor object for the geometry.
     * \param aPose The pose of the geometry (relative to the anchor).
     */
    geometry_2D(const std::string& aName = "",
                const shared_ptr< pose_2D<double> >& aAnchor = shared_ptr< pose_2D<double> >(),
                const pose_2D<double>& aPose = pose_2D<double>());
    
    /**
     * Default destructor.
     */
    virtual ~geometry_2D() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const;

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int);

    RK_RTTI_MAKE_CONCRETE_1BASE(geometry_2D,0xC3100000,1,"geometry_2D",named_object)

};


};

};

#endif











/**
 * \file proximity_record_2D.hpp
 *
 * This library declares a class to record the results of a proximity query between 2D shapes.
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

#ifndef REAK_PROXIMITY_RECORD_2D_HPP
#define REAK_PROXIMITY_RECORD_2D_HPP

#include <ReaK/core/base/shared_object.hpp>

#include <ReaK/math/lin_alg/vect_alg.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {

/**
 * This class stores the data which results from a proximity query with 2D shapes.
 */
class proximity_record_2D : public shared_object {
  public:
    
    /** Holds the closest point (in global coordinates) on the first shape involved in the proximity query. */
    vect<double,2> mPoint1;
    /** Holds the closest point (in global coordinates) on the second shape involved in the proximity query. */
    vect<double,2> mPoint2;
    
    /** Holds the distance between the shapes, a negative value denotes penetration. */
    double mDistance;
    
    /** Default constructor. */
    proximity_record_2D() { };
    
    /** Destructor. */
    virtual ~proximity_record_2D() { };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const;

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int);

    RK_RTTI_MAKE_CONCRETE_1BASE(proximity_record_2D,0xC3200001,1,"proximity_record_2D",shared_object)

};



};

};

#endif











/**
 * \file plane_point_mindist.hpp
 *
 * This library implements geometric models used to compute the kinematics of an end-frame which
 * should follow the minimum distance from the base-frame to a given (fixed) plance in 3D space.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date May 2010
 */

/*
 *    Copyright 2011 Sven Mikael Persson
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

#ifndef REAK_PLANE_POINT_MINDIST_HPP
#define REAK_PLANE_POINT_MINDIST_HPP

#include "kte_map.hpp"

#include "kinetostatics/kinetostatics.hpp"

namespace ReaK {

namespace kte {


/**
 * This class implements the geometric calculation of a 3D end-frame which follows a (fixed) plane while
 * maintain the minimum distance to a base-frame.
 */
class plane_point_mindist_3D : public kte_map {
  private:
    shared_pointer< frame_3D<double> >::type mBase; ///< Holds the base-frame, or kinematic input, or free-point.
    shared_pointer< frame_3D<double> >::type mEnd; ///< Holds the end-frame, or kinematic output, or min-dist point to base-frame, on the plane.
    vect<double,3> mNormal; ///< Holds the normal vector of the plane, in global coordinates.
    double mOrigin; ///< Holds the distance to the origin, i.e., mOrigin * mNormal is the vector from the plane to a plane parallel and intersecting the origin.

  public:

    /** Get read-write access to mNormal. */
    vect<double,3>& Normal() { return mNormal; };
    /** Get read-only access to mNormal. */
    const vect<double,3>& Normal() const { return mNormal; };

    /** Get read-write access to mOrigin. */
    double& Origin() { return mOrigin; };
    /** Get read-only access to mOrigin. */
    double Origin() const { return mOrigin; };

    /**
     * Default constructor.
     */
    plane_point_mindist_3D(const std::string& aName = "") : kte_map(aName), mBase(), mEnd(), mNormal(), mOrigin(0.0) { };

    /**
     * Parametrized constructor.
     * \param aName the name of the KTE model.
     * \param aBase the base-frame, or kinematic input, or free-point.
     * \param aEnd the end-frame, or kinematic output, or min-dist point to base-frame, on the plane.
     * \param aNormal the normal vector of the plane, in global coordinates.
     * \param aOrigin the distance to the origin, i.e., aOrigin * aNormal is the vector from the plane to a plane parallel and intersecting the origin.
     */
    plane_point_mindist_3D(const std::string& aName,
                           const shared_pointer< frame_3D<double> >::type& aBase,
                           const shared_pointer< frame_3D<double> >::type& aEnd,
                           const vect<double,3>& aNormal,
                           double aOrigin) :
                           kte_map(aName),
                           mBase(aBase),
                           mEnd(aEnd),
                           mNormal(aNormal),
                           mOrigin(aOrigin) { };

    /**
     * Default destructor.
     */
    virtual ~plane_point_mindist_3D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mBase)
        & RK_SERIAL_SAVE_WITH_NAME(mEnd)
        & RK_SERIAL_SAVE_WITH_NAME(mNormal)
        & RK_SERIAL_SAVE_WITH_NAME(mOrigin);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mBase)
        & RK_SERIAL_LOAD_WITH_NAME(mEnd)
        & RK_SERIAL_LOAD_WITH_NAME(mNormal)
        & RK_SERIAL_LOAD_WITH_NAME(mOrigin);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(plane_point_mindist_3D,0xC2100030,1,"plane_point_mindist_3D",kte_map)


};


};

};


#endif





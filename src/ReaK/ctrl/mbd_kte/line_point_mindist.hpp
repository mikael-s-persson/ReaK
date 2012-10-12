/**
 * \file line_point_mindist.hpp
 *
 * This library implements geometric models used to compute the kinematics of an end-frame which
 * should follow the minimum distance from the base-frame to a given (fixed) line in 2D or 3D space.
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

#ifndef REAK_LINE_POINT_MINDIST_HPP
#define REAK_LINE_POINT_MINDIST_HPP

#include "kte_map.hpp"

#include "kinetostatics/kinetostatics.hpp"

namespace ReaK {

namespace kte {




/**
 * This class implements the geometric calculation of a 2D end-frame which follows a (fixed) line while
 * maintain the minimum distance to a base-frame.
 */
class line_point_mindist_2D : public kte_map {
  private:
    shared_ptr< frame_2D<double> > mBase; ///< Holds the base-frame, or kinematic input, or free-point.
    shared_ptr< frame_2D<double> > mEnd; ///< Holds the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
    vect<double,2> mTangent; ///< Holds the tangent vector of the (fixed) line in global coordinates.
    vect<double,2> mOriginMinDist; ///< Holds the minimum-distance vector from origin to the line, in global coordinates.

  public:
    
    /**
     * Returns a reference to the base-frame, or kinematic input, or free-point.
     * \return A reference to the base-frame, or kinematic input, or free-point.
     */
    shared_ptr< frame_2D<double> >& Base() { return mBase; };
    /**
     * Returns a const-reference to the base-frame, or kinematic input, or free-point.
     * \return A const-reference to the base-frame, or kinematic input, or free-point.
     */
    const shared_ptr< frame_2D<double> >& Base() const { return mBase; };
    
    /**
     * Returns a reference to the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
     * \return A reference to the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
     */
    shared_ptr< frame_2D<double> >& End() { return mEnd; };
    /**
     * Returns a const-reference to the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
     * \return A const-reference to the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
     */
    const shared_ptr< frame_2D<double> >& End() const { return mEnd; };
    
    /** Get read-write access to mTangent. */
    vect<double,2>& Tangent() { return mTangent; };
    /** Get read-only access to mTangent. */
    const vect<double,2>& Tangent() const { return mTangent; };

    /** Get read-write access to mOriginMinDist. */
    vect<double,2>& OriginMinDist() { return mOriginMinDist; };
    /** Get read-only access to mOriginMinDist. */
    const vect<double,2>& OriginMinDist() const { return mOriginMinDist; };

    /**
     * Default constructor.
     */
    line_point_mindist_2D(const std::string& aName = "") : kte_map(aName), mBase(), mEnd(), mTangent(), mOriginMinDist() { };

    /**
     * Parametrized constructor.
     * \param aName the name for the KTE model.
     * \param aBase the base-frame, or kinematic input, or free-point.
     * \param aEnd the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
     * \param aTangent the tangent vector of the (fixed) line in global coordinates.
     * \param aOriginMinDist the minimum-distance vector from origin to the line, in global coordinates.
     */
    line_point_mindist_2D(const std::string& aName,
                          const shared_ptr< frame_2D<double> >& aBase,
                          const shared_ptr< frame_2D<double> >& aEnd,
                          const vect<double,2>& aTangent,
                          const vect<double,2>& aOriginMinDist) :
                          kte_map(aName),
                          mBase(aBase),
                          mEnd(aEnd),
                          mTangent(aTangent),
                          mOriginMinDist(aOriginMinDist) { };

    /**
     * Default destructor.
     */
    virtual ~line_point_mindist_2D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mBase)
        & RK_SERIAL_SAVE_WITH_NAME(mEnd)
        & RK_SERIAL_SAVE_WITH_NAME(mTangent)
        & RK_SERIAL_SAVE_WITH_NAME(mOriginMinDist);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mBase)
        & RK_SERIAL_LOAD_WITH_NAME(mEnd)
        & RK_SERIAL_LOAD_WITH_NAME(mTangent)
        & RK_SERIAL_LOAD_WITH_NAME(mOriginMinDist);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(line_point_mindist_2D,0xC2100031,1,"line_point_mindist_2D",kte_map)


};

/**
 * This class implements the geometric calculation of a 3D end-frame which follows a (fixed) line while
 * maintain the minimum distance to a base-frame.
 */
class line_point_mindist_3D : public kte_map {
  private:
    shared_ptr< frame_3D<double> > mBase; ///< Holds the base-frame, or kinematic input, or free-point.
    shared_ptr< frame_3D<double> > mEnd; ///< Holds the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
    vect<double,3> mTangent; ///< Holds the tangent vector of the (fixed) line in global coordinates.
    vect<double,3> mOriginMinDist; ///< Holds the minimum-distance vector from origin to the line, in global coordinates.

  public:
    
    /**
     * Returns a reference to the base-frame, or kinematic input, or free-point.
     * \return A reference to the base-frame, or kinematic input, or free-point.
     */
    shared_ptr< frame_3D<double> >& Base() { return mBase; };
    /**
     * Returns a const-reference to the base-frame, or kinematic input, or free-point.
     * \return A const-reference to the base-frame, or kinematic input, or free-point.
     */
    const shared_ptr< frame_3D<double> >& Base() const { return mBase; };
    
    /**
     * Returns a reference to the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
     * \return A reference to the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
     */
    shared_ptr< frame_3D<double> >& End() { return mEnd; };
    /**
     * Returns a const-reference to the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
     * \return A const-reference to the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
     */
    const shared_ptr< frame_3D<double> >& End() const { return mEnd; };
    
    /** Get read-write access to mTangent. */
    vect<double,3>& Tangent() { return mTangent; };
    /** Get read-only access to mTangent. */
    const vect<double,3>& Tangent() const { return mTangent; };

    /** Get read-write access to mOriginMinDist. */
    vect<double,3>& OriginMinDist() { return mOriginMinDist; };
    /** Get read-only access to mOriginMinDist. */
    const vect<double,3>& OriginMinDist() const { return mOriginMinDist; };

    /**
     * Default constructor.
     */
    line_point_mindist_3D(const std::string& aName = "") : kte_map(aName), mBase(), mEnd(), mTangent(), mOriginMinDist() { };

    /**
     * Parametrized constructor.
     * \param aName the name for the KTE model.
     * \param aBase the base-frame, or kinematic input, or free-point.
     * \param aEnd the end-frame, or kinematic output, or min-dist point to base-frame, on the line.
     * \param aTangent the tangent vector of the (fixed) line in global coordinates.
     * \param aOriginMinDist the minimum-distance vector from origin to the line, in global coordinates.
     */
    line_point_mindist_3D(const std::string& aName,
                          const shared_ptr< frame_3D<double> >& aBase,
                          const shared_ptr< frame_3D<double> >& aEnd,
                          const vect<double,3>& aTangent,
                          const vect<double,3>& aOriginMinDist) :
                          kte_map(aName),
                          mBase(aBase),
                          mEnd(aEnd),
                          mTangent(aTangent),
                          mOriginMinDist(aOriginMinDist) { };

    /**
     * Default destructor.
     */
    virtual ~line_point_mindist_3D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mBase)
        & RK_SERIAL_SAVE_WITH_NAME(mEnd)
        & RK_SERIAL_SAVE_WITH_NAME(mTangent)
        & RK_SERIAL_SAVE_WITH_NAME(mOriginMinDist);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mBase)
        & RK_SERIAL_LOAD_WITH_NAME(mEnd)
        & RK_SERIAL_LOAD_WITH_NAME(mTangent)
        & RK_SERIAL_LOAD_WITH_NAME(mOriginMinDist);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(line_point_mindist_3D,0xC2100032,1,"line_point_mindist_3D",kte_map)


};


};


};

#endif





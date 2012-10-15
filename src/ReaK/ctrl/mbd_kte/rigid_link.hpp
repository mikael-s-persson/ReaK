/**
 * \file rigid_link.hpp
 *
 * This library declares the KTE models for rigid-links in 1D (generalized coordinate), 2D and 3D space.
 * These models implement a fixed coordinate transformation from a base frame to an end frame as well as
 * the rigid transmission of the force in the reverse direction.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2010
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

#ifndef REAK_RIGID_LINK_HPP
#define REAK_RIGID_LINK_HPP

#include "kte_map.hpp"

#include "kinetostatics/kinetostatics.hpp"

namespace ReaK {

namespace kte {


/**
 * This class implements the KTE model for a rigid-link on a generalized coordinate (linear or angular).
 * It maps the base frame to an end frame via a constant offset value.
 */
class rigid_link_gen : public kte_map {
  private:
    shared_ptr< gen_coord<double> > mBase; ///< Holds the base frame of the rigid-link.
    shared_ptr< gen_coord<double> > mEnd; ///< Holds the end frame of the rigid-link.
    double mOffset; ///< Holds the offset of the rigid-link (or length of the link).

  public:
    
    /**
     * Sets the link's base frame.
     * \param aPtr The new link's base frame.
     */
    void setBaseFrame(const shared_ptr< gen_coord<double> >& aPtr) { mBase = aPtr; };
    /**
     * Returns the link's base frame.
     * \return The link's base frame.
     */
    shared_ptr< gen_coord<double> > BaseFrame() const { return mBase; };
    
    /**
     * Sets the link's output frame.
     * \param aPtr The new link's output frame.
     */
    void setEndFrame(const shared_ptr< gen_coord<double> >& aPtr) { mEnd = aPtr; };
    /**
     * Returns the link's output frame.
     * \return The link's output frame.
     */
    shared_ptr< gen_coord<double> > EndFrame() const { return mEnd; };
    
    /**
     * Sets the link's offset.
     * \param aValue The link's new offset.
     */
    void setOffset(double aValue) { mOffset = aValue; };
    /**
     * Returns the link's offset.
     * \return The link's offset.
     */
    double Offset() const { return mOffset; };
    
    /**
     * Default constructor.
     */
    rigid_link_gen(const std::string& aName = "") : kte_map(aName), mBase(), mEnd(), mOffset(0.0) { };

    /**
     * Parametrized constructor.
     * \param aName the name of the KTE model.
     * \param aBase the base frame of the rigid-link.
     * \param aEnd the end frame of the rigid-link.
     * \param aOffset the offset of the rigid-link (or length of the link).
     */
    rigid_link_gen(const std::string& aName,
		   const shared_ptr< gen_coord<double> >& aBase,
		   const shared_ptr< gen_coord<double> >& aEnd,
		   double aOffset) :
		   kte_map(aName),
		   mBase(aBase),
		   mEnd(aEnd),
		   mOffset(aOffset) { };

    /**
     * Default destructor.
     */
    virtual ~rigid_link_gen() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mBase)
        & RK_SERIAL_SAVE_WITH_NAME(mEnd)
	& RK_SERIAL_SAVE_WITH_NAME(mOffset);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mBase)
        & RK_SERIAL_LOAD_WITH_NAME(mEnd)
	& RK_SERIAL_LOAD_WITH_NAME(mOffset);

    };

    RK_RTTI_MAKE_CONCRETE_1BASE(rigid_link_gen,0xC2100007,1,"rigid_link_gen",kte_map)


};

/**
 * This class implements the KTE model for a rigid-link on a 2D space.
 * It maps the base frame to an end frame via a constant offset value (translation and rotation).
 */
class rigid_link_2D : public kte_map {
  private:
    shared_ptr< frame_2D<double> > mBase; ///< Holds the base frame of the rigid-link.
    shared_ptr< frame_2D<double> > mEnd; ///< Holds the end frame of the rigid-link.
    pose_2D<double> mPoseOffset; ///< Holds the pose offset of the rigid-link (length and twist of the link).

  public:
    
    /**
     * Sets the link's base frame.
     * \param aPtr The new link's base frame.
     */
    void setBaseFrame(const shared_ptr< frame_2D<double> >& aPtr) { mBase = aPtr; };
    /**
     * Returns the link's base frame.
     * \return The link's base frame.
     */
    shared_ptr< frame_2D<double> > BaseFrame() const { return mBase; };
    
    /**
     * Sets the link's output frame.
     * \param aPtr The new link's output frame.
     */
    void setEndFrame(const shared_ptr< frame_2D<double> >& aPtr) { mEnd = aPtr; };
    /**
     * Returns the link's output frame.
     * \return The link's output frame.
     */
    shared_ptr< frame_2D<double> > EndFrame() const { return mEnd; };
    
    /**
     * Sets the link's pose offset (position vector and rotation).
     * \param aValue The link's new pose offset (position vector and rotation).
     */
    void setPoseOffset(const pose_2D<double>& aValue) { mPoseOffset = aValue; };
    /**
     * Returns the link's pose offset (position vector and rotation).
     * \return The link's pose offset (position vector and rotation).
     */
    pose_2D<double> PoseOffset() const { return mPoseOffset; };
    
    /**
     * Default constructor.
     */
    rigid_link_2D(const std::string& aName = "") : kte_map(aName), mBase(), mEnd(), mPoseOffset() { };

    /**
     * Parametrized constructor.
     * \param aName the name of the KTE model.
     * \param aBase the base frame of the rigid-link.
     * \param aEnd the end frame of the rigid-link.
     * \param aPoseOffset the pose offset of the rigid-link (length and twist of the link).
     */
    rigid_link_2D(const std::string& aName,
		  const shared_ptr< frame_2D<double> >& aBase,
		  const shared_ptr< frame_2D<double> >& aEnd,
		  const pose_2D<double>& aPoseOffset) :
		  kte_map(aName),
		  mBase(aBase),
		  mEnd(aEnd),
		  mPoseOffset(aPoseOffset) { };

    /**
     * Default destructor.
     */
    virtual ~rigid_link_2D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mBase)
        & RK_SERIAL_SAVE_WITH_NAME(mEnd)
	& RK_SERIAL_SAVE_WITH_NAME(mPoseOffset);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mBase)
        & RK_SERIAL_LOAD_WITH_NAME(mEnd)
	& RK_SERIAL_LOAD_WITH_NAME(mPoseOffset);

    };

    RK_RTTI_MAKE_CONCRETE_1BASE(rigid_link_2D,0xC2100008,1,"rigid_link_2D",kte_map)

};

/**
 * This class implements the KTE model for a rigid-link on a 3D space.
 * It maps the base frame to an end frame via a constant offset value (translation and rotation).
 */
class rigid_link_3D : public kte_map {
  private:
    shared_ptr< frame_3D<double> > mBase; ///< Holds the base frame of the rigid-link.
    shared_ptr< frame_3D<double> > mEnd; ///< Holds the end frame of the rigid-link.
    pose_3D<double> mPoseOffset; ///< Holds the pose offset of the rigid-link (length and twist of the link).

  public:
    
    /**
     * Sets the link's base frame.
     * \param aPtr The new link's base frame.
     */
    void setBaseFrame(const shared_ptr< frame_3D<double> >& aPtr) { mBase = aPtr; };
    /**
     * Returns the link's base frame.
     * \return The link's base frame.
     */
    shared_ptr< frame_3D<double> > BaseFrame() const { return mBase; };
    
    /**
     * Sets the link's output frame.
     * \param aPtr The new link's output frame.
     */
    void setEndFrame(const shared_ptr< frame_3D<double> >& aPtr) { mEnd = aPtr; };
    /**
     * Returns the link's output frame.
     * \return The link's output frame.
     */
    shared_ptr< frame_3D<double> > EndFrame() const { return mEnd; };
    
    /**
     * Sets the link's pose offset (position vector and rotation).
     * \param aValue The link's new pose offset (position vector and rotation).
     */
    void setPoseOffset(const pose_3D<double>& aValue) { mPoseOffset = aValue; };
    /**
     * Returns the link's pose offset (position vector and rotation).
     * \return The link's pose offset (position vector and rotation).
     */
    pose_3D<double> PoseOffset() const { return mPoseOffset; };
    
    /**
     * Default constructor.
     */
    rigid_link_3D(const std::string& aName = "") : kte_map(aName), mBase(), mEnd(), mPoseOffset() { };

    /**
     * Parametrized constructor.
     * \param aName the name of the KTE model.
     * \param aBase the base frame of the rigid-link.
     * \param aEnd the end frame of the rigid-link.
     * \param aPoseOffset the pose offset of the rigid-link (length and twist of the link).
     */
    rigid_link_3D(const std::string& aName,
		  const shared_ptr< frame_3D<double> >& aBase,
		  const shared_ptr< frame_3D<double> >& aEnd,
		  const pose_3D<double>& aPoseOffset) :
		  kte_map(aName),
		  mBase(aBase),
		  mEnd(aEnd),
		  mPoseOffset(aPoseOffset) { };

    /**
     * Default destructor.
     */
    virtual ~rigid_link_3D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mBase)
        & RK_SERIAL_SAVE_WITH_NAME(mEnd)
        & RK_SERIAL_SAVE_WITH_NAME(mPoseOffset);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mBase)
        & RK_SERIAL_LOAD_WITH_NAME(mEnd)
        & RK_SERIAL_LOAD_WITH_NAME(mPoseOffset);

    };

    RK_RTTI_MAKE_CONCRETE_1BASE(rigid_link_3D,0xC2100009,1,"rigid_link_3D",kte_map)


};

};

};


#endif





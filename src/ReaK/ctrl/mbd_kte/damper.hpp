/**
 * \file damper.hpp
 *
 * This library declares the KTE models for linear dampers, often reference to as dashpots.
 * Here damper classes are available for 2D and 3D point-to-point dashpots as well as a
 * damper between two generalized coordinates. The model of the damper is a basic linear, constant
 * damping coefficient that multiplies the velocity to obtain the force that will oppose the relative
 * velocity between the anchors.
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

#ifndef DAMPER_HPP
#define DAMPER_HPP

#include "kte_map.hpp"

#include "math/kinetostatics.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.KTE */
namespace kte {


/** This class defines a damper acting between two generalized coordinates. */
class damper_gen : public kte_map {
  private:
    boost::shared_ptr< gen_coord<double> > mAnchor1; ///< Holds the first generalized coordinate.
    boost::shared_ptr< gen_coord<double> > mAnchor2; ///< Holds the second generalized coordinate.
    double mDamping; ///< The damping coefficient (in Ns/m or Nms/rad).

  public:

    /** Get read-write access to mDamping. */
    double& Damping() { return mDamping; };
    /** Get read-only access to mDamping. */
    double Damping() const { return mDamping; };

    /**
     * Default constructor.
     */
    damper_gen(const std::string& aName = "") : kte_map(aName), mAnchor1(), mAnchor2(), mDamping(0.0) { };

    /**
     * Parametrized constructor.
     * \param aName name of the KTE model.
     * \param aAnchor1 first attach point of the damper.
     * \param aAnchor2 second attach point of the damper.
     * \param aDamping damping coefficient (in Ns/m for a linear generalized coord. or Nms/rad for an angular generalized coord.).
     */
    damper_gen(const std::string& aName,
	       boost::shared_ptr< gen_coord<double> > aAnchor1,
	       boost::shared_ptr< gen_coord<double> > aAnchor2,
	       double aDamping) :
	       kte_map(aName),
	       mAnchor1(aAnchor1),
	       mAnchor2(aAnchor2),
	       mDamping(aDamping) { };

    /**
     * Default destructor.
     */
    virtual ~damper_gen() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchor1)
        & RK_SERIAL_SAVE_WITH_NAME(mAnchor2)
	& RK_SERIAL_SAVE_WITH_NAME(mDamping);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchor1)
        & RK_SERIAL_LOAD_WITH_NAME(mAnchor2)
	& RK_SERIAL_LOAD_WITH_NAME(mDamping);

    };

    RK_RTTI_MAKE_CONCRETE_1BASE(damper_gen,0xC2100010,1,"damper_gen",kte_map)

};

/** This class defines a damper acting between two 2D frames. */
class damper_2D : public kte_map {
  private:
    boost::shared_ptr< frame_2D<double> > mAnchor1; ///< Holds the first 2D frame.
    boost::shared_ptr< frame_2D<double> > mAnchor2; ///< Holds the second 2D frame.
    double mDamping; ///< The damping coefficient (in Ns/m).

  public:

    /** Get read-write access to mDamping. */
    double& Damping() { return mDamping; };
    /** Get read-only access to mDamping. */
    double Damping() const { return mDamping; };

    /**
     * Default constructor.
     */
    damper_2D(const std::string& aName = "") : kte_map(aName), mAnchor1(), mAnchor2(), mDamping(0.0) { };

    /**
     * Parametrized constructor.
     * \param aName name of the KTE model.
     * \param aAnchor1 first attach point of the damper.
     * \param aAnchor2 second attach point of the damper.
     * \param aDamping damping coefficient (in Ns/m).
     */
    damper_2D(const std::string& aName,
	      boost::shared_ptr< frame_2D<double> > aAnchor1,
	      boost::shared_ptr< frame_2D<double> > aAnchor2,
	      double aDamping) :
	      kte_map(aName),
	      mAnchor1(aAnchor1),
	      mAnchor2(aAnchor2),
	      mDamping(aDamping) { };

    /**
     * Default destructor.
     */
    virtual ~damper_2D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchor1)
        & RK_SERIAL_SAVE_WITH_NAME(mAnchor2)
	& RK_SERIAL_SAVE_WITH_NAME(mDamping);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchor1)
        & RK_SERIAL_LOAD_WITH_NAME(mAnchor2)
	& RK_SERIAL_LOAD_WITH_NAME(mDamping);

    };

    RK_RTTI_MAKE_CONCRETE_1BASE(damper_2D,0xC2100011,1,"damper_2D",kte_map)

};

/** This class defines a damper acting between two 3D frames. */
class damper_3D : public kte_map {
  private:
    boost::shared_ptr< frame_3D<double> > mAnchor1; ///< Holds the first 3D frame.
    boost::shared_ptr< frame_3D<double> > mAnchor2; ///< Holds the second 3D frame.
    double mDamping; ///< The damping coefficient (in Ns/m).

  public:

    /** Get read-write access to mDamping. */
    double& Damping() { return mDamping; };
    /** Get read-only access to mDamping. */
    double Damping() const { return mDamping; };

    /**
     * Default constructor.
     */
    damper_3D(const std::string& aName = "") : kte_map(aName), mAnchor1(), mAnchor2(), mDamping(0.0) { };

    /**
     * Parametrized constructor.
     * \param aName name of the KTE model.
     * \param aAnchor1 first attach point of the damper.
     * \param aAnchor2 second attach point of the damper.
     * \param aDamping damping coefficient (in Ns/m).
     */
    damper_3D(const std::string& aName,
	      boost::shared_ptr< frame_3D<double> > aAnchor1,
	      boost::shared_ptr< frame_3D<double> > aAnchor2,
	      double aDamping) :
	      kte_map(aName),
	      mAnchor1(aAnchor1),
	      mAnchor2(aAnchor2),
	      mDamping(aDamping) { };

    /**
     * Default destructor.
     */
    virtual ~damper_3D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchor1)
        & RK_SERIAL_SAVE_WITH_NAME(mAnchor2)
	& RK_SERIAL_SAVE_WITH_NAME(mDamping);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchor1)
        & RK_SERIAL_LOAD_WITH_NAME(mAnchor2)
	& RK_SERIAL_LOAD_WITH_NAME(mDamping);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(damper_3D,0xC2100012,1,"damper_3D",kte_map)

};


};

};

#endif











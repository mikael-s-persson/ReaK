/**
 * \file inertia.hpp
 *
 * This library declares KTE models for the inertial forces acting on generalized coordinates,
 * 2D and 3D kinetostatic frames.
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

#ifndef REAK_INERTIA_HPP
#define REAK_INERTIA_HPP

#include "base/defs.hpp"

#include "kte_map.hpp"

#include "kinetostatics/kinetostatics.hpp"
#include "lin_alg/mat_alg.hpp"

#include "jacobian_joint_map.hpp"

namespace ReaK {

namespace kte {


/**
 * This class implements inertia on a generalized coordinate.
 */
class inertia_gen : public kte_map {
  private:
    shared_pointer< joint_dependent_gen_coord >::type mCenterOfMass; ///< Holds the center-of-mass generalized coordinate for the inertial element.
    double mMass; ///< Holds the mass of the inertial element (in kg or kgm2).

  public:

    /** Get read-write access to mMass. */
    double& Mass() { return mMass; };
    /** Get read-only access to mMass. */
    double Mass() const { return mMass; };
    
    /** Get read-write access to mCenterOfMass. */
    shared_pointer< joint_dependent_gen_coord >::type& CenterOfMass() { return mCenterOfMass; };
    /** Get read-only access to mCenterOfMass. */
    shared_pointer< joint_dependent_gen_coord >::type CenterOfMass() const { return mCenterOfMass; };

    /**
     * Default constructor.
     */
    inertia_gen(const std::string& aName = "") : kte_map(aName), mCenterOfMass(), mMass(0.0) { };

    /**
     * Parametrized constructor.
     * \param aName name of the KTE model.
     * \param aCenterOfMass the center-of-mass coordinate for the inertial element.
     * \param aMass the mass of the inertial element (in kg or kgm2).
     * \param aUpStreamJoints the jacobian mappings of up-stream joints.
     * \param aUpStream2DJoints the jacobian mappings of up-stream joints.
     * \param aUpStream3DJoints the jacobian mappings of up-stream joints.
     */
    inertia_gen(const std::string& aName,
		const shared_pointer< joint_dependent_gen_coord >::type& aCenterOfMass,
		double aMass) :
		kte_map(aName),
		mCenterOfMass(aCenterOfMass),
		mMass(aMass) { };

    /**
     * Default destructor.
     */
    virtual ~inertia_gen() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type());

    virtual void clearForce();

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mCenterOfMass)
	& RK_SERIAL_SAVE_WITH_NAME(mMass);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mCenterOfMass)
	& RK_SERIAL_LOAD_WITH_NAME(mMass);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(inertia_gen,0xC210000A,1,"inertia_gen",kte_map)


};

/**
 * This class implements inertia on a 2D frame.
 */
class inertia_2D : public kte_map {
  private:
    shared_pointer< joint_dependent_frame_2D >::type mCenterOfMass; ///< Holds the center-of-mass 2D frame for the inertial element.
    double mMass; ///< Holds the mass of the inertial element (in kg).
    double mMomentOfInertia; ///< Holds the moment of inertia of the inertial element (in kgm2).

  public:

    /** Get read-write access to mMass. */
    double& Mass() { return mMass; };
    /** Get read-only access to mMass. */
    double Mass() const { return mMass; };

    /** Get read-write access to mMomentOfInertia. */
    double& MomentOfInertia() { return mMomentOfInertia; };
    /** Get read-only access to mMomentOfInertia. */
    double MomentOfInertia() const { return mMomentOfInertia; };
    
    /** Get read-write access to mCenterOfMass. */
    shared_pointer< joint_dependent_frame_2D >::type& CenterOfMass() { return mCenterOfMass; };
    /** Get read-only access to mCenterOfMass. */
    shared_pointer< joint_dependent_frame_2D >::type CenterOfMass() const { return mCenterOfMass; };

    /**
     * Default constructor.
     */
    inertia_2D(const std::string& aName = "") : kte_map(aName), mCenterOfMass(), mMass(0.0), mMomentOfInertia(0.0) { };

    /**
     * Parametrized constructor.
     * \param aName name of the KTE model.
     * \param aCenterOfMass the center-of-mass 2D frame for the inertial element.
     * \param aMass the mass of the inertial element (in kg).
     * \param aMomentOfInertia the moment-of-inertia of the inertial element (in kgm2).
     * \param aUpStreamJoints the jacobian mappings of up-stream joints.
     * \param aUpStream2DJoints the jacobian mappings of up-stream joints.
     * \param aUpStream3DJoints the jacobian mappings of up-stream joints.
     */
    inertia_2D(const std::string& aName,
	       const shared_pointer< joint_dependent_frame_2D >::type& aCenterOfMass,
	       double aMass,
	       double aMomentOfInertia) :
	       kte_map(aName),
	       mCenterOfMass(aCenterOfMass),
	       mMass(aMass),
	       mMomentOfInertia(aMomentOfInertia) { };

    /**
     * Default destructor.
     */
    virtual ~inertia_2D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mCenterOfMass)
        & RK_SERIAL_SAVE_WITH_NAME(mMass)
	& RK_SERIAL_SAVE_WITH_NAME(mMomentOfInertia);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mCenterOfMass)
        & RK_SERIAL_LOAD_WITH_NAME(mMass)
	& RK_SERIAL_LOAD_WITH_NAME(mMomentOfInertia);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(inertia_2D,0xC210000B,1,"inertia_2D",kte_map)


};

/**
 * This class implements inertia on a 3D frame.
 */
class inertia_3D : public kte_map {
  private:
    shared_pointer< joint_dependent_frame_3D >::type mCenterOfMass; ///< Holds the center-of-mass 3D frame for the inertial element.
    double mMass; ///< Holds the mass of the inertial element (in kg).
    mat<double,mat_structure::symmetric> mInertiaTensor; ///< Holds the inertia tensor of the inertial element (in kgm2).

  public:

    /** Get read-write access to mMass. */
    double& Mass() { return mMass; };
    /** Get read-only access to mMass. */
    double Mass() const { return mMass; };

    /** Get read-write access to mInertiaTensor. */
    mat<double,mat_structure::symmetric>& InertiaTensor() { return mInertiaTensor; };
    /** Get read-only access to mInertiaTensor. */
    const mat<double,mat_structure::symmetric>& InertiaTensor() const { return mInertiaTensor; };
    
    /** Get read-write access to mCenterOfMass. */
    shared_pointer< joint_dependent_frame_3D >::type& CenterOfMass() { return mCenterOfMass; };
    /** Get read-only access to mCenterOfMass. */
    shared_pointer< joint_dependent_frame_3D >::type CenterOfMass() const { return mCenterOfMass; };

    /**
     * Default constructor.
     */
    inertia_3D(const std::string& aName = "") : kte_map(aName), mCenterOfMass(), mMass(0.0), mInertiaTensor(3) { };

    /**
     * Parametrized constructor.
     * \param aName name of the KTE model.
     * \param aCenterOfMass the center-of-mass 3D frame for the inertial element.
     * \param aMass the mass of the inertial element (in kg).
     * \param aInertiaTensor the inertia tensor of the inertial element (in kgm2).
     * \param aUpStreamJoints the jacobian mappings of up-stream joints.
     */
    inertia_3D(const std::string& aName,
	       const shared_pointer< joint_dependent_frame_3D >::type& aCenterOfMass,
	       double aMass,
	       const mat<double,mat_structure::symmetric>& aInertiaTensor) :
               kte_map(aName),
	       mCenterOfMass(aCenterOfMass),
	       mMass(aMass),
	       mInertiaTensor(aInertiaTensor) { };

    /**
     * Default destructor.
     */
    virtual ~inertia_3D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mCenterOfMass)
        & RK_SERIAL_SAVE_WITH_NAME(mMass)
        & RK_SERIAL_SAVE_WITH_NAME(mInertiaTensor);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mCenterOfMass)
        & RK_SERIAL_LOAD_WITH_NAME(mMass)
        & RK_SERIAL_LOAD_WITH_NAME(mInertiaTensor);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(inertia_3D,0xC210000C,1,"inertia_3D",kte_map)

};


};

};


#endif










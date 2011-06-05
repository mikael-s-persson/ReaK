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

#ifndef INERTIA_HPP
#define INERTIA_HPP

#include "kte_map.hpp"

#include "math/kinetostatics.hpp"
#include "math/mat_alg.hpp"

#include "jacobian_joint_map.hpp"

namespace ReaK {

namespace kte {


/**
 * This class implements inertia on a generalized coordinate.
 */
class inertia_gen : public kte_map {
  private:
    boost::shared_ptr< gen_coord<double> > mCenterOfMass; ///< Holds the center-of-mass generalized coordinate for the inertial element.
    double mMass; ///< Holds the mass of the inertial element (in kg or kgm2).

    jacobian_joint_map_gen mUpStreamJoints; ///< Holds the jacobian mappings of up-stream joints.
    jacobian_joint2D_map_gen mUpStream2DJoints; ///< Holds the jacobian mappings of up-stream joints.
    jacobian_joint3D_map_gen mUpStream3DJoints; ///< Holds the jacobian mappings of up-stream joints.

    friend class mass_matrix_calc;

  public:

    /** Get read-write access to mMass. */
    double& Mass() { return mMass; };
    /** Get read-only access to mMass. */
    double Mass() const { return mMass; };

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
		boost::shared_ptr< ReaK::gen_coord<double> > aCenterOfMass,
		double aMass,
                const jacobian_joint_map_gen& aUpStreamJoints = jacobian_joint_map_gen(),
                const jacobian_joint2D_map_gen& aUpStream2DJoints = jacobian_joint2D_map_gen(),
                const jacobian_joint3D_map_gen& aUpStream3DJoints = jacobian_joint3D_map_gen()) :
		kte_map(aName),
		mCenterOfMass(aCenterOfMass),
		mMass(aMass),
		mUpStreamJoints(aUpStreamJoints),
		mUpStream2DJoints(aUpStream2DJoints),
		mUpStream3DJoints(aUpStream3DJoints) { };

    /**
     * Default destructor.
     */
    virtual ~inertia_gen() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mCenterOfMass)
	& RK_SERIAL_SAVE_WITH_NAME(mMass)
        & RK_SERIAL_SAVE_WITH_NAME(mUpStreamJoints)
        & RK_SERIAL_SAVE_WITH_NAME(mUpStream2DJoints)
        & RK_SERIAL_SAVE_WITH_NAME(mUpStream3DJoints);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mCenterOfMass)
	& RK_SERIAL_LOAD_WITH_NAME(mMass)
	& RK_SERIAL_LOAD_WITH_NAME(mUpStreamJoints)
	& RK_SERIAL_LOAD_WITH_NAME(mUpStream2DJoints)
	& RK_SERIAL_LOAD_WITH_NAME(mUpStream3DJoints);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(inertia_gen,0xC210000A,1,"inertia_gen",kte_map)


};

/**
 * This class implements inertia on a 2D frame.
 */
class inertia_2D : public kte_map {
  private:
    boost::shared_ptr< frame_2D<double> > mCenterOfMass; ///< Holds the center-of-mass 2D frame for the inertial element.
    double mMass; ///< Holds the mass of the inertial element (in kg).
    double mMomentOfInertia; ///< Holds the moment of inertia of the inertial element (in kgm2).

    jacobian_joint_map_2D mUpStreamJoints; ///< Holds the jacobian mappings of up-stream joints.
    jacobian_joint2D_map_2D mUpStream2DJoints; ///< Holds the jacobian mappings of up-stream joints.
    jacobian_joint3D_map_2D mUpStream3DJoints; ///< Holds the jacobian mappings of up-stream joints.

    friend class mass_matrix_calc;

  public:

    /** Get read-write access to mMass. */
    double& Mass() { return mMass; };
    /** Get read-only access to mMass. */
    double Mass() const { return mMass; };

    /** Get read-write access to mMomentOfInertia. */
    double& MomentOfInertia() { return mMomentOfInertia; };
    /** Get read-only access to mMomentOfInertia. */
    double MomentOfInertia() const { return mMomentOfInertia; };

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
	       boost::shared_ptr< frame_2D<double> > aCenterOfMass,
	       double aMass,
	       double aMomentOfInertia,
               const jacobian_joint_map_2D& aUpStreamJoints = jacobian_joint_map_2D(),
               const jacobian_joint2D_map_2D& aUpStream2DJoints = jacobian_joint2D_map_2D(),
               const jacobian_joint3D_map_2D& aUpStream3DJoints = jacobian_joint3D_map_2D()) :
	       kte_map(aName),
	       mCenterOfMass(aCenterOfMass),
	       mMass(aMass),
	       mMomentOfInertia(aMomentOfInertia),
	       mUpStreamJoints(aUpStreamJoints),
	       mUpStream2DJoints(aUpStream2DJoints),
	       mUpStream3DJoints(aUpStream3DJoints) { };

    /**
     * Default destructor.
     */
    virtual ~inertia_2D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mCenterOfMass)
        & RK_SERIAL_SAVE_WITH_NAME(mMass)
	& RK_SERIAL_SAVE_WITH_NAME(mMomentOfInertia)
	& RK_SERIAL_SAVE_WITH_NAME(mUpStreamJoints)
        & RK_SERIAL_SAVE_WITH_NAME(mUpStream2DJoints)
        & RK_SERIAL_SAVE_WITH_NAME(mUpStream3DJoints);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mCenterOfMass)
        & RK_SERIAL_LOAD_WITH_NAME(mMass)
	& RK_SERIAL_LOAD_WITH_NAME(mMomentOfInertia)
	& RK_SERIAL_LOAD_WITH_NAME(mUpStreamJoints)
	& RK_SERIAL_LOAD_WITH_NAME(mUpStream2DJoints)
	& RK_SERIAL_LOAD_WITH_NAME(mUpStream3DJoints);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(inertia_2D,0xC210000B,1,"inertia_2D",kte_map)


};

/**
 * This class implements inertia on a 3D frame.
 */
class inertia_3D : public kte_map {
  private:
    boost::shared_ptr< frame_3D<double> > mCenterOfMass; ///< Holds the center-of-mass 3D frame for the inertial element.
    double mMass; ///< Holds the mass of the inertial element (in kg).
    mat<double,mat_structure::symmetric> mInertiaTensor; ///< Holds the inertia tensor of the inertial element (in kgm2).

    jacobian_joint_map_3D mUpStreamJoints; ///< Holds the jacobian mappings of up-stream joints.
    jacobian_joint2D_map_3D mUpStream2DJoints; ///< Holds the jacobian mappings of up-stream joints.
    jacobian_joint3D_map_3D mUpStream3DJoints; ///< Holds the jacobian mappings of up-stream joints.

    friend class mass_matrix_calc;

  public:

    /** Get read-write access to mMass. */
    double& Mass() { return mMass; };
    /** Get read-only access to mMass. */
    double Mass() const { return mMass; };

    /** Get read-write access to mInertiaTensor. */
    mat<double,mat_structure::symmetric>& InertiaTensor() { return mInertiaTensor; };
    /** Get read-only access to mInertiaTensor. */
    const mat<double,mat_structure::symmetric>& InertiaTensor() const { return mInertiaTensor; };

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
	       boost::shared_ptr< frame_3D<double> > aCenterOfMass,
	       double aMass,
	       const mat<double,mat_structure::symmetric>& aInertiaTensor,
               const jacobian_joint_map_3D& aUpStreamJoints = jacobian_joint_map_3D(),
               const jacobian_joint2D_map_3D& aUpStream2DJoints = jacobian_joint2D_map_3D(),
               const jacobian_joint3D_map_3D& aUpStream3DJoints = jacobian_joint3D_map_3D()) :
               kte_map(aName),
	       mCenterOfMass(aCenterOfMass),
	       mMass(aMass),
	       mInertiaTensor(aInertiaTensor),
	       mUpStreamJoints(aUpStreamJoints),
	       mUpStream2DJoints(aUpStream2DJoints),
	       mUpStream3DJoints(aUpStream3DJoints) { };

    /**
     * Default destructor.
     */
    virtual ~inertia_3D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mCenterOfMass)
        & RK_SERIAL_SAVE_WITH_NAME(mMass)
        & RK_SERIAL_SAVE_WITH_NAME(mInertiaTensor)
        & RK_SERIAL_SAVE_WITH_NAME(mUpStreamJoints)
        & RK_SERIAL_SAVE_WITH_NAME(mUpStream2DJoints)
        & RK_SERIAL_SAVE_WITH_NAME(mUpStream3DJoints);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mCenterOfMass)
        & RK_SERIAL_LOAD_WITH_NAME(mMass)
        & RK_SERIAL_LOAD_WITH_NAME(mInertiaTensor)
        & RK_SERIAL_LOAD_WITH_NAME(mUpStreamJoints)
        & RK_SERIAL_LOAD_WITH_NAME(mUpStream2DJoints)
        & RK_SERIAL_LOAD_WITH_NAME(mUpStream3DJoints);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(inertia_3D,0xC210000C,1,"inertia_3D",kte_map)

};


};

};


#endif










/**
 * \file driving_actuator.hpp
 *
 * This library declares KTE models for a driving actuator, that is, an actuator that applies some
 * given force, from a controller for example. This model is needed to not only compute the net force
 * acting on the joint coordinate, that will lead to its motion, but also the reaction force that
 * the actuator will imbue its base frame.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date April 2010
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

#ifndef REAK_DRIVING_ACTUATOR_HPP
#define REAK_DRIVING_ACTUATOR_HPP

#include "force_actuator.hpp"
#include "kte_system_input.hpp"

namespace ReaK {

namespace kte {


/**
 * This class implements the model for an actuator that is driving a generalized coordinate, as for a
 * powered revolute or prismatic joint.
 */
class driving_actuator_gen : public force_actuator_gen, public system_input {
  protected:
    double mDriveForce; ///< Holds the current force applied to the generalized coordinate.

  public:
    /** 
     * Sets the drive force currently applied by this actuator.
     * \param aValue The new drive force applied by this actuator.
     */
    void setDriveForce(double aValue) { mDriveForce = aValue; };
    /**
     * Returns the drive force currently applied by this actuator.
     * \return The drive force currently applied by this actuator.
     */
    double DriveForce() const { return mDriveForce; };
    
    virtual unsigned int getInputCount() const { return 1; };
    
    virtual void setInput(unsigned int i, double aDriveForce) { mDriveForce = aDriveForce; };
    virtual double getInput(unsigned int i) const { return mDriveForce; }; 
    
    /**
     * Default constructor.
     */
    driving_actuator_gen(const std::string& aName = "") : force_actuator_gen(aName), system_input(aName), mDriveForce(0.0) { };

    /**
     * Parametrized constructor.
     * \param aName the name of the KTE model.
     * \param aFrame the generalized coordinate on which the actuator acts.
     * \param aJoint the joint which will react to the actuator's force.
     */
    driving_actuator_gen(const std::string& aName,
                         const shared_ptr< gen_coord<double> >& aFrame,
                         const shared_ptr< reacting_kte_gen >& aJoint) :
                         force_actuator_gen(aName,aFrame,aJoint),
			 system_input(aName),
                         mDriveForce(0.0) { };

    /**
     * Default destructor.
     */
    virtual ~driving_actuator_gen() { };

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      force_actuator_gen::save(A,force_actuator_gen::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mDriveForce);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      force_actuator_gen::load(A,force_actuator_gen::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mDriveForce);
    };

    RK_RTTI_MAKE_CONCRETE_2BASE(driving_actuator_gen,0xC2100023,1,"driving_actuator_gen",force_actuator_gen,system_input)

};

/**
 * This class implements the model for an actuator that is driving a 2D frame, as for a
 * powered planar joint.
 */
class driving_actuator_2D : public force_actuator_2D, public system_input {
  protected:
    vect<double,2> mDriveForce; ///< Holds the force vector applied to the 2D frame.
    double mDriveTorque; ///< Holds the torque applied to the 2D frame.

  public:
    /** 
     * Sets the drive force currently applied by this actuator.
     * \param aValue The new drive force applied by this actuator.
     */
    void setDriveForce(const vect<double,2>& aValue) { mDriveForce = aValue; };
    /**
     * Returns the drive force currently applied by this actuator.
     * \return The drive force currently applied by this actuator.
     */
    vect<double,2> DriveForce() const { return mDriveForce; };
    
    /** 
     * Sets the drive torque currently applied by this actuator.
     * \param aValue The new drive torque applied by this actuator.
     */
    void setDriveTorque(double aValue) { mDriveTorque = aValue; };
    /**
     * Returns the drive torque currently applied by this actuator.
     * \return The drive torque currently applied by this actuator.
     */
    double DriveTorque() const { return mDriveTorque; };
    
    virtual unsigned int getInputCount() const { return 3; };
    
    virtual void setInput(unsigned int i, double aValue) { 
      if(i < 2) 
	mDriveForce[i] = aValue;
      else
	mDriveTorque = aValue;
    };
    virtual double getInput(unsigned int i) const { 
      if(i < 2)
        return mDriveForce[i]; 
      else
	return mDriveTorque;
    }; 
    
    /**
     * Default constructor.
     */
    driving_actuator_2D(const std::string& aName = "") : force_actuator_2D(aName), mDriveForce(0.0,0.0), mDriveTorque(0.0) { };

    /**
     * Parametrized constructor.
     * \param aName the name of the KTE model.
     * \param aFrame the 2D frame on which the actuator acts.
     * \param aJoint the joint which will react to the actuator's force and torque.
     */
    driving_actuator_2D(const std::string& aName,
                         const shared_ptr< frame_2D<double> >& aFrame,
                         const shared_ptr< reacting_kte_2D >& aJoint) :
                         force_actuator_2D(aName,aFrame,aJoint),
                         mDriveForce(0.0,0.0), mDriveTorque(0.0) { };

    /**
     * Default destructor.
     */
    virtual ~driving_actuator_2D() { };

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      force_actuator_2D::save(A,force_actuator_2D::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mDriveForce)
        & RK_SERIAL_SAVE_WITH_NAME(mDriveTorque);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      force_actuator_2D::load(A,force_actuator_2D::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mDriveForce)
        & RK_SERIAL_LOAD_WITH_NAME(mDriveTorque);
    };

    RK_RTTI_MAKE_CONCRETE_2BASE(driving_actuator_2D,0xC2100024,1,"driving_actuator_2D",force_actuator_2D,system_input)

};

/**
 * This class implements the model for an actuator that is driving a 3D frame, as for a
 * powered planar joint.
 */
class driving_actuator_3D : public force_actuator_3D, public system_input {
  protected:
    vect<double,3> mDriveForce; ///< Holds the force vector applied to the 3D frame.
    vect<double,3> mDriveTorque; ///< Holds the torque vector applied to the 3D frame.

  public:
    /** 
     * Sets the drive force currently applied by this actuator.
     * \param aValue The new drive force applied by this actuator.
     */
    void setDriveForce(const vect<double,3>& aValue) { mDriveForce = aValue; };
    /**
     * Returns the drive force currently applied by this actuator.
     * \return The drive force currently applied by this actuator.
     */
    vect<double,3> DriveForce() const { return mDriveForce; };
    
    /** 
     * Sets the drive torque currently applied by this actuator.
     * \param aValue The new drive torque applied by this actuator.
     */
    void setDriveTorque(const vect<double,3>& aValue) { mDriveTorque = aValue; };
    /**
     * Returns the drive torque currently applied by this actuator.
     * \return The drive torque currently applied by this actuator.
     */
    vect<double,3> DriveTorque() const { return mDriveTorque; };
    
    virtual unsigned int getInputCount() const { return 6; };
    
    virtual void setInput(unsigned int i, double aValue) { 
      if(i < 3) 
	mDriveForce[i] = aValue;
      else if(i < 6)
	mDriveTorque[i-3] = aValue;
      else
	mDriveForce[0] = aValue;
    };
    virtual double getInput(unsigned int i) const { 
      if(i < 3) 
	return mDriveForce[i];
      else if(i < 6)
	return mDriveTorque[i-3];
      else
	return mDriveForce[0];
    }; 
    
    /**
     * Default constructor.
     */
    driving_actuator_3D(const std::string& aName = "") : force_actuator_3D(aName),
                                                         mDriveForce(0.0,0.0,0.0),
                                                         mDriveTorque(0.0,0.0,0.0) { };

    /**
     * Parametrized constructor.
     * \param aName the name of the KTE model.
     * \param aFrame the 3D frame on which the actuator acts.
     * \param aJoint the joint which will react to the actuator's force and torque.
     */
    driving_actuator_3D(const std::string& aName,
                         const shared_ptr< frame_3D<double> >& aFrame,
                         const shared_ptr< reacting_kte_3D >& aJoint) :
                         force_actuator_3D(aName,aFrame,aJoint),
                         mDriveForce(0.0,0.0,0.0),
                         mDriveTorque(0.0,0.0,0.0) { };

    /**
     * Default destructor.
     */
    virtual ~driving_actuator_3D() { };

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      force_actuator_3D::save(A,force_actuator_3D::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mDriveForce)
        & RK_SERIAL_SAVE_WITH_NAME(mDriveTorque);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      force_actuator_3D::load(A,force_actuator_3D::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mDriveForce)
        & RK_SERIAL_LOAD_WITH_NAME(mDriveTorque);
    };

    RK_RTTI_MAKE_CONCRETE_2BASE(driving_actuator_3D,0xC2100025,1,"driving_actuator_3D",force_actuator_3D,system_input)

};


};

};


#endif //DRIVING_ACTUATOR_HPP









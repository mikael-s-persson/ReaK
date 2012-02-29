/**
 * \file joint_friction.hpp
 * 
 * This library does not contain any working class. An attempt was made at modeling joints
 * with friction, but this attempt failed. However, dry_revolute_joint.hpp contains a working 
 * version for a revolute joint.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
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

#ifndef REAK_JOINT_FRICTION_HPP
#define REAK_JOINT_FRICTION_HPP

#include "kte_map.hpp"
#include "kinetostatics/kinetostatics.hpp"

namespace ReaK {

namespace kte {


/**
 * NOT A VALID MODEL.
 */
class joint_dry_microslip_gen : public kte_map {
  private:
    shared_ptr<gen_coord<double> > mAnchor;
    double mStictionVelocity;
    double mSlipVelocity;
    double mStictionCoef;
    double mSlipCoef;

  public:
    double& StictionVelocity() { return mStictionVelocity; };
    double StictionVelocity() const { return mStictionVelocity; };

    double& SlipVelocity() { return mSlipVelocity; };
    double SlipVelocity() const { return mSlipVelocity; };

    double& StictionCoefficient() { return mStictionCoef; };
    double StictionCoefficient() const { return mStictionCoef; };

    double& SlipCoefficient() { return mSlipCoef; };
    double SlipCoefficient() const { return mSlipCoef; };

    /**
     * Default constructor.
     */
    joint_dry_microslip_gen(const std::string& aName = "") : kte_map(aName),
                                                             mAnchor(),
                                                             mStictionVelocity(1E-6),
                                                             mSlipVelocity(2E-6),
                                                             mStictionCoef(0.0),
                                                             mSlipCoef(0.0){ };

    /**
     * Parametrized constructor.
     */
    joint_dry_microslip_gen(const std::string& aName,
                            const shared_ptr< gen_coord<double> >& aAnchor,
                            double aStictionVelocity,
                            double aSlipVelocity,
                            double aStictionCoef,
                            double aSlipCoef) :
                            kte_map(aName),
                            mAnchor(aAnchor),
                            mStictionVelocity(aStictionVelocity),
                            mSlipVelocity(aSlipVelocity),
                            mStictionCoef(aStictionCoef),
                            mSlipCoef(aSlipCoef){ };

    /**
     * Default destructor.
     */
    virtual ~joint_dry_microslip_gen() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void clearForce();


    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchor)
        & RK_SERIAL_SAVE_WITH_NAME(mStictionVelocity)
        & RK_SERIAL_SAVE_WITH_NAME(mSlipVelocity)
        & RK_SERIAL_SAVE_WITH_NAME(mStictionCoef)
        & RK_SERIAL_SAVE_WITH_NAME(mSlipCoef);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchor)
        & RK_SERIAL_LOAD_WITH_NAME(mStictionVelocity)
        & RK_SERIAL_LOAD_WITH_NAME(mSlipVelocity)
        & RK_SERIAL_LOAD_WITH_NAME(mStictionCoef)
        & RK_SERIAL_LOAD_WITH_NAME(mSlipCoef);

    };

    RK_RTTI_MAKE_CONCRETE_1BASE(joint_dry_microslip_gen,0xC2100020,1,"joint_dry_microslip_gen",kte_map)

};


/**
 * NOT A VALID MODEL.
 */
class joint_viscosity_gen : public kte_map {
  private:
    shared_ptr<gen_coord<double> > mAnchor;
    double mViscosity;

  public:
    double& Viscosity() { return mViscosity; };
    double Viscosity() const { return mViscosity; };

    /**
     * Default constructor.
     */
    joint_viscosity_gen(const std::string& aName = "") : kte_map(aName),
                                                         mAnchor(),
                                                         mViscosity(0.0) { };

    /**
     * Parametrized constructor.
     */
    joint_viscosity_gen(const std::string& aName,
                        const shared_ptr< gen_coord<double> >& aAnchor,
                        double aViscosity) :
                        kte_map(aName),
                        mAnchor(aAnchor),
                        mViscosity(aViscosity) { };

    /**
     * Default destructor.
     */
    virtual ~joint_viscosity_gen() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchor)
        & RK_SERIAL_SAVE_WITH_NAME(mViscosity);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchor)
        & RK_SERIAL_LOAD_WITH_NAME(mViscosity);

    };

    RK_RTTI_MAKE_CONCRETE_1BASE(joint_viscosity_gen,0xC2100021,1,"joint_viscosity_gen",kte_map)

};


};

};


#endif








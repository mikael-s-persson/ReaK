/**
 * \file vmc_revolute_joint.hpp
 *
 * This library declares KTE models for revolute joints which also incorporate dry-friction modelling
 * in preparation for Virtual Model Control, adding the friction to be compensated.
 * The friction model used here is the micro-slip model of Coulomb friction. The dry-friction
 * model is incorporated into the joint model because of the need to know the normal force on the joint
 * to obtain the friction force. The micro-slip model is describe, mathematically as:
 * \f[\tau_f = \{ \begin{array}{cr}
                -\mu_{stiction} \|\vec{N}\| \frac{\dot{q}}{v_{slip}} & |\dot{q}| \leq v_{slip} \\
                -sign(\dot{q}) \mu_{slip} \|\vec{N}\| & |\dot{q}| > v_{slip}
                \end{array}\f]
 * Where \f$\mu_{stiction}\f$ and \f$\mu_{stiction}\f$ are the stiction and slip friction coefficients
 * respectively, which are, of course, specific to the joint properties, including material, gear-ratio and stages,
 * joint radius, etc. Also here, \f$\|\vec{N}\|\f$ is the norm of the normal force vector, \f$\dot{q}\f$ is the
 * joint's velocity, and finally, \f$v_{slip}\f$ is the slip velocity (the tolerated micro-slip).
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

#ifndef REAK_VMC_REVOLUTE_JOINT_HPP
#define REAK_VMC_REVOLUTE_JOINT_HPP

#include "revolute_joint.hpp"

namespace ReaK {

namespace kte {


/**
 * This class implements the dry-frictional revolute joint model for 2D. Dry friction model uses the
 * micro-slip model of Coulomb friction and the friction is prepared for compensation by the VMC.
 */
class vmc_revolute_joint_2D : public revolute_joint_2D {
  protected:
    double mStictionCoef; ///< Holds the stiction coefficient in Dry Coulomb friction.
    double mSlipCoef; ///< Holds the slip coefficient in Dry Coulomb friction.
    double mSlipVelocity; ///< Holds the micro-slip tolerance in joint speed.

  public:

    /** Get read-write access to mStictionCoef. */
    double& StictionCoefficient() { return mStictionCoef; };
    /** Get read-only access to mStictionCoef. */
    double StictionCoefficient() const { return mStictionCoef; };

    /** Get read-write access to mSlipCoef. */
    double& SlipCoefficient() { return mSlipCoef; };
    /** Get read-only access to mSlipCoef. */
    double SlipCoefficient() const { return mSlipCoef; };

    /** Get read-write access to mSlipVelocity. */
    double& SlipVelocity() { return mSlipVelocity; };
    /** Get read-only access to mSlipVelocity. */
    double SlipVelocity() const { return mSlipVelocity; };

     /**
     * Default constructor.
     */
    vmc_revolute_joint_2D(const std::string& aName = "") : revolute_joint_2D(aName), mStictionCoef(0.5), mSlipCoef(0.3), mSlipVelocity(1E-5) { };

    /**
     * Parametrized constructor.
     * \param aName the name of this KTE model.
     * \param aAngle the generalized coordinate associated with the displacement of this joint.
     * \param aBase the coordinate frame at the base of the joint.
     * \param aEnd the coordinate frame just after the joint transformations are applied.
     * \param aJacobian a pointer to contain the Jacobian frame produced by this joint, default value will disable the Jacobian frame's calculation.
     * \param aStictionCoef the stiction coefficient in Dry Coulomb friction.
     * \param aSlipCoef the slip coefficient in Dry Coulomb friction.
     * \param aSlipVelocity the micro-slip tolerance in joint speed.
     */
    vmc_revolute_joint_2D(const std::string& aName,
                          const shared_pointer< gen_coord<double> >::type& aAngle,
                          const shared_pointer< frame_2D<double> >::type& aBase,
                          const shared_pointer< frame_2D<double> >::type& aEnd,
                          const shared_pointer< jacobian_gen_2D<double> >::type& aJacobian = shared_pointer< jacobian_gen_2D<double> >::type(),
                          double aStictionCoef = 0.5,
			  double aSlipCoef = 0.3,
			  double aSlipVelocity = 1E-5) :
			  revolute_joint_2D(aName,aAngle,aBase,aEnd,aJacobian),
			  mStictionCoef(aStictionCoef),
			  mSlipCoef(aSlipCoef),
			  mSlipVelocity(aSlipVelocity) { };

    /**
     * Default destructor.
     */
    virtual ~vmc_revolute_joint_2D() { };

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type());

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      revolute_joint_2D::save(A,revolute_joint_2D::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mStictionCoef)
        & RK_SERIAL_SAVE_WITH_NAME(mSlipCoef)
        & RK_SERIAL_SAVE_WITH_NAME(mSlipVelocity);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      revolute_joint_2D::load(A,revolute_joint_2D::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mStictionCoef)
        & RK_SERIAL_LOAD_WITH_NAME(mSlipCoef)
        & RK_SERIAL_LOAD_WITH_NAME(mSlipVelocity);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(vmc_revolute_joint_2D,0xC2100026,1,"vmc_revolute_joint_2D",revolute_joint_2D)

};

/**
 * This class implements the dry-frictional revolute joint model for 3D. Dry friction model uses the
 * micro-slip model of Coulomb friction and the friction is prepared for compensation by the VMC.
 */
class vmc_revolute_joint_3D : public revolute_joint_3D {
  protected:
    double mStictionCoef; ///< Holds the stiction coefficient in Dry Coulomb friction.
    double mSlipCoef; ///< Holds the slip coefficient in Dry Coulomb friction.
    double mSlipVelocity; ///< Holds the micro-slip tolerance in joint speed.

  public:

    /** Get read-write access to mStictionCoef. */
    double& StictionCoefficient() { return mStictionCoef; };
    /** Get read-only access to mStictionCoef. */
    double StictionCoefficient() const { return mStictionCoef; };

    /** Get read-write access to mSlipCoef. */
    double& SlipCoefficient() { return mSlipCoef; };
    /** Get read-only access to mSlipCoef. */
    double SlipCoefficient() const { return mSlipCoef; };

    /** Get read-write access to mSlipVelocity. */
    double& SlipVelocity() { return mSlipVelocity; };
    /** Get read-only access to mSlipVelocity. */
    double SlipVelocity() const { return mSlipVelocity; };

    /**
     * Default constructor.
     */
    vmc_revolute_joint_3D(const std::string& aName = "") : revolute_joint_3D(aName), mStictionCoef(0.5), mSlipCoef(0.3), mSlipVelocity(1E-5) { };

    /**
     * Parametrized constructor.
     * \param aName the name of this KTE model.
     * \param aAngle the generalized coordinate associated with the displacement of this joint.
     * \param aBase the coordinate frame at the base of the joint.
     * \param aEnd the coordinate frame just after the joint transformations are applied.
     * \param aJacobian a pointer to contain the Jacobian frame produced by this joint, default value will disable the Jacobian frame's calculation.
     * \param aStictionCoef the stiction coefficient in Dry Coulomb friction.
     * \param aSlipCoef the slip coefficient in Dry Coulomb friction.
     * \param aSlipVelocity the micro-slip tolerance in joint speed.
     */
    vmc_revolute_joint_3D(const std::string& aName,
                          const shared_pointer< gen_coord<double> >::type& aAngle,
			  const vect<double,3>& aAxis,
			  const shared_pointer< frame_3D<double> >::type& aBase,
			  const shared_pointer< frame_3D<double> >::type& aEnd,
			  const shared_pointer< jacobian_gen_3D<double> >::type& aJacobian = shared_pointer< jacobian_gen_3D<double> >::type(),
			  double aStictionCoef = 0.5,
			  double aSlipCoef = 0.3,
			  double aSlipVelocity = 1E-5) :
			  revolute_joint_3D(aName,aAngle,aAxis,aBase,aEnd,aJacobian),
			  mStictionCoef(aStictionCoef),
			  mSlipCoef(aSlipCoef),
			  mSlipVelocity(aSlipVelocity) { };

    /**
     * Default destructor.
     */
    virtual ~vmc_revolute_joint_3D() { };

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type());

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      revolute_joint_3D::save(A,revolute_joint_3D::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mStictionCoef)
        & RK_SERIAL_SAVE_WITH_NAME(mSlipCoef)
        & RK_SERIAL_SAVE_WITH_NAME(mSlipVelocity);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      revolute_joint_3D::load(A,revolute_joint_3D::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mStictionCoef)
        & RK_SERIAL_LOAD_WITH_NAME(mSlipCoef)
        & RK_SERIAL_LOAD_WITH_NAME(mSlipVelocity);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(vmc_revolute_joint_3D,0xC2100027,1,"vmc_revolute_joint_3D",revolute_joint_3D)


};


};

};

#endif //VMC_REVOLUTE_JOINT_HPP










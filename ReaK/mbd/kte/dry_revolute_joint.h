/**
 * \file dry_revolute_joint.h
 *
 * This library declares KTE models for revolute joints which also incorporate dry-friction modelling.
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

#ifndef REAK_MBD_KTE_DRY_REVOLUTE_JOINT_H_
#define REAK_MBD_KTE_DRY_REVOLUTE_JOINT_H_


#include "ReaK/mbd/kte/revolute_joint.h"

namespace ReaK::kte {

/**
 * This class implements the dry-frictional revolute joint model for 2D. Dry friction model uses the
 * micro-slip model of Coulomb friction.
 */
class dry_revolute_joint_2D : public revolute_joint_2D {
 protected:
  double
      mStictionCoef;  ///< Holds the stiction coefficient in Dry Coulomb friction.
  double mSlipCoef;  ///< Holds the slip coefficient in Dry Coulomb friction.
  double mSlipVelocity;  ///< Holds the micro-slip tolerance in joint speed.

 public:
  /**
   * Sets the stiction coefficient of the dry-revolute joint.
   * \param aValue The new stiction coefficient of the dry-revolute joint.
   */
  void setStictionCoefficient(double aValue) { mStictionCoef = aValue; }
  /**
   * Returns the stiction coefficient of the dry-revolute joint.
   * \return The stiction coefficient of the dry-revolute joint.
   */
  double StictionCoefficient() const { return mStictionCoef; }

  /**
   * Sets the slip coefficient of the dry-revolute joint.
   * \param aValue The new slip coefficient of the dry-revolute joint.
   */
  void setSlipCoefficient(double aValue) { mSlipCoef = aValue; }
  /**
   * Returns the slip coefficient of the dry-revolute joint.
   * \return The slip coefficient of the dry-revolute joint.
   */
  double SlipCoefficient() const { return mSlipCoef; }

  /**
   * Sets the slip velocity of the dry-revolute joint.
   * \param aValue The new slip velocity of the dry-revolute joint.
   */
  void setSlipVelocity(double aValue) { mSlipVelocity = aValue; }
  /**
   * Returns the slip velocity of the dry-revolute joint.
   * \return The slip velocity of the dry-revolute joint.
   */
  double SlipVelocity() const { return mSlipVelocity; }

  /**
  * Default constructor.
  */
  explicit dry_revolute_joint_2D(const std::string& aName = "")
      : revolute_joint_2D(aName),
        mStictionCoef(0.5),
        mSlipCoef(0.3),
        mSlipVelocity(1E-5) {}

  /**
   * Parametrized constructor.
   * \param aName the name of this KTE model.
   * \param aAngle the generalized coordinate associated with the displacement of this joint.
   * \param aBase the coordinate frame at the base of the joint.
   * \param aEnd the coordinate frame just after the joint transformations are applied.
   * \param aJacobian a pointer to contain the Jacobian frame produced by this joint, default value will disable the
   * Jacobian frame's calculation.
   * \param aStictionCoef the stiction coefficient in Dry Coulomb friction.
   * \param aSlipCoef the slip coefficient in Dry Coulomb friction.
   * \param aSlipVelocity the micro-slip tolerance in joint speed.
   */
  dry_revolute_joint_2D(
      const std::string& aName,
      const std::shared_ptr<gen_coord<double>>& aAngle,
      const std::shared_ptr<frame_2D<double>>& aBase,
      const std::shared_ptr<frame_2D<double>>& aEnd,
      const std::shared_ptr<jacobian_gen_2D<double>>& aJacobian =
          std::shared_ptr<jacobian_gen_2D<double>>(),
      double aStictionCoef = 0.5, double aSlipCoef = 0.3,
      double aSlipVelocity = 1E-5)
      : revolute_joint_2D(aName, aAngle, aBase, aEnd, aJacobian),
        mStictionCoef(aStictionCoef),
        mSlipCoef(aSlipCoef),
        mSlipVelocity(aSlipVelocity) {}

  /**
   * Default destructor.
   */
  ~dry_revolute_joint_2D() override = default;

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    revolute_joint_2D::save(
        A, revolute_joint_2D::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(mStictionCoef) &
        RK_SERIAL_SAVE_WITH_NAME(mSlipCoef) &
        RK_SERIAL_SAVE_WITH_NAME(mSlipVelocity);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    revolute_joint_2D::load(
        A, revolute_joint_2D::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(mStictionCoef) &
        RK_SERIAL_LOAD_WITH_NAME(mSlipCoef) &
        RK_SERIAL_LOAD_WITH_NAME(mSlipVelocity);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(dry_revolute_joint_2D, 0xC210002A, 1,
                              "dry_revolute_joint_2D", revolute_joint_2D)
};

/**
 * This class implements the dry-frictional revolute joint model for 3D. Dry friction model uses the
 * micro-slip model of Coulomb friction.
 */
class dry_revolute_joint_3D : public revolute_joint_3D {
 protected:
  double
      mStictionCoef;  ///< Holds the stiction coefficient in Dry Coulomb friction.
  double mSlipCoef;  ///< Holds the slip coefficient in Dry Coulomb friction.
  double mSlipVelocity;  ///< Holds the micro-slip tolerance in joint speed.

 public:
  /**
   * Sets the stiction coefficient of the dry-revolute joint.
   * \param aValue The new stiction coefficient of the dry-revolute joint.
   */
  void setStictionCoefficient(double aValue) { mStictionCoef = aValue; }
  /**
   * Returns the stiction coefficient of the dry-revolute joint.
   * \return The stiction coefficient of the dry-revolute joint.
   */
  double StictionCoefficient() const { return mStictionCoef; }

  /**
   * Sets the slip coefficient of the dry-revolute joint.
   * \param aValue The new slip coefficient of the dry-revolute joint.
   */
  void setSlipCoefficient(double aValue) { mSlipCoef = aValue; }
  /**
   * Returns the slip coefficient of the dry-revolute joint.
   * \return The slip coefficient of the dry-revolute joint.
   */
  double SlipCoefficient() const { return mSlipCoef; }

  /**
   * Sets the slip velocity of the dry-revolute joint.
   * \param aValue The new slip velocity of the dry-revolute joint.
   */
  void setSlipVelocity(double aValue) { mSlipVelocity = aValue; }
  /**
   * Returns the slip velocity of the dry-revolute joint.
   * \return The slip velocity of the dry-revolute joint.
   */
  double SlipVelocity() const { return mSlipVelocity; }

  /**
   * Default constructor.
   */
  explicit dry_revolute_joint_3D(const std::string& aName = "")
      : revolute_joint_3D(aName),
        mStictionCoef(0.5),
        mSlipCoef(0.3),
        mSlipVelocity(1E-5) {}

  /**
   * Parametrized constructor.
   * \param aName the name of this KTE model.
   * \param aAngle the generalized coordinate associated with the displacement of this joint.
   * \param aBase the coordinate frame at the base of the joint.
   * \param aEnd the coordinate frame just after the joint transformations are applied.
   * \param aJacobian a pointer to contain the Jacobian frame produced by this joint, default value will disable the
   * Jacobian frame's calculation.
   * \param aStictionCoef the stiction coefficient in Dry Coulomb friction.
   * \param aSlipCoef the slip coefficient in Dry Coulomb friction.
   * \param aSlipVelocity the micro-slip tolerance in joint speed.
   */
  dry_revolute_joint_3D(
      const std::string& aName,
      const std::shared_ptr<gen_coord<double>>& aAngle,
      const vect<double, 3>& aAxis,
      const std::shared_ptr<frame_3D<double>>& aBase,
      const std::shared_ptr<frame_3D<double>>& aEnd,
      const std::shared_ptr<jacobian_gen_3D<double>>& aJacobian =
          std::shared_ptr<jacobian_gen_3D<double>>(),
      double aStictionCoef = 0.5, double aSlipCoef = 0.3,
      double aSlipVelocity = 1E-5)
      : revolute_joint_3D(aName, aAngle, aAxis, aBase, aEnd, aJacobian),
        mStictionCoef(aStictionCoef),
        mSlipCoef(aSlipCoef),
        mSlipVelocity(aSlipVelocity) {}

  /**
   * Default destructor.
   */
  ~dry_revolute_joint_3D() override = default;

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    revolute_joint_3D::save(
        A, revolute_joint_3D::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(mStictionCoef) &
        RK_SERIAL_SAVE_WITH_NAME(mSlipCoef) &
        RK_SERIAL_SAVE_WITH_NAME(mSlipVelocity);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    revolute_joint_3D::load(
        A, revolute_joint_3D::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(mStictionCoef) &
        RK_SERIAL_LOAD_WITH_NAME(mSlipCoef) &
        RK_SERIAL_LOAD_WITH_NAME(mSlipVelocity);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(dry_revolute_joint_3D, 0xC210002B, 1,
                              "dry_revolute_joint_3D", revolute_joint_3D)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_DRY_REVOLUTE_JOINT_H_

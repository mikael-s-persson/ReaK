/**
 * \file frame_3D.h
 *
 * This library implements the frame_3D class template which can be used to represent a 3D kinetostatic
 * frame of reference (for kinematics transformations mainly).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2011 (last revision)
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

#ifndef REAK_MATH_KINETOSTATICS_FRAME_3D_H_
#define REAK_MATH_KINETOSTATICS_FRAME_3D_H_

#include "ReaK/math/kinetostatics/pose_3D.h"

namespace ReaK {

/**
 * This class extends pose_3D to include velocity and acceleration as well as applied forces.
 * \note Linear kinematics are expressed in Parent coordinates while rotation kinematics are
 *       expressed in this coordinate system (local or "body-fixed"). However, all forces
 *       (force and torque) are expressed in local coordinates.
 */
template <class T>
class frame_3D : public pose_3D<T> {
 public:
  using self = frame_3D<T>;
  using base = pose_3D<T>;
  using value_type = T;

  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;

  using position_type = typename base::position_type;
  using vector_type = typename base::vector_type;
  using rotation_type = typename base::rotation_type;

  using linear_vector_type = vect<T, 3>;
  using angular_vector_type = vect<T, 3>;
  using angular_rate_type = vect<T, 4>;

  /// Velocity vector of this kinematic frame, relative-to and expressed-in parent coordinates.
  linear_vector_type Velocity;
  /// Angular velocity of this kinematic frame, relative to parent coordinates,
  /// expressed in this coordinate system (local).
  angular_vector_type AngVelocity;

  /// Acceleration of this kinematic frame, relative-to and expressed-in parent coordinates.
  linear_vector_type Acceleration;
  /// Angular Acceleration of this kinematic frame, relative to parent
  /// coordinates, expressed in this coordinate system (local).
  angular_vector_type AngAcceleration;

  /// Force flowing through this frame, expressed in this coordinate system (local).
  linear_vector_type Force;
  /// Torque flowing through this frame, expressed in this coordinate system (local).
  angular_vector_type Torque;

  /// Quaternion time-derivative, relative to parent coordinate, expressed in this
  /// coordinate system (local).
  angular_rate_type QuatDot;

  /**
   * Default constructor, all is set to zero.
   */
  frame_3D()
      : base(),
        Velocity(),
        AngVelocity(),
        Acceleration(),
        AngAcceleration(),
        Force(),
        Torque(),
        QuatDot() {}

  /**
   * Parametrized constructor, all is set to corresponding parameters.
   */
  frame_3D(const std::weak_ptr<base>& aParent, const position_type& aPosition,
           const rotation_type& aQuat, const linear_vector_type& aVelocity,
           const angular_vector_type& aAngVelocity,
           const linear_vector_type& aAcceleration,
           const angular_vector_type& aAngAcceleration,
           const linear_vector_type& aForce, const angular_vector_type& aTorque)
      : base(aParent, aPosition, aQuat),
        Velocity(aVelocity),
        AngVelocity(aAngVelocity),
        Acceleration(aAcceleration),
        AngAcceleration(aAngAcceleration),
        Force(aForce),
        Torque(aTorque) {
    UpdateQuatDot();
  }

  /**
   * Copy-constructor.
   */
  frame_3D(const self& aFrame)
      : base(aFrame),
        Velocity(aFrame.Velocity),
        AngVelocity(aFrame.AngVelocity),
        Acceleration(aFrame.Acceleration),
        AngAcceleration(aFrame.AngAcceleration),
        Force(aFrame.Force),
        Torque(aFrame.Torque) {
    UpdateQuatDot();
  }

  /**
   * Explicit conversion from static 3D pose.
   */
  explicit frame_3D(const base& Pose_)
      : base(Pose_),
        Velocity(),
        AngVelocity(),
        Acceleration(),
        AngAcceleration(),
        Force(),
        Torque(),
        QuatDot() {}

  /**
   * Default destructor.
   */
  ~frame_3D() override = default;

 protected:
  self getFrameRelativeToImpl(const base* F) const {
    if (!F) {
      return getGlobalFrame();
    }
    // If this is at the global node, F can meet this there.
    if (this->Parent.expired()) {
      if (rtti::rk_is_of_type<self>(F)) {
        return (~(static_cast<const self*>(F)->getGlobalFrame())) * (*this);
      }
      return (~(self(*F).getGlobalFrame())) * (*this);
    }
    // If F is somewhere down "this"'s chain.
    if (this->isParentPoseImpl(F)) {
      std::shared_ptr<base> p = this->Parent.lock();
      if (p.get() == F) {
        return *this;
      }
      if (rtti::rk_is_of_type<self>(p)) {
        return static_cast<const self*>(p.get())->getFrameRelativeToImpl(F) *
               (*this);
      }
      return self(*p).getFrameRelativeToImpl(F) * (*this);
    }
    // If this is somewhere down F's chain.
    if (F->isParentPose(*this)) {
      if (rtti::rk_is_of_type<self>(F)) {
        return ~(static_cast<const self*>(F)->getFrameRelativeToImpl(this));
      }
      return ~(self(*F).getFrameRelativeToImpl(this));
    }
    // Else means F's chain meets "this"'s chain somewhere down, possibly all the way to the global node.
    std::shared_ptr<base> p = this->Parent.lock();
    if (rtti::rk_is_of_type<self>(p)) {
      return static_cast<const self*>(p.get())->getFrameRelativeToImpl(F) *
             (*this);
    }
    return self(*p).getFrameRelativeToImpl(F) * (*this);
  }

 public:
  /**
   * Updates the quaternion time-derivative to correspond to current angular velocity.
   */
  void UpdateQuatDot() { QuatDot = this->Quat.getQuaternionDot(AngVelocity); }

  /**
   * Returns this coordinate frame relative to the global inertial frame.
   * \note no operations are performed on forces.
   */
  self getGlobalFrame() const {
    if (!this->Parent.expired()) {
      self result;
      std::shared_ptr<const self> p =
          rtti::rk_dynamic_ptr_cast<const self>(this->Parent.lock());
      if (p) {
        result = p->getGlobalFrame();
      } else {
        result = self(*(this->Parent.lock())).getGlobalFrame();
      }

      rot_mat_3D<value_type> R(result.Quat.getRotMat());
      result.Position += R * this->Position;
      result.Velocity += R * ((result.AngVelocity % this->Position) + Velocity);
      result.Acceleration +=
          R * ((result.AngVelocity % (result.AngVelocity % this->Position)) +
               value_type(2.0) * (result.AngVelocity % Velocity) +
               (result.AngAcceleration % this->Position) + Acceleration);

      rot_mat_3D<value_type> R2(invert(this->Quat).getRotMat());
      result.Quat *= this->Quat;
      result.AngVelocity = (R2 * result.AngVelocity) + AngVelocity;
      result.AngAcceleration = (R2 * result.AngAcceleration) +
                               ((R2 * result.AngVelocity) % AngVelocity) +
                               AngAcceleration;

      result.UpdateQuatDot();

      return result;
    }
    return *this;
  }

  /**
   * Returns this coordinate frame relative to the frame or pose F.
   * \note No operations are performed on forces. F is tested for being
   *       castablet to a frame or not.
   */
  self getFrameRelativeTo(const base& F) const {
    return getFrameRelativeToImpl(&F);
  }

  /**
   * Returns this coordinate frame relative to the frame or pose F.
   * \note No operations are performed on forces. F is tested for being
   *       castablet to a frame or not.
   */
  self getFrameRelativeTo(const std::shared_ptr<const base>& F) const {
    if (!F) {
      return getGlobalFrame();
    }
    return getFrameRelativeToImpl(F.get());
  }

  /**
   * Adds frame Frame_ before this coordinate frame ("before" is meant in the same sense as for pose_3D::addBefore()).
   * The transformation uses classic "rotating frame" formulae.
   * \note No operations are performed on forces.
   */
  self& addBefore(const self& aFrame) {

    rot_mat_3D<value_type> R(this->Quat.getRotMat());
    this->Position += R * aFrame.Position;
    Velocity += R * ((AngVelocity % aFrame.Position) + aFrame.Velocity);
    Acceleration +=
        R * ((AngVelocity % (AngVelocity % aFrame.Position)) +
             T(2.0) * (AngVelocity % aFrame.Velocity) +
             (AngAcceleration % aFrame.Position) + aFrame.Acceleration);

    rot_mat_3D<value_type> R2(aFrame.Quat.getRotMat());
    this->Quat *= aFrame.Quat;
    AngAcceleration = (AngAcceleration * R2) +
                      ((AngVelocity * R2) % aFrame.AngVelocity) +
                      aFrame.AngAcceleration;
    AngVelocity = (AngVelocity * R2) + aFrame.AngVelocity;

    UpdateQuatDot();

    return *this;
  }

  self& addBefore(const base& aPose) {

    rot_mat_3D<value_type> R(this->Quat.getRotMat());
    this->Position += R * aPose.Position;
    Velocity += R * (AngVelocity % aPose.Position);
    Acceleration += R * ((AngVelocity % (AngVelocity % aPose.Position)) +
                         (AngAcceleration % aPose.Position));

    rot_mat_3D<value_type> R2(aPose.Quat.getRotMat());
    this->Quat *= aPose.Quat;
    AngAcceleration = (AngAcceleration * R2);
    AngVelocity = (AngVelocity * R2);

    UpdateQuatDot();

    return *this;
  }

  /**
   * Adds frame Frame_ after this coordinate frame ("after" is meant in the same sense as for pose_3D::addAfter()).
   * The transformation uses classic "rotating frame" formulae.
   * \note No operations are performed on forces.
   */
  self& addAfter(const self& aFrame) {

    rot_mat_3D<value_type> R(aFrame.Quat.getRotMat());
    Acceleration =
        aFrame.Acceleration +
        R * ((aFrame.AngVelocity % (aFrame.AngVelocity % this->Position)) +
             value_type(2.0) * (aFrame.AngVelocity % Velocity) +
             (aFrame.AngAcceleration % this->Position) + Acceleration);
    Velocity = aFrame.Velocity +
               R * ((aFrame.AngVelocity % this->Position) + Velocity);
    this->Position = aFrame.Position + R * this->Position;

    rot_mat_3D<value_type> R2(this->Quat.getRotMat());
    AngAcceleration += (aFrame.AngAcceleration * R2) +
                       ((aFrame.AngVelocity * R2) % AngVelocity);
    AngVelocity += (aFrame.AngVelocity * R2);
    this->Quat = aFrame.Quat * this->Quat;

    this->Parent = aFrame.Parent;

    UpdateQuatDot();

    return *this;
  }

  self& addAfter(const base& aPose) {

    rot_mat_3D<value_type> R(aPose.Quat.getRotMat());
    Acceleration = R * (Acceleration);
    Velocity = R * (Velocity);
    this->Position = aPose.Position + R * this->Position;

    this->Quat = aPose.Quat * this->Quat;

    this->Parent = aPose.Parent;

    UpdateQuatDot();

    return *this;
  }

  /**
   * Assignment operator.
   */
  self& operator=(const self& F) {

    this->Parent = F.Parent;
    this->Position = F.Position;
    this->Quat = F.Quat;
    Velocity = F.Velocity;
    AngVelocity = F.AngVelocity;
    Acceleration = F.Acceleration;
    AngAcceleration = F.AngAcceleration;

    UpdateQuatDot();
    return *this;
  }

  /**
   * Assignment operator.
   */
  self& operator=(const base& P) {

    this->Parent = P.Parent;
    this->Position = P.Position;
    this->Quat = P.Quat;
    Velocity = linear_vector_type();
    Acceleration = linear_vector_type();
    Force = linear_vector_type();
    AngVelocity = angular_vector_type();
    AngAcceleration = angular_vector_type();
    Torque = angular_vector_type();

    UpdateQuatDot();
    return *this;
  }

  /**
   * Multiplication-assignment operator, equivalent to "this->addBefore( F )".
   * \note No operations are performed on forces.
   */
  self& operator*=(const self& F) { return addBefore(F); }

  /**
   * Multiplication-assignment operator, equivalent to "this->addBefore( P )".
   * \note No operations are performed on forces.
   */
  self& operator*=(const base& P) { return addBefore(P); }

  /**
   * Multiplication operator, equivalent to "result = *this; result->addBefore( F )".
   * \note No operations are performed on forces.
   */
  friend self operator*(self F1, const self& F2) {
    F1 *= F2;
    return F1;
  }

  /**
   * Multiplication operator, equivalent to "result = *this; result->addBefore( P )".
   * \note No operations are performed on forces.
   */
  friend self operator*(self F, const base& P) {
    F *= P;
    return F;
  }

  friend self operator*(const base& P, const self& F) {
    self result(P);
    result *= F;
    return result;
  }

  /**
   * Inversion operator, i.e. "this->addBefore( ~this ) == Parent".
   * \note Forces are negated and rotated.
   */
  self operator~() {
    rot_mat_3D<value_type> R(this->Quat.getRotMat());
    self result;
    result.Quat = invert(this->Quat);
    result.AngVelocity = R * (-AngVelocity);
    result.AngAcceleration = R * (-AngAcceleration);
    result.Position = (-this->Position) * R;
    result.Velocity = -((result.AngVelocity % this->Position) + Velocity) * R;
    result.Acceleration =
        -((result.AngVelocity % (result.AngVelocity % this->Position)) +
          (value_type(2.0) * result.AngVelocity % Velocity) +
          (result.AngAcceleration % this->Position) + Acceleration) *
        R;
    result.Force = -Force * R;
    result.Torque = R * (-Torque);  // action-reaction
    return result;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    pose_3D<T>::save(A, pose_3D<T>::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(Velocity) &
        RK_SERIAL_SAVE_WITH_NAME(AngVelocity) &
        RK_SERIAL_SAVE_WITH_NAME(Acceleration) &
        RK_SERIAL_SAVE_WITH_NAME(AngAcceleration) &
        RK_SERIAL_SAVE_WITH_NAME(Force) & RK_SERIAL_SAVE_WITH_NAME(Torque);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    pose_3D<T>::load(A, pose_3D<T>::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(Velocity) &
        RK_SERIAL_LOAD_WITH_NAME(AngVelocity) &
        RK_SERIAL_LOAD_WITH_NAME(Acceleration) &
        RK_SERIAL_LOAD_WITH_NAME(AngAcceleration) &
        RK_SERIAL_LOAD_WITH_NAME(Force) & RK_SERIAL_LOAD_WITH_NAME(Torque);

    UpdateQuatDot();
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x00000020, 1, "frame_3D", base)
};

template <typename T>
std::ostream& operator<<(std::ostream& out, const frame_3D<T>& g) {
  out << "(Position = " << g.Position << "; Quaternion = " << g.Quat
      << "; Velocity = " << g.Velocity << "; AngVelocity = " << g.AngVelocity
      << "; Acceleration = " << g.Acceleration
      << "; AngAcceleration = " << g.AngAcceleration << "; Force = " << g.Force
      << "; Torque = " << g.Torque << ")";
  return out;
}

}  // namespace ReaK

#endif  // REAK_MATH_KINETOSTATICS_FRAME_3D_H_

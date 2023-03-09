/**
 * \file frame_2D.h
 *
 * This library implements the frame_2D class template which can be used to represent a 2D kinetostatic
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

#ifndef REAK_MATH_KINETOSTATICS_FRAME_2D_H_
#define REAK_MATH_KINETOSTATICS_FRAME_2D_H_

#include "ReaK/math/kinetostatics/pose_2D.h"

namespace ReaK {

/**
 * This class extends pose_2D to include velocity and acceleration as well as applied forces.
 * \note Linear kinematics are expressed in Parent coordinates while rotation kinematics are
 *       expressed in this coordinate system (local or "body-fixed"). However, all forces
 *       (force and torque) are expressed in local coordinates.
 */
template <typename T>
class frame_2D : public pose_2D<T> {
 public:
  using self = frame_2D<T>;
  using base = pose_2D<T>;
  using value_type = T;

  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;

  using position_type = typename base::position_type;
  using vector_type = typename base::vector_type;
  using rotation_type = typename base::rotation_type;

  using linear_vector_type = vect<T, 2>;
  using angular_vector_type = T;

  /// Velocity vector of this kinematic frame, relative-to and expressed-in parent coordinates.
  linear_vector_type Velocity;
  /// Angular velocity of this kinematic frame, relative to parent coordinates.
  angular_vector_type AngVelocity;

  /// Acceleration of this kinematic frame, relative-to and expressed-in parent coordinates.
  linear_vector_type Acceleration;
  /// Angular Acceleration of this kinematic frame, relative to parent coordinates.
  angular_vector_type AngAcceleration;

  /// Force flowing through this frame, expressed in this coordinate system (local).
  linear_vector_type Force;
  /// Torque flowing through this frame.
  angular_vector_type Torque;

  /**
   * Default constructor, all is set to zero.
   */
  frame_2D()
      : pose_2D<T>(),
        Velocity(),
        AngVelocity(0.0),
        Acceleration(),
        AngAcceleration(0.0),
        Force(),
        Torque(0.0) {}

  /**
   * Parametrized constructor, all is set to corresponding parameters.
   */
  frame_2D(const std::weak_ptr<base>& aParent, const position_type& aPosition,
           const rotation_type& aRotation, const linear_vector_type& aVelocity,
           const angular_vector_type& aAngVelocity,
           const linear_vector_type& aAcceleration,
           const angular_vector_type& aAngAcceleration,
           const linear_vector_type& aForce, const angular_vector_type& aTorque)
      : pose_2D<T>(aParent, aPosition, aRotation),
        Velocity(aVelocity),
        AngVelocity(aAngVelocity),
        Acceleration(aAcceleration),
        AngAcceleration(aAngAcceleration),
        Force(aForce),
        Torque(aTorque) {}

  /**
   * Copy-constructor.
   */
  frame_2D(const self& aFrame)
      : pose_2D<T>(aFrame),
        Velocity(aFrame.Velocity),
        AngVelocity(aFrame.AngVelocity),
        Acceleration(aFrame.Acceleration),
        AngAcceleration(aFrame.AngAcceleration),
        Force(aFrame.Force),
        Torque(aFrame.Torque) {}

  /**
   * Explicit conversion from a simple pose, the additional values are set to zero.
   */
  explicit frame_2D(const base& aPose)
      : pose_2D<T>(aPose),
        Velocity(),
        AngVelocity(0.0),
        Acceleration(),
        AngAcceleration(0.0),
        Force(),
        Torque(0.0) {}

  /**
   * Default virtual destructor.
   */
  ~frame_2D() override = default;

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
    // If F is somewhere down "this"'s chain
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

      result.Position += result.Rotation * this->Position;
      result.Velocity +=
          result.Rotation * ((result.AngVelocity % this->Position) + Velocity);
      result.Acceleration +=
          result.Rotation *
          ((-result.AngVelocity * result.AngVelocity) * this->Position +
           (value_type(2.0) * result.AngVelocity) % Velocity +
           result.AngAcceleration % this->Position + Acceleration);

      result.Rotation *= this->Rotation;
      result.AngVelocity += AngVelocity;
      result.AngAcceleration += AngAcceleration;

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
   * Adds frame Frame_ before this coordinate frame ("before" is meant in the same sense as for pose_2D::addBefore()).
   * The transformation uses classic "rotating frame" formulae.
   * \note No operations are performed on forces.
   */
  self& addBefore(const self& aFrame) {

    this->Position += this->Rotation * aFrame.Position;
    Velocity +=
        this->Rotation * ((AngVelocity % aFrame.Position) + aFrame.Velocity);
    Acceleration += this->Rotation *
                    ((-AngVelocity * AngVelocity) * aFrame.Position +
                     (value_type(2.0) * AngVelocity) % aFrame.Velocity +
                     AngAcceleration % aFrame.Position + aFrame.Acceleration);

    this->Rotation *= aFrame.Rotation;
    AngVelocity += aFrame.AngVelocity;
    AngAcceleration += aFrame.AngAcceleration;

    return *this;
  }

  /**
   * Adds frame Frame_ after this coordinate frame ("after" is meant in the same sense as for pose_2D::addAfter()).
   * The transformation uses classic "rotating frame" formulae.
   * \note No operations are performed on forces.
   */
  self& addAfter(const self& aFrame) {

    Acceleration =
        aFrame.Acceleration +
        aFrame.Rotation *
            ((-aFrame.AngVelocity * aFrame.AngVelocity) * this->Position +
             (value_type(2.0) * aFrame.AngVelocity) % Velocity +
             aFrame.AngAcceleration % this->Position + Acceleration);
    Velocity =
        aFrame.Velocity +
        aFrame.Rotation * (aFrame.AngVelocity % this->Position + Velocity);
    this->Position = aFrame.Position + aFrame.Rotation * this->Position;

    this->Rotation *= aFrame.Rotation;
    AngVelocity += aFrame.AngVelocity;
    AngAcceleration += aFrame.AngAcceleration;

    this->Parent = aFrame.Parent;

    return *this;
  }

  /**
   * Assignment operator.
   */
  self& operator=(const self& F) {

    this->Parent = F.Parent;
    this->Position = F.Position;
    this->Rotation = F.Rotation;
    Velocity = F.Velocity;
    AngVelocity = F.AngVelocity;
    Acceleration = F.Acceleration;
    AngAcceleration = F.AngAcceleration;
    Force = F.Force;
    Torque = F.Torque;

    return *this;
  }

  /**
   * Multiplication-assignment operator, equivalent to "this->addBefore( F )".
   * \note No operations are performed on forces.
   */
  self& operator*=(const self& F) {

    this->Position += this->Rotation * F.Position;
    Velocity += this->Rotation * ((AngVelocity % F.Position) + F.Velocity);
    Acceleration +=
        this->Rotation * ((-AngVelocity * AngVelocity) * F.Position +
                          (value_type(2.0) * AngVelocity) % F.Velocity +
                          AngAcceleration % F.Position + F.Acceleration);

    this->Rotation *= F.Rotation;
    AngVelocity += F.AngVelocity;
    AngAcceleration += F.AngAcceleration;

    return *this;
  }

  /**
   * Multiplication-assignment operator, equivalent to "this->addBefore( P )".
   * \note No operations are performed on forces.
   */
  self& operator*=(const base& P) {

    this->Position += this->Rotation * P.Position;
    Velocity += this->Rotation * (AngVelocity % P.Position);
    Acceleration +=
        this->Rotation * ((-AngVelocity * AngVelocity) * P.Position +
                          AngAcceleration % P.Position);

    this->Rotation *= P.Rotation;

    return *this;
  }

  /**
   * Multiplication operator, equivalent to "result = *this; result->addBefore( F )".
   * \note No operations are performed on forces.
   */
  friend self operator*(const self& F1, const self& F2) {
    self result;

    result.Parent = F1.Parent;

    result.Position = F1.Position + F1.Rotation * F2.Position;
    result.Velocity =
        F1.Velocity +
        F1.Rotation * ((F1.AngVelocity % F2.Position) + F2.Velocity);
    result.Acceleration =
        F1.Acceleration +
        F1.Rotation * ((-F1.AngVelocity * F1.AngVelocity) * F2.Position +
                       (value_type(2.0) * F1.AngVelocity) % F2.Velocity +
                       F1.AngAcceleration % F2.Position + F2.Acceleration);

    result.Rotation = F1.Rotation * F2.Rotation;
    result.AngVelocity = F1.AngVelocity + F2.AngVelocity;
    result.AngAcceleration = F1.AngAcceleration + F2.AngAcceleration;

    return result;
  }

  /**
   * Multiplication operator, equivalent to "result = *this; result->addBefore( P )".
   * \note No operations are performed on forces.
   */
  friend self operator*(const self& F, const base& P) {
    self result;

    result.Parent = F.Parent;

    result.Position = F.Position + F.Rotation * P.Position;
    result.Velocity = F.Velocity + F.Rotation * (F.AngVelocity % P.Position);
    result.Acceleration =
        F.Acceleration +
        F.Rotation * ((-F.AngVelocity * F.AngVelocity) * P.Position +
                      F.AngAcceleration % P.Position);

    result.Rotation = F.Rotation * P.Rotation;
    result.AngVelocity = F.AngVelocity;
    result.AngAcceleration = F.AngAcceleration;

    return result;
  }

  /**
   * Multiplication operator, equivalent to "frame_2D<T>(P) * F".
   * \note No operations are performed on forces.
   */
  friend self operator*(const base& P, const self& F) {
    self result;

    result.Parent = P.Parent;

    result.Position = P.Position + P.Rotation * F.Position;
    result.Velocity = P.Rotation * (F.Velocity);
    result.Acceleration = P.Rotation * (F.Acceleration);

    result.Rotation = P.Rotation * F.Rotation;
    result.AngVelocity = F.AngVelocity;
    result.AngAcceleration = F.AngAcceleration;

    return result;
  }

  /**
   * Inversion operator, i.e. "this->addBefore( ~this ) == Parent".
   * \note Forces are negated and rotated.
   */
  self operator~() const {
    return self(this->Parent, (-this->Position) * this->Rotation,
                invert(this->Rotation),
                ((AngVelocity % this->Position) - Velocity) * this->Rotation,
                -AngVelocity,
                ((AngVelocity * AngVelocity * this->Position) +
                 (value_type(2.0) * AngVelocity) % Velocity +
                 AngAcceleration % this->Position - Acceleration) *
                    this->Rotation,
                -AngAcceleration, -Force * this->Rotation, -Torque);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base::save(A, base::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(Velocity) &
        RK_SERIAL_SAVE_WITH_NAME(AngVelocity) &
        RK_SERIAL_SAVE_WITH_NAME(Acceleration) &
        RK_SERIAL_SAVE_WITH_NAME(AngAcceleration) &
        RK_SERIAL_SAVE_WITH_NAME(Force) & RK_SERIAL_SAVE_WITH_NAME(Torque);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    base::load(A, base::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(Velocity) &
        RK_SERIAL_LOAD_WITH_NAME(AngVelocity) &
        RK_SERIAL_LOAD_WITH_NAME(Acceleration) &
        RK_SERIAL_LOAD_WITH_NAME(AngAcceleration) &
        RK_SERIAL_LOAD_WITH_NAME(Force) & RK_SERIAL_LOAD_WITH_NAME(Torque);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x0000001F, 1, "frame_2D", base)
};

template <typename T>
std::ostream& operator<<(std::ostream& out, const frame_2D<T>& g) {
  out << "(Position = " << g.Position << "; Rotation = " << g.Rotation
      << "; Velocity = " << g.Velocity << "; AngVelocity = " << g.AngVelocity
      << "; Acceleration = " << g.Acceleration
      << "; AngAcceleration = " << g.AngAcceleration << "; Force = " << g.Force
      << "; Torque = " << g.Torque << ")";
  return out;
}

}  // namespace ReaK

#endif  // REAK_MATH_KINETOSTATICS_FRAME_2D_H_

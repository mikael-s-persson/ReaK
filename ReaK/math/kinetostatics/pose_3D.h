/**
 * \file pose_3D.h
 *
 * This library provides the class template which can represent a static 3D pose (position and rotation).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2011
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

#ifndef REAK_MATH_KINETOSTATICS_POSE_3D_H_
#define REAK_MATH_KINETOSTATICS_POSE_3D_H_

#include "ReaK/math/kinetostatics/rotations_3D.h"

#include <utility>
#include "ReaK/core/base/shared_object.h"

namespace ReaK {

/// This class represents the pose of a 3D coordinate frame (static).
template <typename T>
class pose_3D : public shared_object {
 public:
  using value_type = T;
  using self = pose_3D<T>;

  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;

  using position_type = vect<T, 3>;
  using vector_type = vect<T, 3>;
  using rotation_type = quaternion<T>;

  /// Holds a weak pointer to the pose relative to which this pose is expressed.
  std::weak_ptr<self> Parent;

  /// Position vector of this coordinate system, expressed in parent coordinates.
  position_type Position;
  /// Rotation quaternion of this coordinate system, expressed in this coordinates (local).
  rotation_type Quat;

  /// Default constructor, all is set to zero.
  pose_3D() noexcept : shared_object(), Parent(), Position(), Quat() {}

  /// Parametrized constructor, all is set to corresponding parameters.
  pose_3D(std::weak_ptr<self> aParent, const position_type& aPosition,
          const rotation_type& aQuat) noexcept
      : shared_object(),
        Parent(std::move(aParent)),
        Position(aPosition),
        Quat(aQuat) {}

  /// Copy-constructor.
  pose_3D(const self& aPose) noexcept = default;

  /// Default virtual destructor.
  ~pose_3D() override = default;

 protected:
  bool isParentPoseImpl(const self* P) const noexcept {
    if (Parent.expired()) {
      return !static_cast<bool>(P);
    }
    if (P == Parent.lock().get()) {
      return true;
    }
    return Parent.lock()->isParentPoseImpl(P);
  }

  self getPoseRelativeToImpl(const self* P) const noexcept {
    if (isParentPoseImpl(P)) {
      if ((Parent.expired()) || (Parent.lock().get() == P)) {
        return *this;
      }
      return Parent.lock()->getPoseRelativeToImpl(P) * (*this);
    }
    if (P->isParentPoseImpl(this)) {
      return ~(P->getPoseRelativeToImpl(this));
    }
    if (Parent.expired()) {
      return (~(P->getGlobalPose())) * (*this);
    }
    return Parent.lock()->getPoseRelativeToImpl(P) * (*this);
  }

 public:
  /// Returns this 2D pose relative to the global (null) coordinate system.
  self getGlobalPose() const noexcept {
    if (!Parent.expired()) {
      self result = Parent.lock()->getGlobalPose();
      result.Position += result.Quat * Position;
      result.Quat *= Quat;
      return result;
    }
    return *this;
  }

  /// Returns true if P is part of the parent chain from this pose.
  bool isParentPose(const self& P) const noexcept {
    return isParentPoseImpl(&P);
  }

  /// Returns true if P is part of the parent chain from this pose.
  bool isParentPose(const std::shared_ptr<const self>& P) const noexcept {
    return isParentPoseImpl(P.get());
  }

  /// Returns this 3D pose relative to pose P.
  self getPoseRelativeTo(const self& P) const noexcept {
    return getPoseRelativeToImpl(&P);
  }

  /// Returns this 3D pose relative to pose P.
  self getPoseRelativeTo(const std::shared_ptr<const self>& P) const noexcept {
    return getPoseRelativeToImpl(P.get());
  }

  /// Returns the free vector V (expressed in this coordinate system) expressed in the parent coordinate system.
  vector_type rotateToParent(const vector_type& V) const noexcept {
    return Quat * V;
  }

  /// Returns the free vector V (expressed in this coordinate system) expressed in the global coordinate system.
  vector_type rotateToGlobal(const vector_type& V) const noexcept {
    return getGlobalPose().Quat * V;
  }

  /// Returns the free vector V (expressed in the parent coordinate system) expressed in this coordinate system.
  vector_type rotateFromParent(const vector_type& V) const noexcept {
    return invert(Quat) * V;
  }

  /// Returns the free vector V (expressed in the global coordinate system) expressed in this coordinate system.
  vector_type rotateFromGlobal(const vector_type& V) const noexcept {
    return invert(getGlobalPose().Quat) * V;
  }

  /// Returns the position vector V (expressed in this coordinate system) expressed in the parent coordinate system.
  position_type transformToParent(const position_type& V) const noexcept {
    return Position + Quat * V;
  }

  /// Returns the position vector V (expressed in this coordinate system) expressed in the global coordinate system.
  position_type transformToGlobal(const position_type& V) const noexcept {
    return getGlobalPose().transformToParent(V);
  }

  /// Returns the position vector V (expressed in the parent coordinate system) expressed in this coordinate system.
  position_type transformFromParent(const position_type& V) const noexcept {
    return invert(Quat) * (V - Position);
  }

  /// Returns the position vector V (expressed in the global coordinate system) expressed in this coordinate system.
  position_type transformFromGlobal(const position_type& V) const noexcept {
    return getGlobalPose().transformFromParent(V);
  }

  /// Adds the coordinate tranform of Pose_ before this coordinate transform.
  /// \pre if "V == this->transformToParent( Pose_.transformToParent( U ) )" before
  /// \post then "V == this->transformToParent( U )" after.
  /// \note ignores the parent of Pose_.
  self& addBefore(const self& aPose) noexcept {
    Position += Quat * aPose.Position;
    Quat *= aPose.Quat;
    return *this;
  }

  /// Adds the coordinate tranform of Pose_ after this coordinate transform.
  /// \pre if "V == Pose_.transfromToParent( this->transformToParent( U ) )" before
  /// \post then "V == this->transformToParent( U )" after.
  /// \note ignores the parent of this coordinate system.
  self& addAfter(const self& aPose) noexcept {
    Position = aPose.Position + (aPose.Quat * Position);
    Quat = aPose.Quat * Quat;
    Parent = aPose.Parent;
    return *this;
  }

  /// Adds a translation V to this transformation, where V is expressed in the local coordinate system.
  self& translateLocal(const vector_type& V) noexcept {
    Position += rotateToParent(V);
    return *this;
  }

  /// Adds a translation V to this transformation, where V is expressed in the global coordinate system.
  self& translateGlobal(const vector_type& V) noexcept {
    Position += rotateToParent(rotateFromGlobal(V));
    return *this;
  }

  /// Adds a rotation Q to this transformation, where Q is expressed in local coordinates.
  self& rotateLocal(const rotation_type& Q) noexcept {
    Quat *= Q;
    return *this;
  }

  /// Adds a rotation Q to this transformation, where Q is expressed in global coordinates.
  self& rotateGlobal(const rotation_type& Q) noexcept {
    axis_angle<T> A(Q);
    A.axis() = rotateFromGlobal(A.axis());
    Quat *= A.getQuaternion();
    return *this;
  }

  /// Assignment operator.
  self& operator=(const self& P) noexcept = default;

  /// Multiplication-assignment operator, equivalent to "this->addBefore( P )".
  self& operator*=(const self& P) noexcept { return addBefore(P); }

  /// Multiplication operator, equivalent to "result = *this; result->addBefore( P )".
  friend self operator*(const self& P1, const self& P2) noexcept {
    return self(P1.Parent, P1.Position + (P1.Quat * P2.Position),
                P1.Quat * P2.Quat);
  }

  /// Inversion operator, i.e. "this->addBefore( ~this ) == Parent".
  self operator~() const noexcept {
    return self(Parent, invert(Quat) * (-Position), invert(Quat));
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    if (Parent.expired()) {
      A& RK_SERIAL_SAVE_WITH_ALIAS("Parent", std::shared_ptr<serializable>());
    } else {
      A& RK_SERIAL_SAVE_WITH_ALIAS("Parent", Parent.lock());
    }
    A& RK_SERIAL_SAVE_WITH_NAME(Position) & RK_SERIAL_SAVE_WITH_NAME(Quat);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    std::shared_ptr<pose_3D<T>> tmp;
    A& RK_SERIAL_LOAD_WITH_ALIAS("Parent", tmp) &
        RK_SERIAL_LOAD_WITH_NAME(Position) & RK_SERIAL_LOAD_WITH_NAME(Quat);
    Parent = tmp;
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x0000001E, 1, "pose_3D", shared_object)
};

template <typename T>
std::ostream& operator<<(std::ostream& out, const pose_3D<T>& g) {
  out << "(Position = " << g.Position << "; Quaternion = " << g.Quat << ")";
  return out;
}

}  // namespace ReaK

#endif  // REAK_MATH_KINETOSTATICS_POSE_3D_H_

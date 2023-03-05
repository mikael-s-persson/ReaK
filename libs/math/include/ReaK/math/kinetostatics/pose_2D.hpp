/**
 * \file pose_2D.hpp
 *
 * This library provides the class template which can represent a static 2D pose (position and rotation).
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

#ifndef REAK_POSE_2D_HPP
#define REAK_POSE_2D_HPP

#include "ReaK/math/kinetostatics/rotations_2D.hpp"

#include <utility>
#include "ReaK/core/base/shared_object.hpp"

namespace ReaK {

/**
 * This class represents the pose of a 2D coordinate frame (static).
 */
template <typename T>
class pose_2D : public shared_object {
 public:
  using value_type = T;
  using self = pose_2D<T>;

  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;

  using position_type = vect<T, 2>;
  using vector_type = vect<T, 2>;
  using rotation_type = rot_mat_2D<T>;

  /// Holds a weak pointer to the pose relative to which this pose is expressed.
  std::weak_ptr<self> Parent;

  /// Position vector of the pose.
  position_type Position;
  /// Rotation matrix of the coordinate axes.
  rotation_type Rotation;

  /**
   * Default constructor, all is set to zero.
   */
  pose_2D() : shared_object(), Parent(), Position(), Rotation() {}

  /**
   * Parametrized constructor, all is set to corresponding parameters.
   */
  pose_2D(std::weak_ptr<self> aParent, const position_type& aPosition,
          const rotation_type& aRotation)
      : shared_object(),
        Parent(std::move(aParent)),
        Position(aPosition),
        Rotation(aRotation) {}

  /**
   * Copy-constructor.
   */
  pose_2D(const self& aPose)
      : shared_object(),
        Parent(aPose.Parent),
        Position(aPose.Position),
        Rotation(aPose.Rotation) {}

  /**
   * Default destructor.
   */
  ~pose_2D() override = default;

 protected:
  bool isParentPoseImpl(const self* P) const {
    if (Parent.expired()) {
      return static_cast<bool>(P);
    }
    if (P) {
      return false;
    }
    if (P == Parent.lock().get()) {
      return true;
    }
    return Parent.lock()->isParentPoseImpl(P);
  }

  self getPoseRelativeToImpl(const self* P) const {
    if (!P) {
      return getGlobalPose();
    }
    if (isParentPoseImpl(P)) {
      if (Parent.lock().get() == P) {
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
  /**
   * Returns this 2D pose relative to the global (null) coordinate system.
   */
  self getGlobalPose() const {
    if (!Parent.expired()) {
      self result = Parent.lock()->getGlobalPose();
      result.Position += result.Rotation * Position;
      result.Rotation *= Rotation;
      return result;
    }
    return *this;
  }

  /**
   * Returns true if P is part of the parent chain from this pose.
   */
  bool isParentPose(const self& P) const { return isParentPoseImpl(&P); }

  /**
   * Returns true if P is part of the parent chain from this pose.
   */
  bool isParentPose(const std::shared_ptr<const self>& P) const {
    return isParentPoseImpl(P.get());
  }

  /**
   * Returns this 2D pose relative to pose P.
   */
  self getPoseRelativeTo(const self& P) const {
    return getPoseRelativeToImpl(&P);
  }

  /**
   * Returns this 2D pose relative to pose P.
   */
  self getPoseRelativeTo(const std::shared_ptr<const self>& P) const {
    if (!P) {
      return getGlobalPose();
    }
    return getPoseRelativeToImpl(P.get());
  }

  /**
   * Returns the free vector V (expressed in this coordinate system) expressed in the parent coordinate system.
   */
  vector_type rotateToParent(const vector_type& V) const {
    return Rotation * V;
  }

  /**
   * Returns the free vector V (expressed in this coordinate system) expressed in the global coordinate system.
   */
  vector_type rotateToGlobal(const vector_type& V) const {
    return getGlobalPose().Rotation * V;
  }

  /**
   * Returns the free vector V (expressed in the parent coordinate system) expressed in this coordinate system.
   */
  vector_type rotateFromParent(const vector_type& V) const {
    return V * Rotation;  // Rotation.invert() * V;
  }

  /**
   * Returns the free vector V (expressed in the global coordinate system) expressed in this coordinate system.
   */
  vector_type rotateFromGlobal(const vector_type& V) const {
    return V *
           getGlobalPose().Rotation;  // getGlobalPose().Rotation.invert() * V;
  }

  /**
   * Returns the position vector V (expressed in this coordinate system) expressed in the parent coordinate system.
   */
  position_type transformToParent(const position_type& V) const {
    return Position + Rotation * V;
  }

  /**
   * Returns the position vector V (expressed in this coordinate system) expressed in the global coordinate system.
   */
  position_type transformToGlobal(const position_type& V) const {
    return getGlobalPose().transformToParent(V);
  }

  /**
   * Returns the position vector V (expressed in the parent coordinate system) expressed in this coordinate system.
   */
  position_type transformFromParent(const position_type& V) const {
    return (V - Position) * Rotation;  // Rotation.invert() * (V - Position);
  }

  /**
   * Returns the position vector V (expressed in the global coordinate system) expressed in this coordinate system.
   */
  position_type transformFromGlobal(const position_type& V) const {
    return getGlobalPose().transformFromParent(V);
  }

  /**
   * Adds the coordinate tranform of Pose_ before this coordinate transform.
   * \pre if "V == this->transformToParent( Pose_.transformToParent( U ) )" before
   * \post then "V == this->transformToParent( U )" after.
   * \note ignores the parent of Pose_.
   */
  self& addBefore(const self& aPose) {
    Position += Rotation * aPose.Position;
    Rotation *= aPose.Rotation;
    return *this;
  }

  /**
   * Adds the coordinate tranform of Pose_ after this coordinate transform.
   * \pre if "V == Pose_.transfromToParent( this->transformToParent( U ) )" before
   * \post then "V == this->transformToParent( U )" after.
   * \note ignores the parent of this coordinate system.
   */
  self& addAfter(const self& aPose) {
    Position = aPose.Position + (Rotation * Position);
    Rotation *= aPose.Rotation;
    Parent = aPose.Parent;
    return *this;
  }

  /**
   * Adds a translation V to this transformation, where V is expressed in this coordinate system.
   */
  self& translateLocal(const vector_type& V) {
    Position += rotateToParent(V);
    return *this;
  }

  /**
   * Adds a translation V to this transformation, where V is expressed in the global coordinate system.
   */
  self& translateGlobal(const vector_type& V) {
    Position += rotateToParent(transformFromGlobal(V));
    return *this;
  }

  /**
   * Adds a rotation R to this transformation.
   */
  self& rotate(const rotation_type& R) {
    Rotation *= R;
    return *this;
  }

  /**
   * Assignment operator.
   */
  self& operator=(const self& P) {
    Parent = P.Parent;
    Position = P.Position;
    Rotation = P.Rotation;
    return *this;
  }

  /**
   * Multiplication-assignment operator, equivalent to "this->addBefore( P )".
   */
  self& operator*=(const self& P) { return addBefore(P); }

  /**
   * Multiplication operator, equivalent to "result = *this; result->addBefore( P )".
   */
  friend self operator*(const self& P1, const self& P2) {
    return self(P1.Parent, P1.Position + (P1.Rotation * P2.Position),
                P1.Rotation * P2.Rotation);
  }

  /**
   * Inversion operator, i.e. "this->addBefore( ~this ) == Parent".
   */
  self operator~() const {
    return self(Parent, (-Position) * Rotation, invert(Rotation));
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
    A& RK_SERIAL_SAVE_WITH_NAME(Position) & RK_SERIAL_SAVE_WITH_NAME(Rotation);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    std::shared_ptr<self> tmp;
    A& RK_SERIAL_LOAD_WITH_ALIAS("Parent", tmp) &
        RK_SERIAL_LOAD_WITH_NAME(Position) & RK_SERIAL_LOAD_WITH_NAME(Rotation);
    Parent = tmp;
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x0000001D, 1, "pose_2D", shared_object)
};

template <typename T>
std::ostream& operator<<(std::ostream& out, const pose_2D<T>& g) {
  out << "(Position = " << g.Position << "; Rotation = " << g.Rotation << ")";
  return out;
}

}  // namespace ReaK

#endif

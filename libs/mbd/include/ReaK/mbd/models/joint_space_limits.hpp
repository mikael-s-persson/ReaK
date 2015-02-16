/**
 * \file joint_space_limits.hpp
 *
 * This library provides classes to help create and manipulate joint-space topologies in over
 * a joint-space with limits (speed, acceleration, and jerk limits).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_JOINT_SPACE_LIMITS_HPP
#define REAK_JOINT_SPACE_LIMITS_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/named_object.hpp>
#include <ReaK/math/lin_alg/vect_alg.hpp>

namespace ReaK {

namespace kte {


/**
 * This class template stores a set of vectors to represent the rate-limits on the joints
 * of a manipulator. Basically, this class is just a POD class, but it also provides functions
 * to construct a rate-limited joint-space from a normal joint-space, or vice-versa. Also,
 * it can act as a mapping between rate-limited joint coordinates and normal joint coordinates.
 * \tparam T The value type of the underlying joint-space.
 */
template <typename T>
struct joint_limits_collection : public named_object {
  /** Holds the speed limit for all generalized coordinates. */
  vect_n<T> gen_speed_limits;
  /** Holds the acceleration limit for all generalized coordinates. */
  vect_n<T> gen_accel_limits;
  /** Holds the jerk limit for all generalized coordinates. */
  vect_n<T> gen_jerk_limits;
  /** Holds the speed limit for all 2D frames (alternating velocity limit and angular velocity limit). */
  vect_n<T> frame2D_speed_limits;
  /** Holds the acceleration limit for all 2D frames (alternating acceleration limit and angular acceleration limit). */
  vect_n<T> frame2D_accel_limits;
  /** Holds the jerk limit for all 2D frames (alternating jerk limit and angular jerk limit). */
  vect_n<T> frame2D_jerk_limits;
  /** Holds the speed limit for all 3D frames (alternating velocity limit and angular velocity limit). */
  vect_n<T> frame3D_speed_limits;
  /** Holds the acceleration limit for all 3D frames (alternating acceleration limit and angular acceleration limit). */
  vect_n<T> frame3D_accel_limits;
  /** Holds the jerk limit for all 3D frames (alternating jerk limit and angular jerk limit). */
  vect_n<T> frame3D_jerk_limits;

  typedef T value_type;
  typedef joint_limits_collection<T> self;

  /**
   * Default constructor.
   */
  joint_limits_collection(const std::string& aName = "") : named_object() {
    this->setName(aName);
  };

/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(gen_speed_limits)
        & RK_SERIAL_SAVE_WITH_NAME(gen_accel_limits)
        & RK_SERIAL_SAVE_WITH_NAME(gen_jerk_limits)
        & RK_SERIAL_SAVE_WITH_NAME(frame2D_speed_limits)
        & RK_SERIAL_SAVE_WITH_NAME(frame2D_accel_limits)
        & RK_SERIAL_SAVE_WITH_NAME(frame2D_jerk_limits)
        & RK_SERIAL_SAVE_WITH_NAME(frame3D_speed_limits)
        & RK_SERIAL_SAVE_WITH_NAME(frame3D_accel_limits)
        & RK_SERIAL_SAVE_WITH_NAME(frame3D_jerk_limits);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(gen_speed_limits)
        & RK_SERIAL_LOAD_WITH_NAME(gen_accel_limits)
        & RK_SERIAL_LOAD_WITH_NAME(gen_jerk_limits)
        & RK_SERIAL_LOAD_WITH_NAME(frame2D_speed_limits)
        & RK_SERIAL_LOAD_WITH_NAME(frame2D_accel_limits)
        & RK_SERIAL_LOAD_WITH_NAME(frame2D_jerk_limits)
        & RK_SERIAL_LOAD_WITH_NAME(frame3D_speed_limits)
        & RK_SERIAL_LOAD_WITH_NAME(frame3D_accel_limits)
        & RK_SERIAL_LOAD_WITH_NAME(frame3D_jerk_limits);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400011,1,"joint_limits_collection",named_object)

};

};

};

#endif









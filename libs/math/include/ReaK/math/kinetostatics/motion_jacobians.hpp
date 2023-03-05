/**
 * \file motion_jacobians.hpp
 *
 * This library provides a number of classes to represent motion jacobians.
 * These classes use the kinetostatic classes and map the required quantities
 * to describe the motion jacobians between frames (and generalized coordinates).
 * These classes are mostly POD types which hold the motion jacobians.
 * Motion jacobians hold the quantities that can map the velocities and accelerations of
 * one frame to the velocities and accelerations of another.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2011
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

#ifndef REAK_MOTION_JACOBIANS_HPP
#define REAK_MOTION_JACOBIANS_HPP

#include "ReaK/math/kinetostatics/kinetostatics.hpp"

#include <type_traits>
#include <utility>

namespace ReaK {

/**
 * This class template represents the jacobians between two generalized coordinates.
 * \tparam T The value-type.
 */
template <typename T>
class jacobian_gen_gen : public shared_object {
 public:
  using value_type = T;
  using self = jacobian_gen_gen<T>;

  /// Holds how much velocity is generated at output from the input velocity.
  value_type qd_qd;
  /// Holds how much acceleration is generated at output from the input velocity.
  value_type qd_qdd;

  /**
   * Parametrized constructor.
   */
  explicit jacobian_gen_gen(const value_type& aQdQd = value_type(),
                            const value_type& aQdQdd = value_type())
      : qd_qd(aQdQd), qd_qdd(aQdQdd) {}

  self get_jac_relative_to(
      const std::shared_ptr<gen_coord<value_type>>& /*unused*/) const {
    return *this;
  }

  template <typename Matrix1, typename Matrix2>
  void write_to_matrices(Matrix1& Jac, Matrix2& JacDot) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    static_assert(is_fully_writable_matrix_v<Matrix2>);
    if ((Jac.get_row_count() != 1) || (Jac.get_col_count() != 1) ||
        (JacDot.get_row_count() != 1) || (JacDot.get_col_count() != 1)) {
      throw std::range_error(
          "Jacobian and JacobianDot matrices have the wrong size!");
    }

    Jac(0, 0) = qd_qd;
    JacDot(0, 0) = qd_qdd;
  }

  template <typename Matrix1>
  void write_to_matrices(Matrix1& Jac) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    if ((Jac.get_row_count() != 1) || (Jac.get_col_count() != 1)) {
      throw std::range_error("Jacobian matrix has the wrong size!");
    }

    Jac(0, 0) = qd_qd;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(qd_qd) & RK_SERIAL_SAVE_WITH_NAME(qd_qdd);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(qd_qd) & RK_SERIAL_LOAD_WITH_NAME(qd_qdd);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x00000022, 1, "jacobian_gen_gen",
                              shared_object)
};

/**
 * This class template represents the jacobians between generalized coordinate and a 2D frame.
 * \tparam T The value-type.
 */
template <typename T>
class jacobian_gen_2D : public shared_object {
 public:
  using value_type = T;
  using self = jacobian_gen_2D<T>;

  /// Holds the frame to which the jacobians are relative to.
  std::weak_ptr<frame_2D<value_type>> Parent;
  /// Holds how much velocity is generated at output from the input velocity.
  vect<value_type, 2> qd_vel;
  /// Holds how much angular velocity is generated at output from the input velocity.
  value_type qd_avel;
  /// Holds how much acceleration is generated at output from the input velocity.
  vect<value_type, 2> qd_acc;
  /// Holds how much angular acceleration is generated at output from the input velocity.
  value_type qd_aacc;

  /**
   * Parametrized constructor.
   */
  explicit jacobian_gen_2D(
      std::weak_ptr<frame_2D<value_type>> aParent =
          std::weak_ptr<frame_2D<value_type>>(),
      const vect<value_type, 2>& aQdVel = (vect<value_type, 2>()),
      const value_type& aQdAVel = value_type(),
      const vect<value_type, 2>& aQdAcc = (vect<value_type, 2>()),
      const value_type& aQdAAcc = value_type())
      : Parent(std::move(aParent)),
        qd_vel(aQdVel),
        qd_avel(aQdAVel),
        qd_acc(aQdAcc),
        qd_aacc(aQdAAcc) {}

  self get_jac_relative_to(
      const std::shared_ptr<frame_2D<value_type>>& aFrame) const {
    frame_2D<value_type> f2 = aFrame->getFrameRelativeTo(Parent.lock());
    vect<value_type, 2> v_tmp = (qd_avel % f2.Position + qd_vel) * f2.Rotation;
    return self(
        aFrame, v_tmp, qd_avel,
        (qd_avel % f2.Velocity + qd_aacc % f2.Position + qd_acc) * f2.Rotation -
            f2.AngVelocity % v_tmp,
        qd_aacc);
  }

  template <typename Matrix1, typename Matrix2>
  void write_to_matrices(Matrix1& Jac, Matrix2& JacDot) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    static_assert(is_fully_writable_matrix_v<Matrix2>);
    if ((Jac.get_row_count() != 3) || (Jac.get_col_count() != 1) ||
        (JacDot.get_row_count() != 3) || (JacDot.get_col_count() != 1)) {
      throw std::range_error(
          "Jacobian and JacobianDot matrices have the wrong size!");
    }

    Jac(0, 0) = qd_vel[0];
    Jac(1, 0) = qd_vel[1];
    Jac(2, 0) = qd_avel;
    JacDot(0, 0) = qd_acc[0];
    JacDot(1, 0) = qd_acc[1];
    JacDot(2, 0) = qd_aacc;
  }

  template <typename Matrix1>
  void write_to_matrices(Matrix1& Jac) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    if ((Jac.get_row_count() != 3) || (Jac.get_col_count() != 1)) {
      throw std::range_error("Jacobian matrix has the wrong size!");
    }

    Jac(0, 0) = qd_vel[0];
    Jac(1, 0) = qd_vel[1];
    Jac(2, 0) = qd_avel;
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
    A& RK_SERIAL_SAVE_WITH_NAME(qd_vel) & RK_SERIAL_SAVE_WITH_NAME(qd_avel) &
        RK_SERIAL_SAVE_WITH_NAME(qd_acc) & RK_SERIAL_SAVE_WITH_NAME(qd_aacc);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    std::shared_ptr<frame_2D<value_type>> tmp;
    A& RK_SERIAL_LOAD_WITH_ALIAS("Parent", tmp) &
        RK_SERIAL_LOAD_WITH_NAME(qd_vel) & RK_SERIAL_LOAD_WITH_NAME(qd_avel) &
        RK_SERIAL_LOAD_WITH_NAME(qd_acc) & RK_SERIAL_LOAD_WITH_NAME(qd_aacc);
    Parent = tmp;
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x00000023, 1, "jacobian_gen_2D",
                              shared_object)
};

/**
 * This class template represents the jacobians between generalized coordinate and a 3D frame.
 * \tparam T The value-type.
 */
template <typename T>
class jacobian_gen_3D : public shared_object {
 public:
  using value_type = T;
  using self = jacobian_gen_3D<T>;

  /// Holds the frame to which the jacobians are relative to.
  std::weak_ptr<frame_3D<value_type>> Parent;
  /// Holds how much velocity is generated at output from the input velocity.
  vect<value_type, 3> qd_vel;
  /// Holds how much angular velocity is generated at output from the input velocity.
  vect<value_type, 3> qd_avel;
  /// Holds how much acceleration is generated at output from the input velocity.
  vect<value_type, 3> qd_acc;
  /// Holds how much angular acceleration is generated at output from the input velocity.
  vect<value_type, 3> qd_aacc;

  /**
   * Parametrized constructor.
   */
  explicit jacobian_gen_3D(
      std::weak_ptr<frame_3D<value_type>> aParent =
          std::weak_ptr<frame_3D<value_type>>(),
      const vect<value_type, 3>& aQdVel = (vect<value_type, 3>()),
      const vect<value_type, 3>& aQdAVel = (vect<value_type, 3>()),
      const vect<value_type, 3>& aQdAcc = (vect<value_type, 3>()),
      const vect<value_type, 3>& aQdAAcc = (vect<value_type, 3>()))
      : Parent(std::move(aParent)),
        qd_vel(aQdVel),
        qd_avel(aQdAVel),
        qd_acc(aQdAcc),
        qd_aacc(aQdAAcc) {}

  self get_jac_relative_to(
      const std::shared_ptr<frame_3D<value_type>>& aFrame) const {

    frame_3D<value_type> f2 = aFrame->getFrameRelativeTo(Parent.lock());
    rot_mat_3D<value_type> R(f2.Quat.getRotMat());

    vect<value_type, 3> w_tmp = qd_avel * R;
    vect<value_type, 3> v_tmp = (qd_avel % f2.Position + qd_vel) * R;
    return self(aFrame, v_tmp, w_tmp,
                (qd_avel % f2.Velocity + qd_aacc % f2.Position + qd_acc) * R -
                    f2.AngVelocity % v_tmp,
                qd_aacc * R - f2.AngVelocity % w_tmp);
  }

  template <typename Matrix1, typename Matrix2>
  void write_to_matrices(Matrix1& Jac, Matrix2& JacDot) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    static_assert(is_fully_writable_matrix_v<Matrix2>);
    if ((Jac.get_row_count() != 6) || (Jac.get_col_count() != 1) ||
        (JacDot.get_row_count() != 6) || (JacDot.get_col_count() != 1)) {
      throw std::range_error(
          "Jacobian and JacobianDot matrices have the wrong size!");
    }

    Jac(0, 0) = qd_vel[0];
    Jac(1, 0) = qd_vel[1];
    Jac(2, 0) = qd_vel[2];
    Jac(3, 0) = qd_avel[0];
    Jac(4, 0) = qd_avel[1];
    Jac(5, 0) = qd_avel[2];
    JacDot(0, 0) = qd_acc[0];
    JacDot(1, 0) = qd_acc[1];
    JacDot(2, 0) = qd_acc[2];
    JacDot(3, 0) = qd_aacc[0];
    JacDot(4, 0) = qd_aacc[1];
    JacDot(5, 0) = qd_aacc[2];
  }

  template <typename Matrix1>
  void write_to_matrices(Matrix1& Jac) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    if ((Jac.get_row_count() != 6) || (Jac.get_col_count() != 1)) {
      throw std::range_error("Jacobian matrix has the wrong size!");
    }

    Jac(0, 0) = qd_vel[0];
    Jac(1, 0) = qd_vel[1];
    Jac(2, 0) = qd_vel[2];
    Jac(3, 0) = qd_avel[0];
    Jac(4, 0) = qd_avel[1];
    Jac(5, 0) = qd_avel[2];
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
    A& RK_SERIAL_SAVE_WITH_NAME(qd_vel) & RK_SERIAL_SAVE_WITH_NAME(qd_avel) &
        RK_SERIAL_SAVE_WITH_NAME(qd_acc) & RK_SERIAL_SAVE_WITH_NAME(qd_aacc);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    std::shared_ptr<frame_3D<value_type>> tmp;
    A& RK_SERIAL_LOAD_WITH_ALIAS("Parent", tmp) &
        RK_SERIAL_LOAD_WITH_NAME(qd_vel) & RK_SERIAL_LOAD_WITH_NAME(qd_avel) &
        RK_SERIAL_LOAD_WITH_NAME(qd_acc) & RK_SERIAL_LOAD_WITH_NAME(qd_aacc);
    Parent = tmp;
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x00000024, 1, "jacobian_gen_3D",
                              shared_object)
};

/**
 * This class template represents the jacobians between a 2D frame and a generalized coordinate.
 * \tparam T The value-type.
 */
template <typename T>
class jacobian_2D_gen : public shared_object {
 public:
  using value_type = T;
  using self = jacobian_2D_gen<T>;

  /// Holds how much velocity is generated at output from the input velocity.
  vect<value_type, 2> vel_qd;
  /// Holds how much velocity is generated at output from the input angular velocity.
  value_type avel_qd;
  /// Holds how much acceleration is generated at output from the input velocity.
  vect<value_type, 2> vel_qdd;
  /// Holds how much acceleration is generated at output from the input angular velocity.
  value_type avel_qdd;

  /**
   * Parametrized constructor.
   */
  explicit jacobian_2D_gen(
      const vect<value_type, 2>& aVelQd = (vect<value_type, 2>()),
      const value_type& aAVelQd = value_type(),
      const vect<value_type, 2>& aVelQdd = (vect<value_type, 2>()),
      const value_type& aAVelQdd = value_type())
      : vel_qd(aVelQd),
        avel_qd(aAVelQd),
        vel_qdd(aVelQdd),
        avel_qdd(aAVelQdd) {}

  self get_jac_relative_to(
      const std::shared_ptr<gen_coord<value_type>>& /*unused*/) const {
    return *this;
  }

  template <typename Matrix1, typename Matrix2>
  void write_to_matrices(Matrix1& Jac, Matrix2& JacDot) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    static_assert(is_fully_writable_matrix_v<Matrix2>);
    if ((Jac.get_row_count() != 1) || (Jac.get_col_count() != 3) ||
        (JacDot.get_row_count() != 1) || (JacDot.get_col_count() != 3)) {
      throw std::range_error(
          "Jacobian and JacobianDot matrices have the wrong size!");
    }

    Jac(0, 0) = vel_qd[0];
    Jac(0, 1) = vel_qd[1];
    Jac(0, 2) = avel_qd;
    JacDot(0, 0) = vel_qdd[0];
    JacDot(0, 1) = vel_qdd[1];
    JacDot(0, 2) = avel_qdd;
  }

  template <typename Matrix1>
  void write_to_matrices(Matrix1& Jac) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    if ((Jac.get_row_count() != 1) || (Jac.get_col_count() != 3)) {
      throw std::range_error("Jacobian matrix has the wrong size!");
    }

    Jac(0, 0) = vel_qd[0];
    Jac(0, 1) = vel_qd[1];
    Jac(0, 2) = avel_qd;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(vel_qd) & RK_SERIAL_SAVE_WITH_NAME(avel_qd) &
        RK_SERIAL_SAVE_WITH_NAME(vel_qdd) & RK_SERIAL_SAVE_WITH_NAME(avel_qdd);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(vel_qd) & RK_SERIAL_LOAD_WITH_NAME(avel_qd) &
        RK_SERIAL_LOAD_WITH_NAME(vel_qdd) & RK_SERIAL_LOAD_WITH_NAME(avel_qdd);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x00000025, 1, "jacobian_2D_gen",
                              shared_object)
};

/**
 * This class template represents the jacobians between a 2D frame and a 2D frame.
 * \tparam T The value-type.
 */
template <typename T>
class jacobian_2D_2D : public shared_object {
 public:
  using value_type = T;
  using self = jacobian_2D_2D<T>;

  /// Holds the frame to which the jacobians are relative to.
  std::weak_ptr<frame_2D<value_type>> Parent;
  /// Holds how much velocity is generated at output from the input velocity.
  vect<vect<value_type, 2>, 2> vel_vel;
  /// Holds how much angular velocity is generated at output from the input velocity.
  vect<value_type, 2> vel_avel;
  /// Holds how much velocity is generated at output from the input angular velocity.
  vect<value_type, 2> avel_vel;
  /// Holds how much angular velocity is generated at output from the input angular velocity.
  value_type avel_avel;
  /// Holds how much acceleration is generated at output from the input velocity.
  vect<vect<value_type, 2>, 2> vel_acc;
  /// Holds how much angular acceleration is generated at output from the input velocity.
  vect<value_type, 2> vel_aacc;
  /// Holds how much acceleration is generated at output from the input angular velocity.
  vect<value_type, 2> avel_acc;
  /// Holds how much angular acceleration is generated at output from the input angular velocity.
  value_type avel_aacc;

  /**
   * Parametrized constructor.
   */
  explicit jacobian_2D_2D(
      std::weak_ptr<frame_2D<value_type>> aParent =
          std::weak_ptr<frame_2D<value_type>>(),
      const vect<vect<value_type, 2>, 2>& aVelVel =
          (vect<vect<value_type, 2>, 2>()),
      const vect<value_type, 2>& aVelAVel = (vect<value_type, 2>()),
      const vect<value_type, 2>& aAVelVel = (vect<value_type, 2>()),
      const value_type& aAVelAVel = value_type(),
      const vect<vect<value_type, 2>, 2>& aVelAcc =
          (vect<vect<value_type, 2>, 2>()),
      const vect<value_type, 2>& aVelAAcc = (vect<value_type, 2>()),
      const vect<value_type, 2>& aAVelAcc = (vect<value_type, 2>()),
      const value_type& aAVelAAcc = value_type())
      : Parent(std::move(aParent)),
        vel_vel(aVelVel),
        vel_avel(aVelAVel),
        avel_vel(aAVelVel),
        avel_avel(aAVelAVel),
        vel_acc(aVelAcc),
        vel_aacc(aVelAAcc),
        avel_acc(aAVelAcc),
        avel_aacc(aAVelAAcc) {}

  self get_jac_relative_to(
      const std::shared_ptr<frame_2D<value_type>>& aFrame) const {

    frame_2D<value_type> f2 = aFrame->getFrameRelativeTo(Parent.lock());

    vect<vect<value_type, 2>, 2> new_vel_vel(
        (vel_avel[0] % f2.Position + vel_vel[0]) * f2.Rotation,
        (vel_avel[1] % f2.Position + vel_vel[1]) * f2.Rotation);

    vect<value_type, 2> new_avel_vel =
        (avel_avel % f2.Position + avel_vel) * f2.Rotation;

    return self(aFrame,
                new_vel_vel,   // VelVel
                vel_avel,      // VelAVel
                new_avel_vel,  // AVelVel
                avel_avel,     // AVelAVel
                vect<vect<value_type, 2>, 2>(
                    (vel_avel[0] % f2.Velocity + vel_aacc[0] % f2.Position +
                     vel_acc[0]) *
                            f2.Rotation -
                        f2.AngVelocity % new_vel_vel[0],
                    (vel_avel[1] % f2.Velocity + vel_aacc[1] % f2.Position +
                     vel_acc[1]) *
                            f2.Rotation -
                        f2.AngVelocity % new_vel_vel[1]),  // VelAcc
                vel_aacc,                                  // VelAAcc
                (avel_avel % f2.Velocity + avel_aacc % f2.Position + avel_acc) *
                        f2.Rotation -
                    f2.AngVelocity % new_avel_vel,  // AVelAcc
                avel_aacc                           // AVelAAcc
    );
  }

  template <typename Matrix1, typename Matrix2>
  void write_to_matrices(Matrix1& Jac, Matrix2& JacDot) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    static_assert(is_fully_writable_matrix_v<Matrix2>);
    if ((Jac.get_row_count() != 3) || (Jac.get_col_count() != 3) ||
        (JacDot.get_row_count() != 3) || (JacDot.get_col_count() != 3)) {
      throw std::range_error(
          "Jacobian and JacobianDot matrices have the wrong size!");
    }

    Jac(0, 0) = vel_vel[0][0];
    Jac(1, 0) = vel_vel[0][1];
    Jac(2, 0) = vel_avel[0];
    Jac(0, 1) = vel_vel[1][0];
    Jac(1, 1) = vel_vel[1][1];
    Jac(2, 1) = vel_avel[1];
    Jac(0, 2) = avel_vel[0];
    Jac(1, 2) = avel_vel[1];
    Jac(2, 2) = avel_avel;
    JacDot(0, 0) = vel_acc[0][0];
    JacDot(1, 0) = vel_acc[0][1];
    JacDot(2, 0) = vel_aacc[0];
    JacDot(0, 1) = vel_acc[1][0];
    JacDot(1, 1) = vel_acc[1][1];
    JacDot(2, 1) = vel_aacc[1];
    JacDot(0, 2) = avel_acc[0];
    JacDot(1, 2) = avel_acc[1];
    JacDot(2, 2) = avel_aacc;
  }

  template <typename Matrix1>
  void write_to_matrices(Matrix1& Jac) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    if ((Jac.get_row_count() != 3) || (Jac.get_col_count() != 3)) {
      throw std::range_error("Jacobian matrix has the wrong size!");
    }

    Jac(0, 0) = vel_vel[0][0];
    Jac(1, 0) = vel_vel[0][1];
    Jac(2, 0) = vel_avel[0];
    Jac(0, 1) = vel_vel[1][0];
    Jac(1, 1) = vel_vel[1][1];
    Jac(2, 1) = vel_avel[1];
    Jac(0, 2) = avel_vel[0];
    Jac(1, 2) = avel_vel[1];
    Jac(2, 2) = avel_avel;
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
    A& RK_SERIAL_SAVE_WITH_NAME(vel_vel) & RK_SERIAL_SAVE_WITH_NAME(vel_avel) &
        RK_SERIAL_SAVE_WITH_NAME(avel_vel) &
        RK_SERIAL_SAVE_WITH_NAME(avel_avel) &
        RK_SERIAL_SAVE_WITH_NAME(vel_acc) & RK_SERIAL_SAVE_WITH_NAME(vel_aacc) &
        RK_SERIAL_SAVE_WITH_NAME(avel_acc) &
        RK_SERIAL_SAVE_WITH_NAME(avel_aacc);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    std::shared_ptr<frame_2D<value_type>> tmp;
    A& RK_SERIAL_LOAD_WITH_ALIAS("Parent", tmp) &
        RK_SERIAL_LOAD_WITH_NAME(vel_vel) & RK_SERIAL_LOAD_WITH_NAME(vel_avel) &
        RK_SERIAL_LOAD_WITH_NAME(avel_vel) &
        RK_SERIAL_LOAD_WITH_NAME(avel_avel) &
        RK_SERIAL_LOAD_WITH_NAME(vel_acc) & RK_SERIAL_LOAD_WITH_NAME(vel_aacc) &
        RK_SERIAL_LOAD_WITH_NAME(avel_acc) &
        RK_SERIAL_LOAD_WITH_NAME(avel_aacc);
    Parent = tmp;
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x00000026, 1, "jacobian_2D_2D",
                              shared_object)
};

/**
 * This class template represents the jacobians between a 2D frame and a 3D frame.
 * \tparam T The value-type.
 */
template <typename T>
class jacobian_2D_3D : public shared_object {
 public:
  using value_type = T;
  using self = jacobian_2D_3D<T>;

  /// Holds the frame to which the jacobians are relative to.
  std::weak_ptr<frame_3D<value_type>> Parent;
  /// Holds how much velocity is generated at output from the input velocity.
  vect<vect<value_type, 3>, 2> vel_vel;
  /// Holds how much angular velocity is generated at output from the input velocity.
  vect<vect<value_type, 3>, 2> vel_avel;
  /// Holds how much velocity is generated at output from the input angular velocity.
  vect<value_type, 3> avel_vel;
  /// Holds how much angular velocity is generated at output from the input angular velocity.
  vect<value_type, 3> avel_avel;
  /// Holds how much acceleration is generated at output from the input velocity.
  vect<vect<value_type, 3>, 2> vel_acc;
  /// Holds how much angular acceleration is generated at output from the input velocity.
  vect<vect<value_type, 3>, 2> vel_aacc;
  /// Holds how much acceleration is generated at output from the input angular velocity.
  vect<value_type, 3> avel_acc;
  /// Holds how much angular acceleration is generated at output from the input angular velocity.
  vect<value_type, 3> avel_aacc;

  /**
   * Parametrized constructor.
   */
  explicit jacobian_2D_3D(
      std::weak_ptr<frame_3D<value_type>> aParent =
          std::weak_ptr<frame_3D<value_type>>(),
      const vect<vect<value_type, 3>, 2>& aVelVel =
          (vect<vect<value_type, 3>, 2>()),
      const vect<vect<value_type, 3>, 2>& aVelAVel =
          (vect<vect<value_type, 3>, 2>()),
      const vect<value_type, 3>& aAVelVel = (vect<value_type, 3>()),
      const vect<value_type, 3>& aAVelAVel = (vect<value_type, 3>()),
      const vect<vect<value_type, 3>, 2>& aVelAcc =
          (vect<vect<value_type, 3>, 2>()),
      const vect<vect<value_type, 3>, 2>& aVelAAcc =
          (vect<vect<value_type, 3>, 2>()),
      const vect<value_type, 3>& aAVelAcc = (vect<value_type, 3>()),
      const vect<value_type, 3>& aAVelAAcc = (vect<value_type, 3>()))
      : Parent(std::move(aParent)),
        vel_vel(aVelVel),
        vel_avel(aVelAVel),
        avel_vel(aAVelVel),
        avel_avel(aAVelAVel),
        vel_acc(aVelAcc),
        vel_aacc(aVelAAcc),
        avel_acc(aAVelAcc),
        avel_aacc(aAVelAAcc) {}

  self get_jac_relative_to(
      const std::shared_ptr<frame_3D<value_type>>& aFrame) const {

    frame_3D<value_type> f2 = aFrame->getFrameRelativeTo(Parent.lock());
    rot_mat_3D<value_type> R(f2.Quat.getRotMat());

    vect<vect<value_type, 3>, 2> new_vel_avel(vel_avel[0] * R, vel_avel[1] * R);
    vect<vect<value_type, 3>, 2> new_vel_aacc(
        vel_aacc[0] * R - f2.AngVelocity % new_vel_avel[0],
        vel_aacc[1] * R - f2.AngVelocity % new_vel_avel[1]);

    vect<value_type, 3> new_avel_avel(avel_avel * R);
    vect<value_type, 3> new_avel_aacc(avel_aacc * R -
                                      f2.AngVelocity % new_avel_avel);

    vect<vect<value_type, 3>, 2> new_vel_vel(
        (vel_avel[0] % f2.Position + vel_vel[0]) * R,
        (vel_avel[1] % f2.Position + vel_vel[1]) * R);
    vect<vect<value_type, 3>, 2> new_vel_acc(
        (vel_avel[0] % f2.Velocity + vel_aacc[0] % f2.Position + vel_acc[0]) *
                R -
            f2.AngVelocity % new_vel_vel[0],
        (vel_avel[1] % f2.Velocity + vel_aacc[1] % f2.Position + vel_acc[1]) *
                R -
            f2.AngVelocity % new_vel_vel[1]);

    vect<value_type, 3> new_avel_vel((avel_avel % f2.Position + avel_vel) * R);
    vect<value_type, 3> new_avel_acc(
        (avel_avel % f2.Velocity + avel_aacc % f2.Position + avel_acc) * R -
        f2.AngVelocity % new_avel_vel);

    return self(aFrame, new_vel_vel, new_vel_avel, new_avel_vel, new_avel_avel,
                new_vel_acc, new_vel_aacc, new_avel_acc, new_avel_aacc);
  }

  template <typename Matrix1, typename Matrix2>
  void write_to_matrices(Matrix1& Jac, Matrix2& JacDot) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    static_assert(is_fully_writable_matrix_v<Matrix2>);
    if ((Jac.get_row_count() != 6) || (Jac.get_col_count() != 3) ||
        (JacDot.get_row_count() != 6) || (JacDot.get_col_count() != 3)) {
      throw std::range_error(
          "Jacobian and JacobianDot matrices have the wrong size!");
    }

    Jac(0, 0) = vel_vel[0][0];
    Jac(1, 0) = vel_vel[0][1];
    Jac(2, 0) = vel_vel[0][2];
    Jac(3, 0) = vel_avel[0][0];
    Jac(4, 0) = vel_avel[0][1];
    Jac(5, 0) = vel_avel[0][2];
    Jac(0, 1) = vel_vel[1][0];
    Jac(1, 1) = vel_vel[1][1];
    Jac(2, 1) = vel_vel[1][2];
    Jac(3, 1) = vel_avel[1][0];
    Jac(4, 1) = vel_avel[1][1];
    Jac(5, 1) = vel_avel[1][2];
    Jac(0, 2) = avel_vel[0];
    Jac(1, 2) = avel_vel[1];
    Jac(2, 2) = avel_vel[2];
    Jac(3, 2) = avel_avel[0];
    Jac(4, 2) = avel_avel[1];
    Jac(5, 2) = avel_avel[2];

    JacDot(0, 0) = vel_acc[0][0];
    JacDot(1, 0) = vel_acc[0][1];
    JacDot(2, 0) = vel_acc[0][2];
    JacDot(3, 0) = vel_aacc[0][0];
    JacDot(4, 0) = vel_aacc[0][1];
    JacDot(5, 0) = vel_aacc[0][2];
    JacDot(0, 1) = vel_acc[1][0];
    JacDot(1, 1) = vel_acc[1][1];
    JacDot(2, 1) = vel_acc[1][2];
    JacDot(3, 1) = vel_aacc[1][0];
    JacDot(4, 1) = vel_aacc[1][1];
    JacDot(5, 1) = vel_aacc[1][2];
    JacDot(0, 2) = avel_acc[0];
    JacDot(1, 2) = avel_acc[1];
    JacDot(2, 2) = avel_acc[2];
    JacDot(3, 2) = avel_aacc[0];
    JacDot(4, 2) = avel_aacc[1];
    JacDot(5, 2) = avel_aacc[2];
  }

  template <typename Matrix1>
  void write_to_matrices(Matrix1& Jac) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    if ((Jac.get_row_count() != 6) || (Jac.get_col_count() != 3)) {
      throw std::range_error("Jacobian matrix has the wrong size!");
    }

    Jac(0, 0) = vel_vel[0][0];
    Jac(1, 0) = vel_vel[0][1];
    Jac(2, 0) = vel_vel[0][2];
    Jac(3, 0) = vel_avel[0][0];
    Jac(4, 0) = vel_avel[0][1];
    Jac(5, 0) = vel_avel[0][2];
    Jac(0, 1) = vel_vel[1][0];
    Jac(1, 1) = vel_vel[1][1];
    Jac(2, 1) = vel_vel[1][2];
    Jac(3, 1) = vel_avel[1][0];
    Jac(4, 1) = vel_avel[1][1];
    Jac(5, 1) = vel_avel[1][2];
    Jac(0, 2) = avel_vel[0];
    Jac(1, 2) = avel_vel[1];
    Jac(2, 2) = avel_vel[2];
    Jac(3, 2) = avel_avel[0];
    Jac(4, 2) = avel_avel[1];
    Jac(5, 2) = avel_avel[2];
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
    A& RK_SERIAL_SAVE_WITH_NAME(vel_vel) & RK_SERIAL_SAVE_WITH_NAME(vel_avel) &
        RK_SERIAL_SAVE_WITH_NAME(avel_vel) &
        RK_SERIAL_SAVE_WITH_NAME(avel_avel) &
        RK_SERIAL_SAVE_WITH_NAME(vel_acc) & RK_SERIAL_SAVE_WITH_NAME(vel_aacc) &
        RK_SERIAL_SAVE_WITH_NAME(avel_acc) &
        RK_SERIAL_SAVE_WITH_NAME(avel_aacc);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    std::shared_ptr<frame_3D<value_type>> tmp;
    A& RK_SERIAL_LOAD_WITH_ALIAS("Parent", tmp) &
        RK_SERIAL_LOAD_WITH_NAME(vel_vel) & RK_SERIAL_LOAD_WITH_NAME(vel_avel) &
        RK_SERIAL_LOAD_WITH_NAME(avel_vel) &
        RK_SERIAL_LOAD_WITH_NAME(avel_avel) &
        RK_SERIAL_LOAD_WITH_NAME(vel_acc) & RK_SERIAL_LOAD_WITH_NAME(vel_aacc) &
        RK_SERIAL_LOAD_WITH_NAME(avel_acc) &
        RK_SERIAL_LOAD_WITH_NAME(avel_aacc);
    Parent = tmp;
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x00000027, 1, "jacobian_2D_3D",
                              shared_object)
};

/**
 * This class template represents the jacobians between a 3D frame and a generalized coordinate.
 * \tparam T The value-type.
 */
template <typename T>
class jacobian_3D_gen : public shared_object {
 public:
  using value_type = T;
  using self = jacobian_3D_gen<T>;

  /// Holds how much velocity is generated at output from the input velocity.
  vect<value_type, 3> vel_qd;
  /// Holds how much velocity is generated at output from the input angular velocity.
  vect<value_type, 3> avel_qd;
  /// Holds how much acceleration is generated at output from the input velocity.
  vect<value_type, 3> vel_qdd;
  /// Holds how much acceleration is generated at output from the input angular velocity.
  vect<value_type, 3> avel_qdd;

  /**
   * Parametrized constructor.
   */
  explicit jacobian_3D_gen(
      const vect<value_type, 3>& aVelQd = (vect<value_type, 3>()),
      const vect<value_type, 3>& aAVelQd = (vect<value_type, 3>()),
      const vect<value_type, 3>& aVelQdd = (vect<value_type, 3>()),
      const vect<value_type, 3>& aAVelQdd = (vect<value_type, 3>()))
      : vel_qd(aVelQd),
        avel_qd(aAVelQd),
        vel_qdd(aVelQdd),
        avel_qdd(aAVelQdd) {}

  self get_jac_relative_to(
      const std::shared_ptr<gen_coord<value_type>>& /*unused*/) const {
    return *this;
  }

  template <typename Matrix1, typename Matrix2>
  void write_to_matrices(Matrix1& Jac, Matrix2& JacDot) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    static_assert(is_fully_writable_matrix_v<Matrix2>);
    if ((Jac.get_row_count() != 1) || (Jac.get_col_count() != 6) ||
        (JacDot.get_row_count() != 1) || (JacDot.get_col_count() != 6)) {
      throw std::range_error(
          "Jacobian and JacobianDot matrices have the wrong size!");
    }

    Jac(0, 0) = vel_qd[0];
    Jac(0, 1) = vel_qd[1];
    Jac(0, 2) = vel_qd[2];
    Jac(0, 3) = avel_qd[0];
    Jac(0, 4) = avel_qd[1];
    Jac(0, 5) = avel_qd[2];

    JacDot(0, 0) = vel_qdd[0];
    JacDot(0, 1) = vel_qdd[1];
    JacDot(0, 2) = vel_qdd[2];
    JacDot(0, 3) = avel_qdd[0];
    JacDot(0, 4) = avel_qdd[1];
    JacDot(0, 5) = avel_qdd[2];
  }

  template <typename Matrix1>
  void write_to_matrices(Matrix1& Jac) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    if ((Jac.get_row_count() != 1) || (Jac.get_col_count() != 6)) {
      throw std::range_error("Jacobian matrix has the wrong size!");
    }

    Jac(0, 0) = vel_qd[0];
    Jac(0, 1) = vel_qd[1];
    Jac(0, 2) = vel_qd[2];
    Jac(0, 3) = avel_qd[0];
    Jac(0, 4) = avel_qd[1];
    Jac(0, 5) = avel_qd[2];
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(vel_qd) & RK_SERIAL_SAVE_WITH_NAME(avel_qd) &
        RK_SERIAL_SAVE_WITH_NAME(vel_qdd) & RK_SERIAL_SAVE_WITH_NAME(avel_qdd);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(vel_qd) & RK_SERIAL_LOAD_WITH_NAME(avel_qd) &
        RK_SERIAL_LOAD_WITH_NAME(vel_qdd) & RK_SERIAL_LOAD_WITH_NAME(avel_qdd);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x00000028, 1, "jacobian_3D_gen",
                              shared_object)
};

/**
 * This class template represents the jacobians between a 3D frame and a 2D frame.
 * \tparam T The value-type.
 */
template <typename T>
class jacobian_3D_2D : public shared_object {
 public:
  using value_type = T;
  using self = jacobian_3D_2D<T>;

  /// Holds the frame to which the jacobians are relative to.
  std::weak_ptr<frame_2D<value_type>> Parent;
  /// Holds how much velocity is generated at output from the input velocity.
  vect<vect<value_type, 2>, 3> vel_vel;
  /// Holds how much angular velocity is generated at output from the input velocity.
  vect<value_type, 3> vel_avel;
  /// Holds how much velocity is generated at output from the input angular velocity.
  vect<vect<value_type, 2>, 3> avel_vel;
  /// Holds how much angular velocity is generated at output from the input angular velocity.
  vect<value_type, 3> avel_avel;
  /// Holds how much acceleration is generated at output from the input velocity.
  vect<vect<value_type, 2>, 3> vel_acc;
  /// Holds how much angular acceleration is generated at output from the input velocity.
  vect<value_type, 3> vel_aacc;
  /// Holds how much acceleration is generated at output from the input angular velocity.
  vect<vect<value_type, 2>, 3> avel_acc;
  /// Holds how much angular acceleration is generated at output from the input angular velocity.
  vect<value_type, 3> avel_aacc;

  /**
   * Parametrized constructor.
   */
  explicit jacobian_3D_2D(
      std::weak_ptr<frame_2D<value_type>> aParent =
          std::weak_ptr<frame_2D<value_type>>(),
      const vect<vect<value_type, 2>, 3>& aVelVel =
          (vect<vect<value_type, 2>, 3>()),
      const vect<value_type, 3>& aVelAVel = (vect<value_type, 3>()),
      const vect<vect<value_type, 2>, 3>& aAVelVel =
          (vect<vect<value_type, 2>, 3>()),
      const vect<value_type, 3>& aAVelAVel = (vect<value_type, 3>()),
      const vect<vect<value_type, 2>, 3>& aVelAcc =
          (vect<vect<value_type, 2>, 3>()),
      const vect<value_type, 3>& aVelAAcc = (vect<value_type, 3>()),
      const vect<vect<value_type, 2>, 3>& aAVelAcc =
          (vect<vect<value_type, 2>, 3>()),
      const vect<value_type, 3>& aAVelAAcc = (vect<value_type, 3>()))
      : Parent(std::move(aParent)),
        vel_vel(aVelVel),
        vel_avel(aVelAVel),
        avel_vel(aAVelVel),
        avel_avel(aAVelAVel),
        vel_acc(aVelAcc),
        vel_aacc(aVelAAcc),
        avel_acc(aAVelAcc),
        avel_aacc(aAVelAAcc) {}

  self get_jac_relative_to(
      const std::shared_ptr<frame_2D<value_type>>& aFrame) const {

    frame_2D<value_type> f2 = aFrame->getFrameRelativeTo(Parent.lock());

    vect<vect<value_type, 2>, 3> new_vel_vel(
        (vel_avel[0] % f2.Position + vel_vel[0]) * f2.Rotation,
        (vel_avel[1] % f2.Position + vel_vel[1]) * f2.Rotation,
        (vel_avel[2] % f2.Position + vel_vel[2]) * f2.Rotation);
    vect<vect<value_type, 2>, 3> new_vel_acc(
        (vel_avel[0] % f2.Velocity + vel_aacc[0] % f2.Position + vel_acc[0]) *
                f2.Rotation -
            f2.AngVelocity % new_vel_vel[0],
        (vel_avel[1] % f2.Velocity + vel_aacc[1] % f2.Position + vel_acc[1]) *
                f2.Rotation -
            f2.AngVelocity % new_vel_vel[1],
        (vel_avel[2] % f2.Velocity + vel_aacc[2] % f2.Position + vel_acc[2]) *
                f2.Rotation -
            f2.AngVelocity % new_vel_vel[2]);
    vect<vect<value_type, 2>, 3> new_avel_vel(
        (avel_avel[0] % f2.Position + avel_vel[0]) * f2.Rotation,
        (avel_avel[1] % f2.Position + avel_vel[1]) * f2.Rotation,
        (avel_avel[2] % f2.Position + avel_vel[2]) * f2.Rotation);
    vect<vect<value_type, 2>, 3> new_avel_acc(
        (avel_avel[0] % f2.Velocity + avel_aacc[0] % f2.Position +
         avel_acc[0]) *
                f2.Rotation -
            f2.AngVelocity % new_avel_vel[0],
        (avel_avel[1] % f2.Velocity + avel_aacc[1] % f2.Position +
         avel_acc[1]) *
                f2.Rotation -
            f2.AngVelocity % new_avel_vel[1],
        (avel_avel[2] % f2.Velocity + avel_aacc[2] % f2.Position +
         avel_acc[2]) *
                f2.Rotation -
            f2.AngVelocity % new_avel_vel[2]);

    return self(aFrame,
                new_vel_vel,   // VelVel
                vel_avel,      // VelAVel
                new_avel_vel,  // AVelVel
                avel_avel,     // AVelAVel
                new_vel_acc,   // VelAcc
                vel_aacc,      // VelAAcc
                new_avel_acc,  // AVelAcc
                avel_aacc      // AVelAAcc
    );
  }

  template <typename Matrix1, typename Matrix2>
  void write_to_matrices(Matrix1& Jac, Matrix2& JacDot) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    static_assert(is_fully_writable_matrix_v<Matrix2>);
    if ((Jac.get_row_count() != 3) || (Jac.get_col_count() != 6) ||
        (JacDot.get_row_count() != 3) || (JacDot.get_col_count() != 6)) {
      throw std::range_error(
          "Jacobian and JacobianDot matrices have the wrong size!");
    }

    Jac(0, 0) = vel_vel[0][0];
    Jac(1, 0) = vel_vel[0][1];
    Jac(2, 0) = vel_avel[0];
    Jac(0, 1) = vel_vel[1][0];
    Jac(1, 1) = vel_vel[1][1];
    Jac(2, 1) = vel_avel[1];
    Jac(0, 2) = vel_vel[2][0];
    Jac(1, 2) = vel_vel[2][1];
    Jac(2, 2) = vel_avel[2];
    Jac(0, 3) = avel_vel[0][0];
    Jac(1, 3) = avel_vel[0][1];
    Jac(2, 3) = avel_avel[0];
    Jac(0, 4) = avel_vel[1][0];
    Jac(1, 4) = avel_vel[1][1];
    Jac(2, 4) = avel_avel[1];
    Jac(0, 5) = avel_vel[2][0];
    Jac(1, 5) = avel_vel[2][1];
    Jac(2, 5) = avel_avel[2];

    JacDot(0, 0) = vel_acc[0][0];
    JacDot(1, 0) = vel_acc[0][1];
    JacDot(2, 0) = vel_aacc[0];
    JacDot(0, 1) = vel_acc[1][0];
    JacDot(1, 1) = vel_acc[1][1];
    JacDot(2, 1) = vel_aacc[1];
    JacDot(0, 2) = vel_acc[2][0];
    JacDot(1, 2) = vel_acc[2][1];
    JacDot(2, 2) = vel_aacc[2];
    JacDot(0, 3) = avel_acc[0][0];
    JacDot(1, 3) = avel_acc[0][1];
    JacDot(2, 3) = avel_aacc[0];
    JacDot(0, 4) = avel_acc[1][0];
    JacDot(1, 4) = avel_acc[1][1];
    JacDot(2, 4) = avel_aacc[1];
    JacDot(0, 5) = avel_acc[2][0];
    JacDot(1, 5) = avel_acc[2][1];
    JacDot(2, 5) = avel_aacc[2];
  }

  template <typename Matrix1>
  void write_to_matrices(Matrix1& Jac) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    if ((Jac.get_row_count() != 3) || (Jac.get_col_count() != 6)) {
      throw std::range_error("Jacobian matrix has the wrong size!");
    }

    Jac(0, 0) = vel_vel[0][0];
    Jac(1, 0) = vel_vel[0][1];
    Jac(2, 0) = vel_avel[0];
    Jac(0, 1) = vel_vel[1][0];
    Jac(1, 1) = vel_vel[1][1];
    Jac(2, 1) = vel_avel[1];
    Jac(0, 2) = vel_vel[2][0];
    Jac(1, 2) = vel_vel[2][1];
    Jac(2, 2) = vel_avel[2];
    Jac(0, 3) = avel_vel[0][0];
    Jac(1, 3) = avel_vel[0][1];
    Jac(2, 3) = avel_avel[0];
    Jac(0, 4) = avel_vel[1][0];
    Jac(1, 4) = avel_vel[1][1];
    Jac(2, 4) = avel_avel[1];
    Jac(0, 5) = avel_vel[2][0];
    Jac(1, 5) = avel_vel[2][1];
    Jac(2, 5) = avel_avel[2];
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
    A& RK_SERIAL_SAVE_WITH_NAME(vel_vel) & RK_SERIAL_SAVE_WITH_NAME(vel_avel) &
        RK_SERIAL_SAVE_WITH_NAME(avel_vel) &
        RK_SERIAL_SAVE_WITH_NAME(avel_avel) &
        RK_SERIAL_SAVE_WITH_NAME(vel_acc) & RK_SERIAL_SAVE_WITH_NAME(vel_aacc) &
        RK_SERIAL_SAVE_WITH_NAME(avel_acc) &
        RK_SERIAL_SAVE_WITH_NAME(avel_aacc);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    std::shared_ptr<frame_2D<value_type>> tmp;
    A& RK_SERIAL_LOAD_WITH_ALIAS("Parent", tmp) &
        RK_SERIAL_LOAD_WITH_NAME(vel_vel) & RK_SERIAL_LOAD_WITH_NAME(vel_avel) &
        RK_SERIAL_LOAD_WITH_NAME(avel_vel) &
        RK_SERIAL_LOAD_WITH_NAME(avel_avel) &
        RK_SERIAL_LOAD_WITH_NAME(vel_acc) & RK_SERIAL_LOAD_WITH_NAME(vel_aacc) &
        RK_SERIAL_LOAD_WITH_NAME(avel_acc) &
        RK_SERIAL_LOAD_WITH_NAME(avel_aacc);
    Parent = tmp;
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x00000029, 1, "jacobian_3D_2D",
                              shared_object)
};

/**
 * This class template represents the jacobians between a 3D frame and a 3D frame.
 * \tparam T The value-type.
 */
template <typename T>
class jacobian_3D_3D : public shared_object {
 public:
  using value_type = T;
  using self = jacobian_3D_3D<T>;

  /// Holds the frame to which the jacobians are relative to.
  std::weak_ptr<frame_3D<value_type>> Parent;
  /// Holds how much velocity is generated at output from the input velocity.
  vect<vect<value_type, 3>, 3> vel_vel;
  /// Holds how much angular velocity is generated at output from the input velocity.
  vect<vect<value_type, 3>, 3> vel_avel;
  /// Holds how much velocity is generated at output from the input angular velocity.
  vect<vect<value_type, 3>, 3> avel_vel;
  /// Holds how much angular velocity is generated at output from the input angular velocity.
  vect<vect<value_type, 3>, 3> avel_avel;
  /// Holds how much acceleration is generated at output from the input velocity.
  vect<vect<value_type, 3>, 3> vel_acc;
  /// Holds how much angular acceleration is generated at output from the input velocity.
  vect<vect<value_type, 3>, 3> vel_aacc;
  /// Holds how much acceleration is generated at output from the input angular velocity.
  vect<vect<value_type, 3>, 3> avel_acc;
  /// Holds how much angular acceleration is generated at output from the input angular velocity.
  vect<vect<value_type, 3>, 3> avel_aacc;

  /**
   * Parametrized constructor.
   */
  explicit jacobian_3D_3D(std::weak_ptr<frame_3D<value_type>> aParent =
                              std::weak_ptr<frame_3D<value_type>>(),
                          const vect<vect<value_type, 3>, 3>& aVelVel =
                              (vect<vect<value_type, 3>, 3>()),
                          const vect<vect<value_type, 3>, 3>& aVelAVel =
                              (vect<vect<value_type, 3>, 3>()),
                          const vect<vect<value_type, 3>, 3>& aAVelVel =
                              (vect<vect<value_type, 3>, 3>()),
                          const vect<vect<value_type, 3>, 3>& aAVelAVel =
                              (vect<vect<value_type, 3>, 3>()),
                          const vect<vect<value_type, 3>, 3>& aVelAcc =
                              (vect<vect<value_type, 3>, 3>()),
                          const vect<vect<value_type, 3>, 3>& aVelAAcc =
                              (vect<vect<value_type, 3>, 3>()),
                          const vect<vect<value_type, 3>, 3>& aAVelAcc =
                              (vect<vect<value_type, 3>, 3>()),
                          const vect<vect<value_type, 3>, 3>& aAVelAAcc =
                              (vect<vect<value_type, 3>, 3>()))
      : Parent(std::move(aParent)),
        vel_vel(aVelVel),
        vel_avel(aVelAVel),
        avel_vel(aAVelVel),
        avel_avel(aAVelAVel),
        vel_acc(aVelAcc),
        vel_aacc(aVelAAcc),
        avel_acc(aAVelAcc),
        avel_aacc(aAVelAAcc) {}

  self get_jac_relative_to(
      const std::shared_ptr<frame_3D<value_type>>& aFrame) const {

    frame_3D<value_type> f2 = aFrame->getFrameRelativeTo(Parent.lock());
    rot_mat_3D<value_type> R(f2.Quat.getRotMat());

    vect<vect<value_type, 3>, 3> new_vel_avel(vel_avel[0] * R, vel_avel[1] * R,
                                              vel_avel[2] * R);
    vect<vect<value_type, 3>, 3> new_vel_aacc(
        vel_aacc[0] * R - f2.AngVelocity % new_vel_avel[0],
        vel_aacc[1] * R - f2.AngVelocity % new_vel_avel[1],
        vel_aacc[2] * R - f2.AngVelocity % new_vel_avel[2]);

    vect<vect<value_type, 3>, 3> new_avel_avel(
        avel_avel[0] * R, avel_avel[1] * R, avel_avel[2] * R);
    vect<vect<value_type, 3>, 3> new_avel_aacc(
        avel_aacc[0] * R - f2.AngVelocity % new_avel_avel[0],
        avel_aacc[1] * R - f2.AngVelocity % new_avel_avel[1],
        avel_aacc[2] * R - f2.AngVelocity % new_avel_avel[2]);

    vect<vect<value_type, 3>, 3> new_vel_vel(
        (vel_avel[0] % f2.Position + vel_vel[0]) * R,
        (vel_avel[1] % f2.Position + vel_vel[1]) * R,
        (vel_avel[2] % f2.Position + vel_vel[2]) * R);
    vect<vect<value_type, 3>, 3> new_vel_acc(
        (vel_avel[0] % f2.Velocity + vel_aacc[0] % f2.Position + vel_acc[0]) *
                R -
            f2.AngVelocity % new_vel_vel[0],
        (vel_avel[1] % f2.Velocity + vel_aacc[1] % f2.Position + vel_acc[1]) *
                R -
            f2.AngVelocity % new_vel_vel[1],
        (vel_avel[2] % f2.Velocity + vel_aacc[2] % f2.Position + vel_acc[2]) *
                R -
            f2.AngVelocity % new_vel_vel[2]);

    vect<vect<value_type, 3>, 3> new_avel_vel(
        (avel_avel[0] % f2.Position + avel_vel[0]) * R,
        (avel_avel[1] % f2.Position + avel_vel[1]) * R,
        (avel_avel[2] % f2.Position + avel_vel[2]) * R);
    vect<vect<value_type, 3>, 3> new_avel_acc(
        (avel_avel[0] % f2.Velocity + avel_aacc[0] % f2.Position +
         avel_acc[0]) *
                R -
            f2.AngVelocity % new_avel_vel[0],
        (avel_avel[1] % f2.Velocity + avel_aacc[1] % f2.Position +
         avel_acc[1]) *
                R -
            f2.AngVelocity % new_avel_vel[1],
        (avel_avel[2] % f2.Velocity + avel_aacc[2] % f2.Position +
         avel_acc[2]) *
                R -
            f2.AngVelocity % new_avel_vel[2]);

    return self(aFrame, new_vel_vel, new_vel_avel, new_avel_vel, new_avel_avel,
                new_vel_acc, new_vel_aacc, new_avel_acc, new_avel_aacc);
  }

  template <typename Matrix1, typename Matrix2>
  void write_to_matrices(Matrix1& Jac, Matrix2& JacDot) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    static_assert(is_fully_writable_matrix_v<Matrix2>);
    if ((Jac.get_row_count() != 6) || (Jac.get_col_count() != 6) ||
        (JacDot.get_row_count() != 6) || (JacDot.get_col_count() != 6)) {
      throw std::range_error(
          "Jacobian and JacobianDot matrices have the wrong size!");
    }

    Jac(0, 0) = vel_vel[0][0];
    Jac(1, 0) = vel_vel[0][1];
    Jac(2, 0) = vel_vel[0][2];
    Jac(3, 0) = vel_avel[0][0];
    Jac(4, 0) = vel_avel[0][1];
    Jac(5, 0) = vel_avel[0][2];
    Jac(0, 1) = vel_vel[1][0];
    Jac(1, 1) = vel_vel[1][1];
    Jac(2, 1) = vel_vel[1][2];
    Jac(3, 1) = vel_avel[1][0];
    Jac(4, 1) = vel_avel[1][1];
    Jac(5, 1) = vel_avel[1][2];
    Jac(0, 2) = vel_vel[2][0];
    Jac(1, 2) = vel_vel[2][1];
    Jac(2, 2) = vel_vel[2][2];
    Jac(3, 2) = vel_avel[2][0];
    Jac(4, 2) = vel_avel[2][1];
    Jac(5, 2) = vel_avel[2][2];
    Jac(0, 3) = avel_vel[0][0];
    Jac(1, 3) = avel_vel[0][1];
    Jac(2, 3) = avel_vel[0][2];
    Jac(3, 3) = avel_avel[0][0];
    Jac(4, 3) = avel_avel[0][1];
    Jac(5, 3) = avel_avel[0][2];
    Jac(0, 4) = avel_vel[1][0];
    Jac(1, 4) = avel_vel[1][1];
    Jac(2, 4) = avel_vel[1][2];
    Jac(3, 4) = avel_avel[1][0];
    Jac(4, 4) = avel_avel[1][1];
    Jac(5, 4) = avel_avel[1][2];
    Jac(0, 5) = avel_vel[2][0];
    Jac(1, 5) = avel_vel[2][1];
    Jac(2, 5) = avel_vel[2][2];
    Jac(3, 5) = avel_avel[2][0];
    Jac(4, 5) = avel_avel[2][1];
    Jac(5, 5) = avel_avel[2][2];

    JacDot(0, 0) = vel_acc[0][0];
    JacDot(1, 0) = vel_acc[0][1];
    JacDot(2, 0) = vel_acc[0][2];
    JacDot(3, 0) = vel_aacc[0][0];
    JacDot(4, 0) = vel_aacc[0][1];
    JacDot(5, 0) = vel_aacc[0][2];
    JacDot(0, 1) = vel_acc[1][0];
    JacDot(1, 1) = vel_acc[1][1];
    JacDot(2, 1) = vel_acc[1][2];
    JacDot(3, 1) = vel_aacc[1][0];
    JacDot(4, 1) = vel_aacc[1][1];
    JacDot(5, 1) = vel_aacc[1][2];
    JacDot(0, 2) = vel_acc[2][0];
    JacDot(1, 2) = vel_acc[2][1];
    JacDot(2, 2) = vel_acc[2][2];
    JacDot(3, 2) = vel_aacc[2][0];
    JacDot(4, 2) = vel_aacc[2][1];
    JacDot(5, 2) = vel_aacc[2][2];
    JacDot(0, 3) = avel_acc[0][0];
    JacDot(1, 3) = avel_acc[0][1];
    JacDot(2, 3) = avel_acc[0][2];
    JacDot(3, 3) = avel_aacc[0][0];
    JacDot(4, 3) = avel_aacc[0][1];
    JacDot(5, 3) = avel_aacc[0][2];
    JacDot(0, 4) = avel_acc[1][0];
    JacDot(1, 4) = avel_acc[1][1];
    JacDot(2, 4) = avel_acc[1][2];
    JacDot(3, 4) = avel_aacc[1][0];
    JacDot(4, 4) = avel_aacc[1][1];
    JacDot(5, 4) = avel_aacc[1][2];
    JacDot(0, 5) = avel_acc[2][0];
    JacDot(1, 5) = avel_acc[2][1];
    JacDot(2, 5) = avel_acc[2][2];
    JacDot(3, 5) = avel_aacc[2][0];
    JacDot(4, 5) = avel_aacc[2][1];
    JacDot(5, 5) = avel_aacc[2][2];
  }

  template <typename Matrix1>
  void write_to_matrices(Matrix1& Jac) const {
    static_assert(is_fully_writable_matrix_v<Matrix1>);
    if ((Jac.get_row_count() != 6) || (Jac.get_col_count() != 6)) {
      throw std::range_error("Jacobian matrix has the wrong size!");
    }

    Jac(0, 0) = vel_vel[0][0];
    Jac(1, 0) = vel_vel[0][1];
    Jac(2, 0) = vel_vel[0][2];
    Jac(3, 0) = vel_avel[0][0];
    Jac(4, 0) = vel_avel[0][1];
    Jac(5, 0) = vel_avel[0][2];
    Jac(0, 1) = vel_vel[1][0];
    Jac(1, 1) = vel_vel[1][1];
    Jac(2, 1) = vel_vel[1][2];
    Jac(3, 1) = vel_avel[1][0];
    Jac(4, 1) = vel_avel[1][1];
    Jac(5, 1) = vel_avel[1][2];
    Jac(0, 2) = vel_vel[2][0];
    Jac(1, 2) = vel_vel[2][1];
    Jac(2, 2) = vel_vel[2][2];
    Jac(3, 2) = vel_avel[2][0];
    Jac(4, 2) = vel_avel[2][1];
    Jac(5, 2) = vel_avel[2][2];
    Jac(0, 3) = avel_vel[0][0];
    Jac(1, 3) = avel_vel[0][1];
    Jac(2, 3) = avel_vel[0][2];
    Jac(3, 3) = avel_avel[0][0];
    Jac(4, 3) = avel_avel[0][1];
    Jac(5, 3) = avel_avel[0][2];
    Jac(0, 4) = avel_vel[1][0];
    Jac(1, 4) = avel_vel[1][1];
    Jac(2, 4) = avel_vel[1][2];
    Jac(3, 4) = avel_avel[1][0];
    Jac(4, 4) = avel_avel[1][1];
    Jac(5, 4) = avel_avel[1][2];
    Jac(0, 5) = avel_vel[2][0];
    Jac(1, 5) = avel_vel[2][1];
    Jac(2, 5) = avel_vel[2][2];
    Jac(3, 5) = avel_avel[2][0];
    Jac(4, 5) = avel_avel[2][1];
    Jac(5, 5) = avel_avel[2][2];
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
    A& RK_SERIAL_SAVE_WITH_NAME(vel_vel) & RK_SERIAL_SAVE_WITH_NAME(vel_avel) &
        RK_SERIAL_SAVE_WITH_NAME(avel_vel) &
        RK_SERIAL_SAVE_WITH_NAME(avel_avel) &
        RK_SERIAL_SAVE_WITH_NAME(vel_acc) & RK_SERIAL_SAVE_WITH_NAME(vel_aacc) &
        RK_SERIAL_SAVE_WITH_NAME(avel_acc) &
        RK_SERIAL_SAVE_WITH_NAME(avel_aacc);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    std::shared_ptr<frame_3D<value_type>> tmp;
    A& RK_SERIAL_LOAD_WITH_ALIAS("Parent", tmp) &
        RK_SERIAL_LOAD_WITH_NAME(vel_vel) & RK_SERIAL_LOAD_WITH_NAME(vel_avel) &
        RK_SERIAL_LOAD_WITH_NAME(avel_vel) &
        RK_SERIAL_LOAD_WITH_NAME(avel_avel) &
        RK_SERIAL_LOAD_WITH_NAME(vel_acc) & RK_SERIAL_LOAD_WITH_NAME(vel_aacc) &
        RK_SERIAL_LOAD_WITH_NAME(avel_acc) &
        RK_SERIAL_LOAD_WITH_NAME(avel_aacc);
    Parent = tmp;
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x0000002A, 1, "jacobian_3D_3D",
                              shared_object)
};

}  // namespace ReaK

#endif

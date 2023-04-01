/**
 * \file rotations_3D.h
 *
 * This library declares all geometric 3D rotation classes for fixed (2,3) and variable dimensions.
 *
 * Note: All matrix memory is organized, by default, such that columns are concatenated. This
 *       was found to be a more efficient representation since columns often have
 *       more significances than rows (representing basis vectors for example).
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

#ifndef REAK_MATH_KINETOSTATICS_ROTATIONS_3D_H_
#define REAK_MATH_KINETOSTATICS_ROTATIONS_3D_H_

#include "ReaK/math/lin_alg/mat_alg_skew_symmetric.h"
#include "ReaK/math/lin_alg/mat_alg_square.h"
#include "ReaK/math/lin_alg/mat_alg_symmetric.h"
#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/vect_alg.h"

#include <cassert>
#include <type_traits>

namespace ReaK {

// Forward declaration
template <class T>
class quaternion;

template <class T>
class unit_quat;

template <class T>
class euler_angles_TB;

template <class T>
class axis_angle;

template <class T>
class trans_mat_3D;

/// This class is a rotation matrix 3 by 3.
template <typename T>
class rot_mat_3D {
 public:
  using self = rot_mat_3D<T>;

  using value_type = T;
  using container_type = void;

  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;

  using col_iterator = void;
  using const_col_iterator = void;
  using row_iterator = void;
  using const_row_iterator = void;

  using size_type = int;
  using difference_type = int;

  static constexpr unsigned int static_row_count = 3;
  static constexpr unsigned int static_col_count = 3;
  static constexpr mat_alignment::tag alignment = mat_alignment::column_major;
  static constexpr mat_structure::tag structure = mat_structure::orthogonal;

 private:
  std::array<value_type, 9> q = {};

  /// Constructor from the components of the rotation matrix.
  rot_mat_3D(const_reference a11, const_reference a12, const_reference a13,
             const_reference a21, const_reference a22, const_reference a23,
             const_reference a31, const_reference a32,
             const_reference a33) noexcept {
    q[0] = a11;
    q[1] = a21;
    q[2] = a31;
    q[3] = a12;
    q[4] = a22;
    q[5] = a32;
    q[6] = a13;
    q[7] = a23;
    q[8] = a33;
  }

 public:
  friend class quaternion<value_type>;
  friend class euler_angles_TB<value_type>;
  friend class axis_angle<value_type>;
  friend class trans_mat_3D<value_type>;

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /// Default constructor, sets the matrix to identity (no rotation).
  rot_mat_3D() noexcept {
    q[0] = 1.0;
    q[1] = 0.0;
    q[2] = 0.0;
    q[3] = 0.0;
    q[4] = 1.0;
    q[5] = 0.0;
    q[6] = 0.0;
    q[7] = 0.0;
    q[8] = 1.0;
  }

  /// Constructor from an array of components.
  explicit rot_mat_3D(const_pointer M) noexcept {
    vect<value_type, 3> v1 = unit(vect<value_type, 3>(M));
    q[0] = v1[0];
    q[1] = v1[1];
    q[2] = v1[2];
    vect<value_type, 3> v2(&M[3]);
    v2 = unit(v2 - (v2 * v1) * v1);
    q[3] = v2[0];
    q[4] = v2[1];
    q[5] = v2[2];
    v2 = v1 % v2;
    q[6] = v2[0];
    q[7] = v2[1];
    q[8] = v2[2];
  }
  explicit rot_mat_3D(pointer M) : rot_mat_3D(static_cast<const_pointer>(M)) {}

  rot_mat_3D(const self&) noexcept = default;

  template <typename Matrix>
  explicit rot_mat_3D(const Matrix& M) {
    static_assert(is_readable_matrix_v<Matrix>);
    if ((M.get_col_count() != 3) || (M.get_row_count() != 3)) {
      throw std::range_error(
          "Right-hand-side of assignment to a 3D rotation matrix is not of "
          "dimension 3x3!");
    }
    vect<value_type, 3> v1(M(0, 0), M(1, 0), M(2, 0));
    q[0] = v1[0];
    q[1] = v1[1];
    q[2] = v1[2];
    vect<value_type, 3> v2(M(0, 1), M(1, 1), M(2, 1));
    v2 = unit(v2 - (v2 * v1) * v1);
    q[3] = v2[0];
    q[4] = v2[1];
    q[5] = v2[2];
    v2 = v1 % v2;
    q[6] = v2[0];
    q[7] = v2[1];
    q[8] = v2[2];
  }

  explicit rot_mat_3D(const quaternion<value_type>& Q) noexcept
      : rot_mat_3D(Q.getRotMat()) {}
  explicit rot_mat_3D(const euler_angles_TB<value_type>& E) noexcept
      : rot_mat_3D(E.getRotMat()) {}
  explicit rot_mat_3D(const axis_angle<value_type>& A) noexcept
      : rot_mat_3D(A.getRotMat()) {}

  // Copy constructor. Default is good. \test PASSED

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Provides a copy of the rotation matrix as an ordinary 3x3 matrix.
  mat<value_type, mat_structure::square> getMat() const {
    return mat<value_type, mat_structure::square>(q[0], q[3], q[6], q[1], q[4],
                                                  q[7], q[2], q[5], q[8]);
  }

  quaternion<value_type> getQuaternion() const noexcept {
    using std::sqrt;
    std::array<value_type, 4> a;
    value_type tra = q[0] + q[4] + q[8];
    if (tra > 0.01) {
      a[0] = value_type(0.5) * sqrt(value_type(1.0) + tra);
      a[1] = value_type(0.25) * (q[5] - q[7]) / a[0];
      a[2] = value_type(0.25) * (q[6] - q[2]) / a[0];
      a[3] = value_type(0.25) * (q[1] - q[3]) / a[0];
    } else if ((q[0] > q[4]) && (q[0] > q[8])) {
      a[1] = value_type(0.5) * sqrt(value_type(1.0) + q[0] - q[4] - q[8]);
      a[0] = value_type(0.25) * (q[7] - q[5]) / a[1];
      a[2] = value_type(0.25) * (q[3] + q[1]) / a[1];
      a[3] = value_type(0.25) * (q[6] + q[2]) / a[1];
    } else if (q[4] > q[8]) {
      a[2] = value_type(0.5) * sqrt(value_type(1.0) + q[4] - q[0] - q[8]);
      a[0] = value_type(0.25) * (q[6] - q[2]) / a[2];
      a[1] = value_type(0.25) * (q[3] + q[1]) / a[2];
      a[3] = value_type(0.25) * (q[7] + q[5]) / a[2];
    } else {
      a[3] = value_type(0.5) * sqrt(value_type(1.0) + q[8] - q[0] - q[4]);
      a[0] = value_type(0.25) * (q[3] - q[1]) / a[3];
      a[1] = value_type(0.25) * (q[6] + q[2]) / a[3];
      a[2] = value_type(0.25) * (q[7] + q[5]) / a[3];
    }
    return quaternion<value_type>(a[0], a[1], a[2], a[3]);
  }

  euler_angles_TB<value_type> getEulerAnglesTB() const noexcept {
    using std::asin;
    using std::atan2;
    using std::cos;
    value_type yaw{};
    value_type pitch{};
    value_type roll{};
    if ((q[2] != value_type(1.0)) && (q[2] != value_type(-1.0))) {
      pitch = asin(-q[2]);
      value_type cp = value_type(1.0) / cos(pitch);
      roll = atan2(cp * q[5], cp * q[8]);
      yaw = atan2(cp * q[1], cp * q[0]);
    } else {
      yaw = value_type(0.0);
      roll = atan2(-q[2] * q[3], -q[2] * q[6]);
      pitch = -q[2] * value_type(1.57079632679489662);
    }
    return euler_angles_TB<value_type>(yaw, pitch, roll);
  }

  axis_angle<value_type> getAxisAngle() const noexcept {
    using std::acos;
    using std::sin;
    value_type tmp(value_type(0.5) * (trace(*this) - value_type(1.0)));
    vect<value_type, 3> axis{};
    value_type angle{};
    if (tmp > value_type(0.0000001)) {
      angle = acos(tmp);
      value_type cosec_a = value_type(0.5) / sin(angle);
      axis[0] = (q[5] - q[7]) * cosec_a;
      axis[1] = (q[6] - q[2]) * cosec_a;
      axis[2] = (q[1] - q[3]) * cosec_a;
    } else {
      axis[0] = value_type(1.0);
      axis[1] = value_type(0.0);
      axis[2] = value_type(0.0);
      angle = value_type(0.0);
    }
    return axis_angle<value_type>(angle, axis);
  }

  /// Array indexing operator, accessor for read only.
  value_type operator[](int i) const noexcept {
    assert(i < 9);
    return q[i];
  }

  /// Array double-indexing operator, ith row and jth column, accessor for read only.
  value_type operator()(int i, int j) const noexcept {
    assert((i < 3) && (j < 3));
    return q[j * 3 + i];
  }

  int get_row_count() const noexcept { return 3; }
  int get_col_count() const noexcept { return 3; }

  /*******************************************************************************
                         Assignment Operators
*******************************************************************************/

  self& operator=(const self& M) noexcept = default;

  template <typename Other>
  self& operator=(const Other& M) {
    return *this = self{M};
  }

  /// Multiplication by a rotation matrix and store.
  self& operator*=(const self& M) noexcept { return *this = *this * M; }

  template <typename Other>
  self& operator*=(const Other& M) {
    return *this = *this * self{M};
  }

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /// Multiplication by a rotation matrix.
  friend self operator*(const self& M1, const self& M2) noexcept {
    return self(M1.q[0] * M2.q[0] + M1.q[3] * M2.q[1] + M1.q[6] * M2.q[2],
                M1.q[0] * M2.q[3] + M1.q[3] * M2.q[4] + M1.q[6] * M2.q[5],
                M1.q[0] * M2.q[6] + M1.q[3] * M2.q[7] + M1.q[6] * M2.q[8],
                M1.q[1] * M2.q[0] + M1.q[4] * M2.q[1] + M1.q[7] * M2.q[2],
                M1.q[1] * M2.q[3] + M1.q[4] * M2.q[4] + M1.q[7] * M2.q[5],
                M1.q[1] * M2.q[6] + M1.q[4] * M2.q[7] + M1.q[7] * M2.q[8],
                M1.q[2] * M2.q[0] + M1.q[5] * M2.q[1] + M1.q[8] * M2.q[2],
                M1.q[2] * M2.q[3] + M1.q[5] * M2.q[4] + M1.q[8] * M2.q[5],
                M1.q[2] * M2.q[6] + M1.q[5] * M2.q[7] + M1.q[8] * M2.q[8]);
  }

  template <typename Other>
  friend self operator*(const self& R, const Other& M) {
    return R * self{M};
  }

  template <typename Other>
  friend self operator*(const Other& M, const self& R) {
    return self{M} * R;
  }

  /// Multiplication with a column vector.
  friend vect<value_type, 3> operator*(const self& R,
                                       const vect<value_type, 3>& V) noexcept {
    return vect<value_type, 3>(R.q[0] * V[0] + R.q[3] * V[1] + R.q[6] * V[2],
                               R.q[1] * V[0] + R.q[4] * V[1] + R.q[7] * V[2],
                               R.q[2] * V[0] + R.q[5] * V[1] + R.q[8] * V[2]);
  }

  friend vect<value_type, 3> operator*(const vect<value_type, 3>& V,
                                       const self& R) noexcept {
    return vect<value_type, 3>(R.q[0] * V[0] + R.q[1] * V[1] + R.q[2] * V[2],
                               R.q[3] * V[0] + R.q[4] * V[1] + R.q[5] * V[2],
                               R.q[6] * V[0] + R.q[7] * V[1] + R.q[8] * V[2]);
  }

  /*******************************************************************************
                           Standard Matrix Methods
  *******************************************************************************/

  /// Produces a transpose matrix which is the inverse rotation.
  friend self transpose(const self& R) noexcept {
    return self(R.q[0], R.q[1], R.q[2], R.q[3], R.q[4], R.q[5], R.q[6], R.q[7],
                R.q[8]);
  }

  /// Produces a cofactor matrix which is the same as the rotation matrix itself.
  friend self cofactor(const self& R) noexcept { return R; }

  /// Invert the transformation.
  friend self invert(const self& R) noexcept {
    return self(R.q[0], R.q[1], R.q[2], R.q[3], R.q[4], R.q[5], R.q[6], R.q[7],
                R.q[8]);
  }

  /// Gets the trace of the matrix.
  friend value_type trace(const self& R) noexcept {
    return R.q[0] + R.q[4] + R.q[8];
  }

  /// Gets the determinant of the matrix.
  friend value_type determinant(const self& /*unused*/) noexcept {
    return value_type(1.0);
  }

  /// Gets the symmetric part of the matrix.
  mat<value_type, mat_structure::symmetric> getSymPart() const {
    return mat<value_type, mat_structure::symmetric>(*this);
  }

  /// Gets the skew-symmetric part of the matrix.
  mat<value_type, mat_structure::skew_symmetric> getSkewSymPart() const {
    return mat<value_type, mat_structure::skew_symmetric>(*this);
  }

  /// Loading a rot_mat_3D value with a name.
  friend serialization::iarchive& operator&(
      serialization::iarchive& in,
      const std::pair<std::string, rot_mat_3D<T>&>& R) {
    return in & RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_r11", R.second.q[0]) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_r21", R.second.q[1]) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_r31", R.second.q[2]) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_r12", R.second.q[3]) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_r22", R.second.q[4]) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_r32", R.second.q[5]) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_r13", R.second.q[6]) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_r23", R.second.q[7]) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_r33", R.second.q[8]);
  }

  /// Loading a rot_mat_3D value.
  friend serialization::iarchive& operator>>(serialization::iarchive& in,
                                             rot_mat_3D<T>& R) {
    return in & RK_SERIAL_LOAD_WITH_NAME(R);
  }

  /// Saving a rot_mat_3D value with a name.
  friend serialization::oarchive& operator&(
      serialization::oarchive& out,
      const std::pair<std::string, const rot_mat_3D<T>&>& R) {
    return out & RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_r11", R.second.q[0]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_r21", R.second.q[1]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_r31", R.second.q[2]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_r12", R.second.q[3]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_r22", R.second.q[4]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_r32", R.second.q[5]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_r13", R.second.q[6]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_r23", R.second.q[7]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_r33", R.second.q[8]);
  }

  /// Saving a rot_mat_3D value.
  friend serialization::oarchive& operator<<(serialization::oarchive& out,
                                             const rot_mat_3D<T>& R) {
    return out & RK_SERIAL_SAVE_WITH_NAME(R);
  }
};

namespace rtti {

template <typename T>
struct get_type_id<rot_mat_3D<T>> {
  static constexpr unsigned int ID = 0x00000018;
  static constexpr auto type_name = std::string_view{"ReaK::rot_mat_3D"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }

  using save_type = const rot_mat_3D<T>&;
  using load_type = rot_mat_3D<T>&;
};
}  // namespace rtti

/// Prints a rotation matrix to a standard output stream (<<) as
/// "((a11; a12; a13); (a21; a22; a23); (a31; a32; a33))".
template <class T>
std::ostream& operator<<(std::ostream& out_stream, const rot_mat_3D<T>& R) {
  return out_stream << "((" << R(0, 0) << "; " << R(0, 1) << "; " << R(0, 2)
                    << "); (" << R(1, 0) << "; " << R(1, 1) << "; " << R(1, 2)
                    << "); (" << R(2, 0) << "; " << R(2, 1) << "; " << R(2, 2)
                    << "))";
}

template <typename T>
struct is_readable_matrix<rot_mat_3D<T>> {
  static constexpr bool value = true;
  using type = is_readable_matrix<rot_mat_3D<T>>;
};

/// This class represents a rotation using quaternions (or Euler-Rodriguez parameters).
/// The convention used is with the leading scalar.
template <typename T>
class quaternion {
 public:
  using self = quaternion<T>;

  using value_type = T;
  using container_type = void;

  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;

  using iterator = void;
  using const_iterator = void;

  using size_type = int;
  using difference_type = int;

 private:
  std::array<value_type, 4> q = {};

  quaternion(const_reference q0, const_reference q1, const_reference q2,
             const_reference q3) noexcept {
    q[0] = q0;
    q[1] = q1;
    q[2] = q2;
    q[3] = q3;
  }

 public:
  friend class rot_mat_3D<value_type>;
  friend class euler_angles_TB<value_type>;
  friend class axis_angle<value_type>;
  friend class trans_mat_3D<value_type>;
  friend class unit_quat<value_type>;

  class xrot {
   private:
    value_type q0;
    value_type qx;

    xrot(const_reference Q0, const_reference QX) noexcept : q0(Q0), qx(QX) {}

   public:
    friend class quaternion<value_type>;  // befriend parent.

    explicit xrot(const_reference ang = value_type(0.0)) noexcept {
      set_angle(ang);
    }

    value_type s() const noexcept { return q0; }
    value_type v() const noexcept { return qx; }

    value_type get_angle() const noexcept {
      using std::atan2;
      return value_type(2.0) * atan2(qx, q0);
    }

    void set_angle(const_reference ang) noexcept {
      using std::cos;
      using std::sin;
      q0 = cos(ang * 0.5);
      qx = sin(ang * 0.5);
    }

    vect<value_type, 3> get_axis() const noexcept {
      return vect<value_type, 3>(value_type(1.0), value_type(0.0),
                                 value_type(0.0));
    }

    quaternion<value_type> getQuaternion() const noexcept {
      return quaternion<value_type>(q0, qx, value_type(0.0), value_type(0.0));
    }

    xrot& operator*=(const xrot& Q2) noexcept {
      value_type tmp = Q2.q0 * q0 - Q2.qx * qx;
      qx = Q2.q0 * qx + Q2.qx * q0;
      q0 = tmp;
      return *this;
    }

    friend xrot operator*(const xrot& Q1, const xrot& Q2) noexcept {
      return xrot(Q2.q0 * Q1.q0 - Q2.qx * Q1.qx, Q2.q0 * Q1.qx + Q2.qx * Q1.q0);
    }

    friend vect<value_type, 3> operator*(
        const xrot& Q, const vect<value_type, 3>& V) noexcept {
      value_type t0 = Q.q0 * Q.qx;
      value_type t3 = -Q.qx * Q.qx;
      return vect<value_type, 3>(
          V[0], value_type(2.0) * (t3 * V[1] - t0 * V[2]) + V[1],
          value_type(2.0) * (t0 * V[1] + t3 * V[2]) + V[2]);
    }
  };

  class yrot {
   private:
    value_type q0;
    value_type qy;

    yrot(const_reference Q0, const_reference QY) noexcept : q0(Q0), qy(QY) {}

   public:
    friend class quaternion<value_type>;  // befriend parent.

    explicit yrot(const_reference ang = value_type(0.0)) noexcept {
      set_angle(ang);
    }

    value_type s() const noexcept { return q0; }
    value_type v() const noexcept { return qy; }

    value_type get_angle() const noexcept {
      using std::atan2;
      return value_type(2.0) * atan2(qy, q0);
    }

    void set_angle(const_reference ang) noexcept {
      using std::cos;
      using std::sin;
      q0 = cos(ang * 0.5);
      qy = sin(ang * 0.5);
    }

    vect<value_type, 3> get_axis() const noexcept {
      return vect<value_type, 3>(value_type(0.0), value_type(1.0),
                                 value_type(0.0));
    }

    quaternion<value_type> getQuaternion() const noexcept {
      return quaternion<value_type>(q0, value_type(0.0), qy, value_type(0.0));
    }

    yrot& operator*=(const yrot& Q2) noexcept {
      value_type tmp = Q2.q0 * q0 - Q2.qy * qy;
      qy = Q2.q0 * qy + Q2.qy * q0;
      q0 = tmp;
      return *this;
    }

    friend yrot operator*(const yrot& Q1, const yrot& Q2) noexcept {
      return yrot(Q2.q0 * Q1.q0 - Q2.qy * Q1.qy, Q2.q0 * Q1.qy + Q2.qy * Q1.q0);
    }

    friend vect<value_type, 3> operator*(
        const yrot& Q, const vect<value_type, 3>& V) noexcept {
      value_type t1 = Q.q0 * Q.qy;
      value_type t6 = -Q.qy * Q.qy;
      return vect<value_type, 3>(
          value_type(2.0) * (t6 * V[0] + t1 * V[2]) + V[0], V[1],
          value_type(2.0) * (-t1 * V[0] + t6 * V[2]) + V[2]);
    }
  };

  class zrot {
   private:
    value_type q0;
    value_type qz;

    zrot(const_reference Q0, const_reference QZ) noexcept : q0(Q0), qz(QZ) {}

   public:
    friend class quaternion<value_type>;  // befriend parent.

    explicit zrot(const_reference ang = value_type(0.0)) noexcept {
      set_angle(ang);
    }

    value_type s() const noexcept { return q0; }
    value_type v() const noexcept { return qz; }

    value_type get_angle() const noexcept {
      using std::atan2;
      return value_type(2.0) * atan2(qz, q0);
    }

    void set_angle(const_reference ang) noexcept {
      using std::cos;
      using std::sin;
      q0 = cos(ang * 0.5);
      qz = sin(ang * 0.5);
    }

    vect<value_type, 3> get_axis() const noexcept {
      return vect<value_type, 3>(value_type(0.0), value_type(0.0),
                                 value_type(1.0));
    }

    quaternion<value_type> getQuaternion() const noexcept {
      return quaternion<value_type>(q0, value_type(0.0), value_type(0.0), qz);
    }

    zrot& operator*=(const zrot& Q2) noexcept {
      value_type tmp = Q2.q0 * q0 - Q2.qz * qz;
      qz = Q2.q0 * qz + Q2.qz * q0;
      q0 = tmp;
      return *this;
    }

    friend zrot operator*(const zrot& Q1, const zrot& Q2) noexcept {
      return zrot(Q2.q0 * Q1.q0 - Q2.qz * Q1.qz, Q2.q0 * Q1.qz + Q2.qz * Q1.q0);
    }

    friend vect<value_type, 3> operator*(
        const zrot& Q, const vect<value_type, 3>& V) noexcept {
      value_type t2 = Q.q0 * Q.qz;
      value_type t8 = -Q.qz * Q.qz;
      return vect<value_type, 3>(
          value_type(2.0) * (t8 * V[0] - t2 * V[1]) + V[0],
          value_type(2.0) * (t2 * V[0] + t8 * V[1]) + V[1], V[2]);
    }
  };

  friend self operator*(const xrot& q1, const yrot& q2) noexcept {
    return self(q2.s() * q1.s(), q2.s() * q1.v(), q2.v() * q1.s(),
                q2.v() * q1.v());
  }

  friend self operator*(const yrot& q1, const xrot& q2) noexcept {
    return self(q2.s() * q1.s(), q2.v() * q1.s(), q2.s() * q1.v(),
                -q2.v() * q1.v());
  }

  friend self operator*(const zrot& q1, const xrot& q2) noexcept {
    return self(q2.s() * q1.s(), q2.v() * q1.s(), q2.v() * q1.v(),
                q2.s() * q1.v());
  }

  friend self operator*(const xrot& q1, const zrot& q2) noexcept {
    return self(q2.s() * q1.s(), q2.s() * q1.v(), -q2.v() * q1.v(),
                q2.v() * q1.s());
  }

  friend self operator*(const yrot& q1, const zrot& q2) noexcept {
    return self(q2.s() * q1.s(), q2.v() * q1.v(), q2.s() * q1.v(),
                q2.v() * q1.s());
  }

  friend self operator*(const zrot& q1, const yrot& q2) noexcept {
    return self(q2.s() * q1.s(), -q2.v() * q1.v(), q2.v() * q1.s(),
                q2.s() * q1.v());
  }

  friend self operator*(const self& q1, const xrot& q2) noexcept {
    return self(q2.s() * q1.q[0] - q2.v() * q1.q[1],
                q2.s() * q1.q[1] + q2.v() * q1.q[0],
                q2.s() * q1.q[2] + q2.v() * q1.q[3],
                q2.s() * q1.q[3] - q2.v() * q1.q[2]);
  }

  friend self operator*(const self& q1, const yrot& q2) noexcept {
    return self(q2.s() * q1.q[0] - q2.v() * q1.q[2],
                q2.s() * q1.q[1] - q2.v() * q1.q[3],
                q2.s() * q1.q[2] + q2.v() * q1.q[0],
                q2.s() * q1.q[3] + q2.v() * q1.q[1]);
  }

  friend self operator*(const self& q1, const zrot& q2) noexcept {
    return self(q2.s() * q1.q[0] - q2.v() * q1.q[3],
                q2.s() * q1.q[1] + q2.v() * q1.q[2],
                q2.s() * q1.q[2] - q2.v() * q1.q[1],
                q2.s() * q1.q[3] + q2.v() * q1.q[0]);
  }

  friend self operator*(const xrot& q1, const self& q2) noexcept {
    return self(q2.q[0] * q1.s() - q2.q[1] * q1.v(),
                q2.q[1] * q1.s() + q2.q[0] * q1.v(),
                q2.q[2] * q1.s() - q2.q[3] * q1.v(),
                q2.q[3] * q1.s() + q2.q[2] * q1.v());
  }

  friend self operator*(const yrot& q1, const self& q2) noexcept {
    return self(q2.q[0] * q1.s() - q2.q[2] * q1.v(),
                q2.q[1] * q1.s() + q2.q[3] * q1.v(),
                q2.q[2] * q1.s() + q2.q[0] * q1.v(),
                q2.q[3] * q1.s() - q2.q[1] * q1.v());
  }

  friend self operator*(const zrot& q1, const self& q2) noexcept {
    return self(q2.q[0] * q1.s() - q2.q[3] * q1.v(),
                q2.q[1] * q1.s() - q2.q[2] * q1.v(),
                q2.q[2] * q1.s() + q2.q[1] * q1.v(),
                q2.q[3] * q1.s() + q2.q[0] * q1.v());
  }

  self& operator*=(const xrot& q2) noexcept { return *this = *this * q2; }

  self& operator*=(const yrot& q2) noexcept { return *this = *this * q2; }

  self& operator*=(const zrot& q2) noexcept { return *this = *this * q2; }

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /// Default Constructor.
  quaternion() noexcept : quaternion(1.0, 0.0, 0.0, 0.0) {}

  quaternion(const self&) noexcept = default;

  template <typename Vector>
  explicit quaternion(const Vector& aV) noexcept {
    static_assert(is_readable_vector_v<Vector>);
    auto v = unit(vect<value_type, 4>(aV[0], aV[1], aV[2], aV[3]));
    q[0] = v[0];
    q[1] = v[1];
    q[2] = v[2];
    q[3] = v[3];
  }

  /// Constructor from a rotation matrix.
  explicit quaternion(const rot_mat_3D<value_type>& R) noexcept
      : quaternion(R.getQuaternion()) {}

  explicit quaternion(const axis_angle<value_type>& A) noexcept
      : quaternion(A.getQuaternion()) {}

  explicit quaternion(const euler_angles_TB<value_type>& E) noexcept
      : quaternion(E.getQuaternion()) {}

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Provides the rotation matrix as an ordinary 3x3 matrix.
  mat<value_type, mat_structure::square> getMat() const {
    return getRotMat().getMat();
  }

  /// Provides the rotation matrix corresponding to the quaternion.
  rot_mat_3D<value_type> getRotMat() const noexcept {
    value_type t01(value_type(2.0) * q[0] * q[1]);
    value_type t02(value_type(2.0) * q[0] * q[2]);
    value_type t03(value_type(2.0) * q[0] * q[3]);
    value_type t11(value_type(2.0) * q[1] * q[1]);
    value_type t12(value_type(2.0) * q[1] * q[2]);
    value_type t13(value_type(2.0) * q[1] * q[3]);
    value_type t22(value_type(2.0) * q[2] * q[2]);
    value_type t23(value_type(2.0) * q[2] * q[3]);
    value_type t33(value_type(2.0) * q[3] * q[3]);
    return rot_mat_3D<value_type>(
        value_type(1.0) - t22 - t33, t12 - t03, t02 + t13, t12 + t03,
        value_type(1.0) - t11 - t33, t23 - t01, t13 - t02, t01 + t23,
        value_type(1.0) - t11 - t22);
  }

  euler_angles_TB<value_type> getEulerAnglesTB() const noexcept {
    using std::asin;
    using std::atan2;
    using std::cos;
    value_type yaw{};
    value_type pitch{};
    value_type roll{};
    pitch = value_type(2.0) * (q[0] * q[2] - q[1] * q[3]);
    if ((pitch != value_type(1.0)) && (pitch != value_type(-1.0))) {
      pitch = asin(pitch);
      value_type cp = value_type(1.0) / cos(pitch);
      roll = atan2(value_type(2.0) * cp * (q[2] * q[3] + q[0] * q[1]),
                   cp * (value_type(1.0) -
                         value_type(2.0) * (q[1] * q[1] + q[2] * q[2])));
      yaw = atan2(value_type(2.0) * cp * (q[1] * q[2] + q[0] * q[3]),
                  cp * (value_type(1.0) -
                        value_type(2.0) * (q[2] * q[2] + q[3] * q[3])));
    } else {
      yaw = value_type(0.0);
      roll = atan2(pitch * value_type(2.0) * (q[1] * q[2] - q[0] * q[3]),
                   pitch * value_type(2.0) * (q[1] * q[3] + q[0] * q[2]));
      pitch *= value_type(1.57079632679489662);
    }
    return euler_angles_TB<value_type>(yaw, pitch, roll);
  }

  axis_angle<value_type> getAxisAngle() const noexcept {
    using std::acos;
    using std::sqrt;
    vect<value_type, 4> v(q[0], q[1], q[2], q[3]);
    v = unit(v);
    value_type tmp(sqrt(v[1] * v[1] + v[2] * v[2] + v[3] * v[3]));
    vect<value_type, 3> axis{};
    value_type angle{};
    if (tmp > value_type(0.0000001)) {
      axis[0] = v[1] / tmp;
      axis[1] = v[2] / tmp;
      axis[2] = v[3] / tmp;
      if (v[0] < value_type(0.0)) {
        angle = value_type(2.0) * acos(-v[0]);
        axis = -axis;
      } else {
        angle = value_type(2.0) * acos(v[0]);
      }
    } else {
      axis[0] = value_type(1.0);
      axis[1] = value_type(0.0);
      axis[2] = value_type(0.0);
      angle = value_type(0.0);
    }
    return axis_angle<value_type>(angle, axis);
  }

  /// Array indexing operator, accessor for read only.
  const_reference operator[](int i) const noexcept {
    assert(i < 4);
    return q[i];
  }

  /*******************************************************************************
                         Assignment Operators
*******************************************************************************/

  self& operator=(const self&) noexcept = default;

  /// Assignment operator from a rotation matrix.
  template <typename Other>
  self& operator=(const Other& R) noexcept {
    return *this = self(R);
  }

  /// Multiply-and-store operator from a quaternion.
  self& operator*=(const self& Q) noexcept { return (*this = *this * Q); }

  /// Multiply-and-store operator from a rotation matrix.
  template <typename Other>
  self& operator*=(const Other& R) noexcept {
    return *this *= self{R};
  }

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /// Multiplication by a quaternion.
  friend self operator*(const self& Q1, const self& Q2) noexcept {
    return self(Q2.q[0] * Q1.q[0] - Q2.q[1] * Q1.q[1] - Q2.q[2] * Q1.q[2] -
                    Q2.q[3] * Q1.q[3],
                Q2.q[0] * Q1.q[1] + Q2.q[3] * Q1.q[2] - Q2.q[2] * Q1.q[3] +
                    Q2.q[1] * Q1.q[0],
                Q2.q[0] * Q1.q[2] - Q2.q[3] * Q1.q[1] + Q2.q[1] * Q1.q[3] +
                    Q2.q[2] * Q1.q[0],
                Q2.q[0] * Q1.q[3] + Q2.q[2] * Q1.q[1] - Q2.q[1] * Q1.q[2] +
                    Q2.q[3] * Q1.q[0]);
  }

  /// Multiplication by a matrix or other rotation representation.
  template <typename Other,
            std::enable_if_t<!std::is_same_v<rot_mat_3D<value_type>, Other>,
                             void*> = nullptr>
  friend auto operator*(const self& Q, const Other& R) {
    if constexpr (is_readable_matrix_v<Other>) {
      return Q.getRotMat() * R;
    } else {
      return Q * self{R};
    }
  }

  /// Multiplication by a matrix or other rotation representation.
  template <typename Other,
            std::enable_if_t<!std::is_same_v<rot_mat_3D<value_type>, Other>,
                             void*> = nullptr>
  friend auto operator*(const Other& R, const self& Q) {
    if constexpr (is_readable_matrix_v<Other>) {
      return R * Q.getRotMat();
    } else {
      return self{R} * Q;
    }
  }

  /// Multiplication by a column vector.
  friend vect<value_type, 3> operator*(const self& Q,
                                       const vect<value_type, 3>& V) noexcept {
    std::array<value_type, 9> t;
    t[0] = Q.q[0] * Q.q[1];
    t[1] = Q.q[0] * Q.q[2];
    t[2] = Q.q[0] * Q.q[3];
    t[3] = -Q.q[1] * Q.q[1];
    t[4] = Q.q[1] * Q.q[2];
    t[5] = Q.q[1] * Q.q[3];
    t[6] = -Q.q[2] * Q.q[2];
    t[7] = Q.q[2] * Q.q[3];
    t[8] = -Q.q[3] * Q.q[3];
    return {value_type(2.0) * ((t[6] + t[8]) * V[0] + (t[4] - t[2]) * V[1] +
                               (t[1] + t[5]) * V[2]) +
                V[0],
            value_type(2.0) * ((t[2] + t[4]) * V[0] + (t[3] + t[8]) * V[1] +
                               (t[7] - t[0]) * V[2]) +
                V[1],
            value_type(2.0) * ((t[5] - t[1]) * V[0] + (t[0] + t[7]) * V[1] +
                               (t[3] + t[6]) * V[2]) +
                V[2]};
  }

  friend vect<value_type, 3> operator*(
      const self& Q, const vect_component<value_type, 0>& x_value) noexcept {
    return {x_value.q - 2.0 * (Q.q[2] * Q.q[2] + Q.q[3] * Q.q[3]) * x_value.q,
            2.0 * (Q.q[0] * Q.q[3] + Q.q[1] * Q.q[2]) * x_value.q,
            2.0 * (Q.q[1] * Q.q[3] - Q.q[0] * Q.q[2]) * x_value.q};
  }

  friend vect<value_type, 3> operator*(
      const self& Q, const vect_component<value_type, 1>& y_value) noexcept {
    return {2.0 * (Q.q[1] * Q.q[2] - Q.q[0] * Q.q[3]) * y_value.q,
            y_value.q - 2.0 * (Q.q[1] * Q.q[1] + Q.q[3] * Q.q[3]) * y_value.q,
            2.0 * (Q.q[0] * Q.q[1] + Q.q[2] * Q.q[3]) * y_value.q};
  }

  friend vect<value_type, 3> operator*(
      const self& Q, const vect_component<value_type, 2>& z_value) noexcept {
    return {2.0 * (Q.q[0] * Q.q[2] + Q.q[1] * Q.q[3]) * z_value.q,
            2.0 * (Q.q[2] * Q.q[3] - Q.q[0] * Q.q[1]) * z_value.q,
            z_value.q - 2.0 * (Q.q[1] * Q.q[1] + Q.q[2] * Q.q[2]) * z_value.q};
  }

  /*******************************************************************************
                           Special Methods
  *******************************************************************************/

  /// Gets the time-derivative of the quaternion that corresponds to the angular velocity Omega.
  vect<value_type, 4> getQuaternionDot(
      const vect<value_type, 3>& Omega) const noexcept {
    return {-value_type(0.5) *
                (q[1] * Omega.q[0] + q[2] * Omega.q[1] + q[3] * Omega.q[2]),
            value_type(0.5) *
                (q[0] * Omega.q[0] - q[3] * Omega.q[1] + q[2] * Omega.q[2]),
            value_type(0.5) *
                (q[0] * Omega.q[1] + q[3] * Omega.q[0] - q[1] * Omega.q[2]),
            value_type(0.5) *
                (q[0] * Omega.q[2] - q[2] * Omega.q[0] + q[1] * Omega.q[1])};
  }

  /// Gets the angular velocity that corresponds to the time-derivative of the quaternion.
  vect<value_type, 3> getOmega(
      const vect<value_type, 4>& QuaternionDot) const noexcept {
    return {value_type(2.0) *
                (-q[1] * QuaternionDot.q[0] + q[0] * QuaternionDot.q[1] +
                 q[3] * QuaternionDot.q[2] - q[2] * QuaternionDot.q[3]),
            value_type(2.0) *
                (-q[2] * QuaternionDot.q[0] - q[3] * QuaternionDot.q[1] +
                 q[0] * QuaternionDot.q[2] + q[1] * QuaternionDot.q[3]),
            value_type(2.0) *
                (-q[3] * QuaternionDot.q[0] + q[2] * QuaternionDot.q[1] -
                 q[1] * QuaternionDot.q[2] + q[0] * QuaternionDot.q[3])};
  }

  /// Gets the 2-time-derivative of the quaternion that corresponds to the angular velocity Omega.
  vect<value_type, 4> getQuaternionDotDot(
      const vect<value_type, 4>& QD, const vect<value_type, 3>& W,
      const vect<value_type, 3>& WD) const noexcept {
    return {-value_type(0.5) *
                (q[1] * WD.q[0] + q[2] * WD.q[1] + q[3] * WD.q[2] +
                 QD.q[1] * W.q[0] + QD.q[2] * W.q[1] + QD.q[3] * W.q[2]),
            value_type(0.5) *
                (q[0] * WD.q[0] - q[3] * WD.q[1] + q[2] * WD.q[2] +
                 QD.q[0] * W.q[0] - QD.q[3] * W.q[1] + QD.q[2] * W.q[2]),
            value_type(0.5) *
                (q[3] * WD.q[0] + q[0] * WD.q[1] - q[1] * WD.q[2] +
                 QD.q[3] * W.q[0] + QD.q[0] * W.q[1] - QD.q[1] * W.q[2]),
            value_type(0.5) *
                (-q[2] * WD.q[0] + q[1] * WD.q[1] + q[0] * WD.q[2] -
                 QD.q[2] * W.q[0] + QD.q[1] * W.q[1] + QD.q[0] * W.q[2])};
  }

  /// Gets the angular acceleration that corresponds to the 2-time-derivative of the quaternion.
  vect<value_type, 3> getOmegaDot(
      const vect<value_type, 4>& QD,
      const vect<value_type, 4>& QDD) const noexcept {
    return {value_type(2.0) *
                (-q[1] * QDD.q[0] + q[0] * QDD.q[1] + q[3] * QDD.q[2] -
                 q[2] * QDD.q[3] - QD.q[1] * QD.q[0] + QD.q[0] * QD.q[1] +
                 QD.q[3] * QD.q[2] - QD.q[2] * QD.q[3]),
            value_type(2.0) *
                (-q[2] * QDD.q[0] - q[3] * QDD.q[1] + q[0] * QDD.q[2] +
                 q[1] * QDD.q[3] - QD.q[2] * QD.q[0] - QD.q[3] * QD.q[1] +
                 QD.q[0] * QD.q[2] + QD.q[1] * QD.q[3]),
            value_type(2.0) *
                (-q[3] * QDD.q[0] + q[2] * QDD.q[1] - q[1] * QDD.q[2] +
                 q[0] * QDD.q[3] - QD.q[3] * QD.q[0] + QD.q[2] * QD.q[1] -
                 QD.q[1] * QD.q[2] + QD.q[0] * QD.q[3])};
  }

  /*******************************************************************************
                           Standard Matrix Methods
  *******************************************************************************/

  /// Produces a transpose quaternion which is the inverse rotation.
  friend self transpose(const self& Q) noexcept {
    return {Q.q[0], -Q.q[1], -Q.q[2], -Q.q[3]};
  }

  /// Produces a cofactor matrix which is the same as the rotation matrix itself.
  friend self cofactor(const self& Q) noexcept { return Q; }

  /// Invert the rotation.
  friend self invert(const self& Q) noexcept {
    return {Q.q[0], -Q.q[1], -Q.q[2], -Q.q[3]};
  }

  /// Gets the trace of the matrix.
  friend value_type trace(const self& Q) noexcept {
    return value_type(4.0) * Q.q[0] * Q.q[0] - value_type(1.0);
  }

  /// Gets the determinant of the matrix.
  friend value_type determinant(const self& Q) noexcept {
    return value_type(1.0);
  }

  /// Gets the symmetric part of the matrix.
  mat<value_type, mat_structure::symmetric> getSymPart() const {
    value_type t11(value_type(2.0) * q[1] * q[1]);
    value_type t12(value_type(2.0) * q[1] * q[2]);
    value_type t13(value_type(2.0) * q[1] * q[3]);
    value_type t22(value_type(2.0) * q[2] * q[2]);
    value_type t23(value_type(2.0) * q[2] * q[3]);
    value_type t33(value_type(2.0) * q[3] * q[3]);
    return {value_type(1.0) - t22 - t33, t12, t13,
            value_type(1.0) - t11 - t33, t23, value_type(1.0) - t11 - t22};
  }

  /// Gets the skew-symmetric part of the matrix.
  mat<value_type, mat_structure::skew_symmetric> getSkewSymPart() const {
    value_type t01(value_type(2.0) * q[0] * q[1]);
    value_type t02(value_type(2.0) * q[0] * q[2]);
    value_type t03(value_type(2.0) * q[0] * q[3]);
    return {-t03, t02, -t01};
  }

  /// Loading a quaternion value with a name.
  friend serialization::iarchive& operator&(
      serialization::iarchive& in,
      const std::pair<std::string, quaternion<T>&>& R) {
    return in & RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_q0", R.second.q[0]) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_q1", R.second.q[1]) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_q2", R.second.q[2]) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_q3", R.second.q[3]);
  }

  /// Loading a quaternion value.
  friend serialization::iarchive& operator>>(serialization::iarchive& in,
                                             quaternion<T>& R) {
    return in & RK_SERIAL_LOAD_WITH_NAME(R);
  }

  /// Saving a quaternion value with a name.
  friend serialization::oarchive& operator&(
      serialization::oarchive& out,
      const std::pair<std::string, const quaternion<T>&>& R) {
    return out & RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_q0", R.second.q[0]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_q1", R.second.q[1]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_q2", R.second.q[2]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_q3", R.second.q[3]);
  }

  /// Saving a quaternion value.
  friend serialization::oarchive& operator<<(serialization::oarchive& out,
                                             const quaternion<T>& R) {
    return out & RK_SERIAL_SAVE_WITH_NAME(R);
  }
};

namespace rtti {

template <typename T>
struct get_type_id<quaternion<T>> {
  static constexpr unsigned int ID = 0x0000001A;
  static constexpr auto type_name = std::string_view{"ReaK::quaternion"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }

  using save_type = const quaternion<T>&;
  using load_type = quaternion<T>&;
};
}  // namespace rtti

/// Prints a quaternion to a standard output stream (<<) as "(q0; q1; q2; q3)".
template <class T>
std::ostream& operator<<(std::ostream& out_stream, const quaternion<T>& Q) {
  return (out_stream << "(" << Q[0] << "; " << Q[1] << "; " << Q[2] << "; "
                     << Q[3] << ")");
}

/// This class represents a rotation using Euler angles (Tait-Bryan), 321-body-fixed, in body frame.
template <class T>
class euler_angles_TB {
 public:
  using self = euler_angles_TB<T>;

  using value_type = T;
  using container_type = void;

  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;

  using size_type = int;
  using difference_type = int;

 private:
  std::array<value_type, 3> q;

 public:
  friend class rot_mat_3D<value_type>;
  friend class quaternion<value_type>;
  friend class axis_angle<value_type>;
  friend class trans_mat_3D<value_type>;

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /// Default Constructor.
  euler_angles_TB() noexcept {
    q[0] = 0.0;
    q[1] = 0.0;
    q[2] = 0.0;
  }

  /// Constructor from three euler angles.
  euler_angles_TB(const_reference Yaw_, const_reference Pitch_,
                  const_reference Roll_) noexcept {
    q[0] = Yaw_;
    q[1] = Pitch_;
    q[2] = Roll_;
  }

  euler_angles_TB(const self&) noexcept = default;

  /// Constructor from a quaternion.
  explicit euler_angles_TB(const quaternion<value_type>& Q) noexcept
      : euler_angles_TB(Q.getEulerAnglesTB()) {}

  explicit euler_angles_TB(const axis_angle<value_type>& A) noexcept
      : euler_angles_TB(A.getEulerAnglesTB()) {}

  /// Constructor from a rotation matrix.
  explicit euler_angles_TB(const rot_mat_3D<value_type>& R) noexcept
      : euler_angles_TB(R.getEulerAnglesTB()) {}

  euler_angles_TB(const rot_mat_3D<value_type>& R,
                  const euler_angles_TB<value_type>& /*Predicted*/) noexcept
      : euler_angles_TB(R.getEulerAnglesTB()) {}

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Get yaw, read-write.
  reference yaw() noexcept { return q[0]; }
  /// Get pitch, read-write.
  reference pitch() noexcept { return q[1]; }
  /// Get roll, read-write.
  reference roll() noexcept { return q[2]; }

  /// Get yaw, read-only.
  const_reference yaw() const noexcept { return q[0]; }
  /// Get pitch, read-only.
  const_reference pitch() const noexcept { return q[1]; }
  /// Get roll, read-only.
  const_reference roll() const noexcept { return q[2]; }

  /// Provides a quaternion corresponding to this rotation.
  quaternion<value_type> getQuaternion() const noexcept {
    using std::cos;
    using std::sin;
    value_type cpsi = cos(value_type(0.5) * q[0]);
    value_type spsi = sin(value_type(0.5) * q[0]);
    value_type ctheta = cos(value_type(0.5) * q[1]);

    value_type stheta = sin(value_type(0.5) * q[1]);
    value_type cphi = cos(value_type(0.5) * q[2]);
    value_type sphi = sin(value_type(0.5) * q[2]);

    return {cphi * ctheta * cpsi + sphi * stheta * spsi,
            sphi * ctheta * cpsi - cphi * stheta * spsi,
            cphi * stheta * cpsi + sphi * ctheta * spsi,
            cphi * ctheta * spsi - sphi * stheta * cpsi};
  }

  axis_angle<value_type> getAxisAngle() const noexcept {
    return axis_angle<value_type>(getQuaternion());
  }

  /// Provides a rotation matrix corresponding to this rotation.
  rot_mat_3D<value_type> getRotMat() const noexcept {
    using std::cos;
    using std::sin;
    value_type s1(sin(q[0]));
    value_type c1(cos(q[0]));
    value_type s2(sin(q[1]));
    value_type c2(cos(q[1]));
    value_type s3(sin(q[2]));
    value_type c3(cos(q[2]));

    return {c1 * c2,
            -(s1 * c3) + (c1 * s2 * s3),
            (s1 * s3) + (c1 * s2 * c3),
            s1 * c2,
            (c1 * c3) + (s1 * s2 * s3),
            -(c1 * s3) + (s1 * s2 * c3),
            -s2,
            c2 * s3,
            c2 * c3};
  }

  /// Provides a rotation matrix as a regular 3x3 matris corresponding to this rotation.
  mat<value_type, mat_structure::square> getMat() const {
    return getRotMat().getMat();
  }

  /*******************************************************************************
                         Assignment Operators
*******************************************************************************/

  self& operator=(const self&) noexcept = default;

  /// Assignment from other representation.
  template <typename Other>
  self& operator=(const Other& R) noexcept {
    return *this = self{R};
  }

  /// Multiply-and-store from another rotation.
  template <typename Other>
  self& operator*=(const Other& R) noexcept {
    return *this = *this * R;
  }

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /// Multiply by a euler angle representation.
  friend auto operator*(const self& E1, const self& E2) noexcept {
    return E1.getRotMat() * E2.getRotMat();
  }

  /// Multiply by a matrix or rotation representation.
  template <typename Other,
            std::enable_if_t<!std::is_same_v<Other, rot_mat_3D<value_type>> &&
                                 !std::is_same_v<Other, quaternion<value_type>>,
                             void*> = nullptr>
  friend auto operator*(const self& E, const Other& R) {
    return E.getRotMat() * R;
  }

  /// Multiply by a matrix or rotation representation.
  template <typename Other,
            std::enable_if_t<!std::is_same_v<Other, rot_mat_3D<value_type>> &&
                                 !std::is_same_v<Other, quaternion<value_type>>,
                             void*> = nullptr>
  friend auto operator*(const Other& R, const self& E) {
    return R * E.getRotMat();
  }

  /// Multiply by a vector, rotating it.
  friend vect<value_type, 3> operator*(const self& E,
                                       const vect<value_type, 3>& V) noexcept {
    return E.getRotMat() * V;
  }

  /*******************************************************************************
                           Standard Matrix Methods
  *******************************************************************************/
  /// Produces a transpose quaternion which is the inverse rotation.
  friend self transpose(const self& E) noexcept {
    using std::asin;
    using std::atan2;
    using std::cos;
    using std::sin;

    value_type s1(sin(E.q[0]));
    value_type c1(cos(E.q[0]));
    value_type s2(sin(E.q[1]));
    value_type c2(cos(E.q[1]));
    value_type s3(sin(E.q[2]));
    value_type c3(cos(E.q[2]));

    value_type R2((s1 * s3) + (c1 * s2 * c3));

    self result;

    if ((R2 != value_type(1.0)) && (R2 != value_type(-1.0))) {
      value_type R0(c1 * c2);
      value_type R1(-(s1 * c3) + (c1 * s2 * s3));
      value_type R5(-(c1 * s3) + (s1 * s2 * c3));
      value_type R8(c2 * c3);
      result.q[1] = asin(-R2);
      value_type cp = value_type(1.0) / cos(result.q[1]);
      result.q[2] = atan2(cp * R5, cp * R8);
      result.q[0] = atan2(cp * R1, cp * R0);
    } else {
      value_type R3(s1 * c2);
      result.q[0] = value_type(0.0);
      result.q[2] = atan2(-R2 * R3, R2 * s2);
      result.q[1] = -R2 * value_type(1.57079632679489662);
    };
    return result;
  }

  /// Produces a cofactor matrix which is the same as the rotation matrix itself.
  friend self cofactor(const self& E) noexcept { return E; }

  /// Invert the rotation.
  friend self invert(const self& E) noexcept { return transpose(E); }

  /// Gets the trace of the matrix.
  friend value_type trace(const self& E) noexcept {
    using std::cos;
    using std::sin;
    value_type t =
        cos(value_type(0.5) * E.q[2]) * cos(value_type(0.5) * E.q[1]) *
            cos(value_type(0.5) * E.q[0]) +
        sin(value_type(0.5) * E.q[2]) * sin(value_type(0.5) * E.q[1]) *
            sin(value_type(0.5) * E.q[0]);
    return value_type(4.0) * t * t - value_type(1.0);
  }

  /// Gets the determinant of the matrix.
  friend value_type determinant(const self& /*unused*/) noexcept {
    return value_type(1.0);
  }

  /// Gets the symmetric part of the matrix.
  mat<value_type, mat_structure::symmetric> getSymPart() const {
    return mat<value_type, mat_structure::symmetric>(this->getMat());
  }

  /// Gets the skew-symmetric part of the matrix.
  mat<value_type, mat_structure::skew_symmetric> getSkewSymPart() const {
    return mat<value_type, mat_structure::skew_symmetric>(this->getMat());
  }

  /// Loading a euler_angles_TB value with a name.
  friend serialization::iarchive& operator&(
      serialization::iarchive& in,
      const std::pair<std::string, euler_angles_TB<T>&>& R) {
    return in & RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_yaw", R.second.q[0]) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_pitch", R.second.q[1]) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_roll", R.second.q[2]);
  }

  /// Loading a euler_angles_TB value.
  friend serialization::iarchive& operator>>(serialization::iarchive& in,
                                             euler_angles_TB<T>& R) {
    return in & RK_SERIAL_LOAD_WITH_NAME(R);
  }

  /// Saving a euler_angles_TB value with a name.
  friend serialization::oarchive& operator&(
      serialization::oarchive& out,
      const std::pair<std::string, const euler_angles_TB<T>&>& R) {
    return out & RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_yaw", R.second.q[0]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_pitch", R.second.q[1]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_roll", R.second.q[2]);
  }

  /// Saving a euler_angles_TB value.
  friend serialization::oarchive& operator<<(serialization::oarchive& out,
                                             const euler_angles_TB<T>& R) {
    return out & RK_SERIAL_SAVE_WITH_NAME(R);
  }
};

namespace rtti {

template <typename T>
struct get_type_id<euler_angles_TB<T>> {
  static constexpr unsigned int ID = 0x0000001B;
  static constexpr auto type_name = std::string_view{"ReaK::euler_angles_TB"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }

  using save_type = const euler_angles_TB<T>&;
  using load_type = euler_angles_TB<T>&;
};

}  // namespace rtti

/// Prints a euler angles to a standard output stream (<<) as "(Yaw = value; Pitch = value; Roll = value)".
template <class T>
std::ostream& operator<<(std::ostream& out_stream,
                         const euler_angles_TB<T>& E) {
  return out_stream << "(Yaw = " << E.yaw() << "; Pitch = " << E.pitch()
                    << "; Roll = " << E.roll() << ")";
}

/**
 * This class is a 3D rotation represented by an axis and angle.
 * \test All tests for this class have been passed!
 */
template <class T>
class axis_angle {
 public:
  using self = axis_angle<T>;

  using value_type = T;
  using container_type = void;

  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;

  using size_type = int;
  using difference_type = int;

 private:
  value_type mAngle;
  vect<value_type, 3> mAxis;

 public:
  friend class trans_mat_3D<value_type>;

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /// Default Constructor.
  axis_angle() noexcept : mAngle(0.0), mAxis(1.0, 0.0, 0.0) {}

  /// Constructor from angle and axis.
  axis_angle(const value_type& aAngle,
             const vect<value_type, 3>& aAxis) noexcept
      : mAngle(aAngle) {
    using std::sqrt;
    value_type tmp = norm_2(aAxis);
    if (tmp > value_type(0.0000001)) {
      mAxis.q[0] = aAxis[0] / tmp;
      mAxis.q[1] = aAxis[1] / tmp;
      mAxis.q[2] = aAxis[2] / tmp;
    } else {
      mAxis.q[0] = value_type(1.0);
      mAxis.q[1] = value_type(0.0);
      mAxis.q[2] = value_type(0.0);
    }
  }

  axis_angle(const self& A) noexcept = default;

  /// Constructor from a quaternion.
  explicit axis_angle(const quaternion<value_type>& Q) noexcept
      : axis_angle(Q.getAxisAngle()) {}

  /// Constructor from a rotation matrix.
  explicit axis_angle(const rot_mat_3D<value_type>& R) noexcept
      : axis_angle(R.getAxisAngle()) {}

  /// Constructor from euler angles.
  explicit axis_angle(const euler_angles_TB<value_type>& E) noexcept
      : axis_angle(E.getAxisAngle()) {}

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Provides the angle, read-write.
  reference angle() noexcept { return mAngle; }

  /// Provides the axis, read-write.
  vect<value_type, 3>& axis() noexcept { return mAxis; }

  /// Provides the angle, read-only.
  const_reference angle() const noexcept { return mAngle; }

  /// Provides the axis, read-only.
  const vect<value_type, 3>& axis() const noexcept { return mAxis; }

  /// Provides a quaternion representation of this rotation.
  quaternion<value_type> getQuaternion() const noexcept {
    using std::cos;
    using std::sin;
    value_type t = norm_2(mAxis);
    if (t == value_type(0.0)) {
      return {value_type(1.0), value_type(0.0), value_type(0.0),
              value_type(0.0)};
    }
    t = sin(value_type(0.5) * mAngle);
    return {cos(value_type(0.5) * mAngle), mAxis.q[0] * t, mAxis.q[1] * t,
            mAxis.q[2] * t};
  }

  /// Provides a euler angles representation of this rotation.
  euler_angles_TB<value_type> getEulerAnglesTB() const noexcept {
    using std::asin;
    using std::atan2;
    using std::cos;
    using std::sin;
    std::array<value_type, 4> quat;
    value_type t = norm_2(mAxis);
    if (t == value_type(0.0)) {
      quat[0] = value_type(1.0);
      quat[1] = value_type(0.0);
      quat[2] = value_type(0.0);
      quat[3] = value_type(0.0);
    } else {
      quat[0] = sin(mAngle / value_type(2.0));
      quat[1] = mAxis.q[0] * quat[0];
      quat[2] = mAxis.q[1] * quat[0];
      quat[3] = mAxis.q[2] * quat[0];
      quat[0] = cos(mAngle / value_type(2.0));
    }
    euler_angles_TB<value_type> result;
    result.q[1] = value_type(2.0) * (quat[0] * quat[2] - quat[1] * quat[3]);
    if ((result.q[1] != value_type(1.0)) && (result.q[1] != value_type(-1.0))) {
      result.q[1] = asin(result.q[1]);
      value_type cp = value_type(1.0) / cos(result.q[1]);
      result.q[2] = atan2(
          value_type(2.0) * cp * (quat[2] * quat[3] + quat[0] * quat[1]),
          cp * (value_type(1.0) -
                value_type(2.0) * (quat[1] * quat[1] + quat[2] * quat[2])));
      result.q[0] = atan2(
          value_type(2.0) * cp * (quat[1] * quat[2] + quat[0] * quat[3]),
          cp * (value_type(1.0) -
                value_type(2.0) * (quat[2] * quat[2] + quat[3] * quat[3])));
    } else {
      result.q[0] = value_type(0.0);
      result.q[2] = atan2(result.q[1] * value_type(2.0) *
                              (quat[1] * quat[2] - quat[0] * quat[3]),
                          result.q[1] * value_type(2.0) *
                              (quat[1] * quat[3] + quat[0] * quat[2]));
      result.q[1] *= value_type(1.57079632679489662);
    }
    return result;
  }

  /// Provides a rotation matrix representation of this rotation.
  rot_mat_3D<value_type> getRotMat() const noexcept {
    using std::cos;
    using std::sin;
    value_type ca(cos(mAngle));
    value_type one_minus_ca(value_type(1.0) - ca);
    value_type t11(ca + one_minus_ca * mAxis.q[0] * mAxis.q[0]);
    value_type t22(ca + one_minus_ca * mAxis.q[1] * mAxis.q[1]);
    value_type t33(ca + one_minus_ca * mAxis.q[2] * mAxis.q[2]);
    value_type t12(one_minus_ca * mAxis.q[0] * mAxis.q[1]);
    value_type t13(one_minus_ca * mAxis.q[0] * mAxis.q[2]);
    value_type t23(one_minus_ca * mAxis.q[1] * mAxis.q[2]);
    value_type sin_a(sin(mAngle));
    value_type t01(sin_a * mAxis.q[0]);
    value_type t02(sin_a * mAxis.q[1]);
    value_type t03(sin_a * mAxis.q[2]);

    return {t11,       t12 - t03, t13 + t02, t12 + t03, t22,
            t23 - t01, t13 - t02, t23 + t01, t33};
  }

  /// Provides a 3x3 matrix representation of this rotation.
  mat<value_type, mat_structure::square> getMat() const {
    return getRotMat().getMat();
  }

  /*******************************************************************************
                         Assignment Operators
*******************************************************************************/

  self& operator=(const self& A) noexcept = default;

  /// Standard assignment operator.
  template <typename Other>
  self& operator=(const Other& R) noexcept {
    return *this = self{R};
  }

  /// Multiply-and-store from a axis / angle.
  self& operator*=(const self& A) noexcept {
    return *this = this->getQuaternion() * A.getQuaternion();
  }

  /// Multiply-and-store from a euler angles.
  self& operator*=(const euler_angles_TB<value_type>& E) noexcept {
    return (*this = (this->getRotMat() * E.getRotMat()));
  }

  /// Multiply-and-store from a rotation matrix.
  self& operator*=(const rot_mat_3D<value_type>& R) noexcept {
    return (*this = (this->getRotMat() * R));
  }

  /// Multiply-and-store from a quaternion.
  self& operator*=(const quaternion<value_type>& Q) noexcept {
    return (*this = (this->getQuaternion() * Q));
  }

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /// Multiplication with an axis / angle representation.
  friend quaternion<value_type> operator*(const self& A1,
                                          const self& A2) noexcept {
    return A1.getQuaternion() * A2.getQuaternion();
  }

  /// Multiply by a matrix.
  template <
      typename Other,
      std::enable_if_t<!std::is_same_v<Other, rot_mat_3D<value_type>> &&
                           !std::is_same_v<Other, quaternion<value_type>> &&
                           !std::is_same_v<Other, euler_angles_TB<value_type>>,
                       void*> = nullptr>
  friend auto operator*(const self& A, const Other& R) {
    return A.getRotMat() * R;
  }

  /// Multiply by a matrix.
  template <
      typename Other,
      std::enable_if_t<!std::is_same_v<Other, rot_mat_3D<value_type>> &&
                           !std::is_same_v<Other, quaternion<value_type>> &&
                           !std::is_same_v<Other, euler_angles_TB<value_type>>,
                       void*> = nullptr>
  friend auto operator*(const Other& R, const self& A) {
    return R * A.getRotMat();
  }

  /// Multiplication with a column vector.
  friend vect<value_type, 3> operator*(const self& A,
                                       const vect<value_type, 3>& V) noexcept {
    return A.getRotMat() * V;
  }

  /*******************************************************************************
                           Standard Matrix Methods
  *******************************************************************************/

  /// Produces a transpose axis/angle which is the inverse rotation.
  friend self transpose(const self& A) noexcept { return {-A.mAngle, A.mAxis}; }

  /// Produces a cofactor matrix which is the same as the rotation itself.
  friend self cofactor(const self& A) noexcept { return A; }

  /// Invert the rotation.
  friend self invert(const self& A) noexcept { return {-A.mAngle, A.mAxis}; }

  /// Gets the trace of the matrix.
  friend value_type trace(const self& A) noexcept {
    using std::cos;
    return value_type(2.0) * cos(A.mAngle) + value_type(1.0);
  }

  /// Gets the determinant of the matrix.
  friend value_type determinant(const self& /*unused*/) noexcept {
    return value_type(1.0);
  }

  /// Gets the symmetric part of the matrix.
  mat<value_type, mat_structure::symmetric> getSymPart() const {
    using std::cos;
    value_type ca(cos(mAngle));
    value_type one_minus_ca(value_type(1.0) - ca);
    value_type t11(ca + one_minus_ca * mAxis.q[0] * mAxis.q[0]);
    value_type t22(ca + one_minus_ca * mAxis.q[1] * mAxis.q[1]);
    value_type t33(ca + one_minus_ca * mAxis.q[2] * mAxis.q[2]);
    value_type t12(one_minus_ca * mAxis.q[0] * mAxis.q[1]);
    value_type t13(one_minus_ca * mAxis.q[0] * mAxis.q[2]);
    value_type t23(one_minus_ca * mAxis.q[1] * mAxis.q[2]);
    return {t11, t12, t13, t22, t23, t33};
  }

  /// Gets the skew-symmetric part of the matrix.
  mat<value_type, mat_structure::skew_symmetric> getSkewSymPart() const {
    using std::sin;
    value_type sin_a(sin(mAngle));
    value_type t01(sin_a * mAxis.q[0]);
    value_type t02(sin_a * mAxis.q[1]);
    value_type t03(sin_a * mAxis.q[2]);
    return mat<value_type, mat_structure::skew_symmetric>(-t03, t02, -t01);
  }

  /// Loading a axis_angle value with a name.
  friend serialization::iarchive& operator&(
      serialization::iarchive& in,
      const std::pair<std::string, axis_angle<T>&>& R) {
    return in & RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_angle", R.second.mAngle) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_axis", R.second.mAxis);
  }

  /// Loading a axis_angle value.
  friend serialization::iarchive& operator>>(serialization::iarchive& in,
                                             axis_angle<T>& R) {
    return in & RK_SERIAL_LOAD_WITH_NAME(R);
  }

  /// Saving a axis_angle value with a name.
  friend serialization::oarchive& operator&(
      serialization::oarchive& out,
      const std::pair<std::string, const axis_angle<T>&>& R) {
    return out &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_angle", R.second.mAngle) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_axis", R.second.mAxis);
  }

  /// Saving a axis_angle value.
  friend serialization::oarchive& operator<<(serialization::oarchive& out,
                                             const axis_angle<T>& R) {
    return out & RK_SERIAL_SAVE_WITH_NAME(R);
  }
};

namespace rtti {

template <typename T>
struct get_type_id<axis_angle<T>> {
  static constexpr unsigned int ID = 0x0000001C;
  static constexpr auto type_name = std::string_view{"ReaK::axis_angle"};
  static construct_ptr CreatePtr() noexcept { return nullptr; };

  using save_type = const axis_angle<T>&;
  using load_type = axis_angle<T>&;
};

}  // namespace rtti

/// Prints a axis / angle to a standard output stream (<<) as "(Angle = value; Axis = (a1; a2; a3))".
template <class T>
std::ostream& operator<<(std::ostream& out_stream, const axis_angle<T>& A) {
  return out_stream << "(Angle = " << A.angle() << "; Axis = " << A.axis()
                    << ")";
}

/// This class is a transformation matrix 4 by 4, i.e. to rotate and translate a 3D vector.
template <class T>
class trans_mat_3D {
 public:
  using self = trans_mat_3D<T>;

  using value_type = T;
  using container_type = void;

  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;

  using col_iterator = void;
  using const_col_iterator = void;
  using row_iterator = void;
  using const_row_iterator = void;

  using size_type = int;
  using difference_type = int;

  static constexpr unsigned int static_row_count = 4;
  static constexpr unsigned int static_col_count = 4;
  static constexpr mat_alignment::tag alignment = mat_alignment::column_major;
  static constexpr mat_structure::tag structure = mat_structure::orthogonal;

  using translation_type = vect<value_type, 3>;

 private:
  std::array<value_type, 16> q = {};

  trans_mat_3D(const_reference a11, const_reference a12, const_reference a13,
               const_reference a14, const_reference a21, const_reference a22,
               const_reference a23, const_reference a24, const_reference a31,
               const_reference a32, const_reference a33,
               const_reference a34) noexcept {
    q[0] = a11;
    q[1] = a21;
    q[2] = a31;
    q[3] = value_type(0.0);
    q[4] = a12;
    q[5] = a22;
    q[6] = a32;
    q[7] = value_type(0.0);
    q[8] = a13;
    q[9] = a23;
    q[10] = a33;
    q[11] = value_type(0.0);
    q[12] = a14;
    q[13] = a24;
    q[14] = a34;
    q[15] = value_type(1.0);
  }

 public:
  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /// Default Constructor.
  trans_mat_3D() noexcept {
    std::fill(q.begin(), q.end(), value_type(0.0));
    q[0] = value_type(1.0);
    q[5] = value_type(1.0);
    q[10] = value_type(1.0);
    q[15] = value_type(1.0);
  }

  /// Constructor from a 4x4 array (16 values).
  explicit trans_mat_3D(const_pointer M) noexcept {
    vect<value_type, 3> v1 = unit(vect<value_type, 3>(M));
    q[0] = v1[0];
    q[1] = v1[1];
    q[2] = v1[2];
    q[3] = value_type(0.0);
    vect<value_type, 3> v2(M[4], M[5], M[6]);
    v2 = unit(v2 - (v2 * v1) * v1);
    q[4] = v2[0];
    q[5] = v2[1];
    q[6] = v2[2];
    q[7] = value_type(0.0);
    v2 = v1 % v2;
    q[8] = v2[0];
    q[9] = v2[1];
    q[10] = v2[2];
    q[11] = value_type(0.0);
    q[12] = M[12];
    q[13] = M[13];
    q[14] = M[14];
    q[15] = value_type(1.0);
  }
  explicit trans_mat_3D(pointer M)
      : trans_mat_3D(static_cast<const_pointer>(M)) {}

  /// Constructor from a regular matrix.
  template <typename Matrix>
  explicit trans_mat_3D(const Matrix& M) {
    static_assert(is_readable_matrix_v<Matrix>);
    if ((M.get_row_count() != 4) || (M.get_col_count() != 4)) {
      throw std::range_error(
          "Matrix for creating the 3D transformation matrix is not of correct "
          "dimensions!");
    }
    vect<value_type, 3> v1 =
        unit(vect<value_type, 3>(M(0, 0), M(1, 0), M(2, 0)));
    q[0] = v1[0];
    q[1] = v1[1];
    q[2] = v1[2];
    q[3] = value_type(0.0);
    vect<value_type, 3> v2(M(0, 1), M(1, 1), M(2, 1));
    v2 = unit(v2 - (v2 * v1) * v1);
    q[4] = v2[0];
    q[5] = v2[1];
    q[6] = v2[2];
    q[7] = value_type(0.0);
    v2 = v1 % v2;
    q[8] = v2[0];
    q[9] = v2[1];
    q[10] = v2[2];
    q[11] = value_type(0.0);
    q[12] = M(0, 3);
    q[13] = M(1, 3);
    q[14] = M(2, 3);
    q[15] = value_type(1.0);
  }

  /// Constructor from a rotation matrix and an optional translation vector V.
  explicit trans_mat_3D(const rot_mat_3D<value_type>& R,
                        const translation_type& V =
                            translation_type(value_type(0.0), value_type(0.0),
                                             value_type(0.0))) noexcept {
    setRotMat(R);
    q[3] = value_type(0.0);
    q[7] = value_type(0.0);
    q[11] = value_type(0.0);
    setTranslation(V);
    q[15] = value_type(1.0);
  }

  /// Constructor from a quaternion representation and an optional translation vector V.
  explicit trans_mat_3D(const quaternion<value_type>& Q,
                        const translation_type& V =
                            translation_type(value_type(0.0), value_type(0.0),
                                             value_type(0.0))) noexcept
      : trans_mat_3D(Q.getRotMat(), V) {}

  /// Constructor from a euler angles TB representation and an optional translation vector V.
  explicit trans_mat_3D(const euler_angles_TB<value_type>& E,
                        const translation_type& V =
                            translation_type(value_type(0.0), value_type(0.0),
                                             value_type(0.0))) noexcept
      : trans_mat_3D(E.getRotMat(), V) {}

  /// Constructor from an axis / angle representation and an optional translation vector V.
  explicit trans_mat_3D(const axis_angle<value_type>& A,
                        const translation_type& V =
                            translation_type(value_type(0.0), value_type(0.0),
                                             value_type(0.0))) noexcept
      : trans_mat_3D(A.getRotMat(), V) {}

  trans_mat_3D(const self&) noexcept = default;

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /// Provides a copy of the transformation matrix as an ordinary 4x4 matrix.
  mat<value_type, mat_structure::square> getMat() const {
    return mat<value_type, mat_structure::square>(
        q[0], q[4], q[8], q[12], q[1], q[5], q[9], q[13], q[2], q[6], q[10],
        q[14], value_type(0.0), value_type(0.0), value_type(0.0),
        value_type(1.0));
  }

  /// Provides the rotation part of the transformation as a rotation matrix.
  rot_mat_3D<value_type> getRotMat() const noexcept {
    return rot_mat_3D<value_type>(q[0], q[4], q[8], q[1], q[5], q[9], q[2],
                                  q[6], q[10]);
  }

  /// Sets the rotation part of the transformation from a rotation matrix.
  void setRotMat(const rot_mat_3D<value_type>& R) noexcept {
    q[0] = R.q[0];
    q[1] = R.q[1];
    q[2] = R.q[2];
    q[4] = R.q[3];
    q[5] = R.q[4];
    q[6] = R.q[5];
    q[8] = R.q[6];
    q[9] = R.q[7];
    q[10] = R.q[8];
  }

  /// Returns the quaternion of the rotation matrix.
  quaternion<value_type> getQuaternion() const noexcept {
    return quaternion<value_type>(getRotMat());
  }

  /// Sets the quaternion of the rotation matrix.
  void setQuaternion(const quaternion<value_type>& Q) noexcept {
    setRotMat(Q.getRotMat());
  }

  /// Returns the euler angles TB of the rotation matrix.
  euler_angles_TB<value_type> getEulerAnglesTB() const noexcept {
    return euler_angles_TB<value_type>(getRotMat());
  }

  /// Sets the euler angles TB of the rotation matrix.
  void setEulerAnglesTB(const euler_angles_TB<value_type>& E) noexcept {
    setRotMat(E.getRotMat());
  }

  /// Returns the axis / angle of the rotation matrix.
  axis_angle<value_type> getAxisAngle() const noexcept {
    return axis_angle<value_type>(getRotMat());
  }

  /// Sets the axis / angle of the rotation matrix.
  void setAxisAngle(const axis_angle<value_type>& A) noexcept {
    setRotMat(A.getRotMat());
  }

  /// Provides the translation part of the transformation matrix as a vector.
  translation_type getTranslation() const noexcept {
    return translation_type(q[12], q[13], q[14]);
  }

  /// Sets the translation part of the transformation matrix to a vector.
  void setTranslation(const translation_type& Translation) noexcept {
    q[12] = Translation[0];
    q[13] = Translation[1];
    q[14] = Translation[2];
  }

  /// Array indexing operator, accessor for read only.
  const_reference operator[](int i) const noexcept {
    assert(i < 16);
    return q[i];
  }

  /// Array double-indexing operator, ith row and jth column, accessor for read only.
  const_reference operator()(int i, int j) const noexcept {
    assert((i < 4) || (j < 4));
    return q[j * 4 + i];
  }

  int get_row_count() const noexcept { return 4; }
  int get_col_count() const noexcept { return 4; }

  /*******************************************************************************
                         Assignment Operators
*******************************************************************************/

  self& operator=(const self&) noexcept = default;

  /// Assignment operator with regular matrix.
  template <typename Matrix>
  self& operator=(const Matrix& M) {
    return *this = self{M};
  }

  /// Assignment operator with rotation matrix.
  self& operator=(const rot_mat_3D<value_type>& R) noexcept {
    return *this = self{R};
  }

  /// Assignment operator with a quaternion representation.
  self& operator=(const quaternion<value_type>& Q) noexcept {
    return *this = self{Q.getRotMat()};
  }

  /// Assignment operator with euler angles TB representation.
  self& operator=(const euler_angles_TB<value_type>& E) noexcept {
    return *this = self{E.getRotMat()};
  }

  /// Assignment operator with an axis / angle representation.
  self& operator=(const axis_angle<value_type>& A) noexcept {
    return *this = self{A.getRotMat()};
  }

  /// Multiply-and-store with a transformation matrix.
  self& operator*=(const self& M) noexcept {
    (*this) = self(q[0] * M.q[0] + q[4] * M.q[1] + q[8] * M.q[2],
                   q[0] * M.q[4] + q[4] * M.q[5] + q[8] * M.q[6],
                   q[0] * M.q[8] + q[4] * M.q[9] + q[8] * M.q[10],
                   q[0] * M.q[12] + q[4] * M.q[13] + q[8] * M.q[14] + q[12],
                   q[1] * M.q[0] + q[5] * M.q[1] + q[9] * M.q[2],
                   q[1] * M.q[4] + q[5] * M.q[5] + q[9] * M.q[6],
                   q[1] * M.q[8] + q[5] * M.q[9] + q[9] * M.q[10],
                   q[1] * M.q[12] + q[5] * M.q[13] + q[9] * M.q[14] + q[13],
                   q[2] * M.q[0] + q[6] * M.q[1] + q[10] * M.q[2],
                   q[2] * M.q[4] + q[6] * M.q[5] + q[10] * M.q[6],
                   q[2] * M.q[8] + q[6] * M.q[9] + q[10] * M.q[10],
                   q[2] * M.q[12] + q[6] * M.q[13] + q[10] * M.q[14] + q[14]);
    return *this;
  }

  /// Multiply-and-store with a rotation matrix.
  self& operator*=(const rot_mat_3D<value_type>& R) noexcept {
    (*this) = self(q[0] * R.q[0] + q[4] * R.q[1] + q[8] * R.q[2],
                   q[0] * R.q[3] + q[4] * R.q[4] + q[8] * R.q[5],
                   q[0] * R.q[6] + q[4] * R.q[7] + q[8] * R.q[8], q[12],
                   q[1] * R.q[0] + q[5] * R.q[1] + q[9] * R.q[2],
                   q[1] * R.q[3] + q[5] * R.q[4] + q[9] * R.q[5],
                   q[1] * R.q[6] + q[5] * R.q[7] + q[9] * R.q[8], q[13],
                   q[2] * R.q[0] + q[6] * R.q[1] + q[10] * R.q[2],
                   q[2] * R.q[3] + q[6] * R.q[4] + q[10] * R.q[5],
                   q[2] * R.q[6] + q[6] * R.q[7] + q[10] * R.q[8], q[14]);
    return *this;
  }

  /// Multiply-and-store with a quaternion representation.
  self& operator*=(const quaternion<value_type>& Q) noexcept {
    return (*this *= Q.getRotMat());
  }

  /// Multiply-and-store with a euler angles TB representation.
  self& operator*=(const euler_angles_TB<value_type>& E) noexcept {
    return (*this *= E.getRotMat());
  }

  /// Multiply-and-store with an axis / angle representation.
  self& operator*=(const axis_angle<value_type>& A) noexcept {
    return (*this *= A.getRotMat());
  }

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /// Multiplication with a transformation matrix.
  friend self operator*(const self& M1, const self& M2) noexcept {
    return {
        M1.q[0] * M2.q[0] + M1.q[4] * M2.q[1] + M1.q[8] * M2.q[2],
        M1.q[0] * M2.q[4] + M1.q[4] * M2.q[5] + M1.q[8] * M2.q[6],
        M1.q[0] * M2.q[8] + M1.q[4] * M2.q[9] + M1.q[8] * M2.q[10],
        M1.q[0] * M2.q[12] + M1.q[4] * M2.q[13] + M1.q[8] * M2.q[14] + M1.q[12],
        M1.q[1] * M2.q[0] + M1.q[5] * M2.q[1] + M1.q[9] * M2.q[2],
        M1.q[1] * M2.q[4] + M1.q[5] * M2.q[5] + M1.q[9] * M2.q[6],
        M1.q[1] * M2.q[8] + M1.q[5] * M2.q[9] + M1.q[9] * M2.q[10],
        M1.q[1] * M2.q[12] + M1.q[5] * M2.q[13] + M1.q[9] * M2.q[14] + M1.q[13],
        M1.q[2] * M2.q[0] + M1.q[6] * M2.q[1] + M1.q[10] * M2.q[2],
        M1.q[2] * M2.q[4] + M1.q[6] * M2.q[5] + M1.q[10] * M2.q[6],
        M1.q[2] * M2.q[8] + M1.q[6] * M2.q[9] + M1.q[10] * M2.q[10],
        M1.q[2] * M2.q[12] + M1.q[6] * M2.q[13] + M1.q[10] * M2.q[14] +
            M1.q[14]};
  }

  /// Multiply by a matrix.
  template <typename Matrix>
  friend auto operator*(const self& M1, const Matrix& M2) {
    static_assert(is_readable_matrix_v<Matrix>);
    return M1.getMat() * M2;
  }

  /// Multiply by a matrix.
  template <typename Matrix>
  friend auto operator*(const Matrix& M1, const self& M2) {
    static_assert(is_readable_matrix_v<Matrix>);
    return M1 * M2.getMat();
  }

  /// Multiplication with a rotation matrix.
  friend self operator*(const self& M,
                        const rot_mat_3D<value_type>& R) noexcept {
    return {M.q[0] * R[0] + M.q[4] * R[1] + M.q[8] * R[2],
            M.q[0] * R[3] + M.q[4] * R[4] + M.q[8] * R[5],
            M.q[0] * R[6] + M.q[4] * R[7] + M.q[8] * R[8],
            M.q[12],
            M.q[1] * R[0] + M.q[5] * R[1] + M.q[9] * R[2],
            M.q[1] * R[3] + M.q[5] * R[4] + M.q[9] * R[5],
            M.q[1] * R[6] + M.q[5] * R[7] + M.q[9] * R[8],
            M.q[13],
            M.q[2] * R[0] + M.q[6] * R[1] + M.q[10] * R[2],
            M.q[2] * R[3] + M.q[6] * R[4] + M.q[10] * R[5],
            M.q[2] * R[6] + M.q[6] * R[7] + M.q[10] * R[8],
            M.q[14]};
  }

  /// Multiplication with a transformation matrix.
  friend self operator*(const rot_mat_3D<value_type>& R,
                        const self& M) noexcept {
    return trans_mat_3D<value_type>(R) * M;
  }

  /// Multiplication with a 3D column vector.
  friend vect<value_type, 3> operator*(const self& M,
                                       const vect<value_type, 3>& V) noexcept {
    return {M.q[0] * V[0] + M.q[4] * V[1] + M.q[8] * V[2] + M.q[12],
            M.q[1] * V[0] + M.q[5] * V[1] + M.q[9] * V[2] + M.q[13],
            M.q[2] * V[0] + M.q[6] * V[1] + M.q[10] * V[2] + M.q[14]};
  }

  /// Multiplication with a 4D column vector.
  friend vect<value_type, 4> operator*(const self& M,
                                       const vect<value_type, 4>& V) noexcept {
    return {M.q[0] * V[0] + M.q[4] * V[1] + M.q[8] * V[2] + M.q[12] * V[3],
            M.q[1] * V[0] + M.q[5] * V[1] + M.q[9] * V[2] + M.q[13] * V[3],
            M.q[2] * V[0] + M.q[6] * V[1] + M.q[10] * V[2] + M.q[14] * V[3],
            V[3]};
  }

  /*******************************************************************************
                           Special Methods
  *******************************************************************************/

  /// Rotate-only a 3D column vector.
  vect<value_type, 3> rotate(const vect<value_type, 3>& V) const noexcept {
    return {q[0] * V[0] + q[4] * V[1] + q[8] * V[2],
            q[1] * V[0] + q[5] * V[1] + q[9] * V[2],
            q[2] * V[0] + q[6] * V[1] + q[10] * V[2]};
  }

  /*******************************************************************************
                           Standard Matrix Methods
  *******************************************************************************/

  /// Creates the transpose matrix.
  /// \note the matrix is no longer a transformation matrix.
  friend mat<value_type, mat_structure::square> transpose(
      const self& M) noexcept {
    return {M.q[0], M.q[1], M.q[2],  0.0, M.q[4],  M.q[5],  M.q[6],  0.0,
            M.q[8], M.q[9], M.q[10], 0.0, M.q[12], M.q[13], M.q[14], 1.0};
  }

  /// Gets the trace of the matrix.
  friend value_type trace(const self& M) noexcept {
    return M.q[0] + M.q[5] + M.q[10] + value_type(1.0);
  }

  /// Gets the determinant of the matrix.
  friend value_type determinant(const self& /*unused*/) noexcept {
    return value_type(1.0);
  }

  /// Invert the transformation.
  friend self invert(const self& M) noexcept {
    return {M.q[0],  M.q[1],
            M.q[2],  -M.q[0] * M.q[12] - M.q[1] * M.q[13] - M.q[2] * M.q[14],
            M.q[4],  M.q[5],
            M.q[6],  -M.q[4] * M.q[12] - M.q[5] * M.q[13] - M.q[6] * M.q[14],
            M.q[8],  M.q[9],
            M.q[10], -M.q[8] * M.q[12] - M.q[9] * M.q[13] - M.q[10] * M.q[14]};
  }

  /// Gets the symmetric part of the matrix.
  mat<value_type, mat_structure::symmetric> getSymPart() const {
    return mat<value_type, mat_structure::symmetric>(getMat());
  }

  /// Gets the skew-symmetric part of the matrix.
  mat<value_type, mat_structure::skew_symmetric> getSkewSymPart() const {
    return mat<value_type, mat_structure::skew_symmetric>(getMat());
  }

  /// Loading a trans_mat_3D value with a name.
  friend serialization::iarchive& operator&(
      serialization::iarchive& in,
      const std::pair<std::string, trans_mat_3D<T>&>& M) {
    in& RK_SERIAL_LOAD_WITH_ALIAS(M.first + "_r11", M.second.q[0]) &
        RK_SERIAL_LOAD_WITH_ALIAS(M.first + "_r21", M.second.q[1]) &
        RK_SERIAL_LOAD_WITH_ALIAS(M.first + "_r31", M.second.q[2]) &
        RK_SERIAL_LOAD_WITH_ALIAS(M.first + "_r12", M.second.q[4]) &
        RK_SERIAL_LOAD_WITH_ALIAS(M.first + "_r22", M.second.q[5]) &
        RK_SERIAL_LOAD_WITH_ALIAS(M.first + "_r32", M.second.q[6]) &
        RK_SERIAL_LOAD_WITH_ALIAS(M.first + "_r13", M.second.q[8]) &
        RK_SERIAL_LOAD_WITH_ALIAS(M.first + "_r23", M.second.q[9]) &
        RK_SERIAL_LOAD_WITH_ALIAS(M.first + "_r33", M.second.q[10]) &
        RK_SERIAL_LOAD_WITH_ALIAS(M.first + "_t_x", M.second.q[12]) &
        RK_SERIAL_LOAD_WITH_ALIAS(M.first + "_t_y", M.second.q[13]) &
        RK_SERIAL_LOAD_WITH_ALIAS(M.first + "_t_z", M.second.q[14]);
    M.second.q[3] = 0.0;
    M.second.q[7] = 0.0;
    M.second.q[11] = 0.0;
    M.second.q[15] = 1.0;
    return in;
  }

  /// Loading a trans_mat_3D value.
  friend serialization::iarchive& operator>>(serialization::iarchive& in,
                                             trans_mat_3D<T>& M) {
    return in & RK_SERIAL_LOAD_WITH_ALIAS("T", M);
  }

  /// Saving a trans_mat_3D value with a name.
  friend serialization::oarchive& operator&(
      serialization::oarchive& out,
      const std::pair<std::string, const trans_mat_3D<T>&>& M) {
    return out & RK_SERIAL_SAVE_WITH_ALIAS(M.first + "_r11", M.second.q[0]) &
           RK_SERIAL_SAVE_WITH_ALIAS(M.first + "_r21", M.second.q[1]) &
           RK_SERIAL_SAVE_WITH_ALIAS(M.first + "_r31", M.second.q[2]) &
           RK_SERIAL_SAVE_WITH_ALIAS(M.first + "_r12", M.second.q[4]) &
           RK_SERIAL_SAVE_WITH_ALIAS(M.first + "_r22", M.second.q[5]) &
           RK_SERIAL_SAVE_WITH_ALIAS(M.first + "_r32", M.second.q[6]) &
           RK_SERIAL_SAVE_WITH_ALIAS(M.first + "_r13", M.second.q[8]) &
           RK_SERIAL_SAVE_WITH_ALIAS(M.first + "_r23", M.second.q[9]) &
           RK_SERIAL_SAVE_WITH_ALIAS(M.first + "_r33", M.second.q[10]) &
           RK_SERIAL_SAVE_WITH_ALIAS(M.first + "_t_x", M.second.q[12]) &
           RK_SERIAL_SAVE_WITH_ALIAS(M.first + "_t_y", M.second.q[13]) &
           RK_SERIAL_SAVE_WITH_ALIAS(M.first + "_t_z", M.second.q[14]);
  }

  /// Saving a trans_mat_3D value.
  friend serialization::oarchive& operator<<(serialization::oarchive& out,
                                             const trans_mat_3D<T>& M) {
    return out & RK_SERIAL_SAVE_WITH_ALIAS("T", M);
  }
};

namespace rtti {

template <typename T>
struct get_type_id<trans_mat_3D<T>> {
  static constexpr unsigned int ID = 0x00000019;
  static constexpr auto type_name = std::string_view{"ReaK::trans_mat_3D"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }

  using save_type = const trans_mat_3D<T>&;
  using load_type = trans_mat_3D<T>&;
};

}  // namespace rtti

/// Prints a 3D transformation matrix to a standard output stream (<<)
/// as "(quaternion = (q0; q1; q2; q3); translation = (tx; ty; tz))".
template <class T>
std::ostream& operator<<(std::ostream& out_stream, const trans_mat_3D<T>& M) {
  out_stream << "(quaternion = " << M.getQuaternion()
             << "; translation = " << M.getTranslation() << ")";
  return out_stream;
}

template <typename T>
struct is_readable_matrix<trans_mat_3D<T>> {
  static constexpr bool value = true;
  using type = is_readable_matrix<trans_mat_3D<T>>;
};

}  // namespace ReaK

#endif  // REAK_MATH_KINETOSTATICS_ROTATIONS_3D_H_

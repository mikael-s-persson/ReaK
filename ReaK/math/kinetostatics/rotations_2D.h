/**
 * \file rotations_2D.h
 *
 * This library declares all geometric 2D rotation classes for fixed (2,3) and variable dimensions.
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

#ifndef REAK_MATH_KINETOSTATICS_ROTATIONS_2D_H_
#define REAK_MATH_KINETOSTATICS_ROTATIONS_2D_H_

#include "ReaK/math/lin_alg/mat_alg_skew_symmetric.h"
#include "ReaK/math/lin_alg/mat_alg_square.h"
#include "ReaK/math/lin_alg/mat_alg_symmetric.h"
#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/vect_alg.h"

#include <cassert>
#include <type_traits>

namespace ReaK {

// forward declaration, for friend declaration.
template <class T>
class trans_mat_2D;

/**
 * This class is a rotation matrix (proper orthogonal) of dimension 2 by 2.
 * \test All unit test for this class have been passed!
 */
template <typename T>
class rot_mat_2D {
 public:
  using self = rot_mat_2D<T>;

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

  static constexpr unsigned int static_row_count = 2;
  static constexpr unsigned int static_col_count = 2;
  static constexpr mat_alignment::tag alignment = mat_alignment::column_major;
  static constexpr mat_structure::tag structure = mat_structure::orthogonal;

 private:
  std::array<T, 2> q = {};

  rot_mat_2D(const_reference cos_a, const_reference sin_a) noexcept {
    q[0] = cos_a;
    q[1] = sin_a;
  }

 public:
  friend class trans_mat_2D<T>;

  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /**
   * Default Constructor with no rotation.
   * \test PASSED
   */
  rot_mat_2D() noexcept {
    q[0] = 1.0;
    q[1] = 0.0;
  }

  /**
   * Explicit constructor of a rotation matrix from a rotation angle.
   * \test PASSED
   */
  explicit rot_mat_2D(const_reference Angle) noexcept {
    using std::cos;
    using std::sin;
    q[0] = cos(Angle);
    q[1] = sin(Angle);
  }

  /**
   * Explicit constructor of a rotation matrix from a cosine and sine of an angle.
   * \test PASSED
   */
  explicit rot_mat_2D(vect<value_type, 2> v) noexcept {
    v = unit(v);
    q[0] = v[0];
    q[1] = v[1];
  }

  template <typename Matrix>
  explicit rot_mat_2D(const Matrix& R) {
    static_assert(is_readable_matrix_v<Matrix>);
    if ((R.get_col_count() != 2) || (R.get_row_count() != 2)) {
      throw std::range_error(
          "Right-hand-side of 2D rotation matrix assignment is not a 2x2 "
          "matrix!");
    }
    vect<value_type, 2> v = unit(vect<value_type, 2>(R(0, 0), R(1, 0)));
    q[0] = v[0];
    q[1] = v[1];
  }

  rot_mat_2D(const self& R) noexcept = default;

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Provides a copy of the rotation matrix as an ordinary 2x2 matrix.
   * \return this rotation matrix as a normal column-major matrix.
   * \test PASSED
   */
  mat<value_type, mat_structure::square> getMat() const {
    return {q[0], -q[1], q[1], q[0]};
  }

  /**
   * Returns the angle (-pi .. pi) of the rotation matrix.
   * \test PASSED
   */
  value_type getAngle() const noexcept {
    using std::atan2;
    return atan2(q[1], q[0]);
  }

  /**
   * Sets the angle (in radians) of the rotation matrix.
   * \test PASSED
   */
  void setAngle(const_reference Angle) noexcept {
    using std::cos;
    using std::sin;
    q[0] = cos(Angle);
    q[1] = sin(Angle);
  }

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  value_type operator[](int i) const {
    if (i >= 4) {
      throw std::range_error("Matrix index out of range.");
    }
    if (i == 3) {
      return q[0];
    }
    if (i == 2) {
      return -q[1];
    }
    return q[i];
  }

  /**
   * Array double-indexing operator, ith row and jth column, accessor for read only.
   * \test PASSED
   */
  value_type operator()(int i, int j) const {
    if ((i >= 2) || (j >= 2)) {
      throw std::range_error("Matrix index out of range.");
    }
    if ((j == 1) && (i == 0)) {
      return -q[1];
    }
    if ((j == 1) && (i == 1)) {
      return q[0];
    }
    return q[i];
  }

  int get_row_count() const noexcept { return 2; }
  int get_col_count() const noexcept { return 2; }

  /*******************************************************************************
                         Assignment Operators
*******************************************************************************/

  self& operator=(const self& R) noexcept = default;

  template <typename Matrix>
  self& operator=(const Matrix& R) {
    static_assert(is_readable_matrix_v<Matrix>);
    if ((R.get_col_count() != 2) || (R.get_row_count() != 2)) {
      throw std::range_error(
          "Right-hand-side of 2D rotation matrix assignment is not a 2x2 "
          "matrix!");
    }
    vect<value_type, 2> v = unit(vect<value_type, 2>(R(0, 0), R(1, 0)));
    q[0] = v[0];
    q[1] = v[1];
    return *this;
  }

  /**
   * Multiply-and-store operator.
   * \test PASSED
   */
  self& operator*=(const self& R) noexcept {
    value_type tmp = q[0] * R.q[0] - q[1] * R.q[1];
    q[1] = q[1] * R.q[0] + q[0] * R.q[1];
    q[0] = tmp;
    return *this;
  }

  template <typename Matrix>
  self& operator*=(const Matrix& R) {
    static_assert(is_readable_matrix_v<Matrix>);
    if ((R.get_col_count() != 2) || (R.get_row_count() != 2)) {
      throw std::range_error(
          "Right-hand-side of 2D rotation matrix assignment is not a 2x2 "
          "matrix!");
    }
    vect<value_type, 2> v = unit(vect<value_type, 2>(R(0, 0), R(1, 0)));
    value_type tmp = q[0] * v[0] - q[1] * v[1];
    q[1] = q[1] * v[0] + q[0] * v[1];
    q[0] = tmp;
    return *this;
  }

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /**
   * Rotation matrix multiplication.
   * \test PASSED
   */
  friend self operator*(const self& R1, const self& R2) noexcept {
    return self(R1.q[0] * R2.q[0] - R1.q[1] * R2.q[1],
                R1.q[1] * R2.q[0] + R1.q[0] * R2.q[1]);
  }

  /**
   * Matrix multiplication.
   * \test PASSED
   */
  template <typename Matrix>
  friend Matrix operator*(const self& R, const Matrix& M) {
    static_assert(is_fully_writable_matrix_v<Matrix>);
    if (M.get_row_count() != 2) {
      throw std::range_error(
          "Matrix M's row count is not 2, cannot perform 2D rotation!");
    }
    Matrix result(M);
    for (unsigned int jj = 0; jj < result.get_col_count(); ++jj) {
      result(0, jj) = R.q[0] * M(0, jj) - R.q[1] * M(1, jj);
      result(1, jj) = R.q[1] * M(0, jj) + R.q[0] * M(1, jj);
    }
    return result;
  }

  /**
   * 2D Rotation matrix times a column 2D vector.
   * \test PASSED
   */
  friend vect<value_type, 2> operator*(const self& R,
                                       const vect<value_type, 2>& V) noexcept {
    return {V[0] * R.q[0] - V[1] * R.q[1], V[0] * R.q[1] + V[1] * R.q[0]};
  }

  /**
   * Row 2D vector times a rotation matrix.
   * \test PASSED
   */
  friend vect<value_type, 2> operator*(const vect<value_type, 2>& V,
                                       const self& R) noexcept {
    return {V[0] * R.q[0] + V[1] * R.q[1], V[1] * R.q[0] - V[0] * R.q[1]};
  }

  /**
   * Matrix multiplication.
   * \test PASSED
   */
  template <typename Matrix>
  friend Matrix operator*(const Matrix& M, const self& R) {
    static_assert(is_fully_writable_matrix_v<Matrix>);
    if (M.get_col_count() != 2) {
      throw std::range_error(
          "Matrix M's column count is not 2, cannot perform 2D rotation!");
    }
    Matrix result(M);
    for (unsigned int i = 0; i < result.get_row_count(); ++i) {
      result(i, 0) = R.q[0] * M(i, 0) + R.q[1] * M(i, 1);
      result(i, 1) = R.q[0] * M(i, 1) - R.q[1] * M(i, 0);
    }
    return result;
  }

  /*******************************************************************************
                           Comparison Operators
  *******************************************************************************/

  /**
   * Equality Comparison operator.
   * \test PASSED
   */
  friend bool operator==(const self& R1, const self& R2) noexcept {
    return ((R1.q[0] == R2.q[0]) && (R1.q[1] == R2.q[1]));
  }

  /**
   * Inequality Comparison operator.
   * \test PASSED
   */
  friend bool operator!=(const self& R1, const self& R2) noexcept {
    return ((R1.q[0] != R2.q[0]) || (R1.q[1] != R2.q[1]));
  }

  /*******************************************************************************
                           Standard Matrix Methods
  *******************************************************************************/

  /**
   * Creates the transpose matrix.
   * \test PASSED
   */
  friend self transpose(const self& R) noexcept { return {R.q[0], -R.q[1]}; }

  /**
   * Creates the transpose matrix.
   * \test PASSED
   */
  friend self transpose_move(const self& R) noexcept {
    return {R.q[0], -R.q[1]};
  }

  /**
   * Gets the trace of the matrix.
   * \test PASSED
   */
  friend value_type trace(const self& R) noexcept {
    return value_type(2.0) * R.q[0];
  }

  /**
   * Gets the determinant of the matrix.
   * \test PASSED
   */
  friend value_type determinant(const self& /*unused*/) noexcept {
    return value_type(1.0);
  }

  /**
   * Gets the inverse of the matrix.
   * \test PASSED
   */
  friend self invert(const self& R) noexcept { return transpose(R); }

  /**
   * Gets the symmetric part of the matrix.
   * \test PASSED
   */
  mat<value_type, mat_structure::symmetric> getSymPart() const {
    return {q[0], value_type(0.0), q[0]};
  }

  /**
   * Gets the skew-symmetric part of the matrix.
   * \test PASSED
   */
  mat<value_type, mat_structure::skew_symmetric> getSkewSymPart() const {
    return mat<value_type, mat_structure::skew_symmetric>(-q[1]);
  }

  /// Loading a rot_mat_2D value with a name.
  friend serialization::iarchive& operator&(
      serialization::iarchive& in,
      const std::pair<std::string, rot_mat_2D<T>&>& R) {
    return in & RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_cos", R.second.q[0]) &
           RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_sin", R.second.q[1]);
  }

  /// Loading a rot_mat_2D value.
  friend serialization::iarchive& operator>>(serialization::iarchive& in,
                                             rot_mat_2D<T>& R) {
    return in & RK_SERIAL_LOAD_WITH_NAME(R);
  }

  /// Saving a rot_mat_2D value with a name.
  friend serialization::oarchive& operator&(
      serialization::oarchive& out,
      const std::pair<std::string, const rot_mat_2D<T>&>& R) {
    return out & RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_cos", R.second.q[0]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_sin", R.second.q[1]);
  }

  /// Saving a rot_mat_2D value.
  friend serialization::oarchive& operator<<(serialization::oarchive& out,
                                             const rot_mat_2D<T>& R) {
    return out & RK_SERIAL_SAVE_WITH_NAME(R);
  }
};

namespace rtti {

template <typename T>
struct get_type_id<rot_mat_2D<T>> {
  static constexpr unsigned int ID = 0x00000016;
  static constexpr auto type_name = std::string_view{"ReaK::rot_mat_2D"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }

  using save_type = const rot_mat_2D<T>&;
  using load_type = rot_mat_2D<T>&;
};

}  // namespace rtti

/**
 * Prints a 2D rotation matrix to a standard output stream (<<) as "(angle = a)".
 * \test PASSED
 */
template <typename T>
std::ostream& operator<<(std::ostream& out_stream, const rot_mat_2D<T>& R) {
  out_stream << "(angle = " << R.getAngle() << ")";
  return out_stream;
}

template <typename T>
struct is_readable_matrix<rot_mat_2D<T>> {
  static constexpr bool value = true;
  using type = is_readable_matrix<rot_mat_2D<T>>;
};

/**
 * This class is a transformation matrix 3 by 3, i.e. to rotate and translate a 2D vector.
 * \test All unit tests for this class have been passed!
 */
template <typename T>
class trans_mat_2D {
 public:
  using self = trans_mat_2D<T>;

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
  static constexpr mat_structure::tag structure = mat_structure::square;

  using translation_type = vect<value_type, 2>;

 private:
  std::array<value_type, 9> q = {};

  trans_mat_2D(const_reference a11, const_reference a12, const_reference a13,
               const_reference a21, const_reference a22,
               const_reference a23) noexcept {
    q[0] = a11;
    q[3] = a12;
    q[6] = a13;
    q[1] = a21;
    q[4] = a22;
    q[7] = a23;
    q[2] = 0.0;
    q[5] = 0.0;
    q[8] = 1.0;
  }

 public:
  /*******************************************************************************
                           Constructors / Destructors
  *******************************************************************************/

  /**
   * Default constructor.
   * \test PASSED
   */
  trans_mat_2D() noexcept {
    q[0] = 1.0;
    q[3] = 0.0;
    q[6] = 0.0;
    q[1] = 0.0;
    q[4] = 1.0;
    q[7] = 0.0;
    q[2] = 0.0;
    q[5] = 0.0;
    q[8] = 1.0;
  }

  /**
   * Constructor from a rotation angle and a translation vector.
   * \test PASSED
   */
  explicit trans_mat_2D(
      const_reference Angle,
      translation_type Translation = translation_type()) noexcept {
    q[4] = (q[0] = cos(Angle));
    q[3] = -(q[1] = sin(Angle));
    q[6] = Translation[0];
    q[7] = Translation[1];
    q[2] = 0.0;
    q[5] = 0.0;
    q[8] = 1.0;
  }

  explicit trans_mat_2D(
      const rot_mat_2D<value_type>& R,
      translation_type Translation = translation_type()) noexcept {
    q[4] = (q[0] = R(0, 0));
    q[3] = -(q[1] = R(1, 0));
    q[6] = Translation[0];
    q[7] = Translation[1];
    q[2] = 0.0;
    q[5] = 0.0;
    q[8] = 1.0;
  }

  template <typename Matrix>
  explicit trans_mat_2D(const Matrix& M) {
    static_assert(is_readable_matrix_v<Matrix>);
    if ((M.get_col_count() != 3) || (M.get_row_count() != 3)) {
      throw std::range_error(
          "Right-hand-side of 2D transformation matrix assignment is not a 3x3 "
          "matrix!");
    }
    translation_type v = unit(translation_type(M(0, 0), M(1, 0)));
    q[0] = v[0];
    q[1] = v[1];
    q[2] = 0.0;
    q[3] = -q[1];
    q[4] = q[0];
    q[5] = 0.0;
    q[6] = M(0, 2);
    q[7] = M(1, 2);
    q[8] = 1.0;
  }

  trans_mat_2D(const self&) noexcept = default;

  /*******************************************************************************
                           Accessors and Methods
  *******************************************************************************/

  /**
   * Provides a copy of the transformation matrix as an ordinary 3x3 matrix.
   * \test PASSED
   */
  mat<value_type, mat_structure::square> getMat() const {
    return {q[0], q[3], q[6], q[1], q[4], q[7], 0.0, 0.0, 1.0};
  }

  /**
   * Provides a copy of the rotation matrix part of the transformation matrix.
   * \test PASSED
   */
  rot_mat_2D<value_type> getRotMat() const noexcept { return {q[0], q[1]}; }

  /**
   * Sets the rotation part of the transformation matrix.
   * \test PASSED
   */
  void setRotMat(const rot_mat_2D<value_type>& R) noexcept {
    q[4] = (q[0] = R.q[0]);
    q[3] = -(q[1] = R.q[1]);
  }

  /**
   * Returns the angle of the rotation matrix.
   * \test PASSED
   */
  value_type getAngle() const noexcept { return atan2(q[1], q[0]); }

  /**
   * Returns the angle of the rotation matrix.
   * \test PASSED
   */
  void setAngle(const_reference Angle) noexcept {
    q[4] = (q[0] = cos(Angle));
    q[3] = -(q[1] = sin(Angle));
  }

  /**
   * Provides a copy of the translation part of the transformation matrix.
   * \test PASSED
   */
  translation_type getTranslation() const noexcept {
    return translation_type(q[6], q[7]);
  }

  /**
   * Sets the translation part of the transformation matrix.
   * \test PASSED
   */
  void setTranslation(const translation_type& Translation) noexcept {
    q[6] = Translation.q[0];
    q[7] = Translation.q[1];
  }

  /**
   * Array indexing operator, accessor for read only.
   * \test PASSED
   */
  const_reference operator[](int i) const {
    if (i >= 9) {
      throw std::range_error("Matrix index out of range.");
    }
    return q[i];
  }

  /**
   * Array double-indexing operator, ith row and jth column, accessor for read only.
   * \test PASSED
   */
  const_reference operator()(int i, int j) const {
    if ((i >= 3) || (j >= 3)) {
      throw std::range_error("Matrix index out of range.");
    }
    return q[j * 3 + i];
  }

  int get_row_count() const noexcept { return 3; }
  int get_col_count() const noexcept { return 3; }

  /*******************************************************************************
                         Assignment Operators
*******************************************************************************/

  self& operator=(const self&) noexcept = default;

  template <typename Matrix>
  self& operator=(const Matrix& M) {
    static_assert(is_readable_matrix_v<Matrix>);
    if ((M.get_col_count() != 3) || (M.get_row_count() != 3)) {
      throw std::range_error(
          "Right-hand-side of 2D transformation matrix assignment is not a 3x3 "
          "matrix!");
    }
    translation_type v = unit(translation_type(M(0, 0), M(1, 0)));
    q[0] = v[0];
    q[1] = v[1];
    q[2] = 0.0;
    q[3] = -q[1];
    q[4] = q[0];
    q[5] = 0.0;
    q[6] = M(0, 2);
    q[7] = M(1, 2);
    q[8] = 1.0;
  }

  /**
   * Multiply-and-store operator.
   * \test PASSED
   */
  self& operator*=(const self& M) noexcept { return (*this = (*this) * M); }

  template <typename Matrix>
  self& operator*=(const Matrix& M) {
    static_assert(is_readable_matrix_v<Matrix>);
    return (*this = (*this) * M);
  }

  /*******************************************************************************
                           Basic Operators
  *******************************************************************************/

  /**
   * Multiplication with rotation matrix.
   */
  friend self operator*(const rot_mat_2D<value_type>& R,
                        const self& M) noexcept {
    return self(R) * M;
  }

  /**
   * Multiplication with rotation matrix.
   */
  friend self operator*(const self& M,
                        const rot_mat_2D<value_type>& R) noexcept {
    return M * self(R);
  }

  /**
   * Matrix multiplication.
   * \test PASSED
   */
  friend self operator*(const self& M1, const self& M2) noexcept {
    return {M1.q[0] * M2.q[0] + M1.q[3] * M2.q[1],
            M1.q[0] * M2.q[3] + M1.q[3] * M2.q[4],
            M1.q[0] * M2.q[6] + M1.q[3] * M2.q[7] + M1.q[6],
            M1.q[1] * M2.q[0] + M1.q[4] * M2.q[1],
            M1.q[1] * M2.q[3] + M1.q[4] * M2.q[4],
            M1.q[1] * M2.q[6] + M1.q[4] * M2.q[7] + M1.q[7]};
  }

  /**
   * Matrix multiplication.
   * \test PASSED
   */
  template <typename Matrix>
  friend Matrix operator*(const self& M1, const Matrix& M2) {
    static_assert(is_fully_writable_matrix_v<Matrix>);
    if (M2.get_row_count() != 3) {
      throw std::range_error(
          "Matrix M's row count is not 3, 2D transformation impossible!");
    }
    Matrix result(M2);
    for (int i = 0; i < 3; ++i) {
      for (int jj = 0; jj < result.get_col_count(); ++jj) {
        result(i, jj) = 0;
        for (int j = 0; j < 3; ++j) {
          result(i, jj) += M1.q[j * 3 + i] * M2(j, jj);
        }
      }
    }
    return result;
  }

  /**
   * Matrix multiplication.
   * \test PASSED
   */
  template <typename Matrix>
  friend Matrix operator*(const Matrix& M1, const self& M2) {
    static_assert(is_fully_writable_matrix_v<Matrix>);
    if (M1.get_col_count() != 3) {
      throw std::range_error(
          "Matrix M1's column count is not 3, 2D transformation impossible!");
    }
    Matrix result(M1.get_row_count(), 3);
    for (int i = 0; i < result.get_row_count(); ++i) {
      for (int jj = 0; jj < 3; ++jj) {
        result(i, jj) = 0;
        for (int j = 0; j < 3; ++j) {
          result(i, jj) += M1(i, j) * M2.q[jj * 3 + j];
        }
      }
    }
    return result;
  }

  /**
   * 2D Transformation matrix times a column 2D vector.
   * \test PASSED
   */
  friend vect<value_type, 2> operator*(const self& M,
                                       const vect<value_type, 2>& V) noexcept {
    return {V[0] * M.q[0] + V[1] * M.q[3] + M.q[6],
            V[0] * M.q[1] + V[1] * M.q[4] + M.q[7]};
  }

  /**
   * 2D Transformation matrix times a column 2D augmented vector.
   * \test PASSED
   */
  friend vect<value_type, 3> operator*(const self& M,
                                       const vect<value_type, 3>& V) noexcept {
    return {V[0] * M.q[0] + V[1] * M.q[3] + V[2] * M.q[6],
            V[0] * M.q[1] + V[1] * M.q[4] + V[2] * M.q[7], V[2] * M.q[8]};
  }

  /*******************************************************************************
                           Comparison Operators
  *******************************************************************************/

  /**
   * Standard equality operator.
   * \test PASSED
   */
  friend bool operator==(const self& M1, const self& M2) noexcept {
    return ((M1.q[0] == M2.q[0]) && (M1.q[1] == M2.q[1]) &&
            (M1.q[6] == M2.q[6]) && (M1.q[7] == M2.q[7]));
  }

  /**
   * Standard inequality operator.
   * \test PASSED
   */
  friend bool operator!=(const self& M1, const self& M2) noexcept {
    return ((M1.q[0] != M2.q[0]) || (M1.q[1] != M2.q[1]) ||
            (M1.q[6] != M2.q[6]) || (M1.q[7] != M2.q[7]));
  }

  /*******************************************************************************
                           Special Methods
  *******************************************************************************/

  /**
   * Rotate-only a 2D vector.
   * \test PASSED
   */
  vect<value_type, 2> rotate(const vect<value_type, 2>& V) const noexcept {
    return {V[0] * q[0] + V[1] * q[3], V[0] * q[1] + V[1] * q[4]};
  }

  /*******************************************************************************
                           Standard Matrix Methods
  *******************************************************************************/

  /**
   * Creates the transpose matrix.
   * \note the matrix is no longer a transformation matrix.
   * \test PASSED
   */
  friend mat<value_type, mat_structure::square> transpose(
      const self& M) noexcept {
    return {M.q[0], M.q[1], 0.0, M.q[3], M.q[4], 0.0, M.q[6], M.q[7], 1.0};
  }

  /**
   * Gets the trace of the matrix.
   * \test PASSED
   */
  friend value_type trace(const self& M) noexcept {
    return M.q[0] + M.q[4] + value_type(1.0);
  }

  /**
   * Gets the determinant of the matrix.
   * \test PASSED
   */
  friend value_type determinant(const self& /*unused*/) noexcept {
    return value_type(1.0);
  }

  /**
   * Invert the transformation.
   * \test PASSED
   */
  friend self invert(const self& M) noexcept {
    return {M.q[0], M.q[1], -M.q[0] * M.q[6] - M.q[1] * M.q[7],
            M.q[3], M.q[4], -M.q[3] * M.q[6] - M.q[4] * M.q[7]};
  }

  /**
   * Gets the symmetric part of the matrix.
   * \test PASSED
   */
  mat<value_type, mat_structure::symmetric> getSymPart() const {
    return {q[0], value_type(0.0),        value_type(0.5) * q[6],
            q[0], value_type(0.5) * q[7], value_type(1.0)};
  }

  /**
   * Gets the skew-symmetric part of the matrix.
   * \test PASSED
   */
  mat<value_type, mat_structure::skew_symmetric> getSkewSymPart() const {
    return {-q[1], value_type(0.5) * q[6], value_type(0.5) * q[7]};
  }

  /// Loading a trans_mat_2D value with a name.
  friend serialization::iarchive& operator&(
      serialization::iarchive& in,
      const std::pair<std::string, trans_mat_2D<T>&>& R) {
    in& RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_cos", R.second.q[0]) &
        RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_sin", R.second.q[1]) &
        RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_t_x", R.second.q[6]) &
        RK_SERIAL_LOAD_WITH_ALIAS(R.first + "_t_y", R.second.q[7]);
    R.second.q[3] = -R.second.q[1];
    R.second.q[4] = R.second.q[0];
    R.second.q[2] = 0.0;
    R.second.q[5] = 0.0;
    R.second.q[8] = 1.0;
    return in;
  }

  /// Loading a trans_mat_2D value.
  friend serialization::iarchive& operator<<(serialization::iarchive& out,
                                             const trans_mat_2D<T>& R) {
    return out & RK_SERIAL_LOAD_WITH_NAME(R);
  }

  /// Saving a trans_mat_2D value with a name.
  friend serialization::oarchive& operator&(
      serialization::oarchive& out,
      const std::pair<std::string, const trans_mat_2D<T>&>& R) {
    return out & RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_cos", R.second.q[0]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_sin", R.second.q[1]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_t_x", R.second.q[6]) &
           RK_SERIAL_SAVE_WITH_ALIAS(R.first + "_t_y", R.second.q[7]);
  }

  /// Saving a trans_mat_2D value.
  friend serialization::oarchive& operator<<(serialization::oarchive& out,
                                             const trans_mat_2D<T>& R) {
    return out & RK_SERIAL_SAVE_WITH_NAME(R);
  }
};

namespace rtti {

template <typename T>
struct get_type_id<trans_mat_2D<T>> {
  static constexpr unsigned int ID = 0x00000017;
  static constexpr auto type_name = std::string_view{"ReaK::trans_mat_2D"};
  static construct_ptr CreatePtr() noexcept { return nullptr; }

  using save_type = const trans_mat_2D<T>&;
  using load_type = trans_mat_2D<T>&;
};

}  // namespace rtti

/** * Prints a 2D rotation matrix to a standard output stream (<<) as "(angle = a; translation = (tx; ty))".
 * \test PASSED
 */
template <class T>
std::ostream& operator<<(std::ostream& out_stream, const trans_mat_2D<T>& M) {
  out_stream << "(angle = " << M.getAngle()
             << "; translation = " << M.getTranslation() << ")";
  return out_stream;
}

template <typename T>
struct is_readable_matrix<trans_mat_2D<T>> {
  static constexpr bool value = true;
  using type = is_readable_matrix<trans_mat_2D<T>>;
};

}  // namespace ReaK

#endif  // REAK_MATH_KINETOSTATICS_ROTATIONS_2D_H_

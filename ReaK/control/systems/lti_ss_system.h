/**
 * \file lti_ss_system.h
 *
 * This library provides a class template which can be used to create a simple continuous-time LTI
 * state-space system, as used in ReaK::ctrl.
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

#ifndef REAK_CONTROL_SYSTEMS_LTI_SS_SYSTEM_H_
#define REAK_CONTROL_SYSTEMS_LTI_SS_SYSTEM_H_

#include "ReaK/core/base/named_object.h"
#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/vect_alg.h"

#include "ReaK/control/systems/linear_ss_system_concept.h"

namespace ReaK::ctrl {

/**
 * This class template can be used to create a simple continuous-time LTI
 * state-space system, as used in ReaK::ctrl. A continuous-time LTI state-space system is
 * basically described by four system matrices (A,B,C,D) which make the linear mapping
 * between the current state and input and the state-derivative and current output.
 * \tparam T The value-type of the system matrices.
 */
template <typename T>
class lti_system_ss : public named_object {
 private:
  mat<T, mat_structure::square> A;
  mat<T, mat_structure::rectangular> B;
  mat<T, mat_structure::rectangular> C;
  mat<T, mat_structure::rectangular> D;

 public:
  using self = lti_system_ss<T>;
  using value_type = T;
  using size_type = std::size_t;

  using point_type = vect_n<T>;
  using point_difference_type = vect_n<T>;
  using point_derivative_type = vect_n<T>;
  using topology = self;

  using time_type = T;
  using time_difference_type = T;

  using input_type = vect_n<T>;
  using output_type = vect_n<T>;

  using matrixA_type = mat<T, mat_structure::square>;
  using matrixB_type = mat<T, mat_structure::rectangular>;
  using matrixC_type = mat<T, mat_structure::rectangular>;
  using matrixD_type = mat<T, mat_structure::rectangular>;

  static constexpr std::size_t dimensions = 0;
  static constexpr std::size_t input_dimensions = 0;
  static constexpr std::size_t output_dimensions = 0;

  /**
   * Parametrized and default constructor.
   * \param aX_size The size of the state-vector.
   * \param aU_size The size of the input-vector.
   * \param aY_size The size of the output-vector.
   */
  lti_system_ss(size_type aX_size = 0, size_type aU_size = 0,
                size_type aY_size = 0, const std::string& aName = "")
      : A(aX_size),
        B(aX_size, aU_size),
        C(aY_size, aX_size),
        D(aY_size, aU_size) {
    setName(aName);
  }

  /**
   * Standard copy-constructor.
   */
  lti_system_ss(const self& rhs) : A(rhs.A), B(rhs.B), C(rhs.C), D(rhs.D) {
    setName(rhs.getName());
  }

  /**
   * Parametrized constructor.
   * \param aA The continuous-time system matrix A.
   * \param aB The continuous-time system matrix B.
   * \param aC The continuous-time system matrix C.
   * \param aD The continuous-time system matrix D.
   */
  template <typename MatrixA, typename MatrixB, typename MatrixC,
            typename MatrixD>
  lti_system_ss(const MatrixA& aA, const MatrixB& aB, const MatrixC& aC,
                const MatrixD& aD, const std::string& aName = "")
      : A(aA), B(aB), C(aC), D(aD) {
    static_assert(is_readable_matrix_v<MatrixA>);
    static_assert(is_readable_matrix_v<MatrixB>);
    static_assert(is_readable_matrix_v<MatrixC>);
    static_assert(is_readable_matrix_v<MatrixD>);
    setName(aName);
  }

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) throw() {
    using std::swap;
    swap(lhs.A, rhs.A);
    swap(lhs.B, rhs.B);
    swap(lhs.C, rhs.C);
    swap(lhs.D, rhs.D);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(self rhs) {
    swap(*this, rhs);
    return *this;
  }

  /**
   * Sets the continuous-time system matrix A of the continuous-time system.
   * \param aA The new continuous-time system matrix A for the system.
   */
  template <typename MatrixA>
  void setA(const MatrixA& aA) {
    static_assert(is_readable_matrix_v<MatrixA>);
    A = aA;
  }

  /**
   * Sets the continuous-time system matrix B of the continuous-time system.
   * \param aB The new continuous-time system matrix B for the system.
   */
  template <typename MatrixB>
  void setB(const MatrixB& aB) {
    static_assert(is_readable_matrix_v<MatrixB>);
    B = aB;
  }

  /**
   * Sets the continuous-time system matrix C of the continuous-time system.
   * \param aC The new continuous-time system matrix C for the system.
   */
  template <typename MatrixC>
  void setC(const MatrixC& aC) {
    static_assert(is_readable_matrix_v<MatrixC>);
    C = aC;
  }

  /**
   * Sets the continuous-time system matrix D of the continuous-time system.
   * \param aD The new continuous-time system matrix D for the system.
   */
  template <typename MatrixD>
  void setD(const MatrixD& aD) {
    static_assert(is_readable_matrix_v<MatrixD>);
    D = aD;
  }

  /**
   * Returns the continuous-time system matrix A of the system.
   * \return The continuous-time system matrix A of the system.
   */
  const matrixA_type& getA() const { return A; }
  /**
   * Returns the continuous-time system matrix B of the system.
   * \return The continuous-time system matrix B of the system.
   */
  const matrixB_type& getB() const { return B; }
  /**
   * Returns the continuous-time system matrix C of the system.
   * \return The continuous-time system matrix C of the system.
   */
  const matrixC_type& getC() const { return C; }
  /**
   * Returns the continuous-time system matrix D of the system.
   * \return The continuous-time system matrix D of the system.
   */
  const matrixD_type& getD() const { return D; }

  /**
   * Returns the state-vector's dimension.
   * \return The state-vector's dimension.
   */
  size_type get_state_dimensions() const { return A.get_col_count(); }
  /**
   * Returns the input-vector's dimension.
   * \return The input-vector's dimension.
   */
  size_type get_input_dimensions() const { return B.get_col_count(); }
  /**
   * Returns the output-vector's dimension.
   * \return The output-vector's dimension.
   */
  size_type get_output_dimensions() const { return C.get_row_count(); }

  /**
   * Fills the given matrices with the continuous-time system matrices.
   * \param aA Stores, as output, the system matrix A.
   * \param aB Stores, as output, the system matrix B.
   * \param aC Stores, as output, the system matrix C.
   * \param aD Stores, as output, the system matrix D.
   */
  template <typename MatrixA, typename MatrixB, typename MatrixC,
            typename MatrixD>
  void get_linear_blocks(MatrixA& aA, MatrixB& aB, MatrixC& aC,
                         MatrixD& aD) const {
    static_assert(is_writable_matrix_v<MatrixA>);
    static_assert(is_writable_matrix_v<MatrixB>);
    static_assert(is_writable_matrix_v<MatrixC>);
    static_assert(is_writable_matrix_v<MatrixD>);
    aA = A;
    aB = B;
    aC = C;
    aD = D;
  }

  /**
   * Fills the given matrices with the continuous-time system matrices.
   * \param aA Stores, as output, the system matrix A.
   * \param aB Stores, as output, the system matrix B.
   * \param aC Stores, as output, the system matrix C.
   * \param aD Stores, as output, the system matrix D.
   */
  template <typename MatrixA, typename MatrixB, typename MatrixC,
            typename MatrixD, typename StateSpaceType>
  void get_linear_blocks(MatrixA& aA, MatrixB& aB, MatrixC& aC, MatrixD& aD,
                         const StateSpaceType&, const time_type& = 0,
                         const point_type& = point_type(),
                         const input_type& = input_type()) const {
    static_assert(is_writable_matrix_v<MatrixA>);
    static_assert(is_writable_matrix_v<MatrixB>);
    static_assert(is_writable_matrix_v<MatrixC>);
    static_assert(is_writable_matrix_v<MatrixD>);
    aA = A;
    aB = B;
    aC = C;
    aD = D;
  }

  /**
   * Returns state-derivative given the current state, input and time.
   * \param p The current state.
   * \param u The current input.
   * \return The state-derivative.
   */
  template <typename StateSpaceType>
  point_derivative_type get_state_derivative(const StateSpaceType&,
                                             const point_type& p,
                                             const input_type& u,
                                             const time_type& = 0) const {
    return A * p + B * u;
  }

  /**
   * Returns output of the system given the current state, input and time.
   * \param p The current state.
   * \param u The current input.
   * \return The current output.
   */
  template <typename StateSpaceType>
  output_type get_output(const StateSpaceType&, const point_type& p,
                         const input_type& u, const time_type& = 0) const {
    return C * p + D * u;
  }

  virtual void save(ReaK::serialization::oarchive& aA, unsigned int) const {
    ReaK::named_object::save(
        aA, ReaK::named_object::getStaticObjectType()->TypeVersion());
    aA& RK_SERIAL_SAVE_WITH_NAME(A) & RK_SERIAL_SAVE_WITH_NAME(B) &
        RK_SERIAL_SAVE_WITH_NAME(C) & RK_SERIAL_SAVE_WITH_NAME(D);
  }
  virtual void load(ReaK::serialization::iarchive& aA, unsigned int) {
    ReaK::named_object::load(
        aA, ReaK::named_object::getStaticObjectType()->TypeVersion());
    aA& RK_SERIAL_LOAD_WITH_NAME(A) & RK_SERIAL_LOAD_WITH_NAME(B) &
        RK_SERIAL_LOAD_WITH_NAME(C) & RK_SERIAL_LOAD_WITH_NAME(D);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2300001, 1, "lti_system_ss",
                              named_object)
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_SYSTEMS_LTI_SS_SYSTEM_H_

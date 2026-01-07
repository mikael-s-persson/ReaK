/**
 * \file discretized_lti_sys.h
 *
 * This library provides a class template which can take a continuous-time linear time-invariant
 * state-space system and turn it into a discrete-time linear time-invariant state-space system.
 * The class template uses a matrix exponential method to compute the system matrix of the
 * discrete-time system, and thus, the resulting system is exact, as far as the matrix exponential
 * is correct (see mat_exp_methods.h).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2011
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

#ifndef REAK_CONTROL_SYSTEMS_DISCRETIZED_LTI_SYS_H_
#define REAK_CONTROL_SYSTEMS_DISCRETIZED_LTI_SYS_H_

#include "ReaK/core/base/named_object.h"
#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/mat_exp_methods.h"
#include "ReaK/math/lin_alg/mat_qr_decomp.h"

#include "ReaK/control/systems/discrete_linear_sss_concept.h"
#include "ReaK/control/systems/linear_ss_system_concept.h"

namespace ReaK::ctrl {

/**
 * This class template can take a continuous-time linear time-invariant
 * state-space system and turn it into a discrete-time linear time-invariant state-space system.
 * The class template uses a matrix exponential method to compute the system matrix of the
 * discrete-time system, and thus, the resulting system is exact, as far as the matrix exponential
 * is correct (see mat_exp_methods.hpp).
 *
 * \tparam LTISystem The linear time-invariance, continuous-time state-space system type, should model
 *LinearSSSystemConcept with LTISystemType.
 */
template <typename LTISystem>
class discretized_lti_sys : public named_object {
 public:
  using self = discretized_lti_sys<LTISystem>;
  using value_type = LTISystem::value_type;
  using size_type = LTISystem::size_type;

  using point_type = typename ss_system_traits<LTISystem>::point_type;
  using point_difference_type =
      typename ss_system_traits<LTISystem>::point_difference_type;

  using time_type = typename ss_system_traits<LTISystem>::time_type;
  using time_difference_type =
      typename ss_system_traits<LTISystem>::time_difference_type;

  using input_type = typename ss_system_traits<LTISystem>::input_type;
  using output_type = typename ss_system_traits<LTISystem>::output_type;

  using matrixA_type =
      typename linear_ss_system_traits<LTISystem>::matrixA_type;
  using matrixB_type =
      typename linear_ss_system_traits<LTISystem>::matrixB_type;
  using matrixC_type =
      typename linear_ss_system_traits<LTISystem>::matrixC_type;
  using matrixD_type =
      typename linear_ss_system_traits<LTISystem>::matrixD_type;

  static constexpr std::size_t dimensions =
      ss_system_traits<LTISystem>::dimensions;
  static constexpr std::size_t input_dimensions =
      ss_system_traits<LTISystem>::input_dimensions;
  static constexpr std::size_t output_dimensions =
      ss_system_traits<LTISystem>::output_dimensions;

  BOOST_CONCEPT_ASSERT((LinearSSSystemConcept<LTISystem, LTISystemType>));

 private:
  time_difference_type dt;

  matrixA_type Ad;
  matrixB_type Bd;
  matrixC_type Cd;
  matrixD_type Dd;

 public:
  /**
   * Returns the dimensions of the state vectors.
   * \return The dimensions of the state vectors.
   */
  size_type get_state_dimensions() const { return Ad.get_row_count(); }

  /**
   * Returns the dimensions of the input vectors.
   * \return The dimensions of the input vectors.
   */
  size_type get_input_dimensions() const { return Bd.get_col_count(); }

  /**
   * Returns the dimensions of the output vectors.
   * \return The dimensions of the output vectors.
   */
  size_type get_output_dimensions() const { return Cd.get_row_count(); }

  /**
   * Parametrized and default constructor.
   * \param aSys The continuous-time system.
   * \param aDt The time-step of the discrete-time system.
   */
  discretized_lti_sys(const LTISystem& aSys = LTISystem(),
                      const time_difference_type& aDt = 1)
      : dt(aDt) {
    setName(aSys.getName());

    aSys.get_linear_blocks(Ad, Bd, Cd, Dd);

    mat<value_type, mat_structure::square> A_aug(
        Ad.get_col_count() + Bd.get_col_count(), value_type(0));
    set_block(A_aug, Ad * dt, 0, 0);
    set_block(A_aug, Bd * dt, 0, Ad.get_col_count());
    mat<value_type, mat_structure::square> A_aug_exp(
        Ad.get_col_count() + Bd.get_col_count(), value_type(0));
    exp_PadeSAS(A_aug, A_aug_exp, QR_linlsqsolver());
    Ad = get_block(A_aug_exp, 0, 0, Ad.get_row_count(), Ad.get_col_count());
    Bd = get_block(A_aug_exp, 0, Ad.get_col_count(), Bd.get_row_count(),
                   Bd.get_col_count());
  }

  /**
   * Standard copy-constructor.
   */
  discretized_lti_sys(const self& rhs)
      : dt(rhs.dt), Ad(rhs.Ad), Bd(rhs.Bd), Cd(rhs.Cd), Dd(rhs.Dd) {
    setName(rhs.getName());
  }

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) throw() {
    using std::swap;
    swap(lhs.dt, rhs.dt);
    swap(lhs.Ad, rhs.Ad);
    swap(lhs.Bd, rhs.Bd);
    swap(lhs.Cd, rhs.Cd);
    swap(lhs.Dd, rhs.Dd);
  }

  /**
   * Standard assignment operator.
   */
  self& operator=(self rhs) {
    swap(*this, rhs);
    return *this;
  }

  /**
   * Fills the given matrices with the discrete-time system's state transition matrices.
   * \param aA Stores, as output, the system matrix A.
   * \param aB Stores, as output, the system matrix B.
   */
  template <WritableMatrix MatrixA, WritableMatrix MatrixB>
  void get_state_transition_blocks(MatrixA& aA, MatrixB& aB) const {
    aA = Ad;
    aB = Bd;
  }

  /**
   * Fills the given matrices with the discrete-time system's state transition matrices.
   * \param aA Stores, as output, the system matrix A.
   * \param aB Stores, as output, the system matrix B.
   */
  template <typename StateSpaceType, WritableMatrix MatrixA,
            WritableMatrix MatrixB>
  void get_state_transition_blocks(MatrixA& aA, MatrixB& aB,
                                   const StateSpaceType&,
                                   [[maybe_unused]] const time_type& t_0 = time_type(),
                                   [[maybe_unused]] const time_type& t_1 = time_type(),
                                   [[maybe_unused]] const point_type& p_0 = point_type(),
                                   [[maybe_unused]] const point_type& p_1 = point_type(),
                                   [[maybe_unused]] const input_type& u_0 = input_type(),
                                   [[maybe_unused]] const input_type& u_1 = input_type()) const {
    aA = Ad;
    aB = Bd;
  }

  /**
   * Fills the given matrices with the discrete-time system's output function matrices.
   * \param aC Stores, as output, the system matrix C.
   * \param aD Stores, as output, the system matrix D.
   */
  template <WritableMatrix MatrixC, WritableMatrix MatrixD>
  void get_output_function_blocks(MatrixC& aC, MatrixD& aD) const {
    aC = Cd;
    aD = Dd;
  }

  /**
   * Fills the given matrices with the discrete-time system's output function matrices.
   * \param aC Stores, as output, the system matrix C.
   * \param aD Stores, as output, the system matrix D.
   */
  template <typename StateSpaceType, WritableMatrix MatrixC,
            WritableMatrix MatrixD>
  void get_output_function_blocks(MatrixC& aC, MatrixD& aD,
                                  const StateSpaceType&,
                                  [[maybe_unused]] const time_type& t = time_type(),
                                  [[maybe_unused]] const point_type& p = point_type(),
                                  [[maybe_unused]] const input_type& u = input_type()) const {
    aC = Cd;
    aD = Dd;
  }

  /**
   * Returns the time-step of the discrete-time system.
   */
  time_difference_type get_time_step() const { return dt; }

  /**
   * Returns next state of the system given the current state, input and time.
   * \param p The current state.
   * \param u The current input.
   * \param t The current time.
   * \return The next state, at t + get_time_step().
   */
  template <typename StateSpaceType>
  point_type get_next_state(const StateSpaceType&, const point_type& p,
                            const input_type& u, [[maybe_unused]] const time_type& t = 0) const {
    using ReaK::from_vect;
    using ReaK::to_vect;
    return from_vect<point_type>(Ad * to_vect<value_type>(p) +
                                 Bd * to_vect<value_type>(u));
  }

  /**
   * Returns output of the system given the current state, input and time.
   * \param p The current state.
   * \param u The current input.
   * \param t The current time.
   * \return The current output.
   */
  template <typename StateSpaceType>
  output_type get_output(const StateSpaceType&, const point_type& p,
                         const input_type& u, [[maybe_unused]] const time_type& t = 0) const {
    return from_vect<output_type>(Cd * to_vect<value_type>(p) +
                                  Dd * to_vect<value_type>(u));
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void save(ReaK::serialization::oarchive& aA, unsigned int) const {
    ReaK::named_object::save(
        aA, ReaK::named_object::getStaticObjectType()->TypeVersion());
    aA& RK_SERIAL_SAVE_WITH_NAME(dt) & RK_SERIAL_SAVE_WITH_NAME(Ad) &
        RK_SERIAL_SAVE_WITH_NAME(Bd) & RK_SERIAL_SAVE_WITH_NAME(Cd) &
        RK_SERIAL_SAVE_WITH_NAME(Dd);
  }
  virtual void load(ReaK::serialization::iarchive& aA, unsigned int) {
    ReaK::named_object::load(
        aA, ReaK::named_object::getStaticObjectType()->TypeVersion());
    aA& RK_SERIAL_LOAD_WITH_NAME(dt) & RK_SERIAL_LOAD_WITH_NAME(Ad) &
        RK_SERIAL_LOAD_WITH_NAME(Bd) & RK_SERIAL_LOAD_WITH_NAME(Cd) &
        RK_SERIAL_LOAD_WITH_NAME(Dd);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2300005, 1, "discretized_lti_sys",
                              named_object)
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_SYSTEMS_DISCRETIZED_LTI_SYS_H_

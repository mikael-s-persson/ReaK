
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

#ifndef REAK_INVARIANT_KALMAN_BUCY_FILTER_HPP
#define REAK_INVARIANT_KALMAN_BUCY_FILTER_HPP

#include "ReaK/math/integrators/integrator.hpp"
#include "ReaK/math/lin_alg/mat_alg.hpp"
#include "ReaK/math/lin_alg/mat_cholesky.hpp"
#include "ReaK/math/lin_alg/vect_concepts.hpp"

#include "ReaK/control/estimators/belief_state_concept.hpp"
#include "ReaK/control/estimators/covariance_concept.hpp"
#include "ReaK/control/systems/invariant_system_concept.hpp"
#include "ReaK/control/systems/linear_ss_system_concept.hpp"

#include <type_traits>

namespace ReaK::ctrl {

namespace detail {

template <typename T, typename InvariantSystem, typename StateSpaceType,
          typename SystemNoiseCovariance, typename MeasurementNoiseCovariance>
struct invariant_kb_system : public state_rate_function<T> {

  using value_type = T;
  using size_type = std::size_t;
  using state_type = typename ss_system_traits<InvariantSystem>::point_type;
  using state_deriv_type =
      typename ss_system_traits<InvariantSystem>::point_derivative_type;
  using input_type = typename ss_system_traits<InvariantSystem>::input_type;
  using output_type = typename ss_system_traits<InvariantSystem>::output_type;

  using matrixA_type =
      typename linear_ss_system_traits<InvariantSystem>::matrixA_type;
  using matrixB_type =
      typename linear_ss_system_traits<InvariantSystem>::matrixB_type;
  using matrixC_type =
      typename linear_ss_system_traits<InvariantSystem>::matrixC_type;
  using matrixD_type =
      typename linear_ss_system_traits<InvariantSystem>::matrixD_type;

  using invariant_error_type =
      typename invariant_system_traits<InvariantSystem>::invariant_error_type;
  using invariant_correction_type = typename invariant_system_traits<
      InvariantSystem>::invariant_correction_type;

  BOOST_CONCEPT_ASSERT(
      (InvariantContinuousSystemConcept<InvariantSystem, StateSpaceType>));
  BOOST_CONCEPT_ASSERT(
      (CovarianceMatrixConcept<SystemNoiseCovariance, input_type>));
  BOOST_CONCEPT_ASSERT(
      (CovarianceMatrixConcept<MeasurementNoiseCovariance, output_type>));

  const InvariantSystem& sys;
  const StateSpaceType& state_space;
  const input_type& u;
  const output_type& z;

  mat<value_type, mat_structure::symmetric> Q;
  mat<value_type, mat_structure::symmetric> R_inv;
  mat<value_type, mat_structure::rectangular> K;
  mat<value_type, mat_structure::square> P;

  matrixA_type A;
  matrixB_type B;
  matrixC_type C;
  matrixD_type D;

  invariant_kb_system(const InvariantSystem& aSys,
                      const StateSpaceType& aStateSpace, const input_type& aU,
                      const output_type& aZ, const SystemNoiseCovariance& aQ,
                      const MeasurementNoiseCovariance& aR)
      : sys(aSys), state_space(aStateSpace), u(aU), z(aZ), Q(aQ.get_matrix()) {
    invert_Cholesky(aR.get_matrix(), R_inv);
    K.set_col_count(R_inv.get_row_count());
    K.set_row_count(Q.get_row_count());
    P.set_row_count(Q.get_row_count());
  }

  virtual void computeStateRate(double aTime,
                                const ReaK::vect_n<value_type>& aState,
                                ReaK::vect_n<value_type>& aStateRate) {
    using ReaK::from_vect;
    using ReaK::to_vect;

    ReaK::vect_n<value_type> x;
    x.resize(aState.size() - Q.get_row_count() * Q.get_row_count());
    for (int i = 0; i < x.size(); ++i) {
      x[i] = aState[i];
    }

    state_type x_state = from_vect<state_type>(x);
    sys.get_linear_blocks(A, B, C, D, state_space, aTime, x_state, u);
    invariant_error_type e =
        sys.get_invariant_error(state_space, x_state, u, z, aTime);

    for (int j = 0; j < Q.get_row_count(); ++j) {
      for (int i = 0; i < Q.get_row_count(); ++i) {
        P(i, j) = aState[x.size() + Q.get_row_count() * j + i];
      }
    }

    K = P * transpose_view(C) * R_inv;

    vect_n<value_type> xd = to_vect<value_type>(sys.apply_correction(
        state_space, x_state,
        sys.get_state_derivative(state_space, x_state, u, aTime),
        from_vect<invariant_correction_type>(K * e), u, aTime));
    P = (A - K * C) * P + B * Q + Q * transpose_view(B) + P * transpose_view(A);

    for (int i = 0; i < xd.size(); ++i) {
      aStateRate[i] = xd[i];
    }

    for (int j = 0; j < Q.get_row_count(); ++j) {
      for (int i = 0; i < Q.get_row_count(); ++i) {
        aStateRate[xd.size() + Q.get_row_count() * j + i] =
            0.5 * (P(i, j) + P(j, i));
      }
    }
  }
};
}  // namespace detail

template <typename InvariantSystem, typename StateSpaceType,
          typename BeliefState, typename InputBelief,
          typename MeasurementBelief>
void invariant_kalman_bucy_filter_step(
    const InvariantSystem& sys,
    integrator<mat_value_type_t<
        typename covariance_mat_traits<typename continuous_belief_state_traits<
            BeliefState>::covariance_type>::matrix_type>>& integ,
    const StateSpaceType& state_space, BeliefState& b_x, const InputBelief& b_u,
    const MeasurementBelief& b_z,
    typename ss_system_traits<InvariantSystem>::time_difference_type dt,
    typename ss_system_traits<InvariantSystem>::time_type t = 0) {
  // here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) prediction
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) prediction
  BOOST_CONCEPT_ASSERT(
      (InvariantContinuousSystemConcept<InvariantSystem, StateSpaceType>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<InputBelief>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<MeasurementBelief>));
  static_assert(is_continuous_belief_state_v<BeliefState>);
  static_assert(belief_state_traits<BeliefState>::representation ==
                belief_representation::gaussian);
  static_assert(belief_state_traits<BeliefState>::distribution ==
                belief_distribution::unimodal);

  using StateType = typename ss_system_traits<InvariantSystem>::point_type;
  using CovType =
      typename continuous_belief_state_traits<BeliefState>::covariance_type;
  using MatType = typename covariance_mat_traits<CovType>::matrix_type;
  using ValueType = mat_value_type_t<MatType>;

  using InputCovType =
      typename continuous_belief_state_traits<InputBelief>::covariance_type;
  using OutputCovType = typename continuous_belief_state_traits<
      MeasurementBelief>::covariance_type;

  using ReaK::from_vect;
  using ReaK::to_vect;

  integ.setTime(t);
  integ.clearStateVector();
  vect_n<ValueType> x = to_vect<ValueType>(b_x.get_mean_state());
  integ.addStateElements(x);
  mat<ValueType, mat_structure::square> P =
      mat<ValueType, mat_structure::square>(b_x.get_covariance().get_matrix());

  for (int j = 0; j < P.get_col_count(); ++j) {
    for (int i = 0; i < P.get_row_count(); ++i) {
      integ.addStateElement(P(i, j));
    }
  }

  auto integ_sys = std::make_shared<detail::invariant_kb_system<
      ValueType, InvariantSystem, StateSpaceType, InputCovType, OutputCovType>>(
      sys, state_space, b_u.get_mean_state(), b_z.get_mean_state(),
      b_u.get_covariance(), b_z.get_covariance());

  integ.setStateRateFunc(integ_sys);

  integ.integrate(t + dt);

  auto it = integ.getStateBegin();
  for (int i = 0; i < x.size(); ++it, ++i) {
    x[i] = *it;
  }

  for (int j = 0; j < P.get_col_count(); ++j) {
    for (int i = 0; i < P.get_row_count(); ++i, ++it) {
      P(i, j) = *it;
    }
  }

  b_x.set_mean_state(from_vect<StateType>(x));
  b_x.set_covariance(CovType(MatType(P)));

  integ.setStateRateFunc(std::shared_ptr<state_rate_function<ValueType>>());
}

}  // namespace ReaK::ctrl

#endif

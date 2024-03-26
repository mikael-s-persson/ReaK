
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

#ifndef REAK_CONTROL_CONTROLLERS_HYBRID_KALMAN_FILTER_H_
#define REAK_CONTROL_CONTROLLERS_HYBRID_KALMAN_FILTER_H_

#include "ReaK/math/integrators/integrator.h"
#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_cholesky.h"
#include "ReaK/math/lin_alg/vect_concepts.h"

#include "ReaK/control/estimators/belief_state_concept.h"
#include "ReaK/control/estimators/covariance_concept.h"
#include "ReaK/control/systems/linear_ss_system_concept.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"

#include <type_traits>

namespace ReaK::ctrl {

namespace detail {

template <typename T, typename LinearSystem, typename StateSpaceType,
          CovarianceMatrix<pp::topology_point_difference_type_t<StateSpaceType>>
              SystemNoiseCovariance>
requires LinearSSSystem<LinearSystem, StateSpaceType,
                        LinearizedSystemType> struct kalman_bucy_predictor
    : public state_rate_function<T> {

  using value_type = T;
  using size_type = std::size_t;
  using state_type = typename ss_system_traits<LinearSystem>::point_type;
  using input_type = typename ss_system_traits<LinearSystem>::input_type;
  using output_type = typename ss_system_traits<LinearSystem>::output_type;

  using matrixA_type =
      typename linear_ss_system_traits<LinearSystem>::matrixA_type;
  using matrixB_type =
      typename linear_ss_system_traits<LinearSystem>::matrixB_type;
  using matrixC_type =
      typename linear_ss_system_traits<LinearSystem>::matrixC_type;
  using matrixD_type =
      typename linear_ss_system_traits<LinearSystem>::matrixD_type;

  const LinearSystem& sys;
  const StateSpaceType& state_space;
  const input_type& u;

  mat<value_type, mat_structure::symmetric> Q;
  mat<value_type, mat_structure::square> P;

  matrixA_type A;
  matrixB_type B;
  matrixC_type C;
  matrixD_type D;

  kalman_bucy_system(const LinearSystem& aSys,
                     const StateSpaceType& aStateSpace, const input_type& aU,
                     const SystemNoiseCovariance& aQ)
      : sys(aSys), state_space(aStateSpace), u(aU), Q(aQ.get_matrix()) {
    P.set_row_count(Q.get_row_count());
  }

  virtual void computeStateRate(double aTime,
                                const ReaK::vect_n<value_type>& aState,
                                ReaK::vect_n<value_type>& aStateRate) {
    state_type x;
    x.resize(Q.get_row_count());
    for (size_type i = 0; i < x.size(); ++i) {
      x[i] = aState[i];
    }

    for (size_type j = 0; j < x.size(); ++j) {
      for (size_type i = 0; i < x.size(); ++i) {
        P(i, j) = aState[x.size() * (j + 1) + i];
      }
    }

    sys.get_linear_blocks(A, B, C, D, state_space, aTime, x, u);

    x = A * x + B * u;
    P = A * P;
    P += B * Q + Q * transpose_view(B) + transpose_view(P);

    for (size_type i = 0; i < x.size(); ++i) {
      aStateRate[i] = x[i];
    }

    for (size_type j = 0; j < x.size(); ++j) {
      for (size_type i = 0; i < x.size(); ++i) {
        aStateRate[x.size() * (j + 1) + i] = 0.5 * (P(i, j) + P(j, i));
      }
    }
  }
};

}  // namespace detail

template <pp::Topology StateSpaceType,
          LinearSSSystem<StateSpaceType, LinearizedSystemType> LinearSystem,
          ContinuousBeliefState BState, ContinuousBeliefState InputBelief,
          ContinuousBeliefState MeasurementBelief>
void hybrid_kalman_filter_step(
    const LinearSystem& sys,
    integrator<mat_value_type_t<
        typename covariance_mat_traits<typename continuous_belief_state_traits<
            BState>::covariance_type>::matrix_type>>& integ,
    const StateSpaceType& state_space, BState& b_x, const InputBelief& b_u,
    const MeasurementBelief& b_z,
    typename ss_system_traits<LinearSystem>::time_difference_type dt,
    typename ss_system_traits<LinearSystem>::time_type t = 0) {
  static_assert(belief_state_traits<BState>::representation ==
                belief_representation::gaussian);
  static_assert(belief_state_traits<BState>::distribution ==
                belief_distribution::unimodal);

  using StateType = typename ss_system_traits<LinearSystem>::point_type;
  using StateDiffType =
      typename ss_system_traits<LinearSystem>::point_difference_type;
  using CovType =
      typename continuous_belief_state_traits<BState>::covariance_type;
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
  ReaK::vect_n<ValueType> x = to_vect<ValueType>(b_x.get_mean_state());
  integ.addStateElements(x);
  mat<ValueType, mat_structure::square> P = b_x.get_covariance().get_matrix();

  for (int j = 0; j < P.get_col_count(); ++j) {
    for (int i = 0; i < P.get_row_count(); ++i) {
      integ.addStateElement(P(i, j));
    }
  }

  integ.setStateRateFunc(
      std::make_shared<
          detail::kalman_bucy_system<ValueType, LinearSystem, StateSpaceType,
                                     InputCovType, OutputCovType>>(
          sys, state_space, to_vect<ValueType>(b_u.get_mean_state()),
          to_vect<ValueType>(b_z.get_mean_state()), b_u.get_covariance(),
          b_z.get_covariance()));

  integ.integrate(t + dt);

  std::vector<ValueType>::const_iterator it = integ.getStateBegin();
  for (int i = 0; i < x.size(); ++it, ++i) {
    x[i] = *it;
  }

  for (int j = 0; j < P.get_col_count(); ++j) {
    for (int i = 0; i < P.get_row_count(); ++i, ++it) {
      P(i, j) = *it;
    }
  }

  typename linear_ss_system_traits<LinearSystem>::matrixA_type A;
  typename linear_ss_system_traits<LinearSystem>::matrixB_type B;
  typename linear_ss_system_traits<LinearSystem>::matrixC_type C;
  typename linear_ss_system_traits<LinearSystem>::matrixD_type D;

  sys.get_linear_blocks(A, B, C, D, state_space, t, from_vect<StateType>(x),
                        b_u.get_mean_state());

  vect_n<ValueType> y = to_vect<ValueType>(b_z.get_mean_state()) - C * x -
                        D * to_vect<ValueType>(b_u.get_mean_state());
  mat<ValueType, mat_structure::rectangular, mat_alignment::column_major> CP =
      C * P;
  mat<ValueType, mat_structure::symmetric> S =
      CP * transpose_view(C) + b_z.get_covariance().get_matrix();
  linsolve_Cholesky(S, CP);
  mat<ValueType, mat_structure::rectangular> K = transpose_view(CP);

  b_x.set_mean_state(state_space.adjust(from_vect<StateType>(x),
                                        from_vect<StateDiffType>(K * y)));
  b_x.set_covariance(CovType(MatType(
      (mat<ValueType, mat_structure::identity>(K.get_row_count()) - K * C) *
      P)));

  integ.setStateRateFunc(std::shared_ptr<state_rate_function<ValueType>>());
}

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_CONTROLLERS_HYBRID_KALMAN_FILTER_H_


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

#ifndef REAK_CONTROL_CONTROLLERS_INVARIANT_AGGREGATE_KALMAN_FILTER_H_
#define REAK_CONTROL_CONTROLLERS_INVARIANT_AGGREGATE_KALMAN_FILTER_H_

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_cholesky.h"
#include "ReaK/math/lin_alg/mat_star_product.h"
#include "ReaK/math/lin_alg/vect_concepts.h"

#include "ReaK/control/estimators/belief_state_concept.h"
#include "ReaK/control/estimators/covariance_concept.h"
#include "ReaK/control/systems/discrete_linear_sss_concept.h"
#include "ReaK/control/systems/invariant_system_concept.h"

#include <type_traits>

namespace ReaK::ctrl {

template <InvariantDiscreteSystem ISystem, ContinuousBeliefState BState,
          typename SystemNoiseCovariance, typename MeasurementNoiseCovariance>
requires CovarianceMatrix<typename discrete_sss_traits<ISystem>::input_type>&&
    CovarianceMatrix<typename discrete_sss_traits<ISystem>::output_type> void
    invariant_aggregate_kf_step(
        const ISystem& sys, BState& b,
        const discrete_sss_traits<ISystem>::input_type& b_u,
        const discrete_sss_traits<ISystem>::output_type& b_z,
        hamiltonian_mat_t<mat_value_type_t<typename covariance_mat_traits<
            typename continuous_belief_state_traits<BState>::covariance_type>::
                                               matrix_type>>& ScSm,
        hamiltonian_mat_t<mat_value_type_t<typename covariance_mat_traits<
            typename continuous_belief_state_traits<BState>::covariance_type>::
                                               matrix_type>>& Sc,
        typename discrete_sss_traits<ISystem>::time_type t = 0) {
  // here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) update
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) update
  using StateType = typename discrete_sss_traits<ISystem>::point_type;
  using InputType = typename discrete_sss_traits<ISystem>::input_type;
  using OutputType = typename discrete_sss_traits<ISystem>::output_type;
  using CovType =
      typename continuous_belief_state_traits<BState>::covariance_type;
  using MatType = typename covariance_mat_traits<CovType>::matrix_type;
  using ValueType = mat_value_type_t<MatType>;

  static_assert(belief_state_traits<BState>::representation ==
                belief_representation::gaussian);
  static_assert(belief_state_traits<BState>::distribution ==
                belief_distribution::unimodal);

  typename discrete_linear_sss_traits<ISystem>::matrixA_type A;
  typename discrete_linear_sss_traits<ISystem>::matrixB_type B;
  typename discrete_linear_sss_traits<ISystem>::matrixC_type C;
  typename discrete_linear_sss_traits<ISystem>::matrixD_type D;

  using HamilMat = typename hamiltonian_mat<ValueType>::type;
  using HamilMatUp = typename hamiltonian_mat<ValueType>::upper;
  using HamilMatLo = typename hamiltonian_mat<ValueType>::lower;
  using HamilMatUL = typename hamiltonian_mat<ValueType>::upper_left;
  using HamilMatUR = typename hamiltonian_mat<ValueType>::upper_right;
  using HamilMatLL = typename hamiltonian_mat<ValueType>::lower_left;
  using HamilMatLR = typename hamiltonian_mat<ValueType>::lower_right;

  using InvFrameType =
      typename invariant_system_traits<ISystem>::invariant_frame_type;
  using InvErrorType =
      typename invariant_system_traits<ISystem>::invariant_error_type;
  using InvCorrType = typename invariant_system_traits<
      InvariantDiscreteSystem>::invariant_correction_type;

  StateType x = b.get_mean_state();
  MatType P = b.get_covariance().get_matrix();
  sys.get_linear_blocks(A, B, C, D, t, x, b_u.get_mean_state());
  int N = A.get_col_count();

  x = sys.get_next_state(x, u, t);
  P = (A * P * transpose_view(A)) + b_u.get_covariance().get_matrix();

  InvErrorType e = sys.get_output_error(
      x, b_u.get_mean_state(), b_z.get_mean_state(), t + sys.get_time_step());
  InvFrameType W = sys.get_invariant_prior_frame(
      b.get_mean_state(), x, b_u.get_mean_state(), t + sys.get_time_step());

  mat<ValueType, mat_structure::rectangular, mat_alignment::column_major> CP =
      C * P;
  mat<ValueType, mat_structure::symmetric> S =
      CP * transpose_view(C) + b_z.get_covariance().get_matrix();
  linsolve_Cholesky(S, CP);
  mat<ValueType, mat_structure::rectangular, mat_alignment::row_major> K =
      transpose_view(CP);

  b.set_mean_state(sys.apply_correction(x, from_vect<InvCorrType>(W * K * e),
                                        b_u.get_mean_state(),
                                        t + sys.get_time_step()));
  W = sys.get_invariant_posterior_frame(x_prior, b.get_mean_state(),
                                        b_u.get_mean_state(),
                                        t + sys.get_time_step()) *
      W;
  InvFrameType Wt = InvFrameType(transpose_view(W));
  b.set_covariance(CovType(MatType(
      W * (mat<ValueType, mat_structure::identity>(K.get_row_count()) - K * C) *
      P * transpose_view(W))));

  // TODO Apply the W transform somehow.

  HamilMat Sc_tmp(
      HamilMatUp(HamilMatUL(A), HamilMatUR(b_u.get_covariance().get_matrix())),
      HamilMatLo(HamilMatLL(mat<ValueType, mat_structure::nil>(N)),
                 HamilMatLR(transpose_view(A))));

  swap(Sc, Sc_tmp);
  HamilMat ScSm_tmp(star_product(
      Sc,
      HamilMatUp(HamilMatUL(mat<ValueType, mat_structure::identity>(N)),
                 HamilMatUR(mat<ValueType, mat_structure::nil>(N))),
      HamilMatLo(HamilMatLL(transpose_view(C) *
                            b_z.get_covariance().get_inverse_matrix() * C),
                 HamilMatLR(mat<ValueType, mat_structure::identity>(N)))));
  swap(ScSm, ScSm_tmp);
}

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_CONTROLLERS_INVARIANT_AGGREGATE_KALMAN_FILTER_H_

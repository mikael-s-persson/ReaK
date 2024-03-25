/**
 * \file unscented_kalman_filter.h
 *
 * This library contains an attempt at implementing the unscented Kalman Filter (UKF), however
 * tests have shown problems with the implementation and since it is not useful for anything
 * right now, it has been left in this malfunctioning state.
 *
 * \todo Revise and fix this implementation.
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

#ifndef REAK_CONTROL_CONTROLLERS_UNSCENTED_KALMAN_FILTER_H_
#define REAK_CONTROL_CONTROLLERS_UNSCENTED_KALMAN_FILTER_H_

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_cholesky.h"
#include "ReaK/math/lin_alg/mat_svd_method.h"
#include "ReaK/math/lin_alg/vect_concepts.h"

#include "ReaK/control/estimators/belief_state_concept.h"
#include "ReaK/control/estimators/covariance_concept.h"
#include "ReaK/control/systems/discrete_sss_concept.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"

#include <type_traits>

namespace ReaK::ctrl {

/** UKF is not usable or even working at all! And does not consider Q as input-noise. */
template <pp::Topology StateSpaceType, DiscreteSSS<StateSpaceType> System,
          ContinuousBeliefState BState, ContinuousBeliefState InputBelief>
void unscented_kalman_predict(
    const System& sys, const StateSpaceType& state_space, BState& b_x,
    const InputBelief& b_u,
    typename discrete_sss_traits<System>::time_type t = 0,
    typename belief_state_traits<BState>::scalar_type alpha = 1E-3,
    typename belief_state_traits<BState>::scalar_type kappa = 1,
    typename belief_state_traits<BState>::scalar_type beta = 2) {
  static_assert(belief_state_traits<BState>::representation ==
                belief_representation::gaussian);
  static_assert(belief_state_traits<BState>::distribution ==
                belief_distribution::unimodal);

  using std::sqrt;

  using StateType = typename discrete_sss_traits<System>::point_type;
  using CovType =
      typename continuous_belief_state_traits<BState>::covariance_type;
  using MatType = typename covariance_mat_traits<CovType>::matrix_type;
  using ValueType = mat_value_type_t<MatType>;
  using SizeType = typename vect_n<ValueType>::size_type;

  using InputType =
      typename continuous_belief_state_traits<InputBelief>::state_type;
  using InputCovType =
      typename continuous_belief_state_traits<InputBelief>::covariance_type;
  using InputMatType =
      typename covariance_mat_traits<InputCovType>::matrix_type;

  StateType x = b_x.get_mean_state();

  MatType P = b_x.get_covariance().get_matrix();
  const InputMatType& Q = b_u.get_covariance().get_matrix();

  SizeType N = P.get_row_count();
  SizeType M = Q.get_row_count();

  mat<ValueType, mat_structure::square> L_p(N + M);
  mat<ValueType, mat_structure::square> P_aug(N + M);
  sub(P_aug)(range(0, N), range(0, N)) = P;
  sub(P_aug)(range(N, N + M), range(N, N + M)) = Q;

  try {
    decompose_Cholesky(P_aug, L_p);
  } catch (singularity_error&) {
    // use SVD instead.
    mat<ValueType, mat_structure::square> svd_U, svd_V;
    mat<ValueType, mat_structure::diagonal> svd_E;
    decompose_SVD(P_aug, svd_U, svd_E, svd_V);
    if (svd_E(0, 0) < 0)
      throw singularity_error(
          "'A-Priori Covariance P, in UKF prediction is singular, beyond "
          "repair!'");
    ValueType min_tolerable_sigma = sqrt(svd_E(0, 0)) * 1E-2;
    for (unsigned int i = 0; i < svd_E.get_row_count(); ++i) {
      if (svd_E(i, i) < min_tolerable_sigma * min_tolerable_sigma) {
        svd_E(i, i) = min_tolerable_sigma;
      } else {
        svd_E(i, i) = sqrt(svd_E(i, i));
      }
    }
    L_p = svd_U * svd_E;
  }

  ValueType lambda = alpha * alpha * (N + M + kappa) - N - M;
  ValueType gamma = sqrt(ValueType(N + M) + lambda);

  vect_n<StateType> X_a(1 + 2 * (N + M));
  X_a[0] = sys.get_next_state(state_space, x, b_u.get_mean_state(), t);
  for (SizeType j = 0; j < N + M; ++j) {
    StateType x_right = x;
    StateType x_left = x;
    for (SizeType i = 0; i < N; ++i) {
      x_right[i] += gamma * L_p(i, j);
      x_left[i] -= gamma * L_p(i, j);
    }
    InputType u_right = b_u.get_mean_state();
    InputType u_left = b_u.get_mean_state();
    for (SizeType i = 0; i < M; ++i) {
      u_right[i] += gamma * L_p(N + i, j);
      u_left[i] -= gamma * L_p(N + i, j);
    }
    X_a[1 + 2 * j] = sys.get_next_state(state_space, x_right, u_right, t);
    X_a[2 + 2 * j] = sys.get_next_state(state_space, x_left, u_left, t);
  }

  gamma = ValueType(1) / (ValueType(N + M) + lambda);
  x = (lambda * gamma) * X_a[0];
  for (SizeType j = 0; j < N + M; ++j) {
    x += (0.5 * gamma) * (X_a[1 + 2 * j] + X_a[2 + 2 * j]);
  }
  for (SizeType j = 0; j < 1 + 2 * (N + M); ++j) {
    X_a[j] -= x;
  }

  ValueType W_c = (lambda * gamma + ValueType(1) - alpha * alpha + beta);
  for (SizeType i = 0; i < N; ++i) {
    for (SizeType j = i; j < N; ++j) {
      P(i, j) = W_c * X_a[0][i] * X_a[0][j];
    }
  }

  W_c = ValueType(0.5) * gamma;
  for (SizeType k = 1; k < 1 + 2 * (N + M); ++k) {
    for (SizeType i = 0; i < N; ++i) {
      for (SizeType j = i; j < N; ++j) {
        P(i, j) += W_c * X_a[k][i] * X_a[k][j];
      }
    }
  }

  b_x.set_mean_state(x);
  b_x.set_covariance(CovType(P));
}

/** UKF is not usable or even working at all! And does not consider Q as input-noise. */
template <pp::Topology StateSpaceType, DiscreteSSS<StateSpaceType> System,
          ContinuousBeliefState BState, ContinuousBeliefState InputBelief,
          ContinuousBeliefState MeasurementBelief>
void unscented_kalman_update(
    const System& sys, const StateSpaceType& state_space, BState& b_x,
    const InputBelief& b_u, const MeasurementBelief& b_z,
    typename discrete_sss_traits<System>::time_type t = 0,
    typename belief_state_traits<BState>::scalar_type alpha = 1E-3,
    typename belief_state_traits<BState>::scalar_type kappa = 1,
    typename belief_state_traits<BState>::scalar_type beta = 2) {
  static_assert(belief_state_traits<BState>::representation ==
                belief_representation::gaussian);
  static_assert(belief_state_traits<BState>::distribution ==
                belief_distribution::unimodal);

  using std::sqrt;

  using StateType = typename discrete_sss_traits<System>::point_type;
  using OutputType = typename discrete_sss_traits<System>::output_type;
  using CovType =
      typename continuous_belief_state_traits<BState>::covariance_type;
  using MatType = typename covariance_mat_traits<CovType>::matrix_type;
  using ValueType = mat_value_type_t<MatType>;
  using SizeType = typename vect_n<ValueType>::size_type;

  using OutputCovType = typename continuous_belief_state_traits<
      MeasurementBelief>::covariance_type;
  using OutputMatType =
      typename covariance_mat_traits<OutputCovType>::matrix_type;

  StateType x = b_x.get_mean_state();

  const MatType& P = b_x.get_covariance().get_matrix();
  const OutputMatType& R = b_z.get_covariance().get_matrix();

  SizeType N = P.get_row_count();
  SizeType M = R.get_row_count();

  mat<ValueType, mat_structure::square> L_p(N + M);
  mat<ValueType, mat_structure::square> P_aug(N + M);
  sub(P_aug)(range(0, N), range(0, N)) = P;
  sub(P_aug)(range(N, N + M), range(N, N + M)) = R;

  try {
    decompose_Cholesky(P_aug, L_p);
  } catch (singularity_error&) {
    // use SVD instead.
    mat<ValueType, mat_structure::square> svd_U, svd_V;
    mat<ValueType, mat_structure::diagonal> svd_E;
    decompose_SVD(P_aug, svd_U, svd_E, svd_V);
    if (svd_E(0, 0) < 0) {
      throw singularity_error(
          "'A-Priori Covariance P, in UKF prediction is singular, beyond "
          "repair!'");
    }
    ValueType min_tolerable_sigma = sqrt(svd_E(0, 0)) * 1E-2;
    for (unsigned int i = 0; i < svd_E.get_row_count(); ++i) {
      if (svd_E(i, i) < min_tolerable_sigma * min_tolerable_sigma) {
        svd_E(i, i) = min_tolerable_sigma;
      } else {
        svd_E(i, i) = sqrt(svd_E(i, i));
      }
    }
    L_p = svd_U * svd_E;
  }

  ValueType lambda = alpha * alpha * (N + M + kappa) - N - M;
  ValueType gamma = sqrt(ValueType(N + M) + lambda);

  vect_n<StateType> X_a(1 + 2 * (N + M));
  vect_n<OutputType> Y_a(1 + 2 * (N + M));
  X_a[0] = x;
  Y_a[0] = sys.get_output(state_space, x, b_u.get_mean_state(), t);
  for (SizeType j = 0; j < N + M; ++j) {
    X_a[1 + 2 * j] = x;
    X_a[2 + 2 * j] = x;
    for (SizeType i = 0; i < N; ++i) {
      X_a[1 + 2 * j][i] += gamma * L_p(i, j);
      X_a[2 + 2 * j][i] -= gamma * L_p(i, j);
    }
    OutputType z_right = b_z.get_mean_state();
    OutputType z_left = b_z.get_mean_state();
    for (SizeType i = 0; i < M; ++i) {
      z_right[i] = gamma * L_p(N + i, j);
      z_left[i] = -gamma * L_p(N + i, j);
    }
    Y_a[1 + 2 * j] =
        sys.get_output(state_space, X_a[1 + 2 * j], b_u.get_mean_state(), t) +
        z_right;
    Y_a[2 + 2 * j] =
        sys.get_output(state_space, X_a[2 + 2 * j], b_u.get_mean_state(), t) +
        z_left;
  }

  gamma = ValueType(1) / (ValueType(N + M) + lambda);
  OutputType z_p = (lambda * gamma) * Y_a[0];
  for (SizeType j = 0; j < N + M; ++j) {
    z_p += (0.5 * gamma) * (Y_a[1 + 2 * j] + Y_a[2 + 2 * j]);
  }
  for (SizeType j = 0; j < 1 + 2 * (N + M); ++j) {
    Y_a[j] -= z_p;
  }

  mat<ValueType, mat_structure::symmetric> P_zz(z_p.size());
  ValueType W_c = (lambda * gamma + ValueType(1) - alpha * alpha + beta);
  for (SizeType i = 0; i < M; ++i) {
    for (SizeType j = i; j < M; ++j) {
      P_zz(i, j) = W_c * Y_a[0][i] * Y_a[0][j];
    }
  }

  W_c = ValueType(0.5) * gamma;
  for (SizeType k = 1; k < 1 + 2 * (N + M); ++k) {
    for (SizeType i = 0; i < M; ++i) {
      for (SizeType j = i; j < M; ++j) {
        P_zz(i, j) += W_c * Y_a[k][i] * Y_a[k][j];
      }
    }
  }

  mat<ValueType, mat_structure::rectangular> P_xz_t(z_p.size(), N);
  for (SizeType k = 1; k < 1 + 2 * (N + M); ++k) {
    for (SizeType i = 0; i < N; ++i) {
      for (SizeType j = i; j < M; ++j) {
        P_xz_t(j, i) += W_c * (X_a[k][i] - x[i]) * Y_a[k][j];
      }
    }
  }

  mat<ValueType, mat_structure::rectangular> Kt(P_xz_t);

  try {
    linsolve_Cholesky(P_zz, Kt);
  } catch (singularity_error&) {
    // use SVD instead.
    mat<ValueType, mat_structure::square> Pzz_pinv(P_zz.get_row_count());
    pseudoinvert_SVD(P_zz, Pzz_pinv);
    Kt = Pzz_pinv * Kt;
    throw singularity_error(
        "'A-Posteriori Measurement Covariance Pzz, in UKF update'");
  }

  b_x.set_mean_state(state_space.adjust(x, (b_z.get_mean_state() - z_p) * Kt));
  b_x.set_covariance(CovType(MatType(P - transpose_view(Kt) * P_xz_t)));
}

/** UKF is not usable or even working at all! And does not consider Q as input-noise. */
template <pp::Topology StateSpaceType, DiscreteSSS<StateSpaceType> System,
          ContinuousBeliefState BState, ContinuousBeliefState InputBelief,
          ContinuousBeliefState MeasurementBelief>
void unscented_kalman_filter_step(
    const System& sys, const StateSpaceType& state_space, BState& b_x,
    const InputBelief& b_u, const MeasurementBelief& b_z,
    typename discrete_sss_traits<System>::time_type t = 0,
    typename belief_state_traits<BState>::scalar_type alpha = 1E-3,
    typename belief_state_traits<BState>::scalar_type kappa = 1,
    typename belief_state_traits<BState>::scalar_type beta = 2) {
  static_assert(belief_state_traits<BState>::representation ==
                belief_representation::gaussian);
  static_assert(belief_state_traits<BState>::distribution ==
                belief_distribution::unimodal);
  unscented_kalman_predict(sys, state_space, b_x, b_u, t, alpha, kappa, beta);
  unscented_kalman_update(sys, state_space, b_x, b_u, b_z, t, alpha, kappa,
                          beta);
}

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_CONTROLLERS_UNSCENTED_KALMAN_FILTER_H_

/**
 * \file invariant_kalman_filter.h
 *
 * This library provides a number of functions and classes to do state estimation
 * using the Invariant Kalman Filter. This filtering technique applies to a
 * gaussian belief state and an invariant state-space system. The Kalman filter is
 * an optimal filter for a linear
 * state-space system where all sources of noise or disturbances are Gaussian
 * (normally distributed). This Kalman filter implementation requires
 * that the system be a linearized system which has an invariant frame which
 * can map the state-space into a frame in which the non-linearities have less or
 * no effect on the covariance transitions during the prediction and measurement
 * update.
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

#ifndef REAK_CONTROL_CONTROLLERS_INVARIANT_KALMAN_FILTER_H_
#define REAK_CONTROL_CONTROLLERS_INVARIANT_KALMAN_FILTER_H_

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_cholesky.h"
#include "ReaK/math/lin_alg/vect_concepts.h"

#include "ReaK/topologies/spaces/metric_space_concept.h"

#include "ReaK/control/estimators/belief_state_concept.h"
#include "ReaK/control/estimators/covariance_concept.h"
#include "ReaK/control/estimators/covariance_matrix.h"
#include "ReaK/control/estimators/gaussian_belief_state.h"
#include "ReaK/control/estimators/kalman_filter.h"
#include "ReaK/control/systems/discrete_linear_sss_concept.h"
#include "ReaK/control/systems/invariant_system_concept.h"

#include <type_traits>
#include <utility>

namespace ReaK::ctrl {

/**
 * This function template performs one prediction step using the Invariant Kalman Filter method.
 * \tparam ISystem An invariant discrete-time state-space system.
 * \tparam StateSpaceType A topology type on which the state-vectors can reside.
 * \tparam BState A belief state with a unimodular gaussian representation.
 * \tparam InputBelief A belief state with a unimodular gaussian representation.
 * \param sys The invariant discrete-time state-space system used in the state estimation.
 * \param state_space The state-space topology on which the state representations lie.
 * \param b_x As input, it stores the belief-state before the prediction step. As output, it stores
 *        the belief-state after the prediction step.
 * \param b_u The input belief to apply to the state-space system to make the transition of the
 *        mean-state, i.e., the current input vector and its covariance.
 * \param t The current time (before the prediction).
 *
 */
template <pp::Topology StateSpaceType,
          InvariantDiscreteSystem<StateSpaceType> ISystem,
          ContinuousBeliefState BState, ContinuousBeliefState InputBelief>
void invariant_kalman_predict(
    const ISystem& sys, const StateSpaceType& state_space, BState& b_x,
    const InputBelief& b_u,
    typename discrete_sss_traits<ISystem>::time_type t = 0) {
  static_assert(belief_state_traits<BState>::representation ==
                belief_representation::gaussian);
  static_assert(belief_state_traits<BState>::distribution ==
                belief_distribution::unimodal);

  using StateType = typename discrete_sss_traits<ISystem>::point_type;
  using CovType =
      typename continuous_belief_state_traits<BState>::covariance_type;
  using MatType = typename covariance_mat_traits<CovType>::matrix_type;
  using InvarFrame =
      typename invariant_system_traits<ISystem>::invariant_frame_type;

  typename discrete_linear_sss_traits<ISystem>::matrixA_type A;
  typename discrete_linear_sss_traits<ISystem>::matrixB_type B;

  StateType x = b_x.get_mean_state();
  MatType P = b_x.get_covariance().get_matrix();

  StateType x_prior =
      sys.get_next_state(state_space, x, b_u.get_mean_state(), t);
  sys.get_state_transition_blocks(A, B, state_space, t, t + sys.get_time_step(),
                                  x, x_prior, b_u.get_mean_state(),
                                  b_u.get_mean_state());
  InvarFrame W = sys.get_invariant_prior_frame(
      state_space, x, x_prior, b_u.get_mean_state(), t + sys.get_time_step());
  P = W *
      ((A * P * transpose_view(A)) +
       B * b_u.get_covariance().get_matrix() * transpose_view(B)) *
      transpose_view(W);
  b_x.set_mean_state(x_prior);
  b_x.set_covariance(CovType(P));
}

template <
    pp::Topology StateSpaceType,
    DiscreteLinearSSS<StateSpaceType, DiscreteLinearizedSystemType> ISystem,
    ContinuousBeliefState BState, ContinuousBeliefState InputBelief>
void invariant_kalman_predict(
    const ISystem& sys, const StateSpaceType& state_space, BState& b_x,
    const InputBelief& b_u,
    typename discrete_sss_traits<ISystem>::time_type t = 0) {
  kalman_predict(sys, state_space, b_x, b_u, t);
}

/**
 * This function template performs one measurement update step using the Invariant Kalman Filter method.
 * \tparam ISystem An invariant discrete-time state-space system.
 * \tparam StateSpaceType A topology type on which the state-vectors can reside.
 * \tparam BState A belief state with a unimodular gaussian representation.
 * \tparam InputBelief A belief state with a unimodular gaussian representation.
 * \tparam MeasurementBelief A belief state with a unimodular gaussian representation.
 * \param sys The invariant discrete-time state-space system used in the state estimation.
 * \param state_space The state-space topology on which the state representations lie.
 * \param b_x As input, it stores the belief-state before the update step. As output, it stores
 *        the belief-state after the update step.
 * \param b_u The input vector to apply to the state-space system to make the transition of the
 *        mean-state, i.e., the current input vector and its covariance.
 * \param b_z The output belief that was measured, i.e. the measurement vector and its covariance.
 * \param t The current time.
 *
 */
template <pp::Topology StateSpaceType,
          InvariantDiscreteSystem<StateSpaceType> ISystem,
          ContinuousBeliefState BState, ContinuousBeliefState InputBelief,
          ContinuousBeliefState MeasurementBelief>
void invariant_kalman_update(
    const ISystem& sys, const StateSpaceType& state_space, BState& b_x,
    const InputBelief& b_u, const MeasurementBelief& b_z,
    typename discrete_sss_traits<ISystem>::time_type t = 0) {
  static_assert(belief_state_traits<BState>::representation ==
                belief_representation::gaussian);
  static_assert(belief_state_traits<BState>::distribution ==
                belief_distribution::unimodal);

  using StateType = typename discrete_sss_traits<ISystem>::point_type;
  using CovType =
      typename continuous_belief_state_traits<BState>::covariance_type;
  using MatType = typename covariance_mat_traits<CovType>::matrix_type;
  using ValueType = mat_value_type_t<MatType>;
  using InvarFrame =
      typename invariant_system_traits<ISystem>::invariant_frame_type;
  using InvarCorr =
      typename invariant_system_traits<ISystem>::invariant_correction_type;

  typename discrete_linear_sss_traits<ISystem>::matrixC_type C;
  typename discrete_linear_sss_traits<ISystem>::matrixD_type D;

  StateType x = b_x.get_mean_state();
  MatType P = b_x.get_covariance().get_matrix();
  sys.get_output_function_blocks(C, D, state_space, t, x, b_u.get_mean_state());

  vect_n<ValueType> e = to_vect<ValueType>(
      sys.get_invariant_error(state_space, x, b_u.get_mean_state(),
                              b_z.get_mean_state(), t + sys.get_time_step()));

  mat<ValueType, mat_structure::rectangular, mat_alignment::column_major> CP =
      C * P;
  mat<ValueType, mat_structure::symmetric> S(CP * transpose_view(C) +
                                             b_z.get_covariance().get_matrix());
  linsolve_Cholesky(S, CP);
  mat<ValueType, mat_structure::rectangular, mat_alignment::row_major> K(
      transpose_view(CP));

  b_x.set_mean_state(
      sys.apply_correction(state_space, x, from_vect<InvarCorr>(K * e),
                           b_u.get_mean_state(), t + sys.get_time_step()));
  InvarFrame W = sys.get_invariant_posterior_frame(
      state_space, x, b_x.get_mean_state(), b_u.get_mean_state(),
      t + sys.get_time_step());
  b_x.set_covariance(CovType(MatType(
      W *
      ((mat<ValueType, mat_structure::identity>(K.get_row_count()) - K * C) *
       P) *
      transpose_view(W))));
}
template <
    pp::Topology StateSpaceType,
    DiscreteLinearSSS<StateSpaceType, DiscreteLinearizedSystemType> ISystem,
    ContinuousBeliefState BState, ContinuousBeliefState InputBelief,
    ContinuousBeliefState MeasurementBelief>
void invariant_kalman_update(
    const ISystem& sys, const StateSpaceType& state_space, BState& b_x,
    const InputBelief& b_u, const MeasurementBelief& b_z,
    typename discrete_sss_traits<ISystem>::time_type t = 0) {
  kalman_update(sys, state_space, b_x, b_u, b_z, t);
}

/**
 * This function template performs one complete estimation step using the (Extended) Kalman
 * Filter method, which includes a prediction and measurement update step. This function is,
 * in general, more efficient than applying the prediction and update separately.
 * \tparam ISystem An invariant discrete-time state-space system.
 * \tparam StateSpaceType A topology type on which the state-vectors can reside.
 * \tparam BState A belief state with a unimodular gaussian representation.
 * \tparam InputBelief A belief state with a unimodular gaussian representation.
 * \tparam MeasurementBelief A belief state with a unimodular gaussian representation.
 * \param sys The invariant discrete-time state-space system used in the state estimation.
 * \param state_space The state-space topology on which the state representations lie.
 * \param b_x As input, it stores the belief-state before the estimation step. As output, it stores
 *        the belief-state after the estimation step.
 * \param b_u The input vector to apply to the state-space system to make the transition of the
 *        mean-state, i.e., the current input vector and its covariance.
 * \param b_z The output belief that was measured, i.e. the measurement vector and its covariance.
 * \param t The current time (before the prediction).
 *
 */
template <pp::Topology StateSpaceType,
          InvariantDiscreteSystem<StateSpaceType> ISystem,
          ContinuousBeliefState BState, ContinuousBeliefState InputBelief,
          ContinuousBeliefState MeasurementBelief>
void invariant_kalman_filter_step(
    const ISystem& sys, const StateSpaceType& state_space, BState& b_x,
    const InputBelief& b_u, const MeasurementBelief& b_z,
    typename discrete_sss_traits<ISystem>::time_type t = 0) {
  static_assert(belief_state_traits<BState>::representation ==
                belief_representation::gaussian);
  static_assert(belief_state_traits<BState>::distribution ==
                belief_distribution::unimodal);

  using StateType = typename discrete_sss_traits<ISystem>::point_type;
  using CovType =
      typename continuous_belief_state_traits<BState>::covariance_type;
  using MatType = typename covariance_mat_traits<CovType>::matrix_type;
  using ValueType = mat_value_type_t<MatType>;
  using InvarFrame =
      typename invariant_system_traits<ISystem>::invariant_frame_type;
  using InvarCorr =
      typename invariant_system_traits<ISystem>::invariant_correction_type;

  typename discrete_linear_sss_traits<ISystem>::matrixA_type A;
  typename discrete_linear_sss_traits<ISystem>::matrixB_type B;
  typename discrete_linear_sss_traits<ISystem>::matrixC_type C;
  typename discrete_linear_sss_traits<ISystem>::matrixD_type D;

  StateType x = b_x.get_mean_state();
  MatType P = b_x.get_covariance().get_matrix();

  StateType x_prior =
      sys.get_next_state(state_space, x, b_u.get_mean_state(), t);
  sys.get_state_transition_blocks(A, B, state_space, t, t + sys.get_time_step(),
                                  x, x_prior, b_u.get_mean_state(),
                                  b_u.get_mean_state());
  InvarFrame W = sys.get_invariant_prior_frame(
      state_space, x, x_prior, b_u.get_mean_state(), t + sys.get_time_step());
  P = W *
      ((A * P * transpose_view(A)) +
       B * b_u.get_covariance().get_matrix() * transpose_view(B)) *
      transpose_view(W);

  sys.get_output_function_blocks(C, D, state_space, t + sys.get_time_step(),
                                 x_prior, b_u.get_mean_state());
  vect_n<ValueType> e = to_vect<ValueType>(
      sys.get_invariant_error(state_space, x_prior, b_u.get_mean_state(),
                              b_z.get_mean_state(), t + sys.get_time_step()));

  mat<ValueType, mat_structure::rectangular, mat_alignment::column_major> CP =
      C * P;
  mat<ValueType, mat_structure::symmetric> S(CP * transpose_view(C) +
                                             b_z.get_covariance().get_matrix());
  linsolve_Cholesky(S, CP);
  mat<ValueType, mat_structure::rectangular, mat_alignment::row_major> K(
      transpose_view(CP));

  b_x.set_mean_state(
      sys.apply_correction(state_space, x_prior, from_vect<InvarCorr>(K * e),
                           b_u.get_mean_state(), t + sys.get_time_step()));
  W = sys.get_invariant_posterior_frame(
      state_space, x_prior, b_x.get_mean_state(), b_u.get_mean_state(),
      t + sys.get_time_step());
  b_x.set_covariance(CovType(MatType(
      W *
      ((mat<ValueType, mat_structure::identity>(K.get_row_count()) - K * C) *
       P) *
      transpose_view(W))));
}

template <
    pp::Topology StateSpaceType,
    DiscreteLinearSSS<StateSpaceType, DiscreteLinearizedSystemType> ISystem,
    ContinuousBeliefState BState, ContinuousBeliefState InputBelief,
    ContinuousBeliefState MeasurementBelief>
void invariant_kalman_filter_step(
    const ISystem& sys, const StateSpaceType& state_space, BState& b_x,
    const InputBelief& b_u, const MeasurementBelief& b_z,
    typename discrete_sss_traits<ISystem>::time_type t = 0) {
  kalman_filter_step(sys, state_space, b_x, b_u, b_z, t);
}

/**
 * This class template can be used as a belief-state predictor (and transfer) that uses the
 * Invariant (Extended) Kalman Filter method. This class template models the BeliefTransferConcept and
 * the BeliefPredictorConcept.
 * \tparam IKFTransferFactory The factory type which can create this kalman predictor.
 */
template <typename IKFTransferFactory>
struct IKF_belief_transfer {
  using self = IKF_belief_transfer<IKFTransferFactory>;
  using state_space_system = typename IKFTransferFactory::state_space_system;
  using state_space_system_ptr = std::shared_ptr<state_space_system>;
  using time_type = typename discrete_sss_traits<state_space_system>::time_type;
  using time_difference_type =
      typename discrete_sss_traits<state_space_system>::time_difference_type;

  using covariance_type = covariance_matrix<vect_n<double>>;

  using input_type =
      typename discrete_sss_traits<state_space_system>::input_type;
  using output_type =
      typename discrete_sss_traits<state_space_system>::output_type;

  using input_belief_type = gaussian_belief_state<input_type, covariance_type>;
  using output_belief_type =
      gaussian_belief_state<output_type, covariance_type>;

  const IKFTransferFactory* factory;

  /**
   * Parametrized constructor.
   * \param aFactory A pointer to the factory object that is creating this object.
   */
  explicit IKF_belief_transfer(const IKFTransferFactory* aFactory)
      : factory(aFactory) {}

  IKF_belief_transfer() : IKF_belief_transfer(nullptr) {}

  /**
   * Returns the time-step of the predictor.
   * \return The time-step of the predictor.
   */
  time_difference_type get_time_step() const {
    return factory->get_time_step();
  }

  /**
   * Returns a reference to the underlying state-space system.
   * \return A reference to the underlying state-space system.
   */
  const state_space_system_ptr& get_ss_system() const {
    return factory->get_state_space_system();
  }

  /**
   * Returns the belief-state at the next time instant.
   * \tparam BeliefSpace The belief-space type on which to operate.
   * \param b_space The belief-space on which the belief-states lie.
   * \param b The current belief-state.
   * \param t The current time.
   * \param u The current input given to the system.
   * \param y The output that was measured at the next time instant.
   * \return the belief-state at the next time instant.
   */
  template <typename BeliefSpace>
  pp::topology_point_type_t<BeliefSpace> get_next_belief(
      const BeliefSpace& b_space, pp::topology_point_type_t<BeliefSpace> b,
      const time_type& t, const input_type& u, const output_type& y) const {
    invariant_kalman_filter_step(
        *(factory->get_state_space_system()), b_space.get_state_topology(), b,
        input_belief_type(
            u, covariance_type(factory->get_input_disturbance_cov())),
        output_belief_type(
            y, covariance_type(factory->get_measurement_noise_cov())),
        t);
    return b;
  }

  /**
   * Returns the prediction belief-state at the next time instant.
   * \tparam BeliefSpace The belief-space type on which to operate.
   * \param b_space The belief-space on which the belief-states lie.
   * \param b The current belief-state.
   * \param t The current time.
   * \param u The current input given to the system.
   * \return the belief-state at the next time instant, predicted by the filter.
   */
  template <typename BeliefSpace>
  pp::topology_point_type_t<BeliefSpace> predict_belief(
      const BeliefSpace& b_space, pp::topology_point_type_t<BeliefSpace> b,
      const time_type& t, const input_type& u) const {
    invariant_kalman_predict(
        *(factory->get_state_space_system()), b_space.get_state_topology(), b,
        input_belief_type(
            u, covariance_type(factory->get_input_disturbance_cov())),
        t);
    return b;
  }

  /**
   * Converts a prediction belief-state into an updated belief-state which assumes the most likely measurement.
   * \tparam BeliefSpace The belief-space type on which to operate.
   * \param b_space The belief-space on which the belief-states lie.
   * \param b The current prediction's belief-state.
   * \param t The current time.
   * \param u The current input given to the system.
   * \return the updated belief-state when assuming the most likely measurement.
   */
  template <typename BeliefSpace>
  pp::topology_point_type_t<BeliefSpace> prediction_to_ML_belief(
      const BeliefSpace& b_space, pp::topology_point_type_t<BeliefSpace> b,
      const time_type& t, const input_type& u) const {
    invariant_kalman_update(
        *(factory->get_state_space_system()), b_space.get_state_topology(), b,
        input_belief_type(
            u, covariance_type(factory->get_input_disturbance_cov())),
        output_belief_type(
            factory->get_state_space_system()->get_output(
                b_space.get_state_topology(), b.get_mean_state(), u, t),
            covariance_type(factory->get_measurement_noise_cov())),
        t);
    return b;
  }

  /**
   * Returns the prediction belief-state at the next time instant, assuming the upcoming measurement to be the most
   * likely one.
   * \tparam BeliefSpace The belief-space type on which to operate.
   * \param b_space The belief-space on which the belief-states lie.
   * \param b The current belief-state.
   * \param t The current time.
   * \param u The current input given to the system.
   * \return the belief-state at the next time instant, predicted by the filter.
   */
  template <typename BeliefSpace>
  pp::topology_point_type_t<BeliefSpace> predict_ML_belief(
      const BeliefSpace& b_space, pp::topology_point_type_t<BeliefSpace> b,
      const time_type& t, const input_type& u) const {
    input_belief_type b_u(
        u, covariance_type(factory->get_input_disturbance_cov()));
    invariant_kalman_predict(*(factory->get_state_space_system()),
                             b_space.get_state_topology(), b, b_u, t);
    invariant_kalman_update(
        *(factory->get_state_space_system()), b_space.get_state_topology(), b,
        b_u,
        output_belief_type(
            factory->get_state_space_system()->get_output(
                b_space.get_state_topology(), b.get_mean_state(), u, t),
            covariance_type(factory->get_measurement_noise_cov())),
        t + factory->get_state_space_system()->get_time_step());
    return b;
  }
};

/**
 * This class is a factory class for invariant Kalman filtering predictors on a belief-space.
 * \tparam InvariantSystem An invariant discrete-time state-space system modeling the InvariantDiscreteSystemConcept.
 */
template <typename InvariantSystem>
class IKF_belief_transfer_factory : public serializable {
 public:
  using self = IKF_belief_transfer_factory<InvariantSystem>;
  using predictor_type = IKF_belief_transfer<self>;

  using state_space_system = InvariantSystem;
  using state_space_system_ptr = std::shared_ptr<state_space_system>;
  using covariance_type = covariance_matrix<vect_n<double>>;
  using matrix_type = covariance_mat_traits<covariance_type>::matrix_type;

  using time_type = typename discrete_sss_traits<state_space_system>::time_type;
  using time_difference_type =
      typename discrete_sss_traits<state_space_system>::time_difference_type;

  using input_type =
      typename discrete_sss_traits<state_space_system>::input_type;

  template <BeliefSpace BSpace>
  struct predictor {
    using type = predictor_type;
  };

 private:
  /// Holds the reference to the system used for the filter.
  state_space_system_ptr sys;
  /// Holds the system's input noise covariance matrix.
  matrix_type Q;
  /// Holds the system's output measurement's covariance matrix.
  matrix_type R;

 public:
  /**
   * Parametrized constructor.
   * \param aSys The reference to the system used for the filter.
   * \param aQ The system's input noise covariance matrix.
   * \param aR The system's output measurement's covariance matrix.
   */
  explicit IKF_belief_transfer_factory(const state_space_system_ptr& aSys,
                                       matrix_type aQ = matrix_type(),
                                       matrix_type aR = matrix_type())
      : sys(aSys), Q(std::move(aQ)), R(std::move(aR)) {}

  IKF_belief_transfer_factory()
      : IKF_belief_transfer_factory(state_space_system_ptr()) {}

  /**
   * Returns the time-step of the discrete-time system.
   * \return The time-step of the discrete-time system.
   */
  time_difference_type get_time_step() const { return sys->get_time_step(); }

  /**
   * Sets the state-space system used by this kalman filter transfer factory.
   * \param aSys The new state-space system, by shared-pointer.
   */
  void set_state_space_system(const state_space_system_ptr& aSys) {
    sys = aSys;
  }
  /**
   * Gets the state-space system used by this kalman filter transfer factory.
   * \param aSys The new state-space system, by shared-pointer.
   */
  const state_space_system_ptr& get_state_space_system() const { return sys; }

  /**
   * Sets the system input disturbance covariance matrix used by this kalman filter transfer factory.
   * \param aQ The new system input disturbance covariance matrix.
   */
  void set_input_disturbance_cov(const matrix_type& aQ) { Q = aQ; }
  /**
   * Returns the system input disturbance covariance matrix used by this kalman filter transfer factory.
   * \return The system input disturbance covariance matrix.
   */
  const matrix_type& get_input_disturbance_cov() const { return Q; }

  /**
   * Sets the system measurement noise covariance matrix used by this kalman filter transfer factory.
   * \param aR The new system measurement noise covariance matrix.
   */
  void set_measurement_noise_cov(const matrix_type& aR) { R = aR; }
  /**
   * Returns the system measurement noise covariance matrix used by this kalman filter transfer factory.
   * \return The system measurement noise covariance matrix.
   */
  const matrix_type& get_measurement_noise_cov() const { return R; }

  template <BeliefSpace BSpace>
  predictor_type create_predictor(
      const BSpace& /*unused*/,
      const typename pp::topology_traits<BSpace>::point_type* /*unused*/,
      const time_type& /*unused*/, const input_type& /*unused*/) const {
    return predictor_type(this);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(sys) & RK_SERIAL_SAVE_WITH_NAME(Q) &
        RK_SERIAL_SAVE_WITH_NAME(R);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(sys) & RK_SERIAL_LOAD_WITH_NAME(Q) &
        RK_SERIAL_LOAD_WITH_NAME(R);
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2320002, 1,
                              "IKF_belief_transfer_factory", serializable)
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_CONTROLLERS_INVARIANT_KALMAN_FILTER_H_

/**
 * \file invariant_symplectic_kalman_filter.hpp
 * 
 * This library provides a number of functions and classes to do state estimation 
 * using the Invariant Symplectic Kalman Filter. This Kalman filtering technique applies to a 
 * gaussian belief state where the covariance is decomposed into a covarying matrix
 * and an informing matrix (i.e. P = X * invert(Y)), and  an invariant state-space system. 
 * The transition between covariances
 * is achieved using the tranformation matrices (which constitutes a symplectic mapping).
 * The matrices can be multiplied together beyond a single 
 * estimation step and can be aggregated over several steps. If the system is non-linear,
 * then it would need the recomputation of those transformation matrices for the individual
 * steps if the mean-states change too much. The estimation functions will output the 
 * transformation matrices, allowing the caller to aggregate them if needed. If the aggregation
 * is not necessary, then there is no need to use this filter technique as it differs in 
 * no way from the regular Kalman filter (see kalman_filter.hpp).
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

#ifndef REAK_INVARIANT_SYMPLECTIC_KALMAN_FILTER_HPP
#define REAK_INVARIANT_SYMPLECTIC_KALMAN_FILTER_HPP

#include "belief_state_concept.hpp"
#include "discrete_linear_sss_concept.hpp"
#include "invariant_system_concept.hpp"
#include <boost/utility/enable_if.hpp>
#include "lin_alg/vect_concepts.hpp"
#include "lin_alg/mat_alg.hpp"

#include <boost/static_assert.hpp>
#include "covariance_concept.hpp"
#include "lin_alg/mat_cholesky.hpp"
#include "lin_alg/mat_qr_decomp.hpp"
#include "lin_alg/mat_views.hpp"

#include "gaussian_belief_state.hpp"
#include "covariance_matrix.hpp"
#include "decomp_covariance_matrix.hpp"

#include "path_planning/metric_space_concept.hpp"

namespace ReaK {

namespace ctrl {






/**
 * This function template performs one prediction step using the Invariant Symplectic Kalman Filter method.
 * \tparam InvariantSystem An invariant discrete-time state-space system modeling the 
 *         InvariantDiscreteSystemConcept.
 * \tparam StateSpaceType A topology type on which the state-vectors can reside, should model
 *         the pp::TopologyConcept.
 * \tparam BeliefState A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam InputBelief A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam InvCovTransMatrix A matrix type to store the invariant covariance transformation matrix.
 * \tparam InvFrameMatrix A matrix type to store the invariant frame transformation.
 * \param sys The invariant discrete-time state-space system used in the state estimation.
 * \param b As input, it stores the belief-state before the prediction step. As output, it stores
 *        the belief-state after the prediction step.
 * \param u The input vector to apply to the state-space system to make the transition of the 
 *        mean-state, i.e., the current input vector.
 * \param Q The input noise covariance matrix. This is the level of uncertainty on the input 
 *        vector components (not the noise on the state transition). This was chosen as the most
 *        common way application (usually system disturbance comes on the input, not on the state).
 * \param Tc Stores, as output, the transformation matrix (a symplectic mapping) which can transfer the 
 *        covarying and informing components of the covariance matrix of the belief-state. However, note 
 *        that this transformation matrix does not include the invariant frame transformation Wp, which 
 *        means that the resulting covariance in the output belief-state is actually P = diag(Wp,Wp) * Tc * P_prev,
 *        where P is, of course, the two components of the decomposed covariance matrix stacked vertically.
 * \param Wp The invariant frame transformation matrix.
 * \param t The current time (before the prediction).
 * 
 */
template <typename InvariantSystem,  
          typename StateSpaceType,
          typename BeliefState, 
          typename InputBelief,
          typename InvCovTransMatrix,
          typename InvFrameMatrix>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal) &&
                             is_fully_writable_matrix<InvCovTransMatrix>::value &&
                             is_writable_matrix<InvFrameMatrix>::value,
void >::type invariant_symplectic_kf_predict(const InvariantSystem& sys,
                                             const StateSpaceType& state_space,
                                             BeliefState& b_x,
                                             const InputBelief& b_u,
                                             InvCovTransMatrix& Tc,
                                             InvFrameMatrix& Wp,
                                             typename discrete_sss_traits<InvariantSystem>::time_type t = 0) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) update
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) update
  
  typedef typename discrete_sss_traits<InvariantSystem>::point_type StateType;
  typedef typename discrete_sss_traits<InvariantSystem>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  typedef typename invariant_system_traits<InvariantSystem>::invariant_error_type ErrorType;
  
  BOOST_CONCEPT_ASSERT((pp::TopologyConcept< StateSpaceType >));
  BOOST_CONCEPT_ASSERT((InvariantDiscreteSystemConcept<InvariantSystem, StateSpaceType>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<InputBelief>));
  BOOST_CONCEPT_ASSERT((WritableMatrixConcept<InvCovTransMatrix>));
  BOOST_CONCEPT_ASSERT((WritableMatrixConcept<InvFrameMatrix>));
  BOOST_CONCEPT_ASSERT((DecomposedCovarianceConcept<CovType>));
  
  typedef typename decomp_covariance_mat_traits< CovType >::matrix_block_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename mat_traits<MatType>::size_type SizeType;
  
  typename discrete_linear_sss_traits<InvariantSystem>::matrixA_type A;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixB_type B;
  
  StateType x = b_x.get_mean_state();
  const MatType& X = b_x.get_covariance().get_covarying_block();
  const MatType& Y = b_x.get_covariance().get_informing_inv_block(); 
  SizeType N = X.get_col_count();
  
  x = sys.get_next_state(state_space, x, b_u.get_mean_state(), t);
  sys.get_state_transition_blocks(A, B, state_space, t, t + sys.get_time_step(), b_x.get_mean_state(), x, b_u.get_mean_state(), b_u.get_mean_state());
  Wp = sys.get_invariant_prior_frame(state_space, b_x.get_mean_state(), x, b_u.get_mean_state(), t + sys.get_time_step());
  
  Tc.set_row_count(2 * N);
  Tc.set_col_count(2 * N);
  
  set_block(Tc, A, 0, 0);
  typename discrete_linear_sss_traits<InvariantSystem>::matrixA_type A_inv; 
  pseudoinvert_QR(A,A_inv);
  mat_sub_block< InvCovTransMatrix > T_lr(Tc,N,N,N,N);
  T_lr = transpose_view(A_inv);
  mat_sub_block< InvCovTransMatrix > T_ur(Tc,N,N,0,N);
  T_ur = (B * b_u.get_covariance().get_matrix() * transpose_view(B) ) * T_lr;
  set_block(Tc, mat<ValueType,mat_structure::nil>(N), N, 0);
  
  b_x.set_covariance( CovType( MatType( Wp * ( A * X + T_ur * Y ) ), MatType( Wp * T_lr * Y ) ) );
  b_x.set_mean_state(x);
};







/**
 * This function template performs one measurement update step using the Invariant Symplectic Kalman Filter method.
 * \tparam InvariantSystem An invariant discrete-time state-space system modeling the 
 *         InvariantDiscreteSystemConcept.
 * \tparam StateSpaceType A topology type on which the state-vectors can reside, should model
 *         the pp::TopologyConcept.
 * \tparam BeliefState A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam InputBelief A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam MeasurementBelief A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam InvCovTransMatrix A matrix type to store the invariant covariance transformation matrix.
 * \tparam InvFrameMatrix A matrix type to store the invariant frame transformation.
 * \param sys The invariant discrete-time state-space system used in the state estimation.
 * \param state_space The state-space topology on which the state representations lie.
 * \param b_x As input, it stores the belief-state before the update step. As output, it stores
 *        the belief-state after the update step.
 * \param b_u The input vector to apply to the state-space system to make the transition of the 
 *        mean-state, i.e., the current input vector and its covariance.
 * \param b_z The output belief that was measured, i.e. the measurement vector and its covariance.
 * \param Tm Stores, as output, the transformation matrix (a symplectic mapping) which can transfer the 
 *        covarying and informing components of the covariance matrix of the belief-state. However, note 
 *        that this transformation matrix does not include the invariant frame transformation Wu, which 
 *        means that the resulting covariance in the output belief-state is actually P = diag(Wu,Wu) * Tm * P_prev,
 *        where P is, of course, the two components of the decomposed covariance matrix stacked vertically.
 * \param Wu The invariant frame transformation matrix.
 * \param t The current time.
 * 
 */
template <typename InvariantSystem, 
          typename StateSpaceType,
          typename BeliefState, 
          typename InputBelief, 
          typename MeasurementBelief,
          typename InvCovTransMatrix,
          typename InvFrameMatrix>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal) &&
                             is_fully_writable_matrix<InvCovTransMatrix>::value &&
                             is_writable_matrix<InvFrameMatrix>::value,
void >::type invariant_symplectic_kf_update(const InvariantSystem& sys,
                                            const StateSpaceType& state_space,
                                            BeliefState& b_x,
                                            const InputBelief& b_u,
                                            const MeasurementBelief& b_z,
                                            InvCovTransMatrix& Tm,
                                            InvFrameMatrix& Wu,
                                            typename discrete_sss_traits<InvariantSystem>::time_type t = 0) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) update
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) update
  
  typedef typename discrete_sss_traits<InvariantSystem>::point_type StateType;
  typedef typename discrete_sss_traits<InvariantSystem>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  typedef typename invariant_system_traits<InvariantSystem>::invariant_error_type ErrorType;
  typedef typename invariant_system_traits<InvariantSystem>::invariant_correction_type CorrType;
  
  BOOST_CONCEPT_ASSERT((pp::TopologyConcept< StateSpaceType >));
  BOOST_CONCEPT_ASSERT((InvariantDiscreteSystemConcept<InvariantSystem, StateSpaceType>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<InputBelief>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<MeasurementBelief>));
  BOOST_CONCEPT_ASSERT((WritableMatrixConcept<InvCovTransMatrix>));
  BOOST_CONCEPT_ASSERT((WritableMatrixConcept<InvFrameMatrix>));
  BOOST_CONCEPT_ASSERT((DecomposedCovarianceConcept<CovType>));
  
  typedef typename decomp_covariance_mat_traits< CovType >::matrix_block_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename mat_traits<MatType>::size_type SizeType;
  
  typename discrete_linear_sss_traits<InvariantSystem>::matrixC_type C;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixD_type D;
  
  StateType x = b_x.get_mean_state();
  const MatType& X = b_x.get_covariance().get_covarying_block();
  const MatType& Y = b_x.get_covariance().get_informing_inv_block(); 
  SizeType N = X.get_col_count();
  
  Tm.set_row_count(2 * N);
  Tm.set_col_count(2 * N);
  
  sys.get_output_function_blocks(C, D, state_space, t, x, b_u.get_mean_state());
  vect_n<ValueType> e = to_vect<ValueType>(sys.get_output_error(state_space, x, b_u.get_mean_state(), b_z.get_mean_state(), t));
  
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > Ct(transpose_view(C));
  mat< ValueType, mat_structure::symmetric > M = Ct * b_z.get_covariance().get_inverse_matrix() * C;
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > YC;
  linlsq_QR(Y,YC,Ct);
  mat< ValueType, mat_structure::symmetric > S = C * X * YC + b_z.get_covariance().get_matrix();
  YC = transpose_view(X * YC);
  linsolve_Cholesky(S,YC);
  mat< ValueType, mat_structure::rectangular, mat_alignment::row_major > K = transpose_view(YC);
   
  b_x.set_mean_state( sys.apply_correction(state_space, x, from_vect<CorrType>(K * e), b_u.get_mean_state(), t) );
  Wu = sys.get_invariant_posterior_frame(state_space, x, b_x.get_mean_state(), b_u.get_mean_state(), t);
  
  set_block(Tm, M, N, 0);
  set_block(Tm, mat< ValueType, mat_structure::identity>(N), 0, 0);
  set_block(Tm, mat< ValueType, mat_structure::identity>(N), N, N);
    
  b_x.set_covariance( CovType( MatType( Wu * X ), MatType( Wu * ( Y + M * X ) ) ) );
};





/**
 * This function template performs one complete estimation step using the Invariant Symplectic Kalman 
 * Filter method, which includes a prediction and measurement update step. This function is, 
 * in general, more efficient than applying the prediction and update separately.
 * \tparam InvariantSystem An invariant discrete-time state-space system modeling the 
 *         InvariantDiscreteSystemConcept.
 * \tparam StateSpaceType A topology type on which the state-vectors can reside, should model
 *         the pp::TopologyConcept.
 * \tparam BeliefState A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam InputBelief A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam MeasurementBelief A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam InvCovTransMatrix A matrix type to store the invariant covariance transformation matrix.
 * \tparam InvFrameMatrix A matrix type to store the invariant frame transformation.
 * \param sys The invariant discrete-time state-space system used in the state estimation.
 * \param b As input, it stores the belief-state before the estimation step. As output, it stores
 *        the belief-state after the estimation step.
 * \param u The input vector to apply to the state-space system to make the transition of the 
 *        mean-state, i.e., the current input vector.
 * \param z The output vector to that was measured.
 * \param Q The input noise covariance matrix. This is the level of uncertainty on the input 
 *        vector components (not the noise on the state transition). This was chosen as the most
 *        common way application (usually system disturbance comes on the input, not on the state).
 * \param R The output noise covariance matrix. This is the level of uncertainty on the output 
 *        vector components coming from the measurements.
 * \param T Stores, as output, the transformation matrix (a symplectic mapping) which can transfer the 
 *        covarying and informing components of the covariance matrix of the belief-state. However, note 
 *        that this transformation matrix does not include the invariant frame transformation W, which 
 *        means that the resulting covariance in the output belief-state is actually P = diag(W,W) * T * P_prev,
 *        where P is, of course, the two components of the decomposed covariance matrix stacked vertically.
 * \param W The invariant frame transformation matrix.
 * \param t The current time (before the prediction).
 * 
 */
template <typename InvariantSystem, 
          typename StateSpaceType,
          typename BeliefState, 
          typename InputBelief, 
          typename MeasurementBelief,
          typename InvCovTransMatrix,
          typename InvFrameMatrix>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal) &&
                             is_fully_writable_matrix<InvCovTransMatrix>::value &&
                             is_writable_matrix<InvFrameMatrix>::value,
void >::type invariant_symplectic_kf_step(const InvariantSystem& sys,
                                          const StateSpaceType& state_space,
                                          BeliefState& b_x,
                                          const InputBelief& b_u,
                                          const MeasurementBelief& b_z,
                                          InvCovTransMatrix& T,
                                          InvFrameMatrix& W,
                                          typename discrete_sss_traits<InvariantSystem>::time_type t = 0) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) update
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) update
  
  typedef typename discrete_sss_traits<InvariantSystem>::point_type StateType;
  typedef typename discrete_sss_traits<InvariantSystem>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  typedef typename invariant_system_traits<InvariantSystem>::invariant_error_type ErrorType;
  typedef typename invariant_system_traits<InvariantSystem>::invariant_correction_type CorrType;
  
  BOOST_CONCEPT_ASSERT((pp::TopologyConcept< StateSpaceType >));
  BOOST_CONCEPT_ASSERT((InvariantDiscreteSystemConcept<InvariantSystem, StateSpaceType>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<InputBelief>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<MeasurementBelief>));
  BOOST_CONCEPT_ASSERT((WritableMatrixConcept<InvCovTransMatrix>));
  BOOST_CONCEPT_ASSERT((WritableMatrixConcept<InvFrameMatrix>));
  BOOST_CONCEPT_ASSERT((DecomposedCovarianceConcept<CovType>));
  
  typedef typename decomp_covariance_mat_traits< CovType >::matrix_block_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename mat_traits<MatType>::size_type SizeType;
  
  typename discrete_linear_sss_traits<InvariantSystem>::matrixA_type A;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixB_type B;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixC_type C;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixD_type D;
  
  StateType x = b_x.get_mean_state();
  MatType X = b_x.get_covariance().get_covarying_block();
  MatType Y = b_x.get_covariance().get_informing_inv_block(); 
  SizeType N = X.get_col_count();
  
  x = sys.get_next_state(state_space, x, b_u.get_mean_state(), t);
  sys.get_state_transition_blocks(A, B, state_space, t, t + sys.get_time_step(), b_x.get_mean_state(), x, b_u.get_mean_state(), b_u.get_mean_state());
  W = sys.get_invariant_prior_frame(state_space, b_x.get_mean_state(), x, b_u.get_mean_state(), t + sys.get_time_step());
  
  T.set_row_count(2 * N);
  T.set_col_count(2 * N);
  
  set_block(T, A, 0, 0);
  typename discrete_linear_sss_traits<InvariantSystem>::matrixA_type A_inv; 
  pseudoinvert_QR(A,A_inv);
  mat_sub_block< InvCovTransMatrix > T_lr(T,N,N,N,N);
  T_lr = transpose_view(A_inv);
  mat_sub_block< InvCovTransMatrix > T_ur(T,N,N,0,N);
  T_ur = (B * b_u.get_covariance().get_matrix() * transpose_view(B) ) * T_lr;
  set_block(T, mat<ValueType,mat_structure::nil>(N), N, 0);
  
  X = A * X + T_ur * Y;
  Y = T_lr * Y;
  
  sys.get_output_function_blocks(C, D, state_space, t + sys.get_time_step(), x, b_u.get_mean_state());
  vect_n<ValueType> e = to_vect<ValueType>(sys.get_output_error(state_space, x, b_u.get_mean_state(), b_z.get_mean_state(), t + sys.get_time_step()));
  
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > Ct(transpose_view(C));
  mat< ValueType, mat_structure::symmetric > M = Ct * b_z.get_covariance().get_inverse_matrix() * C;
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > YC;
  linlsq_QR(Y,YC,Ct);
  mat< ValueType, mat_structure::symmetric > S = C * X * YC + b_z.get_covariance().get_matrix();
  YC = transpose_view(X * YC);
  linsolve_Cholesky(S,YC);
  mat< ValueType, mat_structure::rectangular, mat_alignment::row_major > K(transpose_view(YC));
   
  b_x.set_mean_state( sys.apply_correction(state_space, x, from_vect<CorrType>(W * K * e), b_u.get_mean_state(), t + sys.get_time_step()) );
  W = sys.get_invariant_posterior_frame(state_space, x, b_x.get_mean_state(), b_u.get_mean_state(), t + sys.get_time_step()) * W;
  
  set_block(T, M * A, N, 0);
  T_lr += M * T_ur;
    
  b_x.set_covariance( CovType( MatType( W * X ), MatType( W * ( Y + M * X ) ) ) );
};








/**
 * This class template can be used as a belief-state predictor (and transfer) that uses the 
 * Invariant Symplectic Kalman Filter method. This class template models the BeliefTransferConcept and 
 * the BeliefPredictorConcept.
 * \tparam InvariantSystem An invariant discrete-time state-space system modeling the 
 *         InvariantDiscreteSystemConcept.
 * \tparam BeliefState A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam SystemNoiseCovar A covariance matrix type modeling the CovarianceMatrixConcept.
 * \tparam MeasurementCovar A covariance matrix type modeling the CovarianceMatrixConcept.
 */
template <typename InvariantSystem,
          typename BeliefState = gaussian_belief_state< typename discrete_sss_traits<InvariantSystem>::point_type, decomp_covariance_matrix< typename discrete_sss_traits<InvariantSystem>::point_type > >,
          typename SystemNoiseCovar = covariance_matrix< typename discrete_sss_traits< InvariantSystem >::input_type >,
          typename MeasurementCovar = covariance_matrix< typename discrete_sss_traits< InvariantSystem >::output_type > >
struct ISKF_belief_transfer {
  typedef ISKF_belief_transfer<InvariantSystem, BeliefState> self;
  typedef BeliefState belief_state;
  typedef InvariantSystem state_space_system;
  typedef shared_ptr< state_space_system > state_space_system_ptr;
  typedef typename discrete_sss_traits< state_space_system >::time_type time_type;
  typedef typename discrete_sss_traits< state_space_system >::time_difference_type time_difference_type;

  typedef typename belief_state_traits< belief_state >::state_type state_type;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type covariance_type;
  typedef typename decomp_covariance_mat_traits< covariance_type >::matrix_block_type matrix_type;
  typedef typename mat_traits< matrix_type >::value_type mat_value_type;
  typedef mat_value_type value_type;
  typedef typename mat_traits< matrix_type >::size_type mat_size_type;

  typedef typename discrete_sss_traits< state_space_system >::input_type input_type;
  typedef typename discrete_sss_traits< state_space_system >::output_type output_type;
  
  typedef gaussian_belief_state<input_type, SystemNoiseCovar> input_belief_type;
  typedef gaussian_belief_state<output_type, MeasurementCovar> output_belief_type;
  
  BOOST_CONCEPT_ASSERT((InvariantDiscreteSystemConcept<InvariantSystem>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  BOOST_CONCEPT_ASSERT((CovarianceMatrixConcept<SystemNoiseCovar,input_type>));
  BOOST_CONCEPT_ASSERT((CovarianceMatrixConcept<MeasurementCovar,output_type>));
  BOOST_CONCEPT_ASSERT((DecomposedCovarianceConcept<covariance_type>));
  
  state_space_system_ptr sys; ///< Holds the reference to the system used for the filter.
  SystemNoiseCovar Q; ///< Holds the system's input noise covariance matrix.
  MeasurementCovar R; ///< Holds the system's output measurement's covariance matrix.
  value_type reupdate_threshold; ///< The threshold at which the state change is considered too high and the state transition matrices are recomputed.
  state_type initial_state; ///< The initial mean-state at which the predictor is linearized (if non-linear).
  state_type predicted_state; ///< The predicted mean-state at which the updator is linearized (if non-linear).
  mat< mat_value_type, mat_structure::square > Tc; ///< Holds the prediction invariant covariance transformation matrix.
  mat< mat_value_type, mat_structure::square > Wp; ///< Holds the prediction invariant frame transformation.
  mat< mat_value_type, mat_structure::square > Tm; ///< Holds the updating invariant covariance transformation matrix.
  mat< mat_value_type, mat_structure::square > Wu; ///< Holds the updating invariant frame transformation.

  /**
   * Parametrized constructor.
   * \param aSys The reference to the system used for the filter.
   * \param state_space The state-space in which to operate.
   * \param aQ The system's input noise covariance matrix.
   * \param aR The system's output measurement's covariance matrix.
   * \param aReupdateThreshold The threshold at which the state change is considered too high and the state transition matrices are recomputed.
   * \param aInitialState The initial mean-state at which the predictor is linearized (if non-linear).
   * \param aInitialInput The input at the time of the initial state.
   * \param aInitialTime The time of the initial state.
   */
  template <typename StateSpaceType>
  ISKF_belief_transfer(const state_space_system_ptr& aSys, 
                       const StateSpaceType& state_space,
                       const SystemNoiseCovar& aQ,
                       const MeasurementCovar& aR,
                       const value_type& aReupdateThreshold,
                       const state_type& aInitialState,
                       const input_type& aInitialInput,
                       const time_type& aInitialTime) : 
                       sys(aSys), Q(aQ), R(aR),
                       reupdate_threshold(aReupdateThreshold),
                       initial_state(aInitialState) {
    belief_state b = belief_state(initial_state, 
                                  covariance_type(state_space.difference(aInitialState, aInitialState).size(), 
                                                  covariance_initial_level::no_info));
    input_belief_type b_u(aInitialInput, Q);
    invariant_symplectic_kf_predict(*sys,
                                    b,
                                    b_u,
                                    Tc,
                                    Wp,
                                    aInitialTime);
    predicted_state = b.get_mean_state();
    invariant_symplectic_kf_update(*sys,
                                   b,
                                   b_u,
                                   output_belief_type( 
                                     sys->get_output(predicted_state,
                                                     aInitialInput,
                                                     aInitialTime),
                                     R),
                                   Tm,
                                   Wu,
                                   aInitialTime);
  };
  
  /**
   * Returns the time-step of the predictor.
   * \return The time-step of the predictor.
   */
  time_difference_type get_time_step() const { return sys->get_time_step(); };

  /**
   * Returns a reference to the underlying state-space system.
   * \return A reference to the underlying state-space system.
   */
  const state_space_system& get_ss_system() const { return *sys; };

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
  belief_state get_next_belief(const BeliefSpace& b_space, belief_state b, const time_type& t, const input_type& u, const input_type& y) {
    initial_state = b.get_mean_state();
    input_belief_type b_u(u, Q);
    invariant_symplectic_kf_predict(*sys,b_space.get_state_topology(),b,b_u,Tc,Wp,t);
    predicted_state = b.get_mean_state();
    invariant_symplectic_kf_update(*sys,b_space.get_state_topology(),b,b_u,output_belief_type(y,R),Tm,Wu,t);
    return b;
  };
  
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
  belief_state predict_belief(const BeliefSpace& b_space, belief_state b, const time_type& t, const input_type& u) {
    if( norm( b_space.get_state_topology().difference(b.get_mean_state(), initial_state) ) > reupdate_threshold ) {
      initial_state = b.get_mean_state();
      input_belief_type b_u(u, Q);
      invariant_symplectic_kf_predict(*sys,b_space.get_state_topology(),b,b_u,Tc,Wp,t);
    } else {
      state_type x = sys->get_next_state(b_space.get_state_topology(), b.get_mean_state(), u, t);
      Wp = sys->get_invariant_prior_frame(b_space.get_state_topology(), b.get_mean_state(), x, u, t);
      b.set_mean_state( x );
      mat_size_type N = Tc.get_row_count() / 2;
      mat< mat_value_type, mat_structure::rectangular > P_tmp = 
        Tc * ( b.get_covariance().get_covarying_block() |
               b.get_covariance().get_informing_inv_block() );
      b.set_covariance( covariance_type( matrix_type( Wp * sub(P_tmp)(range(0,N-1),range(0,N-1)) ), 
                                         matrix_type( Wp * sub(P_tmp)(range(N,2*N-1),range(0,N-1)) ) ) ); 
    };
    return b;
  };
  
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
  belief_state prediction_to_ML_belief(const BeliefSpace& b_space, belief_state b, const time_type& t, const input_type& u) {
    if( norm( b_space.get_state_topology().difference(b.get_mean_state(), predicted_state) ) > reupdate_threshold ) {
      predicted_state = b.get_mean_state();
      input_belief_type b_u(u, Q);
      invariant_symplectic_kf_update(*sys,b_space.get_state_topology(),b,b_u,output_belief_type(sys->get_output(b_space.get_state_topology(), predicted_state, b_u, t + sys->get_time_step()),R),Tm,Wu,t);
    } else {
      mat_size_type N = Tm.get_row_count() / 2;
      Wu = mat< mat_value_type, mat_structure::identity >(N);
      mat< mat_value_type, mat_structure::rectangular > P_tmp = 
        Tm * ( b.get_covariance().get_covarying_block() |
               b.get_covariance().get_informing_inv_block() );
      b.set_covariance( covariance_type( matrix_type( sub(P_tmp)(range(0,N-1),range(0,N-1)) ), 
                                         matrix_type( sub(P_tmp)(range(N,2*N-1),range(0,N-1)) ) ) ); 
    };
    return b;
  };
  
  /**
   * Returns the prediction belief-state at the next time instant, assuming the upcoming measurement to be the most likely one.
   * \tparam BeliefSpace The belief-space type on which to operate.
   * \param b_space The belief-space on which the belief-states lie.
   * \param b The current belief-state.
   * \param t The current time.
   * \param u The current input given to the system.
   * \return the belief-state at the next time instant, predicted by the filter.
   */
  template <typename BeliefSpace>
  belief_state predict_ML_belief(const BeliefSpace& b_space, belief_state b, const time_type& t, const input_type& u) {
    return prediction_to_ML_belief(b_space, predict_belief(b_space, b, t, u),t,u);
  };
  
};










};

};


#endif




















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

#include <ReaK/math/lin_alg/vect_concepts.hpp>
#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/math/lin_alg/mat_cholesky.hpp>
#include <ReaK/math/lin_alg/mat_qr_decomp.hpp>
#include <ReaK/math/lin_alg/mat_views.hpp>

#include <ReaK/topologies/spaces/metric_space_concept.hpp>

#include "belief_state_concept.hpp"
#include <ReaK/control/systems/discrete_linear_sss_concept.hpp>
#include <ReaK/control/systems/invariant_system_concept.hpp>
#include "covariance_concept.hpp"
#include "gaussian_belief_state.hpp"
#include "covariance_matrix.hpp"
#include "decomp_covariance_matrix.hpp"

#include <boost/static_assert.hpp>
#include <boost/utility/enable_if.hpp>

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
 * \tparam ISKFTransferFactory The factory type which can create this kalman predictor.
 * \tparam BeliefSpace The belief-space type on which to operate.
 */
template <typename ISKFTransferFactory, typename BeliefSpace>
struct ISKF_belief_transfer {
  typedef ISKF_belief_transfer<ISKFTransferFactory, BeliefSpace> self;
  typedef BeliefSpace belief_space;
  typedef typename pp::topology_traits<BeliefSpace>::point_type belief_state;
  
  typedef typename ISKFTransferFactory::state_space_system state_space_system;
  typedef shared_ptr< state_space_system > state_space_system_ptr;
  typedef typename discrete_sss_traits< state_space_system >::time_type time_type;
  typedef typename discrete_sss_traits< state_space_system >::time_difference_type time_difference_type;

  typedef typename belief_state_traits< belief_state >::state_type state_type;
  
  typedef double value_type;
  typedef typename continuous_belief_state_traits< belief_state >::covariance_type covariance_type;
  typedef typename decomp_covariance_mat_traits< covariance_type >::matrix_block_type matrix_type;
  
  typedef typename discrete_sss_traits< state_space_system >::input_type input_type;
  typedef typename discrete_sss_traits< state_space_system >::output_type output_type;
  
  typedef covariance_matrix< vect_n<double> > io_covariance_type;
  typedef gaussian_belief_state<input_type,  io_covariance_type> input_belief_type;
  typedef gaussian_belief_state<output_type, io_covariance_type> output_belief_type;
  
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<belief_state>));
  BOOST_CONCEPT_ASSERT((DecomposedCovarianceConcept<covariance_type, vect_n<double> >));
  
  const ISKFTransferFactory* factory;
  
  state_type initial_state; ///< The initial mean-state at which the predictor is linearized (if non-linear).
  state_type predicted_state; ///< The predicted mean-state at which the updator is linearized (if non-linear). 
  mat< double, mat_structure::square > Tc; ///< Holds the prediction invariant covariance transformation matrix.
  mat< double, mat_structure::square > Wp; ///< Holds the prediction invariant frame transformation.
  mat< double, mat_structure::square > Tm; ///< Holds the updating invariant covariance transformation matrix.
  mat< double, mat_structure::square > Wu; ///< Holds the updating invariant frame transformation.

  ISKF_belief_transfer() : factory(NULL) { };
  
  /**
   * Parametrized constructor.
   * \param aFactory A pointer to the factory object that is creating this object.
   * \param b_space The belief-space in which to operate.
   * \param b_x The initial belief-state at which the predictor is linearized (if non-linear).
   * \param t The time of the initial state.
   * \param u The input at the time of the initial state.
   */
  ISKF_belief_transfer(const ISKFTransferFactory* aFactory,
                       const belief_space& b_space,
                       const belief_state* b_x,
                       const time_type& t,
                       const input_type& u) : 
                       factory(aFactory),
                       initial_state(b_x->get_mean_state()) {
    belief_state b = *b_x;
    input_belief_type b_u(u, io_covariance_type(factory->get_input_disturbance_cov()));
    invariant_symplectic_kf_predict(*(factory->get_state_space_system()), 
                                    b_space.get_state_topology(),
                                    b, b_u, Tc, Wp, t);
    predicted_state = b.get_mean_state();
    output_belief_type b_y(factory->get_state_space_system()->get_output(b_space.get_state_topology(), predicted_state, u, t + factory->get_state_space_system()->get_time_step()), io_covariance_type(factory->get_measurement_noise_cov()));
    invariant_symplectic_kf_update(*(factory->get_state_space_system()), 
                                   b_space.get_state_topology(),
                                   b, b_u, b_y, Tm, Wu, t + factory->get_state_space_system()->get_time_step());
  };
  
  /**
   * Returns the time-step of the predictor.
   * \return The time-step of the predictor.
   */
  time_difference_type get_time_step() const { return factory->get_time_step(); };
  
  /**
   * Returns a reference to the underlying state-space system.
   * \return A reference to the underlying state-space system.
   */
  const state_space_system_ptr& get_ss_system() const { return factory->get_state_space_system(); };
  
  /**
   * Returns the belief-state at the next time instant.
   * \param b_space The belief-space on which the belief-states lie.
   * \param b The current belief-state.
   * \param t The current time.
   * \param u The current input given to the system.
   * \param y The output that was measured at the next time instant.
   * \return the belief-state at the next time instant.
   */
  belief_state get_next_belief(const belief_space& b_space, belief_state b, const time_type& t, const input_type& u, const output_type& y) {
    initial_state = b.get_mean_state();
    input_belief_type b_u(u, io_covariance_type(factory->get_input_disturbance_cov()));
    invariant_symplectic_kf_predict(*(factory->get_state_space_system()), b_space.get_state_topology(), b, b_u, Tc, Wp, t);
    predicted_state = b.get_mean_state();
    output_belief_type b_y(y, io_covariance_type(factory->get_measurement_noise_cov()))
    invariant_symplectic_kf_update(*(factory->get_state_space_system()), b_space.get_state_topology(), b, b_u, b_y, Tm, Wu, t + factory->get_state_space_system()->get_time_step());
    return b;
  };
  
  /**
   * Returns the prediction belief-state at the next time instant.
   * \param b_space The belief-space on which the belief-states lie.
   * \param b The current belief-state.
   * \param t The current time.
   * \param u The current input given to the system.
   * \return the belief-state at the next time instant, predicted by the filter.
   */
  belief_state predict_belief(const belief_space& b_space, belief_state b, const time_type& t, const input_type& u) {
    double dist = get(pp::distance_metric, b_space.get_state_topology())(b.get_mean_state(), initial_state, b_space.get_state_topology());
    if( dist > factory->get_reupdate_threshold() ) {
      initial_state = b.get_mean_state();
      input_belief_type b_u(u, io_covariance_type(factory->get_input_disturbance_cov()));
      invariant_symplectic_kf_predict(*(factory->get_state_space_system()), b_space.get_state_topology(), b, b_u, Tc, Wp, t);
    } else {
      state_type x = factory->get_state_space_system()->get_next_state(b_space.get_state_topology(), b.get_mean_state(), u, t);
      Wp = factory->get_state_space_system()->get_invariant_prior_frame(b_space.get_state_topology(), b.get_mean_state(), x, u, t);
      b.set_mean_state( x );
      std::size_t N = Tc.get_row_count() / 2;
      mat< double, mat_structure::rectangular > P_tmp = 
        Tc * ( b.get_covariance().get_covarying_block() |
               b.get_covariance().get_informing_inv_block() );
      b.set_covariance( covariance_type( matrix_type( Wp * sub(P_tmp)(range(0,N),range(0,N)) ), 
                                         matrix_type( Wp * sub(P_tmp)(range(N,2*N),range(0,N)) ) ) ); 
    };
    return b;
  };
  
  /**
   * Converts a prediction belief-state into an updated belief-state which assumes the most likely measurement.
   * \param b_space The belief-space on which the belief-states lie.
   * \param b The current prediction's belief-state.
   * \param t The current time.
   * \param u The current input given to the system.
   * \return the updated belief-state when assuming the most likely measurement.
   */
  belief_state prediction_to_ML_belief(const belief_space& b_space, belief_state b, const time_type& t, const input_type& u) {
    double dist = get(pp::distance_metric, b_space.get_state_topology())(b.get_mean_state(), predicted_state, b_space.get_state_topology());
    if( dist > factory->get_reupdate_threshold() ) {
      predicted_state = b.get_mean_state();
      input_belief_type b_u(u, io_covariance_type(factory->get_input_disturbance_cov()));
      output_belief_type b_y(factory->get_state_space_system()->get_output(b_space.get_state_topology(), predicted_state, u, t + factory->get_state_space_system()->get_time_step()), io_covariance_type(factory->get_measurement_noise_cov()));
      invariant_symplectic_kf_update(*(factory->get_state_space_system()), b_space.get_state_topology(), b, b_u, b_y, Tm, Wu, t + factory->get_state_space_system()->get_time_step());
    } else {
      std::size_t N = Tm.get_row_count() / 2;
      Wu = mat_ident<double>(N);
      mat< double, mat_structure::rectangular > P_tmp = 
        Tm * ( b.get_covariance().get_covarying_block() |
               b.get_covariance().get_informing_inv_block() );
      b.set_covariance( covariance_type( matrix_type( sub(P_tmp)(range(0,N),range(0,N)) ), 
                                         matrix_type( sub(P_tmp)(range(N,2*N),range(0,N)) ) ) ); 
    };
    return b;
  };
  
  /**
   * Returns the prediction belief-state at the next time instant, assuming the upcoming measurement to be the most likely one.
   * \param b_space The belief-space on which the belief-states lie.
   * \param b The current belief-state.
   * \param t The current time.
   * \param u The current input given to the system.
   * \return the belief-state at the next time instant, predicted by the filter.
   */
  belief_state predict_ML_belief(const belief_space& b_space, belief_state b, const time_type& t, const input_type& u) {
    return prediction_to_ML_belief(b_space, predict_belief(b_space, b, t, u), t, u);
  };
  
};





/**
 * This class is a factory class for invariant symplectic Kalman filtering predictors on a belief-space.
 * \tparam LinearSystem A discrete state-space system modeling the DiscreteLinearSSSConcept 
 *         at least as a DiscreteLinearizedSystemType.
 */
template <typename InvariantSystem>
class ISKF_belief_transfer_factory : public serializable {
  public:
    typedef ISKF_belief_transfer_factory<InvariantSystem> self;
    
    typedef InvariantSystem state_space_system;
    typedef shared_ptr< state_space_system > state_space_system_ptr;
    typedef covariance_matrix< vect_n<double> > covariance_type;
    typedef covariance_mat_traits< covariance_type >::matrix_type matrix_type;
    
    typedef typename discrete_sss_traits< state_space_system >::time_type time_type;
    typedef typename discrete_sss_traits< state_space_system >::time_difference_type time_difference_type;
    
    typedef typename discrete_sss_traits< state_space_system >::input_type input_type;
    
    template <typename BeliefSpace>
    struct predictor {
      typedef ISKF_belief_transfer<self, BeliefSpace> type;
    };
    
  private:
    state_space_system_ptr sys; ///< Holds the reference to the system used for the filter.
    matrix_type Q; ///< Holds the system's input noise covariance matrix.
    matrix_type R; ///< Holds the system's output measurement's covariance matrix.
    double reupdate_threshold; ///< The threshold at which the state change is considered too high and the state transition matrices are recomputed.
    
  public:
    
    
    /**
     * Parametrized constructor.
     * \param aSys The reference to the system used for the filter.
     * \param aQ The system's input noise covariance matrix.
     * \param aR The system's output measurement's covariance matrix.
     * \param aReupdateThreshold The threshold at which the state change is considered too high and the state transition matrices are recomputed.
     */
    ISKF_belief_transfer_factory(const state_space_system_ptr& aSys = state_space_system_ptr(), 
                                 const matrix_type& aQ = matrix_type(),
                                 const matrix_type& aR = matrix_type(), 
                                 double aReupdateThreshold = 1e-6) : 
                                 sys(aSys), Q(aQ), R(aR), reupdate_threshold(aReupdateThreshold) { };
    
    /**
     * Returns the time-step of the discrete-time system.
     * \return The time-step of the discrete-time system.
     */
    time_difference_type get_time_step() const { return sys->get_time_step(); };
    
    /**
     * Sets the state-space system used by this kalman filter transfer factory.
     * \param aSys The new state-space system, by shared-pointer.
     */
    void set_state_space_system(const state_space_system_ptr& aSys) { sys = aSys; };
    /**
     * Gets the state-space system used by this kalman filter transfer factory.
     * \param aSys The new state-space system, by shared-pointer.
     */
    const state_space_system_ptr& get_state_space_system() const { return sys; };
    
    /**
     * Sets the system input disturbance covariance matrix used by this kalman filter transfer factory.
     * \param aQ The new system input disturbance covariance matrix.
     */
    void set_input_disturbance_cov(const matrix_type& aQ) { Q = aQ; };
    /**
     * Returns the system input disturbance covariance matrix used by this kalman filter transfer factory.
     * \return The system input disturbance covariance matrix.
     */
    const matrix_type& get_input_disturbance_cov() const { return Q; };
    
    /**
     * Sets the system measurement noise covariance matrix used by this kalman filter transfer factory.
     * \param aR The new system measurement noise covariance matrix.
     */
    void set_measurement_noise_cov(const matrix_type& aR) { R = aR; };
    /**
     * Returns the system measurement noise covariance matrix used by this kalman filter transfer factory.
     * \return The system measurement noise covariance matrix.
     */
    const matrix_type& get_measurement_noise_cov() const { return R; };
    
    /**
     * Sets the re-updating threshold (w.r.t. state-space metric) used by kalman filter predictors.
     * \param aReupdateThreshold The new re-updating threshold (w.r.t. state-space metric).
     */
    void set_reupdate_threshold(double aReupdateThreshold) const { reupdate_threshold = aReupdateThreshold; };
    /**
     * Returns the re-updating threshold (w.r.t. state-space metric) used by kalman filter predictors.
     * \return The re-updating threshold (w.r.t. state-space metric).
     */
    double get_reupdate_threshold() const { return reupdate_threshold; };
    
    template <typename BeliefSpace>
    ISKF_belief_transfer<self, BeliefSpace> 
        create_predictor(const BeliefSpace& b_space, 
                         const typename pp::topology_traits<BeliefSpace>::point_type* pb, 
                         const time_type& t, const input_type& u) const {
      return ISKF_belief_transfer<self, BeliefSpace>(this, b_space, pb, t, u);
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_NAME(sys)
        & RK_SERIAL_SAVE_WITH_NAME(Q)
        & RK_SERIAL_SAVE_WITH_NAME(R)
        & RK_SERIAL_SAVE_WITH_NAME(reupdate_threshold);
    };
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_NAME(sys)
        & RK_SERIAL_LOAD_WITH_NAME(Q)
        & RK_SERIAL_LOAD_WITH_NAME(R)
        & RK_SERIAL_LOAD_WITH_NAME(reupdate_threshold);
    };
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2320004,1,"ISKF_belief_transfer_factory",serializable)
};




};

};


#endif




















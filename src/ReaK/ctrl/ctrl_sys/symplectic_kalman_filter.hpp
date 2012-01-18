/**
 * \file symplectic_kalman_filter.hpp
 * 
 * This library provides a number of functions and classes to do state estimation 
 * using the Symplectic Kalman Filter. This Kalman filtering technique applies to a 
 * gaussian belief state where the covariance is decomposed into a covarying matrix
 * and an informing matrix (i.e. P = X * invert(Y)). The transition between covariances
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

#ifndef REAK_SYMPLECTIC_KALMAN_FILTER_HPP
#define REAK_SYMPLECTIC_KALMAN_FILTER_HPP

#include "belief_state_concept.hpp"
#include "discrete_linear_sss_concept.hpp"
#include <boost/utility/enable_if.hpp>
#include "lin_alg/vect_concepts.hpp"
#include "lin_alg/mat_alg.hpp"

#include <boost/static_assert.hpp>
#include "covariance_concept.hpp"
#include "lin_alg/mat_cholesky.hpp"
#include "lin_alg/mat_qr_decomp.hpp"
#include "lin_alg/mat_views.hpp"

namespace ReaK {

namespace ctrl {




/**
 * This function template performs one prediction step using the Symplectic Kalman Filter method.
 * \tparam LinearSystem A discrete state-space system modeling the DiscreteLinearSSSConcept 
 *         at least as a DiscreteLinearizedSystemType.
 * \tparam BeliefState A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation and a covariance matrix modeling the DecomposedCovarianceConcept.
 * \tparam SystemNoiseCovariance A covariance matrix type modeling the CovarianceMatrixConcept.
 * \tparam PredictionCovTransMatrix A matrix type to store the covariance transformation matrix.
 * \param sys The discrete state-space system used in the state estimation.
 * \param b As input, it stores the belief-state before the prediction step. As output, it stores
 *        the belief-state after the prediction step.
 * \param u The input vector to apply to the state-space system to make the transition of the 
 *        mean-state, i.e., the current input vector.
 * \param Q The input noise covariance matrix. This is the level of uncertainty on the input 
 *        vector components (not the noise on the state transition). This was chosen as the most
 *        common way application (usually system disturbance comes on the input, not on the state).
 * \param Tc Stores, as output, the transformation matrix (a symplectic mapping) which can transfer the 
 *        covarying and informing components of the covariance matrix of the belief-state.
 * \param t The current time (before the prediction).
 * 
 */
template <typename LinearSystem, 
          typename BeliefState, 
	  typename SystemNoiseCovariance,
	  typename PredictionCovTransMatrix>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal) &&
                             is_fully_writable_matrix<PredictionCovTransMatrix>::value,
void >::type symplectic_kalman_predict(const LinearSystem& sys,
				       BeliefState& b,
				       const typename discrete_sss_traits<LinearSystem>::input_type& u,
				       const SystemNoiseCovariance& Q,
				       PredictionCovTransMatrix& Tc,
				       typename discrete_sss_traits<LinearSystem>::time_type t = 0) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) prediction
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) prediction
  
  typedef typename discrete_sss_traits<LinearSystem>::point_type StateType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  
  BOOST_CONCEPT_ASSERT((DiscreteLinearSSSConcept< LinearSystem, DiscreteLinearizedSystemType >));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  BOOST_CONCEPT_ASSERT((CovarianceMatrixConcept<SystemNoiseCovariance>));
  BOOST_CONCEPT_ASSERT((WritableMatrixConcept<PredictionCovTransMatrix>));
  BOOST_CONCEPT_ASSERT((DecomposedCovarianceConcept<CovType>));
  
  typedef typename decomp_covariance_mat_traits< CovType >::matrix_block_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename mat_traits<MatType>::size_type SizeType;
  
  typename discrete_linear_sss_traits<LinearSystem>::matrixA_type A;
  typename discrete_linear_sss_traits<LinearSystem>::matrixB_type B;
  StateType x = b.get_mean_state();
  
  const MatType& X = b.get_covariance().get_covarying_block();
  const MatType& Y = b.get_covariance().get_informing_inv_block(); 
  SizeType N = X.get_row_count();
  
  Tc.set_row_count(2 * N);
  Tc.set_col_count(2 * N);

  b.set_mean_state( sys.get_next_state(x,u,t) );
  sys.get_state_transition_blocks(A, B, t, t + sys.get_time_step(), x, b.get_mean_state(), u, u);
  
  set_block(Tc, mat<ValueType, mat_structure::nil>(N), N, 0);
  typename discrete_linear_sss_traits<LinearSystem>::matrixA_type A_inv; 
  pseudoinvert_QR(A,A_inv);
  set_block(Tc, transpose_view(A_inv), N, N);
  set_block(Tc, (B * Q.get_matrix() * transpose_view(B)) * mat_const_sub_block< PredictionCovTransMatrix >(Tc, N, N, N, N), 0, N);
  set_block(Tc, A, 0, 0);
  
  b.set_covariance( CovType( MatType( mat_const_sub_block< PredictionCovTransMatrix >(Tc, N, 2*N, 0, 0) * mat_const_ref_vert_cat< MatType, MatType >(X,Y) ),
                             MatType( mat_const_sub_block< PredictionCovTransMatrix >(Tc, N, N, N, N) * Y ) ) );
};


/**
 * This function template performs one measurement update step using the Symplectic Kalman Filter method.
 * \tparam LinearSystem A discrete state-space system modeling the DiscreteLinearSSSConcept 
 *         at least as a DiscreteLinearizedSystemType.
 * \tparam BeliefState A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation and a covariance matrix modeling the DecomposedCovarianceConcept.
 * \tparam MeasurementNoiseCovariance A covariance matrix type modeling the CovarianceMatrixConcept.
 * \tparam UpdateCovTransMatrix A matrix type to store the covariance transformation matrix.
 * \param sys The discrete state-space system used in the state estimation.
 * \param b As input, it stores the belief-state before the update step. As output, it stores
 *        the belief-state after the update step.
 * \param u The input vector to apply to the state-space system to make the transition of the 
 *        mean-state, i.e., the current input vector.
 * \param z The output vector to that was measured.
 * \param R The output noise covariance matrix. This is the level of uncertainty on the output 
 *        vector components coming from the measurements.
 * \param Tm Stores, as output, the transformation matrix (a symplectic mapping) which can transfer the 
 *        covarying and informing components of the covariance matrix of the belief-state.
 * \param t The current time.
 * 
 */
template <typename LinearSystem, 
          typename BeliefState, 
	  typename MeasurementNoiseCovariance,
	  typename UpdateCovTransMatrix>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal) &&
                             is_fully_writable_matrix<UpdateCovTransMatrix>::value,
void >::type symplectic_kalman_update(const LinearSystem& sys,
				      BeliefState& b,
				      const typename discrete_sss_traits<LinearSystem>::input_type& u,
				      const typename discrete_sss_traits<LinearSystem>::output_type& z,
				      const MeasurementNoiseCovariance& R,
				      UpdateCovTransMatrix& Tm,
				      typename discrete_sss_traits<LinearSystem>::time_type t = 0) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) update
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) update
  
  typedef typename discrete_sss_traits<LinearSystem>::point_type StateType;
  typedef typename discrete_sss_traits<LinearSystem>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  
  BOOST_CONCEPT_ASSERT((DiscreteLinearSSSConcept< LinearSystem, DiscreteLinearizedSystemType >));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  BOOST_CONCEPT_ASSERT((CovarianceMatrixConcept<MeasurementNoiseCovariance>));
  BOOST_CONCEPT_ASSERT((WritableMatrixConcept<UpdateCovTransMatrix>));
  BOOST_CONCEPT_ASSERT((DecomposedCovarianceConcept<CovType>));
  
  typedef typename decomp_covariance_mat_traits< CovType >::matrix_block_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename mat_traits<MatType>::size_type SizeType;
  
  typename discrete_linear_sss_traits<LinearSystem>::matrixC_type C;
  typename discrete_linear_sss_traits<LinearSystem>::matrixD_type D;
  
  StateType x = b.get_mean_state();
  const MatType& X = b.get_covariance().get_covarying_block();
  const MatType& Y = b.get_covariance().get_informing_inv_block(); 
  sys.get_output_function_blocks(C, D, t, x, u);
  
  OutputType y = z - sys.get_output(x,u,t);
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > YC(transpose_view(C));
  mat< ValueType, mat_structure::symmetric > M = YC * R.get_inverse_matrix() * C;
  linsolve_QR(Y,YC);
  mat< ValueType, mat_structure::symmetric > S = C * X * YC + R.get_matrix();
  YC = transpose_view(X * YC);
  linsolve_Cholesky(S,YC);
  mat< ValueType, mat_structure::rectangular, mat_alignment::row_major > K(transpose_view(YC));
   
  b.set_mean_state( x + K * y );
  b.set_covariance( CovType( X, MatType( Y + M * X ) ) );
  
  SizeType N = X.get_row_count();
  Tm.set_row_count(2 * N);
  Tm.set_col_count(2 * N);
  
  set_block(Tm, mat<ValueType, mat_structure::identity>(N), 0, 0);
  set_block(Tm, mat<ValueType, mat_structure::nil>(N), 0, N);
  set_block(Tm, M, N, 0);
  set_block(Tm, mat<ValueType, mat_structure::identity>(N), N, N);
};



/**
 * This function template performs one complete estimation step using the Symplectic Kalman 
 * Filter method, which includes a prediction and measurement update step. This function is, 
 * in general, more efficient than applying the prediction and update separately.
 * \tparam LinearSystem A discrete state-space system modeling the DiscreteLinearSSSConcept 
 *         at least as a DiscreteLinearizedSystemType.
 * \tparam BeliefState A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation and a covariance matrix modeling the DecomposedCovarianceConcept.
 * \tparam SystemNoiseCovariance A covariance matrix type modeling the CovarianceMatrixConcept.
 * \tparam MeasurementNoiseCovariance A covariance matrix type modeling the CovarianceMatrixConcept.
 * \tparam CovTransMatrix A matrix type to store the covariance transformation matrix.
 * \param sys The discrete state-space system used in the state estimation.
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
 *        covarying and informing components of the covariance matrix of the belief-state.
 * \param t The current time (before the prediction).
 * 
 */
template <typename LinearSystem, 
          typename BeliefState, 
	  typename SystemNoiseCovariance,
	  typename MeasurementNoiseCovariance,
	  typename CovTransMatrix>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal) &&
                             is_fully_writable_matrix<CovTransMatrix>::value,
void >::type symplectic_kalman_filter_step(const LinearSystem& sys,
					   BeliefState& b,
					   const discrete_sss_traits<LinearSystem>::input_type& u,
					   const discrete_sss_traits<LinearSystem>::output_type& z,
					   const SystemNoiseCovariance& Q,
					   const MeasurementNoiseCovariance& R,
					   CovTransMatrix& T,
					   typename discrete_sss_traits<LinearSystem>::time_type t = 0) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) update
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) update
  
  typedef typename discrete_sss_traits<LinearSystem>::point_type StateType;
  typedef typename discrete_sss_traits<LinearSystem>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  
  BOOST_CONCEPT_ASSERT((DiscreteLinearSSSConcept< LinearSystem, DiscreteLinearizedSystemType >));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  BOOST_CONCEPT_ASSERT((CovarianceMatrixConcept<SystemNoiseCovariance>));
  BOOST_CONCEPT_ASSERT((CovarianceMatrixConcept<MeasurementNoiseCovariance>));
  BOOST_CONCEPT_ASSERT((WritableMatrixConcept<CovTransMatrix>));
  BOOST_CONCEPT_ASSERT((DecomposedCovarianceConcept<CovType>));
  
  typedef typename decomp_covariance_mat_traits< CovType >::matrix_block_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename mat_traits<MatType>::size_type SizeType;
  
  typename discrete_linear_sss_traits<LinearSystem>::matrixA_type A;
  typename discrete_linear_sss_traits<LinearSystem>::matrixB_type B;
  typename discrete_linear_sss_traits<LinearSystem>::matrixC_type C;
  typename discrete_linear_sss_traits<LinearSystem>::matrixD_type D;
  
  StateType x = b.get_mean_state();
  const MatType& X = b.get_covariance().get_covarying_block();
  const MatType& Y = b.get_covariance().get_informing_inv_block(); 
  
  x = sys.get_next_state(x,u,t);
  sys.get_state_transition_blocks(A, B, t, t + sys.get_time_step(), b.get_mean_state(), x, u, u);
  
  SizeType N = X.get_row_count();
  T.set_row_count(2 * N);
  T.set_col_count(2 * N);
  
  set_block(T, A, 0, 0);
  typename discrete_linear_sss_traits<LinearSystem>::matrixA_type A_inv; 
  pseudoinvert_QR(A,A_inv);
  mat_sub_block< CovTransMatrix > T_lr(T,N,N,N,N);
  T_lr = transpose_view(A_inv);
  mat_sub_block< CovTransMatrix > T_ur(T,N,N,0,N);
  T_ur = (B * Q.get_matrix() * transpose_view(B)) * T_lr;
  set_block(T, mat<ValueType,mat_structure::nil>(N), N, 0);
  
  X = A * X + T_ur * Y;
  Y = T_lr * Y;
  
  t += sys.get_time_step();
  
  sys.get_output_function_blocks(C, D, t, x, u);
  OutputType y = z - sys.get_output(x,u,t);
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > YC(transpose_view(C));
  mat< ValueType, mat_structure::symmetric > M = YC * R.get_inverse_matrix() * C;
  linsolve_QR(Y,YC);
  mat< ValueType, mat_structure::symmetric > S = C * X * YC + R.get_matrix();
  YC = transpose(X * YC);
  linsolve_Cholesky(S,YC);
  mat< ValueType, mat_structure::rectangular, mat_alignment::row_major > K(transpose_view(YC));
   
  b.set_mean_state( x + K * y );
  
  set_block(T, M * A, N, 0);
  T_lr += M * T_ur;
  
  b.set_covariance( CovType( X, MatType( Y + M * X ) ) );
};







/**
 * This class template can be used as a belief-state predictor (and transfer) that uses the 
 * Symplectic Kalman Filter method. This class template models the BeliefTransferConcept and 
 * the BeliefPredictorConcept.
 * \tparam LinearSystem A discrete state-space system modeling the DiscreteLinearSSSConcept 
 *         at least as a DiscreteLinearizedSystemType.
 * \tparam BeliefState A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam SystemNoiseCovar A covariance matrix type modeling the CovarianceMatrixConcept.
 * \tparam MeasurementCovar A covariance matrix type modeling the CovarianceMatrixConcept.
 */
template <typename LinearSystem,
          typename BeliefState = gaussian_belief_state< decomp_covariance_matrix< typename discrete_sss_traits<LinearSystem>::point_type > >,
          typename SystemNoiseCovar = covariance_matrix< typename discrete_sss_traits< LinearSystem >::input_type >,
          typename MeasurementCovar = covariance_matrix< typename discrete_sss_traits< LinearSystem >::output_type > >
struct SKF_belief_transfer {
  typedef SKF_belief_transfer<LinearSystem, BeliefState> self;
  typedef BeliefState belief_state;
  typedef LinearSystem state_space_system;
  typedef typename discrete_sss_traits< state_space_system >::time_type time_type;
  typedef typename discrete_sss_traits< state_space_system >::time_difference_type time_difference_type;

  typedef typename belief_state_traits< belief_state >::state_type state_type;
  typedef typename state_vector_traits< state_type >::value_type value_type;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type covariance_type;
  typedef typename decomp_covariance_mat_traits< covariance_type >::matrix_block_type matrix_type;
  typedef typename mat_traits< matrix_type >::value_type mat_value_type;
  typedef typename mat_traits< matrix_type >::size_type mat_size_type;

  typedef typename discrete_sss_traits< state_space_system >::input_type input_type;
  typedef typename discrete_sss_traits< state_space_system >::output_type output_type;

  BOOST_CONCEPT_ASSERT((DiscreteLinearSSSConcept< LinearSystem, DiscreteLinearizedSystemType >));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  BOOST_CONCEPT_ASSERT((CovarianceMatrixConcept<SystemNoiseCovar>));
  BOOST_CONCEPT_ASSERT((CovarianceMatrixConcept<MeasurementCovar>));
  BOOST_CONCEPT_ASSERT((DecomposedCovarianceConcept<covariance_type>));
  
  const LinearSystem* sys; ///< Holds the reference to the system used for the filter.
  SystemNoiseCovar Q; ///< Holds the system's input noise covariance matrix.
  MeasurementCovar R; ///< Holds the system's output measurement's covariance matrix.
  value_type reupdate_threshold; ///< The threshold at which the state change is considered too high and the state transition matrices are recomputed.
  state_type initial_state; ///< The initial mean-state at which the predictor is linearized (if non-linear).
  state_type predicted_state; ///< The predicted mean-state at which the updator is linearized (if non-linear). 
  mat< mat_value_type, mat_structure::square > Tc; ///< Holds the prediction covariance transformation matrix.
  mat< mat_value_type, mat_structure::square > Tm; ///< Holds the updating covariance transformation matrix.

  /**
   * Parametrized constructor.
   * \param aSys The reference to the system used for the filter.
   * \param aQ The system's input noise covariance matrix.
   * \param aR The system's output measurement's covariance matrix.
   * \param aReupdateThreshold The threshold at which the state change is considered too high and the state transition matrices are recomputed.
   * \param aInitialState The initial mean-state at which the predictor is linearized (if non-linear).
   * \param aInitialInput The input at the time of the initial state.
   * \param aInitialTime The time of the initial state.
   */
  SKF_belief_transfer(const LinearSystem& aSys, 
                      const SystemNoiseCovar& aQ,
                      const MeasurementCovar& aR,
                      const value_type& aReupdateThreshold,
                      const state_type& aInitialState,
                      const input_type& aInitialInput,
                      const time_type& aInitialTime) : 
                      sys(&aSys), Q(aQ), R(aR),
                      reupdate_threshold(aReupdateThreshold),
                      initial_state(aInitialState) {
    symplectic_kalman_predict(*sys,
                              b,
                              aInitialInput,
                              Q,
                              Tc,
                              aInitialTime);
    predicted_state = b.get_mean_state();
    symplectic_kalman_update(*sys,
                             b,
                             aInitialInput,
                             sys->get_output(predicted_state,
                                             aInitialInput,
                                             aInitialTime),
                             R,
                             Tm,
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
   * \param b The current belief-state.
   * \param t The current time.
   * \param u The current input given to the system.
   * \param y The output that was measured at the next time instant.
   * \return the belief-state at the next time instant.
   */
  belief_state get_next_belief(belief_state b, const time_type& t, const input_type& u, const input_type& y) {
    initial_state = b.get_mean_state();
    symplectic_kalman_predict(*sys,b,u,Q,Tc,t);
    predicted_state = b.get_mean_state();
    symplectic_kalman_update(*sys,b,u,y,R,Tm,t);
    return b;
  };
  
  /**
   * Returns the prediction belief-state at the next time instant.
   * \param b The current belief-state.
   * \param t The current time.
   * \param u The current input given to the system.
   * \return the belief-state at the next time instant, predicted by the filter.
   */
  belief_state predict_belief(belief_state b, const time_type& t, const input_type& u) {
    if( norm( diff(b.get_mean_state(), initial_state) ) > reupdate_threshold ) {
      initial_state = b.get_mean_state();
      symplectic_kalman_predict(*sys,b,u,Q,Tc,t);
    } else {
      b.set_mean_state( sys->get_next_state(b.get_mean_state(), u, t) );
      mat_size_type N = Tc.get_row_count() / 2;
      mat< mat_value_type, mat_structure::rectangular > P_tmp = 
        Tc * ( b.get_covariance().get_covarying_block() |
               b.get_covariance().get_informing_inv_block() );
      b.set_covariance( covariance_type( matrix_type( sub(P_tmp)(range(0,N-1),range(0,N-1)) ), 
                                         matrix_type( sub(P_tmp)(range(N,2*N-1),range(0,N-1)) ) ) ); 
    };
    return b;
  };
  
  /**
   * Converts a prediction belief-state into an updated belief-state which assumes the most likely measurement.
   * \param b The current prediction's belief-state.
   * \param t The current time.
   * \param u The current input given to the system.
   * \return the updated belief-state when assuming the most likely measurement.
   */
  belief_state prediction_to_ML_belief(belief_state b, const time_type& t, const input_type& u) {
    if( norm( diff(b.get_mean_state(), predicted_state) ) > reupdate_threshold ) {
      predicted_state = b.get_mean_state();
      symplectic_kalman_update(*sys,b,u,sys->get_output(predicted_state,u,t + sys->get_time_step()),R,Tm,t);
    } else {
      mat_size_type N = Tm.get_row_count() / 2;
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
   * \param b The current belief-state.
   * \param t The current time.
   * \param u The current input given to the system.
   * \return the belief-state at the next time instant, predicted by the filter.
   */
  belief_state predict_ML_belief(belief_state b, const time_type& t, const input_type& u) {
    return prediction_to_ML_belief(predict_belief(b, t, u),t,u);
  };
  
};









};

};

#endif












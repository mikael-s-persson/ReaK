/**
 * \file aggregate_kalman_filter.hpp
 * 
 * This library provides a number of functions and classes to do state estimation 
 * using the Aggregate Kalman Filter. This Kalman filtering technique applies to a 
 * gaussian belief state where the covariance is decomposed into a covarying matrix
 * and an informing matrix (i.e. P = X * invert(Y)). The transition between covariances
 * is achieved using the aggregation of hamiltonian matrices (akin to scattering matrices)
 * via the Redeffer Star-product. This aggregation of matrices extends beyond a single 
 * estimation step and can be aggregated over several steps. If the system is non-linear,
 * then it would need the recomputation of those hamiltonian matrices for the individual
 * steps if the mean-states change too much. The estimation functions will output the 
 * hamiltonian matrices, allowing the caller to aggregate them if needed. If the aggregation
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

#ifndef REAK_AGGREGATE_KALMAN_FILTER_HPP
#define REAK_AGGREGATE_KALMAN_FILTER_HPP

#include "belief_state_concept.hpp"
#include "discrete_linear_sss_concept.hpp"
#include <boost/utility/enable_if.hpp>
#include <lin_alg/vect_concepts.hpp>
#include <lin_alg/mat_alg.hpp>
#include <lin_alg/mat_cholesky.hpp>

#include <boost/static_assert.hpp>
#include "covariance_concept.hpp"

#include "covariance_matrix.hpp"
#include "decomp_covariance_matrix.hpp"
#include "gaussian_belief_state.hpp"

#include "lin_alg/mat_star_product.hpp"

namespace ReaK {

namespace ctrl {




/**
 * This function template performs one prediction step using the Aggregate Kalman Filter method.
 * \tparam LinearSystem A discrete state-space system modeling the DiscreteLinearSSSConcept 
 *         at least as a DiscreteLinearizedSystemType.
 * \tparam StateSpaceType A topology type on which the state-vectors can reside, should model
 *         the pp::TopologyConcept.
 * \tparam BeliefState A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam InputBelief A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \param sys The discrete state-space system used in the state estimation.
 * \param state_space The state-space topology on which the state representations lie.
 * \param b_x As input, it stores the belief-state before the prediction step. As output, it stores
 *        the belief-state after the prediction step.
 * \param b_u The input belief to apply to the state-space system to make the transition of the 
 *        mean-state, i.e., the current input vector and its covariance.
 * \param Sc Stores, as output, the scattering matrix (a hamiltonian matrix) which can scatter the 
 *        covarying and informing components of the covariance matrix of the belief-state.
 * \param t The current time (before the prediction).
 * 
 */
template <typename LinearSystem,  
          typename StateSpaceType,
          typename BeliefState, 
	  typename InputBelief>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal),
void >::type aggregate_kalman_predict(const LinearSystem& sys,
			              const StateSpaceType& state_space,
			              BeliefState& b_x,
			              const InputBelief& b_u,
			              typename hamiltonian_mat< typename mat_traits< typename covariance_mat_traits< typename continuous_belief_state_traits<BeliefState>::covariance_type >::matrix_type >::value_type >::type& Sc,
				      typename discrete_sss_traits<LinearSystem>::time_type t = 0) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) update
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) update
  typedef typename pp::topology_traits<StateSpaceType>::point_type StateType;
  typedef typename discrete_sss_traits<LinearSystem>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  typedef typename covariance_mat_traits< CovType >::matrix_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename mat_traits<MatType>::size_type SizeType;
  
  BOOST_CONCEPT_ASSERT((pp::TopologyConcept< StateSpaceType >));
  BOOST_CONCEPT_ASSERT((DiscreteLinearSSSConcept< LinearSystem, StateSpaceType, DiscreteLinearizedSystemType >));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<InputBelief>));
  
  typedef typename hamiltonian_mat< ValueType >::type HamilMat;
  typedef typename hamiltonian_mat< ValueType >::upper HamilMatUp;
  typedef typename hamiltonian_mat< ValueType >::lower HamilMatLo;
  typedef typename hamiltonian_mat< ValueType >::upper_left HamilMatUL;
  typedef typename hamiltonian_mat< ValueType >::upper_right HamilMatUR;
  typedef typename hamiltonian_mat< ValueType >::lower_left HamilMatLL;
  typedef typename hamiltonian_mat< ValueType >::lower_right HamilMatLR;
  
  using std::swap;
  
  typename discrete_linear_sss_traits<LinearSystem>::matrixA_type A;
  typename discrete_linear_sss_traits<LinearSystem>::matrixB_type B;
  
  StateType x = b_x.get_mean_state();
  MatType P = b_x.get_covariance().get_matrix();

  StateType x_prior = sys.get_next_state(state_space, x, b_u.get_mean_state(), t);
  sys.get_state_transition_blocks(A, B, state_space, t, t + sys.get_time_step(), x, x_prior, b_u.get_mean_state(), b_u.get_mean_state());
  SizeType N = A.get_row_count();
  HamilMatUR Q_tmp(B * b_u.get_covariance().get_matrix() * transpose_view(B));
  HamilMat Sc_tmp(HamilMatUp(HamilMatUL(A),Q_tmp),
                  HamilMatLo(HamilMatLL(mat<ValueType,mat_structure::nil>(N)),HamilMatLR(transpose_view(A))));
  P = A * P * transpose_view(A) + Q_tmp;
  b_x.set_mean_state( x_prior );
  b_x.set_covariance( CovType( P ) );
  swap(Sc,Sc_tmp);
};




/**
 * This function template performs one measurement update step using the Aggregate Kalman Filter method.
 * \tparam LinearSystem A discrete state-space system modeling the DiscreteLinearSSSConcept 
 *         at least as a DiscreteLinearizedSystemType.
 * \tparam StateSpaceType A topology type on which the state-vectors can reside, should model
 *         the pp::TopologyConcept.
 * \tparam BeliefState A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam InputBelief A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam MeasurementBelief A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \param sys The discrete state-space system used in the state estimation.
 * \param state_space The state-space topology on which the state representations lie.
 * \param b_x As input, it stores the belief-state before the update step. As output, it stores
 *        the belief-state after the update step.
 * \param b_u The input vector to apply to the state-space system to make the transition of the 
 *        mean-state, i.e., the current input vector and its covariance.
 * \param b_z The output belief that was measured, i.e. the measurement vector and its covariance.
 * \param Sm Stores, as output, the scattering matrix (a hamiltonian matrix) which can scatter the 
 *        covarying and informing components of the covariance matrix of the belief-state.
 * \param t The current time.
 * 
 */
template <typename LinearSystem, 
          typename StateSpaceType,
          typename BeliefState, 
	  typename InputBelief, 
	  typename MeasurementBelief>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal),
void >::type aggregate_kalman_update(const LinearSystem& sys,
				     const StateSpaceType& state_space,
				     BeliefState& b_x,
				     const InputBelief& b_u,
				     const MeasurementBelief& b_z,
				     typename hamiltonian_mat< typename mat_traits< typename covariance_mat_traits< typename continuous_belief_state_traits<BeliefState>::covariance_type >::matrix_type >::value_type >::type& Sm,
				     typename discrete_sss_traits<LinearSystem>::time_type t = 0) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) update
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) update
  BOOST_CONCEPT_ASSERT((pp::TopologyConcept< StateSpaceType >));
  BOOST_CONCEPT_ASSERT((DiscreteLinearSSSConcept< LinearSystem, StateSpaceType, DiscreteLinearizedSystemType >));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<InputBelief>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<MeasurementBelief>));

  typedef typename discrete_sss_traits<LinearSystem>::point_type StateType;
  typedef typename discrete_sss_traits<LinearSystem>::point_difference_type StateDiffType;
  typedef typename discrete_sss_traits<LinearSystem>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  typedef typename covariance_mat_traits< CovType >::matrix_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename mat_traits<MatType>::size_type SizeType;
  
  typedef typename hamiltonian_mat< ValueType >::type HamilMat;
  typedef typename hamiltonian_mat< ValueType >::upper HamilMatUp;
  typedef typename hamiltonian_mat< ValueType >::lower HamilMatLo;
  typedef typename hamiltonian_mat< ValueType >::upper_left HamilMatUL;
  typedef typename hamiltonian_mat< ValueType >::upper_right HamilMatUR;
  typedef typename hamiltonian_mat< ValueType >::lower_left HamilMatLL;
  typedef typename hamiltonian_mat< ValueType >::lower_right HamilMatLR;
  
  using std::swap;
  using ReaK::to_vect; using ReaK::from_vect;
  
  typename discrete_linear_sss_traits<LinearSystem>::matrixC_type C;
  typename discrete_linear_sss_traits<LinearSystem>::matrixD_type D;
  
  StateType x = b_x.get_mean_state();
  MatType P = b_x.get_covariance().get_matrix();
  sys.get_output_function_blocks(C, D, state_space, t, x, b_u.get_mean_state());
  SizeType N = C.get_col_count();

  OutputType y = b_z - sys.get_output(state_space, x, b_u.get_mean_state(), t);
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > CP = C * P;
  mat< ValueType, mat_structure::symmetric > S = CP * transpose_view(C) + b_z.get_covariance().get_matrix();
  linsolve_Cholesky(S,CP);
  mat< ValueType, mat_structure::rectangular, mat_alignment::row_major > K( transpose_view(CP) );
   
  b_x.set_mean_state( state_space.adjust( x, from_vect<StateDiffType>( K * to_vect<ValueType>(y) ) ) );
  b_x.set_covariance( CovType( MatType( (mat< ValueType, mat_structure::identity>(K.get_row_count()) - K * C) * P ) ) );

  HamilMat Sm_tmp(HamilMatUp(HamilMatUL(mat<ValueType,mat_structure::identity>(N)),HamilMatUR(mat<ValueType,mat_structure::nil>(N))),HamilMatLo(HamilMatLL( transpose_view(C) * b_z.get_covariance().get_inverse_matrix() * C ),HamilMatLR(mat<ValueType,mat_structure::identity>(N))));
  swap(Sm,Sm_tmp);
};




/**
 * This function template performs one complete estimation step using the Aggregate Kalman 
 * Filter method, which includes a prediction and measurement update step. This function is, 
 * in general, more efficient than applying the prediction and update separately.
 * \tparam LinearSystem A discrete state-space system modeling the DiscreteLinearSSSConcept 
 *         at least as a DiscreteLinearizedSystemType.
 * \tparam StateSpaceType A topology type on which the state-vectors can reside, should model
 *         the pp::TopologyConcept.
 * \tparam BeliefState A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam InputBelief A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam MeasurementBelief A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \param sys The discrete state-space system used in the state estimation.
 * \param state_space The state-space topology on which the state representations lie.
 * \param b_x As input, it stores the belief-state before the estimation step. As output, it stores
 *        the belief-state after the estimation step.
 * \param b_u The input vector to apply to the state-space system to make the transition of the 
 *        mean-state, i.e., the current input vector and its covariance.
 * \param b_z The output belief that was measured, i.e. the measurement vector and its covariance.
 * \param Sc Stores, as output, the scattering matrix (a hamiltonian matrix) which can scatter the 
 *        covarying and informing components of the covariance matrix of the belief-state. This 
 *        scattering matrix is associated to the prediction step.
 * \param Sm Stores, as output, the scattering matrix (a hamiltonian matrix) which can scatter the 
 *        covarying and informing components of the covariance matrix of the belief-state. This 
 *        scattering matrix is associated to the measurement update step.
 * \param t The current time (before the prediction).
 * 
 */
template <typename LinearSystem, 
          typename StateSpaceType,
          typename BeliefState, 
	  typename InputBelief, 
	  typename MeasurementBelief>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal),
void >::type aggregate_kalman_filter_step(const LinearSystem& sys,
			                  const StateSpaceType& state_space,
					  BeliefState& b_x,
					  const InputBelief& b_u,
					  const MeasurementBelief& b_z,
					  typename hamiltonian_mat< typename mat_traits< typename covariance_mat_traits< typename continuous_belief_state_traits<BeliefState>::covariance_type >::matrix_type >::value_type >::type& Sc,
					  typename hamiltonian_mat< typename mat_traits< typename covariance_mat_traits< typename continuous_belief_state_traits<BeliefState>::covariance_type >::matrix_type >::value_type >::type& Sm,
					  typename discrete_sss_traits<LinearSystem>::time_type t = 0) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) update
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) update
  BOOST_CONCEPT_ASSERT((pp::TopologyConcept< StateSpaceType >));
  BOOST_CONCEPT_ASSERT((DiscreteLinearSSSConcept< LinearSystem, StateSpaceType, DiscreteLinearizedSystemType >));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<InputBelief>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<MeasurementBelief>));
  
  typedef typename discrete_sss_traits<LinearSystem>::point_type StateType;
  typedef typename discrete_sss_traits<LinearSystem>::point_difference_type StateDiffType;
  typedef typename discrete_sss_traits<LinearSystem>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  typedef typename covariance_mat_traits< CovType >::matrix_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename mat_traits<MatType>::size_type SizeType;
  
  typedef typename hamiltonian_mat< ValueType >::type HamilMat;
  typedef typename hamiltonian_mat< ValueType >::upper HamilMatUp;
  typedef typename hamiltonian_mat< ValueType >::lower HamilMatLo;
  typedef typename hamiltonian_mat< ValueType >::upper_left HamilMatUL;
  typedef typename hamiltonian_mat< ValueType >::upper_right HamilMatUR;
  typedef typename hamiltonian_mat< ValueType >::lower_left HamilMatLL;
  typedef typename hamiltonian_mat< ValueType >::lower_right HamilMatLR;
  
  using std::swap;
  using ReaK::to_vect; using ReaK::from_vect;
  
  typename discrete_linear_sss_traits<LinearSystem>::matrixA_type A;
  typename discrete_linear_sss_traits<LinearSystem>::matrixB_type B;
  typename discrete_linear_sss_traits<LinearSystem>::matrixC_type C;
  typename discrete_linear_sss_traits<LinearSystem>::matrixD_type D;

  StateType x = b_x.get_mean_state();
  MatType P = b_x.get_covariance().get_matrix();

  StateType x_prior = sys.get_next_state(state_space, x, b_u.get_mean_state(), t);
  sys.get_state_transition_blocks(A, B, state_space, t, t + sys.get_time_step(), x, x_prior, b_u.get_mean_state(), b_u.get_mean_state());
  SizeType N = A.get_row_count();
  HamilMatUR Q_tmp(B * b_u.get_covariance().get_matrix() * transpose_view(B));
  HamilMat Sc_tmp(HamilMatUp(HamilMatUL(A),Q_tmp),
                  HamilMatLo(HamilMatLL(mat<ValueType,mat_structure::nil>(N)),HamilMatLR(transpose_view(A))));
  swap(Sc,Sc_tmp);
  P = A * P * transpose_view(A) + Q_tmp;
  
  sys.get_output_function_blocks(C, D, state_space, t + sys.get_time_step(), x_prior, b_u.get_mean_state());
  OutputType y = b_z.get_mean_state() - sys.get_output(state_space, x_prior, b_u.get_mean_state(), t + sys.get_time_step());
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > CP = C * P;
  mat< ValueType, mat_structure::symmetric > S = CP * transpose_view(C) + b_z.get_covariance().get_matrix();
  linsolve_Cholesky(S,CP);
  mat< ValueType, mat_structure::rectangular, mat_alignment::row_major > K( transpose_view(CP) );
   
  b_x.set_mean_state( state_space.adjust(x_prior, from_vect<StateDiffType>( K * to_vect<ValueType>(y) ) ) );
  b_x.set_covariance( CovType( MatType( (mat< ValueType, mat_structure::identity>(K.get_row_count()) - K * C) * P ) ) );

  HamilMat Sm_tmp(HamilMatUp(HamilMatUL(mat<ValueType,mat_structure::identity>(N)),HamilMatUR(mat<ValueType,mat_structure::nil>(N))),HamilMatLo(HamilMatLL( transpose_view(C) * b_z.get_covariance().get_inverse_matrix() * C ),HamilMatLR(mat<ValueType,mat_structure::identity>(N))));
  swap(Sm,Sm_tmp);
};





/**
 * This class template can be used as a belief-state predictor (and transfer) that uses the 
 * Aggregate Kalman Filter method. This class template models the BeliefTransferConcept and 
 * the BeliefPredictorConcept.
 * \tparam LinearSystem A discrete state-space system modeling the DiscreteLinearSSSConcept 
 *         at least as a DiscreteLinearizedSystemType.
 * \tparam BeliefState A belief state type modeling the ContinuousBeliefStateConcept with
 *         a unimodular gaussian representation.
 * \tparam SystemNoiseCovar A covariance matrix type modeling the CovarianceMatrixConcept.
 * \tparam MeasurementCovar A covariance matrix type modeling the CovarianceMatrixConcept.
 */
template <typename LinearSystem,
          typename BeliefState = gaussian_belief_state<typename discrete_sss_traits<LinearSystem>::point_type, decomp_covariance_matrix< typename discrete_sss_traits<LinearSystem>::point_type > >,
          typename SystemNoiseCovar = covariance_matrix< typename discrete_sss_traits< LinearSystem >::input_type >,
          typename MeasurementCovar = covariance_matrix< typename discrete_sss_traits< LinearSystem >::output_type > >
struct AKF_belief_transfer {
  
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  
  typedef AKF_belief_transfer<LinearSystem, BeliefState, SystemNoiseCovar, MeasurementCovar> self;
  typedef BeliefState belief_state;
  typedef LinearSystem state_space_system;
  typedef shared_ptr< LinearSystem > state_space_system_ptr;
  typedef typename discrete_sss_traits< state_space_system >::time_type time_type;
  typedef typename discrete_sss_traits< state_space_system >::time_difference_type time_difference_type;

  typedef typename belief_state_traits< belief_state >::state_type state_type;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type covariance_type;
  typedef typename covariance_mat_traits< covariance_type >::matrix_type matrix_type;
  typedef typename mat_traits< matrix_type >::value_type value_type;
  typedef typename mat_traits< matrix_type >::size_type mat_size_type;

  typedef typename discrete_sss_traits< state_space_system >::input_type input_type;
  typedef typename discrete_sss_traits< state_space_system >::output_type output_type;
  
  typedef gaussian_belief_state<input_type, SystemNoiseCovar> input_belief_type;
  typedef gaussian_belief_state<output_type, MeasurementCovar> output_belief_type;

  state_space_system_ptr sys; ///< Holds the reference to the system used for the filter.
  SystemNoiseCovar Q; ///< Holds the system's input noise covariance matrix.
  MeasurementCovar R; ///< Holds the system's output measurement's covariance matrix.
  value_type reupdate_threshold; ///< The threshold at which the state change is considered too high and the state transition matrices are recomputed.
  state_type initial_state; ///< The initial mean-state at which the predictor is linearized (if non-linear).
  state_type predicted_state; ///< The predicted mean-state at which the updator is linearized (if non-linear). 
  typename hamiltonian_mat< mat_value_type >::type Sc; ///< Holds the prediction scattering matrix.
  typename hamiltonian_mat< mat_value_type >::type Sm; ///< Holds the updating scattering matrix.

  /**
   * Parametrized constructor.
   * \param aSys The pointer to the system used for the filter.
   * \param state_space The state-space in which to operate.
   * \param aQ The system's input noise covariance matrix.
   * \param aR The system's output measurement's covariance matrix.
   * \param aReupdateThreshold The threshold at which the state change is considered too high and the state transition matrices are recomputed.
   * \param aInitialState The initial mean-state at which the predictor is linearized (if non-linear).
   * \param aInitialInput The input at the time of the initial state.
   * \param aInitialTime The time of the initial state.
   */
  template <typename StateSpaceType>
  AKF_belief_transfer(const state_space_system_ptr& aSys, 
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
    aggregate_kalman_predict(*sys,
			     state_space,
                             b,
                             b_u,
                             Sc,
                             aInitialTime);
    predicted_state = b.get_mean_state();
    aggregate_kalman_update(*sys,
			    state_space,
                            b,
                            b_u,
                            sys->get_output(state_space,
			                    predicted_state,
                                            aInitialInput,
                                            aInitialTime + sys->get_time_step()),
                            R,
                            Sm,
                            aInitialTime + sys->get_time_step());
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
    output_belief_type b_z(y, R);
    aggregate_kalman_predict(*sys,b_space.get_state_topology(),b,b_u,Sc,t);
    predicted_state = b.get_mean_state();
    aggregate_kalman_update(*sys,b_space.get_state_topology(),b,b_u,b_z,Sm,t);
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
      aggregate_kalman_predict(*sys,b_space.get_state_topology(),b,b_u,Sc,t);
    } else {
      b.set_mean_state( sys->get_next_state(b_space.get_state_topology(), b.get_mean_state(), u, t));
      mat_size_type N = Sc.get_row_count() / 2;
      typename hamiltonian_mat< mat_value_type >::type P_h =
        star_product( ( mat<value_type,mat_structure::identity>(N) & b.get_covariance().get_matrix()             |
                        mat<value_type,mat_structure::nil>(N)      & mat<value_type,mat_structure::identity>(N)  ) , 
                      Sc );
      b.set_covariance( covariance_type( matrix_type( sub(P_h)(range(0,N-1),range(N,2*N-1)) ) ) ); 
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
      output_belief_type b_z(sys->get_output(b_space.get_state_topology(),predicted_state,u,t + sys->get_time_step()), R);
      aggregate_kalman_update(*sys,b_space.get_state_topology(),b,b_u,b_z,Sm,t);
    } else {
      mat_size_type N = Sm.get_row_count() / 2;
      typename hamiltonian_mat< mat_value_type >::type P_h =
        star_product( ( mat<value_type,mat_structure::identity>(N) & b.get_covariance().get_matrix()             |
                        mat<value_type,mat_structure::nil>(N)      & mat<value_type,mat_structure::identity>(N)  ) , 
                      Sm );
      b.set_covariance( covariance_type( matrix_type( sub(P_h)(range(0,N-1),range(N,2*N-1)) ) ) );
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













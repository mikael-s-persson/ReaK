
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

#ifndef INVARIANT_SYMPLECTIC_KALMAN_FILTER_HPP
#define INVARIANT_SYMPLECTIC_KALMAN_FILTER_HPP

#include "belief_state_concept.hpp"
#include "discrete_linear_sss_concept.hpp"
#include "invariant_system_concept.hpp"
#include <boost/utility/enable_if.hpp>
#include "math/vect_concepts.hpp"
#include "math/mat_alg.hpp"

#include <boost/static_assert.hpp>
#include "covariance_concept.hpp"
#include "math/mat_cholesky.hpp"
#include "math/mat_qr_decomp.hpp"
#include "math/mat_views.hpp"

namespace ReaK {

namespace ctrl {






template <typename InvariantSystem, 
          typename BeliefState, 
	  typename SystemNoiseCovariance,
	  typename InvCovTransMatrix,
	  typename InvFrameMatrix>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal) &&
                             is_fully_writable_matrix<InvCovTransMatrix>::value &&
                             is_writable_matrix<InvFrameMatrix>::value,
void >::type invariant_symplectic_kf_predict(const InvariantSystem& sys,
					     BeliefState& b,
					     const discrete_sss_traits<InvariantSystem>::input_type& u,
					     const SystemNoiseCovariance& Q,
					     InvCovTransMatrix& Tc,
					     InvFrameMatrix& Wp,
					     typename discrete_sss_traits<InvariantSystem>::time_type t = 0) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) update
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) update
  boost::function_requires< DiscreteLinearSSSConcept< InvariantSystem, DiscreteLinearizedSystemType > >();
  boost::function_requires< InvariantDiscreteSystemConcept<InvariantSystem> >();
  boost::function_requires< ContinuousBeliefStateConcept<BeliefState> >();

  typedef typename discrete_sss_traits<InvariantSystem>::point_type StateType;
  typedef typename discrete_sss_traits<InvariantSystem>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  typedef typename invariant_system_traits<InvariantSystem>::invariant_error_type ErrorType;
  
  boost::function_requires< DecomposedCovarianceConcept<CovType> >();
  
  typedef typename decomp_covariance_mat_traits< CovType >::matrix_block_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename mat_traits<MatType>::size_type SizeType;
  
  typename discrete_linear_sss_traits<InvariantSystem>::matrixA_type A;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixB_type B;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixC_type C;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixD_type D;
  
  StateType x = b.get_mean_state();
  const MatType& X = b.get_covariance().get_covarying_block();
  const MatType& Y = b.get_covariance().get_informing_inv_block(); 
  sys.get_linear_blocks(A, B, C, D, t, x, u);
  
  x = sys.get_next_state(x,u,t);
  Wp = sys.get_invariant_prior_frame(b.get_mean_state(), x, u, t + sys.get_time_step());
  
  SizeType N = A.get_col_count();
  Tc.set_row_count(2 * N);
  Tc.set_col_count(2 * N);
  
  set_block(Tc, A, 0, 0);
  typename discrete_linear_sss_traits<InvariantSystem>::matrixA_type A_inv; 
  pseudoinvert_QR(A,A_inv);
  mat_sub_block< InvCovTransMatrix > T_lr(Tc,N,N,N,N);
  T_lr = transpose_move(A_inv);
  mat_sub_block< InvCovTransMatrix > T_ur(Tc,N,N,0,N);
  T_ur = (B * Q.get_matrix() * transpose(B) ) * T_lr;
  set_block(Tc, mat<ValueType,mat_structure::nil>(N), N, 0);
  
  b.set_covariance( CovType( MatType( Wp * ( A * X + T_ur * Y ) ), MatType( Wp * T_lr * Y ) ) );
  b.set_mean_state(x);
};







template <typename InvariantSystem, 
          typename BeliefState, 
	  typename MeasurementNoiseCovariance,
	  typename InvCovTransMatrix,
	  typename InvFrameMatrix>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal) &&
                             is_fully_writable_matrix<InvCovTransMatrix>::value &&
                             is_writable_matrix<InvFrameMatrix>::value,
void >::type invariant_symplectic_kf_update(const InvariantSystem& sys,
					    BeliefState& b,
					    const discrete_sss_traits<InvariantSystem>::input_type& u,
					    const discrete_sss_traits<InvariantSystem>::output_type& z,
					    const MeasurementNoiseCovariance& R,
					    InvCovTransMatrix& Tm,
					    InvFrameMatrix& Wu,
					    typename discrete_sss_traits<InvariantSystem>::time_type t = 0) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) update
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) update
  boost::function_requires< DiscreteLinearSSSConcept< InvariantSystem, DiscreteLinearizedSystemType > >();
  boost::function_requires< InvariantDiscreteSystemConcept<InvariantSystem> >();
  boost::function_requires< ContinuousBeliefStateConcept<BeliefState> >();

  typedef typename discrete_sss_traits<InvariantSystem>::point_type StateType;
  typedef typename discrete_sss_traits<InvariantSystem>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  typedef typename invariant_system_traits<InvariantSystem>::invariant_error_type ErrorType;
  
  boost::function_requires< DecomposedCovarianceConcept<CovType> >();
  
  typedef typename decomp_covariance_mat_traits< CovType >::matrix_block_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename mat_traits<MatType>::size_type SizeType;
  
  typename discrete_linear_sss_traits<InvariantSystem>::matrixA_type A;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixB_type B;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixC_type C;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixD_type D;
  
  StateType x = b.get_mean_state();
  const MatType& X = b.get_covariance().get_covarying_block();
  const MatType& Y = b.get_covariance().get_informing_inv_block(); 
  sys.get_linear_blocks(A, B, C, D, t, x, u);
  
  SizeType N = A.get_col_count();
  Tm.set_row_count(2 * N);
  Tm.set_col_count(2 * N);
  
  ErrorType e = sys.get_output_error(x, u, z, t);
  
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > Ct = transpose(C);
  mat< ValueType, mat_structure::symmetric > M = Ct * R.get_inverse_matrix() * C;
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > YC;
  linlsq_QR(Y,YC,Ct);
  mat< ValueType, mat_structure::symmetric > S = C * X * YC + R.get_matrix();
  YC = transpose(X * YC);
  linsolve_Cholesky(S,YC);
  mat< ValueType, mat_structure::rectangular, mat_alignment::row_major > K = transpose_move(YC);
   
  b.set_mean_state( sys.apply_correction(x, K * e, u, t) );
  Wu = sys.get_invariant_posterior_frame(x, b.get_mean_state(), u, t);
  
  set_block(Tm, M, N, 0);
  set_block(Tm, mat< ValueType, mat_structure::identity>(N), 0, 0);
  set_block(Tm, mat< ValueType, mat_structure::identity>(N), N, N);
    
  b.set_covariance( CovType( MatType( Wu * X ), MatType( Wu * ( Y + M * X ) ) ) );
};





template <typename InvariantSystem, 
          typename BeliefState, 
	  typename SystemNoiseCovariance,
	  typename MeasurementNoiseCovariance,
	  typename InvCovTransMatrix,
	  typename InvFrameMatrix>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal) &&
                             is_fully_writable_matrix<InvCovTransMatrix>::value &&
                             is_writable_matrix<InvFrameMatrix>::value,
void >::type invariant_symplectic_kf_step(const InvariantSystem& sys,
					  BeliefState& b,
					  const discrete_sss_traits<InvariantSystem>::input_type& u,
					  const discrete_sss_traits<InvariantSystem>::output_type& z,
					  const SystemNoiseCovariance& Q,
					  const MeasurementNoiseCovariance& R,
					  InvCovTransMatrix& T,
					  InvFrameMatrix& W,
					  typename discrete_sss_traits<InvariantSystem>::time_type t = 0) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) update
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) update
  boost::function_requires< DiscreteLinearSSSConcept< InvariantSystem, DiscreteLinearizedSystemType > >();
  boost::function_requires< InvariantDiscreteSystemConcept<InvariantSystem> >();
  boost::function_requires< ContinuousBeliefStateConcept<BeliefState> >();

  typedef typename discrete_sss_traits<InvariantSystem>::point_type StateType;
  typedef typename discrete_sss_traits<InvariantSystem>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  typedef typename invariant_system_traits<InvariantSystem>::invariant_error_type ErrorType;
  
  boost::function_requires< DecomposedCovarianceConcept<CovType> >();
  
  typedef typename decomp_covariance_mat_traits< CovType >::matrix_block_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename mat_traits<MatType>::size_type SizeType;
  
  typename discrete_linear_sss_traits<InvariantSystem>::matrixA_type A;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixB_type B;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixC_type C;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixD_type D;
  
  StateType x = b.get_mean_state();
  MatType X = b.get_covariance().get_covarying_block();
  MatType Y = b.get_covariance().get_informing_inv_block(); 
  sys.get_linear_blocks(A, B, C, D, t, x, u);
  
  x = sys.get_next_state(x,u,t);
  W = sys.get_invariant_prior_frame(b.get_mean_state(), x, u, t + sys.get_time_step());
  
  SizeType N = A.get_col_count();
  T.set_row_count(2 * N);
  T.set_col_count(2 * N);
  
  set_block(T, A, 0, 0);
  typename discrete_linear_sss_traits<InvariantSystem>::matrixA_type A_inv; 
  pseudoinvert_QR(A,A_inv);
  mat_sub_block< InvCovTransMatrix > T_lr(T,N,N,N,N);
  T_lr = transpose_move(A_inv);
  mat_sub_block< InvCovTransMatrix > T_ur(T,N,N,0,N);
  T_ur = (B * Q.get_matrix() * transpose(B) ) * T_lr;
  set_block(T, mat<ValueType,mat_structure::nil>(N), N, 0);
  
  X = A * X + T_ur * Y;
  Y = T_lr * Y;
  
  ErrorType e = sys.get_output_error(x, u, z, t + sys.get_time_step());
  
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > Ct = transpose(C);
  mat< ValueType, mat_structure::symmetric > M = Ct * R.get_inverse_matrix() * C;
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > YC;
  linlsq_QR(Y,YC,Ct);
  mat< ValueType, mat_structure::symmetric > S = C * X * YC + R.get_matrix();
  YC = transpose(X * YC);
  linsolve_Cholesky(S,YC);
  mat< ValueType, mat_structure::rectangular, mat_alignment::row_major > K = transpose_move(YC);
   
  b.set_mean_state( sys.apply_correction(x, W * K * e, u, t + sys.get_time_step()) );
  W = sys.get_invariant_posterior_frame(x, b.get_mean_state(), u, t + sys.get_time_step()) * W;
  
  set_block(T, M * A, N, 0);
  T_lr += M * T_ur;
    
  b.set_covariance( CovType( MatType( W * X ), MatType( W * ( Y + M * X ) ) ) );
};









template <typename LinearSystem,
          typename BeliefState = gaussian_belief_state< decomp_covariance_matrix< typename discrete_sss_traits<LinearSystem>::point_type > >,
          typename SystemNoiseCovar = covariance_matrix< typename discrete_sss_traits< LinearSystem >::input_type >,
          typename MeasurementCovar = covariance_matrix< typename discrete_sss_traits< LinearSystem >::output_type > >
struct ISKF_belief_transfer {
  typedef ISKF_belief_transfer<LinearSystem, BeliefState> self;
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

  const LinearSystem* sys;
  SystemNoiseCovar Q;
  MeasurementCovar R;
  value_type reupdate_threshold;
  state_type initial_state;
  state_type predicted_state;
  mat< mat_value_type, mat_structure::square > Tc;
  mat< mat_value_type, mat_structure::square > Wp;
  mat< mat_value_type, mat_structure::square > Tm;
  mat< mat_value_type, mat_structure::square > Wu;

  ISKF_belief_transfer(const LinearSystem& aSys, 
                       const SystemNoiseCovar& aQ,
                       const MeasurementCovar& aR,
                       const value_type& aReupdateThreshold,
                       const state_type& aInitialState,
                       const input_type& aInitialInput,
                       const time_type& aInitialTime) : 
                       sys(&aSys), Q(aQ), R(aR),
                       reupdate_threshold(aReupdateThreshold),
                       initial_state(aInitialState) {
    invariant_symplectic_kf_predict(*sys,
                                    b,
                                    aInitialInput,
                                    Q,
                                    Tc,
                                    Wp,
                                    aInitialTime);
    predicted_state = b.get_mean_state();
    invariant_symplectic_kf_update(*sys,
                                   b,
                                   aInitialInput,
                                   sys->get_output(predicted_state,
                                                   aInitialInput,
                                                   aInitialTime),
                                   R,
                                   Tm,
                                   Wu,
                                   aInitialTime);
  };
  
  time_difference_type get_time_step() const { return sys->get_time_step(); };

  const state_space_system& get_ss_system() const { return *sys; };

  belief_state get_next_belief(belief_state b, const time_type& t, const input_type& u, const input_type& y) {
    initial_state = b.get_mean_state();
    invariant_symplectic_kf_predict(*sys,b,u,Q,Tc,Wp,t);
    predicted_state = b.get_mean_state();
    invariant_symplectic_kf_update(*sys,b,u,y,R,Tm,Wu,t);
    return b;
  };
  
  belief_state predict_belief(belief_state b, const time_type& t, const input_type& u) {
    if( norm( diff(b.get_mean_state(), initial_state) ) > reupdate_threshold ) {
      initial_state = b.get_mean_state();
      invariant_symplectic_kf_predict(*sys,b,u,Q,Tc,Wp,t);
    } else {
      state_type x = sys->get_next_state(b.get_mean_state(), u, t);
      Wp = sys->get_invariant_prior_frame(b.get_mean_state(),x,u,t);
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
  
  belief_state prediction_to_ML_belief(belief_state b, const time_type& t, const input_type& u) {
    if( norm( diff(b.get_mean_state(), predicted_state) ) > reupdate_threshold ) {
      predicted_state = b.get_mean_state();
      invariant_symplectic_kf_update(*sys,b,u,sys->get_output(predicted_state,u,t + sys->get_time_step()),R,Tm,Wu,t);
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
  
  belief_state predict_ML_belief(belief_state b, const time_type& t, const input_type& u) {
    return prediction_to_ML_belief(predict_belief(b, t, u),t,u);
  };
  
};










};

};


#endif




















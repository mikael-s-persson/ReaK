
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

#ifndef SYMPLECTIC_KALMAN_FILTER_HPP
#define SYMPLECTIC_KALMAN_FILTER_HPP

#include "belief_state_concept.hpp"
#include "discrete_linear_sss_concept.hpp"
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
  boost::function_requires< DiscreteLinearSSSConcept< LinearSystem, DiscreteLinearizedSystemType > >();
  boost::function_requires< ContinuousBeliefStateConcept<BeliefState> >();
  
  typedef typename discrete_sss_traits<LinearSystem>::point_type StateType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  
  boost::function_requires< DecomposedCovarianceConcept<CovType> >();
  
  typedef typename decomp_covariance_mat_traits< CovType >::matrix_block_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename mat_traits<MatType>::size_type SizeType;
  
  typename discrete_linear_sss_traits<LinearSystem>::matrixA_type A;
  typename discrete_linear_sss_traits<LinearSystem>::matrixB_type B;
  typename discrete_linear_sss_traits<LinearSystem>::matrixC_type C;
  typename discrete_linear_sss_traits<LinearSystem>::matrixD_type D;
  StateType x = b.get_mean_state();
  sys.get_linear_blocks(A, B, C, D, t, x, u);
  SizeType N = A.get_col_count();
  
  b.set_mean_state(A * x + B * u);
  const MatType& X = b.get_covariance().get_covarying_block();
  const MatType& Y = b.get_covariance().get_informing_inv_block(); 
  
  Tc.set_row_count(2 * N);
  Tc.set_col_count(2 * N);
  
  set_block(Tc, mat<ValueType, mat_structure::nil>(N), N, 0);
  typename discrete_linear_sss_traits<LinearSystem>::matrixA_type A_inv; 
  pseudoinvert_QR(A,A_inv);
  set_block(Tc, transpose_move(A_inv), N, N);
  set_block(Tc, Q.get_matrix() * mat_const_sub_block< PredictionCovTransMatrix >(Tc, N, N, N, N), 0, N);
  set_block(Tc, A, 0, 0);
  
  b.set_covariance( CovType( MatType( mat_const_sub_block< PredictionCovTransMatrix >(Tc, N, 2*N, 0, 0) * mat_const_ref_vert_cat< MatType, MatType >(X,Y) ),
                             MatType( mat_const_sub_block< PredictionCovTransMatrix >(Tc, N, N, N, N) * Y ) ) );
};


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
  boost::function_requires< DiscreteLinearSSSConcept< LinearSystem, DiscreteLinearizedSystemType > >();
  boost::function_requires< ContinuousBeliefStateConcept<BeliefState> >();
  
  typedef typename discrete_sss_traits<LinearSystem>::point_type StateType;
  typedef typename discrete_sss_traits<LinearSystem>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  
  boost::function_requires< DecomposedCovarianceConcept<CovType> >();
  
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
  sys.get_linear_blocks(A, B, C, D, t, x, u);
  
  OutputType y = z - C * x - D * u;
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > YC = transpose(C);
  mat< ValueType, mat_structure::symmetric > M = YC * R.get_inverse_matrix() * C;
  linsolve_QR(Y,YC);
  mat< ValueType, mat_structure::symmetric > S = C * X * YC + R.get_matrix();
  YC = transpose(X * YC);
  linsolve_Cholesky(S,YC);
  mat< ValueType, mat_structure::rectangular, mat_alignment::row_major > K = transpose_move(YC);
   
  b.set_mean_state( x + K * y );
  b.set_covariance( CovType( X, MatType( Y + M * X ) ) );
  
  SizeType N = A.get_col_count();
  Tm.set_row_count(2 * N);
  Tm.set_col_count(2 * N);
  
  set_block(Tm, mat<ValueType, mat_structure::identity>(N), 0, 0);
  set_block(Tm, mat<ValueType, mat_structure::nil>(N), 0, N);
  set_block(Tm, M, N, 0);
  set_block(Tm, mat<ValueType, mat_structure::identity>(N), N, N);
};



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
  boost::function_requires< DiscreteLinearSSSConcept< LinearSystem, DiscreteLinearizedSystemType > >();
  boost::function_requires< ContinuousBeliefStateConcept<BeliefState> >();

  typedef typename discrete_sss_traits<LinearSystem>::point_type StateType;
  typedef typename discrete_sss_traits<LinearSystem>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  
  boost::function_requires< DecomposedCovarianceConcept<CovType> >();
  
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
  sys.get_linear_blocks(A, B, C, D, t, x, u);
  
  x = A * x + B * u;
  
  SizeType N = A.get_col_count();
  T.set_row_count(2 * N);
  T.set_col_count(2 * N);
  
  set_block(T, A, 0, 0);
  typename discrete_linear_sss_traits<LinearSystem>::matrixA_type A_inv; 
  pseudoinvert_QR(A,A_inv);
  mat_sub_block< CovTransMatrix > T_lr(T,N,N,N,N);
  T_lr = transpose_move(A_inv);
  mat_sub_block< CovTransMatrix > T_ur(T,N,N,0,N);
  T_ur = Q.get_matrix() * T_lr;
  set_block(T, mat<ValueType,mat_structure::nil>(N), N, 0);
  
  X = A * X + T_ur * Y;
  Y = T_lr * Y;
  
  OutputType y = z - C * x - D * u;
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > YC = transpose(C);
  mat< ValueType, mat_structure::symmetric > M = YC * R.get_inverse_matrix() * C;
  linsolve_QR(Y,YC);
  mat< ValueType, mat_structure::symmetric > S = C * X * YC + R.get_matrix();
  YC = transpose(X * YC);
  linsolve_Cholesky(S,YC);
  mat< ValueType, mat_structure::rectangular, mat_alignment::row_major > K = transpose_move(YC);
   
  b.set_mean_state( x + K * y );
  
  set_block(T, M * A, N, 0);
  T_lr += M * T_ur;
  
  b.set_covariance( CovType( X, MatType( Y + M * X ) ) );
};








};

};

#endif












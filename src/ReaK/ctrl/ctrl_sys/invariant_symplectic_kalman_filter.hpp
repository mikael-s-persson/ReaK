
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
	  typename MeasurementNoiseCovariance,
	  typename CovTransMatrix>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal) &&
                             is_fully_writable_matrix<CovTransMatrix>::value,
void >::type invariant_symplectic_kf_step(const InvariantSystem& sys,
					  BeliefState& b,
					  const discrete_sss_traits<InvariantSystem>::input_type& u,
					  const discrete_sss_traits<InvariantSystem>::output_type& z,
					  const SystemNoiseCovariance& Q,
					  const MeasurementNoiseCovariance& R,
					  CovTransMatrix& T,
					  typename discrete_sss_traits<InvariantSystem>::time_type t = 0) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) update
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) update
  boost::function_requires< DiscreteLinearSSSConcept< InvariantSystem, DiscreteLinearizedSystemType > >();
  boost::function_requires< InvariantSystemConcept<InvariantSystem> >();
  boost::function_requires< ContinuousBeliefStateConcept<BeliefState> >();

  typedef typename discrete_sss_traits<InvariantSystem>::point_type StateType;
  typedef typename discrete_sss_traits<InvariantSystem>::point_difference_type StateDiffType;
  typedef typename discrete_sss_traits<InvariantSystem>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  
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
  
  SizeType N = A.get_col_count();
  T.set_row_count(2 * N);
  T.set_col_count(2 * N);
  
  set_block(T, A, 0, 0);
  typename discrete_linear_sss_traits<InvariantSystem>::matrixA_type A_inv; 
  pseudoinvert_QR(A,A_inv);
  mat_sub_block< CovTransMatrix > T_lr(T,N,N,N,N);
  T_lr = transpose_move(A_inv);
  mat_sub_block< CovTransMatrix > T_ur(T,N,N,0,N);
  T_ur = Q.get_matrix() * T_lr;
  set_block(T, mat<ValueType,mat_structure::nil>(N), N, 0);
  
  X = A * X + T_ur * Y;
  Y = T_lr * Y;
  
  sys.get_output_error(e, t, x, u, z);
  sys.get_invariant_frame(W, t, x);
  
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > YC = transpose(C);
  mat< ValueType, mat_structure::symmetric > M = YC * R.get_inverse_matrix() * C;
  linsolve_QR(Y,YC);
  mat< ValueType, mat_structure::symmetric > S = C * X * YC + R.get_matrix();
  YC = transpose(X * YC);
  linsolve_Cholesky(S,YC);
  mat< ValueType, mat_structure::rectangular, mat_alignment::row_major > K = transpose_move(YC);
   
  b.set_mean_state( x + W * K * e );
  
  set_block(T, M * A, N, 0);
  T_lr += M * T_ur;
  
  b.set_covariance( CovType( X, MatType( Y + M * X ) ) );
};




};

};


#endif




















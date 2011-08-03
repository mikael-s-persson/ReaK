
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

#ifndef INVARIANT_AGGREGATE_KALMAN_FILTER_HPP
#define INVARIANT_AGGREGATE_KALMAN_FILTER_HPP

#include "belief_state_concept.hpp"
#include "discrete_linear_sss_concept.hpp"
#include "invariant_system_concept.hpp"

#include <boost/utility/enable_if.hpp>
#include <math/vect_concepts.hpp>
#include <math/mat_alg.hpp>
#include <math/mat_cholesky.hpp>

#include <boost/static_assert.hpp>
#include "covariance_concept.hpp"

#include "math/mat_star_product.hpp"


namespace ReaK {

namespace ctrl {


template <typename InvariantSystem, 
          typename BeliefState, 
	  typename SystemNoiseCovariance,
	  typename MeasurementNoiseCovariance>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal),
void >::type invariant_aggregate_kf_step(const InvariantSystem& sys,
				         BeliefState& b,
					 const discrete_sss_traits<InvariantSystem>::input_type& u,
					 const discrete_sss_traits<InvariantSystem>::output_type& z,
					 const SystemNoiseCovariance& Q,
					 const MeasurementNoiseCovariance& R,
					 typename hamiltonian_mat< typename mat_traits< typename covariance_mat_traits< typename continuous_belief_state_traits<BeliefState>::covariance_type >::matrix_type >::value_type >::type& ScSm,
					 typename hamiltonian_mat< typename mat_traits< typename covariance_mat_traits< typename continuous_belief_state_traits<BeliefState>::covariance_type >::matrix_type >::value_type >::type& Sc,
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
  typedef typename covariance_mat_traits< CovType >::matrix_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename mat_traits<MatType>::size_type SizeType;
  
  typename discrete_linear_sss_traits<InvariantSystem>::matrixA_type A;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixB_type B;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixC_type C;
  typename discrete_linear_sss_traits<InvariantSystem>::matrixD_type D;
  
  typedef typename hamiltonian_mat< ValueType >::type HamilMat;
  typedef typename hamiltonian_mat< ValueType >::upper HamilMatUp;
  typedef typename hamiltonian_mat< ValueType >::lower HamilMatLo;
  typedef typename hamiltonian_mat< ValueType >::upper_left HamilMatUL;
  typedef typename hamiltonian_mat< ValueType >::upper_right HamilMatUR;
  typedef typename hamiltonian_mat< ValueType >::lower_left HamilMatLL;
  typedef typename hamiltonian_mat< ValueType >::lower_right HamilMatLR;
  
  typedef typename invariant_system_traits<InvariantSystem>::invariant_frame_type InvFrameType;
  typedef typename invariant_system_traits<InvariantSystem>::invariant_error_type InvErrorType;
  
  StateType x = b.get_mean_state();
  MatType P = b.get_covariance().get_matrix();
  sys.get_linear_blocks(A, B, C, D, t, x, u);
  SizeType N = A.get_col_count();
  
  x = sys.get_next_state(x,u,t);
  P = ( A * P * transpose(A)) + Q.get_matrix();
  
  InvErrorType e = sys.get_output_error(x, u, z, t + sys.get_time_step());
  InvFrameType W = sys.get_invariant_prior_frame(b.get_mean_state(), x, u, t + sys.get_time_step());
  
  mat< ValueType, mat_structure::rectangular, mat_alignment::column_major > CP = C * P;
  mat< ValueType, mat_structure::symmetric > S = CP * transpose(C) + R.get_matrix();
  linsolve_Cholesky(S,CP);
  mat< ValueType, mat_structure::rectangular, mat_alignment::row_major > K = transpose_move(CP);
   
  b.set_mean_state( sys.apply_correction(x, W * K * e, u, t + sys.get_time_step()) );
  W = sys.get_invariant_posterior_frame(x_prior, b.get_mean_state(), u, t + sys.get_time_step()) * W;
  InvFrameType Wt = transpose(W);
  b.set_covariance( CovType( MatType( W * (mat< ValueType, mat_structure::identity>(K.get_row_count()) - K * C) * P * Wt ) ) );
  
  //TODO Apply the W transform somehow.
  
  HamilMat Sc_tmp(HamilMatUp(HamilMatUL(A),HamilMatUR(Q.get_matrix())),HamilMatLo(HamilMatLL(mat<ValueType,mat_structure::nil>(N)),HamilMatLR(transpose_move(A))));
  
  swap(Sc,Sc_tmp);
  HamilMat ScSm_tmp(star_product(Sc,HamilMatUp(HamilMatUL(mat<ValueType,mat_structure::identity>(N)),HamilMatUR(mat<ValueType,mat_structure::nil>(N))),HamilMatLo(HamilMatLL( transpose(C) * R.get_inverse_matrix() * C ),HamilMatLR(mat<ValueType,mat_structure::identity>(N)))));
  swap(ScSm,ScSm_tmp);
};








};

};


#endif







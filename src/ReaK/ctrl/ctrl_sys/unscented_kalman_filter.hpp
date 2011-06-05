
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

#ifndef UNSCENTED_KALMAN_FILTER_HPP
#define UNSCENTED_KALMAN_FILTER_HPP

#include "belief_state_concept.hpp"
#include "discrete_sss_concept.hpp"
#include <boost/utility/enable_if.hpp>
#include <math/vect_concepts.hpp>
#include <math/mat_alg.hpp>
#include <math/mat_cholesky.hpp>

#include <boost/static_assert.hpp>
#include "covariance_concept.hpp"


namespace ReaK {

namespace ctrl {



template <typename System, 
          typename BeliefState, 
	  typename SystemNoiseCovariance>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal),
void >::type unscented_kalman_predict(const System& sys,
				      BeliefState& b,
				      const discrete_sss_traits<System>::input_type& u,
				      const SystemNoiseCovariance& Q,
				      typename discrete_sss_traits<System>::time_type t = 0,
				      typename belief_state_traits<BeliefState>::scalar_type alpha = 1E-3,
				      typename belief_state_traits<BeliefState>::scalar_type kappa = 1,
				      typename belief_state_traits<BeliefState>::scalar_type beta = 2) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) prediction
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) prediction
  boost::function_requires< DiscreteSSSConcept< System > >();
  boost::function_requires< ContinuousBeliefStateConcept<BeliefState> >();
  
  using std::sqrt;
  
  typedef typename discrete_sss_traits<System>::point_type StateType;
  typedef typename discrete_sss_traits<System>::point_difference_type StateDiffType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  typedef typename covariance_mat_traits< CovType >::matrix_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename vect_n<ValueType>::size_type SizeType;
  
  StateType x = b.get_mean_state();
    
  MatType P = b.get_covariance().get_matrix();
  
  mat<ValueType, mat_structure::square> L_p(P.get_row_count());
  decompose_Cholesky(P,L_p);
  
  SizeType N = P.get_row_count();
  
  ValueType lambda = alpha * alpha * (N + kappa) - N;
  ValueType gamma = sqrt(ValueType(N) + lambda);
  
  vect_n< StateType > X_a(1 + 2*N);
  X_a[0] = sys.get_next_state(x, u, t);
  for(SizeType j = 0; j < N; ++j) {
    for(SizeType i = 0; i < N; ++i) {
      X_a[1+2*j][i] = x[i] + gamma * L_p(i,j);
      X_a[2+2*j][i] = x[i] - gamma * L_p(i,j);
    };
    X_a[1+2*j] = sys.get_next_state(X_a[1+2*j], u, t);
    X_a[2+2*j] = sys.get_next_state(X_a[2+2*j], u, t);
  };
  
  gamma = ValueType(1) / (ValueType(N) + lambda);
  x = (lambda * gamma) * X_a[0];
  for(SizeType j = 0; j < N; ++j)
    x += (0.5 * gamma) * (X_a[1+2*j] + X_a[2+2*j]);
  for(SizeType j = 0; j < 1 + 2*N; ++j)
    X_a[j] -= x;
  
  ValueType W_c = (lambda * gamma + ValueType(1) + alpha*alpha + beta);
  for(SizeType i = 0; i < N; ++i)
    for(SizeType j = i; j < N; ++j)
      P(i,j) = W_c * X_a[0][i] * X_a[0][j];
  
  W_c = ValueType(0.5) * gamma;
  for(SizeType k = 1; k < 1 + 2*N; ++k)
    for(SizeType i = 0; i < N; ++i)
      for(SizeType j = i; j < N; ++j)
        P(i,j) += W_c * X_a[k][i] * X_a[k][j];
  
  b.set_mean_state(x);
  b.set_covariance( CovType( P + Q.get_matrix() ) );
};


template <typename System, 
          typename BeliefState, 
	  typename MeasurementNoiseCovariance>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal),
void >::type unscented_kalman_update(const System& sys,
				     BeliefState& b,
				     const discrete_sss_traits<LinearSystem>::input_type& u,
				     const discrete_sss_traits<LinearSystem>::output_type& z,
				     const MeasurementNoiseCovariance& R,
				     typename discrete_sss_traits<LinearSystem>::time_type t = 0,
				     typename belief_state_traits<BeliefState>::scalar_type alpha = 1E-3,
				     typename belief_state_traits<BeliefState>::scalar_type kappa = 1,
				     typename belief_state_traits<BeliefState>::scalar_type beta = 2) {
  //here the requirement is that the system models a linear system which is at worse a linearized system
  // - if the system is LTI or LTV, then this will result in a basic Kalman Filter (KF) update
  // - if the system is linearized, then this will result in an Extended Kalman Filter (EKF) update
  boost::function_requires< DiscreteSSSConcept< System > >();
  boost::function_requires< ContinuousBeliefStateConcept<BeliefState> >();
  
  using std::sqrt;
  
  typedef typename discrete_sss_traits<System>::point_type StateType;
  typedef typename discrete_sss_traits<System>::point_difference_type StateDiffType;
  typedef typename discrete_sss_traits<System>::output_type OutputType;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type CovType;
  typedef typename covariance_mat_traits< CovType >::matrix_type MatType;
  typedef typename mat_traits<MatType>::value_type ValueType;
  typedef typename vect_n<ValueType>::size_type SizeType;
  
  StateType x = b.get_mean_state();
  const MatType& P = b.get_covariance().get_matrix();
  
  mat<ValueType, mat_structure::square> L_p(P.get_row_count());
  decompose_Cholesky(P,L_p);
  
  SizeType N = P.get_row_count();
  
  ValueType lambda = alpha * alpha * (N + kappa) - N;
  ValueType gamma = sqrt(ValueType(N) + lambda);
  
  vect_n< StateType > X_a(1 + 2*N);
  X_a[0] = x;
  for(SizeType j = 0; j < N; ++j) {
    for(SizeType i = 0; i < N; ++i) {
      X_a[1+2*j][i] = x[i] + gamma * L_p(i,j);
      X_a[2+2*j][i] = x[i] - gamma * L_p(i,j);
    };
  };
  
  vect_n< OutputType > Y_a(1 + 2*N);
  Y_a[0] = sys.get_output(X_a[0],u,t);
  for(SizeType j = 0; j < N; ++j) {
    Y_a[1+2*j] = sys.get_output(X_a[1+2*j],u,t);
    Y_a[2+2*j] = sys.get_output(X_a[2+2*j],u,t);
  };
  
  gamma = ValueType(1) / (ValueType(N) + lambda);
  OutputType z_p = (lambda * gamma) * Y_a[0];
  for(SizeType j = 0; j < N; ++j)
    z_p += (0.5 * gamma) * (Y_a[1+2*j] + Y_a[2+2*j]);
  for(SizeType j = 0; j < 1 + 2*N; ++j)
    Y_a[j] -= z_p;
  
  const typename covariance_mat_traits<MeasurementNoiseCovariance>::matrix_type& R_mat = R.get_matrix();
  mat<ValueType, mat_structure::symmetric> P_zz(z.size());
  ValueType W_c = (lambda * gamma + ValueType(1) + alpha*alpha + beta);
  for(SizeType i = 0; i < z.size(); ++i)
    for(SizeType j = i; j < z.size(); ++j)
      P_zz(i,j) = W_c * Y_a[0][i] * Y_a[0][j] + R_mat(i,j);
  
  W_c = ValueType(0.5) * gamma;
  for(SizeType k = 1; k < 1 + 2*N; ++k)
    for(SizeType i = 0; i < z.size(); ++i)
      for(SizeType j = i; j < z.size(); ++j)
        P_zz(i,j) += W_c * Y_a[k][i] * Y_a[k][j];
  
  mat<ValueType, mat_structure::rectangular> P_xz_t(z.size(), N);
  W_c = (lambda * gamma + ValueType(1) + alpha*alpha + beta);
  for(SizeType i = 0; i < N; ++i)
    for(SizeType j = i; j < z.size(); ++j)
      P_xz_t(j,i) = W_c * (X_a[0][i] - x[i]) * Y_a[0][j];
  
  W_c = ValueType(0.5) * gamma;
  for(SizeType k = 1; k < 1 + 2*N; ++k)
    for(SizeType i = 0; i < N; ++i)
      for(SizeType j = i; j < z.size(); ++j)
        P_xz_t(j,i) += W_c * (X_a[k][i] - x[i]) * Y_a[k][j];
  
  mat<ValueType, mat_structure::rectangular> Kt(P_xz_t);
  linsolve_Cholesky(P_zz,Kt);
  
  b.set_mean_state( x + (z - z_p) * Kt );
  b.set_covariance( CovType( P - transpose_move(Kt) * P_xz_t ) );
};



template <typename System, 
          typename BeliefState, 
	  typename SystemNoiseCovariance,
	  typename MeasurementNoiseCovariance>
typename boost::enable_if_c< is_continuous_belief_state<BeliefState>::value &&
                             (belief_state_traits<BeliefState>::representation == belief_representation::gaussian) &&
                             (belief_state_traits<BeliefState>::distribution == belief_distribution::unimodal),
void >::type unscented_kalman_filter_step(const System& sys,
					  BeliefState& b,
					  const discrete_sss_traits<System>::input_type& u,
					  const discrete_sss_traits<System>::output_type& z,
					  const SystemNoiseCovariance& Q,
					  const MeasurementNoiseCovariance& R,
					  typename discrete_sss_traits<System>::time_type t = 0,
					  typename belief_state_traits<BeliefState>::scalar_type alpha = 1E-3,
					  typename belief_state_traits<BeliefState>::scalar_type kappa = 1,
					  typename belief_state_traits<BeliefState>::scalar_type beta = 2) {
  unscented_kalman_predict(sys,b,u,Q,t,alpha,kappa,beta);
  unscented_kalman_update(sys,b,u,z,R,t,alpha,kappa,beta);
};








};

};

#endif















/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_TARGET_PREDICT_DATA_HPP
#define REAK_TARGET_PREDICT_DATA_HPP

#include "CRS_planner_data.hpp"

#include "ss_systems/satellite_invar_models.hpp"
#include "ctrl_sys/invariant_kalman_filter.hpp"
#include "ctrl_sys/gaussian_belief_space.hpp"
#include "ctrl_sys/covar_topology.hpp"
#include "ctrl_sys/belief_state_predictor.hpp"
#include "ctrl_sys/maximum_likelihood_mapping.hpp"
#include "topologies/vector_topology.hpp"
#include "topologies/temporal_space.hpp"
#include "interpolation/constant_trajectory.hpp"
#include "path_planning/transformed_trajectory.hpp"

class SoTimerSensor;


struct satellite_predict_data {
  
  typedef ReaK::ctrl::satellite3D_inv_dt_system system_type;
  typedef system_type::state_space_type state_space_type;
  
  typedef typename ReaK::ctrl::discrete_sss_traits< system_type >::input_type input_type;
  typedef ReaK::pp::constant_trajectory< ReaK::pp::vector_topology< input_type > > input_traj_type;
  
  typedef ReaK::ctrl::IKF_belief_transfer_factory< system_type > pred_factory_type;
  
  typedef ReaK::ctrl::covariance_matrix< ReaK::vect_n<double> > covar_type;
  typedef ReaK::ctrl::covar_topology< covar_type > covar_space_type;
  typedef ReaK::ctrl::gaussian_belief_space<state_space_type, covar_space_type> belief_space_type;
  
  typedef ReaK::ctrl::belief_predicted_trajectory<belief_space_type, pred_factory_type, input_traj_type> belief_pred_traj_type;
  
  typedef ReaK::pp::temporal_space<state_space_type, ReaK::pp::time_poisson_topology, ReaK::pp::time_distance_only> temp_state_space_type;
  typedef ReaK::pp::transformed_trajectory<temp_state_space_type, belief_pred_traj_type, ReaK::ctrl::maximum_likelihood_map> ML_traj_type;
  
  ReaK::shared_ptr< belief_pred_traj_type > predictor;
  ReaK::shared_ptr< ML_traj_type > trajectory;
  SoTimerSensor* animation_timer;
  volatile bool enabled;
  
};




#endif


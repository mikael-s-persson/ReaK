
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


#include "lin_alg/vect_alg.hpp"
#include "ss_systems/airship_assembled_models.hpp"

#include "ctrl_sys/invariant_kalman_filter.hpp"
#include "ctrl_sys/tsos_aug_inv_kalman_filter.hpp"

#include "ctrl_sys/gaussian_belief_state.hpp"
#include "ctrl_sys/covariance_matrix.hpp"

#include "serialization/archiver_factory.hpp"
#include "recorders/data_record_po.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;





int main(int argc, char** argv) {
  
  using namespace ReaK;
  using namespace ctrl;
  
  
  // ----------- General airship configs -----------
  
  airship_parameter_pack a_params;
  a_params.mass = 1.0;
  a_params.J = mat<double,mat_structure::symmetric>(1.0,0.0,0.0,1.0,0.0,1.0);
  a_params.added_mass_factor = 0.5;
  
  a_params.use_hot_del_q_terms = false;
  a_params.use_momentum_transfer_terms = true;
  
  a_params.gravity_acc_vect    = vect<double,3>(0.0,0.0,-9.81);
  a_params.magnetic_field_vect = vect<double,3>(0.0,0.0,0.0);
  
  a_params.IMU_position = vect<double,3>(0.0,0.0,0.0);            // relative to body-fixed frame.
  a_params.IMU_orientation = quaternion<double>(vect<double,4>(1.0,0.0,0.0,0.0)); // relative to body-fixed frame.
  
  
  // ----------- Configs of sonars -----------
  
  sonars_in_room_output_model sonars(6);
  sonars.lower_corner = vect<double,3>(0.0,0.0,0.0); // lower-corner of room, in meters.
  sonars.upper_corner = vect<double,3>(0.0,0.0,0.0); // upper-corner of room, in meters.
  
  sonars.sonar_pos[0] = vect<double,3>(0.0,0.0,0.0); // relative position of sonar
  sonars.sonar_dir[0] = vect<double,3>(0.0,0.0,-1.0); // relative direction of sonar
  
  sonars.sonar_pos[1] = vect<double,3>(0.0,0.0,0.0); // relative position of sonar
  sonars.sonar_dir[1] = vect<double,3>(0.0,0.0,-1.0); // relative direction of sonar
  
  sonars.sonar_pos[2] = vect<double,3>(0.0,0.0,0.0); // relative position of sonar
  sonars.sonar_dir[2] = vect<double,3>(0.0,0.0,-1.0); // relative direction of sonar
  
  sonars.sonar_pos[3] = vect<double,3>(0.0,0.0,0.0); // relative position of sonar
  sonars.sonar_dir[3] = vect<double,3>(0.0,0.0,-1.0); // relative direction of sonar
  
  sonars.sonar_pos[4] = vect<double,3>(0.0,0.0,0.0); // relative position of sonar
  sonars.sonar_dir[4] = vect<double,3>(0.0,0.0,-1.0); // relative direction of sonar
  
  sonars.sonar_pos[5] = vect<double,3>(0.0,0.0,0.0); // relative position of sonar
  sonars.sonar_dir[5] = vect<double,3>(0.0,0.0,-1.0); // relative direction of sonar
  
  
  
  
  
  tryphon_sgam_me sys(
    "tryphon_sgam_me", 
    a_params,
    arithmetic_tuple<
      satellite_state_model, 
      near_buoyancy_state_model, 
      eccentricity_state_model>(),
    arithmetic_tuple<airship3D_6dof_thrusters>(),
    arithmetic_tuple< 
      sonars_in_room_output_model,
      sat_gyros_output_model,
      sat_accelerometer_output_model,
      sat_magnetometer_output_model
    >(sonars, sat_gyros_output_model(), sat_accelerometer_output_model(), sat_magnetometer_output_model()),
    0.01);
  
  
  
  
  return 0;
};



















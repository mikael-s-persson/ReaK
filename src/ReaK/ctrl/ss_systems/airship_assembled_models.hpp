/**
 * \file airship_assembled_models.hpp
 * 
 * This library contains a number of building blocks for invariantized discrete-time state-space systems 
 * to describe the dynamics of an airship (UAV). These are simplified models, 
 * with forces of limited complexity applied, just free-floating dynamics with 6 dof actuation forces 
 * and some simple augmented states for drag and imbalances. These systems
 * benefit from a special integration method called the "momentum-conserving trapezoidal method" (TRAPM),
 * which is an invariant variational method that guarantees conservation of angular momentum 
 * when no actuation is applied, i.e., it is an efficient and highly stable method.
 * 
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date May 2014
 */

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

#ifndef REAK_AIRSHIP_ASSEMBLED_MODELS_HPP
#define REAK_AIRSHIP_ASSEMBLED_MODELS_HPP

#include "base/named_object.hpp"

#include "ss_systems/state_space_system_tuple.hpp"
#include "ss_systems/airship_basic_mixins.hpp"
#include "ss_systems/airship_drag_mixins.hpp"
#include "ss_systems/airship_thruster_mixins.hpp"
#include "ss_systems/airship_IMU_mixins.hpp"
#include "ss_systems/airship_sonar_mixins.hpp"


namespace ReaK {

namespace ctrl {


typedef state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<satellite_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model> 
  > airship_6dof_pq;

typedef state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<satellite_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model> 
  > airship_6dof_pqg;

typedef state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    gyros_bias_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model> 
  > airship_6dof_pqg_g;

typedef state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model, 
    gyros_bias_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model> 
  > airship_6dof_pqg_meg;

typedef state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    linear_drag_state_model,
    torsional_drag_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model> 
  > airship_6dof_pqg_med;

typedef state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    linear_drag_state_model,
    torsional_drag_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  > airship_6dof_pqgam_med;

typedef state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    gyros_bias_state_model,
    accelerometer_bias_state_model,
    magnetometer_bias_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  > airship_6dof_pqgam_megam;

typedef state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    linear_drag_state_model,
    torsional_drag_state_model,
    gyros_bias_state_model,
    accelerometer_bias_state_model,
    magnetometer_bias_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  > airship_6dof_pqgam_medgam;



typedef state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple< 
    sonars_in_room_output_model,
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  > tryphon_sgam_me;

typedef state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    linear_drag_state_model,
    torsional_drag_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple< 
    sonars_in_room_output_model,
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  > tryphon_sgam_med;

typedef state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    gyros_bias_state_model,
    accelerometer_bias_state_model,
    magnetometer_bias_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple< 
    sonars_in_room_output_model,
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  > tryphon_sgam_megam;

typedef state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    linear_drag_state_model,
    torsional_drag_state_model,
    gyros_bias_state_model,
    accelerometer_bias_state_model,
    magnetometer_bias_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple< 
    sonars_in_room_output_model,
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  > tryphon_sgam_medgam;

typedef state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    gyros_bias_state_model,
    accelerometer_bias_state_model,
    magnetometer_bias_state_model,
    room_orientation_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple< 
    sonars_in_room_output_model,
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  > tryphon_sgam_megamr;

typedef state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    linear_drag_state_model,
    torsional_drag_state_model,
    gyros_bias_state_model,
    accelerometer_bias_state_model,
    magnetometer_bias_state_model,
    room_orientation_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple< 
    sonars_in_room_output_model,
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  > tryphon_sgam_medgamr;


#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE


extern template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<satellite_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model> 
  >;

extern template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<satellite_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model> 
  >;

extern template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    gyros_bias_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model> 
  >;

extern template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model, 
    gyros_bias_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model> 
  >;

extern template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    linear_drag_state_model,
    torsional_drag_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model> 
  >;

extern template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    linear_drag_state_model,
    torsional_drag_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  >;

extern template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    gyros_bias_state_model,
    accelerometer_bias_state_model,
    magnetometer_bias_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  >;

extern template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    linear_drag_state_model,
    torsional_drag_state_model,
    gyros_bias_state_model,
    accelerometer_bias_state_model,
    magnetometer_bias_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  >;


extern template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model>, 
  arithmetic_tuple<tryphon_n_thrusters>, 
  arithmetic_tuple< 
    sonars_in_room_output_model,
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  >;

extern template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    linear_drag_state_model,
    torsional_drag_state_model>, 
  arithmetic_tuple<tryphon_n_thrusters>, 
  arithmetic_tuple< 
    sonars_in_room_output_model,
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  >;

extern template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    gyros_bias_state_model,
    accelerometer_bias_state_model,
    magnetometer_bias_state_model>, 
  arithmetic_tuple<tryphon_n_thrusters>, 
  arithmetic_tuple< 
    sonars_in_room_output_model,
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  >;

extern template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    linear_drag_state_model,
    torsional_drag_state_model,
    gyros_bias_state_model,
    accelerometer_bias_state_model,
    magnetometer_bias_state_model>, 
  arithmetic_tuple<tryphon_n_thrusters>, 
  arithmetic_tuple< 
    sonars_in_room_output_model,
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  >;

extern template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    gyros_bias_state_model,
    accelerometer_bias_state_model,
    magnetometer_bias_state_model,
    room_orientation_state_model>, 
  arithmetic_tuple<tryphon_n_thrusters>, 
  arithmetic_tuple< 
    sonars_in_room_output_model,
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  >;

extern template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<
    satellite_state_model, 
    near_buoyancy_state_model, 
    eccentricity_state_model,
    linear_drag_state_model,
    torsional_drag_state_model,
    gyros_bias_state_model,
    accelerometer_bias_state_model,
    magnetometer_bias_state_model,
    room_orientation_state_model>, 
  arithmetic_tuple<tryphon_n_thrusters>, 
  arithmetic_tuple< 
    sonars_in_room_output_model,
    sat_gyros_output_model,
    sat_accelerometer_output_model,
    sat_magnetometer_output_model> 
  >;


#endif


};

};

#endif





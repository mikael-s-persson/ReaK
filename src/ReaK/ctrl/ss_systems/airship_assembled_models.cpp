
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

#include <ReaK/core/base/defs.hpp>

#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

#include <ReaK/ctrl/ss_systems/airship_assembled_models.hpp>

namespace ReaK {

namespace ctrl {

template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<satellite_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model> 
  >;

template class state_space_system_tuple<
  airship_parameter_pack, 
  arithmetic_tuple<satellite_state_model>, 
  arithmetic_tuple<airship3D_6dof_thrusters>, 
  arithmetic_tuple<
    sat_position_output_model, 
    sat_quaternion_output_model, 
    sat_gyros_output_model> 
  >;

template class state_space_system_tuple<
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

template class state_space_system_tuple<
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

template class state_space_system_tuple<
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

template class state_space_system_tuple<
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

template class state_space_system_tuple<
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

template class state_space_system_tuple<
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


template class state_space_system_tuple<
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

template class state_space_system_tuple<
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

template class state_space_system_tuple<
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

template class state_space_system_tuple<
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

template class state_space_system_tuple<
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

template class state_space_system_tuple<
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

};

};

#else

namespace ReaK {

namespace ctrl {

void dummy_airship_assembled_models_externs_1_symbol() { };

};

};

#endif



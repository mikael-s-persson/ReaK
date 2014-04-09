
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

#include "satellite_modeling_options.hpp"

#include "serialization/archiver_factory.hpp"

#include <limits>

namespace ReaK {

namespace ctrl {



shared_ptr< satellite_model_options::temp_state_space_type > satellite_model_options::get_temporal_state_space(double aStartTime, double aEndTime) const {
#define RK_D_INF std::numeric_limits<double>::infinity()
  return shared_ptr< temp_state_space_type >(new temp_state_space_type(
    "satellite3D_temporal_space", 
    pp::make_se3_space(
      "satellite3D_state_space",
      vect<double,3>(-RK_D_INF, -RK_D_INF, -RK_D_INF),
      vect<double,3>( RK_D_INF,  RK_D_INF,  RK_D_INF),
      RK_D_INF, RK_D_INF),
    pp::time_poisson_topology("satellite3D_time_space", time_step, (aEndTime - aStartTime) * 0.5)));
};

shared_ptr< satellite_model_options::state_space_type > satellite_model_options::get_state_space() const {
  return shared_ptr< state_space_type >(new state_space_type(pp::make_se3_space(
    "satellite3D_state_space",
    vect<double,3>(-RK_D_INF, -RK_D_INF, -RK_D_INF),
    vect<double,3>( RK_D_INF,  RK_D_INF,  RK_D_INF),
    RK_D_INF, RK_D_INF)));
#undef RK_D_INF
};


shared_ptr< satellite_model_options::system_base_type > satellite_model_options::get_base_sat_system() const {
  switch(system_kind & 7) {
    case satellite_model_options::invariant:
      return shared_ptr< system_base_type >(new satellite3D_inv_dt_system(
        "satellite3D_inv", mass, inertia_tensor, time_step));
    case satellite_model_options::invar_mom2:
      return shared_ptr< system_base_type >(new satellite3D_imdt_sys(
        "satellite3D_invmid", mass, inertia_tensor, time_step, 2));
    case satellite_model_options::invar_mom:
    default:
      return shared_ptr< system_base_type >(new satellite3D_imdt_sys(
        "satellite3D_invmom", mass, inertia_tensor, time_step));
  };
};

shared_ptr< satellite_model_options::system_gyro_type > satellite_model_options::get_gyro_sat_system() const {
  switch(system_kind & 7) {
    case satellite_model_options::invariant:
      return shared_ptr< system_gyro_type >(new satellite3D_gyro_inv_dt_system(
        "satellite3D_inv_with_gyros", mass, inertia_tensor, time_step));
    case satellite_model_options::invar_mom2:
      return shared_ptr< system_gyro_type >(new satellite3D_gyro_imdt_sys(
        "satellite3D_invmid_with_gyros", mass, inertia_tensor, time_step, 2));
    case satellite_model_options::invar_mom:
    default:
      return shared_ptr< system_gyro_type >(new satellite3D_gyro_imdt_sys(
        "satellite3D_invmom_with_gyros", mass, inertia_tensor, time_step));
  };
};

shared_ptr< satellite_model_options::system_IMU_type > satellite_model_options::get_IMU_sat_system() const {
  switch(system_kind & 7) {
    case satellite_model_options::invar_mom2:
      return shared_ptr< system_IMU_type >(new satellite3D_IMU_imdt_sys(
        "satellite3D_invmid_with_IMU", mass, inertia_tensor, time_step,
        IMU_orientation, IMU_location, earth_orientation, mag_field_direction, 2));
    case satellite_model_options::invar_mom:
    default:
      return shared_ptr< system_IMU_type >(new satellite3D_IMU_imdt_sys(
        "satellite3D_invmom_with_IMU", mass, inertia_tensor, time_step,
        IMU_orientation, IMU_location, earth_orientation, mag_field_direction));
  };
};


satellite_model_options::state_belief_type satellite_model_options::get_init_state_belief(double aCovDiag) const {
  state_type x_init;
  set_frame_3D(x_init, initial_motion);
  return state_belief_type(x_init, covar_type(covar_type::matrix_type(mat<double,mat_structure::diagonal>(12, aCovDiag))));
};

satellite_model_options::input_belief_type satellite_model_options::get_zero_input_belief() const {
  return input_belief_type(input_type(vect_n<double>(6, 0.0)), covar_type(covar_type::matrix_type(input_disturbance)));
};

satellite_model_options::output_belief_type satellite_model_options::get_zero_output_belief() const {
  vect_n<double> z(get_measurement_count(), 0.0);
  z[3] = 1.0;
  if( artificial_noise.get_row_count() == measurement_noise.get_row_count() )
    return output_belief_type(output_type(z), 
                              covar_type(covar_type::matrix_type(measurement_noise + artificial_noise)));
  else
    return output_belief_type(output_type(z), 
                              covar_type(covar_type::matrix_type(measurement_noise)));
};

void satellite_model_options::imbue_names_for_generated_meas(recorder::data_stream_options& data_opt) const {
  data_opt.names.clear();
  
  data_opt
    .add_name("time").add_name("p_x").add_name("p_y").add_name("p_z")
    .add_name("q_0").add_name("q_1").add_name("q_2").add_name("q_3");
  
  switch(system_kind & 24) {
    case 8:
      data_opt
        .add_name("w_x").add_name("w_y").add_name("w_z");
      break;
    case 24:
      data_opt
        .add_name("w_x").add_name("w_y").add_name("w_z")
        .add_name("acc_x").add_name("acc_y").add_name("acc_z")
        .add_name("mag_x").add_name("mag_y").add_name("mag_z");
      break;
    default:
      break;
  };
  
  data_opt
    .add_name("f_x").add_name("f_y").add_name("f_z")
    .add_name("t_x").add_name("t_y").add_name("t_z")
    .add_name("p_x_true").add_name("p_y_true").add_name("p_z_true")
    .add_name("q_0_true").add_name("q_1_true").add_name("q_2_true").add_name("q_3_true")
    .add_name("v_x_true").add_name("v_y_true").add_name("v_z_true")
    .add_name("w_x_true").add_name("w_y_true").add_name("w_z_true");
};

void satellite_model_options::imbue_names_for_meas_stddevs(recorder::data_stream_options& data_opt) const {
  data_opt.names.clear();
  
  data_opt
    .add_name("ep_x").add_name("ep_y").add_name("ep_z")
    .add_name("ea_x").add_name("ea_y").add_name("ea_z").add_name("ep_m").add_name("ea_m");
  
  switch(system_kind & 24) {
    case 8:
      data_opt
        .add_name("ew_x").add_name("ew_y").add_name("ew_z").add_name("ew_m");
      break;
    case 24:
      data_opt
        .add_name("ew_x").add_name("ew_y").add_name("ew_z").add_name("ew_m")
        .add_name("eacc_x").add_name("eacc_y").add_name("eacc_z").add_name("eacc_m")
        .add_name("emag_x").add_name("emag_y").add_name("emag_z").add_name("emag_m");
      break;
    default:
      break;
  };
};

void satellite_model_options::imbue_names_for_state_estimates(recorder::data_stream_options& data_opt) const {
  data_opt.names.clear();
  data_opt
    .add_name("time").add_name("pos_x").add_name("pos_y").add_name("pos_z")
    .add_name("q0").add_name("q1").add_name("q2").add_name("q3")
    .add_name("vel_x").add_name("vel_y").add_name("vel_z")
    .add_name("avel_x").add_name("avel_y").add_name("avel_z")
    .add_name("ep_x").add_name("ep_y").add_name("ep_z")
    .add_name("ea_x").add_name("ea_y").add_name("ea_z")
    .add_name("ev_x").add_name("ev_y").add_name("ev_z")
    .add_name("ew_x").add_name("ew_y").add_name("ew_z")
    .add_name("P_xx").add_name("P_yy").add_name("P_zz")
    .add_name("P_aax").add_name("P_aay").add_name("P_aaz")
    .add_name("P_vvx").add_name("P_vvy").add_name("P_vvz")
    .add_name("P_wwx").add_name("P_wwy").add_name("P_wwz");
};

void satellite_model_options::imbue_names_for_state_estimates_stddevs(recorder::data_stream_options& data_opt) const {
  data_opt.names.clear();
  data_opt
    .add_name("ep_x").add_name("ep_y").add_name("ep_z")
    .add_name("ea_x").add_name("ea_y").add_name("ea_z")
    .add_name("ev_x").add_name("ev_y").add_name("ev_z")
    .add_name("ew_x").add_name("ew_y").add_name("ew_z")
    .add_name("ep_m").add_name("ea_m").add_name("ev_m").add_name("ew_m")
    .add_name("P_xx").add_name("P_yy").add_name("P_zz")
    .add_name("P_aax").add_name("P_aay").add_name("P_aaz")
    .add_name("P_vvx").add_name("P_vvy").add_name("P_vvz")
    .add_name("P_wwx").add_name("P_wwy").add_name("P_wwz");
};



void satellite_model_options::load_all_configs_impl(serialization::iarchive& in) {
  in
    & RK_SERIAL_LOAD_WITH_NAME(mass)
    & RK_SERIAL_LOAD_WITH_NAME(inertia_tensor)
    & RK_SERIAL_LOAD_WITH_NAME(IMU_orientation)
    & RK_SERIAL_LOAD_WITH_NAME(IMU_location)
    & RK_SERIAL_LOAD_WITH_NAME(earth_orientation)
    & RK_SERIAL_LOAD_WITH_NAME(mag_field_direction)
    & RK_SERIAL_LOAD_WITH_NAME(input_disturbance)
    & RK_SERIAL_LOAD_WITH_NAME(measurement_noise)
    & RK_SERIAL_LOAD_WITH_NAME(artificial_noise)
    & RK_SERIAL_LOAD_WITH_NAME(initial_motion)
    & RK_SERIAL_LOAD_WITH_NAME(system_kind);
};

void satellite_model_options::save_all_configs_impl(serialization::oarchive& out) const {
  out
    & RK_SERIAL_SAVE_WITH_NAME(mass)
    & RK_SERIAL_SAVE_WITH_NAME(inertia_tensor)
    & RK_SERIAL_SAVE_WITH_NAME(IMU_orientation)
    & RK_SERIAL_SAVE_WITH_NAME(IMU_location)
    & RK_SERIAL_SAVE_WITH_NAME(earth_orientation)
    & RK_SERIAL_SAVE_WITH_NAME(mag_field_direction)
    & RK_SERIAL_SAVE_WITH_NAME(input_disturbance)
    & RK_SERIAL_SAVE_WITH_NAME(measurement_noise)
    & RK_SERIAL_SAVE_WITH_NAME(artificial_noise)
    & RK_SERIAL_SAVE_WITH_NAME(initial_motion)
    & RK_SERIAL_SAVE_WITH_NAME(system_kind);
};


void satellite_model_options::load_all_configs(const std::string& aFileName) {
  load_all_configs_impl(*(serialization::open_iarchive(aFileName)));
};

void satellite_model_options::save_all_configs(const std::string& aFileName) const {
  save_all_configs_impl(*(serialization::open_oarchive(aFileName)));
};



void satellite_model_options::load_mass_configs(const std::string& aFileName) {
  *(serialization::open_iarchive(aFileName))
    & RK_SERIAL_LOAD_WITH_NAME(mass)
    & RK_SERIAL_LOAD_WITH_NAME(inertia_tensor);
};

void satellite_model_options::save_mass_configs(const std::string& aFileName) const {
  *(serialization::open_oarchive(aFileName))
    & RK_SERIAL_SAVE_WITH_NAME(mass)
    & RK_SERIAL_SAVE_WITH_NAME(inertia_tensor);
};


void satellite_model_options::load_IMU_configs(const std::string& aFileName) {
  *(serialization::open_iarchive(aFileName))
    & RK_SERIAL_LOAD_WITH_NAME(IMU_orientation)
    & RK_SERIAL_LOAD_WITH_NAME(IMU_location)
    & RK_SERIAL_LOAD_WITH_NAME(earth_orientation)
    & RK_SERIAL_LOAD_WITH_NAME(mag_field_direction);
};

void satellite_model_options::save_IMU_configs(const std::string& aFileName) const {
  *(serialization::open_oarchive(aFileName))
    & RK_SERIAL_SAVE_WITH_NAME(IMU_orientation)
    & RK_SERIAL_SAVE_WITH_NAME(IMU_location)
    & RK_SERIAL_SAVE_WITH_NAME(earth_orientation)
    & RK_SERIAL_SAVE_WITH_NAME(mag_field_direction);
};


void satellite_model_options::load_input_disturbance(const std::string& aFileName) {
  *(serialization::open_iarchive(aFileName))
    & RK_SERIAL_LOAD_WITH_NAME(input_disturbance);
};

void satellite_model_options::save_input_disturbance(const std::string& aFileName) const {
  *(serialization::open_oarchive(aFileName))
    & RK_SERIAL_SAVE_WITH_NAME(input_disturbance);
};

void satellite_model_options::load_measurement_noise(const std::string& aFileName) {
  *(serialization::open_iarchive(aFileName))
    & RK_SERIAL_LOAD_WITH_NAME(measurement_noise);
};

void satellite_model_options::save_measurement_noise(const std::string& aFileName) const {
  *(serialization::open_oarchive(aFileName))
    & RK_SERIAL_SAVE_WITH_NAME(measurement_noise);
};

void satellite_model_options::load_artificial_noise(const std::string& aFileName) {
  *(serialization::open_iarchive(aFileName))
    & RK_SERIAL_LOAD_WITH_NAME(artificial_noise);
};

void satellite_model_options::save_artificial_noise(const std::string& aFileName) const {
  *(serialization::open_oarchive(aFileName))
    & RK_SERIAL_SAVE_WITH_NAME(artificial_noise);
};

void satellite_model_options::load_initial_motion(const std::string& aFileName) {
  *(serialization::open_iarchive(aFileName))
    & RK_SERIAL_LOAD_WITH_NAME(initial_motion);
};

void satellite_model_options::save_initial_motion(const std::string& aFileName) const {
  *(serialization::open_oarchive(aFileName))
    & RK_SERIAL_SAVE_WITH_NAME(initial_motion);
};




void satellite_predictor_options::load_all_configs_impl(serialization::iarchive& in) {
  satellite_model_options::load_all_configs_impl(in);
  load_predict_configs_impl(in);
};

void satellite_predictor_options::save_all_configs_impl(serialization::oarchive& out) const {
  satellite_model_options::save_all_configs_impl(out);
  save_predict_configs_impl(out);
};


void satellite_predictor_options::load_predict_configs_impl(serialization::iarchive& in) {
  unsigned int pred_assume = 0;
  in
    & RK_SERIAL_LOAD_WITH_NAME(predict_time_horizon)
    & RK_SERIAL_LOAD_WITH_NAME(predict_Pnorm_threshold)
    & RK_SERIAL_LOAD_WITH_ALIAS("predict_assumption",pred_assume);
  switch(pred_assume) {
    case 2:
      predict_assumption = satellite_predictor_options::full_certainty;
      break;
    case 1:
      predict_assumption = satellite_predictor_options::most_likely_measurements;
      break;
    case 0:
    default:
      predict_assumption = satellite_predictor_options::no_measurements;
      break;
  };
};

void satellite_predictor_options::save_predict_configs_impl(serialization::oarchive& out) const {
  unsigned int pred_assume = predict_assumption;
  out
    & RK_SERIAL_SAVE_WITH_NAME(initial_motion)
    & RK_SERIAL_SAVE_WITH_NAME(predict_Pnorm_threshold)
    & RK_SERIAL_SAVE_WITH_ALIAS("predict_assumption",pred_assume);
};



void satellite_predictor_options::load_prediction_configs(const std::string& aFileName) {
  load_predict_configs_impl(*(serialization::open_iarchive(aFileName)));
};

void satellite_predictor_options::save_prediction_configs(const std::string& aFileName) const {
  save_predict_configs_impl(*(serialization::open_oarchive(aFileName)));
};



};

};









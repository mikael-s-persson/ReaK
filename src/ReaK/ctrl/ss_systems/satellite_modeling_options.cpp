
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


namespace ReaK {

namespace ctrl {


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









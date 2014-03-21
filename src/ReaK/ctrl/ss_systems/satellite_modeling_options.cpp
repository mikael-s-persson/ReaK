
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
  


};

};









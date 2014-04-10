/**
 * \file data_record_po.hpp
 *
 * This library declares utility functions for creating and dealing with program-options related 
 * to a data recording or extraction stream (see data_record.hpp). Here, "data" is meant as
 * columns of floating-point records of data, such as simulation results for example.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2014
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

#ifndef REAK_SATELLITE_MODELING_PO_HPP
#define REAK_SATELLITE_MODELING_PO_HPP

#include "satellite_modeling_options.hpp"
#include "satellite_invar_models.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>


namespace ReaK {

namespace ctrl {

/**
 * This function constructs a Boost.Program-Options descriptor for satellite modeling options.
 * This function can either construct options with or without predictor options.
 */
boost::program_options::options_description get_satellite_model_options_po_desc(bool aPredictOpt = false) {
  using boost::program_options::options_description;
  using boost::program_options::value;
  
  options_description result;
  
  options_description model_options("Satellite Modeling options");
  model_options.add_options()
    ("sat-config-file",value< std::string >(), "specify the filename for all the satellite modeling options")
    ("init-motion",    value< std::string >()->default_value("models/satellite3D_init.rkx"), "specify the filename for the satellite's initial motion (state) (default is 'models/satellite3D_init.rkx')")
    ("inertia",        value< std::string >()->default_value("models/satellite3D_inertia.rkx"), "specify the filename for the satellite's inertial data (default is 'models/satellite3D_inertia.rkx')")
    ("Q-matrix",       value< std::string >()->default_value("models/satellite3D_Q.rkx"), "specify the filename for the satellite's input disturbance covariance matrix (default is 'models/satellite3D_Q.rkx')")
    ("R-matrix",       value< std::string >()->default_value("models/satellite3D_R.rkx"), "specify the filename for the satellite's measurement noise covariance matrix (default is 'models/satellite3D_R.rkx')")
    ("R-added",        value< std::string >(), "specify the filename for the satellite's artificial measurement noise covariance matrix")
    ("IMU-config",     value< std::string >()->default_value("models/satellite3D_IMU_config.rkx"), "specify the filename for the satellite's IMU configuration data, specifying its placement on the satellite and the inertial / magnetic-field frame it is relative to (default is 'models/satellite3D_IMU_config.rkx')")
    ("generate-mdl-files", "if set, the output will be the generation of all the modeling files (with default values) into the specified file-names")
    ("system-output",  value< std::string >(), "specify the filename-stem for the output of the satellite system, when 'generate-files' is set")
    ("time-step", value< double >()->default_value(0.01), "time-step of the satellite system (default is 0.01)")
    ("gyro",      "if set, a set of gyros is added to the model (angular velocity measurements). This requires the 'R-matrix' file to contain a 9x9 matrix.")
    ("IMU",       "if set, a set of gyros is added to the model (angular velocity, magnetic field, and accelerometer measurements).\
 This requires the 'R-matrix' file to contain a 15x15 matrix. This option also automatically implies the 'midpoint' option.\
 This option will trigger the use of the 'IMU-config' file to obtain the information necessary about the IMU and the Earth's inertial frame.")
//     ("mekf",      "if set, results for the multiplicative extended Kalman filter (MEKF) will be generated.")
    ("iekf",      "if set, results for the invariant extended Kalman filter (IEKF) will be generated.")
    ("imkf",      "if set, results for the invariant momentum-tracking Kalman filter (IMKF) will be generated.")
    ("imkfv2",    "if set, results for the invariant midpoint Kalman filter (IMKFv2) will be generated.")
    ("imkf-em",   "if set, results for the invariant momentum-tracking Kalman filter (IMKF) with adaptive mass-eccentricity parameters will be generated.")
    ("imkf-emd",  "if set, results for the invariant momentum-tracking Kalman filter (IMKF) with adaptive mass-eccentricity-drag parameters will be generated.")
  ;
  result.add(model_options);
  
  if(aPredictOpt) {
    options_description pred_options("Prediction options");
    pred_options.add_options()
      ("pred-config-file", value< std::string >(), "the file containing all the predictor configurations")
      ("time-horizon",    value< double >()->default_value(100.0), "time-horizon of the satellite state predictor (default is 100.0)")
      ("Pnorm-threshold", value< double >()->default_value(1.0),   "threshold on the norm of the state prediction covariance matrix below which predictions can start (default is 1.0)")
      ("pred-assumption", value< int >()->default_value(0),        "prediction assumption to be used (0: no-measurements, 1: most-likely-measurements, 2: full-certainty")
    ;
    result.add(pred_options);
  };
  
  return result;
};


namespace detail {


void fill_satellite_model_options_from_po(satellite_model_options& result, boost::program_options::variables_map& vm) {
  
  if(vm.count("sat-config-file")) {
    result.load_all_configs(vm["sat-config-file"].as<std::string>());
  } else {
    result.load_mass_configs(vm["inertia"].as<std::string>());
    result.load_IMU_configs(vm["IMU-config"].as<std::string>());
    result.load_input_disturbance(vm["Q-matrix"].as<std::string>());
    result.load_measurement_noise(vm["R-matrix"].as<std::string>());
    if(vm.count("R-added"))
      result.load_artificial_noise(vm["R-added"].as<std::string>());
    result.load_initial_motion(vm["init-motion"].as<std::string>());
  };
  
  result.time_step = vm["time-step"].as<double>();
  
  result.system_kind = 0;
  if(vm.count("gyro"))
    result.system_kind |= satellite_model_options::gyro_measures;
  if(vm.count("IMU"))
    result.system_kind |= satellite_model_options::IMU_measures;
  
  if(vm.count("iekf"))
    result.system_kind |= satellite_model_options::invariant;
  else if(vm.count("imkfv2"))
    result.system_kind |= satellite_model_options::invar_mom2;
  else if(vm.count("imkf-em"))
    result.system_kind |= satellite_model_options::invar_mom_em;
  else if(vm.count("imkf-emd"))
    result.system_kind |= satellite_model_options::invar_mom_emd;
  else 
    result.system_kind |= satellite_model_options::invar_mom;
  
};


void save_satellite_model_options_to_files(satellite_model_options& result, boost::program_options::variables_map& vm) {
  
  if(result.input_disturbance.get_row_count() != 6)
    result.input_disturbance = mat<double,mat_structure::diagonal>(6, true);
  if(result.measurement_noise.get_row_count() != result.get_meas_error_count())
    result.measurement_noise = mat<double,mat_structure::diagonal>(result.get_meas_error_count(), true);
  if(result.artificial_noise.get_row_count() != result.get_meas_error_count())
    result.artificial_noise = mat<double,mat_structure::diagonal>(result.get_meas_error_count(), true);
  
  if(vm.count("sat-config-file")) {
    boost::filesystem::create_directories(boost::filesystem::path(vm["sat-config-file"].as<std::string>()).parent_path());
    result.save_all_configs(vm["sat-config-file"].as<std::string>());
  } else {
    boost::filesystem::create_directories(boost::filesystem::path(vm["inertia"].as<std::string>()).parent_path());
    result.save_mass_configs(vm["inertia"].as<std::string>());
    
    boost::filesystem::create_directories(boost::filesystem::path(vm["IMU-config"].as<std::string>()).parent_path());
    result.save_IMU_configs(vm["IMU-config"].as<std::string>());
    
    boost::filesystem::create_directories(boost::filesystem::path(vm["Q-matrix"].as<std::string>()).parent_path());
    result.save_input_disturbance(vm["Q-matrix"].as<std::string>());
    
    boost::filesystem::create_directories(boost::filesystem::path(vm["R-matrix"].as<std::string>()).parent_path());
    result.save_measurement_noise(vm["R-matrix"].as<std::string>());
    
    if(vm.count("R-added")) {
      boost::filesystem::create_directories(boost::filesystem::path(vm["R-added"].as<std::string>()).parent_path());
      result.save_artificial_noise(vm["R-added"].as<std::string>());
    };
    
    boost::filesystem::create_directories(boost::filesystem::path(vm["init-motion"].as<std::string>()).parent_path());
    result.save_initial_motion(vm["init-motion"].as<std::string>());
  };
  if(vm.count("system-output")) {
    
  };
  
};


};


/**
 * This function constructs a satellite-modeling options object from a Boost.Program-Options variable-map.
 */
satellite_model_options get_satellite_model_options_from_po(boost::program_options::variables_map& vm) {
  satellite_model_options result;
  if(vm.count("generate-mdl-files")) {
    detail::save_satellite_model_options_to_files(result, vm);
    return result;
  };
  detail::fill_satellite_model_options_from_po(result, vm);
  return result;
};

/**
 * This function constructs a satellite-predictor options object from a Boost.Program-Options variable-map.
 */
satellite_predictor_options get_satellite_predictor_options_from_po(boost::program_options::variables_map& vm) {
  satellite_predictor_options result;
  if(vm.count("generate-mdl-files")) {
    detail::save_satellite_model_options_to_files(result, vm);
    if(vm.count("pred-config-file")) {
      boost::filesystem::create_directories(boost::filesystem::path(vm["pred-config-file"].as<std::string>()).parent_path());
      result.save_prediction_configs(vm["pred-config-file"].as<std::string>());
    };
    return result;
  };
  
  detail::fill_satellite_model_options_from_po(result, vm);
  
  if(vm.count("pred-config-file")) {
    result.load_prediction_configs(vm["pred-config-file"].as<std::string>());
  } else {
    result.predict_time_horizon = vm["time-horizon"].as<double>();
    result.predict_Pnorm_threshold = vm["Pnorm-threshold"].as<double>();
    switch(vm["pred-assumption"].as<int>()) {
      case 1:
        result.predict_assumption = satellite_predictor_options::most_likely_measurements;
        break;
      case 2:
        result.predict_assumption = satellite_predictor_options::full_certainty;
        break;
      case 0:
      default:
        result.predict_assumption = satellite_predictor_options::no_measurements;
        break;
    };
  };
  
  return result;
};


};

};

#endif // RK_DATA_RECORD_PO_HPP













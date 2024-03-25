
/*
 *    Copyright 2023 Sven Mikael Persson
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
#include "ReaK/control/systems/satellite_modeling_po.h"
#include "ReaK/control/systems/satellite_modeling_options.h"

#include <filesystem>

#include "absl/flags/flag.h"

ABSL_FLAG(std::string, sat_config_file, "",
          "Specify the filename for all the satellite modeling options.");
ABSL_FLAG(std::string, sat_init_motion, "models/satellite3D_init.rkx",
          "Specify the filename for the satellite's initial motion (state) "
          "(default is 'models/satellite3D_init.rkx').");
ABSL_FLAG(std::string, sat_inertia, "models/satellite3D_inertia.rkx",
          "Specify the filename for the satellite's inertial data (default is "
          "'models/satellite3D_inertia.rkx').");
ABSL_FLAG(std::string, sat_Q_matrix, "models/satellite3D_Q.rkx",
          "Specify the filename for the satellite's input disturbance "
          "covariance matrix (default is 'models/satellite3D_Q.rkx').");
ABSL_FLAG(std::string, sat_R_matrix, "models/satellite3D_R.rkx",
          "Specify the filename for the satellite's measurement noise "
          "covariance matrix (default is 'models/satellite3D_R.rkx').");
ABSL_FLAG(std::string, sat_R_added, "",
          "Specify the filename for the satellite's artificial measurement "
          "noise covariance matrix.");
ABSL_FLAG(std::string, sat_IMU_config, "models/satellite3D_IMU_config.rkx",
          "Specify the filename for the satellite's IMU configuration data, "
          "specifying its placement on the satellite and the inertial / "
          "magnetic-field frame it is relative to (default is "
          "'models/satellite3D_IMU_config.rkx').");
ABSL_FLAG(bool, sat_generate_mdl_files, false,
          "If set, the output will be the generation of all the modeling files "
          "(with default values) into the specified file-names.");
ABSL_FLAG(std::string, sat_system_output, "",
          "Specify the filename-stem for the output of the satellite system, "
          "when 'generate-files' is set.");
ABSL_FLAG(double, sat_time_step, 0.01,
          "Time-step of the satellite system (default is 0.01).");
ABSL_FLAG(bool, sat_gyro, false,
          "If set, a set of gyros is added to the model (angular velocity "
          "measurements). This requires the 'R-matrix' file to contain a 9x9 "
          "matrix.");
ABSL_FLAG(bool, sat_IMU, false,
          "If set, a set of gyros is added to the model (angular velocity, "
          "magnetic field, and accelerometer measurements). This requires the "
          "'R-matrix' file to contain a 15x15 matrix. This option also "
          "automatically implies the 'midpoint' option. This option will "
          "trigger the use of the 'IMU-config' file to obtain the information "
          "necessary about the IMU and the Earth's inertial frame.");
ABSL_FLAG(bool, sat_iekf, false,
          "If set, results for the invariant extended Kalman filter (IEKF) "
          "will be generated.");
ABSL_FLAG(bool, sat_imkf, false,
          "If set, results for the invariant momentum-tracking Kalman filter "
          "(IMKF) will be generated.");
ABSL_FLAG(bool, sat_imkfv2, false,
          "If set, results for the invariant midpoint Kalman filter (IMKFv2) "
          "will be generated.");
ABSL_FLAG(
    bool, sat_imkf_em, false,
    "If set, results for the invariant momentum-tracking Kalman filter (IMKF) "
    "with adaptive mass-eccentricity parameters will be generated.");
ABSL_FLAG(
    bool, sat_imkf_emd, false,
    "If set, results for the invariant momentum-tracking Kalman filter (IMKF) "
    "with adaptive mass-eccentricity-drag parameters will be generated.");
ABSL_FLAG(bool, sat_imkf_emdJ, false,
          "If set, results for the invariant momentum-tracking Kalman filter "
          "(IMKF) with adaptive mass-eccentricity-drag-inertia parameters will "
          "be generated.");
ABSL_FLAG(
    bool, sat_tsosakf, false,
    "If set, the two-stage online-steady augmented Kalman filter (TSOSAKF) "
    "will be used for the augmented systems (parameter identification).");
ABSL_FLAG(std::string, sat_Pa_matrix, "models/satellite3D_Pa.rkx",
          "specify the filename for the satellite's steady augmented parameter "
          "covariance matrix (default is 'models/satellite3D_Pa.rkx').");

ABSL_FLAG(std::string, sat_pred_config_file, "",
          "The file containing all the predictor configurations.");
ABSL_FLAG(double, sat_time_horizon, 100.0,
          "Time-horizon of the satellite state predictor (default is 100.0).");
ABSL_FLAG(double, sat_Pnorm_threshold, 1.0,
          "Threshold on the norm of the state prediction covariance matrix "
          "below which predictions can start (default is 1.0).");
ABSL_FLAG(int, sat_pred_assumption, 0,
          "Prediction assumption to be used (0: no-measurements, 1: "
          "most-likely-measurements, 2: full-certainty.");

namespace ReaK::ctrl {

namespace detail {

void fill_satellite_model_options_from_flags(satellite_model_options& result) {

  if (!absl::GetFlag(FLAGS_sat_config_file).empty()) {
    result.load_all_configs(absl::GetFlag(FLAGS_sat_config_file));
  } else {
    result.load_mass_configs(absl::GetFlag(FLAGS_sat_inertia));
    result.load_IMU_configs(absl::GetFlag(FLAGS_sat_IMU_config));
    result.load_input_disturbance(absl::GetFlag(FLAGS_sat_Q_matrix));
    result.load_measurement_noise(absl::GetFlag(FLAGS_sat_R_matrix));
    if (!absl::GetFlag(FLAGS_sat_R_added).empty()) {
      result.load_artificial_noise(absl::GetFlag(FLAGS_sat_R_added));
    }
    result.load_initial_motion(absl::GetFlag(FLAGS_sat_init_motion));
  }

  result.time_step = absl::GetFlag(FLAGS_sat_time_step);

  result.sys_output_stem_name = absl::GetFlag(FLAGS_sat_system_output);

  result.system_kind = {};
  if (absl::GetFlag(FLAGS_sat_gyro)) {
    result.system_kind.measures[static_cast<std::size_t>(
        satellite_model_options::available_measurements::gyro_measures)] = true;
  }
  if (absl::GetFlag(FLAGS_sat_IMU)) {
    result.system_kind.measures[static_cast<std::size_t>(
        satellite_model_options::available_measurements::gyro_measures)] = true;
    result.system_kind.measures[static_cast<std::size_t>(
        satellite_model_options::available_measurements::IMU_measures)] = true;
  }

  if (absl::GetFlag(FLAGS_sat_iekf)) {
    result.system_kind.model = satellite_model_options::model_kind::invariant;
  } else if (absl::GetFlag(FLAGS_sat_imkfv2)) {
    result.system_kind.model = satellite_model_options::model_kind::invar_mom2;
  } else if (absl::GetFlag(FLAGS_sat_imkf_em)) {
    result.system_kind.model =
        satellite_model_options::model_kind::invar_mom_em;
  } else if (absl::GetFlag(FLAGS_sat_imkf_emd)) {
    result.system_kind.model =
        satellite_model_options::model_kind::invar_mom_emd;
  } else if (absl::GetFlag(FLAGS_sat_imkf_emdJ)) {
    result.system_kind.model =
        satellite_model_options::model_kind::invar_mom_emdJ;
  } else {
    result.system_kind.model = satellite_model_options::model_kind::invar_mom;
  }

  if (absl::GetFlag(FLAGS_sat_tsosakf)) {
    result.load_steady_param_covariance(absl::GetFlag(FLAGS_sat_Pa_matrix));
    result.system_kind.filters[static_cast<std::size_t>(
        satellite_model_options::filtering_options::TSOSAKF)] = true;
  }
}

void save_satellite_model_options_to_files(satellite_model_options& result) {

  result.time_step = absl::GetFlag(FLAGS_sat_time_step);

  result.system_kind = {};
  if (absl::GetFlag(FLAGS_sat_gyro)) {
    result.system_kind.measures[static_cast<std::size_t>(
        satellite_model_options::available_measurements::gyro_measures)] = true;
  }
  if (absl::GetFlag(FLAGS_sat_IMU)) {
    result.system_kind.measures[static_cast<std::size_t>(
        satellite_model_options::available_measurements::gyro_measures)] = true;
    result.system_kind.measures[static_cast<std::size_t>(
        satellite_model_options::available_measurements::IMU_measures)] = true;
  }

  if (absl::GetFlag(FLAGS_sat_iekf)) {
    result.system_kind.model = satellite_model_options::model_kind::invariant;
  } else if (absl::GetFlag(FLAGS_sat_imkfv2)) {
    result.system_kind.model = satellite_model_options::model_kind::invar_mom2;
  } else if (absl::GetFlag(FLAGS_sat_imkf_em)) {
    result.system_kind.model =
        satellite_model_options::model_kind::invar_mom_em;
  } else if (absl::GetFlag(FLAGS_sat_imkf_emd)) {
    result.system_kind.model =
        satellite_model_options::model_kind::invar_mom_emd;
  } else if (absl::GetFlag(FLAGS_sat_imkf_emdJ)) {
    result.system_kind.model =
        satellite_model_options::model_kind::invar_mom_emdJ;
  } else {
    result.system_kind.model = satellite_model_options::model_kind::invar_mom;
  }

  if (absl::GetFlag(FLAGS_sat_tsosakf)) {
    result.system_kind.filters[static_cast<std::size_t>(
        satellite_model_options::filtering_options::TSOSAKF)] = true;
  }

  if (result.input_disturbance.get_row_count() != 6) {
    result.input_disturbance = mat<double, mat_structure::diagonal>(6, true);
  }
  if (result.measurement_noise.get_row_count() !=
      result.get_meas_error_count()) {
    result.measurement_noise = mat<double, mat_structure::diagonal>(
        result.get_meas_error_count(), true);
  }
  if (result.artificial_noise.get_row_count() !=
      result.get_meas_error_count()) {
    result.artificial_noise = mat<double, mat_structure::diagonal>(
        result.get_meas_error_count(), true);
  }
  if (result.steady_param_covariance.get_row_count() !=
      result.get_total_inv_state_count() -
          result.get_actual_inv_state_count()) {
    result.steady_param_covariance = mat<double, mat_structure::diagonal>(
        result.get_total_inv_state_count() -
            result.get_actual_inv_state_count(),
        true);
  }

  if (!absl::GetFlag(FLAGS_sat_config_file).empty()) {
    std::filesystem::create_directories(
        std::filesystem::path(absl::GetFlag(FLAGS_sat_config_file))
            .parent_path());
    result.save_all_configs(absl::GetFlag(FLAGS_sat_config_file));
  } else {
    std::filesystem::create_directories(
        std::filesystem::path(absl::GetFlag(FLAGS_sat_inertia)).parent_path());
    result.save_mass_configs(absl::GetFlag(FLAGS_sat_inertia));

    std::filesystem::create_directories(
        std::filesystem::path(absl::GetFlag(FLAGS_sat_IMU_config))
            .parent_path());
    result.save_IMU_configs(absl::GetFlag(FLAGS_sat_IMU_config));

    std::filesystem::create_directories(
        std::filesystem::path(absl::GetFlag(FLAGS_sat_Q_matrix)).parent_path());
    result.save_input_disturbance(absl::GetFlag(FLAGS_sat_Q_matrix));

    std::filesystem::create_directories(
        std::filesystem::path(absl::GetFlag(FLAGS_sat_R_matrix)).parent_path());
    result.save_measurement_noise(absl::GetFlag(FLAGS_sat_R_matrix));

    if (!absl::GetFlag(FLAGS_sat_R_added).empty()) {
      std::filesystem::create_directories(
          std::filesystem::path(absl::GetFlag(FLAGS_sat_R_added))
              .parent_path());
      result.save_artificial_noise(absl::GetFlag(FLAGS_sat_R_added));
    }

    if (absl::GetFlag(FLAGS_sat_tsosakf)) {
      std::filesystem::create_directories(
          std::filesystem::path(absl::GetFlag(FLAGS_sat_Pa_matrix))
              .parent_path());
      result.save_steady_param_covariance(absl::GetFlag(FLAGS_sat_Pa_matrix));
    }

    std::filesystem::create_directories(
        std::filesystem::path(absl::GetFlag(FLAGS_sat_init_motion))
            .parent_path());
    result.save_initial_motion(absl::GetFlag(FLAGS_sat_init_motion));
  }
}

}  // namespace detail

/**
 * This function constructs a satellite-modeling options object from a Boost.Program-Options variable-map.
 */
satellite_model_options get_satellite_model_options_from_flags() {
  satellite_model_options result;
  if (absl::GetFlag(FLAGS_sat_generate_mdl_files)) {
    detail::save_satellite_model_options_to_files(result);
    return result;
  }
  detail::fill_satellite_model_options_from_flags(result);
  return result;
}

/**
 * This function constructs a satellite-predictor options object from a Boost.Program-Options variable-map.
 */
satellite_predictor_options get_satellite_predictor_options_from_flags() {
  satellite_predictor_options result;
  if (absl::GetFlag(FLAGS_sat_generate_mdl_files)) {
    detail::save_satellite_model_options_to_files(result);
    if (!absl::GetFlag(FLAGS_sat_pred_config_file).empty()) {
      std::filesystem::create_directories(
          std::filesystem::path(absl::GetFlag(FLAGS_sat_pred_config_file))
              .parent_path());
      result.save_prediction_configs(absl::GetFlag(FLAGS_sat_pred_config_file));
    }
    return result;
  }

  detail::fill_satellite_model_options_from_flags(result);

  if (!absl::GetFlag(FLAGS_sat_pred_config_file).empty()) {
    result.load_prediction_configs(absl::GetFlag(FLAGS_sat_pred_config_file));
  } else {
    result.predict_time_horizon = absl::GetFlag(FLAGS_sat_time_horizon);
    result.predict_Pnorm_threshold = absl::GetFlag(FLAGS_sat_Pnorm_threshold);
    switch (absl::GetFlag(FLAGS_sat_pred_assumption)) {
      case 1:
        result.predict_assumption =
            satellite_predictor_options::most_likely_measurements;
        break;
      case 2:
        result.predict_assumption = satellite_predictor_options::full_certainty;
        break;
      case 0:
      default:
        result.predict_assumption =
            satellite_predictor_options::no_measurements;
        break;
    }
  }

  return result;
}

}  // namespace ReaK::ctrl

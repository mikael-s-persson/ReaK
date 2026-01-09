/**
 * \file satellite_modeling_options.h
 *
 * This library declares utility classes for building option-sets related
 * to satellite system models (see satellite_invar_models.hpp), and using those options
 * in a factory-function to create an abstract systems and filters.
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

#ifndef REAK_CONTROL_SYSTEMS_SATELLITE_MODELING_OPTIONS_H_
#define REAK_CONTROL_SYSTEMS_SATELLITE_MODELING_OPTIONS_H_

#include "ReaK/control/systems/near_buoyant_airship_models.h"
#include "ReaK/control/systems/satellite_invar_models.h"

#include "ReaK/control/estimators/belief_state_predictor.h"
#include "ReaK/control/estimators/covar_topology.h"
#include "ReaK/control/estimators/covariance_matrix.h"
#include "ReaK/control/estimators/gaussian_belief_space.h"
#include "ReaK/control/estimators/invariant_kalman_filter.h"
#include "ReaK/control/estimators/maximum_likelihood_mapping.h"
#include "ReaK/control/systems/discrete_sss_concept.h"
#include "ReaK/topologies/interpolation/constant_trajectory.h"
#include "ReaK/topologies/interpolation/transformed_trajectory.h"
#include "ReaK/topologies/spaces/temporal_space.h"
#include "ReaK/topologies/spaces/time_poisson_topology.h"
#include "ReaK/topologies/spaces/vector_topology.h"

#include "ReaK/math/lin_alg/mat_alg_diagonal.h"
#include "ReaK/math/lin_alg/mat_alg_symmetric.h"

#include "ReaK/core/recorders/data_record_options.h"
#include "ReaK/core/serialization/archiver.h"

#include <bitset>
#include <string>

namespace ReaK::ctrl {

/**
 * This class stores a number of options related to the modeling of satellite
 * systems (see satellite_invar_models.hpp).
 * \note This class is mainly intended to be used with satellite_modeling_po (program-options).
 */
struct satellite_model_options {

  using system_base_type = satellite3D_inv_dt_system;
  using system_gyro_type = satellite3D_gyro_inv_dt_system;
  using system_IMU_type = satellite3D_IMU_imdt_sys;

  using system_em_type = airship3D_imdt_em_sys;
  using system_emd_type = airship3D_imdt_emd_sys;
  using system_gyro_emd_type = airship3D_gyro_imdt_emd_sys;
  using system_emdJ_type = airship3D_imdt_emdJ_sys;
  using system_gyro_emdJ_type = airship3D_gyro_imdt_emdJ_sys;

  using state_type = discrete_sss_traits<system_base_type>::point_type;
  using input_type = discrete_sss_traits<system_base_type>::input_type;
  using output_type = discrete_sss_traits<system_base_type>::output_type;

  using state_em_type = discrete_sss_traits<system_em_type>::point_type;
  using state_emd_type = discrete_sss_traits<system_emd_type>::point_type;
  using state_emdJ_type = discrete_sss_traits<system_emdJ_type>::point_type;

  using covar_type = covariance_matrix<vect_n<double>>;
  using covar_space_type = covar_topology<covar_type>;

  using state_belief_type = gaussian_belief_state<state_type, covar_type>;
  using input_belief_type = gaussian_belief_state<input_type, covar_type>;
  using output_belief_type = gaussian_belief_state<output_type, covar_type>;

  using state_em_belief_type = gaussian_belief_state<state_em_type, covar_type>;
  using state_emd_belief_type =
      gaussian_belief_state<state_emd_type, covar_type>;
  using state_emdJ_belief_type =
      gaussian_belief_state<state_emdJ_type, covar_type>;

  using state_space_type = system_base_type::state_space_type;
  using temp_state_space_type =
      pp::temporal_space<state_space_type, pp::time_poisson_topology,
                         pp::time_distance_only>;
  using belief_space_type =
      gaussian_belief_space<state_space_type, covar_space_type>;

  using state_space_em_type = system_em_type::state_space_type;
  using temp_state_space_em_type =
      pp::temporal_space<state_space_em_type, pp::time_poisson_topology,
                         pp::time_distance_only>;
  using belief_space_em_type =
      gaussian_belief_space<state_space_em_type, covar_space_type>;

  using state_space_emd_type = system_emd_type::state_space_type;
  using temp_state_space_emd_type =
      pp::temporal_space<state_space_emd_type, pp::time_poisson_topology,
                         pp::time_distance_only>;
  using belief_space_emd_type =
      gaussian_belief_space<state_space_emd_type, covar_space_type>;

  using state_space_emdJ_type = system_emdJ_type::state_space_type;
  using temp_state_space_emdJ_type =
      pp::temporal_space<state_space_emdJ_type, pp::time_poisson_topology,
                         pp::time_distance_only>;
  using belief_space_emdJ_type =
      gaussian_belief_space<state_space_emdJ_type, covar_space_type>;

  /// Stores the time-step of the system (if discrete-time).
  double time_step{0.01};

  /// Stores the mass of the satellite (in kg).
  double mass{1.0};

  /// Stores the inertia-tensor of the satellite (in body-fixed frame, in kg-m2).
  mat<double, mat_structure::symmetric> inertia_tensor;

  /// Stores the input-disturbance covariance matrix (Q), disturbance on the net wrench on the system.
  mat<double, mat_structure::diagonal> input_disturbance;

  /// Stores the measurement-noise covariance matrix (R).
  mat<double, mat_structure::diagonal> measurement_noise;

  /// Stores the artificial measurement-noise covariance matrix (R). This is sometimes used to add noise to measurements
  /// (for testing the filtering methods).
  mat<double, mat_structure::diagonal> artificial_noise;

  /// Stores the steady-state parameter identification covariance matrix (Pa). This is to be used in combination with
  /// the TSOSAKF methods.
  mat<double, mat_structure::diagonal> steady_param_covariance;

  /// Stores the orientation of the IMU with respect to the body-fixed frame of the satellite.
  unit_quat<double> IMU_orientation;
  /// Stores the position of the IMU with respect to the body-fixed frame of the satellite.
  vect<double, 3> IMU_location;
  /// Stores the orientation of the 'Earth' frame with respect to the global inertial frame used.
  unit_quat<double> earth_orientation;
  /// Stores the expected direction of the Earth's magnetic field lines at the location of the satellite, and relative
  /// to the Earth-frame.
  vect<double, 3> mag_field_direction;

  /// Stores the initial motion of the satellite.
  frame_3D<double> initial_motion;

  /// Specify the filename-stem for the output of the satellite system, when 'generate-files' is set.
  std::string sys_output_stem_name;

  enum class available_measurements : std::uint8_t {
    pose_measures,
    gyro_measures,
    IMU_measures
  };

  enum class model_kind : std::uint8_t {
    invariant,
    invar_mom,
    invar_mom2,
    invar_mom_em,
    invar_mom_emd,
    invar_mom_emdJ
  };

  enum class filtering_options : std::uint8_t { plain_KF, TSOSAKF };

  struct system_kind_record {
    model_kind model = model_kind::invariant;
    std::bitset<8> measures{};
    std::bitset<8> filters{};
  };

  /// Stores the kind of system used to model the satellite (OR-combination of 'available_measurements' 'model_kind').
  system_kind_record system_kind;

  bool has_gyro_measures() const {
    return system_kind.measures[static_cast<std::size_t>(
        available_measurements::gyro_measures)];
  }

  bool has_IMU_measures() const {
    return system_kind.measures[static_cast<std::size_t>(
        available_measurements::IMU_measures)];
  }

  bool has_TSOSAKF_filter() const {
    return system_kind
        .filters[static_cast<std::size_t>(filtering_options::TSOSAKF)];
  }

  std::size_t get_measurement_count() const {
    std::size_t result = 7;
    if (has_gyro_measures()) {
      result = 10;
    }
    if (has_IMU_measures()) {
      result = 16;
    }
    return result;
  }

  std::size_t get_meas_error_count() const {
    std::size_t result = 6;
    if (has_gyro_measures()) {
      result = 9;
    }
    if (has_IMU_measures()) {
      result = 15;
    }
    return result;
  }

  std::size_t get_total_inv_state_count() const {
    std::size_t result = 0;
    switch (system_kind.model) {
      case model_kind::invar_mom_em:
        result = 16;
        break;
      case model_kind::invar_mom_emd:
        result = 18;
        break;
      case model_kind::invar_mom_emdJ:
        result = 24;
        break;
      case model_kind::invariant:
      case model_kind::invar_mom2:
      case model_kind::invar_mom:
      default:
        result = 12;
        break;
    }
    return result;
  }

  std::size_t get_actual_inv_state_count() const { return 12; }

  std::string get_kf_accronym() const {
    std::string result;

    std::string kf_method = "imkf";
    if (has_TSOSAKF_filter()) {
      kf_method = "tsosakf";
    }

    switch (system_kind.model) {
      case model_kind::invariant:
        result = "iekf";
        break;
      case model_kind::invar_mom2:
        result = "imkfv2";
        break;
      case model_kind::invar_mom_em:
        result = kf_method + "_em";
        break;
      case model_kind::invar_mom_emd:
        result = kf_method + "_emd";
        break;
      case model_kind::invar_mom_emdJ:
        result = kf_method + "_emdJ";
        break;
      case model_kind::invar_mom:
      default:
        result = "imkf";
        break;
    }

    if (has_IMU_measures()) {
      result += "_IMU";
    } else if (has_gyro_measures()) {
      result += "_gyro";
    }
    return result;
  }

  std::string get_sys_abbreviation() const {
    std::string result;
    switch (system_kind.model) {
      case model_kind::invariant:
        result = "inv";
        break;
      case model_kind::invar_mom2:
        result = "invmid";
        break;
      case model_kind::invar_mom_em:
        result = "invmid_em";
        break;
      case model_kind::invar_mom_emd:
        result = "invmid_emd";
        break;
      case model_kind::invar_mom_emdJ:
        result = "invmid_emdJ";
        break;
      case model_kind::invar_mom:
      default:
        result = "invmom";
        break;
    }

    if (has_IMU_measures()) {
      result += "_IMU";
    } else if (has_gyro_measures()) {
      result += "_gyro";
    }
    return result;
  }

 protected:
  virtual void load_all_configs_impl(serialization::iarchive& in);
  virtual void save_all_configs_impl(serialization::oarchive& out) const;

 public:
  /**
   * Constructs a base satellite system.
   * \return A newly created base satellite system (as shared-pointer).
   */
  std::shared_ptr<system_base_type> get_base_sat_system() const;

  /**
   * Constructs a satellite system (with gyro).
   * \return A newly created gyro satellite system (as shared-pointer).
   */
  std::shared_ptr<system_gyro_type> get_gyro_sat_system() const;

  /**
   * Constructs a satellite system (with IMU).
   * \return A newly created IMU satellite system (as shared-pointer).
   */
  std::shared_ptr<system_IMU_type> get_IMU_sat_system() const;

  /**
   * Constructs a eccentricity-imbalance airship system.
   * \return A newly created eccentricity-imbalance airship system (as shared-pointer).
   */
  std::shared_ptr<system_em_type> get_em_airship_system() const;

  /**
   * Constructs a eccentricity-imbalance-drag airship system.
   * \return A newly created eccentricity-imbalance-drag airship system (as shared-pointer).
   */
  std::shared_ptr<system_emd_type> get_emd_airship_system() const;

  /**
   * Constructs a eccentricity-imbalance-drag airship system with gyro measurements.
   * \return A newly created eccentricity-imbalance-drag airship system with gyro measurements (as shared-pointer).
   */
  std::shared_ptr<system_gyro_emd_type> get_gyro_emd_airship_system() const;

  /**
   * Constructs a eccentricity-imbalance-drag-inertia airship system.
   * \return A newly created eccentricity-imbalance-drag-inertia airship system (as shared-pointer).
   */
  std::shared_ptr<system_emdJ_type> get_emdJ_airship_system() const;

  /**
   * Constructs a eccentricity-imbalance-drag-inertia airship system with gyro measurements.
   * \return A newly created eccentricity-imbalance-drag-inertia airship system with gyro measurements (as
   * shared-pointer).
   */
  std::shared_ptr<system_gyro_emdJ_type> get_gyro_emdJ_airship_system() const;

  /**
   * Create a belief point for the state, with mean set to initial-motion.
   * \param aCovDiag The initial diagonal values for the covariant matrix (should be high).
   * \return A belief point for the state, with mean set to initial-motion and high covariance.
   */
  state_belief_type get_init_state_belief(double aCovDiag = 10.0) const;

  /**
   * Create a belief point for the state (augmented with imbalance and eccentricity parameters),
   * with mean set to initial-motion.
   * \param aCovDiag The initial diagonal values for the covariant matrix (should be high).
   * \return A belief point for the state, with mean set to initial-motion and high covariance.
   */
  state_em_belief_type get_init_state_em_belief(double aCovDiag = 10.0) const;

  /**
   * Create a belief point for the state (augmented with imbalance, eccentricity and drag parameters),
   * with mean set to initial-motion.
   * \param aCovDiag The initial diagonal values for the covariant matrix (should be high).
   * \return A belief point for the state, with mean set to initial-motion and high covariance.
   */
  state_emd_belief_type get_init_state_emd_belief(double aCovDiag = 10.0) const;

  /**
   * Create a belief point for the state (augmented with imbalance, eccentricity, drag and inertia parameters),
   * with mean set to initial-motion.
   * \param aCovDiag The initial diagonal values for the covariant matrix (should be high).
   * \return A belief point for the state, with mean set to initial-motion and high covariance.
   */
  state_emdJ_belief_type get_init_state_emdJ_belief(
      double aCovDiag = 10.0) const;

  /**
   * Create a belief point for the input, assuming zero-mean.
   * \return A belief point for the input, with zero mean and input-disturbance covariance.
   */
  input_belief_type get_zero_input_belief() const;

  /**
   * Create a belief point for the measurement vector, assuming zero-mean.
   * \return A belief point for the measurement vector, with zero mean and measurement-noise covariance.
   */
  output_belief_type get_zero_output_belief() const;

  void imbue_names_for_received_meas(
      recorder::data_stream_options& data_opt) const;

  void imbue_names_for_generated_meas(
      recorder::data_stream_options& data_opt) const;
  void imbue_names_for_meas_stddevs(
      recorder::data_stream_options& data_opt) const;

  void imbue_names_for_state_estimates(
      recorder::data_stream_options& data_opt) const;
  void imbue_names_for_state_estimates_stddevs(
      recorder::data_stream_options& data_opt) const;

  /**
   * Default constructor.
   */
  satellite_model_options()
      : inertia_tensor(1.0, 0.0, 0.0, 1.0, 0.0, 1.0),
        input_disturbance(6, 0.0),
        measurement_noise(6, 0.0) {}

  virtual ~satellite_model_options() = default;

  /**
   * Loads all the configurations from the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive from which to load all the configurations.
   */
  void load_all_configs(const std::string& file_name);

  /**
   * Saves all the configurations to the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive to which to save all the configurations.
   */
  void save_all_configs(const std::string& file_name) const;

  /**
   * Loads the mass configurations (mass and inertia-tensor) from the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive from which to load the mass information.
   */
  void load_mass_configs(const std::string& file_name);

  /**
   * Saves the mass configurations (mass and inertia-tensor) to the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive to which to save the mass information.
   */
  void save_mass_configs(const std::string& file_name) const;

  /**
   * Loads the IMU configurations (orientation, location, Earth frame, and mag-field) from the
   * given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive from which to load the configs.
   */
  void load_IMU_configs(const std::string& file_name);

  /**
   * Saves the IMU configurations (orientation, location, Earth frame, and mag-field) to the
   * given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive to which to save the configs.
   */
  void save_IMU_configs(const std::string& file_name) const;

  /**
   * Loads the input-disturbance covariance matrix from the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive from which to load the matrix.
   */
  void load_input_disturbance(const std::string& file_name);

  /**
   * Saves the input-disturbance covariance matrix to the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive to which to save the matrix.
   */
  void save_input_disturbance(const std::string& file_name) const;

  /**
   * Loads the measurement-noise covariance matrix from the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive from which to load the matrix.
   */
  void load_measurement_noise(const std::string& file_name);

  /**
   * Saves the measurement-noise covariance matrix to the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive to which to save the matrix.
   */
  void save_measurement_noise(const std::string& file_name) const;

  /**
   * Loads the artificial measurement-noise covariance matrix from the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive from which to load the matrix.
   */
  void load_artificial_noise(const std::string& file_name);

  /**
   * Saves the artificial measurement-noise covariance matrix to the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive to which to save the matrix.
   */
  void save_artificial_noise(const std::string& file_name) const;

  /**
   * Loads the steady-state parameter covariance matrix from the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive from which to load the matrix.
   */
  void load_steady_param_covariance(const std::string& file_name);

  /**
   * Saves the steady-state parameter covariance matrix to the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive to which to save the matrix.
   */
  void save_steady_param_covariance(const std::string& file_name) const;

  /**
   * Loads the initial motion of the satellite from the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive from which to load the initial motion.
   */
  void load_initial_motion(const std::string& file_name);

  /**
   * Saves the initial motion of the satellite to the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive to which to save the initial motion.
   */
  void save_initial_motion(const std::string& file_name) const;
};

/**
 * This class stores a number of options related to the state prediction of satellite
 * systems (see satellite_invar_models.hpp).
 * \note This class is mainly intended to be used with satellite_predictor_po (program-options).
 */
struct satellite_predictor_options : satellite_model_options {

  using system_base_type = satellite_model_options::system_base_type;
  using system_gyro_type = satellite_model_options::system_gyro_type;
  using system_IMU_type = satellite_model_options::system_IMU_type;

  using input_type = satellite_model_options::input_type;
  using state_space_type = satellite_model_options::state_space_type;
  using temp_state_space_type = satellite_model_options::temp_state_space_type;
  using belief_space_type = satellite_model_options::belief_space_type;

  using input_cst_traj_type =
      pp::constant_trajectory<pp::vector_topology<input_type>>;
  using pred_factory_type = IKF_belief_transfer_factory<system_base_type>;
  using belief_pred_traj_type =
      belief_predicted_trajectory<belief_space_type, pred_factory_type,
                                  input_cst_traj_type>;
  using ML_pred_traj_type =
      pp::transformed_trajectory<temp_state_space_type, belief_pred_traj_type,
                                 maximum_likelihood_map>;

  /// Stores the prediction time horizon, for when using the satellite model for state predictions.
  double predict_time_horizon{100.0};

  /// Stores the threshold on the norm (Frobenius) of the state covariance matrix P, this threshold is used to start the
  /// prediction after sufficient online estimation.
  double predict_Pnorm_threshold{1.0};

  /// Distinguishes different types of prediction assumptions.
  enum predict_assumption_type {
    /// Assume that no measurements are made from the start of the predictions onward, i.e., pure a priori predictions.
    /// This gives the worst case uncertainty evaluations.
    no_measurements = 0,
    /// Assume that every measurement from the start of the predictions onward is equal to the most likely (mode)
    /// measurement given the a priori state prediction. This gives the best case uncertainty evaluations.
    most_likely_measurements,
    /// Discard the covariance calculations from the start of the predictions onward, i.e., do only state predictions,
    /// with no uncertainty evaluations.
    full_certainty
  };

  /// Stores the type of prediction assumption to be applied when computing the state predictions.
  predict_assumption_type predict_assumption{no_measurements};

 protected:
  void load_all_configs_impl(serialization::iarchive& in) override;
  void save_all_configs_impl(serialization::oarchive& out) const override;

  void load_predict_configs_impl(serialization::iarchive& in);
  void save_predict_configs_impl(serialization::oarchive& out) const;

 public:
  /**
   * Default constructor.
   */
  satellite_predictor_options()

      = default;

  /**
   * Loads the prediction configurations from the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive from which to load the prediction configurations.
   */
  void load_prediction_configs(const std::string& file_name);

  /**
   * Saves the prediction configurations to the given file-name (ReaK archive).
   * \param file_name The file-name of the ReaK archive to which to save the prediction configurations.
   */
  void save_prediction_configs(const std::string& file_name) const;
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_SYSTEMS_SATELLITE_MODELING_OPTIONS_H_

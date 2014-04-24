/**
 * \file satellite_modeling_options.hpp
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

#ifndef REAK_SATELLITE_MODELING_OPTIONS_HPP
#define REAK_SATELLITE_MODELING_OPTIONS_HPP

#include "satellite_invar_models.hpp"
#include "near_buoyant_airship_models.hpp"

#include "ctrl_sys/discrete_sss_concept.hpp"
#include "ctrl_sys/covar_topology.hpp"
#include "ctrl_sys/covariance_matrix.hpp"
#include "ctrl_sys/gaussian_belief_space.hpp"
#include "topologies/temporal_space.hpp"
#include "topologies/vector_topology.hpp"
#include "topologies/time_poisson_topology.hpp"
#include "ctrl_sys/invariant_kalman_filter.hpp"
#include "ctrl_sys/belief_state_predictor.hpp"
#include "ctrl_sys/maximum_likelihood_mapping.hpp"
#include "interpolation/constant_trajectory.hpp"
#include "path_planning/transformed_trajectory.hpp"

#include "lin_alg/mat_alg_symmetric.hpp"
#include "lin_alg/mat_alg_diagonal.hpp"

#include "serialization/archiver.hpp"
#include "recorders/data_record_options.hpp"

#include <string>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK's Data Recorders and Extractors */
namespace ctrl {


/**
 * This class stores a number of options related to the modeling of satellite 
 * systems (see satellite_invar_models.hpp).
 * \note This class is mainly intended to be used with satellite_modeling_po (program-options).
 */
struct satellite_model_options {
  
  
  typedef satellite3D_inv_dt_system system_base_type;
  typedef satellite3D_gyro_inv_dt_system system_gyro_type;
  typedef satellite3D_IMU_imdt_sys system_IMU_type;
  
  typedef airship3D_imdt_em_sys system_em_type;
  typedef airship3D_imdt_emd_sys system_emd_type;
  typedef airship3D_gyro_imdt_emd_sys system_gyro_emd_type;
  typedef airship3D_imdt_emdJ_sys system_emdJ_type;
  typedef airship3D_gyro_imdt_emdJ_sys system_gyro_emdJ_type;
  
  typedef discrete_sss_traits< system_base_type >::point_type state_type;
  typedef discrete_sss_traits< system_base_type >::input_type input_type;
  typedef discrete_sss_traits< system_base_type >::output_type output_type;
  
  typedef discrete_sss_traits< system_em_type >::point_type state_em_type;
  typedef discrete_sss_traits< system_emd_type >::point_type state_emd_type;
  typedef discrete_sss_traits< system_emdJ_type >::point_type state_emdJ_type;
  
  typedef covariance_matrix< vect_n<double> > covar_type;
  typedef covar_topology< covar_type > covar_space_type;
  
  typedef gaussian_belief_state< state_type,  covar_type > state_belief_type;
  typedef gaussian_belief_state< input_type,  covar_type > input_belief_type;
  typedef gaussian_belief_state< output_type, covar_type > output_belief_type;
  
  typedef gaussian_belief_state< state_em_type, covar_type > state_em_belief_type;
  typedef gaussian_belief_state< state_emd_type, covar_type > state_emd_belief_type;
  typedef gaussian_belief_state< state_emdJ_type, covar_type > state_emdJ_belief_type;
  
  
  typedef system_base_type::state_space_type state_space_type;
  typedef pp::temporal_space<state_space_type, pp::time_poisson_topology, pp::time_distance_only> temp_state_space_type;
  typedef gaussian_belief_space<state_space_type, covar_space_type> belief_space_type;
  
  typedef system_em_type::state_space_type state_space_em_type;
  typedef pp::temporal_space<state_space_em_type, pp::time_poisson_topology, pp::time_distance_only> temp_state_space_em_type;
  typedef gaussian_belief_space<state_space_em_type, covar_space_type> belief_space_em_type;
  
  typedef system_emd_type::state_space_type state_space_emd_type;
  typedef pp::temporal_space<state_space_emd_type, pp::time_poisson_topology, pp::time_distance_only> temp_state_space_emd_type;
  typedef gaussian_belief_space<state_space_emd_type, covar_space_type> belief_space_emd_type;
  
  typedef system_emdJ_type::state_space_type state_space_emdJ_type;
  typedef pp::temporal_space<state_space_emdJ_type, pp::time_poisson_topology, pp::time_distance_only> temp_state_space_emdJ_type;
  typedef gaussian_belief_space<state_space_emdJ_type, covar_space_type> belief_space_emdJ_type;
  
  
  /// Stores the time-step of the system (if discrete-time).
  double time_step;
  
  /// Stores the mass of the satellite (in kg).
  double mass;
  
  /// Stores the inertia-tensor of the satellite (in body-fixed frame, in kg-m2).
  mat<double,mat_structure::symmetric> inertia_tensor;
  
  /// Stores the input-disturbance covariance matrix (Q), disturbance on the net wrench on the system.
  mat<double,mat_structure::diagonal> input_disturbance;
  
  /// Stores the measurement-noise covariance matrix (R).
  mat<double,mat_structure::diagonal> measurement_noise;
  
  /// Stores the artificial measurement-noise covariance matrix (R). This is sometimes used to add noise to measurements (for testing the filtering methods).
  mat<double,mat_structure::diagonal> artificial_noise;
  
  /// Stores the steady-state parameter identification covariance matrix (Pa). This is to be used in combination with the TSOSAKF methods.
  mat<double,mat_structure::diagonal> steady_param_covariance;
  
  /// Stores the orientation of the IMU with respect to the body-fixed frame of the satellite.
  unit_quat<double> IMU_orientation;
  /// Stores the position of the IMU with respect to the body-fixed frame of the satellite.
  vect<double,3>    IMU_location;
  /// Stores the orientation of the 'Earth' frame with respect to the global inertial frame used.
  unit_quat<double> earth_orientation;
  /// Stores the expected direction of the Earth's magnetic field lines at the location of the satellite, and relative to the Earth-frame.
  vect<double,3>    mag_field_direction;
  
  /// Stores the initial motion of the satellite.
  frame_3D<double> initial_motion;
  
  enum available_measurements {
    pose_measures = 0,
    gyro_measures = 16,
    IMU_measures = 48
  };
  
  enum model_kind {
    invariant = 0,
    invar_mom,
    invar_mom2,
    invar_mom_em,
    invar_mom_emd,
    invar_mom_emdJ
  };
  
  enum filtering_options {
    plain_KF = 0,
    TSOSAKF = 1024
  };
  
  /// Stores the kind of system used to model the satellite (OR-combination of 'available_measurements' 'model_kind').
  int system_kind;
  
  std::size_t get_measurement_count() const {
    switch(system_kind & 48) {
      case 16:
        return 10;
      case 48:
        return 16;
      default:
        return 7;
    };
  };
  
  std::size_t get_meas_error_count() const {
    switch(system_kind & 48) {
      case 16:
        return 9;
      case 48:
        return 15;
      default:
        return 6;
    };
  };
  
  std::size_t get_total_inv_state_count() const {
    std::size_t result = 0;
    switch(system_kind & 15) {
      case invar_mom_em:
        result = 16;
        break;
      case invar_mom_emd:
        result = 18;
        break;
      case invar_mom_emdJ:
        result = 24;
        break;
      case invariant:
      case invar_mom2:
      case invar_mom:
      default:
        result = 12;
        break;
    };
    return result;
  };
  
  std::size_t get_actual_inv_state_count() const {
    return 12;
  };
  
  std::string get_kf_accronym() const {
    std::string result;
    
    std::string kf_method = "imkf";
    if(system_kind & TSOSAKF)
      kf_method = "tsosakf";
    
    switch(system_kind & 15) {
      case invariant:
        result = "iekf";
        break;
      case invar_mom2:
        result = "imkfv2";
        break;
      case invar_mom_em:
        result = kf_method + "_em";
        break;
      case invar_mom_emd:
        result = kf_method + "_emd";
        break;
      case invar_mom_emdJ:
        result = kf_method + "_emdJ";
        break;
      case invar_mom:
      default:
        result = "imkf";
        break;
    };
    switch(system_kind & 48) {
      case 16:
        result += "_gyro";
        break;
      case 48:
        result += "_IMU";
        break;
      default:
        break;
    };
    return result;
  };
  
  std::string get_sys_abbreviation() const {
    std::string result;
    switch(system_kind & 15) {
      case invariant:
        result = "inv";
        break;
      case invar_mom2:
        result = "invmid";
        break;
      case invar_mom_em:
        result = "invmid_em";
        break;
      case invar_mom_emd:
        result = "invmid_emd";
        break;
      case invar_mom_emdJ:
        result = "invmid_emdJ";
        break;
      case invar_mom:
      default:
        result = "invmom";
        break;
    };
    switch(system_kind & 48) {
      case 16:
        result += "_gyro";
        break;
      case 48:
        result += "_IMU";
        break;
      default:
        break;
    };
    return result;
  };
  
  
protected:
  
  virtual void load_all_configs_impl(serialization::iarchive& in);
  virtual void save_all_configs_impl(serialization::oarchive& out) const;
  
public:
  
  /**
   * Constructs a base satellite system.
   * \return A newly created base satellite system (as shared-pointer).
   */
  shared_ptr< system_base_type > get_base_sat_system() const;
  
  /**
   * Constructs a satellite system (with gyro).
   * \return A newly created gyro satellite system (as shared-pointer).
   */
  shared_ptr< system_gyro_type > get_gyro_sat_system() const;
  
  /**
   * Constructs a satellite system (with IMU).
   * \return A newly created IMU satellite system (as shared-pointer).
   */
  shared_ptr< system_IMU_type > get_IMU_sat_system() const;
  
  /**
   * Constructs a eccentricity-imbalance airship system.
   * \return A newly created eccentricity-imbalance airship system (as shared-pointer).
   */
  shared_ptr< system_em_type > get_em_airship_system() const;
  
  /**
   * Constructs a eccentricity-imbalance-drag airship system.
   * \return A newly created eccentricity-imbalance-drag airship system (as shared-pointer).
   */
  shared_ptr< system_emd_type > get_emd_airship_system() const;
  
  /**
   * Constructs a eccentricity-imbalance-drag airship system with gyro measurements.
   * \return A newly created eccentricity-imbalance-drag airship system with gyro measurements (as shared-pointer).
   */
  shared_ptr< system_gyro_emd_type > get_gyro_emd_airship_system() const;
  
  /**
   * Constructs a eccentricity-imbalance-drag-inertia airship system.
   * \return A newly created eccentricity-imbalance-drag-inertia airship system (as shared-pointer).
   */
  shared_ptr< system_emdJ_type > get_emdJ_airship_system() const;
  
  /**
   * Constructs a eccentricity-imbalance-drag-inertia airship system with gyro measurements.
   * \return A newly created eccentricity-imbalance-drag-inertia airship system with gyro measurements (as shared-pointer).
   */
  shared_ptr< system_gyro_emdJ_type > get_gyro_emdJ_airship_system() const;
  
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
  state_emdJ_belief_type get_init_state_emdJ_belief(double aCovDiag = 10.0) const;
  
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
  
  void imbue_names_for_received_meas(recorder::data_stream_options& data_opt) const;
  
  void imbue_names_for_generated_meas(recorder::data_stream_options& data_opt) const;
  void imbue_names_for_meas_stddevs(recorder::data_stream_options& data_opt) const;
  
  void imbue_names_for_state_estimates(recorder::data_stream_options& data_opt) const;
  void imbue_names_for_state_estimates_stddevs(recorder::data_stream_options& data_opt) const;
  
  /**
   * Default constructor.
   */
  satellite_model_options() : time_step(0.01), mass(1.0), 
                              inertia_tensor(1.0, 0.0, 0.0, 1.0, 0.0, 1.0),
                              input_disturbance(6, 0.0),
                              measurement_noise(6, 0.0),
                              IMU_orientation(), 
                              IMU_location(), 
                              earth_orientation(), 
                              mag_field_direction(),
                              initial_motion(),
                              system_kind(invariant | pose_measures) { };
  
  virtual ~satellite_model_options() { };
  
  /**
   * Loads all the configurations from the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive from which to load all the configurations.
   */
  void load_all_configs(const std::string& aFileName);
  
  /**
   * Saves all the configurations to the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive to which to save all the configurations.
   */
  void save_all_configs(const std::string& aFileName) const;
  
  
  /**
   * Loads the mass configurations (mass and inertia-tensor) from the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive from which to load the mass information.
   */
  void load_mass_configs(const std::string& aFileName);
  
  /**
   * Saves the mass configurations (mass and inertia-tensor) to the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive to which to save the mass information.
   */
  void save_mass_configs(const std::string& aFileName) const;
  
  /**
   * Loads the IMU configurations (orientation, location, Earth frame, and mag-field) from the 
   * given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive from which to load the configs.
   */
  void load_IMU_configs(const std::string& aFileName);
  
  /**
   * Saves the IMU configurations (orientation, location, Earth frame, and mag-field) to the 
   * given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive to which to save the configs.
   */
  void save_IMU_configs(const std::string& aFileName) const;
  
  /**
   * Loads the input-disturbance covariance matrix from the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive from which to load the matrix.
   */
  void load_input_disturbance(const std::string& aFileName);
  
  /**
   * Saves the input-disturbance covariance matrix to the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive to which to save the matrix.
   */
  void save_input_disturbance(const std::string& aFileName) const;
  
  /**
   * Loads the measurement-noise covariance matrix from the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive from which to load the matrix.
   */
  void load_measurement_noise(const std::string& aFileName);
  
  /**
   * Saves the measurement-noise covariance matrix to the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive to which to save the matrix.
   */
  void save_measurement_noise(const std::string& aFileName) const;
  
  /**
   * Loads the artificial measurement-noise covariance matrix from the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive from which to load the matrix.
   */
  void load_artificial_noise(const std::string& aFileName);
  
  /**
   * Saves the artificial measurement-noise covariance matrix to the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive to which to save the matrix.
   */
  void save_artificial_noise(const std::string& aFileName) const;
  
  /**
   * Loads the steady-state parameter covariance matrix from the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive from which to load the matrix.
   */
  void load_steady_param_covariance(const std::string& aFileName);
  
  /**
   * Saves the steady-state parameter covariance matrix to the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive to which to save the matrix.
   */
  void save_steady_param_covariance(const std::string& aFileName) const;
  
  /**
   * Loads the initial motion of the satellite from the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive from which to load the initial motion.
   */
  void load_initial_motion(const std::string& aFileName);
  
  /**
   * Saves the initial motion of the satellite to the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive to which to save the initial motion.
   */
  void save_initial_motion(const std::string& aFileName) const;
  
  
  
};



/**
 * This class stores a number of options related to the state prediction of satellite 
 * systems (see satellite_invar_models.hpp).
 * \note This class is mainly intended to be used with satellite_predictor_po (program-options).
 */
struct satellite_predictor_options : satellite_model_options {
  
  typedef satellite_model_options::system_base_type system_base_type;
  typedef satellite_model_options::system_gyro_type system_gyro_type;
  typedef satellite_model_options::system_IMU_type system_IMU_type;
  
  typedef satellite_model_options::input_type input_type;
  typedef satellite_model_options::state_space_type state_space_type;
  typedef satellite_model_options::temp_state_space_type temp_state_space_type;
  typedef satellite_model_options::belief_space_type belief_space_type;
  
  typedef pp::constant_trajectory< pp::vector_topology< input_type > > input_cst_traj_type;
  typedef IKF_belief_transfer_factory< system_base_type > pred_factory_type;
  typedef belief_predicted_trajectory< belief_space_type, pred_factory_type, input_cst_traj_type > belief_pred_traj_type;
  typedef pp::transformed_trajectory< temp_state_space_type, belief_pred_traj_type, maximum_likelihood_map> ML_pred_traj_type;
  
  /// Stores the prediction time horizon, for when using the satellite model for state predictions.
  double predict_time_horizon;
  
  /// Stores the threshold on the norm (Frobenius) of the state covariance matrix P, this threshold is used to start the prediction after sufficient online estimation.
  double predict_Pnorm_threshold;
  
  /// Distinguishes different types of prediction assumptions.
  enum predict_assumption_type {
    /// Assume that no measurements are made from the start of the predictions onward, i.e., pure a priori predictions. This gives the worst case uncertainty evaluations.
    no_measurements = 0,
    /// Assume that every measurement from the start of the predictions onward is equal to the most likely (mode) measurement given the a priori state prediction. This gives the best case uncertainty evaluations.
    most_likely_measurements,
    /// Discard the covariance calculations from the start of the predictions onward, i.e., do only state predictions, with no uncertainty evaluations.
    full_certainty
  };
  
  /// Stores the type of prediction assumption to be applied when computing the state predictions.
  predict_assumption_type predict_assumption;
  
protected:
  
  void load_all_configs_impl(serialization::iarchive& in);
  void save_all_configs_impl(serialization::oarchive& out) const;
  
  void load_predict_configs_impl(serialization::iarchive& in);
  void save_predict_configs_impl(serialization::oarchive& out) const;
  
public:
  
  /**
   * Default constructor.
   */
  satellite_predictor_options() : predict_time_horizon(100.0), 
                                  predict_Pnorm_threshold(1.0), 
                                  predict_assumption(no_measurements) { };
  
  
  /**
   * Loads the prediction configurations from the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive from which to load the prediction configurations.
   */
  void load_prediction_configs(const std::string& aFileName);
  
  /**
   * Saves the prediction configurations to the given file-name (ReaK archive).
   * \param aFileName The file-name of the ReaK archive to which to save the prediction configurations.
   */
  void save_prediction_configs(const std::string& aFileName) const;
  
};


};

};

#endif // REAK_SATELLITE_MODELING_OPTIONS_HPP













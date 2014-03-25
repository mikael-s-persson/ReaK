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

#include "lin_alg/mat_alg_symmetric.hpp"
#include "lin_alg/mat_alg_diagonal.hpp"

#include "serialization/archiver.hpp"

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
    gyro_measures = 8,
    IMU_measures = 24
  };
  
  enum model_kind {
    invariant = 0,
    invar_mom,
    invar_mom2
  };
  
  /// Stores the kind of system used to model the satellite (OR-combination of 'available_measurements' 'model_kind').
  int system_kind;
  
  std::size_t get_measurement_count() const {
    switch(system_kind & 24) {
      case 8:
        return 10;
      case 24:
        return 16;
      default:
        return 7;
    };
  };
  
  std::size_t get_meas_error_count() const {
    switch(system_kind & 24) {
      case 8:
        return 9;
      case 24:
        return 15;
      default:
        return 6;
    };
  };
  
protected:
  
  virtual void load_all_configs_impl(serialization::iarchive& in);
  virtual void save_all_configs_impl(serialization::oarchive& out) const;
  
public:
  
  
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













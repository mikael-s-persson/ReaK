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
  
  /**
   * Default constructor.
   */
  satellite_model_options() : time_step(0.01), mass(1.0), 
                              inertia_tensor(1.0, 0.0, 0.0, 1.0, 0.0, 1.0),
                              input_disturbance(6, 0.0),
                              measurement_noise(6, 0.0),
                              IMU_orientation(), 
                              IMU_location(), 
                              Earth_orientation(), 
                              mag_field_direction(),
                              initial_motion() { };
  
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



};

};

#endif // REAK_SATELLITE_MODELING_OPTIONS_HPP













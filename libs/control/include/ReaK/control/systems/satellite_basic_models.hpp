/**
 * \file satellite_basic_models.hpp
 *
 * This library contains a number of discrete-time state-space systems to describe the dynamics
 * of a free-floating satellite. These are basic and simplified models, with no complex forces
 * applied, just free-floating dynamics with 6 dof actuation forces applied to it. These systems
 * benefit from a special integration method called the "momentum-conserving trapezoidal method" (TRAPM),
 * which is an invariant variational method that guarantees conservation of angular momentum
 * when no actuation is applied, i.e., it is an efficient and highly stable method.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date October 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_SATELLITE_BASIC_MODELS_HPP
#define REAK_SATELLITE_BASIC_MODELS_HPP

#include "ReaK/core/base/named_object.hpp"

#include "ReaK/control/estimators/covar_topology.hpp"
#include "ReaK/control/estimators/covariance_matrix.hpp"
#include "ReaK/control/estimators/gaussian_belief_space.hpp"
#include "ReaK/control/systems/invariant_system_concept.hpp"
#include "ReaK/topologies/spaces/se3_topologies.hpp"
#include "ReaK/topologies/spaces/temporal_space.hpp"
#include "ReaK/topologies/spaces/time_poisson_topology.hpp"

#include "ReaK/math/lin_alg/mat_alg.hpp"

#include <type_traits>

namespace ReaK::ctrl {

/**
 * This class implements a basic linearized discrete-time state-space system for
 * simple free-floating dynamics characteristic of satellites in orbit. This system
 * benefits from a special integration method called the "momentum-conserving trapezoidal method" (TRAPM),
 * which is an invariant variational method that guarantees conservation of angular momentum
 * when no actuation is applied, i.e., it is an efficient and highly stable method.
 * Also, this system operates within a first-order (once-differentiable) SE(3) topology.
 * \note THIS CLASS SHOULD NOT BE USED FOR THE LINEARIZATION / KALMAN-FILTERING.
 */
class satellite3D_lin_dt_system : public named_object {
 public:
  using state_space_type = pp::se3_1st_order_topology_t<double>;

  using point_type = pp::topology_traits<state_space_type>::point_type;
  using point_difference_type =
      pp::topology_traits<state_space_type>::point_difference_type;
  using point_derivative_type =
      pp::topology_traits<state_space_type>::point_difference_type;

  using time_type = double;
  using time_difference_type = double;

  using input_type = vect_n<double>;
  using output_type = vect_n<double>;

  static constexpr std::size_t dimensions = 13;
  static constexpr std::size_t input_dimensions = 6;
  static constexpr std::size_t output_dimensions = 7;

  using matrixA_type = mat<double, mat_structure::square>;
  using matrixB_type = mat<double, mat_structure::rectangular>;
  using matrixC_type = mat<double, mat_structure::rectangular>;
  using matrixD_type = mat<double, mat_structure::rectangular>;

  struct zero_input_trajectory {
    auto get_point(time_type /*unused*/) const {
      return input_type(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
  };

  using covar_type = covariance_matrix<vect_n<double>>;
  using covar_space_type = covar_topology<covar_type>;
  using temporal_state_space_type =
      pp::temporal_space<state_space_type, pp::time_poisson_topology,
                         pp::time_distance_only>;
  using belief_space_type =
      gaussian_belief_space<state_space_type, covar_space_type>;
  using temporal_belief_space_type =
      pp::temporal_space<belief_space_type, pp::time_poisson_topology,
                         pp::time_distance_only>;
  using state_belief_type = gaussian_belief_state<point_type, covar_type>;
  using input_belief_type = gaussian_belief_state<input_type, covar_type>;
  using output_belief_type = gaussian_belief_state<output_type, covar_type>;

  virtual std::shared_ptr<temporal_state_space_type> get_temporal_state_space(
      double aStartTime = 0.0, double aEndTime = 1.0) const;
  virtual std::shared_ptr<state_space_type> get_state_space() const;

  virtual std::shared_ptr<temporal_belief_space_type> get_temporal_belief_space(
      double aStartTime = 0.0, double aEndTime = 1.0) const;
  virtual std::shared_ptr<belief_space_type> get_belief_space() const;

  virtual state_belief_type get_zero_state_belief(
      double aCovValue = 10.0) const;
  virtual input_belief_type get_zero_input_belief(double aCovValue = 1.0) const;
  virtual output_belief_type get_zero_output_belief(
      double aCovValue = 1.0) const;

 protected:
  double mMass;
  mat<double, mat_structure::symmetric> mInertiaMoment;
  mat<double, mat_structure::symmetric> mInertiaMomentInv;
  time_difference_type mDt;

 public:
  /**
   * Returns the dimensions of the states of the system.
   * \return The dimensions of the states of the system.
   */
  virtual std::size_t get_state_dimensions() const { return 13; }

  /**
   * Returns the dimensions of the input of the system.
   * \return The dimensions of the input of the system.
   */
  virtual std::size_t get_input_dimensions() const { return 6; }

  /**
   * Returns the dimensions of the output of the system.
   * \return The dimensions of the output of the system.
   */
  virtual std::size_t get_output_dimensions() const { return 7; }

  /**
   * Constructor.
   * \param aName The name for this object.
   * \param aMass The mass of the satellite.
   * \param aInertiaMoment The inertia tensor of the satellite.
   * \param aDt The time-step for this discrete-time system.
   */
  explicit satellite3D_lin_dt_system(
      const std::string& aName, double aMass = 1.0,
      const mat<double, mat_structure::symmetric>& aInertiaMoment =
          (mat<double, mat_structure::symmetric>(
              mat<double, mat_structure::identity>(3))),
      double aDt = 0.001);

  satellite3D_lin_dt_system() : satellite3D_lin_dt_system("") {}

  /**
   * This function returns the time-step for this discrete-time system.
   * \return The time-step for this discrete-time system.
   */
  time_difference_type get_time_step() const { return mDt; }

  /**
   * This function sets the time-step for this discrete-time system.
   * \param aDt The new time-step for this discrete-time system.
   */
  virtual void set_time_step(time_difference_type aDt) { mDt = aDt; }

  /**
   * This function computes the next state of the system, i.e., the state at one time-step after the current time.
   * \param space The state-space within which the states reside.
   * \param x The current state of the system.
   * \param u The current input being applied to the system.
   * \param t The current time.
   * \return The state after one time-step beyond the given current state of the system.
   */
  virtual point_type get_next_state(const state_space_type& space,
                                    const point_type& x, const input_type& u,
                                    const time_type& t = 0.0) const;

  /**
   * This function computes the linearization of the state-transitions of the system.
   * In other words, it populates the system matrices with the values appropriate for
   * the given state-transition.
   * \param A Holds, as output, the state-to-state jacobian matrix of the state-transition of the system.
   * \param B Holds, as output, the input-to-state jacobian matrix of the state-transition of the system.
   * \param space The state-space within which the states reside.
   * \param t_0 The time before the state-transition occurred.
   * \param t_1 The time after the state-transition occurred.
   * \param p_0 The state before the state-transition occurred.
   * \param p_1 The state after the state-transition occurred.
   * \param u_0 The input before the state-transition occurred.
   * \param u_1 The input after the state-transition occurred.
   */
  virtual void get_state_transition_blocks(
      matrixA_type& A, matrixB_type& B, const state_space_type& space,
      const time_type& t_0, const time_type& t_1, const point_type& p_0,
      const point_type& p_1, const input_type& u_0,
      const input_type& u_1) const;

  /**
   * This function computes the output of the system corresponding to the current state.
   * \param space The state-space within which the states reside.
   * \param x The current state of the system.
   * \param u The current input being applied to the system.
   * \param t The current time.
   * \return The output for the given current state of the system.
   */
  virtual output_type get_output(const state_space_type& space,
                                 const point_type& x, const input_type& u,
                                 const time_type& t = 0.0) const;

  /**
   * This function computes the linearization of the output-function of the system.
   * In other words, it populates the system matrices with the values appropriate at
   * the given state.
   * \param C Holds, as output, the state-to-output jacobian matrix of the output-function of the system.
   * \param D Holds, as output, the input-to-output jacobian matrix of the output-function of the system.
   * \param space The state-space within which the states reside.
   * \param t The current time.
   * \param p The current state of the system.
   * \param u The input at the current time.
   */
  virtual void get_output_function_blocks(matrixC_type& C, matrixD_type& D,
                                          const state_space_type& space,
                                          const time_type& t,
                                          const point_type& p,
                                          const input_type& u) const;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override;
  void load(ReaK::serialization::iarchive& A, unsigned int /*unused*/) override;

  RK_RTTI_MAKE_CONCRETE_1BASE(satellite3D_lin_dt_system, 0xC2310013, 1,
                              "satellite3D_lin_dt_system", named_object)
};

/**
 * This class implements a basic linearized discrete-time state-space system for
 * simple free-floating dynamics characteristic of satellites in orbit. In this class,
 * the measurements (output) from gyros is incorporated into the system, that is, measurements
 * of the angular velocity is available. This system benefits from a special integration method
 * called the "momentum-conserving trapezoidal method" (TRAPM), which is an invariant variational
 * method that guarantees conservation of angular momentum when no actuation is applied, i.e.,
 * it is an efficient and highly stable method. Also, this system operates within a first-order
 * (once-differentiable) SE(3) topology.
 * \note THIS CLASS SHOULD NOT BE USED FOR THE LINEARIZATION / KALMAN-FILTERING.
 */
class satellite3D_gyro_lin_dt_system : public satellite3D_lin_dt_system {
 public:
  using state_space_type = satellite3D_lin_dt_system::state_space_type;

  using covar_type = satellite3D_lin_dt_system::covar_type;
  using covar_space_type = satellite3D_lin_dt_system::covar_space_type;
  using temporal_state_space_type =
      satellite3D_lin_dt_system::temporal_state_space_type;
  using belief_space_type = satellite3D_lin_dt_system::belief_space_type;
  using temporal_belief_space_type =
      satellite3D_lin_dt_system::temporal_belief_space_type;
  using state_belief_type = satellite3D_lin_dt_system::state_belief_type;
  using input_belief_type = satellite3D_lin_dt_system::input_belief_type;
  using output_belief_type = satellite3D_lin_dt_system::output_belief_type;

  using point_type = satellite3D_lin_dt_system::point_type;
  using point_difference_type =
      satellite3D_lin_dt_system::point_difference_type;
  using point_derivative_type =
      satellite3D_lin_dt_system::point_derivative_type;

  using time_type = satellite3D_lin_dt_system::time_type;
  using time_difference_type = satellite3D_lin_dt_system::time_difference_type;

  using input_type = satellite3D_lin_dt_system::input_type;
  using output_type = satellite3D_lin_dt_system::output_type;

  static constexpr std::size_t dimensions = 13;
  static constexpr std::size_t input_dimensions = 6;
  static constexpr std::size_t output_dimensions = 10;

  using matrixA_type = satellite3D_lin_dt_system::matrixA_type;
  using matrixB_type = satellite3D_lin_dt_system::matrixB_type;
  using matrixC_type = satellite3D_lin_dt_system::matrixC_type;
  using matrixD_type = satellite3D_lin_dt_system::matrixD_type;

  using zero_input_trajectory =
      satellite3D_lin_dt_system::zero_input_trajectory;

  output_belief_type get_zero_output_belief(
      double aCovValue = 1.0) const override;

  /**
   * Returns the dimensions of the output of the system.
   * \return The dimensions of the output of the system.
   */
  std::size_t get_output_dimensions() const override { return 10; }

  /**
   * Constructor.
   * \param aName The name for this object.
   * \param aMass The mass of the satellite.
   * \param aInertiaMoment The inertia tensor of the satellite.
   * \param aDt The time-step for this discrete-time system.
   */
  explicit satellite3D_gyro_lin_dt_system(
      const std::string& aName, double aMass = 1.0,
      const mat<double, mat_structure::symmetric>& aInertiaMoment =
          (mat<double, mat_structure::symmetric>(
              mat<double, mat_structure::identity>(3))),
      double aDt = 0.001);

  satellite3D_gyro_lin_dt_system() : satellite3D_gyro_lin_dt_system("") {}

  output_type get_output(const state_space_type& /*space*/, const point_type& x,
                         const input_type& u,
                         const time_type& t = 0.0) const override;

  void get_output_function_blocks(matrixC_type& C, matrixD_type& D,
                                  const state_space_type& /*space*/,
                                  const time_type& /*t*/,
                                  const point_type& /*p*/,
                                  const input_type& /*u*/) const override;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    satellite3D_lin_dt_system::save(
        A, satellite3D_lin_dt_system::getStaticObjectType()->TypeVersion());
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    satellite3D_lin_dt_system::load(
        A, satellite3D_lin_dt_system::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(satellite3D_gyro_lin_dt_system, 0xC2310018, 1,
                              "satellite3D_gyro_lin_dt_system",
                              satellite3D_lin_dt_system)
};

/**
 * This class implements an invariantized discrete-time state-space system for
 * simple free-floating dynamics characteristic of satellites in orbit. This system
 * benefits from a special integration method called the "momentum-conserving trapezoidal method" (TRAPM),
 * which is an invariant variational method that guarantees conservation of angular momentum
 * when no actuation is applied, i.e., it is an efficient and highly stable method.
 * Also, this system operates within a first-order (once-differentiable) SE(3) topology.
 */
class satellite3D_inv_dt_system : public satellite3D_lin_dt_system {
 public:
  using state_space_type = satellite3D_lin_dt_system::state_space_type;

  using covar_type = satellite3D_lin_dt_system::covar_type;
  using covar_space_type = satellite3D_lin_dt_system::covar_space_type;
  using temporal_state_space_type =
      satellite3D_lin_dt_system::temporal_state_space_type;
  using belief_space_type = satellite3D_lin_dt_system::belief_space_type;
  using temporal_belief_space_type =
      satellite3D_lin_dt_system::temporal_belief_space_type;
  using state_belief_type = satellite3D_lin_dt_system::state_belief_type;
  using input_belief_type = satellite3D_lin_dt_system::input_belief_type;
  using output_belief_type = satellite3D_lin_dt_system::output_belief_type;

  using point_type = satellite3D_lin_dt_system::point_type;
  using point_difference_type =
      satellite3D_lin_dt_system::point_difference_type;
  using point_derivative_type =
      satellite3D_lin_dt_system::point_derivative_type;

  using time_type = satellite3D_lin_dt_system::time_type;
  using time_difference_type = satellite3D_lin_dt_system::time_difference_type;

  using input_type = satellite3D_lin_dt_system::input_type;
  using output_type = satellite3D_lin_dt_system::output_type;

  using invariant_error_type = vect_n<double>;
  using invariant_correction_type = vect_n<double>;
  using invariant_frame_type = mat<double, mat_structure::square>;

  static constexpr std::size_t dimensions = 13;
  static constexpr std::size_t input_dimensions = 6;
  static constexpr std::size_t output_dimensions = 7;
  static constexpr std::size_t invariant_error_dimensions = 6;
  static constexpr std::size_t invariant_correction_dimensions = 12;

  using matrixA_type = satellite3D_lin_dt_system::matrixA_type;
  using matrixB_type = satellite3D_lin_dt_system::matrixB_type;
  using matrixC_type = satellite3D_lin_dt_system::matrixC_type;
  using matrixD_type = satellite3D_lin_dt_system::matrixD_type;

  using zero_input_trajectory =
      satellite3D_lin_dt_system::zero_input_trajectory;

  std::shared_ptr<temporal_belief_space_type> get_temporal_belief_space(
      double aStartTime = 0.0, double aEndTime = 1.0) const override;
  std::shared_ptr<belief_space_type> get_belief_space() const override;

  state_belief_type get_zero_state_belief(
      double aCovValue = 10.0) const override;
  output_belief_type get_zero_output_belief(
      double aCovValue = 1.0) const override;

  /**
   * Returns the dimensions of the states of the system.
   * \return The dimensions of the states of the system.
   */
  std::size_t get_state_dimensions() const override { return 13; }

  /**
   * Returns the dimensions of the input of the system.
   * \return The dimensions of the input of the system.
   */
  std::size_t get_input_dimensions() const override { return 6; }

  /**
   * Returns the dimensions of the output of the system.
   * \return The dimensions of the output of the system.
   */
  std::size_t get_output_dimensions() const override { return 7; }

  /**
   * Returns the dimensions of the invariant errors of the system.
   * \return The dimensions of the invariant errors of the system.
   */
  virtual std::size_t get_invariant_error_dimensions() const { return 6; }

  /**
   * Returns the dimensions of the corrections to the states of the system.
   * \return The dimensions of the corrections to the states of the system.
   */
  virtual std::size_t get_correction_dimensions() const { return 12; }

  /**
   * Constructor.
   * \param aName The name for this object.
   * \param aMass The mass of the satellite.
   * \param aInertiaMoment The inertia tensor of the satellite.
   * \param aDt The time-step for this discrete-time system.
   */
  explicit satellite3D_inv_dt_system(
      const std::string& aName, double aMass = 1.0,
      const mat<double, mat_structure::symmetric>& aInertiaMoment =
          (mat<double, mat_structure::symmetric>(
              mat<double, mat_structure::identity>(3))),
      double aDt = 0.001);

  satellite3D_inv_dt_system() : satellite3D_inv_dt_system("") {}

  void get_state_transition_blocks(
      matrixA_type& A, matrixB_type& B, const state_space_type& /*space*/,
      const time_type& t_0, const time_type& /*t_1*/, const point_type& p_0,
      const point_type& /*p_1*/, const input_type& /*u_0*/,
      const input_type& /*u_1*/) const override;

  void get_output_function_blocks(matrixC_type& C, matrixD_type& D,
                                  const state_space_type& /*space*/,
                                  const time_type& /*t*/,
                                  const point_type& /*p*/,
                                  const input_type& /*u*/) const override;

  /**
   * This function computes the invariant output-error of the system corresponding to the current state and the given
   * output.
   * \param space The state-space within which the states reside.
   * \param x The current state of the system.
   * \param u The current input being applied to the system.
   * \param y The output against which to compute the invariant error.
   * \param t The current time.
   * \return The invariant output-error for the given state and output.
   */
  virtual invariant_error_type get_invariant_error(
      const state_space_type& space, const point_type& x, const input_type& u,
      const output_type& y, const time_type& t) const;

  /**
   * This function computes a state corresponding to the given state corrected by a given invariant term.
   * \param space The state-space within which the states reside.
   * \param x The current state of the system.
   * \param c The invariant correction term to apply to the state.
   * \param u The current input being applied to the system.
   * \param t The current time.
   * \return The corrected state of the system.
   */
  virtual point_type apply_correction(const state_space_type& space,
                                      const point_type& x,
                                      const invariant_correction_type& c,
                                      const input_type& u,
                                      const time_type& t) const;

  /**
   * This function computes the invariant frame transition matrix for the prior stage,
   * i.e., during a state transition from x_0 to x_1, what invariant frame transition matrix
   * describes the shift from one frame to the other.
   * \param space The state-space within which the states reside.
   * \param x_0 The state of the system before the state-transition.
   * \param x_1 The state of the system after the state-transition.
   * \param u The input being applied to the system before the state-transition.
   * \param t The time before the state-transition.
   * \return The invariant frame transition matrix for the prior stage.
   */
  virtual invariant_frame_type get_invariant_prior_frame(
      const state_space_type& space, const point_type& x_0,
      const point_type& x_1, const input_type& u, const time_type& t) const {
    return invariant_frame_type(
        mat<double, mat_structure::identity>(invariant_correction_dimensions));
  }

  /**
   * This function computes the invariant frame transition matrix for the posterior stage,
   * i.e., during a state correction from x_0 to x_1, what invariant frame transition matrix
   * describes the shift from one frame to the other.
   * \param space The state-space within which the states reside.
   * \param x_0 The state of the system before the correction.
   * \param x_1 The state of the system after the correction.
   * \param u The current input being applied to the system.
   * \param t The current time.
   * \return The invariant frame transition matrix for the posterior stage.
   */
  virtual invariant_frame_type get_invariant_posterior_frame(
      const state_space_type& space, const point_type& x_0,
      const point_type& x_1, const input_type& u, const time_type& t) const {
    return invariant_frame_type(
        mat<double, mat_structure::identity>(invariant_correction_dimensions));
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    satellite3D_lin_dt_system::save(
        A, satellite3D_lin_dt_system::getStaticObjectType()->TypeVersion());
  };
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    satellite3D_lin_dt_system::load(
        A, satellite3D_lin_dt_system::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(satellite3D_inv_dt_system, 0xC2310014, 1,
                              "satellite3D_inv_dt_system",
                              satellite3D_lin_dt_system)
};

template <>
struct is_invariant_system<satellite3D_inv_dt_system> : std::true_type {};

/**
 * This class implements an invariantized discrete-time state-space system for
 * simple free-floating dynamics characteristic of satellites in orbit. In this class,
 * the measurements (output) from gyros is incorporated into the system, that is, measurements
 * of the angular velocity is available. This system benefits from a special integration method
 * called the "momentum-conserving trapezoidal method" (TRAPM), which is an invariant variational
 * method that guarantees conservation of angular momentum when no actuation is applied, i.e.,
 * it is an efficient and highly stable method. Also, this system operates within a first-order
 * (once-differentiable) SE(3) topology.
 */
class satellite3D_gyro_inv_dt_system : public satellite3D_inv_dt_system {
 public:
  using state_space_type = satellite3D_inv_dt_system::state_space_type;

  using covar_type = satellite3D_inv_dt_system::covar_type;
  using covar_space_type = satellite3D_inv_dt_system::covar_space_type;
  using temporal_state_space_type =
      satellite3D_inv_dt_system::temporal_state_space_type;
  using belief_space_type = satellite3D_inv_dt_system::belief_space_type;
  using temporal_belief_space_type =
      satellite3D_inv_dt_system::temporal_belief_space_type;
  using state_belief_type = satellite3D_inv_dt_system::state_belief_type;
  using input_belief_type = satellite3D_inv_dt_system::input_belief_type;
  using output_belief_type = satellite3D_inv_dt_system::output_belief_type;

  using point_type = satellite3D_inv_dt_system::point_type;
  using point_difference_type =
      satellite3D_inv_dt_system::point_difference_type;
  using point_derivative_type =
      satellite3D_inv_dt_system::point_derivative_type;

  using time_type = satellite3D_inv_dt_system::time_type;
  using time_difference_type = satellite3D_inv_dt_system::time_difference_type;

  using input_type = satellite3D_inv_dt_system::input_type;
  using output_type = satellite3D_inv_dt_system::output_type;

  using invariant_error_type = satellite3D_inv_dt_system::invariant_error_type;
  using invariant_correction_type =
      satellite3D_inv_dt_system::invariant_correction_type;
  using invariant_frame_type = satellite3D_inv_dt_system::invariant_frame_type;

  static constexpr std::size_t dimensions = 13;
  static constexpr std::size_t input_dimensions = 6;
  static constexpr std::size_t output_dimensions = 10;
  static constexpr std::size_t invariant_error_dimensions = 9;
  static constexpr std::size_t invariant_correction_dimensions = 12;

  using matrixA_type = satellite3D_inv_dt_system::matrixA_type;
  using matrixB_type = satellite3D_inv_dt_system::matrixB_type;
  using matrixC_type = satellite3D_inv_dt_system::matrixC_type;
  using matrixD_type = satellite3D_inv_dt_system::matrixD_type;

  using zero_input_trajectory =
      satellite3D_inv_dt_system::zero_input_trajectory;

  output_belief_type get_zero_output_belief(
      double aCovValue = 1.0) const override;

  /**
   * Returns the dimensions of the output of the system.
   * \return The dimensions of the output of the system.
   */
  std::size_t get_output_dimensions() const override { return 10; }

  /**
   * Returns the dimensions of the invariant errors of the system.
   * \return The dimensions of the invariant errors of the system.
   */
  std::size_t get_invariant_error_dimensions() const override { return 9; }

  /**
   * Constructor.
   * \param aName The name for this object.
   * \param aMass The mass of the satellite.
   * \param aInertiaMoment The inertia tensor of the satellite.
   * \param aDt The time-step for this discrete-time system.
   */
  explicit satellite3D_gyro_inv_dt_system(
      const std::string& aName, double aMass = 1.0,
      const mat<double, mat_structure::symmetric>& aInertiaMoment =
          (mat<double, mat_structure::symmetric>(
              mat<double, mat_structure::identity>(3))),
      double aDt = 0.001);

  satellite3D_gyro_inv_dt_system() : satellite3D_gyro_inv_dt_system("") {}

  output_type get_output(const state_space_type& /*space*/, const point_type& x,
                         const input_type& u,
                         const time_type& t = 0.0) const override;

  void get_output_function_blocks(matrixC_type& C, matrixD_type& D,
                                  const state_space_type& /*unused*/,
                                  const time_type& /*unused*/,
                                  const point_type& /*unused*/,
                                  const input_type& /*unused*/) const override;

  invariant_error_type get_invariant_error(const state_space_type& /*space*/,
                                           const point_type& x,
                                           const input_type& u,
                                           const output_type& y,
                                           const time_type& t) const override;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    satellite3D_inv_dt_system::save(
        A, satellite3D_inv_dt_system::getStaticObjectType()->TypeVersion());
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    satellite3D_inv_dt_system::load(
        A, satellite3D_inv_dt_system::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(satellite3D_gyro_inv_dt_system, 0xC2310019, 1,
                              "satellite3D_gyro_inv_dt_system",
                              satellite3D_inv_dt_system)
};

template <>
struct is_invariant_system<satellite3D_gyro_inv_dt_system> : std::true_type {};

}  // namespace ReaK::ctrl

#endif

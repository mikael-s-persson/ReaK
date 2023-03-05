/**
 * \file airship_IMU_mixins.hpp
 *
 * This library contains a number of building blocks for invariantized discrete-time state-space systems
 * to describe the dynamics of an airship (UAV). These are simplified models,
 * with forces of limited complexity applied, just free-floating dynamics with 6 dof actuation forces
 * and some simple augmented states for drag and imbalances. These systems
 * benefit from a special integration method called the "momentum-conserving trapezoidal method" (TRAPM),
 * which is an invariant variational method that guarantees conservation of angular momentum
 * when no actuation is applied, i.e., it is an efficient and highly stable method.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date May 2014
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

#ifndef REAK_AIRSHIP_IMU_MIXINS_HPP
#define REAK_AIRSHIP_IMU_MIXINS_HPP

#include "ReaK/core/base/named_object.hpp"

#include "ReaK/control/systems/airship_basic_mixins.hpp"
#include "ReaK/control/systems/state_space_system_tuple.hpp"

#include "ReaK/topologies/spaces/hyperball_topology.hpp"

namespace ReaK::ctrl {

class gyros_bias_state_model : public named_object {
 public:
  using state_space_type = pp::hyperball_topology<vect<double, 3>>;

  using point_type = pp::topology_point_type_t<state_space_type>;
  using point_difference_type =
      pp::topology_point_difference_type_t<state_space_type>;
  using point_derivative_type =
      pp::topology_point_difference_type_t<state_space_type>;

  using time_type = double;
  using time_difference_type = double;

 private:
  std::size_t state_start_index;
  std::size_t inv_corr_start_index;

 public:
  state_space_type create_state_space() const {
    return pp::hyperball_topology<vect<double, 3>>(
        "gyros_bias_param_space", vect<double, 3>(0.0, 0.0, 0.0),
        std::numeric_limits<double>::infinity());
  }

  void get_zero_state(point_type& x) const {
    x = vect<double, 3>(0.0, 0.0, 0.0);
  }

  std::size_t get_state_start_index() const { return state_start_index; }
  std::size_t get_inv_corr_start_index() const { return inv_corr_start_index; }

  gyros_bias_state_model() = default;

  void construct_all_dimensions(std::size_t& state_dim,
                                std::size_t& inv_corr_dim,
                                std::size_t& actual_dim) {
    state_start_index = state_dim;
    state_dim += 3;
    inv_corr_start_index = inv_corr_dim;
    inv_corr_dim += 3;
    RK_UNUSED(actual_dim);
  }

  template <typename FlyWeight, typename StateSpaceType, typename InputType>
  void add_to_fly_weight_params(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      const InputType& /*unused*/, time_difference_type dt, time_type t) const {
  }

  template <typename FlyWeight, typename StateSpaceType, typename InputType>
  void add_state_difference(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      pp::topology_point_difference_type_t<StateSpaceType>& dx,
      const InputType& /*unused*/, time_difference_type dt, time_type t) const {
  }

  template <typename MatrixA, typename MatrixB, typename FlyWeight,
            typename StateSpaceType, typename InputType>
  void add_state_transition_blocks(
      MatrixA& A, MatrixB& B, const FlyWeight& params,
      const StateSpaceType& space, time_type t_0, time_type t_1,
      const pp::topology_point_type_t<StateSpaceType>& p_0,
      const pp::topology_point_type_t<StateSpaceType>& p_1,
      const InputType& u_0, const InputType& u_1) const {
    const std::pair<std::size_t, std::size_t> gb_r(inv_corr_start_index,
                                                   inv_corr_start_index + 3);
    sub(A)(gb_r, gb_r) += mat_ident<double>(3);
  }

  template <typename FlyWeight, typename StateSpaceType, typename InvCorrType,
            typename InputType>
  void apply_correction_to_state(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      pp::topology_point_type_t<StateSpaceType>& x_c, const InvCorrType& c,
      const InputType& u, const time_type& t) const {
    const point_type& gb =
        params.get_state_models()
            .template get_state_for_system<gyros_bias_state_model>(x);
    point_type& gb_c =
        params.get_state_models()
            .template get_state_for_system<gyros_bias_state_model>(x_c);

    gb_c = gb + vect<double, 3>(c[inv_corr_start_index],
                                c[inv_corr_start_index + 1],
                                c[inv_corr_start_index + 2]);
  }

  template <typename FlyWeight, typename StateSpaceType, typename InputType,
            typename InvarFrameType>
  void set_invariant_frame_blocks(
      const FlyWeight& params, const StateSpaceType& space,
      InvarFrameType& invar_frame,
      const pp::topology_point_type_t<StateSpaceType>& x_0,
      const pp::topology_point_type_t<StateSpaceType>& x_1, const InputType& u,
      const time_type& t) const {
    /* identity is OK */
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(gyros_bias_state_model, 0xC231002B, 1,
                              "gyros_bias_state_model", named_object)
};

class sat_gyros_output_model : public named_object {
 public:
  using output_type = vect_n<double>;
  using invariant_error_type = vect_n<double>;

  using time_type = double;
  using time_difference_type = double;

 private:
  std::size_t start_index{0};
  std::size_t inv_start_index{0};

 public:
  sat_gyros_output_model() = default;

  void construct_output_dimensions(std::size_t& cur_dim,
                                   std::size_t& cur_inv_dim) {
    start_index = cur_dim;
    cur_dim += 3;
    inv_start_index = cur_inv_dim;
    cur_inv_dim += 3;
  }

  template <typename FlyWeight, typename StateSpaceType>
  void add_output_from_state_for_gb(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x, output_type& y,
      bool is_inv_err) const {
    if constexpr (FlyWeight::state_models_type::template has_system_v<
                      gyros_bias_state_model>) {
      // if there is an estimate of the gyro-bias, use it:
      if (!is_inv_err) {
        y[range(start_index, start_index + 3)] +=
            params.get_state_models()
                .template get_state_for_system<gyros_bias_state_model>(x);
      } else {
        y[range(inv_start_index, inv_start_index + 3)] -=
            params.get_state_models()
                .template get_state_for_system<gyros_bias_state_model>(x);
      }
    }
  }

  template <typename FlyWeight, typename StateSpaceType>
  void set_output_from_state(const FlyWeight& params,
                             const StateSpaceType& space,
                             const pp::topology_point_type_t<StateSpaceType>& x,
                             output_type& y, time_type t) const {
    const auto& x_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(x);

    y[range(start_index, start_index + 3)] = get_ang_velocity(x_se3);

    add_output_from_state_for_gb(params, space, x, y, false);
  }

  template <typename FlyWeight, typename StateSpaceType>
  void set_inv_err_from_output(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x, const output_type& y,
      invariant_error_type& e, time_type t) const {
    const auto& x_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(x);

    e[range(inv_start_index, inv_start_index + 3)] =
        vect<double, 3>(y[start_index], y[start_index + 1],
                        y[start_index + 2]) -
        get_ang_velocity(x_se3);

    add_output_from_state_for_gb(params, space, x, e, true);
  }

  template <typename MatrixC, typename FlyWeight>
  void add_output_block_for_gb(MatrixC& C, const FlyWeight& params) const {
    if constexpr (FlyWeight::state_models_type::template has_system_v<
                      gyros_bias_state_model>) {
      // if there is an estimate of the gyro-bias, use it:
      const std::size_t gb_index =
          params.get_state_models()
              .template get_system<gyros_bias_state_model>()
              .get_inv_corr_start_index();
      const std::pair<std::size_t, std::size_t> gb_r(gb_index, gb_index + 3);
      const std::pair<std::size_t, std::size_t> wm_r(inv_start_index,
                                                     inv_start_index + 3);
      sub(C)(wm_r, gb_r) += mat_ident<double>(3);
    }
  }

  template <typename MatrixC, typename MatrixD, typename FlyWeight,
            typename StateSpaceType, typename InputType>
  void add_output_function_blocks(
      MatrixC& C, MatrixD& D, const FlyWeight& params,
      const StateSpaceType& space, time_type t,
      const pp::topology_point_type_t<StateSpaceType>& p,
      const InputType& u) const {
    const std::size_t sat3d_state_index =
        params.get_state_models()
            .template get_system<satellite_state_model>()
            .get_inv_corr_start_index();
    const std::pair<std::size_t, std::size_t> w_r(sat3d_state_index + 9,
                                                  sat3d_state_index + 12);
    const std::pair<std::size_t, std::size_t> wm_r(inv_start_index,
                                                   inv_start_index + 3);

    sub(C)(wm_r, w_r) += mat_ident<double>(
        3);  // TODO Add a frame transition ? (in invariant posterior frame)

    add_output_block_for_gb(C, params);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(sat_gyros_output_model, 0xC231002C, 1,
                              "sat_gyros_output_model", named_object)
};

class accelerometer_bias_state_model : public named_object {
 public:
  using state_space_type = pp::hyperball_topology<vect<double, 3>>;

  using point_type = pp::topology_point_type_t<state_space_type>;
  using point_difference_type =
      pp::topology_point_difference_type_t<state_space_type>;
  using point_derivative_type =
      pp::topology_point_difference_type_t<state_space_type>;

  using time_type = double;
  using time_difference_type = double;

 private:
  std::size_t state_start_index;
  std::size_t inv_corr_start_index;

 public:
  state_space_type create_state_space() const {
    return pp::hyperball_topology<vect<double, 3>>(
        "accel_bias_param_space", vect<double, 3>(0.0, 0.0, 0.0),
        std::numeric_limits<double>::infinity());
  }

  void get_zero_state(point_type& x) const {
    x = vect<double, 3>(0.0, 0.0, 0.0);
  }

  std::size_t get_state_start_index() const { return state_start_index; }
  std::size_t get_inv_corr_start_index() const { return inv_corr_start_index; }

  accelerometer_bias_state_model() = default;

  void construct_all_dimensions(std::size_t& state_dim,
                                std::size_t& inv_corr_dim,
                                std::size_t& actual_dim) {
    state_start_index = state_dim;
    state_dim += 3;
    inv_corr_start_index = inv_corr_dim;
    inv_corr_dim += 3;
    RK_UNUSED(actual_dim);
  }

  template <typename FlyWeight, typename StateSpaceType, typename InputType>
  void add_to_fly_weight_params(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      const InputType& /*unused*/, time_difference_type dt, time_type t) const {
  }

  template <typename FlyWeight, typename StateSpaceType, typename InputType>
  void add_state_difference(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      pp::topology_point_difference_type_t<StateSpaceType>& dx,
      const InputType& /*unused*/, time_difference_type dt, time_type t) const {
  }

  template <typename MatrixA, typename MatrixB, typename FlyWeight,
            typename StateSpaceType, typename InputType>
  void add_state_transition_blocks(
      MatrixA& A, MatrixB& B, const FlyWeight& params,
      const StateSpaceType& space, time_type t_0, time_type t_1,
      const pp::topology_point_type_t<StateSpaceType>& p_0,
      const pp::topology_point_type_t<StateSpaceType>& p_1,
      const InputType& u_0, const InputType& u_1) const {
    const std::pair<std::size_t, std::size_t> ab_r(inv_corr_start_index,
                                                   inv_corr_start_index + 3);
    sub(A)(ab_r, ab_r) += mat_ident<double>(3);
  }

  template <typename FlyWeight, typename StateSpaceType, typename InvCorrType,
            typename InputType>
  void apply_correction_to_state(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      pp::topology_point_type_t<StateSpaceType>& x_c, const InvCorrType& c,
      const InputType& u, const time_type& t) const {
    const point_type& ab =
        params.get_state_models()
            .template get_state_for_system<accelerometer_bias_state_model>(x);
    point_type& ab_c =
        params.get_state_models()
            .template get_state_for_system<accelerometer_bias_state_model>(x_c);

    ab_c = ab + vect<double, 3>(c[inv_corr_start_index],
                                c[inv_corr_start_index + 1],
                                c[inv_corr_start_index + 2]);
  }

  template <typename FlyWeight, typename StateSpaceType, typename InputType,
            typename InvarFrameType>
  void set_invariant_frame_blocks(
      const FlyWeight& params, const StateSpaceType& space,
      InvarFrameType& invar_frame,
      const pp::topology_point_type_t<StateSpaceType>& x_0,
      const pp::topology_point_type_t<StateSpaceType>& x_1, const InputType& u,
      const time_type& t) const {
    /* identity is OK */
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(accelerometer_bias_state_model, 0xC231002D, 1,
                              "accelerometer_bias_state_model", named_object)
};

class sat_accelerometer_output_model : public named_object {
 public:
  using output_type = vect_n<double>;
  using invariant_error_type = vect_n<double>;

  using time_type = double;
  using time_difference_type = double;

 private:
  std::size_t start_index{0};
  std::size_t inv_start_index{0};

 public:
  sat_accelerometer_output_model() = default;
  ;

  void construct_output_dimensions(std::size_t& cur_dim,
                                   std::size_t& cur_inv_dim) {
    start_index = cur_dim;
    cur_dim += 3;
    inv_start_index = cur_inv_dim;
    cur_inv_dim += 3;
  }

  template <typename FlyWeight, typename StateSpaceType>
  void add_output_from_state_for_ab(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x, output_type& y,
      bool is_inv_err) const {
    if constexpr (FlyWeight::state_models_type::template has_system_v<
                      accelerometer_bias_state_model>) {
      // if there is an estimate of the gyro-bias, use it:
      if (!is_inv_err) {
        y[range(start_index, start_index + 3)] +=
            params.get_state_models()
                .template get_state_for_system<accelerometer_bias_state_model>(
                    x);
      } else {
        y[range(inv_start_index, inv_start_index + 3)] -=
            params.get_state_models()
                .template get_state_for_system<accelerometer_bias_state_model>(
                    x);
      }
    }
  }

  template <typename FlyWeight, typename StateSpaceType>
  void set_output_from_state(const FlyWeight& params,
                             const StateSpaceType& space,
                             const pp::topology_point_type_t<StateSpaceType>& x,
                             output_type& y, time_type t) const {
    const auto& x_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(x);
    auto& sys_params = params.get_system_parameters();

    y[range(start_index, start_index + 3)] =
        invert(get_quaternion(x_se3).as_rotation()) *
        sys_params.gravity_acc_vect;

    // TODO maybe add the acceleration (change of velocity), but how?

    add_output_from_state_for_ab(params, space, x, y, false);
  }

  template <typename FlyWeight, typename StateSpaceType>
  void set_inv_err_from_output(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x, const output_type& y,
      invariant_error_type& e, time_type t) const {
    const auto& x_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(x);
    auto& sys_params = params.get_system_parameters();

    e[range(inv_start_index, inv_start_index + 3)] =
        vect<double, 3>(y[start_index], y[start_index + 1],
                        y[start_index + 2]) -
        invert(get_quaternion(x_se3).as_rotation()) *
            sys_params.gravity_acc_vect;

    // TODO maybe add the acceleration (change of velocity), but how?

    add_output_from_state_for_ab(params, space, x, e, true);
  }

  template <typename MatrixC, typename FlyWeight>
  void add_output_block_for_ab(MatrixC& C, const FlyWeight& params) const {
    if constexpr (FlyWeight::state_models_type::template has_system_v<
                      accelerometer_bias_state_model>) {
      // if there is an estimate of the accel-bias, use it:
      const std::size_t accel_bias_index =
          params.get_state_models()
              .template get_system<accelerometer_bias_state_model>()
              .get_inv_corr_start_index();
      const std::pair<std::size_t, std::size_t> ab_r(accel_bias_index,
                                                     accel_bias_index + 3);
      const std::pair<std::size_t, std::size_t> am_r(inv_start_index,
                                                     inv_start_index + 3);
      sub(C)(am_r, ab_r) += mat_ident<double>(3);
    }
  }

  template <typename MatrixC, typename MatrixD, typename FlyWeight,
            typename StateSpaceType, typename InputType>
  void add_output_function_blocks(
      MatrixC& C, MatrixD& D, const FlyWeight& params,
      const StateSpaceType& space, time_type t,
      const pp::topology_point_type_t<StateSpaceType>& p,
      const InputType& u) const {
    const auto& x_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(p);
    auto& sys_params = params.get_system_parameters();

    const std::size_t sat3d_state_index =
        params.get_state_models()
            .template get_system<satellite_state_model>()
            .get_inv_corr_start_index();
    const std::pair<std::size_t, std::size_t> q_r(sat3d_state_index + 6,
                                                  sat3d_state_index + 9);
    const std::pair<std::size_t, std::size_t> am_r(inv_start_index,
                                                   inv_start_index + 3);

    vect<double, 3> local_g = invert(get_quaternion(x_se3).as_rotation()) *
                              sys_params.gravity_acc_vect;

    sub(C)(am_r, q_r) += mat<double, mat_structure::skew_symmetric>(-local_g);

    add_output_block_for_ab(C, params);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(sat_accelerometer_output_model, 0xC231002E, 1,
                              "sat_accelerometer_output_model", named_object)
};

class magnetometer_bias_state_model : public named_object {
 public:
  using state_space_type = pp::hyperball_topology<vect<double, 3>>;

  using point_type = pp::topology_point_type_t<state_space_type>;
  using point_difference_type =
      pp::topology_point_difference_type_t<state_space_type>;
  using point_derivative_type =
      pp::topology_point_difference_type_t<state_space_type>;

  using time_type = double;
  using time_difference_type = double;

 private:
  std::size_t state_start_index;
  std::size_t inv_corr_start_index;

 public:
  state_space_type create_state_space() const {
    return pp::hyperball_topology<vect<double, 3>>(
        "magnetometer_bias_param_space", vect<double, 3>(0.0, 0.0, 0.0),
        std::numeric_limits<double>::infinity());
  }

  void get_zero_state(point_type& x) const {
    x = vect<double, 3>(0.0, 0.0, 0.0);
  }

  std::size_t get_state_start_index() const { return state_start_index; }
  std::size_t get_inv_corr_start_index() const { return inv_corr_start_index; }

  magnetometer_bias_state_model() = default;

  void construct_all_dimensions(std::size_t& state_dim,
                                std::size_t& inv_corr_dim,
                                std::size_t& actual_dim) {
    state_start_index = state_dim;
    state_dim += 3;
    inv_corr_start_index = inv_corr_dim;
    inv_corr_dim += 3;
    RK_UNUSED(actual_dim);
  }

  template <typename FlyWeight, typename StateSpaceType, typename InputType>
  void add_to_fly_weight_params(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      const InputType& /*unused*/, time_difference_type dt, time_type t) const {
  }

  template <typename FlyWeight, typename StateSpaceType, typename InputType>
  void add_state_difference(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      pp::topology_point_difference_type_t<StateSpaceType>& dx,
      const InputType& /*unused*/, time_difference_type dt, time_type t) const {
  }

  template <typename MatrixA, typename MatrixB, typename FlyWeight,
            typename StateSpaceType, typename InputType>
  void add_state_transition_blocks(
      MatrixA& A, MatrixB& B, const FlyWeight& params,
      const StateSpaceType& space, time_type t_0, time_type t_1,
      const pp::topology_point_type_t<StateSpaceType>& p_0,
      const pp::topology_point_type_t<StateSpaceType>& p_1,
      const InputType& u_0, const InputType& u_1) const {
    const std::pair<std::size_t, std::size_t> mb_r(inv_corr_start_index,
                                                   inv_corr_start_index + 3);
    sub(A)(mb_r, mb_r) += mat_ident<double>(3);
  }

  template <typename FlyWeight, typename StateSpaceType, typename InvCorrType,
            typename InputType>
  void apply_correction_to_state(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      pp::topology_point_type_t<StateSpaceType>& x_c, const InvCorrType& c,
      const InputType& u, const time_type& t) const {
    const point_type& mb =
        params.get_state_models()
            .template get_state_for_system<magnetometer_bias_state_model>(x);
    point_type& mb_c =
        params.get_state_models()
            .template get_state_for_system<magnetometer_bias_state_model>(x_c);

    mb_c = mb + vect<double, 3>(c[inv_corr_start_index],
                                c[inv_corr_start_index + 1],
                                c[inv_corr_start_index + 2]);
  }

  template <typename FlyWeight, typename StateSpaceType, typename InputType,
            typename InvarFrameType>
  void set_invariant_frame_blocks(
      const FlyWeight& params, const StateSpaceType& space,
      InvarFrameType& invar_frame,
      const pp::topology_point_type_t<StateSpaceType>& x_0,
      const pp::topology_point_type_t<StateSpaceType>& x_1, const InputType& u,
      const time_type& t) const {
    /* identity is OK */
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(magnetometer_bias_state_model, 0xC231002F, 1,
                              "magnetometer_bias_state_model", named_object)
};

class sat_magnetometer_output_model : public named_object {
 public:
  using output_type = vect_n<double>;
  using invariant_error_type = vect_n<double>;

  using time_type = double;
  using time_difference_type = double;

 private:
  std::size_t start_index{0};
  std::size_t inv_start_index{0};

 public:
  sat_magnetometer_output_model() = default;

  void construct_output_dimensions(std::size_t& cur_dim,
                                   std::size_t& cur_inv_dim) {
    start_index = cur_dim;
    cur_dim += 3;
    inv_start_index = cur_inv_dim;
    cur_inv_dim += 3;
  }

  template <typename FlyWeight, typename StateSpaceType>
  void add_output_from_state_for_mb(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x, output_type& y,
      bool is_inv_err) const {
    if constexpr (FlyWeight::state_models_type::template has_system_v<
                      magnetometer_bias_state_model>) {
      // if there is an estimate of the gyro-bias, use it:
      if (!is_inv_err) {
        y[range(start_index, start_index + 3)] +=
            params.get_state_models()
                .template get_state_for_system<magnetometer_bias_state_model>(
                    x);
      } else {
        y[range(inv_start_index, inv_start_index + 3)] -=
            params.get_state_models()
                .template get_state_for_system<magnetometer_bias_state_model>(
                    x);
      }
    }
  }

  template <typename FlyWeight, typename StateSpaceType>
  void set_output_from_state(const FlyWeight& params,
                             const StateSpaceType& space,
                             const pp::topology_point_type_t<StateSpaceType>& x,
                             output_type& y, time_type t) const {
    const auto& x_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(x);
    auto& sys_params = params.get_system_parameters();

    y[range(start_index, start_index + 3)] =
        invert(get_quaternion(x_se3).as_rotation()) *
        sys_params.magnetic_field_vect;

    add_output_from_state_for_mb(params, space, x, y, false);
  }

  template <typename FlyWeight, typename StateSpaceType>
  void set_inv_err_from_output(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x, const output_type& y,
      invariant_error_type& e, time_type t) const {
    const auto& x_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(x);
    auto& sys_params = params.get_system_parameters();

    e[range(inv_start_index, inv_start_index + 3)] =
        vect<double, 3>(y[start_index], y[start_index + 1],
                        y[start_index + 2]) -
        invert(get_quaternion(x_se3).as_rotation()) *
            sys_params.magnetic_field_vect;

    add_output_from_state_for_mb(params, space, x, e, true);
  }

  template <typename MatrixC, typename FlyWeight>
  void add_output_block_for_mb(MatrixC& C, const FlyWeight& params) const {
    if constexpr (FlyWeight::state_models_type::template has_system_v<
                      magnetometer_bias_state_model>) {
      // if there is a mag-bias estimate, use it:
      const std::size_t mag_bias_index =
          params.get_state_models()
              .template get_system<magnetometer_bias_state_model>()
              .get_inv_corr_start_index();
      const std::pair<std::size_t, std::size_t> mb_r(mag_bias_index,
                                                     mag_bias_index + 3);
      const std::pair<std::size_t, std::size_t> mm_r(inv_start_index,
                                                     inv_start_index + 3);
      sub(C)(mm_r, mb_r) += mat_ident<double>(3);
    }
  }

  template <typename MatrixC, typename MatrixD, typename FlyWeight,
            typename StateSpaceType, typename InputType>
  void add_output_function_blocks(
      MatrixC& C, MatrixD& D, const FlyWeight& params,
      const StateSpaceType& space, time_type t,
      const pp::topology_point_type_t<StateSpaceType>& p,
      const InputType& u) const {
    const auto& x_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(p);
    auto& sys_params = params.get_system_parameters();

    const std::size_t sat3d_state_index =
        params.get_state_models()
            .template get_system<satellite_state_model>()
            .get_inv_corr_start_index();
    const std::pair<std::size_t, std::size_t> q_r(sat3d_state_index + 6,
                                                  sat3d_state_index + 9);
    const std::pair<std::size_t, std::size_t> mm_r(inv_start_index,
                                                   inv_start_index + 3);

    vect<double, 3> local_m = invert(get_quaternion(x_se3).as_rotation()) *
                              sys_params.magnetic_field_vect;

    sub(C)(mm_r, q_r) += mat<double, mat_structure::skew_symmetric>(-local_m);

    add_output_block_for_mb(C, params);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(sat_magnetometer_output_model, 0xC2310030, 1,
                              "sat_magnetometer_output_model", named_object)
};

}  // namespace ReaK::ctrl

#endif

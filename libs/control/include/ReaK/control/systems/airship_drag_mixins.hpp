/**
 * \file airship_drag_mixins.hpp
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

#ifndef REAK_AIRSHIP_DRAG_MIXINS_HPP
#define REAK_AIRSHIP_DRAG_MIXINS_HPP

#include <ReaK/core/base/named_object.hpp>

#include "airship_basic_mixins.hpp"
#include "state_space_system_tuple.hpp"

#include <ReaK/topologies/spaces/line_topology.hpp>

namespace ReaK::ctrl {

class linear_drag_state_model : public named_object {
 public:
  using state_space_type = pp::line_segment_topology<double>;

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
    return pp::line_segment_topology<double>(
        "linear_drag_param_space", 0.0,
        std::numeric_limits<double>::infinity());
  }

  void get_zero_state(point_type& x) const { x = 0.0; }

  std::size_t get_state_start_index() const { return state_start_index; }
  std::size_t get_inv_corr_start_index() const { return inv_corr_start_index; }

  linear_drag_state_model() = default;

  void construct_all_dimensions(std::size_t& state_dim,
                                std::size_t& inv_corr_dim,
                                std::size_t& actual_dim) {
    state_start_index = state_dim;
    state_dim += 1;
    inv_corr_start_index = inv_corr_dim;
    inv_corr_dim += 1;
    RK_UNUSED(actual_dim);
  }

  template <typename FlyWeight, typename StateSpaceType, typename InputType>
  void add_to_fly_weight_params(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      const InputType& /*unused*/, time_difference_type dt,
      time_type t) const {};

  template <typename FlyWeight, typename StateSpaceType, typename InputType>
  void add_state_difference(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      pp::topology_point_difference_type_t<StateSpaceType>& dx,
      const InputType& /*unused*/, time_difference_type dt, time_type t) const {
    const auto& x_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(x);
    auto& dx_se3 =
        params.get_state_models()
            .template get_state_diff_for_system<satellite_state_model>(dx);
    auto& sys_params = params.get_system_parameters();

    const point_type d_f =
        params.get_state_models()
            .template get_state_for_system<linear_drag_state_model>(x);

    vect<double, 3> fd =
        (-d_f * norm_2(get_velocity(x_se3))) * get_velocity(x_se3);

    // velocity:
    fd *= (dt / (sys_params.effective_mass + sys_params.added_mass));
    get<1>(get<0>(dx_se3)) += fd;
    // position:
    get<0>(get<0>(dx_se3)) += (0.5 * dt) * fd;
  }

  template <typename MatrixA, typename MatrixB, typename FlyWeight,
            typename StateSpaceType, typename InputType>
  void add_state_transition_blocks(
      MatrixA& A, MatrixB& B, const FlyWeight& params,
      const StateSpaceType& space, time_type t_0, time_type t_1,
      const pp::topology_point_type_t<StateSpaceType>& p_0,
      const pp::topology_point_type_t<StateSpaceType>& p_1,
      const InputType& u_0, const InputType& u_1) const {
    const auto& x0_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(p_0);
    const auto& x1_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(p_1);
    auto& sys_params = params.get_system_parameters();

    const point_type d_f =
        params.get_state_models()
            .template get_state_for_system<linear_drag_state_model>(p_0);
    const std::size_t sat3d_state_index =
        params.get_state_models()
            .template get_system<satellite_state_model>()
            .get_inv_corr_start_index();

    const double dt = t_1 - t_0;

    const std::pair<std::size_t, std::size_t> p_r(sat3d_state_index,
                                                  sat3d_state_index + 3);
    const std::pair<std::size_t, std::size_t> v_r(sat3d_state_index + 3,
                                                  sat3d_state_index + 6);

    const std::size_t fd_r = inv_corr_start_index;

    vect<double, 3> v_avg = 0.5 * (get_velocity(x0_se3) + get_velocity(x1_se3));

    mat<double, mat_structure::symmetric> v_avg_outer(
        v_avg[0] * v_avg[0], v_avg[0] * v_avg[1], v_avg[0] * v_avg[2],
        v_avg[1] * v_avg[1], v_avg[1] * v_avg[2], v_avg[2] * v_avg[2]);
    double v_avg_mag = norm_2(v_avg);
    // v-v block:
    if (v_avg_mag > 1e-4) {
      mat<double, mat_structure::square> delv(
          (d_f * dt / (sys_params.effective_mass + sys_params.added_mass)) *
          ((1.0 / v_avg_mag) * v_avg_outer + v_avg_mag * mat_ident<double>(3)));
      sub(A)(v_r, v_r) -= delv;
      sub(A)(p_r, v_r) -= (0.5 * dt) * delv;
    }

    slice(A)(v_r, fd_r) -=
        (dt * v_avg_mag / (sys_params.effective_mass + sys_params.added_mass)) *
        v_avg;
    slice(A)(p_r, fd_r) -=
        (0.5 * dt * dt * v_avg_mag /
         (sys_params.effective_mass + sys_params.added_mass)) *
        v_avg;

    A(fd_r, fd_r) += 1.0;
  }

  template <typename FlyWeight, typename StateSpaceType, typename InvCorrType,
            typename InputType>
  void apply_correction_to_state(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      pp::topology_point_type_t<StateSpaceType>& x_c, const InvCorrType& c,
      const InputType& u, const time_type& t) const {
    const point_type& d_f =
        params.get_state_models()
            .template get_state_for_system<linear_drag_state_model>(x);
    point_type& d_f_c =
        params.get_state_models()
            .template get_state_for_system<linear_drag_state_model>(x_c);

    d_f_c = d_f + c[inv_corr_start_index];
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

  RK_RTTI_MAKE_CONCRETE_1BASE(linear_drag_state_model, 0xC2310024, 1,
                              "linear_drag_state_model", named_object)
};

class torsional_drag_state_model : public named_object {
 public:
  using state_space_type = pp::line_segment_topology<double>;

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
    return pp::line_segment_topology<double>(
        "torsional_drag_param_space", 0.0,
        std::numeric_limits<double>::infinity());
  }

  void get_zero_state(point_type& x) const { x = 0.0; }

  std::size_t get_state_start_index() const { return state_start_index; }
  std::size_t get_inv_corr_start_index() const { return inv_corr_start_index; }

  torsional_drag_state_model() = default;

  void construct_all_dimensions(std::size_t& state_dim,
                                std::size_t& inv_corr_dim,
                                std::size_t& actual_dim) {
    state_start_index = state_dim;
    state_dim += 1;
    inv_corr_start_index = inv_corr_dim;
    inv_corr_dim += 1;
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
    const auto& x_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(x);
    auto& dx_se3 =
        params.get_state_models()
            .template get_state_diff_for_system<satellite_state_model>(dx);
    auto& sys_params = params.get_system_parameters();

    const point_type d_t =
        params.get_state_models()
            .template get_state_for_system<torsional_drag_state_model>(x);

    vect<double, 3> td =
        (-d_t * norm_2(get_ang_velocity(x_se3))) * get_ang_velocity(x_se3);

    // ang-velocity:
    td = dt * sys_params.effective_J_inv * td;
    get<1>(get<1>(dx_se3)) += td;
    // quat-diff:
    get<0>(get<1>(dx_se3)) += (0.5 * dt) * td;
  }

  template <typename MatrixA, typename MatrixB, typename FlyWeight,
            typename StateSpaceType, typename InputType>
  void add_state_transition_blocks(
      MatrixA& A, MatrixB& B, const FlyWeight& params,
      const StateSpaceType& space, time_type t_0, time_type t_1,
      const pp::topology_point_type_t<StateSpaceType>& p_0,
      const pp::topology_point_type_t<StateSpaceType>& p_1,
      const InputType& u_0, const InputType& u_1) const {
    const auto& x0_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(p_0);
    const auto& x1_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(p_1);
    auto& sys_params = params.get_system_parameters();

    const point_type d_t =
        params.get_state_models()
            .template get_state_for_system<torsional_drag_state_model>(p_0);
    const std::size_t sat3d_state_index =
        params.get_state_models()
            .template get_system<satellite_state_model>()
            .get_inv_corr_start_index();

    const double dt = t_1 - t_0;

    const std::pair q_r(sat3d_state_index + 6, sat3d_state_index + 9);
    const std::pair w_r(sat3d_state_index + 9, sat3d_state_index + 12);

    const std::size_t td_r = inv_corr_start_index;

    mat<double, mat_structure::square> R_0_1(
        (invert(get_quaternion(x1_se3).as_rotation()) *
         get_quaternion(x0_se3).as_rotation())
            .getMat());

    // get the average ang-velocity in the destination frame.
    vect<double, 3> w_avg =
        0.5 * (R_0_1 * get_ang_velocity(x0_se3) + get_ang_velocity(x1_se3));

    mat<double, mat_structure::symmetric> w_avg_outer(
        w_avg[0] * w_avg[0], w_avg[0] * w_avg[1], w_avg[0] * w_avg[2],
        w_avg[1] * w_avg[1], w_avg[1] * w_avg[2], w_avg[2] * w_avg[2]);
    double w_avg_mag = norm_2(w_avg);
    // v-v block:
    if (w_avg_mag > 1e-4) {
      mat<double, mat_structure::square> delw(
          (d_t * dt) * sys_params.effective_J_inv *
          ((1.0 / w_avg_mag) * w_avg_outer + w_avg_mag * mat_ident<double>(3)));
      sub(A)(w_r, w_r) -= delw;
      sub(A)(q_r, w_r) -= (0.5 * dt) * delw;
    }
    w_avg = sys_params.effective_J_inv * w_avg;
    slice(A)(w_r, td_r) -= (dt * w_avg_mag) * w_avg;
    slice(A)(q_r, td_r) -= (0.5 * dt * dt * w_avg_mag) * w_avg;

    A(td_r, td_r) += 1.0;
  }

  template <typename FlyWeight, typename StateSpaceType, typename InvCorrType,
            typename InputType>
  void apply_correction_to_state(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      pp::topology_point_type_t<StateSpaceType>& x_c, const InvCorrType& c,
      const InputType& u, const time_type& t) const {
    const point_type& d_t =
        params.get_state_models()
            .template get_state_for_system<torsional_drag_state_model>(x);
    point_type& d_t_c =
        params.get_state_models()
            .template get_state_for_system<torsional_drag_state_model>(x_c);

    d_t_c = d_t + c[inv_corr_start_index];
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

  RK_RTTI_MAKE_CONCRETE_1BASE(torsional_drag_state_model, 0xC2310025, 1,
                              "torsional_drag_state_model", named_object)
};

}  // namespace ReaK::ctrl

#endif

/**
 * \file airship_thruster_mixins.h
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

#ifndef REAK_CONTROL_SYSTEMS_AIRSHIP_THRUSTER_MIXINS_H_
#define REAK_CONTROL_SYSTEMS_AIRSHIP_THRUSTER_MIXINS_H_

#include "ReaK/core/base/named_object.h"

#include "ReaK/control/systems/airship_basic_mixins.h"
#include "ReaK/control/systems/state_space_system_tuple.h"

namespace ReaK::ctrl {

class airship3D_6dof_thrusters : public named_object {
 public:
  using input_type = vect_n<double>;

  using time_type = double;
  using time_difference_type = double;

 private:
  std::size_t start_index{0};

 public:
  airship3D_6dof_thrusters() = default;

  void construct_input_dimensions(std::size_t& cur_dim) {
    start_index = cur_dim;
    cur_dim += 6;
  }

  template <typename FlyWeight, typename StateSpaceType>
  void add_state_difference(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      pp::topology_point_difference_type_t<StateSpaceType>& dx,
      const input_type& u, time_difference_type dt, time_type t) const {
    using ReaK::get;  // for ADL

    const auto& x_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(x);
    auto& dx_se3 =
        params.get_state_models()
            .template get_state_diff_for_system<satellite_state_model>(dx);
    auto& sys_params = params.get_system_parameters();

    vect<double, 3> f(u[start_index], u[start_index + 1], u[start_index + 2]);
    vect<double, 3> tau(u[start_index + 3], u[start_index + 4],
                        u[start_index + 5]);

    // velocity:
    f = get_quaternion(x_se3).as_rotation() * f;
    f *= (dt / (sys_params.effective_mass + sys_params.added_mass));
    get<1>(get<0>(dx_se3)) += f;
    // position:
    get<0>(get<0>(dx_se3)) += (0.5 * dt) * f;

    // ang-velocity
    tau = sys_params.effective_J_inv * tau;
    tau *= dt;
    get<1>(get<1>(dx_se3)) += tau;
    // quaternion-diff (Lie alg.):
    get<0>(get<1>(dx_se3)) += (0.5 * dt) * tau;
  }

  template <typename MatrixA, typename MatrixB, typename FlyWeight,
            typename StateSpaceType>
  void add_state_transition_blocks(
      MatrixA& A, MatrixB& B, const FlyWeight& params,
      const StateSpaceType& space, time_type t_0, time_type t_1,
      const pp::topology_point_type_t<StateSpaceType>& p_0,
      const pp::topology_point_type_t<StateSpaceType>& p_1,
      const input_type& u_0, const input_type& u_1) const {
    const auto& x0_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(p_0);
    const auto& x1_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(p_1);
    auto& sys_params = params.get_system_parameters();

    const std::size_t sat3d_state_index =
        params.get_state_models()
            .template get_system<satellite_state_model>()
            .get_inv_corr_start_index();

    const double dt = t_1 - t_0;

    const std::pair p_r(sat3d_state_index, sat3d_state_index + 3);
    const std::pair v_r(sat3d_state_index + 3, sat3d_state_index + 6);
    const std::pair q_r(sat3d_state_index + 6, sat3d_state_index + 9);
    const std::pair w_r(sat3d_state_index + 9, sat3d_state_index + 12);

    const std::pair f_r(start_index, start_index + 3);
    const std::pair t_r(start_index + 3, start_index + 6);

    // (p,v)-f block:
    mat<double, mat_structure::square> R_0(
        get_quaternion(x0_se3).as_rotation().getMat());
    mat<double, mat_structure::square> delv(
        (dt / (sys_params.effective_mass + sys_params.added_mass)) * R_0);
    sub(B)(v_r, f_r) += delv;
    delv *= 0.5 * dt;
    sub(B)(p_r, f_r) += delv;

    // (q,w)-t block:
    mat<double, mat_structure::square> R_0_1(
        (invert(get_quaternion(x1_se3).as_rotation()) *
         get_quaternion(x0_se3).as_rotation())
            .getMat());
    mat<double, mat_structure::square> delw(
        dt * (R_0_1 * sys_params.effective_J_inv));
    sub(B)(w_r, t_r) += delw;
    delw *= 0.5 * dt;
    sub(B)(q_r, t_r) += delw;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::get_static_object_type()->version());
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::get_static_object_type()->version());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(airship3D_6dof_thrusters, 0xC2310027, 1,
                              "airship3D_6dof_thrusters", named_object)
};

class tryphon_n_thrusters : public named_object {
 public:
  using input_type = vect_n<double>;

  using time_type = double;
  using time_difference_type = double;

 private:
  std::size_t start_index{0};

 public:
  std::vector<vect<double, 3>> thruster_pos;
  std::vector<vect<double, 3>> thruster_dir;

  explicit tryphon_n_thrusters(std::size_t N)
      : thruster_pos(N), thruster_dir(N) {}

  tryphon_n_thrusters() : tryphon_n_thrusters(8) {}

  void construct_input_dimensions(std::size_t& cur_dim) {
    start_index = cur_dim;
    cur_dim += thruster_pos.size();
  }

  template <typename FlyWeight, typename StateSpaceType>
  void add_state_difference(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      pp::topology_point_difference_type_t<StateSpaceType>& dx,
      const input_type& u, time_difference_type dt, time_type t) const {
    using ReaK::get;  // for ADL

    const auto& x_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(x);
    auto& dx_se3 =
        params.get_state_models()
            .template get_state_diff_for_system<satellite_state_model>(dx);
    typename FlyWeight::system_param_type& sys_params =
        params.get_system_parameters();

    vect<double, 3> f;
    vect<double, 3> tau;
    for (std::size_t i = 0; i < 8; ++i) {
      f += u[start_index + i] * thruster_dir[i];
      tau += u[start_index + i] * (thruster_pos[i] % thruster_dir[i]);
    }

    // velocity:
    f = get_quaternion(x_se3).as_rotation() * f;
    f *= (dt / (sys_params.effective_mass + sys_params.added_mass));
    get<1>(get<0>(dx_se3)) += f;
    // position:
    get<0>(get<0>(dx_se3)) += (0.5 * dt) * f;

    // ang-velocity
    tau = sys_params.effective_J_inv * tau;
    tau *= dt;
    get<1>(get<1>(dx_se3)) += tau;
    // quaternion-diff (Lie alg.):
    get<0>(get<1>(dx_se3)) += (0.5 * dt) * tau;
  }

  template <typename MatrixA, typename MatrixB, typename FlyWeight,
            typename StateSpaceType>
  void add_state_transition_blocks(
      MatrixA& A, MatrixB& B, const FlyWeight& params,
      const StateSpaceType& space, time_type t_0, time_type t_1,
      const pp::topology_point_type_t<StateSpaceType>& p_0,
      const pp::topology_point_type_t<StateSpaceType>& p_1,
      const input_type& u_0, const input_type& u_1) const {
    const auto& x0_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(p_0);
    const auto& x1_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(p_1);
    auto& sys_params = params.get_system_parameters();

    const std::size_t sat3d_state_index =
        params.get_state_models()
            .template get_system<satellite_state_model>()
            .get_inv_corr_start_index();

    const double dt = t_1 - t_0;

    const std::pair p_r(sat3d_state_index, sat3d_state_index + 3);
    const std::pair v_r(sat3d_state_index + 3, sat3d_state_index + 6);
    const std::pair q_r(sat3d_state_index + 6, sat3d_state_index + 9);
    const std::pair w_r(sat3d_state_index + 9, sat3d_state_index + 12);

    mat<double, mat_structure::square> R_0(
        get_quaternion(x0_se3).as_rotation().getMat());
    mat<double, mat_structure::square> R_0_1(
        (invert(get_quaternion(x1_se3).as_rotation()) *
         get_quaternion(x0_se3).as_rotation())
            .getMat());

    for (std::size_t i = 0; i < 8; ++i) {
      // (p,v)-f block:
      vect<double, 3> f_jac =
          (dt / (sys_params.effective_mass + sys_params.added_mass)) *
          (R_0 * thruster_dir[i]);
      slice(B)(v_r, start_index + i) += f_jac;
      f_jac *= 0.5 * dt;
      slice(B)(p_r, start_index + i) += f_jac;

      // (q,w)-t block:
      vect<double, 3> tau_jac =
          dt * (R_0_1 * (sys_params.effective_J_inv *
                         (thruster_pos[i] % thruster_dir[i])));
      slice(B)(w_r, start_index + i) += tau_jac;
      tau_jac *= 0.5 * dt;
      slice(B)(q_r, start_index + i) += tau_jac;
    }
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(thruster_pos) &
        RK_SERIAL_SAVE_WITH_NAME(thruster_dir);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(thruster_pos) &
        RK_SERIAL_LOAD_WITH_NAME(thruster_dir);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(tryphon_n_thrusters, 0xC2310028, 1,
                              "tryphon_n_thrusters", named_object)
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_SYSTEMS_AIRSHIP_THRUSTER_MIXINS_H_

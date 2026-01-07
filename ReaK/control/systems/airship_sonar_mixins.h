/**
 * \file airship_sonar_mixins.h
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

#ifndef REAK_CONTROL_SYSTEMS_AIRSHIP_SONAR_MIXINS_H_
#define REAK_CONTROL_SYSTEMS_AIRSHIP_SONAR_MIXINS_H_

#include <cmath>

#include "ReaK/core/base/named_object.h"

#include "ReaK/control/systems/airship_basic_mixins.h"
#include "ReaK/control/systems/state_space_system_tuple.h"

#include "ReaK/topologies/spaces/line_topology.h"

namespace ReaK::ctrl {

class room_orientation_state_model : public named_object {
 public:
  using state_space_type = pp::line_segment_topology<double>;

  using point_type = pp::topology_traits<state_space_type>::point_type;
  using point_difference_type =
      pp::topology_traits<state_space_type>::point_difference_type;
  using point_derivative_type =
      pp::topology_traits<state_space_type>::point_difference_type;

  using time_type = double;
  using time_difference_type = double;

 private:
  std::size_t state_start_index;
  std::size_t inv_corr_start_index;

 public:
  state_space_type create_state_space() const {
    return pp::line_segment_topology<double>(
        "room_angle_param_space", -std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity());
  }

  void get_zero_state(point_type& x) const { x = 0.0; }

  std::size_t get_state_start_index() const { return state_start_index; }
  std::size_t get_inv_corr_start_index() const { return inv_corr_start_index; }

  room_orientation_state_model() = default;

  void construct_all_dimensions(std::size_t& state_dim,
                                std::size_t& inv_corr_dim,
                                [[maybe_unused]] std::size_t& actual_dim) {
    state_start_index = state_dim;
    state_dim += 1;
    inv_corr_start_index = inv_corr_dim;
    inv_corr_dim += 1;
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
    const std::size_t ro_r = inv_corr_start_index;
    A(ro_r, ro_r) += 1.0;
  }

  template <typename FlyWeight, typename StateSpaceType, typename InvCorrType,
            typename InputType>
  void apply_correction_to_state(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      pp::topology_point_type_t<StateSpaceType>& x_c, const InvCorrType& c,
      const InputType& u, const time_type& t) const {
    const point_type& ro =
        params.get_state_models()
            .template get_state_for_system<room_orientation_state_model>(x);
    point_type& ro_c =
        params.get_state_models()
            .template get_state_for_system<room_orientation_state_model>(x_c);

    ro_c = ro + c[inv_corr_start_index];
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

  RK_RTTI_MAKE_CONCRETE_1BASE(room_orientation_state_model, 0xC2310031, 1,
                              "room_orientation_state_model", named_object)
};

class sonars_in_room_output_model : public named_object {
 public:
  using output_type = vect_n<double>;
  using invariant_error_type = vect_n<double>;

  using time_type = double;
  using time_difference_type = double;

 private:
  std::size_t start_index{0};
  std::size_t inv_start_index{0};

  void get_sonar_distance_to_room(const vect<double, 3>& spos_gbl,
                                  const vect<double, 3>& sdir_gbl, double& y,
                                  int& surface_id) const {

    y = std::numeric_limits<double>::infinity();
    surface_id = -1;

    // NOTE: the negative sign in the expressions below is correct.

    if (std::abs(sdir_gbl[0]) > 1e-4) {
      double tmp = -(spos_gbl[0] - lower_corner[0]) / sdir_gbl[0];
      if ((tmp > 0.0) && (tmp < y)) {
        y = tmp;
        surface_id = 0;
      }
      tmp = (upper_corner[0] - spos_gbl[0]) / sdir_gbl[0];
      if ((tmp > 0.0) && (tmp < y)) {
        y = tmp;
        surface_id = 1;
      }
    }

    if (std::abs(sdir_gbl[1]) > 1e-4) {
      double tmp = -(spos_gbl[1] - lower_corner[1]) / sdir_gbl[1];
      if ((tmp > 0.0) && (tmp < y)) {
        y = tmp;
        surface_id = 2;
      }
      tmp = (upper_corner[1] - spos_gbl[1]) / sdir_gbl[1];
      if ((tmp > 0.0) && (tmp < y)) {
        y = tmp;
        surface_id = 3;
      }
    }

    if (std::abs(sdir_gbl[2]) > 1e-4) {
      double tmp = -(spos_gbl[2] - lower_corner[2]) / sdir_gbl[2];
      if ((tmp > 0.0) && (tmp < y)) {
        y = tmp;
        surface_id = 4;
      }
      tmp = (upper_corner[2] - spos_gbl[2]) / sdir_gbl[2];
      if ((tmp > 0.0) && (tmp < y)) {
        y = tmp;
        surface_id = 5;
      }
    }
    if (y == std::numeric_limits<double>::infinity()) {
      y = 0.0;  // prevent impossible distances (outside the box) from wreaking havoc in the system.
    }
  }

 public:
  std::vector<vect<double, 3>> sonar_pos;
  std::vector<vect<double, 3>> sonar_dir;

  vect<double, 3> lower_corner;
  vect<double, 3> upper_corner;

  explicit sonars_in_room_output_model(std::size_t N)
      : sonar_pos(N),
        sonar_dir(N),
        lower_corner(0.0, 0.0, 0.0),
        upper_corner(1.0, 1.0, 1.0) {}

  sonars_in_room_output_model() : sonars_in_room_output_model(6) {}

  void construct_output_dimensions(std::size_t& cur_dim,
                                   std::size_t& cur_inv_dim) {
    start_index = cur_dim;
    cur_dim += sonar_pos.size();
    inv_start_index = cur_inv_dim;
    cur_inv_dim += sonar_pos.size();
  }

  template <typename FlyWeight, typename StateSpaceType>
  void accum_rotation_from_ro(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      mat<double, mat_structure::square>& R) const {
    if constexpr (FlyWeight::state_models_type::template has_system_v<
                      room_orientation_state_model>) {
      double ro_angle =
          params.get_state_models()
              .template get_state_for_system<room_orientation_state_model>(x);
      R = quaternion<double>::zrot(ro_angle).getQuaternion().getMat() * R;
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

    mat<double, mat_structure::square> R(
        get_quaternion(x_se3).as_rotation().getMat());

    accum_rotation_from_ro(params, space, x, R);

    int id = 0;
    for (std::size_t i = 0; i < sonar_pos.size(); ++i) {
      vect<double, 3> spos_gbl = R * sonar_pos[i] + get_position(x_se3);
      vect<double, 3> sdir_gbl = R * sonar_dir[i];
      get_sonar_distance_to_room(spos_gbl, sdir_gbl, y[start_index + i], id);
    }
  }

  template <typename FlyWeight, typename StateSpaceType>
  void set_inv_err_from_output(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x, const output_type& y,
      invariant_error_type& e, time_type t) const {
    const auto& x_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(x);

    mat<double, mat_structure::square> R(
        get_quaternion(x_se3).as_rotation().getMat());

    accum_rotation_from_ro(params, space, x, R);

    int id = 0;
    for (std::size_t i = 0; i < sonar_pos.size(); ++i) {
      vect<double, 3> spos_gbl = R * sonar_pos[i] + get_position(x_se3);
      vect<double, 3> sdir_gbl = R * sonar_dir[i];
      get_sonar_distance_to_room(spos_gbl, sdir_gbl, e[inv_start_index + i],
                                 id);
      e[inv_start_index + i] = y[start_index + i] - e[inv_start_index + i];
    }
  }

  template <typename MatrixC, typename FlyWeight>
  void accum_output_del_from_ro(MatrixC& C, const FlyWeight& params,
                                const vect<double, 3>& x_pos,
                                const mat<double, mat_structure::square>& R,
                                const vect<double, 3>& spos_gbl,
                                const vect<double, 3>& sdir_gbl, std::size_t i,
                                std::size_t coord) const {
    if constexpr (FlyWeight::state_models_type::template has_system_v<
                      room_orientation_state_model>) {
      const std::size_t ro_r =
          params.get_state_models()
              .template get_system<room_orientation_state_model>()
              .get_inv_corr_start_index();

      vect<double, 3> spos_ro_arm = vect_k % (spos_gbl - x_pos);
      C(inv_start_index + i, ro_r) -=
          (1.0 / sdir_gbl[coord]) * spos_ro_arm[coord];
    }
  }

  template <typename MatrixC, typename MatrixD, typename FlyWeight,
            typename StateSpaceType, typename InputType>
  void add_output_function_blocks(
      MatrixC& C, MatrixD& D, const FlyWeight& params,
      const StateSpaceType& space, time_type t,
      const pp::topology_point_type_t<StateSpaceType>& x,
      const InputType& u) const {
    const auto& x_se3 =
        params.get_state_models()
            .template get_state_for_system<satellite_state_model>(x);

    const std::size_t sat3d_state_index =
        params.get_state_models()
            .template get_system<satellite_state_model>()
            .get_inv_corr_start_index();
    const std::pair q_r(sat3d_state_index + 6, sat3d_state_index + 9);
    const std::pair mm_r(inv_start_index, inv_start_index + 3);

    mat<double, mat_structure::square> R(
        get_quaternion(x_se3).as_rotation().getMat());

    accum_rotation_from_ro(params, space, x, R);

    for (std::size_t i = 0; i < sonar_pos.size(); ++i) {
      double dist = NAN;
      int id = 0;
      vect<double, 3> spos_gbl = R * sonar_pos[i] + get_position(x_se3);
      vect<double, 3> sdir_gbl = R * sonar_dir[i];
      get_sonar_distance_to_room(spos_gbl, sdir_gbl, dist, id);

      mat<double, mat_structure::square> Rp(
          R * mat<double, mat_structure::skew_symmetric>(sonar_pos[i]));

      switch (id) {
        case 0: {  // lower-bound on x
          C(inv_start_index + i, sat3d_state_index) -= 1.0 / sdir_gbl[0];
          slice(C)(inv_start_index + i, q_r) +=
              (1.0 / sdir_gbl[0]) * slice(Rp)(0, range(0, 3));
          accum_output_del_from_ro(C, params, get_position(x_se3), R, spos_gbl,
                                   sdir_gbl, i, 0);
          break;
        }
        case 1: {  // upper-bound on x
          C(inv_start_index + i, sat3d_state_index) -= 1.0 / sdir_gbl[0];
          slice(C)(inv_start_index + i, q_r) +=
              (1.0 / sdir_gbl[0]) * slice(Rp)(0, range(0, 3));
          accum_output_del_from_ro(C, params, get_position(x_se3), R, spos_gbl,
                                   sdir_gbl, i, 0);
          break;
        }
        case 2: {  // lower-bound on y
          C(inv_start_index + i, sat3d_state_index + 1) -= 1.0 / sdir_gbl[1];
          slice(C)(inv_start_index + i, q_r) +=
              (1.0 / sdir_gbl[1]) * slice(Rp)(1, range(0, 3));
          accum_output_del_from_ro(C, params, get_position(x_se3), R, spos_gbl,
                                   sdir_gbl, i, 1);
          break;
        }
        case 3: {  // upper-bound on y
          C(inv_start_index + i, sat3d_state_index + 1) -= 1.0 / sdir_gbl[1];
          slice(C)(inv_start_index + i, q_r) +=
              (1.0 / sdir_gbl[1]) * slice(Rp)(1, range(0, 3));
          accum_output_del_from_ro(C, params, get_position(x_se3), R, spos_gbl,
                                   sdir_gbl, i, 1);
          break;
        }
        case 4: {  // lower-bound on z
          C(inv_start_index + i, sat3d_state_index + 2) -= 1.0 / sdir_gbl[2];
          slice(C)(inv_start_index + i, q_r) +=
              (1.0 / sdir_gbl[2]) * slice(Rp)(2, range(0, 3));
          // NOTE: no point in doing it for a z-component because it's a planar rotation.
          // accum_output_del_from_ro(C, params, get_position(x_se3), R, spos_gbl, sdir_gbl, i, 1);
          break;
        }
        case 5: {  // upper-bound on z
          C(inv_start_index + i, sat3d_state_index + 2) -= 1.0 / sdir_gbl[2];
          slice(C)(inv_start_index + i, q_r) +=
              (1.0 / sdir_gbl[2]) * slice(Rp)(2, range(0, 3));
          // NOTE: no point in doing it for a z-component because it's a planar rotation.
          // accum_output_del_from_ro(C, params, get_position(x_se3), R, spos_gbl, sdir_gbl, i, 2);
          break;
        }
        default:
          break;
      }
    }
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(sonar_pos) &
        RK_SERIAL_SAVE_WITH_NAME(sonar_dir) &
        RK_SERIAL_SAVE_WITH_NAME(lower_corner) &
        RK_SERIAL_SAVE_WITH_NAME(upper_corner);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(sonar_pos) &
        RK_SERIAL_LOAD_WITH_NAME(sonar_dir) &
        RK_SERIAL_LOAD_WITH_NAME(lower_corner) &
        RK_SERIAL_LOAD_WITH_NAME(upper_corner);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(sonars_in_room_output_model, 0xC2310032, 1,
                              "sonars_in_room_output_model", named_object)
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_SYSTEMS_AIRSHIP_SONAR_MIXINS_H_

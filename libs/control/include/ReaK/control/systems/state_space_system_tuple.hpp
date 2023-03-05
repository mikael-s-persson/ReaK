/**
 * \file state_space_system_tuple.hpp
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

#ifndef REAK_STATE_SPACE_SYSTEM_TUPLE_HPP
#define REAK_STATE_SPACE_SYSTEM_TUPLE_HPP

#include "ReaK/core/base/named_object.hpp"

#include "ReaK/control/systems/augmented_sss_concept.hpp"
#include "ReaK/control/systems/invariant_system_concept.hpp"

#include "ReaK/topologies/spaces/metric_space_tuple.hpp"
#include "ReaK/topologies/spaces/temporal_space.hpp"
#include "ReaK/topologies/spaces/time_poisson_topology.hpp"
#include "ReaK/topologies/spaces/tuple_distance_metrics.hpp"

#include "ReaK/control/estimators/covar_topology.hpp"
#include "ReaK/control/estimators/covariance_matrix.hpp"
#include "ReaK/control/estimators/gaussian_belief_space.hpp"
#include "ReaK/control/systems/sss_exceptions.hpp"
#include "ReaK/math/lin_alg/arithmetic_tuple.hpp"
#include "ReaK/math/lin_alg/mat_alg.hpp"

#include <type_traits>
#include <utility>

namespace ReaK::ctrl {

template <typename InputTuple>
class ss_system_input_tuple : public named_object {
 public:
  using self = ss_system_input_tuple<InputTuple>;

  using time_type = double;
  using time_difference_type = double;

  using input_type = vect_n<double>;
  using covar_type = covariance_matrix<vect_n<double>>;
  using input_belief_type = gaussian_belief_state<input_type, covar_type>;

  static constexpr std::size_t input_dimensions = 0;

 private:
  InputTuple data;
  std::size_t total_input_dim;

  void construct_input_dimensions() {
    total_input_dim = 0;
    tuple_for_each(data, [&](auto& data_elem) {
      data_elem.construct_input_dimensions(total_input_dim);
    });
  }

 public:
  template <typename InSystem>
  const InSystem& get_input_system() const {
    return get_by_type<InSystem>(data);
  }

  template <typename InSystem>
  InSystem& get_input_system() {
    return get_by_type<InSystem>(data);
  }

  template <typename InSystem>
  static constexpr bool has_input_system =
      arithmetic_tuple_has_type_v<InSystem, InputTuple>;

  /**
   * Returns a belief point for zero input values and a given uniform covariance value.
   * \param aCovValue A uniform covariance value to give to all the input component beliefs.
   * \return A belief point for zero input values and a given uniform covariance value.
   */
  input_belief_type get_zero_input_belief(double aCovValue = 1.0) const {
    return {
        input_type(vect_n<double>(total_input_dim, 0.0)),
        covar_type(covar_type::matrix_type(
            mat<double, mat_structure::diagonal>(total_input_dim, aCovValue)))};
  }

  /**
   * Returns the dimensions of the input of the system.
   * \return The dimensions of the input of the system.
   */
  std::size_t get_input_dimensions() const { return total_input_dim; }

  explicit ss_system_input_tuple(InputTuple aData) : data(std::move(aData)) {
    construct_input_dimensions();
  }

  template <typename... Args>
  explicit ss_system_input_tuple(Args&&... args)
      : data(std::forward<Args>(args)...) {
    construct_input_dimensions();
  }

  ss_system_input_tuple() : ss_system_input_tuple(InputTuple()) {}

  /**
   * Fills in the state-difference object with the effects of the given input on the state
   * over the given time difference.
   * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
   * \tparam StateSpaceType The type of the state-space (topology) used by the parent ss-system.
   * \param params The fly-weight parameters that describe the dynamics of the system.
   * \param space The space object in which to operate.
   * \param x The current state object for the system.
   * \param dx The state-difference accumulated (as output) for the different effects.
   * \param u The current input vector for the system.
   * \param dt The time period over which to compute the effect of the input vector.
   * \param t The current time corresponding to the current state of the system.
   */
  template <typename FlyWeight, typename StateSpaceType>
  void add_state_difference(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x,
      pp::topology_point_difference_type_t<StateSpaceType>& dx,
      const input_type& u, time_difference_type dt, time_type t = 0.0) const {
    tuple_for_each(data, [&](const auto& data_elem) {
      data_elem.add_state_difference(params, space, x, dx, u, dt, t);
    });
  }

  /**
   * Fills in the state-difference object with the effects of the given input on the state
   * over the given time difference.
   * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
   * \tparam StateSpaceType The type of the state-space (topology) used by the parent ss-system.
   * \param A Holds, as output, the state-to-state jacobian matrix of the state-transition of the system.
   * \param B Holds, as output, the input-to-state jacobian matrix of the state-transition of the system.
   * \param params The fly-weight parameters that describe the dynamics of the system.
   * \param space The space object in which to operate.
   * \param t_0 The previous time corresponding to the previous state of the system.
   * \param t_1 The next time corresponding to the next state of the system.
   * \param p_0 The previous state object for the system.
   * \param p_1 The next state object for the system.
   * \param u_0 The previous input vector for the system.
   * \param u_1 The next input vector for the system.
   */
  template <typename MatrixA, typename MatrixB, typename FlyWeight,
            typename StateSpaceType>
  void add_state_transition_blocks(
      MatrixA& A, MatrixB& B, const FlyWeight& params,
      const StateSpaceType& space, time_type t_0, time_type t_1,
      const pp::topology_point_type_t<StateSpaceType>& p_0,
      const pp::topology_point_type_t<StateSpaceType>& p_1,
      const input_type& u_0, const input_type& u_1) const {
    tuple_for_each(data, [&](const auto& data_elem) {
      data_elem.add_state_transition_blocks(A, B, params, space, t_0, t_1, p_0,
                                            p_1, u_0, u_1);
    });
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(data);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(data);
    construct_input_dimensions();
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2300018, 1, "ss_system_input_tuple",
                              named_object)
};

template <typename OutputTuple>
class ss_system_output_tuple : public named_object {
 public:
  using self = ss_system_output_tuple<OutputTuple>;

  using time_type = double;
  using time_difference_type = double;

  using output_type = vect_n<double>;
  using invariant_error_type = vect_n<double>;
  using covar_type = covariance_matrix<vect_n<double>>;
  using output_belief_type = gaussian_belief_state<output_type, covar_type>;

  static constexpr std::size_t output_dimensions = 0;
  static constexpr std::size_t invariant_error_dimensions = 0;

 private:
  OutputTuple data;
  std::size_t total_output_dim;
  std::size_t total_inv_err_dim;

  void construct_output_dimensions() {
    total_output_dim = 0;
    total_inv_err_dim = 0;
    tuple_for_each(data, [&](auto& data_elem) {
      data_elem.construct_output_dimensions(total_output_dim,
                                            total_inv_err_dim);
    });
  }

 public:
  template <typename OutSystem>
  const OutSystem& get_output_system() const {
    return get_by_type<OutSystem>(data);
  }

  template <typename OutSystem>
  OutSystem& get_output_system() {
    return get_by_type<OutSystem>(data);
  }

  template <typename OutSystem>
  static constexpr bool has_output_system =
      arithmetic_tuple_has_type_v<OutSystem, OutputTuple>;

  /**
   * Returns a belief point for zero input values and a given uniform covariance value.
   * \param aCovValue A uniform covariance value to give to all the input component beliefs.
   * \return A belief point for zero input values and a given uniform covariance value.
   */
  output_belief_type get_zero_output_belief(double aCovValue = 1.0) const {
    return {
        output_type(vect_n<double>(total_output_dim, 0.0)),
        covar_type(covar_type::matrix_type(mat<double, mat_structure::diagonal>(
            total_output_dim, aCovValue)))};
  }

  /**
   * Returns the dimensions of the output of the system.
   * \return The dimensions of the output of the system.
   */
  std::size_t get_output_dimensions() const { return total_output_dim; }

  /**
   * Returns the dimensions of the invariant output error of the system.
   * \return The dimensions of the invariant output error of the system.
   */
  std::size_t get_invariant_error_dimensions() const {
    return total_inv_err_dim;
  }

  explicit ss_system_output_tuple(OutputTuple aData) : data(std::move(aData)) {
    construct_output_dimensions();
  }

  template <typename... Args>
  explicit ss_system_output_tuple(Args&&... args)
      : data(std::forward<Args>(args)...) {
    construct_output_dimensions();
  }

  ss_system_output_tuple() : ss_system_output_tuple(OutputTuple()) {}

  /**
   * Fills in the state-difference object with the effects of the given input on the state
   * over the given time difference.
   * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
   * \tparam StateSpaceType The type of the state-space (topology) used by the parent ss-system.
   * \param params The fly-weight parameters that describe the dynamics of the system.
   * \param space The space object in which to operate.
   * \param x The current state object for the system.
   * \param y The current output vector for the system (accumulated as output).
   * \param t The current time corresponding to the current state of the system.
   */
  template <typename FlyWeight, typename StateSpaceType>
  void set_output_from_state(const FlyWeight& params,
                             const StateSpaceType& space,
                             const pp::topology_point_type_t<StateSpaceType>& x,
                             output_type& y, time_type t) const {
    tuple_for_each(data, [&](const auto& data_elem) {
      data_elem.set_output_from_state(params, space, x, y, t);
    });
  }

  /**
   * Fills in the invariant-output-error object with the error between the given output and
   * the current state.
   * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
   * \tparam StateSpaceType The type of the state-space (topology) used by the parent ss-system.
   * \param params The fly-weight parameters that describe the dynamics of the system.
   * \param space The space object in which to operate.
   * \param x The current state object for the system.
   * \param y The current output vector obtained from a measurement on the real system.
   * \param e The current invariant-output-error vector for the system (accumulated as output).
   * \param t The current time corresponding to the current state of the system.
   */
  template <typename FlyWeight, typename StateSpaceType>
  void set_inv_err_from_output(
      const FlyWeight& params, const StateSpaceType& space,
      const pp::topology_point_type_t<StateSpaceType>& x, const output_type& y,
      invariant_error_type& e, time_type t) const {
    tuple_for_each(data, [&](const auto& data_elem) {
      data_elem.set_inv_err_from_output(params, space, x, y, e, t);
    });
  }

  /**
   * Fills in the state-difference object with the effects of the given input on the state
   * over the given time difference.
   * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
   * \tparam StateSpaceType The type of the state-space (topology) used by the parent ss-system.
   * \param C Holds, as output, the state-to-output jacobian matrix of the output-function of the system.
   * \param D Holds, as output, the input-to-input jacobian matrix of the output-function of the system.
   * \param params The fly-weight parameters that describe the dynamics of the system.
   * \param space The space object in which to operate.
   * \param t The current time corresponding to the current state of the system.
   * \param p The current state object for the system.
   * \param u The current input vector for the system.
   */
  template <typename MatrixC, typename MatrixD, typename FlyWeight,
            typename StateSpaceType, typename InputType>
  void add_output_function_blocks(
      MatrixC& C, MatrixD& D, const FlyWeight& params,
      const StateSpaceType& space, time_type t,
      const pp::topology_point_type_t<StateSpaceType>& p,
      const InputType& u) const {
    tuple_for_each(data, [&](const auto& data_elem) {
      data_elem.add_output_function_blocks(C, D, params, space, t, p, u);
    });
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(data);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(data);
    construct_output_dimensions();
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2300019, 1, "ss_system_output_tuple",
                              named_object)
};

namespace detail {

template <std::size_t Size, typename SystemModelsTuple>
struct gather_system_state_space_tuple_impl {
  // static_assert(false);
};

template <std::size_t Size, typename... Systems>
struct gather_system_state_space_tuple_impl<Size, std::tuple<Systems...>> {
  using type = arithmetic_tuple<typename Systems::state_space_type...>;
};

template <std::size_t Size, typename... Systems>
struct gather_system_state_space_tuple_impl<Size,
                                            arithmetic_tuple<Systems...>> {
  using type = arithmetic_tuple<typename Systems::state_space_type...>;
};

}  // namespace detail

template <typename StateSysTuple>
class ss_system_state_tuple : public named_object {
 public:
  using self = ss_system_state_tuple<StateSysTuple>;

  using state_space_type = pp::metric_space_tuple<
      typename detail::gather_system_state_space_tuple_impl<
          arithmetic_tuple_size_v<StateSysTuple>, StateSysTuple>::type,
      pp::manhattan_tuple_distance>;

  using time_type = double;
  using time_difference_type = double;

  using point_type = typename pp::topology_traits<state_space_type>::point_type;
  using point_difference_type =
      typename pp::topology_traits<state_space_type>::point_difference_type;
  using point_derivative_type =
      typename pp::topology_traits<state_space_type>::point_difference_type;

  using invariant_correction_type = vect_n<double>;

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

  static constexpr std::size_t dimensions = 0;
  static constexpr std::size_t invariant_correction_dimensions = 0;
  static constexpr std::size_t actual_state_dimensions = 0;

 private:
  StateSysTuple systems;
  std::size_t total_state_dim;
  std::size_t total_inv_corr_dim;
  std::size_t total_actual_state_dim;

  static constexpr unsigned int system_tuple_size =
      arithmetic_tuple_size<StateSysTuple>::value;

  void construct_all_dimensions() {
    total_state_dim = 0;
    total_inv_corr_dim = 0;
    total_actual_state_dim = 0;
    tuple_for_each(systems, [&](auto& system_elem) {
      system_elem.construct_all_dimensions(total_state_dim, total_inv_corr_dim,
                                           total_actual_state_dim);
    });
  }

 public:
  state_space_type create_state_space() const {
    state_space_type space;
    tuple_for_each(space, systems,
                   [&](auto& space_elem, const auto& system_elem) {
                     space_elem = system_elem.create_state_space();
                   });
    return space;
  }

  template <typename System>
  static const typename System::point_type& get_state_for_system(
      const point_type& x) {
    using ReaK::get;
    return get<arithmetic_tuple_index_of_v<System, StateSysTuple>>(x);
  }

  template <typename System>
  static typename System::point_type& get_state_for_system(point_type& x) {
    using ReaK::get;
    return get<arithmetic_tuple_index_of_v<System, StateSysTuple>>(x);
  }

  template <typename System>
  static const typename System::point_difference_type&
  get_state_diff_for_system(const point_difference_type& x) {
    using ReaK::get;
    return get<arithmetic_tuple_index_of_v<System, StateSysTuple>>(x);
  }

  template <typename System>
  static typename System::point_difference_type& get_state_diff_for_system(
      point_difference_type& x) {
    using ReaK::get;
    return get<arithmetic_tuple_index_of_v<System, StateSysTuple>>(x);
  }

  template <typename System>
  const System& get_system() const {
    return get_by_type<System>(systems);
  }

  template <typename System>
  System& get_system() {
    return get_by_type<System>(systems);
  }

  template <typename System>
  static constexpr bool has_system_v =
      arithmetic_tuple_has_type_v<System, StateSysTuple>;

  /**
   * Returns a belief point for zero input values and a given uniform covariance value.
   * \param aCovValue A uniform covariance value to give to all the input component beliefs.
   * \return A belief point for zero input values and a given uniform covariance value.
   */
  state_belief_type get_zero_state_belief(double aCovValue = 1.0) const {
    point_type x_init;
    tuple_for_each(x_init, systems, [&](auto& x_elem, const auto& system_elem) {
      system_elem.get_zero_state(x_elem);
    });
    return state_belief_type(
        x_init,
        covar_type(covar_type::matrix_type(mat<double, mat_structure::diagonal>(
            total_inv_corr_dim, aCovValue))));
  }

  /**
   * Returns the dimensions of the state of the system.
   * \return The dimensions of the state of the system.
   */
  std::size_t get_state_dimensions() const { return total_state_dim; }

  /**
   * Returns the dimensions of the invariant correction of the system.
   * \return The dimensions of the invariant correction of the system.
   */
  std::size_t get_correction_dimensions() const { return total_inv_corr_dim; }

  /**
   * Returns the dimensions of the actual state of the system.
   * \return The dimensions of the actual state of the system.
   */
  std::size_t get_actual_state_dimensions() const {
    return total_actual_state_dim;
  }

  explicit ss_system_state_tuple(StateSysTuple aSystems)
      : systems(std::move(aSystems)) {
    construct_all_dimensions();
  }

  template <typename... Args>
  explicit ss_system_state_tuple(Args&&... args)
      : systems(std::forward<Args>(args)...) {
    construct_all_dimensions();
  }

  ss_system_state_tuple() : ss_system_state_tuple(StateSysTuple()) {}

  template <typename FlyWeight, typename InputType>
  void add_to_fly_weight_params(const FlyWeight& params,
                                const state_space_type& space,
                                const point_type& x, const InputType& u,
                                time_difference_type dt, time_type t) const {
    tuple_for_each(systems, [&](const auto& system_elem) {
      system_elem.add_to_fly_weight_params(params, space, x, u, dt, t);
    });
  }

  /**
   * Fills in the state-difference object with the effects of the given input on the state
   * over the given time difference.
   * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
   * \param params The fly-weight parameters that describe the dynamics of the system.
   * \param space The space object in which to operate.
   * \param x The current state object for the system.
   * \param dx The state-difference accumulated (as output) for the different effects.
   * \param u The current input vector for the system.
   * \param dt The time period over which to compute the effect of the input vector.
   * \param t The current time corresponding to the current state of the system.
   */
  template <typename FlyWeight, typename InputType>
  void add_state_difference(const FlyWeight& params,
                            const state_space_type& space, const point_type& x,
                            point_difference_type& dx, const InputType& u,
                            time_difference_type dt, time_type t = 0.0) const {
    tuple_for_each(systems, [&](const auto& system_elem) {
      system_elem.add_state_difference(params, space, x, dx, u, dt, t);
    });
  }

  /**
   * Fills in the state-difference object with the effects of the given input on the state
   * over the given time difference.
   * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
   * \param A Holds, as output, the state-to-state jacobian matrix of the state-transition of the system.
   * \param B Holds, as output, the input-to-state jacobian matrix of the state-transition of the system.
   * \param params The fly-weight parameters that describe the dynamics of the system.
   * \param space The space object in which to operate.
   * \param t_0 The previous time corresponding to the previous state of the system.
   * \param t_1 The next time corresponding to the next state of the system.
   * \param p_0 The previous state object for the system.
   * \param p_1 The next state object for the system.
   * \param u_0 The previous input vector for the system.
   * \param u_1 The next input vector for the system.
   */
  template <typename MatrixA, typename MatrixB, typename FlyWeight,
            typename InputType>
  void add_state_transition_blocks(MatrixA& A, MatrixB& B,
                                   const FlyWeight& params,
                                   const state_space_type& space, time_type t_0,
                                   time_type t_1, const point_type& p_0,
                                   const point_type& p_1, const InputType& u_0,
                                   const InputType& u_1) const {
    tuple_for_each(systems, [&](const auto& system_elem) {
      system_elem.add_state_transition_blocks(A, B, params, space, t_0, t_1,
                                              p_0, p_1, u_0, u_1);
    });
  }

  /**
   * This function computes a state corresponding to the given state corrected by a given invariant term.
   * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
   * \param params The fly-weight parameters that describe the dynamics of the system.
   * \param space The state-space within which the states reside.
   * \param x The current state of the system.
   * \param x_c The corrected state of the system, as output.
   * \param c The invariant correction term to apply to the state.
   * \param u The current input being applied to the system.
   * \param t The current time.
   * \return The corrected state of the system.
   */
  template <typename FlyWeight, typename InputType>
  void apply_correction_to_state(const FlyWeight& params,
                                 const state_space_type& space,
                                 const point_type& x, point_type& x_c,
                                 const invariant_correction_type& c,
                                 const InputType& u, const time_type& t) const {
    tuple_for_each(systems, [&](const auto& system_elem) {
      system_elem.apply_correction_to_state(params, space, x, x_c, c, u, t);
    });
  }

  /**
   * This function computes the invariant frame transition matrix for a
   * state transition from x_0 to x_1, i.e., what invariant frame transition matrix
   * describes the shift from one frame to the other.
   * \tparam FlyWeight The type of a set of records of dynamic parameters for the system.
   * \param params The fly-weight parameters that describe the dynamics of the system.
   * \param space The state-space within which the states reside.
   * \param invar_frame The matrix into which the frame transition is accumulated (it is initialized to identity).
   * \param x_0 The state of the system before the state-transition.
   * \param x_1 The state of the system after the state-transition.
   * \param u The input being applied to the system before the state-transition.
   * \param t The time before the state-transition.
   * \return The invariant frame transition matrix for the prior stage.
   */
  template <typename FlyWeight, typename InputType, typename InvarFrameType>
  void set_invariant_frame_blocks(const FlyWeight& params,
                                  const state_space_type& space,
                                  InvarFrameType& invar_frame,
                                  const point_type& x_0, const point_type& x_1,
                                  const InputType& u,
                                  const time_type& t) const {
    tuple_for_each(systems, [&](const auto& system_elem) {
      system_elem.set_invariant_frame_blocks(params, space, invar_frame, x_0,
                                             x_1, u, t);
    });
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(systems);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(systems);
    construct_all_dimensions();
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC230001A, 1, "ss_system_state_tuple",
                              named_object)
};

template <typename SystemParamPack, typename StateModelsTuple,
          typename InputModelsTuple, typename OutputModelsTuple>
class state_space_system_tuple : public named_object {
 public:
  using self = state_space_system_tuple<SystemParamPack, StateModelsTuple,
                                        InputModelsTuple, OutputModelsTuple>;

  using system_param_type = SystemParamPack;

  using state_models_type = ss_system_state_tuple<StateModelsTuple>;
  using input_models_type = ss_system_input_tuple<InputModelsTuple>;
  using output_models_type = ss_system_output_tuple<OutputModelsTuple>;

  using state_space_type = typename state_models_type::state_space_type;

  using point_type = typename pp::topology_traits<state_space_type>::point_type;
  using point_difference_type =
      typename pp::topology_traits<state_space_type>::point_difference_type;
  using point_derivative_type =
      typename pp::topology_traits<state_space_type>::point_difference_type;

  using time_type = double;
  using time_difference_type = double;

  using input_type = vect_n<double>;
  using output_type = vect_n<double>;

  using invariant_error_type = vect_n<double>;
  using invariant_correction_type = vect_n<double>;
  using invariant_frame_type = mat<double, mat_structure::square>;

  static constexpr std::size_t dimensions = 0;
  static constexpr std::size_t input_dimensions = 0;
  static constexpr std::size_t output_dimensions = 0;
  static constexpr std::size_t invariant_error_dimensions = 0;
  static constexpr std::size_t invariant_correction_dimensions = 0;
  static constexpr std::size_t actual_state_dimensions = 0;

  using matrixA_type = mat<double, mat_structure::square>;
  using matrixB_type = mat<double, mat_structure::rectangular>;
  using matrixC_type = mat<double, mat_structure::rectangular>;
  using matrixD_type = mat<double, mat_structure::rectangular>;

  struct zero_input_trajectory {
    const self* parent;
    explicit zero_input_trajectory(const self* aParent) : parent(aParent) {}
    input_type get_point(time_type /*unused*/) const {
      return input_type(parent->get_input_dimensions(), 0.0);
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

 protected:
  double dt;

  mutable system_param_type sys_params;
  state_models_type state_models;
  input_models_type input_models;
  output_models_type output_models;

 public:
  std::shared_ptr<temporal_state_space_type> get_temporal_state_space(
      double aStartTime = 0.0, double aEndTime = 1.0) const {
    return std::make_shared<temporal_state_space_type>(
        "temporal_space", state_models.create_state_space(),
        pp::time_poisson_topology("time_space", dt,
                                  (aEndTime - aStartTime) * 0.5));
  }

  std::shared_ptr<state_space_type> get_state_space() const {
    return std::make_shared<state_space_type>(
        state_models.create_state_space());
  }

  std::shared_ptr<temporal_belief_space_type> get_temporal_belief_space(
      double aStartTime = 0.0, double aEndTime = 1.0) const {
    return std::make_shared<temporal_belief_space_type>(
        "temporal_state_belief_space",
        belief_space_type(get_state_space(),
                          std::make_shared<covar_space_type>(
                              state_models.get_correction_dimensions()),
                          "state_belief_space"),
        pp::time_poisson_topology("time_space", dt,
                                  (aEndTime - aStartTime) * 0.5));
  }

  std::shared_ptr<belief_space_type> get_belief_space() const {
    return std::make_shared<belief_space_type>(
        get_state_space(),
        std::make_shared<covar_space_type>(
            state_models.get_correction_dimensions()),
        "state_belief_space");
  }

  /**
   * Returns a zero-valued state-belief for this system, with the given uniform covariance on the diagonal.
   * \param aCovValue A uniform covariance value to initialize the diagonal of the belief's covariance matrix.
   * \return A zero-valued state-belief for this system.
   */
  state_belief_type get_zero_state_belief(double aCovValue = 10.0) const {
    return state_models.get_zero_state_belief(aCovValue);
  }

  /**
   * Returns a zero-valued input-belief for this system, with the given uniform covariance on the diagonal.
   * \param aCovValue A uniform covariance value to initialize the diagonal of the belief's covariance matrix.
   * \return A zero-valued input-belief for this system.
   */
  input_belief_type get_zero_input_belief(double aCovValue = 1.0) const {
    return input_models.get_zero_input_belief(aCovValue);
  }

  /**
   * Returns a zero-valued output-belief for this system, with the given uniform covariance on the diagonal.
   * \param aCovValue A uniform covariance value to initialize the diagonal of the belief's covariance matrix.
   * \return A zero-valued output-belief for this system.
   */
  output_belief_type get_zero_output_belief(double aCovValue = 1.0) const {
    return output_models.get_zero_output_belief(aCovValue);
  }

  /**
   * Returns the dimensions of the states of the system.
   * \return The dimensions of the states of the system.
   */
  std::size_t get_state_dimensions() const {
    return state_models.get_state_dimensions();
  }

  /**
   * Returns the dimensions of the input of the system.
   * \return The dimensions of the input of the system.
   */
  std::size_t get_input_dimensions() const {
    return input_models.get_input_dimensions();
  }

  /**
   * Returns the dimensions of the output of the system.
   * \return The dimensions of the output of the system.
   */
  std::size_t get_output_dimensions() const {
    return output_models.get_output_dimensions();
  }

  /**
   * Returns the dimensions of the invariant errors of the system.
   * \return The dimensions of the invariant errors of the system.
   */
  std::size_t get_invariant_error_dimensions() const {
    return output_models.get_invariant_error_dimensions();
  }

  /**
   * Returns the dimensions of the corrections to the states of the system.
   * \return The dimensions of the corrections to the states of the system.
   */
  std::size_t get_correction_dimensions() const {
    return state_models.get_correction_dimensions();
  }

  /**
   * Returns the dimensions of the actual states of the system.
   * \return The dimensions of the actual states of the system.
   */
  std::size_t get_actual_state_dimensions() const {
    return state_models.get_actual_state_dimensions();
  }

  /**
   * Constructor.
   * \param aName The name for this object..
   * \param aDt The time-step for this discrete-time system.
   */
  explicit state_space_system_tuple(const std::string& aName, double aDt = 0.01)
      : dt(aDt) {
    setName(aName);
  }

  state_space_system_tuple() : state_space_system_tuple("") {}

  /**
   * Parametrized constructor.
   * \param aName The name for this object..
   * \param aDt The time-step for this discrete-time system.
   */
  state_space_system_tuple(
      const std::string& aName, system_param_type aSysParams,
      const StateModelsTuple& aStateModels = StateModelsTuple(),
      const InputModelsTuple& aInputModels = InputModelsTuple(),
      const OutputModelsTuple& aOutputModels = OutputModelsTuple(),
      double aDt = 0.01)
      : dt(aDt),
        sys_params(std::move(aSysParams)),
        state_models(aStateModels),
        input_models(aInputModels),
        output_models(aOutputModels) {
    setName(aName);
  }

  /**
   * This function returns a reference to the system parameter-pack used by this system.
   * \return A reference to the system parameter-pack used by this system.
   */
  system_param_type& get_system_parameters() const { return sys_params; }

  state_models_type& get_state_models() { return state_models; }
  const state_models_type& get_state_models() const { return state_models; }

  input_models_type& get_input_models() { return input_models; }
  const input_models_type& get_input_models() const { return input_models; }

  output_models_type& get_output_models() { return output_models; }
  const output_models_type& get_output_models() const { return output_models; }

  /**
   * This function returns the time-step for this discrete-time system.
   * \return The time-step for this discrete-time system.
   */
  time_difference_type get_time_step() const { return dt; }

  /**
   * This function sets the time-step for this discrete-time system.
   * \param aDt The new time-step for this discrete-time system.
   */
  void set_time_step(time_difference_type aDt) { dt = aDt; }

  /**
   * This function computes the next state of the system, i.e., the state at one time-step after the current time.
   * \param space The state-space within which the states reside.
   * \param x The current state of the system.
   * \param u The current input being applied to the system.
   * \param t The current time.
   * \return The state after one time-step beyond the given current state of the system.
   */
  point_type get_next_state(const state_space_type& space, const point_type& x,
                            const input_type& u,
                            const time_type& t = 0.0) const {
    sys_params.reset_parameters();
    state_models.add_to_fly_weight_params(*this, space, x, u, dt, t);
    sys_params.finalize_parameters();
    point_difference_type dx = space.difference(x, x);
    input_models.add_state_difference(*this, space, x, dx, u, dt, t);
    state_models.add_state_difference(*this, space, x, dx, u, dt, t);
    return space.adjust(x, dx);
  }

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
  void get_state_transition_blocks(matrixA_type& A, matrixB_type& B,
                                   const state_space_type& space,
                                   const time_type& t_0, const time_type& t_1,
                                   const point_type& p_0, const point_type& p_1,
                                   const input_type& u_0,
                                   const input_type& u_1) const {
    sys_params.reset_parameters();
    state_models.add_to_fly_weight_params(*this, space, p_0, u_0, dt, t_0);
    sys_params.finalize_parameters();
    A = mat_nil<double>(get_correction_dimensions(),
                        get_correction_dimensions());
    B = mat_nil<double>(get_correction_dimensions(), get_input_dimensions());
    input_models.add_state_transition_blocks(A, B, *this, space, t_0, t_1, p_0,
                                             p_1, u_0, u_1);
    state_models.add_state_transition_blocks(A, B, *this, space, t_0, t_1, p_0,
                                             p_1, u_0, u_1);
  }

  /**
   * This function computes the output of the system corresponding to the current state.
   * \param space The state-space within which the states reside.
   * \param x The current state of the system.
   * \param u The current input being applied to the system.
   * \param t The current time.
   * \return The output for the given current state of the system.
   */
  output_type get_output(const state_space_type& space, const point_type& x,
                         const input_type& u, const time_type& t = 0.0) const {
    sys_params.reset_parameters();
    state_models.add_to_fly_weight_params(*this, space, x, u, dt, t);
    sys_params.finalize_parameters();
    output_type y(get_output_dimensions(), 0.0);
    output_models.set_output_from_state(*this, space, x, y,
                                        t);  // TODO should have u.
    return y;
  }

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
  void get_output_function_blocks(matrixC_type& C, matrixD_type& D,
                                  const state_space_type& space,
                                  const time_type& t, const point_type& x,
                                  const input_type& u) const {
    sys_params.reset_parameters();
    state_models.add_to_fly_weight_params(*this, space, x, u, dt, t);
    sys_params.finalize_parameters();
    C = mat_nil<double>(get_invariant_error_dimensions(),
                        get_correction_dimensions());
    D = mat_nil<double>(get_invariant_error_dimensions(),
                        get_input_dimensions());
    output_models.add_output_function_blocks(C, D, *this, space, t, x, u);
  }

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
  invariant_error_type get_invariant_error(const state_space_type& space,
                                           const point_type& x,
                                           const input_type& u,
                                           const output_type& y,
                                           const time_type& t) const {
    sys_params.reset_parameters();
    state_models.add_to_fly_weight_params(*this, space, x, u, dt, t);
    sys_params.finalize_parameters();
    invariant_error_type e(get_invariant_error_dimensions(), 0.0);
    output_models.set_inv_err_from_output(*this, space, x, y, e,
                                          t);  // TODO should have u.
    return e;
  }

  /**
   * This function computes a state corresponding to the given state corrected by a given invariant term.
   * \param space The state-space within which the states reside.
   * \param x The current state of the system.
   * \param c The invariant correction term to apply to the state.
   * \param u The current input being applied to the system.
   * \param t The current time.
   * \return The corrected state of the system.
   */
  point_type apply_correction(const state_space_type& space,
                              const point_type& x,
                              const invariant_correction_type& c,
                              const input_type& u, const time_type& t) const {
    point_type result;
    state_models.apply_correction_to_state(*this, space, x, result, c, u, t);
    return result;
  }

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
  invariant_frame_type get_invariant_prior_frame(const state_space_type& space,
                                                 const point_type& x_0,
                                                 const point_type& x_1,
                                                 const input_type& u,
                                                 const time_type& t) const {
    invariant_frame_type result(mat_ident<double>(get_correction_dimensions()));
    state_models.set_invariant_frame_blocks(*this, space, result, x_0, x_1, u,
                                            t);
    return result;
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
  invariant_frame_type get_invariant_posterior_frame(
      const state_space_type& space, const point_type& x_0,
      const point_type& x_1, const input_type& u, const time_type& t) const {
    return get_invariant_prior_frame(space, x_0, x_1, u, t);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(dt) & RK_SERIAL_SAVE_WITH_NAME(sys_params) &
        RK_SERIAL_SAVE_WITH_NAME(state_models) &
        RK_SERIAL_SAVE_WITH_NAME(input_models) &
        RK_SERIAL_SAVE_WITH_NAME(output_models);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(dt) & RK_SERIAL_LOAD_WITH_NAME(sys_params) &
        RK_SERIAL_LOAD_WITH_NAME(state_models) &
        RK_SERIAL_LOAD_WITH_NAME(input_models) &
        RK_SERIAL_LOAD_WITH_NAME(output_models);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC230001B, 1, "state_space_system_tuple",
                              named_object)
};

template <typename SystemParamPack, typename StateModelsTuple,
          typename InputModelsTuple, typename OutputModelsTuple>
struct is_invariant_system<state_space_system_tuple<
    SystemParamPack, StateModelsTuple, InputModelsTuple, OutputModelsTuple>>
    : std::true_type {};

template <typename SystemParamPack, typename StateModelsTuple,
          typename InputModelsTuple, typename OutputModelsTuple>
struct is_augmented_ss_system<state_space_system_tuple<
    SystemParamPack, StateModelsTuple, InputModelsTuple, OutputModelsTuple>>
    : std::true_type {};

}  // namespace ReaK::ctrl

#endif

/**
 * \file discrete_linear_sss_concept.h
 *
 * This library defines a number of traits and concept classes related to the definition
 * of a discrete-time linear state-space system, with different system types such as a
 * linear-time-invariant (LTI), linear-time-varying (LTV), and linearized. The main
 * characteristic of such systems is that they can provide system matrices (A,B,C,D),
 * whether they are independent of state and time (LTI), independent of state (LTV),
 * or a linear approximation based on the state and time (Linearized).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2011
 */

/*
 *    Copyright 2011 Sven Mikael Persson
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

#ifndef REAK_CONTROL_SYSTEMS_DISCRETE_LINEAR_SSS_CONCEPT_H_
#define REAK_CONTROL_SYSTEMS_DISCRETE_LINEAR_SSS_CONCEPT_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/math/lin_alg/arithmetic_tuple.h"

#include "ReaK/math/lin_alg/mat_traits.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"

#include "ReaK/control/systems/discrete_sss_concept.h"
#include "ReaK/control/systems/linear_ss_system_concept.h"

#include <concepts>

namespace ReaK::ctrl {

/**
 * This traits class defines the traits of a discrete-time linear state-space system.
 * This traits class only includes the traits that are not already included in the
 * discrete_sss_traits class.
 */
template <typename T>
struct discrete_linear_sss_traits {
  using matrixA_type = typename T::matrixA_type;
  using matrixB_type = typename T::matrixB_type;
  using matrixC_type = typename T::matrixC_type;
  using matrixD_type = typename T::matrixD_type;
};

/**
 * This concept class defines no requirements for a system to be able to provide
 * linear(ized) system matrices.
 */
struct DiscreteNonLinearSystemType {
  template <typename System, typename StateSpaceType, typename Point,
            typename Input, typename Time, typename A_t, typename B_t,
            typename C_t, typename D_t>
  static void constraints(const System& /*unused*/,
                          const StateSpaceType& /*unused*/,
                          const Point& /*unused*/, const Input& /*unused*/,
                          const Time& /*unused*/, A_t& /*unused*/,
                          B_t& /*unused*/, C_t& /*unused*/, D_t& /*unused*/) {}
};

/**
 * This concept class defines the requirement for a system to be able to provide
 * linear system matrices for an LTI system, that is, the matrices are independent
 * of time or state or input.
 *
 * Valid expression:
 *
 * sys.get_state_transition_blocks(A,B);  The system's state transition matrices can be obtained without providing a
 *time or state.
 *
 * sys.get_output_function_blocks(C,D);  The system's output function matrices can be obtained without providing a time
 *or state.
 */
struct DiscreteLTISystemType {
  template <typename System, typename StateSpaceType, typename Point,
            typename Input, typename Time, typename A_t, typename B_t,
            typename C_t, typename D_t>
  static void constraints(const System& sys, const StateSpaceType& state_space,
                          const Point& /*unused*/, const Input& /*unused*/,
                          const Time& /*unused*/, A_t& A, B_t& B, C_t& C,
                          D_t& D) {
    sys.get_state_transition_blocks(A, B, state_space);
    sys.get_output_function_blocks(C, D, state_space);
  }
};

/**
 * This concept class defines the requirement for a system to be able to provide
 * linear system matrices for an LTV system, that is, the matrices are independent
 * of state or input, but dependent on time. Note that an LTI system is a subset of
 * an LTV system.
 *
 * Valid expression:
 *
 * sys.get_state_transition_blocks(A,B,t_0,t_1);  The system's state transition matrices can be obtained without
 *providing states, but providing the time before and after the step.
 *
 * sys.get_output_function_blocks(C,D,t_0);  The system's output function matrices can be obtained without providing a
 *state.
 */
struct DiscreteLTVSystemType {
  template <typename System, typename StateSpaceType, typename Point,
            typename Input, typename Time, typename A_t, typename B_t,
            typename C_t, typename D_t>
  static void constraints(const System& sys, const StateSpaceType& state_space,
                          const Point& /*unused*/, const Input& /*unused*/,
                          const Time& t, A_t& A, B_t& B, C_t& C, D_t& D) {
    sys.get_state_transition_blocks(A, B, state_space, t, t);
    sys.get_output_function_blocks(C, D, state_space, t);
  }
};

/**
 * This concept class defines the requirement for a system to be able to provide
 * linear system matrices for a linearized system, that is, the matrices are dependent
 * on time, state and-or input, but dependent on time. Note that an LTI system and LTV
 * system are subsets of a linearized system.
 *
 * Valid expression:
 *
 * sys.get_state_transition_blocks(A,B,t_0,t_1,p_0,p_1,u_0,u_1);  The system's state transition matrices can be obtained
 *by providing the time, state and input, before and after the step is taken.
 *
 * sys.get_output_function_blocks(C,D,t_0,p_0,u_0);  The system's output function matrices can be obtained from the
 *current time, state and input.
 */
struct DiscreteLinearizedSystemType {
  template <typename System, typename StateSpaceType, typename Point,
            typename Input, typename Time, typename A_t, typename B_t,
            typename C_t, typename D_t>
  static void constraints(const System& sys, const StateSpaceType& state_space,
                          const Point& p, const Input& u, const Time& t, A_t& A,
                          B_t& B, C_t& C, D_t& D) {
    sys.get_state_transition_blocks(A, B, state_space, t, t, p, p, u, u);
    sys.get_output_function_blocks(C, D, state_space, t, p, u);
  }
};

/**
 * This concept class defines the requirements for a discrete-time state-space system to be
 * considered a linear system, as used in ReaK::ctrl. This concept class depends on a helping
 * concept class which is given as SystemType which defines what kind of linear system it is
 * (DiscreteLTISystemType, DiscreteLTVSystemType or DiscreteLinearizedSystemType).
 *
 * Required concepts:
 *
 * The state-space system should model DiscreteSSSConcept.
 *
 * Valid expressions:
 *
 * sys_type.constraints(sys, state_space, p, u, t, A, B, C, D);  The system should comply to the constraints of the
 *SystemType concept class.
 *
 * p = A * p + B * u;  The next state can be obtained by linear transformation of the state and input using the system
 *matrices.
 *
 * y = C * p + D * u;  The output can be obtained by linear transformation of the state and input using the system
 *matrices.
 *
 * \tparam SystemType The concept class that tests the system-type, can be either DiscreteLTISystemType,
 *DiscreteLTVSystemType or DiscreteLinearizedSystemType.
 */
template <typename T, typename StateSpace,
          typename SystemType = DiscreteLTISystemType>
concept DiscreteLinearSSS = DiscreteSSS<T, StateSpace>&& requires(
    const T& sys, const StateSpace& space,
    typename discrete_sss_traits<T>::point_type p,
    typename discrete_sss_traits<T>::input_type u,
    typename discrete_sss_traits<T>::time_type t,
    typename discrete_sss_traits<T>::point_difference_type dp,
    typename discrete_sss_traits<T>::output_type y,
    typename discrete_linear_sss_traits<T>::matrixA_type A,
    typename discrete_linear_sss_traits<T>::matrixB_type B,
    typename discrete_linear_sss_traits<T>::matrixC_type C,
    typename discrete_linear_sss_traits<T>::matrixD_type D) {
  SystemType::constraints(sys, space, p, u, t, A, B, C, D);
  dp = ReaK::from_vect<decltype(dp)>(
      A * to_vect<mat_value_type_t<decltype(A)>>(dp) +
      B * to_vect<mat_value_type_t<decltype(B)>>(u));
  y = ReaK::from_vect<decltype(y)>(
      C * to_vect<mat_value_type_t<decltype(C)>>(dp) +
      D * to_vect<mat_value_type_t<decltype(D)>>(u));
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_SYSTEMS_DISCRETE_LINEAR_SSS_CONCEPT_H_

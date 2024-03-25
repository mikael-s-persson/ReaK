/**
 * \file belief_state_concept.h
 *
 * This library defines a number of concept and traits class templates related to
 * the definition of a belief state and the transitions between them. Belief states
 * are the foundation of estimation theory. A belief state represents the estimate
 * of a state of a system by representing the state as a probability distribution of
 * the state (e.g. a gaussian belief state would hold the mean state and the covariance
 * of that estimate).
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

#ifndef REAK_CONTROL_CONTROLLERS_BELIEF_STATE_CONCEPT_H_
#define REAK_CONTROL_CONTROLLERS_BELIEF_STATE_CONCEPT_H_

#include "ReaK/control/estimators/covariance_concept.h"
#include "ReaK/control/systems/discrete_sss_concept.h"
#include "ReaK/control/systems/state_space_sys_concept.h"
#include "ReaK/control/systems/state_vector_concept.h"
#include "ReaK/topologies/spaces/metric_space_concept.h"
#include "boost/concept_check.hpp"

namespace ReaK::ctrl {

/** This namespace includes a tag that tells whether a belief representation is unimodal or multimodal */
namespace belief_distribution {
/** This tag tells whether a belief representation is unimodal or multimodal */
enum tag { unimodal = 1, multimodal };
}  // namespace belief_distribution

/** This namespace includes a tag that tells whether a belief representation is gaussian or point-based (and possibly
 * other in the future) */
namespace belief_representation {
/** This tag tells whether a belief representation is gaussian or point-based (and possibly other in the future) */
enum tag { gaussian = 1, point_based };
}  // namespace belief_representation

/**
 * This traits class defines the traits that describe a belief-state.
 * \tparam BeliefState The belief-state type for which the traits are sought.
 */
template <typename Belief>
struct belief_state_traits {
  /** The state-vector type, should model StateVector. */
  using state_type = typename Belief::state_type;
  /** The scalar type which is arithmetically compatible with the state vector. */
  using scalar_type = typename Belief::scalar_type;

  /** This constant defines the distribution tag of the belief-state. */
  static constexpr belief_distribution::tag distribution = Belief::distribution;
  /** This constant defines the representation tag of the belief-state. */
  static constexpr belief_representation::tag representation =
      Belief::representation;
};

/**
 * This concept class checks that a belief-state models the required concept as used in ReaK::ctrl.
 *
 * Valid expressions:
 *
 * v = b.get_most_likely_state();  The most likely state can be obtained from the belief-state.
 *
 * \tparam BeliefState The belief-state type to be tested against the belief-state concept.
 */
template <typename Belief>
concept BeliefState = requires(Belief b) {
  { b.get_most_likely_state() } -> StateVector;
};

/**
 * This traits class defines the traits that describe a belief-space (topology of belief-states).
 */
template <typename BSpace>
struct belief_space_traits {
  /** The state-vector topology type, should model pp::Topology. */
  using state_topology = typename BSpace::state_topology;
};

/**
 * This concept class checks that a belief-space models the required concept as used in ReaK::ctrl.
 *
 * Required Concepts:
 *
 * The belief-space type should model the pp::Topology.
 *
 * The associated state-topology should model the pp::Topology.
 *
 * The point-type should model the BeliefState.
 *
 * Valid expressions:
 *
 * s_space = b_space.get_state_topology();  The state-topology can be obtained from the belief-space.
 */
template <typename BSpace>
concept BeliefSpace = pp::Topology<BSpace>&& requires(BSpace b_space) {
  { b_space.get_state_topology() } -> pp::Topology;
  { b_space.origin() } -> BeliefState;
};

/**
 * This class template defines the traits that a continuous belief-state should have. A
 * continuous belief state is characterized by the fact that the state variables are continuous
 * (as opposed to discrete).
 */
template <typename ContBelief>
struct continuous_belief_state_traits {

  /** The state-vector type, should model StateVectorConcept. */
  using state_type = typename ContBelief::state_type;
  /** The state-difference type, as defined by state_vector_traits. */
  using state_difference_type = typename ContBelief::state_difference_type;
  /** The size type, as defined by state_vector_traits. */
  using size_type = typename ContBelief::size_type;
  /** The scalar type which is arithmetically compatible with the state vector. */
  using scalar_type = typename ContBelief::scalar_type;

  /** The covariance matrix type which can represent the covariance of the belief-state, should model
   * CovarianceMatrixConcept. */
  using covariance_type = typename ContBelief::covariance_type;
};

/**
 * This concept class checks that a continuous belief-state models the required concept as used in ReaK::ctrl.
 *
 * Require concepts:
 *
 * state-vector type associated to the belief-state should model StateVectorConcept.
 *
 * covariance type associated to the belief-state should model CovarianceMatrixConcept.
 *
 * the belief-state should model BeliefStateConcept.
 *
 * Valid expressions (in addition to those of BeliefStateConcept):
 *
 * v = b.get_mean_state();  The mean-state can be obtained from the belief-state.
 *
 * c = b.get_covariance();  The covariance matrix can be obtained from the belief-state.
 *
 * b.set_mean_state(v);  The mean-state can be set for the belief-state.
 *
 * b.set_covariance(c);  The covariance matrix can be set for the belief-state.
 */
template <typename ContBelief>
concept ContinuousBeliefState =
    BeliefState<ContBelief>&& requires(ContBelief b) {
  { b.get_mean_state() } -> StateVector;
  {
    b.get_covariance()
    } -> CovarianceMatrix<typename continuous_belief_state_traits<
        ContBelief>::state_difference_type>;
  b.set_mean_state(b.get_mean_state());
  b.set_covariance(b.get_covariance());
};

/**
 * This traits class defines the traits that describe a continuous belief-space (topology of continuous belief-states).
 */
template <typename ContBeliefSpace>
struct continuous_belief_space_traits {
  /** The state-vector topology type, should model pp::Topology. */
  using covariance_topology = typename ContBeliefSpace::covariance_topology;
};

/**
 * This concept class checks that a belief-space models the required concept as used in ReaK::ctrl.
 *
 * Required Concepts:
 *
 * The belief-space type should also model the BeliefSpaceConcept.
 *
 * The associated covariance-topology should model the pp::TopologyConcept.
 *
 * The point-type should model the ContinuousBeliefStateConcept.
 *
 * Valid expressions:
 *
 * c_space = b_space.get_covariance_topology();  The state-topology can be obtained from the belief-space.
 */
template <typename ContBeliefSpace>
concept ContinuousBeliefSpace =
    BeliefSpace<ContBeliefSpace>&& requires(ContBeliefSpace b_space) {
  { b_space.get_covariance_topology() } -> pp::Topology;
  { b_space.origin() } -> ContinuousBeliefState;
};

/**
 * This traits class defines the traits that a belief transfer function should have. A
 * belief transfer function can be used to model the transmission of a belief-state from
 * time to time.
 */
template <typename BTransfer>
struct belief_transfer_traits {
  using state_space_system = int;
  using time_type = int;
  using time_difference_type = int;
};
template <typename BTransfer>
concept HasAllBeliefTransferTraits = requires {
  typename BTransfer::state_space_system;
  typename BTransfer::time_type;
  typename BTransfer::time_difference_type;
};
template <HasAllBeliefTransferTraits BTransfer>
struct belief_transfer_traits<BTransfer> {
  /** The state-space system type associated to the transfer function */
  using state_space_system = typename BTransfer::state_space_system;
  /** The time type associated to the transfer function */
  using time_type = typename BTransfer::time_type;
  /** The time-difference type associated to the transfer function */
  using time_difference_type = typename BTransfer::time_difference_type;
};

/**
 * This concept class defines the requirements for a type to be a belief transfer function type. A
 * belief transfer function can be used to model the transmission of a belief-state from
 * time to time.
 *
 * Required concepts:
 *
 * The associated belief state should model BeliefState.
 *
 * The associated belief-space should model the pp::Topology.
 *
 * Valid expressions:
 *
 * dt = f.get_time_step();  The time-step of the transfer function can be obtained.
 *
 * sys = *f.get_ss_system();  The state-space system can be obtained from the transfer function.
 *
 * b = f.get_next_belief(belief_space, b, t, u, y);  The next belief state (b) can be obtained from the current
 *belief-state (b), current time (t), current input (u), and next measurement (y).
 */
template <typename BTransfer, typename BSpaceType>
concept BeliefTransfer = BeliefSpace<BSpaceType>&& requires(BTransfer f) {
  {
    f.get_time_step()
    } -> std::convertible_to<
        typename belief_transfer_traits<BTransfer>::time_difference_type>;
}
&&requires(BTransfer f, BSpaceType b_space) {
  {
    *f.get_ss_system()
    } -> DiscreteSSS<std::decay_t<decltype(b_space.get_state_topology())>>;
}
&&requires(BTransfer f, BSpaceType b_space,
           typename belief_transfer_traits<BTransfer>::time_type t,
           typename ss_system_traits<typename belief_transfer_traits<
               BTransfer>::state_space_system>::input_type u,
           typename ss_system_traits<typename belief_transfer_traits<
               BTransfer>::state_space_system>::output_type y) {
  { f.get_next_belief(b_space, b_space.origin(), t, u, y) } -> BeliefState;
};

/**
 * This concept class defines the requirements for a type to be a belief predictor function type. A
 * belief predictor function can be used to model the prediction of a belief-state from
 * time to time.
 *
 * Required concepts:
 *
 * The belief predictor should also model BeliefTransfer.
 *
 * Valid expressions:
 *
 * b = f.predict_belief(b, t, u);  The next predicted belief state (b) can be obtained from the current belief-state
 *(b), current time (t), and current input (u).
 *
 * b = f.prediction_to_ML_belief(b, t, u);  The most-likely updated belief-state can be obtained from the predicted
 *belief-state (b), at time (t) with the current input (u).
 *
 * b = f.predict_ML_belief(b, t, u);  The next most-likely updated belief-state (b) can be obtained from the current
 *belief-state (b), current time (t), and current input (u).
 */
template <typename BPredictor, typename BSpaceType>
concept BeliefPredictor = BeliefTransfer<BPredictor, BSpaceType>&& requires(
    BPredictor f, BSpaceType b_space,
    typename belief_transfer_traits<BPredictor>::time_type t,
    typename ss_system_traits<typename belief_transfer_traits<
        BPredictor>::state_space_system>::input_type u) {
  { f.predict_belief(b_space, b_space.origin(), t, u) } -> BeliefState;
  { f.prediction_to_ML_belief(b_space, b_space.origin(), t, u) } -> BeliefState;
  { f.predict_ML_belief(b_space, b_space.origin(), t, u) } -> BeliefState;
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_CONTROLLERS_BELIEF_STATE_CONCEPT_H_

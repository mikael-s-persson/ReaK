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
#include "ReaK/control/systems/state_space_sys_concept.h"
#include "ReaK/control/systems/state_vector_concept.h"
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
template <typename BeliefState>
struct belief_state_traits {

  /** The state-vector type, should model StateVectorConcept. */
  using state_type = typename BeliefState::state_type;
  /** The scalar type which is arithmetically compatible with the state vector. */
  using scalar_type = typename BeliefState::scalar_type;

  /** This constant defines the distribution tag of the belief-state. */
  static constexpr belief_distribution::tag distribution =
      BeliefState::distribution;
  /** This constant defines the representation tag of the belief-state. */
  static constexpr belief_representation::tag representation =
      BeliefState::representation;
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
template <typename BeliefState>
struct BeliefStateConcept {
  BeliefState b;
  typename belief_state_traits<BeliefState>::state_type v;

  BOOST_CONCEPT_USAGE(BeliefStateConcept) { v = b.get_most_likely_state(); }
};

/**
 * This meta-function evaluates if a belief-state models the BeliefStateConcept. This meta-function
 * does not attempt to instantiate the BeliefStateConcept class because that would break Sfinae
 * functionality. Instead, if one creates a new belief-state type, one should specialize this
 * meta-function to return true.
 * \tparam BeliefState The belief-state type to be tested.
 */
template <typename BeliefState>
struct is_belief_state {
  static constexpr bool value = false;
  using type = is_belief_state<BeliefState>;
};

template <typename BeliefState>
static constexpr bool is_belief_state_v = is_belief_state<BeliefState>::value;

/**
 * This traits class defines the traits that describe a belief-space (topology of belief-states).
 * \tparam BeliefSpace The belief-space type for which the traits are sought.
 */
template <typename BeliefSpace>
struct belief_space_traits {

  /** The state-vector topology type, should model pp::TopologyConcept. */
  using state_topology = typename BeliefSpace::state_topology;
};

/**
 * This concept class checks that a belief-space models the required concept as used in ReaK::ctrl.
 *
 * Required Concepts:
 *
 * The belief-space type should model the pp::TopologyConcept.
 *
 * The associated state-topology should model the pp::TopologyConcept.
 *
 * The point-type should model the BeliefStateConcept.
 *
 * Valid expressions:
 *
 * s_space = b_space.get_state_topology();  The state-topology can be obtained from the belief-space.
 *
 * \tparam BeliefSpace The belief-space type to be tested against the belief-space concept.
 */
template <typename BeliefSpace>
struct BeliefSpaceConcept : pp::TopologyConcept<BeliefSpace> {

  using state_space_type =
      typename belief_space_traits<BeliefSpace>::state_topology;
  using belief_state_type =
      typename pp::topology_traits<BeliefSpace>::point_type;

  BOOST_CONCEPT_ASSERT((pp::TopologyConcept<state_space_type>));
  BOOST_CONCEPT_ASSERT((BeliefStateConcept<belief_state_type>));

  BeliefSpace b_space;

  BOOST_CONCEPT_USAGE(BeliefSpaceConcept) {
    const state_space_type& s_space = b_space.get_state_topology();
    RK_UNUSED(s_space);
  }
};

/**
 * This class template defines the traits that a continuous belief-state should have. A
 * continuous belief state is characterized by the fact that the state variables are continuous
 * (as opposed to discrete).
 * \tparam ContBeliefState
 */
template <typename ContBeliefState>
struct continuous_belief_state_traits {

  /** The state-vector type, should model StateVectorConcept. */
  using state_type = typename ContBeliefState::state_type;
  /** The state-difference type, as defined by state_vector_traits. */
  using state_difference_type = typename ContBeliefState::state_difference_type;
  /** The size type, as defined by state_vector_traits. */
  using size_type = typename ContBeliefState::size_type;
  /** The scalar type which is arithmetically compatible with the state vector. */
  using scalar_type = typename ContBeliefState::scalar_type;

  /** The covariance matrix type which can represent the covariance of the belief-state, should model
   * CovarianceMatrixConcept. */
  using covariance_type = typename ContBeliefState::covariance_type;
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
 *
 * \tparam ContBeliefState The continuous belief-state type to be tested against the continuous belief-state concept.
 */
template <typename ContBeliefState>
struct ContinuousBeliefStateConcept : BeliefStateConcept<ContBeliefState> {
  using state_difference_type = typename continuous_belief_state_traits<
      ContBeliefState>::state_difference_type;
  using covariance_type =
      typename continuous_belief_state_traits<ContBeliefState>::covariance_type;

  BOOST_CONCEPT_ASSERT(
      (CovarianceMatrixConcept<covariance_type, state_difference_type>));

  typename continuous_belief_state_traits<ContBeliefState>::covariance_type c;

  BOOST_CONCEPT_USAGE(ContinuousBeliefStateConcept) {
    this->v = this->b.get_mean_state();
    this->c = this->b.get_covariance();

    this->b.set_mean_state(this->v);
    this->b.set_covariance(this->c);
  }
};

/**
 * This meta-function evaluates if a belief-state models the ContinuousBeliefStateConcept. This meta-function
 * does not attempt to instantiate the ContinuousBeliefStateConcept class because that would break Sfinae
 * functionality. Instead, if one creates a new continuous belief-state type, one should specialize this
 * meta-function to return true.
 * \tparam BeliefState The continuous belief-state type to be tested.
 */
template <typename BeliefState>
struct is_continuous_belief_state {
  static constexpr bool value = false;
  using type = is_continuous_belief_state<BeliefState>;
};

template <typename BeliefState>
static constexpr bool is_continuous_belief_state_v =
    is_continuous_belief_state<BeliefState>::value;

/**
 * This traits class defines the traits that describe a continuous belief-space (topology of continuous belief-states).
 * \tparam ContBeliefSpace The continuous belief-space type for which the traits are sought.
 */
template <typename ContBeliefSpace>
struct continuous_belief_space_traits {

  /** The state-vector topology type, should model pp::TopologyConcept. */
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
 *
 * \tparam BeliefSpace The belief-space type to be tested against the belief-space concept.
 */
template <typename ContBeliefSpace>
struct ContinuousBeliefSpaceConcept : BeliefSpaceConcept<ContBeliefSpace> {

  using covariance_space_type = typename continuous_belief_space_traits<
      ContBeliefSpace>::covariance_topology;
  using belief_state_type =
      typename pp::topology_traits<ContBeliefSpace>::point_type;

  BOOST_CONCEPT_ASSERT((pp::TopologyConcept<covariance_space_type>));
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<belief_state_type>));

  BOOST_CONCEPT_USAGE(ContinuousBeliefSpaceConcept) {
    const covariance_space_type& c_space =
        this->b_space.get_covariance_topology();
    RK_UNUSED(c_space);
  }
};

/**
 * This traits class defines the traits that a belief transfer function should have. A
 * belief transfer function can be used to model the transmission of a belief-state from
 * time to time.
 * \tparam BeliefTransfer The belief transfer function for which the traits are sought.
 */
template <typename BeliefTransfer>
struct belief_transfer_traits {
  /** The state-space system type associated to the transfer function */
  using state_space_system = typename BeliefTransfer::state_space_system;
  /** The time type associated to the transfer function */
  using time_type = typename BeliefTransfer::time_type;
  /** The time-difference type associated to the transfer function */
  using time_difference_type = typename BeliefTransfer::time_difference_type;
};

/**
 * This concept class defines the requirements for a type to be a belief transfer function type. A
 * belief transfer function can be used to model the transmission of a belief-state from
 * time to time.
 *
 * Required concepts:
 *
 * The associated belief state should model BeliefStateConcept.
 *
 * The associated belief-space should model the pp::TopologyConcept.
 *
 * Valid expressions:
 *
 * dt = f.get_time_step();  The time-step of the transfer function can be obtained.
 *
 * sys = *f.get_ss_system();  The state-space system can be obtained from the transfer function.
 *
 * b = f.get_next_belief(belief_space, b, t, u, y);  The next belief state (b) can be obtained from the current
 *belief-state (b), current time (t), current input (u), and next measurement (y).
 *
 * \tparam BeliefTransfer The belief transfer function type to be checked.
 * \tparam BeliefSpaceType The belief topology type on which the beliefs lie.
 */
template <typename BeliefTransfer, typename BeliefSpaceType>
struct BeliefTransferConcept {
  using BeliefState = typename pp::topology_traits<BeliefSpaceType>::point_type;
  using StateSpaceSystem =
      typename belief_transfer_traits<BeliefTransfer>::state_space_system;
  using TimeType = typename belief_transfer_traits<BeliefTransfer>::time_type;
  using TimeDiffType =
      typename belief_transfer_traits<BeliefTransfer>::time_difference_type;

  BOOST_CONCEPT_ASSERT((BeliefStateConcept<BeliefState>));
  BOOST_CONCEPT_ASSERT((BeliefSpaceConcept<BeliefSpaceType>));

  BeliefTransfer f;
  typename ss_system_traits<StateSpaceSystem>::input_type u;
  typename ss_system_traits<StateSpaceSystem>::output_type y;
  BeliefState b;
  StateSpaceSystem sys;
  BeliefSpaceType belief_space;
  TimeType t;
  TimeDiffType dt;

  BOOST_CONCEPT_USAGE(BeliefTransferConcept) {
    dt = f.get_time_step();
    sys = *f.get_ss_system();
    b = f.get_next_belief(belief_space, b, t, u, y);
  }
};

/**
 * This meta-function evaluates if a belief transfer models the BeliefTransferConcept. This meta-function
 * does not attempt to instantiate the BeliefTransferConcept class because that would break Sfinae
 * functionality. Instead, if one creates a new belief transfer type, one should specialize this
 * meta-function to return true.
 * \tparam BeliefTransfer The belief transfer function type to be tested.
 */
template <typename BeliefTransfer>
struct is_belief_transfer {
  static constexpr bool value = false;
  using type = is_belief_transfer<BeliefTransfer>;
};

template <typename BeliefTransfer>
static constexpr bool is_belief_transfer_v =
    is_belief_transfer<BeliefTransfer>::value;

/**
 * This concept class defines the requirements for a type to be a belief predictor function type. A
 * belief predictor function can be used to model the prediction of a belief-state from
 * time to time.
 *
 * Required concepts:
 *
 * The associated belief state should model BeliefStateConcept.
 *
 * The associated state-space system should model SSSystemConcept.
 *
 * The belief predictor should also model BeliefTransferConcept.
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
 *
 * \tparam BeliefPredictor The belief predictor function type to be checked.
 * \tparam BeliefSpaceType The belief topology type on which the beliefs lie.
 */
template <typename BeliefPredictor, typename BeliefSpaceType>
struct BeliefPredictorConcept
    : BeliefTransferConcept<BeliefPredictor, BeliefSpaceType> {

  BOOST_CONCEPT_USAGE(BeliefPredictorConcept) {
    this->b =
        this->f.predict_belief(this->belief_space, this->b, this->t, this->u);
    this->b = this->f.prediction_to_ML_belief(this->belief_space, this->b,
                                              this->t, this->u);
    this->b = this->f.predict_ML_belief(this->belief_space, this->b, this->t,
                                        this->u);
  }
};

/**
 * This meta-function evaluates if a belief predictor models the BeliefPredictorConcept. This meta-function
 * does not attempt to instantiate the BeliefPredictorConcept class because that would break Sfinae
 * functionality. Instead, if one creates a new belief predictor type, one should specialize this
 * meta-function to return true.
 * \tparam BeliefTransfer The belief predictor function type to be tested.
 */
template <typename BeliefPredictor>
struct is_belief_predictor {
  static constexpr bool value = false;
  using type = is_belief_predictor<BeliefPredictor>;
};

template <typename BeliefPredictor>
static constexpr bool is_belief_predictor_v =
    is_belief_predictor<BeliefPredictor>::value;

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_CONTROLLERS_BELIEF_STATE_CONCEPT_H_

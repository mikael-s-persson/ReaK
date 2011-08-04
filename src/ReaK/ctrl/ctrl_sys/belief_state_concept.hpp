
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

#ifndef BELIEF_STATE_CONCEPT_HPP
#define BELIEF_STATE_CONCEPT_HPP

#include <boost/concept_check.hpp>
#include "state_vector_concept.hpp"
#include "covariance_concept.hpp"
#include "state_space_sys_concept.hpp"

namespace ReaK {

namespace ctrl {


namespace belief_distribution {
  enum tag {
    unimodal = 1,
    multimodal
  };
};

namespace belief_representation {
  enum tag {
    gaussian = 1,
    point_based
  };
};



template <typename BeliefState>
struct belief_state_traits {
  
  typedef typename BeliefState::state_type state_type;
  typedef typename BeliefState::scalar_type scalar_type;
  typedef typename BeliefState::random_sampler_type random_sampler_type;
  typedef typename BeliefState::pdf_type pdf_type;
  
  BOOST_STATIC_CONSTANT(belief_distribution::tag, distribution = BeliefState::distribution);
  BOOST_STATIC_CONSTANT(belief_representation::tag, representation = BeliefState::representation);
  
};



template <typename BeliefState>
struct BeliefStateConcept {
  BeliefState b;
  typename belief_state_traits<BeliefState>::state_type v;
  typename belief_state_traits<BeliefState>::scalar_type s;
  typename belief_state_traits<BeliefState>::pdf_type p;
  typename belief_state_traits<BeliefState>::random_sampler_type rs;
  
  void constraints() {
    boost::function_requires< StateVectorConcept< typename belief_state_traits<BeliefState>::state_type > >();
    
    p = b.get_pdf();
    s = p(v);
    v = b.get_most_likely_state();
    rs = b.get_random_sampler();
    v = rs();
  };
  
};


template <typename BeliefState>
struct is_belief_state {
  BOOST_STATIC_CONSTANT(bool, value = false);
  typedef is_belief_state<BeliefState> type;
};


template <typename ContBeliefState>
struct continuous_belief_state_traits {
  
  typedef typename ContBeliefState::state_type state_type;
  typedef typename ContBeliefState::state_difference_type state_difference_type;
  typedef typename ContBeliefState::size_type size_type;
  typedef typename ContBeliefState::scalar_type scalar_type;
  typedef typename ContBeliefState::random_sampler_type random_sampler_type;
  
  typedef typename ContBeliefState::covariance_type covariance_type;
  
};


template <typename ContBeliefState>
struct ContinuousBeliefStateConcept {
  ContBeliefState b;
  typename continuous_belief_state_traits<ContBeliefState>::state_type v;
  typename continuous_belief_state_traits<ContBeliefState>::scalar_type s;
  typename continuous_belief_state_traits<ContBeliefState>::covariance_type c;
  
  void constraints() {
    boost::function_requires< BeliefStateConcept<ContBeliefState> >();
    boost::function_requires< CovarianceMatrixConcept< typename continuous_belief_state_traits<ContBeliefState>::covariance_type > >();
    v = b.get_mean_state();
    c = b.get_covariance();
    
    s = diff(v,v) * ( c.get_inverse_matrix() * diff(v,v) );
    
    b.set_mean_state(v);
    b.set_covariance(c);
  };
  
};


template <typename BeliefState>
struct is_continuous_belief_state {
  BOOST_STATIC_CONSTANT(bool, value = false);
  typedef is_continuous_belief_state<BeliefState> type;
};





template <typename BeliefTransfer>
struct belief_transfer_traits {
  typedef typename BeliefTransfer::belief_state belief_state;
  typedef typename BeliefTransfer::state_space_system state_space_system;
  typedef typename BeliefTransfer::time_type time_type;
  typedef typename BeliefTransfer::time_difference_type time_difference_type;
};


template <typename BeliefTransfer>
struct BeliefTransferConcept {
  BeliefTransfer f;
  typedef typename belief_transfer_traits< BeliefTransfer >::belief_state BeliefState;
  typedef typename belief_transfer_traits< BeliefTransfer >::state_space_system StateSpaceSystem;
  typedef typename belief_transfer_traits< BeliefTransfer >::time_type TimeType;
  typedef typename belief_transfer_traits< BeliefTransfer >::time_difference_type TimeDiffType;
  typename ss_system_traits<SSSystem>::input_type u;
  typename ss_system_traits<SSSystem>::output_type y;
  BeliefState b;
  StateSpaceSystem sys;
  TimeType t;
  TimeDiffType dt;

  void constraints() {
    boost::function_requires< BeliefStateConcept< BeliefState > >();
    boost::function_requires< SSSystemConcept< StateSpaceSystem > >();

    dt = f.get_time_step();
    sys = f.get_ss_system();
    b = f.get_next_belief(b, t, u, y);
  };
};

template <typename BeliefTransfer>
struct is_belief_transfer {
  BOOST_STATIC_CONSTANT(bool, value = false);
  typedef is_belief_transfer<BeliefTransfer> type;
};







template <typename BeliefPredictor>
struct BeliefPredictorConcept {
  BeliefPredictor f;
  typedef typename belief_transfer_traits< BeliefPredictor >::belief_state BeliefState;
  typedef typename belief_transfer_traits< BeliefPredictor >::state_space_system StateSpaceSystem;
  typedef typename belief_transfer_traits< BeliefPredictor >::time_type TimeType;
  typedef typename belief_transfer_traits< BeliefPredictor >::time_difference_type TimeDiffType;
  typename ss_system_traits<StateSpaceSystem>::input_type u;
  typename ss_system_traits<StateSpaceSystem>::output_type y;
  BeliefState b;
  StateSpaceSystem sys;
  TimeType t;
  TimeDiffType dt;

  void constraints() {
    boost::function_requires< BeliefStateConcept< BeliefState > >();
    boost::function_requires< SSSystemConcept< StateSpaceSystem > >();
    boost::function_requires< BeliefTransferConcept< BeliefPredictor > >();

    b = f.predict_belief(b, t, u);
    b = f.prediction_to_ML_belief(b, t, u);
    b = f.predict_ML_belief(b, t, u);
  };
};

template <typename BeliefPredictor>
struct is_belief_predictor {
  BOOST_STATIC_CONSTANT(bool, value = false);
  typedef is_belief_predictor<BeliefPredictor> type;
};




};

};

#endif













